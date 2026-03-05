"""
compute.py — thin wrapper that calls the existing flow_computations pipeline
and returns JSON-serialisable results.

This module is imported by app.py, which already inserts the flow_computations/
directory into sys.path before importing this module.

Caching
-------
The expensive part of each run is buildmatrix + inv(pkernel).  Because the
solution is a linear Stokes system, pcoeff ∝ amu (viscosity).  So if only
viscosity, plot depth, or grid resolution changes we can reuse the cached
pcoeff (scaled by amu_new/amu_base) and call only outputgrids — skipping the
~13 s matrix inversion entirely.

Cache key = (model, basal_bc, flux_slab, flux_width, flux_alpha,
             no_flux_for_slabtails,
             euler_lat, euler_lon, euler_rate, trench_vew, trench_vns)
"""
import sys
import time
import math
import pathlib
import numpy as np
from numpy.linalg import inv

# Ensure flow_computations is importable even when this module is used standalone.
_FC_DIR = pathlib.Path(__file__).parent.parent / 'flow_computations'
if str(_FC_DIR) not in sys.path:
    sys.path.insert(0, str(_FC_DIR))

from functions import (
    readdomains, readgrid, readbounds, organizebounds,
    pressurepoints, buildmatrix, buildvector, outputgrids, outputDP,
    partition_polygon_points,
)
from euler_pole import EulerPole, sph2cart as _ep_sph2cart, local_coords as _ep_local_coords

# ── Fixed physical constants (match global_pressure_withPressurePlot.py) ──────
_RAD_KM      = 6378.0
_ALITH       = 80.0e3        # lithospheric thickness, m
_AH1         = 660.0e3       # total thickness, m
_SHIFT_EDGES = 0
_EPSLRC      = 100.0e3
_EPSLIT      = 1.0e3
_EPS_FACT    = 0.01
_EPSDP_FACT  = 2 * _EPS_FACT
_ISEG_MIN    = 0
_FLUX_VEL_CONST = 50.0       # mm/yr, only used when flux_slab == 1
_PRESS_DEPTH = 330.0e3       # depth at which pressure field is evaluated, m (default)
_DIP_DEPTH   = 330.0e3       # depth for DP calculation, m

# ── Module-level cache ────────────────────────────────────────────────────────
_cache     = None   # set after first successful solve
_poly_cache = {}    # (id(bound_ind), spacing) → (polygons, polygon_points)

# ── Earth / real-plate models (velocities come entirely from input files) ─────
_EARTH_MODELS = {'Slab2.0Final_NoJapTail_nnr_FS'}

# ── Simple parametric model ───────────────────────────────────────────────────
_SIMPLE_MODEL = 'Simple'


def generate_simple_inp_files(width_deg: float, length_deg: float) -> None:
    """
    Write Subbon_Simple.inp and Subfil_Simple.inp for a rectangular SP geometry.

    The subducting plate is centred at (lon=180°, lat=0°) with:
      - E–W extent = width_deg  → trench at 180-W/2, far edge at 180+W/2
      - N–S extent = length_deg → north at +L/2, south at -L/2
    The plate moves eastward at ~50 mm/yr (Euler pole at 90°N, 180°E, 0.45°/Myr).
    Segment count: ~1 per 10° on each side.
    Segment order: right-edge (S→N) → top (E→W) → trench (N→S) → bottom (W→E).
    """
    # width_deg  = slab length (along-trench, N–S)
    # length_deg = trench-perpendicular extent (E–W)
    lon_t = 180.0 - length_deg / 2.0
    lon_f = 180.0 + length_deg / 2.0
    lat_n =         width_deg  / 2.0
    lat_s =        -width_deg  / 2.0

    n_ns = max(2, math.ceil(width_deg  / 2.5))  # segments along N–S edges
    n_ew = max(2, math.ceil(length_deg / 2.5))  # segments along E–W edges

    def fmt(*vals):
        return '\t'.join(str(v) for v in vals)

    lines = []
    seg = 1

    # Right edge: n_ns segments south → north at lon_f  (N–S = width_deg)
    for k in range(n_ns):
        la = lat_s + k       * width_deg / n_ns
        lb = lat_s + (k + 1) * width_deg / n_ns
        lines.append(fmt(0, lon_f, la, lon_f, lb, 2, 1, f'bound {seg}'))
        seg += 1

    # Top edge: n_ew segments east → west at lat_n  (E–W = length_deg)
    for k in range(n_ew):
        loa = lon_f - k       * length_deg / n_ew
        lob = lon_f - (k + 1) * length_deg / n_ew
        lines.append(fmt(0, loa, lat_n, lob, lat_n, 2, 1, f'bound {seg}'))
        seg += 1

    # Trench: n_ns segments north → south at lon_t (type 1)  (N–S = width_deg)
    for k in range(n_ns):
        la = lat_n - k       * width_deg / n_ns
        lb = lat_n - (k + 1) * width_deg / n_ns
        lines.append(fmt(1, lon_t, la, lon_t, lb, 2, 1, 50., 0., 'r', 1, f'bound {seg}'))
        seg += 1

    # Bottom edge: n_ew segments west → east at lat_s  (E–W = length_deg)
    for k in range(n_ew):
        loa = lon_t + k       * length_deg / n_ew
        lob = lon_t + (k + 1) * length_deg / n_ew
        lines.append(fmt(0, loa, lat_s, lob, lat_s, 2, 1, f'bound {seg}'))
        seg += 1

    total = seg - 1
    bound_list = ','.join(str(i) for i in range(1, total + 1))

    subfil_lines = [
        '100\t0.\t90.\t0.0\t-\t-\t-\t*1\t0',
        f'100\t180.\t90.\t0.45\t-\t-\t-\t*2\t{bound_list}',
    ]

    inp_dir = _FC_DIR / 'inputs'
    (inp_dir / f'Subbon_{_SIMPLE_MODEL}.inp').write_text('\n'.join(lines) + '\n')
    (inp_dir / f'Subfil_{_SIMPLE_MODEL}.inp').write_text('\n'.join(subfil_lines) + '\n')


def _to_float(v):
    """Convert numpy scalar to plain Python float for JSON serialisation."""
    return float(v)


def _seg_midpoint(lona_i, lata_i, lonb_i, latb_i):
    """Midpoint lon (0-360) / lat, handling the dateline correctly."""
    dlon = lonb_i - lona_i
    if dlon > 180:  dlon -= 360
    if dlon < -180: dlon += 360
    return (lona_i + dlon / 2) % 360, (lata_i + latb_i) / 2


def get_geometry(model: str,
                 width_deg: float = None, length_deg: float = None) -> dict:
    """
    Read boundary and domain files for a model.
    Returns boundaries, plate velocity vectors, and trench velocity vectors
    at segment midpoints — all without running the matrix solve.
    For model='Simple', generates the input files from width_deg / length_deg first.
    """
    if model == _SIMPLE_MODEL:
        generate_simple_inp_files(
            width_deg  if width_deg  is not None else 60.0,
            length_deg if length_deg is not None else 45.0,
        )
    inp_dir   = _FC_DIR / 'inputs'
    f_bounds  = str(inp_dir / f'Subbon_{model}.inp')
    f_domains = str(inp_dir / f'Subfil_{model}.inp')

    (ndomain,
     pole_top_lon, pole_top_lat, pole_top_rate,
     pole_bott_lon, pole_bott_lat, pole_bott_rate,
     rigid_vew, rigid_vns,
     domain_bounds) = readdomains(f_domains)

    (num_bounds,
     iwall, lona, lata, lonb, latb, bound_ind,
     idl, idr, vt_ew, vt_ns, polarity, large_wall_inds) = readbounds(f_bounds)

    boundaries = []
    trench_velocities = []

    for i in range(num_bounds):
        boundaries.append({
            'type':         int(iwall[i]),
            'lona':         _to_float(lona[i]),
            'lata':         _to_float(lata[i]),
            'lonb':         _to_float(lonb[i]),
            'latb':         _to_float(latb[i]),
            'left_domain':  int(idl[i]),
            'right_domain': int(idr[i]),
            'polarity':     int(polarity[i]),
        })

        if iwall[i] == 1:
            mid_lon, mid_lat = _seg_midpoint(lona[i], lata[i], lonb[i], latb[i])
            trench_velocities.append({
                'lon': _to_float(mid_lon),
                'lat': _to_float(mid_lat),
                'vew': _to_float(vt_ew[i]),  # mm/yr
                'vns': _to_float(vt_ns[i]),  # mm/yr
            })

    # ── Plate velocities on an interior grid ────────────────────────────────
    # Use partition_polygon_points to assign each coarse grid point to a domain,
    # then compute the Euler-pole velocity there.
    GRID_DEG = 15.0
    lons_c = np.arange(0.0, 360.0, GRID_DEG)
    lats_c = np.arange(-82.5, 83.0, GRID_DEG)

    _, polygon_points = partition_polygon_points(
        lons_c, lats_c,
        list(range(1, num_bounds + 1)),   # bound_ind: 1-based sequential
        lona, lata, lonb, latb,
        domain_bounds, _RAD_KM,
    )

    plate_velocities = []
    for ii in range(len(lats_c)):
        for jj in range(len(lons_c)):
            domain = int(polygon_points[ii, jj])
            if domain == 0:
                continue
            dom = domain - 1
            lat, lon = float(lats_c[ii]), float(lons_c[jj])
            if ndomain[dom] == 500:
                vew = float(rigid_vew[dom])
                vns = float(rigid_vns[dom])
            else:
                pole = EulerPole(lat=pole_top_lat[dom],
                                 lon=pole_top_lon[dom],
                                 rate=pole_top_rate[dom])
                vew, vns, _ = pole.velocity_components(lat, lon)
                vew, vns = float(vew), float(vns)
            plate_velocities.append({
                'lon':    lon,
                'lat':    lat,
                'vew':    vew,    # mm/yr
                'vns':    vns,    # mm/yr
                'domain': domain, # 1-indexed
            })

    is_earth = model in _EARTH_MODELS
    sp_euler_pole     = None
    trench_vel_default = None
    sp_domain_idx     = None

    if not is_earth:
        for i in range(num_bounds):
            if iwall[i] == 1:
                sp_domain_idx = int(idl[i])   # 1-indexed
                sp_dom = sp_domain_idx - 1
                sp_euler_pole = {
                    'lat':  float(pole_top_lat[sp_dom]),
                    'lon':  float(pole_top_lon[sp_dom]),
                    'rate': float(pole_top_rate[sp_dom]),
                }
                break
        slab_vews = [float(vt_ew[i]) for i in range(num_bounds) if iwall[i] == 1]
        slab_vnss = [float(vt_ns[i]) for i in range(num_bounds) if iwall[i] == 1]
        if slab_vews:
            trench_vel_default = {
                'vew': sum(slab_vews) / len(slab_vews),
                'vns': sum(slab_vnss) / len(slab_vnss),
            }

    return {
        'boundaries':        boundaries,
        'plate_velocities':  plate_velocities,
        'trench_velocities': trench_velocities,
        'is_earth':          is_earth,
        'sp_euler_pole':     sp_euler_pole,
        'trench_vel_default': trench_vel_default,
        'sp_domain_idx':     sp_domain_idx,
    }


def _safe(v):
    """Return None for non-finite floats so JSON stays valid (no bare Infinity/NaN)."""
    f = float(v)
    return f if math.isfinite(f) else None


def _serialise_outputs(P_out, vel_ew, vel_ns, DP, iwall,
                       lons_grd, lats_grd, orig_bounds, t0):
    """Pack outputgrids / outputDP results into a JSON-serialisable dict."""
    lons_1d = lons_grd[0, :].tolist()
    lats_1d = lats_grd[:, 0].tolist()

    dp_list = [
        {
            'lon':    _to_float(DP[i, 0]),
            'lat':    _to_float(DP[i, 1]),
            'dp_MPa': _safe(DP[i, 4]),
        }
        for i in range(len(DP))
        if iwall[i] == 1
    ]

    step_lon = max(1, len(lons_1d) // 20)
    step_lat = max(1, len(lats_1d) // 10)
    vel_list = [
        {
            'lon': lons_1d[j],
            'lat': lats_1d[i],
            'vew': _safe(vel_ew[i, j]),
            'vns': _safe(vel_ns[i, j]),
        }
        for i in range(0, len(lats_1d), step_lat)
        for j in range(0, len(lons_1d), step_lon)
        if math.isfinite(vel_ew[i, j]) and math.isfinite(vel_ns[i, j])
    ]

    vel_max_cmyr = float(np.nanmax(np.hypot(vel_ew, vel_ns))) if vel_ew.size else 0.0
    if not math.isfinite(vel_max_cmyr):
        vel_max_cmyr = 0.0

    return {
        'lons':         lons_1d,
        'lats':         lats_1d,
        'pressure':     [[_safe(p) for p in row] for row in P_out.tolist()],
        'dp':           dp_list,
        'velocity':     vel_list,
        'boundaries':   orig_bounds,
        'vel_max_cmyr': vel_max_cmyr,
        'timing_s':     time.perf_counter() - t0,
    }


_EP_EARTH_RADIUS_KM = 6371.009  # must match EulerPole.earth_radius


def _euler_vel_arr(pole_lat, pole_lon, pole_rate, lats, lons):
    """
    Vectorized Euler-pole velocity.  lats/lons are 1-D arrays.
    Returns (east_cmyr, north_cmyr) arrays (mm/yr → cm/yr via ×0.1).
    """
    ox, oy, oz = _ep_sph2cart(pole_lat, pole_lon, np.radians(-pole_rate))
    omega = np.array([ox, oy, oz])                      # (3,)
    x, y, z = _ep_sph2cart(lats, lons, _EP_EARTH_RADIUS_KM)
    r = np.column_stack([x, y, z])                      # (N, 3)
    v = np.cross(omega, r)                              # (N, 3)
    east, north, _ = _ep_local_coords(lats, lons, v[:, 0], v[:, 1], v[:, 2])
    return east * 0.1, north * 0.1                      # cm/yr


def _fast_outputgrids(spacing, lona, lata, lonb, latb, lono, lato,
                      iwall, gam, alpha, amu, a_asth,
                      n_segs, num_segs, pcoeff,
                      rad_km, domain_bounds, bound_ind,
                      pole_top_lon, pole_top_lat, pole_top_rate,
                      vt_ew, vt_ns, alith, press_depth,
                      coefftr1, coefftr2,
                      pole_bott_lon, pole_bott_lat, pole_bott_rate,
                      rigid_vew, rigid_vns, ah1, ndomain):
    """
    Fast replacement for functions.outputgrids using NumPy vectorisation.

    Instead of calling scalar findpressure_edge / findpressure_wall once per
    (lat, lon, segment) triple, this function iterates over *segments* and
    evaluates each Green's-function kernel on the entire observation grid at
    once.  This trades O(nlat × nlon × n_segs) Python function calls for
    O(n_segs) NumPy ufunc calls, giving a 20-100× speedup on the grid phase.

    Returns the same 18-element tuple as outputgrids().
    """
    # ── Grid setup ────────────────────────────────────────────────────────────
    nlat = int(round(180.0 / spacing)) + 1
    nlon = int(round(360.0 / spacing)) + 1
    lats = np.linspace(-90.0, 90.0, nlat)
    lons = np.linspace(0.0, 360.0, nlon)
    # lat_grd[i, j] = lats[i],  lon_grd[i, j] = lons[j]  (matches original)
    lat_grd, lon_grd = np.meshgrid(lats, lons, indexing='ij')

    rad          = rad_km * 1.0e3
    rad_pressloc = rad_km - press_depth / 1.0e3

    # ── Pre-convert observation grid to radians (reused across all segments) ──
    lonobs_r = np.deg2rad(lon_grd)   # (nlat, nlon)
    latobs_r = np.deg2rad(lat_grd)   # (nlat, nlon)

    # ── Segment midpoints (vectorised over all segments at once) ───────────────
    lona_arr = np.asarray(lona, dtype=float)
    lata_arr = np.asarray(lata, dtype=float)
    lonb_arr = np.asarray(lonb, dtype=float)
    latb_arr = np.asarray(latb, dtype=float)
    gam_arr  = np.asarray(gam,  dtype=float)
    alp_arr  = np.asarray(alpha, dtype=float)

    lonA_r = np.deg2rad(lona_arr); latA_r = np.deg2rad(lata_arr)
    lonB_r = np.deg2rad(lonb_arr); latB_r = np.deg2rad(latb_arr)
    dLon   = lonB_r - lonA_r
    Bx     = np.cos(latB_r) * np.cos(dLon)
    By     = np.cos(latB_r) * np.sin(dLon)
    latC_r = np.arctan2(
        np.sin(latA_r) + np.sin(latB_r),
        np.sqrt((np.cos(latA_r) + Bx) ** 2 + By ** 2),
    )
    lonC_r = lonA_r + np.arctan2(By, np.cos(latA_r) + Bx)
    lonC_r = (lonC_r + 3.0 * np.pi) % (2.0 * np.pi) - np.pi   # → [-π, π]
    lonmid_arr = np.rad2deg(lonC_r)   # [-180, 180], same as midpoint()
    latmid_arr = np.rad2deg(latC_r)

    # Haversine helper: scalar source point → 2-D grid distance (km)
    def _hav(lon_s_r, lat_s_r):
        dlo = lonobs_r - lon_s_r
        dla = latobs_r - lat_s_r
        aa  = (np.sin(dla / 2) ** 2
               + np.cos(lat_s_r) * np.cos(latobs_r) * np.sin(dlo / 2) ** 2)
        return 2.0 * np.arcsin(np.sqrt(np.clip(aa, 0.0, 1.0))) * rad_km

    # ── Pressure accumulation: loop over segments, NumPy over grid ────────────
    P_grd     = np.zeros((nlat, nlon))
    Pwall_grd = np.zeros((nlat, nlon))
    Pedge_grd = np.zeros((nlat, nlon))

    for iset in range(n_segs):
        gm       = gam_arr[iset]
        pc       = float(pcoeff[iset])
        lonmid_r = lonC_r[iset];    latmid_r = latC_r[iset]
        lonmid   = lonmid_arr[iset]; latmid  = latmid_arr[iset]
        lonaa_r  = lonA_r[iset];    lataa_r  = latA_r[iset]
        lonbb_r  = lonB_r[iset];    latbb_r  = latB_r[iset]

        if iset >= num_segs:
            # ── Slab-wall kernel (replicates findpressure_wall) ───────────────
            alp    = (alp_arr[iset] + 2.0 * np.pi) % (2.0 * np.pi)
            A_pt_w = gm / (4.0 * rad)

            angle  = _hav(lonmid_r, latmid_r) / rad_km   # radians

            # azimuth transform — matches the scalar formula exactly
            lon_rad  = (np.deg2rad(lon_grd - lonmid - 180.0) + 2.0 * np.pi) % (2.0 * np.pi)
            dlat_seg = np.deg2rad(90.0 - latmid)
            lons_tr  = np.arctan2(
                np.sin(lon_rad),
                np.tan(latobs_r) * np.sin(dlat_seg) + np.cos(lon_rad) * np.cos(dlat_seg),
            )
            pt_azim = (lons_tr + alp + 2.0 * np.pi) % (2.0 * np.pi)

            dista = _hav(lonaa_r, lataa_r) * 1.0e3   # m
            distb = _hav(lonbb_r, latbb_r) * 1.0e3   # m

            coshlam = 0.5 / gm * (dista + distb)
            sinhlam = np.sqrt(np.maximum(coshlam ** 2 - 1.0, 0.0))
            cossig  = np.clip(0.5 / gm * (dista - distb), -1.0, 1.0)
            sinsig  = np.sqrt(np.maximum(1.0 - cossig ** 2, 0.0))
            sinsig  = np.where(pt_azim <= np.pi, -sinsig, sinsig)

            t       = coshlam - sinhlam
            P_plane = (gm * t * sinsig
                       - (1.0 / 3.0) * gm * t ** 3
                       * (3.0 * sinsig * cossig ** 2 - sinsig ** 3))

            cos_a    = np.where(np.cos(angle) >= 1.0, np.cos(1e-10), np.cos(angle))
            P_sphere = A_pt_w * gm * np.sin(-pt_azim) * np.sin(angle) / (1.0 - cos_a)

            angle_s = np.where(angle == 0.0, 1e-10, angle)
            P_xy    = A_pt_w * gm * 2.0 / angle_s * np.sin(-pt_azim)

            contrib    = (P_plane - P_xy + P_sphere) * pc
            P_grd     += contrib
            Pwall_grd += contrib

        else:
            # ── Plate-edge kernel (replicates findpressure_edge) ──────────────
            A_pt_e   = -gm / np.pi
            shift_pt = (-2.0 * gm / np.pi) * (1.0 + np.log(rad * np.sqrt(2.0) / gm))
            Ppt_avg  = A_pt_e * (np.log(2.0) - 1.0)

            angle = _hav(lonmid_r, latmid_r) / rad_km   # radians

            dista = _hav(lonaa_r, lataa_r) * 1.0e3
            distb = _hav(lonbb_r, latbb_r) * 1.0e3

            coshlam = 0.5 / gm * (dista + distb)
            sinhlam = np.sqrt(np.maximum(coshlam ** 2 - 1.0, 0.0))
            cossig  = np.clip(0.5 / gm * (dista - distb), -1.0, 1.0)
            sinsig  = np.sqrt(np.maximum(1.0 - cossig ** 2, 0.0))

            xf = gm * coshlam * cossig
            yf = gm * sinhlam * sinsig
            yf = np.where(yf == 0.0, gm * 1e-10, yf)
            x  = xf / gm;  y = yf / gm

            P_plane = (gm / np.pi) * (
                y * (np.arctan((x - 1.0) / y) - np.arctan((x + 1.0) / y))
                + 0.5 * (x - 1.0) * np.log((x - 1.0) ** 2 + y ** 2)
                - 0.5 * (x + 1.0) * np.log((x + 1.0) ** 2 + y ** 2)
            )
            P_plane -= A_pt_e * angle ** 2 / 4.0
            P_plane -= shift_pt

            cos_a    = np.where(np.cos(angle) >= 1.0, np.cos(1e-10), np.cos(angle))
            P_sphere = A_pt_e * np.log(1.0 - cos_a)

            angle_s = np.where(angle == 0.0, 1e-10, angle)
            P_xy    = 2.0 * A_pt_e * (np.log(angle_s) - np.log(np.sqrt(2.0)))
            P_xy   -= A_pt_e * angle ** 2 / 4.0

            contrib    = (P_plane - P_xy + P_sphere - Ppt_avg) * pc
            P_grd     += contrib
            Pedge_grd += contrib

    # ── Convert to MPa ────────────────────────────────────────────────────────
    P_grd     /= 1.0e6
    Pwall_grd /= 1.0e6
    Pedge_grd /= 1.0e6

    # ── Pressure derivatives (fully vectorised, no Python loops) ──────────────
    def _gradients(F):
        """Finite-difference gradients of F (MPa), returned as Pa/degree."""
        G = F * 1.0e6
        dlon = np.empty_like(G)
        dlon[:, 1:-1] = (G[:, 2:]  - G[:, :-2]) / (2.0 * spacing)
        dlon[:, 0]    = (G[:, 1]   - G[:, 0])   / spacing
        dlon[:, -1]   = (G[:, -1]  - G[:, -2])  / spacing
        dlat = np.empty_like(G)
        dlat[1:-1, :] = (G[2:, :]  - G[:-2, :]) / (2.0 * spacing)
        dlat[0,    :] = (G[1,  :]  - G[0,  :])  / spacing
        dlat[-1,   :] = (G[-1, :]  - G[-2, :])  / spacing
        return dlon, dlat

    # Pa/degree → same units as original (replicates the scalar loop formula)
    _fac    = 360.0 / (2.0 * np.pi * rad_km * 1.0e3) * 1.0e-3
    cos_lat = np.cos(np.deg2rad(lat_grd))

    dPdl,  dPdb  = _gradients(P_grd)
    dPwdl, dPwdb = _gradients(Pwall_grd)
    dPedl, dPedb = _gradients(Pedge_grd)

    dPdlon_grd     = cos_lat * _fac * dPdl
    dPdlat_grd     =           _fac * dPdb
    dPwalldlon_grd = cos_lat * _fac * dPwdl
    dPwalldlat_grd =           _fac * dPwdb
    dPedgedlon_grd = cos_lat * _fac * dPedl
    dPedgedlat_grd =           _fac * dPedb

    # ── Domain assignment (cached by spacing; geometry is fixed per _cache) ──
    _pk = spacing
    if _pk not in _poly_cache:
        _poly_cache[_pk] = partition_polygon_points(
            lons, lats, bound_ind,
            lona, lata, lonb, latb,
            domain_bounds, rad_km,
        )
    polygons, polygon_points = _poly_cache[_pk]

    # ── Plate/base velocities (vectorised, grouped by domain) ─────────────────
    poly_idx     = polygon_points.astype(int)   # (nlat, nlon)
    plate_vel_ew = np.zeros((nlat, nlon))
    plate_vel_ns = np.zeros((nlat, nlon))
    base_vel_ew  = np.zeros((nlat, nlon))
    base_vel_ns  = np.zeros((nlat, nlon))

    for dom_idx in np.unique(poly_idx):
        mask = (poly_idx == dom_idx)
        dom  = dom_idx - 1
        nd   = ndomain[dom]
        lats_m = lat_grd[mask]
        lons_m = lon_grd[mask]
        if nd == 500:
            plate_vel_ew[mask] = rigid_vew[dom] * 0.1
            plate_vel_ns[mask] = rigid_vns[dom] * 0.1
        else:
            ew, ns = _euler_vel_arr(
                pole_top_lat[dom], pole_top_lon[dom], pole_top_rate[dom],
                lats_m, lons_m,
            )
            plate_vel_ew[mask] = ew
            plate_vel_ns[mask] = ns
        if nd in (200, 300):
            ew_b, ns_b = _euler_vel_arr(
                pole_bott_lat[dom], pole_bott_lon[dom], pole_bott_rate[dom],
                lats_m, lons_m,
            )
            base_vel_ew[mask] = ew_b
            base_vel_ns[mask] = ns_b

    # ── Asthenospheric velocities (fully vectorised, no Python loops) ──────────
    nd_grd  = np.array(ndomain, dtype=int)[poly_idx - 1]   # (nlat, nlon)
    slab_bc = np.isin(nd_grd, [200, 400])                   # True → slab-tail coeff
    thick   = np.where(slab_bc, (ah1 - 2.0*alith)*1.0e-3, (ah1 - alith)*1.0e-3)
    coeff_g = np.where(slab_bc, coefftr2, coefftr1)
    free_b  = np.isin(nd_grd, [100, 400])                   # True → free base

    vel_convert = 100.0 * 60.0 * 60.0 * 24.0 * 365.0
    # f_grd = column-averaged pressure-driven coefficient (used as base for depth scaling)
    f_free  = coeff_g / 3.0  * (1.0 + thick / (8.0 * rad_km)) * vel_convert
    f_fixed = coeff_g / 12.0 * (1.0 + thick / rad_km)          * vel_convert
    f_grd   = np.where(free_b, f_free, f_fixed)

    # ── Depth-dependent velocity profile ──────────────────────────────────────
    # ζ = normalised depth within the asthenosphere [0 = top, 1 = bottom]
    # press_depth and alith are in metres; thick is in km.
    z_asth_km = max(0.0, (press_depth - alith) * 1.0e-3)   # km below lithosphere
    zeta      = np.clip(z_asth_km / thick, 0.0, 1.0)        # (nlat, nlon)

    # Pressure-driven shape function (normalised so column-avg = 1):
    #   free base  (v_plate top, dv/dz=0 bottom): 3ζ(2−ζ)/2
    #   fixed base (v=0 at both boundaries):       6ζ(1−ζ)
    dw_p = np.where(free_b, 1.5 * zeta * (2.0 - zeta), 6.0 * zeta * (1.0 - zeta))

    # Plate/base-driven component at depth ζ:
    #   free base:  v_plate uniform across column (Couette with free bottom)
    #               spherical correction: (1 − z/R) vs column avg (1 − H/2R)
    plate_sfac = 1.0 - z_asth_km / rad_km          # at depth z
    pfac_avg   = 1.0 - thick / (2.0 * rad_km)      # column average (original pfac)
    fv   = 1.0 + thick / (6.0 * rad_km)
    fb   = 1.0 - thick / (6.0 * rad_km)

    def _asthen_vel(dP_ew, dP_ns):
        # Pressure-driven: column-avg × depth shape function
        vp_ew = -dP_ew * 1e3 * f_grd * dw_p
        vp_ns = -dP_ns * 1e3 * f_grd * dw_p
        # Plate/base-driven: depth-specific Couette profile
        p_ew = np.where(free_b,
            plate_vel_ew * plate_sfac,                               # free: uniform
            plate_vel_ew * (1.0 - zeta) * fv + base_vel_ew * zeta * fb)  # fixed: linear
        p_ns = np.where(free_b,
            plate_vel_ns * plate_sfac,
            plate_vel_ns * (1.0 - zeta) * fv + base_vel_ns * zeta * fb)
        return vp_ew + p_ew, vp_ns + p_ns

    avgvel_asthen_ew, avgvel_asthen_ns = _asthen_vel(dPdlon_grd, dPdlat_grd)
    pdrivenvel_wall_ew = -dPwalldlon_grd * 1e3 * f_grd * dw_p
    pdrivenvel_wall_ns = -dPwalldlat_grd * 1e3 * f_grd * dw_p
    pdrivenvel_edge_ew = -dPedgedlon_grd * 1e3 * f_grd * dw_p
    pdrivenvel_edge_ns = -dPedgedlat_grd * 1e3 * f_grd * dw_p

    # ── Trench velocities ─────────────────────────────────────────────────────
    trench_vels = np.empty((0, 4), float)
    num_t = 0
    for k in range(len(lato)):
        if iwall[k] == 1:
            if num_t % 10 == 0:
                trench_vels = np.vstack(
                    [trench_vels, [lono[k], lato[k], vt_ew[k] * 0.1, vt_ns[k] * 0.1]]
                )
            num_t += 1

    # ── Depth correction (same scaling as original) ───────────────────────────
    P_zfactor = rad_km / rad_pressloc

    return (
        P_grd      * P_zfactor,
        Pwall_grd  * P_zfactor,
        Pedge_grd  * P_zfactor,
        dPdlon_grd * P_zfactor,
        dPdlat_grd * P_zfactor,
        polygon_points,
        plate_vel_ew,
        plate_vel_ns,
        trench_vels,
        avgvel_asthen_ew,
        avgvel_asthen_ns,
        pdrivenvel_wall_ew,
        pdrivenvel_wall_ns,
        pdrivenvel_edge_ew,
        pdrivenvel_edge_ns,
        lon_grd,
        lat_grd,
        polygons,
    )


# ── Vectorized helpers for _fast_buildmatrix ──────────────────────────────────
_PI = math.pi


def _hav_km(lon1, lat1, lon2, lat2, rad_km):
    """Haversine distance in km; inputs may be scalars or arrays (degrees)."""
    lon1r = np.radians(lon1); lat1r = np.radians(lat1)
    lon2r = np.radians(lon2); lat2r = np.radians(lat2)
    dlat  = lat2r - lat1r;    dlon  = lon2r - lon1r
    a = (np.sin(dlat / 2.0) ** 2
         + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2.0) ** 2)
    return 2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))) * rad_km


def _project_vec(lon1_rad, lat1_rad, azim, len_ratio):
    """Vectorized project_to_point. All inputs may be arrays. Returns (lon_rad, lat_rad)."""
    lat2 = np.arcsin(
        np.sin(lat1_rad) * np.cos(len_ratio)
        + np.cos(lat1_rad) * np.sin(len_ratio) * np.cos(azim)
    )
    lon2 = lon1_rad + np.arctan2(
        np.sin(azim) * np.sin(len_ratio) * np.cos(lat1_rad),
        np.cos(len_ratio) - np.sin(lat1_rad) * np.sin(lat2),
    )
    lon2 = (lon2 + 2.0 * _PI) % (2.0 * _PI)
    return lon2, lat2


def _midpoint_arr(latA, lonA, latB, lonB):
    """Vectorized midpoint; inputs may be scalars or arrays (degrees). Returns (lat_deg, lon_deg)."""
    lonA_r = np.radians(lonA); lonB_r = np.radians(lonB)
    latA_r = np.radians(latA); latB_r = np.radians(latB)
    dLon   = lonB_r - lonA_r
    Bx     = np.cos(latB_r) * np.cos(dLon)
    By     = np.cos(latB_r) * np.sin(dLon)
    latC   = np.arctan2(np.sin(latA_r) + np.sin(latB_r),
                        np.sqrt((np.cos(latA_r) + Bx) ** 2 + By ** 2))
    lonC   = lonA_r + np.arctan2(By, np.cos(latA_r) + Bx)
    lonC   = (lonC + 3.0 * _PI) % (2.0 * _PI) - _PI
    return np.degrees(latC), np.degrees(lonC)


def _findpressure_edge_vec(lonobs, latobs, lonaa, lataa, lonbb, latbb, gm, rad_km):
    """Vectorized findpressure_edge. lonobs/latobs may be arrays; all other args scalars."""
    latmid, lonmid = _midpoint_arr(lataa, lonaa, latbb, lonbb)
    dist   = _hav_km(lonmid, latmid, lonobs, latobs, rad_km)
    angle  = dist / rad_km
    rad    = rad_km * 1.0e3

    A_pt     = -gm / _PI
    shift_pt = (-2.0 * gm / _PI) * (1.0 + math.log(rad * math.sqrt(2.0) / gm))
    Ppt_avg  = A_pt * (math.log(2.0) - 1.0)

    dista = _hav_km(lonaa, lataa, lonobs, latobs, rad_km) * 1.0e3
    distb = _hav_km(lonbb, latbb, lonobs, latobs, rad_km) * 1.0e3

    coshlam = 0.5 / gm * (dista + distb)
    sinhlam = np.sqrt(np.maximum(coshlam ** 2 - 1.0, 0.0))
    cossig  = np.clip(0.5 / gm * (dista - distb), -1.0, 1.0)
    sinsig  = np.sqrt(1.0 - cossig ** 2)

    xf = gm * coshlam * cossig
    yf = gm * sinhlam * sinsig
    yf = np.where(yf == 0.0, gm * 1.0e-10, yf)
    x  = xf / gm;  y = yf / gm

    P_plane = (gm / _PI) * (
        y * (np.arctan((x - 1.0) / y) - np.arctan((x + 1.0) / y))
        + 0.5 * (x - 1.0) * np.log((x - 1.0) ** 2 + y ** 2)
        - 0.5 * (x + 1.0) * np.log((x + 1.0) ** 2 + y ** 2)
    )
    P_plane -= A_pt * angle ** 2 / 4.0
    P_plane -= shift_pt

    angle_s  = np.where(np.cos(angle) == 1.0, 1.0e-10, angle)
    P_sphere = A_pt * np.log(1.0 - np.cos(angle_s))
    P_xy     = 2.0 * A_pt * (np.log(angle_s) - math.log(math.sqrt(2.0)))
    P_xy    -= A_pt * angle ** 2 / 4.0

    return P_plane - P_xy + P_sphere - Ppt_avg


def _findpressure_wall_vec(lonobs, latobs, lonaa, lataa, lonbb, latbb, gm, alp, rad_km):
    """Vectorized findpressure_wall. lonobs/latobs may be arrays; all other args scalars."""
    latmid, lonmid = _midpoint_arr(lataa, lonaa, latbb, lonbb)
    angle  = _hav_km(float(lonmid), float(latmid), lonobs, latobs, rad_km) / rad_km
    alp_n  = (alp + 2.0 * _PI) % (2.0 * _PI)
    rad    = rad_km * 1.0e3
    A_pt   = gm / (4.0 * rad)

    lat_rad  = np.radians(latobs)
    lon_rad  = (np.radians(lonobs - float(lonmid) - 180.0) + 2.0 * _PI) % (2.0 * _PI)
    dlat_rad = math.radians(90.0 - float(latmid))
    lons_t   = np.arctan2(
        np.sin(lon_rad),
        np.tan(lat_rad) * math.sin(dlat_rad) + np.cos(lon_rad) * math.cos(dlat_rad),
    )
    pt_azim = ((lons_t + alp_n) + 2.0 * _PI) % (2.0 * _PI)

    dista  = _hav_km(lonaa, lataa, lonobs, latobs, rad_km) * 1.0e3
    distb  = _hav_km(lonbb, latbb, lonobs, latobs, rad_km) * 1.0e3
    coshlam = 0.5 / gm * (dista + distb)
    sinhlam = np.sqrt(np.maximum(coshlam ** 2 - 1.0, 0.0))
    cossig  = np.clip(0.5 / gm * (dista - distb), -1.0, 1.0)
    sinsig  = np.sqrt(1.0 - cossig ** 2)
    sinsig  = np.where(pt_azim <= _PI, -sinsig, sinsig)

    t       = coshlam - sinhlam
    P_plane = (gm * t * sinsig
               - (1.0 / 3.0) * gm * t ** 3 * (3.0 * sinsig * cossig ** 2 - sinsig ** 3))

    angle_s  = np.where(np.cos(angle) == 1.0, 1.0e-10, angle)
    P_sphere = A_pt * gm * np.sin(-pt_azim) * (np.sin(angle_s) / (1.0 - np.cos(angle_s)))
    P_xy     = A_pt * gm * 2.0 * (1.0 / angle_s) * np.sin(-pt_azim)

    return P_plane - P_xy + P_sphere


def _calcvel_wall_vec(lonobsa, latobsa, lonobsb, latobsb,
                      azimuth_j, lonaa, lataa, lonbb, latbb,
                      gm_src, alp_src, dsegtr, elit_j, ebig_j, side, rad_km):
    """
    Vectorized calcvel_wall.
    lonobsa/b, latobsa/b, azimuth_j, elit_j, ebig_j: arrays (n_segs,)
    All source-segment args (lonaa/b, lataa/b, gm_src, alp_src): scalars.
    side: scalar (1=right, 0=left).
    """
    # Observation segment midpoints
    latobs_mid, lonobs_mid = _midpoint_arr(latobsa, lonobsa, latobsb, lonobsb)
    # Source midpoint (scalar)
    lat_src_mid, lon_src_mid = _midpoint_arr(lataa, lonaa, latbb, lonbb)

    angle = _hav_km(lonobs_mid, latobs_mid, float(lon_src_mid), float(lat_src_mid), rad_km) / rad_km

    thresh = 0.5 * (12.5 * 0.5 * dsegtr * 1.0e3 + 1500.0e3)
    near   = (angle * rad_km * 1.0e3) <= thresh

    # ── Near (planar) branch ───────────────────────────────────────────────────
    gamma_p = 1.0e3 * _hav_km(lonobsa, latobsa, lonobsb, latobsb, rad_km) / 2.0

    da_a = _hav_km(lonaa, lataa, lonobsa, latobsa, rad_km) * 1.0e3
    db_a = _hav_km(lonbb, latbb, lonobsa, latobsa, rad_km) * 1.0e3
    coshlama = 0.5 / gm_src * (da_a + db_a)
    sinhlama = np.sqrt(np.maximum(coshlama ** 2 - 1.0, 0.0))
    cossiga  = np.clip(0.5 / gm_src * (da_a - db_a), -1.0, 1.0)
    sinsiga  = np.sqrt(1.0 - cossiga ** 2)

    da_b = _hav_km(lonaa, lataa, lonobsb, latobsb, rad_km) * 1.0e3
    db_b = _hav_km(lonbb, latbb, lonobsb, latobsb, rad_km) * 1.0e3
    coshlamb = 0.5 / gm_src * (da_b + db_b)
    sinhlamb = np.sqrt(np.maximum(coshlamb ** 2 - 1.0, 0.0))
    cossigb  = np.clip(0.5 / gm_src * (da_b - db_b), -1.0, 1.0)
    sinsigb  = np.sqrt(1.0 - cossigb ** 2)

    explama  = coshlama - sinhlama
    explamb  = coshlamb - sinhlamb
    safe_gp  = np.where(gamma_p == 0.0, 1.0e-30, gamma_p)
    fvel_plane = (
        gm_src * (explama * cossiga - explamb * cossigb) / (2.0 * safe_gp)
        - (1.0 / 3.0) * gm_src * explama ** 3 * (cossiga ** 3 - 3.0 * cossiga * sinsiga ** 2) / (2.0 * safe_gp)
        + (1.0 / 3.0) * gm_src * explamb ** 3 * (cossigb ** 3 - 3.0 * cossigb * sinsigb ** 2) / (2.0 * safe_gp)
    )

    # ── Far (spherical) branch ─────────────────────────────────────────────────
    lr_close = elit_j / (rad_km * 1.0e3)
    ebig_f   = elit_j + 1.0e3
    lr_far   = ebig_f / (rad_km * 1.0e3)
    pt_azim  = azimuth_j + (_PI / 2.0 if side == 1 else -_PI / 2.0)

    lon_close, lat_close = _project_vec(np.radians(lonobs_mid), np.radians(latobs_mid), pt_azim, lr_close)
    lon_far,   lat_far   = _project_vec(np.radians(lonobs_mid), np.radians(latobs_mid), pt_azim, lr_far)

    Pclose = _findpressure_wall_vec(
        np.degrees(lon_close), np.degrees(lat_close),
        lonaa, lataa, lonbb, latbb, gm_src, alp_src, rad_km,
    )
    Pfar = _findpressure_wall_vec(
        np.degrees(lon_far), np.degrees(lat_far),
        lonaa, lataa, lonbb, latbb, gm_src, alp_src, rad_km,
    )
    if side == 1:
        fvel_sphere = (Pfar - Pclose) / (ebig_f - elit_j)
    else:
        fvel_sphere = (Pclose - Pfar) / (ebig_f - elit_j)

    return np.where(near, fvel_plane, fvel_sphere)


def _fast_buildmatrix(lona, lata, lonb, latb, gam, alpha, lono, lato,
                      iwall, idl, idr, n_segs, num_segs,
                      coeff1, coeff2, coefftr1, coefftr2,
                      ndomain, epslit, dsegtr, rad_km, alith, ah1, eps_fact):
    """
    Vectorized replacement for functions.buildmatrix.

    Keeps the outer for-iset loop (source segments); the inner for-jobs loop
    (observation segments) is replaced by NumPy operations over all n_segs
    entries simultaneously.  Gives a ~n_segs× speedup on the matrix build.
    """
    gam_j   = np.asarray(gam,   dtype=float)
    alpha_j = np.asarray(alpha, dtype=float)
    iwall_j = np.asarray(iwall, dtype=int)
    lono_j  = np.asarray(lono,  dtype=float)
    lato_j  = np.asarray(lato,  dtype=float)
    lona_j  = np.asarray(lona,  dtype=float)
    lata_j  = np.asarray(lata,  dtype=float)
    lonb_j  = np.asarray(lonb,  dtype=float)
    latb_j  = np.asarray(latb,  dtype=float)
    idl_j   = np.asarray(idl,   dtype=int)
    idr_j   = np.asarray(idr,   dtype=int)

    elit_j    = eps_fact * gam_j
    ebig_j    = elit_j + 1.0e3
    len_lit_j = elit_j / (rad_km * 1.0e3)
    len_big_j = ebig_j / (rad_km * 1.0e3)

    lona_rad_j = np.radians(lona_j);  lata_rad_j = np.radians(lata_j)
    lonb_rad_j = np.radians(lonb_j);  latb_rad_j = np.radians(latb_j)

    # per-jobs ndomain lookups
    nd_idl = np.array([ndomain[idl_j[j] - 1] for j in range(n_segs)], dtype=int)
    nd_idr = np.array([ndomain[idr_j[j] - 1] for j in range(n_segs)], dtype=int)

    # per-jobs secondary geometry coefficients (depend on iwall[jobs])
    asth_t  = (ah1 - alith)     * 1.0e-3
    asth_t2 = (ah1 - 2 * alith) * 1.0e-3
    is_wall = (iwall_j == 1)

    co1_j   = np.where(is_wall, coefftr1, coeff1)
    co2_j   = np.where(is_wall, coefftr2, coeff2)
    scfb_j  = np.where(is_wall,
        (2.0 * (8 * rad_km + asth_t))  / (8.0 * (2 * rad_km - asth_t)),
        1.0 + asth_t  / (8 * rad_km))
    scfb2_j = np.where(is_wall,
        (2.0 * (8 * rad_km + asth_t2)) / (8.0 * (2 * rad_km - asth_t2)),
        1.0 + asth_t2 / (8 * rad_km))
    scxb_j  = np.where(is_wall,
        (2 * rad_km + 2 * asth_t)  / (2 * rad_km - asth_t),
        1.0 + asth_t  / rad_km)
    scxb2_j = np.where(is_wall,
        (2 * rad_km + 2 * asth_t2) / (2 * rad_km - asth_t2),
        1.0 + asth_t2 / rad_km)

    jobs_idx = np.arange(n_segs)
    pkernel  = np.zeros((n_segs, n_segs))

    for iset in range(n_segs):
        lonaa = float(lona[iset]);  lataa = float(lata[iset])
        lonbb = float(lonb[iset]);  latbb = float(latb[iset])
        gm_s  = float(gam[iset])
        alp_s = float(alpha[iset])

        # ── Compute deriv1, deriv2 for all jobs simultaneously ────────────────
        if iwall[iset] == 1 and iset >= num_segs:          # wall source
            obs_azim_r = alpha_j + _PI / 2.0
            la_r, lo_r = _project_vec(lona_rad_j, lata_rad_j, obs_azim_r, len_lit_j)
            lb_r, lo_rb = _project_vec(lonb_rad_j, latb_rad_j, obs_azim_r, len_lit_j)
            deriv1 = _calcvel_wall_vec(
                np.degrees(la_r), np.degrees(lo_r),
                np.degrees(lb_r), np.degrees(lo_rb),
                alpha_j, lonaa, lataa, lonbb, latbb,
                gm_s, alp_s, dsegtr, elit_j, ebig_j, 1, rad_km,
            )
            obs_azim_l = alpha_j - _PI / 2.0
            la_l, lo_l = _project_vec(lona_rad_j, lata_rad_j, obs_azim_l, len_lit_j)
            lb_l, lo_lb = _project_vec(lonb_rad_j, latb_rad_j, obs_azim_l, len_lit_j)
            deriv2 = _calcvel_wall_vec(
                np.degrees(la_l), np.degrees(lo_l),
                np.degrees(lb_l), np.degrees(lo_lb),
                alpha_j, lonaa, lataa, lonbb, latbb,
                gm_s, alp_s, dsegtr, elit_j, ebig_j, 0, rad_km,
            )

        elif iwall[iset] == 1 or iwall[iset] == 0:         # edge source
            lono_r = np.radians(lono_j);  lato_r = np.radians(lato_j)

            obs_azim_r = alpha_j + _PI / 2.0
            lon_rf, lat_rf = _project_vec(lono_r, lato_r, obs_azim_r, len_big_j)
            Pr2 = _findpressure_edge_vec(np.degrees(lon_rf), np.degrees(lat_rf),
                                         lonaa, lataa, lonbb, latbb, gm_s, rad_km)
            lon_rn, lat_rn = _project_vec(lono_r, lato_r, obs_azim_r, len_lit_j)
            Pr1 = _findpressure_edge_vec(np.degrees(lon_rn), np.degrees(lat_rn),
                                         lonaa, lataa, lonbb, latbb, gm_s, rad_km)
            deriv1 = (Pr2 - Pr1) / (ebig_j - elit_j)

            obs_azim_l = alpha_j - _PI / 2.0
            lon_lf, lat_lf = _project_vec(lono_r, lato_r, obs_azim_l, len_big_j)
            Pl2 = _findpressure_edge_vec(np.degrees(lon_lf), np.degrees(lat_lf),
                                         lonaa, lataa, lonbb, latbb, gm_s, rad_km)
            lon_ln, lat_ln = _project_vec(lono_r, lato_r, obs_azim_l, len_lit_j)
            Pl1 = _findpressure_edge_vec(np.degrees(lon_ln), np.degrees(lat_ln),
                                         lonaa, lataa, lonbb, latbb, gm_s, rad_km)
            deriv2 = (Pl1 - Pl2) / (ebig_j - elit_j)

        else:                                               # strike-slip
            deriv1 = np.zeros(n_segs)
            deriv2 = np.zeros(n_segs)

        # ── Domain weighting (vectorised over jobs) ────────────────────────────
        def _weight(nd):
            return np.where((nd == 100) | (nd == 500), co1_j * scfb_j  / 3.0,
                   np.where( nd == 400,                co2_j * scfb2_j / 3.0,
                   np.where( nd == 200,                co2_j * scxb2_j / 12.0,
                   np.where( nd == 300,                co1_j * scxb_j  / 12.0,
                             np.zeros(n_segs)))))

        d1 = deriv1 * _weight(nd_idr)
        d2 = deriv2 * _weight(nd_idl)

        # ── Fill matrix column ─────────────────────────────────────────────────
        col = np.where(iwall_j == 2, 0.0,
              np.where(iwall_j == 0, d1 - d2,
              np.where((iwall_j == 1) & (jobs_idx < num_segs), d1,
                       d2)))
        pkernel[:, iset] = col

    return pkernel


def _fast_outputDP(lona, lata, lonb, latb, lono, lato, iwall, gam, alpha,
                   n_segs, num_segs, pcoeff, rad_km,
                   lon_subslab, lat_subslab, lon_wedge, lat_wedge,
                   polarity, vtopl, vtopr, vt, dip_depth):
    """
    Vectorized replacement for functions.outputDP.

    The original has O(n_obs × n_segs) scalar Python calls.  Here we loop over
    source segments (n_segs iterations) and evaluate each Green's-function
    kernel on all observation points at once via _findpressure_{edge,wall}_vec,
    giving O(n_segs) NumPy calls instead.
    """
    n_obs = len(lon_wedge)
    rad_diploc = rad_km - dip_depth / 1.0e3
    z_factor   = rad_km / rad_diploc

    lon_ss = np.asarray(lon_subslab, dtype=float)
    lat_ss = np.asarray(lat_subslab, dtype=float)
    lon_w  = np.asarray(lon_wedge,   dtype=float)
    lat_w  = np.asarray(lat_wedge,   dtype=float)

    sumpress_ss = np.zeros(n_obs)
    sumpress_w  = np.zeros(n_obs)

    for iset in range(n_segs):
        lonaa = float(lona[iset]); lataa = float(lata[iset])
        lonbb = float(lonb[iset]); latbb = float(latb[iset])
        gm    = float(gam[iset])
        pc    = float(pcoeff[iset])

        if iset >= num_segs:   # slab-wall kernel
            alp = float(alpha[iset])
            sumpress_ss += _findpressure_wall_vec(
                lon_ss, lat_ss, lonaa, lataa, lonbb, latbb, gm, alp, rad_km) * pc
            sumpress_w  += _findpressure_wall_vec(
                lon_w,  lat_w,  lonaa, lataa, lonbb, latbb, gm, alp, rad_km) * pc
        else:                  # plate-edge kernel
            sumpress_ss += _findpressure_edge_vec(
                lon_ss, lat_ss, lonaa, lataa, lonbb, latbb, gm, rad_km) * pc
            sumpress_w  += _findpressure_edge_vec(
                lon_w,  lat_w,  lonaa, lataa, lonbb, latbb, gm, rad_km) * pc

    ss_z = sumpress_ss * z_factor
    w_z  = sumpress_w  * z_factor

    vel_term = 1.0e-3 / (365.0 * 24.0 * 60.0 * 60.0)
    DP = np.zeros((n_obs, 7))
    DP[:, 0] = np.asarray(lono, dtype=float)[:n_obs]
    DP[:, 1] = np.asarray(lato, dtype=float)[:n_obs]
    DP[:, 2] = ss_z / 1.0e6
    DP[:, 3] = w_z  / 1.0e6
    DP[:, 4] = (ss_z - w_z) / 1.0e6
    DP[:, 5] = np.abs(np.asarray(vtopl, dtype=float)[:n_obs]
                      - np.asarray(vtopr, dtype=float)[:n_obs]) / vel_term
    pol = np.asarray(polarity, dtype=float)[:n_obs]
    vt_arr = np.asarray(vt, dtype=float)[:n_obs]
    DP[:, 6] = np.where(pol == 1, -vt_arr / vel_term, vt_arr / vel_term)
    return DP


def _run_outputgrids(amu, press_depth, grid_spacing, pcoeff, geom):
    """Call outputgrids + outputDP with the given pcoeff and geometry cache."""
    coefftr1 = (_AH1 - _ALITH)    ** 2 / amu
    coefftr2 = (_AH1 - 2 * _ALITH) ** 2 / amu

    (P_out,
     _Pwall, _Pedge,
     _dPlon, _dPlat,
     _polygon_pts,
     _pvel_ew, _pvel_ns,
     _trench_vels,
     vel_ew, vel_ns,
     _pdwall_ew, _pdwall_ns,
     _pdedge_ew, _pdedge_ns,
     lons_grd, lats_grd,
     _polygons) = _fast_outputgrids(
        grid_spacing,
        geom['lona'], geom['lata'], geom['lonb'], geom['latb'],
        geom['lono'], geom['lato'], geom['iwall'], geom['gam'], geom['alpha'],
        amu, _AH1 - _ALITH,
        geom['n_segs'], geom['num_segs'], pcoeff,
        _RAD_KM, geom['domain_bounds'], geom['bound_ind'],
        geom['pole_top_lon'], geom['pole_top_lat'], geom['pole_top_rate'],
        geom['vt_ew'], geom['vt_ns'], _ALITH, press_depth,
        coefftr1, coefftr2,
        geom['pole_bott_lon'], geom['pole_bott_lat'], geom['pole_bott_rate'],
        geom['rigid_vew'], geom['rigid_vns'], _AH1, geom['ndomain'],
    )

    DP = _fast_outputDP(
        geom['lona'], geom['lata'], geom['lonb'], geom['latb'],
        geom['lono'], geom['lato'], geom['iwall'], geom['gam'], geom['alpha'],
        geom['n_segs'], geom['num_segs'], pcoeff, _RAD_KM,
        geom['lon_subslab'], geom['lat_subslab'],
        geom['lon_wedge'], geom['lat_wedge'],
        geom['polarity'], geom['vtopl'], geom['vtopr'], geom['vt'], _DIP_DEPTH,
    )
    return P_out, vel_ew, vel_ns, DP, lons_grd, lats_grd


def solve_only(payload: dict) -> dict:
    """
    Phase 1: read inputs, build matrix, invert.
    Returns immediately from cache if the heavy params haven't changed.
    Result contains only {cached, timing_s} — no grid data yet.
    """
    global _cache
    t0 = time.perf_counter()

    model                 = payload.get('model',  'LargeSP_RetreatingTrench')
    amu                   = float(payload.get('viscosity',              3e20))
    basal_bc              =       payload.get('basal_bc',             'free')
    flux_slab             = int(  payload.get('flux_slab',               0))
    flux_width            = float(payload.get('flux_width',        500_000.0))
    flux_alpha            = float(payload.get('flux_alpha',              0.0))
    no_flux_for_slabtails = int(  payload.get('no_flux_for_slabtails',   1))
    euler_lat  = payload.get('euler_lat',  None)
    euler_lon  = payload.get('euler_lon',  None)
    euler_rate = payload.get('euler_rate', None)
    trench_vew = payload.get('trench_vew', None)
    trench_vns = payload.get('trench_vns', None)

    width_deg  = float(payload.get('width_deg',  60.0)) if model == _SIMPLE_MODEL else None
    length_deg = float(payload.get('length_deg', 45.0)) if model == _SIMPLE_MODEL else None

    cache_key = (
        model, basal_bc, flux_slab, flux_width, flux_alpha,
        no_flux_for_slabtails,
        round(float(euler_lat),  4) if euler_lat  is not None else None,
        round(float(euler_lon),  4) if euler_lon  is not None else None,
        round(float(euler_rate), 4) if euler_rate is not None else None,
        round(float(trench_vew), 4) if trench_vew is not None else None,
        round(float(trench_vns), 4) if trench_vns is not None else None,
        round(width_deg,  2) if width_deg  is not None else None,
        round(length_deg, 2) if length_deg is not None else None,
    )

    if _cache is not None and _cache['key'] == cache_key:
        return {'cached': True, 'timing_s': time.perf_counter() - t0}

    _poly_cache.clear()  # geometry changed — invalidate polygon cache

    # ── Generate Simple model input files if needed ──────────────────────────
    if model == _SIMPLE_MODEL:
        generate_simple_inp_files(width_deg, length_deg)

    # ── Full matrix solve ────────────────────────────────────────────────────
    coeff1   = (_AH1 - _ALITH)    ** 3 / (amu * _AH1)
    coeff2   = (_AH1 - 2 * _ALITH) ** 3 / (amu * _AH1)
    coefftr1 = (_AH1 - _ALITH)    ** 2 / amu
    coefftr2 = (_AH1 - 2 * _ALITH) ** 2 / amu

    inp_dir   = _FC_DIR / 'inputs'
    f_bounds  = str(inp_dir / f'Subbon_{model}.inp')
    f_domains = str(inp_dir / f'Subfil_{model}.inp')
    f_grid    = str(inp_dir / 'Subgrd_Fast.inp')

    (ndomain,
     pole_top_lon, pole_top_lat, pole_top_rate,
     pole_bott_lon, pole_bott_lat, pole_bott_rate,
     rigid_vew, rigid_vns,
     domain_bounds) = readdomains(f_domains)

    if basal_bc == 'no_slip':
        ndomain = [300 if v == 100 else v for v in ndomain]

    grid_spacing_default, _prof, dsegtr, dseged = readgrid(f_grid)

    (num_bounds,
     iwall, lona, lata, lonb, latb, bound_ind,
     idl, idr, vt_ew, vt_ns, polarity, large_wall_inds) = readbounds(f_bounds)

    # Apply Euler pole override to SP domain (non-Earth models only)
    if all(v is not None for v in [euler_lat, euler_lon, euler_rate]):
        sp_dom = next((idl[i]-1 for i in range(num_bounds) if iwall[i] == 1), None)
        if sp_dom is not None:
            pole_top_lat[sp_dom]  = float(euler_lat)
            pole_top_lon[sp_dom]  = float(euler_lon)
            pole_top_rate[sp_dom] = float(euler_rate)
            pole_bott_lat[sp_dom]  = float(euler_lat)
            pole_bott_lon[sp_dom]  = float(euler_lon)
            pole_bott_rate[sp_dom] = float(euler_rate)

    if trench_vew is not None and trench_vns is not None:
        for i in range(num_bounds):
            if iwall[i] == 1:
                vt_ew[i] = float(trench_vew)
                vt_ns[i] = float(trench_vns)

    orig_bounds = [
        {
            'type':        int(iwall[i]),
            'lona':        _to_float(lona[i]),
            'lata':        _to_float(lata[i]),
            'lonb':        _to_float(lonb[i]),
            'latb':        _to_float(latb[i]),
            'left_domain': int(idl[i]),
            'right_domain':int(idr[i]),
            'polarity':    int(polarity[i]),
        }
        for i in range(num_bounds)
    ]

    (n_segs, num_segs,
     iwall, idl, idr,
     lona, lata, lonb, latb,
     bound_ind, large_wall_inds,
     vt_ew, vt_ns, polarity,
     num_wall_segs) = organizebounds(
        num_bounds, iwall, idl, idr,
        lona, lata, lonb, latb,
        bound_ind, large_wall_inds,
        vt_ew, vt_ns,
        dsegtr, dseged, polarity,
        _RAD_KM, _ISEG_MIN,
    )

    (lono, lato, gam, alpha,
     vtopl, vtopr, vbotl, vbotr, vt,
     lon_subslab, lat_subslab,
     lon_wedge, lat_wedge) = pressurepoints(
        lona, lata, lonb, latb,
        vt_ew, vt_ns,
        iwall, idl, idr, n_segs,
        pole_top_lon, pole_top_lat, pole_top_rate,
        pole_bott_lon, pole_bott_lat, pole_bott_rate,
        rigid_vew, rigid_vns,
        ndomain, _EPSLRC, _RAD_KM, _ALITH,
        _SHIFT_EDGES, polarity, _EPSDP_FACT,
    )

    pkernel = _fast_buildmatrix(
        lona, lata, lonb, latb,
        gam, alpha, lono, lato,
        iwall, idl, idr, n_segs, num_segs,
        coeff1, coeff2, coefftr1, coefftr2,
        ndomain, _EPSLIT, dsegtr,
        _RAD_KM, _ALITH, _AH1, _EPS_FACT,
    )

    vector = buildvector(
        iwall, alpha, ndomain, idl, idr,
        vtopl, vtopr, vbotl, vbotr, vt,
        n_segs, num_segs,
        flux_slab, _FLUX_VEL_CONST,
        flux_width, polarity, flux_alpha,
        no_flux_for_slabtails,
        _RAD_KM, _ALITH, _AH1,
    )

    pcoeff = inv(pkernel).dot(vector)

    geom = dict(
        lona=lona, lata=lata, lonb=lonb, latb=latb,
        lono=lono, lato=lato, iwall=iwall, gam=gam, alpha=alpha,
        n_segs=n_segs, num_segs=num_segs,
        domain_bounds=domain_bounds, bound_ind=bound_ind,
        pole_top_lon=pole_top_lon, pole_top_lat=pole_top_lat,
        pole_top_rate=pole_top_rate,
        vt_ew=vt_ew, vt_ns=vt_ns,
        pole_bott_lon=pole_bott_lon, pole_bott_lat=pole_bott_lat,
        pole_bott_rate=pole_bott_rate,
        rigid_vew=rigid_vew, rigid_vns=rigid_vns,
        ndomain=ndomain,
        lon_subslab=lon_subslab, lat_subslab=lat_subslab,
        lon_wedge=lon_wedge, lat_wedge=lat_wedge,
        polarity=polarity, vtopl=vtopl, vtopr=vtopr, vt=vt,
    )
    _cache = dict(
        key=cache_key,
        amu_base=amu,
        pcoeff_base=pcoeff.copy(),
        geom=geom,
        orig_bounds=orig_bounds,
        grid_spacing_default=grid_spacing_default,
    )

    return {'cached': False, 'timing_s': time.perf_counter() - t0}


def grid_only(payload: dict) -> dict:
    """
    Phase 2: run outputgrids on the cached solve state.
    Must be called after solve_only has populated _cache.
    """
    if _cache is None:
        raise RuntimeError('No cached solve — call /compute/solve first')

    t0          = time.perf_counter()
    amu         = float(payload.get('viscosity',    _cache['amu_base']))
    press_depth = float(payload.get('press_depth',  _PRESS_DEPTH))
    gsd         =       payload.get('grid_spacing_deg', None)
    grid_spacing = float(gsd) if gsd is not None else _cache['grid_spacing_default']

    pcoeff = _cache['pcoeff_base'] * (amu / _cache['amu_base'])

    # Restore polygon cache from solve state so partition_polygon_points is
    # only called once per spacing, even across subprocess spawns.
    global _poly_cache
    _poly_cache = dict(_cache.get('poly', {}))

    P_out, vel_ew, vel_ns, DP, lons_grd, lats_grd = _run_outputgrids(
        amu, press_depth, grid_spacing, pcoeff, _cache['geom'])

    # Persist any newly computed polygon entries back into the cache so
    # the main process can pass them to future subprocesses.
    _cache['poly'] = dict(_poly_cache)

    return _serialise_outputs(
        P_out, vel_ew, vel_ns, DP, _cache['geom']['iwall'],
        lons_grd, lats_grd, _cache['orig_bounds'], t0)


def run_computation(payload: dict) -> dict:
    """Single-call wrapper — solve then grid (kept for compatibility)."""
    t0 = time.perf_counter()
    solve_result = solve_only(payload)
    grid_result  = grid_only(payload)
    grid_result['cached']   = solve_result['cached']
    grid_result['timing_s'] = time.perf_counter() - t0
    return grid_result
