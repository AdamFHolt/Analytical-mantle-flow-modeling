"""
compute.py — thin wrapper that calls the existing flow_computations pipeline
and returns JSON-serialisable results.

This module is imported by app.py, which already inserts the flow_computations/
directory into sys.path before importing this module.
"""
import sys
import time
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
)

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
_PRESS_DEPTH = 330.0e3       # depth at which pressure field is evaluated, m
_DIP_DEPTH   = 330.0e3       # depth for DP calculation, m


def _to_float(v):
    """Convert numpy scalar to plain Python float for JSON serialisation."""
    return float(v)


def get_geometry(model: str) -> dict:
    """
    Read only the boundary and domain files for a model — no computation.
    Returns boundaries suitable for immediate map display.
    """
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

    boundaries = [
        {
            'type':         int(iwall[i]),
            'lona':         _to_float(lona[i]),
            'lata':         _to_float(lata[i]),
            'lonb':         _to_float(lonb[i]),
            'latb':         _to_float(latb[i]),
            'left_domain':  int(idl[i]),
            'right_domain': int(idr[i]),
        }
        for i in range(num_bounds)
    ]
    return {'boundaries': boundaries}


def run_computation(payload: dict) -> dict:
    """
    Run the BEM pressure-computation pipeline for a named plate model.

    Parameters (from JSON payload)
    --------------------------------
    model               : str   plate model name (default: LargeSP_RetreatingTrench)
    viscosity           : float asthenospheric viscosity, Pa·s  (default: 3e20)
    flux_slab           : int   0=none, 1=constant-vel flux, 2=convergence-rate flux
    flux_width          : float total slab-flux channel width, m  (default: 500 000)
    flux_alpha          : float fraction of flux on SP side       (default: 0)
    no_flux_for_slabtails: int  1 = suppress flux where slab tail present
    """
    t0 = time.perf_counter()

    model                 = payload.get('model',  'LargeSP_RetreatingTrench')
    amu                   = float(payload.get('viscosity',              3e20))
    basal_bc              =       payload.get('basal_bc',             'free')
    flux_slab             = int(  payload.get('flux_slab',               2))
    flux_width            = float(payload.get('flux_width',        500_000.0))
    flux_alpha            = float(payload.get('flux_alpha',              0.0))
    no_flux_for_slabtails = int(  payload.get('no_flux_for_slabtails',   1))

    # ── Derived coefficients ────────────────────────────────────────────────
    coeff1   = (_AH1 - _ALITH)    ** 3 / (amu * _AH1)
    coeff2   = (_AH1 - 2 * _ALITH) ** 3 / (amu * _AH1)
    coefftr1 = (_AH1 - _ALITH)    ** 2 / amu
    coefftr2 = (_AH1 - 2 * _ALITH) ** 2 / amu

    # ── Input file paths (absolute, relative to flow_computations/) ─────────
    inp_dir   = _FC_DIR / 'inputs'
    f_bounds  = str(inp_dir / f'Subbon_{model}.inp')
    f_domains = str(inp_dir / f'Subfil_{model}.inp')
    f_grid    = str(inp_dir / 'Subgrd_Fast.inp')

    # ── Read inputs ─────────────────────────────────────────────────────────
    (ndomain,
     pole_top_lon, pole_top_lat, pole_top_rate,
     pole_bott_lon, pole_bott_lat, pole_bott_rate,
     rigid_vew, rigid_vns,
     domain_bounds) = readdomains(f_domains)

    # Remap domain type codes for basal boundary condition:
    #   100 = free-slip base (default), 300 = no-slip base.
    if basal_bc == 'no_slip':
        ndomain = [300 if v == 100 else v for v in ndomain]

    grid_spacing, _prof, dsegtr, dseged = readgrid(f_grid)

    (num_bounds,
     iwall, lona, lata, lonb, latb, bound_ind,
     idl, idr, vt_ew, vt_ns, polarity, large_wall_inds) = readbounds(f_bounds)

    # Save original (unsegmented) boundaries for the frontend.
    orig_bounds = [
        {
            'type':        int(iwall[i]),
            'lona':        _to_float(lona[i]),
            'lata':        _to_float(lata[i]),
            'lonb':        _to_float(lonb[i]),
            'latb':        _to_float(latb[i]),
            'left_domain': int(idl[i]),
            'right_domain':int(idr[i]),
        }
        for i in range(num_bounds)
    ]

    # ── Segment boundaries and double-up slab walls ─────────────────────────
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

    # ── Pressure-point setup ────────────────────────────────────────────────
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

    # ── Matrix build ────────────────────────────────────────────────────────
    pkernel = buildmatrix(
        lona, lata, lonb, latb,
        gam, alpha, lono, lato,
        iwall, idl, idr, n_segs, num_segs,
        coeff1, coeff2, coefftr1, coefftr2,
        ndomain, _EPSLIT, dsegtr,
        _RAD_KM, _ALITH, _AH1, _EPS_FACT,
    )

    # ── RHS vector ──────────────────────────────────────────────────────────
    vector = buildvector(
        iwall, alpha, ndomain, idl, idr,
        vtopl, vtopr, vbotl, vbotr, vt,
        n_segs, num_segs,
        flux_slab, _FLUX_VEL_CONST,
        flux_width, polarity, flux_alpha,
        no_flux_for_slabtails,
        _RAD_KM, _ALITH, _AH1,
    )

    # ── Inversion ───────────────────────────────────────────────────────────
    pcoeff = inv(pkernel).dot(vector)

    # ── Grid output ─────────────────────────────────────────────────────────
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
     _polygons) = outputgrids(
        grid_spacing,
        lona, lata, lonb, latb,
        lono, lato, iwall, gam, alpha,
        amu, _AH1 - _ALITH,
        n_segs, num_segs, pcoeff,
        _RAD_KM, domain_bounds, bound_ind,
        pole_top_lon, pole_top_lat, pole_top_rate,
        vt_ew, vt_ns, _ALITH, _PRESS_DEPTH,
        coefftr1, coefftr2,
        pole_bott_lon, pole_bott_lat, pole_bott_rate,
        rigid_vew, rigid_vns, _AH1, ndomain,
    )

    # ── DP output ───────────────────────────────────────────────────────────
    DP = outputDP(
        lona, lata, lonb, latb,
        lono, lato, iwall, gam, alpha,
        n_segs, num_segs, pcoeff, _RAD_KM,
        lon_subslab, lat_subslab,
        lon_wedge, lat_wedge,
        polarity, vtopl, vtopr, vt, _DIP_DEPTH,
    )

    # ── Serialise results ───────────────────────────────────────────────────
    lons_1d = lons_grd[0, :].tolist()
    lats_1d = lats_grd[:, 0].tolist()

    dp_list = [
        {
            'lon':    _to_float(DP[i, 0]),
            'lat':    _to_float(DP[i, 1]),
            'dp_MPa': _to_float(DP[i, 4]),
        }
        for i in range(len(DP))
        if iwall[i] == 1
    ]

    # Subsample velocity field to ~20×10 arrows for visualisation.
    step_lon = max(1, len(lons_1d) // 20)
    step_lat = max(1, len(lats_1d) // 10)
    vel_list = [
        {
            'lon': lons_1d[j],
            'lat': lats_1d[i],
            'vew': _to_float(vel_ew[i, j]),
            'vns': _to_float(vel_ns[i, j]),
        }
        for i in range(0, len(lats_1d), step_lat)
        for j in range(0, len(lons_1d), step_lon)
    ]

    return {
        'lons':       lons_1d,
        'lats':       lats_1d,
        'pressure':   P_out.tolist(),   # 2-D list [nlat][nlon], MPa
        'dp':         dp_list,
        'velocity':   vel_list,
        'boundaries': orig_bounds,
        'timing_s':   time.perf_counter() - t0,
    }
