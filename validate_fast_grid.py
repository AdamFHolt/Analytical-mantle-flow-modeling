#!/usr/bin/env python3
"""
validate_fast_grid.py
Validates that _fast_outputgrids in webapp/compute.py produces results
numerically identical to the reference functions.outputgrids.

Run from the repo root:
    conda run -n mantle-flow-modeling python validate_fast_grid.py
"""
import sys, time
import pathlib
import numpy as np

ROOT = pathlib.Path(__file__).parent
FC   = ROOT / 'flow_computations'
sys.path.insert(0, str(FC))
sys.path.insert(0, str(ROOT / 'webapp'))

from functions import (
    readdomains, readgrid, readbounds, organizebounds,
    pressurepoints, buildmatrix, buildvector, outputgrids,
)
from numpy.linalg import inv
import compute as _cmod

# ── Parameters (match compute.py defaults) ────────────────────────────────────
MODEL        = 'Simple'
SPACING      = 10.0          # coarse grid — fast reference run
AMU          = 3e20
BASAL_BC     = 'free'
FLUX_SLAB    = 2
FLUX_WIDTH   = 500_000.0
FLUX_ALPHA   = 0.0
NO_FLUX_TAIL = 1
EULER_LAT, EULER_LON, EULER_RATE = 90.0, 180.0, 0.45
TRENCH_VEW,  TRENCH_VNS          = 0.0, 0.0
WIDTH_DEG,   LENGTH_DEG          = 45.0, 60.0

ALITH = 80.0e3
AH1   = 660.0e3
RAD_KM     = 6378.0
EPS_FACT   = 0.01
EPSDP_FACT = 2 * EPS_FACT
EPSLRC     = 100.0e3
EPSLIT     = 1.0e3
ISEG_MIN   = 0
DIP_DEPTH  = 330.0e3

# ── Build inputs ──────────────────────────────────────────────────────────────
print('Generating Simple model input files...')
_cmod.generate_simple_inp_files(WIDTH_DEG, LENGTH_DEG)

inp_dir   = FC / 'inputs'
f_bounds  = str(inp_dir / f'Subbon_{MODEL}.inp')
f_domains = str(inp_dir / f'Subfil_{MODEL}.inp')
f_grid    = str(inp_dir / 'Subgrd_Fast.inp')

coeff1   = (AH1 - ALITH)    ** 3 / (AMU * AH1)
coeff2   = (AH1 - 2*ALITH)  ** 3 / (AMU * AH1)
coefftr1 = (AH1 - ALITH)    ** 2 / AMU
coefftr2 = (AH1 - 2*ALITH)  ** 2 / AMU

(ndomain,
 pole_top_lon, pole_top_lat, pole_top_rate,
 pole_bott_lon, pole_bott_lat, pole_bott_rate,
 rigid_vew, rigid_vns,
 domain_bounds) = readdomains(f_domains)

if BASAL_BC == 'no_slip':
    ndomain = [300 if v == 100 else v for v in ndomain]

grid_spacing_default, _prof, dsegtr, dseged = readgrid(f_grid)

(num_bounds,
 iwall, lona, lata, lonb, latb, bound_ind,
 idl, idr, vt_ew, vt_ns, polarity, large_wall_inds) = readbounds(f_bounds)

# Apply Euler pole override
sp_dom = next((idl[i]-1 for i in range(num_bounds) if iwall[i] == 1), None)
if sp_dom is not None:
    pole_top_lat[sp_dom] = pole_bott_lat[sp_dom] = EULER_LAT
    pole_top_lon[sp_dom] = pole_bott_lon[sp_dom] = EULER_LON
    pole_top_rate[sp_dom] = pole_bott_rate[sp_dom] = EULER_RATE
for i in range(num_bounds):
    if iwall[i] == 1:
        vt_ew[i] = TRENCH_VEW
        vt_ns[i] = TRENCH_VNS

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
    RAD_KM, ISEG_MIN,
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
    ndomain, EPSLRC, RAD_KM, ALITH,
    0, polarity, EPSDP_FACT,
)

print(f'n_segs={n_segs}, num_segs={num_segs}, spacing={SPACING}°')
print('Building matrix and solving...')
t0 = time.perf_counter()
pkernel = buildmatrix(
    lona, lata, lonb, latb,
    gam, alpha, lono, lato,
    iwall, idl, idr, n_segs, num_segs,
    coeff1, coeff2, coefftr1, coefftr2,
    ndomain, EPSLIT, dsegtr,
    RAD_KM, ALITH, AH1, EPS_FACT,
)
vector = buildvector(
    iwall, alpha, ndomain, idl, idr,
    vtopl, vtopr, vbotl, vbotr, vt,
    n_segs, num_segs,
    FLUX_SLAB, 50.0,
    FLUX_WIDTH, polarity, FLUX_ALPHA,
    NO_FLUX_TAIL,
    RAD_KM, ALITH, AH1,
)
pcoeff = inv(pkernel).dot(vector)
print(f'  solve took {time.perf_counter()-t0:.2f}s')

# ── Common kwargs for both outputgrids calls ──────────────────────────────────
kwargs = dict(
    spacing=SPACING,
    lona=lona, lata=lata, lonb=lonb, latb=latb,
    lono=lono, lato=lato, iwall=iwall, gam=gam, alpha=alpha,
    amu=AMU, a=AH1-ALITH, n_segs=n_segs, num_segs=num_segs, pcoeff=pcoeff,
    rad_km=RAD_KM, domain_bounds=domain_bounds, bound_ind=bound_ind,
    pole_top_lon=pole_top_lon, pole_top_lat=pole_top_lat,
    pole_top_rate=pole_top_rate,
    vt_ew=vt_ew, vt_ns=vt_ns, alith=ALITH, press_depth=330e3,
    coefftr1=coefftr1, coefftr2=coefftr2,
    pole_bott_lon=pole_bott_lon, pole_bott_lat=pole_bott_lat,
    pole_bott_rate=pole_bott_rate,
    rigid_vew=rigid_vew, rigid_vns=rigid_vns,
    ah1=AH1, ndomain=ndomain,
)

# ── Run reference outputgrids ─────────────────────────────────────────────────
print('\nRunning reference outputgrids...')
t0 = time.perf_counter()
ref = outputgrids(**kwargs)
t_ref = time.perf_counter() - t0
print(f'  reference: {t_ref:.2f}s')

# ── Run _fast_outputgrids ─────────────────────────────────────────────────────
print('Running _fast_outputgrids...')
t0 = time.perf_counter()
fast = _cmod._fast_outputgrids(
    SPACING,
    lona, lata, lonb, latb,
    lono, lato, iwall, gam, alpha,
    AMU, AH1-ALITH,
    n_segs, num_segs, pcoeff,
    RAD_KM, domain_bounds, bound_ind,
    pole_top_lon, pole_top_lat, pole_top_rate,
    vt_ew, vt_ns, ALITH, 330e3,
    coefftr1, coefftr2,
    pole_bott_lon, pole_bott_lat, pole_bott_rate,
    rigid_vew, rigid_vns, AH1, ndomain,
)
t_fast = time.perf_counter() - t0
print(f'  fast:      {t_fast:.2f}s  ({t_ref/t_fast:.1f}× speedup)')

# ── Compare outputs ───────────────────────────────────────────────────────────
NAMES = [
    'P_grd', 'Pwall_grd', 'Pedge_grd',
    'dPdlon_grd', 'dPdlat_grd',
    'polygon_points',
    'plate_vel_ew', 'plate_vel_ns',
    'trench_vels',
    'avgvel_asthen_ew', 'avgvel_asthen_ns',
    'pdrivenvel_wall_ew', 'pdrivenvel_wall_ns',
    'pdrivenvel_edge_ew', 'pdrivenvel_edge_ns',
    'lon_grd', 'lat_grd',
    'polygons',
]

print('\n── Field comparison ─────────────────────────────────────────────────')
print(f'{"Field":<24}  {"max|Δ|":>12}  {"max|ref|":>12}  {"rel err":>10}  status')
print('-' * 75)

all_ok = True
for name, r, f in zip(NAMES, ref, fast):
    if name == 'polygons':
        continue  # Shapely objects — skip
    r_arr = np.asarray(r, dtype=float)
    f_arr = np.asarray(f, dtype=float)
    if r_arr.shape != f_arr.shape:
        print(f'{name:<24}  shape mismatch: ref={r_arr.shape} fast={f_arr.shape}')
        all_ok = False
        continue
    mask = np.isfinite(r_arr) & np.isfinite(f_arr)
    if not mask.any():
        print(f'{name:<24}  (all NaN/Inf — skipped)')
        continue
    abs_err  = np.abs(r_arr[mask] - f_arr[mask]).max()
    ref_max  = np.abs(r_arr[mask]).max()
    rel_err  = abs_err / ref_max if ref_max > 0 else 0.0
    ok = rel_err < 1e-4 or abs_err < 1e-8
    status = 'OK' if ok else 'FAIL'
    if not ok:
        all_ok = False
    print(f'{name:<24}  {abs_err:12.4e}  {ref_max:12.4e}  {rel_err:10.2e}  {status}')

print()
if all_ok:
    print('VALIDATION PASSED — _fast_outputgrids matches reference.')
else:
    print('VALIDATION FAILED — see FAIL rows above.')
