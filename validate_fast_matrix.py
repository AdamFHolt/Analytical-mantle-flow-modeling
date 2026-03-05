#!/usr/bin/env python3
"""
validate_fast_matrix.py
Validates that _fast_buildmatrix in webapp/compute.py produces a matrix
numerically identical to the reference functions.buildmatrix.

Run from the repo root:
    conda run -n mantle-flow-modeling python validate_fast_matrix.py
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
    pressurepoints, buildmatrix,
)
import compute as _cmod

# ── Parameters (match compute.py defaults) ────────────────────────────────────
MODEL        = 'Simple'
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

print(f'n_segs={n_segs}, num_segs={num_segs}')

# Common kwargs
bm_args = (
    lona, lata, lonb, latb,
    gam, alpha, lono, lato,
    iwall, idl, idr, n_segs, num_segs,
    coeff1, coeff2, coefftr1, coefftr2,
    ndomain, EPSLIT, dsegtr,
    RAD_KM, ALITH, AH1, EPS_FACT,
)

# ── Reference buildmatrix ─────────────────────────────────────────────────────
print('\nRunning reference buildmatrix...')
t0 = time.perf_counter()
ref = buildmatrix(*bm_args)
t_ref = time.perf_counter() - t0
print(f'  reference: {t_ref:.2f}s')

# ── Fast buildmatrix ──────────────────────────────────────────────────────────
print('Running _fast_buildmatrix...')
t0 = time.perf_counter()
fast = _cmod._fast_buildmatrix(*bm_args)
t_fast = time.perf_counter() - t0
print(f'  fast:      {t_fast:.2f}s  ({t_ref/t_fast:.1f}x speedup)')

# ── Compare ───────────────────────────────────────────────────────────────────
ref  = np.asarray(ref,  dtype=float)
fast = np.asarray(fast, dtype=float)

mask    = np.isfinite(ref) & np.isfinite(fast)
abs_err = np.abs(ref[mask] - fast[mask]).max()
ref_max = np.abs(ref[mask]).max()
rel_err = abs_err / ref_max if ref_max > 0 else 0.0

print(f'\nMatrix shape: {ref.shape}')
print(f'max |Δ|:  {abs_err:.4e}')
print(f'max |ref|:{ref_max:.4e}')
print(f'rel err:  {rel_err:.2e}')

ok = rel_err < 1e-4 or abs_err < 1e-8
print()
if ok:
    print('VALIDATION PASSED — _fast_buildmatrix matches reference.')
else:
    print('VALIDATION FAILED')
    # Show worst offenders
    diff = np.abs(ref - fast)
    idx  = np.unravel_index(np.argmax(diff), diff.shape)
    print(f'  worst cell: ({idx[0]}, {idx[1]})  ref={ref[idx]:.6e}  fast={fast[idx]:.6e}')
