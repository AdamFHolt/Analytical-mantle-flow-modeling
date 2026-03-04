# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment setup

```bash
conda env create -f environment.yml
conda activate mantle-flow-modeling
python -m ipykernel install --user --name mantle-flow-modeling --display-name "Python (mantle-flow-modeling)"
```

Key dependencies: Python 3.10, numpy, scipy, matplotlib/basemap, xarray, shapely, geographiclib, GMT, ImageMagick.

## Running the workflow

### Notebook (recommended)
```bash
jupyter lab
# Open notebooks/workflow_examples.ipynb with "Python (mantle-flow-modeling)" kernel
# RUN_HEAVY = False → use precomputed outputs
# RUN_HEAVY = True  → rerun all scripts (uses Subgrd_Fast.inp)
```

### CLI driver (from repo root)
```bash
./run_example_workflow.py --mode quick   # validate precomputed artifacts only
./run_example_workflow.py --mode full    # run full pressure + age + dip pipeline
```

### Core scripts (run from `flow_computations/`)

**1. Pressure field computation:**
```bash
python global_pressure_withPressurePlot.py Slab2.0Final_NoJapTailNoPhil_nnr_FS 3.0e20 2 550e3 0 1 Subgrd_Fast.inp 4.0e20
```

**2. Subducting plate ages:**
```bash
./get_SPages.py Slab2.0Final_NoJapTail_nnr_FS
```

**3. Modeled vs observed dip comparison:**
```bash
python plot_DipComparison_varyDPfactor.py Slab2.0Final_NoJapTail_nnr_FS 3.0e20 2 500e3 0 1 Subgrd_Fast.inp 2 12 4 5
```

Override default dip catalogue: `export DIPS_OBS_TXT=/absolute/path/to/AllDips.txt`

### Observed dip extraction (from `dip_observations/`)
```bash
./many_dip_extractions.sh
./many_dip_extractions_2Contours.sh
```

## Architecture

### Computational pipeline

`global_pressure_withPressurePlot.py` is the main solver. It:
1. Reads plate boundary geometry (`inputs/Subbon_<model>.inp`) and domain/BC file (`inputs/Subfil_<model>.inp`)
2. Segments boundaries into small elements using `functions.organizebounds()`
3. Solves for mantle pressure everywhere via a linear system (BEM-style): `buildmatrix()` → `buildvector()` → matrix inversion
4. Outputs `text_files/<run_name>/DP.txt` (pressure discontinuity at each slab wall) and `Pcoeff.txt`
5. Plots pressure/velocity field to `plots/pressure_fields/`

`plot_DipComparison_varyDPfactor.py` then combines `DP.txt` with plate ages (from `ages/`) and observed dips (from `dip_observations/dip_catalogues/`) to find the best-fit viscosity and produce comparison plots in `plots/dip_comparisons/`.

### Input file types (in `flow_computations/inputs/`)

- `Subbon_<model>.inp` — plate boundary segments (type, endpoints, domain indices, wall velocities, dip polarity, subduction index)
- `Subfil_<model>.inp` — plate domains (Euler pole parameters, domain type codes 100/200/300/400 controlling BCs, encircling segment list)
- `Subgrd.inp` / `Subgrd_Fast.inp` — output grid spacing and boundary segment lengths. `Subgrd_Fast.inp` uses coarser resolution for faster runs.

### Key modules

- `flow_computations/functions.py` — all core physics: `readdomains`, `readbounds`, `organizebounds`, `pressurepoints`, `buildmatrix`, `buildvector`, `outputDP`, `outputgrids`; also oceanic buoyancy calculations
- `flow_computations/plotting.py` — pressure/velocity plots via Basemap (non-interactive, `Agg` backend)
- `flow_computations/euler_pole.py` — Euler pole → velocity vector conversion
- `dip_observations/functions.py` — dip extraction utilities

### Output locations

- `flow_computations/text_files/<run_name>/DP.txt` — pressure discontinuity per slab segment
- `flow_computations/text_files/<run_name>/Pcoeff.txt` — pressure coefficients
- `flow_computations/ages/<model>.txt` — subducting plate ages
- `flow_computations/plots/pressure_fields/` — pressure field maps
- `flow_computations/plots/dip_comparisons/` — model vs observed dip plots

### Run name convention

Output directories and filenames encode all run parameters, e.g.:
`Slab2.0Final_NoJapTail_nnr_FS.3e+20_VcSlabFlux_width500000.0_alpha0.0.NoTailFlux`

### Plate models available

- `Slab2.0Final_NoJapTail_nnr_FS` — reference Earth model (no Japan slab tail, no-net-rotation frame)
- `Slab2.0Final_NoJapTailNoPhil_nnr_FS` — same but also excluding Philippine slab
- `LargeSP_RetreatingTrench`, `SmallSP_RetreatingTrench`, `LargeSP_RetreatingTrenchSlabGap` — idealized models for paper figures

### Plot conversion

Scripts call ImageMagick (`magick` when available, fallback to `convert`) to convert EPS/PS output to PNG.
