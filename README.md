# Analytical mantle flow modeling

Code and data accompanying Holt and Royden (2020, *Geochemistry, Geophysics, Geosystems*):

- Holt, A. F., and Royden, L. H. (2020), *Subduction dynamics and mantle pressure: (ii) Towards a Global Understanding of Slab Dip and Upper Mantle Circulation*, doi:10.1029/2019GC008771
- Royden, L. H., and Holt, A. F. (2020), *Subduction dynamics and mantle pressure: (i) An Analytical Framework Relating Subduction Geometry, Plate Motion, and Asthenospheric Pressure*, doi:10.1029/2020GC009032

## What this repository contains

- `dip_observations/`: workflows and Slab2.0 inputs used to build the observed slab-dip catalogues.
- `flow_computations/`: analytical mantle pressure solver, dip prediction workflow, model inputs, and plotting scripts.
- `Holt-Royden_2020_G3.pdf`: manuscript PDF.

## Scientific workflow in this repo

1. Derive observed dips from Slab2.0 surfaces (`dip_observations/`).
2. Build/choose a global plate-slab geometry (`flow_computations/inputs/`).
3. Solve for mantle pressure and slab pressure discontinuity (`global_pressure_withPressurePlot.py`).
4. Convert pressure discontinuity + plate age into predicted slab dips and compare to observations (`plot_DipComparison_varyDPfactor.py`).

## Reproducibility status (audited February 15, 2026)

### What is reproducible immediately

The repo already includes key generated products, including:

- Observed dip catalogue: `dip_observations/dip_catalogues/Slab2_const-depth/AllDips.txt`
- Model outputs: `flow_computations/text_files/*/DP.txt` and `flow_computations/text_files/*/Pcoeff.txt`
- Main figures/plots in `flow_computations/plots/` and `dip_observations/plots/`

So you can inspect and analyze results directly without rerunning heavy computations.

### What currently blocks full rerun on a clean modern machine

The scripts are now Python 3 syntax compatible, but full reruns are still not plug-and-play on a clean machine because:

- several legacy geoscience dependencies are required (Basemap, SciPy stack, NetCDF tools, GMT, etc.);
- some platform-specific solver/tool issues may still require local tuning.

### Fix applied in this repository

To improve rerun portability, `flow_computations/plot_DipComparison_varyDPfactor.py` was updated to:

- use a repository-relative default observed dip file (`dip_observations/dip_catalogues/Slab2_const-depth/AllDips.txt`) instead of a machine-specific absolute path;
- honor the grid-file argument from the command line, instead of overriding it internally.

## Suggested environment for full reruns

For strict reproducibility of original behavior, use Python 3 with the following geoscience stack:

- Python 3.10+ (scripts migrated to Python 3 syntax)
- `numpy`, `scipy`, `matplotlib`, `basemap`, `netCDF4`, `xarray`, `shapely`, `geographiclib`, `pandas`, `geopy`, `gpxpy`
- GMT (for dip extraction scripts using `grdtrack`)
- ImageMagick `convert` (used by plotting scripts to make PNG outputs)
- `gawk`, standard Unix tools

A pinned Conda environment is provided at `environment.yml`.

Create and activate it:

```bash
conda env create -f environment.yml
conda activate mantle-flow-modeling
python -m ipykernel install --user --name mantle-flow-modeling --display-name "Python (mantle-flow-modeling)"
```

Launch Jupyter:

```bash
jupyter lab
```

## Main driver (Python 3)

Run from repository root:

```bash
./run_example_workflow.py --mode quick
```

`quick` mode validates that key precomputed outputs and inputs are present.

To rerun the full paper-style modeling chain (requires dependencies above):

```bash
./run_example_workflow.py --mode full
```

## Notebook workflow (simple -> global)

Use:

- `notebooks/workflow_examples.ipynb`

The notebook:

1. Runs/loads two simple idealized examples (`LargeSP_RetreatingTrench`, `LargeSP_RetreatingTrenchSlabGap`).
2. Plots and summarizes DP outputs.
3. Builds up to a global case (`Slab2.0Final_NoJapTail_nnr_FS`) with pressure, age, and dip-comparison steps.

Set `RUN_HEAVY = False` in the notebook to use precomputed shipped outputs, or `RUN_HEAVY = True` to rerun scripts.

## Notebook workflow (publication-style reproduction)

Use:

- `notebooks/publication_style_reproduction.ipynb`

This notebook is figure-oriented and deterministic:

1. Validates a manifest of key figure/data assets.
2. Computes reproducibility summary statistics from shipped outputs.
3. Displays publication-style figure panels directly.
4. Optionally reruns core generation scripts if `RUN_REGENERATE = True`.

## Running the core computations

Run from `flow_computations/`.

### 1) Pressure/discontinuity computation

```bash
python global_pressure_withPressurePlot.py Slab2.0Final_NoJapTailNoPhil_nnr_FS 3.0e20 2 550e3 0 1 Subgrd.inp 4.0e20
```

Arguments:

1. Plate model name (for `inputs/Subbon_<model>.inp` and `inputs/Subfil_<model>.inp`)
2. Asthenosphere viscosity used for inversion (Pa s)
3. Slab flux mode (`0` none, `1` constant 50 mm/yr, `2` convergence-scaled)
4. Flux width (m)
5. Flux partition (`0` overriding-plate side, `1` subducting-plate side)
6. `no_flux_for_slabtails` (`0` allow flux at slab tails, `1` suppress)
7. Grid file in `inputs/` (for example `Subgrd.inp`)
8. Viscosity used for plotting scale (Pa s)

### 2) Subducting-plate ages (for dip prediction)

```bash
./get_SPages.py Slab2.0Final_NoJapTail_nnr_FS
```

### 3) Dip comparison and viscosity-factor search

```bash
python plot_DipComparison_varyDPfactor.py Slab2.0Final_NoJapTail_nnr_FS 3.0e20 2 500e3 0 1 Subgrd.inp 2 12 4 5
```

If your observed dip file is elsewhere, set:

```bash
export DIPS_OBS_TXT=/absolute/path/to/AllDips.txt
```

## Running observed-dip extraction

Run from `dip_observations/`:

```bash
./many_dip_extractions.sh
./many_dip_extractions_2Contours.sh
```

Outputs are written to `dip_catalogues/Slab2_const-depth/` and `dip_catalogues/Slab2_two-depths/`.

## Notes

- The repository includes many precomputed outputs from the publication workflow.
- If you want, a next step is to add a pinned Conda environment (or Docker image) to make full reruns one-command reproducible.
