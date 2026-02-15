# Analytical mantle flow modeling

Code and data accompanying:

- Holt, A. F., and Royden, L. H. (2020), *Subduction dynamics and mantle pressure: (ii) Towards a Global Understanding of Slab Dip and Upper Mantle Circulation*, doi:10.1029/2019GC008771
- Royden, L. H., and Holt, A. F. (2020), *Subduction dynamics and mantle pressure: (i) An Analytical Framework Relating Subduction Geometry, Plate Motion, and Asthenospheric Pressure*, doi:10.1029/2020GC009032

## Repository layout

- `flow_computations/`: pressure solver, dip-prediction scripts, model inputs, plots, and text outputs.
- `dip_observations/`: Slab2-based observed dip catalogues and extraction scripts.
- `notebooks/workflow_examples.ipynb`: educational notebook from simple to global examples.
- `environment.yml`: pinned Conda environment.
- `Holt-Royden_2020_G3.pdf`: manuscript PDF.

## Quick start

```bash
conda env create -f environment.yml
conda activate mantle-flow-modeling
python -m ipykernel install --user --name mantle-flow-modeling --display-name "Python (mantle-flow-modeling)"
jupyter lab
```

Open `notebooks/workflow_examples.ipynb` and select the `Python (mantle-flow-modeling)` kernel.

## Main workflow options

### Notebook workflow (recommended)

Use `notebooks/workflow_examples.ipynb`.

- `RUN_HEAVY = False`: use shipped precomputed outputs.
- `RUN_HEAVY = True`: rerun model scripts and regenerate outputs.

### CLI driver

From repository root:

```bash
./run_example_workflow.py --mode quick
./run_example_workflow.py --mode full
```

- `quick`: validates required precomputed artifacts.
- `full`: runs pressure, age, and dip-comparison pipeline.

## Core scripts (manual run)

Run from `flow_computations/`.

### 1) Pressure and DP computation

```bash
python global_pressure_withPressurePlot.py Slab2.0Final_NoJapTailNoPhil_nnr_FS 3.0e20 2 550e3 0 1 Subgrd.inp 4.0e20
```

Arguments:

1. Plate model name (for `inputs/Subbon_<model>.inp` and `inputs/Subfil_<model>.inp`)
2. Asthenosphere viscosity for inversion (Pa s)
3. Slab flux mode (`0` none, `1` constant 50 mm/yr, `2` convergence-scaled)
4. Flux width (m)
5. Flux partition (`0` overriding-plate side, `1` subducting-plate side)
6. `no_flux_for_slabtails` (`0` allow, `1` suppress)
7. Grid file in `inputs/` (e.g., `Subgrd.inp`)
8. Plot viscosity scale (Pa s)

### 2) Subducting-plate ages

```bash
./get_SPages.py Slab2.0Final_NoJapTail_nnr_FS
```

### 3) Modeled vs observed dip comparison

```bash
python plot_DipComparison_varyDPfactor.py Slab2.0Final_NoJapTail_nnr_FS 3.0e20 2 500e3 0 1 Subgrd.inp 2 12 4 5
```

By default this uses:

`dip_observations/dip_catalogues/Slab2_const-depth/AllDips.txt`

Override with:

```bash
export DIPS_OBS_TXT=/absolute/path/to/AllDips.txt
```

## Observed dip extraction

Run from `dip_observations/`:

```bash
./many_dip_extractions.sh
./many_dip_extractions_2Contours.sh
```

Outputs are written to:

- `dip_catalogues/Slab2_const-depth/`
- `dip_catalogues/Slab2_two-depths/`

## Notes

- The repository already includes key precomputed outputs in `flow_computations/text_files/` and `flow_computations/plots/`.
- Plot conversion uses ImageMagick (`magick` when available, fallback to `convert`).
