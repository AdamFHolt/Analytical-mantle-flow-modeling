### This directory contains scripts necessary to compute upper mantle pressure and flow driven by subduction in a spherical shell (Holt and Royden, 2019)

### RUNNING PRESSURE COMPUTATIONS

The main driver script is global_pressure_withPressurePlot.py, which computes pressure for an input plate and slab model (location in inputs/) and plots the pressure field (into plots/). 

Inside inputs/ are the slab geometries.  For example, the reference model geometry is 'Slab2.0FinalNo_JapTail_nnr_FS' and so there are two files corresponding to this model. Subbon_Slab2.0Final_NoJapTail_nnr_FS.inp, which contains the plate boundary geometry (and slab wall velocities) and Subbon_Slab2.0Final_NoJapTail_nnr_FS.inp, which contains the plate domain boundary conditions (e.g. the plate velocities). One can run "./plot_boundary_locations.py Slab2.0FinalNo_JapTail_nnr_FS" to check/plot the plate boundaries of this particular file.

In addition to these geometry files, a Subgrd* file is required, which contains the boundary segment lengths and the interval over which to grid the pressure field for output. 

As an example, to run a model using this reference geometry, make a plot of the pressure and velocity field, and calculate the pressure discontinuity (DP) at all slab walls, one would execute:

python global_pressure_withPressurePlot.py Slab2.0Final_NoJapTailNoPhil_nnr_FS 3.0e20 2 550e3 0 1 Subgrd.inp 4.0e20

Therefore, there are 8 parameters that correspond to the following:

$1 = The plate model. The code for the following two files in inputs/, Subfil_{plate model}.inp and Subbon_{plate_model}.inp
$2 = Asthenosphere viscosity. This viscosity is used to compute the pressure field and resulting DP values.
$3 = Whether to flux material into the lower mantle: 0 = no fluxing, 1 = flux material with const. velocity at all slabs (50 mm/yr), 2 = flux material at convergence velocity.
$4 = Flux width: total width of the layers fluxing material into/out of the lower mantle. (Positive value = downward flux, as considered in study).
$5 = How flux is partitioned between upper and lower plate side of slab: 0 = all flux on overriding plate side, 1 = all on subducting plate side.
$6 = The name of the grid file located in inputs/, e.g. Subgrid.inp
$7 = Asthenospheric viscosity for plot. If you want to plot the pressure field using a different viscosity to that which you originally run the model (and output the DP for).

Hence, the above example is a down-flux model with a flux with of 550 km, a flux velocity that is the convergence rate. The scripts outputs the model DP (here, for an asthenospheric viscosity of 3e20 Pas) into text_files/{plate model} and a pressure plot (scaled to an asthenospheric viscosity of 4e20 Pas) into plots/{plate model}

### ANALYZING DIP ANGLES

Using the outputted DP file from the previous step, text_files/{plate model}/DP.txt, one can use the plot_DipComparison_varyDPfactor.py script to make a plot that compares the model dips (calculated using DP and plate age) to observed dips (in ../dip_observations/dip_catalogues/Slab2_const-depth). Prior to running this script, we need to run get_SPages.py - this takes the input boundary file (inputs/Subbon_{plate_model}.inp) and calculates the age of the subducting plate at all of the subduction segments. These ages are used to compute plate buoyancy which, when combined with model DP, gives model dip. For, this reference geometry, ages are extracted as:

./get_SPages.py Slab2.0Final_NoJapTail_nnr_FS

This uses the Muller et al. (2016) grid to extract ages, and places an age file into ages/Slab2.0Final_NoJapTail_nnr_FS.txt
For the particular model ran above, the dip comparison can then be ran as:

python plot_DipComparison_varyDPfactor.py Slab2.0Final_NoJapTail_nnr_FS 3.0e20 2 500e3 0 1 Subgrd.inp 2 12 4 5

Note firstly that the first 7 input parameters need to be the same as those used to run the initial model (see above). The other four parameters correspond to the following:

$8  = subduction index (column 11 of inputs/Subbon_{plate_model}.inp file) of the Pacific subduction zone. (So that the subduction zone can be split into three sections for averaging)
$9  = subduction index of the S. America subduction zone. (So that the subduction zone can be split into two sections for averaging) 
$10 = subduction index of any complex subduction zone we want to exclude from the analysis (here, index 4 and 5 correspond to the Solomon and New Hebrides trenches)
$11 = subduction index of any other complex subduction zone we want to exclude from the analysis (here, index 4 and 5 correspond to the Solomon and New Hebrides trenches)

Using the DP file initially outputted by the model run, this script loops through a range of viscosities, searches for the one that gives the best fit (lowest RMS misfit) between model and observed dips, and places a observed-modeled dip comparison plot for this viscosity into plots/dip_comparisons 


















