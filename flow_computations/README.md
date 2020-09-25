### This directory contains scripts necessary to compute upper mantle pressure and flow driven by subduction in a spherical shell (Holt and Royden, 2020, G-cubed)

### RUNNING PRESSURE COMPUTATIONS

The main driver script is global_pressure_withPressurePlot.py, which computes pressure for an input plate and slab model (location in inputs/) and plots the pressure field (into plots/). 

Inside inputs/ are the slab geometries (see end of README for file structure details).  For example, the reference model geometry is 'Slab2.0FinalNo_JapTail_nnr_FS' and so there are two files corresponding to this model. Subbon_Slab2.0Final_NoJapTail_nnr_FS.inp, which contains the plate boundary geometry (and slab wall velocities) and Subbon_Slab2.0Final_NoJapTail_nnr_FS.inp, which contains the plate domain boundary conditions (e.g. the plate velocities). One can run "./plot_boundary_locations.py Slab2.0FinalNo_JapTail_nnr_FS" to check/plot the plate boundaries of this particular file. In addition to these geometry files, a Subgrd* file is required, which contains the boundary segment lengths and the interval over which to grid the pressure field for output. See the last section of this README for specific details about the structure of these three types of input files. 

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

### INPUT FILE STRUCTURES

Running a model requires a plate boundary file (e.g. Subbon_Slab2.0Final_NoJapTail_nnr_FS.inp), a plate domain file (Subfil_Slab2.0Final_NoJapTail_nnr_FS.inp) and a grid file (e.g. Subgrd.inp). Within the inputs/ directory are a range of input files to produce both Earth-like models (e.g. paper reference model) and highly idealized models (e.g. paper Figures 2, 3) Describe below is the structure of each of these files.

Boundary file (e.g. Subbon_Slab2.0Final_NoJapTail_nnr_FS.inp):

	Each row corresponds to a plate boundary segment, and the column values corresponid to the following;

	col. 1 = What type of boundary? 0 = plate 'edge' (i.e. no wall), 1 = subudction zone (i.e. slab wall)
	col. 2 = longitude of end point #1
	col. 3 = latitude of end point #1
	col. 4 = longitude of end point #2
	col. 5 = latitude of end point #2
	col. 6 = The indice of the domain (see corresponding plate domain file) that is to the left of the boundary when looking from point #1 to point #2
	col. 7 = The indice of the domain (see corresponding plate domain file) that is to the right of the boundary when looking from point #1 to point #2
	(Additional columns only needed if boundary is a slab wall, i.e. col. 1 = 1)
	col. 8 = East-west velocity of the slab wall, mm/yr, where East is positive 
	col. 9 = North-south velocity of the slab wall, mm/yr, where North is positive 
	col. 10 = Dip direction of the slab whenlooking from point #1 to point #2. r = right, l = left.
	col. 11 = Subduction index (the number of the individual subduction zone - e.g. the Aleutian subd. zone - that the segment is a part of. Only used for subsequent dip angle post-processing analysis)
	col. 8 (edge boundary) or col. 12 (wall boundary) = index of segment (used in the Plate domain file, see below)


	One can run "./plot_boundary_locations.py Slab2.0FinalNo_JapTail_nnr_FS nnr" to check/plot the plate boundaries of this particular file. Furthermore, this script interpolates trench velocities (here, in the no-net rotation ref. frame) to the segment lats/lons, so that columns 8 and 9 (trench velocity) can be filled in for slab walls.

Plate domain file (e.g. Subfil_Slab2.0Final_NoJapTail_nnr_FS.inp):

	This file specific the nature of the domains that the plate boundaries encircle. Each row corresponds to a distinct 'domain', with the following properties.

	col. 1 = The code for the type of domain (i.e. the type of basal and upper boundary conditions). See below for full list of options.
	col. 2 = If upper boundary velocity prescribed: Longitude of Euler pole
	col. 3 = If upper boundary velocity prescribed: Latiude of Euler pole
	col. 4 = If upper boundary velocity prescribed: Rotation rate of Euler pole [degrees/Myr]
	col. 5 = If lower boundary velocity prescribed (e.g. slab tail): Longitude of Euler pole
	col. 6 = If lower boundary velocity prescribed (e.g. slab tail): Latiude of Euler pole
	col. 7 = If lower boundary velocity prescribed (e.g. slab tail): Rotation rate of Euler pole [degrees/Myr]
	col. 8 = Index of domain (i.e. what is specified in columns 6 and 7 of boundary file) 
	col. 9 = Comma-separated list of the indices of the boundary segments that encircle the domain. These can be listed in either clock- or anti-clockwise direction, but need to be in the correct order.
	col. 10 = Name of the domain, e.g. India. (Just for book-keeping)

	Different domain options (i.e. what is entered into col. 1):

	100 = fixed upper plate, free slip base (i.e. asthenosphere thickness = 580 km)
	200 = fixed upper plate, fixed velocity slab tail (i.e. asth. thickness = 500 km)
	300 = fixed upper plate, fixed velocity base (i.e. asth. thickness = 580 km)
	400 = fixed upper plate, free slip slab tail (i.e. asth. thickness = 500 km)

Grid file (e.g. Subgrid.inp):

	This file should not change between model runs / plate geometries and just contains a few parameters:

	row 1, col 1 = The lat/lon interval used to output global pressure/velocity grids for storing and/or plotting. Note that it takes a long time to output over a fine grid interval (e.g. < 2 degrees)
	row 1, col 2 = The lat/lon interval used to output a pressure/velocity profile (can be much lower than that specified for the full, global grid.)
	row 2, col 1 = Length of slab wall segments, km (i.e. how the wall plate boundaries specified in the Plate boundary file are broken up)
	row 2, col 2 = Length of plate edge segments, km (i.e. how the wall plate boundaries specified in the Plate boundary file are broken up)





