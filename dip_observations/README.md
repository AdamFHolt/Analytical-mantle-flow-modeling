### This directory contains scripts necessary to compute Earth dip angles using Slab2.0 (Hayes et al. 2018), contained here in Slab2Distribute_Mar2018

Requirements: GMT, Python, plotting and geographic Python modules (see first few lines of .py scripts).

In order to re-compute the dip angle 'observations' in Holt and Royden (2020), run the following wrapper scripts:

	./many_dip_extractions_2Contours.sh
	./many_dip_extractions.sh

- These scripts place dip angle catalogues into dip_catalogues/Slab2_two-depths and dip_catalogues/Slab2_const-depth, respectively.

- The dips used in the main manuscript are dip_catalogues/Slab2_const-depth/AllDips.txt (dips extracted from one depth on Slab2 slab surface), and these are compared to dip_catalogues/Slab2_two-depths/AllDips.txt (dips computed with two contours - separated by 50km depth, on Slab2 slab surface).

- For these two dip calculation methods, there is an AllDips.txt and an AllDips.full.txt catalogue. The former computed dips at a distance interval of 250 km along each trench (and is used for the comparison and in the main paper) and the latter is at a much smaller interval that is specified as a constant longitude or latitude interval (and so the corresponding distance is often unequal between trenches - hence, these files are just for testing purposes).

- In order to make a comparison plot of the two dip catalogues, and the dips of Lallemand et al.(2005), run 'compare_all_dips.py' to produce Figure S3.
 
