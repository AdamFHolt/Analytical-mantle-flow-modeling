### This directory contains scripts necessary to compute Earth dip angles using Slab2.0 (Hayes et al. 2018).

Requirements: GMT, Python, plotting and geographic Python modules (see first few lines of .py scripts).

In order to re-compute the dip angle 'observations' in Holt and Royden (2019), run the following wrappers:

./many_dip_extractions_2Contours.sh  
./many_dip_extractions.sh

-These script place dip angle catalogues into dip_catalogues/Slab2_two-depths and dip_catalogues/Slab2_const-depth, respectively.
-The dips used in the main manuscript are dip_catalogues/Slab2_const-depth/AllDips.txt (dips extracted from one depth on Slab2 slab surface), ...
... and these are compared to dip_catalogues/Slab2_two-depths/AllDips.txt (dips computed with two contours - separated by 50km depth, on Slab2 slab surface).
-For these two dip calculation methods, there is an AllDips.txt and an AllDips.full.txt catalogue. The former computed dips every 250 km along a trench ...
... (and is used for the comparison) and the latter is at a much smaller interval that is specified as a longitude or latitude (and so is sometimes unequal between trenches).
-In order to make a comparison plot of the two dip catalogues, and the dips of Lallemand, run 'compare_all_dips.py' to produce Figure S3.
 
