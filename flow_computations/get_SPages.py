#!/usr/bin/python
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from functions import haversine, pressurepoints, project_to_point
from functions import readbounds, organizebounds, readgrid, readdomains
from netCDF4 import Dataset
from plotting import shiftedColorMap
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.cm as cm
import subprocess
import sys, os, math

plt.ioff()

plates=str(sys.argv[1])	         
rad_km = 6378.;

#### MAKE SURE FOLLOWING PARAMS CONSISTENT WITH GLOBAL_PRESSURE*py SCRIPT ######
max_ocean_age = 80.
shift_edges=0   	
epslrc = 100.e3 	
epsdp_fact =  1  	# (doesn't get used here...)
rad_km = 6378.;
alith = 80.0e3; 
#################################################################################

infile_grid='inputs/Subgrd.inp' 
infile_bounds  =''.join(['inputs/Subbon_',plates,'.inp']);
infile_domains  =''.join(['inputs/Subfil_',plates,'.inp']);
age_save_name = ''.join(['ages/',str(plates),'.txt'])

ndomain, pole_top_lon, pole_top_lat, pole_top_rate, pole_bott_lon, pole_bott_lat, pole_bott_rate, rigid_vew, rigid_vns, domain_bounds  = readdomains(infile_domains)
grid_spacing, prof_spacing, dsegtr, dseged = readgrid(infile_grid) 
num_bounds,iwall,lona,lata,lonb,latb,bound_ind,idl,idr,vt_ew,vt_ns,polarity,large_wall_inds = readbounds(infile_bounds)

iseg_min=0
n_segs,num_segs,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,polarity,num_wall_segs = \
	organizebounds(num_bounds,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,dsegtr,dseged,polarity,rad_km,iseg_min)
print "input done.\n---"

print "setting up pressure inversion points..."
lono,lato,gam,alpha,vtopl,vtopr,vbotl,vbotr,vt,lon_subslab,lat_subslab,lon_wedge,lat_wedge =  \
	pressurepoints(lona,lata,lonb,latb,vt_ew,vt_ns,iwall,idl,idr,n_segs,pole_top_lon,pole_top_lat,pole_top_rate,pole_bott_lon, \
		pole_bott_lat,pole_bott_rate,rigid_vew,rigid_vns,ndomain,epslrc,rad_km,alith,shift_edges,polarity,epsdp_fact)
print "pressure points set up.\n---"

# get subducting plate ages #############################################################

#read grid for plotting
age_grd_name='ages/Muller_etal_2016_AREPS_Agegrids/just-present/agegrid_0.low_res.nc'
age_grd = Dataset(age_grd_name, mode='r')
lons = age_grd.variables['x'][:]
lats = age_grd.variables['y'][:]
ages = np.array(age_grd.variables['z'][:])
age_grd.close()

age_txt_name='ages/Muller_etal_2016_AREPS_Agegrids/just-present/agegrid_0.low_res.txt'
age_txt = np.loadtxt(age_txt_name)

num_DPs = len(lono)
age_len_ratio = 225.e3/(rad_km * 1e3)
sp_ages = np.zeros((num_DPs,3))
for i in range(num_DPs):

	# get age coordinates
	if iwall[i] == 1:
		if polarity[i] == 1: # dipping to left
			age_azim = alpha[i] + (np.pi/2.)
		else:				 # dipping to right 
			age_azim = alpha[i] - (np.pi/2.)
		lon_age_rad , lat_age_rad = project_to_point(np.deg2rad(lono[i]), np.deg2rad(lato[i]), age_azim, age_len_ratio)
		if lon_age_rad < 0.:
			lon_age_rad = (2.*np.pi) + lon_age_rad
		sp_ages[i,0] = math.degrees(lon_age_rad)
		sp_ages[i,1] = math.degrees(lat_age_rad)

		dist_max = 1.e9
		for j in range(0,len(age_txt)):
			dist_pt = haversine(sp_ages[i,0], sp_ages[i,1], age_txt[j,0], age_txt[j,1],6378.)
			if dist_pt < dist_max:
				dist_max = dist_pt
				closest_ind = j
		print i
	
		sp_ages[i,2] = age_txt[closest_ind,2]

print "getting rid of nans..."
for i in range(len(sp_ages)):

	if np.isnan(sp_ages[i,2]):

		dist_max = 1.e9
		for k in range(0,len(sp_ages)):
			dist_pt = haversine(sp_ages[i,0], sp_ages[i,1], sp_ages[k,0], sp_ages[k,1],6378.)
			if dist_pt < dist_max and np.isnan(sp_ages[k,2]) == 0 and iwall[k] == 1:
				dist_max = dist_pt
				closest_ind = k

		sp_ages[i,2] = sp_ages[closest_ind,2]

print "got sp ages, saving..."
np.savetxt(age_save_name,sp_ages,fmt='%.4f')	

print "making a test plot..."
fig = plt.figure()

ax = fig.add_subplot(221)
map = Basemap(projection='hammer',lon_0=180,resolution='l')
map.drawmeridians(np.arange(0,360,45),linewidth=0.3)
map.drawparallels(np.arange(-90,90,45),linewidth=0.3)	
# segments
for i in range (num_DPs):
	if iwall[i] == 1:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown',zorder=0)	
	else:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.3,color='black',zorder=0)

# plot age grid
lon, lat = np.meshgrid(lons, lats)
xi, yi = map(lon, lat)
age_plot = map.pcolor(xi,yi,ages,cmap='RdBu',vmin=0,vmax=160,zorder=1)

# plot age points
for j in range(num_DPs):
	xo, yo = map(sp_ages[j,0],sp_ages[j,1])
	if np.isnan(sp_ages[j,2]) == 1 and iwall[j]==1:
		map.scatter(xo, yo, c='r',s=2.5,lw=0,zorder=2)
	elif iwall[j] == 1:
		ages_plot = map.scatter(xo, yo, c=sp_ages[j,2],s=5,cmap='RdBu',vmin=0,vmax=160,edgecolor='black',linewidth=0.2,zorder=2)
cbar = map.colorbar(ages_plot,location='bottom',pad="5%",size="4%")
cbar.set_ticks(np.array([0,40,80,120,160,180]))
cbar.set_label('age, Ma',size=9)
cbar.ax.tick_params(labelsize=8)

test_plot_name=''.join(['plots/ages_and_boundary_tests/ages/',str(plates),'.pdf'])
plt.savefig(test_plot_name, bbox_inches='tight', format='pdf')
