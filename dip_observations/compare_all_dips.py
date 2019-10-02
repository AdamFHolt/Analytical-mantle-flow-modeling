#!/usr/bin/python
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
from functions import haversine, pressurepoints, shiftedColorMap
from functions import readbounds, organizebounds, readgrid, readdomains
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.cm as cm
import subprocess
import sys, os, math
import glob

plt.ioff()
plot_full_dips=0
rad_km = 6378.;

dips_slab2_text='dip_catalogues/Slab2_const-depth/AllDips.txt'
dips_slab2_full_text='dip_catalogues/Slab2_const-depth/AllDips.full.txt'
dips_slab2cont_text='dip_catalogues/Slab2_two-depths/AllDips_2Contours.txt'
dips_slab2cont_full_text='dip_catalogues/Slab2_two-depths/AllDips_2Contours.full.txt'
dips_lall_txt = 'dip_catalogues/Lallemand/dips.penetration_indicated.txt'

plates='Slab2.0Final_NoJapTail_nnr_FS'	           
infile_grid='model_input/Subgrd.inp' 
infile_bounds  =''.join(['model_input/Subbon_',plates,'.inp']);
infile_domains  =''.join(['model_input/Subfil_',plates,'.inp']);

plot_name = ''.join(['plots/compare_all_dips.pdf']);
plot_name_png = ''.join(['plots/compare_all_dips.png']);

grdfiles = []
for file in glob.glob("/home/aholt/Dropbox/datasets/Slab2Distribute_Mar2018/*dep_*.grd"):
    grdfiles.append(file)

ndomain, pole_top_lon, pole_top_lat, pole_top_rate, pole_bott_lon, pole_bott_lat, pole_bott_rate, domain_bounds  = readdomains(infile_domains)

grid_spacing, prof_spacing, dsegtr, dseged = readgrid(infile_grid)

num_bounds,iwall,lona,lata,lonb,latb,bound_ind,idl,idr,vt_ew,vt_ns,polarity,large_wall_inds = readbounds(infile_bounds)

n_segs,num_segs,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,polarity,num_wall_segs = \
	organizebounds(num_bounds,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,dsegtr,dseged,polarity,rad_km)

lono,lato,gam,alpha,vtopl,vtopr,vbotl,vbotr,vt,lon_subslab,lat_subslab,lon_wedge,lat_wedge,lon_age,lat_age,lon_age2,lat_age2 =  \
	pressurepoints(lona,lata,lonb,latb,vt_ew,vt_ns,iwall,idl,idr,n_segs,pole_top_lon,pole_top_lat,pole_top_rate,pole_bott_lon, \
		pole_bott_lat,pole_bott_rate,ndomain,1.e3,rad_km,0,1.e3,polarity,150e3,50e3)

dips_lall = np.loadtxt(dips_lall_txt) 
dips_slab2 = np.loadtxt(dips_slab2_text)	
dips_slab2_full = np.loadtxt(dips_slab2_full_text)
dips_slab2cont = np.loadtxt(dips_slab2cont_text)
dips_slab2cont_full = np.loadtxt(dips_slab2cont_full_text)

num_obs = len(dips_lall)
slab2_dips_LallLocations = np.zeros((num_obs,4))
for i in range(num_obs):

	lon_obs = dips_lall[i,0]
	lat_obs = dips_lall[i,1]

	dist_max = 5000.0; dist_ind = 1e6
	for j in range(len(dips_slab2_full)):

		lon_calc = dips_slab2_full[j,0]
		lat_calc = dips_slab2_full[j,1]

		dist = haversine(lon_obs, lat_obs, lon_calc, lat_calc,6378.)
		if dist < dist_max:
			dist_max = dist;
			dist_ind = j;

	if dist_ind != 1e6:
		slab2_dips_LallLocations[i,0] = dips_slab2_full[dist_ind,0]
		slab2_dips_LallLocations[i,1] = dips_slab2_full[dist_ind,1]
		slab2_dips_LallLocations[i,2] = dips_slab2_full[dist_ind,2]
	slab2_dips_LallLocations[i,3] = dist_max


slab2cont_dips_LallLocations = np.zeros((num_obs,4))
for i in range(num_obs):

	lon_obs = dips_lall[i,0]
	lat_obs = dips_lall[i,1]

	dist_max = 5000.0; dist_ind = 1e6
	for j in range(len(dips_slab2cont_full)):

		lon_calc = dips_slab2cont_full[j,0]
		lat_calc = dips_slab2cont_full[j,1]

		dist = haversine(lon_obs, lat_obs, lon_calc, lat_calc,6378.)
		if dist < dist_max:
			dist_max = dist;
			dist_ind = j;

	if dist_ind != 1e6:
		slab2cont_dips_LallLocations[i,0] = dips_slab2cont_full[dist_ind,0]
		slab2cont_dips_LallLocations[i,1] = dips_slab2cont_full[dist_ind,1]
		slab2cont_dips_LallLocations[i,2] = dips_slab2cont_full[dist_ind,2]
	slab2cont_dips_LallLocations[i,3] = dist_max


num_obs = len(dips_slab2cont)
slab2_dips_ContLocations = np.zeros((num_obs,4))
for i in range(num_obs):

	lon_obs = dips_slab2cont[i,0]
	lat_obs = dips_slab2cont[i,1]

	dist_max = 5000.0; dist_ind = 1e6
	for j in range(len(dips_slab2_full)):

		lon_calc = dips_slab2_full[j,0]
		lat_calc = dips_slab2_full[j,1]

		dist = haversine(lon_obs, lat_obs, lon_calc, lat_calc,6378.)
		if dist < dist_max:
			dist_max = dist;
			dist_ind = j;

	if dist_ind != 1e6:
		slab2_dips_ContLocations[i,0] = dips_slab2_full[dist_ind,0]
		slab2_dips_ContLocations[i,1] = dips_slab2_full[dist_ind,1]
		slab2_dips_ContLocations[i,2] = dips_slab2_full[dist_ind,2]
	slab2_dips_ContLocations[i,3] = dist_max


fig = plt.figure(num=None, figsize=(12, 6), facecolor='w', edgecolor='k')

ax = fig.add_subplot(231)
map = Basemap(projection='hammer',lon_0=180,resolution='l')
map.drawmeridians(np.arange(0,360,30),linewidth=0.1,linecolor='gray')
map.drawparallels(np.arange(-90,90,30),linewidth=0.1,linecolor='gray')
orig_dip_cmap = cm.get_cmap('PuOr')
shifted_dip_cmap = shiftedColorMap(orig_dip_cmap, midpoint=60./90., name='shifted')
# segments
for i in range (0,len(lata)):
	if iwall[i] == 1:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown',zorder=0)	
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='red',zorder=0)
	else:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.3,color='black',zorder=0)
for j in range(num_obs):
	xo, yo = map(dips_lall[j,0],dips_lall[j,1])
	dips_plot = map.scatter(xo, yo, c=dips_lall[j,2],s=25,cmap=shifted_dip_cmap,vmin=30,vmax=90,edgecolor='black',linewidth=0.2)
cbar = map.colorbar(dips_plot,location='bottom',pad="5%",size="4%")
cbar.set_ticks(np.array([30,45,60,75,90]))
cbar.set_label('dip$\mathregular{_{Lallemand}}$ [$^\circ$]',size=9)
cbar.ax.tick_params(labelsize=8)


ax = fig.add_subplot(232)
map = Basemap(projection='hammer',lon_0=180,resolution='l')
map.drawmeridians(np.arange(0,360,30),linewidth=0.1,linecolor='gray')
map.drawparallels(np.arange(-90,90,30),linewidth=0.1,linecolor='gray')
orig_dip_cmap = cm.get_cmap('PuOr')
shifted_dip_cmap = shiftedColorMap(orig_dip_cmap, midpoint=60./90., name='shifted')
# segments
for i in range (0,len(lata)):
	if iwall[i] == 1:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown',zorder=0)	
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='red',zorder=0)
	else:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.3,color='black',zorder=0)
for j in range(len(dips_slab2)):
	xo, yo = map(dips_slab2[j,0],dips_slab2[j,1])
	dips_plot = map.scatter(xo, yo, c=dips_slab2[j,2],s=25,cmap=shifted_dip_cmap,vmin=30,vmax=90,edgecolor='black',linewidth=0.2)
cbar = map.colorbar(dips_plot,location='bottom',pad="5%",size="4%")
cbar.set_ticks(np.array([30,45,60,75,90]))
cbar.set_label('dip$\mathregular{_{Slab2.0\ (one\ depth)}}$  [$^\circ$]',size=9)
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(233)
map = Basemap(projection='hammer',lon_0=180,resolution='l')
map.drawmeridians(np.arange(0,360,30),linewidth=0.1,linecolor='gray')
map.drawparallels(np.arange(-90,90,30),linewidth=0.1,linecolor='gray')
orig_dip_cmap = cm.get_cmap('PuOr')
shifted_dip_cmap = shiftedColorMap(orig_dip_cmap, midpoint=60./90., name='shifted')
# segments
for i in range (0,len(lata)):
	if iwall[i] == 1:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown',zorder=0)	
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='red',zorder=0)
	else:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.3,color='black',zorder=0)
for j in range(len(dips_slab2cont)):
	xo, yo = map(dips_slab2cont[j,0],dips_slab2cont[j,1])
	dips_plot = map.scatter(xo, yo, c=dips_slab2cont[j,2],s=25,cmap=shifted_dip_cmap,vmin=30,vmax=90,edgecolor='black',linewidth=0.2)
cbar = map.colorbar(dips_plot,location='bottom',pad="5%",size="4%")
cbar.set_ticks(np.array([30,45,60,75,90]))
cbar.set_label('dip$\mathregular{_{Slab2.0\ (two\ depths)}}$  [$^\circ$]',size=9)
cbar.ax.tick_params(labelsize=8)


#### SCATTER PLOTS ####
def fixed_aspect_ratio(ratio):
    '''
    Set a fixed aspect ratio on matplotlib plots 
    regardless of axis units
    '''
    xvals,yvals = plt.gca().axes.get_xlim(),plt.gca().axes.get_ylim()

    xrange = xvals[1]-xvals[0]
    yrange = yvals[1]-yvals[0]
    plt.gca().set_aspect(ratio*(xrange/yrange), adjustable='box')

ax = fig.add_subplot(234)
for k in range(len(dips_lall)):
	if slab2_dips_LallLocations[k,3] < 200.:
		plt.scatter(dips_lall[k,2],slab2_dips_LallLocations[k,2],cmap='Paired',vmin=0,vmax=13,s=22.5,edgecolor='black',lw=0.25,zorder=3)
plt.xlabel('dip$\mathregular{_{Lallemand}}$  [$^\circ$]',size=9)
plt.ylabel('dip$\mathregular{_{Slab2.0\ (one\ depth)}}$  [$^\circ$]',size=9)
plt.xlim(10,  100)
plt.ylim(10,  100)
plt.plot([10, 100], [10, 100], color='bisque', linewidth=3, zorder=0)
ax.tick_params(axis='x', labelsize=7.5)
ax.tick_params(axis='y', labelsize=7.5)
fixed_aspect_ratio(1)
coeff1 = np.corrcoef(dips_lall[:,2],slab2_dips_LallLocations[:,2])[1,0]
coeff_string = ''.join(['$\mathregular{R_{Pearson}}$ = ',str(round(coeff1, 3)),' (n = ',str(len(dips_lall)),')'])
ax.text(0.05,0.91,coeff_string,size=7, color="black",transform = ax.transAxes)

ax = fig.add_subplot(235)
for k in range(len(dips_lall)):
	if slab2cont_dips_LallLocations[k,3] < 200.:
		plt.scatter(dips_lall[k,2],slab2cont_dips_LallLocations[k,2],cmap='Paired',vmin=0,vmax=13,s=22.5,edgecolor='black',lw=0.25,zorder=3)
plt.xlabel('dip$\mathregular{_{Lallemand}}$  [$^\circ$]',size=9)
plt.ylabel('dip$\mathregular{_{Slab2.0\ (two\ depths)}}$  [$^\circ$]',size=9)
plt.xlim(10,  100)
plt.ylim(10,  100)
plt.plot([10, 100], [10, 100], color='bisque', linewidth=3, zorder=0)
ax.tick_params(axis='x', labelsize=7.5)
ax.tick_params(axis='y', labelsize=7.5)
fixed_aspect_ratio(1)
coeff2 = np.corrcoef(dips_lall[:,2],slab2cont_dips_LallLocations[:,2])[1,0]
coeff_string = ''.join(['$\mathregular{R_{Pearson}}$ = ',str(round(coeff2, 3)),' (n = ',str(len(dips_lall)),')'])
ax.text(0.05,0.91,coeff_string,size=7, color="black",transform = ax.transAxes)

ax = fig.add_subplot(236)
for k in range(len(dips_slab2cont)):
	if slab2_dips_ContLocations[k,3] < 200.:
		plt.scatter(dips_slab2cont[k,2],slab2_dips_ContLocations[k,2],cmap='Paired',vmin=0,vmax=13,s=22.5,edgecolor='black',lw=0.25,zorder=3)
plt.xlabel('dip$\mathregular{_{Slab2.0\ (two\ depths)}}$  [$^\circ$]',size=9)
plt.ylabel('dip$\mathregular{_{Slab2.0\ (one\ depth)}}$  [$^\circ$]',size=9)
plt.xlim(10,  100)
plt.ylim(10,  100)
plt.plot([10, 100], [10, 100], color='bisque', linewidth=3, zorder=0)
ax.tick_params(axis='x', labelsize=7.5)
ax.tick_params(axis='y', labelsize=7.5)
fixed_aspect_ratio(1)
coeff3 = np.corrcoef(dips_slab2cont[:,2],slab2_dips_ContLocations[:,2])[1,0]
coeff_string = ''.join(['$\mathregular{R_{Pearson}}$ = ',str(round(coeff3, 3)),' (n = ',str(len(dips_slab2cont)),')'])
ax.text(0.05,0.91,coeff_string,size=7, color="black",transform = ax.transAxes)

bash_command = ''.join(['convert -density 400 -flatten ',plot_name,' ',plot_name_png]);
plt.savefig(plot_name, bbox_inches='tight', format='pdf')
process = subprocess.Popen(['/bin/bash','-c',bash_command])
process.wait()
os.remove(plot_name)
