#!/usr/bin/python

from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from functions import haversine, get_oceanic_buoyancy_wCrust, pressurepoints
from functions import readbounds, organizebounds, readgrid, readdomains
from functions import findpressure_wall, findpressure_edge, project_to_point
from netCDF4 import Dataset
from plotting_spherical import shiftedColorMap
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.cm as cm
import subprocess
import sys, os, math

plt.ioff()

plates=str(sys.argv[1])	            # options: 'Pacific', 'PacPSPAus', 'LotsOfPlates', etc.
amu = float(sys.argv[2])	        # asthenospheric viscosity, Pa.s
flux_slab = int(sys.argv[3])        # 0 = no slab fluxing, 1 = slab fluxing with constant velocity, 2 = flux velocity is convergence rate
flux_width = float(sys.argv[4])     # width of layer fluxing into/out of LM. (m)
flux_alpha = float(sys.argv[5]);    # how is flux into LM partioned between SP and OP side? 1 = all on SP, 0 = all OP
no_flux_for_slabtails = int(sys.argv[6])            # if =1, this means that no lower mantle flux can occur where there is a slab tail
infile_grid = str(sys.argv[7])        		# if=1, no lower mantle flux can occur where there is a slab tail
# some big subd indices
Pac_index = int(sys.argv[8])
SAm_index = int(sys.argv[9])
# some big subd indices to not include
subd_index1 = int(sys.argv[10])
subd_index2 = int(sys.argv[11])
flux_vel_const = 50.    # only used for flux_slab = 1, [mm/yr]

taper_length = 0.
rad_km = 6378.;
if no_flux_for_slabtails == 0:
	tail_flux_string = ''
else:
	tail_flux_string = '.NoTailFlux'

# break up Pacific segments into three sections
Pac_northbound_lat=44
Pac_southbound_lat=32
SAm_bound_lat=-17

infile_grid='../input_grids/Subgrd.inp' 
infile_bounds  =''.join(['../input_bounds/Subbon_',plates,'.inp']);
infile_domains  =''.join(['../input_bounds/Subfil_',plates,'.inp']);

dips_obs_txt = '/home/aholt/Dropbox/Analytical/python/sphere/forward_flux_tests/calculate_dips/dips_from_slab2.0/AllDips.with_penetration.txt'
sp_ages_name = ''.join(['ages/',str(plates),'.txt'])
sp_ages = np.loadtxt(sp_ages_name)


## plot name given options specified
if flux_slab == 1:
	plot_string = ''.join([plates,'.',str(amu),'ConstFlux_width',str(flux_width),'_v',str(flux_vel_const),'_alpha',str(flux_alpha),tail_flux_string,'C']);
	text_string = ''.join(['../text_files/',plates,'.',str(amu),'ConstFlux_width',str(flux_width),'_v',str(flux_vel_const),'_alpha',str(flux_alpha),tail_flux_string,'C']);
elif flux_slab == 2:
	plot_string = ''.join([plates,'.',str(amu),'_VcSlabFlux_width',str(flux_width),'_alpha',str(flux_alpha),tail_flux_string,'C']);
	text_string = ''.join(['../text_files/',plates,'.',str(amu),'_VcSlabFlux_width',str(flux_width),'_alpha',str(flux_alpha),tail_flux_string,'C']);
else:
	plot_string = ''.join([plates,'.',str(amu),'noslabflux',tail_flux_string,'C']);
	text_string = ''.join(['../text_files/',plates,'.',str(amu),'noslabflux',tail_flux_string,'C']);
DPs_name = ''.join([text_string,'/DP.txt'])

print "reading input and segmenting boundaries..."
ndomain, pole_top_lon, pole_top_lat, pole_top_rate, pole_bott_lon, pole_bott_lat, pole_bott_rate, rigid_vew, rigid_vns, domain_bounds  = readdomains(infile_domains)
grid_spacing, prof_spacing, dsegtr, dseged = readgrid(infile_grid) 
num_bounds,iwall,lona,lata,lonb,latb,bound_ind,idl,idr,vt_ew,vt_ns,polarity,large_wall_inds = readbounds(infile_bounds)
iseg_min=0
n_segs,num_segs,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,polarity,num_wall_segs = \
	organizebounds(num_bounds,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,dsegtr,dseged,polarity,rad_km,iseg_min)
shift_edges=0; epslrc = 100.e3; eps_fact = 0.01  
epsdp_fact = 2*eps_fact; alith = 80.0e3; ah1=660.*1.e3   
lono,lato,gam,alpha,vtopl,vtopr,vbotl,vbotr,vt,lon_subslab,lat_subslab,lon_wedge,lat_wedge =  \
	pressurepoints(lona,lata,lonb,latb,vt_ew,vt_ns,iwall,idl,idr,n_segs,pole_top_lon,pole_top_lat,pole_top_rate,pole_bott_lon, \
		pole_bott_lat,pole_bott_rate,rigid_vew,rigid_vns,ndomain,epslrc,rad_km,alith,shift_edges,polarity,epsdp_fact)

DPs_mod = np.loadtxt(DPs_name)
dips_obs_tmp = np.loadtxt(dips_obs_txt) 	
num_obs = len(dips_obs_tmp)
num_DPs = len(DPs_mod)

synthetic_dips_tmp = np.zeros((num_obs,8))
misfit_dips = np.zeros((num_obs,3))
seg_dist = np.zeros((num_obs))
is_trench_fixed = np.zeros((num_obs))
is_just_upper_mant = np.zeros((num_obs))
num_good_subds = 0
for i in range(num_obs):

	lon_obs = dips_obs_tmp[i,0]
	lat_obs = dips_obs_tmp[i,1]

	dist_max = 1000.0;
	for j in range(num_DPs):

		lon_mod = DPs_mod[j,0]
		lat_mod = DPs_mod[j,1]

		dist = haversine(lon_obs, lat_obs, lon_mod, lat_mod,6378.)
		if dist < dist_max and iwall[j] == 1:
			dist_max = dist;
			dist_ind = j;

	lon_mod = DPs_mod[dist_ind,0]
	lat_mod = DPs_mod[dist_ind,1]
	DP = DPs_mod[dist_ind,4] 
	age = sp_ages[dist_ind,2]
	vc = DPs_mod[dist_ind,5]
	vt = DPs_mod[dist_ind,6]

	crust_density  = 3450	# rho_eclogite [kg/m3], C-T Lee, W-P Chen, 2007, EPSL
	crust_thick = 7.5e3 	# m

	# fill array...
	synthetic_dips_tmp[i,0] = lon_mod
	synthetic_dips_tmp[i,1] = lat_mod
	synthetic_dips_tmp[i,2] = DP #dip_mod
	synthetic_dips_tmp[i,3] = get_oceanic_buoyancy_wCrust(age,crust_thick,crust_density,T0=273.,DT=1300.,k=1.e-6,rho0=3300.,alpha=3.e-5)
	synthetic_dips_tmp[i,4] = large_wall_inds[dist_ind] # segment index!
	synthetic_dips_tmp[i,5] = age
	synthetic_dips_tmp[i,6] = vc
	synthetic_dips_tmp[i,7] = vt

	seg_dist[i] = dist_max 					# distance between compared segments
	if np.abs(vt) < 0.0025:
		is_trench_fixed[i] = 1			    # stationary trench?
	is_just_upper_mant[i] = dips_obs_tmp[i,3] 	# penetration? 

	if synthetic_dips_tmp[i,4] != subd_index1 and synthetic_dips_tmp[i,4] != subd_index2:
		num_good_subds = num_good_subds + 1

# get rid of subzones with index = subd_index1 or subd_index2
synthetic_dips = np.zeros((num_good_subds,8)); ind = 0
dips_obs = np.zeros((num_good_subds,3))
for i in range(num_obs):
	if synthetic_dips_tmp[i,4] != subd_index1 and synthetic_dips_tmp[i,4] != subd_index2:
		synthetic_dips[ind,0] = synthetic_dips_tmp[i,0]
		synthetic_dips[ind,1] = synthetic_dips_tmp[i,1]		
		synthetic_dips[ind,2] = synthetic_dips_tmp[i,2]
		synthetic_dips[ind,3] = synthetic_dips_tmp[i,3]
		synthetic_dips[ind,4] = synthetic_dips_tmp[i,4]
		synthetic_dips[ind,5] = synthetic_dips_tmp[i,5]
		synthetic_dips[ind,6] = synthetic_dips_tmp[i,6]
		synthetic_dips[ind,7] = synthetic_dips_tmp[i,7]
		dips_obs[ind,0] = dips_obs_tmp[i,0]
		dips_obs[ind,1] = dips_obs_tmp[i,1]
		dips_obs[ind,2] = dips_obs_tmp[i,2]
		ind = ind + 1

# # sort by subduction index
original_segment_order = np.arange(num_good_subds)
original_segment_order = original_segment_order[synthetic_dips[:,4].argsort()] 

dips_obs = dips_obs[synthetic_dips[:,4].argsort()]
is_trench_fixed = is_trench_fixed[synthetic_dips[:,4].argsort()]
is_just_upper_mant = is_just_upper_mant[synthetic_dips[:,4].argsort()]
seg_dist = seg_dist[synthetic_dips[:,4].argsort()]
synthetic_dips = synthetic_dips[synthetic_dips[:,4].argsort()]

after_Pac_index = 100
for i in range(len(synthetic_dips)):
	if synthetic_dips[i,4] > Pac_index and synthetic_dips[i,4] < after_Pac_index:
		after_Pac_index = int(synthetic_dips[i,4])

after_SAm_index = 100
for i in range(len(synthetic_dips)):
	if synthetic_dips[i,4] > SAm_index and synthetic_dips[i,4] < after_SAm_index:
		after_SAm_index = int(synthetic_dips[i,4])
	
# simplify subduction indices
last_subd_ind=1; replace_subd_ind=1;
for i in range(len(synthetic_dips)):

	if synthetic_dips[i,4] == last_subd_ind:
		synthetic_dips[i,4] = replace_subd_ind		
		if last_subd_ind == Pac_index:
			Pac_index_new = replace_subd_ind
		if last_subd_ind == SAm_index:
			SAm_index_new = replace_subd_ind
	else:
		last_subd_ind = synthetic_dips[i,4]
		if last_subd_ind == after_Pac_index:
			replace_subd_ind =  synthetic_dips[i-1,4] + 3
			synthetic_dips[i,4] = replace_subd_ind
		elif last_subd_ind == after_SAm_index:
			replace_subd_ind =  synthetic_dips[i-1,4] + 2
			synthetic_dips[i,4] = replace_subd_ind
		else:
			replace_subd_ind =  synthetic_dips[i-1,4] + 1
			synthetic_dips[i,4] = replace_subd_ind

# break up Pacific segments into three sections
for i in range(len(synthetic_dips)):
	if synthetic_dips[i,4] == Pac_index_new:
		if synthetic_dips[i,1] >= Pac_southbound_lat and synthetic_dips[i,1] < Pac_northbound_lat:
			synthetic_dips[i,4] = Pac_index_new + 1;
		elif synthetic_dips[i,1] < Pac_southbound_lat :
			synthetic_dips[i,4] = Pac_index_new + 2;

# break up South America into two sections
for i in range(len(synthetic_dips)):
	if synthetic_dips[i,4] == SAm_index_new:
		if synthetic_dips[i,1] >= SAm_bound_lat:
			synthetic_dips[i,4] = SAm_index_new + 1;


num_subd_zones = int(np.amax(synthetic_dips[:,4]))
rms_lowest = 100.; rms_lowest_shifted = 100.;
visc_min = 18.5;	visc_max = 21.5
viscosities = np.linspace(visc_min,visc_max,161)
factors = (10.0**viscosities)/amu
for k in range(len(factors)):

	test_factor = factors[k]
	rms_sum_all = 0.; rms_sum_all_shifted = 0.; mean_sum_all = 0;
	n_all = 0; n_all_shifted = 0;
	dip_array = np.zeros((num_good_subds)) # the dips predicted for the factor value
	dip_array_shifted = np.zeros((num_good_subds))
	dip_bounds = np.zeros((num_good_subds,2)) # for the error bars

	# for averaging
	dip_tots = np.zeros((num_subd_zones))
	dip_nums = np.zeros((num_subd_zones))
	dip_averages = np.zeros((num_subd_zones,2)) # for the subd zone averages
	dip_tots_shifted = np.zeros((num_subd_zones))
	dip_nums_shifted = np.zeros((num_subd_zones))
	dip_averages_shifted = np.zeros((num_subd_zones,2)) # for the subd zone averages	
	dip_obs_tots = np.zeros((num_subd_zones))
	dip_obs_nums = np.zeros((num_subd_zones))

	for i in range(num_good_subds):

		DP_mod       	= synthetic_dips[i,2] * test_factor
		plate_buoy   	= synthetic_dips[i,3]
		dip_mod      	= np.rad2deg(np.arccos((DP_mod * 1.e6)/(9.81*plate_buoy))) 
		dip_mod_shifted = np.rad2deg(np.arccos((DP_mod * 1.e6)/(9.81*plate_buoy))) - 5.0

		dip_array[i] = dip_mod
		dip_bounds[i,0] = np.rad2deg(np.arccos(((DP_mod+5) * 1.e6)/(9.81*plate_buoy))) 
		dip_bounds[i,1] = np.rad2deg(np.arccos(((DP_mod-5) * 1.e6)/(9.81*plate_buoy))) 
		dip_array_shifted[i] = dip_mod_shifted

		if np.isnan(dip_mod) == True: 
			dip_mod = 0.
		misfit_dip   = dip_mod - dips_obs[i,2]

		if np.isnan(dip_mod_shifted) == True: 
			dip_mod_shifted = 0.
		misfit_dip_shifted   = dip_mod_shifted - dips_obs[i,2]

		# averages
		dip_tots[int(synthetic_dips[i,4]-1)] = dip_tots[int(synthetic_dips[i,4]-1)] + dip_mod
		dip_nums[int(synthetic_dips[i,4]-1)] = dip_nums[int(synthetic_dips[i,4]-1)] + 1
		dip_tots_shifted[int(synthetic_dips[i,4]-1)] = dip_tots_shifted[int(synthetic_dips[i,4]-1)] + dip_mod_shifted
		dip_nums_shifted[int(synthetic_dips[i,4]-1)] = dip_nums_shifted[int(synthetic_dips[i,4]-1)] + 1
		dip_obs_tots[int(synthetic_dips[i,4]-1)] = dip_obs_tots[int(synthetic_dips[i,4]-1)] + dips_obs[i,2]
		dip_obs_nums[int(synthetic_dips[i,4]-1)] = dip_obs_nums[int(synthetic_dips[i,4]-1)] + 1

		if seg_dist[i] < 300 and is_trench_fixed[i] != 1:
			rms_sum_all = rms_sum_all + (misfit_dip)**2
			mean_sum_all = mean_sum_all + abs(misfit_dip)
			n_all = n_all + 1

		if seg_dist[i] < 300 and is_trench_fixed[i] != 1:
			rms_sum_all_shifted = rms_sum_all_shifted + (misfit_dip_shifted)**2
			n_all_shifted = n_all_shifted + 1	

	rms_all  		= round(np.sqrt(rms_sum_all/n_all),2)
	mean_all 		= round(mean_sum_all/n_all,2)
	rms_all_shifted = round(np.sqrt(rms_sum_all_shifted/n_all_shifted),2)	
	dip_averages[:,0] = dip_tots/dip_nums
	dip_averages[:,1] = np.arange((num_subd_zones))+1
	dip_averages_shifted[:,0] = dip_tots_shifted/dip_nums_shifted
	dip_averages_shifted[:,1] = np.arange((num_subd_zones))+1

	dip_obs_averages = dip_obs_tots/dip_obs_nums

	if rms_all < rms_lowest:
		rms_lowest     = rms_all
		best_factor    = test_factor
		best_rms_all   = rms_all
		best_mean_all  = mean_all
		best_dips      = dip_array
		best_dipbounds = dip_bounds
		best_avg_dips  = dip_averages

	if rms_all_shifted < rms_lowest_shifted:
		rms_lowest_shifted    = rms_all_shifted
		best_factor_shifted   = test_factor
		best_rms_all_shifted  = rms_all_shifted
		best_dips_shifted     = dip_array_shifted
		best_avg_dips_shifted = dip_averages_shifted

		
n_avg = 0; rms_sum_avg = 0; rms_sum_avg_shifted = 0; mean_sum_avg = 0
for f in range(len(dip_averages)):
	misfit_dip = best_avg_dips[f,0] - dip_obs_averages[f]
	rms_sum_avg = rms_sum_avg + (misfit_dip)**2
	mean_sum_avg = mean_sum_avg + abs(misfit_dip)
	misfit_dip_shifted = best_avg_dips_shifted[f,0] - dip_obs_averages[f]
	rms_sum_avg_shifted = rms_sum_avg_shifted + (misfit_dip_shifted)**2
	n_avg = n_avg + 1
rms_avg = round(np.sqrt(rms_sum_avg/n_avg),2)
mean_avg = round(mean_sum_avg/n_avg,2)
rms_avg_shifted = round(np.sqrt(rms_sum_avg_shifted/n_avg),2)



colormap_loc = '/home/aholt/progs/python/ScientificColourMaps5'
cmap_name='batlow'
cm_data = np.loadtxt(''.join([colormap_loc,'/',cmap_name,'/',cmap_name,'.txt']))
CBname_map = LinearSegmentedColormap.from_list(cmap_name, cm_data)



print "best dip match (%.2f) for a pressure factor of %.4f (%.2f)" % (rms_lowest,best_factor,np.log10(amu*best_factor))
print "best avg dip match (%.2f) for a pressure factor of %.4f (%.2f)" % (rms_avg,best_factor,np.log10(amu*best_factor))


best_factor_rounded = str(round(best_factor, 3))
plot_name = ''.join(['plots/best_visc_search/pdfs/wMisfits/',plot_string,'.fact',str(best_factor_rounded),'.wMisfits.pdf'])
plot_name_png = ''.join(['plots/best_visc_search/pngs/wMisfits/',plot_string,'.fact',str(best_factor_rounded),'.wMisfits.png'])
plot_name_eps = ''.join(['plots/best_visc_search/eps/wMisfits/',plot_string,'.fact',str(best_factor_rounded),'.wMisfits.eps'])

fig = plt.figure()

#### PANEL 1 OBSERVED DIPS ####
ax = fig.add_subplot(221)
map = Basemap(projection='hammer',lon_0=180,resolution='c')
map.drawcoastlines(linewidth=0.2,color='darkgray',zorder=1)
map.fillcontinents(color='darkgray',zorder=1)
# segments
for i in range (0,len(lata)):
	if iwall[i] == 1:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown',zorder=1)	
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='red',zorder=1)
	else:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.3,color='black',zorder=1)
for j in range(num_good_subds):
	for k in range(num_good_subds):
		if original_segment_order[k] == j:
			x = dips_obs[k,0]
			y = dips_obs[k,1]
			dip_point = dips_obs[k,2]
			xo, yo = map(x,y)
			dips_plot = map.scatter(xo, yo, c=dip_point,s=20,cmap=cm.get_cmap('inferno'),vmin=30,vmax=90,edgecolor='black',linewidth=0.2,zorder=2)
			# dips_plot = map.scatter(xo, yo, c=dip_point,s=20,cmap=CBname_map,vmin=30,vmax=90,edgecolor='black',linewidth=0.2,zorder=2)

cbar = map.colorbar(dips_plot,location='bottom',pad="5%",size="4%")
cbar.set_ticks(np.array([30,45,60,75,90,100]))
cbar.set_label('dip$\mathregular{_{observed}}$ [$^\circ$]',size=9)
cbar.ax.tick_params(labelsize=8)


#### PANEL 2 MODELED DIPS ####
ax = fig.add_subplot(222)
map = Basemap(projection='hammer',lon_0=180,resolution='c')
map.fillcontinents(color='darkgray',zorder=1)
map.drawcoastlines(linewidth=0.2,color='darkgray',zorder=1)
# segments
for i in range (0,len(lata)):
	if iwall[i] == 1:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown',zorder=1)	
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='red',zorder=1)
	else:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.3,color='black',zorder=1)	

for j in range(num_good_subds):
	for k in range(num_good_subds):
		if original_segment_order[k] == j:
			x = synthetic_dips[k,0]
			y = synthetic_dips[k,1]
			xo, yo = map(x,y)
			if seg_dist[k] < 300 and np.isnan(best_dips[k]) == 0: 
				if is_trench_fixed[k] == 0:  
					if best_dips[k] > 90:
						dips_plot = map.scatter(xo, yo, c=best_dips[k],s=20,cmap=cm.get_cmap('inferno'),vmin=30,vmax=90,linewidth=0.2,edgecolor='red',zorder=2)
						# dips_plot = map.scatter(xo, yo, c=best_dips[k],s=20,cmap=CBname_map,vmin=30,vmax=90,linewidth=0.2,edgecolor='red',zorder=2)
					else:
						dips_plot = map.scatter(xo, yo, c=best_dips[k],s=20,cmap=cm.get_cmap('inferno'),vmin=30,vmax=90,linewidth=0.2,edgecolor='black',zorder=2)
						# dips_plot = map.scatter(xo, yo, c=best_dips[k],s=20,cmap=CBname_map,vmin=30,vmax=90,linewidth=0.2,edgecolor='black',zorder=2)

			elif seg_dist[k] < 300 and np.isnan(best_dips[k]) == 1 and is_trench_fixed[k] == 0:
				dips_plot = map.scatter(xo, yo, s=5,color='red',zorder=2)


cbar2 = map.colorbar(dips_plot,location='bottom',pad="5%",size="4%")
cbar2.set_ticks(np.array([30,45,60,75,90,100]))
cbar2.set_label('dip$\mathregular{_{modeled}}$  [$^\circ$]',size=9)
cbar2.ax.tick_params(labelsize=8)


#### PANEL 3 MODEL-OBSERVED DIP DIFFERENCE ####
ax = fig.add_subplot(223)
map = Basemap(projection='hammer',lon_0=180,resolution='c')
map.drawcoastlines(linewidth=0.2,color='darkgray',zorder=1)
map.fillcontinents(color='darkgray',zorder=1)

indice_colors = synthetic_dips[:,4]
for j in range(len(indice_colors)):
	if indice_colors[j] > 6:
		indice_colors[j] = indice_colors[j] + 1
# segments
for i in range (0,len(lata)):
	if iwall[i] == 1:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown',zorder=1)	
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='red',zorder=1)
	else:
		if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
			pass;
		elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
			pass;
		else:
			map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.3,color='black',zorder=1)

for k in range(num_good_subds):
	for j in range(num_good_subds):
		if original_segment_order[j] == k:
			xo, yo = map(synthetic_dips[j,0],synthetic_dips[j,1]) 
			if seg_dist[j] < 300 and np.isnan(best_dips[j]) == 0 and is_trench_fixed[j] == 0: 
		            if best_dips[j] > 90:
		                misfits_plot = map.scatter(xo, yo, c=best_dips[j] - dips_obs[j,2],s=20,cmap='RdBu',vmin=-40,vmax=40,linewidth=0.2,edgecolor='red',zorder=2)
		            else:
		                misfits_plot = map.scatter(xo, yo, c=best_dips[j] - dips_obs[j,2],s=20,cmap='RdBu',vmin=-40,vmax=40,linewidth=0.2,edgecolor='black',zorder=2)
				# inds_plot = map.scatter(xo, yo, c=indice_colors[j],s=20,cmap='Paired',vmin=0,vmax=num_subd_zones,linewidth=0,zorder=3)
			elif seg_dist[j] < 300 and np.isnan(best_dips[j]) == 1 and is_trench_fixed[j] == 0:
				# inds_plot = map.scatter(xo, yo, c=indice_colors[j],s=20,cmap='Paired',vmin=0,vmax=num_subd_zones,linewidth=0,zorder=3)
				pass

cbar3 = map.colorbar(misfits_plot,location='bottom',pad="5%",size="4%")
cbar3.set_ticks(np.array([-40,-20,0,20,40]))
cbar3.set_label('dip$\mathregular{_{mod.}}$ - dip$\mathregular{_{obs.}}$  [$^\circ$]',size=9)
cbar3.ax.tick_params(labelsize=8)
ax.text(-0.04,0.0,"underpredicts",size=7, color="black",transform = ax.transAxes)
ax.text(0.825,0.0,"overpredicts",size=7, color="black",transform = ax.transAxes)

#### PANEL 4 SCATTER PLOTS ####
def fixed_aspect_ratio(ratio):
	'''
	Set a fixed aspect ratio on matplotlib plots 
	regardless of axis units
	'''
	xvals,yvals = plt.gca().axes.get_xlim(),plt.gca().axes.get_ylim()

	xrange = xvals[1]-xvals[0]
	yrange = yvals[1]-yvals[0]
	plt.gca().set_aspect(ratio*(xrange/yrange), adjustable='box')

ax = fig.add_subplot(224)
num_good_points = 0
for k in range(len(dips_obs)):
	if seg_dist[k] < 300 and np.isnan(best_dips[k]) == 0 and is_trench_fixed[k] == 0: 
		if synthetic_dips[k,4] == Pac_index_new:
			plt.scatter(dips_obs[k,2],best_dips[k],c=indice_colors[k],marker='s',cmap='Paired',vmin=0,vmax=num_subd_zones,s=8,edgecolor='black',lw=0,zorder=6,alpha=0.6)
		elif synthetic_dips[k,4] == Pac_index_new+1:
			plt.scatter(dips_obs[k,2],best_dips[k],c=indice_colors[k],marker='D',cmap='Paired',vmin=0,vmax=num_subd_zones,s=8,edgecolor='black',lw=0,zorder=6,alpha=0.6)
		elif synthetic_dips[k,4] == Pac_index_new+2:
			plt.scatter(dips_obs[k,2],best_dips[k],c=indice_colors[k],marker='^',cmap='Paired',vmin=0,vmax=num_subd_zones,s=8,edgecolor='black',lw=0,zorder=6,alpha=0.6)
		elif synthetic_dips[k,4] == SAm_index_new:
			plt.scatter(dips_obs[k,2],best_dips[k],c=indice_colors[k],marker='s',cmap='Paired',vmin=0,vmax=num_subd_zones,s=8,edgecolor='black',lw=0,zorder=6,alpha=0.6)
		elif synthetic_dips[k,4] == SAm_index_new+1:
			plt.scatter(dips_obs[k,2],best_dips[k],c=indice_colors[k],marker='D',cmap='Paired',vmin=0,vmax=num_subd_zones,s=8,edgecolor='black',lw=0,zorder=6,alpha=0.6)
		else:
			plt.scatter(dips_obs[k,2],best_dips[k],c=indice_colors[k],cmap='Paired',vmin=0,vmax=num_subd_zones,s=8,edgecolor='black',lw=0,zorder=6,alpha=0.6)
		plt.errorbar(dips_obs[k,2],best_dips[k], xerr=2.5, ecolor='lightgray', elinewidth=0.8, capsize=0, zorder=1)
		asymmetric_error = [[5],[0]]
#		plt.errorbar(dips_obs[k,2],best_dips[k], yerr=asymmetric_error, ecolor='lightgray', elinewidth=0.8, capsize=0, zorder=1)

		num_good_points = num_good_points + 1

	elif seg_dist[k] < 300 and np.isnan(best_dips[k]) == 1 and is_trench_fixed[k] == 0:
		plt.scatter(dips_obs[k,2],0.,c=indice_colors[k],cmap='Paired',vmin=0,vmax=num_subd_zones,s=8,lw=0,alpha=0.6)

good_dip_obs = np.zeros((num_good_points))
good_dip_mod = np.zeros((num_good_points))
ind = 0
for k in range(len(dips_obs)):
	if seg_dist[k] < 300 and np.isnan(best_dips[k]) == 0 and is_trench_fixed[k] == 0: 
		good_dip_mod[ind] = best_dips[k]
		good_dip_obs[ind] = dips_obs[k,2]
		ind = ind + 1
coeff_all = np.corrcoef(good_dip_mod,good_dip_obs)[1,0]

indice_colors_avg=best_avg_dips[:,1]
for j in range(len(indice_colors_avg)):
	if indice_colors_avg[j] > 6:
		indice_colors_avg[j] = indice_colors_avg[j] + 1

for j in range(len(best_avg_dips)): 
	if best_avg_dips[j,1] == Pac_index_new:
		plt.scatter(dip_obs_averages[j],best_avg_dips[j,0],c=indice_colors_avg[j],marker='s',cmap='Paired',vmin=0,vmax=num_subd_zones,s=int(dip_nums[j]*5),edgecolor='black',lw=0.3,zorder=7) # scaled by trench length
	elif best_avg_dips[j,1] == Pac_index_new+1:
		plt.scatter(dip_obs_averages[j],best_avg_dips[j,0],c=indice_colors_avg[j],marker='D',cmap='Paired',vmin=0,vmax=num_subd_zones,s=int(dip_nums[j]*5),edgecolor='black',lw=0.3,zorder=7) # scaled by trench length
	elif best_avg_dips[j,1] == Pac_index_new+2:
		plt.scatter(dip_obs_averages[j],best_avg_dips[j,0],c=indice_colors_avg[j],marker='^',cmap='Paired',vmin=0,vmax=num_subd_zones,s=int(dip_nums[j]*5),edgecolor='black',lw=0.3,zorder=7) # scaled by trench length
	elif best_avg_dips[j,1] == SAm_index_new+1: # extra one because of the above +1 shift for indices > 6
		plt.scatter(dip_obs_averages[j],best_avg_dips[j,0],c=indice_colors_avg[j],marker='s',cmap='Paired',vmin=0,vmax=num_subd_zones,s=int(dip_nums[j]*5),edgecolor='black',lw=0.3,zorder=7) # scaled by trench length	
	elif best_avg_dips[j,1] == SAm_index_new+2: # extra one because of the above +1 shift for indices > 6
		plt.scatter(dip_obs_averages[j],best_avg_dips[j,0],c=indice_colors_avg[j],marker='D',cmap='Paired',vmin=0,vmax=num_subd_zones,s=int(dip_nums[j]*5),edgecolor='black',lw=0.3,zorder=7) # scaled by trench length
	else:
		plt.scatter(dip_obs_averages[j],best_avg_dips[j,0],c=indice_colors_avg[j],cmap='Paired',vmin=0,vmax=num_subd_zones,s=int(dip_nums[j]*5),edgecolor='black',lw=0.3,zorder=7) # scaled by trench length


plt.xlabel('dip$\mathregular{_{observed}}$  [$^\circ$]',size=9)
plt.ylabel('dip$\mathregular{_{modeled}}$  [$^\circ$]',size=9)
plt.xlim(10,  100)
plt.ylim(10,  100)
plt.plot([10, 100], [10, 100], color='bisque', linewidth=6, zorder=0)

ax.tick_params(axis='x', labelsize=7.5)
ax.tick_params(axis='y', labelsize=7.5)
rms_all_string = ''.join(['$\mathregular{RMS}$ = ',str(best_rms_all),'$^\circ$, $\mathregular{RMS_{avg}}$ = ',str(rms_avg),'$^\circ$'])
ax.text(0.3,0.13,rms_all_string,size=5.5, color="black",transform = ax.transAxes)
mean_all_string_shifted = ''.join(['$\mathregular{mean}$ = ',str(best_mean_all),'$^\circ$, $\mathregular{mean_{avg}}$ = ',str(mean_avg),'$^\circ$'])
ax.text(0.3,0.07,mean_all_string_shifted,size=5.5, color="black",transform = ax.transAxes)

fixed_aspect_ratio(1)

bash_command = ''.join(['convert -density 400 -flatten ',plot_name,' ',plot_name_png]);
print plot_name
plt.savefig(plot_name, bbox_inches='tight', format='pdf')
plt.savefig(plot_name_eps, bbox_inches='tight', format='eps')
process = subprocess.Popen(['/bin/bash','-c',bash_command])
process.wait()

