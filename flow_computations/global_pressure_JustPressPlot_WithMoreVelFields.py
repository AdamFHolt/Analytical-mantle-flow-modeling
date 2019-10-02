#!/usr/bin/python
import numpy as np
import time, sys, os
from functions import readdomains, readgrid, readbounds , outputDP, organizebounds
from functions import pressurepoints, buildmatrix, outputgrids_more, buildvector, outputPprof
from plotting_spherical import plot_pressure_components_MoreVelFields, plot_plate_indices
from numpy.linalg import inv
import xarray as xr
t1 = time.time()

## command line args
plates=str(sys.argv[1])	            		# name of plate model (to read Subbon* and Subfil* files) 
amu = float(sys.argv[2])               		# asthenospheric viscosity for which computation was ran, Pa.s
flux_slab = int(sys.argv[3])        		# 0 = no slab fluxing, 1 = slab fluxing with constant velocity, 2 = flux velocity is convergence rate
flux_width = float(sys.argv[4])     		# total width of layer(s) fluxing into/out of LM. (m)
flux_alpha = float(sys.argv[5]);    		# how is flux into LM partioned between SP and OP side? 1 = all on SP, 0 = all OP
no_flux_for_slabtails = int(sys.argv[6]) 	# if=1, no lower mantle flux can occur where there is a slab tail
infile_grid = str(sys.argv[7])        		# if=1, no lower mantle flux can occur where there is a slab tail
amu_plot = float(sys.argv[8])        	# asthenospheric viscosity to plot pressure for, Pa.s

## calculation parameters
max_ocean_age = 80.
shift_edges=0   		# 1 = shift non-walls by epslrc, 0 = don't bother
epslrc = 100.e3 		# distance to relocate non-wall boundaries (for shift_edges=1)
epslit = 1.e3;  		# distance for pressure derivatives (edges)
eps_fact = 0.01  		# factor that controls how far out from boundaries, pressure gradients are set
epsdp_fact = 2*eps_fact # factor that controls how far out from walls, DPs are calculated
rad_km = 6378.;
flux_vel_const = 50. 	# only used for flux_slab = 1, [mm/yr]
iseg_min = 0		 	# minimum number of segments per boundary (0 = no minimum)
alith = 80.0e3; 		# lithospheric thickness, m
ah1=660.*1.e3   		# total thickness

## pressure-to-velocity coefficietns
coeff1=   ((ah1-alith)**3/(amu*ah1)) 
coeff2=   ((ah1-2.*alith)**3/(amu*ah1)) 
coefftr1= ((ah1-alith)**2/amu)      
coefftr2= ((ah1-2.*alith)**2/amu) 

if no_flux_for_slabtails == 0:
	tail_flux_string = ''
else:
	tail_flux_string = '.NoTailFlux'

## input files
infile_bounds  =''.join(['input_bounds/Subbon_',plates,'.inp']);
infile_domains =''.join(['input_bounds/Subfil_',plates,'.inp']);
infile_grid    =''.join(['input_grids/',infile_grid]);

## plot name given options specified
if flux_slab == 1:
	plot_string = ''.join(['plots/',plates,'.',str(amu),'ConstFlux_width',str(flux_width),'_v',str(flux_vel_const),'_alpha',str(flux_alpha),tail_flux_string,'.',str(amu_plot),'_MoreVelFields']);
	text_string = ''.join(['text_files/',plates,'.',str(amu),'ConstFlux_width',str(flux_width),'_v',str(flux_vel_const),'_alpha',str(flux_alpha),tail_flux_string,'C']);
elif flux_slab == 2:
	plot_string = ''.join(['plots/',plates,'.',str(amu),'_VcSlabFlux_width',str(flux_width),'_alpha',str(flux_alpha),tail_flux_string,'.',str(amu_plot),'_MoreVelFields']);
	text_string = ''.join(['text_files/',plates,'.',str(amu),'_VcSlabFlux_width',str(flux_width),'_alpha',str(flux_alpha),tail_flux_string,'C']);
else:
	plot_string = ''.join(['plots/',plates,'.',str(amu),'noslabflux',tail_flux_string,'.',str(amu_plot),'_MoreVelFields']);
	text_string = ''.join(['text_files/',plates,'.',str(amu),'noslabflux',tail_flux_string,'C']);

if not os.path.exists(text_string):
	os.mkdir(text_string)

print "reading input and segmenting boundaries..."
ndomain, pole_top_lon, pole_top_lat, pole_top_rate, pole_bott_lon, pole_bott_lat, pole_bott_rate, rigid_vew, rigid_vns, domain_bounds  = readdomains(infile_domains)
grid_spacing, prof_spacing, dsegtr, dseged = readgrid(infile_grid) 

# calculate distances to switch between planar and spherical solutions
thresh_dist_wall1 = 12.5 * (0.5*dsegtr*1.e3) 	# .5% error for segment vs point
thresh_dist_wall2 = 1500e3 					 	# .5% error for plane vs. sphere
thresh_dist_wall =  0.5 * (thresh_dist_wall1 + thresh_dist_wall2)  # [m]
taper_length=0. 			# length over which to taper between planar and spherical solutions

# read boundary information
num_bounds,iwall,lona,lata,lonb,latb,bound_ind,idl,idr,vt_ew,vt_ns,polarity,large_wall_inds = readbounds(infile_bounds)

# segment boundaries, and double up wall segments
n_segs,num_segs,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,polarity,num_wall_segs = \
	organizebounds(num_bounds,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,dsegtr,dseged,polarity,rad_km,iseg_min)
print "input done.\n---"

print "setting up pressure inversion points..."
lono,lato,gam,alpha,vtopl,vtopr,vbotl,vbotr,vt,lon_subslab,lat_subslab,lon_wedge,lat_wedge =  \
	pressurepoints(lona,lata,lonb,latb,vt_ew,vt_ns,iwall,idl,idr,n_segs,pole_top_lon,pole_top_lat,pole_top_rate,pole_bott_lon, \
		pole_bott_lat,pole_bott_rate,rigid_vew,rigid_vns,ndomain,epslrc,rad_km,alith,shift_edges,polarity,epsdp_fact)
print "pressure points set up.\n---"

print "loading precomputed Pcoeff and DP from %s/" % text_string 
DP_out= ''.join([text_string,'/DP.txt'])
Pcoeff_out=''.join([text_string,'/Pcoeff.txt'])
DP = np.loadtxt((DP_out))
pcoeff = np.loadtxt((Pcoeff_out))	

print "outputting solution on a grid..."
press_depth = 330.e3
P_out,Pwall_out,Pedge_out,dPdlon_out,dPdlat_out,polygon_points,plate_vel_ew,plate_vel_ns,trench_vels,avgvel_asthen_ew,avgvel_asthen_ns, pdrivenvel_wall_ew, pdrivenvel_wall_ns, pdrivenvel_edge_ew, pdrivenvel_edge_ns,lons_out,lats_out,polygons = \
	outputgrids_more(grid_spacing,lona,lata,lonb,latb,lono,lato,iwall,gam,alpha,amu,ah1-alith,n_segs,num_segs,pcoeff,rad_km,domain_bounds,bound_ind,pole_top_lon,pole_top_lat,\
		pole_top_rate,vt_ew,vt_ns,alith,press_depth,coefftr1,coefftr2,pole_bott_lon,pole_bott_lat,pole_bott_rate,rigid_vew,rigid_vns,ah1,ndomain)

print "outputting across-slab DP..."
dip_depth = 330.e3
DP = outputDP(lona,lata,lonb,latb,lono,lato,iwall,gam,alpha,n_segs,num_segs,pcoeff,rad_km,lon_subslab,\
	lat_subslab,lon_wedge,lat_wedge,polarity,vtopl,vtopr,vt,dip_depth)

DP[:,4] = DP[:,4] * (amu_plot/amu)
P_out = P_out * (amu_plot/amu)
Pwall_out = Pwall_out * (amu_plot/amu)
Pedge_out = Pedge_out * (amu_plot/amu)
avgvel_asthen_ew = avgvel_asthen_ew * (amu/amu_plot)
avgvel_asthen_ns = avgvel_asthen_ns * (amu/amu_plot)

print "plotting..."
plot_name   = ''.join([plot_string,'.pressure']);
plot_pressure_components_MoreVelFields(lons_out,lats_out,P_out,Pwall_out,Pedge_out,lona,lata,lonb,latb,lono,lato,iwall,DP,vt_ew,vt_ns,polygon_points,\
	avgvel_asthen_ew,avgvel_asthen_ns,pdrivenvel_wall_ew, pdrivenvel_wall_ns, pdrivenvel_edge_ew, pdrivenvel_edge_ns,plot_name)
