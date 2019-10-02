 #!/usr/bin/python

import numpy as np
import math
import scipy.special
import scipy.integrate
from geographiclib.geodesic import Geodesic
from euler_pole import EulerPole
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt

def readdomains(infile):
	
	f=open(infile,"r").readlines()

	num_domains = sum(1 for line in f)
	ndomain = [0] * num_domains
	pole_top_lon  = np.zeros((num_domains))		# deg
	pole_top_lat  = np.zeros((num_domains))		# deg
	pole_top_rate = np.zeros((num_domains))		# deg/Myr 
	pole_bott_lon  = np.zeros((num_domains))
	pole_bott_lat  = np.zeros((num_domains))
	pole_bott_rate = np.zeros((num_domains))
	domain_bounds = []

	for i in range(num_domains):
		ndomain[i] = int(f[i].split('\t')[0])

		if f[i].split('\t')[1] != '-':
			lon = f[i].split('\t')[1]
			if lon < 0:
				pole_top_lon[i] = 360. + lon
			else:
				pole_top_lon[i] = lon
		if f[i].split('\t')[2] != '-':
			pole_top_lat[i] = f[i].split('\t')[2]
		if f[i].split('\t')[3] != '-':
			pole_top_rate[i] = f[i].split('\t')[3]
		if f[i].split('\t')[4] != '-':
			pole_bott_lon[i] = f[i].split('\t')[4]
		if f[i].split('\t')[5] != '-':
			pole_bott_lat[i] = f[i].split('\t')[5]
		if f[i].split('\t')[6] != '-':
			pole_bott_rate[i] = f[i].split('\t')[6]
		if f[i].split('\t')[8] != '-':
			boundaries = [x.strip() for x in f[i].split('\t')[8].split(',')]
			domain_bounds.append(boundaries) # get indices of segments that surround domain

		if ndomain[i] == -100:	# (plate motion on the bottom)
			pole_bott_lon[i] = pole_top_lon[i]
			pole_bott_lat[i] = pole_top_lat[i]
			pole_bott_rate[i] = pole_top_rate[i]

	return (ndomain, pole_top_lon, pole_top_lat, pole_top_rate, pole_bott_lon, pole_bott_lat, pole_bott_rate, domain_bounds)

def readgrid(infile):
	
	# read profile/grid/spacing info
	g=open(infile,"r").readlines()
	grid_spacing=float(g[0].split('\t')[0])
	prof_spacing=float(g[0].split('\t')[1])
	dsegtr=float(g[1].split('\t')[0])	# trench segent length
	dseged=float(g[1].split('\t')[1])	# edge segment length
	return (grid_spacing, prof_spacing, dsegtr, dseged)

def readbounds(infile):
	
	# read boundary info
	b=open(infile,"r").readlines()
	num_bounds = sum(1 for line in b)
	iwall = [0] * num_bounds;
	polarity = np.zeros((num_bounds))
	idl = [0] * num_bounds; idr = [0] * num_bounds; 
	lona = np.zeros((num_bounds)); lata = np.zeros((num_bounds));
	lonb = np.zeros((num_bounds)); latb = np.zeros((num_bounds)); 
	bound_ind = np.zeros((num_bounds));
	vt_ew = np.zeros((num_bounds)); 
	vt_ns = np.zeros((num_bounds)); 
	large_wall_inds = np.zeros((num_bounds))

	for i in range(num_bounds):

		iwall[i] = int(b[i].split('\t')[0])
		lona_tmp = b[i].split('\t')[1]
		if lona_tmp < 0:
			lona[i] = 360.0 + lona_tmp
		else:
			lona[i] = lona_tmp
		lata[i] = b[i].split('\t')[2]
		lonb_tmp = b[i].split('\t')[3]
		if lonb_tmp < 0:
			lonb[i] = 360.0 + lonb_tmp
		else:
			lonb[i] = lonb_tmp
		latb[i] = b[i].split('\t')[4]
		bound_ind[i] = i + 1;

		if iwall[i] < 5: # edges / walls
			idl[i] = int(b[i].split('\t')[5]) # plate to the left
			idr[i] = int(b[i].split('\t')[6]) # plate to the right
		if iwall[i] == 1: # wall
			vt_ew[i] = b[i].split('\t')[7]
			vt_ns[i] = b[i].split('\t')[8]
			if b[i].split('\t')[9] == 'l': # subducting to left of segment
				polarity[i] = 1;
			elif b[i].split('\t')[9] == 'r': # subducting to the right of segment
				polarity[i] = 2;
			large_wall_inds[i] = int(b[i].split('\t')[10])

	return (num_bounds,iwall,lona,lata,lonb,latb,bound_ind,idl,idr,vt_ew,vt_ns,polarity,large_wall_inds)


def haversine(lon1, lat1, lon2, lat2, rad_km):
	#Great circle distance (km's) between two points (given in degrees)
	lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
	dlon = lon2 - lon1 
	dlat = lat2 - lat1 
	a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
	c = 2 * math.asin(math.sqrt(a)) 
	return c * rad_km


def midpoint(latA, lonA, latB, lonB):
	lonA = math.radians(lonA); lonB = math.radians(lonB)
	latA = math.radians(latA); latB = math.radians(latB)
	dLon = lonB - lonA
	Bx = math.cos(latB) * math.cos(dLon)
	By = math.cos(latB) * math.sin(dLon)
	latC = math.atan2(math.sin(latA) + math.sin(latB),math.sqrt((math.cos(latA) + Bx) * (math.cos(latA) + Bx) + By * By))
	lonC = lonA + math.atan2(By, math.cos(latA) + Bx)
	lonC = (lonC + 3 * math.pi) % (2 * math.pi) - math.pi
	return math.degrees(latC), math.degrees(lonC)

def project_to_point(lon1, lat1, azim, len_ratio):
	#Point at a given great circle angle and azimuth fron (lon1,lat1)
	lat2 = math.asin( math.sin(lat1)*math.cos(len_ratio) + math.cos(lat1)*math.sin(len_ratio)*math.cos(azim))
	lon2 = lon1 + math.atan2(math.sin(azim)*math.sin(len_ratio)*math.cos(lat1),math.cos(len_ratio)-math.sin(lat1)*math.sin(lat2))
	lon2 = (lon2 + 2.*np.pi) % (2.*np.pi);
	return (lon2,lat2)

def organizebounds(num_bounds,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,dsegtr,dseged,polarity,rad_km):

	# divide trenches and edges into segments
	lona_temp  = []; lata_temp  = []; 
	lonb_temp  = []; latb_temp  = []; iwall_temp = [];
	idl_temp   = []; idr_temp   = []; bound_ind_temp = []
	vt_ew_temp = []; vt_ns_temp = []; 
	polarity_temp = [];
	large_wall_inds_temp = [];

	num_segs = 0; num_wall_segs = 0
	print "------------------------"
	print "%.0f original boundaries" % num_bounds
	for i in range(num_bounds):

		# great circle distance [km]
		length = haversine(lona[i],lata[i],lonb[i],latb[i],rad_km)
		print "orig. length of segment %.0f = %.8f km" % (i+1,length)
		if iwall[i] == 1:
			iseg = int(.999 * length/dsegtr) + 1
		elif (idl[i] != idr[i]) and iwall[i] == 0:
			iseg = int(.999 * length/dseged) + 1
		else:
			iseg = 1  	# strike-slip (iwall=2)

		# https://geographiclib.sourceforge.io/html/python/
		geod = Geodesic(rad_km * 1e3, 0) # sphere
		gd = geod.Inverse(lata[i], lona[i], latb[i], lonb[i]) # total great-circle distance
		line = geod.Line(gd['lat1'], gd['lon1'], gd['azi1'])
		for iset in range(0,iseg):
			pointa = line.Position((gd['s12'] / iseg) * iset)
			pointb = line.Position((gd['s12'] / iseg) * (iset+1))
			lona2 = pointa['lon2']
			if lona2 < 0.:
				lona2 = 360. + lona2
			lona_temp.append(lona2)
			lata_temp.append(pointa['lat2'])
			lonb2 = pointb['lon2']
			if lonb2 < 0.:
				lonb2 = 360. + lonb2
			lonb_temp.append(lonb2)
			latb_temp.append(pointb['lat2'])

			iwall_temp.append(iwall[i]);
			idl_temp.append(idl[i]);
			idr_temp.append(idr[i]);
			bound_ind_temp.append(bound_ind[i])
			large_wall_inds_temp.append(large_wall_inds[i])
			vt_ew_temp.append(vt_ew[i])
			vt_ns_temp.append(vt_ns[i])
			polarity_temp.append(polarity[i])

		num_segs += iseg;
		if iwall[i] == 1:
			num_wall_segs += iseg;

	# # Double up wall boundaries
	n_segs = num_segs + num_wall_segs
	lona = np.zeros((n_segs)); lata = np.zeros((n_segs))
	lonb = np.zeros((n_segs)); latb = np.zeros((n_segs))
	vt_ew = np.zeros((n_segs)); vt_ns = np.zeros((n_segs)); 
	polarity = np.zeros((n_segs)); 
	iwall = [0] * n_segs; bound_ind = [0] * n_segs
	idl = [0] * n_segs; idr = [0] * n_segs; large_wall_inds = [0] * n_segs
	nwall = 0;
	for i in range(num_segs):
		lona[i] = lona_temp[i]
		lata[i] = lata_temp[i]
		lonb[i] = lonb_temp[i]
		latb[i] = latb_temp[i]
		iwall[i] = iwall_temp[i]
		idl[i] = idl_temp[i]
		idr[i] = idr_temp[i]
		bound_ind[i] = bound_ind_temp[i]
		large_wall_inds[i] = large_wall_inds_temp[i]
		vt_ew[i] = vt_ew_temp[i]
		vt_ns[i] = vt_ns_temp[i]
		polarity[i] = polarity_temp[i]
		if iwall_temp[i] == 1:
			lona[num_segs+nwall] = lona_temp[i]
			lata[num_segs+nwall] = lata_temp[i]
			lonb[num_segs+nwall] = lonb_temp[i]
			latb[num_segs+nwall] = latb_temp[i]
			iwall[num_segs+nwall] = iwall_temp[i]
			idl[num_segs+nwall] = idl_temp[i]
			idr[num_segs+nwall] = idr_temp[i]
			bound_ind[num_segs+nwall] = bound_ind_temp[i]
			large_wall_inds[num_segs+nwall] = large_wall_inds_temp[i]
			vt_ew[num_segs+nwall] = vt_ew_temp[i]
			vt_ns[num_segs+nwall] = vt_ns_temp[i]
			polarity[num_segs+nwall] = polarity_temp[i]
			nwall += 1;

	return (n_segs,num_segs,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,polarity,num_wall_segs)


def pressurepoints(lona,lata,lonb,latb,vt_ew,vt_ns,iwall,idl,idr,n_segs,pole_top_lon,pole_top_lat,pole_top_rate, \
	pole_bot_lon,pole_bot_lat,pole_bot_rate,ndomain,epslrc,rad_km,shift_edges,DP_dist,polarity,age_dist1,age_dist2):

	# Set up the points about which to conduct the inversion for mantle pressure.
	# Also find velocities of plates and slabs that are orthogonal to boundaries.

	vel_term = 1.e-3/(365. * 24. * 60. * 60.) # mm/yr -> m/s
	vtopl = np.zeros((n_segs));
	vbotl = np.zeros((n_segs));
	vtopr = np.zeros((n_segs));
	vbotr = np.zeros((n_segs));
	vt = np.zeros((n_segs));
	lato = np.zeros((n_segs));
	lono = np.zeros((n_segs));
	alpha = np.zeros((n_segs));
	gam = np.zeros((n_segs));

	lon_subslab = np.zeros((n_segs));
	lat_subslab = np.zeros((n_segs));
	lon_wedge 	= np.zeros((n_segs));
	lat_wedge	= np.zeros((n_segs));

	lon_age		= np.zeros((n_segs));
	lat_age		= np.zeros((n_segs));
	lon_age2	= np.zeros((n_segs));
	lat_age2	= np.zeros((n_segs));

	for i in range(n_segs):

		lato[i], lono[i] = midpoint(lata[i],lona[i],latb[i],lonb[i])
		if lono[i] < 0.:
			lono[i] = 360. + lono[i]
		lono_rad=math.radians(lono[i]);	lato_rad=math.radians(lato[i])
		lona_rad=math.radians(lona[i]);	lata_rad=math.radians(lata[i])
		lonb_rad=math.radians(lonb[i]);	latb_rad=math.radians(latb[i])
		
		gam[i] = (1.e3 * haversine(lona[i],lata[i],lonb[i],latb[i],rad_km))/2.  # halfwidth [m]
		if gam[i] <= 0:
			gam[i] = 1.e3

		alpha[i] = math.atan2(np.sin(lonb_rad-lona_rad)*np.cos(latb_rad),\
			np.cos(lata_rad)*np.sin(latb_rad) - np.sin(lata_rad)*np.cos(latb_rad)*np.cos(lonb_rad-lona_rad)); # -pi to pi (from north)

		idleft = ndomain[idl[i]-1];   # plate to left
		idright = ndomain[idr[i]-1];  # plate to right

		# calculate locations for calculating DP and for extracting SP age

		if polarity[i] == 1: # dipping to left
			reloc_azim_wedge   = alpha[i] - (np.pi/2.)
			reloc_azim_subslab = alpha[i] + (np.pi/2.)
		else: # dipping to right (and edges which are meaningless for this)
			reloc_azim_wedge   = alpha[i] + (np.pi/2.)
			reloc_azim_subslab = alpha[i] - (np.pi/2.)

		DP_len_ratio = DP_dist/(rad_km * 1e3)
		lon_wedge_rad , lat_wedge_rad = project_to_point(lono_rad, lato_rad, reloc_azim_wedge, DP_len_ratio)
		if lon_wedge_rad < 0.:
			lon_wedge_rad = (2.*np.pi) + lon_wedge_rad
		lon_wedge[i] = math.degrees(lon_wedge_rad)
		lat_wedge[i] = math.degrees(lat_wedge_rad)

		lon_subslab_rad , lat_subslab_rad = project_to_point(lono_rad, lato_rad, reloc_azim_subslab, DP_len_ratio)
		if lon_subslab_rad < 0.:
			lon_subslab_rad = (2.*np.pi) + lon_subslab_rad
		lon_subslab[i] = math.degrees(lon_subslab_rad)
		lat_subslab[i] = math.degrees(lat_subslab_rad)

		age_len_ratio1 = age_dist1/(rad_km * 1e3)
		lon_age_rad , lat_age_rad = project_to_point(lono_rad, lato_rad, reloc_azim_subslab, age_len_ratio1)
		if lon_age_rad < 0.:
			lon_age_rad = (2.*np.pi) + lon_age_rad
		lon_age[i] = math.degrees(lon_age_rad)
		lat_age[i] = math.degrees(lat_age_rad)

		age_len_ratio2 = age_dist2/(rad_km * 1e3)
		lon_age2_rad , lat_age2_rad = project_to_point(lono_rad, lato_rad, reloc_azim_subslab, age_len_ratio2)
		if lon_age2_rad < 0.:
			lon_age2_rad = (2.*np.pi) + lon_age2_rad
		lon_age2[i] = math.degrees(lon_age2_rad)
		lat_age2[i] = math.degrees(lat_age2_rad)

		# relocate non-walls by epslrc toward less "fixed" side of the boundary
		if iwall[i] == 0 and shift_edges == 1:
			len_ratio = epslrc/(rad_km * 1e3)
			if abs(idleft) < abs(idright):
				reloc_azim = alpha[i] - (np.pi/2.)
				lono_rad_new , lato_rad_new = project_to_point(lono_rad, lato_rad, reloc_azim, len_ratio)
				if lono_rad_new < 0.:
					lono_rad_new = (2.*np.pi) + lono_rad_new
				lato[i] = math.degrees(lato_rad_new)
				lono[i] = math.degrees(lono_rad_new)
			elif abs(idleft) > abs(idright):
				reloc_azim = alpha[i] + (np.pi/2.)
				lono_rad_new , lato_rad_new = project_to_point(lono_rad, lato_rad, reloc_azim, len_ratio)
				if lono_rad_new < 0.:
					lono_rad_new = (2.*np.pi) + lono_rad_new
				lato[i] = math.degrees(lato_rad_new)
				lono[i] = math.degrees(lono_rad_new)
		
		# find velocities that are orthogonal to the plate boundaries.
		# positive velocities point to right side of the boundary (idr direction)
		if iwall[i] != 2: # (i.e. not strike-slip)

			# get unit vector pointing towards right hand side
			vproj_azi = alpha[i] + (np.pi/2.)
			unit_vect_right = np.array([np.sin(vproj_azi) , np.cos(vproj_azi)]);  # E-W (E positive), N-S (N positive)

			# on left side
			if idleft != 0:
				if idleft == 200 or idleft == 100 or idleft == 300:
					pole = EulerPole(lat=pole_top_lat[idl[i]-1], lon=pole_top_lon[idl[i]-1], rate=pole_top_rate[idl[i]-1])
					(vtopl_azi, vtopl_mag) = pole.velocity(lato[i],lono[i]) # degrees, mm/yr
					vtopl_vect = vtopl_mag * np.array([np.sin(np.deg2rad(vtopl_azi)),np.cos(np.deg2rad(vtopl_azi))]) # E positive, N positive
					vtopl[i] = np.dot(unit_vect_right,vtopl_vect)

				if idleft == 200 or idleft == -100 or idleft == 300:
					pole = EulerPole(lat=pole_bot_lat[idl[i]-1], lon=pole_bot_lon[idl[i]-1], rate=pole_bot_rate[idl[i]-1])
					(vbotl_azi, vbotl_mag) = pole.velocity(lato[i],lono[i]) # degrees, mm/yr
					vbotl_vect =  vtopl_mag * np.array([np.sin(np.deg2rad(vbotl_azi)),np.cos(np.deg2rad(vbotl_azi))])
					vbotl[i] = np.dot(unit_vect_right,vbotl_vect)

			# on right side 
			if idright != 0:
				if idright == 200 or idright == 100 or idright == 300:
					pole = EulerPole(lat=pole_top_lat[idr[i]-1], lon=pole_top_lon[idr[i]-1], rate=pole_top_rate[idr[i]-1])
					(vtopr_azi, vtopr_mag) = pole.velocity(lato[i],lono[i])	# degrees, mm/yr
					vtopr_vect =  vtopr_mag * np.array([np.sin(np.deg2rad(vtopr_azi)),np.cos(np.deg2rad(vtopr_azi))]) 	
					vtopr[i] = np.dot(unit_vect_right,vtopr_vect)

				if idright == 200 or idright == -100 or idright == 300:
					pole = EulerPole(lat=pole_bot_lat[idr[i]-1], lon=pole_bot_lon[idr[i]-1], rate=pole_bot_rate[idr[i]-1])
					(vbotr_azi, vbotr_mag) = pole.velocity(lato[i],lono[i]) # degrees, mm/yr
					vbotr_vect =  vbotr_mag * np.array([np.sin(np.deg2rad(vbotr_azi)),np.cos(np.deg2rad(vbotr_azi))])  
					vbotr[i] = np.dot(unit_vect_right,vbotr_vect)

		# # trench/wall velocities are orthogonal to walls, and also point towards right hand side
		if iwall[i] == 1:
			vproj_azi = alpha[i] + (np.pi/2.)
			unit_vect_right = np.array([np.sin(vproj_azi) , np.cos(vproj_azi)]);  # E-W (E positive), N-S (N positive)
			vt_vect = np.array([vt_ew[i] , vt_ns[i]])  # mm/yr
			vt[i] = np.dot(unit_vect_right,vt_vect)

	return (lono,lato,gam,alpha,vtopl*vel_term,vtopr*vel_term,vbotl*vel_term,vbotr*vel_term,vt*vel_term,lon_subslab,lat_subslab,lon_wedge,lat_wedge,lon_age,lat_age,lon_age2,lat_age2)

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
	'''
	Function to offset the "center" of a colormap. Useful for
	data with a negative min and positive max and you want the
	middle of the colormap's dynamic range to be at zero
	'''
	cdict = {
		'red': [],
		'green': [],
		'blue': [],
		'alpha': []
	}

	# regular index to compute the colors
	reg_index = np.linspace(start, stop, 257)

	# shifted index to match the data
	shift_index = np.hstack([
		np.linspace(0.0, midpoint, 128, endpoint=False), 
		np.linspace(midpoint, 1.0, 129, endpoint=True)
	])

	for ri, si in zip(reg_index, shift_index):
		r, g, b, a = cmap(ri)

		cdict['red'].append((si, r, r))
		cdict['green'].append((si, g, g))
		cdict['blue'].append((si, b, b))
		cdict['alpha'].append((si, a, a))

	newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
	plt.register_cmap(cmap=newcmap)

	return newcmap
