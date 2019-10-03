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

def get_oceanic_buoyancy(age,T0,DT,k,rho0,alpha):
	T1=T0+DT;    # K
	z=np.arange(0,400.0,0.5)*1e3 # km
	T_erf= T1 - DT * scipy.special.erfc(z/(2*np.sqrt(k*age*1e6*365*24*60*60)));
	B_erf= (T1 - T_erf) * rho0 * alpha; # kg/m3
	B_erf_int=scipy.integrate.simps(y=B_erf, x=z, even='avg') # kg/m2
	return B_erf_int 

def get_oceanic_buoyancy_wCrust(age,crust_thick,crust_density,T0,DT,k,rho0,alpha):

	T1=T0+DT;    # K
	z=np.arange(0,400.0,0.5)*1e3 # m
	T_erf= T1 - DT * scipy.special.erfc(z/(2*np.sqrt(k*age*1e6*365*24*60*60)));
	B_erf= (T1 - T_erf) * rho0 * alpha; # rho - rho0, kg/m3
	B_erf_int=scipy.integrate.simps(y=B_erf, x=z, even='avg') + (crust_thick * (crust_density-rho0)) # kg/m2
	return B_erf_int 

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
	rigid_vew  = np.zeros((num_domains)) # only for ndomain=500
	rigid_vns  = np.zeros((num_domains)) # only for ndomain=500
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

		if ndomain[i] == 500:	# (fixed velocity block)
			rigid_vew[i] = f[i].split('\t')[1]
			rigid_vns[i] = f[i].split('\t')[2]

	return (ndomain, pole_top_lon, pole_top_lat, pole_top_rate, pole_bott_lon, pole_bott_lat, pole_bott_rate, rigid_vew, rigid_vns, domain_bounds)

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

def organizebounds(num_bounds,iwall,idl,idr,lona,lata,lonb,latb,bound_ind,large_wall_inds,vt_ew,vt_ns,dsegtr,dseged,polarity,rad_km,iseg_min):

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

		if iseg < iseg_min and iwall[i] == 1:
			iseg = iseg_min

		# see https://geographiclib.sourceforge.io/html/python/
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
	pole_bot_lon,pole_bot_lat,pole_bot_rate,rigid_vew,rigid_vns,ndomain,epslrc,rad_km,alith,shift_edges,polarity,epsdp_fact):

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

		DP_dist = epsdp_fact * gam[i]
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

			# for reference:
			# 100:	fixed plate, free base
			# 200:	fixed plate, fixed slab tail
			# 300:	fixed plate, fixed base
			# 400:	fixed plate, free slab stail
			# 500:  fixed plate, free base (fixed plate velocity for microplates)

			# on left side
			if idleft != 0:
				if idleft == 200 or idleft == 100 or idleft == 300 or idleft == 400:
					pole = EulerPole(lat=pole_top_lat[idl[i]-1], lon=pole_top_lon[idl[i]-1], rate=pole_top_rate[idl[i]-1])
					(vtopl_azi, vtopl_mag) = pole.velocity(lato[i],lono[i]) # degrees, mm/yr
					vtopl_vect = vtopl_mag * np.array([np.sin(np.deg2rad(vtopl_azi)),np.cos(np.deg2rad(vtopl_azi))]) # E positive, N positive
					vtopl[i] = np.dot(unit_vect_right,vtopl_vect)
				elif idleft == 500:
					vtopl_vect = np.array([rigid_vew[idl[i]-1],rigid_vns[idl[i]-1]]) # E positive, N positive
					vtopl[i] = np.dot(unit_vect_right,vtopl_vect)

				if idleft == 200 or idleft == 300:
					pole = EulerPole(lat=pole_bot_lat[idl[i]-1], lon=pole_bot_lon[idl[i]-1], rate=pole_bot_rate[idl[i]-1])
					(vbotl_azi, vbotl_mag) = pole.velocity(lato[i],lono[i]) # degrees, mm/yr
					vbotl_vect =  vtopl_mag * np.array([np.sin(np.deg2rad(vbotl_azi)),np.cos(np.deg2rad(vbotl_azi))])
					vbotl[i] = np.dot(unit_vect_right,vbotl_vect)

			# on right side 
			if idright != 0:
				if idright == 200 or idright == 100 or idright == 300 or idright == 400:
					pole = EulerPole(lat=pole_top_lat[idr[i]-1], lon=pole_top_lon[idr[i]-1], rate=pole_top_rate[idr[i]-1])
					(vtopr_azi, vtopr_mag) = pole.velocity(lato[i],lono[i])	# degrees, mm/yr
					vtopr_vect =  vtopr_mag * np.array([np.sin(np.deg2rad(vtopr_azi)),np.cos(np.deg2rad(vtopr_azi))]) 	
					vtopr[i] = np.dot(unit_vect_right,vtopr_vect)
				elif idright == 500:
					vtopr_vect = np.array([rigid_vew[idr[i]-1],rigid_vns[idr[i]-1]]) # E positive, N positive
					vtopr[i] = np.dot(unit_vect_right,vtopr_vect)

				if idright == 200 or idright == 300:
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

	return (lono,lato,gam,alpha,vtopl*vel_term,vtopr*vel_term,vbotl*vel_term,vbotr*vel_term,vt*vel_term,\
		lon_subslab,lat_subslab,lon_wedge,lat_wedge)


def buildmatrix(lona,lata,lonb,latb,gam,alpha,lono,lato,iwall,idl,idr,n_segs,num_segs,coeff1,coeff2,\
	coefftr1,coefftr2,ndomain,epslit,dsegtr,rad_km,alith,ah1,eps_fact):

	pkernel = np.zeros((n_segs,n_segs))
	for iset in range(n_segs):  # segments driving the velocity/pressure field

		lonaa = lona[iset];  lataa = lata[iset]
		lonbb = lonb[iset];  latbb = latb[iset]
		gm  = gam[iset]; 	# halfwidth [m]
		alp = alpha[iset] 	# azimuth (-pi to pi)

		for jobs in range(n_segs):  # segments where velocity/pressure is calculated

			gamma = gam[jobs]; azimuth = alpha[jobs]
			elit = eps_fact * gamma
			ebig = elit + 1.e3
			len_ratio_lit = elit/(rad_km*1e3)
			len_ratio_big = ebig/(rad_km*1e3)

			if iwall[iset] == 1 and iset >= num_segs: # wall

				lata_rad = math.radians(lata[jobs]); lona_rad = math.radians(lona[jobs]);
				latb_rad = math.radians(latb[jobs]); lonb_rad = math.radians(lonb[jobs]);

				# right side of boundary
				side = 1;	obs_azim = azimuth + (np.pi/2.)
				lonobsa, latobsa = project_to_point(lona_rad, lata_rad, obs_azim, len_ratio_lit)			
				lonobsb, latobsb = project_to_point(lonb_rad, latb_rad, obs_azim, len_ratio_lit)		
				deriv1 = calcvel_wall(math.degrees(lonobsa),math.degrees(latobsa),math.degrees(lonobsb),math.degrees(latobsb),azimuth,gamma,\
					lonaa,lataa,lonbb,latbb,gm,alp,iset,dsegtr,elit,ebig,side,rad_km)

				# left side of boundary
				side = 0;	obs_azim = azimuth - (np.pi/2.)
				lonobsa, latobsa = project_to_point(lona_rad, lata_rad, obs_azim, len_ratio_lit)			
				lonobsb, latobsb = project_to_point(lonb_rad, latb_rad, obs_azim, len_ratio_lit)	
				deriv2 = calcvel_wall(math.degrees(lonobsa),math.degrees(latobsa),math.degrees(lonobsb),math.degrees(latobsb),azimuth,gamma,\
					lonaa,lataa,lonbb,latbb,gm,alp,iset,dsegtr,elit,ebig,side,rad_km)

			elif iwall[iset] == 1 or iwall[iset] == 0: 	# edge component

				lato_rad = math.radians(lato[jobs]); lono_rad = math.radians(lono[jobs])
				
				# right side (further from boundary)
				obs_azim = azimuth + (np.pi/2.)
				lonobs , latobs = project_to_point(lono_rad, lato_rad, obs_azim, len_ratio_big)			
				Pr2 = findpressure_edge(math.degrees(lonobs),math.degrees(latobs),lonaa,lataa,lonbb,latbb,gm,rad_km)
				# right side (closer to boundary)
				lonobs , latobs = project_to_point(lono_rad, lato_rad, obs_azim, len_ratio_lit)			
				Pr1 = findpressure_edge(math.degrees(lonobs),math.degrees(latobs),lonaa,lataa,lonbb,latbb,gm,rad_km)
				deriv1 = (Pr2 - Pr1)/(ebig - elit) 	# right direction positive
				
				# left side (further from boundary)
				obs_azim = azimuth - (np.pi/2.)
				lonobs , latobs = project_to_point(lono_rad, lato_rad, obs_azim, len_ratio_big)			
				Pl2 = findpressure_edge(math.degrees(lonobs),math.degrees(latobs),lonaa,lataa,lonbb,latbb,gm,rad_km)
				# left side (closer to boundary)
				lonobs , latobs = project_to_point(lono_rad, lato_rad, obs_azim, len_ratio_lit)			
				Pl1 = findpressure_edge(math.degrees(lonobs),math.degrees(latobs),lonaa,lataa,lonbb,latbb,gm,rad_km)
				deriv2 = (Pl1 - Pl2)/(ebig-elit) 	# right direction positive

			else: # iwall[iset] == 2 (stike-slip)

				deriv1 = 0.; deriv2 = 0.

			# now weight pressure derivatives to velocity dependence on pressure...
			# derivative according to the domain boundary conditions (see p.44 in my book and Wiki's paper)
			if iwall[jobs] == 1:
				asth_thick = (ah1-alith)*1.e-3
				co1 = coefftr1 							# ((a)**2)/mu * (R/(R-a/2)) 	
				asth_thick2 = (ah1-2.*alith)*1.e-3
				co2 = coefftr2 							# ((a')**2)/mu * (R/(R-a'/2))
				# secondary, geometrical coefficients:
				second_coeff_freebase   = (2*(8*rad_km+asth_thick)) /(8*(2*rad_km-asth_thick))    	# 2(8R+a)  /8(2R-a)	
				second_coeff_freebase2  = (2*(8*rad_km+asth_thick2))/(8*(2*rad_km-asth_thick2))   	# 2(8R+a') /8(2R-a')	
				second_coeff_fixedbase  = (2*rad_km+2*asth_thick)   /(2*rad_km-asth_thick)			# (2R+2a)  /(2R-a)
				second_coeff_fixedbase2 = (2*rad_km+2*asth_thick2)  /(2*rad_km-asth_thick2)	 		# (2R+2a') /(2R-a')														
											
			else:
				asth_thick = (ah1-alith)*1.e-3
				co1 = coeff1 							# ((a)**3)/h*mu * (R/(R-a/2))
				asth_thick2 = (ah1-2.*alith)*1.e-3
				co2 = coeff2 							# ((a')**3)/h*mu * (R/(R-a/2))
 				# secondary, geometrical coefficients:
				second_coeff_freebase   = (1 + (asth_thick  /(8*rad_km)))  	# (1 + a/8R)
				second_coeff_freebase2  = (1 + (asth_thick2 /(8*rad_km)))  	# (1 + a'/8R)
				second_coeff_fixedbase  = (1 + (asth_thick  / rad_km))	# (1 + a/R)			
				second_coeff_fixedbase2 = (1 + (asth_thick2 / rad_km)) 	# (1 + a/R)		  	

			# ~~~~~~~~~~~~~~~~~~~~~~~
			# for reference:
			# 100:	fixed plate, free base
			# 200:	fixed plate, fixed slab tail
			# 300:	fixed plate, fixed base
			# 400:	fixed plate, free slab stail
			# 500:  fixed plate, free base (fixed plate velocity for microplates)
			# ~~~~~~~~~~~~~~~~~~~~~~~

			# idr: right side of boundary
			if ndomain[idr[jobs]-1] == 100 or ndomain[idr[jobs]-1] == 500:	# fixed plate, free base
				deriv1 = (deriv1/3.) * co1 * second_coeff_freebase;
			elif ndomain[idr[jobs]-1] == 400:  					# fixed plate, free slab stail
				deriv1 = (deriv1/3.) * co2 * second_coeff_freebase2;
			elif ndomain[idr[jobs]-1] == 200:					# fixed plate, fixed slab tail
				deriv1 = (deriv1/12.) * co2 * second_coeff_fixedbase2;
			elif ndomain[idr[jobs]-1] == 300:					# fixed plate, fixed base
				deriv1 = (deriv1/12.) * co1 * second_coeff_fixedbase;

			# idl: left side of boundary
			if ndomain[idl[jobs]-1] == 100 or ndomain[idl[jobs]-1] == 500:	# fixed plate, free base
				deriv2 = (deriv2/3.) * co1 * second_coeff_freebase;
			elif ndomain[idl[jobs]-1] == 400:  					# fixed plate, free slab stail
				deriv2 = (deriv2/3.) * co2 * second_coeff_freebase2;	
			elif ndomain[idl[jobs]-1] == 200:					# fixed plate, fixed slab tail
				deriv2 = (deriv2/12.) * co2 * second_coeff_fixedbase2;
			elif ndomain[idl[jobs]-1] == 300:					# fixed plate, fixed base
				deriv2 = (deriv2/12.) * co1 * second_coeff_fixedbase;

			# fill matrix elements
			if iwall[jobs] == 2:
				pkernel[jobs,iset] = 0
			elif iwall[jobs] == 0:
				pkernel[jobs,iset] = deriv1 - deriv2;
			elif jobs < num_segs and iwall[jobs] == 1:
				pkernel[jobs,iset] = deriv1
			elif iwall[jobs] == 1:
				pkernel[jobs,iset] = deriv2

	return pkernel

def calcvel_wall(lonobsa,latobsa,lonobsb,latobsb,azimuth,gamma,lonaa,lataa,lonbb,latbb,gm,alp,iset,dsegtr,elit,ebig,side,rad_km):

	# finds avg vel through lonobsa,latobsa to lonobsb,latobsb due to wall segment from lonaa,lataa to lonbb,latbb.
	# positive direction is towards the right side of the boundary (a to b)
	# first find distance between source segment and local segment 
	latobs_mid, lonobs_mid = midpoint(latobsa,lonobsa,latobsb,lonobsb) # observation points
	lat_mid, lon_mid = midpoint(lataa,lonaa,latbb,lonbb)			   # segment points
	angle = haversine(lonobs_mid,latobs_mid,lon_mid,lat_mid,rad_km) / rad_km;

	## distance at which to switch between planar and spherical solutions for computing flux at wall (function of segment length)
	thresh_dist_wall1 = 12.5 * (0.5*dsegtr*1.e3)    # distance at which less than .5% error for segment vs point [m]
	thresh_dist_wall2 = 1500e3                      # distance at which less than .5% error for plane vs. sphere [m]
	thresh_dist_wall =  0.5 * (thresh_dist_wall1 + thresh_dist_wall2)  

	A_seg = 1.; A_pt  = A_seg * (gm/(4. * rad_km * 1.e3))
	if angle * (rad_km * 1.e3) <= (thresh_dist_wall): # segment close by: need plane solution

		gamma = (1.e3 * haversine(lonobsa,latobsa,lonobsb,latobsb,rad_km))/2.

		# end point a
		dista = haversine(lonaa,lataa,lonobsa,latobsa,rad_km) * 1.e3;
		distb = haversine(lonbb,latbb,lonobsa,latobsa,rad_km) * 1.e3;
		coshlama = (.5/gm) * (dista + distb)
		sinhlama = np.sqrt(coshlama**2 -1.)
		cossiga  = (.5/gm) * (dista - distb)
		if cossiga > 1.:
			cossiga = 1.;
		elif cossiga < -1.:
			cossiga = -1;
		sinsiga=np.sqrt(1.-cossiga**2)

		# end point b
		dista = haversine(lonaa,lataa,lonobsb,latobsb,rad_km) * 1.e3;
		distb = haversine(lonbb,latbb,lonobsb,latobsb,rad_km) * 1.e3;
		coshlamb = (.5/gm) * (dista + distb)
		sinhlamb = np.sqrt(coshlamb**2 -1.)
		cossigb  = (.5/gm) * (dista - distb)
		if cossigb > 1.:
			cossigb = 1.;
		elif cossigb < -1.:
			cossigb = -1.;
		sinsigb=np.sqrt(1.-cossigb**2)

		explama = coshlama - sinhlama
		explamb = coshlamb - sinhlamb
		fvel_plane = gm*A_seg*(explama*cossiga-explamb*cossigb)/(2 * gamma) \
			-1./3.*gm*A_seg*(coshlama-sinhlama)**3*(cossiga**3-3.*cossiga*sinsiga**2)/(2 * gamma) \
			+1./3.*gm*A_seg*(coshlamb-sinhlamb)**3*(cossigb**3-3.*cossigb*sinsigb**2)/(2 * gamma)

	if angle * (rad_km * 1.e3) > (thresh_dist_wall): 

		len_ratio_close = elit/(rad_km*1.e3)
		ebig = elit + 1.e3
		len_ratio_far = ebig/(rad_km*1.e3)

		if side == 1:
			pt_azim = azimuth + (np.pi/2.) # points towards the right hand side
		else:
			pt_azim = azimuth - (np.pi/2.) # points towards the left  hand side

		# near to wall
		lonobs_close, latobs_close = project_to_point(math.radians(lonobs_mid), math.radians(latobs_mid), pt_azim, len_ratio_close)
		Pclose = findpressure_wall(math.degrees(lonobs_close),math.degrees(latobs_close),lonaa,lataa,lonbb,latbb,gm,alp,rad_km)
		# further from wall
		lonobs, latobs = project_to_point(math.radians(lonobs_mid), math.radians(latobs_mid), pt_azim, len_ratio_far)		
		Pfar = findpressure_wall(math.degrees(lonobs),math.degrees(latobs),lonaa,lataa,lonbb,latbb,gm,alp,rad_km)

		if side == 1: 	# (right side of boundary)
			fvel_sphere = (Pfar - Pclose)/(ebig - elit)
		else: 			# (left side of boundary)
			fvel_sphere = (Pclose - Pfar)/(ebig - elit)

	if angle * (rad_km * 1.e3) <= (thresh_dist_wall): 	# pure planar solution
		return fvel_plane
	elif angle * (rad_km * 1.e3) > (thresh_dist_wall):	# pure spherical solution
		return fvel_sphere



def findpressure_wall(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,alp,rad_km):

	(latmid,lonmid) = midpoint(lataa, lonaa, latbb, lonbb);
	angle =  haversine(lonmid,latmid,lonobs,latobs,rad_km) / rad_km; 
	alp = (alp + 2.*np.pi) % (2.*np.pi); # convert from -pi to pi to 0 to 2pi

	rad = rad_km * 1.e3
	A_seg = 1.; A_pt  = A_seg * (gm/(4. * rad))

	# calculate azimuth about segment location (by shifting mid-point to north pole, then azim = longitude)
	# see: gis.stackexchange.com/questions/10808/manually-transforming-rotated-lat-lon-to-regular-lat-lon
	lat_rad = math.radians(latobs)
	lon_rad = (math.radians(lonobs - lonmid - 180) + 2.*np.pi) % (2.*np.pi);
	dlat_rad = math.radians(90. - latmid)  # shift mid-point to north-pole (then azim = longitude)
	lons_transformed = math.atan2(np.sin(lon_rad), np.tan(lat_rad)*np.sin(dlat_rad) + np.cos(lon_rad)*np.cos(dlat_rad))  # radians
	pt_azim = ((lons_transformed + alp) + 2.*np.pi) % (2.*np.pi)

	# sign convention of sin(sigma) gives positive pressure on right hand side (aligns with velocity convention)
	dista = haversine(lonaa,lataa,lonobs,latobs,rad_km) * 1.e3; 
	distb = haversine(lonbb,latbb,lonobs,latobs,rad_km) * 1.e3; 
	coshlam = (.5/gm) * (dista + distb)
	sinhlam = np.sqrt(coshlam**2 -1.)
	cossig = (.5/gm) * (dista - distb)
	if cossig > 1.:
		cossig = 1.;
	elif cossig < -1.:
		cossig = -1;
	sinsig = np.sqrt(1. - (cossig**2))
	if pt_azim <= np.pi:
		sinsig = -1.0 * sinsig
	P_plane = A_seg * ( (gm * (coshlam - sinhlam) * sinsig) - (1./3.*gm*(coshlam-sinhlam)**3 * (3.*sinsig*cossig**2-sinsig**3)) )

	# Ppt, sphere
	if np.cos(angle) == 1.0:
		angle = 1e-10
	P_sphere = A_pt * gm * np.sin(-1.*pt_azim) * ((np.sin(angle))/(1.0 - np.cos(angle)))

	# Ppt, plane
	P_xy = A_pt * gm * 2. * (1./angle) * np.sin(-1.*pt_azim)

	return (P_plane - P_xy + P_sphere)


# def findpressure_edge(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,rad_km):

# 	(latmid,lonmid) = midpoint(lataa,lonaa,latbb,lonbb);
# 	lonobs_rad = math.radians(lonobs); latobs_rad = math.radians(latobs)
# 	lonmid_rad = math.radians(lonmid); latmid_rad = math.radians(latmid)
# 	dist = haversine(lonmid,latmid,lonobs,latobs,rad_km)
# 	angle =  dist / rad_km; 
# 	rad = rad_km * 1e3

# 	A_seg = 1.0;
# 	A_pt = (-1.0*A_seg*gm)/np.pi
# 	shift_pt = ((-2.0*A_seg*gm)/(np.pi))*(1.0 + np.log(rad*np.sqrt(2.)/gm)) # see "shift-for-pedge" document 
# 	Ppt_avg = A_pt *  (np.log(2) - 1) + shift_pt 	# average pressure 

# 	# Pplane
# 	dista = haversine(lonaa,lataa,lonobs,latobs,rad_km) * 1.e3; 
# 	distb = haversine(lonbb,latbb,lonobs,latobs,rad_km) * 1.e3; 
# 	coshlam = (.5/gm) * (dista + distb)
# 	sinhlam = np.sqrt(coshlam**2 - 1.)

# 	cossig = (.5/gm) * (dista - distb)
# 	if cossig > 1.:
# 		cossig = 1.;
# 	elif cossig < -1.:
# 		cossig = -1;
# 	sinsig = np.sqrt(1. - (cossig**2))
# 	xf = gm * coshlam * cossig; 
# 	yf = gm * sinhlam * sinsig;
# 	if yf == 0: yf = gm * 1.e-10
# 	x = xf/gm; y = yf/gm;

# 	# below, additional (last) term to make Laplacian = 0
# 	P_plane = ((A_seg * gm)/np.pi) * (  y*(np.arctan((x-1)/y)-np.arctan((x+1)/y)) + (.5*(x-1))*np.log(((x-1)**2 + y**2)) - (.5*(x+1))*np.log(((x+1)**2 + y**2)) ) - (((A_pt*angle**2))/4.) 
# 	P_plane = P_plane - Ppt_avg  # first , second to make average P = 0

# 	# Ppt, sphere
# 	if np.cos(angle) == 1.0:
# 		angle = 1e-10
# 	P_sphere = A_pt * np.log((1.0 - np.cos(angle))) 
# 	P_sphere = P_sphere + shift_pt - Ppt_avg

# 	# Ppt, plane
# 	P_xy = 2 * A_pt *  (np.log(angle) - np.log(np.sqrt(2.))) - (((A_pt*angle**2))/4.)
# 	P_xy = P_xy  + shift_pt - Ppt_avg

# 	return (P_plane - P_xy + P_sphere)

def findpressure_edge(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,rad_km):

	(latmid,lonmid) = midpoint(lataa,lonaa,latbb,lonbb);
	lonobs_rad = math.radians(lonobs); latobs_rad = math.radians(latobs)
	lonmid_rad = math.radians(lonmid); latmid_rad = math.radians(latmid)
	dist = haversine(lonmid,latmid,lonobs,latobs,rad_km)
	angle =  dist / rad_km; 
	rad = rad_km * 1e3

	A_seg = 1.0;
	A_pt = (-1.0*A_seg*gm)/np.pi
	shift_pt = ((-2.0*A_seg*gm)/(np.pi))*(1.0 + np.log(rad*np.sqrt(2.)/gm)) # see "shift-for-pedge" document 
	Ppt_avg = A_pt *  (np.log(2) - 1) # average pressure 

	# Pplane
	dista = haversine(lonaa,lataa,lonobs,latobs,rad_km) * 1.e3; 
	distb = haversine(lonbb,latbb,lonobs,latobs,rad_km) * 1.e3; 
	coshlam = (.5/gm) * (dista + distb)
	sinhlam = np.sqrt(coshlam**2 - 1.)

	cossig = (.5/gm) * (dista - distb)
	if cossig > 1.:
		cossig = 1.;
	elif cossig < -1.:
		cossig = -1;
	sinsig = np.sqrt(1. - (cossig**2))
	xf = gm * coshlam * cossig; 
	yf = gm * sinhlam * sinsig;
	if yf == 0: yf = gm * 1.e-10
	x = xf/gm; y = yf/gm;

	P_plane = ((A_seg * gm)/np.pi) * (  y*(np.arctan((x-1)/y)-np.arctan((x+1)/y)) + (.5*(x-1))*np.log(((x-1)**2 + y**2)) - (.5*(x+1))*np.log(((x+1)**2 + y**2)) ) 
	P_plane = P_plane - (((A_pt*angle**2))/4.) # additional term to force Laplacian = 0
	P_plane = P_plane - shift_pt

	# Ppt, sphere
	if np.cos(angle) == 1.0:
		angle = 1e-10
	P_sphere = A_pt * np.log((1.0 - np.cos(angle))) 

	# Ppt, plane
	P_xy = 2 * A_pt *  (np.log(angle) - np.log(np.sqrt(2.)))
	P_xy = P_xy - (((A_pt*angle**2))/4.) # additional term to force Laplacian = 0

	return (P_plane - P_xy + P_sphere - Ppt_avg)


def buildvector(iwall,alpha,ndomain,idl,idr,vtopl,vtopr,vbotl,vbotr,vt,n_segs,num_segs,flux_slab,flux_vel,flux_width,polarity,flux_alpha,no_flux_for_slabtails,rad_km,alith,ah1):

	vel_term = 1.e-3/(365. * 24. * 60. * 60.) # mm/yr -> m/s

	# Set up the vector (RHS) with plate/trench velocity information.
	vector  = np.zeros((n_segs));

	for jobs in range(n_segs): # get total flux due to motion of surface plates

		if flux_slab == 2:
			trench_flux_vel = (flux_width/(ah1 - alith)) * np.abs(-1.*vtopl[jobs] + vtopr[jobs])
		elif flux_slab == 1:
			trench_flux_vel = (flux_width/(ah1 - alith)) * (flux_vel * vel_term)
		else:
			trench_flux_vel = 0.

		if iwall[jobs] == 1 and polarity[jobs] == 1: 	# dipping to left
			alph = 1.0 - flux_alpha
		elif iwall[jobs] == 1 and polarity[jobs] == 2: 	# dipping to right
			alph = flux_alpha

		# for reference:
		# 100:	fixed plate, free base
		# 200:	fixed plate, fixed slab tail
		# 300:	fixed plate, fixed base
		# 400:	fixed plate, free slab stail
		# 500:  fixed plate, free base (fixed plate velocity for microplates)

		# slabs walls and fill vector
		if iwall[jobs] == 1:


			if  jobs < num_segs:  # top/right side of slab wall

				# adjust asthenosphere thickness for presence of slab tail?
				if ndomain[idr[jobs]-1] == 200 or ndomain[idr[jobs]-1] == 400:
					asth_thick = (ah1-2.*alith)*1.e-3
				else:
					asth_thick = (ah1-alith)*1.e-3

				# fixed or free asthenosphere base?
				if ndomain[idr[jobs]-1] == 100 or ndomain[idr[jobs]-1] == 400 or ndomain[idr[jobs]-1] == 500: 	# fixed plate, free base
					right_top_coeff  = 1.0
					right_bott_coeff = 0
				else:															# fixed plate, fixed base
					right_top_coeff  = (1./6.)*((6.*rad_km + asth_thick)/(2.*rad_km - asth_thick))
					right_bott_coeff = (1./6.)*((6.*rad_km - asth_thick)/(2.*rad_km - asth_thick))
				vp_component_right = vtopr[jobs]*right_top_coeff + vbotr[jobs]*right_bott_coeff

				# fill vector
				if (ndomain[idr[jobs]-1] == 200 or ndomain[idr[jobs]-1] == 400) and no_flux_for_slabtails == 1:
					vector[jobs] = -1.*vt[jobs] + vp_component_right
				else:
					vector[jobs] = -1.*(vt[jobs] - (1-alph)*trench_flux_vel) + vp_component_right

			elif jobs >= num_segs: # bott/left side of slab wall

				# adjust asthenosphere thickness for presence of slab tail?
				if ndomain[idl[jobs]-1] == 200 or ndomain[idl[jobs]-1] == 400:
					asth_thick = (ah1-2.*alith)*1.e-3
				else:
					asth_thick = (ah1-alith)*1.e-3

				# fixed or free asthenosphere base?
				if ndomain[idl[jobs]-1] == 100 or ndomain[idl[jobs]-1] == 400 or ndomain[idl[jobs]-1] == 500: 	# fixed plate, free base
					left_top_coeff  = 1.0
					left_bott_coeff = 0
				else:															# fixed plate, fixed base
					left_top_coeff  = (1./6.)*((6.*rad_km + asth_thick)/(2.*rad_km - asth_thick))
					left_bott_coeff = (1./6.)*((6.*rad_km - asth_thick)/(2.*rad_km - asth_thick))
				vp_component_left = vtopl[jobs]*left_top_coeff + vbotl[jobs]*left_bott_coeff

				# fill vector
				if (ndomain[idl[jobs]-1] == 200 or ndomain[idl[jobs]-1] == 400) and no_flux_for_slabtails == 1:
					vector[jobs] = -1.*vt[jobs] + vp_component_left
				else:
					vp_component = vtopl[jobs]*left_top_coeff + vbotl[jobs]*left_bott_coeff
					vector[jobs] = -1.*(vt[jobs] +  alph*trench_flux_vel) + vp_component_left

		# edge boundaries			
		elif ndomain[idl[jobs]-1] != 0 and ndomain[idr[jobs]-1] != 0 and iwall[jobs] != 2: # no free domains/strike-slip

			# asthenosphere thicknesses on both sides of boundary
			if ndomain[idl[jobs]-1] == 200 or ndomain[idl[jobs]-1] == 400:
				asth_thick_left = (ah1-2.*alith)*1.e-3
				alithb_left = 80.e3
			else:
				asth_thick_left = (ah1-alith)*1.e-3	
				alithb_left = 0.

			if ndomain[idr[jobs]-1] == 200 or ndomain[idr[jobs]-1] == 400:
				asth_thick_right = (ah1-2.*alith)*1.e-3
				alithb_right = 80.e3
			else:
				asth_thick_right = (ah1-alith)*1.e-3
				alithb_right = 0.

			# coeffs for right side of boundary
			if ndomain[idr[jobs]-1] == 100 or ndomain[idr[jobs]-1] == 400 or ndomain[idr[jobs]-1] == 500: 	# fixed plate, free base
				right_top_coeff  = 1.0 - ((asth_thick_right**2)/(2.*rad_km*ah1*1.e-3))
				right_bott_coeff = 0
			else:															# fixed plate, fixed base
				right_top_coeff  = (alith*1.e-3        + asth_thick_right*(0.5 + (asth_thick_right/(12.*rad_km))))/(ah1*1.e-3) 
				right_bott_coeff = (alithb_right*1.e-3 + asth_thick_right*(0.5 - (asth_thick_right/(12.*rad_km))))/(ah1*1.e-3)  

			# coeffs for left side of boundary
			if ndomain[idl[jobs]-1] == 100 or ndomain[idl[jobs]-1] == 400 or ndomain[idl[jobs]-1] == 500: 	# fixed plate, free base
				left_top_coeff  = 1.0 - ((asth_thick_left**2)/(2.*rad_km*ah1*1.e-3))
				left_bott_coeff = 0
			else:															# fixed plate, fixed base
				left_top_coeff  = (alith*1.e-3       + asth_thick_left*(0.5 + (asth_thick_left/(12.*rad_km))))/(ah1*1.e-3)     
				left_bott_coeff = (alithb_left*1.e-3 + asth_thick_left*(0.5 - (asth_thick_left/(12.*rad_km))))/(ah1*1.e-3) 


			vector[jobs] = (vtopr[jobs]*right_top_coeff + vbotr[jobs]*right_bott_coeff) - (vtopl[jobs]*left_top_coeff + vbotl[jobs]*left_bott_coeff) 

	return vector


def outputgrids(spacing,lona,lata,lonb,latb,lono,lato,iwall,gam,alpha,amu,a,n_segs,num_segs,pcoeff,\
		rad_km,domain_bounds,bound_ind,pole_top_lon,pole_top_lat,pole_top_rate,vt_ew,vt_ns,alith,press_depth,\
		coefftr1,coefftr2,pole_bott_lon,pole_bott_lat,pole_bott_rate,rigid_vew,rigid_vns,ah1,ndomain):

	## grids
	lats = np.linspace(-90,90,(180/spacing)+1)
	lons = np.linspace(0,360,(360/spacing)+1)
	lat_grd = np.zeros((len(lats),len(lons)))
	lon_grd = np.zeros((len(lats),len(lons)))
	## pressure and pressure derivatives
	P_grd = np.zeros((len(lats),len(lons)))
	Pwall_grd = np.zeros((len(lats),len(lons)))
	Pedge_grd = np.zeros((len(lats),len(lons)))
	dPdlat_grd = np.zeros((len(lats),len(lons)))
	dPdlon_grd = np.zeros((len(lats),len(lons)))
	dPedgedlat_grd = np.zeros((len(lats),len(lons)))
	dPedgedlon_grd = np.zeros((len(lats),len(lons)))
	dPwalldlat_grd = np.zeros((len(lats),len(lons)))
	dPwalldlon_grd = np.zeros((len(lats),len(lons)))
	## plate velocities
	plate_vel_ew = np.zeros((len(lats),len(lons)))
	plate_vel_ns = np.zeros((len(lats),len(lons)))
	base_vel_ew = np.zeros((len(lats),len(lons)))
	base_vel_ns = np.zeros((len(lats),len(lons)))
	## asthenospheric velocities
	avgvel_asthen_ew = np.zeros((len(lats),len(lons)))
	avgvel_asthen_ns = np.zeros((len(lats),len(lons)))
	pdrivenvel_wall_ew = np.zeros((len(lats),len(lons)))
	pdrivenvel_wall_ns = np.zeros((len(lats),len(lons)))
	pdrivenvel_edge_ew = np.zeros((len(lats),len(lons)))
	pdrivenvel_edge_ns = np.zeros((len(lats),len(lons)))

	## adjusted radius
	rad_pressloc = rad_km - (press_depth/1.e3)

	for j in range(0,len(lons)):
		lat_grd[:,j] = lats;
	for i in range(0,len(lats)):
		lon_grd[i,:] = lons;

	# generate P points on a grid 
	for j in range(0,len(lons)):
		for i in range(0,len(lats)):
		
			latobs = lats[i]
			lonobs = lons[j]
			sumpress = 0.
			sumpress_walls = 0.
			sumpress_edges = 0.

			for iset in range(n_segs):
				lonaa = lona[iset]
				lataa = lata[iset]
				lonbb = lonb[iset]
				latbb = latb[iset]
				gm  = gam[iset]
				alp = alpha[iset]

				if iset >= num_segs:
					sumpress = sumpress + findpressure_wall(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,alp,rad_km) * pcoeff[iset]
					sumpress_walls = sumpress_walls + findpressure_wall(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,alp,rad_km) * pcoeff[iset]
				else:
					sumpress = sumpress + findpressure_edge(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,rad_km) * pcoeff[iset]
					sumpress_edges = sumpress_edges + findpressure_edge(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,rad_km) * pcoeff[iset]

			P_grd[i,j] = sumpress/1.e6
			Pedge_grd[i,j] = sumpress_edges/1.e6
			Pwall_grd[i,j] = sumpress_walls/1.e6

	# # compute dP/dX (east-west), dP/dy (north-south) points on same grid 
	for j in range(0,len(lons)):
		for i in range(0,len(lats)):

			if j == 0:
				dPdlon = 	 1.e6 * (P_grd[i,j+1] - P_grd[i,j])/(spacing) 
				dPwalldlon = 1.e6 * (Pwall_grd[i,j+1] - Pwall_grd[i,j])/(spacing) 
				dPedgedlon = 1.e6 * (Pedge_grd[i,j+1] - Pedge_grd[i,j])/(spacing) 
			elif j == len(lons) - 1:
				dPdlon = 	 1.e6 * (P_grd[i,j] - P_grd[i,j-1])/(spacing) 
				dPwalldlon = 1.e6 * (Pwall_grd[i,j] - Pwall_grd[i,j-1])/(spacing) 
				dPedgedlon = 1.e6 * (Pedge_grd[i,j] - Pedge_grd[i,j-1])/(spacing) 
			else:
				dPdlon = 	 1.e6 * (P_grd[i,j+1]-P_grd[i,j-1])/(2.*spacing) 
				dPwalldlon = 1.e6 * (Pwall_grd[i,j+1]-Pwall_grd[i,j-1])/(2.*spacing) 
				dPedgedlon = 1.e6 * (Pedge_grd[i,j+1]-Pedge_grd[i,j-1])/(2.*spacing) 

			if i == 0:
				dPdlat = 	 1.e6 * (P_grd[i+1,j] - P_grd[i,j])/(spacing) 
				dPwalldlat = 1.e6 * (Pwall_grd[i+1,j] - Pwall_grd[i,j])/(spacing) 
				dPedgedlat = 1.e6 * (Pedge_grd[i+1,j] - Pedge_grd[i,j])/(spacing) 
			elif i == len(lats) - 1:
				dPdlat = 	 1.e6 * (P_grd[i,j] - P_grd[i-1,j])/(spacing)  
				dPwalldlat = 1.e6 * (Pwall_grd[i,j] - Pwall_grd[i-1,j])/(spacing)  
				dPedgedlat = 1.e6 * (Pedge_grd[i,j] - Pedge_grd[i-1,j])/(spacing)  
			else:
				dPdlat = 	 1.e6 * (P_grd[i+1,j]-P_grd[i-1,j])/(2.*spacing)
				dPwalldlat = 1.e6 * (Pwall_grd[i+1,j]-Pwall_grd[i-1,j])/(2.*spacing)
				dPedgedlat = 1.e6 * (Pedge_grd[i+1,j]-Pedge_grd[i-1,j])/(2.*spacing)

			dPdlon_grd[i,j] = np.cos(math.radians(lat_grd[i,j])) * 360. * (dPdlon/(2. * np.pi * rad_km * 1e3)) * 1e-3
			dPdlat_grd[i,j] = 360. * (dPdlat/(2. * np.pi * rad_km * 1e3)) * 1e-3; # Pa/degrees -> Pa/m -> MPa/km\

			dPwalldlon_grd[i,j] = np.cos(math.radians(lat_grd[i,j])) * 360. * (dPwalldlon/(2. * np.pi * rad_km * 1e3)) * 1e-3
			dPwalldlat_grd[i,j] = 360. * (dPwalldlat/(2. * np.pi * rad_km * 1e3)) * 1e-3; # Pa/degrees -> Pa/m -> MPa/km\

			dPedgedlon_grd[i,j] = np.cos(math.radians(lat_grd[i,j])) * 360. * (dPedgedlon/(2. * np.pi * rad_km * 1e3)) * 1e-3
			dPedgedlat_grd[i,j] = 360. * (dPedgedlat/(2. * np.pi * rad_km * 1e3)) * 1e-3; # Pa/degrees -> Pa/m -> MPa/km

	polygons, polygon_points = partition_polygon_points(lons,lats,bound_ind,lona,lata,lonb,latb,domain_bounds,rad_km)

	# calculate plate velocities for each grid point
	for j in range(0,len(lons)):
		for i in range(0,len(lats)):

			domain = int(polygon_points[i,j])

			if ndomain[domain-1] == 500:
				plate_vel_ew[i,j]  = rigid_vew[domain-1] * .1;  # cm/yr
				plate_vel_ns[i,j]  = rigid_vns[domain-1] * .1;
			else:
				pole = EulerPole(lat=pole_top_lat[domain-1], lon=pole_top_lon[domain-1], rate=pole_top_rate[domain-1])
				(vel_azi, vel_mag) = pole.velocity(lats[i],lons[j]) # degrees, mm/yr
				plate_vel_ew[i,j]  = vel_mag * np.sin(np.deg2rad(vel_azi)) * .1;  # cm/yr
				plate_vel_ns[i,j]  = vel_mag * np.cos(np.deg2rad(vel_azi)) * .1;

			if ndomain[domain-1] == 200 or ndomain[domain-1] == 300:
				pole_bott = EulerPole(lat=pole_bott_lat[domain-1], lon=pole_bott_lon[domain-1], rate=pole_bott_rate[domain-1])
				(vbott_azi, vbott_mag) = pole_bott.velocity(lats[i],lons[j]) # degrees, mm/yr
				base_vel_ew[i,j]  = vbott_mag * np.sin(np.deg2rad(vbott_azi)) * .1;  # cm/yr
				base_vel_ns[i,j]  = vbott_mag * np.cos(np.deg2rad(vbott_azi)) * .1;		

	# for reference:
	# 100:	fixed plate, free base
	# 200:	fixed plate, fixed slab tail
	# 300:	fixed plate, fixed base
	# 400:	fixed plate, free slab stail
	# 500:  fixed plate, free base (fixed plate velocity for microplates)

	# average asthenospheric velocities 
	vel_convert = 100.*60.*60.*24.*365.
	for j in range(0,len(lons)):
		for i in range(0,len(lats)):

			domain = int(polygon_points[i,j])
			if ndomain[domain-1] == 200 or ndomain[domain-1] == 400:
				asth_thick = (ah1-2.*alith)*1.e-3
				coeff = coefftr2
			else:
				asth_thick = (ah1-alith)*1.e-3
				coeff = coefftr1

			if ndomain[domain-1] == 100 or ndomain[domain-1] == 400:	# free base
				avgvel_asthen_ew[i,j] = (-1. * dPdlon_grd[i,j] * 1e3) * (coeff/3.) * (1.0+(asth_thick/(8*rad_km))) * vel_convert +  plate_vel_ew[i,j] * (1.0 - (asth_thick/(2.0 * rad_km)))
				avgvel_asthen_ns[i,j] = (-1. * dPdlat_grd[i,j] * 1e3) * (coeff/3.) * (1.0+(asth_thick/(8*rad_km))) * vel_convert +  plate_vel_ns[i,j] * (1.0 - (asth_thick/(2.0 * rad_km)))
			else:														# fixed base
				avgvel_asthen_ew[i,j] = (-1. * dPdlon_grd[i,j] * 1e3) * (coeff/12) * (1.0+(asth_thick/(rad_km))) * vel_convert \
					+ 0.5*plate_vel_ew[i,j]*(1.0 + (asth_thick/(6.0 * rad_km))) + 0.5*base_vel_ew[i,j]*(1.0 - (asth_thick/(6.0 * rad_km))) 
				avgvel_asthen_ns[i,j] = (-1. * dPdlat_grd[i,j] * 1e3) * (coeff/12) * (1.0+(asth_thick/(rad_km))) * vel_convert \
					+ 0.5*plate_vel_ns[i,j]*(1.0 + (asth_thick/(6.0 * rad_km))) + 0.5*base_vel_ns[i,j]*(1.0 - (asth_thick/(6.0 * rad_km))) 


			if ndomain[domain-1] == 100 or ndomain[domain-1] == 400:	# free base
				pdrivenvel_wall_ew[i,j] = (-1. * dPwalldlon_grd[i,j] * 1e3) * (coeff/3.) * (1.0+(asth_thick/(8*rad_km))) * vel_convert
				pdrivenvel_wall_ns[i,j] = (-1. * dPwalldlat_grd[i,j] * 1e3) * (coeff/3.) * (1.0+(asth_thick/(8*rad_km))) * vel_convert
				pdrivenvel_edge_ew[i,j] = (-1. * dPedgedlon_grd[i,j] * 1e3) * (coeff/3.) * (1.0+(asth_thick/(8*rad_km))) * vel_convert
				pdrivenvel_edge_ns[i,j] = (-1. * dPedgedlat_grd[i,j] * 1e3) * (coeff/3.) * (1.0+(asth_thick/(8*rad_km))) * vel_convert
			else:														# fixed base
				pdrivenvel_wall_ew[i,j] = (-1. * dPwalldlon_grd[i,j] * 1e3) * (coeff/12) * (1.0+(asth_thick/(rad_km))) * vel_convert
				pdrivenvel_wall_ns[i,j] = (-1. * dPwalldlat_grd[i,j] * 1e3) * (coeff/12) * (1.0+(asth_thick/(rad_km))) * vel_convert
				pdrivenvel_edge_ew[i,j] = (-1. * dPedgedlon_grd[i,j] * 1e3) * (coeff/12) * (1.0+(asth_thick/(rad_km))) * vel_convert
				pdrivenvel_edge_ns[i,j] = (-1. * dPedgedlat_grd[i,j] * 1e3) * (coeff/12) * (1.0+(asth_thick/(rad_km))) * vel_convert    


	# calculate trench velocites (for potential plotting)
	trench_vels = np.empty((0,4), float); num = 0
	for k in range(0,len(lato)):
		if iwall[k] == 1:
			if num % 10 == 0:
				trench_vel_ew = vt_ew[k] * .1;  # mm/yr -> cm/yr
				trench_vel_ns = vt_ns[k] * .1;  # mm/yr -> cm/yr
				trench_vels = np.append(trench_vels, np.array([[lono[k],lato[k],trench_vel_ew,trench_vel_ns]]), axis=0)
			num = num + 1;

	P_zfactor = rad_km/rad_pressloc
	return (P_grd*P_zfactor, Pwall_grd*P_zfactor, Pedge_grd*P_zfactor, dPdlon_grd*P_zfactor, dPdlat_grd*P_zfactor, polygon_points, plate_vel_ew, plate_vel_ns, trench_vels, avgvel_asthen_ew, avgvel_asthen_ns, \
	 pdrivenvel_wall_ew, pdrivenvel_wall_ns, pdrivenvel_edge_ew, pdrivenvel_edge_ns, lon_grd, lat_grd, polygons)


def partition_polygon_points(lons,lats,bound_ind,lona,lata,lonb,latb,domain_bounds,rad_km):

	# figure out which polygon each grid point is in
	## which plate grid points are in:
	polygon_points = np.zeros((len(lats),len(lons)))
	bound_inds = np.array(bound_ind);
	latas = np.array(lata); lonas = np.array(lona);
	latbs = np.array(latb); lonbs = np.array(lonb);
	polygons = []
	for k in range(0,len(domain_bounds)): # loop over each domain
		# create polygon using the correct boundaries
		if len(domain_bounds[k]) == 1 and int(domain_bounds[k][0]) == int(0):
			polygon_points[polygon_points == 0] = k + 1
		else:
			# get polygon point for correct boundaries
			poly_lats = []; poly_lons = [];
			for d in range(0,len(domain_bounds[k])):

				pt_inds = np.where(bound_inds == int(domain_bounds[k][d]))
				latas_pts   = latas[pt_inds]; lonas_pts   = lonas[pt_inds]
				latbs_pts   = latbs[pt_inds]; lonbs_pts   = lonbs[pt_inds]
				num_pts   = np.size(pt_inds)-1

				if len(poly_lats) > 0: # see if segment is reversed
					if haversine(lonbs_pts[num_pts],latbs_pts[num_pts],last_lon,last_lat,rad_km) < haversine(lonas_pts[0],latas_pts[0],last_lon,last_lat,rad_km):
						latas_pts = latas_pts[::-1]
						lonas_pts = lonas_pts[::-1]
						latbs_pts = latbs_pts[::-1]
						lonbs_pts = lonbs_pts[::-1]
						latas_pts, latbs_pts =  latbs_pts, latas_pts
						lonas_pts, lonbs_pts =  lonbs_pts, lonas_pts

				poly_lats = np.concatenate((poly_lats,latas_pts),axis=0)
				poly_lons = np.concatenate((poly_lons,lonas_pts),axis=0)
				last_lat = latbs_pts[num_pts]; last_lon = lonbs_pts[num_pts]
				poly_lats = np.concatenate((poly_lats,[last_lat]),axis=0)
				poly_lons = np.concatenate((poly_lons,[last_lon]),axis=0)


			# check if any segments cross the meridian (lon = 0 = 360)
			orig_crosses_meridian = 0
			for a in range(0,len(poly_lons)-1):
				if (330 < poly_lons[a] <= 360) and (0 < poly_lons[a+1] <= 30):
					orig_crosses_meridian = 1
				elif (330 < poly_lons[a+1] <= 360) and (0 < poly_lons[a] <= 30):
					orig_crosses_meridian = 1
			# print "original segment crosses meridian = %.0f" % orig_crosses_meridian


			# if meridian crossing happens, try to fix by rotating polygon longitudinally
			if orig_crosses_meridian == 1:

				poly_lons_rot = np.zeros((len(poly_lons)))
				poly_lats_rot = np.zeros((len(poly_lats)))

				crosses_meridian = 1; dlat = -90
				while crosses_meridian == 1:

					dlat = dlat + 45

					for a in range(len(poly_lons_rot)):
						lat_rad = math.radians(poly_lats[a])
						lon_rad = math.radians(poly_lons[a])  
						dlon_rad = math.radians(0.)
						dlat_rad = math.radians(dlat) 
						lons_transformed = math.atan2(np.sin(lon_rad), np.tan(lat_rad)*np.sin(dlat_rad) + np.cos(lon_rad)*np.cos(dlat_rad))  
						lats_transformed = math.asin(np.cos(dlat_rad)*np.sin(lat_rad) - np.cos(lon_rad)*np.sin(dlat_rad)*np.cos(lat_rad)) - dlon_rad
						if math.degrees(lons_transformed) < 0:
							poly_lons_rot[a] = 360. + math.degrees(lons_transformed)
						else:
							poly_lons_rot[a] = math.degrees(lons_transformed)
						poly_lats_rot[a] = math.degrees(lats_transformed)

					crosses_meridian = 0
					for a in range(0,len(poly_lons_rot)-1):
						if (330 < poly_lons_rot[a] <= 360) and (0 < poly_lons_rot[a+1] <= 30):
							crosses_meridian = 1
						elif (330 < poly_lons_rot[a+1] <= 360) and (0 < poly_lons_rot[a] <= 30):
							crosses_meridian = 1

					if dlat > 90:
						break

			# last resort: try to fix by rotating polygon longitudinally
			if orig_crosses_meridian == 1 and crosses_meridian == 1:
				poly_lons_rot = np.zeros((len(poly_lons)))
				poly_lats_rot = np.zeros((len(poly_lats)))

				crosses_meridian = 1; dlon = -180
				while crosses_meridian == 1:

					dlon = dlon + 45

					for a in range(len(poly_lons_rot)):
						lat_rad = math.radians(poly_lats[a])
						lon_rad = math.radians(poly_lons[a])  
						dlon_rad = math.radians(dlon)
						dlat_rad = math.radians(0.) 
						lons_transformed = math.atan2(np.sin(lon_rad), np.tan(lat_rad)*np.sin(dlat_rad) + np.cos(lon_rad)*np.cos(dlat_rad))  
						lats_transformed = math.asin(np.cos(dlat_rad)*np.sin(lat_rad) - np.cos(lon_rad)*np.sin(dlat_rad)*np.cos(lat_rad)) - dlon_rad
						if math.degrees(lons_transformed) < 0:
							poly_lons_rot[a] = 360. + math.degrees(lons_transformed)
						else:
							poly_lons_rot[a] = math.degrees(lons_transformed)
						poly_lats_rot[a] = math.degrees(lats_transformed)

					crosses_meridian = 0
					for a in range(0,len(poly_lons_rot)-1):
						if (330 < poly_lons_rot[a] <= 360) and (0 < poly_lons_rot[a+1] <= 30):
							crosses_meridian = 1
						elif (330 < poly_lons_rot[a+1] <= 360) and (0 < poly_lons_rot[a] <= 30):
							crosses_meridian = 1

					if dlon > 180:
						break

			if orig_crosses_meridian == 0: 
				polygon = Polygon(np.column_stack((poly_lons, poly_lats))) # create polygon
			else:
				# print "polygon %.0f rotated longitudinally by %.0f, latitudinally by %.0f to avoid crosssing meridian" % (k+1,math.degrees(dlon_rad),math.degrees(dlat_rad))
				polygon = Polygon(np.column_stack((poly_lons_rot, poly_lats_rot)))
			polygons.append(polygon)

			# check which points are inside the polygon
			if orig_crosses_meridian == 1:
				for j in range(0,len(lons)):
					for i in range(0,len(lats)):
						lon_rad = math.radians(lons[j])
						lat_rad = math.radians(lats[i])
						rot_lon = math.atan2(np.sin(lon_rad), np.tan(lat_rad)*np.sin(dlat_rad) + np.cos(lon_rad)*np.cos(dlat_rad))  
						rot_lat = math.asin(np.cos(dlat_rad)*np.sin(lat_rad) - np.cos(lon_rad)*np.sin(dlat_rad)*np.cos(lat_rad)) - dlon_rad
						if math.degrees(rot_lon) < 0:
							rot_lon = 360. + math.degrees(rot_lon)
						else:
							rot_lon = math.degrees(rot_lon)
						rot_lat = math.degrees(rot_lat)
						point = Point((rot_lon,rot_lat))
						if str(polygon.contains(point)) == "True":
							polygon_points[i,j] = k + 1;
			else:	
				for j in range(0,len(lons)):
					for i in range(0,len(lats)):
						point = Point((lons[j],lats[i]))
						if str(polygon.contains(point)) == "True":
							polygon_points[i,j] = k + 1;

	return (polygons, polygon_points)

def outputDP(lona,lata,lonb,latb,lono,lato,iwall,gam,alpha,n_segs,num_segs,pcoeff,rad_km,lon_subslab,\
        lat_subslab,lon_wedge,lat_wedge,polarity,vtopl,vtopr,vt,dip_depth):

	DP = np.zeros((len(lon_wedge),7))
	vel_term = 1.e-3/(365. * 24. * 60. * 60.) # mm/yr -> m/s
	rad_diploc = rad_km - (dip_depth/1.e3)

	for i in range(0,len(DP)):
		
		latss = lat_subslab[i]
		lonss = lon_subslab[i]
		sumpress_subslab = 0.

		latw = lat_wedge[i]
		lonw = lon_wedge[i]
		sumpress_wedge = 0.

		# get across-slab pressure differences
		for iset in range(n_segs):
			lonaa = lona[iset]
			lataa = lata[iset]
			lonbb = lonb[iset]
			latbb = latb[iset]
			gm  = gam[iset]
			alp = alpha[iset]

			if iset >= num_segs:
				sumpress_subslab  = sumpress_subslab  + findpressure_wall(lonss,latss,lonaa,lataa,lonbb,latbb,gm,alp,rad_km) * pcoeff[iset]
				sumpress_wedge	  = sumpress_wedge    + findpressure_wall(lonw,latw,lonaa,lataa,lonbb,latbb,gm,alp,rad_km) * pcoeff[iset]
			else:
				sumpress_subslab  = sumpress_subslab  + findpressure_edge(lonss,latss,lonaa,lataa,lonbb,latbb,gm,rad_km) * pcoeff[iset]
				sumpress_wedge    = sumpress_wedge    + findpressure_edge(lonw,latw,lonaa,lataa,lonbb,latbb,gm,rad_km) * pcoeff[iset]

		# convert pressures to correct depth  
		sumpress_subslab_correctz = sumpress_subslab * ((rad_km)/(rad_diploc))
		sumpress_wedge_correctz   = sumpress_wedge   * ((rad_km)/(rad_diploc))

		DP[i,0] = lono[i]
		DP[i,1] = lato[i]
		DP[i,2] = sumpress_subslab_correctz/1.e6
		DP[i,3]	= sumpress_wedge_correctz/1.e6
		DP[i,4] = (sumpress_subslab_correctz - sumpress_wedge_correctz)/1.e6
		DP[i,5] = abs(vtopl[i] - vtopr[i]) * (1./vel_term)
		if polarity[i] == 1:	# subducting to left of segment
			DP[i,6] = -1.0 * vt[i] * (1./vel_term)
		else: 					# subducting to right of segment
			DP[i,6] = vt[i] * (1./vel_term)

	return DP

def outputPprof(prof_spacing,lon1,lat1,lon2,lat2,lona,lata,lonb,latb,gam,alpha,ah1,alith,amu,n_segs,num_segs,pcoeff,rad_km):


	# get points along profiles
	length_km = haversine(lon1,lat1,lon2,lat2,rad_km)
	length_deg = np.rad2deg(length_km/rad_km)
	iprof = int(.999 * length_deg/prof_spacing) + 1
	geod = Geodesic(rad_km * 1e3, 0) # sphere
	gd = geod.Inverse(lat1, lon1, lat2, lon2) # total great-circle distance
	line = geod.Line(gd['lat1'], gd['lon1'], gd['azi1'])

	P_prof = np.zeros((iprof,4))
	coords = np.zeros((iprof,3))
	for iset in range(0,iprof):
		pointa = line.Position((gd['s12'] / iprof) * iset)
		lona2 = pointa['lon2']
		if lona2 < 0.:
			lona2 = 360. + lona2
		coords[iset,0] = lona2
		coords[iset,1] = pointa['lat2']

	for i in range(len(coords)-1):
		dist = haversine(coords[i,0],coords[i,1],coords[i+1,0],coords[i+1,1],rad_km)
		coords[i+1,2] = coords[i,2] + dist

	### PRESSURES
	for i in range(len(coords)):

		lonobs = coords[i,0]
		latobs = coords[i,1]
		P_prof[i,0] = coords[i,2]

		sumpress = 0.
		sumpress_wall = 0.
		sumpress_edge = 0.
		for iset in range(n_segs):
			lonaa = lona[iset]
			lataa = lata[iset]
			lonbb = lonb[iset]
			latbb = latb[iset]
			gm  = gam[iset]
			alp = alpha[iset]

			if iset >= num_segs:	
				sumpress = sumpress + findpressure_wall(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,alp,rad_km) * pcoeff[iset]
				sumpress_wall = sumpress_wall + findpressure_wall(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,alp,rad_km) * pcoeff[iset]
			else:
				sumpress = sumpress + findpressure_edge(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,rad_km) * pcoeff[iset]
				sumpress_edge = sumpress_edge + findpressure_edge(lonobs,latobs,lonaa,lataa,lonbb,latbb,gm,rad_km) * pcoeff[iset]

		P_prof[i,1] = sumpress/1.e6
		P_prof[i,2] = sumpress_wall/1.e6
		P_prof[i,3] = sumpress_edge/1.e6

	
	return P_prof
