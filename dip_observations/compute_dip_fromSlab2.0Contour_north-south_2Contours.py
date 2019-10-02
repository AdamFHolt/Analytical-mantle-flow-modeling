#!/usr/bin/python 

from pandas import *
import sys
import numpy as np
import gpxpy.geo
import geopy.distance
import math

lat_min=float(sys.argv[1])     
lat_max=float(sys.argv[2])     
trench_file=sys.argv[3]
deep_file=sys.argv[4]
shall_file=sys.argv[5]
dip_dir=int(sys.argv[6])
smooth=int(sys.argv[7])
box_width=int(sys.argv[8])
dip_int=float(sys.argv[9]) 
z_deep=float(sys.argv[10]) 
z_shall=float(sys.argv[11]) 

trench_loc=np.loadtxt(trench_file)
shall_slab=np.loadtxt(shall_file)
deep_slab=np.loadtxt(deep_file)

lats=np.arange(lat_min,lat_max,(lat_max-lat_min)/300.);
dips_tmp=np.zeros([len(lats),6])

# normalize box width so it's an ~equivalent degree distance
box_width=int(round(box_width*(10./(lat_max-lat_min))))

def haversine(lon1, lat1, lon2, lat2):
	#Great circle distance (km's) between two points (given in degrees)
	lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
	dlon = lon2 - lon1 
	dlat = lat2 - lat1 
	a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
	c = 2 * math.asin(math.sqrt(a)) 
	return c * 6371.0

net_distance=0; num_bad = 0
for n in range(len(lats)):

		dlat = 90.0; ind_trench = 100;
		for i in range(len(trench_loc)):
			dlat_tmp = np.abs(lats[n] - trench_loc[i,1]);
			if dlat_tmp < dlat:
				dlat = dlat_tmp;
				ind_trench = i;

		# trench angles (clockwise from N)
		if ind_trench != (len(trench_loc)-1):
			lat1=np.deg2rad(trench_loc[ind_trench,1]); lon1=np.deg2rad(trench_loc[ind_trench,0]); lat2=np.deg2rad(trench_loc[ind_trench+1,1]); lon2=np.deg2rad(trench_loc[ind_trench+1,0]);
		else:
			lat1=np.deg2rad(trench_loc[ind_trench-1,1]); lon1=np.deg2rad(trench_loc[ind_trench-1,0]); lat2=np.deg2rad(trench_loc[ind_trench,1]); lon2=np.deg2rad(trench_loc[ind_trench,0]);
		trench_az = np.rad2deg(math.atan2( np.sin(lon2-lon1)*np.cos(lat2) , np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1) ));

		if dip_dir == 0: # left dipping
			search_az =  trench_az - 90;
		else: # right dipping
			search_az =  trench_az + 90;

		# search for point at two different depths (along azimuth), "deep" and "shallow"
		dang1 = 360.0;
		for j in range(len(shall_slab)):
			lat1=np.deg2rad(trench_loc[ind_trench,1]);lon1=np.deg2rad(trench_loc[ind_trench,0]);lat2=np.deg2rad(shall_slab[j,1]);lon2=np.deg2rad(shall_slab[j,0]);
			ang = np.rad2deg(math.atan2( np.sin(lon2-lon1)*np.cos(lat2) , np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1) ));
			if abs(search_az - ang) < dang1:
				dang1 = abs(search_az - ang);
				ind_shall = j
		dang2 = 360.0;
		for j in range(len(deep_slab)):
			lat1=np.deg2rad(trench_loc[ind_trench,1]);lon1=np.deg2rad(trench_loc[ind_trench,0]);lat2=np.deg2rad(deep_slab[j,1]);lon2=np.deg2rad(deep_slab[j,0]);
			ang = np.rad2deg(math.atan2( np.sin(lon2-lon1)*np.cos(lat2) , np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1) ));
			if abs(search_az - ang) < dang2:
				dang2 = abs(search_az - ang);
				ind_deep = j


		horiz_dist = haversine(shall_slab[ind_shall,0],shall_slab[ind_shall,1],deep_slab[ind_deep,0],deep_slab[ind_deep,1])
		slab_dip = np.rad2deg(np.arctan((z_shall - z_deep)/horiz_dist))

		###          (trench longitude , trench latitude , dip angle)
		dips_tmp[n,:] = (trench_loc[ind_trench,0], trench_loc[ind_trench,1], slab_dip, dang1, dang2, dlat);

		if dang1 >= 3 or dang2 >= 3 or dlat >= 0.1:
			num_bad = num_bad + 1

# get rid of bad points
dips=np.zeros([len(lats)-num_bad,3]); ind=0
for n in range(len(lats)):

	if dips_tmp[n,3] >= 3 or dips_tmp[n,4] >= 3 or dips_tmp[n,5] >= 0.1:
		pass
	else:
		dips[ind,0] = dips_tmp[n,0]
		dips[ind,1] = dips_tmp[n,1]
		dips[ind,2] = dips_tmp[n,2]
		ind = ind + 1


if smooth==1:

	def smooth(y, box_pts):
		box = np.ones(box_pts)/box_pts
		y_smooth = np.convolve(y, box, mode='same')
		return y_smooth

	smoothed_dips=smooth(dips[:,2],box_width)
	dips[box_width:(np.size(dips,0)-box_width),2]=smoothed_dips[box_width:(np.size(dips,0)-box_width)]


# get dips to be every dip_int km 
new_dips = np.empty((0,3), int)
ind_interest=0
for k in range(10000):
	net_dist = 0
	min_error = 1e9
	for i in range(ind_interest,len(dips)-1):
		net_dist = net_dist + haversine(dips[i,0],dips[i,1],dips[i+1,0],dips[i+1,1])
		if np.abs(net_dist - dip_int) < min_error:
			min_error = np.abs(net_dist - dip_int) 
			ind_next  = i+1
			dist_best = net_dist

	if ind_next == len(dips) - 1:
		if dist_best > (dip_int/5.):
			new_dips = np.vstack((  new_dips,  dips[ind_next,:]  ))
		break
	else:
		new_dips = np.vstack((  new_dips,  dips[ind_next,:]  ))

	ind_interest = ind_next

np.savetxt('tmp.txt',new_dips,fmt='%.4f')
np.savetxt('tmp_full.txt',dips,fmt='%.4f')
