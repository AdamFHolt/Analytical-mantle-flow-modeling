#!/usr/bin/python 

from pandas import *
import sys
import numpy as np
import gpxpy.geo
import geopy.distance
import math
from functions import haversine

dips_name=sys.argv[1]
dips_full_name=sys.argv[2]

dips=np.loadtxt(dips_name)
dips_full=np.loadtxt(dips_full_name)

lon1=283.90; lat1=-14.37	# SAm 1
lon2=286.76; lat2=-33.01	# SAm 2
lon3=143.07; lat3=25.33		# Mariana
rad_km=6378.
thresh_dist=200.

# every-so-often case
num = 0
for n in range(len(dips)):
	lon = dips[n,0]; lat = dips[n,1]
	dist1 = haversine(lon,lat,lon1,lat1,rad_km); 
	dist2 = haversine(lon,lat,lon2,lat2,rad_km); 
	dist3 = haversine(lon,lat,lon3,lat3,rad_km); 
	if dist1 < thresh_dist or dist2 < thresh_dist or dist3 < thresh_dist:
		pass
	else:
		num = num + 1

dips_new=np.zeros((num,3))
num = 0
for i in range(len(dips)):
	lon = dips[i,0]; lat = dips[i,1]
	dist1 = haversine(lon,lat,lon1,lat1,rad_km); 
	dist2 = haversine(lon,lat,lon2,lat2,rad_km); 
	dist3 = haversine(lon,lat,lon3,lat3,rad_km); 
	if dist1 < thresh_dist or dist2 < thresh_dist or dist3 < thresh_dist:
		pass
	else:
		dips_new[num,:] = dips[i,:]
		num = num + 1

## full case
num = 0
for n in range(len(dips_full)):
	lon = dips_full[n,0]; lat = dips_full[n,1]
	dist1 = haversine(lon,lat,lon1,lat1,rad_km); 
	dist2 = haversine(lon,lat,lon2,lat2,rad_km); 
	dist3 = haversine(lon,lat,lon3,lat3,rad_km); 
	if dist1 < thresh_dist or dist2 < thresh_dist or dist3 < thresh_dist:
		pass
	else:
		num = num + 1

dips_full_new=np.zeros((num,3))
num = 0
for i in range(len(dips_full)):
	lon = dips_full[i,0]; lat = dips_full[i,1]
	dist1 = haversine(lon,lat,lon1,lat1,rad_km); 
	dist2 = haversine(lon,lat,lon2,lat2,rad_km); 
	dist3 = haversine(lon,lat,lon3,lat3,rad_km); 
	if dist1 < thresh_dist or dist2 < thresh_dist or dist3 < thresh_dist:
		pass
	else:
		dips_full_new[num,:] = dips_full[i,:]
		num = num + 1

np.savetxt(dips_name,dips_new,fmt='%.4f')
np.savetxt(dips_full_name,dips_full_new,fmt='%.4f')
