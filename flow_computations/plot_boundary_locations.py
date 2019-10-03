#!/usr/bin/python
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import math, sys, os
import glob


plate_model=str(sys.argv[1])	            
refframe=str(sys.argv[2])

trenchfile =''.join(['trench_motions/tnew.',refframe,'.dat']);

# read boundary info
infile  =''.join(['inputs/Subbon_',plate_model,'.inp']);
b=open(infile,"r").readlines()

num_bounds = sum(1 for line in b)
bounds = np.zeros([num_bounds,7])
for i in range(num_bounds):
	bounds[i,0] = int(b[i].split('\t')[0])  # iwall
	bounds[i,1] = float(b[i].split('\t')[1]) # lona
	bounds[i,2] = float(b[i].split('\t')[2]) # lata
	bounds[i,3] = float(b[i].split('\t')[3]) # lonb
	bounds[i,4] = float(b[i].split('\t')[4]) # latb
	if bounds[i,0] == 1:
		bounds[i,5] = float(b[i].split('\t')[7]) # vt_ew
		bounds[i,6] = float(b[i].split('\t')[8]) # vt_ns


t=open(trenchfile,"r").readlines()
num_vts = sum(1 for line in t)
trenches = np.zeros([num_vts,4])
for v in range(num_vts):
	trenches[v,0] = float(t[v].split(' ')[0]) # lon
	trenches[v,1] = float(t[v].split(' ')[1]) # lat
	trenches[v,2] = float(t[v].split(' ')[3]) # v_ew
	trenches[v,3] = float(t[v].split(' ')[4]) # v_ns

def midpoint(latA, lonA, latB, lonB):
	lonA = math.radians(lonA); lonB = math.radians(lonB)
	latA = math.radians(latA); latB = math.radians(latB)
	dLon = lonB - lonA
	Bx = math.cos(latB) * math.cos(dLon); By = math.cos(latB) * math.sin(dLon)
	latC = math.atan2(math.sin(latA) + math.sin(latB),math.sqrt((math.cos(latA) + Bx) * (math.cos(latA) + Bx) + By * By))
	lonC = lonA + math.atan2(By, math.cos(latA) + Bx)
	lonC = (lonC + 3 * math.pi) % (2 * math.pi) - math.pi
	return math.degrees(latC), math.degrees(lonC)

def haversine(lon1, lat1, lon2, lat2):
	lon1 = math.radians(lon1); lat1 = math.radians(lat1)
	lon2 = math.radians(lon2); lat2 = math.radians(lat2)
	dlon = lon2 - lon1; dlat = lat2 - lat1 
	a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
	c = 2 * math.asin(math.sqrt(a)) 
	return c 

for j in range(num_bounds):

	if bounds[j,1] < 0:
		bounds[j,1] = 360. + bounds[j,1]
	if bounds[j,3] < 0:
		bounds[j,3] = 360. + bounds[j,3]

odir_name =''.join(['plots/ages_and_boundary_tests/bounds/',plate_model,'/']);
if not os.path.exists(odir_name):
	os.mkdir(odir_name)

for i in range(num_bounds):

	plt.close(); plt.clf()

	print "segment %.0f" % (i+1)
	iwall = int(bounds[i,0])
	lona = bounds[i,1]
	lata = bounds[i,2]
	lonb = bounds[i,3]
	latb = bounds[i,4]
	lato, lono = midpoint(lata,lona,latb,lonb)

	dist_a_max = 1000.0;
	dist_b_max = 1000.0;
	dist_o_max = 1000.0;

	for v in range(num_vts):
		lon_vt = trenches[v,0]
		lat_vt = trenches[v,1]
		dista = haversine(lona, lata, lon_vt, lat_vt)
		if dista < dist_a_max:
			dist_a_max = dista;
			a_ind = v;
			lon_vta = lon_vt; lat_vta = lat_vt
		distb = haversine(lonb, latb, lon_vt, lat_vt)
		if distb < dist_b_max:
			dist_b_max = distb;
			b_ind = v;
			lon_vtb = lon_vt; lat_vtb = lat_vt
		disto = haversine(lono, lato, lon_vt, lat_vt)
		if disto < dist_o_max:
			dist_o_max = disto;
			o_ind = v;
			lon_vto = lon_vt; lat_vto = lat_vt

	vew_a = trenches[a_ind,2];  vns_a = trenches[a_ind,3]
	vew_b = trenches[b_ind,2];  vns_b = trenches[b_ind,3]
	vew_o = trenches[o_ind,2];  vns_o = trenches[o_ind,3]
	vew_avg = (vew_a + vew_o + vew_b)/3.
	vns_avg = (vns_a + vns_o + vns_b)/3.

	vew_a = round(vew_a, 2); 	vns_a = round(vns_a, 2)
	vew_b = round(vew_b, 2); 	vns_b = round(vns_b, 2)
	vew_o = round(vew_o, 2); 	vns_o = round(vns_o, 2)
	vew_avg = round(vew_avg,2); 
	vns_avg = round(vns_avg,2); 

	fig_name =''.join([odir_name,refframe,'_bound',str(i+1),'.png']);

	fig = plt.figure()
	map = Basemap(projection='hammer',lon_0=180,resolution='l')
	map.drawmeridians(np.arange(0,360,30))
	map.drawparallels(np.arange(-90,90,30))	
	map.drawcoastlines(linewidth=0.2)

	# plate bounds
	for d in range(num_bounds):
		if (330 < bounds[d,1] <= 360) and (0 < bounds[d,3] <= 30):
			pass;
		elif (330 < bounds[d,3] <= 360) and (0 < bounds[d,1] <= 30):
			pass;
		else: 
			if int(bounds[d,0]) == 1 and np.sqrt(bounds[d,5]**2 + bounds[d,6]**2) < 0.025:
				map.drawgreatcircle( bounds[d,1],bounds[d,2],bounds[d,3],bounds[d,4],linewidth=1,color='yellow')
			elif int(bounds[d,0]) == 1:
				map.drawgreatcircle( bounds[d,1],bounds[d,2],bounds[d,3],bounds[d,4],linewidth=1,color='pink')
			else:
				map.drawgreatcircle( bounds[d,1],bounds[d,2],bounds[d,3],bounds[d,4],linewidth=1,color='gray')
	if iwall == 1:
		map.drawgreatcircle( lona,lata,lonb,latb,linewidth=2,color='red')
	else:
		map.drawgreatcircle( lona,lata,lonb,latb,linewidth=2,color='black')
	
	x, y = map(lona,lata)
	map.plot(x, y, 'go', markersize=3.5)
	x, y = map(lonb,latb)
	map.plot(x, y, 'bo', markersize=3.5)
	x, y = map(lono,lato)
	map.plot(x, y, 'ro', markersize=3.5)

	x, y = map(lon_vta,lat_vta)
	map.plot(x, y, 'go', markersize=2)
	x, y = map(lon_vtb,lat_vtb)
	map.plot(x, y, 'bo', markersize=2)
	x, y = map(lon_vto,lat_vto)
	map.plot(x, y, 'ro', markersize=2)

	plt.annotate(''.join(['vew = ',str(vew_a),', vns = ',str(vns_a),' mm/yr ']), xy=(0, 1.25), xycoords='axes fraction',color='green',size=10)
	plt.annotate(''.join(['vew = ',str(vew_b),', vns = ',str(vns_b),' mm/yr ']), xy=(0, 1.175), xycoords='axes fraction',color='blue',size=10)
	plt.annotate(''.join(['vew = ',str(vew_o),', vns = ',str(vns_o),' mm/yr ']), xy=(0, 1.1), xycoords='axes fraction',color='red',size=10)
	plt.annotate(''.join(['vew = ',str(vew_avg),', vns = ',str(vns_avg),' mm/yr ']), xy=(0, 1.025), xycoords='axes fraction',color='black',size=10)
	plt.annotate(''.join(['(',str(lona),',',str(lata),')  to  (',str(lonb),',',str(latb),')']), xy=(0.45, 1.1), xycoords='axes fraction',color='black',size=10)
	plt.annotate(''.join(['segment ',str(i + 1),':']), xy=(0.55, 1.2), xycoords='axes fraction',color='black',size=12)

	plt.savefig(fig_name, bbox_inches='tight', dpi = 200)






