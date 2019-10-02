#!/usr/bin/python
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import subprocess
import os
import matplotlib.gridspec as gridspec
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

plt.ioff()

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
	'''
	Function to offset the "center" of a colormap
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


def plot_pressure_components(lons_out,lats_out,P_out,Pwall_out,Pedge_out,lona,lata,lonb,latb,lono,lato,iwall,DP,vt_ew,vt_ns,polygon_points,avgvel_asthen_ew,avgvel_asthen_ns,plot_name):


	press_plot_name=''.join([plot_name,'.pdf']);
	press_plot_name_png=''.join([plot_name,'.png']);
	press_plot_name_eps=''.join([plot_name,'.eps']);

	j_split = 10000;
	lons_plot = np.zeros((len(lons_out[0,:])))
	for j in range(0,len(lons_out[0,:])):
		if lons_out[0,j] > 180.:
			lons_plot[j] = lons_out[0,j] - 360.
			if j < j_split:
				j_split = j
		else:
			lons_plot[j] = lons_out[0,j]
	lons_plot=np.concatenate(( lons_plot[j_split:len(lons_plot)], lons_plot[0:j_split] ), axis=0)
	lats_plot=lats_out[:,0]		

	fig = plt.figure()
	ax = fig.add_subplot(311)
	map = Basemap(projection='hammer',lon_0=180,resolution='l')
	map.drawmeridians(np.arange(0,360,45),linewidth=0.15)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.15)	

	# Panel 1: plot full pressure
	x, y = map(lons_out, lats_out)
	cs = map.contourf(x,y,P_out,levels=np.linspace(-45,45,91),cmap=cm.get_cmap('bwr'),extend="both")
	cbar1 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar1.set_ticks(np.array([-40,-20,0,20,40]))
	cbar1.ax.tick_params(labelsize=6)
	# plot plate bounds
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
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red',zorder=1)
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black',zorder=1)
	
	# avg asthen vel
	vel_plot_thresh = 12.5
	avgvel_asthen_ew_plot=np.concatenate(( avgvel_asthen_ew[:,j_split:len(lons_plot)], avgvel_asthen_ew[:,0:j_split]), axis=1) 
	avgvel_asthen_ns_plot=np.concatenate(( avgvel_asthen_ns[:,j_split:len(lons_plot)], avgvel_asthen_ns[:,0:j_split]), axis=1) 
	avgvel_asthen_ew_proj,avgvel_asthen_ns_proj,xx,yy = map.transform_vector(avgvel_asthen_ew_plot, avgvel_asthen_ns_plot, lons_plot, lats_plot, 24, 12, returnxy=True)
	for i in range(0,avgvel_asthen_ew_proj.shape[0]):
		for j in range(0,avgvel_asthen_ew_proj.shape[1]):
			if np.sqrt(avgvel_asthen_ew_proj[i,j]**2 + avgvel_asthen_ns_proj[i,j]**2) > vel_plot_thresh or (np.sqrt(avgvel_asthen_ew_proj[i,j]**2 + avgvel_asthen_ns_proj[i,j]**2)) < 0.15:
				avgvel_asthen_ew_proj[i,j] = float('nan'); avgvel_asthen_ns_proj[i,j] = float('nan')
	J = map.quiver(xx,yy,avgvel_asthen_ew_proj,avgvel_asthen_ns_proj,scale=200, width=0.0025,zorder=2,color='slategray')
	plt.quiverkey(J, 0.96, 0.0, 10, '10 cm/yr', labelpos='W',fontproperties={'size': '7'})

	# plot DP points
	iwall = np.array(iwall); 
	lono = np.array(lono); lato = np.array(lato);
	wall_inds = np.where(iwall == 1)
	lono = lono[wall_inds]; lato = lato[wall_inds]
	vt_ew_temp = vt_ew[wall_inds]; vt_ns_temp = vt_ns[wall_inds]
	DP_walls = -1.0 * DP[:,4][wall_inds]
	for j in range(0,len(lono)):
		if np.sqrt(vt_ew_temp[j]**2 + vt_ns_temp[j]**2) > 0.05:
			xo, yo = map(lono[j],lato[j])
			if DP_walls[j] > 0:
				dps = map.scatter(xo, yo, c=DP_walls[j],s=17.5,cmap=cm.get_cmap('CMRmap_r'),vmin=-60,vmax=0,edgecolors='white',linewidth = 0.25,zorder=3)
			else:
				dps = map.scatter(xo, yo, c=DP_walls[j],s=17.5,cmap=cm.get_cmap('CMRmap_r'),vmin=-60,vmax=0,lw = 0,zorder=3)
			cbar2 = map.colorbar(dps,location='right',pad="23%",size="2%")
			cbar2.set_ticks(np.array([-60,-40,-20, 0]))
			cbar2.ax.tick_params(labelsize=6)
			cbar2.set_label("$\Delta$P  [MPa]",size=8)

	# Panel 2: plot Pedge pressure component
	ax2 = fig.add_subplot(312)
	map.drawmeridians(np.arange(0,360,45),linewidth=0.15)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.15)	
	cs = map.contourf(x,y,Pedge_out,levels=np.linspace(-45,45,91),cmap=cm.get_cmap('bwr'),extend="both")
	cbar3 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar3.set_ticks(np.array([-40,-20,0,20,40]))
	cbar3.ax.tick_params(labelsize=6)
	# plot plate bounds
	for i in range (0,len(lata)):
		if iwall[i] == 1:

			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown')	
				else:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red')
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black')
	# this colorbar just to keep subplots in proportion (delete w/illustrator in postprocessing)
	cbar4 = map.colorbar(cs,location='right',pad="23%",size="2%")
	cbar4.set_ticks(np.array([-40,-20,0,20,40]))
	cbar4.ax.tick_params(labelsize=6)
	cbar4.set_label("P  [MPa]",size=8)

	# Panel 3: plot Pwall pressure component
	ax3 = fig.add_subplot(313)
	map.drawmeridians(np.arange(0,360,45),linewidth=0.15)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.15)	
	cs = map.contourf(x,y,Pwall_out,levels=np.linspace(-45,45,91),cmap=cm.get_cmap('bwr'),extend="both")
	cbar5 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar5.set_ticks(np.array([-40,-20,0,20,40]))
	cbar5.ax.tick_params(labelsize=6)
	# plot plate bounds
	for i in range (0,len(lata)):
		if iwall[i] == 1:

			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown')	
				else:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red')
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black')
	# this colorbar just to keep subplots in proportion (delete w/illustrator in postprocessing)
	cbar6 = map.colorbar(cs,location='right',pad="23%",size="2%")
	cbar6.set_ticks(np.array([-40,-20,0,20,40]))
	cbar6.ax.tick_params(labelsize=6)
	cbar6.set_label("P  [MPa]",size=8)

	# finalize plot
	bash_command = ''.join(['convert -density 400 -flatten ',press_plot_name,' ',press_plot_name_png]);
	plt.savefig(press_plot_name, bbox_inches='tight', format='pdf')
	plt.savefig(press_plot_name_eps, bbox_inches='tight', format='eps')
	process = subprocess.Popen(['/bin/bash','-c',bash_command])
	process.wait()


def plot_pressure_simple(lons_out,lats_out,P_out,Pwall_out,Pedge_out,avgvel_Pdriven_ew,avgvel_Pdriven_ns,avgvel_total_ew,avgvel_total_ns,plate_vel_ew,plate_vel_ns,\
	trench_vels,pole_top_lon,pole_top_lat,pole_top_rate,lona,lata,lonb,latb,lono,lato,iwall,DP,vt_ew,vt_ns,polygon_points,age_dips,flux_vel,plot_name):


	press_plot_name=''.join([plot_name,'.pdf']);
	press_plot_name_png=''.join([plot_name,'.png']);
	press_plot_name_eps=''.join([plot_name,'.eps']);

	j_split = 10000;
	lons_plot = np.zeros((len(lons_out[0,:])))
	for j in range(0,len(lons_out[0,:])):
		if lons_out[0,j] > 180.:
			lons_plot[j] = lons_out[0,j] - 360.
			if j < j_split:
				j_split = j
		else:
			lons_plot[j] = lons_out[0,j]
	lons_plot=np.concatenate(( lons_plot[j_split:len(lons_plot)], lons_plot[0:j_split] ), axis=0)
	lats_plot=lats_out[:,0]		


	fig = plt.figure()
	ax = fig.add_subplot(311)
	map = Basemap(projection='hammer',lon_0=180,resolution='c')
	map.drawcoastlines(linewidth=0.3,color='darkgray',zorder=1)
	x, y = map(lons_out, lats_out)
	cs = map.contourf(x,y,P_out,levels=np.linspace(-45,45,101),cmap=cm.get_cmap('bwr'),extend="both")
	cbar2 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar2.set_ticks(np.array([-40,-20,0,20,40]))
	cbar2.ax.tick_params(labelsize=6)
	# plot plate bounds
	for i in range (0,len(lata)):
		if iwall[i] == 1:

			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8, color='brown')	
				else:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.9,color='red')
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black')

	# finalize plot
	bash_command = ''.join(['convert -density 400 -flatten ',press_plot_name,' ',press_plot_name_png]);
	plt.savefig(press_plot_name, bbox_inches='tight', format='pdf')
	plt.savefig(press_plot_name_eps, bbox_inches='tight', format='eps')
	process = subprocess.Popen(['/bin/bash','-c',bash_command])
	process.wait()


def plot_plate_indices(lons_out,lats_out,lona,lata,lonb,latb,lono,polygon_points,plot_name):

    index_plot_name=''.join([plot_name,'.pdf']);

    plt.clf()
    fig = plt.figure()
    map = Basemap(projection='hammer',lon_0=180,resolution='l')
    map.drawmeridians(np.arange(0,360,45),linewidth=0.3)
    map.drawparallels(np.arange(-90,90,45),linewidth=0.3)
    x, y = map(lons_out, lats_out)
    cs = map.contourf(x,y,polygon_points,levels=np.linspace(0,18,19))
    cbar = map.colorbar(cs,location='right',pad="5%")
    cbar.set_label('plate index')
    # plot plate bounds
    for i in range (0,len(lata)):
        if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
            pass;
        elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
            pass;
        else:
            map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.5,color='black')

	# finalize plot
    plt.savefig(index_plot_name, bbox_inches='tight', format='pdf')

def plot_pressure_components_MoreVelFields(lons_out,lats_out,P_out,Pwall_out,Pedge_out,lona,lata,lonb,latb,lono,lato,iwall,DP,vt_ew,vt_ns,polygon_points,\
	avgvel_asthen_ew,avgvel_asthen_ns,pdrivenvel_wall_ew, pdrivenvel_wall_ns, pdrivenvel_edge_ew, pdrivenvel_edge_ns,plot_name):


	press_plot_name=''.join([plot_name,'.pdf']);
	press_plot_name_png=''.join([plot_name,'.png']);
	press_plot_name_eps=''.join([plot_name,'.eps']);

	j_split = 10000;
	lons_plot = np.zeros((len(lons_out[0,:])))
	for j in range(0,len(lons_out[0,:])):
		if lons_out[0,j] > 180.:
			lons_plot[j] = lons_out[0,j] - 360.
			if j < j_split:
				j_split = j
		else:
			lons_plot[j] = lons_out[0,j]
	lons_plot=np.concatenate(( lons_plot[j_split:len(lons_plot)], lons_plot[0:j_split] ), axis=0)
	lats_plot=lats_out[:,0]		

	fig = plt.figure()
	ax = fig.add_subplot(311)
	map = Basemap(projection='hammer',lon_0=180,resolution='l')
	map.drawmeridians(np.arange(0,360,45),linewidth=0.1)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.1)	

	# Panel 1: plot full pressure
	x, y = map(lons_out, lats_out)
	cs = map.contourf(x,y,P_out,levels=np.linspace(-60,60,121),cmap=cm.get_cmap('bwr'),extend="both")
        cbar1 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar1.set_ticks(np.array([-60,-40,-20,0,20,40,60]))
        cbar1.ax.tick_params(labelsize=6)
	# plot plate bounds
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
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red',zorder=1)
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black',zorder=1)
	
	# avg asthen vel
	vel_plot_thresh = 11
	avgvel_asthen_ew_plot=np.concatenate(( avgvel_asthen_ew[:,j_split:len(lons_plot)], avgvel_asthen_ew[:,0:j_split]), axis=1) 
	avgvel_asthen_ns_plot=np.concatenate(( avgvel_asthen_ns[:,j_split:len(lons_plot)], avgvel_asthen_ns[:,0:j_split]), axis=1) 
	avgvel_asthen_ew_proj,avgvel_asthen_ns_proj,xx,yy = map.transform_vector(avgvel_asthen_ew_plot, avgvel_asthen_ns_plot, lons_plot, lats_plot, 24, 12, returnxy=True)
	for i in range(0,avgvel_asthen_ew_proj.shape[0]):
		for j in range(0,avgvel_asthen_ew_proj.shape[1]):
			if np.sqrt(avgvel_asthen_ew_proj[i,j]**2 + avgvel_asthen_ns_proj[i,j]**2) > vel_plot_thresh or (np.sqrt(avgvel_asthen_ew_proj[i,j]**2 + avgvel_asthen_ns_proj[i,j]**2)) < 0.15:
				avgvel_asthen_ew_proj[i,j] = float('nan'); avgvel_asthen_ns_proj[i,j] = float('nan')
	J = map.quiver(xx,yy,avgvel_asthen_ew_proj,avgvel_asthen_ns_proj,scale=200, width=0.0025,zorder=2,color='slategray')
	plt.quiverkey(J, 0.96, 0.0, 10, '10 cm/yr', labelpos='W',fontproperties={'size': '7'})

	# plot DP points
	iwall = np.array(iwall); 
	lono = np.array(lono); lato = np.array(lato);
	wall_inds = np.where(iwall == 1)
	lono = lono[wall_inds]; lato = lato[wall_inds]
	vt_ew_temp = vt_ew[wall_inds]; vt_ns_temp = vt_ns[wall_inds]
	DP_walls = -1.0 * DP[:,4][wall_inds]
	for j in range(0,len(lono)):
		if np.sqrt(vt_ew_temp[j]**2 + vt_ns_temp[j]**2) > 0.05:
			xo, yo = map(lono[j],lato[j])
			if DP_walls[j] > 0:
				dps = map.scatter(xo, yo, c=DP_walls[j],s=17.5,cmap=cm.get_cmap('CMRmap_r'),vmin=-80,vmax=0,edgecolors='white',linewidth = 0.25,zorder=3)
			else:
				dps = map.scatter(xo, yo, c=DP_walls[j],s=17.5,cmap=cm.get_cmap('CMRmap_r'),vmin=-80,vmax=0,lw = 0,zorder=3)
			cbar2 = map.colorbar(dps,location='right',pad="23%",size="2%",extend='both')
			cbar2.set_ticks(np.array([-80,-60,-40,-20, 0]))
			cbar2.ax.tick_params(labelsize=6)
			cbar2.set_label("$\Delta$P  [MPa]",size=8)

	# Panel 2: plot Pedge pressure component
	ax2 = fig.add_subplot(312)
	map.drawmeridians(np.arange(0,360,45),linewidth=0.1)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.1)	
	cs = map.contourf(x,y,Pedge_out,levels=np.linspace(-60,60,121),cmap=cm.get_cmap('bwr'),extend="both")
	cbar3 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar3.set_ticks(np.array([-60,-40,-20,0,20,40,60]))
	cbar3.ax.tick_params(labelsize=6)
	# plot plate bounds
	for i in range (0,len(lata)):
		if iwall[i] == 1:

			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown')	
				else:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red')
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black')
	# this colorbar just to keep subplots in proportion (delete w/illustrator in postprocessing)
	cbar4 = map.colorbar(cs,location='right',pad="23%",size="2%")
	cbar4.set_ticks(np.array([-40,-20,0,20,40]))
	cbar4.ax.tick_params(labelsize=6)
	cbar4.set_label("P  [MPa]",size=8)
	# vel vectors
	pdrivenvel_edge_ew_plot=np.concatenate(( pdrivenvel_edge_ew[:,j_split:len(lons_plot)], pdrivenvel_edge_ew[:,0:j_split]), axis=1) 
	pdrivenvel_edge_ns_plot=np.concatenate(( pdrivenvel_edge_ns[:,j_split:len(lons_plot)], pdrivenvel_edge_ns[:,0:j_split]), axis=1) 
	pdrivenvel_edge_ew_proj,pdrivenvel_edge_ns_proj,xx,yy = map.transform_vector(pdrivenvel_edge_ew_plot, pdrivenvel_edge_ns_plot, lons_plot, lats_plot, 24, 12, returnxy=True)
	for i in range(0,pdrivenvel_edge_ew_proj.shape[0]):
		for j in range(0,pdrivenvel_edge_ew_proj.shape[1]):
			if np.sqrt(pdrivenvel_edge_ew_proj[i,j]**2 + pdrivenvel_edge_ns_proj[i,j]**2) > vel_plot_thresh or (np.sqrt(pdrivenvel_edge_ew_proj[i,j]**2 + pdrivenvel_edge_ns_proj[i,j]**2)) < 0.15:
				pdrivenvel_edge_ew_proj[i,j] = float('nan'); pdrivenvel_edge_ns_proj[i,j] = float('nan')
	J = map.quiver(xx,yy,pdrivenvel_edge_ew_proj,pdrivenvel_edge_ns_proj,scale=200, width=0.0025,zorder=2,color='slategray')
	plt.quiverkey(J, 0.96, 0.0, 10, '10 cm/yr', labelpos='W',fontproperties={'size': '7'})


	# Panel 3: plot Pwall pressure component
	ax3 = fig.add_subplot(313)
	map.drawmeridians(np.arange(0,360,45),linewidth=0.1)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.1)	
	cs = map.contourf(x,y,Pwall_out,levels=np.linspace(-60,60,121),cmap=cm.get_cmap('bwr'),extend="both")
	cbar5 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar5.set_ticks(np.array([-60,-40,-20,0,20,40,60]))
	cbar5.ax.tick_params(labelsize=6)
	# plot plate bounds
	for i in range (0,len(lata)):
		if iwall[i] == 1:

			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown')	
				else:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red')
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black')
	# this colorbar just to keep subplots in proportion (delete w/illustrator in postprocessing)
	cbar6 = map.colorbar(cs,location='right',pad="23%",size="2%")
	cbar6.set_ticks(np.array([-60,-40,-20,0,20,40,60]))
	cbar6.ax.tick_params(labelsize=6)
	cbar6.set_label("P  [MPa]",size=8)
	# vel vectors
	pdrivenvel_wall_ew_plot=np.concatenate(( pdrivenvel_wall_ew[:,j_split:len(lons_plot)], pdrivenvel_wall_ew[:,0:j_split]), axis=1) 
	pdrivenvel_wall_ns_plot=np.concatenate(( pdrivenvel_wall_ns[:,j_split:len(lons_plot)], pdrivenvel_wall_ns[:,0:j_split]), axis=1) 
	pdrivenvel_wall_ew_proj,pdrivenvel_wall_ns_proj,xx,yy = map.transform_vector(pdrivenvel_wall_ew_plot, pdrivenvel_wall_ns_plot, lons_plot, lats_plot, 24, 12, returnxy=True)
	for i in range(0,pdrivenvel_wall_ew_proj.shape[0]):
		for j in range(0,pdrivenvel_wall_ew_proj.shape[1]):
			if np.sqrt(pdrivenvel_wall_ew_proj[i,j]**2 + pdrivenvel_wall_ns_proj[i,j]**2) > vel_plot_thresh or (np.sqrt(pdrivenvel_wall_ew_proj[i,j]**2 + pdrivenvel_wall_ns_proj[i,j]**2)) < 0.15:
				pdrivenvel_wall_ew_proj[i,j] = float('nan'); pdrivenvel_wall_ns_proj[i,j] = float('nan')
	J = map.quiver(xx,yy,pdrivenvel_wall_ew_proj,pdrivenvel_wall_ns_proj,scale=200, width=0.0025,zorder=2,color='slategray')
	plt.quiverkey(J, 0.96, 0.0, 10, '10 cm/yr', labelpos='W',fontproperties={'size': '7'})

	# finalize plot
	bash_command = ''.join(['convert -density 400 -flatten ',press_plot_name,' ',press_plot_name_png]);
	plt.savefig(press_plot_name, bbox_inches='tight', format='pdf')
	plt.savefig(press_plot_name_eps, bbox_inches='tight', format='eps')
	process = subprocess.Popen(['/bin/bash','-c',bash_command])
	process.wait()



def plot_pressure_components_WPacPressPoint(lons_out,lats_out,P_out,Pwall_out,Pedge_out,lona,lata,lonb,latb,lono,lato,iwall,DP,vt_ew,vt_ns,polygon_points,avgvel_asthen_ew,avgvel_asthen_ns,plot_name,lon_pac,lat_pac,Press_pac):


	press_plot_name=''.join([plot_name,'.pdf']);
	press_plot_name_png=''.join([plot_name,'.png']);
	press_plot_name_eps=''.join([plot_name,'.eps']);

	j_split = 10000;
	lons_plot = np.zeros((len(lons_out[0,:])))
	for j in range(0,len(lons_out[0,:])):
		if lons_out[0,j] > 180.:
			lons_plot[j] = lons_out[0,j] - 360.
			if j < j_split:
				j_split = j
		else:
			lons_plot[j] = lons_out[0,j]
	lons_plot=np.concatenate(( lons_plot[j_split:len(lons_plot)], lons_plot[0:j_split] ), axis=0)
	lats_plot=lats_out[:,0]		

	fig = plt.figure()
	ax = fig.add_subplot(311)
	map = Basemap(projection='hammer',lon_0=180,resolution='l')
	map.drawmeridians(np.arange(0,360,45),linewidth=0.15)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.15)	

	# Panel 1: plot full pressure
	x, y = map(lons_out, lats_out)
	cs = map.contourf(x,y,P_out,levels=np.linspace(-45,45,91),cmap=cm.get_cmap('bwr'),extend="both")
	cbar1 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar1.set_ticks(np.array([-40,-20,0,20,40]))
	cbar1.ax.tick_params(labelsize=6)
	# plot plate bounds
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
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red',zorder=1)
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black',zorder=1)
	
	# avg asthen vel
	vel_plot_thresh = 12.5
	avgvel_asthen_ew_plot=np.concatenate(( avgvel_asthen_ew[:,j_split:len(lons_plot)], avgvel_asthen_ew[:,0:j_split]), axis=1) 
	avgvel_asthen_ns_plot=np.concatenate(( avgvel_asthen_ns[:,j_split:len(lons_plot)], avgvel_asthen_ns[:,0:j_split]), axis=1) 
	avgvel_asthen_ew_proj,avgvel_asthen_ns_proj,xx,yy = map.transform_vector(avgvel_asthen_ew_plot, avgvel_asthen_ns_plot, lons_plot, lats_plot, 24, 12, returnxy=True)
	for i in range(0,avgvel_asthen_ew_proj.shape[0]):
		for j in range(0,avgvel_asthen_ew_proj.shape[1]):
			if np.sqrt(avgvel_asthen_ew_proj[i,j]**2 + avgvel_asthen_ns_proj[i,j]**2) > vel_plot_thresh or (np.sqrt(avgvel_asthen_ew_proj[i,j]**2 + avgvel_asthen_ns_proj[i,j]**2)) < 0.15:
				avgvel_asthen_ew_proj[i,j] = float('nan'); avgvel_asthen_ns_proj[i,j] = float('nan')
	J = map.quiver(xx,yy,avgvel_asthen_ew_proj,avgvel_asthen_ns_proj,scale=200, width=0.0025,zorder=2,color='slategray')
	plt.quiverkey(J, 0.96, 0.0, 10, '10 cm/yr', labelpos='W',fontproperties={'size': '7'})

	# plot DP points
	iwall = np.array(iwall); 
	lono = np.array(lono); lato = np.array(lato);
	wall_inds = np.where(iwall == 1)
	lono = lono[wall_inds]; lato = lato[wall_inds]
	vt_ew_temp = vt_ew[wall_inds]; vt_ns_temp = vt_ns[wall_inds]
	DP_walls = -1.0 * DP[:,4][wall_inds]
	for j in range(0,len(lono)):
		if np.sqrt(vt_ew_temp[j]**2 + vt_ns_temp[j]**2) > 0.05:
			xo, yo = map(lono[j],lato[j])
			if DP_walls[j] > 0:
				dps = map.scatter(xo, yo, c=DP_walls[j],s=17.5,cmap=cm.get_cmap('CMRmap_r'),vmin=-60,vmax=0,edgecolors='white',linewidth = 0.25,zorder=3)
			else:
				dps = map.scatter(xo, yo, c=DP_walls[j],s=17.5,cmap=cm.get_cmap('CMRmap_r'),vmin=-60,vmax=0,lw = 0,zorder=3)
			cbar2 = map.colorbar(dps,location='right',pad="23%",size="2%")
			cbar2.set_ticks(np.array([-60,-40,-20, 0]))
			cbar2.ax.tick_params(labelsize=6)
			cbar2.set_label("$\Delta$P  [MPa]",size=8)
	xpac,ypac = map(lon_pac, lat_pac)
	map.plot(xpac, ypac, marker='*', color='green', markersize=15, zorder=4, linewidth=0, markeredgewidth=0.0)
	press_string = ''.join([str(str(round(Press_pac,2))),' MPa']);
	plt.annotate(press_string, xy=(0.5, 0.5), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=5.5)


	# Panel 2: plot Pedge pressure component
	ax2 = fig.add_subplot(312)
	map.drawmeridians(np.arange(0,360,45),linewidth=0.15)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.15)	
	cs = map.contourf(x,y,Pedge_out,levels=np.linspace(-45,45,91),cmap=cm.get_cmap('bwr'),extend="both")
	cbar3 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar3.set_ticks(np.array([-40,-20,0,20,40]))
	cbar3.ax.tick_params(labelsize=6)
	# plot plate bounds
	for i in range (0,len(lata)):
		if iwall[i] == 1:

			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown')	
				else:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red')
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black')
	# this colorbar just to keep subplots in proportion (delete w/illustrator in postprocessing)
	cbar4 = map.colorbar(cs,location='right',pad="23%",size="2%")
	cbar4.set_ticks(np.array([-40,-20,0,20,40]))
	cbar4.ax.tick_params(labelsize=6)
	cbar4.set_label("P  [MPa]",size=8)

	# Panel 3: plot Pwall pressure component
	ax3 = fig.add_subplot(313)
	map.drawmeridians(np.arange(0,360,45),linewidth=0.15)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.15)	
	cs = map.contourf(x,y,Pwall_out,levels=np.linspace(-45,45,91),cmap=cm.get_cmap('bwr'),extend="both")
	cbar5 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar5.set_ticks(np.array([-40,-20,0,20,40]))
	cbar5.ax.tick_params(labelsize=6)
	# plot plate bounds
	for i in range (0,len(lata)):
		if iwall[i] == 1:

			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown')	
				else:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red')
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black')
	# this colorbar just to keep subplots in proportion (delete w/illustrator in postprocessing)
	cbar6 = map.colorbar(cs,location='right',pad="23%",size="2%")
	cbar6.set_ticks(np.array([-40,-20,0,20,40]))
	cbar6.ax.tick_params(labelsize=6)
	cbar6.set_label("P  [MPa]",size=8)

	# finalize plot
	bash_command = ''.join(['convert -density 400 -flatten ',press_plot_name,' ',press_plot_name_png]);
	plt.savefig(press_plot_name, bbox_inches='tight', format='pdf')
	plt.savefig(press_plot_name_eps, bbox_inches='tight', format='eps')
	process = subprocess.Popen(['/bin/bash','-c',bash_command])
	process.wait()


def plot_pressure_components_MoreVelFields_wildColors(lons_out,lats_out,P_out,Pwall_out,Pedge_out,lona,lata,lonb,latb,lono,lato,iwall,DP,vt_ew,vt_ns,polygon_points,\
	avgvel_asthen_ew,avgvel_asthen_ns,pdrivenvel_wall_ew, pdrivenvel_wall_ns, pdrivenvel_edge_ew, pdrivenvel_edge_ns,plot_name):


	press_plot_name=''.join([plot_name,'.pdf']);
	press_plot_name_png=''.join([plot_name,'.png']);
	press_plot_name_eps=''.join([plot_name,'.eps']);

	j_split = 10000;
	lons_plot = np.zeros((len(lons_out[0,:])))
	for j in range(0,len(lons_out[0,:])):
		if lons_out[0,j] > 180.:
			lons_plot[j] = lons_out[0,j] - 360.
			if j < j_split:
				j_split = j
		else:
			lons_plot[j] = lons_out[0,j]
	lons_plot=np.concatenate(( lons_plot[j_split:len(lons_plot)], lons_plot[0:j_split] ), axis=0)
	lats_plot=lats_out[:,0]		

	fig = plt.figure()
	ax = fig.add_subplot(311)
	map = Basemap(projection='hammer',lon_0=180,resolution='l')
	map.drawmeridians(np.arange(0,360,45),linewidth=0.1)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.1)	

	# Panel 1: plot full pressure
	x, y = map(lons_out, lats_out)
	cs = map.contourf(x,y,P_out,levels=np.linspace(-60,60,121),cmap=cm.get_cmap('jet'),extend="both")
        cbar1 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar1.set_ticks(np.array([-60,-40,-20,0,20,40,60]))
        cbar1.ax.tick_params(labelsize=6)
	# plot plate bounds
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
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red',zorder=1)
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black',zorder=1)
	
	# avg asthen vel
	vel_plot_thresh = 11
	avgvel_asthen_ew_plot=np.concatenate(( avgvel_asthen_ew[:,j_split:len(lons_plot)], avgvel_asthen_ew[:,0:j_split]), axis=1) 
	avgvel_asthen_ns_plot=np.concatenate(( avgvel_asthen_ns[:,j_split:len(lons_plot)], avgvel_asthen_ns[:,0:j_split]), axis=1) 
	avgvel_asthen_ew_proj,avgvel_asthen_ns_proj,xx,yy = map.transform_vector(avgvel_asthen_ew_plot, avgvel_asthen_ns_plot, lons_plot, lats_plot, 24, 12, returnxy=True)
	for i in range(0,avgvel_asthen_ew_proj.shape[0]):
		for j in range(0,avgvel_asthen_ew_proj.shape[1]):
			if np.sqrt(avgvel_asthen_ew_proj[i,j]**2 + avgvel_asthen_ns_proj[i,j]**2) > vel_plot_thresh or (np.sqrt(avgvel_asthen_ew_proj[i,j]**2 + avgvel_asthen_ns_proj[i,j]**2)) < 0.15:
				avgvel_asthen_ew_proj[i,j] = float('nan'); avgvel_asthen_ns_proj[i,j] = float('nan')
	J = map.quiver(xx,yy,avgvel_asthen_ew_proj,avgvel_asthen_ns_proj,scale=200, width=0.0025,zorder=2,color='slategray')
	plt.quiverkey(J, 0.96, 0.0, 10, '10 cm/yr', labelpos='W',fontproperties={'size': '7'})

	# plot DP points
	iwall = np.array(iwall); 
	lono = np.array(lono); lato = np.array(lato);
	wall_inds = np.where(iwall == 1)
	lono = lono[wall_inds]; lato = lato[wall_inds]
	vt_ew_temp = vt_ew[wall_inds]; vt_ns_temp = vt_ns[wall_inds]
	DP_walls = -1.0 * DP[:,4][wall_inds]
	for j in range(0,len(lono)):
		if np.sqrt(vt_ew_temp[j]**2 + vt_ns_temp[j]**2) > 0.05:
			xo, yo = map(lono[j],lato[j])
			if DP_walls[j] > 0:
				dps = map.scatter(xo, yo, c=DP_walls[j],s=17.5,cmap=cm.get_cmap('CMRmap_r'),vmin=-80,vmax=0,edgecolors='white',linewidth = 0.25,zorder=3)
			else:
				dps = map.scatter(xo, yo, c=DP_walls[j],s=17.5,cmap=cm.get_cmap('CMRmap_r'),vmin=-80,vmax=0,lw = 0,zorder=3)
			cbar2 = map.colorbar(dps,location='right',pad="23%",size="2%",extend='both')
			cbar2.set_ticks(np.array([-80,-60,-40,-20, 0]))
			cbar2.ax.tick_params(labelsize=6)
			cbar2.set_label("$\Delta$P  [MPa]",size=8)

	# Panel 2: plot Pedge pressure component
	ax2 = fig.add_subplot(312)
	map.drawmeridians(np.arange(0,360,45),linewidth=0.1)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.1)	
	cs = map.contourf(x,y,Pedge_out,levels=np.linspace(-60,60,121),cmap=cm.get_cmap('jet'),extend="both")
	cbar3 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar3.set_ticks(np.array([-60,-40,-20,0,20,40,60]))
	cbar3.ax.tick_params(labelsize=6)
	# plot plate bounds
	for i in range (0,len(lata)):
		if iwall[i] == 1:

			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown')	
				else:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red')
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black')
	# this colorbar just to keep subplots in proportion (delete w/illustrator in postprocessing)
	cbar4 = map.colorbar(cs,location='right',pad="23%",size="2%")
	cbar4.set_ticks(np.array([-40,-20,0,20,40]))
	cbar4.ax.tick_params(labelsize=6)
	cbar4.set_label("P  [MPa]",size=8)
	# vel vectors
	pdrivenvel_edge_ew_plot=np.concatenate(( pdrivenvel_edge_ew[:,j_split:len(lons_plot)], pdrivenvel_edge_ew[:,0:j_split]), axis=1) 
	pdrivenvel_edge_ns_plot=np.concatenate(( pdrivenvel_edge_ns[:,j_split:len(lons_plot)], pdrivenvel_edge_ns[:,0:j_split]), axis=1) 
	pdrivenvel_edge_ew_proj,pdrivenvel_edge_ns_proj,xx,yy = map.transform_vector(pdrivenvel_edge_ew_plot, pdrivenvel_edge_ns_plot, lons_plot, lats_plot, 24, 12, returnxy=True)
	for i in range(0,pdrivenvel_edge_ew_proj.shape[0]):
		for j in range(0,pdrivenvel_edge_ew_proj.shape[1]):
			if np.sqrt(pdrivenvel_edge_ew_proj[i,j]**2 + pdrivenvel_edge_ns_proj[i,j]**2) > vel_plot_thresh or (np.sqrt(pdrivenvel_edge_ew_proj[i,j]**2 + pdrivenvel_edge_ns_proj[i,j]**2)) < 0.15:
				pdrivenvel_edge_ew_proj[i,j] = float('nan'); pdrivenvel_edge_ns_proj[i,j] = float('nan')
	J = map.quiver(xx,yy,pdrivenvel_edge_ew_proj,pdrivenvel_edge_ns_proj,scale=200, width=0.0025,zorder=2,color='slategray')
	plt.quiverkey(J, 0.96, 0.0, 10, '10 cm/yr', labelpos='W',fontproperties={'size': '7'})


	# Panel 3: plot Pwall pressure component
	ax3 = fig.add_subplot(313)
	map.drawmeridians(np.arange(0,360,45),linewidth=0.1)
	map.drawparallels(np.arange(-90,90,45),linewidth=0.1)	
	cs = map.contourf(x,y,Pwall_out,levels=np.linspace(-60,60,121),cmap=cm.get_cmap('jet'),extend="both")
	cbar5 = map.colorbar(cs,location='right',pad="-6%",size="2%")
	cbar5.set_ticks(np.array([-60,-40,-20,0,20,40,60]))
	cbar5.ax.tick_params(labelsize=6)
	# plot plate bounds
	for i in range (0,len(lata)):
		if iwall[i] == 1:

			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				if np.sqrt(vt_ew[i]**2 + vt_ns[i]**2) < 0.05:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6, color='brown')	
				else:
					map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.8,color='red')
		else:
			if (330 < lona[i] <= 360) and (0 < lonb[i] <= 30):
				pass;
			elif (330 < lonb[i] <= 360) and (0 < lona[i] <= 30):
				pass;
			else:
				map.drawgreatcircle( lona[i],lata[i],lonb[i],latb[i],linewidth=0.6,color='black')
	# this colorbar just to keep subplots in proportion (delete w/illustrator in postprocessing)
	cbar6 = map.colorbar(cs,location='right',pad="23%",size="2%")
	cbar6.set_ticks(np.array([-60,-40,-20,0,20,40,60]))
	cbar6.ax.tick_params(labelsize=6)
	cbar6.set_label("P  [MPa]",size=8)
	# vel vectors
	pdrivenvel_wall_ew_plot=np.concatenate(( pdrivenvel_wall_ew[:,j_split:len(lons_plot)], pdrivenvel_wall_ew[:,0:j_split]), axis=1) 
	pdrivenvel_wall_ns_plot=np.concatenate(( pdrivenvel_wall_ns[:,j_split:len(lons_plot)], pdrivenvel_wall_ns[:,0:j_split]), axis=1) 
	pdrivenvel_wall_ew_proj,pdrivenvel_wall_ns_proj,xx,yy = map.transform_vector(pdrivenvel_wall_ew_plot, pdrivenvel_wall_ns_plot, lons_plot, lats_plot, 24, 12, returnxy=True)
	for i in range(0,pdrivenvel_wall_ew_proj.shape[0]):
		for j in range(0,pdrivenvel_wall_ew_proj.shape[1]):
			if np.sqrt(pdrivenvel_wall_ew_proj[i,j]**2 + pdrivenvel_wall_ns_proj[i,j]**2) > vel_plot_thresh or (np.sqrt(pdrivenvel_wall_ew_proj[i,j]**2 + pdrivenvel_wall_ns_proj[i,j]**2)) < 0.15:
				pdrivenvel_wall_ew_proj[i,j] = float('nan'); pdrivenvel_wall_ns_proj[i,j] = float('nan')
	J = map.quiver(xx,yy,pdrivenvel_wall_ew_proj,pdrivenvel_wall_ns_proj,scale=200, width=0.0025,zorder=2,color='slategray')
	plt.quiverkey(J, 0.96, 0.0, 10, '10 cm/yr', labelpos='W',fontproperties={'size': '7'})

	# finalize plot
	bash_command = ''.join(['convert -density 400 -flatten ',press_plot_name,' ',press_plot_name_png]);
	plt.savefig(press_plot_name, bbox_inches='tight', format='pdf')
	plt.savefig(press_plot_name_eps, bbox_inches='tight', format='eps')
	process = subprocess.Popen(['/bin/bash','-c',bash_command])
	process.wait()

