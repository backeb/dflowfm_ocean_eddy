# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:37:21 2020
Updated on 

@purpose: 
    plot eddy evolution on map and position of max ssh

@author: 
    backeb
"""

"import libraries"
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

"define some variables for plotting"
expt  = 'expt14'
start = 0
stop  = 120

if expt == 'expt14':
    #set filename
    fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-'+expt+'\oceaneddymankmx0_map.nc'
    days2plt = np.linspace(start, stop, num = 4, endpoint = True)
    
    ugrid_all = get_netdata(file_nc=fname)
    ds = xr.open_dataset(fname)
    
    x = np.empty(ds.time.size)
    y = np.empty(ds.time.size)
    
    for i in np.arange(0, ds.time.size, 1):
        ssh = ds.s1.data[i,:]
        maxi = np.argmax(ssh)
        x[i], y[i] = ugrid_all.verts[maxi,:,:].mean(axis = 0)
    
    cnt = 0
    fig, axs = plt.subplots(1, 4, sharey=True, figsize=(15, 5))
    
    for i in days2plt.astype(int):
    
        #plot water level on map
        ssh = get_ncmodeldata(file_nc=fname, 
                              varname='s1', 
                              timestep=i)
        
        pc = plot_netmapdata(ugrid_all.verts, values=ssh[0,:], ax=axs[cnt], linewidth=0.5, cmap="jet")
        pc.set_clim([0, 0.025])
        axs[cnt].plot(x[:i],y[:i],'r.', markersize = 2)
        #axs[cnt].set_title(np.datetime_as_string(ds.time.data[i], unit = 'h'))
        axs[cnt].set_title('t = '+str(i)+' days')
        axs[cnt].set_aspect('equal')
        xticks = np.linspace(np.min(ds.NetNode_x.data),
                         np.max(ds.NetNode_x.data),
                         num = 5,
                         endpoint = True)
        axs[cnt].set_xticks(xticks)
        axs[cnt].set_xlabel('%s (%s)'%(ds.NetNode_x.long_name, ds.NetNode_x.units))
        if cnt == 0:
            yticks = np.linspace(np.min(ds.NetNode_y.data),
                         np.max(ds.NetNode_y.data),
                         num = 5,
                         endpoint = True)
            axs[cnt].set_yticks(yticks)
            axs[cnt].set_ylabel('%s (%s)'%(ds.NetNode_y.long_name, ds.NetNode_y.units))
    
        cnt = cnt + 1
    
    p0 = axs[0].get_position().get_points().flatten()
    p1 = axs[-1].get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], 0.1, p1[2]-p0[0], 0.015])
    plt.colorbar(pc, cax=ax_cbar, orientation='horizontal')
else:
    #set filename
    fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-'+expt+'\oceaneddymankmx0_map.nc'
    days2plt = np.linspace(start, stop, num = 4, endpoint = True)
    
    ugrid_all = get_netdata(file_nc=fname)
    ds = xr.open_dataset(fname)
    
    x = np.empty(ds.time.size)
    y = np.empty(ds.time.size)
    
    for i in np.arange(0, ds.time.size, 1):
        ssh = ds.mesh2d_s1.data[i,:]
        maxi = np.argmax(ssh)
        x[i], y[i] = ugrid_all.verts[maxi,:,:].mean(axis = 0)
    
    cnt = 0
    fig, axs = plt.subplots(1, 4, sharey=True, figsize=(15, 5))
    
    for i in days2plt.astype(int):
    
        #plot water level on map
        ssh = get_ncmodeldata(file_nc=fname, 
                              varname='mesh2d_s1', 
                              timestep=i)
        
        pc = plot_netmapdata(ugrid_all.verts, values=ssh[0,:], ax=axs[cnt], linewidth=0.5, cmap="jet")
        pc.set_clim([0, 0.025])
        axs[cnt].plot(x[:i],y[:i],'r.', markersize = 2)
        #axs[cnt].set_title(np.datetime_as_string(ds.time.data[i], unit = 'h'))
        axs[cnt].set_title('t = '+str(i)+' days')
        axs[cnt].set_aspect('equal')
        xticks = np.linspace(np.min(ds.mesh2d_node_x.data),
                         np.max(ds.mesh2d_node_x.data),
                         num = 5,
                         endpoint = True)
        axs[cnt].set_xticks(xticks)
        axs[cnt].set_xlabel('%s (%s)'%(ds.mesh2d_node_x.long_name, ds.mesh2d_node_x.units))
        if cnt == 0:
            yticks = np.linspace(np.min(ds.mesh2d_node_y.data),
                         np.max(ds.mesh2d_node_y.data),
                         num = 5,
                         endpoint = True)
            axs[cnt].set_yticks(yticks)
            axs[cnt].set_ylabel('%s (%s)'%(ds.mesh2d_node_y.long_name, ds.mesh2d_node_y.units))
    
        cnt = cnt + 1
    
    p0 = axs[0].get_position().get_points().flatten()
    p1 = axs[-1].get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], 0.1, p1[2]-p0[0], 0.015])
    plt.colorbar(pc, cax=ax_cbar, orientation='horizontal')