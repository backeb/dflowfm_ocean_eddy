# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:37:21 2020
Updated on 

@purpose: 
    plot eddy evolution on map

@author: 
    backeb
"""

"import libraries"
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

#set filename
fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt2\oceaneddymankmx0_map.nc'

ugrid_all = get_netdata(file_nc=fname)
ds = xr.open_dataset(fname)

ssh0 = get_ncmodeldata(file_nc=fname, varname='mesh2d_s1', timestep=0)
maxi = np.argmax(ssh0)
x0, y0 = ugrid_all.verts[maxi,:,:].mean(axis = 0)


cnt = 0
fig, axs = plt.subplots(1, 4, sharey=True, figsize=(15, 5))

for i in [30, 60, 90, 120]:

    #plot water level on map
    ssh = get_ncmodeldata(file_nc=fname, 
                          varname='mesh2d_s1', 
                          timestep=i)
    
    # find location of max ssh to add to plot
    maxi = np.argmax(ssh)
    
    #    fig, axs = plt.subplots(2, int(len(tsteps2plot)/2), sharey=True)
    pc = plot_netmapdata(ugrid_all.verts, values=ssh[0,:], ax=axs[cnt], linewidth=0.5, cmap="jet")
    pc.set_clim([0, 0.025])
    x, y = ugrid_all.verts[maxi,:,:].mean(axis = 0)
    axs[cnt].plot((x0,x),(y0, y),'wx-')
    #ax.set_title('%s (%s)'%(ssh.var_varname, ssh.var_ncvarobject.units))
    axs[cnt].set_title('t = '+str(i)+' days')
    axs[cnt].set_aspect('equal')
    xticks = np.linspace(np.min(ds.mesh2d_node_x.data),
                     np.max(ds.mesh2d_node_x.data),
                     num = 5,
                     endpoint = True)
    xticklabels = (xticks - np.median(xticks))/1000
    axs[cnt].set_xticks(xticks)
    axs[cnt].set_xticklabels(xticklabels.astype(int))
    axs[cnt].set_xlabel('Distance (km)')
    if cnt == 0:
        yticks = np.linspace(np.min(ds.mesh2d_node_y.data),
                     np.max(ds.mesh2d_node_y.data),
                     num = 10,
                     endpoint = True)
        yticklabels = (yticks - np.median(yticks))/1000
        axs[cnt].set_yticks(yticks)
        axs[cnt].set_yticklabels(yticklabels.astype(int))
        axs[cnt].set_ylabel('Distance (km)')

    cnt = cnt + 1

p0 = axs[0].get_position().get_points().flatten()
p1 = axs[-1].get_position().get_points().flatten()
ax_cbar = fig.add_axes([p0[0], 0.1, p1[2]-p0[0], 0.015])
plt.colorbar(pc, cax=ax_cbar, orientation='horizontal')