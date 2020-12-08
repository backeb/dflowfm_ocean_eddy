# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 10:59:39 2020

@purpose:
    plot position of max ssh on model grid

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
fname1 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt1\oceaneddymankmx0_map.nc'
fname4 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt4\oceaneddymankmx0_map.nc'

ds1 = xr.open_dataset(fname1)
ugrid_all1 = get_netdata(file_nc=fname1)
ds4 = xr.open_dataset(fname4)
ugrid_all4 = get_netdata(file_nc=fname4)

x1 = np.empty(ds1.time.size)
y1 = np.empty(ds1.time.size)
x4 = np.empty(ds4.time.size)
y4 = np.empty(ds4.time.size)


for i in np.arange(0, ds1.time.size, 1):
    ssh1 = ds1.mesh2d_s1.data[i,:]
    ssh4 = ds4.mesh2d_s1.data[i,:]
    maxi1 = np.argmax(ssh1)
    maxi4 = np.argmax(ssh4)
    x1[i], y1[i] = ugrid_all1.verts[maxi1,:,:].mean(axis = 0)
    x4[i], y4[i] = ugrid_all4.verts[maxi4,:,:].mean(axis = 0)

fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all1.verts, 
                     values=None, 
                     ax=None, 
                     linewidth=0.1, 
                     color="gray", 
                     facecolor="None")
ax.set_aspect('equal')
ax.plot(x1, y1, '-', label = 'Expt 1')
ax.plot(x4, y4, '-', label = 'Expt 4')
ax.set_title('Maximum sea surface height location')
xticks = np.linspace(np.min(ds1.mesh2d_node_x.data),
                     np.max(ds1.mesh2d_node_x.data),
                     num = 10,
                     endpoint = True)
xticklabels = (xticks - np.median(xticks))/1000
yticks = np.linspace(np.min(ds1.mesh2d_node_y.data),
                     np.max(ds1.mesh2d_node_y.data),
                     num = 10,
                     endpoint = True)
yticklabels = (yticks - np.median(yticks))/1000
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels.astype(int))
ax.set_xlabel('Distance (km)')
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels.astype(int))
ax.set_ylabel('Distance (km)')
ax.legend()