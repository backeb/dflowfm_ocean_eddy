# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 10:59:39 2020

@author: backeber
"""

"import libraries"
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

#set filename
fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt2\oceaneddymankmx0_map.nc'

ds = xr.open_dataset(fname)
ugrid_all = get_netdata(file_nc=fname)

x = np.empty(ds.time.size)
y = np.empty(ds.time.size)


for i in np.arange(0, ds.time.size, 1):
    ssh = ds.mesh2d_s1.data[i,:]
    maxi = np.argmax(ssh)
    x[i], y[i] = ugrid_all.verts[maxi,:,:].mean(axis = 0)

fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, 
                     values=None, 
                     ax=None, 
                     linewidth=0.1, 
                     color="gray", 
                     facecolor="None")
ax.set_aspect('equal')
ax.plot(x, y,'r-',label = 'Expt 1')
ax.set_title('Maximum sea surface height location')
xticks = np.linspace(np.min(ds.mesh2d_node_x.data),
                     np.max(ds.mesh2d_node_x.data),
                     num = 10,
                     endpoint = True)
xticklabels = (xticks - np.median(xticks))/1000
yticks = np.linspace(np.min(ds.mesh2d_node_y.data),
                     np.max(ds.mesh2d_node_y.data),
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