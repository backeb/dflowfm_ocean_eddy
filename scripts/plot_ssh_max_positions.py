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
fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0\oceaneddymankmx0_map.nc'

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
ax.set_xticks(np.arange(0,220000,20000))
ax.set_xticklabels(np.arange(-50,60,10))
ax.set_xlabel('Distance (km)')
ax.set_yticks(np.arange(0,220000,20000))
ax.set_yticklabels(np.arange(-50,60,10))
ax.set_ylabel('Distance (km)')
ax.legend()