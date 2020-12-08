# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:52:46 2020

Plot model grid using dfm_tools (https://github.com/openearth/dfm_tools)

@author: backeb

"""

"import libraries"
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

"set filename - note use *_map.nc to plot grid"
fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt1\oceaneddymankmx0_map.nc'

ugrid_all = get_netdata(file_nc=fname)
ds = xr.open_dataset(fname)

fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, 
                     values=None, 
                     ax=None, 
                     linewidth=0.5, 
                     color="crimson", 
                     facecolor="None")
ax.set_aspect('equal')

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
