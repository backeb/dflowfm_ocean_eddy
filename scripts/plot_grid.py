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

expt = 'expt10'

"set filename - note use *_map.nc to plot grid"
fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-'+expt+'\oceaneddymankmx0_map.nc'

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
                      num = 5,
                      endpoint = True)
yticks = np.linspace(np.min(ds.mesh2d_node_y.data),
                      np.max(ds.mesh2d_node_y.data),
                      num = 5,
                      endpoint = True)
ax.set_xticks(xticks)
ax.set_xlabel('%s (%s)'%(ds.mesh2d_node_x.long_name, ds.mesh2d_node_x.units))
ax.set_yticks(yticks)
ax.set_ylabel('%s (%s)'%(ds.mesh2d_node_y.long_name, ds.mesh2d_node_y.units))
