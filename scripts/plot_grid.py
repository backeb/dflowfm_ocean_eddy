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

expt = 'expt00'

"set filename - note use *_map.nc to plot grid"
fname = 'C:\\Users\\backeber\\OneDrive - Stichting Deltares\\Desktop\\Project-D-HYDRO-Phase-4\\dflowfm\\dflowfm_serial\\DFM_OUTPUT_oceaneddy_'+expt+'\oceaneddy_'+expt+'_map.nc'

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

xticks = np.linspace(np.min(ds.NetNode_x.data),
                      np.max(ds.NetNode_x.data),
                      num = 5,
                      endpoint = True)
yticks = np.linspace(np.min(ds.NetNode_y.data),
                      np.max(ds.NetNode_y.data),
                      num = 5,
                      endpoint = True)
ax.set_xticks(xticks)
ax.set_xlabel('%s (%s)'%(ds.NetNode_x.long_name, ds.NetNode_x.units))
ax.set_yticks(yticks)
ax.set_ylabel('%s (%s)'%(ds.NetNode_y.long_name, ds.NetNode_y.units))
