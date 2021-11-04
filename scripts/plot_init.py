# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 13:42:43 2020

@purpose:
    plot initial condition

@author: backeb
"""

# import libraries
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# set filename
fname = 'C:\\Users\\backeber\\OneDrive - Stichting Deltares\\Desktop\\Project-D-HYDRO-Phase-4\\dflowfm\\dflowfm_serial\\oceaneddy_expt00_20010101_000000_rst.nc'
ds = xr.open_dataset(fname)
# plt.scatter(ds.FlowElem_xzw.data,ds.FlowElem_yzw.data,1,ds.s1.data[0,:])

ugrid_all = get_netdata(file_nc='C:\\Users\\backeber\\OneDrive - Stichting Deltares\\Desktop\\Project-D-HYDRO-Phase-4\\dflowfm\\dflowfm_serial\\DFM_OUTPUT_oceaneddy_expt00\\oceaneddy_expt00_map.nc')
#ugrid_all = get_netdata(file_nc=fname)
ds = xr.open_dataset(fname)

ssh = get_ncmodeldata(file_nc=fname, varname='s1', timestep=0)

# find location of max ssh to add to plot
maxi = np.argmax(ssh)

fig, ax = plt.subplots()#(figsize=(15, 10))

pc = plot_netmapdata(ugrid_all.verts, values=ssh[0,:], ax=ax, linewidth=0.5, cmap="jet")
#pc.set_clim([0, 0.05])
x, y = ugrid_all.verts[maxi,:,:].mean(axis = 0)
ax.plot(x, y,'wx')
ax.set_title('%s (%s)'%(ds.s1.long_name, ds.s1.units))
ax.set_aspect('equal')
fig.colorbar(pc, ax=ax)
xticks = np.linspace(np.min(ds.FlowElem_xzw.data),
                      np.max(ds.FlowElem_xzw.data),
                      num = 5,
                      endpoint = True)
yticks = np.linspace(np.min(ds.FlowElem_yzw.data),
                      np.max(ds.FlowElem_yzw.data),
                      num = 5,
                      endpoint = True)

ax.set_xticks(xticks)
ax.set_xlabel('%s (%s)'%(ds.FlowElem_xzw.long_name, ds.FlowElem_xzw.units))
ax.set_yticks(yticks)
ax.set_ylabel('%s (%s)'%(ds.FlowElem_yzw.long_name, ds.FlowElem_yzw.units))