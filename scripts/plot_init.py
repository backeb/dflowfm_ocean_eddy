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
#from scipy.interpolate import griddata

expt = 'expt12'

# set filename
fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-'+expt+'\oceaneddymankmx0_map.nc'

ugrid_all = get_netdata(file_nc=fname)
ds = xr.open_dataset(fname)

ssh = get_ncmodeldata(file_nc=fname, varname='mesh2d_s1', timestep=0)

# find location of max ssh to add to plot
maxi = np.argmax(ssh)

fig, ax = plt.subplots()#(figsize=(15, 10))

pc = plot_netmapdata(ugrid_all.verts, values=ssh[0,:], ax=ax, linewidth=0.5, cmap="jet")
pc.set_clim([0, 0.05])
x, y = ugrid_all.verts[maxi,:,:].mean(axis = 0)
ax.plot(x, y,'wx')
ax.set_title('%s (%s)'%(ds.mesh2d_s1.long_name, ds.mesh2d_s1.units))
ax.set_aspect('equal')
fig.colorbar(pc, ax=ax)
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


'trying with xarray and regridding / reshaping'
# ds = xr.open_dataset(fname)
# ssh = ds.mesh2d_s1[0,:].data

# reg_x_vec = np.linspace(0,1,100)
# reg_y_vec = np.linspace(0,1,100)
# x_grid,y_grid = np.meshgrid(reg_x_vec,reg_y_vec)

# value_grid = griddata((np.linspace(0,1,100),np.linspace(0,1,100)),
#                       ssh,
#                       (x_grid,y_grid),
#                       method='nearest')
