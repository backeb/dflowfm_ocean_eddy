# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 13:42:43 2020

@purpose:
    plot initial condition

@author: backeb
"""

# import libraries
#from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
#from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# set filename
fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0\oceaneddymankmx0_map.nc'
ds = xr.open_dataset(fname)
ssh = ds.mesh2d_s1[0,:].data

reg_x_vec = np.linspace(0,1,100)
reg_y_vec = np.linspace(0,1,100)
x_grid,y_grid = np.meshgrid(reg_x_vec,reg_y_vec)

value_grid = griddata((np.linspace(0,1,100),np.linspace(0,1,100)),
                      ssh,
                      (x_grid,y_grid),
                      method='nearest')


# ugrid_all = get_netdata(file_nc=fname)
# ssh = get_ncmodeldata(file_nc=fname, varname='mesh2d_s1', timestep=0)

# fig, ax = plt.subplots()#(figsize=(15, 10))

# pc = plot_netmapdata(ugrid_all.verts, values=ssh[0,:], ax=ax, linewidth=0.5, cmap="jet")
# pc.set_clim([0, 0.05])
# ax.set_title('Initial condition (sea surface height in m)')
# ax.set_aspect('equal')
# fig.colorbar(pc, ax=ax)


