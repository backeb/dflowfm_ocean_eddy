# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 08:51:32 2020

@purpose:
    extract max ssh at each time step and plot as time series

@author: backeb
"""

"import libraries"
# from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
# from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

#set filename
fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt2\oceaneddymankmx0_map.nc'

ds = xr.open_dataset(fname)
ssh = ds.mesh2d_s1.data
ssh_max = np.empty(ssh.shape[0])

for i in np.arange(0,len(ssh_max),1):
    ssh_max[i] = np.nanmax(ssh[i,:])
    
fig, ax = plt.subplots(figsize=(15, 5))
ax.plot(ssh_max,'.',label = 'Expt 1')
ax.set_title('Sea surface height maximum over time')
ax.set_xlabel('timestep')
ax.set_ylabel('Sea surface height (m)')
ax.legend()

