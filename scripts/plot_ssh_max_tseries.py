# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 08:51:32 2020

@purpose:
    extract max ssh at each time step and plot as time series

@author: 
    backeb
"""

"import libraries"
# from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
# from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

#set filename
fname1 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt1\oceaneddymankmx0_map.nc'
fname4 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt4\oceaneddymankmx0_map.nc'

ds1 = xr.open_dataset(fname1)
ssh1 = ds1.mesh2d_s1.data
ssh_max1 = np.empty(ssh1.shape[0])

ds4 = xr.open_dataset(fname4)
ssh4 = ds4.mesh2d_s1.data
ssh_max4 = np.empty(ssh4.shape[0])

for i in np.arange(0,len(ssh_max1),1):
    ssh_max1[i] = np.nanmax(ssh1[i,:])
    ssh_max4[i] = np.nanmax(ssh4[i,:])
    
fig, ax = plt.subplots(figsize=(15, 5))
ax.plot(ssh_max1,'.',label = 'Expt 1')
ax.plot(ssh_max4,'.',label = 'Expt 4')
ax.set_title('Sea surface height maximum over time')
ax.set_xlabel('timestep')
ax.set_ylabel('Sea surface height (m)')
ax.legend()

