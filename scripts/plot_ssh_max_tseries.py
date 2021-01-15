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
fname11 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt11\oceaneddymankmx0_map.nc' # beta plane 100x100_dx2000
fname12 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt12\oceaneddymankmx0_map.nc' # sferic 100x100_dx2000
fname13 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt13\oceaneddymankmx0_map.nc' # sferic 100x100_dx2000 oceaneddysizefrac=0.025 
fname14 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt14\oceaneddymankmx0_map.nc' # init from map, use release v2021.03, oceaneddysizefrac=0.05 

ds11 = xr.open_dataset(fname11)
ssh11 = ds11.mesh2d_s1.data
ssh_max11 = np.empty(ssh11.shape[0])

ds12 = xr.open_dataset(fname12)
ssh12 = ds12.mesh2d_s1.data
ssh_max12 = np.empty(ssh12.shape[0])

ds13 = xr.open_dataset(fname13)
ssh13 = ds13.mesh2d_s1.data
ssh_max13 = np.empty(ssh13.shape[0])

ds14 = xr.open_dataset(fname14)
ssh14 = ds14.s1.data
ssh_max14 = np.empty(ssh14.shape[0])


for i in np.arange(0,len(ssh_max11),1):
    ssh_max11[i] = np.nanmax(ssh11[i,:])
    ssh_max12[i] = np.nanmax(ssh12[i,:])
    ssh_max13[i] = np.nanmax(ssh13[i,:])
    ssh_max14[i] = np.nanmax(ssh14[i,:])
    
fig, ax = plt.subplots(figsize=(15, 5))
ax.plot(ssh_max11,'.',label = 'Expt 11')
ax.plot(ssh_max12,'.',label = 'Expt 12')
ax.plot(ssh_max13,'.',label = 'Expt 13')
ax.plot(ssh_max14,'.',label = 'Expt 14')
ax.set_ylim(bottom = 0)
ax.set_title('Sea surface height maximum over time')
ax.set_xlabel('days')
ax.set_ylabel('Sea surface height (m)')
ax.legend()

