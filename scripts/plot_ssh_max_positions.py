# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 10:59:39 2020

@purpose:
    plot position of max ssh on model grid

@author: 
    backeb
"""

"import libraries"
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

#set filename
fname11 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt11\oceaneddymankmx0_map.nc' # beta plane 100x100_dx2000
fname12 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt12\oceaneddymankmx0_map.nc' # sferic 100x100_dx2000
fname13 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt13\oceaneddymankmx0_map.nc' # sferic 100x100_dx2000 oceaneddysizefrac=0.025 
fname14 = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-expt14\oceaneddymankmx0_map.nc' # init from map, use release v2021.03, oceaneddysizefrac=0.05 

ds11 = xr.open_dataset(fname11)
ugrid_all11 = get_netdata(file_nc=fname11)
ds12 = xr.open_dataset(fname12)
ugrid_all12 = get_netdata(file_nc=fname12)
ds13 = xr.open_dataset(fname13)
ugrid_all13 = get_netdata(file_nc=fname13)
ds14 = xr.open_dataset(fname14)
ugrid_all14 = get_netdata(file_nc=fname14)


x11 = np.empty(ds11.time.size)
y11 = np.empty(ds11.time.size)
x12 = np.empty(ds12.time.size)
y12 = np.empty(ds12.time.size)
x13 = np.empty(ds13.time.size)
y13 = np.empty(ds13.time.size)
x14 = np.empty(ds14.time.size)
y14 = np.empty(ds14.time.size)

for i in np.arange(0, ds11.time.size, 1):
    ssh11 = ds11.mesh2d_s1.data[i,:]
    ssh12 = ds12.mesh2d_s1.data[i,:]
    ssh13 = ds13.mesh2d_s1.data[i,:]
    ssh14 = ds14.s1.data[i,:]
    maxi11 = np.argmax(ssh11)
    maxi12 = np.argmax(ssh12)
    maxi13 = np.argmax(ssh13)
    maxi14 = np.argmax(ssh14)
    x11[i], y11[i] = ugrid_all11.verts[maxi11,:,:].mean(axis = 0)
    x12[i], y12[i] = ugrid_all12.verts[maxi12,:,:].mean(axis = 0)
    x13[i], y13[i] = ugrid_all13.verts[maxi13,:,:].mean(axis = 0)
    x14[i], y14[i] = ugrid_all14.verts[maxi14,:,:].mean(axis = 0)
    
fig, axs = plt.subplots(1, 2, figsize=(15, 5))
pc = plot_netmapdata(ugrid_all11.verts, 
                     values=None, 
                     ax=axs[0],
                     linewidth=0.1, 
                     color="gray", 
                     facecolor="None")
axs[0].set_aspect('equal')
axs[0].plot(x11, y11, '-', label = 'Expt 11')
axs[0].set_title('Beta plane expt')
xticks = np.linspace(np.min(ds11.mesh2d_node_x.data),
                     np.max(ds11.mesh2d_node_x.data),
                     num = 5,
                     endpoint = True)
yticks = np.linspace(np.min(ds11.mesh2d_node_y.data),
                     np.max(ds11.mesh2d_node_y.data),
                     num = 5,
                     endpoint = True)
axs[0].set_xticks(xticks)
axs[0].set_xlabel('%s (%s)'%(ds11.mesh2d_node_x.long_name, ds11.mesh2d_node_x.units))
axs[0].set_yticks(yticks)
axs[0].set_ylabel('%s (%s)'%(ds11.mesh2d_node_y.long_name, ds11.mesh2d_node_y.units))
axs[0].legend()

pc = plot_netmapdata(ugrid_all12.verts, 
                     values=None, 
                     ax=axs[1], 
                     linewidth=0.1, 
                     color="gray", 
                     facecolor="None")
axs[1].set_aspect('equal')
axs[1].plot(x12, y12, '-', label = 'Expt 12')
axs[1].plot(x13, y13, '-', label = 'Expt 13')
axs[1].plot(x14, y14, '-', label = 'Expt 14')
axs[1].set_title('Spherical grid expts')
xticks = np.linspace(np.min(ds12.mesh2d_node_x.data),
                     np.max(ds12.mesh2d_node_x.data),
                     num = 5,
                     endpoint = True)
yticks = np.linspace(np.min(ds12.mesh2d_node_y.data),
                     np.max(ds12.mesh2d_node_y.data),
                     num = 5,
                     endpoint = True)
axs[1].set_xticks(xticks)
axs[1].set_xlabel('%s (%s)'%(ds12.mesh2d_node_x.long_name, ds12.mesh2d_node_x.units))
axs[1].set_yticks(yticks)
axs[1].set_ylabel('%s (%s)'%(ds12.mesh2d_node_y.long_name, ds12.mesh2d_node_y.units))
axs[1].legend()
