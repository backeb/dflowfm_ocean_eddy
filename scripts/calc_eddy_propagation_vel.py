# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 10:56:28 2020

@author: backeber
"""

"import libraries"
from dfm_tools.get_nc import get_netdata
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from geopy.distance import geodesic

"define some variables for plotting"
expt  = 'expt13'

#set filename
fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0-'+expt+'\oceaneddymankmx0_map.nc'

ugrid_all = get_netdata(file_nc=fname)
ds = xr.open_dataset(fname)

x = np.empty(ds.time.size)
y = np.empty(ds.time.size)

for i in np.arange(0, ds.time.size, 1):
    ssh = ds.mesh2d_s1.data[i,:]
    maxi = np.argmax(ssh)
    x[i], y[i] = ugrid_all.verts[maxi,:,:].mean(axis = 0)

if expt == 'expt11': # not convinced about this 
    d = x**2 + y**2
    d = np.sqrt(d)/100/10000
else:
    d = np.empty(x.size-1)    
    for i in np.arange(0, x.size-1, 1):
        d[i] = geodesic((x[i], y[i]), (x[i+1], y[i+1])).km
