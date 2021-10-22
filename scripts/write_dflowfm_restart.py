# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:52:46 2020

read restart file and write ocean eddy initial conditions to restart file

@author: backeb
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np




ds = xr.open_dataset('oceaneddy_expt00_20010101_000000_rst.nc')

ssh = ds.s1.data

gridx = ds.FlowElem_xzw.data
gridy = ds.FlowElem_yzw.data



from scipy.interpolate import RegularGridInterpolator#, griddata
# x and y are both 1D from a regularxgrid
# bedlev_corr_grid is 2D with the dimensions y,x (or the other way round, but then switch dimensions below
#create function with relation of lat/lon/watercolumn_grid
f_amp = RegularGridInterpolator((x,y),bedlev_corr_grid, bounds_error=False, fill_value=np.nan) 

pli_coords = np.stack([gridx,gridy])

elev_sel_ext = f_amp(pli_coords) #interpolate the function to the points