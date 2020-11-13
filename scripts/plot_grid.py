# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:52:46 2020
Updated on Tue Nov 10 15:16:00 2020

Plot model grid using dfm_tools (https://github.com/openearth/dfm_tools)

@author: backeb
"""

"import libraries"
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import matplotlib.pyplot as plt

"set filename - note use *_map.nc to plot grid"
fname = 'C:\oceaneddy\DFM_OUTPUT_oceaneddy\oceaneddy_map.nc'

ds = get_netdata(file_nc=fname)
fig, ax = plt.subplots()
pc = plot_netmapdata(ds.verts, 
                     values=None, 
                     ax=None, 
                     linewidth=0.5, 
                     color="crimson", 
                     facecolor="None")
ax.set_aspect('equal')

