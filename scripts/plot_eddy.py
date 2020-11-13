# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:37:21 2020
Updated on 

Purpose: plot eddy on map

@author: backeb
"""

"import libraries"
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import matplotlib.pyplot as plt

"set filename - note use *_map.nc to plot grid"
fname = 'C:\oceaneddy\DFM_OUTPUT_oceaneddy\oceaneddy_map.nc'

ugrid_all = get_netdata(file_nc=fname)

#plot water level on map
ssh = get_ncmodeldata(file_nc=fname, 
                      varname='mesh2d_s1', 
                      timestep=0)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, 
                     values=ssh[0,:], 
                     ax=None, 
                     linewidth=0.5, 
                     cmap="jet")
#pc.set_clim([-0.5,1])
fig.colorbar(pc, ax=ax)
#ax.set_title('%s (%s)'%(ssh.var_varname, ssh.var_ncvarobject.units))
ax.set_xlim(60000, 140000)
ax.set_ylim(60000, 140000)
ax.set_title('Baroclinic Vortex: t = 0 days')
ax.set_aspect('equal')

ssh = get_ncmodeldata(file_nc=fname, 
                      varname='mesh2d_s1', 
                      timestep=240)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, 
                     values=ssh[0,:], 
                     ax=None, 
                     linewidth=0.5, 
                     cmap="jet")
#pc.set_clim([-0.5,1])
fig.colorbar(pc, ax=ax)
#ax.set_title('%s (%s)'%(ssh.var_varname, ssh.var_ncvarobject.units))
ax.set_xlim(60000, 140000)
ax.set_ylim(60000, 140000)
ax.set_title('Baroclinic Vortex: t = 240 days')
ax.set_aspect('equal')
