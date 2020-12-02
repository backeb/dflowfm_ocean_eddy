# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:37:21 2020
Updated on 

Purpose: plot eddy evolution on map

@author: backeb
"""

"import libraries"
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import matplotlib.pyplot as plt

#set filename
fname = 'c:\oceaneddy\DFM_OUTPUT_oceaneddymankmx0\oceaneddymankmx0_map.nc'

cnt = 0
fig, axs = plt.subplots(1, 4, sharey=True, figsize=(15, 10))

for i in [30, 60, 90, 120]:

#for i in tsteps2plot:

    ugrid_all = get_netdata(file_nc=fname)

    #plot water level on map
    ssh = get_ncmodeldata(file_nc=fname, 
                          varname='mesh2d_s1', 
                          timestep=i)

    #    fig, axs = plt.subplots(2, int(len(tsteps2plot)/2), sharey=True)
    pc = plot_netmapdata(ugrid_all.verts, values=ssh[0,:], ax=axs[cnt], linewidth=0.5, cmap="jet")
    pc.set_clim([0, 0.025])
    #ax.set_title('%s (%s)'%(ssh.var_varname, ssh.var_ncvarobject.units))
    axs[cnt].set_title('t = '+str(i)+' days')
    axs[cnt].set_aspect('equal')

    cnt = cnt + 1

p0 = axs[0].get_position().get_points().flatten()
p1 = axs[-1].get_position().get_points().flatten()
ax_cbar = fig.add_axes([p0[0], 0.3, p1[2]-p0[0], 0.015])
plt.colorbar(pc, cax=ax_cbar, orientation='horizontal')