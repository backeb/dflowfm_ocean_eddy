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

#set filename
fname = 'C:\oceaneddy\DFM_OUTPUT_oceaneddy\oceaneddy_map.nc'

for i in [0, 30, 60, 90]:

    ugrid_all = get_netdata(file_nc=fname)

    ssh = get_ncmodeldata(file_nc=fname,
                        varname='mesh2d_s1', 
                        timestep=i)

    pc = plot_netmapdata(ugrid_all.verts, values=ssh[0,:], ax=axs[cnt], linewidth=0.5, cmap="jet")
    pc.set_clim([0, 0.025])
    axs[cnt].set_title('t = '+str(i)+' days')
    axs[cnt].set_aspect('equal')
    cnt = cnt + 1

