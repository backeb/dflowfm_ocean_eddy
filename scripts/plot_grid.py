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

cnt = 0
fig, axs = plt.subplots(1, 4, sharey=True, figsize=(15, 10))

for i in [0, 30, 60, 90]:

    ugrid_all = get_netdata(file_nc=fname)

    #plot water level on map
    ssh = get_ncmodeldata(file_nc=fname,
                        varname='mesh2d_s1', 
                        timestep=i)

    pc = plot_netmapdata(ugrid_all.verts, values=ssh[0,:], ax=axs[cnt], linewidth=0.5, cmap="jet")
    pc.set_clim([0, 0.025])
    axs[cnt].set_title('t = '+str(i)+' days')
    axs[cnt].set_aspect('equal')
    cnt = cnt + 1

