# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:37:21 2020
Updated on 

Purpose: plot eddy on map

@author: backeb
"""
#
# user defined variables
#
expt  = 'expt00'
var2plot = 'ucx'

#
# import libraries
#
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np


fname = 'C:\\Users\\backeber\\OneDrive - Stichting Deltares\\Desktop\\Project-D-HYDRO-Phase-4\\dflowfm\\dflowfm_serial\\DFM_OUTPUT_oceaneddy_'+expt+'\\oceaneddy_'+expt+'_map.nc'
ugrid_all = get_netdata(file_nc=fname)
ds = xr.open_dataset(fname)

for i in np.arange(0, 20, 1):#len(ds.time.data), 1):
    
    fig, axs = plt.subplots(1, 1, figsize=(5, 5))
    
    ssh = get_ncmodeldata(file_nc=fname, 
                          varname=var2plot, 
                          timestep=i)
    
    pc = plot_netmapdata(ugrid_all.verts, values=ssh[0,:], ax=axs, linewidth=0.5, cmap="jet")
    #pc.set_clim([0, 0.25])
    #axs.plot(x[:i],y[:i],'r.', markersize = 2)
    #axs.set_title(np.datetime_as_string(ds.time.data[i], unit = 'h'))
    axs.set_title('t = '+str(i)+' hours')
    axs.set_aspect('equal')
    xticks = np.linspace(np.min(ds.NetNode_x.data),
                      np.max(ds.NetNode_x.data),
                      num = 5,
                      endpoint = True)
    axs.set_xticks(xticks)
    axs.set_xlabel('%s (%s)'%(ds.NetNode_x.long_name, ds.NetNode_x.units))
    yticks = np.linspace(np.min(ds.NetNode_y.data),
                         np.max(ds.NetNode_y.data),
                         num = 5,
                         endpoint = True)
    axs.set_yticks(yticks)
    axs.set_ylabel('%s (%s)'%(ds.NetNode_y.long_name, ds.NetNode_y.units))
    
    # p0 = axs.get_position().get_points().flatten()
    # p1 = axs[-1].get_position().get_points().flatten()
    # ax_cbar = fig.add_axes([p0[0], 0.1, p1[2]-p0[0], 0.015])
    plt.colorbar(pc, orientation='horizontal')
    sname = 'C:\\Users\\backeber\\OneDrive - Stichting Deltares\\Desktop\\Project-D-HYDRO-Phase-4\\figures\\all_timesteps\\'+f'{i:04}'+'.png'
    plt.savefig(sname)
    plt.close("all")