# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 13:16:18 2021

@author: backeber
"""
#%matplotlib auto
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

cd C:\Users\backeber\OneDrive - Stichting Deltares\Desktop\Project-D-HYDRO-Phase-4\dflowfm\dflowfm_serial\map_init_old_netcdf_dx30e3 
ds = xr.open_dataset("oceaneddy_init_map.nc")
plt.scatter(ds.FlowElem_xcc,ds.FlowElem_ycc,1,ds.s1[0,:])
plt.scatter(ds.FlowElem_xcc,ds.FlowElem_ycc,1,ds.ucx[0,:])
plt.scatter(ds.FlowElem_xcc,ds.FlowElem_ycc,1,ds.ucy[0,:])
plt.scatter(ds.FlowLink_xu,ds.FlowLink_yu,1,ds.unorm[0,:])

cd C:\Users\backeber\OneDrive - Stichting Deltares\Desktop\Project-D-HYDRO-Phase-4\dflowfm\dflowfm_serial\map_init_old_netcdf 
ds2 = xr.open_dataset("oceaneddy_init_map.nc")
plt.scatter(ds2.FlowElem_xcc,ds2.FlowElem_ycc,1,ds2.s1[0,:])
plt.scatter(ds2.FlowElem_xcc,ds2.FlowElem_ycc,1,ds2.ucx[0,:])
plt.scatter(ds2.FlowElem_xcc,ds2.FlowElem_ycc,1,ds2.ucy[0,:])
plt.scatter(ds2.FlowLink_xu,ds2.FlowLink_yu,1,ds2.unorm[0,:])

cd C:\Users\backeber\OneDrive - Stichting Deltares\Desktop\Project-D-HYDRO-Phase-4\dflowfm\dflowfm_serial\init_expt01
ds = xr.open_dataset("oceaneddy_init_rst.nc")
plt.scatter(ds.FlowElem_xzw,ds.FlowElem_yzw,1,ds.s1[0,:])
plt.scatter(ds.FlowElem_xzw,ds.FlowElem_yzw,1,ds.ucx[0,:])
plt.scatter(ds.FlowElem_xzw,ds.FlowElem_yzw,1,ds.ucy[0,:])
plt.scatter(ds.FlowLink_xu,ds.FlowLink_yu,1,ds.unorm[0,:])

cd C:\Users\backeber\OneDrive - Stichting Deltares\Desktop\Project-D-HYDRO-Phase-4\dflowfm\dflowfm_serial\init_expt02
ds = xr.open_dataset("oceaneddy_init_rst.nc")
plt.scatter(ds.FlowElem_xzw,ds.FlowElem_yzw,1,ds.s1[0,:])
plt.scatter(ds.FlowElem_xzw,ds.FlowElem_yzw,1,ds.ucx[0,:])
plt.scatter(ds.FlowElem_xzw,ds.FlowElem_yzw,1,ds.ucy[0,:])
plt.scatter(ds.FlowLink_xu,ds.FlowLink_yu,1,ds.unorm[0,:])

cd C:\Users\backeber\OneDrive - Stichting Deltares\Desktop\Project-D-HYDRO-Phase-4\dflowfm\dflowfm_serial\init_expt03
ds = xr.open_dataset("oceaneddy_init_rst.nc")
plt.scatter(ds.FlowElem_xzw,ds.FlowElem_yzw,1,ds.s1[0,:])
plt.scatter(ds.FlowElem_xzw,ds.FlowElem_yzw,1,ds.ucx[0,:])
plt.scatter(ds.FlowElem_xzw,ds.FlowElem_yzw,1,ds.ucy[0,:])
plt.scatter(ds.FlowLink_xu,ds.FlowLink_yu,1,ds.unorm[0,:])

cd C:\Users\backeber\OneDrive - Stichting Deltares\Desktop\Project-D-HYDRO-Phase-4\dflowfm\dflowfm_serial\init_expt04
ds = xr.open_dataset("oceaneddy_init_map.nc")
plt.scatter(ds.FlowElem_xcc,ds.FlowElem_ycc,1,ds.s1[0,:])
plt.scatter(ds.FlowElem_xcc,ds.FlowElem_ycc,1,ds.ucx[0,:])
plt.scatter(ds.FlowElem_xcc,ds.FlowElem_ycc,1,ds.ucy[0,:])
plt.scatter(ds.FlowLink_xu,ds.FlowLink_yu,1,ds.unorm[0,:])
