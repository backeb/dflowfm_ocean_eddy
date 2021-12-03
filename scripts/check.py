# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 13:16:18 2021

@author: backeber
"""
#%matplotlib auto
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

fpath = 'C:\\Users\\backeber\OneDrive - Stichting Deltares\\Desktop\\Project-D-HYDRO-Phase-4\\dflowfm\\dflowfm_serial\\restart_template_for_make_dflowfm_vortex_v3\\'

rst = xr.open_dataset(fpath+'oceaneddy_init_rst.nc')
plt.figure(),plt.scatter(rst.FlowElem_xzw, rst.FlowElem_yzw, 5, rst.s1)
plt.figure(),plt.scatter(rst.FlowElem_xzw, rst.FlowElem_yzw, 5, rst.ucx)
plt.figure(),plt.scatter(rst.FlowElem_xzw, rst.FlowElem_yzw, 5, rst.ucy)
plt.figure(),plt.scatter(rst.FlowLink_xu, rst.FlowLink_yu, 5, rst.unorm)
plt.figure(),plt.scatter(rst.FlowLink_xu, rst.FlowLink_yu, 5, rst.u0)

#mp = xr.open_dataset(fpath+'../DFM_OUTPUT_oceaneddy_expt00/oceaneddy_expt00_map.nc')
#mp = xr.open_dataset(fpath+'oceaneddy_expt00_map.nc')
plt.figure(),plt.scatter(mp.mesh2d_face_x, mp.mesh2d_face_y, 5, mp.mesh2d_s1[0,:])
plt.figure(),plt.scatter(mp.mesh2d_face_x, mp.mesh2d_face_y, 5, mp.mesh2d_ucx[0,:])
plt.figure(),plt.scatter(mp.mesh2d_face_x, mp.mesh2d_face_y, 5, mp.mesh2d_ucy[0,:])
plt.figure(),plt.scatter(mp.mesh2d_edge_x, mp.mesh2d_edge_y, 5, mp.mesh2d_u1[0,:])

ds1 = xr.open_dataset(f"C:/Users/backeber/OneDrive - Stichting Deltares/Desktop/Project-D-HYDRO-Phase-4/dflowfm/dflowfm_serial/restart_template_for_make_dflowfm_vortex/oceaneddy_expt00_map.nc")
ds2 = xr.open_dataset(f"C:/Users/backeber/OneDrive - Stichting Deltares/Desktop/Project-D-HYDRO-Phase-4/dflowfm/dflowfm_serial/restart_template_for_make_dflowfm_vortex_v3/oceaneddy_expt00_map.nc")

mp = xr.open_dataset(f'C:/Users/backeber/OneDrive - Stichting Deltares/Desktop/Project-D-HYDRO-Phase-4/dflowfm/dflowfm_serial/restart_template_for_make_dflowfm_vortex_v3/oceaneddy_init_map.nc')
plt.figure(),plt.scatter(mp.FlowElem_xcc, mp.FlowElem_ycc, 5, mp.s1[0,:])
plt.figure(),plt.scatter(mp.FlowElem_xcc, mp.FlowElem_ycc, 5, mp.ucx[0,:])
plt.figure(),plt.scatter(mp.FlowElem_xcc, mp.FlowElem_ycc, 5, mp.ucy[0,:])
plt.figure(),plt.scatter(mp.FlowLink_xu, mp.FlowLink_yu, 5, mp.unorm[0,:])
