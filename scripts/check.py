# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 13:16:18 2021

@author: backeber
"""
%matplotlib auto
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

cd f'C:\Users\backeber\OneDrive - Stichting Deltares\Desktop\Project-D-HYDRO-Phase-4\dflowfm\dflowfm_serial\restart_template_for_make_dflowfm_vortex'

rst = xr.open_dataset('oceaneddy_init_rst.nc')
plt.figure(),plt.scatter(rst.FlowElem_xzw, rst.FlowElem_yzw, 5, rst.s1)
plt.figure(),plt.scatter(rst.FlowElem_xzw, rst.FlowElem_yzw, 5, rst.ucx)
plt.figure(),plt.scatter(rst.FlowElem_xzw, rst.FlowElem_yzw, 5, rst.ucy)
plt.figure(),plt.scatter(rst.FlowLink_xu, rst.FlowLink_yu, 5, rst.unorm)
plt.figure(),plt.scatter(rst.FlowLink_xu, rst.FlowLink_yu, 5, rst.u0)

mp = xr.open_dataset('../DFM_OUTPUT_oceaneddy_expt00/oceaneddy_expt00_map.nc')
plt.figure(),plt.scatter(mp.mesh2d_face_x, mp.mesh2d_face_y, 5, mp.mesh2d_s1[0,:])
plt.figure(),plt.scatter(mp.mesh2d_face_x, mp.mesh2d_face_y, 5, mp.mesh2d_ucx[0,:])
plt.figure(),plt.scatter(mp.mesh2d_face_x, mp.mesh2d_face_y, 5, mp.mesh2d_ucy[0,:])
plt.figure(),plt.scatter(mp.mesh2d_edge_x, mp.mesh2d_edge_y, 5, mp.mesh2d_u1[0,:])

