# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 10:40:12 2021

create the initial condition for the dflowfm ocean eddy validation case
based on make_vortex.m and barocvortex.m provided by P. Penven

dependencies:
    numpy
    progressbar
    croco_vgrid, zlevs (provided by P. Penven)
    tpx_tools (provided by P. Penven)
    interp_Cgrid (provided by P. Penven)
    xarray
    scipy.interpolate
    os

Ref:    Penven, P., L. Debreu, P. Marchesiello et J.C. McWilliams,   
        Application of the ROMS embedding procedure for the Central
        California Upwelling System,  Ocean Modelling, 2006.

@author: backeb
"""
#
# user defined stuff
#
expt = 'expt04'

# 
# import libraries
#
import numpy as np
from croco_vgrid import zlevs
import matplotlib.pyplot as plt
import xarray
from scipy.interpolate import RegularGridInterpolator#, griddata
import os

#
# Parameters for the grid
#
dx = 1000#30e3;                  # Horizontal resolution m-direction
xmax = 900e3;               # Half the domain length
H0 = 5000;                  # Depth
H = 2500;                   # Level of no-motion
theta = 38.5;               # Latitude (beta-plane)
R = 6367442.76;             # Earth radius
Pa = 1013e2;                # Atmospheric pressure
rho0 = 1024.4;              # Mean ocean density
umax = 1;                   # Max velocity (<0 cyclonic in the northern hemisphere)
radius=np.sqrt(2)*60e3;     # Vortex radius (np.exp(r2/radius2))
g=9.81;                     # Gravity acceleration
N2=(0.003)**2;              # Brunt-Vaissala frequency

#
# Vertical grid parameters
#
N=10;
theta_s=1;
theta_b=0;
hc=H0;
vtransform =  1.; # s-coordinate type (1: old- ; 2: new- coordinates)

#
################### END USERS DEFINED VARIABLES #######################
#

#
# Horzontal Grid
#
x = np.arange(-xmax-dx/2, xmax+dx, dx)
y=x;
dy=dx;
[X,Y]=np.meshgrid(x,y);

#
# Topo
#
h0=H0+0*X;

#
# Coriolis term (beta plane)
#
deg2rad=np.pi/180;
omega=2*np.pi/(24*3600);
f0=2*omega*np.sin(deg2rad*theta);
beta=2*omega/R*np.cos(deg2rad*theta);
f=f0+beta*Y;

#
# Compute zeta,ubar,vbar,u,v,t for the vortex
#

#
# P1 : pressure at z=0
# P1=P0*np.exp(-r2/radius2)
# P0 such as u(r0)=umax
# => P0=rho0 f0 umax radius sqrt(e/2)
#
P0=rho0*f0*umax*radius*np.sqrt(np.exp(1)/2);
P1=Pa+P0*np.exp(-(np.power(X,2)+np.power(Y,2))/radius**2);

#
# Calculate rho at z=0
#
a=-P0*(1-np.exp(-H))/(g*(H-1+np.exp(-H)));
rho1=rho0+a*np.exp(-(np.power(X,2)+np.power(Y,2))/radius**2);

#
# Surface elevation 
#
zeta=(P1-Pa)/(g*rho1);

#
# Vertical grid 
#
zw=zlevs(h0,zeta,theta_s,theta_b,hc,N,'w',vtransform);
zr=zlevs(h0,zeta,theta_s,theta_b,hc,N,'r',vtransform);

#
# Density
#
#
M, L = np.shape(X)
xr = np.reshape(X, (1, M, L))
xr = np.tile(xr, (N, 1, 1))
M, L = np.shape(Y)
yr = np.reshape(Y, (1, M, L))
yr = np.tile(yr, (N, 1, 1))

#
rho=rho0*(1-N2*zr/g);
rhodyn=-P0*(1-np.exp(-zr-H))*np.exp(-(np.power(xr,2)+np.power(yr,2))/radius**2)/(g*(H-1+np.exp(-H)));
## TODO ERROR RuntimeWarning: overflow encountered in exp
## see https://stackoverflow.com/questions/40726490/overflow-error-in-pythons-numpy-exp-function
## https://stackoverflow.com/questions/9559346/deal-with-overflow-in-exp-using-numpy
rho[zr>-H]=rho[zr>-H]+rhodyn[zr>-H];

#
# Temperature
#
# rho=rho0+R0-TCOEF*T
#
R0=30;
TCOEF=0.28;
t=(-rho+1000+R0)/TCOEF;

#
# U and V
#
#a=2*P0/(f0*rho0*radius^2);
#zu=0.5*(zr(:,:,1:end-1)+zr(:,:,2:end));
#xu=0.5*(xr(:,:,1:end-1)+xr(:,:,2:end));
#yu=0.5*(yr(:,:,1:end-1)+yr(:,:,2:end));
#F=(H-1+zu+np.exp(-zu-H))./(H-1+np.exp(-H));
#F(zu<-H)=0;
#u=a.*F.*yu.*np.exp(-(xu.^2+yu.^2)/radius^2);
#zv=0.5*(zr(:,1:end-1,:)+zr(:,2:end,:));
#xv=0.5*(xr(:,1:end-1,:)+xr(:,2:end,:));
#yv=0.5*(yr(:,1:end-1,:)+yr(:,2:end,:));
#F=(H-1+zv+np.exp(-zv-H))./(H-1+np.exp(-H));
#F(zv<-H)=0;
#v=-a.*F.*xv.*np.exp(-(xv.^2+yv.^2)/radius^2);

#
# Compute geostrophic velocities Vg
#
F=(H-1+zr+np.exp(-zr-H))/(H-1+np.exp(-H));
## TODO ERROR RuntimeWarning: overflow encountered in exp
## see https://stackoverflow.com/questions/40726490/overflow-error-in-pythons-numpy-exp-function
## https://stackoverflow.com/questions/9559346/deal-with-overflow-in-exp-using-numpy
F[zr<-H]=0;
r=np.sqrt(np.power(xr,2)+np.power(yr,2));
Vg=-(2*r*P0*F/(rho0*f0*radius*radius))*np.exp(-np.power(r,2)/radius**2);

#
# Compute gradient wind velocities Vgr : Vgr/R^2 + f Vgr = f Vg 
#
a=1+4*Vg/(f0*r);
indx=np.where(a<0);
if len(indx)>0:
    print(str(len(indx))+' points with no gradient wind solution (use geostrophy) ')
a[a<0]=1;
Vgr=2*Vg/(1+np.sqrt(a));
#Vgr=Vg;

#
# Project on the grid
#
ur=-Vgr*yr/r;
vr=Vgr*xr/r;
#u=0.5*(ur[:,:,:-1]+ur[:,:,1:]);
#v=0.5*(vr[:,:-1,:]+vr[:,1:,:]);
u=0.5*(ur[:,:,:]+ur[:,:,:]);
v=0.5*(vr[:,:,:]+vr[:,:,:]);

#
# Barotropic speeds
#
dz=zw[1:,:,:]-zw[:-1,:,:];
#dzu=0.5*(dz[:,:,:-1]+dz[:,:,1:]);
#dzv=0.5*(dz[:,:-1,:]+dz[:,1:,:]);
dzu=0.5*(dz[:,:,:]+dz[:,:,:]);
dzv=0.5*(dz[:,:,:]+dz[:,:,:]);
hu=np.squeeze(np.sum(dzu*u,axis=0));
hv=np.squeeze(np.sum(dzv*v,axis=0));
D_u=np.squeeze(np.sum(dzu,axis=0));
D_v=np.squeeze(np.sum(dzv,axis=0));
ubar=np.squeeze(hu/D_u);
vbar=np.squeeze(hv/D_v);

#
# create function and interpolate to dflowfm meshgrid
#
fpath = 'C:\\Users\\backeber\\OneDrive - Stichting Deltares\\Desktop\\Project-D-HYDRO-Phase-4\\dflowfm\\dflowfm_serial\\dfm_old_nctmplt\\'
ds = xarray.open_dataset(fpath+'oceaneddy_expt00_map.nc')
lon_minmax = [np.min(ds.FlowElem_xcc.data), np.max(ds.FlowElem_xcc.data)]
lat_minmax = [np.min(ds.FlowElem_ycc.data), np.max(ds.FlowElem_ycc.data)]
dlon = (lon_minmax[1]-lon_minmax[0])/(len(x))
dlat = (lat_minmax[1]-lat_minmax[0])/(len(x))
lon = np.linspace(lon_minmax[0]-dlon/2, lon_minmax[1]+dlon/2, len(x))
lat = np.linspace(lat_minmax[0]-dlat/2, lat_minmax[1]+dlat/2, len(y))
dlon = np.mean(np.diff(lon))
dlat = np.mean(np.diff(lat))


#
# interpolate zeta to dflowfm grid
#
f_amp = RegularGridInterpolator((lat,lon), zeta, bounds_error=False, fill_value=np.nan) 
#f_amp = RegularGridInterpolator((lon,lat), zeta, bounds_error=False, fill_value=np.nan) 
pli_coords = np.moveaxis(np.stack([ds.FlowElem_ycc.data,ds.FlowElem_xcc.data]),0,-1) # comparable to transpose...
s = f_amp(pli_coords)[np.newaxis]
s = np.append(s, s, axis=0)
#plt.scatter(ds.FlowElem_xcc,ds.FlowElem_ycc,1,s,cmap='jet')

#
# interpolate barotropic velocities to dflowfm grid
#
f_amp = RegularGridInterpolator((lat,lon), ubar, bounds_error=False, fill_value=np.nan)
pli_coords = np.moveaxis(np.stack([ds.FlowElem_ycc.data,ds.FlowElem_xcc.data]),0,-1) # comparable to transpose...
ucx = f_amp(pli_coords)[np.newaxis]
ucx = np.append(ucx, ucx, axis=0)
f_amp = RegularGridInterpolator((lat,lon), vbar)
pli_coords = np.moveaxis(np.stack([ds.FlowElem_ycc.data,ds.FlowElem_xcc.data]),0,-1) # comparable to transpose...
ucy = f_amp(pli_coords)[np.newaxis]
ucy = np.append(ucy, ucy, axis=0)

#
# write to dflowfm netcdf map file
#
wds = xarray.Dataset()
wds['mesh2d_enc_x'] = ds.mesh2d_enc_x
wds['mesh2d_enc_y'] = ds.mesh2d_enc_y
wds['mesh2d_enc_node_count'] = ds.mesh2d_enc_node_count
wds['mesh2d_enc_part_node_count'] = ds.mesh2d_enc_part_node_count
wds['mesh2d_enc_interior_ring'] = ds.mesh2d_enc_interior_ring
wds['mesh2d_enclosure_container'] = ds.mesh2d_enclosure_container
wds['Mesh2D'] = ds.Mesh2D
wds['wgs84'] = ds.wgs84
wds['NetNode_z'] = ds.NetNode_z
wds['NetLink'] = ds.NetLink
wds['NetLinkType'] = ds.NetLinkType
wds['NetElemNode'] = ds.NetElemNode
wds['NetElemLink'] = ds.NetElemLink
wds['NetLinkContour_x'] = ds.NetLinkContour_x
wds['NetLinkContour_y'] = ds.NetLinkContour_y
wds['NetLink_xu'] = ds.NetLink_xu
wds['NetLink_yu'] = ds.NetLink_yu
wds['BndLink'] = ds.BndLink
wds['FlowElem_zcc'] = ds.FlowElem_zcc
wds['FlowElem_bac'] = ds.FlowElem_bac
wds['FlowElem_xzw'] = ds.FlowElem_xzw
wds['FlowElem_yzw'] = ds.FlowElem_yzw
wds['FlowElemContour_x'] = ds.FlowElemContour_x
wds['FlowElemContour_y'] = ds.FlowElemContour_y
wds['FlowElem_bl'] = ds.FlowElem_bl
wds['ElemLink'] = ds.ElemLink
wds['FlowLink'] = ds.FlowLink
wds['FlowLinkType'] = ds.FlowLinkType
wds['timestep'] = ds.timestep
wds['s1'] = xarray.DataArray(data=s, name=ds.s1.name, dims=ds.s1.dims, attrs=ds.s1.attrs) #wds['s1'] = ds.s1
wds['s0'] = xarray.DataArray(data=s, name=ds.s0.name, dims=ds.s0.dims, attrs=ds.s0.attrs) #wds['s1'] = ds.s1
wds['waterdepth'] = ds.waterdepth
wds['numlimdt'] = ds.numlimdt
wds['taus'] = ds.taus
wds['unorm'] = ds.unorm
wds['u0'] = ds.u0
wds['q1'] = ds.q1
wds['viu'] = ds.viu
wds['diu'] = ds.diu
wds['ucx'] = xarray.DataArray(data=ucx, name=ds.ucx.name, dims=ds.ucx.dims, attrs=ds.ucx.attrs)
wds['ucy'] = xarray.DataArray(data=ucy, name=ds.ucy.name, dims=ds.ucy.dims, attrs=ds.ucy.attrs)
wds['czs'] = ds.czs
wds['czu'] = ds.czu
wds.attrs=ds.attrs

fname = "oceaneddy_init_map.nc"
try:
    os.remove(fpath+'..\\init_'+expt+'\\'+fname)
except OSError:
    pass

wds.to_netcdf(fpath+'..\\init_'+expt+'\\'+fname, 'w', 'NETCDF4')



##%%
##
## plot vortex
##
#fig, ax = plt.subplots(figsize=(15, 10))
#pc = plt.contourf(lon, lat, zeta, vmin = 0, vmax = 1, cmap="jet")
##c = plt.contourf(lon, lat, ubar, vmin = -0.4, vmax = 0.4, cmap="jet")
#ax.set_title('%s (%s)'%("Initial condition: zeta", "m"))
#ax.set_aspect('equal')
#fig.colorbar(pc, ax=ax)
#xticks = np.linspace(np.min(lon),
#                     np.max(lon),
#                     num = 5,
#                     endpoint = True)
#yticks = np.linspace(np.min(lat),
#                     np.max(lat),
#                     num = 5,
#                     endpoint = True)
#ax.set_xticks(xticks)
#ax.set_xlabel('%s (%s)'%("Longitude", "degrees E"))
#ax.set_yticks(yticks)
#ax.set_ylabel('%s (%s)'%("Latitude", "degrees N"))
#ax = plt.quiver(lon, lat, u[-1, :, :], v[-1, :, :], color = "w", scale = 50)
#

