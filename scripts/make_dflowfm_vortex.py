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
#%%
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
dx = 30e3;                  # Horizontal resolution m-direction
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
fpath = 'C:\\Users\\backeber\\OneDrive - Stichting Deltares\\Desktop\\Project-D-HYDRO-Phase-4\\dflowfm\\dflowfm_serial\\restart_template_for_make_dflowfm_vortex\\'
ds = xarray.open_dataset(fpath+'oceaneddy_expt00_20010101_000000_rst.nc')
lon_minmax = [np.min(ds.FlowElem_xzw.data), np.max(ds.FlowElem_xzw.data)]
lat_minmax = [np.min(ds.FlowElem_yzw.data), np.max(ds.FlowElem_yzw.data)]
dlon = (lon_minmax[1]-lon_minmax[0])/(len(x))
dlat = (lat_minmax[1]-lat_minmax[0])/(len(x))
lon = np.linspace(lon_minmax[0]-dlon/2, lon_minmax[1]+dlon/2, len(x))
lat = np.linspace(lat_minmax[0]-dlat/2, lat_minmax[1]+dlon/2, len(y))
dlon = np.mean(np.diff(lon))
dlat = np.mean(np.diff(lat))


#
# interpolate zeta to dflowfm grid
#
f_amp = RegularGridInterpolator((lon,lat), zeta)#, bounds_error=False, fill_value=np.nan) 
#f_amp = RegularGridInterpolator((lon,lat), zeta, bounds_error=False, fill_value=np.nan) 
pli_coords = np.moveaxis(np.stack([ds.FlowElem_xzw.data,ds.FlowElem_yzw.data]),0,-1) # comparable to transpose...
s = f_amp(pli_coords)[np.newaxis]
#plt.scatter(ds.FlowElem_xzw,ds.FlowElem_yzw,1,data,cmap='jet')

#
# interpolate barotropic velocities to dflowfm grid
#
f_amp = RegularGridInterpolator((lon,lat), ubar)
pli_coords = np.moveaxis(np.stack([ds.FlowElem_xzw.data,ds.FlowElem_yzw.data]),0,-1) # comparable to transpose...
ucx = f_amp(pli_coords)[np.newaxis]
f_amp = RegularGridInterpolator((lon,lat), vbar)
pli_coords = np.moveaxis(np.stack([ds.FlowElem_xzw.data,ds.FlowElem_yzw.data]),0,-1) # comparable to transpose...
ucy = f_amp(pli_coords)[np.newaxis]

#%%
# interpolate unorm to dflowfm grid
# TODO add check re direction of faces

unorm = np.array([])
gds = xarray.open_dataset(fpath+'oceaneddy_expt00_map.nc')

# using nearest neighbour search
from sklearn.neighbors import BallTree
# Setup Balltree using df as reference dataset
# Use Haversine calculate distance between points on the earth from lat/long
# haversine - https://pypi.org/project/haversine/ 
# these are the lon lat points where we HAVE ucx and ucy
tree = BallTree(np.deg2rad(np.swapaxes([ds.FlowElem_xzw.data, ds.FlowElem_yzw.data],0,1)), metric='haversine')

# first figure out if its a vertical or horizontal edge
# then interpolate to that location and append to unorm
for i in np.arange(0, np.shape(gds.mesh2d_edge_faces.data)[0], 1):
    if (np.diff(gds.mesh2d_face_x.data[gds.mesh2d_edge_faces.data[i,:].astype(int)-1]) == 0):
        # we want to find ds.FlowElem_xzw and ds.FlowElem_yzw that are closest to gds.mesh2d_edge_x.data[i] and gds.mesh2d_edge_y.data[i]
        # use k = 3 for 3 closest neighbors
        distances, indices = tree.query(np.deg2rad(np.c_[gds.mesh2d_edge_x.data[i], gds.mesh2d_edge_y.data[i]]), k = 3)
        unorm = np.append(unorm, ucx[0,indices[np.where(distances == distances.min())]])
        #unorm = np.append(unorm, np.interp(gds.mesh2d_edge_x.data[i], ds.FlowElem_xzw.data, np.squeeze(ucx)))
        #unorm = np.append(unorm, ucx[0,i])
    elif (np.diff(gds.mesh2d_face_x.data[gds.mesh2d_edge_faces.data[i,:].astype(int)-1]) != 0):
        distances, indices = tree.query(np.deg2rad(np.c_[gds.mesh2d_edge_x.data[i], gds.mesh2d_edge_y.data[i]]), k = 3)
        unorm = np.append(unorm, ucy[0,indices[np.where(distances == distances.min())]])
        #unorm = np.append(unorm, np.interp(gds.mesh2d_edge_x.data[i], ds.FlowElem_xzw.data, np.squeeze(ucy)))
        #unorm = np.append(unorm, ucy[0,i])

#plt.scatter(gds.mesh2d_edge_x, gds.mesh2d_edge_y, 1, unorm)    

#f_amp = RegularGridInterpolator((gds.mesh2d_edge_x.data,gds.mesh2d_edge_y.data), unorm)
#pli_coords = np.moveaxis(np.stack([ds.FlowLink_xu.data,ds.FlowLink_yu.data]),0,-1)
#unorm2 = f_amp(pli_coords)[np.newaxis]
#plt.scatter(ds.FlowLink_xu, ds.FlowLink_yu, 1, unorm2)

#%% some tests and checks
gds.mesh2d_node_x.data[gds.mesh2d_edge_faces.data[:100,:].astype(int)-1] 
# the above returns the nodes of two edge faces at i
np.diff(gds.mesh2d_node_x.data[gds.mesh2d_edge_faces.data[:100,:].astype(int)-1])
# the above returns the difference between the nodes of the two edge faces at i
# when =0 we have a vertical edge, therefore should specify ucx
# when !=0 we have a horizontal edge, 
# BUT there are negative and positive differences 
# AND there are differences > abs(0.1)
gds.mesh2d_face_x.data[gds.mesh2d_edge_faces.data[:100,:].astype(int)-1]
np.diff(gds.mesh2d_face_x.data[gds.mesh2d_edge_faces.data[:100,:].astype(int)-1])
# the above gets rid of differences > abs(0.1)
# AND Arthur said to use faces instead of nodes

#%% try nearest neighbour search
from sklearn.neighbors import BallTree
# Setup Balltree using df as reference dataset
# Use Haversine calculate distance between points on the earth from lat/long
# haversine - https://pypi.org/project/haversine/ 
# these are the lon lat points where we HAVE ucx and ucy
tree = BallTree(np.deg2rad(np.swapaxes([ds.FlowElem_xzw.data, ds.FlowElem_yzw.data],0,1)), metric='haversine')
# we want to find ds.FlowElem_xzw and ds.FlowElem_yzw 
# that are closest to gds.mesh2d_edge_x.data[i] and gds.mesh2d_edge_y.data[i]
# use k = 3 for 3 closest neighbors
distances, indices = tree.query(np.deg2rad(np.c_[gds.mesh2d_edge_x.data[i], gds.mesh2d_edge_y.data[i]]), k = 3)


#%%
# write to dflowfm netcdf restart file
#
wds = xarray.Dataset()
wds['timestep'] = ds.timestep
wds['s1'] = xarray.DataArray(data=s, name=ds.s1.name, dims=ds.s1.dims, attrs=ds.s1.attrs)
wds['s0'] = xarray.DataArray(data=s, name=ds.s0.name, dims=ds.s0.dims, attrs=ds.s0.attrs)
wds['taus'] = ds.taus # taucurrent in flow element center
wds['czs'] = ds.czs # Chezy roughness in flow element center
wds['FlowElem_bl'] = ds.FlowElem_bl # bed level at flow element circumcenter
#wds['ucx'] = ds.ucx # eastward_sea_water_velocity
wds['ucx'] = xarray.DataArray(data=ucx, name=ds.ucx.name, dims=ds.ucx.dims, attrs=ds.ucx.attrs)
#wds['ucy'] = ds.ucy # northward_sea_water_velocity
wds['ucy'] = xarray.DataArray(data=ucy, name=ds.ucy.name, dims=ds.ucy.dims, attrs=ds.ucy.attrs)
#wds['unorm'] = ds.unorm # normal component of sea_water_speed
wds['unorm'] = xarray.DataArray(data=unorm2, name=ds.unorm.name, dims=ds.unorm.dims, attrs=ds.unorm.attrs)
#wds['u0'] = ds.u0 # normal component of velocity through flow link at previous ti...
wds['u0'] = xarray.DataArray(data=unorm2, name=ds.u0.name, dims=ds.u0.dims, attrs=ds.u0.attrs)
wds['q1'] = ds.q1 # discharge
wds['qa'] = ds.qa # discharge used in advection
wds['squ'] = ds.squ # cell center outcoming flux
wds['FlowElem_xzw'] = ds.FlowElem_xzw # longitude
wds['FlowElem_yzw'] = ds.FlowElem_yzw # latitude
wds.attrs=ds.attrs

fname = "oceaneddy_init_rst.nc"
try:
    os.remove(fpath+fname)
except OSError:
    pass
print("writing datastack to netcdf4 :: "+fpath+fname)
wds.to_netcdf(fpath+fname, 'w', 'NETCDF4')


##%%
##
## plot vortex
##
#fig, ax = plt.subplots(figsize=(15, 10))
##pc = plt.contourf(X, Y, zeta, vmin = 0, vmax = 1, cmap="jet")
##pc = plt.pcolormesh(lon, lat, zeta, vmin = 0, vmax = 1, cmap="jet")
#pc = plt.pcolormesh(lon2, lat2, unorm, cmap="jet")
#ax.set_title('%s (%s)'%("Initial condition: sea level", "m"))
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
##ax = plt.quiver(lon, lat, ur[-1, :, :], vr[-1, :, :], color = "w", scale = 50)
#
##spd = np.sqrt(ur**2+vr**2)
##fig, ax = plt.subplots(figsize=(15, 10))
###pc = plt.contourf(X, Y, zeta, vmin = 0, vmax = 1, cmap="jet")
##pc = plt.pcolormesh(X, Y, spd[-1,:,:], vmin = 0, vmax = np.max(spd), cmap="jet")
##ax.set_title('%s (%s)'%("Initial condition: sea level", "m"))
##ax.set_aspect('equal')
##fig.colorbar(pc, ax=ax)
##xticks = np.linspace(np.min(X),
##                     np.max(X),
##                     num = 5,
##                     endpoint = True)
##yticks = np.linspace(np.min(Y),
##                     np.max(Y),
##                     num = 5,
##                     endpoint = True)
##ax.set_xticks(xticks)
##ax.set_xlabel('%s (%s)'%("Distance", "m"))
##ax.set_yticks(yticks)
##ax.set_ylabel('%s (%s)'%("Distance", "m"))
#
