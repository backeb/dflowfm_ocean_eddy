# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 10:40:12 2021

create the initial condition for the dflowfm ocean eddy validation case
based on make_vortex.m and barocvortex.m provided by P. Penven

dependencies:
    numpy
    croco_vgrid, zlevs (provided by P. Penven)
    tpx_tools (provided by P. Penven)
    interp_Cgrid (provided by P. Penven)

Ref:    Penven, P., L. Debreu, P. Marchesiello et J.C. McWilliams,   
        Application of the ROMS embedding procedure for the Central
        California Upwelling System,  Ocean Modelling, 2006.

@author: backeb
"""
# 
# import libraries
#
import numpy as np
from croco_vgrid import zlevs
import matplotlib.pyplot as plt


#
# Parameters for the grid
#
dx = 30e3;                  # Horizontal resolution
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
u=0.5*(ur[:,:,:-1]+ur[:,:,1:]);
v=0.5*(vr[:,:-1,:]+vr[:,1:,:]);
#
# Barotropic speeds
#
dz=zw[1:,:,:]-zw[:-1,:,:];
dzu=0.5*(dz[:,:,:-1]+dz[:,:,1:]);
dzv=0.5*(dz[:,:-1,:]+dz[:,1:,:]);
hu=np.squeeze(np.sum(dzu*u));
hv=np.squeeze(np.sum(dzv*v));
D_u=np.squeeze(np.sum(dzu));
D_v=np.squeeze(np.sum(dzv));
ubar=np.squeeze(hu/D_u);
vbar=np.squeeze(hv/D_v);

#
# write to dflowfm netcdf restart file
#
#TODO


#
# plot vortex
#
fig, ax = plt.subplots(figsize=(15, 10))
#pc = plt.contourf(X, Y, zeta, vmin = 0, vmax = 1, cmap="jet")
pc = plt.pcolormesh(X, Y, zeta, vmin = 0, vmax = 1, cmap="jet")
ax.set_title('%s (%s)'%("Initial condition: sea level", "m"))
ax.set_aspect('equal')
fig.colorbar(pc, ax=ax)
xticks = np.linspace(np.min(X),
                     np.max(X),
                     num = 5,
                     endpoint = True)
yticks = np.linspace(np.min(Y),
                     np.max(Y),
                     num = 5,
                     endpoint = True)
ax.set_xticks(xticks)
ax.set_xlabel('%s (%s)'%("Distance", "m"))
ax.set_yticks(yticks)
ax.set_ylabel('%s (%s)'%("Distance", "m"))
ax = plt.quiver(X, Y, ur[-1, :, :], vr[-1, :, :], color = "w", scale = 50)

#spd = np.sqrt(ur**2+vr**2)
#fig, ax = plt.subplots(figsize=(15, 10))
##pc = plt.contourf(X, Y, zeta, vmin = 0, vmax = 1, cmap="jet")
#pc = plt.pcolormesh(X, Y, spd[-1,:,:], vmin = 0, vmax = np.max(spd), cmap="jet")
#ax.set_title('%s (%s)'%("Initial condition: sea level", "m"))
#ax.set_aspect('equal')
#fig.colorbar(pc, ax=ax)
#xticks = np.linspace(np.min(X),
#                     np.max(X),
#                     num = 5,
#                     endpoint = True)
#yticks = np.linspace(np.min(Y),
#                     np.max(Y),
#                     num = 5,
#                     endpoint = True)
#ax.set_xticks(xticks)
#ax.set_xlabel('%s (%s)'%("Distance", "m"))
#ax.set_yticks(yticks)
#ax.set_ylabel('%s (%s)'%("Distance", "m"))

