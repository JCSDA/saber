#!/usr/bin/env python3

import argparse
import os
import sys
from sys import exit
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import subprocess
import math

# Parser
parser = argparse.ArgumentParser()
parser.add_argument('build_dir', help='Build directory')
parser.add_argument('mode', help='Plot mode: hdiag_sc1 / hdiag_sc2 / hdiag_sc3 / hdiag_sc3-sc4 / fit_iso / fit_ani / fit_multi / rh_sampling / rh_map / nicas_iso / nicas_ani / nicas_mask')
args = parser.parse_args()

####################################################################################################
# READ NETCDF FILES ################################################################################
####################################################################################################

# Load HDIAG sampling (isotropic)
f_sampling = Dataset(args.build_dir + '/saber/test/testdata/doc_6/1-1_sampling_grids_local_000001-000001.nc', 'r', format='NETCDF4')
lon_iso = f_sampling['lon'][0,:,:,:]
lat_iso = f_sampling['lat'][0,:,:,:]
lon_local = f_sampling['lon_local'][0,:,:]
lat_local = f_sampling['lat_local'][0,:,:]
nc4_iso = lon_iso.shape[0]
nc3_iso = lon_iso.shape[1]
nc1 = lon_iso.shape[2]
nc2a = lon_local.shape[0]

# Load HDIAG sampling (anisotropic)
f_sampling = Dataset(args.build_dir + '/saber/test/testdata/doc_7/1-1_sampling_grids_local_000001-000001.nc', 'r', format='NETCDF4')
lon_ani = f_sampling['lon'][0,:,:,:]
lat_ani = f_sampling['lat'][0,:,:,:]
nc3_ani = lon_ani.shape[1]
nc4_ani = lon_ani.shape[0]

# Load HDIAG profile (isotropic)
f_profile = Dataset(args.build_dir + '/saber/test/testdata/doc_6/1-1_diag_plot.nc', 'r', format='NETCDF4')
disth_iso = f_profile['stream_function']['disth'][:]*6371.229
raw_hor_iso = f_profile['stream_function']['cor1']['raw_hor'][0,:,:]
fit_hor_iso = f_profile['stream_function']['cor1']['fit_hor'][0,:,:]

# Load HDIAG profile (anisotropic)
f_profile = Dataset(args.build_dir + '/saber/test/testdata/doc_7/1-1_diag_plot.nc', 'r', format='NETCDF4')
disth_ani = f_profile['stream_function']['disth'][:]*6371.229
angular_sector = f_profile['stream_function']['as'][:]
raw_hor_ani = f_profile['stream_function']['cor1']['raw_hor'][0,:,:]
fit_hor_ani = f_profile['stream_function']['cor1']['fit_hor'][0,:,:]
angular_sector_cyclic = []
raw_hor_ani_cyclic = np.empty([0,nc3_ani])
fit_hor_ani_cyclic = np.empty([0,nc3_ani])
for ic4 in range(0, nc4_ani):
    for ic3 in range(0, nc3_ani):
        if ma.is_masked(raw_hor_ani[ic4,ic3]):
            ic4m1 = ic4-1
            ic4p1 = ic4+1
            if ic4 == 0:
                ic4m1 = nc4_ani-1
            elif ic4 == nc4_ani-1:
                ic4p1 = 0
            if ma.is_masked(raw_hor_ani[ic4m1,ic3]) or ma.is_masked(raw_hor_ani[ic4p1,ic3]):
                print("Too many missing values")
            else:
                raw_hor_ani[ic4,ic3] = 0.5*(raw_hor_ani[ic4m1,ic3]+raw_hor_ani[ic4p1,ic3])
for i in range(0, nc4_ani):
    if angular_sector[i] >= 0:
        angular_sector_cyclic.append(angular_sector[i])
        raw_hor_ani_cyclic = np.append(raw_hor_ani_cyclic, [raw_hor_ani[i,:]], axis=0)
        fit_hor_ani_cyclic = np.append(fit_hor_ani_cyclic, [fit_hor_ani[i,:]], axis=0)
for i in range(0, nc4_ani):
      angular_sector_cyclic.append(angular_sector[i]+math.pi)
      raw_hor_ani_cyclic = np.append(raw_hor_ani_cyclic, [raw_hor_ani[i,:]], axis=0)
      fit_hor_ani_cyclic = np.append(fit_hor_ani_cyclic, [fit_hor_ani[i,:]], axis=0)
for i in range(0, nc4_ani):
    if angular_sector[i] <= 0:
        angular_sector_cyclic.append(angular_sector[i]+2*math.pi)
        raw_hor_ani_cyclic = np.append(raw_hor_ani_cyclic, [raw_hor_ani[i,:]], axis=0)
        fit_hor_ani_cyclic = np.append(fit_hor_ani_cyclic, [fit_hor_ani[i,:]], axis=0)

# Load HDIAG profile (isotropic)
f_profile = Dataset(args.build_dir + '/saber/test/testdata/doc_9/1-1_diag_plot.nc', 'r', format='NETCDF4')
fit_hor_multi = f_profile['stream_function']['cor1']['fit_hor'][0,:,:]
fit_hor_multi_detail = f_profile['stream_function']['cor1']['fit_detail_hor'][:,0,:,:]

# Load HDIAG map (isotropic)
f_map = Dataset(args.build_dir + '/saber/test/testdata/doc_6/cor_rh.nc', 'r', format='NETCDF4')
lon_hdiag = f_map['lon'][0,:]
lat_hdiag = f_map['lat'][:,0]
rh_hdiag = f_map['stream_function'][0,:,:]/1000.0
rh_hdiag_cyclic, lon_hdiag_cyclic = add_cyclic_point(rh_hdiag, coord=lon_hdiag)

# Load HDIAG map (anisotropic)
f_map = Dataset(args.build_dir + '/saber/test/testdata/doc_7/cor_rh1.nc', 'r', format='NETCDF4')
rh1_hdiag = f_map['stream_function'][0,:,:]/1000.0
f_map = Dataset(args.build_dir + '/saber/test/testdata/doc_7/cor_rh2.nc', 'r', format='NETCDF4')
rh2_hdiag = f_map['stream_function'][0,:,:]/1000.0
f_map = Dataset(args.build_dir + '/saber/test/testdata/doc_7/cor_rhc.nc', 'r', format='NETCDF4')
rhc_hdiag = f_map['stream_function'][0,:,:]

# Load NICAS diracs (isotropic)
f_diracs = Dataset(args.build_dir + '/saber/test/testdata/doc_6/dirac_nicas.nc', 'r', format='NETCDF4')
lon_diracs = f_diracs['lon'][0,:]
lat_diracs = f_diracs['lat'][:,0]
diracs_iso = f_diracs['stream_function'][0,:,:]
diracs_iso_cyclic, lon_diracs_cyclic = add_cyclic_point(diracs_iso, coord=lon_diracs)

# Load NICAS steps (isotropic)
f_steps = Dataset(args.build_dir + '/saber/test/testdata/doc_6/1-1_nicas_steps_local_000001-000001.nc', 'r', format='NETCDF4')
steps_iso = f_steps['stream_function']

# Load NICAS diracs (anisotropic)
f_diracs = Dataset(args.build_dir + '/saber/test/testdata/doc_7/dirac_nicas.nc', 'r', format='NETCDF4')
diracs_ani = f_diracs['stream_function'][0,:,:]
diracs_ani_cyclic, lon_diracs_cyclic = add_cyclic_point(diracs_ani, coord=lon_diracs)

# Load NICAS steps (anisotropic)
f_steps = Dataset(args.build_dir + '/saber/test/testdata/doc_7/1-1_nicas_steps_local_000001-000001.nc', 'r', format='NETCDF4')
steps_ani = f_steps['stream_function']

# Load NICAS diracs (isotropic - low resolution)
f_nicas_sampling = Dataset(args.build_dir + '/saber/test/testdata/doc_8/1-1_nicas_grids_local_000001-000001.nc', 'r', format='NETCDF4')
lon_sa_lr = f_nicas_sampling['stream_function']['lon_sa'][:]*180.0/math.pi
lat_sa_lr = f_nicas_sampling['stream_function']['lat_sa'][:]*180.0/math.pi

# Load NICAS diracs (mask)
f_diracs = Dataset(args.build_dir + '/saber/test/testdata/doc_10/dirac_nicas.nc', 'r', format='NETCDF4')
diracs_mask_mask = f_diracs['stream_function'][0,:,:]
diracs_mask_cyclic, lon_diracs_cyclic = add_cyclic_point(diracs_mask_mask, coord=lon_diracs)

####################################################################################################
# PREPROCESSING ####################################################################################
####################################################################################################

# Bounds
lon_plot = -45.0
lat_plot = -0.9
rh_min = 2000.0
rh_max = 4300.0
steps_min = 0.0
steps_max = 1.1
lon_steps = steps_iso['steps_5']['lon'][:]
lat_steps = steps_iso['steps_5']['lat'][:]
lon_min = np.min(lon_steps)
lon_max = np.max(lon_steps)
lat_min = np.min(lat_steps)
lat_max = np.max(lat_steps)

# Sizes
nx = lon_hdiag.shape[0]
ny = lat_hdiag.shape[0]

# rh on Sc2 sampling
rh_c2a = []
rh1_c2a = []
rh2_c2a = []
rhc_c2a = []
for ic2a in range(0, nc2a):
    ix_c2a = -1
    for ix in range(0, nx):
        lon_tmp = lon_hdiag[ix]
        if lon_tmp > 180.0:
            lon_tmp -= 360.0
        if abs(lon_tmp-lon_local[ic2a,0])<1.0e-6:
            ix_c2a = ix
            break
    if ix_c2a==-1:
        print('Cannot find ix for ' + str(lon_local[ic2a,0]))
        print(lon_hdiag)
        exit()
    iy_c2a = -1
    for iy in range(0, ny):
        if abs(lat_hdiag[iy]-lat_local[ic2a,0])<1.0e-6:
            iy_c2a = iy
            break
    if iy_c2a==-1:
        print('Cannot find iy for ' + str(lat_local[ic2a,0]))
        print(lat_hdiag)
        exit()
    rh_c2a.append(rh_hdiag[iy_c2a,ix_c2a])
    rh1_c2a.append(rh1_hdiag[iy_c2a,ix_c2a])
    rh2_c2a.append(rh2_hdiag[iy_c2a,ix_c2a])
    rhc_c2a.append(rhc_hdiag[iy_c2a,ix_c2a])

# Plot indices
ic2a_plot = [38, 150, -1]
for ic2a in range(0, nc2a):
    if (abs(lon_local[ic2a,0]-lon_plot)<1.0e-1) and (abs(lat_local[ic2a,0]-lat_plot)<1.0e-1):
        print(str(lon_local[ic2a,0]) + " / " + str(lat_local[ic2a,0]))
        ic2a_plot[2] = ic2a
        break
if ic2a_plot[2]==-1:
    print('Cannot find ic2a_plot[2]')
    exit()

# Eigen decomposition
D11 = rh1_c2a[ic2a_plot[2]]**2
D22 = rh2_c2a[ic2a_plot[2]]**2
D12 = rhc_c2a[ic2a_plot[2]]
tr = D11+D22
det = D11*D22*(1.0-D12**2)
l1 = 0.5*tr+math.sqrt(0.25*tr**2-det)
l2 = 0.5*tr-math.sqrt(0.25*tr**2-det)
v11 = l1-D22
v21 = math.sqrt(D11*D22)*D12
n1 = math.sqrt(v11**2+v21**2)
v11 = v11/n1
v21 = v21/n1
v12 = l2-D22
v22 = math.sqrt(D11*D22)*D12
n2 = math.sqrt(v12**2+v22**2)
v12 = v12/n2
v22 = v22/n2

####################################################################################################
# PLOTS ############################################################################################
####################################################################################################

print('####################################################################################################')
print('# Plot location: ' + str(lon_local[ic2a_plot[2],0]) + ' / ' + str(lat_local[ic2a_plot[2],0]))
print('####################################################################################################')

if args.mode == 'hdiag_sc1' or args.mode == 'all':
    # Sc1 sampling
    filename = 'hdiag_sc1'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
    plt.plot(lon_iso[0,0,:], lat_iso[0,0,:], color='lightblue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'hdiag_sc2' or args.mode == 'all':
    # Sc2 sampling, general
    filename = 'hdiag_sc2'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
    plt.plot(lon_iso[0,0,:], lat_iso[0,0,:], color='lightblue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
    plt.plot(lon_local[:,0], lat_local[:,0], color='red', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

    # Sc2 sampling, point A
    filename = 'hdiag_sc2_a'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ic2a = ic2a_plot[0]
    plt.plot(lon_iso[0,0,:], lat_iso[0,0,:], color='lightblue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
    plt.plot(lon_local[ic2a,1:], lat_local[ic2a,1:], color='blue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
    plt.plot(lon_local[ic2a,0], lat_local[ic2a,0], color='red', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

    # Sc2 sampling, point B
    filename = 'hdiag_sc2_b'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    plt.plot(lon_iso[0,0,:], lat_iso[0,0,:], color='lightblue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
    ic2a = ic2a_plot[1]
    plt.plot(lon_local[ic2a,1:], lat_local[ic2a,1:], color='blue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
    plt.plot(lon_local[ic2a,0], lat_local[ic2a,0], color='red', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

    # Sc2 sampling, point C
    filename = 'hdiag_sc2_c'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    plt.plot(lon_iso[0,0,:], lat_iso[0,0,:], color='lightblue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
    ic2a = ic2a_plot[2]
    plt.plot(lon_local[ic2a,1:], lat_local[ic2a,1:], color='blue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
    plt.plot(lon_local[ic2a,0], lat_local[ic2a,0], color='red', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'hdiag_sc3' or args.mode == 'all':
    # Sc3 sampling (isotropic)
    nc1sub = sum(~lon_local[ic2a_plot[2],:].mask)-1
    ic1_vec = []
    for ic1sub in range(0, nc1sub):
        lon_sub = lon_local[ic2a_plot[2],ic1sub+1]
        lat_sub = lat_local[ic2a_plot[2],ic1sub+1]
        for ic1 in range(0, nc1):
            if (lon_iso[0,0,ic1]==lon_sub) and (lat_iso[0,0,ic1]==lat_sub):
                break
        ic1_vec.append(ic1)
    delta = 60.0
    lon_min = lon_local[ic2a_plot[2],0]-delta
    lon_max = lon_local[ic2a_plot[2],0]+delta
    lat_min = lat_local[ic2a_plot[2],0]-delta
    lat_max = lat_local[ic2a_plot[2],0]+delta
    for ic3 in range(0, min(nc3_iso, 10)):
        filename = 'hdiag_sc3_' + str(ic3)
        print('Working on ' + filename)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        plt.plot(lon_iso[0,0,:], lat_iso[0,0,:], color='lightblue', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
        for ic1sub in range(0,nc1sub-1):
            lon_sub = lon_local[ic2a_plot[2],ic1sub+1]
            lat_sub = lat_local[ic2a_plot[2],ic1sub+1]
            plt.plot(lon_iso[0,ic3,ic1_vec[ic1sub]], lat_iso[0,ic3,ic1_vec[ic1sub]], color='green', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
            plt.plot([lon_sub,lon_iso[0,ic3,ic1_vec[ic1sub]]], [lat_sub,lat_iso[0,ic3,ic1_vec[ic1sub]]], color='green', linewidth=1, transform=ccrs.PlateCarree())
        plt.plot(lon_local[ic2a_plot[2],1:], lat_local[ic2a_plot[2],1:], color='blue', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
        plt.plot(lon_local[ic2a_plot[2],0], lat_local[ic2a_plot[2],0], color='red', linewidth=0, markersize=16, marker='.', transform=ccrs.PlateCarree())
        ax.set_extent([lon_min,lon_max,lat_min,lat_max], ccrs.PlateCarree())
        plt.savefig(filename + '.jpg', format='jpg', dpi=300)
        plt.close()
        subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'hdiag_sc3-sc4' or args.mode == 'all':
    # Sc3 sampling (anisotropic)
    nc1sub = sum(~lon_local[ic2a_plot[2],:].mask)-1
    ic1_vec = []
    for ic1sub in range(0, nc1sub):
        lon_sub = lon_local[ic2a_plot[2],ic1sub+1]
        lat_sub = lat_local[ic2a_plot[2],ic1sub+1]
        for ic1 in range(0, nc1):
            if (lon_iso[0,0,ic1]==lon_sub) and (lat_iso[0,0,ic1]==lat_sub):
                break
        ic1_vec.append(ic1)
    delta = 60.0
    lon_min = lon_local[ic2a_plot[2],0]-delta
    lon_max = lon_local[ic2a_plot[2],0]+delta
    lat_min = lat_local[ic2a_plot[2],0]-delta
    lat_max = lat_local[ic2a_plot[2],0]+delta
    for ic3 in range(0, min(nc3_iso, 10)):
        if (ic3==3 or ic3==8):
            for ic4 in range(0, nc4_ani):
                filename = 'hdiag_sc3-sc4_' + str(ic3) + '-' + str(ic4)
                print('Working on ' + filename)
                ax = plt.axes(projection=ccrs.PlateCarree())
                ax.coastlines()
                plt.plot(lon_iso[0,0,:], lat_iso[0,0,:], color='lightblue', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
                for ic1sub in range(0,nc1sub):
                    lon_sub = lon_local[ic2a_plot[2],ic1sub+1]
                    lat_sub = lat_local[ic2a_plot[2],ic1sub+1]
                    plt.plot(lon_ani[ic4,ic3,ic1_vec[ic1sub]], lat_ani[ic4,ic3,ic1_vec[ic1sub]], color='green', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
                    plt.plot([lon_sub,lon_ani[ic4,ic3,ic1_vec[ic1sub]]], [lat_sub,lat_ani[ic4,ic3,ic1_vec[ic1sub]]], color='green', linewidth=1, transform=ccrs.PlateCarree())
                plt.plot(lon_local[ic2a_plot[2],1:], lat_local[ic2a_plot[2],1:], color='blue', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
                plt.plot(lon_local[ic2a_plot[2],0], lat_local[ic2a_plot[2],0], color='red', linewidth=0, markersize=16, marker='.', transform=ccrs.PlateCarree())
                ax.set_extent([lon_min,lon_max,lat_min,lat_max], ccrs.PlateCarree())
                plt.savefig(filename + '.jpg', format='jpg', dpi=300)
                plt.close()
                subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'fit_iso' or args.mode == 'all':
    # Raw diagnostic (isotropic)
    filename = 'fit_iso_raw'
    print('Working on ' + filename)
    ax = plt.axes()
    ax.axhline(y=0, color='gray', linestyle='--')
    ax.axhline(y=1, color='gray', linestyle='--')
    ax.set_xlim(0, 4200.0)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel('Horizontal distance (km)')
    ax.set_ylabel('Correlation')
    ax.plot(disth_iso, raw_hor_iso[0,:], 'k', linewidth=2)
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

    # Fit (isotropic)
    filename = 'fit_iso_fit'
    print('Working on ' + filename)
    ax = plt.axes()
    ax.axhline(y=0, color='gray', linestyle='--')
    ax.axhline(y=1, color='gray', linestyle='--')
    ax.set_xlim(0, 4200.0)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel('Horizontal distance (km)')
    ax.set_ylabel('Correlation')
    ax.plot(disth_iso, raw_hor_iso[0,:], 'k', linewidth=2)
    ax.plot(disth_iso, fit_hor_iso[0,:], 'r', linewidth=2)
    ax.annotate('', xy=(rh_c2a[ic2a_plot[2]], 0.0), xycoords='data', xytext=(-10, 0.0), textcoords='data', arrowprops=dict(arrowstyle='-|>', color='blue', linewidth=4))
    ax.annotate('Length-scale: ' + str(int(rh_c2a[ic2a_plot[2]])) + ' km', xy=(2200.0, -0.07), textcoords='data', color='blue')
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'fit_ani' or args.mode == 'all':
    # Raw diagnostic (anisotropic)
    filename = 'fit_ani_raw'
    print('Working on ' + filename)
    ax = plt.axes(projection='polar')
    r, theta = np.meshgrid(disth_ani, angular_sector_cyclic)
    levels = np.linspace(-1.0, 1.0, 21)
    ax.contourf(theta, r, raw_hor_ani_cyclic, levels=levels, cmap='coolwarm')
    ax.contour(theta, r, raw_hor_ani_cyclic, levels=levels, linewidths=[0.5], colors=['k'])
    sm = cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=-1.0, vmax=1.0))
    sm.set_array([])
    plt.colorbar(sm, orientation='vertical', pad=0.04)
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

    # Fit (anisotropic)
    hfac = 0.8
    filename = 'fit_ani_fit'
    print('Working on ' + filename)
    ax = plt.axes(projection='polar')
    r, theta = np.meshgrid(disth_ani, angular_sector_cyclic)
    levels = np.linspace(-1.0, 1.0, 21)
    ax.contourf(theta, r, fit_hor_ani_cyclic, levels=levels, cmap='coolwarm')
    ax.contour(theta, r, fit_hor_ani_cyclic, levels=levels, linewidths=0.5, colors='k')
    t = np.linspace(0, 2*math.pi, 100)
    angle = np.linspace(math.atan2(-v12, v22), math.atan2(-v12, v22), len(t))
    circle = np.linspace(rh_c2a[ic2a_plot[2]], rh_c2a[ic2a_plot[2]], len(t))*hfac
    ell = np.sqrt(l1*l2/(l1*np.sin(t-angle)**2+l2*np.cos(t-angle)**2))*hfac
    plt.plot(t, circle, 'blue', linestyle='--')
    plt.plot(t, ell, 'green')
    ax.annotate('Length-scale', xy=(np.radians(115.0), 3000), textcoords='data', color='blue')
    ax.annotate('Tensor', xy=(np.radians(110.0), 1900), textcoords='data', color='green')
    sm = cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=-1.0, vmax=1.0))
    sm.set_array([])
    plt.colorbar(sm, orientation='vertical', pad=0.04)
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'fit_multi' or args.mode == 'all':
    # Fit (multi-component)
    filename = 'fit_multi_fit_1'
    print('Working on ' + filename)
    ax = plt.axes()
    ax.axhline(y=0, color='gray', linestyle='--')
    ax.axhline(y=1, color='gray', linestyle='--')
    ax.set_xlim(0, 4200.0)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel('Horizontal distance (km)')
    ax.set_ylabel('Correlation')
    ax.plot(disth_iso, raw_hor_iso[0,:], 'k', linewidth=2)
    ax.plot(disth_iso, fit_hor_multi_detail[0,0,:], 'g', linestyle='-.', linewidth=2)
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

    # Fit (multi-component)
    filename = 'fit_multi_fit_2'
    print('Working on ' + filename)
    ax = plt.axes()
    ax.axhline(y=0, color='gray', linestyle='--')
    ax.axhline(y=1, color='gray', linestyle='--')
    ax.set_xlim(0, 4200.0)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel('Horizontal distance (km)')
    ax.set_ylabel('Correlation')
    ax.plot(disth_iso, raw_hor_iso[0,:], 'k', linewidth=2)
    ax.plot(disth_iso, fit_hor_multi_detail[0,0,:], 'g', linestyle='-.', linewidth=2)
    ax.plot(disth_iso, fit_hor_multi_detail[1,0,:]-fit_hor_multi_detail[0,0,:], 'g', linestyle='--', linewidth=2)
    ax.plot(disth_iso, fit_hor_multi[0,:], 'r', linewidth=2)
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])


if args.mode == 'rh_sampling' or args.mode == 'all':
    # Horizontal radius sampling
    filename = 'rh_sampling_iso'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
    cmap = cm.get_cmap('jet')
    for ic2a in range(0, nc2a):
        rh_norm = (rh_c2a[ic2a]-rh_min)/(rh_max-rh_min)
        plt.plot(lon_local[ic2a,0], lat_local[ic2a,0], color=cmap(rh_norm), linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree())
    sm = cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=rh_min, vmax=rh_max))
    sm.set_array([])
    plt.colorbar(sm, orientation='horizontal', pad=0.06)
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'rh_map' or args.mode == 'all':
    # Horizontal radius map
    filename = 'rh_map'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
    levels = np.linspace(rh_min, rh_max, 100)
    ax.contourf(lon_hdiag_cyclic, lat_hdiag, rh_hdiag_cyclic, levels=levels, cmap='jet', transform=ccrs.PlateCarree())
    sm = cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=rh_min, vmax=rh_max))
    sm.set_array([])
    plt.colorbar(sm, orientation='horizontal', pad=0.06)
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'nicas_iso' or args.mode == 'all':
    # NICAS diracs (isotropic)
    filename = 'nicas_iso_diracs'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
    levels = np.linspace(-1.0, 1.0, 100)
    ax.contourf(lon_diracs_cyclic, lat_diracs, diracs_iso_cyclic, levels=levels, cmap='jet', transform=ccrs.PlateCarree())
    sm = cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=-1.0, vmax=1.0))
    sm.set_array([])
    plt.colorbar(sm, orientation='horizontal', pad=0.06)
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

    # NICAS steps (isotropic)
    for i in range(1, 8):
        steps = steps_iso['steps_' + str(i)]
        lon_steps = steps['lon'][:]
        lat_steps = steps['lat'][:]
        values_steps = steps['values'][:]
        n = len(lon_steps)

        filename = 'nicas_iso_steps_' + str(i)
        print('Working on ' + filename)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_extent([lon_min,lon_max,lat_min,lat_max], ccrs.PlateCarree())
        if i<6:
            cmap = cm.get_cmap('rainbow')
            for j in range(0, n):
                value = (values_steps[j]-steps_min)/(steps_max-steps_min)
                plt.plot(lon_steps[j], lat_steps[j], color=cmap(value), linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
        else:
            levels = np.linspace(steps_min, steps_max, 11)
            ax.tricontourf(lon_steps, lat_steps, values_steps, levels=levels, cmap='rainbow')
        sm = cm.ScalarMappable(cmap='rainbow', norm=plt.Normalize(vmin=steps_min, vmax=steps_max))
        sm.set_array([])
        plt.colorbar(sm, orientation='vertical', pad=0.04)
        plt.savefig(filename + '.jpg', format='jpg', dpi=300)
        plt.close()
        subprocess.run(['mogrify', '-trim', filename + '.jpg'])

    # NICAS diracs + rh (isotropic)
    filename = 'nicas_iso_diracs_rh'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
    levels = np.linspace(rh_min, rh_max, 100)
    ax.contourf(lon_hdiag_cyclic, lat_hdiag, rh_hdiag_cyclic, levels=levels, cmap='jet', transform=ccrs.PlateCarree())
    sm = cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=rh_min, vmax=rh_max))
    sm.set_array([])
    plt.colorbar(sm, orientation='horizontal', pad=0.06)
    levels = np.linspace(0.2, 1.0, 5)
    ax.contour(lon_diracs_cyclic, lat_diracs, diracs_iso_cyclic, levels=levels, linewidths=1.0, colors='k', transform=ccrs.PlateCarree())
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'nicas_ani' or args.mode == 'all':
    # NICAS diracs (anisotropic)
    filename = 'nicas_ani_diracs'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
    levels = np.linspace(-1.0, 1.0, 100)
    ax.contourf(lon_diracs_cyclic, lat_diracs, diracs_ani_cyclic, levels=levels, cmap='jet', transform=ccrs.PlateCarree())
    sm = cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=-1.0, vmax=1.0))
    sm.set_array([])
    plt.colorbar(sm, orientation='horizontal', pad=0.06)
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

    # NICAS steps (anisotropic)
    for i in range(1, 8):
        steps = steps_ani['steps_' + str(i)]
        lon_steps = steps['lon'][:]
        lat_steps = steps['lat'][:]
        values_steps = steps['values'][:]
        n = len(lon_steps)

        filename = 'nicas_ani_steps_' + str(i)
        print('Working on ' + filename)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_extent([lon_min,lon_max,lat_min,lat_max], ccrs.PlateCarree())
        if i<6:
            cmap = cm.get_cmap('rainbow')
            for j in range(0, n):
                value = (values_steps[j]-steps_min)/(steps_max-steps_min)
                plt.plot(lon_steps[j], lat_steps[j], color=cmap(value), linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
        else:
            levels = np.linspace(steps_min, steps_max, 11)
            ax.tricontourf(lon_steps, lat_steps, values_steps, levels=levels, cmap='rainbow')
        sm = cm.ScalarMappable(cmap='rainbow', norm=plt.Normalize(vmin=steps_min, vmax=steps_max))
        sm.set_array([])
        plt.colorbar(sm, orientation='vertical', pad=0.04)
        plt.savefig(filename + '.jpg', format='jpg', dpi=300)
        plt.close()
        subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'nicas_iso-lr' or args.mode == 'all':
    # NICAS sampling + rh (isotropic- low resolution)
    filename = 'nicas_iso-lr_sampling_rh'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
    levels = np.linspace(rh_min, rh_max, 100)
    ax.contourf(lon_hdiag_cyclic, lat_hdiag, rh_hdiag_cyclic, levels=levels, cmap='jet', transform=ccrs.PlateCarree())
    sm = cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=rh_min, vmax=rh_max))
    sm.set_array([])
    plt.colorbar(sm, orientation='horizontal', pad=0.06)
    plt.plot(lon_sa_lr, lat_sa_lr, color='k', linewidth=0, markersize=3, marker='.', transform=ccrs.PlateCarree())
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'nicas_mask' or args.mode == 'all':
    # NICAS diracs (mask)
    filename = 'nicas_mask_diracs'
    print('Working on ' + filename)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
    levels = np.linspace(-1.0, 1.0, 100)
    ax.contourf(lon_diracs_cyclic, lat_diracs, diracs_mask_cyclic, levels=levels, cmap='jet', transform=ccrs.PlateCarree())
    sm = cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=-1.0, vmax=1.0))
    sm.set_array([])
    plt.colorbar(sm, orientation='horizontal', pad=0.06)
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

if args.mode == 'resolution' or args.mode == 'all':
    # Resolution definition
    filename = 'resolution'
    print('Working on ' + filename)
    fig, ax = plt.subplots(ncols=4, figsize=(18,4))
    fig.subplots_adjust(wspace=0.3)
    resols = [4,6,8,10]
    for ir in range(0, len(resols)):
        # Define function
        x = np.linspace(0.0, 1.0, resols[ir])
        y = np.zeros((resols[ir]))
        for i in range(0, resols[ir]):
            if x[i] < 0.5:
                y[i] = -8.0*x[i]**5+8.0*x[i]**4+5.0*x[i]**3-20.0/3.0*x[i]**2+1.0
            elif x[i] < 1.0:
                y[i] = 8.0/3.0*x[i]**5-8.0*x[i]**4+5.0*x[i]**3+20.0/3.0*x[i]**2-10.0*x[i]+4.0-1.0/(3.0*x[i])
            else:
                y[i] = 0.0
        ax[ir].axhline(y=0, color='gray', linestyle='--')
        ax[ir].axhline(y=1, color='gray', linestyle='--')
        ax[ir].set_xlim(0.0, 1.1)
        ax[ir].set_ylim(0.0, 1.1)
        ax[ir].set_xlabel('Normalized horizontal distance')
        ax[ir].set_ylabel('GC99 at resolution ' + str(resols[ir]))
        ax[ir].plot(x, y, 'r', linewidth=2, marker='o')
    plt.savefig(filename + '.jpg', format='jpg', dpi=300)
    plt.close()
    subprocess.run(['mogrify', '-trim', filename + '.jpg'])

print('####################################################################################################')
