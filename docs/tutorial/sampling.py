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
parser.add_argument('data_dir', help='Data directory')
parser.add_argument('avg_center', help='Averaging center index')
args = parser.parse_args()

####################################################################################################
# READ NETCDF FILES ################################################################################
####################################################################################################

# Load HDIAG sampling
f_sampling = Dataset(args.data_dir + '/1-1_sampling_grids_local_000001-000001.nc', 'r', format='NETCDF4')
lon = f_sampling['lon'][0,:,:,:]
lat = f_sampling['lat'][0,:,:,:]
inhomogeneous = ('lon_local' in f_sampling.variables)
if inhomogeneous:
  lon_local = f_sampling['lon_local'][0,:,:]
  lat_local = f_sampling['lat_local'][0,:,:]
nc4 = lon.shape[0]
nc3 = lon.shape[1]
nc1 = lon.shape[2]
if inhomogeneous:
  nc2a = lon_local.shape[0]

####################################################################################################
# PLOTS ############################################################################################
####################################################################################################

print('####################################################################################################')
ic2a = int(args.avg_center)
if inhomogeneous:
  print('# Plot location: ' + str(lon_local[ic2a,0]) + ' / ' + str(lat_local[ic2a,0]))
  print('####################################################################################################')

# Sc1 sampling
filename = args.data_dir + '/hdiag_sc1.jpg'
print('Output file: ' + filename)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
plt.plot(lon[0,0,:], lat[0,0,:], color='lightblue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree())
plt.savefig(filename, format='jpg', dpi=300)
plt.close()
subprocess.run(['mogrify', '-trim', filename])

if inhomogeneous:
  # Sc2 sampling, general 
  filename = args.data_dir + '/hdiag_sc2_global.jpg'
  print('Output file: ' + filename) 
  ax = plt.axes(projection=ccrs.PlateCarree()) 
  ax.coastlines() 
  ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree()) 
  plt.plot(lon[0,0,:], lat[0,0,:], color='lightblue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree()) 
  plt.plot(lon_local[:,0], lat_local[:,0], color='red', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree()) 
  plt.savefig(filename, format='jpg', dpi=300) 
  plt.close() 
  subprocess.run(['mogrify', '-trim', filename]) 
 
  # Sc2 sampling, local 
  filename = args.data_dir + '/hdiag_sc2_local.jpg'
  print('Output file: ' + filename)
  ax = plt.axes(projection=ccrs.PlateCarree()) 
  ax.coastlines() 
  plt.plot(lon[0,0,:], lat[0,0,:], color='lightblue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree()) 
  plt.plot(lon_local[ic2a,1:], lat_local[ic2a,1:], color='blue', linewidth=0, markersize=4, marker='.', transform=ccrs.PlateCarree()) 
  plt.plot(lon_local[ic2a,0], lat_local[ic2a,0], color='red', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree()) 
  plt.savefig(filename, format='jpg', dpi=300) 
  plt.close() 
  subprocess.run(['mogrify', '-trim', filename]) 
 
  # Sc3 sampling (isotropic) 
  nc1sub = sum(~lon_local[ic2a,:].mask)-1 
  ic1_vec = [] 
  for ic1sub in range(0, nc1sub): 
    lon_sub = lon_local[ic2a,ic1sub+1] 
    lat_sub = lat_local[ic2a,ic1sub+1] 
    for ic1 in range(0, nc1): 
      if (lon[0,0,ic1]==lon_sub) and (lat[0,0,ic1]==lat_sub): 
        break 
    ic1_vec.append(ic1) 
  delta = 60.0 
  lon_min = lon_local[ic2a,0]-delta 
  lon_max = lon_local[ic2a,0]+delta 
  lat_min = lat_local[ic2a,0]-delta 
  lat_max = lat_local[ic2a,0]+delta 
  for ic3 in range(0, min(nc3, 10)): 
    filename = args.data_dir + '/hdiag_sc3_' + str(ic3) + '.jpg'
    print('Output file: ' + filename) 
    ax = plt.axes(projection=ccrs.PlateCarree()) 
    ax.coastlines() 
    plt.plot(lon[0,0,:], lat[0,0,:], color='lightblue', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree()) 
    for ic1sub in range(0,nc1sub-1): 
      lon_sub = lon_local[ic2a,ic1sub+1] 
      lat_sub = lat_local[ic2a,ic1sub+1] 
      plt.plot(lon[0,ic3,ic1_vec[ic1sub]], lat[0,ic3,ic1_vec[ic1sub]], color='green', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree()) 
      plt.plot([lon_sub,lon[0,ic3,ic1_vec[ic1sub]]], [lat_sub,lat[0,ic3,ic1_vec[ic1sub]]], color='green', linewidth=1, transform=ccrs.PlateCarree()) 
    plt.plot(lon_local[ic2a,1:], lat_local[ic2a,1:], color='blue', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree()) 
    plt.plot(lon_local[ic2a,0], lat_local[ic2a,0], color='red', linewidth=0, markersize=16, marker='.', transform=ccrs.PlateCarree()) 
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], ccrs.PlateCarree()) 
    plt.savefig(filename, format='jpg', dpi=300) 
    plt.close() 
    subprocess.run(['mogrify', '-trim', filename]) 
 
  if nc4 > 1:
    # Sc3 sampling (anisotropic) 
    nc1sub = sum(~lon_local[ic2a,:].mask)-1 
    ic1_vec = [] 
    for ic1sub in range(0, nc1sub): 
      lon_sub = lon_local[ic2a,ic1sub+1] 
      lat_sub = lat_local[ic2a,ic1sub+1] 
      for ic1 in range(0, nc1): 
        if (lon[0,0,ic1]==lon_sub) and (lat[0,0,ic1]==lat_sub): 
          break 
      ic1_vec.append(ic1) 
    delta = 60.0 
    lon_min = lon_local[ic2a,0]-delta 
    lon_max = lon_local[ic2a,0]+delta 
    lat_min = lat_local[ic2a,0]-delta 
    lat_max = lat_local[ic2a,0]+delta 
    for ic3 in range(0, min(nc3, 10)): 
      if (ic3==3 or ic3==8): 
        for ic4 in range(0, nc4): 
          filename = args.data_dir + '/hdiag_sc3-sc4_' + str(ic3) + '-' + str(ic4)  + '.jpg'
          print('Output file: ' + filename) 
          ax = plt.axes(projection=ccrs.PlateCarree()) 
          ax.coastlines() 
          plt.plot(lon[0,0,:], lat[0,0,:], color='lightblue', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree()) 
          for ic1sub in range(0,nc1sub): 
            lon_sub = lon_local[ic2a,ic1sub+1] 
            lat_sub = lat_local[ic2a,ic1sub+1] 
            plt.plot(lon[ic4,ic3,ic1_vec[ic1sub]], lat[ic4,ic3,ic1_vec[ic1sub]], color='green', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree()) 
            plt.plot([lon_sub,lon[ic4,ic3,ic1_vec[ic1sub]]], [lat_sub,lat[ic4,ic3,ic1_vec[ic1sub]]], color='green', linewidth=1, transform=ccrs.PlateCarree()) 
          plt.plot(lon_local[ic2a,1:], lat_local[ic2a,1:], color='blue', linewidth=0, markersize=8, marker='.', transform=ccrs.PlateCarree()) 
          plt.plot(lon_local[ic2a,0], lat_local[ic2a,0], color='red', linewidth=0, markersize=16, marker='.', transform=ccrs.PlateCarree()) 
          ax.set_extent([lon_min,lon_max,lat_min,lat_max], ccrs.PlateCarree()) 
          plt.savefig(filename, format='jpg', dpi=300) 
          plt.close() 
          subprocess.run(['mogrify', '-trim', filename])

print('####################################################################################################')
