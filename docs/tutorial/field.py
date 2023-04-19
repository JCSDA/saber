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
parser.add_argument('ncfile', help='NetCDF file')
args = parser.parse_args()

# Hard-coded parameters
variable = 'stream_function'
level = 0

# Load data
f = Dataset(args.ncfile, 'r', format='NETCDF4')
lon = f['lon'][0,:]
lat = f['lat'][:,0]
var = f[variable][level,:,:]
var_cyclic, lon_cyclic = add_cyclic_point(var, coord=lon)

# Levels
var_min = np.min(var)
var_max = np.max(var)
levels = np.linspace(var_min, var_max, 100)

# Plot
filename = args.ncfile + ".jpg"
print('####################################################################################################')
print('Output file: ' + filename)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
ax.contourf(lon_cyclic, lat, var_cyclic, levels=levels, cmap='jet', transform=ccrs.PlateCarree())
sm = cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=var_min, vmax=var_max))
sm.set_array([])
plt.colorbar(sm, orientation='horizontal', pad=0.06)
plt.savefig(filename, format='jpg', dpi=300)
plt.close()
subprocess.run(['mogrify', '-trim', filename])
print('####################################################################################################')
