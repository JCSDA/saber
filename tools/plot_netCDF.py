#!/usr/bin/env python3
'''
(C) Copyright 2024 UCAR
(C) Crown Copyright 2024 Met Office

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 '''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from netCDF4 import Dataset as NetCDFFile

import cartopy.crs as ccrs
import pathlib
import argparse

'''
Python script for reading and plotting data from a global NetCDF File.
Suited to plot QUENCH fields or outputs of oops util::writeFieldSets.
(This script expects file have arrays of lats/lons in list of variables)
Does not currently support plotting for Gaussian grids.

For a list of variables in the file use the command:

  ncdump -h <FILE>

This will also print information on the sizes and shapes of data in the file.

Call this method with:

  python3 plot_netCDF.py <netCDF_file.nc> <variable> --<kwargs>

Get more information on options with:

  python3 plot_netCDF.py --help

In summary, the command line arguments specify which variable, level, etc will
plotted. See the notes in the `help` argument.
'''

parser = argparse.ArgumentParser()

parser.add_argument("file", help="relative path to netCDF file")
parser.add_argument("variable", help="variable in netCDF file to plot")
parser.add_argument("--level", type=int, default=0, help="vertical level to plot (default=0)")
parser.add_argument("--center-on-zero", action=argparse.BooleanOptionalAction, help="flag for centering variable plot color scale at zero with normalized bounds between [-1, 1]")
parser.add_argument("--ax-title", help="plot title, default is `<variable> - level #`")
parser.add_argument("--cb-title", default='', help="add given title for colorbar")
parser.add_argument("--output-prefix", default='', help="add prefix to output filename")
parser.add_argument("--contours", action=argparse.BooleanOptionalAction)
parser.add_argument("--ncontours", default=25, help="number of contour levels")
parser.add_argument("--show-point-positions", action=argparse.BooleanOptionalAction)
parser.add_argument("--save", action=argparse.BooleanOptionalAction, default=False, help="save plot figure as .jpg")

#%% =========================================================
# Read input parameters

args = parser.parse_args()

ncFilePath = pathlib.Path(args.file)

variable = args.variable
level = args.level
center_on_zero = (args.center_on_zero is not None)
ax_title = f'{variable} - level {level}' if not args.ax_title else args.ax_title
cb_title = args.cb_title
prefix = args.output_prefix  # For output filename
use_contours = (args.contours is not None)
n_contours = int(args.ncontours)
show_point_pos = (args.show_point_positions is not None)
save = args.save

print('##########################################################################################')
print("center_on_zero:", center_on_zero)
print("contours:", use_contours)

# Read NetCDF File:
nc = NetCDFFile(ncFilePath, 'r', format='NETCDF4')

# get lists of lats & lons from file (these are 1D arrays)
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]


# Get variable from file
# for <variable> get all data for <level> (1D array)
var = nc[variable][:,level]

# Levels
var_min = np.min(var)
var_max = np.max(var)
bound = max(np.abs(var_min), np.abs(var_max))
if bound == 0:
    print(f"Bounds for {variable} are zero")
    print("Be ready for trouble...")
    print("Do NOT trust plot produced")
if center_on_zero:
    if bound == 0:
        bound = 1
    var_min = - bound
    var_max = bound

cmap = 'coolwarm' if center_on_zero else 'plasma'
norm = plt.Normalize(vmin=var_min, vmax=var_max)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)

if use_contours:
    try:
        # convert from [0, 360) -> (-180, 180]
        for i, l in enumerate(lon):
          if l > 180:
            lon[i] -= 360
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([-180.0,180.0,-90.0,90.0], ccrs.PlateCarree())
        print("n_contours = ", n_contours)
        plt.tricontourf(lon, lat, var, n_contours,
                              transform=ccrs.PlateCarree(),
                              cmap=cmap,
                              norm=norm)
    except ValueError:
        use_contours = False

if not use_contours:
  ax = plt.axes(projection=ccrs.PlateCarree())
  plt.scatter(
      x=lon,
      y=lat,
      c=var,
      s=8,
      alpha=1,
      transform=ccrs.PlateCarree(),
      cmap=cmap,
      norm=norm
  )

if show_point_pos:
    ax.plot(lon, lat, 'ko', transform=ccrs.PlateCarree(), ms=1)

ax.coastlines()
ax.gridlines()

ax.set_title(ax_title)
cb = plt.colorbar(sm, orientation='horizontal', pad=0.06, ax=ax)
if len(cb_title):
    cb.ax.set_xlabel(cb_title)

if(save):
  filename = ncFilePath.with_suffix(".jpg")
  filename = pathlib.Path("./") / f"{prefix}{variable}_level{level}_{filename.name}"
  plt.savefig(filename, format='jpg', dpi=300)
  print(f"Figure saved as: {filename}")
plt.show()
print('##########################################################################################')
