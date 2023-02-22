#!/usr/bin/env python3
# This script can be used to generate input covariance files for GSI-GEOS with unit variance and
# homogeneous length-scales, which is useful to optimize the scaling factor between GSI and BUMP.

# Modules
import os
import sys
import netCDF4 as nc
import numpy as np

# Input file
inputField = '../test/testdata/gsi-coeffs-gmao-global-l72x72y46.nc4'

# Output files
outputField = []
outputField.append('../test/testdata/gsi-coeffs-gmao-global-l72x72y46_opt1.nc4')
outputField.append('../test/testdata/gsi-coeffs-gmao-global-l72x72y46_opt2.nc4')
outputField.append('../test/testdata/gsi-coeffs-gmao-global-l72x72y46_opt3.nc4')

# Variance fields that should be reinitialized to one
stddev = []
stddev.append(1.0)
stddev.append(1.0)
stddev.append(1.0)
variance = ['ps','sf','vp','t','q','qi','ql','qr','qs','oz','cw','sst']

# Optimization horizontal length-scale
Lh = []
Lh.append(1000000.0)
Lh.append(1500000.0)
Lh.append(-1.0)
horizontal = ['hps','hsf','hvp','ht','hq','hqi','hql','hqr','hqs','hoz','hcw','hsst']

# Optimization vertical length-scale
Lv = []
Lv.append(0.3)
Lv.append(0.5)
Lv.append(-1.0)
vertical = ['vsf','vvp','vt','vq','vqi','vql','vqr','vqs','voz','vcw']

# Create optimization files
for i in range(0,len(outputField)):
  with nc.Dataset(inputField, "r") as src, nc.Dataset(outputField[i], "w") as dst:
    # Copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)

    # Copy dimensions
    for name, dimension in src.dimensions.items():
      dst.createDimension(
        name, (len(dimension) if not dimension.isunlimited() else None))

    # Copy all file data, 
    for name, variable in src.variables.items():
      x = dst.createVariable(name, variable.datatype, variable.dimensions)
      if (name in variance) and (stddev[i] > 0.0):
        dst[name][:] = stddev[i]
      elif (name in horizontal) and (Lh[i] > 0.0):
        dst[name][:] = Lh[i]
      elif (name in vertical) and (Lv[i] > 0.0):
        dst[name][:] = Lv[i]
      else:
        dst[name][:] = src[name][:]

      # Copy variable attributes all at once via dictionary
      dst[name].setncatts(src[name].__dict__)
