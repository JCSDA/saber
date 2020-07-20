#!/usr/bin/env python3

import argparse
import Ngl, Nio
import numpy as np
import os

def lct_cor(testdata, test, mpi, omp, suffix):
   # Loop over files
   first = True
   for impi in range(0, int(mpi)):
      # Open file
      f = Nio.open_file(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + mpi.zfill(6) + "-" + format(impi+1, '06d') + ".nc")




test_2-1_lct_cor_000002-000001_01_01_01_01.nc


      # Get data
      nc1a = f.variables["raw"][:].shape[1]
      if first:
         nl0 = f.variables["raw"][:].shape[0]
         nl0r = f.variables["raw"][:].shape[2]
         nc3 = f.variables["raw"][:].shape[3]
         raw = np.zeros((nl0, 0, nl0r, nc3))
         fit = np.zeros((nl0, 0, nl0r, nc3))
         lon = np.zeros((0, nc3))
         lat = np.zeros((0, nc3))
         first = False
      for ic1a in range(0, nc1a):
         raw = np.append(raw, [f.variables["raw"][:,ic1a,:,:].data], axis=1)
         fit = np.append(fit, [f.variables["fit"][:,ic1a,:,:].data], axis=1)
         lon = np.append(lon, [f.variables["lon"][ic1a,:].data], axis=0)
         lat = np.append(lon, [f.variables["lon"][ic1a,:].data], axis=0)
