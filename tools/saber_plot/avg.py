#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import os

def avg(testdata, test, mpi, omp, suffix, testfig):
   """! Plot script for the "average" files produced by BUMP"""

   # Open file
   f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Get vertical unit
   vunit = f["vunit"][:]

   # Get number of levels
   nl0 = vunit.shape[0]

   for group in f.groups:
      # Get horizontal distance
      disth = f.groups[group]["disth"][:]

      # Get number of horizontal classes
      nc3 = disth.shape[0]

      # Get number of reduced levels
      nl0r = f.groups[group]["cor_hist"][:,:,:,:].shape[1]

      # Get number of bins
      nbins = f.groups[group]["cor_hist"][:,:,:,:].shape[3]

      # Read variable and bins
      l0rl0_to_l0 = f.groups[group]["l0rl0_to_l0"][:,:]
      m11_hist = f.groups[group]["m11_hist"][:,:,:,:]
      m11_bins = f.groups[group]["m11_bins"][:,:,:,:]
      m11m11_hist = f.groups[group]["m11m11_hist"][:,:,:,:]
      m11m11_bins = f.groups[group]["m11m11_bins"][:,:,:,:]
      m2m2_hist = f.groups[group]["m2m2_hist"][:,:,:,:]
      m2m2_bins = f.groups[group]["m2m2_bins"][:,:,:,:]
      m22_hist = f.groups[group]["m22_hist"][:,:,:,:]
      m22_bins = f.groups[group]["m22_bins"][:,:,:,:]
      cor_hist = f.groups[group]["cor_hist"][:,:,:,:]
      cor_bins = f.groups[group]["cor_bins"][:,:,:,:]

      # Plots
      for il0 in range(0, nl0):
         for jl0r in range(0, nl0r):
            for jc3 in range(0, nc3):
               if np.any(cor_hist[il0,jl0r,jc3,:] > 0.0):
                  # Plots
                  fig, ax = plt.subplots(ncols=5)
                  fig.subplots_adjust(hspace=2.0)

                  # Covariance
                  i = 0
                  ax[i].set_title("Covariance")
                  x = 0.5*(m11_bins[il0,jl0r,jc3,0:nbins]+m11_bins[il0,jl0r,jc3,1:nbins+1])
                  ax[i].plot(x, m11_hist[il0,jl0r,jc3,:], 'r-')

                  # Covariance squared
                  i = 1
                  ax[i].set_title("Covariance squared")
                  x = 0.5*(m11m11_bins[il0,jl0r,jc3,0:nbins]+m11m11_bins[il0,jl0r,jc3,1:nbins+1])
                  ax[i].plot(x, m11m11_hist[il0,jl0r,jc3,:], 'r-')

                  # Variance product
                  i = 2
                  ax[i].set_title("Variance product")
                  x = 0.5*(m2m2_bins[il0,jl0r,jc3,0:nbins]+m2m2_bins[il0,jl0r,jc3,1:nbins+1])
                  ax[i].plot(x, m2m2_hist[il0,jl0r,jc3,:], 'r-')

                  # Fourth-order moment
                  i = 3
                  ax[i].set_title("Fourth-order moment")
                  x = 0.5*(m22_bins[il0,jl0r,jc3,0:nbins]+m22_bins[il0,jl0r,jc3,1:nbins+1])
                  ax[i].plot(x, m22_hist[il0,jl0r,jc3,:], 'r-')

                  # Correlation
                  i = 4
                  ax[i].set_title("Correlation")
                  x = 0.5*(cor_bins[il0,jl0r,jc3,0:nbins]+cor_bins[il0,jl0r,jc3,1:nbins+1])
                  ax[i].plot(x, cor_hist[il0,jl0r,jc3,:], 'r-')

                  # Add title
                  fig.suptitle("Between levels " + '%.2e'%vunit[il0] + " and " + '%.2e'%vunit[l0rl0_to_l0[il0,jl0r]-1] + " for distance class " + '%.2e'%disth[jc3])

                  # Save and close figure
                  plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + group + "_" + str(il0+1) + "-" + str(l0rl0_to_l0[il0,jl0r]) + "-" + str(jc3+1) + ".jpg", format="jpg", dpi=300)
                  plt.close()
