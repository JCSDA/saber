#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import os

def randomization(testdata, test, mpi, omp, suffix, testfig):
   """! Plot script for the "randomization" files produced by BUMP"""

   # Open file
   f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Get factors
   nefac = f["nefac"][:]

   # Get number of tests
   nfac = nefac.shape[0]

   # Get number of factors
   nfac = nefac.shape[0]
   x = np.linspace(1, nfac, nfac)

   # Get data, compute mean and standard deviation
   mse_mean = np.zeros((0, nfac))
   mse_mean = np.append(mse_mean, [np.mean(f["mse"][:,:], axis=1)], axis=0)
   mse_mean = np.append(mse_mean, [np.mean(f["mse_th"][:,:], axis=1)], axis=0)
   mse_std = np.zeros((0, nfac))
   mse_std = np.append(mse_std, [np.std(f["mse"][:,:], axis=1)], axis=0)
   mse_std = np.append(mse_std, [np.std(f["mse_th"][:,:], axis=1)], axis=0)

   # Plots
   fig, ax = plt.subplots()
   ax.set_xlim([0.5,float(nfac)+0.5])
   ax.set_xlabel("Number of members")
   ax.set_ylabel("Mean Square Error")
   ax.plot(x, mse_mean[0,:], 'k-')
   ax.plot(x, mse_mean[0,:]-mse_std[0,:], 'k--')
   ax.plot(x, mse_mean[0,:]+mse_std[0,:], 'k--')
   ax.plot(x, mse_mean[1,:], 'r-')
   ax.plot(x, mse_mean[1,:]-mse_std[1,:], 'r--')
   ax.plot(x, mse_mean[1,:]+mse_std[1,:], 'r--')

   # Save and close figure
   plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + ".jpg", format="jpg", dpi=300)
   plt.close()
