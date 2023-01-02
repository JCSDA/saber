#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import os

def bump_normality(args, suffix):
   """! Plot script for the "normality" files produced by BUMP"""

   # Make directory
   if not os.path.exists(args.test):
      os.mkdir(args.test)

   # Loop over files
   first = True
   for impi in range(0, int(args.mpi)):
      # Open file
      f = Dataset(args.testdata + "/" + args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_local_" + args.mpi.zfill(6) + "-" + format(impi+1, '06d') + ".nc", "r", format="NETCDF4")

      # Get data
      if "ens_norm" in f.variables:
         nloc = f["ens_norm"][:].shape[0]
         if first:
            ne = f["ens_norm"][:].shape[1]
            ens_norm = np.zeros((0, ne))
            ens_step = np.zeros((0, ne-1))
            first = False
         for iloc in range(0, nloc):
            ens_norm = np.append(ens_norm, [f["ens_norm"][iloc,:].data], axis=0)
            ens_step = np.append(ens_step, [f["ens_step"][iloc,:].data], axis=0)

   # X axis
   x_norm = np.linspace(1, ne, ne)
   x_step = np.linspace(1, ne-1, ne-1)

   # Plots
   fig, ax = plt.subplots(nrows=2)
   fig.subplots_adjust(hspace=0.5)
   ax[0].set_title("Normalized perturbation amplitude")
   ax[0].set_xlim([1.0,float(ne)])
   ax[0].set_ylim([-0.1,1.1])
   ax[0].set_xlabel("Ordered member")
   ax[0].plot(x_norm, ens_norm[0,:])
   ax[1].set_title("Normalized perturbation amplitude step")
   ax[1].set_xlim([1.0,float(ne-1)])
   ax[1].set_ylim([-0.1,1.1])
   ax[1].set_xlabel("Ordered member")
   ax[1].plot(x_step, ens_step[0,:])

   # Save and close figure
   plotpath = args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + ".jpg"
   plt.savefig(plotpath, format="jpg", dpi=300)
   plt.close()
   print(" -> plot produced: " + plotpath)
