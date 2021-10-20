#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import os

def diag(testdata, test, mpi, omp, suffix, testfig):
   """! Plot script for the "diagnostic" files produced by BUMP"""

   # Open file
   f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Get _FillValue
   _FillValue = f.__dict__["_FillValue"]

   # Color levels
   levels = np.linspace(-1.0, 1.0, 101)

   # Fit profiles list
   diag_list = ["coef_ens","fit_rh","fit_D11","fit_D22","fit_D12","fit_pkh","fit_rv","fit_pkv","fit_nlv"]

   # Fit profiles
   for diag in diag_list:
      fig, ax = plt.subplots()
      fig.subplots_adjust(right=0.8)
      cmap = matplotlib.cm.get_cmap('Spectral')
      ax.set_title(diag)
      valid = False

      for group in f.groups:
         if group != "cov_full_vertical":
            # Get vertical unit
            vunit = f.groups[group]["vunit"][:]
            vunitmin = np.min(vunit)
            vunitmax = np.max(vunit)
            if vunitmin < vunitmax:
               ax.set_ylim([vunitmin,vunitmax])

            # Get number of levels
            nl0 = vunit.shape[0]

            # Profiles only
            if nl0 > 1:
               for subgroup in f.groups[group].groups:
                  if (diag in f.groups[group].groups[subgroup].variables):
                     ax.plot(f.groups[group].groups[subgroup][diag][:], vunit, label=group + " - " + subgroup)
                     valid = True

      if (valid):
         # Single legend
         handles, labels = ax.get_legend_handles_labels()
         fig.legend(handles, labels, loc='center right')

         # Save and close figure
         plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + diag + ".jpg", format="jpg", dpi=300)
         plt.close()
      else:
         # Just close the figure
         plt.close()

   # Raw and fitted fields
   for group in f.groups:
      if group != "cov_full_vertical":
         # Get vertical unit
         vunit = f.groups[group]["vunit"][:]

         # Get number of levels
         nl0 = vunit.shape[0]

         # Horizontal distance
         disth = f.groups[group]["disth"][:]
         nc3 = len(disth)

         # Angular sector
         as_full = f.groups[group]["as"][:]
         nc4 = len(as_full)
         as_full = np.append(as_full, as_full+math.pi)
         as_full = np.append(as_full, as_full[0]+2*math.pi)
         print(as_full)

         # Cylindrical coordinates
         r, theta = np.meshgrid(disth, as_full)

         for subgroup in f.groups[group].groups:
            # One horizontal plot per level
            for il0 in range(0, nl0):
               # Raw data
               raw_hor = np.zeros((2*nc4+1,nc3))
               raw_hor[0:nc4,:] = f.groups[group].groups[subgroup].variables["raw_hor"][il0,:,:]
               raw_hor[nc4:2*nc4,:] = raw_hor[0:nc4,:]
               raw_hor[2*nc4,:] = raw_hor[0,:]

               # Fit data
               fit_hor = np.zeros((2*nc4+1,nc3))
               fit_hor[0:nc4,:] = f.groups[group].groups[subgroup].variables["fit_hor"][il0,:,:]
               fit_hor[nc4:2*nc4,:] = fit_hor[0:nc4,:]
               fit_hor[2*nc4,:] = fit_hor[0,:]

               # Polar plot
               fig, ax = plt.subplots(ncols=2, subplot_kw=dict(projection='polar'))
               fig.subplots_adjust(wspace=0.4, right=0.8)
               ax[0].set_title("Raw")
               ax[0].contourf(theta, r, raw_hor, levels=levels, cmap="bwr")
               ax[1].set_title("Fit")
               im = ax[1].contourf(theta, r, fit_hor, levels=levels, cmap="bwr")

               # Colorbar
               cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
               fig.colorbar(im, cax=cbar_ax)

               # Save and close figure
               plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + group + "_" + subgroup + "_" + str(il0) + ".jpg", format="jpg", dpi=300)
               plt.close()

            # Vertical plot

            # Raw data
            raw_zs = f.groups[group].groups[subgroup].variables["raw_zs"][:,:]

            # Fit data
            fit_zs = f.groups[group].groups[subgroup].variables["fit_zs"][:,:]

            # Plot
            fig, ax = plt.subplots(ncols=2)
            fig.subplots_adjust(wspace=0.4, right=0.8)
            ax[0].set_title("Raw")
            ax[0].matshow(raw_zs, vmin=-1, vmax=1, cmap="bwr")
            ax[1].set_title("Fit")
            im = ax[1].matshow(fit_zs, vmin=-1, vmax=1, cmap="bwr")

            # Colorbar
            cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            fig.colorbar(im, cax=cbar_ax)

            # Save and close figure
            plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + group + "_" + subgroup + "_ver.jpg", format="jpg", dpi=300)
            plt.close()
