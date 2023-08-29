#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import os

def bump_diag(args, suffix):
   """! Plot script for the "diagnostic" files produced by BUMP"""

   # Make directory
   if not os.path.exists(args.test):
      os.mkdir(args.test)

   # Open file
   f = Dataset(args.testdata + "/" + args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Color levels
   levels = np.linspace(-1.0, 1.0, 101)

   # Fit profiles list
   diag_list = ["coef_ens_l0","coef_sta_l0","rh_l0","pkh_l0","D11_l0","D22_l0","D12_l0","rv_l0","pkv_l0","nlv_l0"]

   # Fit profiles
   for diag in diag_list:
      fig, ax = plt.subplots()
      fig.subplots_adjust(right=0.8)
      cmap = matplotlib.cm.get_cmap('Spectral')
      ax.set_title(diag)
      valid = False

      for group in f.groups:
         if group != "cov_full_vertical":
            # Get vertical coordinate
            vert_coord = f.groups[group]["vert_coord"][:]
            vert_coordmin = np.min(vert_coord)
            vert_coordmax = np.max(vert_coord)
            if vert_coordmin < vert_coordmax:
               ax.set_ylim([vert_coordmin,vert_coordmax])

            # Get number of levels
            nl0 = vert_coord.shape[0]

            # Profiles only
            if nl0 > 1:
               for subgroup in f.groups[group].groups:
                  if (diag in f.groups[group].groups[subgroup].variables):
                     ax.plot(f.groups[group].groups[subgroup][diag][:], vert_coord, label=group + " - " + subgroup)
                     valid = True

      if (valid):
         # Single legend
         handles, labels = ax.get_legend_handles_labels()
         fig.legend(handles, labels, loc='center right')

         # Save and close figure
         plotpath = args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_" + diag + ".jpg"
         plt.savefig(plotpath, format="jpg", dpi=300)
         plt.close()
         print(" -> plot produced: " + plotpath)
      else:
         # Just close the figure
         plt.close()

   # Raw and fitted fields
   for group in f.groups:
      if group != "cov_full_vertical":
         # Get vertical coordinate
         vert_coord = f.groups[group]["vert_coord"][:]

         # Get number of levels
         nl0 = vert_coord.shape[0]

         # Horizontal distance
         disth = f.groups[group]["disth"][:]
         nc3 = len(disth)

         # Angular sector
         as_full = f.groups[group]["as"][:]
         nc4 = len(as_full)
         if nc4 > 1:
            as_full = np.append(as_full, as_full+math.pi)
            as_full = np.append(as_full, as_full[0]+2*math.pi)

            # Cylindrical coordinates
            r, theta = np.meshgrid(disth, as_full)

         for subgroup in f.groups[group].groups:
            # One horizontal plot per level
            for il0 in range(0, nl0):
               # Raw data
               if nc4 == 1:
                  raw_hor = np.zeros((nc3))
                  raw_hor = f.groups[group].groups[subgroup].variables["raw_hor"][il0,0,:]
               else:
                  raw_hor = np.zeros((2*nc4+1,nc3))
                  raw_hor[0:nc4,:] = f.groups[group].groups[subgroup].variables["raw_hor"][il0,:,:]
                  raw_hor[nc4:2*nc4,:] = raw_hor[0:nc4,:]
                  raw_hor[2*nc4,:] = raw_hor[0,:]

               # Fit data
               if nc4 == 1:
                  fit_hor = np.zeros((nc3))
                  if "fit_hor" in f.groups[group].groups[subgroup].variables:
                     fit_hor = f.groups[group].groups[subgroup].variables["fit_hor"][il0,0,:]
               else:
                  fit_hor = np.zeros((2*nc4+1,nc3))
                  if "fit_hor" in f.groups[group].groups[subgroup].variables:
                     fit_hor[0:nc4,:] = f.groups[group].groups[subgroup].variables["fit_hor"][il0,:,:]
                     fit_hor[nc4:2*nc4,:] = fit_hor[0:nc4,:]
                     fit_hor[2*nc4,:] = fit_hor[0,:]

               # Polar plot
               if nc4 == 1:
                  fig, ax = plt.subplots(ncols=2)
               else:
                  fig, ax = plt.subplots(ncols=2, subplot_kw=dict(projection='polar'))
               fig.subplots_adjust(wspace=0.4, right=0.8)
               ax[0].set_title("Raw")
               if nc4 == 1:
                  ax[0].plot(disth, raw_hor)
               else:
                  ax[0].contourf(theta, r, raw_hor, levels=levels, cmap="bwr")
               ax[1].set_title("Fit")
               if nc4 == 1:
                  ax[1].plot(disth, fit_hor)
               else:
                  im = ax[1].contourf(theta, r, fit_hor, levels=levels, cmap="bwr")

               if nc4 > 1:
                  # Colorbar
                  cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                  fig.colorbar(im, cax=cbar_ax)

               # Save and close figure
               plotpath = args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_" + group + "_" + subgroup + str(il0) + ".jpg"
               plt.savefig(plotpath, format="jpg", dpi=300)
               plt.close()
               print(" -> plot produced: " + plotpath)

            # Vertical plot

            # Raw data
            raw_zs = f.groups[group].groups[subgroup].variables["raw_zs"][:,:]

            # Fit data
            fit_zs = np.zeros((nl0, nl0))
            if "fit_zs" in f.groups[group].groups[subgroup].variables:
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
            plotpath = args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_" + group + "_" + subgroup + "_ver.jpg"
            plt.savefig(plotpath, format="jpg", dpi=300)
            plt.close()
            print(" -> plot produced: " + plotpath)
