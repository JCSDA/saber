#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import Ngl
import numpy as np
import os

def lct_cor(testdata, test, mpi, omp, suffix, testfig):
   # Processor, index and level to plot
   iproc_plot = 1
   ic1a_plot = 1
   il0_plot = 1

   # Open file
   f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + mpi.zfill(6) + "-" + format(iproc_plot, '06d') + ".nc", "r", format="NETCDF4")

   # Get _FillValue
   _FillValue = f.__dict__["_FillValue"]

   # Get vertical unit
   vunit = f["vunit"][:]

   # Get number of levels
   nl0 = vunit.shape[0]

   for group in f.groups:
      # Get lon/lat/levels
      lon = f.groups[group]["lon"][:,:]
      lat = f.groups[group]["lat"][:,:]
      l0rl0_to_l0 = f.groups[group]["l0rl0_to_l0"][:,:]

      # Contour resources
      cres = Ngl.Resources()
      cres.nglDraw = False
      cres.nglFrame = False
      cres.nglMaximize = True
      cres.mpOutlineOn = False
      cres.mpLandFillColor = -1
      cres.mpProjection = "CylindricalEquidistant"
      cres.mpLimitMode = "LatLon"
      cres.mpMinLatF = 0.0
      cres.mpMaxLatF = 90.0
      cres.mpGridSpacingF = 45.0
      cres.cnFillOn = True
      cres.cnFillMode = "AreaFill"
      cres.cnConstFEnableFill = True
      cres.trGridType = "TriangularMesh"
      cres.cnMonoFillPattern = True
      cres.cnMonoFillColor = False
      cres.lbLabelBarOn = False
      cres.cnInfoLabelOn = False
      cres.cnLineLabelsOn = False
      cres.cnLinesOn = False
      cres.cnNoDataLabelOn = False
      cres.cnMissingValPerimOn = True
      cres.cnMissingValPerimColor = "black"
      cres.cnLevelSelectionMode = "ExplicitLevels"
      cres.sfXArray = lon[ic1a_plot,:]
      cres.sfYArray = lat[ic1a_plot,:]

      # Panel resources
      pnlres = Ngl.Resources()
      pnlres.nglFrame = False
      pnlres.nglPanelLabelBar = True

      for fieldname in ["raw", "fit", "fit_filt"]:
         if fieldname in f.groups[group].variables: 
            # Get field
            field = f.groups[group][fieldname][:,:,:,:]

            # Get dimensions
            nl0r = field.shape[2]
            nc3 = field.shape[3]

            # Open workstation
            wks_type = "png"
            wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + group + "_" + fieldname)
            Ngl.define_colormap(wks, "BlWhRe")

            # Resources
            if np.max(np.abs(field)) > 0.0:
               cres.cnLevels = Ngl.fspan(-np.max(np.abs(field)),np.max(np.abs(field)),21)
            else:
               cres.cnLevels = [-1.0,0.0,1.0]

            # Plots
            plot = []
            for il0r in range(0, nl0r):
               if np.any(field[il0_plot,ic1a_plot,il0r,:] != _FillValue):
                  plot.append(Ngl.contour_map(wks, field[il0_plot,ic1a_plot,il0r,:], cres))

            # Panel
            Ngl.panel(wks, plot, [nl0r,1], pnlres)

            # Advance frame
            Ngl.frame(wks)

            # Delete frame
            Ngl.delete_wks(wks)
