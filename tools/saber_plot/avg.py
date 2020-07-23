#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import Ngl
import numpy as np
import os

def avg(testdata, test, mpi, omp, suffix, testfig):
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

      # XY resources
      xyres = Ngl.Resources()
      xyres.nglFrame = False
      xyres.nglDraw = False
      xyres.xyMarkLineMode = "Markers"
      xyres.xyMarkerColors = 0
      xyres.tiYAxisString = "Distribution"

      # Bar resources
      gsres = Ngl.Resources()
      gsres.gsFillColor = "red"

      # Title resources
      txres = Ngl.Resources()
      txres.txFontHeightF = 0.01
      txres.txJust = "CenterCenter"

      # Panel resources
      pnlres = Ngl.Resources()
      pnlres.nglFrame = False
      pnlres.nglPanelXWhiteSpacePercent = 5.0
      pnlres.nglPanelYWhiteSpacePercent = 5.0

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
                  # Open workstation
                  wks_type = "png"
                  wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + group + "_" + str(il0+1) + "-" + str(l0rl0_to_l0[il0,jl0r]) + "-" + str(jc3+1))

                  # Plots
                  plot = []

                  # Covariance
                  xyres.tiXAxisString = "Covariance"
                  x = 0.5*(m11_bins[il0,jl0r,jc3,0:nbins]+m11_bins[il0,jl0r,jc3,1:nbins+1])
                  p = Ngl.xy(wks, x, m11_hist[il0,jl0r,jc3,:], xyres)
                  for ibin in range(0, nbins):
                     xm = m11_bins[il0,jl0r,jc3,ibin]
                     xp = m11_bins[il0,jl0r,jc3,ibin+1]
                     y = m11_hist[il0,jl0r,jc3,ibin]
                     Ngl.add_polygon(wks,p,[xm,xp,xp,xm,xm],[0.0,0.0,y,y,0.0],gsres)
                     Ngl.add_polyline(wks,p,[xm,xp,xp,xm,xm],[0.0,0.0,y,y,0.0])
                  plot.append(p)

                  # Covariance squared
                  xyres.tiXAxisString = "Covariance squared"
                  x = 0.5*(m11m11_bins[il0,jl0r,jc3,0:nbins]+m11m11_bins[il0,jl0r,jc3,1:nbins+1])
                  p = Ngl.xy(wks, x, m11m11_hist[il0,jl0r,jc3,:], xyres)
                  for ibin in range(0, nbins):
                     xm = m11m11_bins[il0,jl0r,jc3,ibin]
                     xp = m11m11_bins[il0,jl0r,jc3,ibin+1]
                     y = m11m11_hist[il0,jl0r,jc3,ibin]
                     Ngl.add_polygon(wks,p,[xm,xp,xp,xm,xm],[0.0,0.0,y,y,0.0],gsres)
                     Ngl.add_polyline(wks,p,[xm,xp,xp,xm,xm],[0.0,0.0,y,y,0.0])
                  plot.append(p)

                  # Variance product
                  xyres.tiXAxisString = "Variance product"
                  x = 0.5*(m2m2_bins[il0,jl0r,jc3,0:nbins]+m2m2_bins[il0,jl0r,jc3,1:nbins+1])
                  p = Ngl.xy(wks, x, m2m2_hist[il0,jl0r,jc3,:], xyres)
                  for ibin in range(0, nbins):
                     xm = m2m2_bins[il0,jl0r,jc3,ibin]
                     xp = m2m2_bins[il0,jl0r,jc3,ibin+1]
                     y = m2m2_hist[il0,jl0r,jc3,ibin]
                     Ngl.add_polygon(wks,p,[xm,xp,xp,xm,xm],[0.0,0.0,y,y,0.0],gsres)
                     Ngl.add_polyline(wks,p,[xm,xp,xp,xm,xm],[0.0,0.0,y,y,0.0])
                  plot.append(p)

                  # Fourth-order moment
                  xyres.tiXAxisString = "Fourth-order moment"
                  x = 0.5*(m22_bins[il0,jl0r,jc3,0:nbins]+m22_bins[il0,jl0r,jc3,1:nbins+1])
                  p = Ngl.xy(wks, x, m22_hist[il0,jl0r,jc3,:], xyres)
                  for ibin in range(0, nbins):
                     xm = m22_bins[il0,jl0r,jc3,ibin]
                     xp = m22_bins[il0,jl0r,jc3,ibin+1]
                     y = m22_hist[il0,jl0r,jc3,ibin]
                     Ngl.add_polygon(wks,p,[xm,xp,xp,xm,xm],[0.0,0.0,y,y,0.0],gsres)
                     Ngl.add_polyline(wks,p,[xm,xp,xp,xm,xm],[0.0,0.0,y,y,0.0])
                  plot.append(p)

                  # Correlation
                  xyres.tiXAxisString = "Correlation"
                  x = 0.5*(cor_bins[il0,jl0r,jc3,0:nbins]+cor_bins[il0,jl0r,jc3,1:nbins+1])
                  p = Ngl.xy(wks, x, cor_hist[il0,jl0r,jc3,:], xyres)
                  for ibin in range(0, nbins):
                     xm = cor_bins[il0,jl0r,jc3,ibin]
                     xp = cor_bins[il0,jl0r,jc3,ibin+1]
                     y = cor_hist[il0,jl0r,jc3,ibin]
                     Ngl.add_polygon(wks,p,[xm,xp,xp,xm,xm],[0.0,0.0,y,y,0.0],gsres)
                     Ngl.add_polyline(wks,p,[xm,xp,xp,xm,xm],[0.0,0.0,y,y,0.0])
                  plot.append(p)

                  # Add title
                  Ngl.text_ndc(wks, "Between levels " + '%.2e'%vunit[il0] + " and " + '%.2e'%vunit[l0rl0_to_l0[il0,jl0r]-1] + " for distance class " + '%.2e'%disth[jc3], 0.5, 0.61, txres)

                  # Panel
                  Ngl.panel(wks, plot, [1,5], pnlres)

                  # Advance frame
                  Ngl.frame(wks)

                  # Delete frame
                  Ngl.delete_wks(wks)
