#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import Ngl
import numpy as np
import numpy.ma as ma
import os

def diag(testdata, test, mpi, omp, suffix, testfig):
   # Open file
   f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Get _FillValue
   _FillValue = f.__dict__["_FillValue"]

   # Get vertical unit
   vunit = f["vunit"][:]

   # Get number of levels
   nl0 = vunit.shape[0]

   # Profiles only
   if nl0 == 1:
      return

   # Number of subgroups
   ngr = 0
   for group in f.groups:
      for subgroup in f.groups[group].groups:
         ngr += 1

   # Initialize arrays
   if ngr > 1:
      colors = Ngl.fspan(2, 189, ngr).astype(int)
   else:
      colors = np.full(ngr, 2)
   coef_ens = np.zeros((0, nl0))
   fit_rh = np.zeros((0, nl0))
   fit_rv = np.zeros((0, nl0))
   empty_vector = np.full((1, nl0), _FillValue)

   # Get data
   for group in f.groups:
      for subgroup in f.groups[group].groups:
         coef_ens = np.append(coef_ens, [f.groups[group].groups[subgroup]["coef_ens"][:]], axis=0)
         if "fit_rh" in f.groups[group].groups[subgroup].variables:
            fit_rh = np.append(fit_rh, [f.groups[group].groups[subgroup]["fit_rh"][:]], axis=0)
         else:
           fit_rh = np.append(fit_rh, empty_vector, axis=0)
         if "fit_rv" in f.groups[group].groups[subgroup].variables:
            fit_rv = np.append(fit_rv, [f.groups[group].groups[subgroup]["fit_rv"][:]], axis=0)
         else:
            fit_rv = np.append(fit_rv, empty_vector, axis=0)

   # XY resources
   xyres = Ngl.Resources()
   xyres.nglFrame = False
   xyres.nglDraw = False
   xyres.xyMarkLineMode = "Lines"
   xyres.xyLineColors = colors
   xyres.xyLineThicknessF = 4.0
   xyres.trYMinF = min(vunit)
   xyres.trYMaxF = max(vunit)

   # Legend resources
   lnres = Ngl.Resources()
   lnres.gsLineThicknessF = 4.0

   # Text resources
   txres = Ngl.Resources()
   txres.txFontHeightF = 0.01
   txres.txJust = "CenterLeft"

   # Panel resources
   pnlres = Ngl.Resources()
   pnlres.nglFrame = False
   pnlres.nglPanelXWhiteSpacePercent = 5.0
   pnlres.nglPanelRight  = 0.75

   # Open workstation
   wks_type = "png"
   wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_profiles")

   # Plots
   plot = []

   # coef_ens
   xyres.tiXAxisString = "coef_ens"
   plot.append(Ngl.xy(wks, ma.masked_values(coef_ens, _FillValue), vunit, xyres))

   # fit_rh
   xyres.tiXAxisString = "fit_rh"
   plot.append(Ngl.xy(wks, ma.masked_values(fit_rh, _FillValue), vunit, xyres))

   # fit_rv
   xyres.tiXAxisString = "fit_rv"
   plot.append(Ngl.xy(wks, ma.masked_values(fit_rv, _FillValue), vunit, xyres))

   # Legend
   dx = 0.04
   xleg = np.full(ngr, 0.78)
   if ngr > 1:
     yleg = Ngl.fspan(0.62, 0.42, ngr)
   else:
     yleg = np.full(ngr, 0.62)
   xtxt = xleg+dx+0.015
   ytxt = yleg
   i = 0
   for group in f.groups:
      for subgroup in f.groups[group].groups:
         lnres.gsLineColor = colors[i]
         Ngl.polyline_ndc(wks, [xleg[i], xleg[i]+dx], [yleg[i],yleg[i]], lnres)
         Ngl.text_ndc(wks, group + " / " + subgroup, xtxt[i], ytxt[i], txres)
         i = i+1

   # Panel
   Ngl.panel(wks, plot, [1, 3], pnlres)

   # Advance frame
   Ngl.frame(wks)

   # Delete frame
   Ngl.delete_wks(wks)
