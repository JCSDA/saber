#!/usr/bin/env python3

import argparse
import Ngl, Nio
import numpy as np
import numpy.ma as ma
import os

def diag(testdata, test, mpi, omp, suffix):
   # Open file
   f = Nio.open_file(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc")

   # Make output directory
   testfig = testdata + "/" + test + "/fig"
   if not os.path.exists(testfig):
      os.mkdir(testfig)

   # Loop over variables to get roots
   root = []
   for var in f.variables:
      if var.find("cor_")==0 or var.find("loc_")==0:
         # Correlation of localization
         var_root = var.replace("_coef_ens", "").replace("_coef_sta", "").replace("_fit_rh", "").replace("_fit_rv", "").replace("_l0rl0_to_l0", "").replace("_raw", "").replace("_fit", "").replace("_valid", "").replace("_zs", "")

         if not var_root in root:
            root.append(var_root)

   # Get vertical unit
   vunit = f.variables["vunit"][:]

   # Get number of levels
   nl0 = vunit.shape[0]

   # Initialize arrays
   if (len(root)>1):
      colors = Ngl.fspan(2, 189, len(root)).astype(int)
   else:
      colors = np.full(len(root), 2)
   coef_ens = np.zeros((0, nl0))
   fit_rh = np.zeros((0, nl0))
   fit_rv = np.zeros((0, nl0))
   _FillValue = f.variables[root[0] + "_coef_ens"].attributes["_FillValue"]
   empty_vector = np.full((1, nl0), _FillValue)

   # Get data
   for var_root in root:
      coef_ens = np.append(coef_ens, [f.variables[var_root + "_coef_ens"][:]], axis=0)
      if var_root + "_fit_rh" in f.variables:
         fit_rh = np.append(fit_rh, [f.variables[var_root + "_fit_rh"][:]], axis=0)
      else:
         fit_rh = np.append(fit_rh, empty_vector, axis=0)
      if var_root + "_fit_rv" in f.variables:
         fit_rv = np.append(fit_rv, [f.variables[var_root + "_fit_rv"][:]], axis=0)
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
   xleg = np.full(len(root), 0.78)
   if (len(root)>1):
     yleg = Ngl.fspan(0.62, 0.42, len(root))
   else:
     yleg = np.full(len(root), 0.62)
   xtxt = xleg+dx+0.015
   ytxt = yleg
   i = 0
   for var_root in root:
      lnres.gsLineColor = colors[i]
      Ngl.polyline_ndc(wks, [xleg[i], xleg[i]+dx], [yleg[i],yleg[i]], lnres)
      Ngl.text_ndc(wks, var_root, xtxt[i], ytxt[i], txres)
      i = i+1

   # Panel
   Ngl.panel(wks, plot, [1, 3], pnlres)

   # Advance frame
   Ngl.frame(wks)

   # Delete frame
   Ngl.delete_wks(wks)
