#!/usr/bin/env python3

import argparse
import Ngl, Nio
import numpy as np
import os

def normality(testdata, test, mpi, omp, suffix):
   # Loop over files
   first = True
   for impi in range(0, int(mpi)):
      # Open file
      f = Nio.open_file(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + mpi.zfill(6) + "-" + format(impi+1, '06d') + ".nc")

      # Get data
      if "ens_norm" in f.variables:
         nloc = f.variables["ens_norm"][:].shape[0]
         if first:
            ne = f.variables["ens_norm"][:].shape[1]
            ens_norm = np.zeros((0, ne))
            ens_step = np.zeros((0, ne-1))
            first = False
         for iloc in range(0, nloc):
            ens_norm = np.append(ens_norm, [f.variables["ens_norm"][iloc,:].data], axis=0)
            ens_step = np.append(ens_step, [f.variables["ens_step"][iloc,:].data], axis=0)

   # X axis
   x = Ngl.fspan(1, ne, ne)

   # XY resources
   xyres = Ngl.Resources()
   xyres.nglFrame = False
   xyres.nglDraw = False
   xyres.xyMarkLineMode = "Lines"
   xyres.xyLineColor = "gray"
   xyres.xyLineThicknessF = 4.0
   xyres.tiXAxisString = "Ordered member"
   xyres.trXMinF = 1.0
   xyres.trYMinF = 0.0
   xyres.trYMaxF = 1.0
   xyres.tmYLMode = "Explicit"
   xyres.tmYLValues = Ngl.fspan(0.0, 1.0, 6)
   xyres.tmYLLabels = xyres.tmYLValues
   xyres.tmYLMinorOn = True
   xyres.tmYLMinorValues = Ngl.fspan(-0.1, 1.1, 25)

   # Panel resources
   pnlres = Ngl.Resources()
   pnlres.nglFrame = False
   pnlres.nglPanelYWhiteSpacePercent = 5.0

   # Make output directory
   testfig = testdata + "/" + test + "/fig"
   if not os.path.exists(testfig):
      os.mkdir(testfig)

   # Open workstation
   wks_type = "png"
   wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix)

   # Specific resources
   xyres.tiYAxisString = "Normalized perturbation amplitude"
   xyres.trXMaxF = float(ne)

   # Plot
   plot = []
   plot.append(Ngl.xy(wks, x, ens_norm, xyres))

   # Specific resources
   xyres.tiYAxisString = "Normalized perturbation amplitude step"
   xyres.trXMaxF = float(ne-1)

   # Plot
   plot.append(Ngl.xy(wks, x[0:ne-1], ens_step, xyres))

   # Panel
   Ngl.panel(wks, plot, [1,2], pnlres)

   # Advance frame
   Ngl.frame(wks)

   # Delete frame
   Ngl.delete_wks(wks)
