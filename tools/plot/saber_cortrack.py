#!/usr/bin/env python3

import argparse
import Ngl,Nio
import os

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("testdata", help="Test data directory")
parser.add_argument("test", help="Test name")
parser.add_argument("fname", help="File name (no extension)")
args = parser.parse_args()

# Open file
f = Nio.open_file(args.testdata + "/" + args.test + "/" + args.fname + ".nc")

# Get lon/lat
lon = f.variables["lon"][:]
lat = f.variables["lat"][:]

# Variables
var_list = ["cor_00","cor_06","tracker_00_0","tracker_00_1","tracker_06_0"]

# Get number of levels
nl0 = f.variables[var_list[0]][:,:].shape[0]

# Plot resources
res = Ngl.Resources()
res.nglDraw = False
res.nglFrame = False
res.nglMaximize = True
res.cnFillOn = True
res.cnFillMode = "AreaFill"
res.trGridType = "TriangularMesh"
res.cnMonoFillPattern = True
res.cnMonoFillColor = False
res.cnFillPalette = "BlWhRe"
res.lbLabelBarOn = False
res.cnInfoLabelOn = False
res.cnLineLabelsOn = False
res.cnLinesOn = False
res.cnNoDataLabelOn = False
res.cnMissingValFillColor = "black"
res.cnMissingValFillPattern = 1
res.cnMissingValPerimOn = True
res.cnMissingValPerimColor = "black"
res.cnLevelSelectionMode = "ManualLevels"
res.cnMaxLevelValF = 1.0
res.cnMinLevelValF = -1.0
res.cnLevelSpacingF = 0.1
res.mpOutlineOn = False
res.mpLandFillColor = -1
res.mpProjection = "CylindricalEquidistant"
res.mpLimitMode = "LatLon"
res.mpMinLatF = 0.0
res.mpMaxLatF = 90.0
res.sfXArray = lon
res.sfYArray = lat

# Panel resources
pnlres = Ngl.Resources()
pnlres.nglFrame = False
pnlres.nglPanelLabelBar = True

# Make output directory
testfig = args.testdata + "/" + args.test + "/fig"
os.mkdir(testfig)

for var in var_list:
   # Read variable
   field = f.variables[var][:,:]

   # Open workstation
   wks_type = "png"
   wks = Ngl.open_wks(wks_type,testfig + "/" + args.fname + "_" + var)

   # Plots
   plot = []
   for il0 in range(0,nl0):
      p = Ngl.contour_map(wks,field[il0,:],res)
      plot.append(p)

   # Panel
   Ngl.panel(wks,plot[0:nl0],[nl0,1],pnlres)

   # Advance frame
   Ngl.frame(wks)

# Close Ngl
Ngl.end()
