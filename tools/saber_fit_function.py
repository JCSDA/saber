#!/usr/bin/env python3

import argparse
import os
import math
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import scipy.integrate as integrate

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("srcdir", help="SABER source directory")
args = parser.parse_args()

# General parameters
nnd = 51
dnd = 1.0/float(nnd-1)
nd = np.linspace(0, (nnd-1)*dnd, nnd)
epsabs_hor = 1.0e-2
epsabs_ver = 1.0e-4

# Parameters
scalethmin = 0.2
scalethmax = 0.8
dscaleth = 0.1
nscaleth = int((scalethmax-scalethmin)/dscaleth)+1
scaleth = np.linspace(scalethmin, scalethmax, nscaleth)

pkmin = -2.0
pkmax = 4.0
dpk = 0.2
npk = int((pkmax-pkmin)/dpk)+1
pk = np.linspace(pkmin, pkmax, npk)

nlmin = 0.0
nlmax = 6.0
dnl = 1.0
nnl = int((nlmax-nlmin)/dnl)+1
nl = np.linspace(nlmin, nlmax, nnl)

run_horizontal = True
run_vertical = True

# Distance
axis = np.zeros(nnd)
for ind in range(0,nnd):
   axis[ind] = float(ind)/float(nnd-1)

# Functions
def S(r,pk):
   if np.abs(r) <= 0.5:
      if pk > 0.0:
         return (1.0-2.0*np.abs(r))/(1.0+2.0*np.abs(r)*(pk+pk**4))
      else:
         return 1.0-(2.0*np.abs(r))**(1.0+pk**2)
   else:
      return 0.0

def S_hor(x,y,pk):
   r = np.sqrt(x**2+y**2)
   return S(r,pk)

def S_ver(z,pk,nl):
   return S(z,pk)*(1.0-(2.0*np.abs(z)*nl)**2)*(1.0-2.0*np.abs(z))

def GC99(r):
   if r<0.5:
      value = 1.0-r
      value = 1.0+8.0/5.0*r*value
      value = 1.0-3.0/4.0*r*value
      value = 1.0-20.0/3.0*r**2*value
   else:
      if r<1.0:
         value = 1.0-r/3.0
         value = 1.0-8.0/5.0*r*value
         value = 1.0+3.0/4.0*r*value
         value = 1.0-2.0/3.0*r*value
         value = 1.0-5.0/2.0*r*value
         value = 1.0-12.0*r*value
         value = -value/(3.0*r)
      else:
        value = 0.0
   return value

# Initialize arrays
f_sqrt_hor = np.zeros((nnd,npk))
f_int_hor = np.zeros((nnd,npk))
scaleh = np.zeros((nscaleth,npk))
scalehdef = np.zeros(nscaleth)
cost_hor = np.zeros(npk)
f_sqrt_ver = np.zeros((nnd,npk,nnl))
f_int_ver = np.zeros((nnd,npk,nnl))
valid_ver = np.zeros((npk,nnl), dtype=bool)
scalev = np.zeros((nscaleth,npk,nnl))
scalevdef = np.zeros(nscaleth)
cost_ver = np.zeros(npk)
f_sqrt_ver_plot = np.zeros((nnd,npk,nnl))
f_int_ver_plot = np.zeros((nnd,npk,nnl))
f_gc99 = np.zeros(nnd)
scaled_axis = np.zeros((nnd,npk))

# GC99 function
for ind in range(0, nnd):
   f_gc99[ind] = GC99(axis[ind])

if run_horizontal:
   for ipk in range(0, npk):
      for ind in range(0, nnd):
         print("horizontal: " + str(ipk) + " : " + str(ind))

         # Square-root function
         f_sqrt_hor[ind,ipk] = S(axis[ind],pk[ipk])

         # Horizontal integration (2D)
         f = lambda  y, x: S_hor(x,y,pk[ipk])*S_hor(axis[ind]-x,y,pk[ipk])
         fint = integrate.dblquad(f, -0.5, 0.5, lambda x: -0.5, lambda x: 0.5, epsabs = epsabs_hor)
         f_int_hor[ind,ipk] = fint[0]
         if ind == 0:
            norm = f_int_hor[ind,ipk]
         f_int_hor[ind,ipk] = f_int_hor[ind,ipk]/norm

      # Cost
      cost_hor[ipk] = sum((f_int_hor[:,ipk]-f_gc99)**2)

      # Scale at scaleth
      for iscaleth in range(0, nscaleth):
         scaleh[iscaleth,ipk] = 1.0
         for ind in range(0, nnd-1):
            if f_int_hor[ind,ipk]>scaleth[iscaleth] and f_int_hor[ind+1,ipk]<scaleth[iscaleth]:
               A = (f_int_hor[ind,ipk]-f_int_hor[ind+1,ipk])/(axis[ind]-axis[ind+1])
               B = f_int_hor[ind,ipk]-A*axis[ind]
               scaleh[iscaleth,ipk] = (scaleth[iscaleth]-B)/A
               break

   if True:
      # Plot curves
      fig, ax = plt.subplots(ncols=2, figsize=(14,7))
      ax[0].set_xlim([0,1.0])
      ax[0].set_ylim([0,1.1])
      ax[0].set_title("Square-root function")
      ax[0].axhline(y=0, color="k")
      ax[0].axvline(x=0, color="k")
      ax[0].plot(axis, f_sqrt_hor)
      ax[1].set_xlim([0,1.0])
      ax[1].set_ylim([0,1.1])
      ax[1].set_title("Convolution function")
      ax[1].axhline(y=0, color="k")
      ax[1].axvline(x=0, color="k")
      ax[1].plot(axis, f_int_hor)
      plt.savefig("fit_hor.jpg", format="jpg", dpi=300)
      plt.close()

      for iscaleth in range(0, nscaleth):
         for ipk in range(0, npk):
            scaled_axis[:,ipk] = axis/scaleh[iscaleth,ipk]
         fig, ax = plt.subplots()
         ax.set_xlim([0,1.0/min(scaleh[iscaleth,:])])
         ax.set_ylim([0,1.1])
         ax.set_title("Scaled convolution function: " + str(scaleth[iscaleth]))
         ax.axhline(y=0, color="k")
         ax.axvline(x=0, color="k")
         ax.plot(scaled_axis, f_int_hor)
         plt.savefig("fit_hor_" + str(iscaleth) + ".jpg", format="jpg", dpi=300)
         plt.close()

if run_vertical:
   for inl in range(0, nnl):
      for ipk in range(0, npk):
         for ind in range(0, nnd):  
            print("vertical: " + str(inl) + " / " + str(ipk) + " : " + str(ind))

            # Square-root function
            f_sqrt_ver[ind,ipk,inl] = S_ver(axis[ind],pk[ipk],nl[inl])/S_ver(0,pk[ipk],nl[inl])

            # Vertical integration (1D)
            f = lambda  z: S_ver(z,pk[ipk],nl[inl])*S_ver(axis[ind]-z,pk[ipk],nl[inl])
            fint = integrate.quad(f, -0.5, 0.5, epsabs = epsabs_ver)
            f_int_ver[ind,ipk,inl] = fint[0]
            if ind == 0:
               norm = f_int_ver[ind,ipk,inl]
            f_int_ver[ind,ipk,inl] = f_int_ver[ind,ipk,inl]/norm

         # Check conditions
         if inl == 0:
            # No envelop
            valid_ver[ipk,inl] = True
         else:
            # Should have significant negative values
            imin = np.argmin(f_int_ver[:,ipk,inl])
            if f_int_ver[imin,ipk,inl]<0.0:
               # Should have small second maximum
               if max(f_int_ver[imin:,ipk,inl])<0.1:
                  valid_ver[ipk,inl] = True

         # Cost
         if inl == 0:
            cost_ver[ipk] = sum((f_int_ver[:,ipk,inl]-f_gc99)**2)

         # Scale at scaleth
         for iscaleth in range(0, nscaleth):
            scalev[iscaleth,ipk,inl] = 1.0
            for ind in range(0, nnd-1):
               if f_int_ver[ind,ipk,inl]>scaleth[iscaleth] and f_int_ver[ind+1,ipk,inl]<scaleth[iscaleth]:
                  A = (f_int_ver[ind,ipk,inl]-f_int_ver[ind+1,ipk,inl])/(axis[ind]-axis[ind+1])
                  B = f_int_ver[ind,ipk,inl]-A*axis[ind]
                  scalev[iscaleth,ipk,inl] = (scaleth[iscaleth]-B)/A
                  break

      if True:
         for ipk in range(0, npk):
            # Remove invalid curves
            if valid_ver[ipk,inl]:
               f_sqrt_ver_plot[:,ipk,inl] = f_sqrt_ver[:,ipk,inl]
               f_int_ver_plot[:,ipk,inl] = f_int_ver[:,ipk,inl]
            else:
               f_sqrt_ver_plot[:,ipk,inl] = np.nan
               f_int_ver_plot[:,ipk,inl] = np.nan

         # Plot curves
         fig, ax = plt.subplots(ncols=2, figsize=(14,7))
         ax[0].set_xlim([0,1.0])
         ax[0].set_ylim([-0.5,1.1])
         ax[0].set_title("Square-root function")
         ax[0].axhline(y=0, color="k")
         ax[0].axvline(x=0, color="k")
         ax[0].plot(axis, f_sqrt_ver_plot[:,:,inl])
         ax[1].set_xlim([0,1.0])
         ax[1].set_ylim([-0.5,1.1])
         ax[1].set_title("Convolution function")
         ax[1].axhline(y=0, color="k")
         ax[1].axvline(x=0, color="k")
         ax[1].plot(axis, f_int_ver_plot[:,:,inl])
         plt.savefig("fit_ver_" + str(inl) + ".jpg", format="jpg", dpi=300)
         plt.close()

         if False:
            for iscaleth in range(0, nscaleth):
               for ipk in range(0, npk):
                  scaled_axis[:,ipk] = axis/scalev[iscaleth,ipk,inl]
               fig, ax = plt.subplots()
               ax.set_xlim([0,1.0/min(scalev[iscaleth,:,inl])])
               ax.set_ylim([-0.5,1.1])
               ax.set_title("Scaled convolution function: " + str(scaleth[iscaleth]))
               ax.axhline(y=0, color="k")
               ax.axvline(x=0, color="k")
               ax.plot(scaled_axis, f_int_ver_plot[:,:,inl])
               plt.savefig("fit_ver_" + str(inl) + "_" + str(iscaleth) + ".jpg", format="jpg", dpi=300)
               plt.close()

# Get default peaknesses
ipkhdef = np.argmin(cost_hor)
ipkvdef = np.argmin(cost_ver)

if run_horizontal and run_vertical:
   # Create file
   ncfile = Dataset(args.srcdir + "/src/saber/bump/tools_gc99.fypp.nc",mode="w",format='NETCDF4_CLASSIC') 

   # Create dimensions
   nnd_id = ncfile.createDimension('nnd', nnd)
   npk_id = ncfile.createDimension('npk', npk)
   nnl_id = ncfile.createDimension('nnl', nnl)
   nscaleth_id = ncfile.createDimension('nscaleth', nscaleth)

   # Create attributes
   ncfile.ipkhdef = np.argmin(cost_hor)+1
   ncfile.ipkvdef = np.argmin(cost_ver)+1
   ncfile.ndmin = min(nd)
   ncfile.ndmax = max(nd)
   ncfile.pkmin = pkmin
   ncfile.pkmax = pkmax
   ncfile.nlmin = nlmin
   ncfile.nlmax = nlmax
   ncfile.dnd = dnd
   ncfile.dpk = dpk
   ncfile.dnl = dnl
   ncfile.scalethmin = scalethmin
   ncfile.scalethmax = scalethmax

   # Create variables
   pk_id = ncfile.createVariable('pk', np.float64, ('npk'))
   nl_id = ncfile.createVariable('nl', np.float64, ('nnl'))
   scaleth_id = ncfile.createVariable('scaleth', np.float64, ('nscaleth'))
   scaleh_id = ncfile.createVariable('scaleh', np.float64, ('nscaleth', 'npk'))
   func_hor_id = ncfile.createVariable('func_hor', np.float64, ('nnd', 'npk'))
   scalev_id = ncfile.createVariable('scalev', np.float64, ('nscaleth', 'npk', 'nnl'))
   func_ver_id = ncfile.createVariable('func_ver', np.float64, ('nnd', 'npk', 'nnl'))

   # Write variables
   pk_id[:] = pk
   nl_id[:] = nl
   scaleth_id[:] = scaleth
   scaleh_id[:,:] = scaleh
   func_hor_id[:,:] = f_int_hor
   scalev_id[:,:,:] = scalev
   func_ver_id[:,:,:] = f_int_ver

   # Close file
   ncfile.close()
