#!/usr/bin/env python3

import argparse
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("srcdir", help="SABER source directory")
args = parser.parse_args()

# General parameters
nnd = 31
dnd = 1.0/float(nnd-1)
nd = np.linspace(0, (nnd-1)*dnd, nnd)
epsabs = 1.0e-1

# Parameters
newfunc = False
f0 = 1.0e-3
if newfunc:
   pkmin = -6.0
   pkmax = 4.0
   dpk = 0.4
else:
   pkmin = 0.0
   pkmax = 4.0
   dpk = 0.2
npk = int((pkmax-pkmin)/dpk)+1
pk = np.linspace(pkmin, pkmax, npk)

run_horizontal = True
run_vertical = True

# Distance
axis = np.zeros(nnd)
for ind in range(0,nnd):
   axis[ind] = float(ind)/float(nnd-1)

# Functions
def S(r,f0,pk):
   if np.abs(r) <= 0.5:
      if pk > 0.0:
         if newfunc:
            return np.exp(np.log(f0+np.exp(-pk))*np.sqrt(2.0*np.abs(r)))-np.exp(-pk)
         else:
            return (1.0/(1.0+(pk+pk**4)*2.0*np.abs(r))-1.0/(1.0+(pk+pk**4)))/(1.0-1.0/(1.0+(pk+pk**4)))
      else:
         if newfunc:
            return 1.0-np.sqrt(2.0*np.abs(r))**(1+pk**2)
         else:
            return 1.0-(2.0*np.abs(r))
   else:
      return 0.0

def S_hor(x,y,f0,pk):
   r = np.sqrt(x**2+y**2)
   return S(r,f0,pk)

def S_ver(z,f0,pk):
   return S(z,f0,pk)

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
f_sqrt_ver = np.zeros((nnd,npk))
f_int_ver = np.zeros((nnd,npk))
f_gc99 = np.zeros((nnd))
cost_hor = np.zeros((npk))
cost_ver = np.zeros((npk))

# GC99 function
for ind in range(0, nnd):
   f_gc99[ind] = GC99(axis[ind])

if run_horizontal:
   for ipk in range(0, npk):
      for ind in range(0, nnd):
         print("horizontal: " + str(ipk) + " : " + str(ind))

         # Square-root function
         f_sqrt_hor[ind,ipk] = S(axis[ind],f0,pk[ipk])/S(0,f0,pk[ipk])

         # Horizontal integration (2D)
         f = lambda  y, x: S_hor(x,y,f0,pk[ipk])*S_hor(axis[ind]-x,y,f0,pk[ipk])
         fint = integrate.dblquad(f, -0.5, 0.5, lambda x: -0.5, lambda x: 0.5, epsabs = epsabs)
         f_int_hor[ind,ipk] = fint[0]
         if ind == 0:
            norm = f_int_hor[ind,ipk]
         f_int_hor[ind,ipk] = f_int_hor[ind,ipk]/norm

      # Cost
      cost_hor[ipk] = sum((f_int_hor[:,ipk]-f_gc99)**2)

   if True:
      # Plot curves
      fig, ax = plt.subplots(ncols=2)
      ax[0].set_xlim([0,1])
      ax[0].set_ylim([0,1.1])
      ax[0].set_title("Square-root function")
      ax[0].axhline(y=0, color="k")
      ax[0].axvline(x=0, color="k")
      ax[0].plot(axis, f_sqrt_hor)
      ax[1].set_xlim([0,1])
      ax[1].set_ylim([0,1.1])
      ax[1].set_title("Convolution function")
      ax[1].axhline(y=0, color="k")
      ax[1].axvline(x=0, color="k")
      ax[1].plot(axis, f_int_hor)
      plt.savefig("fit_hor.jpg", format="jpg", dpi=300)
      plt.close()

if run_vertical:
   for ipk in range(0, npk):
      for ind in range(0, nnd):  
         print("vertical: " + str(ipk) + " : " + str(ind))

         # Square-root function
         f_sqrt_ver[ind,ipk] = S_ver(axis[ind],f0,pk[ipk])/S_ver(0,f0,pk[ipk])

          # Vertical integration (1D)
         f = lambda  z: S_ver(z,f0,pk[ipk])*S_ver(axis[ind]-z,f0,pk[ipk])
         fint = integrate.quad(f, -0.5, 0.5, epsabs = epsabs)
         f_int_ver[ind,ipk] = fint[0]
         if ind == 0:
            norm = f_int_ver[ind,ipk]
         f_int_ver[ind,ipk] = f_int_ver[ind,ipk]/norm

      # Cost
      cost_ver[ipk] = sum((f_int_ver[:,ipk]-f_gc99)**2)

   if True:
      # Plot curves
      fig, ax = plt.subplots(ncols=2)
      ax[0].set_xlim([0,1])
      ax[0].set_ylim([0,1.1])
      ax[0].set_title("Square-root function")
      ax[0].axhline(y=0, color="k")
      ax[0].axvline(x=0, color="k")
      ax[0].plot(axis, f_sqrt_ver)
      ax[1].set_xlim([0,1])
      ax[1].set_ylim([0,1.1])
      ax[1].set_title("Convolution function")
      ax[1].axhline(y=0, color="k")
      ax[1].axvline(x=0, color="k")
      ax[1].plot(axis, f_int_ver)
      plt.savefig("fit_ver.jpg", format="jpg", dpi=300)
      plt.close()

# Get default peaknesses
pkhdef = pk[np.argmin(cost_hor)]
pkvdef = pk[np.argmin(cost_ver)]

if run_horizontal and run_vertical:
   # Open file
   file = open(args.srcdir + "/src/saber/bump/tools_gc99.fypp", "w")
   
   # Write file 
   file.write("#:include '../instrumentation.fypp'\n")
   file.write("#:include '../generics.fypp'\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Module: tools_gc99\n")
   file.write("!> Gaspari and Cohn (1999)-inspired functions and their square-roots\n")
   file.write("! Author: Benjamin Menetrier\n")
   file.write("! Licensing: this code is distributed under the CeCILL-C license\n")
   file.write("! Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT\n")
   file.write("! WARNING: this module is generated by the python script\n")
   file.write("!            tools/saber_fit_function.py\n")
   file.write("!          to modify this module, update and rerun the python script\n")
   file.write("! 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("module tools_gc99\n")
   file.write("\n")
   file.write("use tools_const, only: zero,half,one,two\n")
   file.write("use tools_kinds, only: kind_real\n")
   file.write("use tools_repro, only: rth,eq,inf,sup\n")
   file.write("use type_mpl, only: mpl_type\n")
   file.write("@:use_probe()\n")
   file.write("\n")
   file.write("implicit none\n")
   file.write("\n")
   file.write("integer,parameter :: nnd = " + str(nnd) + "\n")
   file.write("integer,parameter :: npk = " + str(npk) + "\n")
   if newfunc:
      file.write("real(kind_real),parameter :: f0 = %.8f_kind_real\n" % (f0))
   file.write("real(kind_real),parameter :: ndmin = %.8f_kind_real\n" % (min(nd)))
   file.write("real(kind_real),parameter :: ndmax = %.8f_kind_real\n" % (max(nd)))
   file.write("real(kind_real),parameter :: pkmin = %.8f_kind_real\n" % (pkmin))
   file.write("real(kind_real),parameter :: pkmax = %.8f_kind_real\n" % (pkmax))
   file.write("real(kind_real),parameter :: dnd = %.8f_kind_real\n" % (dnd))
   file.write("real(kind_real),parameter :: dpk = %.8f_kind_real\n" % (dpk))
   file.write("real(kind_real),parameter :: pkhdef = %.8f_kind_real\n" % (pkhdef))
   file.write("real(kind_real),parameter :: pkvdef = %.8f_kind_real\n" % (pkvdef))
   file.write("real(kind_real),parameter :: func_hor(nnd,npk) = reshape((/ &\n")
   for ipk in range(0, npk):
      for ind in range(0, nnd):
         if ipk != npk-1 or ind != nnd-1:
            suffix = ","
         else:
            suffix = "/),"
         file.write(" & %.8f_kind_real" % (f_int_hor[ind,ipk]) + suffix + " &\n")
   file.write(" & (/nnd,npk/))\n")
   file.write("real(kind_real),parameter :: func_ver(nnd,npk) = reshape((/ &\n")
   for ipk in range(0, npk):
      for ind in range(0, nnd):
         if ipk != npk-1 or ind != nnd-1:
            suffix = ","
         else:
            suffix = "/),"
         file.write(" & %.8f_kind_real" % (f_int_ver[ind,ipk]) + suffix + " &\n")
   file.write(" & (/nnd,npk/))\n")
   file.write("\n")
   file.write("interface fit_func\n")
   file.write("   module procedure gc99_fit_func\n")
   file.write("end interface\n")
   file.write("interface fit_func_sqrt\n")
   file.write("   module procedure gc99_fit_func_sqrt\n")
   file.write("end interface\n")
   file.write("\n")
   file.write("public ! TODO: make useless public stuff private\n")
   file.write("\n")
   file.write("contains\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Function: gc99_fit_func\n")
   file.write("!> Fit function\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("function gc99_fit_func(mpl,dir,nd,pk) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("character(len=*),intent(in) :: dir  !< Direction\n")
   file.write("real(kind_real),intent(in) :: nd    !< Normalized distance\n")
   file.write("real(kind_real),intent(in) :: pk    !< Peakness\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Local variables\n")
   file.write("integer :: indm,indp,ipkm,ipkp\n")
   file.write("real(kind_real) :: bnd,bpk,rndm,rndp,rpkm,rpkp\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_func)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Check bounds\n")
   file.write("if (nd<zero) call mpl%abort('${subr}$','negative normalized distance')\n")
   file.write("if ((pk<pkmin).or.(pk>pkmax)) call mpl%abort('${subr}$','peakness out of bounds')\n")
   file.write("\n")
   file.write("if (eq(nd,zero)) then\n")
   file.write("   ! Origin\n")
   file.write("   value = one\n")
   file.write("elseif (sup(nd,one)) then\n")
   file.write("   ! Out of support\n")
   file.write("   value = zero\n")
   file.write("else\n")
   file.write("   ! Bounded values\n")
   file.write("   bnd = max(ndmin,min(nd,ndmax))\n")
   file.write("   bpk = max(pkmin,min(pk,pkmax))\n")
   file.write("\n")
   file.write("   ! Indices\n")
   file.write("   indm = floor(bnd/dnd)+1\n")
   file.write("   indp = indm+1\n")
   file.write("   ipkm = floor((bpk-pkmin)/dpk)+1\n")
   file.write("   ipkp = ipkm+1\n")
   file.write("\n")
   file.write("   ! Coefficients\n")
   file.write("   rndm = real(indp-1,kind_real)-bnd/dnd\n")
   file.write("   rndp = (one-rndm)\n")
   file.write("   rpkm = real(ipkp-1,kind_real)-(bpk-pkmin)/dpk\n")
   file.write("   rpkp = (one-rpkm)\n")
   file.write("\n")
   file.write("   ! Interpolated value\n")
   file.write("   if (dir=='hor') then\n")
   file.write("      ! Horizontal fit function\n")
   file.write("      value = rndm*rpkm*func_hor(indm,ipkm) &\n")
   file.write(" & +rndp*rpkm*func_hor(indp,ipkm) &\n")
   file.write(" & +rndm*rpkp*func_hor(indm,ipkp) &\n")
   file.write(" & +rndp*rpkp*func_hor(indp,ipkp)\n")
   file.write("   elseif (dir=='ver') then\n")
   file.write("      ! Vertical fit function\n")
   file.write("      value = rndm*rpkm*func_ver(indm,ipkm) &\n")
   file.write(" & +rndp*rpkm*func_ver(indp,ipkm) &\n")
   file.write(" & +rndm*rpkp*func_ver(indm,ipkp) &\n")
   file.write(" & +rndp*rpkp*func_ver(indp,ipkp)\n")
   file.write("   else\n")
   file.write("      call mpl%abort('${subr}$','wrong direction: '//dir)\n")
   file.write("   end if\n")
   file.write("end if\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end function gc99_fit_func\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Function: gc99_fit_func_sqrt\n")
   file.write("!> Fit function function square-root\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("function gc99_fit_func_sqrt(mpl,nd,pk) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("real(kind_real),intent(in) :: nd    !< Normalized distance\n")
   file.write("real(kind_real),intent(in) :: pk    !< Peakness\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_func_sqrt)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Check bounds\n")
   file.write("if (nd<zero) call mpl%abort('${subr}$','negative normalized distance')\n")
   file.write("\n")
   file.write("if (eq(nd,zero)) then\n")
   file.write("   ! Origin\n")
   file.write("   value = one\n")
   file.write("elseif (sup(nd,half)) then\n")
   file.write("   ! Out of support\n")
   file.write("   value = zero\n")
   file.write("else\n")
   file.write("   if (pk>zero) then\n")
   if newfunc:
      file.write("      value = exp(log(f0+exp(-pk))*sqrt(two*nd))-exp(-pk)\n")
   else:
      file.write("      value = (one/(one+pk**4*two*nd)-one/(one+pk**4))/(one-one/(one+pk**4))\n")
   file.write("   else\n")
   if newfunc:
      file.write("      value = one-sqrt(two*nd)**(one+pk**2)\n")
   else:
      file.write("      value = one-(two*nd)**(one+pk**2)\n")
   file.write("   end if\n")
   file.write("end if\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end function gc99_fit_func_sqrt\n")
   file.write("\n")
   file.write("end module tools_gc99")
   
   # Close file
   file.close()
