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
nnd = 51
dnd = 1.0/float(nnd-1)
nd = np.linspace(0, (nnd-1)*dnd, nnd)
epsabs_hor = 1.0e-2
epsabs_ver = 1.0e-4
support_th = 1.0e-2

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
support_hor = np.zeros(npk)
cost_hor = np.zeros(npk)
f_sqrt_ver = np.zeros((nnd,npk,nnl))
f_int_ver = np.zeros((nnd,npk,nnl))
valid_ver = np.zeros((npk,nnl), dtype=bool)
scalev = np.zeros((nscaleth,npk,nnl))
scalevdef = np.zeros(nscaleth)
support_ver = np.zeros((npk,nnl))
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

      # Real support radius
      for ind in range(0,nnd):
         if all(np.abs(f_int_hor[ind:,ipk]) < support_th):
            support_hor[ipk] = axis[ind]
            break

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

   print('support_hor: ' + str(support_hor))

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

         # Real support radius
         for ind in range(0,nnd):
            if all(np.abs(f_int_ver[ind:,ipk,inl]) < support_th):
               support_ver[ipk,inl] = axis[ind]
               break

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

      print('support_ver: ' + str(support_ver[:,inl]))

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
   file.write("use tools_const, only: zero,half,one,two,four,eight\n")
   file.write("use tools_kinds, only: kind_real\n")
   file.write("use tools_repro, only: rth,eq,inf,sup\n")
   file.write("use type_mpl, only: mpl_type\n")
   file.write("@:use_probe()\n")
   file.write("\n")
   file.write("implicit none\n")
   file.write("\n")
   file.write("integer,parameter :: nnd = " + str(nnd) + "\n")
   file.write("integer,parameter :: npk = " + str(npk) + "\n")
   file.write("integer,parameter :: nnl = " + str(nnl) + "\n")
   file.write("integer,parameter :: nscaleth = " + str(nscaleth) + "\n")
   file.write("integer,parameter :: ipkhdef = " + str(np.argmin(cost_hor)+1) + "\n")
   file.write("integer,parameter :: ipkvdef = " + str(np.argmin(cost_ver)+1) + "\n")
   file.write("real(kind_real),parameter :: ndmin = %.8f_kind_real\n" % (min(nd)))
   file.write("real(kind_real),parameter :: ndmax = %.8f_kind_real\n" % (max(nd)))
   file.write("real(kind_real),parameter :: pkmin = %.8f_kind_real\n" % (pkmin))
   file.write("real(kind_real),parameter :: pkmax = %.8f_kind_real\n" % (pkmax))
   file.write("real(kind_real),parameter :: nlmin = %.8f_kind_real\n" % (nlmin))
   file.write("real(kind_real),parameter :: nlmax = %.8f_kind_real\n" % (nlmax))
   file.write("real(kind_real),parameter :: dnd = %.8f_kind_real\n" % (dnd))
   file.write("real(kind_real),parameter :: dpk = %.8f_kind_real\n" % (dpk))
   file.write("real(kind_real),parameter :: dnl = %.8f_kind_real\n" % (dnl))
   file.write("real(kind_real),parameter :: scalethmin = %.8f_kind_real\n" % (scalethmin))
   file.write("real(kind_real),parameter :: scalethmax = %.8f_kind_real\n" % (scalethmax))
   file.write("real(kind_real),parameter :: pk(npk) = (/ &\n")
   for ipk in range(0, npk):
      if ipk != npk-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (pk[ipk]) + suffix + "\n")
   file.write("real(kind_real),parameter :: nl(nnl) = (/ &\n")
   for inl in range(0, nnl):
      if inl != nnl-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (nl[inl]) + suffix + "\n")
   file.write("real(kind_real),parameter :: scaleth(nscaleth) = (/ &\n")
   for iscaleth in range(0, nscaleth):
      if iscaleth != nscaleth-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (scaleth[iscaleth]) + suffix + "\n")
   file.write("real(kind_real),parameter :: scaleh(nscaleth,npk) = reshape((/ &\n")
   for ipk in range(0, npk):
      for iscaleth in range(0, nscaleth):
         if ipk != npk-1 or iscaleth != nscaleth-1:
            suffix = ","
         else:
            suffix = "/),"
         file.write(" & %.8f_kind_real" % (scaleh[iscaleth,ipk]) + suffix + " &\n")
   file.write(" & (/nscaleth,npk/))\n")
   file.write("real(kind_real),parameter :: func_hor(nnd,npk) = reshape((/ &\n")
   for ipk in range(0, npk):
      for ind in range(0, nnd):
         if ipk != npk-1 or ind != nnd-1:
            suffix = ","
         else:
            suffix = "/),"
         file.write(" & %.8f_kind_real" % (f_int_hor[ind,ipk]) + suffix + " &\n")
   file.write(" & (/nnd,npk/))\n")
   file.write("real(kind_real),parameter :: support_hor(npk) = (/ &\n")
   for ipk in range(0, npk):
      if ipk != npk-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (support_hor[ipk]) + suffix + "\n")
   file.write("real(kind_real),parameter :: scalev(nscaleth,npk,nnl) = reshape((/ &\n")
   for inl in range(0, nnl):
      for ipk in range(0, npk):
         for iscaleth in range(0, nscaleth):
            if inl != nnl-1 or ipk != npk-1 or iscaleth != nscaleth-1:
               suffix = ","
            else:
               suffix = "/),"
            file.write(" & %.8f_kind_real" % (scalev[iscaleth,ipk,inl]) + suffix + " &\n")
   file.write(" & (/nscaleth,npk,nnl/))\n")
   file.write("real(kind_real),parameter :: func_ver(nnd,npk,nnl) = reshape((/ &\n")
   for inl in range(0, nnl):
      for ipk in range(0, npk):
         for ind in range(0, nnd):
            if inl != nnl-1 or ipk != npk-1 or ind != nnd-1:
               suffix = ","
            else:
               suffix = "/),"
            file.write(" & %.8f_kind_real" % (f_int_ver[ind,ipk,inl]) + suffix + " &\n")
   file.write(" & (/nnd,npk,nnl/))\n")
   file.write("real(kind_real),parameter :: support_ver(npk,nnl) = reshape((/ &\n")
   for inl in range(0, nnl):
      for ipk in range(0, npk):
         if inl != nnl-1 or ipk != npk-1:
            suffix = ","
         else:
            suffix = "/),"
         file.write(" & %.8f_kind_real" % (support_ver[ipk,inl]) + suffix + " &\n")
   file.write(" & (/npk,nnl/))\n")
   file.write("\n")
   file.write("interface fit_func\n")
   file.write("   module procedure gc99_fit_func\n")
   file.write("end interface\n")
   file.write("interface fit_func_sqrt\n")
   file.write("   module procedure gc99_fit_func_sqrt\n")
   file.write("end interface\n")
   file.write("\n")
   file.write("private\n")
   file.write("public :: npk,nnl,nscaleth,scaleth,ipkhdef,ipkvdef,pkmin,pkmax,nlmin,nlmax,scalethmax,pk,nl,scaleh,scalev,support_hor,support_ver\n")
   file.write("public :: fit_func,fit_func_sqrt\n")
   file.write("\n")
   file.write("contains\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Function: gc99_fit_func\n")
   file.write("!> Fit function\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("function gc99_fit_func(mpl,dir,nd,pk,nl) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("character(len=*),intent(in) :: dir  !< Direction\n")
   file.write("real(kind_real),intent(in) :: nd    !< Normalized distance\n")
   file.write("real(kind_real),intent(in) :: pk    !< Peakness\n")
   file.write("real(kind_real),intent(in) :: nl    !< Negative lobe parameter\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Local variables\n")
   file.write("integer :: indm,indp,ipkm,ipkp,inlm,inlp\n")
   file.write("real(kind_real) :: bnd,bpk,bnl,rndm,rndp,rpkm,rpkp,rnlm,rnlp\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_func)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Check bounds\n")
   file.write("if (inf(nd,zero)) call mpl%abort('${subr}$','negative normalized distance')\n")
   file.write("if (inf(pk,pkmin).or.sup(pk,pkmax)) call mpl%abort('${subr}$','peakness out of bounds')\n")
   file.write("if (inf(nl,nlmin).or.sup(nl,nlmax)) call mpl%abort('${subr}$','negative lobe out of bounds')\n")
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
   file.write("   bnl = max(nlmin,min(nl,nlmax))\n")
   file.write("\n")
   file.write("   ! Indices\n")
   file.write("   indm = floor(bnd/dnd)+1\n")
   file.write("   if (indm==nnd) then\n")
   file.write("      indp = indm\n")
   file.write("   else\n")
   file.write("      indp = indm+1\n")
   file.write("   end if\n")
   file.write("   ipkm = floor((bpk-pkmin)/dpk)+1\n")
   file.write("   if (ipkm==npk) then\n")
   file.write("      ipkp = ipkm\n")
   file.write("   else\n")
   file.write("      ipkp = ipkm+1\n")
   file.write("   end if\n")
   file.write("   inlm = floor((bnl-nlmin)/dnl)+1\n")
   file.write("   if (inlm==nnl) then\n")
   file.write("      inlp = inlm\n")
   file.write("   else\n")
   file.write("      inlp = inlm+1\n")
   file.write("   end if\n")
   file.write("\n")
   file.write("   ! Coefficients\n")
   file.write("   if (indm==nnd) then\n")
   file.write("      rndm = one\n")
   file.write("   else\n")
   file.write("      rndm = real(indp-1,kind_real)-bnd/dnd\n")
   file.write("   end if\n")
   file.write("   rndp = (one-rndm)\n")
   file.write("   if (ipkm==npk) then\n")
   file.write("      rpkm = one\n")
   file.write("   else\n")
   file.write("      rpkm = real(ipkp-1,kind_real)-(bpk-pkmin)/dpk\n")
   file.write("   end if\n")
   file.write("   rpkp = (one-rpkm)\n")
   file.write("   if (inlm==nnl) then\n")
   file.write("      rnlm = one\n")
   file.write("   else\n")
   file.write("      rnlm = real(inlp-1,kind_real)-(bnl-nlmin)/dnl\n")
   file.write("   end if\n")
   file.write("   rnlp = (one-rnlm)\n")
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
   file.write("      value = rndm*rpkm*rnlm*func_ver(indm,ipkm,inlm) &\n")
   file.write(" & +rndp*rpkm*rnlm*func_ver(indp,ipkm,inlm) &\n")
   file.write(" & +rndm*rpkp*rnlm*func_ver(indm,ipkp,inlm) &\n")
   file.write(" & +rndp*rpkp*rnlm*func_ver(indp,ipkp,inlm) &\n")
   file.write(" & +rndm*rpkm*rnlp*func_ver(indm,ipkm,inlp) &\n")
   file.write(" & +rndp*rpkm*rnlp*func_ver(indp,ipkm,inlp) &\n")
   file.write(" & +rndm*rpkp*rnlp*func_ver(indm,ipkp,inlp) &\n")
   file.write(" & +rndp*rpkp*rnlp*func_ver(indp,ipkp,inlp)\n")
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
   file.write("function gc99_fit_func_sqrt(mpl,nd,pk,nl) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("real(kind_real),intent(in) :: nd    !< Normalized distance\n")
   file.write("real(kind_real),intent(in) :: pk    !< Peakness\n")
   file.write("real(kind_real),intent(in) :: nl    !< Negative lobe parameter\n")
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
   file.write("if (inf(nd,zero)) call mpl%abort('${subr}$','negative normalized distance')\n")
   file.write("\n")
   file.write("if (eq(nd,zero)) then\n")
   file.write("   ! Origin\n")
   file.write("   value = one\n")
   file.write("elseif (sup(nd,half)) then\n")
   file.write("   ! Out of support\n")
   file.write("   value = zero\n")
   file.write("else\n")
   file.write("   if (pk>zero) then\n")
   file.write("      value = (one-two*nd)/(one+two*nd*(pk+pk**4))\n")
   file.write("   else\n")
   file.write("      value = one-(two*nd)**(one+pk**2)\n")
   file.write("   end if\n")
   file.write("   if (sup(nl,zero)) value = value*(one-(two*nd*nl)**2)*(one-two*nd)\n")
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
