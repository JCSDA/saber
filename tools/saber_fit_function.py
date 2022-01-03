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
pkhdef = pk[ipkhdef]
ipkvdef = np.argmin(cost_ver)
pkvdef = pk[ipkvdef]

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
   file.write("use tools_netcdf, only: open_file,inquire_dim_size,get_att,inquire_var,get_var,close_file\n")
   file.write("use tools_repro, only: rth,eq,inf,sup\n")
   file.write("use type_mpl, only: mpl_type\n")
   file.write("@:use_probe()\n")
   file.write("\n")
   file.write("implicit none\n")
   file.write("\n")
   file.write("! Public parameters\n")
   file.write("logical :: fit_allocated = .false.\n")
   file.write("integer,parameter :: nnd = " + str(nnd) + "\n")
   file.write("integer,parameter :: npk = " + str(npk) + "\n")
   file.write("integer,parameter :: nnl = " + str(nnl) + "\n")
   file.write("integer,parameter :: nscaleth = " + str(nscaleth) + "\n")
   file.write("integer,parameter :: ipkhdef = " + str(ipkhdef+1) + "\n")
   file.write("integer,parameter :: ipkvdef = " + str(ipkvdef+1) + "\n")
   file.write("real(kind_real),parameter :: ndmin = %.8f_kind_real\n" % (min(nd)))
   file.write("real(kind_real),parameter :: ndmax = %.8f_kind_real\n" % (max(nd)))
   file.write("real(kind_real),parameter :: pkmin = %.8f_kind_real\n" % (pkmin))
   file.write("real(kind_real),parameter :: pkmax = %.8f_kind_real\n" % (pkmax))
   file.write("real(kind_real),parameter :: pkhdef = %.8f_kind_real\n" % (pkhdef))
   file.write("real(kind_real),parameter :: pkvdef = %.8f_kind_real\n" % (pkvdef))
   file.write("real(kind_real),parameter :: nlmin = %.8f_kind_real\n" % (nlmin))
   file.write("real(kind_real),parameter :: nlmax = %.8f_kind_real\n" % (nlmax))
   file.write("real(kind_real),parameter :: dnd = %.8f_kind_real\n" % (dnd))
   file.write("real(kind_real),parameter :: dpk = %.8f_kind_real\n" % (dpk))
   file.write("real(kind_real),parameter :: dnl = %.8f_kind_real\n" % (dnl))
   file.write("real(kind_real),parameter :: scalethmin = %.8f_kind_real\n" % (scalethmin))
   file.write("real(kind_real),parameter :: scalethmax = %.8f_kind_real\n" % (scalethmax))
   file.write("real(kind_real),allocatable :: pk(:)\n")
   file.write("real(kind_real),allocatable :: nl(:)\n")
   file.write("real(kind_real),allocatable :: scaleth(:)\n")
   file.write("real(kind_real),allocatable :: scaleh(:,:)\n")
   file.write("real(kind_real),allocatable :: func_hor(:,:)\n")
   file.write("real(kind_real),allocatable :: scalev(:,:,:)\n")
   file.write("real(kind_real),allocatable :: func_ver(:,:,:)\n")
   file.write("\n")
   file.write("interface fit_setup\n")
   file.write("   module procedure gc99_fit_setup\n")
   file.write("end interface\n")
   file.write("interface fit_dealloc\n")
   file.write("   module procedure gc99_fit_dealloc\n")
   file.write("end interface\n")
   file.write("interface fit_func\n")
   file.write("   module procedure gc99_fit_func\n")
   file.write("end interface\n")
   file.write("interface fit_func_sqrt\n")
   file.write("   module procedure gc99_fit_func_sqrt\n")
   file.write("end interface\n")
   file.write("\n")
   file.write("private\n")
   file.write("public :: npk,nnl,nscaleth,scaleth,ipkhdef,ipkvdef,pkmin,pkmax,pkhdef,pkvdef,nlmin,nlmax,scalethmin,scalethmax,dpk,dnl\n")
   file.write("public :: pk,nl,scaleh,scalev\n")
   file.write("public :: fit_setup,fit_dealloc,fit_func,fit_func_sqrt\n")
   file.write("\n")
   file.write("contains\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Subroutine: gc99_fit_setup\n")
   file.write("!> Fit setup\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("subroutine gc99_fit_setup(mpl)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("\n")
   file.write("! Local variables\n")
   file.write("integer :: ncid,pk_id,nl_id,scaleth_id,scaleh_id,func_hor_id,scalev_id,func_ver_id\n")
   file.write("character(len=1024) :: filename \n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_setup)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("if (.not.fit_allocated) then\n")
   file.write("   if (mpl%main) then\n")
   file.write("      ! Get file name\n")
   file.write("      filename = '${_FILE_}$.nc'\n")
   file.write("\n")
   file.write("      ! Open file\n")
   file.write("      ncid = open_file(mpl,filename,0,.true.)\n")
   file.write("   end if\n")
   file.write("\n")
   file.write("   ! Allocation\n")
   file.write("   allocate(pk(npk))\n")
   file.write("   allocate(nl(nnl))\n")
   file.write("   allocate(scaleth(nscaleth))\n")
   file.write("   allocate(scaleh(npk,nscaleth))\n")
   file.write("   allocate(func_hor(npk,nnd))\n")
   file.write("   allocate(scalev(nnl,npk,nscaleth))\n")
   file.write("   allocate(func_ver(nnl,npk,nnd))\n")
   file.write("\n")
   file.write("   if (mpl%main) then\n")
   file.write("      ! Inquire variable\n")
   file.write("      pk_id = inquire_var(mpl,ncid,'pk')\n")
   file.write("      nl_id = inquire_var(mpl,ncid,'nl')\n")
   file.write("      scaleth_id = inquire_var(mpl,ncid,'scaleth')\n")
   file.write("      scaleh_id = inquire_var(mpl,ncid,'scaleh')\n")
   file.write("      func_hor_id = inquire_var(mpl,ncid,'func_hor')\n")
   file.write("      scalev_id = inquire_var(mpl,ncid,'scalev')\n")
   file.write("      func_ver_id = inquire_var(mpl,ncid,'func_ver')\n")
   file.write("\n")
   file.write("      ! Read variable\n")
   file.write("      call get_var(mpl,ncid,pk_id,pk)\n")
   file.write("      call get_var(mpl,ncid,nl_id,nl)\n")
   file.write("      call get_var(mpl,ncid,scaleth_id,scaleth)\n")
   file.write("      call get_var(mpl,ncid,scaleh_id,scaleh)\n")
   file.write("      call get_var(mpl,ncid,func_hor_id,func_hor)\n")
   file.write("      call get_var(mpl,ncid,scalev_id,scalev)\n")
   file.write("      call get_var(mpl,ncid,func_ver_id,func_ver)\n")
   file.write("\n")
   file.write("      ! Close file\n")
   file.write("      call close_file(mpl,ncid)\n")
   file.write("   end if\n")
   file.write("\n")
   file.write("   ! Broadcast variables\n")
   file.write("   call mpl%f_comm%broadcast(pk,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(nl,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(scaleth,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(scaleh,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(func_hor,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(scalev,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(func_ver,mpl%rootproc-1)\n")
   file.write("\n")
   file.write("   ! Set flag\n")
   file.write("   fit_allocated = .true.\n")
   file.write("end if\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end subroutine gc99_fit_setup\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Subroutine: gc99_fit_dealloc\n")
   file.write("!> Fit setup\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("subroutine gc99_fit_dealloc()\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_dealloc)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Release memory\n")
   file.write("if (allocated(pk)) deallocate(pk)\n")
   file.write("if (allocated(nl)) deallocate(nl)\n")
   file.write("if (allocated(scaleth)) deallocate(scaleth)\n")
   file.write("if (allocated(scaleh)) deallocate(scaleh)\n")
   file.write("if (allocated(func_hor)) deallocate(func_hor)\n")
   file.write("if (allocated(scalev)) deallocate(scalev)\n")
   file.write("if (allocated(func_ver)) deallocate(func_ver)\n")
   file.write("\n")
   file.write("! Reset flag\n")
   file.write("fit_allocated = .false.\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end subroutine gc99_fit_dealloc\n")
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
   file.write("      value = rndm*rpkm*func_hor(ipkm,indm) &\n")
   file.write(" & +rndp*rpkm*func_hor(ipkm,indp) &\n")
   file.write(" & +rndm*rpkp*func_hor(ipkp,indm) &\n")
   file.write(" & +rndp*rpkp*func_hor(ipkp,indp)\n")
   file.write("   elseif (dir=='ver') then\n")
   file.write("      ! Vertical fit function\n")
   file.write("      value = rndm*rpkm*rnlm*func_ver(inlm,ipkm,indm) &\n")
   file.write(" & +rndp*rpkm*rnlm*func_ver(inlm,ipkm,indp) &\n")
   file.write(" & +rndm*rpkp*rnlm*func_ver(inlm,ipkp,indm) &\n")
   file.write(" & +rndp*rpkp*rnlm*func_ver(inlm,ipkp,indp) &\n")
   file.write(" & +rndm*rpkm*rnlp*func_ver(inlp,ipkm,indm) &\n")
   file.write(" & +rndp*rpkm*rnlp*func_ver(inlp,ipkm,indp) &\n")
   file.write(" & +rndm*rpkp*rnlp*func_ver(inlp,ipkp,indm) &\n")
   file.write(" & +rndp*rpkp*rnlp*func_ver(inlp,ipkp,indp)\n")
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

   # Create NetCDF file
   ncfile = Dataset(args.srcdir + "/src/saber/bump/tools_gc99.fypp.nc",mode="w",format='NETCDF4_CLASSIC') 

   # Create dimensions
   nnd_id = ncfile.createDimension('nnd', nnd)
   npk_id = ncfile.createDimension('npk', npk)
   nnl_id = ncfile.createDimension('nnl', nnl)
   nscaleth_id = ncfile.createDimension('nscaleth', nscaleth)

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
