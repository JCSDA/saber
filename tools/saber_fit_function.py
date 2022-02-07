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
nncomp = 7
scalethmin = 0.2
scalethmax = 0.9
dscaleth = 0.1
nscaleth = int((scalethmax-scalethmin)/dscaleth+1.0e-6)+1
scaleth = np.linspace(scalethmin, scalethmax, nscaleth)
dsf = 0.01
nsf = int(1.0/dsf+1.0e-6)
run_horizontal = True
run_vertical = True

# Distance
axis = np.zeros(nnd)
for ind in range(0,nnd):
   axis[ind] = float(ind)/float(nnd-1)

# Functions
def S(r):
   if np.abs(r) <= 0.5:
      return 1.0-(2.0*np.abs(r))
   else:
      return 0.0

def S_hor(x,y,):
   r = np.sqrt(x**2+y**2)
   return S(r)

def S_ver(z):
   return S(z)

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
f_sqrt_hor = np.zeros((nnd))
f_sqrt_hor_comp = np.zeros((nnd))
f_int_hor = np.zeros((nnd))
f_int_hor_comp = np.zeros((nnd))
scaleh = np.zeros((nncomp,nscaleth))
f_sqrt_ver = np.zeros((nnd))
f_sqrt_ver_comp = np.zeros((nnd))
f_int_ver = np.zeros((nnd))
f_int_ver_comp = np.zeros((nnd))
scalev = np.zeros((nncomp,nscaleth))
scaled_axis = np.zeros((nnd))
support_factor = np.zeros((nncomp))
f_gc99 = np.zeros((nnd))

if run_horizontal:
   for ind in range(0, nnd):
      # Square-root function
      f_sqrt_hor[ind] = S(axis[ind])

      # Horizontal integration (2D)
      f = lambda  y, x: S_hor(x,y)*S_hor(axis[ind]-x,y)
      fint = integrate.dblquad(f, -0.5, 0.5, lambda x: -0.5, lambda x: 0.5, epsabs = epsabs_hor)
      f_int_hor[ind] = fint[0]
      if ind == 0:
         norm = f_int_hor[ind]
      f_int_hor[ind] = f_int_hor[ind]/norm

   # Loop over nncomp
   for incomp in range(0, nncomp):
      # Composed functions
      f_sqrt_hor_comp_detail = np.zeros((nnd,incomp+1))
      for ind in range(0, nnd-1):
         f_sqrt_hor_comp[ind] = 0.0
         f_sqrt_hor_comp_detail[ind,:] = 0.0
         f_int_hor_comp[ind] = 0.0
         for icomp in range(0, incomp+1):
            if (2**icomp*ind < nnd):
               f_sqrt_hor_comp[ind] = f_sqrt_hor_comp[ind]+f_sqrt_hor[2**icomp*ind]
               for jcomp in range(0, incomp+1):
                  if (jcomp >= icomp):
                     f_sqrt_hor_comp_detail[ind,jcomp] = f_sqrt_hor_comp_detail[ind,jcomp]+f_sqrt_hor[2**icomp*ind]
               f_int_hor_comp[ind] = f_int_hor_comp[ind]+f_int_hor[2**icomp*ind]
         f_sqrt_hor_comp[ind] = f_sqrt_hor_comp[ind]/(incomp+1)
         f_sqrt_hor_comp_detail[ind,:] = f_sqrt_hor_comp_detail[ind,:]/(incomp+1)
         f_int_hor_comp[ind] = f_int_hor_comp[ind]/(incomp+1)

      # Scale at scaleth
      for iscaleth in range(0, nscaleth):
         scaleh[incomp,iscaleth] = 1.0
         for ind in range(0, nnd-1):
            if f_int_hor_comp[ind]>scaleth[iscaleth] and f_int_hor_comp[ind+1]<scaleth[iscaleth]:
               A = (f_int_hor_comp[ind]-f_int_hor_comp[ind+1])/(axis[ind]-axis[ind+1])
               B = f_int_hor_comp[ind]-A*axis[ind]
               scaleh[incomp,iscaleth] = (scaleth[iscaleth]-B)/A
               break

      # Support factor for the best GC99 fit
      err_min = 1.0e12
      support_factor[incomp] = -1.0
      for isf in range(0, nsf):
         sf = (isf+1)*dsf
         err = 0.0
         for ind in range(0, nnd-1):
            err = err+(f_int_hor_comp[ind]-GC99(axis[ind]/sf))**2
         if (err < err_min):
            err_min = err
            support_factor[incomp] = sf
      for ind in range(0, nnd):
         f_gc99[ind] = GC99(axis[ind]/support_factor[incomp])
      print("For " + str(incomp+1) + " component(s), support factor is: " + str(support_factor[incomp]) + " and smallest component is " + str(1.0/2**icomp))

      if True:
         # Plot curves
         fig, ax = plt.subplots(ncols=2, figsize=(14,7))
         ax[0].set_xlim([0,1.0])
         ax[0].set_ylim([0,1.1])
         ax[0].set_title("Square-root function")
         ax[0].axhline(y=0, color="k")
         ax[0].axvline(x=0, color="k")
         ax[0].plot(axis, f_sqrt_hor_comp_detail)
         ax[0].plot(axis, f_sqrt_hor_comp, 'k')
         ax[1].set_xlim([0,1.0])
         ax[1].set_ylim([0,1.1])
         ax[1].set_title("Convolution function")
         ax[1].axhline(y=0, color="k")
         ax[1].axvline(x=0, color="k")
         ax[1].plot(axis, f_int_hor_comp, 'b', axis, f_gc99, 'r')
         plt.savefig("fit_hor_" + str(incomp+1) + ".jpg", format="jpg", dpi=300)
         plt.close()

if run_vertical:
   for ind in range(0, nnd):
      # Square-root function
      f_sqrt_ver[ind] = S_ver(axis[ind])/S_ver(0)

      # Vertical integration (1D)
      f = lambda  z: S_ver(z)*S_ver(axis[ind]-z)
      fint = integrate.quad(f, -0.5, 0.5, epsabs = epsabs_ver)
      f_int_ver[ind] = fint[0]
      if ind == 0:
         norm = f_int_ver[ind]
      f_int_ver[ind] = f_int_ver[ind]/norm

   # Loop over nncomp
   for incomp in range(0, nncomp):
      # Composed functions
      f_sqrt_ver_comp_detail = np.zeros((nnd,incomp+1))
      for ind in range(0, nnd-1):
         f_sqrt_ver_comp[ind] = 0.0
         f_sqrt_ver_comp_detail[ind,:] = 0.0
         f_int_ver_comp[ind] = 0.0
         for icomp in range(0, incomp+1):
            if (2**icomp*ind < nnd):
               f_sqrt_ver_comp[ind] = f_sqrt_ver_comp[ind]+f_sqrt_ver[2**icomp*ind]
               for jcomp in range(0, incomp+1):
                  if (jcomp >= icomp):
                     f_sqrt_ver_comp_detail[ind,jcomp] = f_sqrt_ver_comp_detail[ind,jcomp]+f_sqrt_ver[2**icomp*ind]
               f_int_ver_comp[ind] = f_int_ver_comp[ind]+f_int_ver[2**icomp*ind]
         f_sqrt_ver_comp[ind] = f_sqrt_ver_comp[ind]/(incomp+1)
         f_sqrt_ver_comp_detail[ind,:] = f_sqrt_ver_comp_detail[ind,:]/(incomp+1)
         f_int_ver_comp[ind] = f_int_ver_comp[ind]/(incomp+1)

      # Scale at scaleth
      for iscaleth in range(0, nscaleth):
         scalev[incomp,iscaleth] = 1.0
         for ind in range(0, nnd-1):
            if f_int_ver_comp[ind]>scaleth[iscaleth] and f_int_ver_comp[ind+1]<scaleth[iscaleth]:
               A = (f_int_ver_comp[ind]-f_int_ver_comp[ind+1])/(axis[ind]-axis[ind+1])
               B = f_int_ver_comp[ind]-A*axis[ind]
               scalev[incomp,iscaleth] = (scaleth[iscaleth]-B)/A
               break

      # Support factor for the best GC99 fit
      err_min = 1.0e12
      support_factor[incomp] = -1.0
      for isf in range(0, nsf):
         sf = (isf+1)*dsf
         err = 0.0
         for ind in range(0, nnd-1):
            err = err+(f_int_ver_comp[ind]-GC99(axis[ind]/sf))**2
         if (err < err_min):
            err_min = err
            support_factor[incomp] = sf
      for ind in range(0, nnd):
         f_gc99[ind] = GC99(axis[ind]/support_factor[incomp])
      print("For " + str(incomp+1) + " component(s), support factor is: " + str(support_factor[incomp]) + " and smallest component is " + str(1.0/2**icomp))

      if True:
         # Plot curves
         fig, ax = plt.subplots(ncols=2, figsize=(14,7))
         ax[0].set_xlim([0,1.0])
         ax[0].set_ylim([0,1.1])
         ax[0].set_title("Square-root function")
         ax[0].axhline(y=0, color="k")
         ax[0].axvline(x=0, color="k")
         ax[0].plot(axis, f_sqrt_ver_comp_detail)
         ax[0].plot(axis, f_sqrt_ver_comp, 'k')
         ax[1].set_xlim([0,1.0])
         ax[1].set_ylim([0,1.1])
         ax[1].set_title("Convolution function")
         ax[1].axhline(y=0, color="k")
         ax[1].axvline(x=0, color="k")
         ax[1].plot(axis, f_int_ver_comp, 'b', axis, f_gc99, 'r')
         plt.savefig("fit_ver_" + str(incomp+1) + ".jpg", format="jpg", dpi=300)
         plt.close()

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
   file.write("use tools_repro, only: rth,eq,inf,infeq,sup\n")
   file.write("use type_mpl, only: mpl_type\n")
   file.write("@:use_probe()\n")
   file.write("\n")
   file.write("implicit none\n")
   file.write("\n")
   file.write("! Public parameters\n")
   file.write("logical :: fit_allocated = .false.\n")
   file.write("integer,parameter :: nnd = " + str(nnd) + "\n")
   file.write("integer,parameter :: nncomp = " + str(nncomp) + "\n")
   file.write("integer,parameter :: nscaleth = " + str(nscaleth) + "\n")
   file.write("real(kind_real),parameter :: ndmin = %.8f_kind_real\n" % (min(nd)))
   file.write("real(kind_real),parameter :: ndmax = %.8f_kind_real\n" % (max(nd)))
   file.write("real(kind_real),parameter :: dnd = %.8f_kind_real\n" % (dnd))
   file.write("real(kind_real),parameter :: scalethmin = %.8f_kind_real\n" % (scalethmin))
   file.write("real(kind_real),parameter :: scalethmax = %.8f_kind_real\n" % (scalethmax))
   file.write("real(kind_real),allocatable :: scaleth(:)\n")
   file.write("real(kind_real),allocatable :: scaleh(:,:)\n")
   file.write("real(kind_real),allocatable :: func_hor(:)\n")
   file.write("real(kind_real),allocatable :: scalev(:,:)\n")
   file.write("real(kind_real),allocatable :: func_ver(:)\n")
   file.write("real(kind_real),allocatable :: support_factor(:)\n")
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
   file.write("public :: nncomp,nscaleth,scaleth,scalethmin,scalethmax\n")
   file.write("public :: scaleh,scalev,support_factor\n")
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
   file.write("integer :: ncid,scaleth_id,scaleh_id,func_hor_id,scalev_id,func_ver_id,support_factor_id\n")
   file.write("character(len=1024) :: filename\n")
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
   file.write("   allocate(scaleth(nscaleth))\n")
   file.write("   allocate(scaleh(nscaleth,nncomp))\n")
   file.write("   allocate(func_hor(nnd))\n")
   file.write("   allocate(scalev(nscaleth,nncomp))\n")
   file.write("   allocate(func_ver(nnd))\n")
   file.write("   allocate(support_factor(nncomp))\n")
   file.write("\n")
   file.write("   if (mpl%main) then\n")
   file.write("      ! Inquire variable\n")
   file.write("      scaleth_id = inquire_var(mpl,ncid,'scaleth')\n")
   file.write("      scaleh_id = inquire_var(mpl,ncid,'scaleh')\n")
   file.write("      func_hor_id = inquire_var(mpl,ncid,'func_hor')\n")
   file.write("      scalev_id = inquire_var(mpl,ncid,'scalev')\n")
   file.write("      func_ver_id = inquire_var(mpl,ncid,'func_ver')\n")
   file.write("      support_factor_id = inquire_var(mpl,ncid,'support_factor')\n")
   file.write("\n")
   file.write("      ! Read variable\n")
   file.write("      call get_var(mpl,ncid,scaleth_id,scaleth)\n")
   file.write("      call get_var(mpl,ncid,scaleh_id,scaleh)\n")
   file.write("      call get_var(mpl,ncid,func_hor_id,func_hor)\n")
   file.write("      call get_var(mpl,ncid,scalev_id,scalev)\n")
   file.write("      call get_var(mpl,ncid,func_ver_id,func_ver)\n")
   file.write("      call get_var(mpl,ncid,support_factor_id,support_factor)\n")
   file.write("\n")
   file.write("      ! Close file\n")
   file.write("      call close_file(mpl,ncid)\n")
   file.write("   end if\n")
   file.write("\n")
   file.write("   ! Broadcast variables\n")
   file.write("   call mpl%f_comm%broadcast(scaleth,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(scaleh,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(func_hor,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(scalev,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(func_ver,mpl%rootproc-1)\n")
   file.write("   call mpl%f_comm%broadcast(support_factor,mpl%rootproc-1)\n")
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
   file.write("if (allocated(scaleth)) deallocate(scaleth)\n")
   file.write("if (allocated(scaleh)) deallocate(scaleh)\n")
   file.write("if (allocated(func_hor)) deallocate(func_hor)\n")
   file.write("if (allocated(scalev)) deallocate(scalev)\n")
   file.write("if (allocated(func_ver)) deallocate(func_ver)\n")
   file.write("if (allocated(support_factor)) deallocate(support_factor)\n")
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
   file.write("function gc99_fit_func(mpl,dir,nd,ncomp) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("character(len=*),intent(in) :: dir  !< Direction\n")
   file.write("real(kind_real),intent(in) :: nd    !< Normalized distance\n")
   file.write("integer,intent(in) :: ncomp         !< Number of components\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Local variables\n")
   file.write("integer :: icomp,indm,indp\n")
   file.write("real(kind_real) :: lnd,bnd,rndm,rndp\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_func)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Check bounds\n")
   file.write("if (inf(nd,zero)) call mpl%abort('${subr}$','negative normalized distance')\n")
   file.write("\n")
   file.write("! Initialization\n")
   file.write("value = zero\n")
   file.write("\n")
   file.write("do icomp=1,ncomp\n")
   file.write("   ! Local normalized distance\n")
   file.write("   lnd = real(2**(icomp-1),kind_real)*nd\n")
   file.write("\n")
   file.write("   if (eq(lnd,zero)) then\n")
   file.write("      ! Origin\n")
   file.write("      value = value+one\n")
   file.write("   elseif (infeq(lnd,one)) then\n")
   file.write("      ! Indices\n")
   file.write("      indm = floor(lnd/dnd)+1\n")
   file.write("      if (indm==nnd) then\n")
   file.write("         indp = indm\n")
   file.write("      else\n")
   file.write("         indp = indm+1\n")
   file.write("      end if\n")
   file.write("\n")
   file.write("      ! Coefficients\n")
   file.write("      if (indm==nnd) then\n")
   file.write("         rndm = one\n")
   file.write("      else\n")
   file.write("         rndm = real(indp-1,kind_real)-lnd/dnd\n")
   file.write("      end if\n")
   file.write("      rndp = (one-rndm)\n")
   file.write("\n")
   file.write("      ! Interpolated value\n")
   file.write("      if (dir=='hor') then\n")
   file.write("         ! Horizontal fit function\n")
   file.write("         value = value+rndm*func_hor(indm)+rndp*func_hor(indp)\n")
   file.write("      elseif (dir=='ver') then\n")
   file.write("         ! Vertical fit function\n")
   file.write("         value = value+rndm*func_ver(indm)+rndp*func_ver(indp)\n")
   file.write("      else\n")
   file.write("         call mpl%abort('${subr}$','wrong direction: '//dir)\n")
   file.write("      end if\n")
   file.write("   end if\n")
   file.write("end do\n")
   file.write("\n")
   file.write("! Normalization\n")
   file.write("value = value/real(ncomp,kind_real)\n")
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
   file.write("function gc99_fit_func_sqrt(mpl,nd,ncomp) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("real(kind_real),intent(in) :: nd    !< Normalized distance\n")
   file.write("integer,intent(in) :: ncomp         !< Number of components\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Local variables\n")
   file.write("integer :: icomp\n")
   file.write("real(kind_real) :: lnd\n")
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
   file.write("! Initialization\n")
   file.write("value = zero\n")
   file.write("\n")
   file.write("do icomp=1,ncomp\n")
   file.write("   ! Local normalized distance\n")
   file.write("   lnd = real(2**(icomp-1),kind_real)*nd\n")
   file.write("\n")
   file.write("   if (eq(lnd,zero)) then\n")
   file.write("      ! Origin\n")
   file.write("      value = value+one\n")
   file.write("   elseif (infeq(lnd,half)) then\n")
   file.write("      ! Out of support\n")
   file.write("      value = value+one-(two*lnd)\n")
   file.write("   end if\n")
   file.write("end do\n")
   file.write("\n")
   file.write("! Normalization\n")
   file.write("value = value/real(ncomp,kind_real)\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end function gc99_fit_func_sqrt\n")
   file.write("\n")
   file.write("end module tools_gc99\n")

   # Close file
   file.close()

   # Create NetCDF file
   ncfile = Dataset(args.srcdir + "/src/saber/bump/tools_gc99.fypp.nc",mode="w",format='NETCDF4_CLASSIC')

   # Create dimensions
   nnd_id = ncfile.createDimension('nnd', nnd)
   nncomp_id = ncfile.createDimension('nncomp', nncomp)
   nscaleth_id = ncfile.createDimension('nscaleth', nscaleth)

   # Create variables
   scaleth_id = ncfile.createVariable('scaleth', np.float64, ('nscaleth'))
   scaleh_id = ncfile.createVariable('scaleh', np.float64, ('nncomp','nscaleth'))
   func_hor_id = ncfile.createVariable('func_hor', np.float64, ('nnd'))
   scalev_id = ncfile.createVariable('scalev', np.float64, ('nncomp','nscaleth'))
   func_ver_id = ncfile.createVariable('func_ver', np.float64, ('nnd'))
   support_factor_id = ncfile.createVariable('support_factor', np.float64, ('nncomp'))

   # Write variables
   scaleth_id[:] = scaleth
   scaleh_id[:,:] = scaleh
   func_hor_id[:] = f_int_hor
   scalev_id[:,:] = scalev
   func_ver_id[:] = f_int_ver
   support_factor_id[:] = support_factor

   # Close file
   ncfile.close()
