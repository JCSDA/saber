#!/usr/bin/env python3

import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import scipy.integrate as integrate

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("srcdir", help="SABER source directory")
args = parser.parse_args()

# General parameters
nch = 31
dch = 1.0/float(nch-1)
ch = np.linspace(0, (nch-1)*dch, nch)
ncv = 31
dcv = 1.0/float(ncv-1)
cv = np.linspace(0, (ncv-1)*dcv, ncv)

epsabs = 1.0e-1

# Parameters
f0 = 1.0e-3
npkh = 21
dpkh = 0.25
pkh = np.linspace(0, (npkh-1)*dpkh, npkh)
npkv = 21
dpkv = 0.25
pkv = np.linspace(0, (npkv-1)*dpkv, npkv)
tau = 0.5
nnld = 11
dnld = 0.02
nld = np.linspace(0, (nnld-1)*dnld, nnld)
nnlr = 14
dnlr = 0.1
nlr = np.linspace(0, (nnlr-1)*dnlr, nnlr)+np.ones(nnlr)
nnldt = 50
dnldt = -0.01
nldt = np.linspace(dnldt, nnldt*dnldt, nnldt)
run_horizontal = True
run_vertical = True

def S_hor(x,y,f0,pkh):
   # Horizontal value
   r = np.sqrt(x**2+y**2)
   if r <= 0.5:
      if pkh > 0.0:
         value = np.exp(np.log(f0+np.exp(-pkh))*np.sqrt(2.0*r))-np.exp(-pkh)
      else:
         value = 1.0-2.0*r
   else:
      value = 0.0
   return value

def S_ver(z,f0,pkv,tau,nld):
   # Vertical value
   value = 0.0
   if np.abs(z) <= 0.5*(1.0-tau):
      if pkv > 0:
         value += (1.0+nld)*((np.exp(-pkv*np.sqrt(2.0*np.abs(z)/(1.0-tau)))-np.exp(-pkv))/(1.0-np.exp(-pkv)))
      else:
         value += (1.0+nld)*(1.0-2.0*np.abs(z)/(1.0-tau))
   if np.abs(z) <= 0.5:
      value -= nld*(1.0-2.0*np.abs(z))
   return value

# Distances
haxis = np.zeros(nch)
for ic in range(0,nch):
   haxis[ic] = float(ic)/float(nch-1)
vaxis = np.zeros(ncv)
for ic in range(0,ncv):
   vaxis[ic] = float(ic)/float(ncv-1)

# Initialize arrays
f_sqrt_hor = np.zeros((nch,npkh))
f_int_hor = np.zeros((nch,npkh))
f_sqrt_ver = np.zeros((ncv,npkv,nnld))
f_int_ver = np.zeros((ncv,npkv,nnld))
ver_minval = np.ones((npkv,nnld))
ver_minval[:,:] = 1000.0
ver_ratio = np.zeros((npkv,nnld))
ver_ratio[:,:] = 1000.0
to_inld = -np.ones((nnlr,nnldt), dtype=np.int8)
to_ipkv = -np.ones((nnlr,nnldt), dtype=np.int8)
f_sqrt_ver_plot = np.zeros((ncv,nnlr,nnldt))
f_int_ver_plot = np.zeros((ncv,nnlr,nnldt))

if run_horizontal:
   for ipkh in range(0, npkh):
      for ic in range(0, nch):
          print("horizontal: " + str(ipkh) + " : " + str(ic))

          # Square-root function
          f_sqrt_hor[ic,ipkh] = S_hor(haxis[ic],0,f0,pkh[ipkh])/S_hor(0,0,f0,pkh[ipkh])

          # Horizontal integration (2D)
          f = lambda  y, x: S_hor(x,y,f0,pkh[ipkh])*S_hor(haxis[ic]-x,y,f0,pkh[ipkh])
          fint = integrate.dblquad(f, -0.5, 0.5, lambda x: -0.5, lambda x: 0.5, epsabs = epsabs)
          f_int_hor[ic,ipkh] = fint[0]
          if ic == 0:
             norm = f_int_hor[ic,ipkh]
          f_int_hor[ic,ipkh] = f_int_hor[ic,ipkh]/norm

   # Plot curves
   fig, ax = plt.subplots(ncols=2)
   ax[0].set_xlim([0,1])
   ax[0].set_ylim([0,1.1])
   ax[0].set_title("Square-root function")
   ax[0].axhline(y=0, color="k")
   ax[0].axvline(x=0, color="k")
   ax[0].plot(haxis, f_sqrt_hor)
   ax[1].set_xlim([0,1])
   ax[1].set_ylim([0,1.1])
   ax[1].set_title("Convolution function")
   ax[1].axhline(y=0, color="k")
   ax[1].axvline(x=0, color="k")
   ax[1].plot(haxis, f_int_hor)
   plt.savefig("fit_hor.jpg", format="jpg", dpi=300)
   plt.close()

if run_vertical:
   for ipkv in range(0, npkv):
      for inld in range(0, nnld):
         for icv in range(0, ncv):
             print("vertical:   " + str(ipkv) + " / " + str(inld) + " : " + str(icv))
   
             # Square-root function
             f_sqrt_ver[icv,ipkv,inld] = S_ver(vaxis[icv],f0,pkv[ipkv],tau,nld[inld])/S_ver(0,f0,pkv[ipkv],tau,nld[inld])
   
             # Vertical integration (1D)
             f = lambda  z: S_ver(z,f0,pkv[ipkv],tau,nld[inld])*S_ver(vaxis[icv]-z,f0,pkv[ipkv],tau,nld[inld])
             fint = integrate.quad(f, -0.5, 0.5, epsabs = epsabs)
             f_int_ver[icv,ipkv,inld] = fint[0]
             if icv == 0:
                norm = f_int_ver[icv,ipkv,inld]
             f_int_ver[icv,ipkv,inld] = f_int_ver[icv,ipkv,inld]/norm
   
         if np.any(f_int_ver[:,ipkv,inld] < 0):
            # Find minimum
            ver_minval[ipkv,inld] = min(f_int_ver[:,ipkv,inld])
            ver_minarg = np.argmin(f_int_ver[:,ipkv,inld])

            # Find following maximum
            ver_maxval = max(f_int_ver[ver_minarg:,ipkv,inld])
            if ver_maxval > 0.05:
               f_sqrt_ver[:,ipkv,inld] = np.nan
               f_int_ver[:,ipkv,inld] = np.nan
            else:
               # Find where the function crosses the X axis
               for icv in range(0, ncv-1):
                  if f_int_ver[icv,ipkv,inld] > 0 and f_int_ver[icv+1,ipkv,inld] < 0:
                     ver_zero = (f_int_ver[icv+1,ipkv,inld]*vaxis[icv]-f_int_ver[icv,ipkv,inld]*vaxis[icv+1])/(f_int_ver[icv+1,ipkv,inld]-f_int_ver[icv,ipkv,inld])
                     ver_ratio[ipkv,inld] = vaxis[ver_minarg]/ver_zero
                     break

         else:
            f_sqrt_ver[:,ipkv,inld] = np.nan
            f_int_ver[:,ipkv,inld] = np.nan
 
      # Plot curves
      fig, ax = plt.subplots(ncols=2)
      ax[0].set_xlim([0,1])
      ax[0].set_ylim([-0.5,1.1])
      ax[0].set_title("Square-root function")
      ax[0].axhline(y=0, color="k")
      ax[0].axvline(x=0, color="k")
      ax[0].plot(vaxis, f_sqrt_ver[:,ipkv,:])
      ax[1].set_xlim([0,1])
      ax[1].set_ylim([-0.5,1.1])
      ax[1].set_title("Convolution function")
      ax[1].axhline(y=0, color="k")
      ax[1].axvline(x=0, color="k")
      ax[1].plot(vaxis, f_int_ver[:,ipkv,:])
      plt.savefig("fit_ver_" + str(ipkv) + ".jpg", format="jpg", dpi=300)
      plt.close()

   # Transform arrays
   for inlr in range(0, nnlr):
      for inldt in range(0, nnldt):
          cost = (ver_minval[:,:]-nldt[inldt])**2+(ver_ratio[:,:]-nlr[inlr])**2
          tmp = np.unravel_index(np.argmin(cost, axis=None), cost.shape)
          if abs(ver_minval[tmp[0],tmp[1]]-nldt[inldt]) <= 0.5*abs(dnldt) and abs(ver_ratio[tmp[0],tmp[1]]-nlr[inlr]) <= 0.5*dnlr:
             to_ipkv[inlr,inldt] = tmp[0]
             to_inld[inlr,inldt] = tmp[1]

   # Valid cases
   nvfc = 0
   for inlr in range(0, nnlr):
      for inldt in range(0, nnldt):
         ipkv = to_ipkv[inlr,inldt]
         inld = to_inld[inlr,inldt]
         if inld > 0 and ipkv > 0:
            nvfc += 1
            f_sqrt_ver_plot[:,inlr,inldt] = f_sqrt_ver[:,ipkv,inld]
            f_int_ver_plot[:,inlr,inldt] = f_int_ver[:,ipkv,inld]
         else:
            f_sqrt_ver_plot[:,inlr,inldt] = np.nan
            f_int_ver_plot[:,inlr,inldt] = np.nan

   for inldt in range(0, nnldt):
      # Plot curves
      fig, ax = plt.subplots(ncols=2)
      ax[0].set_xlim([0,1])
      ax[0].set_ylim([-0.5,1.1])
      ax[0].set_title("Square-root function")
      ax[0].axhline(y=0, color="k")
      ax[0].axvline(x=0, color="k")
      ax[0].plot(vaxis, f_sqrt_ver_plot[:,:,inldt])
      ax[1].set_xlim([0,1])
      ax[1].set_ylim([-0.5,1.1])
      ax[1].set_title("Convolution function")
      ax[1].axhline(y=0, color="k")
      ax[1].axvline(x=0, color="k")
      ax[1].plot(vaxis, f_int_ver_plot[:,:,inldt])
      plt.savefig("fit_ver_trans_" + str(inldt) + ".jpg", format="jpg", dpi=300)
      plt.close()

   # Final reformatting into a vertical functions catalog
   vfc_to_pkv = np.zeros((nvfc))
   vfc_to_nld = np.zeros((nvfc))
   vfc_minval = np.zeros((nvfc))
   vfc_ratio = np.zeros((nvfc))
   f_int_vfc = np.zeros((ncv,nvfc))
   ivfc = 0
   for inlr in range(0, nnlr):
      for inldt in range(0, nnldt):
         ipkv = to_ipkv[inlr,inldt]
         inld = to_inld[inlr,inldt]
         if inld > 0 and ipkv > 0:
            vfc_to_pkv[ivfc] = pkv[ipkv]
            vfc_to_nld[ivfc] = nld[inld]
            vfc_minval[ivfc] = ver_minval[ipkv,inld]
            vfc_ratio[ivfc] = ver_ratio[ipkv,inld]
            f_int_vfc[:,ivfc] = f_int_ver[:,ipkv,inld]
            ivfc += 1

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
   file.write("integer,parameter :: nch = " + str(nch) + "\n")
   file.write("integer,parameter :: ncv = " + str(ncv) + "\n")
   file.write("integer,parameter :: npkh = " + str(npkh) + "\n")
   file.write("integer,parameter :: nvfc = " + str(nvfc) + "\n")
   file.write("real(kind_real),parameter :: f0 = %.8f_kind_real\n" % (f0))
   file.write("real(kind_real),parameter :: tau = %.8f_kind_real\n" % (tau))
   file.write("real(kind_real),parameter :: chmin = %.8f_kind_real\n" % (min(ch)))
   file.write("real(kind_real),parameter :: chmax = %.8f_kind_real\n" % (max(ch)))
   file.write("real(kind_real),parameter :: pkhmin = %.8f_kind_real\n" % (min(pkh)))
   file.write("real(kind_real),parameter :: pkhmax = %.8f_kind_real\n" % (max(pkh)))
   file.write("real(kind_real),parameter :: cvmin = %.8f_kind_real\n" % (min(cv)))
   file.write("real(kind_real),parameter :: cvmax = %.8f_kind_real\n" % (max(cv)))
   file.write("real(kind_real),parameter :: pkvmin = %.8f_kind_real\n" % (min(pkv)))
   file.write("real(kind_real),parameter :: pkvmax = %.8f_kind_real\n" % (max(pkv)))
   file.write("real(kind_real),parameter :: dch = %.8f_kind_real\n" % (dch))
   file.write("real(kind_real),parameter :: dpkh = %.8f_kind_real\n" % (dpkh))
   file.write("real(kind_real),parameter :: dcv = %.8f_kind_real\n" % (dcv))
   file.write("real(kind_real),parameter :: dpkv = %.8f_kind_real\n" % (dpkv))
   file.write("real(kind_real),parameter :: ch(nch) = (/ &\n")
   for ich in range(0, nch):
      if ich != nch-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (ch[ich]) + suffix + "\n")
   file.write("real(kind_real),parameter :: func_hor(nch,npkh) = reshape((/ &\n")
   for ipkh in range(0, npkh):
      for ic in range(0, nch):
         if ipkh != npkh-1 or ic != nch-1:
            suffix = ","
         else:
            suffix = "/),"
         file.write(" & %.8f_kind_real" % (f_int_hor[ic,ipkh]) + suffix + " &\n")
   file.write(" & (/ncv,npkh/))\n")
   file.write("real(kind_real),parameter :: cv(ncv) = (/ &\n")
   for icv in range(0, ncv):
      if icv != ncv-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (cv[icv]) + suffix + "\n")
   file.write("real(kind_real),parameter :: vfc_to_pkv(nvfc) = (/ &\n")
   for ivfc in range(0, nvfc):
      if ivfc != nvfc-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (vfc_to_pkv[ivfc]) + suffix + "\n")
   file.write("real(kind_real),parameter :: vfc_to_nld(nvfc) = (/ &\n")
   for ivfc in range(0, nvfc):
      if ivfc != nvfc-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (vfc_to_nld[ivfc]) + suffix + "\n")
   file.write("real(kind_real),parameter :: vfc_minval(nvfc) = (/ &\n")
   for ivfc in range(0, nvfc):
      if ivfc != nvfc-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (vfc_minval[ivfc]) + suffix + "\n")
   file.write("real(kind_real),parameter :: vfc_ratio(nvfc) = (/ &\n")
   for ivfc in range(0, nvfc):
      if ivfc != nvfc-1:
         suffix = ", &"
      else:
         suffix = "/)"
      file.write(" & %.8f_kind_real" % (vfc_ratio[ivfc]) + suffix + "\n")
   file.write("real(kind_real),parameter :: func_vfc(ncv,nvfc) = reshape((/ &\n")
   for ivfc in range(0, nvfc):
      for icv in range(0, ncv):
            if ivfc != nvfc-1 or icv != ncv-1:
               suffix = ","
            else:
               suffix = "/),"
            file.write(" & %.8f_kind_real" % (f_int_vfc[icv,ivfc]) + suffix + " &\n")
   file.write(" & (/ncv,nvfc/))\n")
   file.write("\n")
   file.write("interface fit_func_hor\n")
   file.write("   module procedure gc99_fit_func_hor\n")
   file.write("end interface\n")
   file.write("interface fit_func_ver\n")
   file.write("   module procedure gc99_fit_func_ver\n")
   file.write("end interface\n")
   file.write("interface fit_func_hor_sqrt\n")
   file.write("   module procedure gc99_fit_func_hor_sqrt\n")
   file.write("end interface\n")
   file.write("interface fit_func_ver_sqrt\n")
   file.write("   module procedure gc99_fit_func_ver_sqrt\n")
   file.write("end interface\n")
   file.write("\n")
   file.write("public ! TODO: make useless public stuff private\n")
   file.write("\n")
   file.write("contains\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Function: gc99_fit_func_hor\n")
   file.write("!> Fit function, horizontal component\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("function gc99_fit_func_hor(mpl,ch,pkh) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("real(kind_real),intent(in) :: ch    !< Horizontal normalized distance\n")
   file.write("real(kind_real),intent(in) :: pkh   !< Horizontal peakness\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Local variables\n")
   file.write("integer :: ichm,ichp,ipkhm,ipkhp\n")
   file.write("real(kind_real) :: bch,bpkh\n")
   file.write("real(kind_real) :: rchm,rchp,rpkhm,rpkhp\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_func_hor)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Check bounds\n")
   file.write("if (ch<zero) call mpl%abort('${subr}$','negative horizontal normalized distance')\n")
   file.write("if (pkh<zero) call mpl%abort('${subr}$','negative horizontal peakness')\n")
   file.write("\n")
   file.write("if (eq(ch,zero)) then\n")
   file.write("   ! Normalized origin\n")
   file.write("   value = one\n")
   file.write("elseif (sup(ch,one)) then\n")
   file.write("   ! Out of support\n")
   file.write("   value = zero\n")
   file.write("else\n")
   file.write("   ! Bounded values\n")
   file.write("   bch = max(chmin,min(ch,chmax))\n")
   file.write("   bpkh = max(pkhmin,min(pkh,pkhmax))\n")
   file.write("\n")
   file.write("   ! Indices\n")
   file.write("   ichm = floor(bch/dch)+1\n")
   file.write("   ichp = ichm+1\n")
   file.write("   ipkhm = floor(bpkh/dpkh)+1\n")
   file.write("   ipkhp = ipkhm+1\n")
   file.write("\n")
   file.write("   ! Coefficients\n")
   file.write("   rchm = real(ichp-1,kind_real)-bch/dch\n")
   file.write("   rchp = (one-rchm)\n")
   file.write("   rpkhm = real(ipkhp-1,kind_real)-bpkh/dpkh\n")
   file.write("   rpkhp = (one-rpkhm)\n")
   file.write("\n")
   file.write("   ! Interpolated value\n")
   file.write("   value = rchm*rpkhm*func_hor(ichm,ipkhm) &\n")
   file.write(" & +rchp*rpkhm*func_hor(ichp,ipkhm) &\n")
   file.write(" & +rchm*rpkhp*func_hor(ichm,ipkhp) &\n")
   file.write(" & +rchp*rpkhp*func_hor(ichp,ipkhp)\n")
   file.write("end if\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end function gc99_fit_func_hor\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Function: gc99_fit_func_ver\n")
   file.write("!> Fit function, vertical component\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("function gc99_fit_func_ver(mpl,cv,vfc) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("real(kind_real),intent(in) :: cv    !< Vertical normalized distance\n")
   file.write("integer,intent(in) :: vfc           !< Vertical functions catalog index\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Local variables\n")
   file.write("integer :: icvm,icvp\n")
   file.write("real(kind_real) :: bcv\n")
   file.write("real(kind_real) :: rcvm,rcvp\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_func_ver)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Check bounds\n")
   file.write("if (cv<zero) call mpl%abort('${subr}$','negative vertical normalized distance')\n")
   file.write("if ((vfc<1).or.(vfc>nvfc)) call mpl%abort('${subr}$','vertical functions catalog index out of bounds')\n")
   file.write("\n")
   file.write("if (eq(cv,zero)) then\n")
   file.write("   ! Normalized origin\n")
   file.write("   value = one\n")
   file.write("elseif (sup(cv,one)) then\n")
   file.write("   ! Out of support\n")
   file.write("   value = zero\n")
   file.write("else\n")
   file.write("   ! Bounded values\n")
   file.write("   bcv = max(cvmin,min(cv,cvmax))\n")
   file.write("\n")
   file.write("   ! Indices\n")
   file.write("   icvm = floor(bcv/dcv)+1\n")
   file.write("   icvp = icvm+1\n")
   file.write("\n")
   file.write("   ! Coefficients\n")
   file.write("   rcvm = real(icvp-1,kind_real)-bcv/dcv\n")
   file.write("   rcvp = (one-rcvm)\n")
   file.write("\n")
   file.write("   ! Interpolated value\n")
   file.write("   value = rcvm*func_vfc(icvm,vfc)+rcvp*func_vfc(icvp,vfc)\n")
   file.write("end if\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end function gc99_fit_func_ver\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Function: gc99_fit_func_hor_sqrt\n")
   file.write("!> Fit function function square-root, horizontal component\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("function gc99_fit_func_hor_sqrt(mpl,ch,pkh) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("real(kind_real),intent(in) :: ch    !< Horizontal normalized distance\n")
   file.write("real(kind_real),intent(in) :: pkh   !< Horizontal peakness\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_func_hor_sqrt)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Check bounds\n")
   file.write("if (ch<zero) call mpl%abort('${subr}$','negative horizontal normalized distance')\n")
   file.write("if (pkh<zero) call mpl%abort('${subr}$','negative horizontal peakness')\n")
   file.write("\n")
   file.write("! Horizontal contribution\n")
   file.write("if (eq(ch,zero)) then\n")
   file.write("   ! Normalized origin\n")
   file.write("   value = one\n")
   file.write("elseif (sup(ch,half)) then\n")
   file.write("   ! Out of support\n")
   file.write("   value = zero\n")
   file.write("else\n")
   file.write("   ! Horizontal component\n")
   file.write("   if (pkh>zero) then\n")
   file.write("      value = exp(log(f0+exp(-pkh))*two*sqrt(ch))-exp(-pkh)\n")
   file.write("   else\n")
   file.write("      value = one-two*ch\n")
   file.write("   end if\n")
   file.write("end if\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end function gc99_fit_func_hor_sqrt\n")
   file.write("\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("! Function: gc99_fit_func_ver_sqrt\n")
   file.write("!> Fit function function square-root, vertical component\n")
   file.write("!----------------------------------------------------------------------\n")
   file.write("function gc99_fit_func_ver_sqrt(mpl,cv,pkv,nld) result(value)\n")
   file.write("\n")
   file.write("! Passed variables\n")
   file.write("type(mpl_type),intent(inout) :: mpl !< MPI data\n")
   file.write("real(kind_real),intent(in) :: cv    !< Vertical normalized distance\n")
   file.write("real(kind_real),intent(in) :: pkv   !< Vertical peakness\n")
   file.write("real(kind_real),intent(in) :: nld   !< Negative lobes depth\n")
   file.write("\n")
   file.write("! Returned variable\n")
   file.write("real(kind_real) :: value\n")
   file.write("\n")
   file.write("! Set name\n")
   file.write("@:set_name(gc99_fit_func_ver_sqrt)\n")
   file.write("\n")
   file.write("! Probe in\n")
   file.write("@:probe_in()\n")
   file.write("\n")
   file.write("! Check bounds\n")
   file.write("if (cv<zero) call mpl%abort('${subr}$','negative vertical normalized distance')\n")
   file.write("if (pkv<zero) call mpl%abort('${subr}$','negative vertical peakness')\n")
   file.write("if (nld<zero) call mpl%abort('${subr}$','negative negative lobes depth')\n")
   file.write("\n")
   file.write("! Horizontal contribution\n")
   file.write("if (eq(cv,zero)) then\n")
   file.write("   ! Normalized origin\n")
   file.write("   value = one\n")
   file.write("elseif (sup(cv,half)) then\n")
   file.write("   ! Out of support\n")
   file.write("   value = zero\n")
   file.write("else\n")
   file.write("   ! Vertical component\n")
   file.write("   value = zero\n")
   file.write("   if (cv<half*(one-tau)) then\n")
   file.write("      if (pkv>zero) then\n")
   file.write("         value = value+(one+nld)*((exp(log(f0+exp(-pkv))*two*sqrt(cv)/(one-tau))-exp(-pkv))/(one-exp(-pkv)))\n")
   file.write("      else\n")
   file.write("         value = value+(one+nld)*(one-two*cv/(one-tau))\n")
   file.write("      end if\n")
   file.write("   end if\n")
   file.write("   if (cv<half) value = value-nld*(one-two*cv)\n")
   file.write("end if\n")
   file.write("\n")
   file.write("! Probe out\n")
   file.write("@:probe_out()\n")
   file.write("\n")
   file.write("end function gc99_fit_func_ver_sqrt\n")
   file.write("\n")
   file.write("end module tools_gc99")
   
   # Close file
   file.close()
