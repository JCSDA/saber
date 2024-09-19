#!/usr/bin/env python3

import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("srcdir", help="SABER source directory")
args = parser.parse_args()

# General parameters
nc = 41
cor = np.linspace(0.0, 1.0, nc)
nne = 40
ne = np.logspace(1, 4, nne, dtype=int)
nsample = 100000

# Initialize array
cor_var = np.zeros((nc, nne))

# Initialize random number generator
rng = np.random.default_rng(9)

# Generate normal distribution
z_base = rng.standard_normal((nsample))

for ine in range(0, nne):
   for ic in range(0, nc):
      if cor[ic] < 1.0:
         # Rescale and offset distribution
         z = z_base*np.sqrt(1.0/float(ne[ine]-3))+0.5*np.log((1.0+cor[ic])/(1.0-cor[ic]))

         # Apply inverse Fisher transformation
         rho = np.tanh(z)

         # Compute correlation variance
         cor_var[ic, ine] = np.var(rho)
      else:
         cor_var[ic, ine] = 0.0
      print(str(ne[ine]) + " members for cor = " + str(np.round(cor[ic],3)) + " => var = " + np.format_float_scientific(cor_var[ic,ine], unique=False, precision=4))

# Plot result
fig, ax = plt.subplots(figsize=(7,7))
ax.set_xlim([0,1.0])
ax.set_ylim([0,0.4])
ax.set_xlabel("Correlation")
ax.set_ylabel("Standard-deviation")
ax.plot(cor, np.sqrt(cor_var))
plt.savefig("cor_var.jpg", format="jpg", dpi=300)
plt.close()

# Open file
file = open(args.srcdir + "/src/saber/bump/tools_cor_var.fypp", "w")

# Write file
file.write("#:include 'instrumentation.fypp'\n")
file.write("#:include 'generics.fypp'\n")
file.write("!----------------------------------------------------------------------\n")
file.write("! Module: tools_cor_var\n")
file.write("!> Correlation estimator look-up table\n")
file.write("! Author: Benjamin Menetrier\n")
file.write("! Licensing: this code is distributed under the CeCILL-C license\n")
file.write("! Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT\n")
file.write("! WARNING: this module is generated by the python script\n")
file.write("!            tools/saber_correlation_variance.py\n")
file.write("!          to modify this module, update and rerun the python script\n")
file.write("!----------------------------------------------------------------------\n")
file.write("module tools_cor_var\n")
file.write("\n")
file.write("use tools_const, only: zero,one\n")
file.write("use tools_kinds, only: kind_real\n")
file.write("use tools_repro, only: eq,inf,infeq,sup\n")
file.write("use type_mpl, only: mpl_type\n")
file.write("@:use_probe()\n")
file.write("\n")
file.write("implicit none\n")
file.write("\n")
file.write("! Public parameters\n")
file.write("integer,parameter :: nc = " + str(nc) + "\n")
file.write("integer,parameter :: nne = " + str(nne) + "\n")
file.write("real(kind_real),parameter :: cor(nc) = (/ &\n")
for ic in range(0, nc):
   if ic != nc-1:
      suffix = ", &"
   else:
      suffix = "/)"
   file.write(" & %.8f_kind_real" % (cor[ic]) + suffix + "\n")
file.write("integer,parameter :: ne(nne) = (/ &\n")
for ine in range(0, nne):
   if ine != nne-1:
      suffix = ", &"
   else:
      suffix = "/)"
   file.write(str(ne[ine]) + suffix + "\n")
file.write("real(kind_real),parameter :: cor_var(nc,nne) = reshape((/ &\n")
for ine in range(0, nne):
   for ic in range(0, nc):
      if ic != nc-1 or ine != nne-1:
         suffix = ", &"
      else:
         suffix = "/), (/nc,nne/))"
      file.write(" & %.16f_kind_real" % (cor_var[ic,ine]) + suffix + "\n")
file.write("\n")
file.write("interface eval\n")
file.write("   module procedure cor_var_eval\n")
file.write("end interface\n")
file.write("\n")
file.write("private\n")
file.write("public :: eval\n")
file.write("\n")
file.write("contains\n")
file.write("\n")
file.write("!----------------------------------------------------------------------\n")
file.write("! Function: cor_var_eval\n")
file.write("!> Correlation variance\n")
file.write("!----------------------------------------------------------------------\n")
file.write("function cor_var_eval(mpl,cor_in,ne_in) result(var)\n")
file.write("\n")
file.write("! Passed variables\n")
file.write("type(mpl_type),intent(inout) :: mpl  !< MPI data\n")
file.write("real(kind_real),intent(in) :: cor_in !< Input correlation\n")
file.write("integer,intent(in) :: ne_in          !< Input ensemble size\n")
file.write("\n")
file.write("! Returned variable\n")
file.write("real(kind_real) :: var\n")
file.write("\n")
file.write("! Local variables\n")
file.write("integer :: ic,icm,icp,ine,inem,inep\n")
file.write("real(kind_real) :: rcm,rcp,rnem,rnep\n")
file.write("\n")
file.write("! Set name\n")
file.write("@:set_name(cor_var_eval)\n")
file.write("\n")
file.write("! Probe in\n")
file.write("@:probe_in()\n")
file.write("\n")
file.write("! Check bounds\n")
file.write("if (sup(abs(cor_in),one)) call mpl%abort('${subr}$','absolute correlation over one')\n")
file.write("if (inf(ne_in,minval(ne))) call mpl%abort('${subr}$','ensemble size too small')\n")
file.write("if (sup(ne_in,maxval(ne))) call mpl%abort('${subr}$','ensemble size too large')\n")
file.write("\n")
file.write("! Indices and coefficients\n")
file.write("do ic=1,nc-1\n")
file.write("   if (infeq(cor(ic),cor_in).and.inf(cor_in,cor(ic+1))) then\n")
file.write("      icm = ic\n")
file.write("      icp = ic+1\n")
file.write("      rcm = (cor(icp)-cor_in)/(cor(icp)-cor(icm))\n")
file.write("      rcp = (cor_in-cor(icm))/(cor(icp)-cor(icm))\n")
file.write("   end if\n") 
file.write("end do\n")
file.write("if (eq(cor_in,cor(nc))) then\n")
file.write("   icm = nc-1\n")
file.write("   icp = nc\n")
file.write("   rcm = zero\n")
file.write("   rcp = one\n")
file.write("end if\n")
file.write("do ine=1,nne-1\n")
file.write("   if ((ne(ine)<=ne_in).and.(ne_in<ne(ine+1))) then\n")
file.write("      inem = ine\n")
file.write("      inep = ine+1\n")
file.write("      rnem = real(ne(inep)-ne_in,kind_real)/real(ne(inep)-ne(inem),kind_real)\n")
file.write("      rnep = real(ne_in-ne(inem),kind_real)/real(ne(inep)-ne(inem),kind_real)\n")
file.write("   end if\n") 
file.write("end do\n")
file.write("if (ne_in==ne(nne)) then\n")
file.write("   inem = nne-1\n")
file.write("   inep = nne\n")
file.write("   rnem = zero\n")
file.write("   rnep = one\n")
file.write("end if\n") 
file.write("\n")
file.write("! Interpolated value\n")
file.write("var = rcm*rnem*cor_var(icm,inem) &\n")
file.write(" & + rcm*rnep*cor_var(icm,inep) &\n")
file.write(" & + rcp*rnem*cor_var(icp,inem) &\n")
file.write(" & + rcp*rnep*cor_var(icp,inep)\n")
file.write("\n")
file.write("! Probe out\n")
file.write("@:probe_out()\n")
file.write("\n")
file.write("end function cor_var_eval\n")
file.write("\n")
file.write("end module tools_cor_var\n")

# Close file
file.close()