# (C) Copyright 2019 UCAR
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

if( HAVE_OMP )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
else( )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-openmp")
endif( )

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops -finline-functions ")

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wextra -Wall -pedantic -ftrapv -fall-intrinsics -fcheck=all -fimplicit-none -ffpe-trap=invalid,zero,underflow,overflow,denormal -fcheck-array-temporaries -finit-derived -finit-integer=-999 -finit-real=snan -fbacktrace" )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2 -funroll-all-loops -finline-functions" )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

####################################################################

# Meaning of flags
# ----------------
# -fstack-arrays     : Allocate automatic arrays on the stack (needs large stacksize!!!)
# -funroll-all-loops : Unroll all loops
# -fcheck=bounds     : Bounds checking
  
