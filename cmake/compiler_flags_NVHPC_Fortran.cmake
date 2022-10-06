# (C) Copyright 2019 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ktrap=fp")
if( HAVE_OMP )
    if ( CMAKE_USE_GPU )
        set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp=gpu -gpu=cuda11.7 -Minfo -Mpreprocess" )
    else()
        set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp=multiproc -Minfo -Mpreprocess" )
    endif()
else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -nomp")
endif()

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 --fast" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Mbounds -Mchkstk -traceback" )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2 --fast" )

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
