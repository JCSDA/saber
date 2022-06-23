!-------------------------------------------------------------------------------
! (C) Crown Copyright 2022 Met Office.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!-------------------------------------------------------------------------------
! Derived types module for CVT
!-------------------------------------------------------------------------------
module cvt_derivedtypes_mod

use iso_c_binding, only: c_int, c_int32_t
USE netcdf, only: &
    nf90_max_name

implicit none

! ---------------------------------------------------------------------------
! 1. Derived types:
!    NetCDF coordinate attributes
! ---------------------------------------------------------------------------
type cvt_coordinate_type
  character(len=nf90_max_name)  :: CoName
  character(len=nf90_max_name)  :: long_name
  character(len=nf90_max_name)  :: units
  character(len=nf90_max_name)  :: positive
  integer                       :: CoSize
  integer(kind=c_int32_t)       :: dimid
  integer(kind=c_int)           :: varid
  real                          :: origin
  real                          :: delta
  real, allocatable             :: CoValues(:)
end type cvt_coordinate_type

end module cvt_derivedtypes_mod
