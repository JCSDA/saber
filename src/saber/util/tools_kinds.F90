!----------------------------------------------------------------------
! Module: tools_kinds
!> Kinds definition
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_kinds

use iso_c_binding
use netcdf, only: nf90_double

implicit none

! Kinds
integer,parameter :: kind_int = c_int           ! Integer kind
integer,parameter :: kind_short = c_short       ! Short integer kind
integer,parameter :: kind_real = c_double       ! Real kind

! NetCDF kinds
integer,parameter :: nc_kind_real = nf90_double ! NetCDF real kind

! Huge
real(kind_real),parameter :: huge_real = huge(0.0_kind_real)

private
public kind_int,kind_short,kind_real,nc_kind_real,huge_real

end module tools_kinds
