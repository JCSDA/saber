! (C) Copyright 2022 United States Government as represented by the Administrator of the National
!     Aeronautics and Space Administration
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module gsi_grid_mod

! netcdf
use netcdf

! atlas
use atlas_module,                   only: atlas_field, atlas_fieldset, atlas_real

! fckit
use fckit_mpi_module,               only: fckit_mpi_comm
use fckit_configuration_module,     only: fckit_configuration

! oops
use kinds,                          only: kind_real

! saber
use gsi_utils_mod,                  only: nccheck

implicit none
private
public gsi_grid

! Fortran class header
type :: gsi_grid
  integer :: stub
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: print
    procedure, public :: get_levels
    procedure, public :: set_atlas_lonlat
end type gsi_grid

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, conf, comm)

! Arguments
class(gsi_grid),           intent(inout) :: self
type(fckit_configuration), intent(in)    :: conf
type(fckit_mpi_comm),      intent(in)    :: comm

call abor1_ftn("gsi_grid_mod.create This is a stub. To use this code recompile with GSI B")

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(gsi_grid), intent(inout) :: self

call abor1_ftn("gsi_grid_mod.delete This is a stub. To use this code recompile with GSI B")

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine print(self)

! Arguments
class(gsi_grid), intent(in) :: self

call abor1_ftn("gsi_grid_mod.print This is a stub. To use this code recompile with GSI B")

end subroutine print

! --------------------------------------------------------------------------------------------------

subroutine get_levels(self, levels)

! Arguments
class(gsi_grid), intent(in)    :: self
integer,         intent(inout) :: levels

call abor1_ftn("gsi_grid_mod.get_levels This is a stub. To use this code recompile with GSI B")

end subroutine get_levels

! --------------------------------------------------------------------------------------------------

subroutine set_atlas_lonlat(self, grid_fieldset)

!Arguments
class(gsi_grid),      intent(inout) :: self
type(atlas_fieldset), intent(inout) :: grid_fieldset

call abor1_ftn("gsi_grid_mod.set_atlas_lonlat This is a stub. To use this code recompile with GSI B")

end subroutine set_atlas_lonlat

! --------------------------------------------------------------------------------------------------

end module gsi_grid_mod
