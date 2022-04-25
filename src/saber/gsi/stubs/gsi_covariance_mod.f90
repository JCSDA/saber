! (C) Copyright 2022 United States Government as represented by the Administrator of the National
!     Aeronautics and Space Administration
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module gsi_covariance_mod

! atlas
use atlas_module,                   only: atlas_fieldset, atlas_field

! fckit
use fckit_mpi_module,               only: fckit_mpi_comm
use fckit_configuration_module,     only: fckit_configuration

! oops
use kinds,                          only: kind_real

implicit none
private
public gsi_covariance


! Fortran class header
type :: gsi_covariance
  integer :: stub
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: randomize
    procedure, public :: multiply
    procedure, public :: multiply_ad
end type gsi_covariance

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, comm, config)

! Arguments
class(gsi_covariance),     intent(inout) :: self
type(fckit_mpi_comm),      intent(in)    :: comm
type(fckit_configuration), intent(in)    :: config
!type(atlas_fieldset),      intent(in)    :: background   ! Uncomment once background available as Atlas fieldset
!type(atlas_fieldset),      intent(in)    :: first_guess

call abor1_ftn("gsi_covariance_mod.create This is a stub. To use this code recompile with GSI B")

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(gsi_covariance) :: self

call abor1_ftn("gsi_covariance_mod.delete This is a stub. To use this code recompile with GSI B")

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine randomize(self, fields)

!Arguments
class(gsi_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

call abor1_ftn("gsi_covariance_mod.randomize This is a stub. To use this code recompile with GSI B")

end subroutine randomize

! --------------------------------------------------------------------------------------------------

subroutine multiply(self, fields)

!Arguments
class(gsi_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

call abor1_ftn("gsi_covariance_mod.multiply This is a stub. To use this code recompile with GSI B")

end subroutine multiply

! --------------------------------------------------------------------------------------------------

subroutine multiply_ad(self, fields)

!Arguments
class(gsi_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

call abor1_ftn("gsi_covariance_mod.multiply_ad This is a stub. To use this code recompile with GSI B")

end subroutine multiply_ad

! --------------------------------------------------------------------------------------------------

end module gsi_covariance_mod
