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
use random_mod

! saber
use gsi_grid_mod,                   only: gsi_grid


implicit none
private
public gsi_covariance


! Fortran class header
type :: gsi_covariance
  type(gsi_grid) :: grid
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

! Create the grid
! ---------------
call self%grid%create(config, comm)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(gsi_covariance) :: self

! Delete the grid
! ---------------
call self%grid%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine randomize(self, fields)

!Arguments
class(gsi_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

! Locals
type(atlas_field) :: afield
real(kind=kind_real), pointer :: psi(:,:), chi(:,:), t(:,:), q(:,:), qi(:,:), ql(:,:), o3(:,:)
real(kind=kind_real), pointer :: ps(:)

integer, parameter :: rseed = 3

! Get Atlas field
afield = fields%field('stream_function')
call afield%data(psi)

afield = fields%field('velocity_potential')
call afield%data(chi)

afield = fields%field('air_temperature')
call afield%data(t)

afield = fields%field('surface_pressure')
call afield%data(ps)

afield = fields%field('specific_humidity')
call afield%data(q)

afield = fields%field('cloud_liquid_ice')
call afield%data(qi)

afield = fields%field('cloud_liquid_water')
call afield%data(ql)

afield = fields%field('ozone_mass_mixing_ratio')
call afield%data(o3)


! Set fields to random numbers
call normal_distribution(psi, 0.0_kind_real, 1.0_kind_real, rseed)


end subroutine randomize

! --------------------------------------------------------------------------------------------------

subroutine multiply(self, fields)

!Arguments
class(gsi_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

! To do list for this method
! 1. Convert fields (Atlas fieldsets) to GSI bundle
! 2. Call GSI covariance operator (sqrt version)
! 3. Convert back to Atlas Fields

end subroutine multiply

! --------------------------------------------------------------------------------------------------

subroutine multiply_ad(self, fields)

!Arguments
class(gsi_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

! To do list for this method
! 1. Convert fields (Atlas fieldsets) to GSI bundle
! 2. Call GSI covariance operator adjoint (sqrt version)
! 3. Convert back to Atlas Fields

end subroutine multiply_ad

! --------------------------------------------------------------------------------------------------

end module gsi_covariance_mod
