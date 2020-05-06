! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module oobump_mod

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding
use kinds
use missing_values_mod
use type_bump, only: bump_type

implicit none

private
public :: oobump_type
public :: oobump_registry
public :: oobump_create, oobump_delete, oobump_get_cv_size, oobump_add_member, oobump_remove_member, oobump_run_drivers, &
        & oobump_multiply_vbal, oobump_multiply_vbal_inv, oobump_multiply_vbal_ad, oobump_multiply_vbal_inv_ad, &
        & oobump_multiply_stddev, oobump_multiply_stddev_inv, &
        & oobump_multiply_nicas, oobump_multiply_nicas_sqrt, oobump_multiply_nicas_sqrt_ad, &
        & oobump_randomize_nicas, oobump_get_param, oobump_set_param
! ------------------------------------------------------------------------------
type oobump_type
   type(bump_type) :: bump !> Instances of BUMP
contains
   final :: dummy
end type oobump_type

#define LISTED_TYPE oobump_type

!> Linked list interface - defines registry_t type
#include "saber/util/linkedList_i.f"

!> Global registry
type(registry_t) :: oobump_registry
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Linked list implementation
#include "saber/util/linkedList_c.f"
!-------------------------------------------------------------------------------
!> Create OOBUMP
subroutine oobump_create(self, f_comm, afunctionspace, afieldset, fconf, fgrid)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self                !< OOBUMP
type(fckit_mpi_comm),intent(in) :: f_comm               !< FCKIT MPI communicator wrapper
type(atlas_functionspace), intent(in) :: afunctionspace !< ATLAS function space
type(atlas_fieldset), intent(in) :: afieldset           !< ATLAS fieldset
type(fckit_configuration),intent(in) :: fconf           !< FCKIT configuration
type(fckit_configuration),intent(in) :: fgrid           !< FCKIT grid configuration

! Local variables
real(kind_real) :: msvalr

! Initialize namelist
call self%bump%nam%init(f_comm%size())

! Read configuration
call self%bump%nam%from_conf(fconf)

! Read grid configuration
call self%bump%nam%from_conf(fgrid)

! Get missing value
call fconf%get_or_die('msvalr', msvalr)

! Setup BUMP
call self%bump%setup(f_comm,afunctionspace,afieldset,msvalr=msvalr)

end subroutine oobump_create
!-------------------------------------------------------------------------------
!> Delete OOBUMP
subroutine oobump_delete(self)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP

! Release memory      
call self%bump%dealloc

end subroutine oobump_delete
!-------------------------------------------------------------------------------
!> Get control variable size
subroutine oobump_get_cv_size(self,n)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP
integer, intent(out) :: n                !< Control variable size

! Get control variable size
call self%bump%get_cv_size(n)

end subroutine oobump_get_cv_size
!-------------------------------------------------------------------------------
!> Add ensemble member
subroutine oobump_add_member(self,afieldset,ie,iens)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset
integer, intent(in) :: ie                        !< Ensemble member index
integer, intent(in) :: iens                      !< Ensemble index

! Add member
call self%bump%add_member(afieldset,ie,iens)

end subroutine oobump_add_member
!-------------------------------------------------------------------------------
!> Remove ensemble member
subroutine oobump_remove_member(self,afieldset,ie,iens)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset
integer, intent(in) :: ie                        !< Ensemble member index
integer, intent(in) :: iens                      !< Ensemble index

! Remove member
call self%bump%remove_member(afieldset,ie,iens)

end subroutine oobump_remove_member
!-------------------------------------------------------------------------------
!> Run OOBUMP drivers
subroutine oobump_run_drivers(self)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP

! Run BUMP drivers
call self%bump%run_drivers

end subroutine oobump_run_drivers
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator
subroutine oobump_multiply_vbal(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply vertical balance
call self%bump%apply_vbal(afieldset)

end subroutine oobump_multiply_vbal
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator inverse
subroutine oobump_multiply_vbal_inv(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply vertical balance, inverse
call self%bump%apply_vbal_inv(afieldset)

end subroutine oobump_multiply_vbal_inv
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint
subroutine oobump_multiply_vbal_ad(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply vertical balance, adjoint
call self%bump%apply_vbal_ad(afieldset)

end subroutine oobump_multiply_vbal_ad
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint inverse
subroutine oobump_multiply_vbal_inv_ad(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply vertical balance, inverse adjoint
call self%bump%apply_vbal_inv_ad(afieldset)

end subroutine oobump_multiply_vbal_inv_ad
!-------------------------------------------------------------------------------
!> Multiplication by BUMP standard-deviation
subroutine oobump_multiply_stddev(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply vertical balance
call self%bump%apply_stddev(afieldset)

end subroutine oobump_multiply_stddev
!-------------------------------------------------------------------------------
!> Multiplication by BUMP standard-deviation, inverse
subroutine oobump_multiply_stddev_inv(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply vertical balance
call self%bump%apply_stddev_inv(afieldset)

end subroutine oobump_multiply_stddev_inv
!-------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator
subroutine oobump_multiply_nicas(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply NICAS
call self%bump%apply_nicas(afieldset)

end subroutine oobump_multiply_nicas
!-------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root
subroutine oobump_multiply_nicas_sqrt(self,cv,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
real(kind_real), intent(in) :: cv(:)             !< Control variable
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply NICAS square-root
call self%bump%apply_nicas_sqrt(cv, afieldset)

end subroutine oobump_multiply_nicas_sqrt
!-------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root adjoint
subroutine oobump_multiply_nicas_sqrt_ad(self,afieldset,cv)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset
real(kind_real), intent(inout) :: cv(:)          !< Control variable

! Apply NICAS square-root
call self%bump%apply_nicas_sqrt_ad(afieldset, cv)

end subroutine oobump_multiply_nicas_sqrt_ad
!-------------------------------------------------------------------------------
!> Randomize the BUMP NICAS operator
subroutine oobump_randomize_nicas(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Randomize NICAS
call self%bump%randomize(afieldset)

end subroutine oobump_randomize_nicas
!-------------------------------------------------------------------------------
!> Get BUMP parameter
subroutine oobump_get_param(self,param,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
character(len=*), intent(in) :: param            !< Parameter name
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Get parameter
call self%bump%get_parameter(param,afieldset)

end subroutine oobump_get_param
!-------------------------------------------------------------------------------
!> Set BUMP parameter
subroutine oobump_set_param(self,param,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
character(len=*), intent(in) :: param            !< Parameter name
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Set parameter
call self%bump%set_parameter(param,afieldset)

end subroutine oobump_set_param
!-------------------------------------------------------------------------------
!> Dummy finalization
subroutine dummy(self)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP

end subroutine dummy
!-------------------------------------------------------------------------------
end module oobump_mod
