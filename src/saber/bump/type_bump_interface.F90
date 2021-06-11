!----------------------------------------------------------------------
! Module: type_bump_interface
!> BUMP derived type interface
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright 2015-... UCAR,CERFACS,METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_bump_interface

use atlas_module, only: atlas_functionspace,atlas_fieldset
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding, only: c_int,c_ptr,c_double,c_char
use type_bump, only: bump_type
use type_fieldset, only: fieldset_type

implicit none

private

! BUMP registry
#define LISTED_TYPE bump_type
#include "oops/util/linkedList_i.f"
type(registry_t) :: bump_registry

contains

!----------------------------------------------------------------------
! Linked list implementation
!----------------------------------------------------------------------
#include "oops/util/linkedList_c.f"

!----------------------------------------------------------------------
! Subroutine: bump_create_c
!> Create
!----------------------------------------------------------------------
subroutine bump_create_c(key_bump,c_comm,c_afunctionspace,c_afieldset,c_conf,c_grid) bind(c,name='bump_create_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: key_bump         !< BUMP
type(c_ptr),intent(in),value :: c_comm           !< FCKIT MPI communicator wrapper
type(c_ptr),intent(in),value :: c_afunctionspace !< ATLAS function space
type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset containing geometry elements
type(c_ptr),intent(in),value :: c_conf           !< FCKIT configuration
type(c_ptr),intent(in),value :: c_grid           !< FCKIT grid configuration

! Local variables
type(bump_type),pointer :: bump
type(fckit_mpi_comm) :: f_comm
type(atlas_functionspace) :: f_afunctionspace
type(fieldset_type) :: f_fieldset
type(fckit_configuration) :: f_conf
type(fckit_configuration) :: f_grid

! Interface
call bump_registry%init()
call bump_registry%add(key_bump)
call bump_registry%get(key_bump,bump)
f_comm = fckit_mpi_comm(c_comm)
f_afunctionspace = atlas_functionspace(c_afunctionspace)
f_fieldset = atlas_fieldset(c_afieldset)
f_conf = fckit_configuration(c_conf)
f_grid = fckit_configuration(c_grid)

! Call Fortran
call bump%create(f_comm,f_afunctionspace,f_fieldset,f_conf,f_grid)

end subroutine bump_create_c

!----------------------------------------------------------------------
! Subroutine: bump_run_drivers_c
!> Run drivers
!----------------------------------------------------------------------
subroutine bump_run_drivers_c(key_bump) bind(c,name='bump_run_drivers_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump !< BUMP

! Local variables
type(bump_type),pointer :: bump

! Interface
call bump_registry%get(key_bump,bump)

! Call Fortran
call bump%run_drivers

end subroutine bump_run_drivers_c

!----------------------------------------------------------------------
! Subroutine: bump_add_member_c
!> Add member into bump%ens[1,2]
!----------------------------------------------------------------------
subroutine bump_add_member_c(key_bump,c_afieldset,ie,iens) bind(c,name='bump_add_member_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer
integer(c_int),intent(in) :: ie             !< Ensemble member index
integer(c_int),intent(in) :: iens           !< Ensemble index

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%add_member(f_fieldset,ie,iens)

end subroutine bump_add_member_c

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_c
!> Vertical balance application
!----------------------------------------------------------------------
subroutine bump_apply_vbal_c(key_bump,c_afieldset) bind(c,name='bump_apply_vbal_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%apply_vbal(f_fieldset)

end subroutine bump_apply_vbal_c

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_inv_c
!> Vertical balance application, inverse
!----------------------------------------------------------------------
subroutine bump_apply_vbal_inv_c(key_bump,c_afieldset) bind(c,name='bump_apply_vbal_inv_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%apply_vbal_inv(f_fieldset)

end subroutine bump_apply_vbal_inv_c

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_ad_c
!> Vertical balance application, adjoint
!----------------------------------------------------------------------
subroutine bump_apply_vbal_ad_c(key_bump,c_afieldset) bind(c,name='bump_apply_vbal_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%apply_vbal_ad(f_fieldset)

end subroutine bump_apply_vbal_ad_c

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_inv_ad_c
!> Vertical balance application, inverse adjoint
!----------------------------------------------------------------------
subroutine bump_apply_vbal_inv_ad_c(key_bump,c_afieldset) bind(c,name='bump_apply_vbal_inv_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%apply_vbal_inv_ad(f_fieldset)

end subroutine bump_apply_vbal_inv_ad_c

!----------------------------------------------------------------------
! Subroutine: bump_apply_stddev_c
!> Standard-deviation application
!----------------------------------------------------------------------
subroutine bump_apply_stddev_c(key_bump,c_afieldset) bind(c,name='bump_apply_stddev_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%apply_stddev(f_fieldset)

end subroutine bump_apply_stddev_c

!----------------------------------------------------------------------
! Subroutine: bump_apply_stddev_inv_c
!> Standard-deviation application, inverse
!----------------------------------------------------------------------
subroutine bump_apply_stddev_inv_c(key_bump,c_afieldset) bind(c,name='bump_apply_stddev_inv_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%apply_stddev_inv(f_fieldset)

end subroutine bump_apply_stddev_inv_c

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_c
!> NICAS application
!----------------------------------------------------------------------
subroutine bump_apply_nicas_c(key_bump,c_afieldset) bind(c,name='bump_apply_nicas_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%apply_nicas(f_fieldset)

end subroutine bump_apply_nicas_c

!----------------------------------------------------------------------
! Subroutine: bump_get_cv_size_c
!> Get control variable size
!----------------------------------------------------------------------
subroutine bump_get_cv_size_c(key_bump,n) bind(c,name='bump_get_cv_size_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump !< BUMP
integer(c_int),intent(out) :: n       !< Control variable size

! Local variables
type(bump_type),pointer :: bump

! Interface
call bump_registry%get(key_bump,bump)

! Call Fortran
call bump%get_cv_size(n)

end subroutine bump_get_cv_size_c

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_sqrt
!> NICAS square-root application
!----------------------------------------------------------------------
subroutine bump_apply_nicas_sqrt_c(key_bump,cv,c_afieldset) bind(c,name='bump_apply_nicas_sqrt_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
real(c_double),intent(in) :: cv(:)          !< Control variable
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%apply_nicas_sqrt(cv,f_fieldset)

end subroutine bump_apply_nicas_sqrt_c

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_sqrt_ad_c
!> NICAS square-root adjoint application
!----------------------------------------------------------------------
subroutine bump_apply_nicas_sqrt_ad_c(key_bump,c_afieldset,cv) bind(c,name='bump_apply_nicas_sqrt_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer
real(c_double),intent(inout) :: cv(:)       !< Control variable

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%apply_nicas_sqrt_ad(f_fieldset,cv)

end subroutine bump_apply_nicas_sqrt_ad_c

!----------------------------------------------------------------------
! Subroutine: bump_randomize_c
!> NICAS randomization
!----------------------------------------------------------------------
subroutine bump_randomize_c(key_bump,c_afieldset) bind(c,name='bump_randomize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%randomize(f_fieldset)

end subroutine bump_randomize_c

!----------------------------------------------------------------------
! Subroutine: bump_get_parameter_c
!> Get a parameter
!----------------------------------------------------------------------
subroutine bump_get_parameter_c(key_bump,nstr,cstr,c_afieldset) bind(c,name='bump_get_parameter_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
integer(c_int),intent(in) :: nstr           !< Parameter name size
character(c_char),intent(in) :: cstr(nstr)  !< Parameter name
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
integer :: istr
character(len=nstr) :: param
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
param = ''
do istr=1,nstr
  param = trim(param)//cstr(istr)
end do
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%get_parameter(param,f_fieldset)

end subroutine bump_get_parameter_c

!----------------------------------------------------------------------
! Subroutine: bump_set_parameter_c
!> Set a parameter
!----------------------------------------------------------------------
subroutine bump_set_parameter_c(key_bump,nstr,cstr,c_afieldset) bind(c,name='bump_set_parameter_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump       !< BUMP
integer(c_int),intent(in) :: nstr           !< Parameter name size
character(c_char),intent(in) :: cstr(nstr)  !< Parameter name
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(bump_type),pointer :: bump
integer :: istr
character(len=nstr) :: param
type(fieldset_type) :: f_fieldset

! Interface
call bump_registry%get(key_bump,bump)
param = ''
do istr=1,nstr
  param = trim(param)//cstr(istr)
end do
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call bump%set_parameter(param,f_fieldset)

end subroutine bump_set_parameter_c

!----------------------------------------------------------------------
! Subroutine: bump_partial_dealloc_c
!> Partial deallocation
!----------------------------------------------------------------------
subroutine bump_partial_dealloc_c(key_bump) bind(c,name='bump_partial_dealloc_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_bump !< BUMP

! Local variables
type(bump_type),pointer :: bump

! Interface
call bump_registry%get(key_bump,bump)

! Partially deallocate BUMP
call bump%partial_dealloc

end subroutine bump_partial_dealloc_c

!----------------------------------------------------------------------
! Subroutine: bump_dealloc_c
!> Deallocation
!----------------------------------------------------------------------
subroutine bump_dealloc_c(key_bump) bind(c,name='bump_dealloc_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: key_bump !< BUMP

! Local variables
type(bump_type),pointer :: bump

! Interface
call bump_registry%get(key_bump,bump)

! Call Fortran
call bump%dealloc

! Clean interface
call bump_registry%remove(key_bump)

end subroutine bump_dealloc_c

end module type_bump_interface
