! (C) Copyright 2020- UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! ------------------------------------------------------------------------------
! Module: interpolatorbump_interface
!> BUMP interpolation interface
! ------------------------------------------------------------------------------
module interpolatorbump_interface

use atlas_module, only: atlas_functionspace,atlas_fieldset
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use interpolatorbump_mod, only: bump_interpolator
use iso_c_binding, only: c_int,c_ptr,c_associated
use type_fieldset, only: fieldset_type

private

! BUMP interpolator registry
#define LISTED_TYPE bump_interpolator
#include "../../../oops/src/oops/util/linkedList_i.f"
type(registry_t) :: bump_interpolator_registry

contains

!----------------------------------------------------------------------
! Linked list implementation
!----------------------------------------------------------------------
#include "../../../oops/src/oops/util/linkedList_c.f"

!-------------------------------------------------------------------------------
! Subroutine bint_create_c
!> Create BUMP interpolator (abbreviated as bint)
!-------------------------------------------------------------------------------
subroutine bint_create_c(c_key_bint,c_comm,c_fspace1,c_fspace2,c_masks,c_config) bind(c,name='bint_create_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_bint !< BUMP interpolator
type(c_ptr),value,intent(in) :: c_comm     !< MPI Communicator
type(c_ptr),intent(in),value :: c_fspace1  !< Source grid (atlas functionspace)
type(c_ptr),intent(in),value :: c_fspace2  !< Target grid (atlas functionspace)
type(c_ptr),intent(in),value :: c_masks    !< Masks and other metadata
type(c_ptr),value,intent(in) :: c_config   !< Configuration

! Local variables
type(bump_interpolator),pointer :: bint
type(fckit_mpi_comm) :: f_comm
type(fckit_configuration) :: f_config
type(atlas_functionspace) :: fspace1,fspace2
type(fieldset_type) :: masks

! Interface
f_comm = fckit_mpi_comm(c_comm)
f_config = fckit_configuration(c_config)
call bump_interpolator_registry%init()
call bump_interpolator_registry%add(c_key_bint)
call bump_interpolator_registry%get(c_key_bint,bint)
fspace1 = atlas_functionspace(c_fspace1)
fspace2 = atlas_functionspace(c_fspace2)
if (c_associated(c_masks)) then
   masks = atlas_fieldset(c_masks)
else
   masks = atlas_fieldset()
endif

! Call Fortran
call bint%init(f_comm,afunctionspace_in=fspace1,afunctionspace_out=fspace2,masks=masks,config=f_config)

end subroutine bint_create_c

!-------------------------------------------------------------------------------
! Subroutine: bint_apply_c
!> Apply BUMP interpolator
!-------------------------------------------------------------------------------
subroutine bint_apply_c(c_key_bint,c_infields,c_outfields) bind(c,name='bint_apply_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_bint     !< Key to BUMP interpolator
type(c_ptr),intent(in),value :: c_infields  !< Input fields
type(c_ptr),intent(in),value :: c_outfields !< Output fields

! Local variables
type(bump_interpolator),pointer :: bint
type(fieldset_type) :: infields,outfields

! Interface
call bump_interpolator_registry%get(c_key_bint,bint)
infields = atlas_fieldset(c_infields)
outfields = atlas_fieldset(c_outfields)

! Call Fortran
call bint%apply(infields,outfields)

end subroutine bint_apply_c

!-------------------------------------------------------------------------------
! Subroutine: bint_apply_ad_c
!> Apply BUMP interpolator adjoint
!-------------------------------------------------------------------------------
subroutine bint_apply_ad_c(c_key_bint,c_fields2,c_fields1) bind(c,name='bint_apply_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_bint   !< Key to BUMP interpolator
type(c_ptr),intent(in),value :: c_fields2 !< Input fields
type(c_ptr),intent(in),value :: c_fields1 !< Output fields

! Local variables
type(bump_interpolator),pointer :: bint
type(fieldset_type) :: fields_grid2,fields_grid1

! Interface
call bump_interpolator_registry%get(c_key_bint,bint)
fields_grid2 = atlas_fieldset(c_fields2)
fields_grid1 = atlas_fieldset(c_fields1)

! Call Fortran
call bint%apply_ad(fields_grid2,fields_grid1)

end subroutine bint_apply_ad_c

!----------------------------------------------------------------------
! Subroutine bint_delete_c
!> Delete bump_interpolator
!----------------------------------------------------------------------
subroutine bint_delete_c(c_key_bint) bind(c,name='bint_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_bint !< Key to BUMP interpolator

! Local variables
type(bump_interpolator),pointer :: bint

! Interface
call bump_interpolator_registry%get(c_key_bint,bint)

! Call Fortran
call bint%delete()

! Clean interface
call bump_interpolator_registry%remove(c_key_bint)

end subroutine bint_delete_c

end module interpolatorbump_interface
