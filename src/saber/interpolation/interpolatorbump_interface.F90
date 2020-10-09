! (C) Copyright 2020- UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module interpolatorbump_interface

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding
use missing_values_mod
use bump_interpolation_mod
use type_fieldset

private
! ------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Create Bump Interpolator (abbreviated as bint)
!!
subroutine bint_create_c(c_key_bint, c_comm, c_fspace1, c_fspace2, c_masks, &
                         c_config) bind(c, name='bint_create_f90')
implicit none

! Passed variables
integer(c_int), intent(inout) :: c_key_bint !< bump interpolator
type(c_ptr), value, intent(in) :: c_comm    !< MPI Communicator
type(c_ptr), intent(in),value :: c_fspace1  !< source grid (atlas functionspace)
type(c_ptr), intent(in),value :: c_fspace2  !< target grid (atlas functionspace)
type(c_ptr), intent(in),value :: c_masks    !< masks and other metadata
type(c_ptr), value, intent(in) :: c_config  !< Configuration

! local variables
type(bump_interpolator), pointer :: bint
type(fckit_mpi_comm) :: f_comm
type(fckit_configuration) :: f_config
type(atlas_functionspace) :: fspace1, fspace2
type(fieldset_type) :: masks

f_comm = fckit_mpi_comm(c_comm)
f_config = fckit_configuration(c_config)

call bump_interpolator_registry%init()
call bump_interpolator_registry%add(c_key_bint)
call bump_interpolator_registry%get(c_key_bint, bint)

fspace1 = atlas_functionspace(c_fspace1)
fspace2 = atlas_functionspace(c_fspace2)

if (c_associated(c_masks)) then
   masks = atlas_fieldset(c_masks)
   call bint%init(f_config, f_comm, fspace1, fspace2, masks)
else
   call bint%init(f_config, f_comm, fspace1, fspace2)
endif

end subroutine bint_create_c

!-------------------------------------------------------------------------------
!> Apply Bump Interpolator
!!
subroutine bint_apply_c(c_key_bint, c_infields, c_outfields) bind(c, name='bint_apply_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_bint  !< key to bump interpolator
type(c_ptr), intent(in), value :: c_infields  !< input fields
type(c_ptr), intent(in), value :: c_outfields !< output fields

! Local variables
type(bump_interpolator), pointer :: bint
type(fieldset_type) :: infields, outfields

call bump_interpolator_registry%get(c_key_bint, bint)
infields = atlas_fieldset(c_infields)
outfields = atlas_fieldset(c_outfields)

call bint%apply(infields, outfields)

end subroutine bint_apply_c

!-------------------------------------------------------------------------------
!> Apply Bump Interpolator Adjoint
!!
subroutine bint_apply_ad_c(c_key_bint, c_fields2, c_fields1) bind(c, name='bint_apply_ad_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_bint    !< key to bump interpolator
type(c_ptr), intent(in), value :: c_fields2 !< input fields
type(c_ptr), intent(in), value :: c_fields1 !< output fields

! Local variables
type(bump_interpolator), pointer :: bint
type(fieldset_type) :: fields_grid2, fields_grid1

call bump_interpolator_registry%get(c_key_bint, bint)
fields_grid2 = atlas_fieldset(c_fields2)
fields_grid1 = atlas_fieldset(c_fields1)

call bint%apply_ad(fields_grid2, fields_grid1)

end subroutine bint_apply_ad_c

! ------------------------------------------------------------------------------
!> Delete bump_interpolator
subroutine bint_delete_c(c_key_bint) bind(c, name='bint_delete_f90')

implicit none

! Passed variables
integer(c_int), intent(inout) :: c_key_bint

! Local variables
type(bump_interpolator), pointer :: bint

call bump_interpolator_registry%get(c_key_bint, bint)

call bint%delete()

call bump_interpolator_registry%remove(c_key_bint)

end subroutine bint_delete_c

! ------------------------------------------------------------------------------

end module interpolatorbump_interface
