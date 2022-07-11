! (C) Copyright 2022 United States Government as represented by the Administrator of the National
!     Aeronautics and Space Administration
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module gsi_grid_interface_mod

! iso
use iso_c_binding

! atlas
use atlas_module,               only: atlas_functionspace, atlas_fieldset

! fckit
use fckit_mpi_module,           only: fckit_mpi_comm
use fckit_configuration_module, only: fckit_configuration

! oops
use oops_variables_mod

! saber
use gsi_grid_mod,         only: gsi_grid


implicit none
private
public gsi_grid_registry


! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE gsi_grid

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: gsi_grid_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine gsi_grid_create_cpp(c_self, c_conf, c_comm) &
           bind(c, name='gsi_grid_create_f90')

! Arguments
integer(c_int),     intent(inout) :: c_self
type(c_ptr), value, intent(in)    :: c_conf
type(c_ptr), value, intent(in)    :: c_comm

! Locals
type(gsi_grid), pointer   :: f_self
type(fckit_mpi_comm)      :: f_comm
type(fckit_configuration) :: f_conf

! LinkedList
! ----------
call gsi_grid_registry%init()
call gsi_grid_registry%add(c_self)
call gsi_grid_registry%get(c_self, f_self)

! Fortran APIs
! ------------
f_conf = fckit_configuration(c_conf)
f_comm = fckit_mpi_comm(c_comm)

! Call implementation
! -------------------
call f_self%create(f_conf, f_comm)

end subroutine gsi_grid_create_cpp

! --------------------------------------------------------------------------------------------------

subroutine gsi_grid_delete_cpp(c_self) &
           bind(c, name='gsi_grid_delete_f90')

! Arguments
integer(c_int), intent(inout) :: c_self

! Locals
type(gsi_grid), pointer :: f_self

! LinkedList
! ----------
call gsi_grid_registry%get(c_self, f_self)

! Call implementation
! -------------------
call f_self%delete()

! LinkedList
! ----------
call gsi_grid_registry%remove(c_self)

end subroutine gsi_grid_delete_cpp

! --------------------------------------------------------------------------------------------------

subroutine gsi_grid_print_cpp(c_self) &
           bind(c, name='gsi_grid_print_f90')

! Arguments
integer(c_int), intent(inout) :: c_self

! Locals
type(gsi_grid), pointer :: f_self

! LinkedList
! ----------
call gsi_grid_registry%get(c_self, f_self)

! Call implementation
! -------------------
call f_self%print()

end subroutine gsi_grid_print_cpp

! --------------------------------------------------------------------------------------------------

subroutine gsi_grid_get_levels_cpp(c_self, c_levels) &
           bind(c, name='gsi_grid_get_levels_f90')

! Arguments
integer(c_int), intent(inout) :: c_self
integer(c_int), intent(inout) :: c_levels

! Locals
type(gsi_grid), pointer :: f_self
integer                 :: f_levels

! LinkedList
! ----------
call gsi_grid_registry%get(c_self, f_self)

! Call implementation
! -------------------
call f_self%get_levels(f_levels)

! Convert any precision differeces
! --------------------------------
c_levels = int(f_levels, c_int)

end subroutine gsi_grid_get_levels_cpp

! --------------------------------------------------------------------------------------------------

subroutine gsi_grid_set_atlas_lonlat_cpp(c_self, c_grid_fieldset) &
                                           bind(c,name='gsi_grid_set_atlas_lonlat_f90')

implicit none

!Arguments
integer(c_int),     intent(in) :: c_self
type(c_ptr), value, intent(in) :: c_grid_fieldset

type(gsi_grid), pointer :: f_self
type(atlas_fieldset)    :: f_grid_fieldset

! LinkedList
! ----------
call gsi_grid_registry%get(c_self, f_self)

! Fortran APIs
! ------------
f_grid_fieldset = atlas_fieldset(c_grid_fieldset)

! Call implementation
! -------------------
call f_self%set_atlas_lonlat(f_grid_fieldset)

end subroutine gsi_grid_set_atlas_lonlat_cpp

! --------------------------------------------------------------------------------------------------

end module gsi_grid_interface_mod
