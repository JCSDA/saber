! (C) Copyright 2022 United States Government as represented by the Administrator of the National
!     Aeronautics and Space Administration
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module gsi_covariance_interface_mod

! iso
use iso_c_binding

! oops
use datetime_mod

! atlas
use atlas_module,               only: atlas_functionspace, atlas_fieldset

! fckit
use fckit_mpi_module,           only: fckit_mpi_comm
use fckit_configuration_module, only: fckit_configuration

! saber
use gsi_covariance_mod,         only: gsi_covariance


implicit none
private
public gsi_covariance_registry


! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE gsi_covariance

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: gsi_covariance_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine gsi_covariance_create_cpp(c_self, c_comm, c_conf, c_ntimes, c_bg, c_fg, c_valid_times) &
           bind(c, name='gsi_covariance_create_f90')

! Arguments
integer(c_int),     intent(inout) :: c_self
type(c_ptr), value, intent(in)    :: c_conf
type(c_ptr), value, intent(in)    :: c_comm
integer(c_int),     intent(in)    :: c_ntimes
type(c_ptr),        intent(in)    :: c_bg(c_ntimes)
type(c_ptr),        intent(in)    :: c_fg(c_ntimes)
type(c_ptr),        intent(in)    :: c_valid_times(c_ntimes)

! Locals
type(gsi_covariance), pointer :: f_self
type(fckit_mpi_comm)          :: f_comm
type(fckit_configuration)     :: f_conf
type(atlas_fieldset), dimension(:), allocatable :: f_bg
type(atlas_fieldset), dimension(:), allocatable :: f_fg
type(datetime), dimension(:), allocatable       :: f_valid_times
integer                       :: itime, ntimes

! LinkedList
! ----------
call gsi_covariance_registry%init()
call gsi_covariance_registry%add(c_self)
call gsi_covariance_registry%get(c_self, f_self)

! Fortran APIs
! ------------
f_conf = fckit_configuration(c_conf)
f_comm = fckit_mpi_comm(c_comm)

ntimes = c_ntimes
allocate(f_bg(ntimes), f_fg(ntimes), f_valid_times(ntimes))
do itime = 1, ntimes
  f_bg(itime) = atlas_fieldset(c_bg(itime))
  f_fg(itime) = atlas_fieldset(c_fg(itime))
  call c_f_datetime(c_valid_times(itime), f_valid_times(itime))
enddo

! Call implementation
! -------------------
call f_self%create(f_comm, f_conf, ntimes, f_bg, f_fg, f_valid_times)

deallocate(f_bg, f_fg, f_valid_times)

end subroutine gsi_covariance_create_cpp

! --------------------------------------------------------------------------------------------------

subroutine gsi_covariance_delete_cpp(c_self) &
           bind(c, name='gsi_covariance_delete_f90')

! Arguments
integer(c_int), intent(inout)  :: c_self

! Locals
type(gsi_covariance), pointer :: f_self

! LinkedList
! ----------
call gsi_covariance_registry%get(c_self, f_self)

! Call implementation
! -------------------
call f_self%delete()

! LinkedList
! ----------
call gsi_covariance_registry%remove(c_self)

end subroutine gsi_covariance_delete_cpp

! --------------------------------------------------------------------------------------------------

subroutine gsi_covariance_randomize_cpp(c_self, c_inc) &
           bind(c,name='gsi_covariance_randomize_f90')

implicit none

!Arguments
integer(c_int),     intent(in) :: c_self
type(c_ptr), value, intent(in) :: c_inc

type(gsi_covariance), pointer :: f_self
type(atlas_fieldset)          :: f_inc

! LinkedList
! ----------
call gsi_covariance_registry%get(c_self, f_self)

! Fortran APIs
! ------------
f_inc = atlas_fieldset(c_inc)

! Call implementation
! -------------------
call f_self%randomize(f_inc)

end subroutine gsi_covariance_randomize_cpp

! --------------------------------------------------------------------------------------------------

subroutine gsi_covariance_multiply_cpp(c_self, c_ntimes, c_inc) &
           bind(c,name='gsi_covariance_multiply_f90')

implicit none

!Arguments
integer(c_int),     intent(in) :: c_self
integer(c_int),     intent(in) :: c_ntimes
type(c_ptr),        intent(in) :: c_inc(c_ntimes)

type(gsi_covariance), pointer :: f_self
type(atlas_fieldset), dimension(:), allocatable :: f_inc
integer :: ntimes, itime

! LinkedList
! ----------
call gsi_covariance_registry%get(c_self, f_self)

! Fortran APIs
! ------------
ntimes = c_ntimes
allocate(f_inc(ntimes))
do itime = 1, ntimes
  f_inc(itime) = atlas_fieldset(c_inc(itime))
enddo

! Call implementation
! -------------------
call f_self%multiply(ntimes, f_inc)

deallocate(f_inc)

end subroutine gsi_covariance_multiply_cpp

! --------------------------------------------------------------------------------------------------

end module gsi_covariance_interface_mod
