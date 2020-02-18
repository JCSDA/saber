! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
module type_gaugrid

use tools_kinds, only: kind_real
use tools_const, only: pi, rad2deg
use tools_sp,    only: splat

implicit none

private
real(kind=kind_real):: zero = 0.0_kind_real
real(kind=kind_real):: one  = 1.0_kind_real
real(kind=kind_real):: two  = 2.0_kind_real
real(kind=kind_real):: half = 0.5_kind_real

! Gaussian grid derived type
type gaussian_grid
  integer :: nlat                                      !> Number of longitudes 
  integer :: nlon                                      !> Number of latitudes
  integer :: nlev                                      !> Number of levels
  integer :: nvar                                      !> Number of variables
!  integer :: nts                                       !> Number of timeslots
  character(len=32),allocatable :: vname(:)            !> Name of variables
  real(kind=kind_real),allocatable :: rlats(:)         !> Gaussian latitudes
  real(kind=kind_real),allocatable :: wlats(:)         !> Gaussian weights 
  real(kind=kind_real),allocatable :: rlons(:)         !> Gaussian longitudes 
  real(kind=kind_real),allocatable :: fld(:,:,:,:)     !> Data
contains
  procedure :: create => create_gaugrid
  procedure :: delete => delete_gaugrid 
  procedure :: fld3d_pointer => gaugrid_fld3d_pointer
  procedure :: fld2d_pointer => gaugrid_fld2d_pointer
  procedure :: calc_glb_latlon => gaugrid_calc_glb_latlon
  procedure :: equals
  generic :: assignment(=) => equals 
end type gaussian_grid

public :: gaussian_grid
public :: gaugrid_alloc_coord, gaugrid_dealloc_coord
public :: gaugrid_alloc_field, gaugrid_dealloc_field

contains

!----------------------------------------------------------------------
! Subroutine: create_gaugrid
! Purpose: Create Gaussian grid
!----------------------------------------------------------------------
subroutine create_gaugrid(self)
  implicit none

! Passed variables
  class(gaussian_grid),intent(inout) :: self

! Initialization
  call gaugrid_alloc_coord(self)
  call gaugrid_alloc_field(self)

end subroutine create_gaugrid

!----------------------------------------------------------------------
! Subroutine: delete_gaugrid
! Purpose: Delete Gaussian grid
!----------------------------------------------------------------------
subroutine delete_gaugrid(self)

! Passed variables
  class(gaussian_grid),intent(inout) :: self

  call gaugrid_dealloc_coord(self)
  call gaugrid_dealloc_field(self)

end subroutine delete_gaugrid

!----------------------------------------------------------------------
! Subroutine: gaugrid_alloc_coord
! Purpose: allocate Gaussian grid coordinate
!----------------------------------------------------------------------
subroutine gaugrid_alloc_coord(self)
  implicit none
! Passed variables
  class(gaussian_grid),intent(inout) :: self 

! Initialization
  allocate(self%rlons(self%nlon)); self%rlons=zero
  allocate(self%rlats(self%nlat)); self%rlats=zero
  allocate(self%wlats(self%nlat)); self%wlats=zero
  allocate(self%vname(self%nvar))

end subroutine gaugrid_alloc_coord

!----------------------------------------------------------------------
! Subroutine: gaugrid_dealloc_coord
! Purpose: deallocate Gaussian grid coordinate
!----------------------------------------------------------------------
subroutine gaugrid_dealloc_coord(self)
  implicit none
! Passed variables
  class(gaussian_grid),intent(inout) :: self

  if (allocated(self%rlons))    deallocate(self%rlons)
  if (allocated(self%rlats))    deallocate(self%rlats)
  if (allocated(self%wlats))    deallocate(self%wlats)
  if (allocated(self%vname))    deallocate(self%vname)

end subroutine gaugrid_dealloc_coord

!----------------------------------------------------------------------
! Subroutine: gaugrid_alloc_field
! Purpose: allocate Gaussian grid field
!----------------------------------------------------------------------
subroutine gaugrid_alloc_field(self)
  implicit none
! Passed variables
  class(gaussian_grid),intent(inout) :: self 

! Initialization
  allocate(self%fld(self%nlat,self%nlon,self%nlev,self%nvar))

end subroutine gaugrid_alloc_field

!----------------------------------------------------------------------
! Subroutine: gaugrid_dealloc_field
! Purpose: deallocate Gaussian grid field
!----------------------------------------------------------------------
subroutine gaugrid_dealloc_field(self)
  implicit none
! Passed variables
  class(gaussian_grid),intent(inout) :: self 

  if (allocated(self%fld)) deallocate(self%fld)

end subroutine gaugrid_dealloc_field

!----------------------------------------------------------------------
! Subroutine: gaugrid_calc_glb_latlon
! Purpose: calculate global Gaussian latitudes and longitudes
!----------------------------------------------------------------------
subroutine gaugrid_calc_glb_latlon(self)
  implicit none
! Passed variables
  class(gaussian_grid),intent(inout) :: self

! Local variable
  real(kind=kind_real),allocatable :: slat(:)
  integer :: i,j
  real(kind=kind_real) :: dlon

  ! Gaussian latitudes
  allocate(slat(self%nlat-2))
  call splat(4,self%nlat-2,slat,self%wlats)
  self%rlats = zero
  self%rlats(1) =  -pi*half
  self%rlats(self%nlat) =  pi*half
  do j=1,(self%nlat-2)/2
    self%rlats(self%nlat-j) = asin(slat(j))
    self%rlats(1+j) = -asin(slat(j))
  end do
  deallocate(slat)

  ! equally-spaced longitudes
  dlon = two*pi/real(self%nlon,kind_real)
  do i=1,self%nlon
    self%rlons(i) = real(i-1,kind_real)*dlon
  end do

end subroutine gaugrid_calc_glb_latlon

!----------------------------------------------------------------------
! Subroutine: gaugrid_fld3d_pointer
! Purpose: Set 3D field pointer
!----------------------------------------------------------------------
subroutine gaugrid_fld3d_pointer(self,iv,var,fldpointer)
  implicit none
! Passed variables
  class(gaussian_grid),target, intent(inout) :: self
  integer,                     intent(in)    :: iv
  character(len=*),            intent(in)    :: var
  real(kind_real),     pointer,intent(inout) :: fldpointer(:,:,:)

  self%vname(iv) = trim(var)
  fldpointer => self%fld(1:self%nlat,1:self%nlon,1:self%nlev,iv)

end subroutine gaugrid_fld3d_pointer

!----------------------------------------------------------------------
! Subroutine: gaugrid_fld2d_pointer
! Purpose: Set 2D field pointer
!----------------------------------------------------------------------
subroutine gaugrid_fld2d_pointer(self,iv,var,fldpointer)
  implicit none
! Passed variables
  class(gaussian_grid),target, intent(inout) :: self
  integer,                     intent(in)    :: iv
  character(len=*),            intent(in)    :: var
  real(kind_real),     pointer,intent(inout) :: fldpointer(:,:)

  self%vname(iv) = trim(var)
  fldpointer => self%fld(1:self%nlat,1:self%nlon,1,iv)

end subroutine gaugrid_fld2d_pointer

!----------------------------------------------------------------------
! Subroutine: equals
! Purpose: create new gaussian grid from other
!----------------------------------------------------------------------
subroutine equals(self,rhs)
  implicit none
! Passed variables
  class(gaussian_grid), intent(inout) :: self
  type (gaussian_grid), intent(in)    :: rhs

  call delete_gaugrid(self)
  self%nlat = rhs%nlat; self%nlon = rhs%nlon
  self%nlev = rhs%nlev; self%nvar = rhs%nvar
  call create_gaugrid(self)
  self%vname = rhs%vname
  self%rlats = rhs%rlats; self%rlons = rhs%rlons
  self%wlats = rhs%wlats
  self%fld = rhs%fld

end subroutine equals

end module type_gaugrid
