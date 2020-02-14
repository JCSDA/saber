!----------------------------------------------------------------------
! Module: type_gaugrid
! Purpose: Gaussian grid derived type
! Author: Teppei Kinami
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 019-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
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
  character(len=16),allocatable :: vname(:)            !> Name of variables
  real(kind=kind_real),allocatable :: rlats(:)         !> Gaussian latitudes
  real(kind=kind_real),allocatable :: wlats(:)         !> Gaussian weights 
  real(kind=kind_real),allocatable :: rlons(:)         !> Gaussian longitudes 
  real(kind=kind_real),allocatable :: fld(:,:,:,:)     !> Data
end type gaussian_grid

public :: gaussian_grid
public :: create_gaugrid, delete_gaugrid
public :: gaugrid_alloc_coord, gaugrid_dealloc_coord
public :: gaugrid_alloc_field, gaugrid_dealloc_field
public :: gaugrid_calc_ll_glb

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
  allocate(self%vname(self%nlev))

end subroutine gaugrid_alloc_coord

!----------------------------------------------------------------------
! Subroutine: gaugrid_dealloc_coord
! Purpose: deallocate Gaussian grid coordinate
!----------------------------------------------------------------------
subroutine gaugrid_dealloc_coord(self)
  implicit none
! Passed variables
  class(gaussian_grid),intent(inout) :: self

  if (allocated(self%rlons)) deallocate(self%rlons)
  if (allocated(self%rlats)) deallocate(self%rlats)
  if (allocated(self%wlats)) deallocate(self%wlats)
  if (allocated(self%vname)) deallocate(self%vname)

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
! Subroutine: gaugrid_calc_ll_glb
! Purpose: calculate Gaussian latitudes and longitudes
!----------------------------------------------------------------------
subroutine gaugrid_calc_ll_glb(self)
  implicit none
! Passed variables
  class(gaussian_grid),intent(inout) :: self

! Local variable
  real(kind=kind_real),allocatable :: slat(:)
  integer :: i,j
  real(kind=kind_real) :: dlon

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

  dlon = two*pi/real(self%nlon,kind_real)
  do i=1,self%nlon
    self%rlons(i) = real(i-1,kind_real)*dlon
  end do

end subroutine gaugrid_calc_ll_glb

end module type_gaugrid
