!----------------------------------------------------------------------
! Module: type_gaugrid
! Purpose: Gaussian grid derived type
! Author: Teppei Kinami
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
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
  integer :: nlat                                      !> Number of longitudes (global) 
  integer :: nlon                                      !> Number of latitudes (global)
  integer :: nlev                                      !> Number of levels (global)
  integer :: nvar                                      !> Number of variables
!  integer :: nts                                       !> Number of timeslots
  character(len=16),allocatable :: vname(:)            !> Name of variables
!  integer :: layout(2)                                 !> Layout
!  integer :: hx                                        ! Halo size of x-dim
!  integer :: hy                                        ! Halo size of y-dim
  integer :: lat2                                      !> Number of longitudes
  integer :: lon2                                      !> Number of latitudes
  integer :: lev2                                      !> Number of levels
  real(kind=kind_real),allocatable :: rlats_glb(:)     !> Gaussian latitudes (global)
  real(kind=kind_real),allocatable :: wlats_glb(:)     !> Gaussian weights (global)
  real(kind=kind_real),allocatable :: rlons_glb(:)     !> Gaussian longitudes (global)
  real(kind=kind_real),allocatable :: fld_glb(:,:,:,:) !> Gaussian grid field (global)
  real(kind=kind_real),allocatable :: rlats(:)         !> Gaussian latitudes in subset
  real(kind=kind_real),allocatable :: wlats(:)         !> Gaussian weights in subset
  real(kind=kind_real),allocatable :: rlons(:)         !> Gaussian longitudes in subset
  real(kind=kind_real),allocatable :: fld(:,:,:,:)     !> Gaussian grid field in subset
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
subroutine create_gaugrid(self,nlat,nlon,nlev,nvar,lat2,lon2,lev2)
  implicit none

! Passed variables
  class(gaussian_grid),intent(inout) :: self
  integer,intent(in) :: nlon,nlat,nlev,nvar,lat2,lon2,lev2

! Initialization
  self%nlon = nlon
  self%nlat = nlat
  self%nlev = nlev
  self%nvar = nvar
  self%lon2 = lon2
  self%lat2 = lat2
  self%lev2 = lev2

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
  allocate(self%rlons_glb(self%nlon))
  allocate(self%rlats_glb(self%nlat))
  allocate(self%wlats_glb(self%nlat))
  allocate(self%vname(self%nlev))

  allocate(self%rlons(self%lon2))
  allocate(self%rlats(self%lat2))
  allocate(self%wlats(self%lat2))

end subroutine gaugrid_alloc_coord

!----------------------------------------------------------------------
! Subroutine: gaugrid_dealloc_coord
! Purpose: deallocate Gaussian grid coordinate
!----------------------------------------------------------------------
subroutine gaugrid_dealloc_coord(self)
  implicit none
! Passed variables
  class(gaussian_grid),intent(inout) :: self

  if (allocated(self%rlons_glb)) deallocate(self%rlons_glb)
  if (allocated(self%rlats_glb)) deallocate(self%rlats_glb)
  if (allocated(self%wlats_glb)) deallocate(self%wlats_glb)
  if (allocated(self%vname))     deallocate(self%vname)

  if (allocated(self%rlons))   deallocate(self%rlons)
  if (allocated(self%rlats))   deallocate(self%rlats)
  if (allocated(self%wlats))   deallocate(self%wlats)

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
  allocate(self%fld_glb(self%nlat,self%nlon,self%nlev,self%nvar))
  allocate(self%fld(self%lat2,self%lon2,self%lev2,self%nvar))

end subroutine gaugrid_alloc_field

!----------------------------------------------------------------------
! Subroutine: gaugrid_dealloc_field
! Purpose: deallocate Gaussian grid field
!----------------------------------------------------------------------
subroutine gaugrid_dealloc_field(self)
  implicit none
! Passed variables
  class(gaussian_grid),intent(inout) :: self 

  if (allocated(self%fld_glb)) deallocate(self%fld_glb)
  if (allocated(self%fld))     deallocate(self%fld)

end subroutine gaugrid_dealloc_field

!----------------------------------------------------------------------
! Subroutine: gaugrid_calc_ll_glb
! Purpose: calculate Gaussian latitudes/longitudes
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
  call splat(4,self%nlat-2,slat,self%wlats_glb)
  self%rlats_glb = zero
  self%rlats_glb(1) =  -pi*half
  self%rlats_glb(self%nlat) =  pi*half
  do j=1,(self%nlat-2)/2
    self%rlats_glb(self%nlat-j) = asin(slat(j))
    self%rlats_glb(1+j) = -asin(slat(j))
  end do
  deallocate(slat)

  dlon = two*pi/real(self%nlon,kind_real)
  do i=1,self%nlon
    self%rlons_glb(i) = real(i-1,kind_real)*dlon
  end do

end subroutine gaugrid_calc_ll_glb

end module type_gaugrid
