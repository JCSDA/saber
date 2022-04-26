! (C) Copyright 2022 United States Government as represented by the Administrator of the National
!     Aeronautics and Space Administration
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module gsi_grid_mod

! netcdf
use netcdf

! atlas
use atlas_module,                   only: atlas_field, atlas_fieldset, atlas_real

! fckit
use fckit_mpi_module,               only: fckit_mpi_comm
use fckit_configuration_module,     only: fckit_configuration

! oops
use kinds,                          only: kind_real

! saber
use gsi_utils_mod,                  only: nccheck

! gsibclim
use m_gsibclim,                     only: gsibclim_init

implicit none
private
public gsi_grid

! Fortran class header
type :: gsi_grid
  type(fckit_mpi_comm) :: comm
  character(len=2055) :: filename
  integer :: npx, npy, npz          ! Grid points in global grid
  integer :: layout(2)              ! Number of processors in x (index 1) and y (index 2) directions
  integer :: isc, iec, jsc, jec     ! Start and ending grid points for each processor
  real(kind=kind_real), allocatable :: lats(:), lons(:)
  real(kind=kind_real), allocatable :: grid_lats(:,:), grid_lons(:,:)
  integer :: ngrid ! Number of grid points for each processor
  logical :: debug
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: print
    procedure, public :: get_levels
    procedure, public :: set_atlas_lonlat
end type gsi_grid

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, conf, comm)

! Arguments
class(gsi_grid),           intent(inout) :: self
type(fckit_configuration), intent(in)    :: conf
type(fckit_mpi_comm),      intent(in)    :: comm

! Locals
integer :: ncid, dimid(3), varid(2)
character(len=:), allocatable :: str
integer :: posx, posy, i, j, npx_per_proc, npy_per_proc
logical :: cv

! Create copy of comm
! -------------------
self%comm = comm

! Debug mode
! ----------
call conf%get_or_die("debugging mode", self%debug)

! Read the GSI grid info from file
! --------------------------------
if (comm%rank() == 0) then

  ! Get filename
  call conf%get_or_die("gsi error covariance file", str)
  self%filename = str

  ! Open NetCDF
  call nccheck(nf90_open(trim(self%filename), NF90_NOWRITE, ncid), "nf90_open "//trim(self%filename))

  ! Get grid dimension from file
  call nccheck(nf90_inq_dimid(ncid, "lon", dimid(1)), "nf90_inq_dimid lon")
  call nccheck(nf90_inq_dimid(ncid, "lat", dimid(2)), "nf90_inq_dimid lat")
  call nccheck(nf90_inq_dimid(ncid, "lev", dimid(3)), "nf90_inq_dimid lev")

  call nccheck(nf90_inquire_dimension(ncid, dimid(1), len=self%npx), "nf90_inquire_dimension lon" )
  call nccheck(nf90_inquire_dimension(ncid, dimid(2), len=self%npy), "nf90_inquire_dimension lat" )
  call nccheck(nf90_inquire_dimension(ncid, dimid(3), len=self%npz), "nf90_inquire_dimension lev" )

endif


! Broadcast the dimensions
! ------------------------
call comm%broadcast(self%npx, 0)
call comm%broadcast(self%npy, 0)
call comm%broadcast(self%npz, 0)


! Allocate the lat/lon arrays
! ---------------------------
allocate(self%lons(self%npx))
allocate(self%lats(self%npy))


! Read the latitude and longitude
! -------------------------------
if (comm%rank() == 0) then

  call nccheck(nf90_inq_varid(ncid, "lon", varid(1)), "nf90_inq_varid lon")
  call nccheck(nf90_inq_varid(ncid, "lat", varid(2)), "nf90_inq_varid lat")

  call nccheck(nf90_get_var(ncid, varid(1), self%lons), "nf90_get_var lon" )
  call nccheck(nf90_get_var(ncid, varid(2), self%lats), "nf90_get_var lat" )

  ! Close NetCDF
  call nccheck(nf90_close(ncid), "nf90_close")

end if


! Broadcast the grid
! ------------------
call comm%broadcast(self%lons, 0)
call comm%broadcast(self%lats, 0)


! Domain decomposition
! --------------------
call conf%get_or_die("procesor layout x direction", self%layout(1))
call conf%get_or_die("procesor layout y direction", self%layout(2))

! Check that user choices match comm size
if (.not. self%layout(1)*self%layout(2) == comm%size()) &
  call abor1_ftn("GSI grid: number of processor in layout does not match number in communicator")

! Grid point per processor in each direction
npx_per_proc = floor(real(self%npx, kind_real)/real(self%layout(1), kind_real))
npy_per_proc = floor(real(self%npy, kind_real)/real(self%layout(2), kind_real))

! Start and end points in the x direction
posx = mod(comm%rank(), self%layout(1))
self%isc = posx * npx_per_proc + 1
self%iec = self%isc + npx_per_proc - 1
if (posx == self%layout(1) -1) self%iec = self%npx

! Start and end points in the y direction
do j = 0, self%layout(2)-1
  if ( comm%rank() >= j*self%layout(1) .and. comm%rank() < (j+1)*self%layout(1) ) then
    posy = j
  end if
end do
self%jsc = posy * npy_per_proc + 1
self%jec = self%jsc + npy_per_proc - 1
if (posy == self%layout(2) -1) self%jec = self%npy

self%ngrid = (self%iec-self%isc+1)*(self%jec-self%jsc+1)

! Create arrays of lon/lat to be compatible with interpolation
allocate(self%grid_lons(self%isc:self%iec, self%jsc:self%jec))
allocate(self%grid_lats(self%isc:self%iec, self%jsc:self%jec))

do i = self%isc, self%iec
  self%grid_lons(i,:) = self%lons(i)
enddo
do j = self%jsc, self%jec
  self%grid_lats(:,j) = self%lats(j)
enddo

! This is a test
print*, "Hello gsiblim"
call gsibclim_init(cv)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(gsi_grid), intent(inout) :: self

! Deallocate arrays
deallocate(self%lons)
deallocate(self%lats)
deallocate(self%grid_lons)
deallocate(self%grid_lats)

! Set grid to zero
self%npx = 0
self%npy = 0
self%npz = 0
self%layout = 0
self%isc = 0
self%iec = 0
self%jsc = 0
self%jec = 0

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine print(self)

! Arguments
class(gsi_grid), intent(in) :: self

! Root PE prints grid info
if (self%comm%rank() == 0) then

  write(*,'(A28)')      "ErrorCovarianceGSI GSI grid:"
  write(*,'(A38, I5)')  "  Number of longitudinal grid points: ", self%npx
  write(*,'(A37, I5)')  "  Number of latitudinal grid points: ", self%npy
  write(*,'(A34, I5)')  "  Number of vertical grid points: ", self%npz
  write(*,'(A1)')       " "
  write(*,'(A43, I5)')  "  Number of processors in the x direction: ", self%layout(1)
  write(*,'(A43, I5)')  "  Number of processors in the y direction: ", self%layout(2)

endif

if (self%debug) then
  ! Print index ranges
  write(*,'(A7, I6, A7, I6, A7, I6, A7, I6, A7, I6)')  "  Proc ", self%comm%rank(), &
                        ' isc = ', self%isc, ' iec = ', self%iec, &
                        ' jsc = ', self%jsc, ' jec = ', self%jec

  ! Print latlon
  write(*,'(A10, F10.3, A10, F10.3, A10, F10.3, A10, F10.3)')  &
        "  Lat min ", minval(self%grid_lats), &
        "  Lat max ", maxval(self%grid_lats), &
        "  Lon min ", minval(self%grid_lons), &
        "  Lon max ", maxval(self%grid_lons)
endif

end subroutine print

! --------------------------------------------------------------------------------------------------

subroutine get_levels(self, levels)

! Arguments
class(gsi_grid), intent(in)    :: self
integer,         intent(inout) :: levels

! Get number of levels
! --------------------
levels = self%npz

end subroutine get_levels

! --------------------------------------------------------------------------------------------------

subroutine set_atlas_lonlat(self, grid_fieldset)

!Arguments
class(gsi_grid),      intent(inout) :: self
type(atlas_fieldset), intent(inout) :: grid_fieldset

!Locals
real(kind_real), pointer :: real_ptr(:,:)
type(atlas_field) :: lonlat_field

! Create lon/lat field
lonlat_field = atlas_field(name="lonlat", kind=atlas_real(kind_real), shape=(/2,self%ngrid/))

! Get pointer to the data
call lonlat_field%data(real_ptr)

! Fill lon/lat
real_ptr(1,:) = reshape(self%grid_lons(self%isc:self%iec, &
                                            self%jsc:self%jec), (/self%ngrid/))
real_ptr(2,:) = reshape(self%grid_lats(self%isc:self%iec, &
                                            self%jsc:self%jec), (/self%ngrid/))

! Add field to fieldset
call grid_fieldset%add(lonlat_field)

end subroutine set_atlas_lonlat

! --------------------------------------------------------------------------------------------------

end module gsi_grid_mod
