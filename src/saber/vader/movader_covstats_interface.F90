! * (C) Copyright 2020 Met Office UK
! *
! * This software is licensed under the terms of the Apache Licence Version 2.0
! * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!------------------------------------------------------------------------------

subroutine c_covRegressionMatrices(filename_length, &
 & c_filename, &
 & model_level, bins, values_size, values) &
 & bind(c,name='covRegressionMatrices_f90')

use iso_c_binding, only : c_ptr, c_int, c_float, c_char
use mo_netcdf_mod, only : cvt_nc_read_field_from_file, cvt_nc_err_rpt
use netcdf, only: nf90_max_name
use string_f_c_mod

implicit none

integer(c_int),                intent(in) :: filename_length
character(kind=c_char, len=1), intent(in) :: c_filename(filename_length+1)
integer(c_int),     intent(in) :: model_level
integer(c_int),     intent(in) :: bins
integer(c_int),     intent(in) :: values_size
real(c_float),   intent(inout) :: values(values_size)

character(len=nf90_max_name) :: covariance_file
character(len=nf90_max_name) :: short_name

character(len=:), allocatable :: str

integer(kind=c_int)  :: start_index(3)
integer(kind=c_int)  :: final_index(3)
real(kind=c_float), allocatable :: Field3D(:,:,:)

integer :: i,j,k,n ! loop variables

! read filename for config
call c_f_string(c_filename, covariance_file)
short_name = "RegressionFilter_hP_inc_gP_inc"

call cvt_nc_read_field_from_file(covariance_file, &
                                 short_name, &
                                 Field3D = Field3D, &
                                 start_index = start_index(:), &
                                 final_index = final_index(:))

! maybe apply checks on sizes here.

n = 1
do i = 1, final_index(1)
  do j = 1, final_index(2)
    do k = 1, final_index(3)
      values(n) = Field3D(i,j,k)
      n = n + 1
    end do
  end do
end do


end subroutine c_covRegressionMatrices

!------------------------------------------------------------------------------

subroutine c_covRegressionWeights(filename_length, c_filename, &
 & covGlobalNlats, nBins, weightSize, &
 & startVec, lenVec, covLatitudesVec, regressionWeights1D) &
 & bind(c,name='covRegressionWeights_f90')

use iso_c_binding, only: c_ptr, c_int, c_float, c_char
use mo_cvtcoord_mod, only: cvt_coordinate_type, cvt_initialiseadjustordealloccoord, cvt_create3dcoordinate
use mo_netcdf_mod, only: cvt_nc_read_field_from_file, cvt_nc_err_rpt
use netcdf, only: nf90_max_name
use string_f_c_mod

implicit none

integer(c_int), intent(in) :: filename_length
character(kind=c_char, len=1), intent(in) :: c_filename(filename_length+1)
integer(kind=c_int), intent(in) :: covGlobalNlats
integer(kind=c_int), intent(in) :: nBins
integer(kind=c_int), intent(in) :: weightSize
integer(kind=c_int), intent(inout) :: startVec(nBins)
integer(kind=c_int), intent(inout) :: lenVec(nBins)
real(kind=c_float), intent(inout) :: covLatitudesVec(covGlobalNlats)
real(kind=c_float), intent(inout) :: regressionWeights1D(weightSize)

character(len=nf90_max_name) :: covariance_file
character(len=nf90_max_name) :: short_name

character(len=:), allocatable :: str

integer(kind=c_int)   :: start_index(3)
integer(kind=c_int)   :: final_index(3)
real(kind=c_float), allocatable :: Field3D(:,:,:)

type(cvt_coordinate_type) :: coordinate(3)
type(cvt_coordinate_type) :: coord1
type(cvt_coordinate_type) :: coord2
type(cvt_coordinate_type) :: coord3

integer :: i,j,k,n ! loop variables
integer :: tally

integer :: sh(3)

! read filename for config
covariance_file  = ""
short_name = ""

! read filename for config
call c_f_string(c_filename, covariance_file)

short_name = "WEIGHT_rho_rho"

! 1.1 Initialise 3D coordinate type
CALL cvt_initialiseadjustordealloccoord(coord1, 1)
CALL cvt_initialiseadjustordealloccoord(coord2, 1)
CALL cvt_initialiseadjustordealloccoord(coord3, 1)
CALL cvt_create3dcoordinate(coord1, coord2, coord3, coordinate)

call cvt_nc_read_field_from_file(covariance_file, &
                                 short_name, &
                                 Field3D = Field3D, &
                                 coordinate = coordinate, &
                                 start_index = start_index(:), &
                                 final_index = final_index(:))

do i = 1, 3
! Note that the cov file was generated on a NewDynamics grid
! while we are using an EndGame grid rho is located at different points to
! what is in the cov file!
  if (coordinate(i)%CoName(1:12) == "latitude_rho") then
    covLatitudesVec = coordinate(i)%CoValues
  end if
end do

k = 1
n = 1
do i = 1, nBins
  tally = 0
  do j = 1, covGlobalNlats
    if (Field3D(i,k,j) > 1e-24) then
! Note that unlike the um-jedi version we do not rescale the
! gp interpolation for resolution changes by a constant value here.
! real(covGlobalNlons)/real(globalNlons)
      regressionWeights1D(n) = Field3D(i,k,j)
      tally = tally + 1
      startVec(i) = j
      n = n + 1
    end if
  end do
  lenVec(i) = tally
! we are setting startVec for C++ which indexes from 0 not 1
! so instead of having startVec(i) = startVec(i) - tally + 1
! we have:
  startVec(i) = startVec(i) - tally
end do

end subroutine c_covRegressionWeights

!------------------------------------------------------------------------------

subroutine c_covMuStats(filename_length, &
 & c_filename, &
 & shortname_length, &
 & c_shortname, &
 & modelLevels, muBins, sizeVec, mustats) &
 & bind(c,name='covMuStats_f90')

use iso_c_binding, only: c_ptr, c_int, c_float, c_char
use mo_netcdf_mod, only : cvt_nc_read_field_from_file
use netcdf, only: nf90_max_name
use string_f_c_mod

implicit none

integer(c_int),                intent(in) :: filename_length
character(kind=c_char, len=1), intent(in) :: c_filename(filename_length+1)
integer(c_int),                intent(in) :: shortname_length
character(kind=c_char, len=1), intent(in) :: c_shortname(shortname_length+1)
integer(kind=c_int), intent(in) :: modelLevels
integer(kind=c_int), intent(in) :: muBins
integer(kind=c_int), intent(in) :: sizeVec
real(kind=c_float),  intent(inout) :: mustats(sizeVec)

character(len=nf90_max_name) :: covariance_file
character(len=nf90_max_name) :: short_name

integer(kind=c_int)   :: start_index(3)
integer(kind=c_int)   :: final_index(3)
real(kind=c_float), allocatable :: Field3D(:,:,:)

integer :: i,j,n ! loop variables

! read filename for config
call c_f_string(c_filename, covariance_file)
call c_f_string(c_shortname, short_name)

! note start_index and final_index needed to be called
! to get proper value of Field3D
call cvt_nc_read_field_from_file(covariance_file, &
                                 short_name, &
                                 Field3D = Field3D, &
                                 start_index = start_index(:), &
                                 final_index = final_index(:))

n = 1
do i = 1, modelLevels
  do j = 1, muBins
    mustats(n) = Field3D(j,i,1)
    n = n + 1
  end do
end do

end subroutine c_covMuStats
