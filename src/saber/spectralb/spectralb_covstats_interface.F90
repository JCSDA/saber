! * (C) Copyright 2022 Met Office UK
! *
! * This software is licensed under the terms of the Apache Licence Version 2.0
! * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!------------------------------------------------------------------------------

subroutine c_covSpectralBins(c_conf, &
 & varname_length, c_netcdfvarname, &
 & bins) &
 & bind(c,name='covSpectralBins_f90')

use iso_c_binding, only : c_ptr, c_int, c_float, c_char
use fckit_configuration_module, only: fckit_configuration
use netcdf, only: nf90_max_name
use kinds
use string_f_c_mod
use spectralb_netcdf_mod, only : cvt_nc_read_field_from_file, &
                                 cvt_nc_err_rpt

implicit none

type(c_ptr), value, intent(in) :: c_conf
integer(c_int),     intent(in) :: varname_length
character(kind=c_char, len=1), intent(in) :: c_netcdfvarname(varname_length+1)
integer(c_int),     intent(inout) :: bins

character(len=nf90_max_name) :: covariance_file
character(800)               :: fieldname
character(len=nf90_max_name) :: short_name

character(len=:), allocatable :: str
type(fckit_configuration) :: f_conf

integer(kind=c_int)  :: start_index(3)
integer(kind=c_int)  :: final_index(3)



!------------------------------------------------------------------------------
! read filename for config
!------------------------------------------------------------------------------
covariance_file  = ""
short_name = ""

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("covariance_file", str)
covariance_file = str

call c_f_string(c_netcdfvarname, fieldname)

short_name = TRIM(ADJUSTL(fieldname))

call cvt_nc_read_field_from_file(covariance_file, &
                                 short_name, &
                                 start_index = start_index(:), &
                                 final_index = final_index(:))

bins = final_index(1) - start_index(1)

end subroutine c_covSpectralBins


subroutine c_covSpectralUMatrix(c_conf, &
 & varname_length, c_netcdfvarname, &
 & model_level, bins, values_size, values) &
 & bind(c,name='covSpectralUMatrix_f90')

use iso_c_binding, only : c_ptr, c_int, c_float, c_char
use fckit_configuration_module, only: fckit_configuration
use netcdf, only: nf90_max_name
use kinds
use string_f_c_mod
use spectralb_netcdf_mod, only : cvt_nc_read_field_from_file, &
                                 cvt_nc_err_rpt

implicit none

type(c_ptr), value, intent(in) :: c_conf
integer(c_int),     intent(in) :: varname_length
character(kind=c_char, len=1), intent(in) :: c_netcdfvarname(varname_length+1)
integer(c_int),     intent(in) :: model_level
integer(c_int),     intent(in) :: bins
integer(c_int),     intent(in) :: values_size
real(c_float),   intent(inout) :: values(values_size)

character(len=nf90_max_name) :: covariance_file
character(800)               :: fieldname
character(len=nf90_max_name) :: short_name

character(len=:), allocatable :: str
type(fckit_configuration) :: f_conf

integer(kind=c_int)  :: start_index(3)
integer(kind=c_int)  :: final_index(3)

real(kind=c_float), allocatable :: Field3D(:,:,:)

integer :: i,j,k,n,b ! loop variables

!------------------------------------------------------------------------------
! read filename for config
!------------------------------------------------------------------------------
covariance_file  = ""
short_name = ""

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("covariance_file", str)
covariance_file = str

call c_f_string(c_netcdfvarname, fieldname)

short_name = TRIM(ADJUSTL(fieldname))

call cvt_nc_read_field_from_file(covariance_file, &
                                 short_name, &
                                 Field3D = Field3D, &
                                 start_index = start_index(:), &
                                 final_index = final_index(:))

n = 1
do b = 0, bins -1
  do j = start_index(2), final_index(2)
    do k = start_index(3), final_index(3)
      i = start_index(1) + b
      values(n) = Field3D(i,j,k)
      n = n + 1
    end do
  end do
end do

end subroutine c_covSpectralUMatrix


