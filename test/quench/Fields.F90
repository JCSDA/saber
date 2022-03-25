! (C) Copyright 2022 UCAR.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module Fields

use atlas_module, only: atlas_field,atlas_fieldset,atlas_functionspace_structuredcolumns,atlas_structuredgrid
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use iso_fortran_env, only: output_unit
use netcdf, only: nf90_create,nf90_clobber,nf90_netcdf4,nf90_noerr,nf90_double,nf90_def_dim,nf90_def_var,nf90_put_var,nf90_close, &
 & nf90_strerror
use oops_variables_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Write fields based on structured columns function space to NetCDF file
subroutine fields_write_structuredcolumns_c(c_conf,c_vars,c_globalCoordinates,c_globalData) &
 & bind(c,name='fields_write_structuredcolumns_f90')

implicit none

! Passed variables
type(c_ptr),value,intent(in) :: c_conf              !< Configuration
type(c_ptr),value,intent(in) :: c_vars              !< List of variables
type(c_ptr),intent(in),value :: c_globalCoordinates !< ATLAS fieldset containing coordinates
type(c_ptr),intent(in),value :: c_globalData        !< ATLAS fieldset containing coordinates

! Local variables
integer :: jvar,i,j,k,gidx,nxmax,ny,nz
integer :: ncid,nx_id,ny_id,nz_id,lon_id,lat_id
integer,allocatable :: field_id(:)
real(c_double) :: lonlat(2)
real(c_double),allocatable :: lon(:,:),lat(:,:),field(:,:,:)
real(c_double),pointer :: ptr_lon(:),ptr_lat(:),ptr(:,:)
character(len=1024) :: fieldname,filepath
character(len=:),allocatable :: str
type(atlas_field) :: afield_lon,afield_lat,afield
type(atlas_functionspace_structuredcolumns) :: afunctionspace
type(atlas_fieldset) :: globalCoordinates,globalData
type(atlas_structuredgrid) :: agrid
type(fckit_configuration) :: f_conf
type(oops_variables) :: vars

! Interface
f_conf = fckit_configuration(c_conf)
vars = oops_variables(c_vars)
globalCoordinates = atlas_fieldset(c_globalCoordinates)
globalData = atlas_fieldset(c_globalData)

! Implementation

! Allocation
allocate(field_id(vars%nvars()))

! Get coordinate field
afield_lon = globalCoordinates%field('lon')
afield_lat = globalCoordinates%field('lat')
afield = globalData%field(vars%variable(1))

! Get grid
afunctionspace = afield_lon%functionspace()
agrid = afunctionspace%grid()

! Get grid size
nxmax = agrid%nxmax()
ny = agrid%ny()
nz = afield%levels()

! Get file name
call f_conf%get_or_die('filepath',str)
filepath = str

! Open file
call ncerr(nf90_create(trim(filepath)//'.nc',ior(nf90_clobber,nf90_netcdf4),ncid))

! Define dimensions
call ncerr(nf90_def_dim(ncid,'nx',nxmax,nx_id))
call ncerr(nf90_def_dim(ncid,'ny',ny,ny_id))
call ncerr(nf90_def_dim(ncid,'nz',nz,nz_id))

! Define coordinates
call ncerr(nf90_def_var(ncid,'lon',nf90_double,(/nx_id,ny_id/),lon_id))
call ncerr(nf90_def_var(ncid,'lat',nf90_double,(/nx_id,ny_id/),lat_id))

do jvar=1,vars%nvars()
   ! Get field
   fieldname = vars%variable(jvar)
   afield = globalData%field(fieldname)

   ! Define variables
   call ncerr(nf90_def_var(ncid,fieldname,nf90_double,(/nx_id,ny_id,nz_id/),field_id(jvar)))
end do

! Allocation
allocate(lon(nxmax,ny))
allocate(lat(nxmax,ny))
allocate(field(nxmax,ny,nz))

! Initialization
lon = 0.0_c_double
lat = 0.0_c_double
field = 0.0_c_double

! Copy coordinates
call afield_lon%data(ptr_lon)
call afield_lat%data(ptr_lat)
do j=1,agrid%ny()
   do i=1,agrid%nx(j)
      gidx = agrid%index(i,j)
      lon(i,j) = ptr_lon(gidx)
      lat(i,j) = ptr_lat(gidx)
   end do
end do

! Write coordinates
call ncerr(nf90_put_var(ncid,lon_id,lon))
call ncerr(nf90_put_var(ncid,lat_id,lat))

do jvar=1,vars%nvars()
   ! Get field
   fieldname = vars%variable(jvar)
   afield = globalData%field(fieldname)

   ! Copy data
   call afield%data(ptr)
   do k=1,nz
      do j=1,agrid%ny()
         do i=1,agrid%nx(j)
            gidx = agrid%index(i,j)
            field(i,j,k) = ptr(k,gidx)
         end do
      end do
   end do

   ! Write data
   call ncerr(nf90_put_var(ncid,field_id(jvar),field))
end do

! Close file
call ncerr(nf90_close(ncid))

! Release memory
deallocate(field_id)
deallocate(lon)
deallocate(lat)
deallocate(field)

end subroutine fields_write_structuredcolumns_c
! ------------------------------------------------------------------------------
subroutine ncerr(info)

implicit none

! Passed variables
integer,intent(in) :: info !< Info index

if (info/=nf90_noerr) then
   write(output_unit,'(a)') '!!! NetCDF error in QUENCH: '//trim(nf90_strerror(info))
   call flush(output_unit)
   stop
end if

end subroutine ncerr
! ------------------------------------------------------------------------------
end module Fields
