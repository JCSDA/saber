!-------------------------------------------------------------------------------
! (C) Crown Copyright 2022 Met Office.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!-------------------------------------------------------------------------------
! Purpose : in the training data
!-------------------------------------------------------------------------------

module mo_netcdf_mod

use iso_c_binding, only: c_int, c_float
use fckit_log_module, only : fckit_log
use mo_cvtcoord_mod, only: cvt_coordinate_type

! -----------------------------------------------------------------------------

implicit none
private

public :: cvt_nc_err_rpt, cvt_nc_read_field_from_file

! ------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

subroutine cvt_nc_err_rpt(nf, tag)

use netcdf, only:  &
  nf90_strerror,    &
  nf90_noerr

! subroutine arguments:
integer(kind=c_int),          intent(in) :: nf    ! NetCDF error code
character(len=*),             intent(in) :: tag   !character string description

! Local declarations:
character(len=*), PARAMETER    ::    RoutineName="cvt_nc_err_rpt"
character(len=500)             ::    message
external                       ::    abor1_ftn

!-------------------------------------------------------------------------------

! Exit with error code 2 when a NetCDF failure occurs.
if (nf /= nf90_noerr) then
  WRITE(message,'(A)') RoutineName // " " // trim(adjustl(tag)) //" " // &
    trim(adjustl(nf90_strerror(nf) ) )

  call abor1_ftn(message)
end if

!-------------------------------------------------------------------------------

end subroutine cvt_nc_err_rpt


subroutine cvt_nc_read_field_from_file( &
  InFile,          & ! IN
  short_name,      & ! IN
  Field3D,         & ! out optional
  coordinate,      & ! inout optional
  start_index,     & ! out optional
  final_index,     & ! out optional
  long_name,       & ! out optional
  phi_pole,        & ! out optional
  lambda_pole,     & ! out optional
  vrt_grid,        & ! out optional
  hrz_grid,        & ! out optional
  hbc,             & ! out optional
  units,           & ! out optional
  start_bin,       & ! out optional
  final_bin,       & ! out optional
  TimeEnsBinning,  & ! out optional
  grid_mapping,    & ! out optional
  ens_member,      & ! out optional
  valid_time,      & ! out optional
  title,           & ! out optional, GLOBAL ATTRIBUTE
  institution,     & ! out optional, GLOBAL ATTRIBUTE
  source,          & ! out optional, GLOBAL ATTRIBUTE
  history,         & ! out optional, GLOBAL ATTRIBUTE
  references,      & ! out optional, GLOBAL ATTRIBUTE
  comment,         & ! out optional, GLOBAL ATTRIBUTE
  Conventions)       ! out optional, GLOBAL ATTRIBUTE

! Description:
!
!   NOTE THAT THERE IS A LIMIT ON LENGTH OF THE ----- PATH + FILENAME
!
!   The code expects the array bounds of the field are stored in
!      attributes "start_index" and " final_index".
!
!   If the optional arguments "start index" and "final index" are given
!     then  "start_index" and " final_index" will be output from
!     subroutine.
!
!   The code will read the "long_name" from the file and output it if
!     requested by the optional "long_name" call.
!
!   If optional co-ordinate names are given e.g coord_short
!                                      and they exist in the netCDF file
!     then they are read into subroutine.
!
!   They are output if the optional a
!
!
!

use netcdf, only:         &
  nf90_close,             &
  nf90_get_att,           &
  nf90_get_var,           &
  nf90_global,            &              ! Global attribute
  nf90_inquire_attribute, &
  nf90_inquire_dimension, &
  nf90_inquire_variable,  &
  nf90_inq_attname,       &
  nf90_inq_varid,         &
  nf90_max_name,          &
  nf90_noerr,             &
  nf90_nowrite,           &
  nf90_open,              &
  nf90_strerror

! Parameters
integer,                            PARAMETER :: ndims = 3

! subroutine arguments
character(len=nf90_max_name), intent(in) :: InFile ! file name and path

character(len=nf90_max_name), intent(in) :: short_name ! file name and path

real(kind=c_float), allocatable, optional, intent(out) :: Field3D(:,:,:)

type(cvt_coordinate_type), optional, intent(inout):: coordinate(ndims) ! coordinate info

integer(kind=c_int), optional, intent(out) :: start_index(ndims)
integer(kind=c_int), optional, intent(out) :: final_index(ndims)

character(len=nf90_max_name), optional, intent(out) :: long_name ! file name and path
real, optional, intent(out) :: phi_pole
real, optional, intent(out) :: lambda_pole
! NOTE phi_pole and lambda_pole only exist for the variable "rotated_pole"

character(len=nf90_max_name), optional, intent(out) :: vrt_grid
character(len=nf90_max_name), optional, intent(out) :: hrz_grid
character(len=nf90_max_name), optional, intent(out) :: hbc
character(len=nf90_max_name), optional, intent(out) :: units ! Dimensional units
integer, optional, intent(out) :: start_bin ! First bin no
integer, optional, intent(out) :: final_bin ! Last bin no
character(len=nf90_max_name), optional, intent(out) :: TimeEnsBinning
character(len=nf90_max_name), optional, intent(out) :: grid_mapping
integer(kind=c_int)         , optional, intent(out) :: ens_member
character(len=nf90_max_name), optional, intent(out) :: valid_time


! subroutine arguments : Global attributes
character(len=nf90_max_name), optional, intent(out) :: title
character(len=nf90_max_name), optional, intent(out) :: institution
character(len=nf90_max_name), optional, intent(out) :: source
character(len=nf90_max_name), optional, intent(out) :: history
character(len=nf90_max_name), optional, intent(out) :: references
character(len=nf90_max_name), optional, intent(out) :: comment
character(len=nf90_max_name), optional, intent(out) :: Conventions

! Local constants

CHARACTER (LEN=*), PARAMETER :: RoutineName = "cvt_nc_read_field_from_file"

! Local variables
integer(kind=c_int) :: ncid      ! File ID
integer(kind=c_int) :: varid
integer(kind=c_int) :: att_type
integer(kind=c_int) :: att_len
integer(kind=c_int) :: i ! Loop index
integer(kind=c_int) :: nf

integer(kind=c_int) :: start_index_r(ndims)  ! version read in
integer(kind=c_int) :: final_index_r(ndims)

integer :: start_index_l(3)
integer :: final_index_l(3)

integer(kind=c_int) :: numatts
integer(kind=c_int) :: numDims
integer(kind=c_int), allocatable :: dimid(:)
integer(kind=c_int) :: dim_len(ndims)
integer(kind=c_int) :: local          ! Temporary integer variable


character(len=nf90_max_name) :: temp_name
character(len=nf90_max_name) :: dim_name(ndims)

character(len=nf90_max_name) :: valid_time_local
character(len=80) :: ErrorString
character(len=800) :: output_message

!------------------------------------------------------------------------------

temp_name         = ""
dim_name(1)       = ""
dim_name(2)       = ""
dim_name(3)       = ""
valid_time_local  = ""

if (present(long_name)) long_name = ""

write(output_message,*) 'entering ' // RoutineName // " " // trim(short_name)
call fckit_log % debug(output_message)

nf = nf90_open(InFile, nf90_nowrite, ncid)
call cvt_nc_err_rpt(nf, "open file " // trim(InFile))

nf = nf90_inq_varid     (ncid, trim(short_name), varid)
call cvt_nc_err_rpt(nf, "inq_varid short_name" // trim(short_name))

nf = nf90_inquire_variable(ncid, varid, ndims = numDims, natts = numatts)
call cvt_nc_err_rpt(nf, "inquire variable varid ndims natts" // trim(short_name))

ALLOCATE ( dimid(numDims) )

nf = nf90_inquire_variable(ncid, varid, dimids=dimid)
call cvt_nc_err_rpt(nf, "inquire variable dimids" // trim(short_name))

DO i = 1, numDims
  nf = nf90_inquire_dimension(ncid, dimid(i), name = dim_name(i),len = dim_len(i))
  call cvt_nc_err_rpt(nf, "inquire dimension dimid dim_name" // trim(short_name))
end DO

DEALLOCATE ( dimid )

! The presence of ens_member, validity time, start_bin, final_bin,
! TimeEnsBinning varies from input file:
!
if (present(ens_member)) then
  nf = nf90_inquire_attribute(ncid, varid, "ens_member", &
                     xtype = att_type, len = att_len)

  if (nf /= nf90_noerr) then
    ErrorString = nf90_strerror(nf)
    if ((ErrorString(1:27) == "NetCDF: Attribute not found") .OR. &
        (ErrorString(1:19) == "Attribute not found")) then
      ens_member = -99
    else
      ens_member = -1
    end if
  end if
end if

if (present(valid_time)) then
  nf = nf90_inquire_attribute(ncid, varid, "valid_time", &
                     xtype = att_type, len = att_len)

  ErrorString = nf90_strerror(nf)
  if ((ErrorString(1:27) == "NetCDF: Attribute not found") .OR. &
      (ErrorString(1:19) == "Attribute not found"))             &
    valid_time = "NetCDF: Attribute empty"

end if

if (present(start_bin)) then
  nf = nf90_inquire_attribute(ncid, varid, "start_bin", &
                     xtype = att_type, len = att_len)

  ErrorString = nf90_strerror(nf)
  if ((ErrorString(1:27) == "NetCDF: Attribute not found") .OR. &
      (ErrorString(1:19) == "Attribute not found"))             &
    start_bin = -99
end if

if (present(final_bin)) then
  nf = nf90_inquire_attribute(ncid, varid, "final_bin", &
                     xtype = att_type, len = att_len)

  ErrorString = nf90_strerror(nf)
  if ((ErrorString(1:27) == "NetCDF: Attribute not found") .OR. &
      (ErrorString(1:19) == "Attribute not found"))             &
    final_bin = -99
end if

if (present(TimeEnsBinning)) then
  nf = nf90_inquire_attribute(ncid, varid, "TimeEnsBinning", &
                     xtype = att_type, len = att_len)

  ErrorString = nf90_strerror(nf)
  if ((ErrorString(1:27) ==  "NetCDF: Attribute not found") .OR. &
      (ErrorString(1:19) == "Attribute not found"))              &
    TimeEnsBinning = "NetCDF: Attribute empty"
end if

if (present(grid_mapping)) then
  nf = nf90_inquire_attribute(ncid, varid, "grid_mapping", &
                     xtype = att_type, len = att_len)

  ErrorString = nf90_strerror(nf)
  if ((ErrorString(1:27) ==  "NetCDF: Attribute not found") .OR. &
      (ErrorString(1:19) == "Attribute not found"))              &
    grid_mapping = "NetCDF: Attribute empty"
end if


! Read variable attributes. It is important to set CHARACTER variables to the
! empty string prior to reading to keep out unwanted junk.
DO i = 1, numatts

  nf = nf90_inq_attname(ncid, varid, i, temp_name)
  call cvt_nc_err_rpt(nf, "inq_attname" // trim(short_name))

  nf = nf90_inquire_attribute(ncid, varid, temp_name, &
                          xtype = att_type, len = att_len)
  call cvt_nc_err_rpt(nf, "inquire_attribute" // trim(short_name))

  select case (temp_name)
  case ("ens_member")
    if (present(ens_member)) then
      nf = nf90_get_att(ncid, varid, temp_name, local)
      call cvt_nc_err_rpt(nf,"ens_member" // trim(InFile))
      ens_member = local
      if (ens_member < 0) then
        write(output_message,*) RoutineName, &
            "ens_member attribute value is negative."
        call fckit_log % warning(output_message)
      end if
    end if
  case ("valid_time")
    if (present(valid_time)) then
      nf = nf90_get_att(ncid, varid, temp_name, valid_time_local)
      call cvt_nc_err_rpt(nf,"valid_time" // trim(InFile))
      valid_time = valid_time_local
    end if
  case ("start_index")
    if (present(start_index)) then
      nf = nf90_get_att(ncid, varid, temp_name, start_index_r)
      call cvt_nc_err_rpt(nf, "get_att:start_index_r" // trim(short_name))
    end if
  case ("final_index")
    if (present(final_index)) then
      nf = nf90_get_att(ncid, varid, temp_name, final_index_r)
      call cvt_nc_err_rpt(nf, "get_att:final_index_r" // trim(short_name))
    end if
  case ("long_name")
    if (present(long_name)) then
      long_name = ""
      nf = nf90_get_att(ncid, varid, temp_name, long_name)
      call cvt_nc_err_rpt(nf, "get_att:long_name" // trim(short_name))
    end if
  case ("grid_north_pole_latitude")
    if (present(phi_pole)) then
      nf = nf90_get_att(ncid, varid, temp_name, phi_pole)
      call cvt_nc_err_rpt(nf, "get_att:grid_north_pole_latitude" &
        // trim(short_name))
    end if
  case ("grid_north_pole_longitude")
    if (present(lambda_pole)) then
      nf = nf90_get_att(ncid, varid, temp_name, lambda_pole)
      call cvt_nc_err_rpt(nf, "get_att:grid_north_pole_longitude" &
        // trim(short_name))
    end if
  case ("vrt_grid")
    if (present(vrt_grid)) then
      vrt_grid = ""
      nf = nf90_get_att(ncid, varid, temp_name, vrt_grid)
      call cvt_nc_err_rpt(nf, "get_att:vrt_grid" // trim(short_name))
    end if
  case ("hrz_grid")
    if (present(hrz_grid)) then
      hrz_grid = ""
      nf = nf90_get_att(ncid, varid, temp_name, hrz_grid)
      call cvt_nc_err_rpt(nf, "get_att:hrz_grid" // trim(short_name))
    end if
  case ("hbc")
    if (present(hbc)) then
      hbc = ""
      nf = nf90_get_att(ncid, varid, temp_name, hbc)
      call cvt_nc_err_rpt(nf, "get_att:hbc" // trim(short_name))
    end if
  case ("units")
    if (present(units)) then
      units = ""
      nf = nf90_get_att(ncid, varid, temp_name, units)
      call cvt_nc_err_rpt(nf, "get_att:units " // trim(short_name))
    end if
  case ("start_bin")
    if (present(start_bin)) then
      nf = nf90_get_att(ncid, varid, "start_bin", local)
      call cvt_nc_err_rpt(nf, "get_att:start_bin " // trim(short_name))
      start_bin = local
      if (start_bin < 0) then
        write(output_message,*) RoutineName, &
          "start_bin attribute value is negative."
        call fckit_log % warning(output_message)
      end if
    end if
  case ("final_bin")
    if (present(final_bin)) then
      nf = nf90_get_att(ncid, varid, "final_bin", local)
      call cvt_nc_err_rpt(nf, "get_att:final_bin " // trim(short_name))
      final_bin = local
      if (final_bin < 0) then
        write(output_message,*) RoutineName, &
          "final_bin attribute value is negative."
        call fckit_log % warning(output_message)
      end if
    end if
  case ("TimeEnsBinning")
    if (present(TimeEnsBinning)) then
      TimeEnsBinning = ""
      nf = nf90_get_att(ncid, varid, temp_name, TimeEnsBinning)
      call cvt_nc_err_rpt(nf, "get_att:TimeEnsBinning " // trim(short_name))
    end if
  case ("grid_mapping")
    if (present(grid_mapping)) then
      grid_mapping = ""
      nf = nf90_get_att(ncid, varid, temp_name, grid_mapping)
      call cvt_nc_err_rpt(nf, "get_att:grid_mapping " // trim(short_name))
    end if

  case default
    write(output_message,*) RoutineName, ' not reading attribute ', temp_name
    call fckit_log % info(output_message)
  end select
end DO

if (present(start_index))  start_index = start_index_r
if (present(final_index))  final_index = final_index_r

if (present(Field3D)) then

  start_index_l = start_index_r
  final_index_l = final_index_r

  allocate(Field3D(start_index_l(1):final_index_l(1), &
                   start_index_l(2):final_index_l(2), &
                   start_index_l(3):final_index_l(3) ))


  nf = nf90_get_var(ncid, varid, Field3D)
  call cvt_nc_err_rpt(nf, "get_var:Field3D"//trim(short_name))

end if

! Coordinate information
if (present(coordinate)) then

  ! Store the short names.
  coordinate(:) % CoName = dim_name(:)
  DO i = 1, ndims

    ! Store the coordinate values.
    nf = nf90_inq_varid(ncid, trim(dim_name(i)), coordinate(i) % varid)
    call cvt_nc_err_rpt(nf, "inq_varid:coordinate values" // trim(short_name))
    if (allocated(coordinate(i) % CoValues)) deallocate(coordinate(i) % CoValues)
    allocate(coordinate(i) % CoValues(dim_len(i)))
    nf = nf90_get_var(ncid, coordinate(i) % varid, coordinate(i) % CoValues)
    call cvt_nc_err_rpt(nf, "get_var:coordinate values"// trim(short_name))

    ! Store the long_name, units and positive attributes where present.
    nf = nf90_get_att(ncid, coordinate(i) % varid, "long_name", &
      coordinate(i) % long_name)
    if (nf /= nf90_noerr) coordinate(i) % long_name = ""
    nf = nf90_get_att(ncid, coordinate(i) % varid, "units", &
      coordinate(i) % units)
    if (nf /= nf90_noerr) coordinate(i) % units = ""
    nf = nf90_get_att(ncid, coordinate(i) % varid, "positive", &
      coordinate(i) % positive)
    if (nf /= nf90_noerr) coordinate(i) % positive = ""
  end do
end if

! Global attributes

if (present(title)) then
  nf = nf90_get_att(ncid, nf90_global,"title", title)
  call cvt_nc_err_rpt(nf, "global attribute title " // trim(InFile))
end if

if (present(institution)) then
  nf = nf90_get_att(ncid, nf90_global,"institution", institution)
  call cvt_nc_err_rpt(nf,"global attribute institution" // trim(InFile))
end if

if (present(source)) then
  nf = nf90_get_att(ncid, nf90_global,"source", source)
  call  cvt_nc_err_rpt(nf,"global attribute source" // trim(InFile))
end if

if (present(history)) then
  nf = nf90_get_att(ncid, nf90_global,"history", history)
  call cvt_nc_err_rpt(nf,"global attribute history" // trim(InFile))
end if

if (present(references)) then
  nf = nf90_get_att(ncid, nf90_global,"references", references)
  call cvt_nc_err_rpt(nf,"global attribute references" // trim(InFile))
end if

if (present(comment)) then
  nf = nf90_get_att(ncid, nf90_global,"comment", comment)
  call cvt_nc_err_rpt(nf,"global attribute comment" // trim(InFile))
end if

if (present(Conventions)) then
  nf = nf90_get_att(ncid, nf90_global,"Conventions", Conventions)
  call cvt_nc_err_rpt(nf,"global attribute Conventions" // trim(InFile))
end if

! close file
nf = nf90_close(ncid)
call cvt_nc_err_rpt(nf, "close file " // trim(InFile))

write(output_message,*) 'exiting '//RoutineName//' '//trim(short_name)
call fckit_log % debug(output_message)

end subroutine cvt_nc_read_field_from_file

! ------------------------------------------------------------------------------

end module mo_netcdf_mod

