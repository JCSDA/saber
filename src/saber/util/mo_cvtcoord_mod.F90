!-------------------------------------------------------------------------------
! (C) Crown Copyright 2022 Met Office.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------
! Purpose : in the training data
!-------------------------------------------------------------------------------

module mo_cvtcoord_mod

use iso_c_binding, only : c_int, c_int32_t, c_float
use netcdf, only: nf90_max_name

! -----------------------------------------------------------------------------

implicit none

! ---------------------------------------------------------------------------
! Derived type: NetCDF coordinate attributes
! ---------------------------------------------------------------------------
type cvt_coordinate_type
  character(len=nf90_max_name)  :: CoName
  character(len=nf90_max_name)  :: long_name
  character(len=nf90_max_name)  :: units
  character(len=nf90_max_name)  :: positive
  integer                       :: CoSize
  integer(kind=c_int32_t)       :: dimid
  integer(kind=c_int)           :: varid
  real                          :: origin
  real                          :: delta
  real, allocatable             :: CoValues(:)
end type cvt_coordinate_type

private
public :: cvt_coordinate_type
public :: cvt_create3dcoordinate
public :: cvt_initialiseadjustordealloccoord
! ------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

subroutine cvt_create3dcoordinate( FirCoord1D, SecCoord1D, ThiCoord1D,  &
  CoordOut3D )

! Subroutine arguments:
type(cvt_coordinate_type),             intent(in)  :: FirCoord1D
type(cvt_coordinate_type),             intent(in)  :: SecCoord1D
type(cvt_coordinate_type),             intent(in)  :: ThiCoord1D
type(cvt_coordinate_type),             intent(out) :: CoordOut3D(3)

!------------------------------------------------------------------------------

allocate( CoordOut3D(1) % CoValues(size(FirCoord1D % CoValues) ) )
allocate( CoordOut3D(2) % CoValues(size(SecCoord1D % CoValues) ) )
allocate( CoordOut3D(3) % CoValues(size(ThiCoord1D % CoValues) ) )

CoordOut3D(1) = FirCoord1D
CoordOut3D(2) = SecCoord1D
CoordOut3D(3) = ThiCoord1D

!-------------------------------------------------------------------------------

end subroutine cvt_create3dcoordinate

!------------------------------------------------------------------------------

subroutine cvt_initialiseadjustordealloccoord( CoordOut, &  ! OUT
  CoSize,                                                      &  ! in
  CoordTemplate, CoNa, LongName, CoValues, undo )                 ! in optional

! Subroutine arguments:
type(cvt_coordinate_type),           intent(inout) :: CoordOut
integer(kind=c_int),                    intent(in) :: CoSize

type(cvt_coordinate_type),    intent(in), optional :: CoordTemplate
character(len=nf90_max_name), intent(in), optional :: CoNa
character(len=nf90_max_name), intent(in), optional :: LongName
real(kind=c_float)          , intent(in), optional :: CoValues(CoSize)
logical                     , intent(in), optional :: undo

! Local declarations:
logical :: undo_set ! true if undo was present and .TRUE.

!------------------------------------------------------------------------------
! 1 deallocate structures if undo_set is true
!------------------------------------------------------------------------------
if (present(undo)) then
  undo_set = undo
else
  undo_set = .false.
end if

if (undo_set) then
  deallocate( CoordOut % CoValues )
else

  !------------------------------------------------------------------------------
  ! 2 Setup CoordOut
  !------------------------------------------------------------------------------

  ! 2.1 Copy CoordTemplate structure into CoordOut
  if ( present(CoordTemplate) ) then

    !
    !- Copy CoordTemplate structure into CoordOut
    !
    CoordOut % CoName        = CoordTemplate % CoName
    CoordOut % long_name     = CoordTemplate % long_name
    CoordOut % units         = CoordTemplate % units
    CoordOut % positive      = CoordTemplate % positive
    CoordOut % CoSize        = CoordTemplate % CoSize
    CoordOut % dimid         = CoordTemplate % dimid
    CoordOut % varid         = CoordTemplate % varid
    CoordOut % origin        = CoordTemplate % origin
    CoordOut % delta         = CoordTemplate % delta

    IF (allocated(CoordOut % CoValues)) deallocate(CoordOut % CoValues)

    allocate( CoordOut % CoValues(1:size(CoordTemplate % CoValues) ) )
    CoordOut % CoValues(:) = CoordTemplate % CoValues(:)


  else
    ! 2.2 Set Default CoordOut structure when CoordTemplate not present

    CoordOut % CoName        = ''
    CoordOut % long_name     = ''

    ! default unit is dimensonless which is denoted by '1'
    CoordOut % units         = '1'
    CoordOut % positive      = ''
    CoordOut % CoSize        = CoSize
    CoordOut % dimid         = 0
    CoordOut % varid         = 0
    CoordOut % origin        = 0.0
    CoordOut % delta         = 0.0

    if (allocateD(CoordOut % CoValues))  deallocate(CoordOut % CoValues)

    allocate( CoordOut % CoValues(1:CoSize) )
    CoordOut % CoValues(:)    = 0.0

  end if

  !------------------------------------------------------------------------------
  ! 3  Append where necessary where information is present
  !------------------------------------------------------------------------------
  if (present(CoNa))     CoordOut % CoName     = CoNa
  if (present(LongName)) CoordOut % long_name  = LongName

  if (present(CoValues) ) then

    if ( size(CoordOut % CoValues) /= size(CoValues) ) then

      deallocate( CoordOut % CoValues )
      allocate( CoordOut % CoValues(1:size(CoValues)) )
      CoordOut % CoSize = size(CoValues)

    end if
    CoordOut % CoValues(:) = CoValues(:)

  end if

end if

!-------------------------------------------------------------------------------

end subroutine cvt_initialiseadjustordealloccoord

end module mo_cvtcoord_mod
