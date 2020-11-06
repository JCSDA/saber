!----------------------------------------------------------------------
! Module: tools_repro
!> Reproducibility functions
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_repro

use tools_const, only: pi
use tools_kinds, only: kind_real
use type_mpl, only: mpl_type

implicit none

logical :: repro = .true.                  !< Reproducibility flag
real(kind_real),parameter :: rth = 1.0e-12 !< Reproducibility threshold

private
public :: repro,rth
public :: eq,inf,infeq,sup,supeq,indist,small

contains

!----------------------------------------------------------------------
! Function: eq
!> Equal test for reals
!----------------------------------------------------------------------
function eq(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: eq

if (repro) then
   eq = indist(x,y)
else
   eq = .not.(abs(x-y)>0.0)
end if

end function eq

!----------------------------------------------------------------------
! Function: inf
!> Inferior test for reals
!----------------------------------------------------------------------
function inf(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: inf

inf = (x<y)
if (repro) inf = inf.and.(.not.indist(x,y))

end function inf

!----------------------------------------------------------------------
! Function: infeq
!> Inferior or equal test for reals
!----------------------------------------------------------------------
function infeq(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: infeq

infeq = inf(x,y).or.eq(x,y)

end function infeq

!----------------------------------------------------------------------
! Function: sup
!> Superior test for reals
!----------------------------------------------------------------------
function sup(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: sup

sup = (x>y)
if (repro) sup = sup.and.(.not.indist(x,y))

end function sup

!----------------------------------------------------------------------
! Function: supeq
!> Superior or equal test for reals
!----------------------------------------------------------------------
function supeq(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: supeq

supeq = sup(x,y).or.eq(x,y)

end function supeq

!----------------------------------------------------------------------
! Function: indist
!> Indistiguishability test
!----------------------------------------------------------------------
function indist(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: indist

indist = .false.
if (repro) then
   if ((abs(x)>0.0).or.(abs(y)>0.0)) then
      indist = abs(x-y)<rth*(abs(x+y))
   else
      indist = .true.
   end if
end if

end function indist

!----------------------------------------------------------------------
! Function: small
!> Small value test
!----------------------------------------------------------------------
function small(x,y)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: small

small = .false.
if (repro) small = abs(x)<rth*abs(y)

end function small

end module tools_repro
