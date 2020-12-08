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
function eq(x,y) result(test)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: test

if (repro) then
   test = indist(x,y)
else
   test = .not.(abs(x-y)>0.0)
end if

end function eq

!----------------------------------------------------------------------
! Function: inf
!> Inferior test for reals
!----------------------------------------------------------------------
function inf(x,y) result(test)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: test

test = (x<y)
if (repro) test = test.and.(.not.indist(x,y))

end function inf

!----------------------------------------------------------------------
! Function: infeq
!> Inferior or equal test for reals
!----------------------------------------------------------------------
function infeq(x,y) result(test)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: test

test = inf(x,y).or.eq(x,y)

end function infeq

!----------------------------------------------------------------------
! Function: sup
!> Superior test for reals
!----------------------------------------------------------------------
function sup(x,y) result(test)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: test

test = (x>y)
if (repro) test = test.and.(.not.indist(x,y))

end function sup

!----------------------------------------------------------------------
! Function: supeq
!> Superior or equal test for reals
!----------------------------------------------------------------------
function supeq(x,y) result(test)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: test

test = sup(x,y).or.eq(x,y)

end function supeq

!----------------------------------------------------------------------
! Function: indist
!> Indistiguishability test
!----------------------------------------------------------------------
function indist(x,y) result(test)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: test

test = .false.
if (repro) then
   if ((abs(x)>0.0).or.(abs(y)>0.0)) then
      test = abs(x-y)<rth*(abs(x+y))
   else
      test = .true.
   end if
end if

end function indist

!----------------------------------------------------------------------
! Function: small
!> Small value test
!----------------------------------------------------------------------
function small(x,y) result(test)

implicit none

! Passed variables
real(kind_real),intent(in) :: x !< First real
real(kind_real),intent(in) :: y !< Second real

! Returned variable
logical :: test

test = .false.
if (repro) test = abs(x)<rth*abs(y)

end function small

end module tools_repro
