!----------------------------------------------------------------------
! Module: type_rng
!> Random numbers generator derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_rng

use iso_fortran_env, only: int64
use tools_kinds, only: kind_real
use type_mpl, only: mpl_type
use type_nam, only: nam_type

implicit none

integer,parameter :: default_seed = 140587            ! Default seed
integer(kind=int64),parameter :: a = 1103515245_int64 ! Linear congruential multiplier
integer(kind=int64),parameter :: c = 12345_int64      ! Linear congruential offset
integer(kind=int64),parameter :: m = 2147483648_int64 ! Linear congruential modulo

type rng_type
   integer(kind=int64) :: seed
contains
   procedure :: init => rng_init
   procedure :: reseed => rng_reseed
   procedure :: resync => rng_resync
   procedure :: desync => rng_desync
   procedure :: lcg => rng_lcg
   procedure :: rng_rand_integer_0d
   procedure :: rng_rand_integer_1d
   generic :: rand_integer => rng_rand_integer_0d,rng_rand_integer_1d
   procedure :: rng_rand_real_0d
   procedure :: rng_rand_real_1d
   procedure :: rng_rand_real_2d
   procedure :: rng_rand_real_3d
   procedure :: rng_rand_real_4d
   procedure :: rng_rand_real_5d
   generic :: rand_real => rng_rand_real_0d,rng_rand_real_1d,rng_rand_real_2d,rng_rand_real_3d,rng_rand_real_4d,rng_rand_real_5d
   procedure :: rng_rand_gau_1d
   procedure :: rng_rand_gau_5d
   generic :: rand_gau => rng_rand_gau_1d,rng_rand_gau_5d
end type rng_type

private
public :: rng_type

contains

!----------------------------------------------------------------------
! Subroutine: rng_init
!> Initialize the random number generator
!----------------------------------------------------------------------
subroutine rng_init(rng,mpl,nam)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
type(mpl_type),intent(inout) :: mpl  !< MPI data
type(nam_type),intent(in) :: nam     !< Namelist variables

! Local variable
integer :: seed

! Set seed
if (nam%default_seed) then
   ! Default seed
   seed = default_seed
else
   ! Time-based seed
   call system_clock(count=seed)
end if

! Different seed for each processor
seed = seed+mpl%myproc-1

! Long integer
rng%seed = int(seed,kind=int64)

! Print result
if (nam%default_seed) then
   write(mpl%info,'(a7,a)') '','Linear congruential generator initialized with a default seed'
   call mpl%flush
else
   write(mpl%info,'(a7,a)') '','Linear congruential generator initialized'
   call mpl%flush
end if

end subroutine rng_init

!----------------------------------------------------------------------
! Subroutine: rng_reseed
!> Re-seed the random number generator
!----------------------------------------------------------------------
subroutine rng_reseed(rng,mpl)

implicit none

! Passed variable
class(rng_type),intent(inout) :: rng !< Random number generator
type(mpl_type),intent(inout) :: mpl  !< MPI data

! Local variable
integer :: seed

! Default seed
seed = default_seed

! Different seed for each processor
seed = seed+mpl%myproc-1

! Long integer
rng%seed = int(seed,kind=int64)

end subroutine rng_reseed

!----------------------------------------------------------------------
! Subroutine: rng_resync
!> Resynchronize the random number generator between processors
!----------------------------------------------------------------------
subroutine rng_resync(rng,mpl)

implicit none

! Passed variable
class(rng_type),intent(inout) :: rng !< Random number generator
type(mpl_type),intent(inout) :: mpl  !< MPI data

! Wait
call mpl%f_comm%barrier()

! Broadcast root seed
call mpl%f_comm%broadcast(rng%seed,mpl%rootproc-1)

end subroutine rng_resync

!----------------------------------------------------------------------
! Subroutine: rng_desync
!> Desynchronize the random number generator between processors
!----------------------------------------------------------------------
subroutine rng_desync(rng,mpl)

implicit none

! Passed variable
class(rng_type),intent(inout) :: rng !< Random number generator
type(mpl_type),intent(inout) :: mpl  !< MPI data

! Wait
call mpl%f_comm%barrier()

! Different seed for each processor
rng%seed = rng%seed+int(mpl%myproc-1,kind=int64)

end subroutine rng_desync

!----------------------------------------------------------------------
! Subroutine: rng_lcg
!> Linear congruential generator
!----------------------------------------------------------------------
subroutine rng_lcg(rng,x)

implicit none

! Passed variable
class(rng_type),intent(inout) :: rng !< Random number generator
real(kind_real),intent(out) :: x             !< Random number between 0 and 1

! Update seed
rng%seed = mod(a*rng%seed+c,m)

! Random number
x = real(rng%seed,kind_real)/real(m-1,kind_real)

end subroutine rng_lcg

!----------------------------------------------------------------------
! Subroutine: rng_rand_integer_0d
!> Generate a random integer, 0d
!----------------------------------------------------------------------
subroutine rng_rand_integer_0d(rng,binf,bsup,ir)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
integer,intent(in) :: binf           !< Lower bound
integer,intent(in) :: bsup           !< Upper bound
integer,intent(out) :: ir            !< Random integer

! Local variables
real(kind_real) :: x

! Generate random number between 0 and 1
call rng%lcg(x)

! Adapt range
x = x*real(bsup-binf+1,kind_real)

! Add offset
ir = binf+int(x)

end subroutine rng_rand_integer_0d

!----------------------------------------------------------------------
! Subroutine: rng_rand_integer_1d
!> Generate a random integer, 1d
!----------------------------------------------------------------------
subroutine rng_rand_integer_1d(rng,binf,bsup,ir)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
integer,intent(in) :: binf           !< Lower bound
integer,intent(in) :: bsup           !< Upper bound
integer,intent(out) :: ir(:)         !< Random integer

! Local variables
integer :: i

do i=1,size(ir)
   call rng%rand_integer(binf,bsup,ir(i))
end do

end subroutine rng_rand_integer_1d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_0d
!> Generate a random real, 0d
!----------------------------------------------------------------------
subroutine rng_rand_real_0d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
real(kind_real),intent(in) :: binf   !< Lower bound
real(kind_real),intent(in) :: bsup   !< Upper bound
real(kind_real),intent(out) :: rr    !< Random integer

! Local variables
real(kind_real) :: x

! Generate random number between 0 and 1
call rng%lcg(x)

! Adapt range
x = x*(bsup-binf)

! Add offset
rr = binf+x

end subroutine rng_rand_real_0d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_1d
!> Generate a random real, 1d
!----------------------------------------------------------------------
subroutine rng_rand_real_1d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
real(kind_real),intent(in) :: binf   !< Lower bound
real(kind_real),intent(in) :: bsup   !< Upper bound
real(kind_real),intent(out) :: rr(:) !< Random integer

! Local variables
integer :: i

do i=1,size(rr)
   call rng%rand_real(binf,bsup,rr(i))
end do

end subroutine rng_rand_real_1d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_2d
!> Generate a random real, 2d
!----------------------------------------------------------------------
subroutine rng_rand_real_2d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng    !< Random number generator
real(kind_real),intent(in) :: binf      !< Lower bound
real(kind_real),intent(in) :: bsup      !< Upper bound
real(kind_real),intent(out) :: rr(:,:)  !< Random integer

! Local variables
integer :: i,j

do i=1,size(rr,2)
   do j=1,size(rr,1)
   call rng%rand_real(binf,bsup,rr(j,i))
   end do
end do

end subroutine rng_rand_real_2d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_3d
!> Generate a random real, 3d
!----------------------------------------------------------------------
subroutine rng_rand_real_3d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng     !< Random number generator
real(kind_real),intent(in) :: binf       !< Lower bound
real(kind_real),intent(in) :: bsup       !< Upper bound
real(kind_real),intent(out) :: rr(:,:,:) !< Random integer

! Local variables
integer :: i,j,k

do i=1,size(rr,3)
   do j=1,size(rr,2)
      do k=1,size(rr,1)
         call rng%rand_real(binf,bsup,rr(k,j,i))
      end do
   end do
end do

end subroutine rng_rand_real_3d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_4d
!> Generate a random real, 4d
!----------------------------------------------------------------------
subroutine rng_rand_real_4d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng       !< Random number generator
real(kind_real),intent(in) :: binf         !< Lower bound
real(kind_real),intent(in) :: bsup         !< Upper bound
real(kind_real),intent(out) :: rr(:,:,:,:) !< Random integer

! Local variables
integer :: i,j,k,l

do i=1,size(rr,4)
   do j=1,size(rr,3)
      do k=1,size(rr,2)
         do l=1,size(rr,1)
            call rng%rand_real(binf,bsup,rr(l,k,j,i))
         end do
      end do
   end do
end do

end subroutine rng_rand_real_4d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_5d
!> Generate a random real, 5d
!----------------------------------------------------------------------
subroutine rng_rand_real_5d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng         !< Random number generator
real(kind_real),intent(in) :: binf           !< Lower bound
real(kind_real),intent(in) :: bsup           !< Upper bound
real(kind_real),intent(out) :: rr(:,:,:,:,:) !< Random integer

! Local variables
integer :: i,j,k,l,m

do i=1,size(rr,5)
   do j=1,size(rr,4)
      do k=1,size(rr,3)
         do l=1,size(rr,3)
            do m=1,size(rr,1)
               call rng%rand_real(binf,bsup,rr(m,l,k,j,i))
            end do
         end do
      end do
   end do
end do

end subroutine rng_rand_real_5d

!----------------------------------------------------------------------
! Subroutine: rng_rand_gau_1d
!> Generate random Gaussian deviates, 1d
!----------------------------------------------------------------------
subroutine rng_rand_gau_1d(rng,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng !< Random number generator
real(kind_real),intent(out) :: rr(:) !< Random integer

! Local variables
integer :: i,iset
real(kind_real) :: gasdev,fac,gset,rsq,v1,v2

! Normal distribution
iset = 0
do i=1,size(rr,1)
   if (iset==0) then
      rsq = 0.0
      do while((rsq>=1.0).or.(rsq<=0.0))
         call rng%rand_real(0.0_kind_real,1.0_kind_real,v1)
         v1 = 2.0*v1-1.0
         call rng%rand_real(0.0_kind_real,1.0_kind_real,v2)
         v2 = 2.0*v2-1.0
         rsq = v1**2+v2**2
      end do
      fac = sqrt(-2.0*log(rsq)/rsq)
      gset = v1*fac
      gasdev = v2*fac
      iset = 1
   else
      gasdev = gset
      iset = 0
   end if
   rr(i) = gasdev
end do

end subroutine rng_rand_gau_1d

!----------------------------------------------------------------------
! Subroutine: rng_rand_gau_5d
!> Generate random Gaussian deviates, 5d
!----------------------------------------------------------------------
subroutine rng_rand_gau_5d(rng,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng         !< Random number generator
real(kind_real),intent(out) :: rr(:,:,:,:,:) !< Random integer

! Local variables
integer :: i,j,k,l

do i=1,size(rr,5)
   do j=1,size(rr,4)
      do k=1,size(rr,3)
         do l=1,size(rr,2)
            call rng%rand_gau(rr(:,l,k,j,i))
         end do
      end do
   end do
end do

end subroutine rng_rand_gau_5d

end module type_rng

