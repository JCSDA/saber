!----------------------------------------------------------------------
! Module: type_rng
! Purpose: random numbers generator derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_rng

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_status
use iso_fortran_env, only: int64
use tools_const, only: pi
use tools_func, only: sphere_dist
use tools_kinds, only: kind_real
use tools_qsort, only: qsort
use tools_repro, only: repro,inf,sup,infeq
use type_tree, only: tree_type
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
   procedure :: initialize_sampling => rng_initialize_sampling
end type rng_type

private
public :: rng_type

contains

!----------------------------------------------------------------------
! Subroutine: rng_init
! Purpose: initialize the random number generator
!----------------------------------------------------------------------
subroutine rng_init(rng,mpl,nam)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist variables

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

! Different seed for each task
seed = seed+mpl%myproc

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
! Purpose: re-seed the random number generator
!----------------------------------------------------------------------
subroutine rng_reseed(rng,mpl)

implicit none

! Passed variable
class(rng_type),intent(inout) :: rng ! Random number generator
type(mpl_type),intent(inout) :: mpl  ! MPI data

! Local variable
integer :: seed

! Default seed
seed = default_seed

! Different seed for each task
seed = seed+mpl%myproc

! Long integer
rng%seed = int(seed,kind=int64)

end subroutine rng_reseed

!----------------------------------------------------------------------
! Subroutine: rng_lcg
! Purpose: linear congruential generator
!----------------------------------------------------------------------
subroutine rng_lcg(rng,x)

implicit none

! Passed variable
class(rng_type),intent(inout) :: rng ! Random number generator
real(kind_real),intent(out) :: x             ! Random number between 0 and 1

! Update seed
rng%seed = mod(a*rng%seed+c,m)

! Random number
x = real(rng%seed,kind_real)/real(m-1,kind_real)

end subroutine rng_lcg

!----------------------------------------------------------------------
! Subroutine: rng_rand_integer_0d
! Purpose: generate a random integer, 0d
!----------------------------------------------------------------------
subroutine rng_rand_integer_0d(rng,binf,bsup,ir)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
integer,intent(in) :: binf           ! Lower bound
integer,intent(in) :: bsup           ! Upper bound
integer,intent(out) :: ir            ! Random integer

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
! Purpose: generate a random integer, 1d
!----------------------------------------------------------------------
subroutine rng_rand_integer_1d(rng,binf,bsup,ir)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
integer,intent(in) :: binf           ! Lower bound
integer,intent(in) :: bsup           ! Upper bound
integer,intent(out) :: ir(:)         ! Random integer

! Local variables
integer :: i

do i=1,size(ir)
   call rng%rand_integer(binf,bsup,ir(i))
end do

end subroutine rng_rand_integer_1d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_0d
! Purpose: generate a random real, 0d
!----------------------------------------------------------------------
subroutine rng_rand_real_0d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
real(kind_real),intent(in) :: binf   ! Lower bound
real(kind_real),intent(in) :: bsup   ! Upper bound
real(kind_real),intent(out) :: rr    ! Random integer

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
! Purpose: generate a random real, 1d
!----------------------------------------------------------------------
subroutine rng_rand_real_1d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
real(kind_real),intent(in) :: binf   ! Lower bound
real(kind_real),intent(in) :: bsup   ! Upper bound
real(kind_real),intent(out) :: rr(:) ! Random integer

! Local variables
integer :: i

do i=1,size(rr)
   call rng%rand_real(binf,bsup,rr(i))
end do

end subroutine rng_rand_real_1d

!----------------------------------------------------------------------
! Subroutine: rng_rand_real_2d
! Purpose: generate a random real, 2d
!----------------------------------------------------------------------
subroutine rng_rand_real_2d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng    ! Random number generator
real(kind_real),intent(in) :: binf      ! Lower bound
real(kind_real),intent(in) :: bsup      ! Upper bound
real(kind_real),intent(out) :: rr(:,:)  ! Random integer

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
! Purpose: generate a random real, 3d
!----------------------------------------------------------------------
subroutine rng_rand_real_3d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng     ! Random number generator
real(kind_real),intent(in) :: binf       ! Lower bound
real(kind_real),intent(in) :: bsup       ! Upper bound
real(kind_real),intent(out) :: rr(:,:,:) ! Random integer

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
! Purpose: generate a random real, 4d
!----------------------------------------------------------------------
subroutine rng_rand_real_4d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng       ! Random number generator
real(kind_real),intent(in) :: binf         ! Lower bound
real(kind_real),intent(in) :: bsup         ! Upper bound
real(kind_real),intent(out) :: rr(:,:,:,:) ! Random integer

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
! Purpose: generate a random real, 5d
!----------------------------------------------------------------------
subroutine rng_rand_real_5d(rng,binf,bsup,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng         ! Random number generator
real(kind_real),intent(in) :: binf           ! Lower bound
real(kind_real),intent(in) :: bsup           ! Upper bound
real(kind_real),intent(out) :: rr(:,:,:,:,:) ! Random integer

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
! Purpose: generate random Gaussian deviates, 1d
!----------------------------------------------------------------------
subroutine rng_rand_gau_1d(rng,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng ! Random number generator
real(kind_real),intent(out) :: rr(:) ! Random integer

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
! Purpose: generate random Gaussian deviates, 5d
!----------------------------------------------------------------------
subroutine rng_rand_gau_5d(rng,rr)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng         ! Random number generator
real(kind_real),intent(out) :: rr(:,:,:,:,:) ! Random integer

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

!----------------------------------------------------------------------
! Subroutine: rng_initialize_sampling
! Purpose: intialize sampling
!----------------------------------------------------------------------
subroutine rng_initialize_sampling(rng,mpl,n_loc,lon_loc,lat_loc,mask_loc,rh_loc,loc_to_glb,ntry,nrep,ns2_glb,sam2_glb,fast)

implicit none

! Passed variables
class(rng_type),intent(inout) :: rng          ! Random number generator
type(mpl_type),intent(inout) :: mpl           ! MPI data
integer,intent(in) :: n_loc                   ! Number of points (local)
real(kind_real),intent(in) :: lon_loc(n_loc)  ! Longitudes (local)
real(kind_real),intent(in) :: lat_loc(n_loc)  ! Latitudes (local)
logical,intent(in) :: mask_loc(n_loc)         ! Mask (local)
real(kind_real),intent(in) :: rh_loc(n_loc)   ! Horizontal support radius (local)
integer,intent(in) :: loc_to_glb(n_loc)       ! Local to global index
integer,intent(in) :: ntry                    ! Number of tries
integer,intent(in) :: nrep                    ! Number of replacements
integer,intent(in) :: ns2_glb                 ! Number of samplings points (global)
integer,intent(out) :: sam2_glb(ns2_glb)      ! Horizontal sampling index (global)
logical,intent(in),optional :: fast           ! Fast sampling flag

! Local variables
integer :: n_glb,n_loc_eff,n_glb_eff,n_loc_val,i_loc,i_loc_val,i_glb,is1_loc,is2_glb,iproc,js,irep,irmax,itry,is1_glb,ir
integer :: irval,irvalmin,irvalmax,is2_glb_min,nrep_eff,nn_index(2),proc_to_ns1_loc(mpl%nproc),offset,ns1_loc,ns1_glb
integer :: ns1_glb_val
integer,allocatable :: glb_to_loc(:),glb_to_proc(:),sam1_loc(:),sam1_glb(:),to_valid(:),sam2_glb_tmp(:),order(:)
real(kind_real) :: rhsq_loc,rhsq_glb,d,distmax,distmin,nn_dist(2),cdf_norm,rr
real(kind_real),allocatable :: lon1_glb(:),lat1_glb(:),rh1_glb(:)
real(kind_real),allocatable :: cdf(:)
real(kind_real),allocatable :: lon_rep(:),lat_rep(:),dist(:)
logical :: lfast
logical,allocatable :: mask_glb(:),lmask(:),smask(:),rmask(:)
character(len=1024),parameter :: subr = 'rng_initialize_sampling'
type(fckit_mpi_status) :: status
type(tree_type) :: tree

! Number of effective points
n_loc_eff = count(mask_loc)
call mpl%f_comm%allreduce(n_loc_eff,n_glb_eff,fckit_mpi_sum())

! Check mask size
if (n_glb_eff==0) then
   call mpl%abort(subr,'empty mask in initialize sampling')
elseif (n_glb_eff<ns2_glb) then
   call mpl%abort(subr,'ns2_glb greater than n_glb_eff in initialize_sampling')
elseif (n_glb_eff==ns2_glb) then
   write(mpl%info,'(a)') ' all points are used'
   call mpl%flush

   ! Global size
   call mpl%f_comm%allreduce(n_loc,n_glb,fckit_mpi_sum())

   ! Allocation
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
   allocate(mask_glb(n_glb))

   ! Communication
   call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)
   call mpl%loc_to_glb(n_loc,mask_loc,n_glb,glb_to_proc,glb_to_loc,.false.,mask_glb)  

   if (mpl%main) then
      ! Use all valid points
      is2_glb = 0
      do i_glb=1,n_glb
         if (mask_glb(i_glb)) then
            is2_glb = is2_glb+1
            sam2_glb(is2_glb) = i_glb
         end if
      end do
   end if

   ! Release memory
   deallocate(glb_to_loc)
   deallocate(glb_to_proc)
   deallocate(mask_glb)
else
   ! Define number of subsampling points for the first step
   if (repro) then
      ! First step is skipped
      ns1_loc = n_loc_eff
   else
      ! First subsampling depends on horizontal radius
      rhsq_loc = sum(rh_loc**2,mask=mask_loc)
      call mpl%f_comm%allreduce(rhsq_loc,rhsq_glb,fckit_mpi_sum())
      ns1_loc = int(5.0*real(ns2_glb,kind_real)*rhsq_loc/rhsq_glb)
      ns1_loc = min(ns1_loc,n_loc_eff)
   end if

   ! Communication
   if (mpl%main) then
      ! Receive data on rootproc
      do iproc=1,mpl%nproc
         if (iproc==mpl%rootproc) then
            ! Copy data
            proc_to_ns1_loc(iproc) = ns1_loc
         else
            ! Receive data
            call mpl%f_comm%receive(proc_to_ns1_loc(iproc),iproc-1,mpl%tag,status)
         end if
      end do
   else
      ! Send data to rootproc
      call mpl%f_comm%send(ns1_loc,mpl%rootproc-1,mpl%tag)
   end if
   call mpl%update_tag(1)

   ! Global size
   if (mpl%main) ns1_glb = sum(proc_to_ns1_loc)

   ! Allocation
   allocate(sam1_loc(ns1_loc))
   
   ! First subsampling
   if (n_loc_eff==0) then
      ! No point on this task
      write(mpl%info,'(a)') ' no point on this task'
      call mpl%flush
   elseif (ns1_loc==n_loc_eff) then
      ! All points are used
      write(mpl%info,'(a)') ' all points are used'
      call mpl%flush(.false.)

      is1_loc = 0
      do i_loc=1,n_loc
         if (mask_loc(i_loc)) then
            is1_loc = is1_loc+1
            sam1_loc(is1_loc) = i_loc
         end if
      end do
   else
      ! Define sampling with a cumulative distribution function

      ! Allocation
      n_loc_val = n_loc_eff
      allocate(to_valid(n_loc_val))
      allocate(cdf(n_loc_val))

      ! Initialization
      sam1_loc = mpl%msv%vali
      to_valid = mpl%msv%vali
      i_loc_val = 0
      do i_loc=1,n_loc
         if (mask_loc(i_loc)) then
            i_loc_val = i_loc_val+1
            to_valid(i_loc_val) = i_loc
         end if
      end do
      cdf(1) = 0.0
      do i_loc_val=2,n_loc_val
         i_loc = to_valid(i_loc_val)
         cdf(i_loc_val) = cdf(i_loc_val-1)+1.0/rh_loc(i_loc)**2
      end do
      cdf_norm = 1.0/cdf(n_loc_val)
      cdf = cdf*cdf_norm
      call mpl%prog_init(ns1_loc)

      do is1_loc=1,ns1_loc
         ! Generate random number
         call rng%rand_real(0.0_kind_real,1.0_kind_real,rr)

         ! Dichotomy to find the value
         irval = 1
         irvalmin = 1
         irvalmax = n_loc_val
         do while (irvalmax-irvalmin>1)
            irval = (irvalmin+irvalmax)/2
            if ((cdf(irvalmin)-rr)*(cdf(irval)-rr)>0.0) then
               irvalmin = irval
            else
               irvalmax = irval
            end if
         end do

         ! New sampling point
         ir = to_valid(irval)
         sam1_loc(is1_loc) = ir

         ! Shift valid points array
         if (irval<n_loc_val) then
            cdf(irval:n_loc_val-1) = cdf(irval+1:n_loc_val)
            to_valid(irval:n_loc_val-1) = to_valid(irval+1:n_loc_val)
         end if
         n_loc_val = n_loc_val-1

         ! Renormalize cdf
         cdf_norm = 1.0/cdf(n_loc_val)
         cdf(1:n_loc_val) = cdf(1:n_loc_val)*cdf_norm

         ! Update
         call mpl%prog_print(is1_loc)
      end do
      call mpl%prog_final(.false.)

      ! Release memory
      deallocate(to_valid)
      deallocate(cdf)
   end if

   if (n_loc_eff>0) then
      ! Continue printing
      write(mpl%info,'(a)') ' => '
      call mpl%flush(.false.)
   end if

   ! Communication
   if (mpl%main) then
      ! Allocation
      allocate(sam1_glb(ns1_glb))
      allocate(lon1_glb(ns1_glb))
      allocate(lat1_glb(ns1_glb))
      allocate(rh1_glb(ns1_glb))

      ! Receive data on rootproc
      offset = 0
      do iproc=1,mpl%nproc
         if (proc_to_ns1_loc(iproc)>0) then
            if (iproc==mpl%rootproc) then
               ! Copy data
               sam1_glb(offset+1:offset+proc_to_ns1_loc(iproc)) = loc_to_glb(sam1_loc)
               lon1_glb(offset+1:offset+proc_to_ns1_loc(iproc)) = lon_loc(sam1_loc)
               lat1_glb(offset+1:offset+proc_to_ns1_loc(iproc)) = lat_loc(sam1_loc)
               rh1_glb(offset+1:offset+proc_to_ns1_loc(iproc)) = rh_loc(sam1_loc)   
            else
               ! Receive data
               call mpl%f_comm%receive(sam1_glb(offset+1:offset+proc_to_ns1_loc(iproc)),iproc-1,mpl%tag,status)
               call mpl%f_comm%receive(lon1_glb(offset+1:offset+proc_to_ns1_loc(iproc)),iproc-1,mpl%tag+1,status)
               call mpl%f_comm%receive(lat1_glb(offset+1:offset+proc_to_ns1_loc(iproc)),iproc-1,mpl%tag+2,status)
               call mpl%f_comm%receive(rh1_glb(offset+1:offset+proc_to_ns1_loc(iproc)),iproc-1,mpl%tag+3,status)
            end if

            ! Update
            offset = offset+proc_to_ns1_loc(iproc)
         end if
      end do

      if (repro) then
         ! Reorder data
         allocate(order(ns1_glb))
         call qsort(ns1_glb,sam1_glb,order)
         lon1_glb = lon1_glb(order)
         lat1_glb = lat1_glb(order)
         rh1_glb = rh1_glb(order)
         deallocate(order)
      end if
   else
      if (ns1_loc>0) then
         ! Send data to rootproc
         call mpl%f_comm%send(loc_to_glb(sam1_loc),mpl%rootproc-1,mpl%tag)
         call mpl%f_comm%send(lon_loc(sam1_loc),mpl%rootproc-1,mpl%tag+1)
         call mpl%f_comm%send(lat_loc(sam1_loc),mpl%rootproc-1,mpl%tag+2)
         call mpl%f_comm%send(rh_loc(sam1_loc),mpl%rootproc-1,mpl%tag+3)
      end if
   end if
   call mpl%update_tag(4)

   ! Release memory
   deallocate(sam1_loc)

   if (mpl%main) then
      ! Allocation
      nrep_eff = min(nrep,ns1_glb-ns2_glb)
      allocate(sam2_glb_tmp(ns2_glb+nrep_eff))
      allocate(lmask(ns1_glb))
      allocate(smask(ns1_glb))
      allocate(to_valid(ns1_glb))

      ! Initialization
      sam2_glb_tmp = mpl%msv%vali
      lmask = .true.
      smask = .false.
      to_valid = mpl%msv%vali
      do is1_glb=1,ns1_glb
         to_valid(is1_glb) = is1_glb
      end do
      call mpl%prog_init(ns2_glb+nrep_eff)
      lfast = .false.
      if (present(fast)) lfast = fast

      if (lfast) then
         ! Define sampling with a cumulative distribution function

         ! Allocation
         allocate(cdf(ns1_glb))

         ! Initialization
         ns1_glb_val = ns1_glb
         cdf(1) = 0.0
         do is1_glb=2,ns1_glb
            cdf(is1_glb) = cdf(is1_glb-1)+1.0/rh1_glb(is1_glb)**2
         end do
         cdf_norm = 1.0/cdf(ns1_glb)
         cdf = cdf*cdf_norm

         do is2_glb=1,ns2_glb+nrep_eff
            ! Generate random number
            call rng%rand_real(0.0_kind_real,1.0_kind_real,rr)

            ! Dichotomy to find the value
            irvalmin = 1
            irvalmax = ns1_glb_val
            do while (irvalmax-irvalmin>1)
               irval = (irvalmin+irvalmax)/2
               if ((sup(cdf(irvalmin),rr).and.sup(cdf(irval),rr)).or.(inf(cdf(irvalmin),rr).and.inf(cdf(irval),rr))) then
                  irvalmin = irval
               else
                  irvalmax = irval
               end if
            end do

            ! New sampling point
            ir = to_valid(irval)
            sam2_glb_tmp(is2_glb) = ir

            ! Shift valid points array
            if (irval<ns1_glb_val) then
               cdf(irval:ns1_glb_val-1) = cdf(irval+1:ns1_glb_val)
               to_valid(irval:ns1_glb_val-1) = to_valid(irval+1:ns1_glb_val)
            end if
            ns1_glb_val = ns1_glb_val-1

            ! Renormalize cdf
            cdf_norm = 1.0/cdf(ns1_glb_val)
            cdf(1:ns1_glb_val) = cdf(1:ns1_glb_val)*cdf_norm

            ! Update
            call mpl%prog_print(is2_glb)
         end do
         call mpl%prog_final(.false.)

         ! Release memory
         deallocate(cdf)
      else
         ! Define sampling with a KD-tree

         ! Initialization
         ns1_glb_val = ns1_glb

         do is2_glb=1,ns2_glb+nrep_eff
            if (is2_glb>2) then
               ! Allocation
               call tree%alloc(mpl,ns1_glb,mask=smask)

               ! Initialization
               call tree%init(lon1_glb,lat1_glb)
            end if

            ! Initialization
            distmax = 0.0
            irmax = 0
            irvalmax = 0
            itry = 1

            ! Find a new point
            do itry=1,ntry
               ! Generate a random index among valid points
               call rng%rand_integer(1,ns1_glb_val,irval)
               ir = to_valid(irval)

               ! Check point validity
               if (is2_glb==1) then
                  ! Accept point
                  irvalmax = irval
                  irmax = ir
               else
                  if (is2_glb==2) then
                     ! Compute distance
                     nn_index(1) = sam2_glb_tmp(1)
                     call sphere_dist(lon1_glb(ir),lat1_glb(ir),lon1_glb(nn_index(1)),lat1_glb(nn_index(1)),nn_dist(1))
                  else
                     ! Find nearest neighbor distance
                     call tree%find_nearest_neighbors(lon1_glb(ir),lat1_glb(ir),1,nn_index(1:1),nn_dist(1:1))
                  end if
                  d = nn_dist(1)**2/(rh1_glb(ir)**2+rh1_glb(nn_index(1))**2)

                  ! Check distance
                  if (sup(d,distmax)) then
                     distmax = d
                     irvalmax = irval
                     irmax = ir
                  end if
               end if
            end do

            ! Delete tree
            if (is2_glb>2) call tree%dealloc

            ! Add point to sampling
            if (irmax>0) then
               ! New sampling point
               sam2_glb_tmp(is2_glb) = irmax
               lmask(irmax) = .false.
               smask(irmax) = .true.

               ! Shift valid points array
               if (irvalmax<ns1_glb_val) to_valid(irvalmax:ns1_glb_val-1) = to_valid(irvalmax+1:ns1_glb_val)
               ns1_glb_val = ns1_glb_val-1
            end if

            ! Update
            call mpl%prog_print(is2_glb)
         end do
         call mpl%prog_final(.false.)
      end if

      if (nrep_eff>0) then
         if (n_loc_eff>0) then
            ! Continue printing
            write(mpl%info,'(a)') ' => '
            call mpl%flush(.false.)
         end if

         ! Allocation
         allocate(rmask(ns2_glb+nrep_eff))
         allocate(lon_rep(ns2_glb+nrep_eff))
         allocate(lat_rep(ns2_glb+nrep_eff))
         allocate(dist(ns2_glb+nrep_eff))

         ! Initialization
         rmask = .true.
         do is2_glb=1,ns2_glb+nrep_eff
            lon_rep(is2_glb) = lon1_glb(sam2_glb_tmp(is2_glb))
            lat_rep(is2_glb) = lat1_glb(sam2_glb_tmp(is2_glb))
         end do
         dist = mpl%msv%valr
         call mpl%prog_init(nrep_eff)

         ! Remove closest points
         do irep=1,nrep_eff
            ! Allocation
            call tree%alloc(mpl,ns2_glb+nrep_eff,mask=rmask)

            ! Initialization
            call tree%init(lon_rep,lat_rep)

            ! Get minimum distance
            do is2_glb=1,ns2_glb+nrep_eff
               if (rmask(is2_glb)) then
                  ! Find nearest neighbor distance
                  call tree%find_nearest_neighbors(lon1_glb(sam2_glb_tmp(is2_glb)),lat1_glb(sam2_glb_tmp(is2_glb)), &
                & 2,nn_index,nn_dist)
                  if (nn_index(1)==is2_glb) then
                     dist(is2_glb) = nn_dist(2)
                  elseif (nn_index(2)==is2_glb) then
                     dist(is2_glb) = nn_dist(1)
                  else
                     call mpl%abort(subr,'wrong index in replacement')
                  end if
                  dist(is2_glb) = dist(is2_glb)**2/(rh1_glb(sam2_glb_tmp(nn_index(1)))**2+rh1_glb(sam2_glb_tmp(nn_index(2)))**2)
               end if
            end do

            ! Delete tree
            call tree%dealloc

            ! Remove worst point
            distmin = huge(1.0)
            is2_glb_min = mpl%msv%vali
            do is2_glb=1,ns2_glb+nrep_eff
               if (rmask(is2_glb)) then
                  if (inf(dist(is2_glb),distmin)) then
                     is2_glb_min = is2_glb
                     distmin = dist(is2_glb)
                  end if
               end if
            end do
            rmask(is2_glb_min) = .false.

             ! Update
            call mpl%prog_print(irep)
         end do
         call mpl%prog_final

         ! Copy sam2_glb
         js = 0
         do is2_glb=1,ns2_glb+nrep_eff
            if (rmask(is2_glb)) then
               js = js+1
               sam2_glb(js) = sam2_glb_tmp(is2_glb)
            end if
         end do

         ! Release memory
         deallocate(rmask)
         deallocate(lon_rep)
         deallocate(lat_rep)
         deallocate(dist)
      else
         if (n_loc_eff>0) then
            ! Stop printing
            write(mpl%info,'(a)') ''
            call mpl%flush
         end if

         ! Copy sam2_glb
         sam2_glb = sam2_glb_tmp
      end if

      ! Apply first sampling step
      sam2_glb = sam1_glb(sam2_glb)

      ! Release memory
      deallocate(sam1_glb)
      deallocate(lon1_glb)
      deallocate(lat1_glb)
      deallocate(rh1_glb)
      deallocate(sam2_glb_tmp)
      deallocate(lmask)
      deallocate(smask)
      deallocate(to_valid)
   else
      if (n_loc_eff>0) then
         ! Stop printing
         write(mpl%info,'(a)') ''
         call mpl%flush
      end if
   end if
end if

! Broadcast
call mpl%f_comm%broadcast(sam2_glb,mpl%rootproc-1)

end subroutine rng_initialize_sampling

end module type_rng
