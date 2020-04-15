!----------------------------------------------------------------------
! Module: type_avg
! Purpose: average routines
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_avg

use fckit_mpi_module, only: fckit_mpi_sum
!$ use omp_lib
use tools_func, only: add,divide
use tools_kinds, only: kind_real
use tools_qsort, only: qsort
use type_avg_blk, only: avg_blk_type
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_mom, only: mom_type
use type_mom_blk, only: mom_blk_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type

implicit none

! Averaged statistics derived type
type avg_type
   character(len=1024) :: prefix              ! Prefix
   integer :: ne                              ! Ensemble size
   integer :: nsub                            ! Number of sub-ensembles
   type(avg_blk_type),allocatable :: blk(:,:) ! Averaged statistics blocks
contains
   procedure :: alloc => avg_alloc
   procedure :: dealloc => avg_dealloc
   procedure :: copy => avg_copy
   procedure :: write => avg_write
   procedure :: compute => avg_compute
   procedure :: compute_hyb => avg_compute_hyb
   procedure :: compute_deh => avg_compute_deh
   procedure :: copy_wgt => avg_copy_wgt
   procedure :: compute_bwavg => avg_compute_bwavg
   procedure :: compute_bwavg_hyb => avg_compute_bwavg_hyb
   procedure :: compute_bwavg_deh => avg_compute_bwavg_deh
end type avg_type

private
public :: avg_type

contains

!----------------------------------------------------------------------
! Subroutine: avg_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine avg_alloc(avg,nam,geom,bpar,samp,ne,nsub,prefix)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg  ! Averaged statistics
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(bpar_type),intent(in) :: bpar    ! Block parameters
type(samp_type),intent(in) :: samp    ! Sampling
integer,intent(in) :: ne              ! Ensemble size
integer,intent(in) :: nsub            ! Number of sub-ensembles
character(len=*),intent(in) :: prefix ! Prefix

! Local variables
integer :: ib,ic2a

! Set attributes
avg%prefix = trim(prefix)
avg%ne = ne
avg%nsub = nsub

! Allocation
allocate(avg%blk(0:samp%nc2a,bpar%nbe))
do ib=1,bpar%nbe
   do ic2a=0,samp%nc2a
      call avg%blk(ic2a,ib)%alloc(nam,geom,bpar,ic2a,ib,ne,nsub,prefix)
   end do
end do

end subroutine avg_alloc

!----------------------------------------------------------------------
! Subroutine: avg_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine avg_dealloc(avg)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics

! Local variables
integer :: ib,ic2a

! Allocation
if (allocated(avg%blk)) then
   do ib=1,size(avg%blk,2)
      do ic2a=0,size(avg%blk,1)-1
         call avg%blk(ic2a,ib)%dealloc
      end do
    end do
   deallocate(avg%blk)
end if

end subroutine avg_dealloc

!----------------------------------------------------------------------
! Subroutine: avg_copy
! Purpose: copy
!----------------------------------------------------------------------
subroutine avg_copy(avg_out,avg_in)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_out ! Output averaged statistics
class(avg_type),intent(in) :: avg_in     ! Input averaged statistics

! Local variables
integer :: ib,ic2a

! Copy
do ib=1,size(avg_in%blk,2)
   do ic2a=0,size(avg_in%blk,1)-1
      call avg_out%blk(ic2a,ib)%copy(avg_in%blk(ic2a,ib))
   end do
end do

end subroutine avg_copy

!----------------------------------------------------------------------
! Subroutine: avg_write
! Purpose: write
!----------------------------------------------------------------------
subroutine avg_write(avg,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Diagnostic
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters

! Local variables
integer :: ib
character(len=1024) :: filename

if (mpl%main) then
   filename = trim(nam%prefix)//'_avg'
   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) call avg%blk(0,ib)%write(mpl,nam,geom,bpar,filename)
   end do
end if

end subroutine avg_write

!----------------------------------------------------------------------
! Subroutine: avg_compute
! Purpose: compute averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute(avg,mpl,nam,geom,bpar,samp,mom,ne,prefix)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg  ! Averaged statistics
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(bpar_type),intent(in) :: bpar    ! Block parameters
type(samp_type),intent(in) :: samp    ! Sampling
type(mom_type),intent(in) :: mom      ! Moments
integer,intent(in) :: ne              ! Ensemble size
character(len=*),intent(in) :: prefix ! Prefix

! Local variables
integer :: ib,ic2a
type(mom_blk_type) :: mom_blk

! Allocation
call avg%alloc(nam,geom,bpar,samp,mom%ne,mom%nsub,prefix)

! Compute averaged statistics
write(mpl%info,'(a10,a)') '','Compute averaged statistics'
call mpl%flush
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a13,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)

      ! Global average
      call avg%blk(0,ib)%compute_global(mpl,nam,geom,bpar,samp,mom%blk(ib))

      if (nam%local_diag) then
         ! Moments block extension
         call mom_blk%ext(mpl,geom,bpar,samp,mom%blk(ib))

         ! Local average
         call mpl%prog_init(samp%nc2a)
         do ic2a=1,samp%nc2a
            call avg%blk(ic2a,ib)%compute_local(mpl,nam,geom,bpar,samp,mom_blk)
            call mpl%prog_print(ic2a)
         end do
         call mpl%prog_final

         ! Release memory
         call mom_blk%dealloc
      else
         write(mpl%info,'(a)') ' done'
         call mpl%flush
      end if
   end if
end do

if (mpl%main.and.(nam%avg_nbins>0)) then
   ! Write histograms
   write(mpl%info,'(a10,a)') '','Write histograms'
   call mpl%flush
   call avg%write(mpl,nam,geom,bpar)
end if

! Compute asymptotic statistics
write(mpl%info,'(a10,a)') '','Compute asymptotic statistics:'
call mpl%flush
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a13,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)
      call mpl%prog_init(samp%nc2a+1)
      do ic2a=0,samp%nc2a
         if ((ic2a==0).or.nam%local_diag) call avg%blk(ic2a,ib)%compute_asy(mpl,nam,geom,bpar,ne)
         call mpl%prog_print(ic2a+1)
      end do
      call mpl%prog_final
   end if
end do

end subroutine avg_compute

!----------------------------------------------------------------------
! Subroutine: avg_compute_hyb
! Purpose: compute hybrid averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_hyb(avg_1,mpl,nam,geom,bpar,samp,avg_2)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_1 ! Ensemble 1 averaged statistics
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(samp_type),intent(in) :: samp     ! Sampling
type(avg_type),intent(in) :: avg_2     ! Ensemble 2 averaged statistics

! Local variables
integer :: ib,ic2a

! Compute averaged statistics
write(mpl%info,'(a10,a)') '','Compute averaged statistics'
call mpl%flush
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a13,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)
      call mpl%prog_init(samp%nc2a+1)
      do ic2a=0,samp%nc2a
         if ((ic2a==0).or.nam%local_diag) then
            select case (trim(nam%method))
            case ('hyb-avg')
               ! Static covariance = ensemble covariance
               call avg_1%blk(ic2a,ib)%compute_hyb(mpl,geom,bpar,avg_1%blk(ic2a,ib))
            case ('hyb-rnd')
               ! Static covariance = randomized covariance
               call avg_1%blk(ic2a,ib)%compute_hyb(mpl,geom,bpar,avg_2%blk(ic2a,ib))
            end select
         end if
         call mpl%prog_print(ic2a+1)
      end do
      call mpl%prog_final
   end if
end do

end subroutine avg_compute_hyb

!----------------------------------------------------------------------
! Subroutine: avg_compute_deh
! Purpose: compute dual-ensemble hybrid averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_deh(avg_1,mpl,nam,geom,bpar,samp,mom_1,mom_2)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_1 ! Ensemble 1 averaged statistics
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(samp_type),intent(in) :: samp     ! Sampling
type(mom_type),intent(in) :: mom_1     ! Ensemble 1 moments
type(mom_type),intent(in) :: mom_2     ! Ensemble 2 moments

! Local variables
integer :: ib,ic2a
type(mom_blk_type) :: mom_blk_1,mom_blk_2

! Compute averaged statistics
write(mpl%info,'(a10,a)') '','Compute averaged statistics'
call mpl%flush
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a13,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)

      ! Global average
      call avg_1%blk(0,ib)%compute_deh_global(mpl,nam,geom,bpar,samp,mom_1%blk(ib),mom_2%blk(ib))

      if (nam%local_diag) then
         ! Moments block extension
         call mom_blk_1%ext(mpl,geom,bpar,samp,mom_1%blk(ib))
         call mom_blk_2%ext(mpl,geom,bpar,samp,mom_2%blk(ib))

         ! Local average
         call mpl%prog_init(samp%nc2a)
         do ic2a=1,samp%nc2a
            call avg_1%blk(ic2a,ib)%compute_deh_local(mpl,nam,geom,bpar,samp,mom_blk_1,mom_blk_2)
            call mpl%prog_print(ic2a)
         end do
         call mpl%prog_final

         ! Release memory
         call mom_blk_1%dealloc
         call mom_blk_2%dealloc
      end if
   end if
end do

! Compute asymptotic statistics
write(mpl%info,'(a10,a)') '','Compute asymptotic statistics at low resolution:'
call mpl%flush(.false.)
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a13,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)
      call mpl%prog_init(samp%nc2a+1)
      do ic2a=0,samp%nc2a
         if ((ic2a==0).or.nam%local_diag) call avg_1%blk(ic2a,ib)%compute_asy_deh(mpl,nam,geom,bpar)
         call mpl%prog_print(ic2a+1)
      end do
      call mpl%prog_final
   end if
end do

end subroutine avg_compute_deh

!----------------------------------------------------------------------
! Subroutine: avg_copy_wgt
! Purpose: averaged statistics data copy for weight definition
!----------------------------------------------------------------------
subroutine avg_copy_wgt(avg_out,geom,bpar,avg_in)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg_out ! Output averaged statistics
type(geom_type),intent(in) :: geom       ! Geometry
type(bpar_type),intent(in) :: bpar       ! Block parameters
type(avg_type),intent(in) :: avg_in      ! Input averaged statistics

! Local variables
integer :: ib

if (bpar%diag_block(bpar%nbe)) then
   ! Allocation
   allocate(avg_out%blk(0:0,bpar%nb))

   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) then
         ! Restricted allocation
         allocate(avg_out%blk(0,ib)%m2m2asy(bpar%nc3(ib),bpar%nl0r(ib),geom%nl0))

         ! Restricted copy
         avg_out%blk(0,ib)%m2m2asy = avg_in%blk(0,ib)%m2m2asy
      end if
   end do
end if

end subroutine avg_copy_wgt

!----------------------------------------------------------------------
! Subroutine: avg_compute_bwavg
! Purpose: compute block-averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_bwavg(avg,mpl,nam,geom,bpar,samp,avg_wgt)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
type(samp_type),intent(in) :: samp   ! Sampling
type(avg_type),intent(in) :: avg_wgt ! Averaged statistics for weights

! Local variables
integer :: ib,ic2a,il0,jl0r,jc3
real(kind_real) :: bwgtsq
real(kind_real) :: cor(nam%nc3,bpar%nl0rmax,geom%nl0),nc1a_cor(nam%nc3,bpar%nl0rmax,geom%nl0)
real(kind_real) :: m11asysq(nam%nc3,bpar%nl0rmax,geom%nl0),m11sq(nam%nc3,bpar%nl0rmax,geom%nl0),nc1a(nam%nc3,bpar%nl0rmax,geom%nl0)

! Initialization
write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(bpar%nbe)),':'
call mpl%flush(.false.)
call mpl%prog_init(samp%nc2a+1)

do ic2a=0,samp%nc2a
   ! Copy ensemble size
   avg%blk(ic2a,bpar%nbe)%ne = avg%blk(ic2a,1)%ne
   avg%blk(ic2a,bpar%nbe)%nsub = avg%blk(ic2a,1)%nsub

   if ((ic2a==0).or.nam%local_diag) then
      ! Initialization
      avg%blk(ic2a,bpar%nbe)%cor = 0.0
      cor = 0.0
      avg%blk(ic2a,bpar%nbe)%nc1a_cor = 0.0
      nc1a_cor = 0.0
      avg%blk(ic2a,bpar%nbe)%m11asysq = 0.0
      m11asysq = 0.0
      avg%blk(ic2a,bpar%nbe)%m11sq = 0.0
      m11sq = 0.0
      avg%blk(ic2a,bpar%nbe)%nc1a = 0.0
      nc1a = 0.0

      ! Block averages
      do ib=1,bpar%nb
         if (bpar%avg_block(ib)) then
            !$omp parallel do schedule(static) private(il0,jl0r,bwgtsq,jc3)
            do il0=1,geom%nl0
               do jl0r=1,bpar%nl0r(ib)
                  ! Weight
                  if (avg_wgt%blk(0,ib)%m2m2asy(1,jl0r,il0)>0.0) then
                     bwgtsq = 1.0/avg_wgt%blk(0,ib)%m2m2asy(1,jl0r,il0)
                  else
                     bwgtsq = 0.0
                  end if

                  ! Compute sum
                  do jc3=1,nam%nc3
                     call add(mpl,avg%blk(ic2a,ib)%cor(jc3,jl0r,il0),avg%blk(ic2a,bpar%nbe)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
                     call add(mpl,avg%blk(ic2a,ib)%nc1a_cor(jc3,jl0r,il0),avg%blk(ic2a,bpar%nbe)%nc1a_cor(jc3,jl0r,il0), &
                   & nc1a_cor(jc3,jl0r,il0))
                     call add(mpl,avg%blk(ic2a,ib)%m11asysq(jc3,jl0r,il0),avg%blk(ic2a,bpar%nbe)%m11asysq(jc3,jl0r,il0), &
                   & m11asysq(jc3,jl0r,il0),bwgtsq)
                     call add(mpl,avg%blk(ic2a,ib)%m11sq(jc3,jl0r,il0),avg%blk(ic2a,bpar%nbe)%m11sq(jc3,jl0r,il0), &
                   & m11sq(jc3,jl0r,il0),bwgtsq)
                     call add(mpl,avg%blk(ic2a,ib)%nc1a(jc3,jl0r,il0),avg%blk(ic2a,bpar%nbe)%nc1a(jc3,jl0r,il0), &
                   & nc1a(jc3,jl0r,il0),bwgtsq)
                  end do
               end do
            end do
            !$omp end parallel do
         end if
      end do

      ! Normalization
      !$omp parallel do schedule(static) private(il0,jl0r,jc3)
      do il0=1,geom%nl0
         do jl0r=1,bpar%nl0r(ib)
            do jc3=1,nam%nc3
               call divide(mpl,avg%blk(ic2a,bpar%nbe)%cor(jc3,jl0r,il0),cor(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2a,bpar%nbe)%nc1a_cor(jc3,jl0r,il0),nc1a_cor(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2a,bpar%nbe)%m11asysq(jc3,jl0r,il0),m11asysq(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2a,bpar%nbe)%m11sq(jc3,jl0r,il0),m11sq(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2a,bpar%nbe)%nc1a(jc3,jl0r,il0),nc1a(jc3,jl0r,il0))
            end do
         end do
      end do
      !$omp end parallel do
   end if

   ! Update
   call mpl%prog_print(ic2a+1)
end do
call mpl%prog_final

end subroutine avg_compute_bwavg

!----------------------------------------------------------------------
! Subroutine: avg_compute_bwavg_hyb
! Purpose: compute hybrid block-averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_bwavg_hyb(avg,mpl,nam,geom,bpar,samp,avg_wgt)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
type(samp_type),intent(in) :: samp   ! Sampling
type(avg_type),intent(in) :: avg_wgt ! Averaged statistics for weights

! Local variables
integer :: ib,ic2a,il0,jl0r,jc3
real(kind_real) :: bwgtsq
real(kind_real) :: m11sta(nam%nc3,bpar%nl0rmax,geom%nl0),stasq(nam%nc3,bpar%nl0rmax,geom%nl0)

! Initialization
write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(bpar%nbe)),':'
call mpl%flush(.false.)
call mpl%prog_init(samp%nc2a+1)

do ic2a=0,samp%nc2a
   if ((ic2a==0).or.nam%local_diag) then
      ! Initialization
      avg%blk(ic2a,bpar%nbe)%m11sta = 0.0
      m11sta = 0.0
      avg%blk(ic2a,bpar%nbe)%stasq = 0.0
      stasq = 0.0

      ! Block averages
      do ib=1,bpar%nb
         if (bpar%avg_block(ib)) then
            !$omp parallel do schedule(static) private(il0,jl0r,bwgtsq,jc3)
            do il0=1,geom%nl0
               do jl0r=1,bpar%nl0r(ib)
                  ! Weight
                  if (avg_wgt%blk(0,ib)%m2m2asy(1,jl0r,il0)>0.0) then
                     bwgtsq = 1.0/avg_wgt%blk(0,ib)%m2m2asy(1,jl0r,il0)
                  else
                     bwgtsq = 0.0
                  end if

                  ! Compute sum
                  do jc3=1,nam%nc3
                     call add(mpl,avg%blk(ic2a,ib)%m11sta(jc3,jl0r,il0),avg%blk(ic2a,bpar%nbe)%m11sta(jc3,jl0r,il0), &
                   & m11sta(jc3,jl0r,il0),bwgtsq)
                     call add(mpl,avg%blk(ic2a,ib)%stasq(jc3,jl0r,il0),avg%blk(ic2a,bpar%nbe)%stasq(jc3,jl0r,il0), &
                   & stasq(jc3,jl0r,il0),bwgtsq)
                  end do
               end do
            end do
            !$omp end parallel do
         end if
      end do

      ! Normalization
      !$omp parallel do schedule(static) private(il0,jl0r,jc3)
      do il0=1,geom%nl0
         do jl0r=1,bpar%nl0r(ib)
            do jc3=1,nam%nc3
               call divide(mpl,avg%blk(ic2a,bpar%nbe)%m11sta(jc3,jl0r,il0),m11sta(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2a,bpar%nbe)%stasq(jc3,jl0r,il0),stasq(jc3,jl0r,il0))
            end do
         end do
      end do
      !$omp end parallel do
   end if

   ! Update
   call mpl%prog_print(ic2a+1)
end do
call mpl%prog_final

end subroutine avg_compute_bwavg_hyb

!----------------------------------------------------------------------
! Subroutine: avg_compute_bwavg_deh
! Purpose: compute dual-ensemble hybrid block-averaged statistics
!----------------------------------------------------------------------
subroutine avg_compute_bwavg_deh(avg,mpl,nam,geom,bpar,samp,avg_wgt)

implicit none

! Passed variables
class(avg_type),intent(inout) :: avg ! Averaged statistics
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
type(samp_type),intent(in) :: samp   ! Sampling
type(avg_type),intent(in) :: avg_wgt ! Averaged statistics for weights

! Local variables
integer :: ib,ic2a,il0,jl0r,jc3
real(kind_real) :: bwgtsq
real(kind_real) :: m11lrm11(nam%nc3,bpar%nl0rmax,geom%nl0),m11lrm11asy(nam%nc3,bpar%nl0rmax,geom%nl0)

! Initialization
write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(bpar%nbe)),':'
call mpl%flush(.false.)
call mpl%prog_init(samp%nc2a+1)

do ic2a=0,samp%nc2a
   if ((ic2a==0).or.nam%local_diag) then
      ! Initialization
      avg%blk(ic2a,bpar%nbe)%m11lrm11 = 0.0
      m11lrm11 = 0.0
      avg%blk(ic2a,bpar%nbe)%m11lrm11asy = 0.0
      m11lrm11asy = 0.0

      ! Block averages
      do ib=1,bpar%nb
         if (bpar%avg_block(ib)) then
            !$omp parallel do schedule(static) private(il0,jl0r,bwgtsq,jc3)
            do il0=1,geom%nl0
               do jl0r=1,bpar%nl0r(ib)
                  ! Weight
                  if (avg_wgt%blk(0,ib)%m2m2asy(1,jl0r,il0)>0.0) then
                     bwgtsq = 1.0/avg_wgt%blk(0,ib)%m2m2asy(1,jl0r,il0)
                  else
                     bwgtsq = 0.0
                  end if

                  ! Compute sum
                  do jc3=1,nam%nc3
                     call add(mpl,avg%blk(ic2a,ib)%m11lrm11(jc3,jl0r,il0),avg%blk(ic2a,bpar%nbe)%m11lrm11(jc3,jl0r,il0), &
                   & m11lrm11(jc3,jl0r,il0),bwgtsq)
                     call add(mpl,avg%blk(ic2a,ib)%m11lrm11asy(jc3,jl0r,il0),avg%blk(ic2a,bpar%nbe)%m11lrm11asy(jc3,jl0r,il0), &
                   & m11lrm11asy(jc3,jl0r,il0),bwgtsq)
                  end do
               end do
            end do
            !$omp end parallel do
         end if
      end do

      ! Normalization
      !$omp parallel do schedule(static) private(il0,jl0r,jc3)
      do il0=1,geom%nl0
         do jl0r=1,bpar%nl0r(ib)
            do jc3=1,nam%nc3
               call divide(mpl,avg%blk(ic2a,bpar%nbe)%m11lrm11(jc3,jl0r,il0),m11lrm11(jc3,jl0r,il0))
               call divide(mpl,avg%blk(ic2a,bpar%nbe)%m11lrm11asy(jc3,jl0r,il0),m11lrm11asy(jc3,jl0r,il0))
            end do
         end do
      end do
      !$omp end parallel do
   end if

   ! Update
   call mpl%prog_print(ic2a+1)
end do
call mpl%prog_final

end subroutine avg_compute_bwavg_deh

end module type_avg
