!----------------------------------------------------------------------
! Module: type_lct
! Purpose: LCT data derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_lct

use tools_kinds, only: kind_real
use type_bpar, only: bpar_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_lct_blk, only: lct_blk_type
use type_mom, only: mom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type
use type_rng, only: rng_type

implicit none

! LCT data derived type
type lct_type
   type(samp_type) :: samp                  ! Sampling
   type(mom_type) :: mom                    ! Moments
   type(lct_blk_type),allocatable :: blk(:) ! LCT blocks
   logical :: allocated                     ! Allocation flag
contains
   procedure :: alloc => lct_alloc
   procedure :: partial_dealloc => lct_partial_dealloc
   procedure :: dealloc => lct_dealloc
   procedure :: run_lct => lct_run_lct
   procedure :: compute => lct_compute
   procedure :: filter => lct_filter
   procedure :: interp => lct_interp
   procedure :: write => lct_write
   procedure :: write_cor => lct_write_cor
end type lct_type

private
public :: lct_type

contains

!----------------------------------------------------------------------
! Subroutine: lct_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine lct_alloc(lct,nam,geom,bpar)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct ! LCT
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters

! Local variables
integer :: ib

! Allocation
allocate(lct%blk(bpar%nb))
do ib=1,bpar%nb
   if (bpar%diag_block(ib)) call lct%blk(ib)%alloc(nam,geom,bpar,lct%samp,ib)
end do

! Update allocation flag
lct%allocated = .true.

end subroutine lct_alloc

!----------------------------------------------------------------------
! Subroutine: lct_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine lct_partial_dealloc(lct)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct ! LCT

! Local variables
integer :: ib

! Release memory
if (allocated(lct%blk)) then
   do ib=1,size(lct%blk)
      call lct%blk(ib)%partial_dealloc
   end do
end if

end subroutine lct_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: lct_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine lct_dealloc(lct)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct ! LCT

! Local variables
integer :: ib

! Release memory
if (allocated(lct%blk)) then
   do ib=1,size(lct%blk)
      call lct%blk(ib)%dealloc
   end do
   deallocate(lct%blk)
end if

! Update allocation flag
lct%allocated = .false.

end subroutine lct_dealloc

!----------------------------------------------------------------------
! Subroutine: lct_run_lct
! Purpose: LCT driver
!----------------------------------------------------------------------
subroutine lct_run_lct(lct,mpl,rng,nam,geom,bpar,io,ens)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct ! LCT
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(rng_type),intent(inout) :: rng  ! Random number generator
type(nam_type),intent(inout) :: nam  ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters
type(io_type),intent(in) :: io       ! I/O
type(ens_type),intent(inout) :: ens  ! Ensemble

! Set artificially small local radius
nam%local_rad = 1.0e-12

! Setup sampling
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Setup sampling'
call mpl%flush
call lct%samp%setup('lct',mpl,rng,nam,geom,ens)

if (nam%new_mom) then
   ! Compute sample moments
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Compute sample moments'
   call mpl%flush
   call lct%mom%compute(mpl,nam,geom,bpar,lct%samp,ens,'mom')
else
   ! Load sample moments
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Load sample moments'
   call mpl%flush
   call lct%mom%read(mpl,nam,geom,bpar,lct%samp,ens,'mom')
end if

! Compute LCT
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute LCT'
call mpl%flush
call lct%compute(mpl,rng,nam,geom,bpar)

! Filter LCT
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Filter LCT'
call mpl%flush
call lct%filter(mpl,nam,geom,bpar)

! Interpolate LCT
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Interpolate LCT'
call mpl%flush
call lct%interp(mpl,nam,geom,bpar)

if (nam%write_lct) then
   ! Write LCT
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Write LCT'
   call mpl%flush
   call lct%write(mpl,nam,geom,bpar,io)
end if

if (nam%lct_write_cor) then
   ! Write full correlations
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Write full correlations'
   call mpl%flush
   call lct%write_cor(mpl,nam,geom,bpar)
end if

end subroutine lct_run_lct

!----------------------------------------------------------------------
! Subroutine: lct_compute
! Purpose: compute LCT
!----------------------------------------------------------------------
subroutine lct_compute(lct,mpl,rng,nam,geom,bpar)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct ! LCT
type(rng_type),intent(inout) :: rng  ! Random number generator
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(bpar_type),intent(in) :: bpar   ! Block parameters

! Local variables
integer :: ib

! Allocation
call lct%alloc(nam,geom,bpar)

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      call mpl%flush

      ! Compute
      write(mpl%info,'(a10,a)') '','Compute'
      call mpl%flush
      call lct%blk(ib)%compute(mpl,rng,nam,geom,bpar,lct%samp,lct%mom%blk(ib))

      ! Release memory
      call lct%mom%blk(ib)%dealloc
   end if
end do

end subroutine lct_compute

!----------------------------------------------------------------------
! Subroutine: lct_filter
! Purpose: filter LCT
!----------------------------------------------------------------------
subroutine lct_filter(lct,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct    ! LCT
type(mpl_type),intent(inout) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters

! Local variables
integer :: ib

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      call mpl%flush

      ! Filter
      call lct%blk(ib)%filter(mpl,nam,geom,bpar,lct%samp)
   end if
end do

end subroutine lct_filter

!----------------------------------------------------------------------
! Subroutine: lct_interp
! Purpose: interpolate LCT
!----------------------------------------------------------------------
subroutine lct_interp(lct,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct    ! LCT
type(mpl_type),intent(inout) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters

! Local variables
integer :: ib

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
      call mpl%flush

      ! Interpolation
      call lct%blk(ib)%interp(mpl,nam,geom,bpar,lct%samp)
   end if
end do

end subroutine lct_interp

!----------------------------------------------------------------------
! Subroutine: lct_write
! Purpose: write LCT
!----------------------------------------------------------------------
subroutine lct_write(lct,mpl,nam,geom,bpar,io)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct    ! LCT
type(mpl_type),intent(inout) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters
type(io_type),intent(in) :: io          ! I/O

! Local variables
integer :: ib
character(len=1024) :: filename

! Set file name
filename = trim(nam%prefix)//'_lct'

! Write vertical unit
call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      call lct%blk(ib)%write(mpl,nam,geom,bpar,io,filename)
   end if
end do

end subroutine lct_write

!----------------------------------------------------------------------
! Subroutine: lct_write_cor
! Purpose: write full correlation
!----------------------------------------------------------------------
subroutine lct_write_cor(lct,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(lct_type),intent(inout) :: lct    ! LCT
type(mpl_type),intent(inout) :: mpl     ! MPI data
type(nam_type),intent(in) :: nam        ! Namelist
type(geom_type),intent(in) :: geom      ! Geometry
type(bpar_type),intent(in) :: bpar      ! Block parameters

! Local variables
integer :: ib

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      call lct%blk(ib)%write_cor(mpl,nam,geom,bpar,lct%samp)
   end if
end do

end subroutine lct_write_cor

end module type_lct
