!----------------------------------------------------------------------
! Module: type_diag
! Purpose: diagnostic derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_diag

use fckit_mpi_module, only: fckit_mpi_sum
use tools_const, only: reqkm,rad2deg,pi
use tools_fit, only: ver_smooth
use tools_func, only: fit_diag,fit_diag_dble
use tools_kinds, only: kind_real,nc_kind_real
use type_avg, only: avg_type
use type_bpar, only: bpar_type
use type_diag_blk, only: diag_blk_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type
use type_samp, only: samp_type

implicit none

real(kind_real),parameter :: bound = 5.0_kind_real ! Restriction bound applied on local diagnostics with respect to the global diagnostic

! Diagnostic derived type
type diag_type
   character(len=1024) :: prefix               ! Prefix
   integer :: nc2a                             ! Number of local points
   type(diag_blk_type),allocatable :: blk(:,:) ! Diagnostic blocks
contains
   procedure :: alloc => diag_alloc
   procedure :: dealloc => diag_dealloc
   procedure :: write => diag_write
   procedure :: fit_filter => diag_fit_filter
   procedure :: covariance => diag_covariance
   procedure :: correlation => diag_correlation
   procedure :: localization => diag_localization
   procedure :: hybridization => diag_hybridization
   procedure :: dualens => diag_dualens
end type diag_type

private
public :: diag_type

contains

!----------------------------------------------------------------------
! Subroutine: diag_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine diag_alloc(diag,mpl,nam,geom,bpar,samp,prefix,double_fit)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag ! Diagnostic
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(samp_type),intent(in) :: samp     ! Sampling
character(len=*),intent(in) :: prefix  ! Block prefix
logical,intent(in) :: double_fit       ! Double fit

! Local variables
integer :: ib,ic2a

! Number of local points
if (nam%local_diag) then
   diag%nc2a = samp%nc2a
else
   diag%nc2a = 0
end if

! Prefix
diag%prefix = trim(prefix)

! Allocation
allocate(diag%blk(0:diag%nc2a,bpar%nbe))
do ib=1,bpar%nbe
   if (bpar%diag_block(ib)) then
      do ic2a=0,diag%nc2a
         if (bpar%b_to_v1(ib)==0) then
           ! Common block
           call diag%blk(ic2a,ib)%alloc(mpl,nam,geom,bpar,samp,ic2a,ib,prefix,double_fit.and.any(nam%double_fit(1:nam%nv)))
         else
           ! Specific block
           call diag%blk(ic2a,ib)%alloc(mpl,nam,geom,bpar,samp,ic2a,ib,prefix,double_fit.and.nam%double_fit(bpar%b_to_v1(ib)))
         end if
      end do
   end if
end do

end subroutine diag_alloc

!----------------------------------------------------------------------
! Subroutine: diag_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine diag_dealloc(diag)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag ! Diagnostic

! Local variables
integer :: ib,ic2a

! Release memory
if (allocated(diag%blk)) then
   do ib=1,size(diag%blk,2)
     do ic2a=0,size(diag%blk,1)-1
       call diag%blk(ic2a,ib)%dealloc
     end do
   end do
   deallocate(diag%blk)
end if

end subroutine diag_dealloc

!----------------------------------------------------------------------
! Subroutine: diag_write
! Purpose: write
!----------------------------------------------------------------------
subroutine diag_write(diag,mpl,nam,geom,bpar,io,samp)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag ! Diagnostic
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(io_type),intent(in) :: io         ! I/O
type(samp_type),intent(in) :: samp     ! Sampling

! Local variables
integer :: ib,i,ic2,ic0,il0,il0i,iproc,ic2a,ildw,n
real(kind_real),allocatable :: fld_c2a(:,:),fld_c2b(:,:),fld_c0a(:,:)
character(len=2*1024+12) :: filename
character(len=1024),parameter :: subr = 'diag_write'

write(mpl%info,'(a7,a)') '','Write diagnostic'
call mpl%flush

if (mpl%main) then
   filename = trim(nam%prefix)//'_diag'
   do ib=1,bpar%nbe
      if (bpar%diag_block(ib)) call diag%blk(0,ib)%write(mpl,nam,geom,bpar,filename)
   end do
end if

if (nam%local_diag) then
   ! Allocation
   allocate(fld_c2a(samp%nc2a,geom%nl0))
   allocate(fld_c2b(samp%nc2b,geom%nl0))
   allocate(fld_c0a(geom%nc0a,geom%nl0))

   do ib=1,bpar%nbe
      if (bpar%diag_block(ib)) then
         filename = trim(nam%prefix)//'_local_diag_'//trim(diag%prefix)
         call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)
         n = 1
         if (bpar%fit_block(ib).and.(trim(diag%prefix)/='cov')) then
            n = n+2
            if (diag%blk(0,ib)%double_fit) n = n+2
         end if
         do i=1,n
            ! Copy data
            do ic2a=1,samp%nc2a
               if (i==1) then
                  fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%coef_ens
               elseif (i==2) then
                  fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%fit_rh*reqkm
               elseif (i==3) then
                  fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%fit_rv
               elseif (i==4) then
                  fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%fit_rv_rfac
               elseif (i==5) then
                  fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%fit_rv_coef
               end if
            end do
   
            ! Interpolation
            call samp%com_AB%ext(mpl,geom%nl0,fld_c2a,fld_c2b)
            do il0=1,geom%nl0
               il0i = min(il0,geom%nl0i)
               call samp%h(il0i)%apply(mpl,fld_c2b(:,il0),fld_c0a(:,il0))
            end do
   
            ! Write fields
            if (i==1) then
               call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_coef_ens',fld_c0a)
            elseif (i==2) then
               call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_fit_rh',fld_c0a)
            elseif (i==3) then
               call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_fit_rv',fld_c0a)
            elseif (i==4) then
               call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_fit_rv_rfac',fld_c0a)
            elseif (i==5) then
               call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_fit_rv_coef',fld_c0a)
            end if
         end do
      end if
   end do

   ! Release memory
   deallocate(fld_c2a)
   deallocate(fld_c2b)
   deallocate(fld_c0a)
end if

do ildw=1,nam%nldwv
   ic0 = samp%ldwv_to_c0(ildw)
   if (geom%mask_hor_c0(ic0)) then
      iproc = geom%c0_to_proc(ic0)
      if (mpl%myproc==iproc) then
         ! Build file name
         filename = trim(nam%prefix)//'_diag_'//trim(nam%name_ldwv(ildw))

         ! Find diagnostic point
         do ic2a=1,samp%nc2a
            ic2 = samp%c2a_to_c2(ic2a)
            if (samp%c2_to_c0(ic2)==ic0) then
               do ib=1,bpar%nbe
                  if (bpar%diag_block(ib)) call diag%blk(ic2a,ib)%write(mpl,nam,geom,bpar,filename)
               end do
            end if
         end do
      end if
   else
      call mpl%warning(subr,'missing local profile '//trim(nam%name_ldwv(ildw)))
   end if
end do

end subroutine diag_write

!----------------------------------------------------------------------
! Subroutine: diag_fit_filter
! Purpose: filter fit diagnostics
!----------------------------------------------------------------------
subroutine diag_fit_filter(diag,mpl,nam,geom,bpar,samp)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag ! Diagnostic
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(samp_type),intent(in) :: samp     ! Sampling

! Local variables
integer :: ib,il0,ic2a
real(kind_real) :: rmse,norm,rmse_tot,norm_tot
real(kind_real) :: coef_ens_c2a(samp%nc2a,geom%nl0)
real(kind_real),allocatable :: rh_c2a(:,:),rv_c2a(:,:),rv_rfac_c2a(:,:),rv_coef_c2a(:,:)

if (nam%local_diag.and.(nam%diag_rhflt>0.0)) then
   ! Horizontal filtering
   write(mpl%info,'(a7,a)') '','Horizontal filtering'
   call mpl%flush

   do ib=1,bpar%nbe
      if (bpar%diag_block(ib)) then
         ! Allocation
         if (bpar%fit_block(ib)) then
            allocate(rh_c2a(samp%nc2a,geom%nl0))
            allocate(rv_c2a(samp%nc2a,geom%nl0))
            if (diag%blk(0,ib)%double_fit) then
               allocate(rv_rfac_c2a(samp%nc2a,geom%nl0))
               allocate(rv_coef_c2a(samp%nc2a,geom%nl0))
            end if
         end if
      
         do il0=1,geom%nl0
            do ic2a=1,samp%nc2a
               ! Copy data
               coef_ens_c2a(ic2a,il0) = diag%blk(ic2a,ib)%coef_ens(il0)
               if (bpar%fit_block(ib)) then
                  rh_c2a(ic2a,il0) = diag%blk(ic2a,ib)%fit_rh(il0)
                  rv_c2a(ic2a,il0) = diag%blk(ic2a,ib)%fit_rv(il0)
                  if (diag%blk(0,ib)%double_fit) then
                     rv_rfac_c2a(ic2a,il0) = diag%blk(ic2a,ib)%fit_rv_rfac(il0)
                     rv_coef_c2a(ic2a,il0) = diag%blk(ic2a,ib)%fit_rv_coef(il0)
                  end if
               end if
      
               ! Apply bounds relatively to the global value
               if (mpl%msv%isnot(coef_ens_c2a(ic2a,il0)).and.mpl%msv%isnot(diag%blk(0,ib)%coef_ens(il0))) then
                  if ((coef_ens_c2a(ic2a,il0)<diag%blk(0,ib)%coef_ens(il0)/bound) &
                & .or.(coef_ens_c2a(ic2a,il0)>diag%blk(0,ib)%coef_ens(il0)*bound)) coef_ens_c2a(ic2a,il0) = mpl%msv%valr
               end if
               if (bpar%fit_block(ib)) then
                  if (mpl%msv%isnot(rh_c2a(ic2a,il0)).and.mpl%msv%isnot(diag%blk(0,ib)%fit_rh(il0))) then
                     if ((rh_c2a(ic2a,il0)<diag%blk(0,ib)%fit_rh(il0)/bound) &
                   & .or.(rh_c2a(ic2a,il0)>diag%blk(0,ib)%fit_rh(il0)*bound)) rh_c2a(ic2a,il0) = mpl%msv%valr
                  end if
                  if (mpl%msv%isnot(rv_c2a(ic2a,il0)).and.mpl%msv%isnot(diag%blk(0,ib)%fit_rv(il0))) then
                     if ((rv_c2a(ic2a,il0)<diag%blk(0,ib)%fit_rv(il0)/bound) &
                   & .or.(rv_c2a(ic2a,il0)>diag%blk(0,ib)%fit_rv(il0)*bound)) rv_c2a(ic2a,il0) = mpl%msv%valr
                  end if
               end if
            end do
      
            ! Median filter to remove extreme values
            call samp%diag_filter(mpl,nam,'median',nam%diag_rhflt,coef_ens_c2a(:,il0))
            if (bpar%fit_block(ib)) then
               call samp%diag_filter(mpl,nam,'median',nam%diag_rhflt,rh_c2a(:,il0))
               call samp%diag_filter(mpl,nam,'median',nam%diag_rhflt,rv_c2a(:,il0))
               if (diag%blk(0,ib)%double_fit) then
                  call samp%diag_filter(mpl,nam,'median',nam%diag_rhflt,rv_rfac_c2a(:,il0))
                  call samp%diag_filter(mpl,nam,'median',nam%diag_rhflt,rv_coef_c2a(:,il0))
               end if
            end if
      
            ! Average filter to smooth data
            call samp%diag_filter(mpl,nam,'average',nam%diag_rhflt,coef_ens_c2a(:,il0))
            if (bpar%fit_block(ib)) then
               call samp%diag_filter(mpl,nam,'average',nam%diag_rhflt,rh_c2a(:,il0))
               call samp%diag_filter(mpl,nam,'average',nam%diag_rhflt,rv_c2a(:,il0))
               if (diag%blk(0,ib)%double_fit) then
                  call samp%diag_filter(mpl,nam,'average',nam%diag_rhflt,rv_rfac_c2a(:,il0))
                  call samp%diag_filter(mpl,nam,'average',nam%diag_rhflt,rv_coef_c2a(:,il0))
               end if
            end if
      
            ! Fill missing values
            call samp%diag_fill(mpl,nam,coef_ens_c2a(:,il0))
            if (bpar%fit_block(ib)) then
               call samp%diag_fill(mpl,nam,rh_c2a(:,il0))
               call samp%diag_fill(mpl,nam,rv_c2a(:,il0))
               if (diag%blk(0,ib)%double_fit) then
                  call samp%diag_fill(mpl,nam,rv_rfac_c2a(:,il0))
                  call samp%diag_fill(mpl,nam,rv_coef_c2a(:,il0))
               end if
            end if
      
            ! Copy data
            do ic2a=1,samp%nc2a
               diag%blk(ic2a,ib)%coef_ens(il0) = coef_ens_c2a(ic2a,il0)
               if (bpar%fit_block(ib)) then
                  diag%blk(ic2a,ib)%fit_rh(il0) = rh_c2a(ic2a,il0)
                  diag%blk(ic2a,ib)%fit_rv(il0) = rv_c2a(ic2a,il0)
                  if (diag%blk(0,ib)%double_fit) then
                     diag%blk(ic2a,ib)%fit_rv_rfac(il0) = rv_rfac_c2a(ic2a,il0)
                     diag%blk(ic2a,ib)%fit_rv_coef(il0) = rv_coef_c2a(ic2a,il0)
                  end if
               end if
            end do
         end do
      
         ! Release memory
         if (bpar%fit_block(ib)) then
            deallocate(rh_c2a)
            deallocate(rv_c2a)
            if (diag%blk(0,ib)%double_fit) then
               deallocate(rv_rfac_c2a)
               deallocate(rv_coef_c2a)
            end if
         end if
      end if
   end do
end if

if (nam%diag_rvflt>0.0) then
   ! Vertical filtering
   write(mpl%info,'(a7,a)') '','Vertical filtering'
   call mpl%flush

   do ib=1,bpar%nbe
      if (bpar%diag_block(ib)) then
         ! Vertical filtering
         do ic2a=0,diag%nc2a
            call ver_smooth(mpl,geom%nl0,geom%vunitavg,nam%diag_rvflt,diag%blk(ic2a,ib)%coef_ens)
            if (bpar%fit_block(ib)) then
               call ver_smooth(mpl,geom%nl0,geom%vunitavg,nam%diag_rvflt,diag%blk(ic2a,ib)%fit_rh)
               call ver_smooth(mpl,geom%nl0,geom%vunitavg,nam%diag_rvflt,diag%blk(ic2a,ib)%fit_rv)
               if (diag%blk(0,ib)%double_fit) then
                  call ver_smooth(mpl,geom%nl0,geom%vunitavg,nam%diag_rvflt,diag%blk(ic2a,ib)%fit_rv_rfac)
                  call ver_smooth(mpl,geom%nl0,geom%vunitavg,nam%diag_rvflt,diag%blk(ic2a,ib)%fit_rv_coef)
               end if
            end if
         end do
      end if
   end do
end if

do ib=1,bpar%nbe
   if (bpar%diag_block(ib).and.bpar%fit_block(ib)) then
      ! Rebuild fit
      do ic2a=0,diag%nc2a
         if (diag%blk(0,ib)%double_fit) then
            call fit_diag_dble(mpl,nam%fit_type,nam%nc3,bpar%nl0r(ib),geom%nl0,bpar%l0rl0b_to_l0(:,:,ib),geom%disth, &
          & diag%blk(ic2a,ib)%distv,diag%blk(ic2a,ib)%coef_ens,diag%blk(ic2a,ib)%fit_rh,diag%blk(ic2a,ib)%fit_rv, &
          & diag%blk(ic2a,ib)%fit_rv_rfac,diag%blk(ic2a,ib)%fit_rv_coef,diag%blk(ic2a,ib)%fit)
         else
            call fit_diag(mpl,nam%fit_type,nam%nc3,bpar%nl0r(ib),geom%nl0,bpar%l0rl0b_to_l0(:,:,ib),geom%disth, &
          & diag%blk(ic2a,ib)%distv,diag%blk(ic2a,ib)%coef_ens,diag%blk(ic2a,ib)%fit_rh,diag%blk(ic2a,ib)%fit_rv, &
          & diag%blk(ic2a,ib)%fit)
         end if
      end do
   
      ! Compute RMSE
      rmse = 0.0
      norm = 0.0
      do ic2a=0,diag%nc2a
         rmse = rmse+sum(abs(diag%blk(ic2a,ib)%fit-diag%blk(ic2a,ib)%raw),mask=mpl%msv%isnot(diag%blk(ic2a,ib)%raw))
         norm = norm+real(count(mpl%msv%isnot(diag%blk(ic2a,ib)%raw)),kind_real)
      end do
      call mpl%f_comm%allreduce(rmse,rmse_tot,fckit_mpi_sum())
      call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
      if (norm_tot>0.0) rmse_tot = sqrt(rmse_tot/norm_tot)
      write(mpl%info,'(a10,a,a,a,e15.8,a,i8,a)') '','Fit RMSE for block ',trim(bpar%blockname(ib)),': ',rmse_tot, &
    & ' for ',int(norm_tot),' diagnostic points'
      call mpl%flush
   end if
end do

end subroutine diag_fit_filter

!----------------------------------------------------------------------
! Subroutine: diag_covariance
! Purpose: compute covariance
!----------------------------------------------------------------------
subroutine diag_covariance(diag,mpl,nam,geom,bpar,io,samp,avg,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag ! Diagnostic
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(io_type),intent(in) :: io         ! I/O
type(samp_type),intent(in) :: samp     ! Sampling
type(avg_type),intent(in) :: avg       ! Averaged statistics
character(len=*),intent(in) :: prefix  ! Diagnostic prefix

! Local variables
integer :: ib,ic2a,il0

! Allocation
call diag%alloc(mpl,nam,geom,bpar,samp,prefix,.false.)

do ib=1,bpar%nb
   if (bpar%diag_block(ib)) then
      write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib))
      call mpl%flush

      ! Copy covariance
      do ic2a=0,diag%nc2a
         diag%blk(ic2a,ib)%raw = avg%blk(ic2a,ib)%m11
         diag%blk(ic2a,ib)%valid = avg%blk(ic2a,ib)%nc1a
      end do

      ! Copy (potentially filtered) variance
      do ic2a=0,diag%nc2a
         if (nam%var_filter) then
            diag%blk(ic2a,ib)%coef_ens = sum(avg%blk(ic2a,ib)%m2flt,dim=2)/real(avg%nsub,kind_real)
         else
            diag%blk(ic2a,ib)%coef_ens = sum(avg%blk(ic2a,ib)%m2,dim=2)/real(avg%nsub,kind_real)
         end if
      end do

      ! Print results
      do il0=1,geom%nl0
         if (mpl%msv%isnot(diag%blk(0,ib)%raw(1,bpar%il0rz(il0,ib),il0))) then
            write(mpl%info,'(a13,a,i3,a,a,e9.2,a)') '','Level: ',nam%levs(il0),' ~> cov. at class zero: ',trim(mpl%peach), &
          & diag%blk(0,ib)%coef_ens(il0),trim(mpl%black)
            call mpl%flush
         end if
      end do
   end if
end do

! Write
if (nam%write_hdiag) call diag%write(mpl,nam,geom,bpar,io,samp)

end subroutine diag_covariance

!----------------------------------------------------------------------
! Subroutine: diag_correlation
! Purpose: compute correlation
!----------------------------------------------------------------------
subroutine diag_correlation(diag,mpl,rng,nam,geom,bpar,io,samp,avg,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag ! Diagnostic
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(io_type),intent(in) :: io         ! I/O
type(samp_type),intent(in) :: samp     ! Sampling
type(avg_type),intent(in) :: avg       ! Averaged statistics
character(len=*),intent(in) :: prefix  ! Diagnostic prefix

! Local variables
integer :: ib,ic2a,il0

! Allocation
call diag%alloc(mpl,nam,geom,bpar,samp,prefix,.true.)

do ib=1,bpar%nbe
   if (bpar%diag_block(ib)) then
      ! Initialization
      write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)
      call mpl%prog_init(diag%nc2a+1)

      do ic2a=0,diag%nc2a
         ! Copy correlation
         diag%blk(ic2a,ib)%raw = avg%blk(ic2a,ib)%cor
         diag%blk(ic2a,ib)%valid = avg%blk(ic2a,ib)%nc1a_cor

         ! Set variance to 1
         diag%blk(ic2a,ib)%coef_ens = 1.0

         ! Fitting
         if (bpar%fit_block(ib)) call diag%blk(ic2a,ib)%fitting(mpl,rng,nam,geom,bpar,samp,.false.)

         ! Update
         call mpl%prog_print(ic2a+1)
      end do
      call mpl%prog_final

      ! Print results
      do il0=1,geom%nl0
         if (bpar%fit_block(ib)) then
            if (mpl%msv%isnot(diag%blk(0,ib)%fit_rh(il0))) then
               write(mpl%info,'(a13,a,i3,a,a,f10.2,a,f10.2,a)') '','Level: ',nam%levs(il0),' ~> cor. support radii: ', &
             & trim(mpl%aqua),diag%blk(0,ib)%fit_rh(il0)*reqkm,trim(mpl%black)//' km  / '//trim(mpl%aqua), &
             & diag%blk(0,ib)%fit_rv(il0),trim(mpl%black)//' vert. unit'
               call mpl%flush
               if (diag%blk(0,ib)%double_fit) then
                  write(mpl%info,'(a27,a,a,f10.2,a,f10.2,a)') '','cor. double fit:    ',trim(mpl%aqua), &
                & diag%blk(0,ib)%fit_rv_rfac(il0),trim(mpl%black)//' / '//trim(mpl%aqua),diag%blk(0,ib)%fit_rv_coef(il0), &
                & trim(mpl%black)
                  call mpl%flush
               end if
            end if
         end if
      end do
   end if
end do

! Filtering
call diag%fit_filter(mpl,nam,geom,bpar,samp)

! Write
if (nam%write_hdiag) call diag%write(mpl,nam,geom,bpar,io,samp)

end subroutine diag_correlation

!----------------------------------------------------------------------
! Subroutine: diag_localization
! Purpose: compute diagnostic localization
!----------------------------------------------------------------------
subroutine diag_localization(diag,mpl,rng,nam,geom,bpar,io,samp,avg,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag ! Diagnostic
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(io_type),intent(in) :: io         ! I/O
type(samp_type),intent(in) :: samp     ! Sampling
type(avg_type),intent(in) :: avg       ! Averaged statistics
character(len=*),intent(in) :: prefix  ! Block prefix

! Local variables
integer :: ib,ic2a,il0

! Allocation
call diag%alloc(mpl,nam,geom,bpar,samp,prefix,.false.)

do ib=1,bpar%nbe
   if (bpar%diag_block(ib)) then
      ! Initialization
      write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)
      call mpl%prog_init(diag%nc2a+1)

      do ic2a=0,diag%nc2a
         ! Compute localization
         call diag%blk(ic2a,ib)%localization(mpl,geom,bpar,avg%blk(ic2a,ib))

         ! Normalization
         call diag%blk(ic2a,ib)%normalization(geom,bpar)

         ! Fitting
         if (bpar%fit_block(ib)) call diag%blk(ic2a,ib)%fitting(mpl,rng,nam,geom,bpar,samp)

         ! Update
         call mpl%prog_print(ic2a+1)
      end do
      call mpl%prog_final

      ! Print results
      do il0=1,geom%nl0
         select case (trim(nam%method))
         case ('loc','hyb-avg','hyb-rnd','dual-ens')
            if (mpl%msv%isnot(diag%blk(0,ib)%coef_ens(il0))) then
               write(mpl%info,'(a13,a,i3,a,a,f10.2,a)') '','Level: ',nam%levs(il0),' ~> loc. at class zero: ', &
             & trim(mpl%peach),diag%blk(0,ib)%coef_ens(il0),trim(mpl%black)
               call mpl%flush
            end if
         end select
         if (bpar%fit_block(ib)) then
            if (mpl%msv%isnot(diag%blk(0,ib)%fit_rh(il0))) then
               write(mpl%info,'(a27,a,a,f10.2,a,f10.2,a)') '','loc. support radii: ',trim(mpl%aqua), &
             & diag%blk(0,ib)%fit_rh(il0)*reqkm,trim(mpl%black)//' km  / '//trim(mpl%aqua),diag%blk(0,ib)%fit_rv(il0), &
             & trim(mpl%black)//' vert. unit'
               call mpl%flush
            end if
         end if
      end do
   end if
end do

! Filtering
call diag%fit_filter(mpl,nam,geom,bpar,samp)

! Write
if (nam%write_hdiag) call diag%write(mpl,nam,geom,bpar,io,samp)

end subroutine diag_localization

!----------------------------------------------------------------------
! Subroutine: diag_hybridization
! Purpose: compute diagnostic hybridization
!----------------------------------------------------------------------
subroutine diag_hybridization(diag,mpl,rng,nam,geom,bpar,io,samp,avg,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag ! Diagnostic (localization)
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar     ! Block parameters
type(io_type),intent(in) :: io         ! I/O
type(samp_type),intent(in) :: samp     ! Sampling
type(avg_type),intent(in) :: avg       ! Averaged statistics
character(len=*),intent(in) :: prefix  ! Diagnostic prefix

! Local variables
integer :: ib,ic2a,il0

! Allocation
call diag%alloc(mpl,nam,geom,bpar,samp,prefix,.false.)

do ib=1,bpar%nbe
   if (bpar%diag_block(ib)) then
      ! Initialization
      write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)
      call mpl%prog_init(diag%nc2a+1)

      do ic2a=0,diag%nc2a
         ! Compute hybridization
         call diag%blk(ic2a,ib)%hybridization(mpl,geom,bpar,avg%blk(ic2a,ib))

         ! Normalization
         call diag%blk(ic2a,ib)%normalization(geom,bpar)

         ! Fitting
         if (bpar%fit_block(ib)) call diag%blk(ic2a,ib)%fitting(mpl,rng,nam,geom,bpar,samp)

         ! Update
         call mpl%prog_print(ic2a+1)
      end do
      call mpl%prog_final

      ! Print results
      do il0=1,geom%nl0
         if (mpl%msv%isnot(diag%blk(0,ib)%coef_ens(il0))) then
            write(mpl%info,'(a13,a,i3,a4,a20,a,f10.2,a)') '','Level: ',nam%levs(il0),' ~> ','loc. at class zero: ', &
          & trim(mpl%peach),diag%blk(0,ib)%coef_ens(il0),trim(mpl%black)
            call mpl%flush
         end if
         if (bpar%fit_block(ib)) then
            if (mpl%msv%isnot(diag%blk(0,ib)%fit_rh(il0))) then
               write(mpl%info,'(a27,a,a,f10.2,a,f10.2,a)') '','loc. support radii: ',trim(mpl%aqua), &
             & diag%blk(0,ib)%fit_rh(il0)*reqkm,trim(mpl%black)//' km  / '//trim(mpl%aqua),diag%blk(0,ib)%fit_rv(il0), &
             & trim(mpl%black)//' vert. unit'
               call mpl%flush
            end if
         end if
      end do
      if (mpl%msv%isnot(diag%blk(0,ib)%coef_sta)) write(mpl%info,'(a13,a,a,f4.2,a)') '', &
    & 'Static coeff.:                          ',trim(mpl%purple),diag%blk(0,ib)%coef_sta,trim(mpl%black)
      call mpl%flush
   end if
end do

! Filtering
call diag%fit_filter(mpl,nam,geom,bpar,samp)

! Write
if (nam%write_hdiag) call diag%write(mpl,nam,geom,bpar,io,samp)

end subroutine diag_hybridization

!----------------------------------------------------------------------
! Subroutine: diag_dualens
! Purpose: compute diagnostic dualens
!----------------------------------------------------------------------
subroutine diag_dualens(diag,mpl,rng,nam,geom,bpar,io,samp,avg,avg_lr,diag_lr,prefix,prefix_lr)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag   ! Diagnostic (localization)
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(rng_type),intent(inout) :: rng      ! Random number generator
type(nam_type),intent(in) :: nam         ! Namelist
type(geom_type),intent(in) :: geom       ! Geometry
type(bpar_type),intent(in) :: bpar       ! Block parameters
type(io_type),intent(in) :: io           ! I/O
type(samp_type),intent(in) :: samp       ! Sampling
type(avg_type),intent(in) :: avg         ! Averaged statistics
type(avg_type),intent(in) :: avg_lr      ! LR averaged statistics
type(diag_type),intent(inout) :: diag_lr ! Diagnostic (LR localization)
character(len=*),intent(in) :: prefix    ! Diagnostic prefix
character(len=*),intent(in) :: prefix_lr ! LR diagnostic prefix

! Local variables
integer :: ib,ic2a,il0

! Allocation
call diag%alloc(mpl,nam,geom,bpar,samp,prefix,.false.)
call diag_lr%alloc(mpl,nam,geom,bpar,samp,prefix_lr,.false.)

do ib=1,bpar%nbe
   if (bpar%diag_block(ib)) then
      ! Initialization
      write(mpl%info,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush(.false.)
      call mpl%prog_init(diag%nc2a+1)

      do ic2a=0,diag%nc2a
         ! Compute dualens
         call diag%blk(ic2a,ib)%dualens(mpl,geom,bpar,avg%blk(ic2a,ib),avg_lr%blk(ic2a,ib),diag_lr%blk(ic2a,ib))

         ! Normalization
         call diag%blk(ic2a,ib)%normalization(geom,bpar)
         call diag_lr%blk(ic2a,ib)%normalization(geom,bpar)

         ! Fitting
         if (bpar%fit_block(ib)) then
            call diag%blk(ic2a,ib)%fitting(mpl,rng,nam,geom,bpar,samp)
            call diag_lr%blk(ic2a,ib)%fitting(mpl,rng,nam,geom,bpar,samp)
         end if

         ! Update
         call mpl%prog_print(ic2a+1)
      end do
      call mpl%prog_final

      ! Print results
      do il0=1,geom%nl0
         if (mpl%msv%isnot(diag%blk(0,ib)%coef_ens(il0))) then
            write(mpl%info,'(a13,a,i3,a4,a21,a,f10.2,a)') '','Level: ',nam%levs(il0),' ~> ','loc. at class zero (HR): ', &
          & trim(mpl%peach),diag%blk(0,ib)%coef_ens(il0),trim(mpl%black)
            call mpl%flush
         end if
         if (mpl%msv%isnot(diag%blk(0,ib)%coef_ens(il0))) then
            write(mpl%info,'(a27,a,a,f10.2,a)') '','loc. at class zero (LR): ',trim(mpl%peach), &
          & diag_lr%blk(0,ib)%coef_ens(il0),trim(mpl%black)
            call mpl%flush
         end if
         if (bpar%fit_block(ib)) then
            if (mpl%msv%isnot(diag%blk(0,ib)%fit_rh(il0))) then
               write(mpl%info,'(a27,a,a,f10.2,a,f10.2,a)') '','loc. support radii (HR): ',trim(mpl%aqua), &
             & diag%blk(0,ib)%fit_rh(il0)*reqkm,trim(mpl%black)//' km  / '//trim(mpl%aqua),diag_lr%blk(0,ib)%fit_rv(il0), &
             & trim(mpl%black)//' vert. unit'
               call mpl%flush
            end if
            if (mpl%msv%isnot(diag_lr%blk(0,ib)%fit_rh(il0))) then
               write(mpl%info,'(a27,a,a,f10.2,a,f10.2,a)') '','loc. support radii (LR): ',trim(mpl%aqua), &
             & diag_lr%blk(0,ib)%fit_rh(il0)*reqkm,trim(mpl%black)//' km  / '//trim(mpl%aqua),diag_lr%blk(0,ib)%fit_rv(il0), &
             & trim(mpl%black)//' vert. unit'
               call mpl%flush
            end if
         end if
      end do
   end if
end do

! Filtering
call diag%fit_filter(mpl,nam,geom,bpar,samp)
call diag_lr%fit_filter(mpl,nam,geom,bpar,samp)

! Write
if (nam%write_hdiag) call diag%write(mpl,nam,geom,bpar,io,samp)
if (nam%write_hdiag) call diag_lr%write(mpl,nam,geom,bpar,io,samp)

end subroutine diag_dualens

end module type_diag
