!----------------------------------------------------------------------
! Module: type_lct_blk
! Purpose: LCT data derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_lct_blk

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_status
use netcdf
!$ use omp_lib
use tools_const, only: reqkm,rad2deg
use tools_func, only: gau2gc,fit_lct,lct_d2h,check_cond,Dmin
use tools_kinds, only: kind_real,nc_kind_real
use tools_repro, only: rth,inf,sup
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_minim, only: minim_type
use type_mom_blk, only: mom_blk_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type
use type_samp, only: samp_type

implicit none

! LCT block data derived type
type lct_blk_type
   ! Attributes
   integer :: ib                                    ! Block index
   integer :: nscales                               ! Number of LCT scales
   character(len=1024) :: name                      ! Name

   ! Correlation/variances
   real(kind_real),allocatable :: raw(:,:,:,:)      ! Raw correlations

   ! Diffusion data
   real(kind_real),allocatable :: fit(:,:,:,:)      ! Fitted correlations
   real(kind_real),allocatable :: D(:,:,:,:)        ! Diffusion components
   real(kind_real),allocatable :: coef(:,:,:)       ! Multi-scale coefficients
   real(kind_real),allocatable :: qc_c1a(:,:)       ! Quality control on subset Sc1, halo A
   real(kind_real),allocatable :: qc_c0a(:,:)       ! Quality control on subset Sc0, halo A

   ! Filtered diffusion data
   real(kind_real),allocatable :: fit_filt(:,:,:,:) ! Fitted correlations after filtering
   real(kind_real),allocatable :: D_filt(:,:,:,:)   ! Diffusion components after filtering
   real(kind_real),allocatable :: coef_filt(:,:,:)  ! Multi-scale coefficients after filtering

   ! NICAS-related data
   real(kind_real),allocatable :: H11(:,:,:)        ! Local correlation tensor, component 11
   real(kind_real),allocatable :: H22(:,:,:)        ! Local correlation tensor, component 22
   real(kind_real),allocatable :: H33(:,:,:)        ! Local correlation tensor, component 33
   real(kind_real),allocatable :: H12(:,:,:)        ! Local correlation tensor, component 12

   ! BUMP-output data
   real(kind_real),allocatable :: D11(:,:,:)        ! Daley tensor, component 11
   real(kind_real),allocatable :: D22(:,:,:)        ! Daley tensor, component 22
   real(kind_real),allocatable :: D33(:,:,:)        ! Daley tensor, component 33
   real(kind_real),allocatable :: D12(:,:,:)        ! Daley tensor, component 12
   real(kind_real),allocatable :: Dcoef(:,:,:)      ! Tensor coefficient
   real(kind_real),allocatable :: DLh(:,:,:)        ! Tensor length-scale
contains
   procedure :: alloc => lct_blk_alloc
   procedure :: partial_dealloc => lct_blk_partial_dealloc
   procedure :: dealloc => lct_blk_dealloc
   procedure :: compute => lct_blk_compute
   procedure :: filter => lct_blk_filter
   procedure :: interp => lct_blk_interp
   procedure :: write => lct_blk_write
   procedure :: write_cor => lct_blk_write_cor
end type lct_blk_type

private
public :: lct_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: lct_blk_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine lct_blk_alloc(lct_blk,nam,geom,bpar,samp,ib)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! LCT block
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(samp_type),intent(in) :: samp           ! Sampling
integer,intent(in) :: ib                     ! Block index

! Attributes
lct_blk%ib = ib
lct_blk%nscales = nam%lct_nscales

! Allocation
allocate(lct_blk%raw(nam%nc3,bpar%nl0r(ib),samp%nc1a,geom%nl0))
allocate(lct_blk%fit(nam%nc3,bpar%nl0r(ib),samp%nc1a,geom%nl0))
allocate(lct_blk%D(4,lct_blk%nscales,samp%nc1a,geom%nl0))
allocate(lct_blk%coef(lct_blk%nscales,samp%nc1a,geom%nl0))
allocate(lct_blk%qc_c1a(samp%nc1a,geom%nl0))
allocate(lct_blk%qc_c0a(geom%nc0a,geom%nl0))
if (nam%diag_rhflt>0.0) then
   allocate(lct_blk%fit_filt(nam%nc3,bpar%nl0r(ib),samp%nc1a,geom%nl0))
   allocate(lct_blk%D_filt(4,lct_blk%nscales,samp%nc1a,geom%nl0))
   allocate(lct_blk%coef_filt(lct_blk%nscales,samp%nc1a,geom%nl0))
end if
allocate(lct_blk%D11(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%D22(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%D33(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%D12(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%H11(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%H22(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%H33(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%H12(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%Dcoef(geom%nc0a,geom%nl0,lct_blk%nscales))
allocate(lct_blk%DLh(geom%nc0a,geom%nl0,lct_blk%nscales))

end subroutine lct_blk_alloc

!----------------------------------------------------------------------
! Subroutine: lct_blk_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine lct_blk_partial_dealloc(lct_blk)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! LCT block

! Release memory
if (allocated(lct_blk%raw)) deallocate(lct_blk%raw)
if (allocated(lct_blk%fit)) deallocate(lct_blk%fit)
if (allocated(lct_blk%D)) deallocate(lct_blk%D)
if (allocated(lct_blk%coef)) deallocate(lct_blk%coef)
if (allocated(lct_blk%qc_c1a)) deallocate(lct_blk%qc_c1a)
if (allocated(lct_blk%qc_c0a)) deallocate(lct_blk%qc_c0a)
if (allocated(lct_blk%fit_filt)) deallocate(lct_blk%fit_filt)
if (allocated(lct_blk%D_filt)) deallocate(lct_blk%D_filt)
if (allocated(lct_blk%coef_filt)) deallocate(lct_blk%coef_filt)
if (allocated(lct_blk%H11)) deallocate(lct_blk%H11)
if (allocated(lct_blk%H22)) deallocate(lct_blk%H22)
if (allocated(lct_blk%H33)) deallocate(lct_blk%H33)
if (allocated(lct_blk%H12)) deallocate(lct_blk%H12)

end subroutine lct_blk_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: lct_blk_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine lct_blk_dealloc(lct_blk)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! LCT block

! Release memory
call lct_blk%partial_dealloc
if (allocated(lct_blk%D11)) deallocate(lct_blk%D11)
if (allocated(lct_blk%D22)) deallocate(lct_blk%D22)
if (allocated(lct_blk%D33)) deallocate(lct_blk%D33)
if (allocated(lct_blk%D12)) deallocate(lct_blk%D12)
if (allocated(lct_blk%Dcoef)) deallocate(lct_blk%Dcoef)
if (allocated(lct_blk%DLh)) deallocate(lct_blk%DLh)

end subroutine lct_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: lct_blk_compute
! Purpose: compute raw correlation and fit to get LCT components
!----------------------------------------------------------------------
subroutine lct_blk_compute(lct_blk,mpl,rng,nam,geom,bpar,samp,mom_blk)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! LCT block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(rng_type),intent(inout) :: rng          ! Random number generator
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(samp_type),intent(in) :: samp           ! Sampling
type(mom_blk_type),intent(in) :: mom_blk     ! Moments block

! Local variables
integer :: jsub,il0,jl0r,jl0,jc3,ic1a,ic1,ic0,jc0,iscales,icomp
real(kind_real) :: den,dx,dy,dz,distsq,Dhbar,Dvbar,norm
real(kind_real),allocatable :: norm_raw(:,:),Dh(:),Dv(:)
logical :: valid
logical,allocatable :: Dv_valid(:)
type(minim_type) :: minim

! Associate
associate(ib=>lct_blk%ib)

! Initialization
lct_blk%raw = 0.0

do il0=1,geom%nl0
   write(mpl%info,'(a13,a,i3,a)') '','Level ',nam%levs(il0),':'
   call mpl%flush(.false.)
   call mpl%prog_init(samp%nc1a)

   do ic1a=1,samp%nc1a
      ! Global index
      ic1 = samp%c1a_to_c1(ic1a)

      select case (trim(nam%minim_algo))
      case ('hooke')
         ! Hooke parameters
         minim%hooke_rho = 0.5
         minim%hooke_tol = 1.0e-4
         minim%hooke_itermax = 10
      case ('praxis')
         ! Praxis parameters
         minim%praxis_tol = 1.0
         minim%praxis_itermax = 5
      end select

      ! Allocation
      allocate(norm_raw(nam%nc3,bpar%nl0r(ib)))
      allocate(Dh(nam%nc3))
      allocate(Dv(bpar%nl0r(ib)))
      allocate(Dv_valid(bpar%nl0r(ib)))
      minim%nx = lct_blk%nscales*4
      if (lct_blk%nscales>1) minim%nx = minim%nx+lct_blk%nscales-1
      minim%ny = nam%nc3*bpar%nl0r(ib)
      allocate(minim%x(minim%nx))
      allocate(minim%guess(minim%nx))
      allocate(minim%binf(minim%nx))
      allocate(minim%bsup(minim%nx))
      allocate(minim%obs(minim%ny))
      allocate(minim%dxsq(nam%nc3,bpar%nl0r(ib)))
      allocate(minim%dysq(nam%nc3,bpar%nl0r(ib)))
      allocate(minim%dxdy(nam%nc3,bpar%nl0r(ib)))
      allocate(minim%dzsq(nam%nc3,bpar%nl0r(ib)))
      allocate(minim%dmask(nam%nc3,bpar%nl0r(ib)))

      if (samp%c1l0_log(ic1,il0)) then
         ! Initialization
         norm_raw = 0.0

         ! Compute correlation
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            do jc3=1,nam%nc3
               if (samp%c1ac3l0_log(ic1a,jc3,jl0)) then
                  do jsub=1,mom_blk%nsub
                     den = mom_blk%m2_1(ic1a,il0,jsub)*mom_blk%m2_2(ic1a,jc3,jl0,jsub)
                     if (den>0.0) then
                        lct_blk%raw(jc3,jl0r,ic1a,il0) = lct_blk%raw(jc3,jl0r,ic1a,il0)+ &
                                                       & mom_blk%m11(ic1a,jc3,jl0r,il0,jsub)/sqrt(den)
                        norm_raw(jc3,jl0r) = norm_raw(jc3,jl0r)+1.0
                     end if
                  end do
               end if
            end do
         end do

         ! Normalize
         do jl0r=1,bpar%nl0r(ib)
            do jc3=1,nam%nc3
               if (norm_raw(jc3,jl0r)>0.0) lct_blk%raw(jc3,jl0r,ic1a,il0) = lct_blk%raw(jc3,jl0r,ic1a,il0)/norm_raw(jc3,jl0r)
            end do
         end do

         ! Compute deltas
         ic0 = samp%c1_to_c0(ic1)
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            do jc3=1,nam%nc3
               minim%dmask(jc3,jl0r) = samp%c1ac3l0_log(ic1a,jc3,jl0)
               if (minim%dmask(jc3,jl0r)) then
                  jc0 = samp%c1ac3_to_c0(ic1a,jc3)
                  call geom%compute_deltas(ic0,il0,jc0,jl0,dx,dy,dz)
                  minim%dxsq(jc3,jl0r) = dx**2
                  minim%dysq(jc3,jl0r) = dy**2
                  minim%dxdy(jc3,jl0r) = dx*dy
                  minim%dzsq(jc3,jl0r) = dz**2
               end if
            end do
         end do

         ! Approximate homogeneous horizontal length-scale
         Dh = mpl%msv%valr
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            if (il0==jl0) then
               do jc3=1,nam%nc3
                  if (minim%dmask(jc3,jl0r)) then
                     distsq = minim%dxsq(jc3,jl0r)+minim%dysq(jc3,jl0r)
                     if (sup(lct_blk%raw(jc3,jl0r,ic1a,il0),nam%lct_cor_min).and.inf(lct_blk%raw(jc3,jl0r,ic1a,il0), &
                   & 1.0_kind_real).and.(distsq>0.0)) Dh(jc3) = -distsq/(2.0*log(lct_blk%raw(jc3,jl0r,ic1a,il0)))
                  end if
               end do
            end if
         end do
         Dhbar = mpl%msv%valr
         if (count(mpl%msv%isnot(Dh))>0) Dhbar = sum(Dh,mask=mpl%msv%isnot(Dh))/real(count(mpl%msv%isnot(Dh)),kind_real)

         ! Approximate homogeneous vertical length-scale
         Dv = mpl%msv%valr
         Dv_valid = .false.
         jc3 = 1
         do jl0r=1,bpar%nl0r(ib)
            if (minim%dmask(jc3,jl0r)) then
               distsq = minim%dzsq(jc3,jl0r)
               if (inf(abs(lct_blk%raw(jc3,jl0r,ic1a,il0)),1.0_kind_real).and.(distsq>0.0)) then
                  if (lct_blk%raw(jc3,jl0r,ic1a,il0)>0.0) Dv(jl0r) = -distsq/(2.0*log(lct_blk%raw(jc3,jl0r,ic1a,il0)))
                  Dv_valid(jl0r) = .true.
               end if
            end if
         end do
         Dvbar = mpl%msv%valr
         if (bpar%nl0r(ib)>1) then
            if (count(Dv_valid)>0) then
               if (count(mpl%msv%isnot(Dv))>0) then
                  Dvbar = sum(Dv,mask=mpl%msv%isnot(Dv))/real(count(mpl%msv%isnot(Dv)),kind_real)
               else
                  Dvbar = 1.0/rth
               end if
            end if
         else
             Dvbar = 0.0
         end if

         if (mpl%msv%isnot(Dhbar).and.(mpl%msv%isnot(Dvbar))) then
            ! Define norm and bounds
            do iscales=1,lct_blk%nscales
               minim%guess((iscales-1)*4+1:(iscales-1)*4+3) = (/Dhbar,Dhbar,Dvbar/)*nam%lct_scale_ratio**(iscales-1)
               if (lct_blk%nscales==1) then
                  minim%binf((iscales-1)*4+1:(iscales-1)*4+3) = (/1.0/nam%lct_scale_ratio,1.0/nam%lct_scale_ratio, &
                                                              & 1.0/nam%lct_scale_ratio/)*minim%guess(1:3)
                  minim%bsup((iscales-1)*4+1:(iscales-1)*4+3) = (/nam%lct_scale_ratio,nam%lct_scale_ratio,nam%lct_scale_ratio/) &
                                                              & *minim%guess(1:3)
               else
                  minim%binf((iscales-1)*4+1:(iscales-1)*4+3) = (/1.0/sqrt(nam%lct_scale_ratio),1.0/sqrt(nam%lct_scale_ratio), &
                                                              & 1.0/sqrt(nam%lct_scale_ratio)/) &
                                                              & *minim%guess(1:3)*nam%lct_scale_ratio**(iscales-1)
                  minim%bsup((iscales-1)*4+1:(iscales-1)*4+3) = (/sqrt(nam%lct_scale_ratio),sqrt(nam%lct_scale_ratio), &
                                                              & sqrt(nam%lct_scale_ratio)/)*minim%guess(1:3) &
                                                              & *nam%lct_scale_ratio**(iscales-1)
               end if
               minim%guess((iscales-1)*4+4) = 0.0
               if (nam%lct_diag(iscales)) then
                  ! Diagonal tensor
                  minim%binf((iscales-1)*4+4) = -1.0e-12
                  minim%bsup((iscales-1)*4+4) = 1.0e-12
               else
                  ! Non-diagonal tensor
                  minim%binf((iscales-1)*4+4) = -1.0
                  minim%bsup((iscales-1)*4+4) = 1.0
               end if
            end do
            do iscales=1,lct_blk%nscales-1
               minim%guess(lct_blk%nscales*4+1) = 1.0/real(lct_blk%nscales,kind_real)
               minim%binf(lct_blk%nscales*4+1) = 0.1
               minim%bsup(lct_blk%nscales*4+1) = 1.0
            end do

            ! Fill minim
            minim%obs = pack(lct_blk%raw(:,:,ic1a,il0),.true.)
            minim%cost_function = 'fit_lct'
            minim%algo = trim(nam%minim_algo)
            minim%nc3 = nam%nc3
            minim%nl0 = bpar%nl0r(ib)
            minim%nscales = lct_blk%nscales

            ! Compute fit
            call minim%compute(mpl,rng)

            ! Copy parameters
            do iscales=1,lct_blk%nscales
               do icomp=1,4
                   lct_blk%D(icomp,iscales,ic1a,il0) = minim%x((iscales-1)*4+icomp)
               end do
            end do
            if (lct_blk%nscales>1) then
               lct_blk%coef(1:lct_blk%nscales-1,ic1a,il0) = minim%x(lct_blk%nscales*4+1:lct_blk%nscales*4+lct_blk%nscales-1)
               lct_blk%coef(lct_blk%nscales,ic1a,il0) = 1.0-sum(lct_blk%coef(1:lct_blk%nscales-1,ic1a,il0))
            else
               lct_blk%coef(1,ic1a,il0) = 1.0
            end if

            ! Set vertical value at zero for the 2D case
            if (bpar%nl0r(ib)==1) lct_blk%D(3,:,ic1a,il0) = 0.0

            ! Check tensor validity
            valid = .true.
            do iscales=1,lct_blk%nscales
               if (valid) then
                  call check_cond(lct_blk%D(1,iscales,ic1a,il0),lct_blk%D(2,iscales,ic1a,il0),lct_blk%D(4,iscales,ic1a,il0),valid)
                  if (bpar%nl0r(ib)>1) valid = valid.and.(lct_blk%D(3,iscales,ic1a,il0)>0.0)
                  valid = valid.and.(lct_blk%coef(iscales,ic1a,il0)>0.0)
                  if (lct_blk%nscales>1) valid = valid.and.(inf(lct_blk%coef(iscales,ic1a,il0),1.0_kind_real))
               end if
            end do
            if (valid) then
               do iscales=1,lct_blk%nscales
                  if (nam%lct_diag(iscales)) lct_blk%D(4,iscales,ic1a,il0) = 0.0
               end do

               ! Rebuild fit
               call fit_lct(mpl,nam%nc3,bpar%nl0r(ib),minim%dxsq,minim%dysq,minim%dxdy,minim%dzsq,minim%dmask,lct_blk%nscales, &
             & lct_blk%D(:,:,ic1a,il0),lct_blk%coef(:,ic1a,il0),lct_blk%fit(:,:,ic1a,il0))

               ! Quality control
               lct_blk%qc_c1a(ic1a,il0) = 0.0
               norm = 0.0
               do jl0r=1,bpar%nl0r(ib)
                  jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
                  do jc3=1,nam%nc3
                     if (samp%c1ac3l0_log(ic1a,jc3,jl0)) then
                        if (mpl%msv%isnot(lct_blk%fit(jc3,jl0r,ic1a,il0)).and.lct_blk%raw(jc3,jl0r,ic1a,il0)>nam%lct_qc_th) then
                           lct_blk%qc_c1a(ic1a,il0) = lct_blk%qc_c1a(ic1a,il0)+ &
                                                    & (lct_blk%fit(jc3,jl0r,ic1a,il0)-lct_blk%raw(jc3,jl0r,ic1a,il0))**2
                           norm = norm+1.0
                        end if
                     end if
                  end do
               end do
               lct_blk%qc_c1a(ic1a,il0) = sqrt(lct_blk%qc_c1a(ic1a,il0)/norm)
               if (lct_blk%qc_c1a(ic1a,il0)>nam%lct_qc_max) then
                  ! Missing values
                  lct_blk%D(:,:,ic1a,il0) = mpl%msv%valr
                  lct_blk%coef(:,ic1a,il0) = mpl%msv%valr
                  lct_blk%fit(:,:,ic1a,il0) = mpl%msv%valr
               end if
            else
               ! Missing values
               lct_blk%D(:,:,ic1a,il0) = mpl%msv%valr
               lct_blk%coef(:,ic1a,il0) = mpl%msv%valr
               lct_blk%fit(:,:,ic1a,il0) = mpl%msv%valr
               lct_blk%qc_c1a(ic1a,il0) = mpl%msv%valr
            end if
         else
            ! Missing values
            lct_blk%D(:,:,ic1a,il0) = mpl%msv%valr
            lct_blk%coef(:,ic1a,il0) = mpl%msv%valr
            lct_blk%fit(:,:,ic1a,il0) = mpl%msv%valr
            lct_blk%qc_c1a(ic1a,il0) = mpl%msv%valr
         end if
      else
         ! Missing values
         lct_blk%D(:,:,ic1a,il0) = mpl%msv%valr
         lct_blk%coef(:,ic1a,il0) = mpl%msv%valr
         lct_blk%fit(:,:,ic1a,il0) = mpl%msv%valr
         lct_blk%qc_c1a(ic1a,il0) = mpl%msv%valr
      end if

      ! Release memory
      deallocate(norm_raw)
      deallocate(Dh)
      deallocate(Dv)
      deallocate(Dv_valid)
      deallocate(minim%x)
      deallocate(minim%guess)
      deallocate(minim%binf)
      deallocate(minim%bsup)
      deallocate(minim%obs)
      deallocate(minim%dxsq)
      deallocate(minim%dysq)
      deallocate(minim%dxdy)
      deallocate(minim%dzsq)
      deallocate(minim%dmask)

      ! Update
      call mpl%prog_print(ic1a)
   end do
   call mpl%prog_final
end do

! End associate
end associate

end subroutine lct_blk_compute

!----------------------------------------------------------------------
! Subroutine: lct_blk_filter
! Purpose: filter LCT
!----------------------------------------------------------------------
subroutine lct_blk_filter(lct_blk,mpl,nam,geom,bpar,samp)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! Averaged statistics block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(samp_type),intent(in) :: samp           ! Sampling

! Local variables
integer :: il0,jl0,jl0r,ic1a,ic1,ic0,jc3,jc0,icomp,iscales
real(kind_real) :: fld_c1a(samp%nc1a),dx,dy,dz
real(kind_real),allocatable :: fld_filt_c1a(:),dxsq(:,:),dysq(:,:),dxdy(:,:),dzsq(:,:)
logical :: valid,mask_c1a(samp%nc1a,geom%nl0)
logical,allocatable :: dmask(:,:)

! Associate
associate(ib=>lct_blk%ib)

! Define mask
do il0=1,geom%nl0
   do ic1a=1,samp%nc1a
      ic1 = samp%c1a_to_c1(ic1a)
      mask_c1a(ic1a,il0) = samp%c1l0_log(ic1,il0)
   end do
end do

! Allocation
if (nam%diag_rhflt>0.0) allocate(fld_filt_c1a(samp%nc1a))

do il0=1,geom%nl0
   do iscales=1,lct_blk%nscales
      do icomp=1,4+1
         ! Copy
         if (icomp<=4) then
            fld_c1a = lct_blk%D(icomp,iscales,:,il0)
         else
            fld_c1a = lct_blk%coef(iscales,:,il0)
         end if

         if (nam%diag_rhflt>0.0) then
            ! Copy
            fld_filt_c1a = fld_c1a

            ! Filter
            call samp%diag_filter(mpl,nam,'median',nam%diag_rhflt,fld_filt_c1a)
            call samp%diag_filter(mpl,nam,'gc99',nam%diag_rhflt,fld_filt_c1a)
         end if

         ! Fill missing values
         call samp%diag_fill(mpl,nam,fld_c1a)
         if (nam%diag_rhflt>0.0) call samp%diag_fill(mpl,nam,fld_filt_c1a)

         ! Copy
         if (icomp<=4) then
            lct_blk%D(icomp,iscales,:,il0) = fld_c1a
            if (nam%diag_rhflt>0.0) lct_blk%D_filt(icomp,iscales,:,il0) = fld_filt_c1a
         else
            lct_blk%coef(iscales,:,il0) = fld_c1a
            if (nam%diag_rhflt>0.0) lct_blk%coef_filt(iscales,:,il0) = fld_filt_c1a
         end if
      end do
   end do

   ! Fill missing values for quality control
   call samp%diag_fill(mpl,nam,lct_blk%qc_c1a(:,il0))
end do

if (nam%diag_rhflt>0.0) then
   ! Allocation
   allocate(dxsq(nam%nc3,bpar%nl0r(ib)))
   allocate(dysq(nam%nc3,bpar%nl0r(ib)))
   allocate(dxdy(nam%nc3,bpar%nl0r(ib)))
   allocate(dzsq(nam%nc3,bpar%nl0r(ib)))
   allocate(dmask(nam%nc3,bpar%nl0r(ib)))

   do il0=1,geom%nl0
      do ic1a=1,samp%nc1a
         ! Global index
         ic1 = samp%c1a_to_c1(ic1a)

         if (samp%c1l0_log(ic1,il0)) then
            ! Check tensor validity
            valid = .true.
            do iscales=1,lct_blk%nscales
               if (valid) then
                  call check_cond(lct_blk%D_filt(1,iscales,ic1a,il0),lct_blk%D_filt(2,iscales,ic1a,il0), &
                & lct_blk%D_filt(4,iscales,ic1a,il0),valid)
                  if (bpar%nl0r(ib)>1) valid = valid.and.(lct_blk%D_filt(3,iscales,ic1a,il0)>0.0)
                  valid = valid.and.(lct_blk%coef_filt(iscales,ic1a,il0)>0.0)
                  if (lct_blk%nscales>1) valid = valid.and.(lct_blk%coef_filt(iscales,ic1a,il0)<1.0)
               end if
            end do
            if (valid) then
               if (nam%lct_write_cor) then
                  ! Compute deltas
                  ic0 = samp%c1_to_c0(ic1)
                  do jl0r=1,bpar%nl0r(ib)
                     jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
                     do jc3=1,nam%nc3
                        dmask(jc3,jl0r) = samp%c1l0_log(ic1,il0).and.samp%c1ac3l0_log(ic1a,jc3,jl0)
                        if (dmask(jc3,jl0r)) then
                           jc0 = samp%c1ac3_to_c0(ic1a,jc3)
                           call geom%compute_deltas(ic0,il0,jc0,jl0,dx,dy,dz)
                           dxsq(jc3,jl0r) = dx**2
                           dysq(jc3,jl0r) = dy**2
                           dxdy(jc3,jl0r) = dx*dy
                           dzsq(jc3,jl0r) = dz**2
                        end if
                     end do
                  end do

                  ! Rebuild fit for full correlation output
                  call fit_lct(mpl,nam%nc3,bpar%nl0r(ib),dxsq,dysq,dxdy,dzsq,dmask,lct_blk%nscales, &
                & lct_blk%D_filt(:,:,ic1a,il0),lct_blk%coef_filt(:,ic1a,il0),lct_blk%fit_filt(:,:,ic1a,il0))
               end if
            else
               ! Missing values
               lct_blk%D_filt(:,:,ic1a,il0) = mpl%msv%valr
               lct_blk%coef_filt(:,ic1a,il0) = mpl%msv%valr
               if (nam%lct_write_cor) lct_blk%fit_filt(:,:,ic1a,il0) = mpl%msv%valr
            end if
         else
            ! Missing values
            lct_blk%D_filt(:,:,ic1a,il0) = mpl%msv%valr
            lct_blk%coef_filt(:,ic1a,il0) = mpl%msv%valr
            if (nam%lct_write_cor) lct_blk%fit_filt(:,:,ic1a,il0) = mpl%msv%valr
         end if
      end do
   end do

   ! Release memory
   deallocate(dxsq)
   deallocate(dysq)
   deallocate(dxdy)
   deallocate(dzsq)
   deallocate(dmask)
end if

! Release memory
if (nam%diag_rhflt>0.0) deallocate(fld_filt_c1a)

! End associate
end associate

end subroutine lct_blk_filter

!----------------------------------------------------------------------
! Subroutine: lct_blk_interp
! Purpose: interpolate LCT
!----------------------------------------------------------------------
subroutine lct_blk_interp(lct_blk,mpl,nam,geom,bpar,samp)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! Averaged statistics block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(samp_type),intent(in) :: samp           ! Sampling

! Local variables
integer :: il0,il0i,ic1a,ic1,icomp,ic0a,iscales
real(kind_real) :: det,Lavg_tot,norm_tot
real(kind_real) :: fld_c1a(samp%nc1a,geom%nl0,2*4+2),fld_c1b(samp%nc2b,geom%nl0),fld(geom%nc0a,geom%nl0,2*4+3)
real(kind_real),allocatable :: D(:,:,:,:),coef(:,:,:)
logical :: mask_c1a(samp%nc1a,geom%nl0)
character(len=1024),parameter :: subr = 'lct_interp'

! Associate
associate(ib=>lct_blk%ib)

! Define mask
do il0=1,geom%nl0
   do ic1a=1,samp%nc1a
      ic1 = samp%c1a_to_c1(ic1a)
      mask_c1a(ic1a,il0) = samp%c1l0_log(ic1,il0)
   end do
end do

! Allocation
allocate(D(4,lct_blk%nscales,samp%nc1a,geom%nl0))
allocate(coef(lct_blk%nscales,samp%nc1a,geom%nl0))

! Initialization
if (nam%diag_rhflt>0.0) then
   D = lct_blk%D_filt
   coef = lct_blk%coef_filt
else
   D = lct_blk%D
   coef = lct_blk%coef
end if

do iscales=1,lct_blk%nscales
   write(mpl%info,'(a10,a,i2)') '','Scale: ',iscales
   call mpl%flush

   ! Initialization
   fld_c1a = mpl%msv%valr
   fld = mpl%msv%valr

   ! Copy and inverse diffusion tensor
   write(mpl%info,'(a13,a)') '','Copy and inverse diffusion tensor'
   call mpl%flush
   do il0=1,geom%nl0
      do ic1a=1,samp%nc1a
         ic1 = samp%c1a_to_c1(ic1a)
         if (mask_c1a(ic1a,il0)) then
            ! Ensure positive-definiteness of D
            D(1,iscales,ic1a,il0) = max(Dmin,D(1,iscales,ic1a,il0))
            D(2,iscales,ic1a,il0) = max(Dmin,D(2,iscales,ic1a,il0))
            if (bpar%nl0r(ib)>1) D(3,iscales,ic1a,il0) = max(Dmin,D(3,iscales,ic1a,il0))
            D(4,iscales,ic1a,il0) = max(-1.0_kind_real+Dmin,min(D(4,iscales,ic1a,il0),1.0_kind_real-Dmin))

            ! Copy diffusion tensor
            fld_c1a(ic1a,il0,1) = D(1,iscales,ic1a,il0)
            fld_c1a(ic1a,il0,2) = D(2,iscales,ic1a,il0)
            fld_c1a(ic1a,il0,3) = D(3,iscales,ic1a,il0)
            fld_c1a(ic1a,il0,4) = sqrt(D(1,iscales,ic1a,il0)*D(2,iscales,ic1a,il0))*D(4,iscales,ic1a,il0)

            ! Inverse diffusion tensor
            call lct_d2h(mpl,fld_c1a(ic1a,il0,1),fld_c1a(ic1a,il0,2),fld_c1a(ic1a,il0,3),fld_c1a(ic1a,il0,4), &
          & fld_c1a(ic1a,il0,4+1),fld_c1a(ic1a,il0,4+2),fld_c1a(ic1a,il0,4+3),fld_c1a(ic1a,il0,4+4))

            ! Copy coefficient
            fld_c1a(ic1a,il0,2*4+1) = coef(iscales,ic1a,il0)

            ! Copy quality control
            fld_c1a(ic1a,il0,2*4+2) = lct_blk%qc_c1a(ic1a,il0)
         end if
      end do
   end do

   ! Interpolate components
   write(mpl%info,'(a13,a)') '','Interpolate components'
   call mpl%flush
   do icomp=1,2*4+2
      call samp%com_AB%ext(mpl,geom%nl0,fld_c1a(:,:,icomp),fld_c1b)
      do il0=1,geom%nl0
         il0i = min(il0,geom%nl0i)
         call samp%h(il0i)%apply(mpl,fld_c1b(:,il0),fld(:,il0,icomp))
      end do
   end do

   ! Compute horizontal length-scale and equivalent support radius
   write(mpl%info,'(a13,a)') '','Compute horizontal length-scale and equivalent support radius:'
   call mpl%flush
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         if (geom%mask_c0a(ic0a,il0)) then
            ! Length-scale = D determinant^{1/4}
            det = fld(ic0a,il0,1)*fld(ic0a,il0,2)-fld(ic0a,il0,4)**2
            if (det>0.0) then
               fld(ic0a,il0,2*4+3) = sqrt(sqrt(det))
            else
               call mpl%abort(subr,'non-valid horizontal diffusion tensor determinant, grid c0')
            end if
         end if
      end do
      call mpl%f_comm%allreduce(sum(fld(:,il0,2*4+3),mpl%msv%isnot(fld(:,il0,2*4+3))),Lavg_tot,fckit_mpi_sum())
      call mpl%f_comm%allreduce(real(count(mpl%msv%isnot(fld(:,il0,2*4+3))),kind_real),norm_tot,fckit_mpi_sum())
      if (norm_tot>0.0) then
         write(mpl%info,'(a16,a,i3,a,f10.2,a,f10.2,a)') '','Level',nam%levs(il0),' ~> ', &
       & Lavg_tot/norm_tot*reqkm,' km / ',Lavg_tot/norm_tot*gau2gc*reqkm,' km'
         call mpl%flush
      end if
   end do

   ! Copy output values
   lct_blk%D11(:,:,iscales) = fld(:,:,1)
   lct_blk%D22(:,:,iscales) = fld(:,:,2)
   lct_blk%D33(:,:,iscales) = fld(:,:,3)
   lct_blk%D12(:,:,iscales) = fld(:,:,4)
   lct_blk%H11(:,:,iscales) = fld(:,:,4+1)
   lct_blk%H22(:,:,iscales) = fld(:,:,4+2)
   lct_blk%H33(:,:,iscales) = fld(:,:,4+3)
   lct_blk%H12(:,:,iscales) = fld(:,:,4+4)
   lct_blk%Dcoef(:,:,iscales) = fld(:,:,2*4+1)
   lct_blk%qc_c0a = fld(:,:,2*4+2)
   lct_blk%DLh(:,:,iscales) = fld(:,:,2*4+3)
end do

! Allocation
deallocate(D)
deallocate(coef)

! End associate
end associate

end subroutine lct_blk_interp

!----------------------------------------------------------------------
! Subroutine: lct_blk_write
! Purpose: write LCT
!----------------------------------------------------------------------
subroutine lct_blk_write(lct_blk,mpl,nam,geom,bpar,io,filename)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! Averaged statistics block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(io_type),intent(in) :: io               ! I/O
character(len=*),intent(in) :: filename      ! Filename

! Local variables
integer :: iv,iscales
character(len=1) :: iscaleschar

! Associate
associate(ib=>lct_blk%ib)

iv = bpar%b_to_v2(ib)
do iscales=1,lct_blk%nscales
   ! Write fields
   write(iscaleschar,'(i1)') iscales
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D11_'//iscaleschar,lct_blk%D11(:,:,iscales))
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D22_'//iscaleschar,lct_blk%D22(:,:,iscales))
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D33_'//iscaleschar,lct_blk%D33(:,:,iscales))
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_D12_'//iscaleschar,lct_blk%D12(:,:,iscales))
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_H11_'//iscaleschar,lct_blk%H11(:,:,iscales))
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_H22_'//iscaleschar,lct_blk%H22(:,:,iscales))
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_H33_'//iscaleschar,lct_blk%H33(:,:,iscales))
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_H12_'//iscaleschar,lct_blk%H12(:,:,iscales))
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_coef_'//iscaleschar,lct_blk%Dcoef(:,:,iscales))
   call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_Lh_'//iscaleschar,lct_blk%DLh(:,:,iscales))
end do
call io%fld_write(mpl,nam,geom,filename,trim(nam%varname(iv))//'_qc',lct_blk%qc_c0a)

! End associate
end associate

end subroutine lct_blk_write

!----------------------------------------------------------------------
! Subroutine: lct_blk_write_cor
! Purpose: write full correlations
!----------------------------------------------------------------------
subroutine lct_blk_write_cor(lct_blk,mpl,nam,geom,bpar,samp,filename)

implicit none

! Passed variables
class(lct_blk_type),intent(inout) :: lct_blk ! Averaged statistics block
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(nam_type),intent(in) :: nam             ! Namelist
type(geom_type),intent(in) :: geom           ! Geometry
type(bpar_type),intent(in) :: bpar           ! Block parameters
type(samp_type),intent(in) :: samp           ! Sampling
character(len=*),intent(in) :: filename      ! Filename

! Local variables
integer :: info,ncid,nc3_id,nl0r_id,nc1a_id,nl0_id,lon_id,lat_id,raw_id,fit_id,fit_filt_id
integer :: ic1a,ic3,ic0
real(kind_real) :: lon(nam%nc3,samp%nc1a),lat(nam%nc3,samp%nc1a)
character(len=1024),parameter :: subr = 'lct_blk_write'

! Associate
associate(ib=>lct_blk%ib)

! Lon/lat initialization
do ic1a=1,samp%nc1a
   do ic3=1,nam%nc3
      ic0 = samp%c1ac3_to_c0(ic1a,ic3)
      lon(ic3,ic1a) = geom%lon(ic0)*rad2deg
      lat(ic3,ic1a) = geom%lat(ic0)*rad2deg
   end do
end do

! Check if the file exists
info = nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_noclobber,nf90_64bit_offset),ncid)
if (info==nf90_noerr) then
   ! Write namelist parameters
   call nam%write(mpl,ncid)
else
   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_write,ncid))

   ! Redef mode
   call mpl%ncerr(subr,nf90_redef(ncid))
end if

! Define dimensions and coordinates if necessary
nc3_id = mpl%ncdimcheck(subr,ncid,'nc3',nam%nc3,.true.)
nl0r_id = mpl%ncdimcheck(subr,ncid,'nl0r',bpar%nl0rmax,.true.,.true.)
nc1a_id = mpl%ncdimcheck(subr,ncid,'nc1a',samp%nc1a,.true.)
nl0_id = mpl%ncdimcheck(subr,ncid,'nl0',geom%nl0,.true.,.true.)

! Define variables if necessary
info = nf90_inq_varid(ncid,'lon',lon_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon',nc_kind_real,(/nc3_id,nc1a_id/),lon_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,'lat',lat_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat',nc_kind_real,(/nc3_id,nc1a_id/),lat_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,'raw',raw_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'raw',nc_kind_real,(/nc3_id,nl0r_id,nc1a_id,nl0_id/),raw_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,raw_id,'_FillValue',mpl%msv%valr))
end if
info = nf90_inq_varid(ncid,'fit',fit_id)
if (info/=nf90_noerr) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'fit',nc_kind_real,(/nc3_id,nl0r_id,nc1a_id,nl0_id/),fit_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,fit_id,'_FillValue',mpl%msv%valr))
end if
if (nam%diag_rhflt>0.0) then
   info = nf90_inq_varid(ncid,'fit_filt',fit_filt_id)
   if (info/=nf90_noerr) then
      call mpl%ncerr(subr,nf90_def_var(ncid,'fit_filt',nc_kind_real,(/nc3_id,nl0r_id,nc1a_id,nl0_id/),fit_filt_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,fit_filt_id,'_FillValue',mpl%msv%valr))
   end if
end if

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Write variables
call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,lat))
call mpl%ncerr(subr,nf90_put_var(ncid,raw_id,lct_blk%raw(1:bpar%nc3(ib),1:bpar%nl0r(ib),:,:),(/1,1,1,1/), &
 & (/bpar%nc3(ib),bpar%nl0r(ib),samp%nc1a,geom%nl0/)))
call mpl%ncerr(subr,nf90_put_var(ncid,fit_id,lct_blk%fit(1:bpar%nc3(ib),1:bpar%nl0r(ib),:,:),(/1,1,1,1/), &
 & (/bpar%nc3(ib),bpar%nl0r(ib),samp%nc1a,geom%nl0/)))
if (nam%diag_rhflt>0.0) call mpl%ncerr(subr,nf90_put_var(ncid,fit_filt_id,lct_blk%fit_filt(1:bpar%nc3(ib),1:bpar%nl0r(ib),:,:), &
 & (/1,1,1,1/),(/bpar%nc3(ib),bpar%nl0r(ib),samp%nc1a,geom%nl0/)))

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine lct_blk_write_cor

end module type_lct_blk
