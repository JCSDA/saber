!----------------------------------------------------------------------
! Module: type_vbal_blk
!> Vertical balance block derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_vbal_blk

use tools_kinds, only: kind_real
use tools_func, only: syminv,pseudoinv
use type_bpar, only: bpar_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type

implicit none

! Vertical balance block derived type
type vbal_blk_type
   integer :: iv                                  !< First variable index
   integer :: jv                                  !< Second variable index
   character(len=1024) :: name                    !< Name
   real(kind_real),allocatable :: auto(:,:,:)     !< Auto-covariance
   real(kind_real),allocatable :: cross(:,:,:)    !< Cross-covariance
   real(kind_real),allocatable :: auto_inv(:,:,:) !< Inverse auto-covariance
   real(kind_real),allocatable :: reg(:,:,:)      !< Regression
contains
   procedure :: alloc => vbal_blk_alloc
   procedure :: partial_dealloc => vbal_blk_partial_dealloc
   procedure :: dealloc => vbal_blk_dealloc
   procedure :: compute_covariances => vbal_blk_compute_covariances
   procedure :: compute_regression => vbal_blk_compute_regression
   procedure :: apply => vbal_blk_apply
   procedure :: apply_ad => vbal_blk_apply_ad
end type vbal_blk_type

private
public :: vbal_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: vbal_blk_alloc
!> Allocation
!----------------------------------------------------------------------
subroutine vbal_blk_alloc(vbal_blk,nam,geom,nc2b,iv,jv)

implicit none

! Passed variables
class(vbal_blk_type),intent(inout) :: vbal_blk !< Vertical balance block
type(nam_type),intent(in) :: nam               !< Namelist
type(geom_type),intent(in) :: geom             !< Geometry
integer,intent(in) :: nc2b                     !< Subset Sc2 size, halo B
integer,intent(in) :: iv                       !< First variable index
integer,intent(in) :: jv                       !< Second variable index

! Set attributes
vbal_blk%iv = iv
vbal_blk%jv = jv
vbal_blk%name = trim(nam%variables(jv))//'-'//trim(nam%variables(iv))

! Allocation
allocate(vbal_blk%auto(geom%nl0,geom%nl0,nc2b))
allocate(vbal_blk%cross(geom%nl0,geom%nl0,nc2b))
allocate(vbal_blk%auto_inv(geom%nl0,geom%nl0,nc2b))
allocate(vbal_blk%reg(geom%nl0,geom%nl0,nc2b))

end subroutine vbal_blk_alloc

!----------------------------------------------------------------------
! Subroutine: vbal_blk_partial_dealloc
!> Release memory (partial)
!----------------------------------------------------------------------
subroutine vbal_blk_partial_dealloc(vbal_blk)

implicit none

! Passed variables
class(vbal_blk_type),intent(inout) :: vbal_blk !< Vertical balance block

! Release memory
if (allocated(vbal_blk%auto)) deallocate(vbal_blk%auto)
if (allocated(vbal_blk%cross)) deallocate(vbal_blk%cross)
if (allocated(vbal_blk%auto_inv)) deallocate(vbal_blk%auto_inv)

end subroutine vbal_blk_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: vbal_blk_dealloc
!> Release memory (full)
!----------------------------------------------------------------------
subroutine vbal_blk_dealloc(vbal_blk)

implicit none

! Passed variables
class(vbal_blk_type),intent(inout) :: vbal_blk !< Vertical balance block

! Release memory
call vbal_blk%partial_dealloc
if (allocated(vbal_blk%reg)) deallocate(vbal_blk%reg)

end subroutine vbal_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: vbal_blk_compute_covariances
!> Compute auto- and cross-covariances
!----------------------------------------------------------------------
subroutine vbal_blk_compute_covariances(vbal_blk,mpl,geom,samp,ens,auto,cross)

implicit none

! Passed variables
class(vbal_blk_type),intent(in) :: vbal_blk                                !< Vertical balance block
type(mpl_type),intent(inout) :: mpl                                        !< MPI data
type(geom_type),intent(in) :: geom                                         !< Geometry
type(samp_type), intent(in) :: samp                                        !< Sampling
type(ens_type), intent(in) :: ens                                          !< Ensemble
real(kind_real),intent(out) :: auto(samp%nc1e,geom%nl0,geom%nl0,ens%nsub)  !< Auto-covariance
real(kind_real),intent(out) :: cross(samp%nc1e,geom%nl0,geom%nl0,ens%nsub) !< Cross-covariance

! Local variables
integer :: isub,ie_sub,ie,il0,ic1a,ic0,ic0a,jl0,ic1e,ic1u
real(kind_real) :: fld_c0a_1(geom%nc0a,geom%nl0),fld_c0a_2(geom%nc0a,geom%nl0)
real(kind_real) :: fld_1(samp%nc1a,geom%nl0),fld_2(samp%nc1a,geom%nl0)
real(kind_real) :: fld_ext_1(samp%nc1e,geom%nl0),fld_ext_2(samp%nc1e,geom%nl0)
 
! Initialization
auto = 0.0
cross = 0.0

! Loop on sub-ensembles
do isub=1,ens%nsub
   if (ens%nsub==1) then
      write(mpl%info,'(a10,a)') '','Full ensemble, member:'
      call mpl%flush(.false.)
   else
      write(mpl%info,'(a10,a,i6,a)') '','Sub-ensemble ',isub,', member:'
      call mpl%flush(.false.)
   end if

   ! Compute centered moments
   do ie_sub=1,ens%ne/ens%nsub
      write(mpl%info,'(i6)') ie_sub
      call mpl%flush(.false.)

      ! Full ensemble index
      ie = ie_sub+(isub-1)*ens%ne/ens%nsub

      ! Get perturbation on subset Sc0
      call ens%get_c0(mpl,vbal_blk%iv,geom,'pert',ie,fld_c0a_1)
      call ens%get_c0(mpl,vbal_blk%jv,geom,'pert',ie,fld_c0a_2)

      ! Copy all separations points
      fld_1 = 0.0
      fld_2 = 0.0
      !$omp parallel do schedule(static) private(il0,ic1a,ic0,ic0a)
      do il0=1,geom%nl0
         do ic1a=1,samp%nc1a
            if (samp%smask_c1a(ic1a,il0)) then
               ! Index
               ic0a = samp%c1a_to_c0a(ic1a)

               ! Copy points
               fld_1(ic1a,il0) = fld_c0a_1(ic0a,il0)
               fld_2(ic1a,il0) = fld_c0a_2(ic0a,il0)
            end if
         end do
      end do
      !$omp end parallel do

      ! Halo extension
      call samp%com_AE%ext(mpl,fld_1,fld_ext_1)
      call samp%com_AE%ext(mpl,fld_2,fld_ext_2)

      !$omp parallel do schedule(static) private(il0,jl0,ic1e,ic1u)
      do il0=1,geom%nl0
         do jl0=1,geom%nl0
            do ic1e=1,samp%nc1e
               ! Index
               ic1u = samp%c1e_to_c1u(ic1e)

               ! Auto and cross-covariances
               if (samp%smask_c1u(ic1u,il0).and.samp%smask_c1u(ic1u,jl0)) then
                  auto(ic1e,jl0,il0,isub) = auto(ic1e,jl0,il0,isub)+fld_ext_2(ic1e,il0)*fld_ext_2(ic1e,jl0)
                  cross(ic1e,jl0,il0,isub) = cross(ic1e,jl0,il0,isub)+fld_ext_2(ic1e,il0)*fld_ext_1(ic1e,jl0)
               end if
            end do
         end do
      end do
      !$omp end parallel do
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
end do

end subroutine vbal_blk_compute_covariances

!----------------------------------------------------------------------
! Subroutine: vbal_blk_compute_regression
!> Compute regression
!----------------------------------------------------------------------
subroutine vbal_blk_compute_regression(vbal_blk,mpl,nam,geom,samp,nsub,auto,cross,ic2b)

implicit none

! Passed variables
class(vbal_blk_type),intent(inout) :: vbal_blk                        !< Vertical balance block
type(mpl_type),intent(inout) :: mpl                                   !< MPI data
type(nam_type), intent(in) :: nam                                     !< Namelist
type(geom_type),intent(in) :: geom                                    !< Geometry
type(samp_type), intent(in) :: samp                                   !< Sampling
integer,intent(in) :: nsub                                            !< Number of sub-ensembles
real(kind_real),intent(in) :: auto(samp%nc1e,geom%nl0,geom%nl0,nsub)  !< Auto-covariance
real(kind_real),intent(in) :: cross(samp%nc1e,geom%nl0,geom%nl0,nsub) !< Cross-covariance
integer,intent(in) :: ic2b                                            !< Index

! Local variables
integer :: i,jl0_min,jl0_max,jl0,il0,ic1e,ic1u,nc1a,nc1max,ierr
real(kind_real),allocatable :: list_auto(:),list_cross(:)
logical :: valid

! Initialization
vbal_blk%auto(:,:,ic2b) = 0.0
vbal_blk%cross(:,:,ic2b) = 0.0

! Allocation
nc1max = count(samp%vbal_mask(:,ic2b))
allocate(list_auto(nc1max))
allocate(list_cross(nc1max))

do il0=1,geom%nl0
   ! Indices
   if (nam%vbal_diag_reg((vbal_blk%iv-1)*(vbal_blk%iv-2)/2+vbal_blk%jv)) then
      ! Diagonal regression
      jl0_min = il0
      jl0_max = il0
   else
      jl0_min = 1
      jl0_max = geom%nl0
   end if

   do jl0=jl0_min,jl0_max
      ! Initialization
      list_auto = mpl%msv%valr
      list_cross = mpl%msv%valr

      ! Fill lists
      i = 0
      do ic1e=1,samp%nc1e
         ! Index
         ic1u = samp%c1e_to_c1u(ic1e)

         if (samp%smask_c1u(ic1u,il0).and.samp%smask_c1u(ic1u,jl0).and.samp%vbal_mask(ic1u,ic2b)) then
            ! Update
            i = i+1

            ! Averages over sub-ensembles
            list_auto(i) = sum(auto(ic1e,jl0,il0,:))/real(nsub,kind_real)
            list_cross(i) = sum(cross(ic1e,jl0,il0,:))/real(nsub,kind_real)
         end if
      end do

      ! Average
      nc1a = count(mpl%msv%isnot(list_auto))
      if (nc1a>0) then
         vbal_blk%auto(jl0,il0,ic2b) = sum(list_auto,mask=mpl%msv%isnot(list_auto))/real(nc1a,kind_real)
         vbal_blk%cross(jl0,il0,ic2b) = sum(list_cross,mask=mpl%msv%isnot(list_cross))/real(nc1a,kind_real)
      else
         vbal_blk%auto(jl0,il0,ic2b) = mpl%msv%valr
         vbal_blk%cross(jl0,il0,ic2b) = mpl%msv%valr
      end if
   end do
end do

! Release memory
deallocate(list_auto)
deallocate(list_cross)

! Check if this column is valid (no missing values, positive diagonal for auto-covariance)
valid = mpl%msv%isallnot(vbal_blk%auto(:,:,ic2b)).and.mpl%msv%isallnot(vbal_blk%cross(:,:,ic2b))
do il0=1,geom%nl0
   if (.not.(vbal_blk%auto(il0,il0,ic2b)>0.0)) valid = .false.
end do

if (valid) then
   if (nam%vbal_diag_auto((vbal_blk%iv-1)*(vbal_blk%iv-2)/2+vbal_blk%jv) &
 & .or.nam%vbal_diag_reg((vbal_blk%iv-1)*(vbal_blk%iv-2)/2+vbal_blk%jv)) then
      ! Diagonal inversion
      vbal_blk%auto_inv(:,:,ic2b) = 0.0
      do il0=1,geom%nl0
         vbal_blk%auto_inv(il0,il0,ic2b) = 1.0/vbal_blk%auto(il0,il0,ic2b)
      end do
   else
      ! Inverse the vertical auto-covariance
      if (nam%vbal_pseudo_inv) then
         if (nam%vbal_pseudo_inv_mmax>0) then
            call pseudoinv(mpl,geom%nl0,vbal_blk%auto(:,:,ic2b),vbal_blk%auto_inv(:,:,ic2b),ierr,mmax=nam%vbal_pseudo_inv_mmax)
         elseif (nam%vbal_pseudo_inv_var_th>0.0) then
            call pseudoinv(mpl,geom%nl0,vbal_blk%auto(:,:,ic2b),vbal_blk%auto_inv(:,:,ic2b),ierr,var_th=nam%vbal_pseudo_inv_var_th)
         end if
      else
         call syminv(mpl,geom%nl0,vbal_blk%auto(:,:,ic2b),vbal_blk%auto_inv(:,:,ic2b),ierr)
      end if
      if (ierr/=0) then
         ! Diagonal inversion
         vbal_blk%auto_inv(:,:,ic2b) = 0.0
         do il0=1,geom%nl0
            vbal_blk%auto_inv(il0,il0,ic2b) = 1.0/vbal_blk%auto(il0,il0,ic2b)
         end do
      end if
   end if

   ! Compute the regression
   vbal_blk%reg(:,:,ic2b) = matmul(vbal_blk%cross(:,:,ic2b),vbal_blk%auto_inv(:,:,ic2b))
else
   ! No regression
   vbal_blk%reg(:,:,ic2b) = 0.0
end if

end subroutine vbal_blk_compute_regression

!----------------------------------------------------------------------
! Subroutine: vbal_blk_apply
!> Apply vertical balance block
!----------------------------------------------------------------------
subroutine vbal_blk_apply(vbal_blk,geom,h_n_s,h_c2b,h_S,fld)

implicit none

! Passed variables
class(vbal_blk_type),intent(in) :: vbal_blk              !< Vertical balance block
type(geom_type),intent(in) :: geom                       !< Geometry
integer,intent(in) :: h_n_s(geom%nc0a,geom%nl0i)         !< Number of neighbors for the horizontal interpolation
integer,intent(in) :: h_c2b(3,geom%nc0a,geom%nl0i)       !< Index of neighbors for the horizontal interpolation
real(kind_real),intent(in) :: h_S(3,geom%nc0a,geom%nl0i) !< Weight of neighbors for the horizontal interpolation
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) !< Source/destination vector

! Local variables
integer :: ic0a,il0,jl0,i_s,ic2b
real(kind_real) :: S,fld_tmp(geom%nc0a,geom%nl0)

! Initialization
fld_tmp = 0.0

do ic0a=1,geom%nc0a
   do il0=1,geom%nl0
      do jl0=1,geom%nl0
         ! Loop over neighbors
         do i_s=1,h_n_s(ic0a,min(il0,geom%nl0i))
            ! Find neighbor index and weight
            ic2b = h_c2b(i_s,ic0a,min(il0,geom%nl0i))
            S = h_S(i_s,ic0a,min(il0,geom%nl0i))

            ! Apply regression coefficient weighted by the neighbor weight
            fld_tmp(ic0a,il0) = fld_tmp(ic0a,il0)+S*vbal_blk%reg(il0,jl0,ic2b)*fld(ic0a,jl0)
         end do
      end do
   end do
end do

! Final copy
fld = fld_tmp

end subroutine vbal_blk_apply

!----------------------------------------------------------------------
! Subroutine: vbal_blk_apply_ad
!> Apply adjoint vertical balance block
!----------------------------------------------------------------------
subroutine vbal_blk_apply_ad(vbal_blk,geom,h_n_s,h_c2b,h_S,fld)

implicit none

! Passed variables
class(vbal_blk_type),intent(in) :: vbal_blk              !< Vertical balance block
type(geom_type),intent(in) :: geom                       !< Geometry
integer,intent(in) :: h_n_s(geom%nc0a,geom%nl0i)         !< Number of neighbors for the horizontal interpolation
integer,intent(in) :: h_c2b(3,geom%nc0a,geom%nl0i)       !< Index of neighbors for the horizontal interpolation
real(kind_real),intent(in) :: h_S(3,geom%nc0a,geom%nl0i) !< Weight of neighbors for the horizontal interpolation
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) !< Source/destination vector

! Local variables
integer :: ic0a,il0,jl0,i_s,ic2b
real(kind_real) :: S,fld_tmp(geom%nc0a,geom%nl0)

! Initialization
fld_tmp = 0.0

do ic0a=1,geom%nc0a
   do il0=1,geom%nl0
      do jl0=1,geom%nl0
         ! Loop over neighbors
         do i_s=1,h_n_s(ic0a,min(il0,geom%nl0i))
            ! Find neighbor index and weight
            ic2b = h_c2b(i_s,ic0a,min(il0,geom%nl0i))
            S = h_S(i_s,ic0a,min(il0,geom%nl0i))

            ! Apply regression coefficient weighted by the neighbor weight
            fld_tmp(ic0a,jl0) = fld_tmp(ic0a,jl0)+S*vbal_blk%reg(il0,jl0,ic2b)*fld(ic0a,il0)
         end do
      end do
   end do
end do

! Final copy
fld = fld_tmp

end subroutine vbal_blk_apply_ad

end module type_vbal_blk
