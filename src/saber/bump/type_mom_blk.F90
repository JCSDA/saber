!----------------------------------------------------------------------
! Module: type_mom_blk
!> Moments block derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_mom_blk

use tools_kinds, only: kind_real
use type_bpar, only: bpar_type
use type_geom, only: geom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type

implicit none

! Moments block derived type
type mom_blk_type
   integer :: ib                                  !< Block index
   integer :: ne                                  !< Ensemble size
   integer :: nsub                                !< Number of sub-ensembles
   real(kind_real),allocatable :: m2_1(:,:,:)     !< Variance for variable 1
   real(kind_real),allocatable :: m2_2(:,:,:,:)   !< Variance for variable 2
   real(kind_real),allocatable :: m11(:,:,:,:,:)  !< Covariance
   real(kind_real),allocatable :: m22(:,:,:,:,:)  !< Fourth-order centered moment
contains
   procedure :: alloc => mom_blk_alloc
   procedure :: dealloc => mom_blk_dealloc
   procedure :: ext => mom_blk_ext
end type mom_blk_type

private
public :: mom_blk_type

contains

!----------------------------------------------------------------------
! Subroutine: mom_blk_alloc
!> Allocation
!----------------------------------------------------------------------
subroutine mom_blk_alloc(mom_blk,nc1,geom,bpar,ne,nsub)

implicit none

! Passed variables
class(mom_blk_type),intent(inout) :: mom_blk !< Moments block
integer,intent(in) :: nc1                    !< Subsampling size
type(geom_type),intent(in) :: geom           !< Geometry
type(bpar_type),intent(in) :: bpar           !< Block parameters
integer,intent(in) :: ne                     !< Ensemble size
integer,intent(in) :: nsub                   !< Number of sub-ensembles

! Associate
associate(ib=>mom_blk%ib)

! Attributes
mom_blk%ne = ne
mom_blk%nsub = nsub

! Allocation
allocate(mom_blk%m2_1(nc1,geom%nl0,nsub))
allocate(mom_blk%m2_2(nc1,bpar%nc3(ib),geom%nl0,nsub))
allocate(mom_blk%m11(nc1,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,nsub))
allocate(mom_blk%m22(nc1,bpar%nc3(ib),bpar%nl0r(ib),geom%nl0,nsub))

! End associate
end associate

end subroutine mom_blk_alloc

!----------------------------------------------------------------------
! Subroutine: mom_blk_dealloc
!> Release memory
!----------------------------------------------------------------------
subroutine mom_blk_dealloc(mom_blk)

implicit none

! Passed variables
class(mom_blk_type),intent(inout) :: mom_blk !< Moments block

! Release memory
if (allocated(mom_blk%m2_1)) deallocate(mom_blk%m2_1)
if (allocated(mom_blk%m2_2)) deallocate(mom_blk%m2_2)
if (allocated(mom_blk%m11)) deallocate(mom_blk%m11)
if (allocated(mom_blk%m22)) deallocate(mom_blk%m22)

end subroutine mom_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: mom_blk_ext
!> Halo extension
!----------------------------------------------------------------------
subroutine mom_blk_ext(mom_blk_out,mpl,geom,bpar,samp,mom_blk_in)

implicit none

! Passed variables
class(mom_blk_type),intent(inout) :: mom_blk_out !< Extended moments block
type(mpl_type),intent(inout) :: mpl              !< MPI data
type(geom_type),intent(in) :: geom               !< Geometry
type(bpar_type),intent(in) :: bpar               !< Block parameters
type(samp_type),intent(in) :: samp               !< Sampling
type(mom_blk_type),intent(in) :: mom_blk_in      !< Reduced moments block

! Local variables
integer :: npack,ipack,isub,il0,jc3,jl0r
real(kind_real),allocatable :: sbuf(:,:),rbuf(:,:)

! Associate
associate(ib=>mom_blk_in%ib)

! Allocation
mom_blk_out%ib = ib
call mom_blk_out%alloc(samp%nc1d,geom,bpar,mom_blk_in%ne,mom_blk_in%nsub)
npack = (1+bpar%nc3(ib)*(1+2*bpar%nl0r(ib)))*geom%nl0*mom_blk_out%nsub
allocate(sbuf(samp%nc1a,npack))
allocate(rbuf(samp%nc1d,npack))

! Pack
ipack = 0
do isub=1,mom_blk_out%nsub
   do il0=1,geom%nl0
      ipack = ipack+1
      sbuf(:,ipack) = mom_blk_in%m2_1(:,il0,isub)
      do jc3=1,bpar%nc3(ib)
         ipack = ipack+1
         sbuf(:,ipack) = mom_blk_in%m2_2(:,jc3,il0,isub)
         do jl0r=1,bpar%nl0r(ib)
            ipack = ipack+1
            sbuf(:,ipack) = mom_blk_in%m11(:,jc3,jl0r,il0,isub)
            ipack = ipack+1
            sbuf(:,ipack) = mom_blk_in%m22(:,jc3,jl0r,il0,isub)
         end do
      end do
   end do
end do

! Communication
call samp%com_AD%ext(mpl,sbuf,rbuf)

! Unpack
ipack = 0
do isub=1,mom_blk_out%nsub
   do il0=1,geom%nl0
      ipack = ipack+1
      mom_blk_out%m2_1(:,il0,isub) = rbuf(:,ipack)
      do jc3=1,bpar%nc3(ib)
         ipack = ipack+1
         mom_blk_out%m2_2(:,jc3,il0,isub) = rbuf(:,ipack)
         do jl0r=1,bpar%nl0r(ib)
            ipack = ipack+1
            mom_blk_out%m11(:,jc3,jl0r,il0,isub) = rbuf(:,ipack)
            ipack = ipack+1
            mom_blk_out%m22(:,jc3,jl0r,il0,isub) = rbuf(:,ipack)
         end do
      end do
   end do
end do

! Release memory
deallocate(sbuf)
deallocate(rbuf)

! End associate
end associate

end subroutine mom_blk_ext

end module type_mom_blk
