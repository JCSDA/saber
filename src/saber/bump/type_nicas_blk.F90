!----------------------------------------------------------------------
! Module: type_nicas_blk
! Purpose: NICAS data block derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_nicas_blk

use fckit_mpi_module, only: fckit_mpi_status
use netcdf
!$ use omp_lib
use tools_const, only: pi,req,reqkm,deg2rad,rad2deg
use tools_func, only: gc2gau,lonlatmod,sphere_dist,fit_func
use tools_kinds, only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use tools_repro, only: supeq,sup,inf
use tools_samp, only: initialize_sampling
use type_bpar, only: bpar_type
use type_cmat_blk, only: cmat_blk_type
use type_com, only: com_type,com_ntag
use type_geom, only: geom_type
use type_io, only: io_type
use type_tree, only: tree_type
use type_linop, only: linop_type,linop_ntag
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type
use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_min,fckit_mpi_max

implicit none

real(kind_real),parameter :: sqrt_r = 0.81_kind_real      ! Square-root factor (empirical)
real(kind_real),parameter :: sqrt_r_dble = 0.96_kind_real ! Square-root factor (empirical)
real(kind_real),parameter :: sqrt_rfac = 0.9_kind_real    ! Square-root factor (empirical)
real(kind_real),parameter :: sqrt_coef = 0.54_kind_real   ! Square-root factor (empirical)
real(kind_real),parameter :: S_inf = 1.0e-2_kind_real     ! Minimum value for the convolution coefficients

! Ball data derived type
type balldata_type
   integer :: nbd                        ! Number of values
   integer,allocatable :: bd_to_c1(:)    ! Ball data index to subset Sc1
   integer,allocatable :: bd_to_l1(:)    ! Ball data index to subset Sl1
   real(kind_real),allocatable :: val(:) ! Values
contains
   procedure :: alloc => balldata_alloc
   procedure :: dealloc => balldata_dealloc
   procedure :: pack => balldata_pack
end type balldata_type

! NICAS block derived type
type nicas_blk_type
   ! Block index and name
   integer :: ib                                   ! Block index
   character(len=1024) :: name                     ! Name
   character(len=1024) :: subsamp                  ! Subsampling structure
   logical :: double_fit                           ! Double fit
   logical :: anisotropic                          ! Anisotropic tensor

   ! Specific geometry
   integer :: nc1                                  ! Number of points in subset Sc1
   integer,allocatable :: vbot(:)                  ! Bottom level in grid Gh
   integer,allocatable :: vtop(:)                  ! Top level in grid Gh
   integer,allocatable :: nc2(:)                   ! Number of points in subset Sc2
   integer :: ns                                   ! Number of subgrid nodes

   ! Other data

   ! Support radius
   real(kind_real),allocatable :: rhs_avg(:)       ! Average sampling horizontal support radius at each level

   ! Level selection
   integer :: il0_first                            ! First valid level
   integer :: il0_last                             ! Last valid level
   logical,allocatable :: slev(:)                  ! Selected levels

   ! Parameters/normalization conversion
   integer,allocatable :: s_to_c1(:)               ! Subgrid to subset Sc1
   integer,allocatable :: s_to_l1(:)               ! Subgrid to subset Sl1
   integer,allocatable :: c1_to_c0(:)              ! Subset Sc1 to subset Sc0
   integer,allocatable :: l1_to_l0(:)              ! Subset Sl1 to subset Sl0
   integer,allocatable :: l0_to_l1(:)              ! Subset Sl0 to subset Sl1
   logical,allocatable :: mask_c1(:,:)             ! Mask from subset C1 to subgrid
   logical,allocatable :: mask_c2(:,:)             ! Mask from subset C2 to subgrid
   integer,allocatable :: c1l1_to_s(:,:)           ! Grid Gv to subgrid
   integer,allocatable :: c1_to_proc(:)            ! Subset Sc1 to local task
   integer,allocatable :: s_to_proc(:)             ! Subgrid to local task

   ! MPI distribution
   integer :: nc1a                                 ! Number of points in subset Sc1 on halo A
   integer :: nc1bb                                ! Number of points in subset Sc1 on halo B (extended)
   integer :: nsbb                                 ! Number of points in subgrid on halo B (extended)
   logical,allocatable :: lcheck_sa(:)             ! Detection of halo A on subgrid
   logical,allocatable :: lcheck_sb(:)             ! Detection of halo B on subgrid
   logical,allocatable :: lcheck_sc(:)             ! Detection of halo C on subgrid
   integer,allocatable :: c1a_to_c1(:)             ! Subset Sc1, halo A to global
   integer,allocatable :: c1_to_c1a(:)             ! Subset Sc1, global to halo A
   integer,allocatable :: c1b_to_c1(:)             ! Subset Sc1, halo B to global
   integer,allocatable :: c1_to_c1b(:)             ! Subset Sc1, global to halo B
   integer,allocatable :: c1bl1_to_sb(:,:)         ! Halo B, subset Sc1 to subgrid
   integer,allocatable :: c1a_to_c0a(:)            ! Halo A, subset Sc1 to subset Sc0
   integer,allocatable :: s_to_sa(:)               ! Subgrid, global to halo A
   integer,allocatable :: sb_to_s(:)               ! Subgrid, halo B to global
   integer,allocatable :: sbb_to_s(:)              ! Subgrid, halo B (extended) to global
   integer,allocatable :: c1_to_c1bb(:)            ! Subset Sc1, global to halo B (extended)
   integer,allocatable :: c1bb_to_c1(:)            ! Subset Sc1, halo B (extended) to global
   integer,allocatable :: sc_to_s(:)               ! Subgrid, halo C to global

   ! Convolution parameters
   real(kind_real) :: rhmax                        ! Maximum horizontal support radius
   real(kind_real),allocatable :: rh_c1(:,:)       ! Horizontal support radius on subset Sc1
   real(kind_real),allocatable :: rv_c1(:,:)       ! Vertical support radius on subset Sc1
   real(kind_real),allocatable :: H11_c1(:,:)      ! Local correlation tensor, component 11, on subset Sc1
   real(kind_real),allocatable :: H22_c1(:,:)      ! Local correlation tensor, component 22, on subset Sc1
   real(kind_real),allocatable :: H33_c1(:,:)      ! Local correlation tensor, component 33, on subset Sc1
   real(kind_real),allocatable :: H12_c1(:,:)      ! Local correlation tensor, component 12, on subset Sc1
   type(balldata_type),allocatable :: Hcoef(:)     ! Tensor coefficient on subset Sc1
   type(balldata_type),allocatable :: distnorm(:)  ! Normalized distance
   type(balldata_type),allocatable :: distnormv(:) ! Normalized distance, vertical part
   type(balldata_type),allocatable :: rfac(:)      ! Double-fit radius factor
   type(balldata_type),allocatable :: coef(:)      ! Double-fit coefficient

   ! Extended data for normalization computation
   integer :: nsc_nor                              ! Number of subgrid nodes on halo C (extended for normalization)
   integer,allocatable :: sc_nor_to_s(:)           ! Subgrid, halo C to global (extended for normalization)
   integer,allocatable :: s_to_sc_nor(:)           ! Subgrid, global to halo C (extended for normalization)
   integer,allocatable :: sb_to_sc_nor(:)          ! Subgrid, halo B to halo C (extended for normalization)
   type(linop_type) :: c_nor                       ! Convolution (extended for normalization)

   ! Tree
   type(tree_type) :: tree                         ! Tree

   ! Required data to apply NICAS

   ! Number of points
   integer :: nc0a                                 ! Number of points in subset Sc0 on halo A (required for I/O)
   integer :: nc1b                                 ! Number of points in subset Sc1 on halo B
   integer :: nl1                                  ! Number of levels in subset Sl1
   integer :: nsa                                  ! Number of subgrid nodes on halo A
   integer :: nsb                                  ! Number of subgrid nodes on halo B
   integer :: nsc                                  ! Number of subgrid nodes on halo C
   integer :: nc0d                                 ! Number of points in subset Sc1 on halo D
   integer :: nc0dinv                              ! Number of points in subset Sc1 on halo Dinv

   ! Valid levels
   logical,allocatable :: vlev(:)                  ! Valid levels

   ! Local to global
   integer,allocatable :: sa_to_s(:)               ! Subgrid, halo A to global

   ! Inter-halo conversions
   integer,allocatable :: sa_to_sc(:)              ! Subgrid, halo A to halo C
   integer,allocatable :: sb_to_sc(:)              ! Subgrid, halo B to halo C

   ! Linear operators
   type(linop_type) :: c                           ! Convolution
   type(linop_type),allocatable :: h(:)            ! Horizontal interpolation
   type(linop_type) :: v                           ! Vertical interpolation
   type(linop_type),allocatable :: s(:)            ! Subsample interpolation
   type(linop_type),allocatable :: d(:,:)          ! Advection
   type(linop_type),allocatable :: dinv(:,:)       ! Inverse advection

   ! Copy conversions
   integer,allocatable :: sb_to_c1b(:)             ! Subgrid to subset Sc1 on halo B
   integer,allocatable :: sb_to_l1(:)              ! Subgrid to subset Sl1 on halo B

   ! Normalization
   real(kind_real),allocatable :: norm(:,:)        ! Normalization factor

   ! Localization weights
   real(kind_real),allocatable :: coef_ens(:,:)    ! Ensemble coefficient square-root
   real(kind_real) :: wgt                          ! Main weight

   ! Communications
   type(com_type) :: com_AB                        ! Communication between halos A and B
   type(com_type) :: com_AC                        ! Communication between halos A and C
   type(com_type) :: com_AD                        ! Communication between halos A and D
   type(com_type) :: com_ADinv                     ! Communication between halos A and Dinv

   ! Required data to write grids
   real(kind_real),allocatable :: lon_sa(:)        ! Subgrid, halo A longitudes to subset Sc0, global
   real(kind_real),allocatable :: lat_sa(:)        ! Subgrid, halo A latitudes
   integer,allocatable :: lev_sa(:)                ! Subgrid, halo A levels
   real(kind_real),allocatable :: lon_sb(:)        ! Subgrid, halo A longitudes to subset Sc0, global
   real(kind_real),allocatable :: lat_sb(:)        ! Subgrid, halo A latitudes
   integer,allocatable :: lev_sb(:)                ! Subgrid, halo A levels
   real(kind_real),allocatable :: lon_sc(:)        ! Subgrid, halo A longitudes to subset Sc0, global
   real(kind_real),allocatable :: lat_sc(:)        ! Subgrid, halo A latitudes
   integer,allocatable :: lev_sc(:)                ! Subgrid, halo A levels
contains
   procedure :: partial_dealloc => nicas_blk_partial_dealloc
   procedure :: dealloc => nicas_blk_dealloc
   procedure :: read => nicas_blk_read
   procedure :: write => nicas_blk_write
   procedure :: write_grids => nicas_blk_write_grids
   procedure :: receive => nicas_blk_receive
   procedure :: send => nicas_blk_send
   procedure :: compute_parameters => nicas_blk_compute_parameters
   procedure :: compute_sampling_c1 => nicas_blk_compute_sampling_c1
   procedure :: compute_sampling_v => nicas_blk_compute_sampling_v
   procedure :: compute_sampling_c2 => nicas_blk_compute_sampling_c2
   procedure :: compute_mpi_a => nicas_blk_compute_mpi_a
   procedure :: compute_mpi_ab => nicas_blk_compute_mpi_ab
   procedure :: compute_interp_v => nicas_blk_compute_interp_v
   procedure :: compute_convol => nicas_blk_compute_convol
   procedure :: compute_convol_network => nicas_blk_compute_convol_network
   procedure :: compute_convol_distance => nicas_blk_compute_convol_distance
   procedure :: compute_convol_weights => nicas_blk_compute_convol_weights
   procedure :: compute_mpi_c => nicas_blk_compute_mpi_c
   procedure :: compute_normalization => nicas_blk_compute_normalization
   procedure :: compute_grids => nicas_blk_compute_grids
   procedure :: compute_adv => nicas_blk_compute_adv
   procedure :: apply => nicas_blk_apply
   procedure :: apply_from_sqrt => nicas_blk_apply_from_sqrt
   procedure :: apply_sqrt => nicas_blk_apply_sqrt
   procedure :: apply_sqrt_ad => nicas_blk_apply_sqrt_ad
   procedure :: apply_interp => nicas_blk_apply_interp
   procedure :: apply_interp_ad => nicas_blk_apply_interp_ad
   procedure :: apply_interp_h => nicas_blk_apply_interp_h
   procedure :: apply_interp_h_ad => nicas_blk_apply_interp_h_ad
   procedure :: apply_interp_v => nicas_blk_apply_interp_v
   procedure :: apply_interp_v_ad => nicas_blk_apply_interp_v_ad
   procedure :: apply_interp_s => nicas_blk_apply_interp_s
   procedure :: apply_interp_s_ad => nicas_blk_apply_interp_s_ad
   procedure :: apply_convol => nicas_blk_apply_convol
   procedure :: apply_adv => nicas_blk_apply_adv
   procedure :: apply_adv_ad => nicas_blk_apply_adv_ad
   procedure :: apply_adv_inv => nicas_blk_apply_adv_inv
   procedure :: test_adjoint => nicas_blk_test_adjoint
   procedure :: test_dirac => nicas_blk_test_dirac
end type nicas_blk_type

integer,parameter :: nicas_blk_ntag = 4 ! Number of communication steps to send/receive

private
public :: nicas_blk_type,nicas_blk_ntag

contains

!----------------------------------------------------------------------
! Subroutine: balldata_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine balldata_alloc(balldata)

implicit none

! Passed variables
class(balldata_type),intent(inout) :: balldata ! Ball data

! Allocation
allocate(balldata%bd_to_c1(balldata%nbd))
allocate(balldata%bd_to_l1(balldata%nbd))
allocate(balldata%val(balldata%nbd))

end subroutine balldata_alloc

!----------------------------------------------------------------------
! Subroutine: balldata_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine balldata_dealloc(balldata)

implicit none

! Passed variables
class(balldata_type),intent(inout) :: balldata ! Ball data

! Release memory
if (allocated(balldata%bd_to_c1)) deallocate(balldata%bd_to_c1)
if (allocated(balldata%bd_to_l1)) deallocate(balldata%bd_to_l1)
if (allocated(balldata%val)) deallocate(balldata%val)

end subroutine balldata_dealloc

!----------------------------------------------------------------------
! Subroutine: balldata_pack
! Purpose: pack data into balldata object
!----------------------------------------------------------------------
subroutine balldata_pack(balldata,mpl,nc1,nl1,val)

implicit none

! Passed variables
class(balldata_type),intent(inout) :: balldata ! Ball data
type(mpl_type),intent(inout) :: mpl            ! MPI data
integer,intent(in) :: nc1                      ! Horizontal box size
integer,intent(in) :: nl1                      ! Vertical box size
real(kind_real),intent(in) :: val(nc1,nl1)     ! Box value

! Local variables
integer :: ibd,ic1,il1

! Count non-missing values
balldata%nbd = count(mpl%msv%isnot(val))

! Allocation
call balldata%alloc

! Pack data
ibd = 0
do il1=1,nl1
   do ic1=1,nc1
      if (mpl%msv%isnot(val(ic1,il1))) then
         ibd = ibd+1
         balldata%bd_to_c1(ibd) = ic1
         balldata%bd_to_l1(ibd) = il1
         balldata%val(ibd) = val(ic1,il1)
      end if
   end do
end do

end subroutine balldata_pack

!----------------------------------------------------------------------
! Subroutine: nicas_blk_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine nicas_blk_partial_dealloc(nicas_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block

! Local variables
integer :: isbb

! Release memory
if (allocated(nicas_blk%rhs_avg)) deallocate(nicas_blk%rhs_avg)
if (allocated(nicas_blk%vbot)) deallocate(nicas_blk%vbot)
if (allocated(nicas_blk%vtop)) deallocate(nicas_blk%vtop)
if (allocated(nicas_blk%nc2)) deallocate(nicas_blk%nc2)
if (allocated(nicas_blk%slev)) deallocate(nicas_blk%slev)
if (allocated(nicas_blk%s_to_c1)) deallocate(nicas_blk%s_to_c1)
if (allocated(nicas_blk%s_to_l1)) deallocate(nicas_blk%s_to_l1)
if (allocated(nicas_blk%c1_to_c0)) deallocate(nicas_blk%c1_to_c0)
if (allocated(nicas_blk%l1_to_l0)) deallocate(nicas_blk%l1_to_l0)
if (allocated(nicas_blk%l0_to_l1)) deallocate(nicas_blk%l0_to_l1)
if (allocated(nicas_blk%mask_c1)) deallocate(nicas_blk%mask_c1)
if (allocated(nicas_blk%mask_c2)) deallocate(nicas_blk%mask_c2)
if (allocated(nicas_blk%c1l1_to_s)) deallocate(nicas_blk%c1l1_to_s)
if (allocated(nicas_blk%c1_to_proc)) deallocate(nicas_blk%c1_to_proc)
if (allocated(nicas_blk%s_to_proc)) deallocate(nicas_blk%s_to_proc)
if (allocated(nicas_blk%lcheck_sa)) deallocate(nicas_blk%lcheck_sa)
if (allocated(nicas_blk%lcheck_sb)) deallocate(nicas_blk%lcheck_sb)
if (allocated(nicas_blk%lcheck_sc)) deallocate(nicas_blk%lcheck_sc)
if (allocated(nicas_blk%c1a_to_c1)) deallocate(nicas_blk%c1a_to_c1)
if (allocated(nicas_blk%c1_to_c1a)) deallocate(nicas_blk%c1_to_c1a)
if (allocated(nicas_blk%c1b_to_c1)) deallocate(nicas_blk%c1b_to_c1)
if (allocated(nicas_blk%c1_to_c1b)) deallocate(nicas_blk%c1_to_c1b)
if (allocated(nicas_blk%c1bl1_to_sb)) deallocate(nicas_blk%c1bl1_to_sb)
if (allocated(nicas_blk%c1a_to_c0a)) deallocate(nicas_blk%c1a_to_c0a)
if (allocated(nicas_blk%s_to_sa)) deallocate(nicas_blk%s_to_sa)
if (allocated(nicas_blk%sb_to_s)) deallocate(nicas_blk%sb_to_s)
if (allocated(nicas_blk%sbb_to_s)) deallocate(nicas_blk%sbb_to_s)
if (allocated(nicas_blk%c1_to_c1bb)) deallocate(nicas_blk%c1_to_c1bb)
if (allocated(nicas_blk%c1bb_to_c1)) deallocate(nicas_blk%c1bb_to_c1)
if (allocated(nicas_blk%sc_to_s)) deallocate(nicas_blk%sc_to_s)
if (allocated(nicas_blk%rh_c1)) deallocate(nicas_blk%rh_c1)
if (allocated(nicas_blk%rv_c1)) deallocate(nicas_blk%rv_c1)
if (allocated(nicas_blk%H11_c1)) deallocate(nicas_blk%H11_c1)
if (allocated(nicas_blk%H22_c1)) deallocate(nicas_blk%H22_c1)
if (allocated(nicas_blk%H33_c1)) deallocate(nicas_blk%H33_c1)
if (allocated(nicas_blk%H12_c1)) deallocate(nicas_blk%H12_c1)
if (allocated(nicas_blk%Hcoef)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%Hcoef(isbb)%dealloc
   end do
   deallocate(nicas_blk%Hcoef)
end if

if (allocated(nicas_blk%distnorm)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%distnorm(isbb)%dealloc
   end do
   deallocate(nicas_blk%distnorm)
end if
if (allocated(nicas_blk%distnormv)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%distnormv(isbb)%dealloc
   end do
   deallocate(nicas_blk%distnormv)
end if
if (allocated(nicas_blk%rfac)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%rfac(isbb)%dealloc
   end do
   deallocate(nicas_blk%rfac)
end if
if (allocated(nicas_blk%coef)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%coef(isbb)%dealloc
   end do
   deallocate(nicas_blk%coef)
end if
if (allocated(nicas_blk%sc_nor_to_s)) deallocate(nicas_blk%sc_nor_to_s)
if (allocated(nicas_blk%s_to_sc_nor)) deallocate(nicas_blk%s_to_sc_nor)
if (allocated(nicas_blk%sb_to_sc_nor)) deallocate(nicas_blk%sb_to_sc_nor)
call nicas_blk%c_nor%dealloc
call nicas_blk%tree%dealloc

end subroutine nicas_blk_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: nicas_blk_dealloc
! Purpose: release memory (full)
!----------------------------------------------------------------------
subroutine nicas_blk_dealloc(nicas_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block

! Local variables
integer :: il0,il1,its

! Release memory
call nicas_blk%partial_dealloc
if (allocated(nicas_blk%vlev)) deallocate(nicas_blk%vlev)
if (allocated(nicas_blk%sa_to_s)) deallocate(nicas_blk%sa_to_s)
if (allocated(nicas_blk%sa_to_sc)) deallocate(nicas_blk%sa_to_sc)
if (allocated(nicas_blk%sb_to_sc)) deallocate(nicas_blk%sb_to_sc)
call nicas_blk%c%dealloc
if (allocated(nicas_blk%h)) then
   do il0=1,size(nicas_blk%h)
      call nicas_blk%h(il0)%dealloc
   end do
   deallocate(nicas_blk%h)
end if
call nicas_blk%v%dealloc
if (allocated(nicas_blk%s)) then
   do il1=1,nicas_blk%nl1
     call nicas_blk%s(il1)%dealloc
   end do
   deallocate(nicas_blk%s)
end if
if (allocated(nicas_blk%d)) then
   do its=2,size(nicas_blk%d,2)
      do il0=1,size(nicas_blk%d,1)
        call nicas_blk%d(il0,its)%dealloc
        call nicas_blk%dinv(il0,its)%dealloc
      end do
   end do
   deallocate(nicas_blk%d)
   deallocate(nicas_blk%dinv)
end if
if (allocated(nicas_blk%sb_to_c1b)) deallocate(nicas_blk%sb_to_c1b)
if (allocated(nicas_blk%sb_to_l1)) deallocate(nicas_blk%sb_to_l1)
if (allocated(nicas_blk%norm)) deallocate(nicas_blk%norm)
if (allocated(nicas_blk%coef_ens)) deallocate(nicas_blk%coef_ens)
call nicas_blk%com_AB%dealloc
call nicas_blk%com_AC%dealloc
call nicas_blk%com_AD%dealloc
call nicas_blk%com_ADinv%dealloc
if (allocated(nicas_blk%lon_sa)) deallocate(nicas_blk%lon_sa)
if (allocated(nicas_blk%lat_sa)) deallocate(nicas_blk%lat_sa)
if (allocated(nicas_blk%lev_sa)) deallocate(nicas_blk%lev_sa)
if (allocated(nicas_blk%lon_sb)) deallocate(nicas_blk%lon_sb)
if (allocated(nicas_blk%lat_sb)) deallocate(nicas_blk%lat_sb)
if (allocated(nicas_blk%lev_sb)) deallocate(nicas_blk%lev_sb)
if (allocated(nicas_blk%lon_sc)) deallocate(nicas_blk%lon_sc)
if (allocated(nicas_blk%lat_sc)) deallocate(nicas_blk%lat_sc)
if (allocated(nicas_blk%lev_sc)) deallocate(nicas_blk%lev_sc)


end subroutine nicas_blk_dealloc

!----------------------------------------------------------------------
! Subroutine: nicas_blk_read
! Purpose: read
!----------------------------------------------------------------------
subroutine nicas_blk_read(nicas_blk,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(bpar_type),intent(in) :: bpar               ! Block parameters

! Local variables
integer :: il0i,il1,its,il0,info
integer :: ncid,nl0_id,nc0a_id,nc1b_id,nl1_id,nsa_id,nsb_id,nc0d_id,nc0dinv_id
integer :: vlev_id,sb_to_c1b_id,sb_to_l1_id,sa_to_s_id,sa_to_sc_id,sb_to_sc_id,norm_id,coef_ens_id
integer :: vlev_int(geom%nl0)
character(len=2*1024+1) :: filename
character(len=1024),parameter :: subr = 'nicas_blk_read'

! Associate
associate(ib=>nicas_blk%ib)

! Open file and get dimensions
filename = trim(nam%prefix)//'_'//trim(nicas_blk%name)
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))

! Get dimensions
nl0_id = mpl%ncdimcheck(subr,ncid,'nl0',geom%nl0,.false.)
if (bpar%nicas_block(ib)) then
   info = nf90_inq_dimid(ncid,'nc0a',nc0a_id)
   if (info==nf90_noerr) then
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc0a_id,len=nicas_blk%nc0a))
   else
      nicas_blk%nc0a = 0
   end if
   info = nf90_inq_dimid(ncid,'nc1b',nc1b_id)
   if (info==nf90_noerr) then
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc1b_id,len=nicas_blk%nc1b))
   else
      nicas_blk%nc1b = 0
   end if
   call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nl1',nl1_id))
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nl1_id,len=nicas_blk%nl1))
   info = nf90_inq_dimid(ncid,'nsa',nsa_id)
   if (info==nf90_noerr) then
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nsa_id,len=nicas_blk%nsa))
   else
      nicas_blk%nsa = 0
   end if
   info = nf90_inq_dimid(ncid,'nsb',nsb_id)
   if (info==nf90_noerr) then
      call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nsb_id,len=nicas_blk%nsb))
   else
      nicas_blk%nsb = 0
   end if
   call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'nsc',nicas_blk%nsc))
end if
if ((ib==bpar%nbe).and.(abs(nam%adv_mode)==1)) then
   call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc0d',nc0d_id))
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc0d_id,len=nicas_blk%nc0d))
   call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nc0dinv',nc0dinv_id))
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nc0dinv_id,len=nicas_blk%nc0dinv))
end if

! Allocation
if (bpar%nicas_block(ib)) then
   allocate(nicas_blk%vlev(geom%nl0))
   allocate(nicas_blk%norm(nicas_blk%nc0a,geom%nl0))
   allocate(nicas_blk%coef_ens(nicas_blk%nc0a,geom%nl0))
   if (nicas_blk%nsa>0) then
      allocate(nicas_blk%sa_to_s(nicas_blk%nsa))
      allocate(nicas_blk%sa_to_sc(nicas_blk%nsa))
   end if
   if (nicas_blk%nsb>0) then
      allocate(nicas_blk%sb_to_c1b(nicas_blk%nsb))
      allocate(nicas_blk%sb_to_l1(nicas_blk%nsb))
      allocate(nicas_blk%sb_to_sc(nicas_blk%nsb))
   end if
   allocate(nicas_blk%h(geom%nl0i))
   allocate(nicas_blk%s(nicas_blk%nl1))
end if
if ((ib==bpar%nbe).and.(abs(nam%adv_mode)==1)) then
   allocate(nicas_blk%d(geom%nl0,2:nam%nts))
   allocate(nicas_blk%dinv(geom%nl0,2:nam%nts))
end if

! Get variable id
if (bpar%nicas_block(ib)) then
   call mpl%ncerr(subr,nf90_inq_varid(ncid,'vlev',vlev_id))
   if (nicas_blk%nc0a>0) then
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'norm',norm_id))
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'coef_ens',coef_ens_id))
   end if
   if (nicas_blk%nsa>0) then
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'sa_to_s',sa_to_s_id))
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'sa_to_sc',sa_to_sc_id))
   end if
   if (nicas_blk%nsb>0) then
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'sb_to_c1b',sb_to_c1b_id))
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'sb_to_l1',sb_to_l1_id))
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'sb_to_sc',sb_to_sc_id))
   end if
end if

! Read data
if (bpar%nicas_block(ib)) then
   call mpl%ncerr(subr,nf90_get_var(ncid,vlev_id,vlev_int))
   do il0=1,geom%nl0
      if (vlev_int(il0)==0) then
         nicas_blk%vlev(il0) = .false.
      elseif (vlev_int(il0)==1) then
         nicas_blk%vlev(il0) = .true.
      else
         call mpl%abort(subr,'wrong vlev')
      end if
   end do
   if (nicas_blk%nc0a>0) then
      call mpl%ncerr(subr,nf90_get_var(ncid,norm_id,nicas_blk%norm))
      call mpl%ncerr(subr,nf90_get_var(ncid,coef_ens_id,nicas_blk%coef_ens))
   end if
   if (nicas_blk%nsa>0) then
      call mpl%ncerr(subr,nf90_get_var(ncid,sa_to_s_id,nicas_blk%sa_to_s))
      call mpl%ncerr(subr,nf90_get_var(ncid,sa_to_sc_id,nicas_blk%sa_to_sc))
   end if
   if (nicas_blk%nsb>0) then
      call mpl%ncerr(subr,nf90_get_var(ncid,sb_to_c1b_id,nicas_blk%sb_to_c1b))
      call mpl%ncerr(subr,nf90_get_var(ncid,sb_to_l1_id,nicas_blk%sb_to_l1))
      call mpl%ncerr(subr,nf90_get_var(ncid,sb_to_sc_id,nicas_blk%sb_to_sc))
   end if
   nicas_blk%com_AB%prefix = 'com_AB'
   call nicas_blk%com_AB%read(mpl,ncid)
   nicas_blk%com_AC%prefix = 'com_AC'
   call nicas_blk%com_AC%read(mpl,ncid)
   nicas_blk%c%prefix = 'c'
   call nicas_blk%c%read(mpl,ncid)
   do il0i=1,geom%nl0i
      write(nicas_blk%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
      call nicas_blk%h(il0i)%read(mpl,ncid)
   end do
   nicas_blk%v%prefix = 'v'
   call nicas_blk%v%read(mpl,ncid)
   do il1=1,nicas_blk%nl1
      write(nicas_blk%s(il1)%prefix,'(a,i3.3)') 's_',il1
      call nicas_blk%s(il1)%read(mpl,ncid)
   end do
end if
if ((ib==bpar%nbe).and.(abs(nam%adv_mode)==1)) then
   nicas_blk%com_AD%prefix = 'com_AD'
   call nicas_blk%com_AD%read(mpl,ncid)
   nicas_blk%com_ADinv%prefix = 'com_ADinv'
   call nicas_blk%com_ADinv%read(mpl,ncid)
   do its=2,nam%nts
      do il0=1,geom%nl0
         write(nicas_blk%d(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
         call nicas_blk%d(il0,its)%read(mpl,ncid)
         write(nicas_blk%dinv(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'dinv_',il0,'_',its
         call nicas_blk%dinv(il0,its)%read(mpl,ncid)
      end do
   end do
end if

! Read main weight
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'wgt',nicas_blk%wgt))

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine nicas_blk_read

!----------------------------------------------------------------------
! Subroutine: nicas_blk_write
! Purpose: write
!----------------------------------------------------------------------
subroutine nicas_blk_write(nicas_blk,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(nam_type),intent(in) :: nam              ! Namelist
type(geom_type),intent(in) :: geom            ! Geometry
type(bpar_type),intent(in) :: bpar            ! Block parameters

! Local variables
integer :: il0i,il1,its,il0
integer :: ncid,nl0_id,nc0a_id,nc1b_id,nl1_id,nsa_id,nsb_id,nc0d_id,nc0dinv_id
integer :: vlev_id,sb_to_c1b_id,sb_to_l1_id,sa_to_s_id,sa_to_sc_id,sb_to_sc_id,norm_id,coef_ens_id
integer :: vlev_int(geom%nl0)
character(len=2*1024+1) :: filename
character(len=1024),parameter :: subr = 'nicas_blk_write'

! Associate
associate(ib=>nicas_blk%ib)

! Create file
filename = trim(nam%prefix)//'_'//trim(nicas_blk%name)
call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call nam%write(mpl,ncid)

! Define dimensions
call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
if (bpar%nicas_block(ib)) then
   if (nicas_blk%nc0a>0) call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0a',nicas_blk%nc0a,nc0a_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0i',geom%nl0i))
   if (nicas_blk%nc1b>0) call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1b',nicas_blk%nc1b,nc1b_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nl1',nicas_blk%nl1,nl1_id))
   if (nicas_blk%nsa>0) call mpl%ncerr(subr,nf90_def_dim(ncid,'nsa',nicas_blk%nsa,nsa_id))
   if (nicas_blk%nsb>0) call mpl%ncerr(subr,nf90_def_dim(ncid,'nsb',nicas_blk%nsb,nsb_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'nsc',nicas_blk%nsc))
end if
if ((ib==bpar%nbe).and.nam%adv_diag) then
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0d',nicas_blk%nc0d,nc0d_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0dinv',nicas_blk%nc0dinv,nc0dinv_id))
end if

! Write main weight
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'wgt',nicas_blk%wgt))

! Define variables
if (bpar%nicas_block(ib)) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'vlev',nf90_int,(/nl0_id/),vlev_id))
   if (nicas_blk%nc0a>0) then
      call mpl%ncerr(subr,nf90_def_var(ncid,'norm',nc_kind_real,(/nc0a_id,nl0_id/),norm_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'coef_ens',nc_kind_real,(/nc0a_id,nl0_id/),coef_ens_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,norm_id,'_FillValue',mpl%msv%valr))
      call mpl%ncerr(subr,nf90_put_att(ncid,coef_ens_id,'_FillValue',mpl%msv%valr))
   end if
   if (nicas_blk%nsa>0) then
      call mpl%ncerr(subr,nf90_def_var(ncid,'sa_to_s',nf90_int,(/nsa_id/),sa_to_s_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'sa_to_sc',nf90_int,(/nsa_id/),sa_to_sc_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,sa_to_s_id,'_FillValue',mpl%msv%vali))
      call mpl%ncerr(subr,nf90_put_att(ncid,sa_to_sc_id,'_FillValue',mpl%msv%vali))
   end if
   if (nicas_blk%nsb>0) then
      call mpl%ncerr(subr,nf90_def_var(ncid,'sb_to_c1b',nf90_int,(/nsb_id/),sb_to_c1b_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'sb_to_l1',nf90_int,(/nsb_id/),sb_to_l1_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'sb_to_sc',nf90_int,(/nsb_id/),sb_to_sc_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,sb_to_sc_id,'_FillValue',mpl%msv%vali))
   end if
end if

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Write variables
if (bpar%nicas_block(ib)) then
   do il0=1,geom%nl0
      if (nicas_blk%vlev(il0)) then
         vlev_int(il0) = 1
      else
         vlev_int(il0) = 0
      end if
   end do
   call mpl%ncerr(subr,nf90_put_var(ncid,vlev_id,vlev_int))
   if (nicas_blk%nc0a>0) then
      call mpl%ncerr(subr,nf90_put_var(ncid,norm_id,nicas_blk%norm))
      call mpl%ncerr(subr,nf90_put_var(ncid,coef_ens_id,nicas_blk%coef_ens))
   end if
   if (nicas_blk%nsa>0) then
      call mpl%ncerr(subr,nf90_put_var(ncid,sa_to_s_id,nicas_blk%sa_to_s))
      call mpl%ncerr(subr,nf90_put_var(ncid,sa_to_sc_id,nicas_blk%sa_to_sc))
   end if
   if (nicas_blk%nsb>0) then
      call mpl%ncerr(subr,nf90_put_var(ncid,sb_to_c1b_id,nicas_blk%sb_to_c1b))
      call mpl%ncerr(subr,nf90_put_var(ncid,sb_to_l1_id,nicas_blk%sb_to_l1))
      call mpl%ncerr(subr,nf90_put_var(ncid,sb_to_sc_id,nicas_blk%sb_to_sc))
   end if
   call nicas_blk%com_AB%write(mpl,ncid)
   call nicas_blk%com_AC%write(mpl,ncid)
   call nicas_blk%c%write(mpl,ncid)
   do il0i=1,geom%nl0i
      call nicas_blk%h(il0i)%write(mpl,ncid)
   end do
   call nicas_blk%v%write(mpl,ncid)
   do il1=1,nicas_blk%nl1
      call nicas_blk%s(il1)%write(mpl,ncid)
   end do
end if
if ((ib==bpar%nbe).and.nam%adv_diag) then
   call nicas_blk%com_AD%write(mpl,ncid)
   call nicas_blk%com_ADinv%write(mpl,ncid)
   do its=2,nam%nts
      do il0=1,geom%nl0
         call nicas_blk%d(il0,its)%write(mpl,ncid)
         call nicas_blk%dinv(il0,its)%write(mpl,ncid)
      end do
   end do
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine nicas_blk_write

!----------------------------------------------------------------------
! Subroutine: nicas_blk_write_grids
! Purpose: write NICAS grids
!----------------------------------------------------------------------
subroutine nicas_blk_write_grids(nicas_blk,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(nam_type),intent(in) :: nam              ! Namelist
type(geom_type),intent(in) :: geom            ! Geometry
type(bpar_type),intent(in) :: bpar            ! Block parameters

! Local variables
integer :: ncid,nsa_id,nsb_id,nsc_id,lon_sa_id,lat_sa_id,lev_sa_id,lon_sb_id,lat_sb_id,lev_sb_id,lon_sc_id,lat_sc_id,lev_sc_id
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'nicas_blk_write_grids'

! Associate
associate(ib=>nicas_blk%ib)

! Create file
filename = trim(nam%prefix)//'_'//trim(nicas_blk%name)//'_grids'
call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call nam%write(mpl,ncid)

! Define dimensions
if (bpar%nicas_block(ib)) then
   if (nicas_blk%nsa>0) call mpl%ncerr(subr,nf90_def_dim(ncid,'nsa',nicas_blk%nsa,nsa_id))
   if (nicas_blk%nsb>0) call mpl%ncerr(subr,nf90_def_dim(ncid,'nsb',nicas_blk%nsb,nsb_id))
   if (nicas_blk%nsc>0) call mpl%ncerr(subr,nf90_def_dim(ncid,'nsc',nicas_blk%nsc,nsc_id))
end if

! Define variables
if (bpar%nicas_block(ib)) then
   if (nicas_blk%nsa>0) then
      call mpl%ncerr(subr,nf90_def_var(ncid,'lon_sa',nc_kind_real,(/nsa_id/),lon_sa_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lat_sa',nc_kind_real,(/nsa_id/),lat_sa_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lev_sa',nf90_int,(/nsa_id/),lev_sa_id))
   end if
   if (nicas_blk%nsb>0) then
      call mpl%ncerr(subr,nf90_def_var(ncid,'lon_sb',nc_kind_real,(/nsb_id/),lon_sb_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lat_sb',nc_kind_real,(/nsb_id/),lat_sb_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lev_sb',nf90_int,(/nsb_id/),lev_sb_id))
   end if
   if (nicas_blk%nsc>0) then
      call mpl%ncerr(subr,nf90_def_var(ncid,'lon_sc',nc_kind_real,(/nsc_id/),lon_sc_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lat_sc',nc_kind_real,(/nsc_id/),lat_sc_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lev_sc',nf90_int,(/nsc_id/),lev_sc_id))
   end if
end if

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Write variables
if (bpar%nicas_block(ib)) then
   if (nicas_blk%nsa>0) then
      call mpl%ncerr(subr,nf90_put_var(ncid,lon_sa_id,nicas_blk%lon_sa))
      call mpl%ncerr(subr,nf90_put_var(ncid,lat_sa_id,nicas_blk%lat_sa))
      call mpl%ncerr(subr,nf90_put_var(ncid,lev_sa_id,nicas_blk%lev_sa))
   end if
   if (nicas_blk%nsb>0) then
      call mpl%ncerr(subr,nf90_put_var(ncid,lon_sb_id,nicas_blk%lon_sb))
      call mpl%ncerr(subr,nf90_put_var(ncid,lat_sb_id,nicas_blk%lat_sb))
      call mpl%ncerr(subr,nf90_put_var(ncid,lev_sb_id,nicas_blk%lev_sb))
   end if
   if (nicas_blk%nsc>0) then
      call mpl%ncerr(subr,nf90_put_var(ncid,lon_sc_id,nicas_blk%lon_sc))
      call mpl%ncerr(subr,nf90_put_var(ncid,lat_sc_id,nicas_blk%lat_sc))
      call mpl%ncerr(subr,nf90_put_var(ncid,lev_sc_id,nicas_blk%lev_sc))
   end if
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))
! End associate
end associate

end subroutine nicas_blk_write_grids

!----------------------------------------------------------------------
! Subroutine: nicas_blk_receive
! Purpose: receive
!----------------------------------------------------------------------
subroutine nicas_blk_receive(nicas_blk,mpl,nam,geom,bpar,iproc,tag_start)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(bpar_type),intent(in) :: bpar               ! Block parameters
integer,intent(in) :: iproc                      ! Source task index
integer,intent(in) :: tag_start                  ! MPI tag

! Local variables
integer :: tag,n_dim,n_int,n_real,n_logical,offset_int,offset_real,offset_logical,il0i,il1,its,il0
integer,allocatable :: rbuf_dim(:),rbuf_int(:)
real(kind_real),allocatable :: rbuf_real(:)
logical,allocatable :: rbuf_logical(:),mask_c0a(:,:)
type(fckit_mpi_status) :: status

! Associate
associate(ib=>nicas_blk%ib)

! Initialization
tag = tag_start
n_dim = 0
n_int = 0
n_real = 0
n_logical = 0
offset_int = 0
offset_real = 0
offset_logical = 0

! Define buffer size
if (bpar%nicas_block(ib)) n_dim = n_dim+6
if ((ib==bpar%nbe).and.nam%adv_diag) n_dim = n_dim+2

! Allocation
allocate(rbuf_dim(n_dim))

! Receive buffer
call mpl%f_comm%receive(rbuf_dim,iproc-1,tag,status)
tag = tag+1

! Copy data
if (bpar%nicas_block(ib)) then
   nicas_blk%nc0a = rbuf_dim(1)
   nicas_blk%nc1b = rbuf_dim(2)
   nicas_blk%nl1 = rbuf_dim(3)
   nicas_blk%nsa = rbuf_dim(4)
   nicas_blk%nsb = rbuf_dim(5)
   nicas_blk%nsc = rbuf_dim(6)
end if
if ((ib==bpar%nbe).and.nam%adv_diag) then
   nicas_blk%nc0d = rbuf_dim(7)
   nicas_blk%nc0dinv = rbuf_dim(8)
end if

! Release memory
deallocate(rbuf_dim)

! Allocation
if (bpar%nicas_block(ib)) then
   allocate(nicas_blk%vlev(geom%nl0))
   allocate(nicas_blk%norm(nicas_blk%nc0a,geom%nl0))
   allocate(nicas_blk%coef_ens(nicas_blk%nc0a,geom%nl0))
   if (nicas_blk%nsa>0) then
      allocate(nicas_blk%sa_to_s(nicas_blk%nsa))
      allocate(nicas_blk%sa_to_sc(nicas_blk%nsa))
   end if
   if (nicas_blk%nsb>0) then
      allocate(nicas_blk%sb_to_c1b(nicas_blk%nsb))
      allocate(nicas_blk%sb_to_l1(nicas_blk%nsb))
      allocate(nicas_blk%sb_to_sc(nicas_blk%nsb))
   end if
   if (nam%write_grids) then
      if (nicas_blk%nsa>0) then
         allocate(nicas_blk%lon_sa(nicas_blk%nsa))
         allocate(nicas_blk%lat_sa(nicas_blk%nsa))
         allocate(nicas_blk%lev_sa(nicas_blk%nsa))
      end if
      if (nicas_blk%nsb>0) then
         allocate(nicas_blk%lon_sb(nicas_blk%nsb))
         allocate(nicas_blk%lat_sb(nicas_blk%nsb))
         allocate(nicas_blk%lev_sb(nicas_blk%nsb))
      end if
      if (nicas_blk%nsc>0) then
         allocate(nicas_blk%lon_sc(nicas_blk%nsc))
         allocate(nicas_blk%lat_sc(nicas_blk%nsc))
         allocate(nicas_blk%lev_sc(nicas_blk%nsc))
      end if
   end if
   allocate(nicas_blk%h(geom%nl0i))
   allocate(nicas_blk%s(nicas_blk%nl1))
end if
if ((ib==bpar%nbe).and.(abs(nam%adv_mode)==1)) then
   allocate(nicas_blk%d(geom%nl0,2:nam%nts))
   allocate(nicas_blk%dinv(geom%nl0,2:nam%nts))
end if

! Define buffer sizes
n_real = n_real+1
if (bpar%nicas_block(ib)) then
   n_int = n_int+2*nicas_blk%nsa+3*nicas_blk%nsb
   n_real = n_real+2*nicas_blk%nc0a*geom%nl0
   n_logical = n_logical+geom%nl0
   if (nam%write_grids) then
      n_int = n_int+nicas_blk%nsa+nicas_blk%nsb+nicas_blk%nsc
      n_real = n_real+2*(nicas_blk%nsa+nicas_blk%nsb+nicas_blk%nsc)
   end if
end if

! Allocation
allocate(rbuf_int(n_int))
allocate(rbuf_real(n_real))
allocate(rbuf_logical(n_logical))
if (bpar%nicas_block(ib)) allocate(mask_c0a(nicas_blk%nc0a,geom%nl0))

! Initialization
if (bpar%nicas_block(ib)) mask_c0a = .true.

! Receive buffers
if (n_int>0) call mpl%f_comm%receive(rbuf_int,iproc-1,tag,status)
tag = tag+1
if (n_real>0) call mpl%f_comm%receive(rbuf_real,iproc-1,tag,status)
tag = tag+1
if (n_logical>0) call mpl%f_comm%receive(rbuf_logical,iproc-1,tag,status)
tag = tag+1

! Copy data
nicas_blk%wgt = rbuf_real(1)
offset_real = offset_real+1
if (bpar%nicas_block(ib)) then
   nicas_blk%vlev = rbuf_logical(offset_logical+1:offset_logical+geom%nl0)
   offset_logical = offset_logical+geom%nl0
   if (nicas_blk%nc0a>0) then
      nicas_blk%norm = unpack(rbuf_real(offset_real+1:offset_real+nicas_blk%nc0a*geom%nl0),mask_c0a,nicas_blk%norm)
      offset_real = offset_real+nicas_blk%nc0a*geom%nl0
      nicas_blk%coef_ens = unpack(rbuf_real(offset_real+1:offset_real+nicas_blk%nc0a*geom%nl0),mask_c0a,nicas_blk%coef_ens)
      offset_real = offset_real+nicas_blk%nc0a*geom%nl0
   end if
   if (nicas_blk%nsa>0) then
      nicas_blk%sa_to_s = rbuf_int(offset_int+1:offset_int+nicas_blk%nsa)
      offset_int = offset_int+nicas_blk%nsa
      nicas_blk%sa_to_sc = rbuf_int(offset_int+1:offset_int+nicas_blk%nsa)
      offset_int = offset_int+nicas_blk%nsa
   end if
   if (nicas_blk%nsb>0) then
      nicas_blk%sb_to_c1b = rbuf_int(offset_int+1:offset_int+nicas_blk%nsb)
      offset_int = offset_int+nicas_blk%nsb
      nicas_blk%sb_to_l1 = rbuf_int(offset_int+1:offset_int+nicas_blk%nsb)
      offset_int = offset_int+nicas_blk%nsb
      nicas_blk%sb_to_sc = rbuf_int(offset_int+1:offset_int+nicas_blk%nsb)
      offset_int = offset_int+nicas_blk%nsb
   end if
   if (nam%write_grids) then
      if (nicas_blk%nsa>0) then
         nicas_blk%lon_sa = rbuf_real(offset_real+1:offset_real+nicas_blk%nsa)
         offset_real = offset_real+nicas_blk%nsa
         nicas_blk%lat_sa = rbuf_real(offset_real+1:offset_real+nicas_blk%nsa)
         offset_real = offset_real+nicas_blk%nsa
         nicas_blk%lev_sa = rbuf_int(offset_int+1:offset_int+nicas_blk%nsa)
         offset_int = offset_int+nicas_blk%nsa
      end if
      if (nicas_blk%nsb>0) then
         nicas_blk%lon_sb = rbuf_real(offset_real+1:offset_real+nicas_blk%nsb)
         offset_real = offset_real+nicas_blk%nsb
         nicas_blk%lat_sb = rbuf_real(offset_real+1:offset_real+nicas_blk%nsb)
         offset_real = offset_real+nicas_blk%nsb
         nicas_blk%lev_sb = rbuf_int(offset_int+1:offset_int+nicas_blk%nsb)
         offset_int = offset_int+nicas_blk%nsb
      end if
      if (nicas_blk%nsc>0) then
         nicas_blk%lon_sc = rbuf_real(offset_real+1:offset_real+nicas_blk%nsc)
         offset_real = offset_real+nicas_blk%nsc
         nicas_blk%lat_sc = rbuf_real(offset_real+1:offset_real+nicas_blk%nsc)
         offset_real = offset_real+nicas_blk%nsc
         nicas_blk%lev_sc = rbuf_int(offset_int+1:offset_int+nicas_blk%nsc)
         offset_int = offset_int+nicas_blk%nsc
      end if
   end if
end if

! Receive communications
if (bpar%nicas_block(ib)) then
   nicas_blk%com_AB%prefix = 'com_AB'
   call nicas_blk%com_AB%receive(mpl,iproc,tag)
   tag = tag+com_ntag
   nicas_blk%com_AC%prefix = 'com_AC'
   call nicas_blk%com_AC%receive(mpl,iproc,tag)
   tag = tag+com_ntag
end if

! Receive linear operators
if (bpar%nicas_block(ib)) then
   nicas_blk%c%prefix = 'c'
   call nicas_blk%c%receive(mpl,iproc,tag)
   tag = tag+linop_ntag
   do il0i=1,geom%nl0i
      write(nicas_blk%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
      call nicas_blk%h(il0i)%receive(mpl,iproc,tag)
      tag = tag+linop_ntag
   end do
   nicas_blk%v%prefix = 'v'
   call nicas_blk%v%receive(mpl,iproc,tag)
   tag = tag+linop_ntag
   do il1=1,nicas_blk%nl1
      write(nicas_blk%s(il1)%prefix,'(a,i3.3)') 's_',il1
      call nicas_blk%s(il1)%receive(mpl,iproc,tag)
      tag = tag+linop_ntag
   end do
end if
if ((ib==bpar%nbe).and.nam%adv_diag) then
   nicas_blk%com_AD%prefix = 'com_AD'
   call nicas_blk%com_AD%receive(mpl,iproc,tag)
   tag = tag+com_ntag
   nicas_blk%com_ADinv%prefix = 'com_ADinv'
   call nicas_blk%com_ADinv%receive(mpl,iproc,tag)
   tag = tag+com_ntag
   do its=2,nam%nts
      do il0=1,geom%nl0
         write(nicas_blk%d(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
         call nicas_blk%d(il0,its)%receive(mpl,iproc,tag)
         tag = tag+linop_ntag
         write(nicas_blk%dinv(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'dinv_',il0,'_',its
         call nicas_blk%dinv(il0,its)%receive(mpl,iproc,tag)
         tag = tag+linop_ntag
      end do
   end do
end if

! Release memory
deallocate(rbuf_int)
deallocate(rbuf_real)
deallocate(rbuf_logical)
if (bpar%nicas_block(ib)) deallocate(mask_c0a)

! End associate
end associate

end subroutine nicas_blk_receive

!----------------------------------------------------------------------
! Subroutine: nicas_blk_send
! Purpose: send
!----------------------------------------------------------------------
subroutine nicas_blk_send(nicas_blk,mpl,nam,geom,bpar,iproc,tag_start)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(nam_type),intent(in) :: nam              ! Namelist
type(geom_type),intent(in) :: geom            ! Geometry
type(bpar_type),intent(in) :: bpar            ! Block parameters
integer,intent(in) :: iproc                   ! Destination task index
integer,intent(in) :: tag_start               ! MPI tag

! Local variables
integer :: tag,n_dim,n_int,n_real,n_logical,offset_int,offset_real,offset_logical,il0i,il1,its,il0
integer,allocatable :: sbuf_dim(:),sbuf_int(:)
real(kind_real),allocatable :: sbuf_real(:)
logical,allocatable :: sbuf_logical(:),mask_c0a(:,:)

! Associate
associate(ib=>nicas_blk%ib)

! Initialization
tag = tag_start
n_dim = 0
n_int = 0
n_real = 0
n_logical = 0
offset_int = 0
offset_real = 0
offset_logical = 0

! Define buffer size
if (bpar%nicas_block(ib)) n_dim = n_dim+6
if ((ib==bpar%nbe).and.nam%adv_diag) n_dim = n_dim+2

! Allocation
allocate(sbuf_dim(n_dim))

! Copy data
if (bpar%nicas_block(ib)) then
   sbuf_dim(1) = nicas_blk%nc0a
   sbuf_dim(2) = nicas_blk%nc1b
   sbuf_dim(3) = nicas_blk%nl1
   sbuf_dim(4) = nicas_blk%nsa
   sbuf_dim(5) = nicas_blk%nsb
   sbuf_dim(6) = nicas_blk%nsc
end if
if ((ib==bpar%nbe).and.nam%adv_diag) then
   sbuf_dim(7) = nicas_blk%nc0d
   sbuf_dim(8) = nicas_blk%nc0dinv
end if

! Send buffer
call mpl%f_comm%send(sbuf_dim,iproc-1,tag)
tag = tag+1

! Release memory
deallocate(sbuf_dim)

! Define buffer sizes
n_real = n_real+1
if (bpar%nicas_block(ib)) then
   n_int = n_int+2*nicas_blk%nsa+3*nicas_blk%nsb
   n_real = n_real+2*nicas_blk%nc0a*geom%nl0
   n_logical = n_logical+geom%nl0
   if (nam%write_grids) then
      n_int = n_int+nicas_blk%nsa+nicas_blk%nsb+nicas_blk%nsc
      n_real = n_real+2*(nicas_blk%nsa+nicas_blk%nsb+nicas_blk%nsc)
   end if
end if

! Allocation
allocate(sbuf_int(n_int))
allocate(sbuf_real(n_real))
allocate(sbuf_logical(n_logical))
if (bpar%nicas_block(ib)) allocate(mask_c0a(nicas_blk%nc0a,geom%nl0))

! Initialization
if (bpar%nicas_block(ib)) mask_c0a = .true.

! Copy data
sbuf_real(1) = nicas_blk%wgt
offset_real = offset_real+1
if (bpar%nicas_block(ib)) then
   sbuf_logical(offset_logical+1:offset_logical+geom%nl0) = nicas_blk%vlev
   offset_logical = offset_logical+geom%nl0
   if (nicas_blk%nc0a>0) then
      sbuf_real(offset_real+1:offset_real+nicas_blk%nc0a*geom%nl0) = pack(nicas_blk%norm,mask_c0a)
      offset_real = offset_real+nicas_blk%nc0a*geom%nl0
      sbuf_real(offset_real+1:offset_real+nicas_blk%nc0a*geom%nl0) = pack(nicas_blk%coef_ens,mask_c0a)
      offset_real = offset_real+nicas_blk%nc0a*geom%nl0
   end if
   if (nicas_blk%nsa>0) then
      sbuf_int(offset_int+1:offset_int+nicas_blk%nsa) = nicas_blk%sa_to_s
      offset_int = offset_int+nicas_blk%nsa
      sbuf_int(offset_int+1:offset_int+nicas_blk%nsa) = nicas_blk%sa_to_sc
      offset_int = offset_int+nicas_blk%nsa
   end if
   if (nicas_blk%nsb>0) then
      sbuf_int(offset_int+1:offset_int+nicas_blk%nsb) = nicas_blk%sb_to_c1b
      offset_int = offset_int+nicas_blk%nsb
      sbuf_int(offset_int+1:offset_int+nicas_blk%nsb) = nicas_blk%sb_to_l1
      offset_int = offset_int+nicas_blk%nsb
      sbuf_int(offset_int+1:offset_int+nicas_blk%nsb) = nicas_blk%sb_to_sc
      offset_int = offset_int+nicas_blk%nsb
   end if
   if (nam%write_grids) then
      if (nicas_blk%nsa>0) then
         sbuf_real(offset_real+1:offset_real+nicas_blk%nsa) = nicas_blk%lon_sa
         offset_real = offset_real+nicas_blk%nsa
         sbuf_real(offset_real+1:offset_real+nicas_blk%nsa) = nicas_blk%lat_sa
         offset_real = offset_real+nicas_blk%nsa
         sbuf_int(offset_int+1:offset_int+nicas_blk%nsa) = nicas_blk%lev_sa
         offset_int = offset_int+nicas_blk%nsa
      end if
      if (nicas_blk%nsb>0) then
         sbuf_real(offset_real+1:offset_real+nicas_blk%nsb) = nicas_blk%lon_sb
         offset_real = offset_real+nicas_blk%nsb
         sbuf_real(offset_real+1:offset_real+nicas_blk%nsb) = nicas_blk%lat_sb
         offset_real = offset_real+nicas_blk%nsb
         sbuf_int(offset_int+1:offset_int+nicas_blk%nsb) = nicas_blk%lev_sb
         offset_int = offset_int+nicas_blk%nsb
      end if
      if (nicas_blk%nsc>0) then
         sbuf_real(offset_real+1:offset_real+nicas_blk%nsc) = nicas_blk%lon_sc
         offset_real = offset_real+nicas_blk%nsc
         sbuf_real(offset_real+1:offset_real+nicas_blk%nsc) = nicas_blk%lat_sc
         offset_real = offset_real+nicas_blk%nsc
         sbuf_int(offset_int+1:offset_int+nicas_blk%nsc) = nicas_blk%lev_sc
         offset_int = offset_int+nicas_blk%nsc
      end if
   end if
end if

! Send buffers
if (n_int>0) call mpl%f_comm%send(sbuf_int,iproc-1,tag)
tag = tag+1
if (n_real>0) call mpl%f_comm%send(sbuf_real,iproc-1,tag)
tag = tag+1
if (n_logical>0) call mpl%f_comm%send(sbuf_logical,iproc-1,tag)
tag = tag+1

! Release memory
deallocate(sbuf_int)
deallocate(sbuf_real)
deallocate(sbuf_logical)
if (bpar%nicas_block(ib)) deallocate(mask_c0a)

! Send communications
if (bpar%nicas_block(ib)) then
   call nicas_blk%com_AB%send(mpl,iproc,tag)
   tag = tag+com_ntag
   call nicas_blk%com_AC%send(mpl,iproc,tag)
   tag = tag+com_ntag
end if

! Send linear operators
if (bpar%nicas_block(ib)) then
   call nicas_blk%c%send(mpl,iproc,tag)
   tag = tag+linop_ntag
   do il0i=1,geom%nl0i
      call nicas_blk%h(il0i)%send(mpl,iproc,tag)
      tag = tag+linop_ntag
   end do
   call nicas_blk%v%send(mpl,iproc,tag)
   tag = tag+linop_ntag
   do il1=1,nicas_blk%nl1
      call nicas_blk%s(il1)%send(mpl,iproc,tag)
      tag = tag+linop_ntag
   end do
end if
if ((ib==bpar%nbe).and.nam%adv_diag) then
   call nicas_blk%com_AD%send(mpl,iproc,tag)
   tag = tag+com_ntag
   call nicas_blk%com_ADinv%send(mpl,iproc,tag)
   tag = tag+com_ntag
   do its=2,nam%nts
      do il0=1,geom%nl0
         call nicas_blk%d(il0,its)%send(mpl,iproc,tag)
         tag = tag+linop_ntag
         call nicas_blk%dinv(il0,its)%send(mpl,iproc,tag)
         tag = tag+linop_ntag
      end do
   end do
end if

! End associate
end associate

end subroutine nicas_blk_send

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_parameters
! Purpose: compute NICAS parameters
!----------------------------------------------------------------------
subroutine nicas_blk_compute_parameters(nicas_blk,mpl,rng,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(rng_type),intent(inout) :: rng              ! Random number generator
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       ! C matrix data block

! Local variables
integer :: il0i,il1

! Copy subset Sc0 on halo size A
nicas_blk%nc0a = geom%nc0a

! Compute adaptive sampling, subset Sc1
write(mpl%info,'(a7,a)') '','Compute adaptive sampling, subset Sc1'
call mpl%flush
call nicas_blk%compute_sampling_c1(mpl,rng,nam,geom,cmat_blk)

! Compute MPI distribution, halos A
write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo A'
call mpl%flush
call nicas_blk%compute_mpi_a(mpl,geom)

! Compute adaptive sampling, vertical
write(mpl%info,'(a7,a)') '','Compute adaptive sampling, vertical'
call mpl%flush
call nicas_blk%compute_sampling_v(mpl,nam,geom,cmat_blk)

! Compute adaptive sampling, subset Sc2
write(mpl%info,'(a7,a)') '','Compute adaptive sampling, subset Sc2'
call mpl%flush
call nicas_blk%compute_sampling_c2(mpl,rng,nam,geom,cmat_blk)

! Compute MPI distribution, halos A-B
write(mpl%info,'(a7,a)') '','Compute MPI distribution, halos A-B'
call mpl%flush
call nicas_blk%compute_mpi_ab(mpl,rng,nam,geom)

! Compute vertical interpolation data
write(mpl%info,'(a7,a)') '','Compute vertical interpolation data'
call mpl%flush
call nicas_blk%compute_interp_v(mpl,geom)

! Compute convolution data
write(mpl%info,'(a7,a)') '','Compute convolution data'
call mpl%flush
call nicas_blk%compute_convol(mpl,rng,nam,geom,cmat_blk)

! Compute MPI distribution, halo C
write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo C'
call mpl%flush
call nicas_blk%compute_mpi_c(mpl,geom)

! Compute normalization
write(mpl%info,'(a7,a)') '','Compute normalization'
call mpl%flush
call nicas_blk%compute_normalization(mpl,nam,geom)

if (nam%write_grids) then
   ! Compute grids coordinates
   write(mpl%info,'(a7,a)') '','Compute grids coordinates'
   call mpl%flush
   call nicas_blk%compute_grids(mpl,nam,geom)
end if

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0 =        ',geom%nc0
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0a =       ',nicas_blk%nc0a
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nl0 =        ',geom%nl0
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1 =        ',nicas_blk%nc1
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1a =       ',nicas_blk%nc1a
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1b =       ',nicas_blk%nc1b
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nl1 =        ',nicas_blk%nl1
do il1=1,nicas_blk%nl1
   write(mpl%info,'(a10,a,i3,a,i8)') '','nc2(',il1,') =   ',nicas_blk%nc2(il1)
   call mpl%flush
end do
write(mpl%info,'(a10,a,i8)') '','ns =         ',nicas_blk%ns
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nsa =        ',nicas_blk%nsa
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nsb =        ',nicas_blk%nsb
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nsc =        ',nicas_blk%nsc
do il0i=1,geom%nl0i
   write(mpl%info,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',nicas_blk%h(il0i)%n_s
   call mpl%flush
end do
write(mpl%info,'(a10,a,i8)') '','v%n_s =      ',nicas_blk%v%n_s
call mpl%flush
do il1=1,nicas_blk%nl1
   write(mpl%info,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',nicas_blk%s(il1)%n_s
   call mpl%flush
end do
write(mpl%info,'(a10,a,i9)') '','c%n_s =     ',nicas_blk%c%n_s
call mpl%flush
write(mpl%info,'(a10,a,i9)') '','c_nor%n_s = ',nicas_blk%c_nor%n_s
call mpl%flush

! Release memory (partial)
call nicas_blk%partial_dealloc

end subroutine nicas_blk_compute_parameters

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_sampling_c1
! Purpose: compute NICAS sampling, subset Sc1
!----------------------------------------------------------------------
subroutine nicas_blk_compute_sampling_c1(nicas_blk,mpl,rng,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(rng_type),intent(inout) :: rng              ! Random number generator
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       ! C matrix data block

! Local variables
integer :: il0,ic0,ic0a
real(kind_real) :: rhs_sum(geom%nl0),rvs_sum(geom%nl0),rvs_avg(geom%nl0),norm(geom%nl0)
real(kind_real) :: rhs_minavg
real(kind_real) :: rhs_min(geom%nc0a)
logical :: mask_hor_c0a(geom%nc0a)
character(len=1024),parameter :: subr = 'nicas_blk_compute_sampling_c1'

! Check subsampling method
if ((geom%nl0i>1).and.(trim(nam%subsamp)/='h')) call mpl%abort(subr,'subsamp = h required for variable mask')

! Allocation
allocate(nicas_blk%rhs_avg(geom%nl0))
allocate(nicas_blk%vlev(geom%nl0))

! Reset random numbers seed
if (trim(nam%strategy)=='specific_multivariate') call rng%reseed(mpl)

! Compute support radii (TODO: use draw_type)
norm = 1.0/real(geom%nc0_mask(1:geom%nl0),kind_real)
rhs_sum = sum(cmat_blk%rhs,dim=1,mask=geom%mask_c0a)
call mpl%f_comm%allreduce(rhs_sum,nicas_blk%rhs_avg,fckit_mpi_sum())
nicas_blk%rhs_avg = nicas_blk%rhs_avg*norm
rvs_sum = sum(cmat_blk%rvs,dim=1,mask=geom%mask_c0a)
call mpl%f_comm%allreduce(rvs_sum,rvs_avg,fckit_mpi_sum())
rvs_avg = rvs_avg*norm
write(mpl%info,'(a10,a)') '','Average support radii (H/V): '
call mpl%flush
do il0=1,geom%nl0
   nicas_blk%vlev(il0) = (nicas_blk%rhs_avg(il0)>0.0).or.(rvs_avg(il0)>0.0)
   if (nicas_blk%vlev(il0)) then
      write(mpl%info,'(a13,a,i3,a,f10.2,a,f10.2,a)') '','Level ',nam%levs(il0),': '//trim(mpl%aqua), &
    & nicas_blk%rhs_avg(il0)*reqkm,trim(mpl%black)//' km  / '//trim(mpl%aqua),rvs_avg(il0),trim(mpl%black)//' vert. unit'
      call mpl%flush
   else
      write(mpl%info,'(a13,a,i3,a)') '','Level ',nam%levs(il0),': missing values'
      call mpl%flush
   end if
end do
if (.not.any(nicas_blk%vlev)) call mpl%abort(subr,'no valid level')

if ((trim(nicas_blk%subsamp)=='h').or.(trim(nicas_blk%subsamp)=='hv').or.(trim(nicas_blk%subsamp)=='hvh')) then
   ! Basic horizontal mesh defined with the minimum support radius
   norm(1) = 1.0/real(geom%nc0_mask(0),kind_real)
   rhs_min = huge(1.0)
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      do il0=1,geom%nl0
         if (geom%mask_c0(ic0,il0).and.nicas_blk%vlev(il0)) then
            rhs_min(ic0a) = min(cmat_blk%rhs(ic0a,il0),rhs_min(ic0a))
         end if
      end do
   end do
   call mpl%f_comm%allreduce(sum(rhs_min,mask=geom%mask_hor_c0a),rhs_minavg,fckit_mpi_sum())
   rhs_minavg = rhs_minavg*norm(1)
   nicas_blk%nc1 = floor(2.0*maxval(geom%area)*nam%resol**2/(sqrt(3.0)*rhs_minavg**2))
   write(mpl%info,'(a10,a,i8)') '','Estimated nc1 from horizontal support radius: ',nicas_blk%nc1
   call mpl%flush
   if (nicas_blk%nc1>geom%nc0_mask(0)) then
      call mpl%warning(subr,'required nc1 larger than mask size, resetting to mask size')
      nicas_blk%nc1 = geom%nc0_mask(0)
   end if
   if (nicas_blk%nc1>nam%nc1max) then
      call mpl%warning(subr,'required nc1 larger than nc1max, resetting to nc1max')
      nicas_blk%nc1 = nam%nc1max
   end if
   if (nicas_blk%nc1<3) call mpl%abort(subr,'nicas_blk%nc1 lower than 3')
   write(mpl%info,'(a10,a,i8)') '','Final nc1: ',nicas_blk%nc1
   call mpl%flush
   write(mpl%info,'(a10,a,f5.2)') '','Effective horizontal resolution: ',sqrt(real(nicas_blk%nc1,kind_real)*sqrt(3.0) &
 & *rhs_minavg**2/(2.0*maxval(geom%area)))
   call mpl%flush
else
   ! Use the Sc0 subset
   nicas_blk%nc1 = geom%nc0_mask(0)
end if

! Allocation
allocate(nicas_blk%c1_to_c0(nicas_blk%nc1))

! Compute subset
write(mpl%info,'(a10,a)') '','Compute horizontal subset C1: '
call mpl%flush(.false.)

! Mask initialization
mask_hor_c0a = geom%mask_hor_c0a

if (nam%check_no_point_nicas) then
   ! Mask points on the last MPI task
   if (mpl%myproc==mpl%nproc) mask_hor_c0a = .false.
end if

! Compute subsampling
call initialize_sampling(mpl,rng,geom%nc0a,geom%lon_c0a,geom%lat_c0a,mask_hor_c0a,rhs_min,geom%c0a_to_c0,nam%ntry,nam%nrep, &
 & nicas_blk%nc1,nicas_blk%c1_to_c0,fast=nam%fast_sampling)
nicas_blk%c1_to_proc = geom%c0_to_proc(nicas_blk%c1_to_c0)

end subroutine nicas_blk_compute_sampling_c1

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_sampling_v
! Purpose: compute NICAS sampling, vertical dimension
!----------------------------------------------------------------------
subroutine nicas_blk_compute_sampling_v(nicas_blk,mpl,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       ! C matrix data block

! Local variables
integer :: il0,il0_prev,il1,ic0,ic1,ic0a
real(kind_real) :: distnormmin,distnorm(geom%nc0a),rv
logical :: inside
character(len=1024),parameter :: subr = 'nicas_blk_compute_sampling_v'

! Allocation
allocate(nicas_blk%slev(geom%nl0))

if ((trim(nicas_blk%subsamp)=='hv').or.(trim(nicas_blk%subsamp)=='vh').or.(trim(nicas_blk%subsamp)=='hvh')) then
   ! Vertical sampling
   write(mpl%info,'(a10,a)') '','Compute vertical subset L1'
   call mpl%flush

   ! Initialization
   nicas_blk%il0_first = mpl%msv%vali
   nicas_blk%il0_last = mpl%msv%vali
   do il0=1,geom%nl0
      if (nicas_blk%vlev(il0).and.mpl%msv%is(nicas_blk%il0_first)) nicas_blk%il0_first = il0
   end do
   do il0=geom%nl0,1,-1
      if (nicas_blk%vlev(il0).and.mpl%msv%is(nicas_blk%il0_last)) nicas_blk%il0_last = il0
   end do
   il0_prev = nicas_blk%il0_first
   nicas_blk%slev = .false.

   do il0=1,geom%nl0
      if (nicas_blk%vlev(il0)) then
         ! Look for convolution levels
         if ((il0==nicas_blk%il0_first).or.(il0==nicas_blk%il0_last)) then
            ! Keep first and last levels
            nicas_blk%slev(il0) = .true.
         else
            ! Compute minimum normalized distance with level il0_prev
            distnorm = huge(1.0)
            do ic0a=1,geom%nc0a
               if (geom%mask_c0a(ic0a,il0)) then
                  rv = sqrt(0.5*(cmat_blk%rvs(ic0a,il0)**2+cmat_blk%rvs(ic0a,il0_prev)**2))
                  if (rv>0.0) distnorm(ic0a) = abs(geom%vunit_c0a(ic0a,il0)-geom%vunit_c0a(ic0a,il0_prev))/rv
               end if
            end do
            call mpl%f_comm%allreduce(minval(distnorm),distnormmin,fckit_mpi_min())
            nicas_blk%slev(il0) = distnormmin>1.0/nam%resol
         end if

         ! Update
         if (nicas_blk%slev(il0)) il0_prev = il0
      end if
   end do
else
   ! No vertical sampling
   nicas_blk%il0_first = 1
   nicas_blk%il0_last = geom%nl0
   nicas_blk%slev = .true.
end if

! Count effective levels
nicas_blk%nl1 = count(nicas_blk%slev)
allocate(nicas_blk%l1_to_l0(nicas_blk%nl1))
write(mpl%info,'(a10,a)') '','Effective levels: '
call mpl%flush(.false.)
il1 = 0
do il0=1,geom%nl0
   if (nicas_blk%slev(il0)) then
      write(mpl%info,'(i3,a)') nam%levs(il0),' '
      call mpl%flush(.false.)
      il1 = il1+1
      nicas_blk%l1_to_l0(il1) = il0
   end if
end do
write(mpl%info,'(a)') ''
call mpl%flush

! Find bottom and top for each point of S1
allocate(nicas_blk%vbot(nicas_blk%nc1))
allocate(nicas_blk%vtop(nicas_blk%nc1))
nicas_blk%vbot = mpl%msv%vali
nicas_blk%vtop = mpl%msv%vali
!$omp parallel do schedule(static) private(ic1,ic0,inside,il1,il0)
do ic1=1,nicas_blk%nc1
   ic0 = nicas_blk%c1_to_c0(ic1)
   inside = .false.
   nicas_blk%vtop(ic1) = geom%nl0
   do il1=1,nicas_blk%nl1
      il0 = nicas_blk%l1_to_l0(il1)
      if (.not.inside.and.geom%mask_c0(ic0,il0)) then
         ! Bottom level
         nicas_blk%vbot(ic1) = il0
         inside = .true.
      end if
      if (inside.and.(.not.geom%mask_c0(ic0,il0))) then
         ! Top level
         nicas_blk%vtop(ic1) = il0
         inside = .false.
      end if
   end do
   if (mpl%msv%is(nicas_blk%vbot(ic1))) call mpl%abort(subr,'bottom level not found')
   if (mpl%msv%is(nicas_blk%vtop(ic1))) call mpl%abort(subr,'top level not found')
   if (nicas_blk%vbot(ic1)>nicas_blk%vtop(ic1)) call mpl%abort(subr,'non contiguous mask')
end do
!$omp end parallel do

! Inverse conversion
allocate(nicas_blk%l0_to_l1(geom%nl0))
nicas_blk%l0_to_l1 = mpl%msv%vali
do il1=1,nicas_blk%nl1
   il0 = nicas_blk%l1_to_l0(il1)
   nicas_blk%l0_to_l1(il0) = il1
end do

end subroutine nicas_blk_compute_sampling_v

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_sampling_c2
! Purpose: compute NICAS sampling, subset Sc2
!----------------------------------------------------------------------
subroutine nicas_blk_compute_sampling_c2(nicas_blk,mpl,rng,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(rng_type),intent(inout) :: rng              ! Random number generator
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       ! C matrix data block

! Local variables
integer :: il0,il1,ic1,ic2,is
integer,allocatable :: c2_to_c1(:)
real(kind_real),allocatable :: lon_c1a(:),lat_c1a(:),rhs_c1a(:)
logical,allocatable :: mask_c1a(:)
character(len=1024),parameter :: subr = 'nicas_blk_compute_sampling_c2'

! Allocation
allocate(nicas_blk%nc2(nicas_blk%nl1))
allocate(nicas_blk%mask_c1(nicas_blk%nc1,nicas_blk%nl1))
allocate(nicas_blk%mask_c2(nicas_blk%nc1,nicas_blk%nl1))
allocate(lon_c1a(nicas_blk%nc1a))
allocate(lat_c1a(nicas_blk%nc1a))
allocate(mask_c1a(nicas_blk%nc1a))
allocate(rhs_c1a(nicas_blk%nc1a))

! Initialization
lon_c1a = geom%lon(nicas_blk%c1_to_c0(nicas_blk%c1a_to_c1))
lat_c1a = geom%lat(nicas_blk%c1_to_c0(nicas_blk%c1a_to_c1))

! Vertically dependent horizontal subsampling
if ((trim(nicas_blk%subsamp)=='h').or.(trim(nicas_blk%subsamp)=='vh').or.(trim(nicas_blk%subsamp)=='hvh')) then
   write(mpl%info,'(a10,a)') '','Compute vertically dependent horizontal subsampling: '
   call mpl%flush
end if
do il1=1,nicas_blk%nl1
   ! Initialization
   il0 = nicas_blk%l1_to_l0(il1)
   nicas_blk%mask_c1(:,il1) = geom%mask_c0(nicas_blk%c1_to_c0,il0)

   if ((trim(nicas_blk%subsamp)=='h').or.(trim(nicas_blk%subsamp)=='vh').or.(trim(nicas_blk%subsamp)=='hvh')) then
      write(mpl%info,'(a13,a,i3,a)') '','Level ',il1,':'
      call mpl%flush

      ! Compute nc2
      nicas_blk%nc2(il1) = floor(2.0*geom%area(il0)*nam%resol**2/(sqrt(3.0)*nicas_blk%rhs_avg(il0)**2))
      write(mpl%info,'(a16,a,i8)') '','Estimated nc2 from horizontal support radius: ',nicas_blk%nc2(il1)
      call mpl%flush
      if (nicas_blk%nc2(il1)<3) call mpl%abort(subr,'nicas_blk%nc2 lower than 3')
      nicas_blk%nc2(il1) = min(nicas_blk%nc2(il1),count(nicas_blk%mask_c1(:,il1)))
      write(mpl%info,'(a16,a,i8)') '','Final nc2: ',nicas_blk%nc2(il1)
      call mpl%flush

      ! Compute horizontal subset C2
      write(mpl%info,'(a16,a)') '','Compute horizontal subset C2: '
      call mpl%flush(.false.)

      if (nicas_blk%nc2(il1)<count(nicas_blk%mask_c1(:,il1))) then
         ! Allocation
         allocate(c2_to_c1(nicas_blk%nc2(il1)))

         ! Initialize sampling
         mask_c1a = nicas_blk%mask_c1(nicas_blk%c1a_to_c1,il1)
         rhs_c1a = cmat_blk%rhs(nicas_blk%c1a_to_c0a,il0)
         call initialize_sampling(mpl,rng,nicas_blk%nc1a,lon_c1a,lat_c1a,mask_c1a,rhs_c1a,nicas_blk%c1a_to_c1, &
       & nam%ntry,nam%nrep,nicas_blk%nc2(il1),c2_to_c1,fast=nam%fast_sampling)

         ! Fill C2 mask
         nicas_blk%mask_c2(:,il1) = .false.
         do ic2=1,nicas_blk%nc2(il1)
            ic1 = c2_to_c1(ic2)
            nicas_blk%mask_c2(ic1,il1) = .true.
         end do

         ! Release memory
         deallocate(c2_to_c1)
      else
         if (mpl%main) then
            write(mpl%info,'(a)') ' use all C1 points'
            call mpl%flush
         end if

         ! Fill C2 mask
         nicas_blk%mask_c2(:,il1) = nicas_blk%mask_c1(:,il1)
      end if
   else
      ! No C2 subsampling
      nicas_blk%nc2(il1) = count(nicas_blk%mask_c1(:,il1))
      nicas_blk%mask_c2(:,il1) = nicas_blk%mask_c1(:,il1)
   end if
end do

! Release memory
deallocate(lon_c1a)
deallocate(lat_c1a)
deallocate(mask_c1a)
deallocate(rhs_c1a)

! Size
nicas_blk%ns = sum(nicas_blk%nc2)

! Allocation
allocate(nicas_blk%s_to_c1(nicas_blk%ns))
allocate(nicas_blk%s_to_l1(nicas_blk%ns))
allocate(nicas_blk%s_to_proc(nicas_blk%ns))
allocate(nicas_blk%c1l1_to_s(nicas_blk%nc1,nicas_blk%nl1))

! Subgrid conversions
is = 0
do il1=1,nicas_blk%nl1
   do ic1=1,nicas_blk%nc1
      if (nicas_blk%mask_c2(ic1,il1)) then
         is = is+1
         nicas_blk%s_to_c1(is) = ic1
         nicas_blk%s_to_l1(is) = il1
         nicas_blk%s_to_proc(is) = nicas_blk%c1_to_proc(ic1)
      end if
   end do
end do

! Conversions
nicas_blk%c1l1_to_s = mpl%msv%vali
do is=1,nicas_blk%ns
   ic1 = nicas_blk%s_to_c1(is)
   il1 = nicas_blk%s_to_l1(is)
   nicas_blk%c1l1_to_s(ic1,il1) = is
end do

end subroutine nicas_blk_compute_sampling_c2

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_mpi_a
! Purpose: compute NICAS MPI distribution, halos A
!----------------------------------------------------------------------
subroutine nicas_blk_compute_mpi_a(nicas_blk,mpl,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(geom_type),intent(in) :: geom               ! Geometry

! Local variables
integer :: ic0,ic0a,ic1,ic1a
logical :: lcheck_c1a(nicas_blk%nc1)

! Define halo A
lcheck_c1a = .false.
do ic1=1,nicas_blk%nc1
   ic0 = nicas_blk%c1_to_c0(ic1)
   if (geom%c0_to_proc(ic0)==mpl%myproc) lcheck_c1a(ic1) = .true.
end do
nicas_blk%nc1a = count(lcheck_c1a)

! Allocation
allocate(nicas_blk%c1a_to_c1(nicas_blk%nc1a))
allocate(nicas_blk%c1_to_c1a(nicas_blk%nc1))
allocate(nicas_blk%c1a_to_c0a(nicas_blk%nc1a))

! Global-local conversions for halo A
ic1a = 0
do ic1=1,nicas_blk%nc1
   if (lcheck_c1a(ic1)) then
      ic1a = ic1a+1
      nicas_blk%c1a_to_c1(ic1a) = ic1
   end if
end do
call mpl%glb_to_loc_index(nicas_blk%nc1a,nicas_blk%c1a_to_c1,nicas_blk%nc1,nicas_blk%c1_to_c1a)

! Conversion between subsets Sc0 et Sc1, halo A
do ic1a=1,nicas_blk%nc1a
   ic1 = nicas_blk%c1a_to_c1(ic1a)
   ic0 = nicas_blk%c1_to_c0(ic1)
   ic0a = geom%c0_to_c0a(ic0)
   nicas_blk%c1a_to_c0a(ic1a) = ic0a
end do

end subroutine nicas_blk_compute_mpi_a

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_mpi_ab
! Purpose: compute NICAS MPI distribution, halos A-B
!----------------------------------------------------------------------
subroutine nicas_blk_compute_mpi_ab(nicas_blk,mpl,rng,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(rng_type),intent(inout) :: rng              ! Random number generator
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry

! Local variables
integer :: il0i,ic0,ic1,jc1,ic1b_h,ic1b,il0,il1,isa,isb,i_s,is,js,nc1b_h
integer :: s_to_proc(nicas_blk%ns)
integer,allocatable :: c1b_h_to_c1(:),c1b_h_to_c1b(:),sa_to_sb(:),s_to_sb(:)
real(kind_real) :: lon_c1(nicas_blk%nc1),lat_c1(nicas_blk%nc1)
real(kind_real),allocatable :: lon_c1b_h(:),lat_c1b_h(:)
logical :: mask_c1(nicas_blk%nc1),lcheck_c1b_h(nicas_blk%nc1),lcheck_c1b(nicas_blk%nc1)
logical,allocatable :: mask_c1b_h(:)
character(len=1024),parameter :: subr = 'nicas_blk_compute_mpi_ab'

! Allocation
allocate(nicas_blk%h(geom%nl0i))
allocate(nicas_blk%s(nicas_blk%nl1))
allocate(nicas_blk%lcheck_sa(nicas_blk%ns))
allocate(nicas_blk%lcheck_sb(nicas_blk%ns))

! Compute interpolation
lon_c1 = geom%lon(nicas_blk%c1_to_c0)
lat_c1 = geom%lat(nicas_blk%c1_to_c0)
do il0i=1,geom%nl0i
   mask_c1 = geom%mask_c0(nicas_blk%c1_to_c0,il0i)
   write(nicas_blk%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
   call nicas_blk%h(il0i)%interp(mpl,rng,nam,geom,il0i,nicas_blk%nc1,lon_c1,lat_c1,mask_c1,geom%nc0a, &
    & geom%lon_c0a,geom%lat_c0a,geom%mask_c0a(:,il0i))
end do

! Define halo A
nicas_blk%lcheck_sa = .false.
do is=1,nicas_blk%ns
   ic1 = nicas_blk%s_to_c1(is)
   ic0 = nicas_blk%c1_to_c0(ic1)
   if (geom%c0_to_proc(ic0)==mpl%myproc) nicas_blk%lcheck_sa(is) = .true.
end do
nicas_blk%nsa = count(nicas_blk%lcheck_sa)

! Define halo B (after first horizontal interpolation)
lcheck_c1b_h = .false.
do il0i=1,geom%nl0i
   do i_s=1,nicas_blk%h(il0i)%n_s
      jc1 = nicas_blk%h(il0i)%col(i_s)
      lcheck_c1b_h(jc1) = .true.
   end do
end do
nc1b_h = count(lcheck_c1b_h)

! Allocation
allocate(c1b_h_to_c1(nc1b_h))

! Global-local conversions for halo B (after first horizontal interpolation)
ic1b_h = 0
do ic1=1,nicas_blk%nc1
   if (lcheck_c1b_h(ic1)) then
      ic1b_h = ic1b_h+1
      c1b_h_to_c1(ic1b_h) = ic1
   end if
end do

! Compute interpolation
do il1=1,nicas_blk%nl1
   write(nicas_blk%s(il1)%prefix,'(a,i3.3)') 's_',il1

   ! Allocation
   allocate(lon_c1b_h(nc1b_h))
   allocate(lat_c1b_h(nc1b_h))
   allocate(mask_c1b_h(nc1b_h))

   ! Initialization
   lon_c1b_h = lon_c1(c1b_h_to_c1)
   lat_c1b_h = lat_c1(c1b_h_to_c1)
   mask_c1b_h = nicas_blk%mask_c1(c1b_h_to_c1,il1)

   ! Compute interpolation
   il0 = nicas_blk%l1_to_l0(il1)
   call nicas_blk%s(il1)%interp(mpl,rng,nam,geom,il0,nicas_blk%nc1,lon_c1,lat_c1,nicas_blk%mask_c2(:,il1),nc1b_h, &
 & lon_c1b_h,lat_c1b_h,mask_c1b_h)

   ! Release memory
   deallocate(lon_c1b_h)
   deallocate(lat_c1b_h)
   deallocate(mask_c1b_h)
end do

! Define halo B (required for the second horizontal interpolation)
nicas_blk%lcheck_sb = nicas_blk%lcheck_sa
lcheck_c1b = lcheck_c1b_h
do il1=1,nicas_blk%nl1
   do i_s=1,nicas_blk%s(il1)%n_s
      jc1 = nicas_blk%s(il1)%col(i_s)
      js = nicas_blk%c1l1_to_s(jc1,il1)
      lcheck_c1b(jc1) = .true.
      nicas_blk%lcheck_sb(js) = .true.
   end do
end do
nicas_blk%nc1b = count(lcheck_c1b)
nicas_blk%nsb = count(nicas_blk%lcheck_sb)

! Check halos consistency
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sa(is).and.(.not.nicas_blk%lcheck_sb(is))) call mpl%abort(subr,'point in halo A but not in halo B')
end do

! Allocation
allocate(nicas_blk%sa_to_s(nicas_blk%nsa))
allocate(nicas_blk%s_to_sa(nicas_blk%ns))
allocate(nicas_blk%c1b_to_c1(nicas_blk%nc1b))
allocate(nicas_blk%c1_to_c1b(nicas_blk%nc1))
allocate(nicas_blk%sb_to_s(nicas_blk%nsb))
allocate(s_to_sb(nicas_blk%ns))
allocate(c1b_h_to_c1b(nc1b_h))
allocate(sa_to_sb(nicas_blk%nsa))
allocate(nicas_blk%sb_to_c1b(nicas_blk%nsb))
allocate(nicas_blk%sb_to_l1(nicas_blk%nsb))
allocate(nicas_blk%c1bl1_to_sb(nicas_blk%nc1b,nicas_blk%nl1))

! Global-local conversions for halo A
isa = 0
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sa(is)) then
      isa = isa+1
      nicas_blk%sa_to_s(isa) = is
   end if
end do
call mpl%glb_to_loc_index(nicas_blk%nsa,nicas_blk%sa_to_s,nicas_blk%ns,nicas_blk%s_to_sa)

! Global-local conversions for halo B
nicas_blk%c1_to_c1b = mpl%msv%vali
ic1b = 0
do ic1=1,nicas_blk%nc1
   if (lcheck_c1b(ic1)) then
      ic1b = ic1b+1
      nicas_blk%c1b_to_c1(ic1b) = ic1
      nicas_blk%c1_to_c1b(ic1) = ic1b
   end if
end do
s_to_sb = mpl%msv%vali
isb = 0
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sb(is)) then
      isb = isb+1
      nicas_blk%sb_to_s(isb) = is
      s_to_sb(is) = isb
   end if
end do

! Halos A-B conversions
do isa=1,nicas_blk%nsa
   is = nicas_blk%sa_to_s(isa)
   isb = s_to_sb(is)
   sa_to_sb(isa) = isb
end do

! Halos B_h-B conversions
do ic1b_h=1,nc1b_h
   ic1 = c1b_h_to_c1(ic1b_h)
   ic1b = nicas_blk%c1_to_c1b(ic1)
   c1b_h_to_c1b(ic1b_h) = ic1b
end do

! Local interpolation source
do il0i=1,geom%nl0i
   nicas_blk%h(il0i)%n_src = nicas_blk%nc1b
   do i_s=1,nicas_blk%h(il0i)%n_s
      nicas_blk%h(il0i)%col(i_s) = nicas_blk%c1_to_c1b(nicas_blk%h(il0i)%col(i_s))
   end do
   call nicas_blk%h(il0i)%reorder(mpl)
end do
do il1=1,nicas_blk%nl1
   nicas_blk%s(il1)%n_src = nicas_blk%nc1b
   nicas_blk%s(il1)%n_dst = nicas_blk%nc1b
   do i_s=1,nicas_blk%s(il1)%n_s
      nicas_blk%s(il1)%row(i_s) = c1b_h_to_c1b(nicas_blk%s(il1)%row(i_s))
      nicas_blk%s(il1)%col(i_s) = nicas_blk%c1_to_c1b(nicas_blk%s(il1)%col(i_s))
   end do
   call nicas_blk%s(il1)%reorder(mpl)
end do

! Conversion from subset Sc1 to subgrid
nicas_blk%c1bl1_to_sb = mpl%msv%vali
do isb=1,nicas_blk%nsb
   is = nicas_blk%sb_to_s(isb)
   il1 = nicas_blk%s_to_l1(is)
   ic1 = nicas_blk%s_to_c1(is)
   ic1b = nicas_blk%c1_to_c1b(ic1)
   nicas_blk%sb_to_c1b(isb) = ic1b
   nicas_blk%sb_to_l1(isb) = il1
   nicas_blk%c1bl1_to_sb(ic1b,il1) = isb
end do

! MPI splitting on subgrid
do is=1,nicas_blk%ns
   ic1 = nicas_blk%s_to_c1(is)
   ic0 = nicas_blk%c1_to_c0(ic1)
   s_to_proc(is) = geom%c0_to_proc(ic0)
end do

! Setup communications
call nicas_blk%com_AB%setup(mpl,'com_AB',nicas_blk%ns,nicas_blk%nsa,nicas_blk%nsb,nicas_blk%nsa,nicas_blk%sb_to_s,sa_to_sb, &
 & s_to_proc,nicas_blk%s_to_sa)

! Release memory
deallocate(c1b_h_to_c1)
deallocate(c1b_h_to_c1b)
deallocate(sa_to_sb)
deallocate(s_to_sb)

end subroutine nicas_blk_compute_mpi_ab

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_interp_v
! Purpose: compute vertical interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_compute_interp_v(nicas_blk,mpl,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(geom_type),intent(in) :: geom               ! Geometry

! Local variables
integer :: ic1b,ic1,ic0,jl0,il0,jl1,il0inf,il0sup

! Initialize vertical interpolation
nicas_blk%v%prefix = 'v'
nicas_blk%v%n_src = nicas_blk%nl1
nicas_blk%v%n_dst = geom%nl0

! Count levels
nicas_blk%v%n_s = nicas_blk%nl1
il0inf = nicas_blk%il0_first
do jl0=1,geom%nl0
   if (nicas_blk%slev(jl0)) then
      il0sup = jl0
      do il0=il0inf+1,il0sup-1
         nicas_blk%v%n_s = nicas_blk%v%n_s+2
      end do
      il0inf = jl0
   end if
end do

! Allocation
call nicas_blk%v%alloc(nicas_blk%nc1b)

! Set identity for subsampled levels
do jl1=1,nicas_blk%nl1
   jl0 = nicas_blk%l1_to_l0(jl1)
   nicas_blk%v%row(jl1) = jl0
   nicas_blk%v%col(jl1) = jl0
   do ic1b=1,nicas_blk%nc1b
      nicas_blk%v%Svec(jl1,ic1b) = 1.0
   end do
end do
nicas_blk%v%n_s = nicas_blk%nl1
nicas_blk%v%n_s = nicas_blk%nl1

! Compute linear interpolation for other levels
il0inf = nicas_blk%il0_first
do jl0=1,geom%nl0
   if (nicas_blk%slev(jl0)) then
      il0sup = jl0
      do il0=il0inf+1,il0sup-1
         nicas_blk%v%n_s = nicas_blk%v%n_s+1
         nicas_blk%v%row(nicas_blk%v%n_s) = il0
         nicas_blk%v%col(nicas_blk%v%n_s) = il0inf
         do ic1b=1,nicas_blk%nc1b
            ic1 = nicas_blk%c1b_to_c1(ic1b)
            ic0 = nicas_blk%c1_to_c0(ic1)
            nicas_blk%v%Svec(nicas_blk%v%n_s,ic1b) = abs(geom%vunit_c0(ic0,il0sup)-geom%vunit_c0(ic0,il0)) &
                                                   & /abs(geom%vunit_c0(ic0,il0sup)-geom%vunit_c0(ic0,il0inf))
         end do
         nicas_blk%v%n_s = nicas_blk%v%n_s+1
         nicas_blk%v%row(nicas_blk%v%n_s) = il0
         nicas_blk%v%col(nicas_blk%v%n_s) = il0sup
         do ic1b=1,nicas_blk%nc1b
            ic1 = nicas_blk%c1b_to_c1(ic1b)
            ic0 = nicas_blk%c1_to_c0(ic1)
            nicas_blk%v%Svec(nicas_blk%v%n_s,ic1b) = abs(geom%vunit_c0(ic0,il0)-geom%vunit_c0(ic0,il0inf)) &
                                                   & /abs(geom%vunit_c0(ic0,il0sup)-geom%vunit_c0(ic0,il0inf))
         end do
      end do
      il0inf = jl0
   end if
end do

! Conversion
nicas_blk%v%col = nicas_blk%l0_to_l1(nicas_blk%v%col)

! Reorder
call nicas_blk%v%reorder(mpl)

! Release memory
deallocate(nicas_blk%slev)

end subroutine nicas_blk_compute_interp_v

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol
! Purpose: compute convolution
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol(nicas_blk,mpl,rng,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(rng_type),intent(inout) :: rng              ! Random number generator
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       ! C matrix data block

! Local variables
integer :: n_s_max,ithread,is,ic1,jc1,il1,il0,j,js,isb,ic1b,ic0,ic0a,ic1a,i_s,jc,kc,ks,jbd,jc0,jl0,jl1,ic1bb,isbb
integer :: c_n_s(mpl%nthread)
integer,allocatable :: nn(:),nn_index(:),inec(:),c_ind(:,:)
real(kind_real) :: distvsq,rvsq
real(kind_real) :: lon_c1(nicas_blk%nc1),lat_c1(nicas_blk%nc1)
real(kind_real),allocatable :: rh_c1a(:,:),rv_c1a(:,:),rv_rfac_c1a(:,:),rv_rfac_c1(:,:),rv_coef_c1a(:,:),rv_coef_c1(:,:)
real(kind_real),allocatable :: H11_c1a(:,:),H22_c1a(:,:),H33_c1a(:,:),H12_c1a(:,:),Hcoef_c1a(:,:),Hcoef_c1(:,:)
real(kind_real),allocatable :: distnormv(:,:),rfac(:,:),coef(:,:),Hcoef(:,:)
real(kind_real),allocatable :: c_S(:,:),c_S_conv(:)
logical :: add_op,lcheck_c1bb(nicas_blk%nc1)
type(linop_type) :: ctmp,c(mpl%nthread)

! Associate
associate(ib=>nicas_blk%ib)

! Set double-fit parameter
nicas_blk%double_fit = cmat_blk%double_fit
nicas_blk%anisotropic = cmat_blk%anisotropic

! Allocation
call nicas_blk%tree%alloc(mpl,nicas_blk%nc1)

! Initialization
lon_c1 = geom%lon(nicas_blk%c1_to_c0)
lat_c1 = geom%lat(nicas_blk%c1_to_c0)
call nicas_blk%tree%init(lon_c1,lat_c1)

! Find largest possible radius
call mpl%f_comm%allreduce(maxval(cmat_blk%rh,mask=mpl%msv%isnot(cmat_blk%rh)),nicas_blk%rhmax,fckit_mpi_max())
if (nicas_blk%double_fit) then
   nicas_blk%rhmax = nicas_blk%rhmax*sqrt_r_dble
else
   nicas_blk%rhmax = nicas_blk%rhmax*sqrt_r
end if

if (nam%lsqrt) then
   ! Copy
   nicas_blk%nc1bb = nicas_blk%nc1b
else
   ! Allocation
   allocate(nn(nicas_blk%nc1b))

   do ic1b=1,nicas_blk%nc1b
      ! Indices
      ic1 = nicas_blk%c1b_to_c1(ic1b)
      ic0 = nicas_blk%c1_to_c0(ic1)

      ! Count nearest neighbors
      call nicas_blk%tree%count_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nicas_blk%rhmax,nn(ic1b))
   end do

   ! Initialization
   write(mpl%info,'(a10,a)') '','Define extended halo: '
   call mpl%flush(.false.)
   call mpl%prog_init(nicas_blk%nc1b)
   lcheck_c1bb = .false.

   do ic1b=1,nicas_blk%nc1b
      ! Indices
      ic1 = nicas_blk%c1b_to_c1(ic1b)
      ic0 = nicas_blk%c1_to_c0(ic1)

      ! Allocation
      allocate(nn_index(nn(ic1b)))

      ! Find nearest neighbors
      call nicas_blk%tree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nn(ic1b),nn_index)

      ! Fill mask
      do j=1,nn(ic1b)
         jc1 = nn_index(j)
         lcheck_c1bb(jc1) = .true.
      end do

      ! Release memory
      deallocate(nn_index)

      ! Update
      call mpl%prog_print(ic1b)
   end do
   call mpl%prog_final

   ! Halo size
   nicas_blk%nc1bb = count(lcheck_c1bb)
   write(mpl%info,'(a10,a,i6,a,i6)') '','Halo sizes nc1b / nc1bb: ',nicas_blk%nc1b,' / ',nicas_blk%nc1bb
   call mpl%flush

   ! Release memory
   deallocate(nn)
end if

! Allocation
allocate(nicas_blk%c1bb_to_c1(nicas_blk%nc1bb))
allocate(nicas_blk%c1_to_c1bb(nicas_blk%nc1))

if (nam%lsqrt) then
   ! Copy
   nicas_blk%c1bb_to_c1 = nicas_blk%c1b_to_c1
   nicas_blk%c1_to_c1bb = nicas_blk%c1_to_c1b
   nicas_blk%nsbb = nicas_blk%nsb
else
   ! Global <-> local conversions for fields
   nicas_blk%c1_to_c1bb = mpl%msv%vali
   ic1bb = 0
   do ic1=1,nicas_blk%nc1
      if (lcheck_c1bb(ic1)) then
         ic1bb = ic1bb+1
         nicas_blk%c1bb_to_c1(ic1bb) = ic1
         nicas_blk%c1_to_c1bb(ic1) = ic1bb
      end if
   end do

   ! Count points in extended halo
   nicas_blk%nsbb = 0
   do is=1,nicas_blk%ns
      ic1 = nicas_blk%s_to_c1(is)
      if (lcheck_c1bb(ic1)) nicas_blk%nsbb = nicas_blk%nsbb+1
   end do
   write(mpl%info,'(a10,a,i6,a,i6)') '','Halo sizes nsb / nsbb:   ',nicas_blk%nsb,' / ',nicas_blk%nsbb
   call mpl%flush
end if

! Allocation
allocate(nicas_blk%sbb_to_s(nicas_blk%nsbb))

if (nam%lsqrt) then
   ! Copy
   nicas_blk%sbb_to_s = nicas_blk%sb_to_s
else
   ! Global <-> local conversions for fields
   isbb = 0
   do is=1,nicas_blk%ns
      ic1 = nicas_blk%s_to_c1(is)
      if (lcheck_c1bb(ic1)) then
         isbb = isbb+1
         nicas_blk%sbb_to_s(isbb) = is
      end if
   end do
end if

! Compute horizontal and vertical parameters
write(mpl%info,'(a10,a)') '','Compute horizontal and vertical parameters'
call mpl%flush

! Allocation
allocate(rh_c1a(nicas_blk%nc1a,nicas_blk%nl1))
allocate(rv_c1a(nicas_blk%nc1a,nicas_blk%nl1))
allocate(nicas_blk%rh_c1(nicas_blk%nc1,nicas_blk%nl1))
allocate(nicas_blk%rv_c1(nicas_blk%nc1,nicas_blk%nl1))
if (nicas_blk%double_fit) then
   allocate(rv_rfac_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(rv_coef_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(rv_rfac_c1(nicas_blk%nc1,nicas_blk%nl1))
   allocate(rv_coef_c1(nicas_blk%nc1,nicas_blk%nl1))
   allocate(nicas_blk%rfac(nicas_blk%nsbb))
   allocate(nicas_blk%coef(nicas_blk%nsbb))
   allocate(nicas_blk%distnormv(nicas_blk%nsbb))
end if
if (nicas_blk%anisotropic) then
   allocate(H11_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(H22_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(H33_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(H12_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(Hcoef_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(nicas_blk%H11_c1(nicas_blk%nc1,nicas_blk%nl1))
   allocate(nicas_blk%H22_c1(nicas_blk%nc1,nicas_blk%nl1))
   allocate(nicas_blk%H33_c1(nicas_blk%nc1,nicas_blk%nl1))
   allocate(nicas_blk%H12_c1(nicas_blk%nc1,nicas_blk%nl1))
   allocate(Hcoef_c1(nicas_blk%nc1,nicas_blk%nl1))
   allocate(nicas_blk%Hcoef(nicas_blk%nsbb))
end if

! Copy and rescale
write(mpl%info,'(a13,a)') '','Copy and rescale'
call mpl%flush
do il1=1,nicas_blk%nl1
   do ic1a=1,nicas_blk%nc1a
      ! Indices
      ic0a = nicas_blk%c1a_to_c0a(ic1a)
      il0 = nicas_blk%l1_to_l0(il1)

      if (geom%mask_c0a(ic0a,il0)) then
         ! Copy
         rh_c1a(ic1a,il1) = cmat_blk%rh(ic0a,il0)
         rv_c1a(ic1a,il1) = cmat_blk%rv(ic0a,il0)
         if (nicas_blk%double_fit) then
            rv_rfac_c1a(ic1a,il1) = cmat_blk%rv_rfac(ic0a,il0)
            rv_coef_c1a(ic1a,il1) = cmat_blk%rv_coef(ic0a,il0)
         end if
         if (nicas_blk%anisotropic) then
            H11_c1a(ic1a,il1) = cmat_blk%H11(ic0a,il0)
            H22_c1a(ic1a,il1) = cmat_blk%H22(ic0a,il0)
            H33_c1a(ic1a,il1) = cmat_blk%H33(ic0a,il0)
            H12_c1a(ic1a,il1) = cmat_blk%H12(ic0a,il0)
            Hcoef_c1a(ic1a,il1) = cmat_blk%Hcoef(ic0a,il0)
         end if

         ! Square-root rescaling
         rh_c1a(ic1a,il1) = rh_c1a(ic1a,il1)*sqrt_r
         if (nicas_blk%double_fit) then
            rv_c1a(ic1a,il1) = rv_c1a(ic1a,il1)*sqrt_r_dble
            rv_rfac_c1a(ic1a,il1) = rv_rfac_c1a(ic1a,il1)*sqrt_rfac
            rv_coef_c1a(ic1a,il1) = rv_coef_c1a(ic1a,il1)*sqrt_coef
         else
            rv_c1a(ic1a,il1) = rv_c1a(ic1a,il1)*sqrt_r
         end if
         if (nicas_blk%anisotropic) then
            H11_c1a(ic1a,il1) = H11_c1a(ic1a,il1)/sqrt_r**2
            H22_c1a(ic1a,il1) = H22_c1a(ic1a,il1)/sqrt_r**2
            H33_c1a(ic1a,il1) = H33_c1a(ic1a,il1)/sqrt_r**2
            H12_c1a(ic1a,il1) = H12_c1a(ic1a,il1)/sqrt_r**2
         end if
      else
         ! Missing values
         rh_c1a(ic1a,il1) = mpl%msv%valr
         rv_c1a(ic1a,il1) = mpl%msv%valr
         if (nicas_blk%double_fit) then
            rv_rfac_c1a(ic1a,il1) = mpl%msv%valr
            rv_coef_c1a(ic1a,il1) = mpl%msv%valr
         end if
         if (nicas_blk%anisotropic) then
            H11_c1a(ic1a,il1) = mpl%msv%valr
            H22_c1a(ic1a,il1) = mpl%msv%valr
            H33_c1a(ic1a,il1) = mpl%msv%valr
            H12_c1a(ic1a,il1) = mpl%msv%valr
            Hcoef_c1a(ic1a,il1) = mpl%msv%valr
         end if
      end if
   end do
end do

! Communication
write(mpl%info,'(a13,a)') '','Communication'
call mpl%flush
call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,rh_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & nicas_blk%rh_c1)
call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,rv_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & nicas_blk%rv_c1)
if (nicas_blk%double_fit) then
   call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,rv_rfac_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a, &
 & .true.,rv_rfac_c1)
   call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,rv_coef_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a, &
 & .true.,rv_coef_c1)
end if
if (nicas_blk%anisotropic) then
   call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,H11_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & nicas_blk%H11_c1)
   call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,H22_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & nicas_blk%H22_c1)
   call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,H33_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & nicas_blk%H33_c1)
   call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,H12_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & nicas_blk%H12_c1)
   call mpl%loc_to_glb(nicas_blk%nl1,nicas_blk%nc1a,Hcoef_c1a,nicas_blk%nc1,nicas_blk%c1_to_proc,nicas_blk%c1_to_c1a,.true., &
 & Hcoef_c1)
end if

! Release memory
deallocate(rh_c1a)
deallocate(rv_c1a)
if (nicas_blk%double_fit) then
   deallocate(rv_rfac_c1a)
   deallocate(rv_coef_c1a)
end if
if (nicas_blk%anisotropic) then
   deallocate(H11_c1a)
   deallocate(H22_c1a)
   deallocate(H33_c1a)
   deallocate(H12_c1a)
   deallocate(Hcoef_c1a)
end if

! Allocation
allocate(nicas_blk%distnorm(nicas_blk%nsbb))

! Compute distances
if (nam%network) then
   call nicas_blk%compute_convol_network(mpl,rng,nam,geom)
else
   call nicas_blk%compute_convol_distance(mpl,nam,geom)
end if

! Release memory
deallocate(nicas_blk%rh_c1)
if (nicas_blk%anisotropic) then
   deallocate(nicas_blk%H11_c1)
   deallocate(nicas_blk%H22_c1)
   deallocate(nicas_blk%H33_c1)
   deallocate(nicas_blk%H12_c1)
end if

! Compute ball data
write(mpl%info,'(a13,a)') '','Compute ball data'
call mpl%flush
!$omp parallel do schedule(static) private(isbb,is,ic1,il1,ic0,il0,jbd,jc1,jl1,jc0,jl0,distvsq,rvsq), &
!$omp&                             firstprivate(distnormv,rfac,coef,Hcoef)
do isbb=1,nicas_blk%nsbb
   ! Indices
   is = nicas_blk%sbb_to_s(isbb)
   ic1 = nicas_blk%s_to_c1(is)
   il1 = nicas_blk%s_to_l1(is)
   ic0 = nicas_blk%c1_to_c0(ic1)
   il0 = nicas_blk%l1_to_l0(il1)

   ! Allocation
   if (nicas_blk%double_fit) then
      allocate(distnormv(nicas_blk%nc1,nicas_blk%nl1))
      allocate(rfac(nicas_blk%nc1,nicas_blk%nl1))
      allocate(coef(nicas_blk%nc1,nicas_blk%nl1))
   end if
   if (nicas_blk%anisotropic) allocate(Hcoef(nicas_blk%nc1,nicas_blk%nl1))

   ! Initialization
   if (nicas_blk%double_fit) then
      distnormv = mpl%msv%valr
      rfac = mpl%msv%valr
      coef = mpl%msv%valr
   end if
   if (nicas_blk%anisotropic) Hcoef = mpl%msv%valr

   do jbd=1,nicas_blk%distnorm(isbb)%nbd
      ! Indices
      jc1 = nicas_blk%distnorm(isbb)%bd_to_c1(jbd)
      jl1 = nicas_blk%distnorm(isbb)%bd_to_l1(jbd)
      jc0 = nicas_blk%c1_to_c0(jc1)
      jl0 = nicas_blk%l1_to_l0(jl1)

      if (nicas_blk%double_fit) then
         ! Vertical distance
         distvsq = (geom%vunit_c0(ic0,il0)-geom%vunit_c0(jc0,jl0))**2
         rvsq = 0.5*(nicas_blk%rv_c1(ic1,il1)**2+nicas_blk%rv_c1(jc1,jl1)**2)
         if (rvsq>0.0) then
            distnormv(jc1,jl1) = sqrt(distvsq/rvsq)
         elseif (distvsq>0.0) then
            distnormv(jc1,jl1) = 0.5*huge(0.0)
         end if
         rfac(jc1,jl1) = sqrt(rv_rfac_c1(ic1,il1)*rv_rfac_c1(jc1,jl1))
         coef(jc1,jl1) = sqrt(rv_coef_c1(ic1,il1)*rv_coef_c1(jc1,jl1))
      end if
      if (nicas_blk%anisotropic) Hcoef(jc1,jl1) = sqrt(Hcoef_c1(ic1,il1)*Hcoef_c1(jc1,jl1))
   end do

   ! Pack data
   if (nicas_blk%double_fit) then
      call nicas_blk%distnormv(isbb)%pack(mpl,nicas_blk%nc1,nicas_blk%nl1,distnormv)
      call nicas_blk%rfac(isbb)%pack(mpl,nicas_blk%nc1,nicas_blk%nl1,rfac)
      call nicas_blk%coef(isbb)%pack(mpl,nicas_blk%nc1,nicas_blk%nl1,coef)
   end if
   if (nicas_blk%anisotropic) call nicas_blk%Hcoef(isbb)%pack(mpl,nicas_blk%nc1,nicas_blk%nl1,Hcoef)

   ! Release memory
   if (nicas_blk%double_fit) then
      deallocate(distnormv)
      deallocate(rfac)
      deallocate(coef)
   end if
   if (nicas_blk%anisotropic) deallocate(Hcoef)
end do
!$omp end parallel do

! Release memory
deallocate(nicas_blk%rv_c1)
if (nicas_blk%double_fit) then
   deallocate(rv_rfac_c1)
   deallocate(rv_coef_c1)
end if
if (nicas_blk%anisotropic) deallocate(Hcoef_c1)

! Compute weights
call nicas_blk%compute_convol_weights(mpl,nam,geom,ctmp)

! Release memory
do isbb=1,nicas_blk%nsbb
   call nicas_blk%distnorm(isbb)%dealloc
   if (nicas_blk%double_fit) then
      call nicas_blk%distnormv(isbb)%dealloc
      call nicas_blk%rfac(isbb)%dealloc
      call nicas_blk%coef(isbb)%dealloc
   end if
   if (nicas_blk%anisotropic) call nicas_blk%Hcoef(isbb)%dealloc
end do
deallocate(nicas_blk%distnorm)
if (nicas_blk%double_fit) then
   deallocate(nicas_blk%distnormv)
   deallocate(nicas_blk%rfac)
   deallocate(nicas_blk%coef)
end if
if (nicas_blk%anisotropic) deallocate(nicas_blk%Hcoef)

if (nam%lsqrt) then
   ! Copy
   call nicas_blk%c%copy(ctmp)
else
   ! Compute convolution inverse mapping
   allocate(inec(nicas_blk%ns))
   inec = 0
   do i_s=1,ctmp%n_s
      is = ctmp%col(i_s)
      inec(is) = inec(is)+1
   end do

   allocate(c_ind(maxval(inec),nicas_blk%ns))
   allocate(c_S(maxval(inec),nicas_blk%ns))
   c_ind = mpl%msv%vali
   c_S = mpl%msv%valr
   inec = 0
   do i_s=1,ctmp%n_s
      is = ctmp%col(i_s)
      js = ctmp%row(i_s)
      inec(is) = inec(is)+1
      c_ind(inec(is),is) = js
      c_S(inec(is),is) = ctmp%S(i_s)
   end do

   ! Initialization
   write(mpl%info,'(a10,a)') '','Second pass:     '
   call mpl%flush(.false.)
   call mpl%prog_init(nicas_blk%nsb)
   n_s_max = 100*nint(real(geom%nc0*geom%nl0,kind_real)/real(mpl%nthread*mpl%nproc,kind_real))
   c_n_s = 0
   do ithread=1,mpl%nthread
      c(ithread)%n_s = n_s_max
      call c(ithread)%alloc
      c(ithread)%row = mpl%msv%vali
      c(ithread)%col = mpl%msv%vali
      c(ithread)%S = mpl%msv%valr
   end do

   ! Apply convolution
   !$omp parallel do schedule(static) private(isb,is,ithread,jc,js,kc,ks,add_op) firstprivate(c_S_conv)
   do isb=1,nicas_blk%nsb
      ! Indices
      is = nicas_blk%sb_to_s(isb)
      ithread = 1
!$    ithread = omp_get_thread_num()+1

      ! Allocation
      allocate(c_S_conv(nicas_blk%ns))

      ! Initialization
      c_S_conv = 0.0

      ! Loop twice over points
      do jc=1,inec(is)
         js = c_ind(jc,is)
         do kc=1,inec(js)
            ks = c_ind(kc,js)
            c_S_conv(ks) = c_S_conv(ks)+c_S(jc,is)*c_S(kc,js)
         end do
      end do

      ! Store coefficient for convolution
      do js=1,nicas_blk%ns
         add_op = .false.
         if (nam%mpicom==1) then
            add_op = (nicas_blk%lcheck_sb(js).and.(is<=js)).or.(.not.nicas_blk%lcheck_sb(js))
         elseif (nam%mpicom==2) then
            add_op = nicas_blk%lcheck_sa(is).and.((nicas_blk%lcheck_sa(js).and.(is<=js)) &
                   & .or.(.not.nicas_blk%lcheck_sa(js)))
         end if
         if (add_op) call c(ithread)%add_op(c_n_s(ithread),is,js,c_S_conv(js))
      end do

      ! Release memory
      deallocate(c_S_conv)

      ! Update
      call mpl%prog_print(isb)
   end do
   !$omp end parallel do
   call mpl%prog_final

   ! Gather data from threads
   call nicas_blk%c%gather(mpl,c_n_s,c)

   ! Release memory
   deallocate(inec)
   deallocate(c_ind)
   deallocate(c_S)
   do ithread=1,mpl%nthread
      call c(ithread)%dealloc
   end do
end if

! Set prefix
nicas_blk%c%prefix = 'c'
nicas_blk%c_nor%prefix = 'c_nor'

! End associate
end associate

end subroutine nicas_blk_compute_convol

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol_network
! Purpose: compute convolution with a network approach
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol_network(nicas_blk,mpl,rng,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(rng_type),intent(inout) :: rng              ! Random number generator
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry

! Local variables
integer :: net_nnbmax,is,ic0,ic1,il0,jl1,np,np_new,i,j,k,ip,kc1,jc1,il1,dkl1,kl1,jp,isbb,djl1,inr,jc0,jl0
integer :: net_nnb(nicas_blk%nc1),ic1_loc,nc1_loc(0:mpl%nproc)
integer,allocatable :: net_inb(:,:),plist(:,:),plist_new(:,:)
real(kind_real) :: distnorm_network,disttest
real(kind_real) :: dnb,dx,dy,dz,disthsq,distvsq,rhsq,rvsq,H11,H22,H33,H12
real(kind_real),allocatable :: lon_c1(:),lat_c1(:)
real(kind_real),allocatable :: distnorm(:,:),net_dnb(:,:,:,:)
logical :: init,add_to_front
logical,allocatable :: valid_arc(:,:,:)
type(mesh_type) :: mesh

! Allocation
allocate(lon_c1(nicas_blk%nc1))
allocate(lat_c1(nicas_blk%nc1))
call mesh%alloc(nicas_blk%nc1)

! Initialization
lon_c1 = geom%lon(nicas_blk%c1_to_c0)
lat_c1 = geom%lat(nicas_blk%c1_to_c0)
call mesh%init(mpl,rng,lon_c1,lat_c1,.true.)

! Count neighbors
write(mpl%info,'(a10,a)') '','Count neighbors'
call mpl%flush
net_nnb = 0
do ic1=1,nicas_blk%nc1
   ! Mesh center
   inr = mesh%order_inv(ic1)
   i = mesh%lend(inr)
   init = .true.

   ! Loop over neighbors
   do while ((i/=mesh%lend(inr)).or.init)
      net_nnb(ic1) = net_nnb(ic1)+1
      i = mesh%lptr(i)
      init = .false.
   end do
end do

! Allocation
net_nnbmax = maxval(net_nnb)
allocate(net_inb(nicas_blk%nc1,net_nnbmax))
allocate(net_dnb(nicas_blk%nc1,nicas_blk%nl1,net_nnbmax,3))
allocate(valid_arc(nicas_blk%nc1,nicas_blk%nl1,net_nnbmax))

! MPI splitting
call mpl%split_loop(nicas_blk%nc1,nc1_loc)

! Find mesh neighbors
write(mpl%info,'(a10,a)') '','Find mesh neighbors: '
call mpl%flush(.false.)
call mpl%prog_init(nc1_loc(mpl%myproc))
net_nnb = 0
net_inb = mpl%msv%vali
valid_arc = .false.
do ic1_loc=1,nc1_loc(mpl%myproc)
   ! Index
   ic1 = sum(nc1_loc(0:mpl%myproc-1))+ic1_loc

   ! Mesh center
   inr = mesh%order_inv(ic1)
   i = mesh%lend(inr)
   init = .true.

   ! Loop over neighbors
   do while ((i/=mesh%lend(inr)).or.init)
      ! Get arc
      net_nnb(ic1) = net_nnb(ic1)+1
      net_inb(ic1,net_nnb(ic1)) = mesh%order(abs(mesh%list(i)))
      i = mesh%lptr(i)
      init = .false.

      ! Check arc validity
      ic0 = nicas_blk%c1_to_c0(ic1)
      jc1 = net_inb(ic1,net_nnb(ic1))
      jc0 = nicas_blk%c1_to_c0(jc1)
      do il1=1,nicas_blk%nl1
         il0 = nicas_blk%l1_to_l0(il1)
         valid_arc(ic1,il1,net_nnb(ic1)) = .true.
         if (nam%mask_check) call geom%check_arc(mpl,il0,geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0), &
       & valid_arc(ic1,il1,net_nnb(ic1)))
      end do
   end do

   ! Update
   call mpl%prog_print(ic1_loc)
end do
call mpl%prog_final

! MPI sharing
call mpl%share(nicas_blk%nc1,nc1_loc,net_nnb)
call mpl%share(nicas_blk%nc1,net_nnbmax,nc1_loc,net_inb)
call mpl%share(nicas_blk%nc1,nicas_blk%nl1,net_nnbmax,nc1_loc,valid_arc)

! Compute mesh edges distances
write(mpl%info,'(a10,a)') '','Compute mesh edges distances: '
call mpl%flush(.false.)
call mpl%prog_init(nc1_loc(mpl%myproc))
net_dnb = 1.0
!$omp parallel do schedule(static) private(ic1_loc,ic1,j,ic0,jc1,jc0,dnb,dx,dy,dz,il1,il0,djl1,jl1,jl0,H11,H22,H33), &
!$omp&                             private(H12,disthsq,distvsq,rhsq,rvsq,distnorm_network)
do ic1_loc=1,nc1_loc(mpl%myproc)
   ! Indices
   ic1 = sum(nc1_loc(0:mpl%myproc-1))+ic1_loc
   ic0 = nicas_blk%c1_to_c0(ic1)

   do j=1,net_nnb(ic1)
      ! Indices
      jc1 = net_inb(ic1,j)
      jc0 = nicas_blk%c1_to_c0(jc1)

      if (geom%mask_hor_c0(jc0)) then
         if (nicas_blk%anisotropic) then
            ! Compute longitude/latitude differences
            dx = geom%lon(jc0)-geom%lon(ic0)
            dy = geom%lat(jc0)-geom%lat(ic0)
            call lonlatmod(dx,dy)
            dx = dx*cos(0.5*(geom%lat(ic0)+geom%lat(jc0)))
         else
            ! Compute horizontal distance
            call sphere_dist(geom%lon(ic0),geom%lat(ic0),geom%lon(jc0),geom%lat(jc0),dnb)
         end if

         do il1=1,nicas_blk%nl1
            ! Index
            il0 = nicas_blk%l1_to_l0(il1)

            do djl1=-1,1
               ! Index
               jl1 = max(1,min(il1+djl1,nicas_blk%nl1))
               jl0 = nicas_blk%l1_to_l0(jl1)

               ! Check valid arc for both levels
               if (geom%mask_c0(ic0,il0).and.geom%mask_c0(jc0,jl0).and.valid_arc(ic1,il1,j).and.valid_arc(ic1,jl1,j)) then
                  ! Squared support radii
                  if (nicas_blk%anisotropic) then
                     dz = geom%vunit_c0(ic0,il0)-geom%vunit_c0(jc0,jl0)
                     H11 = 0.5*(nicas_blk%H11_c1(ic1,il1)+nicas_blk%H11_c1(jc1,jl1))
                     H22 = 0.5*(nicas_blk%H22_c1(ic1,il1)+nicas_blk%H22_c1(jc1,jl1))
                     H33 = 0.5*(nicas_blk%H33_c1(ic1,il1)+nicas_blk%H33_c1(jc1,jl1))
                     H12 = 0.5*(nicas_blk%H12_c1(ic1,il1)+nicas_blk%H12_c1(jc1,jl1))
                     net_dnb(ic1,il1,j,djl1+2) = sqrt(H11*dx**2+H22*dy**2+H33*dz**2+2.0*H12*dx*dy)*gc2gau
                  else
                     disthsq = dnb**2
                     distvsq = (geom%vunit_c0(ic0,il0)-geom%vunit_c0(jc0,jl0))**2
                     rhsq = 0.5*(nicas_blk%rh_c1(ic1,il1)**2+nicas_blk%rh_c1(jc1,jl1)**2)
                     rvsq = 0.5*(nicas_blk%rv_c1(ic1,il1)**2+nicas_blk%rv_c1(jc1,jl1)**2)
                     distnorm_network = 0.0
                     if (rhsq>0.0) then
                        distnorm_network = distnorm_network+disthsq/rhsq
                     elseif (disthsq>0.0) then
                        distnorm_network = distnorm_network+0.5*huge(0.0)
                     end if
                     if (rvsq>0.0) then
                        distnorm_network = distnorm_network+distvsq/rvsq
                     elseif (distvsq>0.0) then
                        distnorm_network = distnorm_network+0.5*huge(0.0)
                     end if
                     net_dnb(ic1,il1,j,djl1+2) = sqrt(distnorm_network)
                  end if
               end if
            end do
         end do
      end if
   end do

   ! Update
   call mpl%prog_print(ic1_loc)
end do
!$omp end parallel do
call mpl%prog_final

! MPI sharing
call mpl%share(nicas_blk%nc1,nicas_blk%nl1,net_nnbmax,3,nc1_loc,net_dnb)

! Compute distances
write(mpl%info,'(a10,a)') '','Compute distances: '
call mpl%flush(.false.)
call mpl%prog_init(nicas_blk%nsbb)
!$omp parallel do schedule(static) private(isbb,is,ic1,il1,np,np_new,jc1,jl1,k,kc1,dkl1,kl1,disttest,add_to_front,jp), &
!$omp&                             firstprivate(distnorm,plist,plist_new)
do isbb=1,nicas_blk%nsbb
   ! Indices
   is = nicas_blk%sbb_to_s(isbb)
   ic1 = nicas_blk%s_to_c1(is)
   il1 = nicas_blk%s_to_l1(is)

   ! Allocation
   allocate(distnorm(nicas_blk%nc1,nicas_blk%nl1))
   allocate(plist(nicas_blk%nc1*nicas_blk%nl1,2))
   allocate(plist_new(nicas_blk%nc1*nicas_blk%nl1,2))

   ! Initialize the front
   np = 1
   plist(1,1) = ic1
   plist(1,2) = il1
   distnorm = 1.0
   distnorm(ic1,il1) = 0.0

   do while (np>0)
      ! Propagate the front
      np_new = 0

      do ip=1,np
         ! Indices of the central point
         jc1 = plist(ip,1)
         jl1 = plist(ip,2)

         ! Loop over neighbors
         do k=1,net_nnb(jc1)
            kc1 = net_inb(jc1,k)
            do dkl1=-1,1
               kl1 = max(1,min(jl1+dkl1,nicas_blk%nl1))
               if (nicas_blk%mask_c2(kc1,kl1)) then
                  disttest = distnorm(jc1,jl1)+net_dnb(jc1,jl1,k,dkl1+2)
                  if (inf(disttest,1.0_kind_real)) then
                     ! Point is inside the support
                     if (inf(disttest,distnorm(kc1,kl1))) then
                        ! Update distance
                        distnorm(kc1,kl1) = disttest

                        ! Check if the point should be added to the front (avoid duplicates)
                        add_to_front = .true.
                        do jp=1,np_new
                           if ((plist_new(jp,1)==kc1).and.(plist_new(jp,2)==kl1)) then
                              add_to_front = .false.
                              exit
                           end if
                        end do

                        if (add_to_front) then
                           ! Add point to the front
                           np_new = np_new+1
                           plist_new(np_new,1) = kc1
                           plist_new(np_new,2) = kl1
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end do

      ! Copy new front
      np = np_new
      plist(1:np,:) = plist_new(1:np,:)
   end do

   ! Pack data
   do il1=1,nicas_blk%nl1
      do ic1=1,nicas_blk%nc1
         if (supeq(distnorm(ic1,il1),1.0_kind_real)) distnorm(ic1,il1) = mpl%msv%valr
      end do
   end do
   call nicas_blk%distnorm(isbb)%pack(mpl,nicas_blk%nc1,nicas_blk%nl1,distnorm)

   ! Release memory
   deallocate(distnorm)
   deallocate(plist)
   deallocate(plist_new)

   ! Update
   call mpl%prog_print(isbb)
end do
call mpl%prog_final

! Release memory
deallocate(lon_c1)
deallocate(lat_c1)
deallocate(net_inb)
deallocate(net_dnb)
deallocate(valid_arc)

end subroutine nicas_blk_compute_convol_network

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol_distance
! Purpose: compute convolution with a distance approach
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol_distance(nicas_blk,mpl,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry

! Local variables
integer :: nnmax,is,ic1,jc1,il1,il0,j,js,ic0,jc0,jl0,jl1,ic1bb,isbb
integer :: nn(nicas_blk%nc1bb)
integer,allocatable :: nn_index(:,:)
real(kind_real) :: nnavg,rr,disthsq,distvsq,rhsq,rvsq
real(kind_real) :: dx,dy,dz,H11,H22,H33,H12
real(kind_real),allocatable :: distnorm(:,:),nn_dist(:,:)
logical,allocatable :: valid_arc(:,:,:)

! Count nearest neighbors
write(mpl%info,'(a10,a)') '','Count nearest neighbors: '
call mpl%flush(.false.)
call mpl%prog_init(nicas_blk%nc1bb)
do ic1bb=1,nicas_blk%nc1bb
   ! Indices
   ic1 = nicas_blk%c1bb_to_c1(ic1bb)
   ic0 = nicas_blk%c1_to_c0(ic1)

   ! Research radius
   rr = sqrt(0.5*(maxval(nicas_blk%rh_c1(ic1,:))**2+nicas_blk%rhmax**2))

   ! Count nearest neighbors
   call nicas_blk%tree%count_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),rr,nn(ic1bb))

   ! Update
   call mpl%prog_print(ic1bb)
end do
call mpl%prog_final

! Print results
if (nicas_blk%nc1bb>0) then
   nnmax = maxval(nn)
   nnavg = real(sum(nn),kind_real)/real(nicas_blk%nc1bb,kind_real)
else
   nnmax = 0
   nnavg = 0
end if
write(mpl%info,'(a10,a,f8.1,a,f5.1,a,i6,a,f5.1,a)') '','Average / maximum number of neighbors to find for this task: ',nnavg, &
 & ' (', nnavg/real(nicas_blk%nc1,kind_real)*100.0,'%)  /',nnmax,' (', &
 & real(nnmax,kind_real)/real(nicas_blk%nc1,kind_real)*100.0,'%)'
call mpl%flush

! Allocation
allocate(nn_index(nnmax,nicas_blk%nc1bb))
allocate(nn_dist(nnmax,nicas_blk%nc1bb))
allocate(valid_arc(nnmax,nicas_blk%nc1bb,nicas_blk%nl1))

! Find nearest neighbors
write(mpl%info,'(a10,a)') '','Find nearest neighbors: '
call mpl%flush(.false.)
call mpl%prog_init(nicas_blk%nc1bb)
do ic1bb=1,nicas_blk%nc1bb
   ! Indices
   ic1 = nicas_blk%c1bb_to_c1(ic1bb)
   ic0 = nicas_blk%c1_to_c0(ic1)

   ! Find nearest neighbors
   call nicas_blk%tree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),nn(ic1bb),nn_index(1:nn(ic1bb),ic1bb), &
 & nn_dist(1:nn(ic1bb),ic1bb))

   ! Check arc validity
   do j=1,nn(ic1bb)
      jc1 = nn_index(j,ic1bb)
      jc0 = nicas_blk%c1_to_c0(jc1)
      do il1=1,nicas_blk%nl1
         il0 = nicas_blk%l1_to_l0(il1)
         valid_arc(j,ic1bb,il1) = (geom%mask_c0(ic0,il0).and.geom%mask_c0(jc0,il0))
         if (nam%mask_check.and.valid_arc(j,ic1bb,il1)) call geom%check_arc(mpl,il0,geom%lon(ic0),geom%lat(ic0),geom%lon(jc0), &
       & geom%lat(jc0),valid_arc(j,ic1bb,il1))
      end do
   end do

   ! Update
   call mpl%prog_print(ic1bb)
end do
call mpl%prog_final

! Compute distances
write(mpl%info,'(a10,a)') '','Compute distances: '
call mpl%flush(.false.)
call mpl%prog_init(nicas_blk%nsbb)
!$omp parallel do schedule(static) private(isbb,is,ic1,ic1bb,ic0,il1,il0,j,jl1,jl0,js,jc1,jc0,disthsq,distvsq,rhsq,rvsq), &
!$omp&                             private(dx,dy,dz,H11,H22,H33,H12) firstprivate(distnorm)
do isbb=1,nicas_blk%nsbb
   ! Indices
   is = nicas_blk%sbb_to_s(isbb)
   ic1 = nicas_blk%s_to_c1(is)
   ic1bb = nicas_blk%c1_to_c1bb(ic1)
   ic0 = nicas_blk%c1_to_c0(ic1)
   il1 = nicas_blk%s_to_l1(is)
   il0 = nicas_blk%l1_to_l0(il1)

   ! Allocation
   allocate(distnorm(nicas_blk%nc1,nicas_blk%nl1))

   ! Initialization
   distnorm = mpl%msv%valr

   ! Loop on nearest neighbors
   do j=1,nn(ic1bb)
      jc1 = nn_index(j,ic1bb)
      do jl1=1,nicas_blk%nl1
         if (nicas_blk%mask_c2(jc1,jl1)) then
            jc0 = nicas_blk%c1_to_c0(jc1)
            jl0 = nicas_blk%l1_to_l0(jl1)
            js = nicas_blk%c1l1_to_s(jc1,jl1)

            ! Check valid arc for both levels
            if (valid_arc(j,ic1bb,il1).and.valid_arc(j,ic1bb,jl1)) then
               ! Normalized distance
               if (nicas_blk%anisotropic) then
                  dx = geom%lon(jc0)-geom%lon(ic0)
                  dy = geom%lat(jc0)-geom%lat(ic0)
                  call lonlatmod(dx,dy)
                  dx = dx*cos(0.5*(geom%lat(ic0)+geom%lat(jc0)))
                  dz = geom%vunit_c0(ic0,il0)-geom%vunit_c0(jc0,jl0)
                  H11 = 0.5*(nicas_blk%H11_c1(ic1,il1)+nicas_blk%H11_c1(jc1,jl1))
                  H22 = 0.5*(nicas_blk%H22_c1(ic1,il1)+nicas_blk%H22_c1(jc1,jl1))
                  H33 = 0.5*(nicas_blk%H33_c1(ic1,il1)+nicas_blk%H33_c1(jc1,jl1))
                  H12 = 0.5*(nicas_blk%H12_c1(ic1,il1)+nicas_blk%H12_c1(jc1,jl1))
                  distnorm(jc1,jl1) = sqrt(H11*dx**2+H22*dy**2+H33*dz**2+2.0*H12*dx*dy)*gc2gau
               else
                  disthsq = nn_dist(j,ic1bb)**2
                  distvsq = (geom%vunit_c0(ic0,il0)-geom%vunit_c0(jc0,jl0))**2
                  rhsq = 0.5*(nicas_blk%rh_c1(ic1,il1)**2+nicas_blk%rh_c1(jc1,jl1)**2)
                  rvsq = 0.5*(nicas_blk%rv_c1(ic1,il1)**2+nicas_blk%rv_c1(jc1,jl1)**2)
                  if (rhsq>0.0) then
                     disthsq = disthsq/rhsq
                  elseif (disthsq>0.0) then
                     disthsq = 0.5*huge(0.0)
                  end if
                  if (rvsq>0.0) then
                     distvsq = distvsq/rvsq
                  elseif (distvsq>0.0) then
                     distvsq = 0.5*huge(0.0)
                  end if
                  distnorm(jc1,jl1) = sqrt(disthsq+distvsq)
               end if
            end if
         end if
      end do
   end do

   ! Pack data
   do il1=1,nicas_blk%nl1
      do ic1=1,nicas_blk%nc1
         if (supeq(distnorm(ic1,il1),1.0_kind_real)) distnorm(ic1,il1) = mpl%msv%valr
      end do
   end do
   call nicas_blk%distnorm(isbb)%pack(mpl,nicas_blk%nc1,nicas_blk%nl1,distnorm)

   ! Release memory
   deallocate(distnorm)

   ! Update
   call mpl%prog_print(isbb)
end do
!$omp end parallel do
call mpl%prog_final

! Release memory
deallocate(nn_index)
deallocate(nn_dist)
deallocate(valid_arc)

end subroutine nicas_blk_compute_convol_distance

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol_weights
! Purpose: compute convolution weights
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol_weights(nicas_blk,mpl,nam,geom,ctmp)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(linop_type),intent(inout) :: ctmp           ! Convolution operator

! Local variables
integer :: n_s_max,ithread,is,ic1,jc1,il1,jl1,il0,jbd,js,ic0,ic1bb,isbb
integer :: c_n_s(mpl%nthread),c_nor_n_s(mpl%nthread)
real(kind_real) :: disth,S_test
logical :: add_op
type(linop_type) :: c(mpl%nthread),c_nor(mpl%nthread)

! Allocation
n_s_max = 10*nint(real(geom%nc0*geom%nl0)/real(mpl%nthread*mpl%nproc))
do ithread=1,mpl%nthread
   c(ithread)%n_s = n_s_max
   call c(ithread)%alloc
   c_nor(ithread)%n_s = n_s_max
   call c_nor(ithread)%alloc
end do

! Initialization
write(mpl%info,'(a10,a)') '','Compute weights: '
call mpl%flush(.false.)
call mpl%prog_init(nicas_blk%nsbb)
c_n_s = 0
c_nor_n_s = 0

! Compute weights
!$omp parallel do schedule(static) private(isbb,is,ithread,ic1,ic1bb,ic0,il1,il0,jbd,jc1,jl1,js,disth,S_test,add_op)
do isbb=1,nicas_blk%nsbb
   ! Indices
   is = nicas_blk%sbb_to_s(isbb)
   ithread = 1
!$ ithread = omp_get_thread_num()+1
   ic1 = nicas_blk%s_to_c1(is)
   ic1bb = nicas_blk%c1_to_c1bb(ic1)
   ic0 = nicas_blk%c1_to_c0(ic1)
   il1 = nicas_blk%s_to_l1(is)
   il0 = nicas_blk%l1_to_l0(il1)

   ! Count convolution operations
   do jbd=1,nicas_blk%distnorm(isbb)%nbd
      ! Indices
      jc1 = nicas_blk%distnorm(isbb)%bd_to_c1(jbd)
      jl1 = nicas_blk%distnorm(isbb)%bd_to_l1(jbd)
      js = nicas_blk%c1l1_to_s(jc1,jl1)

      if (nicas_blk%double_fit) then
         ! Double fit function
         disth = sqrt(nicas_blk%distnorm(isbb)%val(jbd)**2-nicas_blk%distnormv(isbb)%val(jbd)**2)
         S_test = fit_func(mpl,nam%fit_type,disth)*((1.0+nicas_blk%coef(isbb)%val(jbd))* &
                & fit_func(mpl,nam%fit_type,nicas_blk%distnormv(isbb)%val(jbd)) &
                & -nicas_blk%coef(isbb)%val(jbd) &
                & *fit_func(mpl,nam%fit_type,nicas_blk%distnormv(isbb)%val(jbd)*nicas_blk%rfac(isbb)%val(jbd)))
      else
         ! Fit function
         S_test = fit_func(mpl,nam%fit_type,nicas_blk%distnorm(isbb)%val(jbd))
      end if
      if (nicas_blk%anisotropic) S_test = S_test*nicas_blk%Hcoef(isbb)%val(jbd)

      if (sup(S_test,S_inf)) then
         ! Store coefficient for convolution
         if (nam%lsqrt) then
            add_op = .false.
            if (nam%mpicom==1) then
               add_op = (nicas_blk%lcheck_sb(js).and.(is<=js)).or.(.not.nicas_blk%lcheck_sb(js))
            elseif (nam%mpicom==2) then
               add_op = nicas_blk%lcheck_sa(is).and.((nicas_blk%lcheck_sa(js).and.(is<=js)) &
                     & .or.(.not.nicas_blk%lcheck_sa(js)))
            end if
         else
            add_op = .true.
         end if
         if (add_op) call c(ithread)%add_op(c_n_s(ithread),is,js,S_test)

         ! Store coefficient for normalization
         add_op = nicas_blk%lcheck_sb(is).and.((nicas_blk%lcheck_sb(js).and.(is<=js)).or.(.not.nicas_blk%lcheck_sb(js)))
         if (add_op) call c_nor(ithread)%add_op(c_nor_n_s(ithread),is,js,S_test)
      end if
   end do

   ! Update
   call mpl%prog_print(isbb)
end do
!$omp end parallel do
call mpl%prog_final

! Gather data from threads
call ctmp%gather(mpl,c_n_s,c)
call nicas_blk%c_nor%gather(mpl,c_nor_n_s,c_nor)

! Release memory
do ithread=1,mpl%nthread
   call c(ithread)%dealloc
   call c_nor(ithread)%dealloc
end do

end subroutine nicas_blk_compute_convol_weights

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_mpi_c
! Purpose: compute NICAS MPI distribution, halo C
!----------------------------------------------------------------------
subroutine nicas_blk_compute_mpi_c(nicas_blk,mpl,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(geom_type),intent(in) :: geom               ! Geometry

! Local variables
integer :: isa,isb,isc,i_s,is,js
integer,allocatable :: s_to_proc(:),s_to_sc(:)
logical,allocatable :: lcheck_sc_nor(:)
character(len=1024),parameter :: subr = 'nicas_blk_compute_mpi_c'

! Allocation
allocate(nicas_blk%lcheck_sc(nicas_blk%ns))
allocate(lcheck_sc_nor(nicas_blk%ns))
allocate(s_to_proc(nicas_blk%ns))

! Define halo C
nicas_blk%lcheck_sc = nicas_blk%lcheck_sb
do i_s=1,nicas_blk%c%n_s
   is = nicas_blk%c%row(i_s)
   js = nicas_blk%c%col(i_s)
   nicas_blk%lcheck_sc(is) = .true.
   nicas_blk%lcheck_sc(js) = .true.
end do
nicas_blk%nsc = count(nicas_blk%lcheck_sc)
lcheck_sc_nor = nicas_blk%lcheck_sb
do i_s=1,nicas_blk%c_nor%n_s
   is = nicas_blk%c_nor%row(i_s)
   js = nicas_blk%c_nor%col(i_s)
   lcheck_sc_nor(is) = .true.
   lcheck_sc_nor(js) = .true.
end do
nicas_blk%nsc_nor = count(lcheck_sc_nor)

! Check halos consistency
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sa(is).and.(.not.nicas_blk%lcheck_sc(is))) call mpl%abort(subr,'point in halo A but not in halo C')
   if (nicas_blk%lcheck_sb(is).and.(.not.nicas_blk%lcheck_sc(is))) call mpl%abort(subr,'point in halo B but not in halo C')
end do

! Allocation
allocate(nicas_blk%sc_to_s(nicas_blk%nsc))
allocate(s_to_sc(nicas_blk%ns))
allocate(nicas_blk%sc_nor_to_s(nicas_blk%nsc_nor))
allocate(nicas_blk%s_to_sc_nor(nicas_blk%ns))
allocate(nicas_blk%sa_to_sc(nicas_blk%nsa))
allocate(nicas_blk%sb_to_sc(nicas_blk%nsb))
allocate(nicas_blk%sb_to_sc_nor(nicas_blk%nsb))

! Global-local conversion for halo C
s_to_sc = mpl%msv%vali
isc = 0
do is=1,nicas_blk%ns
   if (nicas_blk%lcheck_sc(is)) then
      isc = isc+1
      nicas_blk%sc_to_s(isc) = is
      s_to_sc(is) = isc
   end if
end do

! Global-local conversion for halo C (normalization)
nicas_blk%s_to_sc_nor = mpl%msv%vali
isc = 0
do is=1,nicas_blk%ns
   if (lcheck_sc_nor(is)) then
      isc = isc+1
      nicas_blk%sc_nor_to_s(isc) = is
      nicas_blk%s_to_sc_nor(is) = isc
   end if
end do

! Halos A-C conversions
do isa=1,nicas_blk%nsa
   is = nicas_blk%sa_to_s(isa)
   isc = s_to_sc(is)
   nicas_blk%sa_to_sc(isa) = isc
end do
do isb=1,nicas_blk%nsb
   is = nicas_blk%sb_to_s(isb)
   isc = s_to_sc(is)
   nicas_blk%sb_to_sc(isb) = isc
   isc = nicas_blk%s_to_sc_nor(is)
   nicas_blk%sb_to_sc_nor(isb) = isc
end do

! Local convolutions source and destination
nicas_blk%c%n_src = nicas_blk%nsc
nicas_blk%c%n_dst = nicas_blk%nsc
do i_s=1,nicas_blk%c%n_s
   nicas_blk%c%row(i_s) = s_to_sc(nicas_blk%c%row(i_s))
   nicas_blk%c%col(i_s) = s_to_sc(nicas_blk%c%col(i_s))
end do
nicas_blk%c_nor%n_src = nicas_blk%nsc
nicas_blk%c_nor%n_dst = nicas_blk%nsc
do i_s=1,nicas_blk%c_nor%n_s
   nicas_blk%c_nor%row(i_s) = nicas_blk%s_to_sc_nor(nicas_blk%c_nor%row(i_s))
   nicas_blk%c_nor%col(i_s) = nicas_blk%s_to_sc_nor(nicas_blk%c_nor%col(i_s))
end do

! Reorder convolutions
call nicas_blk%c%reorder(mpl)
call nicas_blk%c_nor%reorder(mpl)

! Setup communications
s_to_proc = geom%c0_to_proc(nicas_blk%c1_to_c0(nicas_blk%s_to_c1))
call nicas_blk%com_AC%setup(mpl,'com_AC',nicas_blk%ns,nicas_blk%nsa,nicas_blk%nsc,nicas_blk%nsa,nicas_blk%sc_to_s, &
 & nicas_blk%sa_to_sc,s_to_proc,nicas_blk%s_to_sa)

! Release memory
deallocate(lcheck_sc_nor)
deallocate(s_to_proc)

end subroutine nicas_blk_compute_mpi_c

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_normalization
! Purpose: compute normalization
!----------------------------------------------------------------------
subroutine nicas_blk_compute_normalization(nicas_blk,mpl,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry

! Local variables
integer :: il0i,i_s,ic1,ic1b,jc1b,is,js,isc,jsb,jsc,ic0,ic0a,il0,il1,ih,jv,nlr,ilr,ic,isc_add
integer,allocatable :: ineh(:,:),inev(:),ines(:,:),inec(:),order(:),isc_list(:)
integer,allocatable :: h_col(:,:,:),v_col(:,:),s_col(:,:,:),c_ind(:,:)
real(kind_real) :: S_add
real(kind_real),allocatable :: h_S(:,:,:),v_S(:,:,:),s_S(:,:,:),c_S(:,:)
real(kind_real),allocatable :: list(:),S_list(:),S_list_tmp(:)

! Compute horizontal interpolation inverse mapping
allocate(ineh(geom%nc0a,geom%nl0i))
ineh = 0
do il0i=1,geom%nl0i
   do i_s=1,nicas_blk%h(il0i)%n_s
      ic0a = nicas_blk%h(il0i)%row(i_s)
      ineh(ic0a,il0i) = ineh(ic0a,il0i)+1
   end do
end do
allocate(h_col(maxval(ineh),geom%nc0a,geom%nl0i))
allocate(h_S(maxval(ineh),geom%nc0a,geom%nl0i))
h_col = mpl%msv%vali
h_S = mpl%msv%valr
ineh = 0
do il0i=1,geom%nl0i
   do i_s=1,nicas_blk%h(il0i)%n_s
      ic0a = nicas_blk%h(il0i)%row(i_s)
      ineh(ic0a,il0i) = ineh(ic0a,il0i)+1
      h_col(ineh(ic0a,il0i),ic0a,il0i) = nicas_blk%h(il0i)%col(i_s)
      h_S(ineh(ic0a,il0i),ic0a,il0i) = nicas_blk%h(il0i)%S(i_s)
   end do
end do

! Compute vertical interpolation inverse mapping
allocate(inev(geom%nl0))
inev = 0
do i_s=1,nicas_blk%v%n_s
   il0 = nicas_blk%v%row(i_s)
   inev(il0) = inev(il0)+1
end do
allocate(v_col(maxval(inev),geom%nl0))
allocate(v_S(maxval(inev),nicas_blk%nc1b,geom%nl0))
v_col = mpl%msv%vali
v_S = mpl%msv%valr
inev = 0
do i_s=1,nicas_blk%v%n_s
   il0 = nicas_blk%v%row(i_s)
   inev(il0) = inev(il0)+1
   v_col(inev(il0),il0) = nicas_blk%v%col(i_s)
   do ic1b=1,nicas_blk%nc1b
      v_S(inev(il0),ic1b,il0) = nicas_blk%v%Svec(i_s,ic1b)
   end do
end do

! Compute subsampling interpolation inverse mapping
allocate(ines(nicas_blk%nc1b,nicas_blk%nl1))
ines = 0
do il1=1,nicas_blk%nl1
   do i_s=1,nicas_blk%s(il1)%n_s
      ic1b = nicas_blk%s(il1)%row(i_s)
      ines(ic1b,il1) = ines(ic1b,il1)+1
   end do
end do
allocate(s_col(maxval(ines),nicas_blk%nc1b,nicas_blk%nl1))
allocate(s_S(maxval(ines),nicas_blk%nc1b,nicas_blk%nl1))
s_col = mpl%msv%vali
s_S = mpl%msv%valr
ines = 0
do il1=1,nicas_blk%nl1
   do i_s=1,nicas_blk%s(il1)%n_s
      ic1b = nicas_blk%s(il1)%row(i_s)
      ines(ic1b,il1) = ines(ic1b,il1)+1
      jc1b = nicas_blk%s(il1)%col(i_s)
      jsb = nicas_blk%c1bl1_to_sb(jc1b,il1)
      jsc = nicas_blk%sb_to_sc_nor(jsb)
      s_col(ines(ic1b,il1),ic1b,il1) = jsc
      s_S(ines(ic1b,il1),ic1b,il1) = nicas_blk%s(il1)%S(i_s)
   end do
end do

! Compute convolution inverse mapping
allocate(inec(nicas_blk%nsc_nor))
inec = 0
do i_s=1,nicas_blk%c_nor%n_s
   isc = nicas_blk%c_nor%col(i_s)
   jsc = nicas_blk%c_nor%row(i_s)
   if (jsc/=isc) then
      is = nicas_blk%sc_nor_to_s(isc)
      inec(isc) = inec(isc)+1
      js = nicas_blk%sc_nor_to_s(jsc)
      inec(jsc) = inec(jsc)+1
   end if
end do
allocate(c_ind(maxval(inec),nicas_blk%nsc_nor))
allocate(c_S(maxval(inec),nicas_blk%nsc_nor))
c_ind = mpl%msv%vali
c_S = mpl%msv%valr
inec = 0
do i_s=1,nicas_blk%c_nor%n_s
   isc = nicas_blk%c_nor%col(i_s)
   jsc = nicas_blk%c_nor%row(i_s)
   if (jsc/=isc) then
      is = nicas_blk%sc_nor_to_s(isc)
      inec(isc) = inec(isc)+1
      c_ind(inec(isc),isc) = jsc
      c_S(inec(isc),isc) = nicas_blk%c_nor%S(i_s)
      js = nicas_blk%sc_nor_to_s(jsc)
      inec(jsc) = inec(jsc)+1
      c_ind(inec(jsc),jsc) = isc
      c_S(inec(jsc),jsc) = nicas_blk%c_nor%S(i_s)
   end if
end do

! Re-order indices
do isc=1,nicas_blk%nsc_nor
   ! Allocation
   allocate(order(inec(isc)))
   allocate(list(inec(isc)))

   ! Copy
   list = c_ind(1:inec(isc),isc)

   ! Order
   call qsort(inec(isc),list,order)

   ! Re-order
   c_ind(1:inec(isc),isc) = c_ind(order(1:inec(isc)),isc)
   c_S(1:inec(isc),isc) = c_S(order(1:inec(isc)),isc)

   ! Release memory
   deallocate(order)
   deallocate(list)
end do

! Allocation
allocate(nicas_blk%norm(geom%nc0a,geom%nl0))
nicas_blk%norm = mpl%msv%valr

! Compute normalization weights
do il0=1,geom%nl0
   if (nicas_blk%vlev(il0)) then
      il0i = min(il0,geom%nl0i)
      write(mpl%info,'(a10,a,i3,a)') '','Level ',nam%levs(il0),': '
      call mpl%flush(.false.)
      call mpl%prog_init(geom%nc0a)

      !$omp parallel do schedule(static) private(ic0a,ic0,nlr,isc_add,S_add,ih,ic1b,ic1,jv,il1,is,ilr,ic,isc,jsc), &
      !$omp&                             firstprivate(isc_list,S_list,S_list_tmp)
      do ic0a=1,geom%nc0a
         ! Index
         ic0 = geom%c0a_to_c0(ic0a)

         if (geom%mask_c0(ic0,il0)) then
            ! Allocation
            allocate(isc_list(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
            allocate(S_list(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
            allocate(S_list_tmp(nicas_blk%nsc_nor))

            ! Initialization
            isc_list = 0
            S_list = 0.0

            ! Adjoint interpolation
            nlr = 0
            do ih=1,ineh(ic0a,il0i)
               ic1b = h_col(ih,ic0a,il0i)
               ic1 = nicas_blk%c1b_to_c1(ic1b)
               do jv=1,inev(il0)
                  il1 = v_col(jv,il0)
                  if ((nicas_blk%l1_to_l0(il1)>=nicas_blk%vbot(ic1)).and.(nicas_blk%l1_to_l0(il1)<=nicas_blk%vtop(ic1))) then
                     do is=1,ines(ic1b,il1)
                        isc_add = s_col(is,ic1b,il1)
                        S_add = h_S(ih,ic0a,il0i)*v_S(jv,ic1b,il0)*s_S(is,ic1b,il1)
                        if (nlr==0) then
                           ilr = 1
                           nlr = 1
                        else
                           do ilr=1,nlr
                              if (isc_add==isc_list(ilr)) exit
                           end do
                           if (ilr==nlr+1) nlr = nlr+1
                        end if
                        isc_list(ilr) = isc_add
                        S_list(ilr) = S_list(ilr)+S_add
                     end do
                  end if
               end do
            end do

            ! Initialization
            S_list_tmp = 0.0
            do ilr=1,nlr
               isc = isc_list(ilr)
               S_list_tmp(isc) = S_list(ilr)
            end do

            ! Convolution
            do ilr=1,nlr
               isc = isc_list(ilr)
               do ic=1,inec(isc)
                  jsc = c_ind(ic,isc)
                  S_list_tmp(jsc) = S_list_tmp(jsc)+c_S(ic,isc)*S_list(ilr)
               end do
            end do

            ! Sum of squared values
            nicas_blk%norm(ic0a,il0) = sum(S_list_tmp**2)

            ! Normalization factor
            nicas_blk%norm(ic0a,il0) = 1.0/sqrt(nicas_blk%norm(ic0a,il0))

            ! Update
            call mpl%prog_print(ic0a)

            ! Release memory
            deallocate(isc_list)
            deallocate(S_list)
            deallocate(S_list_tmp)
         end if
      end do
      !$omp end parallel do
      call mpl%prog_final
   else
      ! Not a valid level
      nicas_blk%norm(:,il0) = 1.0
   end if
end do

! Release memory
deallocate(ineh)
deallocate(h_col)
deallocate(h_S)
deallocate(inev)
deallocate(v_col)
deallocate(v_S)
deallocate(ines)
deallocate(s_col)
deallocate(s_S)
deallocate(inec)
deallocate(c_ind)
deallocate(c_S)

end subroutine nicas_blk_compute_normalization

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_grids
! Purpose: compute grids
!----------------------------------------------------------------------
subroutine nicas_blk_compute_grids(nicas_blk,mpl,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry

! Local variables
integer :: isa,isb,isc,is,ic1,il1,ic0,il0

! Allocation
if (nicas_blk%nsa>0) then
   allocate(nicas_blk%lon_sa(nicas_blk%nsa))
   allocate(nicas_blk%lat_sa(nicas_blk%nsa))
   allocate(nicas_blk%lev_sa(nicas_blk%nsa))
end if
if (nicas_blk%nsb>0) then
   allocate(nicas_blk%lon_sb(nicas_blk%nsb))
   allocate(nicas_blk%lat_sb(nicas_blk%nsb))
   allocate(nicas_blk%lev_sb(nicas_blk%nsb))
end if
if (nicas_blk%nsc>0) then
   allocate(nicas_blk%lon_sc(nicas_blk%nsc))
   allocate(nicas_blk%lat_sc(nicas_blk%nsc))
   allocate(nicas_blk%lev_sc(nicas_blk%nsc))
end if

! Copy
if (nicas_blk%nsa>0) then
   do isa=1,nicas_blk%nsa
      is = nicas_blk%sa_to_s(isa)
      ic1 = nicas_blk%s_to_c1(is)
      il1 = nicas_blk%s_to_l1(is)
      ic0 = nicas_blk%c1_to_c0(ic1)
      il0 = nicas_blk%l1_to_l0(il1)
      nicas_blk%lon_sa(isa) = geom%lon(ic0)*rad2deg
      nicas_blk%lat_sa(isa) = geom%lat(ic0)*rad2deg
      nicas_blk%lev_sa(isa) = nam%levs(il0)
   end do
end if
if (nicas_blk%nsb>0) then
   do isb=1,nicas_blk%nsb
      is = nicas_blk%sb_to_s(isb)
      ic1 = nicas_blk%s_to_c1(is)
      il1 = nicas_blk%s_to_l1(is)
      ic0 = nicas_blk%c1_to_c0(ic1)
      il0 = nicas_blk%l1_to_l0(il1)
      nicas_blk%lon_sb(isb) = geom%lon(ic0)*rad2deg
      nicas_blk%lat_sb(isb) = geom%lat(ic0)*rad2deg
      nicas_blk%lev_sb(isb) = nam%levs(il0)
   end do
end if
if (nicas_blk%nsc>0) then
   do isc=1,nicas_blk%nsc
      is = nicas_blk%sc_to_s(isc)
      ic1 = nicas_blk%s_to_c1(is)
      il1 = nicas_blk%s_to_l1(is)
      ic0 = nicas_blk%c1_to_c0(ic1)
      il0 = nicas_blk%l1_to_l0(il1)
      nicas_blk%lon_sc(isc) = geom%lon(ic0)*rad2deg
      nicas_blk%lat_sc(isc) = geom%lat(ic0)*rad2deg
      nicas_blk%lev_sc(isc) = nam%levs(il0)
   end do
end if

end subroutine nicas_blk_compute_grids

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_adv
! Purpose: compute advection
!----------------------------------------------------------------------
subroutine nicas_blk_compute_adv(nicas_blk,mpl,rng,nam,geom,cmat_blk)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(rng_type),intent(inout) :: rng              ! Random number generator
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       ! C matrix data block

! Local variables
integer :: its,il0,ic0,ic0a,i_s,jc0,ic0d,ic0dinv
integer :: c0_to_c0d(geom%nc0),c0_to_c0dinv(geom%nc0),c0a_to_c0d(geom%nc0a),c0a_to_c0dinv(geom%nc0a)
integer,allocatable :: c0d_to_c0(:),c0dinv_to_c0(:)
real(kind_real) :: adv_lon(geom%nc0,geom%nl0),adv_lat(geom%nc0,geom%nl0)
logical :: lcheck_c0d(geom%nc0),lcheck_c0dinv(geom%nc0)

write(mpl%info,'(a7,a)') '','Compute advection'
call mpl%flush

! Allocation
allocate(nicas_blk%d(geom%nl0,2:nam%nts))
allocate(nicas_blk%dinv(geom%nl0,2:nam%nts))

! Compute interpolations
do its=2,nam%nts
   ! Local to global
   call mpl%loc_to_glb(geom%nl0,geom%nc0a,cmat_blk%adv_lon(:,:,its),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,adv_lon) ! TODO: remove that
   call mpl%loc_to_glb(geom%nl0,geom%nc0a,cmat_blk%adv_lat(:,:,its),geom%nc0,geom%c0_to_proc,geom%c0_to_c0a,.true.,adv_lat) ! TODO: remove that

   do il0=1,geom%nl0
      ! Direct      
      write(nicas_blk%d(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
      call nicas_blk%d(il0,its)%interp(mpl,rng,nam,geom,il0,geom%nc0,adv_lon(:,il0),adv_lat(:,il0),geom%mask_c0(:,il0),geom%nc0a, &
    & geom%lon(geom%c0a_to_c0),geom%lat(geom%c0a_to_c0),geom%mask_c0a(:,il0))

      ! Inverse
      write(nicas_blk%dinv(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'dinv_',il0,'_',its
      call nicas_blk%dinv(il0,its)%interp(mpl,rng,nam,geom,il0,geom%nc0,geom%lon,geom%lat,geom%mask_c0(:,il0),geom%nc0a, &
    & cmat_blk%adv_lon(:,il0,its),cmat_blk%adv_lat(:,il0,its),geom%mask_c0a(:,il0))
   end do
end do

! Define halo D/Dinv
lcheck_c0d = .false.
lcheck_c0dinv = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   lcheck_c0d(ic0) = .true.
   lcheck_c0dinv(ic0) = .true.
end do
do its=2,nam%nts
   do il0=1,geom%nl0
      do i_s=1,nicas_blk%d(il0,its)%n_s
         jc0 = nicas_blk%d(il0,its)%col(i_s)
         lcheck_c0d(jc0) = .true.
      end do
      do i_s=1,nicas_blk%dinv(il0,its)%n_s
         jc0 = nicas_blk%dinv(il0,its)%col(i_s)
         lcheck_c0dinv(jc0) = .true.
      end do
   end do
end do
nicas_blk%nc0d = count(lcheck_c0d)
nicas_blk%nc0dinv = count(lcheck_c0dinv)

! Allocation
allocate(c0d_to_c0(nicas_blk%nc0d))
allocate(c0dinv_to_c0(nicas_blk%nc0dinv))

! Global-local conversion for halo D/Dinv
c0_to_c0d = mpl%msv%vali
c0_to_c0dinv = mpl%msv%vali
ic0d = 0
ic0dinv = 0
do ic0=1,geom%nc0
   if (lcheck_c0d(ic0)) then
      ic0d = ic0d+1
      c0d_to_c0(ic0d) = ic0
      c0_to_c0d(ic0) = ic0d
   end if
   if (lcheck_c0dinv(ic0)) then
      ic0dinv = ic0dinv+1
      c0dinv_to_c0(ic0dinv) = ic0
      c0_to_c0dinv(ic0) = ic0dinv
   end if
end do

! Halos A-D/Dinv conversions
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0d = c0_to_c0d(ic0)
   ic0dinv = c0_to_c0dinv(ic0)
   c0a_to_c0d(ic0a) = ic0d
   c0a_to_c0dinv(ic0a) = ic0dinv
end do

! Local interpolation source
do its=2,nam%nts
   do il0=1,geom%nl0
      ! Direct
      nicas_blk%d(il0,its)%n_src = nicas_blk%nc0d
      do i_s=1,nicas_blk%d(il0,its)%n_s
         nicas_blk%d(il0,its)%col(i_s) = c0_to_c0d(nicas_blk%d(il0,its)%col(i_s))
      end do
      call nicas_blk%d(il0,its)%reorder(mpl)

      ! Inverse
      nicas_blk%dinv(il0,its)%n_src = nicas_blk%nc0dinv
      do i_s=1,nicas_blk%dinv(il0,its)%n_s
         nicas_blk%dinv(il0,its)%col(i_s) = c0_to_c0dinv(nicas_blk%dinv(il0,its)%col(i_s))
      end do
      call nicas_blk%dinv(il0,its)%reorder(mpl)
   end do
end do

! Setup communications
call nicas_blk%com_AD%setup(mpl,'com_AD',geom%nc0,geom%nc0a,nicas_blk%nc0d,geom%nc0a,c0d_to_c0,c0a_to_c0d,geom%c0_to_proc, &
 & geom%c0_to_c0a)
call nicas_blk%com_ADinv%setup(mpl,'com_ADinv',geom%nc0,geom%nc0a,nicas_blk%nc0dinv,geom%nc0a,c0dinv_to_c0,c0a_to_c0dinv, &
 & geom%c0_to_proc,geom%c0_to_c0a)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0d =       ',nicas_blk%nc0d
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0dinv =    ',nicas_blk%nc0dinv
call mpl%flush
do its=2,nam%nts
   do il0=1,geom%nl0
      write(mpl%info,'(a10,a,i3,a,i2,a,i8)') '','d(',il0,',',its,')%n_s =    ',nicas_blk%d(il0,its)%n_s
      call mpl%flush
      write(mpl%info,'(a10,a,i3,a,i2,a,i8)') '','dinv(',il0,',',its,')%n_s = ',nicas_blk%dinv(il0,its)%n_s
      call mpl%flush
   end do
end do

! Release memory
deallocate(c0d_to_c0)
deallocate(c0dinv_to_c0)

end subroutine nicas_blk_compute_adv

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply
! Purpose: apply NICAS method
!----------------------------------------------------------------------
subroutine nicas_blk_apply(nicas_blk,mpl,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk            ! NICAS data block
type(mpl_type),intent(inout) :: mpl                      ! MPI data
type(nam_type),intent(in) :: nam                         ! Namelist
type(geom_type),intent(in) :: geom                       ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) ! Field

! Local variables
real(kind_real) :: alpha_a(nicas_blk%nsa),alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

! Adjoint interpolation
call nicas_blk%apply_interp_ad(mpl,geom,fld,alpha_b)

! Communication
if (nam%mpicom==1) then
   ! Initialization
   alpha_c = 0.0

   ! Copy zone B into zone C
   alpha_c(nicas_blk%sb_to_sc) = alpha_b
elseif (nam%mpicom==2) then
   ! Halo reduction from zone B to zone A
   call nicas_blk%com_AB%red(mpl,alpha_b,alpha_a)

   ! Initialization
   alpha_c = 0.0

   ! Copy zone A into zone C
   alpha_c(nicas_blk%sa_to_sc) = alpha_a
end if

! Convolution
call nicas_blk%apply_convol(mpl,alpha_c)

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(mpl,alpha_c,alpha_a)

! Halo extension from zone A to zone B
call nicas_blk%com_AB%ext(mpl,alpha_a,alpha_b)

! Interpolation
call nicas_blk%apply_interp(mpl,geom,alpha_b,fld)

end subroutine nicas_blk_apply

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_from_sqrt
! Purpose: apply NICAS method from its square-root formulation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_from_sqrt(nicas_blk,mpl,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk            ! NICAS data block
type(mpl_type),intent(inout) :: mpl                      ! MPI data
type(geom_type),intent(in) :: geom                       ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) ! Field

! Local variables
real(kind_real) :: alpha(nicas_blk%nsa)

! Apply square-root adjoint
call nicas_blk%apply_sqrt_ad(mpl,geom,fld,alpha)

! Apply square-root
call nicas_blk%apply_sqrt(mpl,geom,alpha,fld)

end subroutine nicas_blk_apply_from_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_sqrt
! Purpose: apply NICAS method square-root
!----------------------------------------------------------------------
subroutine nicas_blk_apply_sqrt(nicas_blk,mpl,geom,alpha,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk          ! NICAS data block
type(mpl_type),intent(inout) :: mpl                    ! MPI data
type(geom_type),intent(in) :: geom                     ! Geometry
real(kind_real),intent(in) :: alpha(nicas_blk%nsa)     ! Subgrid field
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) ! Field

! Local variable
real(kind_real) :: alpha_a(nicas_blk%nsa),alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

! Initialization
alpha_c = 0.0

! Copy zone A into zone C
alpha_c(nicas_blk%sa_to_sc) = alpha

! Convolution
call nicas_blk%apply_convol(mpl,alpha_c)

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(mpl,alpha_c,alpha_a)

! Halo extension from zone A to zone B
call nicas_blk%com_AB%ext(mpl,alpha_a,alpha_b)

! Interpolation
call nicas_blk%apply_interp(mpl,geom,alpha_b,fld)

end subroutine nicas_blk_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_sqrt_ad
! Purpose: apply NICAS method square-root adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_sqrt_ad(nicas_blk,mpl,geom,fld,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk         ! NICAS data block
type(mpl_type),intent(inout) :: mpl                   ! MPI data
type(geom_type),intent(in) :: geom                    ! Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) ! Field
real(kind_real),intent(out) :: alpha(nicas_blk%nsa)   ! Subgrid field

! Local variable
real(kind_real) :: alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

! Adjoint interpolation
call nicas_blk%apply_interp_ad(mpl,geom,fld,alpha_b)

! Halo reduction from zone B to zone A
call nicas_blk%com_AB%red(mpl,alpha_b,alpha)

! Initialization
alpha_c = 0.0

! Copy zone A into zone C
alpha_c(nicas_blk%sa_to_sc) = alpha

! Convolution
call nicas_blk%apply_convol(mpl,alpha_c)

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(mpl,alpha_c,alpha)

end subroutine nicas_blk_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp
! Purpose: apply interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp(nicas_blk,mpl,geom,alpha,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk          ! NICAS data block
type(mpl_type),intent(inout) :: mpl                    ! MPI data
type(geom_type),intent(in) :: geom                     ! Geometry
real(kind_real),intent(in) :: alpha(nicas_blk%nsb)     ! Subgrid field
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) ! Field

! Local variables
real(kind_real) :: gamma(nicas_blk%nc1b,nicas_blk%nl1),delta(nicas_blk%nc1b,geom%nl0)

! Subsampling horizontal interpolation
call nicas_blk%apply_interp_s(mpl,alpha,gamma)

! Vertical interpolation
call nicas_blk%apply_interp_v(mpl,geom,gamma,delta)

! Horizontal interpolation
call nicas_blk%apply_interp_h(mpl,geom,delta,fld)

! Normalization
fld = fld*nicas_blk%norm

end subroutine nicas_blk_apply_interp

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_ad
! Purpose: apply interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_ad(nicas_blk,mpl,geom,fld,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk         ! NICAS data block
type(mpl_type),intent(inout) :: mpl                   ! MPI data
type(geom_type),intent(in) :: geom                    ! Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) ! Field
real(kind_real),intent(out) :: alpha(nicas_blk%nsb)   ! Subgrid field

! Local variables
real(kind_real) :: gamma(nicas_blk%nc1b,nicas_blk%nl1),delta(nicas_blk%nc1b,geom%nl0)
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0)

! Normalization
fld_tmp = fld*nicas_blk%norm

! Horizontal interpolation
call nicas_blk%apply_interp_h_ad(mpl,geom,fld_tmp,delta)

! Vertical interpolation
call nicas_blk%apply_interp_v_ad(mpl,geom,delta,gamma)

! Subsampling horizontal interpolation
call nicas_blk%apply_interp_s_ad(mpl,gamma,alpha)

end subroutine nicas_blk_apply_interp_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_h
! Purpose: apply horizontal interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_h(nicas_blk,mpl,geom,delta,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                ! NICAS data block
type(mpl_type),intent(inout) :: mpl                          ! MPI data
type(geom_type),intent(in) :: geom                           ! Geometry
real(kind_real),intent(in) :: delta(nicas_blk%nc1b,geom%nl0) ! Subset Sc1 field, full levels
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0)       ! Field

! Local variables
integer :: il0

! Horizontal interpolation
!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   if (nicas_blk%vlev(il0)) then
      call nicas_blk%h(min(il0,geom%nl0i))%apply(mpl,delta(:,il0),fld(:,il0),msdst=.false.)
   else
      fld(:,il0) = mpl%msv%valr
   end if
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_h

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_h_ad
! Purpose: apply horizontal interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_h_ad(nicas_blk,mpl,geom,fld,delta)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                 ! NICAS data block
type(mpl_type),intent(inout) :: mpl                           ! MPI data
type(geom_type),intent(in) :: geom                            ! Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0)         ! Field
real(kind_real),intent(out) :: delta(nicas_blk%nc1b,geom%nl0) ! Subset Sc1 field, full levels

! Local variables
integer :: il0

!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   if (nicas_blk%vlev(il0)) then
      call nicas_blk%h(min(il0,geom%nl0i))%apply_ad(mpl,fld(:,il0),delta(:,il0))
   else
      delta(:,il0) = mpl%msv%valr
   end if
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_h_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_v
! Purpose: apply vertical interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_v(nicas_blk,mpl,geom,gamma,delta)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                     ! NICAS data block
type(mpl_type),intent(inout) :: mpl                               ! MPI data
type(geom_type),intent(in) :: geom                                ! Geometry
real(kind_real),intent(in) :: gamma(nicas_blk%nc1b,nicas_blk%nl1) ! Subset Sc1 field, limited levels
real(kind_real),intent(out) :: delta(nicas_blk%nc1b,geom%nl0)     ! Subset Sc1 field, full levels

! Local variables
integer :: ic1b,il0
real(kind_real),allocatable :: gamma_tmp(:),delta_tmp(:)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1b,il0) firstprivate(gamma_tmp,delta_tmp)
do ic1b=1,nicas_blk%nc1b
   ! Allocation
   allocate(gamma_tmp(nicas_blk%nl1))
   allocate(delta_tmp(geom%nl0))

   ! Copy data
   gamma_tmp = gamma(ic1b,:)

   ! Apply interpolation
   call nicas_blk%v%apply(mpl,gamma_tmp,delta_tmp,ivec=ic1b,msdst=.false.)

   ! Copy data
   do il0=1,geom%nl0
      if (nicas_blk%vlev(il0)) then
         delta(ic1b,il0) = delta_tmp(il0)
      else
         delta(ic1b,il0) = mpl%msv%valr
      end if
   end do

   ! Release memory
   deallocate(gamma_tmp)
   deallocate(delta_tmp)
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_v

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_v_ad
! Purpose: apply vertical interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_v_ad(nicas_blk,mpl,geom,delta,gamma)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                      ! NICAS data block
type(mpl_type),intent(inout) :: mpl                                ! MPI data
type(geom_type),intent(in) :: geom                                 ! Geometry
real(kind_real),intent(in) :: delta(nicas_blk%nc1b,geom%nl0)       ! Subset Sc1 field, full levels
real(kind_real),intent(out) :: gamma(nicas_blk%nc1b,nicas_blk%nl1) ! Subset Sc1 field, limited levels

! Local variables
integer :: ic1b
real(kind_real),allocatable :: gamma_tmp(:),delta_tmp(:)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1b) firstprivate(gamma_tmp,delta_tmp)
do ic1b=1,nicas_blk%nc1b
   ! Allocation
   allocate(gamma_tmp(nicas_blk%nl1))
   allocate(delta_tmp(geom%nl0))

   ! Copy data
   delta_tmp = delta(ic1b,:)

   ! Apply interpolation
   call nicas_blk%v%apply_ad(mpl,delta_tmp,gamma_tmp,ivec=ic1b)

   ! Copy data
   gamma(ic1b,:) = gamma_tmp

   ! Release memory
   deallocate(gamma_tmp)
   deallocate(delta_tmp)
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_v_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_s
! Purpose: apply subsampling interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_s(nicas_blk,mpl,alpha,gamma)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                      ! NICAS data block
type(mpl_type),intent(inout) :: mpl                                ! MPI data
real(kind_real),intent(in) :: alpha(nicas_blk%nsb)                 ! Subgrid field
real(kind_real),intent(out) :: gamma(nicas_blk%nc1b,nicas_blk%nl1) ! Subset Sc1 field, limited levels

! Local variables
integer :: isb,il1
real(kind_real) :: beta(nicas_blk%nc1b,nicas_blk%nl1)

! Initialization
beta = 0.0

! Copy
!$omp parallel do schedule(static) private(isb)
do isb=1,nicas_blk%nsb
   beta(nicas_blk%sb_to_c1b(isb),nicas_blk%sb_to_l1(isb)) = alpha(isb)
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,nicas_blk%nl1
   call nicas_blk%s(il1)%apply(mpl,beta(:,il1),gamma(:,il1),msdst=.false.)
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_s

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_interp_s_ad
! Purpose: apply subsampling interpolation adjoint
!----------------------------------------------------------------------
subroutine nicas_blk_apply_interp_s_ad(nicas_blk,mpl,gamma,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                     ! NICAS data block
type(mpl_type),intent(inout) :: mpl                               ! MPI data
real(kind_real),intent(in) :: gamma(nicas_blk%nc1b,nicas_blk%nl1) ! Subset Sc1 field, limited levels
real(kind_real),intent(out) :: alpha(nicas_blk%nsb)               ! Subgrid field

! Local variables
integer :: il1,isb
real(kind_real) :: beta(nicas_blk%nc1b,nicas_blk%nl1)

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,nicas_blk%nl1
   call nicas_blk%s(il1)%apply_ad(mpl,gamma(:,il1),beta(:,il1))
end do
!$omp end parallel do

! Copy
!$omp parallel do schedule(static) private(isb)
do isb=1,nicas_blk%nsb
   alpha(isb) = beta(nicas_blk%sb_to_c1b(isb),nicas_blk%sb_to_l1(isb))
end do
!$omp end parallel do

end subroutine nicas_blk_apply_interp_s_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_convol
! Purpose: apply convolution
!----------------------------------------------------------------------
subroutine nicas_blk_apply_convol(nicas_blk,mpl,alpha)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk         ! NICAS data block
type(mpl_type),intent(inout) :: mpl                   ! MPI data
real(kind_real),intent(inout) :: alpha(nicas_blk%nsc) ! Subgrid field

! Apply linear operator, symmetric
call nicas_blk%c%apply_sym(mpl,alpha)

end subroutine nicas_blk_apply_convol

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_adv
! Purpose: apply advection
!----------------------------------------------------------------------
subroutine nicas_blk_apply_adv(nicas_blk,mpl,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                           ! NICAS data block
type(mpl_type),intent(inout) :: mpl                                     ! MPI data
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Field

! Local variables
integer :: its,iv,il0
real(kind_real) :: fld_d(nicas_blk%nc0d,geom%nl0)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Halo extension from zone A to zone D
      call nicas_blk%com_AD%ext(mpl,geom%nl0,fld(:,:,iv,its),fld_d)

      ! Interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         if (nicas_blk%vlev(il0)) call nicas_blk%d(il0,its)%apply(mpl,fld_d(:,il0),fld(:,il0,iv,its),msdst=.false.)
      end do
      !$omp end parallel do
   end do
end do

end subroutine nicas_blk_apply_adv

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_adv_ad
! Purpose: apply advection
!----------------------------------------------------------------------
subroutine nicas_blk_apply_adv_ad(nicas_blk,mpl,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                           ! NICAS data block
type(mpl_type),intent(inout) :: mpl                                     ! MPI data
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Field

! Local variables
integer :: its,iv,il0
real(kind_real) :: fld_d(nicas_blk%nc0d,geom%nl0)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Adjoint interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         if (nicas_blk%vlev(il0)) call nicas_blk%d(il0,its)%apply_ad(mpl,fld(:,il0,iv,its),fld_d(:,il0))
      end do
      !$omp end parallel do

      ! Halo reduction from zone D to zone A
      call nicas_blk%com_AD%red(mpl,geom%nl0,fld_d,fld(:,:,iv,its))
   end do
end do

end subroutine nicas_blk_apply_adv_ad

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply_adv_inv
! Purpose: apply inverse advection
!----------------------------------------------------------------------
subroutine nicas_blk_apply_adv_inv(nicas_blk,mpl,nam,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk                           ! NICAS data block
type(mpl_type),intent(inout) :: mpl                                     ! MPI data
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Field

! Local variables
integer :: its,iv,il0
real(kind_real) :: fld_dinv(nicas_blk%nc0dinv,geom%nl0)

do its=2,nam%nts
   do iv=1,nam%nv
      ! Halo extension from zone A to zone Dinv
      call nicas_blk%com_ADinv%ext(mpl,geom%nl0,fld(:,:,iv,its),fld_dinv)

      ! Interpolation
      !$omp parallel do schedule(static) private(il0)
      do il0=1,geom%nl0
         if (nicas_blk%vlev(il0)) call nicas_blk%dinv(il0,its)%apply(mpl,fld_dinv(:,il0),fld(:,il0,iv,its),msdst=.false.)
      end do
      !$omp end parallel do
   end do
end do

end subroutine nicas_blk_apply_adv_inv

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_adjoint
! Purpose: test NICAS adjoint accuracy
!----------------------------------------------------------------------
subroutine nicas_blk_test_adjoint(nicas_blk,mpl,rng,nam,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(rng_type),intent(inout) :: rng           ! Random number generator
type(nam_type),intent(in) :: nam              ! Namelist
type(geom_type),intent(in) :: geom            ! Geometry

! Local variables
integer :: isb
real(kind_real) :: sum1,sum2
real(kind_real) :: alpha(nicas_blk%nsb),alpha_save(nicas_blk%nsb),alpha_c(nicas_blk%nsc)
real(kind_real) :: gamma(nicas_blk%nc1b,nicas_blk%nl1),gamma_save(nicas_blk%nc1b,nicas_blk%nl1)
real(kind_real) :: delta(nicas_blk%nc1b,geom%nl0),delta_save(nicas_blk%nc1b,geom%nl0)
real(kind_real) :: fld(geom%nc0a,geom%nl0),fld_save(geom%nc0a,geom%nl0)
real(kind_real) :: fld1(geom%nc0a,geom%nl0),fld1_save(geom%nc0a,geom%nl0)
real(kind_real) :: fld2(geom%nc0a,geom%nl0),fld2_save(geom%nc0a,geom%nl0)
real(kind_real),allocatable :: alpha1(:),alpha1_save(:),alpha2(:),alpha2_save(:)

! Interpolation (subsampling)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha_save)
gamma_save = 0
do isb=1,nicas_blk%nsb
   call rng%rand_real(0.0_kind_real,1.0_kind_real,gamma_save(nicas_blk%sb_to_c1b(isb),nicas_blk%sb_to_l1(isb)))
end do

! Adjoint test
call nicas_blk%apply_interp_s(mpl,alpha_save,gamma)
call nicas_blk%apply_interp_s_ad(mpl,gamma_save,alpha)

! Print result
call mpl%dot_prod(alpha,alpha_save,sum1)
call mpl%dot_prod(gamma,gamma_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (subsampling): ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

! Interpolation (vertical)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,gamma_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,delta_save)

! Adjoint test
call nicas_blk%apply_interp_v(mpl,geom,gamma_save,delta)
call nicas_blk%apply_interp_v_ad(mpl,geom,delta_save,gamma)

! Print result
call mpl%dot_prod(gamma,gamma_save,sum1)
call mpl%dot_prod(delta,delta_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (vertical):    ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

! Interpolation (horizontal)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,delta_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_save)

! Adjoint test
call nicas_blk%apply_interp_h(mpl,geom,delta_save,fld)
call nicas_blk%apply_interp_h_ad(mpl,geom,fld_save,delta)

! Print result
call mpl%dot_prod(delta,delta_save,sum1)
call mpl%dot_prod(fld,fld_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (horizontal):  ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

! Interpolation (total)

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_save)

! Adjoint test
call nicas_blk%apply_interp(mpl,geom,alpha_save,fld)
call nicas_blk%apply_interp_ad(mpl,geom,fld_save,alpha)

! Print result
call mpl%dot_prod(alpha,alpha_save,sum1)
call mpl%dot_prod(fld,fld_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Interpolation adjoint test (total):       ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

! Allocation
allocate(alpha1(nicas_blk%nsc))
allocate(alpha1_save(nicas_blk%nsc))
allocate(alpha2(nicas_blk%nsc))
allocate(alpha2_save(nicas_blk%nsc))

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)
alpha1 = alpha1_save
alpha2 = alpha2_save

! Adjoint test
call nicas_blk%apply_convol(mpl,alpha1)
call nicas_blk%apply_convol(mpl,alpha2)

! Print result
call mpl%dot_prod(alpha1,alpha2_save,sum1)
call mpl%dot_prod(alpha2,alpha1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Convolution adjoint test:                 ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

! Release memory
deallocate(alpha1)
deallocate(alpha1_save)
deallocate(alpha2)
deallocate(alpha2_save)

! Allocation
allocate(alpha1(nicas_blk%nsa))
allocate(alpha1_save(nicas_blk%nsb))
allocate(alpha2(nicas_blk%nsb))
allocate(alpha2_save(nicas_blk%nsa))

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)

! Adjoint test
call nicas_blk%com_AB%red(mpl,alpha1_save,alpha1)
call nicas_blk%com_AB%ext(mpl,alpha2_save,alpha2)

! Print result
call mpl%dot_prod(alpha1,alpha2_save,sum1)
call mpl%dot_prod(alpha2,alpha1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Communication AB adjoint test:            ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

! Release memory
deallocate(alpha1)
deallocate(alpha1_save)
deallocate(alpha2)
deallocate(alpha2_save)

! Allocation
allocate(alpha1(nicas_blk%nsa))
allocate(alpha1_save(nicas_blk%nsc))
allocate(alpha2(nicas_blk%nsc))
allocate(alpha2_save(nicas_blk%nsa))

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)

! Adjoint test
call nicas_blk%com_AC%red(mpl,alpha1_save,alpha1)
call nicas_blk%com_AC%ext(mpl,alpha2_save,alpha2)

! Print result
call mpl%dot_prod(alpha1,alpha2_save,sum1)
call mpl%dot_prod(alpha2,alpha1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Communication AC adjoint test:            ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

! Allocation
deallocate(alpha1)
deallocate(alpha1_save)
deallocate(alpha2)
deallocate(alpha2_save)
allocate(alpha1(nicas_blk%nsa))
allocate(alpha1_save(nicas_blk%nsa))
allocate(alpha2(nicas_blk%nsa))
allocate(alpha2_save(nicas_blk%nsa))

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,alpha2_save)
alpha1 = alpha1_save
alpha2 = alpha2_save

! Adjoint test
alpha_c = 0.0
alpha_c(nicas_blk%sa_to_sc) = alpha1
call nicas_blk%apply_convol(mpl,alpha_c)
call nicas_blk%com_AC%red(mpl,alpha_c,alpha1)
alpha_c = 0.0
alpha_c(nicas_blk%sa_to_sc) = alpha2
call nicas_blk%apply_convol(mpl,alpha_c)
call nicas_blk%com_AC%red(mpl,alpha_c,alpha2)

! Print result
call mpl%dot_prod(alpha1,alpha2_save,sum1)
call mpl%dot_prod(alpha2,alpha1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Convolution / communication adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld2_save)
fld1 = fld1_save
fld2 = fld2_save

! Adjoint test
if (nam%lsqrt) then
   call nicas_blk%apply_from_sqrt(mpl,geom,fld1)
   call nicas_blk%apply_from_sqrt(mpl,geom,fld2)
else
   call nicas_blk%apply(mpl,nam,geom,fld1)
   call nicas_blk%apply(mpl,nam,geom,fld2)
end if

! Print result
call mpl%dot_prod(fld1,fld2_save,sum1)
call mpl%dot_prod(fld2,fld1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','NICAS adjoint test:                       ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

! Release memory
deallocate(alpha1)
deallocate(alpha1_save)
deallocate(alpha2)
deallocate(alpha2_save)

end subroutine nicas_blk_test_adjoint

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_dirac
! Purpose: apply NICAS to diracs
!----------------------------------------------------------------------
subroutine nicas_blk_test_dirac(nicas_blk,mpl,nam,geom,bpar,io)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(nam_type),intent(in) :: nam              ! Namelist
type(geom_type),intent(in) :: geom            ! Geometry
type(bpar_type),intent(in) :: bpar            ! Block parameters
type(io_type),intent(in) :: io                ! I/O

! Local variables
integer :: il0,idir
real(kind_real) :: val,valmin_tot,valmax_tot
real(kind_real) :: fld(geom%nc0a,geom%nl0)
character(len=1024) :: filename

! Associate
associate(ib=>nicas_blk%ib)

! Generate dirac field
fld = 0.0
do idir=1,geom%ndir
   if (geom%iprocdir(idir)==mpl%myproc) fld(geom%ic0adir(idir),geom%il0dir(idir)) = 1.0
end do

! Apply NICAS method
if (nam%lsqrt) then
   call nicas_blk%apply_from_sqrt(mpl,geom,fld)
else
   call nicas_blk%apply(mpl,nam,geom,fld)
end if

! Write field
filename = trim(nam%prefix)//'_dirac'
call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)
call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_dirac',fld)

! Print results
write(mpl%info,'(a7,a)') '','Values at dirac points:'
call mpl%flush
do idir=1,geom%ndir
   if (geom%iprocdir(idir)==mpl%myproc) val = fld(geom%ic0adir(idir),geom%il0dir(idir))
   call mpl%f_comm%broadcast(val,geom%iprocdir(idir)-1)
   if (mpl%msv%isnot(val)) then
      write(mpl%info,'(a10,f6.1,a,f6.1,a,f10.7)') '',geom%londir(idir)*rad2deg,' / ',geom%latdir(idir)*rad2deg,': ',val
      call mpl%flush
   else
      write(mpl%info,'(a10,f6.1,a,f6.1,a)') '',geom%londir(idir)*rad2deg,' / ',geom%latdir(idir)*rad2deg,': missing value'
      call mpl%flush
   end if
end do
write(mpl%info,'(a7,a)') '','Min - max: '
call mpl%flush
do il0=1,geom%nl0
   call mpl%f_comm%allreduce(minval(fld(:,il0),mask=geom%mask_c0a(:,il0)),valmin_tot,fckit_mpi_min())
   call mpl%f_comm%allreduce(maxval(fld(:,il0),mask=geom%mask_c0a(:,il0)),valmax_tot,fckit_mpi_max())
   if (mpl%msv%isnot(valmin_tot).or.mpl%msv%isnot(valmax_tot)) then
      write(mpl%info,'(a10,a,i3,a,f10.7,a,f10.7)') '','Level ',nam%levs(il0),': ',valmin_tot,' - ',valmax_tot
      call mpl%flush
   else
      write(mpl%info,'(a10,a,i3,a)') '','Level ',nam%levs(il0),': missing values'
      call mpl%flush
   end if
end do

! End associate
end associate

end subroutine nicas_blk_test_dirac

end module type_nicas_blk
