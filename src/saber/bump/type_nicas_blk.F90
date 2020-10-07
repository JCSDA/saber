!----------------------------------------------------------------------
! Module: type_nicas_blk
! Purpose: NICAS data block derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_nicas_blk

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_min,fckit_mpi_max
use netcdf
!$ use omp_lib
use tools_const, only: pi,req,reqkm,deg2rad,rad2deg
use tools_func, only: gc2gau,lonlatmod,lonlathash,sphere_dist,fit_func
use tools_kinds, only: kind_real,nc_kind_real,huge_real
use tools_qsort, only: qsort
use tools_repro, only: supeq,sup,inf,eq
use tools_samp, only: initialize_sampling
use type_bpar, only: bpar_type
use type_cmat_blk, only: cmat_blk_type
use type_com, only: com_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_tree, only: tree_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

real(kind_real),parameter :: sqrt_r = 0.725_kind_real ! Square-root factor on support radius (empirical)
real(kind_real),parameter :: sqrt_h = 0.725_kind_real ! Square-root factor on LCT (empirical)
real(kind_real),parameter :: S_inf = 1.0e-2_kind_real ! Minimum value for the convolution coefficients

! Ball data derived type
type balldata_type
   integer :: nbd                        ! Number of values
   integer,allocatable :: bd_to_c1u(:)   ! Ball data index to subset Sc1 in universe
   integer,allocatable :: bd_to_l1(:)    ! Ball data index to subset Sl1
   real(kind_real),allocatable :: val(:) ! Values
contains
   procedure :: alloc => balldata_alloc
   procedure :: dealloc => balldata_dealloc
   procedure :: pack => balldata_pack
end type balldata_type

! NICAS block derived type
type nicas_blk_type
   ! General parameters
   integer :: ib                                   ! Block index
   character(len=1024) :: subsamp                  ! Subsampling structure
   logical :: sqrt_rescaling                       ! Square-root rescaling flag
   logical :: anisotropic                          ! Anisotropic tensor flag
   logical :: verbosity                            ! Verbosity flag
   logical :: var2d                                ! Flag for 2D variables in 3D fields

   ! Subset Sc1
   integer :: nc1                                  ! Number of points in subset Sc1
   integer :: nc1u                                 ! Number of points in subset Sc1, universe
   integer,allocatable :: c1u_to_c1(:)             ! Subset Sc1, universe to global
   integer,allocatable :: c1_to_c1u(:)             ! Subset Sc1, global to universe
   real(kind_real),allocatable :: lon_c1u(:)       ! Longitudes on subset Sc1, universe
   real(kind_real),allocatable :: lat_c1u(:)       ! Latitudes on subset Sc1, universe
   real(kind_real),allocatable :: vunit_c1u(:,:)   ! Latitudes on subset Sc1, universe
   logical,allocatable :: gmask_c1u(:,:)           ! Mask on subset Sc1, universe
   logical,allocatable :: gmask_hor_c1u(:)         ! Union of horizontal masks on subset Sc1, universe
   real(kind_real),allocatable :: lon_c1a(:)       ! Longitudes on subset Sc1, halo A
   real(kind_real),allocatable :: lat_c1a(:)       ! Latitudes on subset Sc1, halo A
   logical,allocatable :: gmask_c1a(:,:)           ! Mask on subset Sc1, halo A
   logical,allocatable :: gmask_hor_c1a(:)         ! Union of horizontal masks on subset Sc1, halo A
   integer :: nc1a                                 ! Number of points in subset Sc1 on halo A
   integer,allocatable :: c1a_to_c1(:)             ! Subset Sc1, halo A to global
   integer,allocatable :: c1u_to_c1a(:)            ! Subset Sc1, universe to halo A
   integer,allocatable :: c1b_to_c1u(:)            ! Subset Sc1, halo B to global
   integer,allocatable :: c1u_to_c1b(:)            ! Subset Sc1, global to halo B
   integer :: nc1bb                                ! Number of points in subset Sc1 on halo B (extended)
   integer,allocatable :: c1u_to_c1bb(:)           ! Subset Sc1, universe to halo B (extended)
   integer,allocatable :: c1bb_to_c1u(:)           ! Subset Sc1, halo B (extended) to universe
   type(com_type) :: com_AU                        ! Communication between halo A and universe on subset Sc1

   ! Link between subset Sc1 and subset Sc0
   integer,allocatable :: c1_to_c0(:)              ! Subset Sc1 to subset Sc0
   integer,allocatable :: c1a_to_c0a(:)            ! Halo A, subset Sc1 to subset Sc0

   ! Subset Sc2 geometry
   integer,allocatable :: nc2(:)                   ! Number of points in subset Sc2
   integer,allocatable :: nc2u(:)                  ! Number of points in subset Sc2, universe
   logical,allocatable :: gmask_c2u(:,:)           ! Mask from subset Sc2 to subgrid

   ! Subgrid geometry
   integer :: nsu                                  ! Number of subgrid nodes, universe
   integer,allocatable :: su_to_s(:)               ! Subgrid, universe to global
   integer,allocatable :: su_to_sb(:)              ! Subgrid, universe to halo B
   integer,allocatable :: sb_to_su(:)              ! Subgrid, halo B to universe
   logical,allocatable :: lcheck_sa(:)             ! Detection of halo A on subgrid
   integer,allocatable :: sa_to_su(:)              ! Subgrid, halo A to universe
   logical,allocatable :: lcheck_sb(:)             ! Detection of halo B on subgrid
   integer :: nsbb                                 ! Number of points in subgrid on halo B (extended)
   integer,allocatable :: sbb_to_su(:)             ! Subgrid, halo B (extended) to universe
   integer,allocatable :: sc_to_su(:)              ! Subgrid, halo C to universe
   integer :: nsc_nor                              ! Number of subgrid nodes on halo C (extended for normalization)
   integer,allocatable :: sa_to_sc_nor(:)          ! Subgrid, halo A to halo C (extended for normalization)
   integer,allocatable :: sb_to_sc_nor(:)          ! Subgrid, halo B to halo C (extended for normalization)

   ! Link between subgrid and subset Sc1
   integer,allocatable :: su_to_c1u(:)             ! Subgrid, universe, to subset Sc1, universe
   integer,allocatable :: su_to_l1(:)              ! Subgrid, universe, to subset Sl1
   integer,allocatable :: c1ul1_to_su(:,:)         ! Grid Gv to subgrid, universe
   integer,allocatable :: c1bl1_to_sb(:,:)         ! Halo B, subset Sc1 to subgrid

   ! Vertical geometry
   integer :: il0_first                            ! First valid level
   integer :: il0_last                             ! Last valid level
   logical,allocatable :: slev(:)                  ! Selected levels
   integer,allocatable :: vbot(:)                  ! Bottom level in grid Gh
   integer,allocatable :: vtop(:)                  ! Top level in grid Gh
   integer,allocatable :: l1_to_l0(:)              ! Subset Sl1 to subset Sl0
   integer,allocatable :: l0_to_l1(:)              ! Subset Sl0 to subset Sl1

   ! Support radius
   real(kind_real),allocatable :: rhs_avg(:)       ! Average sampling horizontal support radius at each level

   ! Convolution parameters
   real(kind_real) :: rhmax                        ! Maximum horizontal support radius
   real(kind_real),allocatable :: rh_c1u(:,:)      ! Horizontal support radius on subset Sc1
   real(kind_real),allocatable :: rv_c1u(:,:)      ! Vertical support radius on subset Sc1
   real(kind_real),allocatable :: H11_c1u(:,:)     ! Local correlation tensor, component 11, on subset Sc1
   real(kind_real),allocatable :: H22_c1u(:,:)     ! Local correlation tensor, component 22, on subset Sc1
   real(kind_real),allocatable :: H33_c1u(:,:)     ! Local correlation tensor, component 33, on subset Sc1
   real(kind_real),allocatable :: H12_c1u(:,:)     ! Local correlation tensor, component 12, on subset Sc1
   type(balldata_type),allocatable :: distnorm(:)  ! Normalized distance
   type(balldata_type),allocatable :: Hcoef(:)     ! Tensor coefficient on subset Sc1

   ! Normalization
   type(linop_type) :: c_nor                       ! Convolution (extended for normalization)
   real(kind_real),allocatable :: inorm_nor(:)     ! Internal normalization factor (extended for normalization)
   type(com_type) :: com_AC_nor                    ! Communication between halos A and C (extended for normalization)
   real(kind_real),allocatable :: smoother_norm(:) ! Smoother normalization

   ! Tree
   type(tree_type) :: tree                         ! Tree

   ! Required data to apply NICAS

   ! Parameters
   integer :: mpicom                               ! Number of communication steps
   integer :: lsqrt                                ! Square-root flag
   integer :: grid_hash                            ! Grid hash

   ! Number of points
   integer :: nc0a                                 ! Number of points in subset Sc0 on halo A (required for I/O)
   integer :: nc1b                                 ! Number of points in subset Sc1 on halo B
   integer :: nl1                                  ! Number of levels in subset Sl1
   integer :: ns                                   ! Number of subgrid nodes
   integer :: nsa                                  ! Number of subgrid nodes on halo A
   integer :: nsb                                  ! Number of subgrid nodes on halo B
   integer :: nsc                                  ! Number of subgrid nodes on halo C

   ! Valid levels
   logical,allocatable :: vlev(:)                  ! Valid levels

   ! Local to global
   integer,allocatable :: sa_to_s(:)               ! Subgrid, halo A to global
   real(kind_real),allocatable :: hash_sa(:)       ! Hash value based on lon/lat/lev

   ! Inter-halo conversions
   integer,allocatable :: sa_to_sc(:)              ! Subgrid, halo A to halo C
   integer,allocatable :: sb_to_sc(:)              ! Subgrid, halo B to halo C

   ! Linear operators
   type(linop_type) :: c                           ! Convolution
   type(linop_type),allocatable :: h(:)            ! Horizontal interpolation
   type(linop_type) :: v                           ! Vertical interpolation
   type(linop_type),allocatable :: s(:)            ! Subsample interpolation

   ! Copy conversions
   integer,allocatable :: sb_to_c1b(:)             ! Subgrid to subset Sc1 on halo B
   integer,allocatable :: sb_to_l1(:)              ! Subgrid to subset Sl1 on halo B

   ! Normalization
   real(kind_real),allocatable :: inorm(:)         ! Internal normalization factor
   real(kind_real),allocatable :: norm(:,:)        ! Normalization factor

   ! Localization weights
   real(kind_real),allocatable :: coef_ens(:,:)    ! Ensemble coefficient square-root
   real(kind_real) :: wgt                          ! Main weight

   ! Communications
   type(com_type) :: com_AB                        ! Communication between halos A and B
   type(com_type) :: com_AC                        ! Communication between halos A and C

   ! Smoother data
   logical :: smoother                             ! Smoother flag
   logical :: horizontal                           ! Horizontal application flag

   ! Required data to write grids
   real(kind_real),allocatable :: lon_sa(:)        ! Subgrid, halo A longitudes
   real(kind_real),allocatable :: lat_sa(:)        ! Subgrid, halo A latitudes
   integer,allocatable :: lev_sa(:)                ! Subgrid, halo A levels
   real(kind_real),allocatable :: lon_sb(:)        ! Subgrid, halo A longitudes
   real(kind_real),allocatable :: lat_sb(:)        ! Subgrid, halo A latitudes
   integer,allocatable :: lev_sb(:)                ! Subgrid, halo A levels
   real(kind_real),allocatable :: lon_sc(:)        ! Subgrid, halo A longitudes
   real(kind_real),allocatable :: lat_sc(:)        ! Subgrid, halo A latitudes
   integer,allocatable :: lev_sc(:)                ! Subgrid, halo A levels
contains
   procedure :: partial_dealloc => nicas_blk_partial_dealloc
   procedure :: dealloc => nicas_blk_dealloc
   procedure :: read => nicas_blk_read
   procedure :: write => nicas_blk_write
   procedure :: write_grids => nicas_blk_write_grids
   procedure :: buffer_size => nicas_blk_buffer_size
   procedure :: serialize => nicas_blk_serialize
   procedure :: deserialize => nicas_blk_deserialize
   procedure :: nicas_blk_compute_parameters
   procedure :: nicas_blk_compute_parameters_horizontal_smoother
   generic :: compute_parameters => nicas_blk_compute_parameters,nicas_blk_compute_parameters_horizontal_smoother
   procedure :: compute_sampling_c1 => nicas_blk_compute_sampling_c1
   procedure :: compute_mpi_a => nicas_blk_compute_mpi_a
   procedure :: compute_sampling_v => nicas_blk_compute_sampling_v
   procedure :: compute_sampling_c2 => nicas_blk_compute_sampling_c2
   procedure :: compute_mpi_ab => nicas_blk_compute_mpi_ab
   procedure :: compute_interp_v => nicas_blk_compute_interp_v
   procedure :: compute_convol => nicas_blk_compute_convol
   procedure :: compute_convol_network => nicas_blk_compute_convol_network
   procedure :: compute_convol_distance => nicas_blk_compute_convol_distance
   procedure :: compute_convol_weights => nicas_blk_compute_convol_weights
   procedure :: compute_mpi_c => nicas_blk_compute_mpi_c
   procedure :: compute_internal_normalization => nicas_blk_compute_internal_normalization
   procedure :: compute_normalization => nicas_blk_compute_normalization
   procedure :: compute_grids => nicas_blk_compute_grids
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
   procedure :: test_adjoint => nicas_blk_test_adjoint
   procedure :: test_dirac => nicas_blk_test_dirac
end type nicas_blk_type

private
public :: nicas_blk_type

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
allocate(balldata%bd_to_c1u(balldata%nbd))
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
if (allocated(balldata%bd_to_c1u)) deallocate(balldata%bd_to_c1u)
if (allocated(balldata%bd_to_l1)) deallocate(balldata%bd_to_l1)
if (allocated(balldata%val)) deallocate(balldata%val)

end subroutine balldata_dealloc

!----------------------------------------------------------------------
! Subroutine: balldata_pack
! Purpose: pack data into balldata object
!----------------------------------------------------------------------
subroutine balldata_pack(balldata,mpl,nc1u,nl1,val)

implicit none

! Passed variables
class(balldata_type),intent(inout) :: balldata ! Ball data
type(mpl_type),intent(inout) :: mpl            ! MPI data
integer,intent(in) :: nc1u                     ! Horizontal box size
integer,intent(in) :: nl1                      ! Vertical box size
real(kind_real),intent(in) :: val(nc1u,nl1)    ! Box value

! Local variables
integer :: ibd,ic1u,il1

! Count non-missing values
balldata%nbd = count(mpl%msv%isnot(val))

! Allocation
call balldata%alloc

! Pack data
ibd = 0
do il1=1,nl1
   do ic1u=1,nc1u
      if (mpl%msv%isnot(val(ic1u,il1))) then
         ibd = ibd+1
         balldata%bd_to_c1u(ibd) = ic1u
         balldata%bd_to_l1(ibd) = il1
         balldata%val(ibd) = val(ic1u,il1)
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
if (allocated(nicas_blk%c1u_to_c1)) deallocate(nicas_blk%c1u_to_c1)
if (allocated(nicas_blk%c1_to_c1u)) deallocate(nicas_blk%c1_to_c1u)
if (allocated(nicas_blk%lon_c1u)) deallocate(nicas_blk%lon_c1u)
if (allocated(nicas_blk%lat_c1u)) deallocate(nicas_blk%lat_c1u)
if (allocated(nicas_blk%vunit_c1u)) deallocate(nicas_blk%vunit_c1u)
if (allocated(nicas_blk%gmask_c1u)) deallocate(nicas_blk%gmask_c1u)
if (allocated(nicas_blk%gmask_hor_c1u)) deallocate(nicas_blk%gmask_hor_c1u)
if (allocated(nicas_blk%lon_c1a)) deallocate(nicas_blk%lon_c1a)
if (allocated(nicas_blk%lat_c1a)) deallocate(nicas_blk%lat_c1a)
if (allocated(nicas_blk%gmask_c1a)) deallocate(nicas_blk%gmask_c1a)
if (allocated(nicas_blk%gmask_hor_c1a)) deallocate(nicas_blk%gmask_hor_c1a)
if (allocated(nicas_blk%c1a_to_c1)) deallocate(nicas_blk%c1a_to_c1)
if (allocated(nicas_blk%c1u_to_c1a)) deallocate(nicas_blk%c1u_to_c1a)
if (allocated(nicas_blk%c1b_to_c1u)) deallocate(nicas_blk%c1b_to_c1u)
if (allocated(nicas_blk%c1u_to_c1b)) deallocate(nicas_blk%c1u_to_c1b)
if (allocated(nicas_blk%c1u_to_c1bb)) deallocate(nicas_blk%c1u_to_c1bb)
if (allocated(nicas_blk%c1bb_to_c1u)) deallocate(nicas_blk%c1bb_to_c1u)
call nicas_blk%com_AU%dealloc
if (allocated(nicas_blk%c1_to_c0)) deallocate(nicas_blk%c1_to_c0)
if (allocated(nicas_blk%c1a_to_c0a)) deallocate(nicas_blk%c1a_to_c0a)
if (allocated(nicas_blk%nc2)) deallocate(nicas_blk%nc2)
if (allocated(nicas_blk%nc2u)) deallocate(nicas_blk%nc2u)
if (allocated(nicas_blk%gmask_c2u)) deallocate(nicas_blk%gmask_c2u)
if (allocated(nicas_blk%su_to_s)) deallocate(nicas_blk%su_to_s)
if (allocated(nicas_blk%su_to_sb)) deallocate(nicas_blk%su_to_sb)
if (allocated(nicas_blk%sb_to_su)) deallocate(nicas_blk%sb_to_su)
if (allocated(nicas_blk%lcheck_sa)) deallocate(nicas_blk%lcheck_sa)
if (allocated(nicas_blk%sa_to_su)) deallocate(nicas_blk%sa_to_su)
if (allocated(nicas_blk%lcheck_sb)) deallocate(nicas_blk%lcheck_sb)
if (allocated(nicas_blk%sbb_to_su)) deallocate(nicas_blk%sbb_to_su)
if (allocated(nicas_blk%sc_to_su)) deallocate(nicas_blk%sc_to_su)
if (allocated(nicas_blk%sa_to_sc_nor)) deallocate(nicas_blk%sa_to_sc_nor)
if (allocated(nicas_blk%sb_to_sc_nor)) deallocate(nicas_blk%sb_to_sc_nor)
if (allocated(nicas_blk%su_to_c1u)) deallocate(nicas_blk%su_to_c1u)
if (allocated(nicas_blk%su_to_l1)) deallocate(nicas_blk%su_to_l1)
if (allocated(nicas_blk%c1ul1_to_su)) deallocate(nicas_blk%c1ul1_to_su)
if (allocated(nicas_blk%c1bl1_to_sb)) deallocate(nicas_blk%c1bl1_to_sb)
if (allocated(nicas_blk%slev)) deallocate(nicas_blk%slev)
if (allocated(nicas_blk%vbot)) deallocate(nicas_blk%vbot)
if (allocated(nicas_blk%vtop)) deallocate(nicas_blk%vtop)
if (allocated(nicas_blk%l1_to_l0)) deallocate(nicas_blk%l1_to_l0)
if (allocated(nicas_blk%l0_to_l1)) deallocate(nicas_blk%l0_to_l1)
if (allocated(nicas_blk%rhs_avg)) deallocate(nicas_blk%rhs_avg)
if (allocated(nicas_blk%rh_c1u)) deallocate(nicas_blk%rh_c1u)
if (allocated(nicas_blk%rv_c1u)) deallocate(nicas_blk%rv_c1u)
if (allocated(nicas_blk%H11_c1u)) deallocate(nicas_blk%H11_c1u)
if (allocated(nicas_blk%H22_c1u)) deallocate(nicas_blk%H22_c1u)
if (allocated(nicas_blk%H33_c1u)) deallocate(nicas_blk%H33_c1u)
if (allocated(nicas_blk%H12_c1u)) deallocate(nicas_blk%H12_c1u)
if (allocated(nicas_blk%distnorm)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%distnorm(isbb)%dealloc
   end do
   deallocate(nicas_blk%distnorm)
end if
if (allocated(nicas_blk%Hcoef)) then
   do isbb=1,nicas_blk%nsbb
      call nicas_blk%Hcoef(isbb)%dealloc
   end do
   deallocate(nicas_blk%Hcoef)
end if
call nicas_blk%c_nor%dealloc
if (allocated(nicas_blk%inorm_nor)) deallocate(nicas_blk%inorm_nor)
call nicas_blk%com_AC_nor%dealloc
call nicas_blk%tree%dealloc
if (allocated(nicas_blk%smoother_norm)) deallocate(nicas_blk%smoother_norm)

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
integer :: il0,il1

! Release memory
call nicas_blk%partial_dealloc
if (allocated(nicas_blk%vlev)) deallocate(nicas_blk%vlev)
if (allocated(nicas_blk%sa_to_s)) deallocate(nicas_blk%sa_to_s)
if (allocated(nicas_blk%hash_sa)) deallocate(nicas_blk%hash_sa)
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
if (allocated(nicas_blk%sb_to_c1b)) deallocate(nicas_blk%sb_to_c1b)
if (allocated(nicas_blk%sb_to_l1)) deallocate(nicas_blk%sb_to_l1)
if (allocated(nicas_blk%inorm)) deallocate(nicas_blk%inorm)
if (allocated(nicas_blk%norm)) deallocate(nicas_blk%norm)
if (allocated(nicas_blk%coef_ens)) deallocate(nicas_blk%coef_ens)
call nicas_blk%com_AB%dealloc
call nicas_blk%com_AC%dealloc
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
subroutine nicas_blk_read(nicas_blk,mpl,nam,geom,bpar,ncid)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(bpar_type),intent(in) :: bpar               ! Block parameters
integer,intent(in) :: ncid                       ! NetCDF file

! Local variables
integer :: il0i,il1,il0
integer :: vlev_id,sb_to_c1b_id,sb_to_l1_id,sa_to_s_id,hash_sa_id,sa_to_sc_id,sb_to_sc_id,inorm_id,norm_id,coef_ens_id
integer :: vlev_int(geom%nl0)
character(len=1024),parameter :: subr = 'nicas_blk_read'

! Associate
associate(ib=>nicas_blk%ib)

! Get or check dimensions
call mpl%nc_dim_check(subr,ncid,'nl0',geom%nl0)
if (bpar%nicas_block(ib)) then
   nicas_blk%nc0a = mpl%nc_dim_inquire(subr,ncid,'nc0a')
   nicas_blk%nc1b = mpl%nc_dim_inquire(subr,ncid,'nc1b')
   nicas_blk%nl1 = mpl%nc_dim_inquire(subr,ncid,'nl1')
   nicas_blk%nsa = mpl%nc_dim_inquire(subr,ncid,'nsa')
   nicas_blk%nsb = mpl%nc_dim_inquire(subr,ncid,'nsb')
   nicas_blk%nsc = mpl%nc_dim_inquire(subr,ncid,'nsc')
   call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'ns',nicas_blk%ns))
end if

! Allocation
if (bpar%nicas_block(ib)) then
   allocate(nicas_blk%vlev(geom%nl0))
   if (.not.nicas_blk%smoother) allocate(nicas_blk%norm(nicas_blk%nc0a,geom%nl0))
   allocate(nicas_blk%coef_ens(nicas_blk%nc0a,geom%nl0))
   if (nicas_blk%nsa>0) then
      allocate(nicas_blk%sa_to_s(nicas_blk%nsa))
      allocate(nicas_blk%hash_sa(nicas_blk%nsa))
      allocate(nicas_blk%sa_to_sc(nicas_blk%nsa))
   end if
   if (nicas_blk%nsb>0) then
      allocate(nicas_blk%sb_to_c1b(nicas_blk%nsb))
      allocate(nicas_blk%sb_to_l1(nicas_blk%nsb))
      allocate(nicas_blk%sb_to_sc(nicas_blk%nsb))
   end if
   if (nicas_blk%nsc>0) then
      if (.not.nicas_blk%smoother) allocate(nicas_blk%inorm(nicas_blk%nsc))
   end if
   allocate(nicas_blk%h(geom%nl0i))
   allocate(nicas_blk%s(nicas_blk%nl1))
end if

! Get variable
if (bpar%nicas_block(ib)) then
   call mpl%ncerr(subr,nf90_inq_varid(ncid,'vlev',vlev_id))
   if (nicas_blk%nc0a>0) then
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'norm',norm_id))
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'coef_ens',coef_ens_id))
   end if
   if (nicas_blk%nsa>0) then
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'sa_to_s',sa_to_s_id))
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'hash_sa',hash_sa_id))
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'sa_to_sc',sa_to_sc_id))
   end if
   if (nicas_blk%nsb>0) then
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'sb_to_c1b',sb_to_c1b_id))
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'sb_to_l1',sb_to_l1_id))
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'sb_to_sc',sb_to_sc_id))
   end if
   if (nicas_blk%nsc>0) then
      call mpl%ncerr(subr,nf90_inq_varid(ncid,'inorm',inorm_id))
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
      if (.not.nicas_blk%smoother) call mpl%ncerr(subr,nf90_get_var(ncid,norm_id,nicas_blk%norm))
      call mpl%ncerr(subr,nf90_get_var(ncid,coef_ens_id,nicas_blk%coef_ens))
   end if
   if (nicas_blk%nsa>0) then
      call mpl%ncerr(subr,nf90_get_var(ncid,sa_to_s_id,nicas_blk%sa_to_s))
      call mpl%ncerr(subr,nf90_get_var(ncid,hash_sa_id,nicas_blk%hash_sa))
      call mpl%ncerr(subr,nf90_get_var(ncid,sa_to_sc_id,nicas_blk%sa_to_sc))
   end if
   if (nicas_blk%nsb>0) then
      call mpl%ncerr(subr,nf90_get_var(ncid,sb_to_c1b_id,nicas_blk%sb_to_c1b))
      call mpl%ncerr(subr,nf90_get_var(ncid,sb_to_l1_id,nicas_blk%sb_to_l1))
      call mpl%ncerr(subr,nf90_get_var(ncid,sb_to_sc_id,nicas_blk%sb_to_sc))
   end if
   if (nicas_blk%nsc>0) then
      if (.not.nicas_blk%smoother) call mpl%ncerr(subr,nf90_get_var(ncid,inorm_id,nicas_blk%inorm))
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

! Read main weight
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'wgt',nicas_blk%wgt))

! End associate
end associate

end subroutine nicas_blk_read

!----------------------------------------------------------------------
! Subroutine: nicas_blk_write
! Purpose: write
!----------------------------------------------------------------------
subroutine nicas_blk_write(nicas_blk,mpl,nam,geom,bpar,ncid)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(nam_type),intent(in) :: nam              ! Namelist
type(geom_type),intent(in) :: geom            ! Geometry
type(bpar_type),intent(in) :: bpar            ! Block parameters
integer,intent(in) :: ncid                    ! NetCDF file

! Local variables
integer :: il0i,il1,il0
integer :: nl0_id,nc0a_id,nc1b_id,nl1_id,nsa_id,nsb_id,nsc_id
integer :: vlev_id,sb_to_c1b_id,sb_to_l1_id,sa_to_s_id,hash_sa_id,sa_to_sc_id,sb_to_sc_id,inorm_id,norm_id,coef_ens_id
integer :: vlev_int(geom%nl0)
character(len=1024),parameter :: subr = 'nicas_blk_write'

! Associate
associate(ib=>nicas_blk%ib)

! Define dimensions
nl0_id = mpl%nc_dim_define_or_get(subr,ncid,'nl0',geom%nl0)
if (bpar%nicas_block(ib)) then
   if (nicas_blk%nc0a>0) nc0a_id = mpl%nc_dim_define_or_get(subr,ncid,'nc0a',nicas_blk%nc0a)
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0i',geom%nl0i))
   if (nicas_blk%nc1b>0) nc1b_id = mpl%nc_dim_define_or_get(subr,ncid,'nc1b',nicas_blk%nc1b)
   nl1_id = mpl%nc_dim_define_or_get(subr,ncid,'nl1',nicas_blk%nl1)
   if (nicas_blk%nsa>0) nsa_id = mpl%nc_dim_define_or_get(subr,ncid,'nsa',nicas_blk%nsa)
   if (nicas_blk%nsb>0) nsb_id = mpl%nc_dim_define_or_get(subr,ncid,'nsb',nicas_blk%nsb)
   if (nicas_blk%nsc>0) nsc_id = mpl%nc_dim_define_or_get(subr,ncid,'nsc',nicas_blk%nsc)
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'ns',nicas_blk%ns))
end if

! Define variables
if (bpar%nicas_block(ib)) then
   vlev_id = mpl%nc_var_define_or_get(subr,ncid,'vlev',nf90_int,(/nl0_id/))
   if (nicas_blk%nc0a>0) then
      if (.not.nicas_blk%smoother) norm_id = mpl%nc_var_define_or_get(subr,ncid,'norm',nc_kind_real,(/nc0a_id,nl0_id/))
      coef_ens_id = mpl%nc_var_define_or_get(subr,ncid,'coef_ens',nc_kind_real,(/nc0a_id,nl0_id/))
   end if
   if (nicas_blk%nsa>0) then
      sa_to_s_id = mpl%nc_var_define_or_get(subr,ncid,'sa_to_s',nf90_int,(/nsa_id/))
      hash_sa_id = mpl%nc_var_define_or_get(subr,ncid,'hash_sa',nc_kind_real,(/nsa_id/))
      sa_to_sc_id = mpl%nc_var_define_or_get(subr,ncid,'sa_to_sc',nf90_int,(/nsa_id/))
   end if
   if (nicas_blk%nsb>0) then
      sb_to_c1b_id = mpl%nc_var_define_or_get(subr,ncid,'sb_to_c1b',nf90_int,(/nsb_id/))
      sb_to_l1_id = mpl%nc_var_define_or_get(subr,ncid,'sb_to_l1',nf90_int,(/nsb_id/))
      sb_to_sc_id = mpl%nc_var_define_or_get(subr,ncid,'sb_to_sc',nf90_int,(/nsb_id/))
   end if
   if (nicas_blk%nsc>0) then
      if (.not.nicas_blk%smoother) inorm_id = mpl%nc_var_define_or_get(subr,ncid,'inorm',nc_kind_real,(/nsc_id/))
   end if
end if

! Write main weight
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'wgt',nicas_blk%wgt))

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
      if (.not.nicas_blk%smoother) call mpl%ncerr(subr,nf90_put_var(ncid,norm_id,nicas_blk%norm))
      call mpl%ncerr(subr,nf90_put_var(ncid,coef_ens_id,nicas_blk%coef_ens))
   end if
   if (nicas_blk%nsa>0) then
      call mpl%ncerr(subr,nf90_put_var(ncid,sa_to_s_id,nicas_blk%sa_to_s))
      call mpl%ncerr(subr,nf90_put_var(ncid,hash_sa_id,nicas_blk%hash_sa))
      call mpl%ncerr(subr,nf90_put_var(ncid,sa_to_sc_id,nicas_blk%sa_to_sc))
   end if
   if (nicas_blk%nsb>0) then
      call mpl%ncerr(subr,nf90_put_var(ncid,sb_to_c1b_id,nicas_blk%sb_to_c1b))
      call mpl%ncerr(subr,nf90_put_var(ncid,sb_to_l1_id,nicas_blk%sb_to_l1))
      call mpl%ncerr(subr,nf90_put_var(ncid,sb_to_sc_id,nicas_blk%sb_to_sc))
   end if
   if (nicas_blk%nsc>0) then
      if (.not.nicas_blk%smoother) call mpl%ncerr(subr,nf90_put_var(ncid,inorm_id,nicas_blk%inorm))
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

! End associate
end associate

end subroutine nicas_blk_write

!----------------------------------------------------------------------
! Subroutine: nicas_blk_write_grids
! Purpose: write NICAS grids
!----------------------------------------------------------------------
subroutine nicas_blk_write_grids(nicas_blk,mpl,ncid)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl           ! MPI data
integer,intent(in) :: ncid                    ! NetCDF file

! Local variables
integer :: nsa_id,nsb_id,nsc_id,lon_sa_id,lat_sa_id,lev_sa_id,lon_sb_id,lat_sb_id,lev_sb_id,lon_sc_id,lat_sc_id,lev_sc_id
character(len=1024),parameter :: subr = 'nicas_blk_write_grids'

! Define dimensions
if (nicas_blk%nsa>0) nsa_id = mpl%nc_dim_define_or_get(subr,ncid,'nsa',nicas_blk%nsa)
if (nicas_blk%nsb>0) nsb_id = mpl%nc_dim_define_or_get(subr,ncid,'nsb',nicas_blk%nsb)
if (nicas_blk%nsc>0) nsc_id = mpl%nc_dim_define_or_get(subr,ncid,'nsc',nicas_blk%nsc)

! Define variables
if (nicas_blk%nsa>0) then
   lon_sa_id = mpl%nc_var_define_or_get(subr,ncid,'lon_sa',nc_kind_real,(/nsa_id/))
   lat_sa_id = mpl%nc_var_define_or_get(subr,ncid,'lat_sa',nc_kind_real,(/nsa_id/))
   lev_sa_id = mpl%nc_var_define_or_get(subr,ncid,'lev_sa',nf90_int,(/nsa_id/))
end if
if (nicas_blk%nsb>0) then
   lon_sb_id = mpl%nc_var_define_or_get(subr,ncid,'lon_sb',nc_kind_real,(/nsb_id/))
   lat_sb_id = mpl%nc_var_define_or_get(subr,ncid,'lat_sb',nc_kind_real,(/nsb_id/))
   lev_sb_id = mpl%nc_var_define_or_get(subr,ncid,'lev_sb',nf90_int,(/nsb_id/))
end if
if (nicas_blk%nsc>0) then
   lon_sc_id = mpl%nc_var_define_or_get(subr,ncid,'lon_sc',nc_kind_real,(/nsc_id/))
   lat_sc_id = mpl%nc_var_define_or_get(subr,ncid,'lat_sc',nc_kind_real,(/nsc_id/))
   lev_sc_id = mpl%nc_var_define_or_get(subr,ncid,'lev_sc',nf90_int,(/nsc_id/))
end if

! Write variables
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

end subroutine nicas_blk_write_grids

!----------------------------------------------------------------------
! Subroutine: nicas_blk_buffer_size
! Purpose: buffer size
!----------------------------------------------------------------------
subroutine nicas_blk_buffer_size(nicas_blk,mpl,nam,geom,bpar,nbufi,nbufr,nbufl)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(nam_type),intent(in) :: nam              ! Namelist
type(geom_type),intent(in) :: geom            ! Geometry
type(bpar_type),intent(in) :: bpar            ! Block parameters
integer,intent(out) :: nbufi                  ! Buffer size (integer)
integer,intent(out) :: nbufr                  ! Buffer size (real)
integer,intent(out) :: nbufl                  ! Buffer size (logical)

! Local variables
integer :: nnbufi,nnbufr
integer :: il0i,il1,il0

! Associate
associate(ib=>nicas_blk%ib)

! Initialization
nbufi = 3
nbufr = 0
nbufl = 0

! Define buffer size
if (bpar%nicas_block(ib)) nbufi = nbufi+6
nbufr = nbufr+1
if (bpar%nicas_block(ib)) then
   nbufi = nbufi+2*nicas_blk%nsa+3*nicas_blk%nsb
   nbufr = nbufr+nicas_blk%nsa+nicas_blk%nc0a*geom%nl0
   if (.not.nicas_blk%smoother) nbufr = nbufr+nicas_blk%nc0a*geom%nl0+nicas_blk%nsc
   nbufl = nbufl+geom%nl0
   if (nam%write_grids) then
      nbufi = nbufi+nicas_blk%nsa+nicas_blk%nsb+nicas_blk%nsc
      nbufr = nbufr+2*(nicas_blk%nsa+nicas_blk%nsb+nicas_blk%nsc)
   end if
end if

! Add communications and linear operators
if (bpar%nicas_block(ib)) then
   call nicas_blk%com_AB%buffer_size(mpl,nnbufi)
   nbufi = nbufi+nnbufi
   call nicas_blk%com_AC%buffer_size(mpl,nnbufi)
   nbufi = nbufi+nnbufi
   call nicas_blk%c%buffer_size(nnbufi,nnbufr)
   nbufi = nbufi+nnbufi
   nbufr = nbufr+nnbufr
   do il0i=1,geom%nl0i
      call nicas_blk%h(il0i)%buffer_size(nnbufi,nnbufr)
      nbufi = nbufi+nnbufi
      nbufr = nbufr+nnbufr
   end do
   call nicas_blk%v%buffer_size(nnbufi,nnbufr)
   nbufi = nbufi+nnbufi
   nbufr = nbufr+nnbufr
   do il1=1,nicas_blk%nl1
      call nicas_blk%s(il1)%buffer_size(nnbufi,nnbufr)
      nbufi = nbufi+nnbufi
      nbufr = nbufr+nnbufr
   end do
end if

! End associate
end associate

end subroutine nicas_blk_buffer_size

!----------------------------------------------------------------------
! Subroutine: nicas_blk_serialize
! Purpose: serialize
!----------------------------------------------------------------------
subroutine nicas_blk_serialize(nicas_blk,mpl,nam,geom,bpar,nbufi,nbufr,nbufl,bufi,bufr,bufl)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(nam_type),intent(in) :: nam              ! Namelist
type(geom_type),intent(in) :: geom            ! Geometry
type(bpar_type),intent(in) :: bpar            ! Block parameters
integer,intent(in) :: nbufi                   ! Buffer size (integer)
integer,intent(in) :: nbufr                   ! Buffer size (real)
integer,intent(in) :: nbufl                   ! Buffer size (logical)
integer,intent(out) :: bufi(nbufi)            ! Buffer (integer)
real(kind_real),intent(out) :: bufr(nbufr)    ! Buffer (real)
logical,intent(out) :: bufl(nbufl)            ! Buffer (logical)

! Local variables
integer :: ibufi,ibufr,ibufl,nnbufi,nnbufr
integer :: il0i,il1,il0
logical,allocatable :: mask_c0a(:,:)
character(len=1024),parameter :: subr = 'nicas_blk_serialize'

! Associate
associate(ib=>nicas_blk%ib)

! Initialization
ibufi = 0
ibufr = 0
ibufl = 0

! Dimensions
bufi(ibufi+1) = nbufi
ibufi = ibufi+1
bufi(ibufi+1) = nbufr
ibufi = ibufi+1
bufi(ibufi+1) = nbufl
ibufi = ibufi+1
if (bpar%nicas_block(ib)) then
   bufi(ibufi+1) = nicas_blk%nc0a
   ibufi = ibufi+1
   bufi(ibufi+1) = nicas_blk%nc1b
   ibufi = ibufi+1
   bufi(ibufi+1) = nicas_blk%nl1
   ibufi = ibufi+1
   bufi(ibufi+1) = nicas_blk%nsa
   ibufi = ibufi+1
   bufi(ibufi+1) = nicas_blk%nsb
   ibufi = ibufi+1
   bufi(ibufi+1) = nicas_blk%nsc
   ibufi = ibufi+1
end if

if (bpar%nicas_block(ib)) then
   if (nicas_blk%nc0a>0) then
      ! Allocation
      allocate(mask_c0a(nicas_blk%nc0a,geom%nl0))

      ! Initialization
      mask_c0a = .true.
   end if
end if

! Data
bufr(ibufr+1) = nicas_blk%wgt
ibufr = ibufr+1
if (bpar%nicas_block(ib)) then
   bufl(ibufl+1:ibufl+geom%nl0) = nicas_blk%vlev
   ibufl = ibufl+geom%nl0
   if (nicas_blk%nc0a>0) then
      if (.not.nicas_blk%smoother) then
         bufr(ibufr+1:ibufr+nicas_blk%nc0a*geom%nl0) = pack(nicas_blk%norm,mask_c0a)
         ibufr = ibufr+nicas_blk%nc0a*geom%nl0
      end if
      bufr(ibufr+1:ibufr+nicas_blk%nc0a*geom%nl0) = pack(nicas_blk%coef_ens,mask_c0a)
      ibufr = ibufr+nicas_blk%nc0a*geom%nl0
   end if
   if (nicas_blk%nsa>0) then
      bufi(ibufi+1:ibufi+nicas_blk%nsa) = nicas_blk%sa_to_s
      ibufi = ibufi+nicas_blk%nsa
      bufr(ibufr+1:ibufr+nicas_blk%nsa) = nicas_blk%hash_sa
      ibufr = ibufr+nicas_blk%nsa
      bufi(ibufi+1:ibufi+nicas_blk%nsa) = nicas_blk%sa_to_sc
      ibufi = ibufi+nicas_blk%nsa
   end if
   if (nicas_blk%nsb>0) then
      bufi(ibufi+1:ibufi+nicas_blk%nsb) = nicas_blk%sb_to_c1b
      ibufi = ibufi+nicas_blk%nsb
      bufi(ibufi+1:ibufi+nicas_blk%nsb) = nicas_blk%sb_to_l1
      ibufi = ibufi+nicas_blk%nsb
      bufi(ibufi+1:ibufi+nicas_blk%nsb) = nicas_blk%sb_to_sc
      ibufi = ibufi+nicas_blk%nsb
   end if
   if (nicas_blk%nsc>0) then
      if (.not.nicas_blk%smoother) then
         bufr(ibufr+1:ibufr+nicas_blk%nsc) = nicas_blk%inorm
         ibufr = ibufr+nicas_blk%nsc
      end if
   end if
   if (nam%write_grids) then
      if (nicas_blk%nsa>0) then
         bufr(ibufr+1:ibufr+nicas_blk%nsa) = nicas_blk%lon_sa
         ibufr = ibufr+nicas_blk%nsa
         bufr(ibufr+1:ibufr+nicas_blk%nsa) = nicas_blk%lat_sa
         ibufr = ibufr+nicas_blk%nsa
         bufi(ibufi+1:ibufi+nicas_blk%nsa) = nicas_blk%lev_sa
         ibufi = ibufi+nicas_blk%nsa
      end if
      if (nicas_blk%nsb>0) then
         bufr(ibufr+1:ibufr+nicas_blk%nsb) = nicas_blk%lon_sb
         ibufr = ibufr+nicas_blk%nsb
         bufr(ibufr+1:ibufr+nicas_blk%nsb) = nicas_blk%lat_sb
         ibufr = ibufr+nicas_blk%nsb
         bufi(ibufi+1:ibufi+nicas_blk%nsb) = nicas_blk%lev_sb
         ibufi = ibufi+nicas_blk%nsb
      end if
      if (nicas_blk%nsc>0) then
         bufr(ibufr+1:ibufr+nicas_blk%nsc) = nicas_blk%lon_sc
         ibufr = ibufr+nicas_blk%nsc
         bufr(ibufr+1:ibufr+nicas_blk%nsc) = nicas_blk%lat_sc
         ibufr = ibufr+nicas_blk%nsc
         bufi(ibufi+1:ibufi+nicas_blk%nsc) = nicas_blk%lev_sc
         ibufi = ibufi+nicas_blk%nsc
      end if
   end if
end if

! Release memory
if (bpar%nicas_block(ib)) then
   if (nicas_blk%nc0a>0) deallocate(mask_c0a)
end if
      
! Communications and linear operators
if (bpar%nicas_block(ib)) then
   call nicas_blk%com_AB%buffer_size(mpl,nnbufi)
   call nicas_blk%com_AB%serialize(mpl,nnbufi,bufi(ibufi+1:ibufi+nnbufi))
   ibufi = ibufi+nnbufi
   call nicas_blk%com_AC%buffer_size(mpl,nnbufi)
   call nicas_blk%com_AC%serialize(mpl,nnbufi,bufi(ibufi+1:ibufi+nnbufi))
   ibufi = ibufi+nnbufi
   call nicas_blk%c%buffer_size(nnbufi,nnbufr)
   call nicas_blk%c%serialize(mpl,nnbufi,nnbufr,bufi(ibufi+1:ibufi+nnbufi),bufr(ibufr+1:ibufr+nnbufr))
   ibufi = ibufi+nnbufi
   ibufr = ibufr+nnbufr
   do il0i=1,geom%nl0i
      call nicas_blk%h(il0i)%buffer_size(nnbufi,nnbufr)
      call nicas_blk%h(il0i)%serialize(mpl,nnbufi,nnbufr,bufi(ibufi+1:ibufi+nnbufi),bufr(ibufr+1:ibufr+nnbufr))
      ibufi = ibufi+nnbufi
      ibufr = ibufr+nnbufr
   end do
   call nicas_blk%v%buffer_size(nnbufi,nnbufr)
   call nicas_blk%v%serialize(mpl,nnbufi,nnbufr,bufi(ibufi+1:ibufi+nnbufi),bufr(ibufr+1:ibufr+nnbufr))
   ibufi = ibufi+nnbufi
   ibufr = ibufr+nnbufr
   do il1=1,nicas_blk%nl1
      call nicas_blk%s(il1)%buffer_size(nnbufi,nnbufr)
      call nicas_blk%s(il1)%serialize(mpl,nnbufi,nnbufr,bufi(ibufi+1:ibufi+nnbufi),bufr(ibufr+1:ibufr+nnbufr))
      ibufi = ibufi+nnbufi
      ibufr = ibufr+nnbufr
   end do
end if

! Check
if (ibufi/=nbufi) call mpl%abort(subr,'inconsistent final offset/buffer size (integer)')
if (ibufr/=nbufr) call mpl%abort(subr,'inconsistent final offset/buffer size (real)')
if (ibufl/=nbufl) call mpl%abort(subr,'inconsistent final offset/buffer size (logical)')

! End associate
end associate

end subroutine nicas_blk_serialize

!----------------------------------------------------------------------
! Subroutine: nicas_blk_deserialize
! Purpose: deserialize
!----------------------------------------------------------------------
subroutine nicas_blk_deserialize(nicas_blk,mpl,nam,geom,bpar,nbufi,nbufr,nbufl,bufi,bufr,bufl)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(bpar_type),intent(in) :: bpar               ! Block parameters
integer,intent(in) :: nbufi                      ! Buffer size (integer)
integer,intent(in) :: nbufr                      ! Buffer size (real)
integer,intent(in) :: nbufl                      ! Buffer size (logical)
integer,intent(in) :: bufi(nbufi)                ! Buffer (integer)
real(kind_real),intent(in) :: bufr(nbufr)        ! Buffer (real)
logical,intent(in) :: bufl(nbufl)                ! Buffer (logical)

! Local variables
integer :: ibufi,ibufr,ibufl,nnbufi,nnbufr
integer :: il0i,il1,il0
logical,allocatable :: mask_c0a(:,:)
character(len=1024),parameter :: subr = 'nicas_blk_deserialize'

! Associate
associate(ib=>nicas_blk%ib)

! Initialization
ibufi = 0
ibufr = 0
ibufl = 0

! Check
if (bufi(ibufi+1)/=nbufi) call mpl%abort(subr,'inconsistent initial value/buffer size (integer)')
ibufi = ibufi+1
if (bufi(ibufi+1)/=nbufr) call mpl%abort(subr,'inconsistent initial value/buffer size (real)')
ibufi = ibufi+1
if (bufi(ibufi+1)/=nbufl) call mpl%abort(subr,'inconsistent initial value/buffer size (logical)')
ibufi = ibufi+1

! Dimensions
if (bpar%nicas_block(ib)) then
   nicas_blk%nc0a = bufi(ibufi+1)
   ibufi = ibufi+1
   nicas_blk%nc1b = bufi(ibufi+1)
   ibufi = ibufi+1
   nicas_blk%nl1 = bufi(ibufi+1)
   ibufi = ibufi+1
   nicas_blk%nsa = bufi(ibufi+1)
   ibufi = ibufi+1
   nicas_blk%nsb = bufi(ibufi+1)
   ibufi = ibufi+1
   nicas_blk%nsc = bufi(ibufi+1)
   ibufi = ibufi+1
end if

! Allocation
if (bpar%nicas_block(ib)) then
   allocate(nicas_blk%vlev(geom%nl0))
   if (.not.nicas_blk%smoother) allocate(nicas_blk%norm(nicas_blk%nc0a,geom%nl0))
   allocate(nicas_blk%coef_ens(nicas_blk%nc0a,geom%nl0))
   if (nicas_blk%nsa>0) then
      allocate(nicas_blk%sa_to_s(nicas_blk%nsa))
      allocate(nicas_blk%hash_sa(nicas_blk%nsa))
      allocate(nicas_blk%sa_to_sc(nicas_blk%nsa))
   end if
   if (nicas_blk%nsb>0) then
      allocate(nicas_blk%sb_to_c1b(nicas_blk%nsb))
      allocate(nicas_blk%sb_to_l1(nicas_blk%nsb))
      allocate(nicas_blk%sb_to_sc(nicas_blk%nsb))
   end if
   if (nicas_blk%nsc>0) then
      if (.not.nicas_blk%smoother) allocate(nicas_blk%inorm(nicas_blk%nsc))
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

if (bpar%nicas_block(ib)) then
   if (nicas_blk%nc0a>0) then
      ! Allocation
      allocate(mask_c0a(nicas_blk%nc0a,geom%nl0))

      ! Initialization
      if (bpar%nicas_block(ib)) mask_c0a = .true.
   end if
end if

! Data
nicas_blk%wgt = bufr(ibufr+1)
ibufr = ibufr+1
if (bpar%nicas_block(ib)) then
   nicas_blk%vlev = bufl(ibufl+1:ibufl+geom%nl0)
   ibufl = ibufl+geom%nl0
   if (nicas_blk%nc0a>0) then
      if (.not.nicas_blk%smoother) then
         nicas_blk%norm = unpack(bufr(ibufr+1:ibufr+nicas_blk%nc0a*geom%nl0),mask_c0a,nicas_blk%norm)
         ibufr = ibufr+nicas_blk%nc0a*geom%nl0
      end if
      nicas_blk%coef_ens = unpack(bufr(ibufr+1:ibufr+nicas_blk%nc0a*geom%nl0),mask_c0a,nicas_blk%coef_ens)
      ibufr = ibufr+nicas_blk%nc0a*geom%nl0
   end if
   if (nicas_blk%nsa>0) then
      nicas_blk%sa_to_s = bufi(ibufi+1:ibufi+nicas_blk%nsa)
      ibufi = ibufi+nicas_blk%nsa
      nicas_blk%hash_sa = bufr(ibufr+1:ibufr+nicas_blk%nsa)
      ibufr = ibufr+nicas_blk%nsa
      nicas_blk%sa_to_sc = bufi(ibufi+1:ibufi+nicas_blk%nsa)
      ibufi = ibufi+nicas_blk%nsa
   end if
   if (nicas_blk%nsb>0) then
      nicas_blk%sb_to_c1b = bufi(ibufi+1:ibufi+nicas_blk%nsb)
      ibufi = ibufi+nicas_blk%nsb
      nicas_blk%sb_to_l1 = bufi(ibufi+1:ibufi+nicas_blk%nsb)
      ibufi = ibufi+nicas_blk%nsb
      nicas_blk%sb_to_sc = bufi(ibufi+1:ibufi+nicas_blk%nsb)
      ibufi = ibufi+nicas_blk%nsb
   end if
   if (nicas_blk%nsc>0) then
      if (.not.nicas_blk%smoother) then
         nicas_blk%inorm = bufr(ibufr+1:ibufr+nicas_blk%nsc)
         ibufr = ibufr+nicas_blk%nsc
      end if
   end if
   if (nam%write_grids) then
      if (nicas_blk%nsa>0) then
         nicas_blk%lon_sa = bufr(ibufr+1:ibufr+nicas_blk%nsa)
         ibufr = ibufr+nicas_blk%nsa
         nicas_blk%lat_sa = bufr(ibufr+1:ibufr+nicas_blk%nsa)
         ibufr = ibufr+nicas_blk%nsa
         nicas_blk%lev_sa = bufi(ibufi+1:ibufi+nicas_blk%nsa)
         ibufi = ibufi+nicas_blk%nsa
      end if
      if (nicas_blk%nsb>0) then
         nicas_blk%lon_sb = bufr(ibufr+1:ibufr+nicas_blk%nsb)
         ibufr = ibufr+nicas_blk%nsb
         nicas_blk%lat_sb = bufr(ibufr+1:ibufr+nicas_blk%nsb)
         ibufr = ibufr+nicas_blk%nsb
         nicas_blk%lev_sb = bufi(ibufi+1:ibufi+nicas_blk%nsb)
         ibufi = ibufi+nicas_blk%nsb
      end if
      if (nicas_blk%nsc>0) then
         nicas_blk%lon_sc = bufr(ibufr+1:ibufr+nicas_blk%nsc)
         ibufr = ibufr+nicas_blk%nsc
         nicas_blk%lat_sc = bufr(ibufr+1:ibufr+nicas_blk%nsc)
         ibufr = ibufr+nicas_blk%nsc
         nicas_blk%lev_sc = bufi(ibufi+1:ibufi+nicas_blk%nsc)
         ibufi = ibufi+nicas_blk%nsc
      end if
   end if
end if

! Release memory
if (bpar%nicas_block(ib)) then
   if (nicas_blk%nc0a>0) deallocate(mask_c0a)
end if
      
! Communications and linear operators
if (bpar%nicas_block(ib)) then
   nicas_blk%com_AB%prefix = 'com_AB'
   nnbufi = bufi(ibufi+1)
   call nicas_blk%com_AB%deserialize(mpl,nnbufi,bufi(ibufi+1:ibufi+nnbufi))
   ibufi = ibufi+nnbufi
   nicas_blk%com_AC%prefix = 'com_AC'
   nnbufi = bufi(ibufi+1)
   call nicas_blk%com_AC%deserialize(mpl,nnbufi,bufi(ibufi+1:ibufi+nnbufi))
   ibufi = ibufi+nnbufi
   nicas_blk%c%prefix = 'c'
   nnbufi = bufi(ibufi+1)
   nnbufr = bufi(ibufi+2)
   call nicas_blk%c%deserialize(mpl,nnbufi,nnbufr,bufi(ibufi+1:ibufi+nnbufi),bufr(ibufr+1:ibufr+nnbufr))
   ibufi = ibufi+nnbufi
   ibufr = ibufr+nnbufr
   do il0i=1,geom%nl0i
      write(nicas_blk%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
      nnbufi = bufi(ibufi+1)
      nnbufr = bufi(ibufi+2)
      call nicas_blk%h(il0i)%deserialize(mpl,nnbufi,nnbufr,bufi(ibufi+1:ibufi+nnbufi),bufr(ibufr+1:ibufr+nnbufr))
      ibufi = ibufi+nnbufi
      ibufr = ibufr+nnbufr
   end do
   nicas_blk%v%prefix = 'v'
   nnbufi = bufi(ibufi+1)
   nnbufr = bufi(ibufi+2)
   call nicas_blk%v%deserialize(mpl,nnbufi,nnbufr,bufi(ibufi+1:ibufi+nnbufi),bufr(ibufr+1:ibufr+nnbufr))
   ibufi = ibufi+nnbufi
   ibufr = ibufr+nnbufr
   do il1=1,nicas_blk%nl1
      write(nicas_blk%s(il1)%prefix,'(a,i3.3)') 's_',il1
      nnbufi = bufi(ibufi+1)
      nnbufr = bufi(ibufi+2)
      call nicas_blk%s(il1)%deserialize(mpl,nnbufi,nnbufr,bufi(ibufi+1:ibufi+nnbufi),bufr(ibufr+1:ibufr+nnbufr))
      ibufi = ibufi+nnbufi
      ibufr = ibufr+nnbufr
   end do
end if

! End associate
end associate

end subroutine nicas_blk_deserialize

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_parameters
! Purpose: compute NICAS parameters
!----------------------------------------------------------------------
subroutine nicas_blk_compute_parameters(nicas_blk,mpl,rng,nam,geom,cmat_blk,sqrt_rescaling)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(rng_type),intent(inout) :: rng              ! Random number generator
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
type(cmat_blk_type),intent(in) :: cmat_blk       ! C matrix data block
logical,intent(in),optional :: sqrt_rescaling    ! Square-root rescaling flag

! Local variables
integer :: il0i,il1
character(len=1024),parameter :: subr = 'nicas_blk_compute_parameters'

! Set square-root rescaling flag
nicas_blk%sqrt_rescaling = .true.
if (present(sqrt_rescaling)) nicas_blk%sqrt_rescaling = sqrt_rescaling

! Copy grid hash and size of subset Sc0 on halo A 
nicas_blk%grid_hash = geom%grid_hash
nicas_blk%nc0a = geom%nc0a

! Check horizontal flag
if (nicas_blk%horizontal) then
   if ((trim(nicas_blk%subsamp)=='hv').or.(trim(nicas_blk%subsamp)=='vh').or.(trim(nicas_blk%subsamp)=='hvh')) &
 & call mpl%abort(subr,'horizontal flag requires h subsamp')
end if

! Compute adaptive sampling, subset Sc1
write(mpl%info,'(a7,a)') '','Compute adaptive sampling, subset Sc1'
if (nicas_blk%verbosity) call mpl%flush
call nicas_blk%compute_sampling_c1(mpl,rng,nam,geom,cmat_blk)

! Compute adaptive sampling, vertical
write(mpl%info,'(a7,a)') '','Compute adaptive sampling, vertical'
if (nicas_blk%verbosity) call mpl%flush
call nicas_blk%compute_sampling_v(mpl,nam,geom,cmat_blk)

! Compute MPI distribution, halos A
write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo A'
if (nicas_blk%verbosity) call mpl%flush
call nicas_blk%compute_mpi_a(mpl,geom)

! Compute adaptive sampling, subset Sc2
write(mpl%info,'(a7,a)') '','Compute adaptive sampling, subset Sc2'
if (nicas_blk%verbosity) call mpl%flush
call nicas_blk%compute_sampling_c2(mpl,rng,nam,geom,cmat_blk)

! Compute MPI distribution, halos A-B
write(mpl%info,'(a7,a)') '','Compute MPI distribution, halos A-B'
if (nicas_blk%verbosity) call mpl%flush
call nicas_blk%compute_mpi_ab(mpl,rng,nam,geom)

! Compute vertical interpolation data
write(mpl%info,'(a7,a)') '','Compute vertical interpolation data'
if (nicas_blk%verbosity) call mpl%flush
call nicas_blk%compute_interp_v(geom)

! Compute convolution data
write(mpl%info,'(a7,a)') '','Compute convolution data'
if (nicas_blk%verbosity) call mpl%flush
call nicas_blk%compute_convol(mpl,rng,nam,geom,cmat_blk)

! Compute MPI distribution, halo C
write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo C'
if (nicas_blk%verbosity) call mpl%flush
call nicas_blk%compute_mpi_c(mpl)

if (.not.nicas_blk%smoother) then
   ! Compute internal normalization
   write(mpl%info,'(a7,a)') '','Compute internal normalization'
   if (nicas_blk%verbosity) call mpl%flush
   call nicas_blk%compute_internal_normalization(mpl)

   ! Compute normalization
   write(mpl%info,'(a7,a)') '','Compute normalization'
   if (nicas_blk%verbosity) call mpl%flush
   call nicas_blk%compute_normalization(mpl,nam,geom)
end if

if (nam%write_grids) then
   ! Compute grids coordinates
   write(mpl%info,'(a7,a)') '','Compute grids coordinates'
   if (nicas_blk%verbosity) call mpl%flush
   call nicas_blk%compute_grids(nam)
end if

! Print results
write(mpl%info,'(a7,a,i6)') '','Parameters for processor #',mpl%myproc
if (nicas_blk%verbosity) call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0 =        ',geom%nc0
if (nicas_blk%verbosity) call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0a =       ',nicas_blk%nc0a
if (nicas_blk%verbosity) call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nl0 =        ',geom%nl0
if (nicas_blk%verbosity) call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1 =        ',nicas_blk%nc1
if (nicas_blk%verbosity) call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1a =       ',nicas_blk%nc1a
if (nicas_blk%verbosity) call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1b =       ',nicas_blk%nc1b
if (nicas_blk%verbosity) call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nl1 =        ',nicas_blk%nl1
do il1=1,nicas_blk%nl1
   write(mpl%info,'(a10,a,i3,a,i8)') '','nc2(',il1,') =   ',nicas_blk%nc2(il1)
   if (nicas_blk%verbosity) call mpl%flush
end do
write(mpl%info,'(a10,a,i8)') '','ns =         ',nicas_blk%ns
if (nicas_blk%verbosity) call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nsa =        ',nicas_blk%nsa
if (nicas_blk%verbosity) call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nsb =        ',nicas_blk%nsb
if (nicas_blk%verbosity) call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nsc =        ',nicas_blk%nsc
do il0i=1,geom%nl0i
   write(mpl%info,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',nicas_blk%h(il0i)%n_s
   if (nicas_blk%verbosity) call mpl%flush
end do
write(mpl%info,'(a10,a,i8)') '','v%n_s =      ',nicas_blk%v%n_s
if (nicas_blk%verbosity) call mpl%flush
do il1=1,nicas_blk%nl1
   write(mpl%info,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',nicas_blk%s(il1)%n_s
   if (nicas_blk%verbosity) call mpl%flush
end do
write(mpl%info,'(a10,a,i9)') '','c%n_s =     ',nicas_blk%c%n_s
if (nicas_blk%verbosity) call mpl%flush
if (.not.nicas_blk%smoother) then
   write(mpl%info,'(a10,a,i9)') '','c_nor%n_s = ',nicas_blk%c_nor%n_s
   if (nicas_blk%verbosity) call mpl%flush
end if

! Release memory (partial)
call nicas_blk%partial_dealloc

end subroutine nicas_blk_compute_parameters

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_parameters_horizontal_smoother
! Purpose: compute NICAS parameters for a horizontal smoother
!----------------------------------------------------------------------
subroutine nicas_blk_compute_parameters_horizontal_smoother(nicas_blk,mpl,rng,nam,geom,rhflt)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(rng_type),intent(inout) :: rng              ! Random number generator
type(nam_type),intent(in) :: nam                 ! Namelist
type(geom_type),intent(in) :: geom               ! Geometry
real(kind_real),intent(in) :: rhflt(geom%nl0)    ! Horizontal support radius profile

! Local variables
integer :: il0
type(cmat_blk_type) :: cmat_blk
type(nam_type) :: nam_smoother

! Initialize local namelist
call nam_smoother%init(mpl%nproc)

! Set local namelist parameters
nam_smoother%levs = nam%levs
if (nam%ntry==nam_smoother%ntry) then
   nam_smoother%ntry = 10
else
   nam_smoother%ntry = nam%ntry
end if
if (eq(nam%resol,nam_smoother%resol)) then
   nam_smoother%resol = 8.0
else
   nam_smoother%resol = nam%resol
end if
if (nam%fast_sampling.eqv.nam_smoother%fast_sampling) then
   nam_smoother%fast_sampling = nam%fast_sampling
else
   nam_smoother%fast_sampling = .false.
end if

! Local cmat_blk allocation
allocate(cmat_blk%coef_ens(geom%nc0a,geom%nl0))
allocate(cmat_blk%coef_sta(geom%nc0a,geom%nl0))
allocate(cmat_blk%rh(geom%nc0a,geom%nl0))
allocate(cmat_blk%rhs(geom%nc0a,geom%nl0))

! Local cmat_blk initialization
cmat_blk%anisotropic = .false.
cmat_blk%coef_ens = 1.0
cmat_blk%coef_sta = 0.0
do il0=1,geom%nl0
   cmat_blk%rh(:,il0) = rhflt(il0)
end do
cmat_blk%rhs = cmat_blk%rh
cmat_blk%wgt = 1.0

! NICAS block initialization
nicas_blk%verbosity = .false.
nicas_blk%smoother = .true.
nicas_blk%horizontal = .true.
nicas_blk%mpicom = 1
nicas_blk%lsqrt = 0
nicas_blk%subsamp = 'h'

! Compute parameters
call nicas_blk%compute_parameters(mpl,rng,nam_smoother,geom,cmat_blk)

! Release memory
call cmat_blk%dealloc

end subroutine nicas_blk_compute_parameters_horizontal_smoother

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
integer :: il0,ic0a,ic0,ic1,ic0u,ic1u,iproc
real(kind_real) :: rhs_sum(geom%nl0),rvs_sum(geom%nl0),rvs_avg(geom%nl0),norm(geom%nl0)
real(kind_real) :: rhs_minavg,rhs_min_norm,rhs_min_norm_tot
real(kind_real) :: rhs_min(geom%nc0a)
logical :: mask_hor_c0a(geom%nc0a)
character(len=1024),parameter :: subr = 'nicas_blk_compute_sampling_c1'

! Check subsampling method
if ((geom%nl0i>1).and.(trim(nicas_blk%subsamp)/='h')) call mpl%abort(subr,'subsamp = h required for variable mask')

! Allocation
allocate(nicas_blk%rhs_avg(geom%nl0))
allocate(nicas_blk%vlev(geom%nl0))

! Reset random numbers seed
if (trim(nam%strategy)=='specific_multivariate') call rng%reseed(mpl)

! Compute support radii
norm = 1.0/real(geom%nc0_gmask(1:geom%nl0),kind_real)
rhs_sum = sum(cmat_blk%rhs,dim=1,mask=geom%gmask_c0a)
call mpl%f_comm%allreduce(rhs_sum,nicas_blk%rhs_avg,fckit_mpi_sum())
nicas_blk%rhs_avg = nicas_blk%rhs_avg*norm
if (nicas_blk%horizontal) then
   rvs_avg = 0.0
else
   rvs_sum = sum(cmat_blk%rvs,dim=1,mask=geom%gmask_c0a)
   call mpl%f_comm%allreduce(rvs_sum,rvs_avg,fckit_mpi_sum())
   rvs_avg = rvs_avg*norm
end if
write(mpl%info,'(a10,a)') '','Average support radii (H/V): '
if (nicas_blk%verbosity) call mpl%flush
do il0=1,geom%nl0
   nicas_blk%vlev(il0) = (nicas_blk%rhs_avg(il0)>0.0).or.(rvs_avg(il0)>0.0)
   if (nicas_blk%vlev(il0)) then
      write(mpl%info,'(a13,a,i3,a,f10.2,a,f10.2,a)') '','Level ',nam%levs(il0),': '//trim(mpl%aqua), &
 & nicas_blk%rhs_avg(il0)*reqkm,trim(mpl%black)//' km  / '//trim(mpl%aqua),rvs_avg(il0),trim(mpl%black)//' vert. unit'
      if (nicas_blk%verbosity) call mpl%flush
   end if
end do
nicas_blk%var2d = (geom%nl0>1).and.(count(nicas_blk%vlev)==1)
if (nicas_blk%var2d) then
   write(mpl%info,'(a10,a)') '','This variable is 2D'
   if (nicas_blk%verbosity) call mpl%flush
end if
if (.not.any(nicas_blk%vlev)) call mpl%abort(subr,'no valid level')

if ((trim(nicas_blk%subsamp)=='h').or.(trim(nicas_blk%subsamp)=='hv').or.(trim(nicas_blk%subsamp)=='hvh')) then
   ! Basic horizontal mesh defined with the minimum support radius
   rhs_min = mpl%msv%valr
   do ic0a=1,geom%nc0a
      if (any(nicas_blk%vlev.and.geom%gmask_c0a(ic0a,:))) rhs_min(ic0a) = minval(cmat_blk%rhs(ic0a,:), &
 & mask=nicas_blk%vlev.and.geom%gmask_c0a(ic0a,:))
   end do
   call mpl%f_comm%allreduce(sum(rhs_min,mask=mpl%msv%isnot(rhs_min)),rhs_minavg,fckit_mpi_sum())
   rhs_min_norm = real(count(mpl%msv%isnot(rhs_min)),kind_real)
   call mpl%f_comm%allreduce(rhs_min_norm,rhs_min_norm_tot,fckit_mpi_sum())
   rhs_minavg = rhs_minavg/rhs_min_norm_tot
   nicas_blk%nc1 = floor(2.0*maxval(geom%area)*nam%resol**2/(sqrt(3.0)*rhs_minavg**2))
   write(mpl%info,'(a10,a,i8)') '','Estimated nc1 from horizontal support radius: ',nicas_blk%nc1
   if (nicas_blk%verbosity) call mpl%flush
   if (nicas_blk%nc1>geom%nc0_gmask(0)) then
      if (nicas_blk%verbosity) call mpl%warning(subr,'required nc1 larger than mask size, resetting to mask size')
      nicas_blk%nc1 = geom%nc0_gmask(0)
   end if
   if (nicas_blk%nc1>nam%nc1max) then
      if (nicas_blk%verbosity) call mpl%warning(subr,'required nc1 larger than nc1max, resetting to nc1max')
      nicas_blk%nc1 = nam%nc1max
   end if
   if (nicas_blk%nc1<3) call mpl%abort(subr,'nicas_blk%nc1 lower than 3')
   write(mpl%info,'(a10,a,i8)') '','Final nc1: ',nicas_blk%nc1
   if (nicas_blk%verbosity) call mpl%flush
   write(mpl%info,'(a10,a,f5.2)') '','Effective horizontal resolution: ',sqrt(real(nicas_blk%nc1,kind_real)*sqrt(3.0) &
 & *rhs_minavg**2/(2.0*maxval(geom%area)))
   if (nicas_blk%verbosity) call mpl%flush
else
   ! Use the Sc0 subset
   nicas_blk%nc1 = geom%nc0_gmask(0)
end if

! Allocation
allocate(nicas_blk%c1_to_c0(nicas_blk%nc1))

! Compute subset
write(mpl%info,'(a10,a)') '','Compute horizontal subset C1: '
if (nicas_blk%verbosity) call mpl%flush(.false.)

! Mask initialization
mask_hor_c0a = geom%gmask_hor_c0a

if (nam%check_no_point_nicas) then
   ! Mask points on the last MPI task
   if (mpl%myproc==mpl%nproc) mask_hor_c0a = .false.
end if

! Compute subsampling
call initialize_sampling(mpl,rng,maxval(geom%area),geom%nc0a,geom%lon_c0a,geom%lat_c0a,mask_hor_c0a,rhs_min,geom%c0a_to_c0, &
 & nam%ntry,nam%nrep,nicas_blk%nc1,nicas_blk%c1_to_c0,fast=nam%fast_sampling,verbosity=nicas_blk%verbosity, &
 & n_uni=geom%nc0u,uni_to_loc=geom%c0u_to_c0a,tree_uni=geom%tree_c0u)

! Count Sc1 point in universe
nicas_blk%nc1u = 0
do ic1=1,nicas_blk%nc1
   ic0 = nicas_blk%c1_to_c0(ic1)
   iproc = geom%c0_to_proc(ic0)
   if (geom%myuniverse(iproc)) nicas_blk%nc1u = nicas_blk%nc1u+1
end do

! Allocation
allocate(nicas_blk%c1u_to_c1(nicas_blk%nc1u))
allocate(nicas_blk%c1_to_c1u(nicas_blk%nc1))
allocate(nicas_blk%lon_c1u(nicas_blk%nc1u))
allocate(nicas_blk%lat_c1u(nicas_blk%nc1u))
allocate(nicas_blk%vunit_c1u(nicas_blk%nc1u,geom%nl0))
allocate(nicas_blk%gmask_c1u(nicas_blk%nc1u,geom%nl0))
allocate(nicas_blk%gmask_hor_c1u(nicas_blk%nc1u))

! Conversions
nicas_blk%c1_to_c1u = mpl%msv%vali
ic1u = 0
do ic1=1,nicas_blk%nc1
   ic0 = nicas_blk%c1_to_c0(ic1)
   iproc = geom%c0_to_proc(ic0)
   if (geom%myuniverse(iproc)) then
      ic1u = ic1u+1
      ic0u = geom%c0_to_c0u(ic0)
      nicas_blk%c1u_to_c1(ic1u) = ic1
      nicas_blk%c1_to_c1u(ic1) = ic1u
      nicas_blk%lon_c1u(ic1u) = geom%lon_c0u(ic0u)
      nicas_blk%lat_c1u(ic1u) = geom%lat_c0u(ic0u)
      nicas_blk%vunit_c1u(ic1u,:) = geom%vunit_c0u(ic0u,:)
      nicas_blk%gmask_c1u(ic1u,:) = geom%gmask_c0u(ic0u,:)
      nicas_blk%gmask_hor_c1u(ic1u) = geom%gmask_hor_c0u(ic0u)
   end if
end do

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
integer :: il0,il0_prev,il1,ic0a
real(kind_real) :: distnormmin,distnorm(geom%nc0a),rv

! Allocation
allocate(nicas_blk%slev(geom%nl0))

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

if ((trim(nicas_blk%subsamp)=='hv').or.(trim(nicas_blk%subsamp)=='vh').or.(trim(nicas_blk%subsamp)=='hvh')) then
   ! Vertical sampling
   write(mpl%info,'(a10,a)') '','Compute vertical subset L1'
   if (nicas_blk%verbosity) call mpl%flush

   do il0=1,geom%nl0
      if (nicas_blk%vlev(il0)) then
         ! Look for convolution levels
         if ((il0==nicas_blk%il0_first).or.(il0==nicas_blk%il0_last)) then
            ! Keep first and last levels
            nicas_blk%slev(il0) = .true.
         else
            ! Compute minimum normalized distance with level il0_prev
            distnorm = huge_real
            do ic0a=1,geom%nc0a
               if (geom%gmask_c0a(ic0a,il0)) then
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
   do il0=nicas_blk%il0_first,nicas_blk%il0_last
      nicas_blk%slev(il0) = .true.
   end do
end if

! Count effective levels
nicas_blk%nl1 = count(nicas_blk%slev)
allocate(nicas_blk%l1_to_l0(nicas_blk%nl1))
write(mpl%info,'(a10,a)') '','Effective levels: '
if (nicas_blk%verbosity) call mpl%flush(.false.)
il1 = 0
do il0=1,geom%nl0
   if (nicas_blk%slev(il0)) then
      write(mpl%info,'(i3,a)') nam%levs(il0),' '
      if (nicas_blk%verbosity) call mpl%flush(.false.)
      il1 = il1+1
      nicas_blk%l1_to_l0(il1) = il0
   end if
end do
write(mpl%info,'(a)') ''
if (nicas_blk%verbosity) call mpl%flush

! Inverse conversion
allocate(nicas_blk%l0_to_l1(geom%nl0))
nicas_blk%l0_to_l1 = mpl%msv%vali
do il1=1,nicas_blk%nl1
   il0 = nicas_blk%l1_to_l0(il1)
   nicas_blk%l0_to_l1(il0) = il1
end do

end subroutine nicas_blk_compute_sampling_v

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
integer :: ic0,ic0a,ic1,ic1a,ic1u,il0,il1,iproc

! Define halo A
nicas_blk%nc1a = 0
do ic1=1,nicas_blk%nc1
   ic0 = nicas_blk%c1_to_c0(ic1)
   if (geom%c0_to_proc(ic0)==mpl%myproc) nicas_blk%nc1a = nicas_blk%nc1a+1
end do

! Allocation
allocate(nicas_blk%c1a_to_c1(nicas_blk%nc1a))
allocate(nicas_blk%c1u_to_c1a(nicas_blk%nc1u))
allocate(nicas_blk%c1a_to_c0a(nicas_blk%nc1a))
allocate(nicas_blk%lon_c1a(nicas_blk%nc1a))
allocate(nicas_blk%lat_c1a(nicas_blk%nc1a))
allocate(nicas_blk%gmask_c1a(nicas_blk%nc1a,nicas_blk%nl1))

! Global-local conversions for halo A
ic1a = 0
do ic1=1,nicas_blk%nc1
   ic0 = nicas_blk%c1_to_c0(ic1)
   iproc = geom%c0_to_proc(ic0)
   if (geom%myuniverse(iproc)) then
      ic1u = nicas_blk%c1_to_c1u(ic1)
      if (iproc==mpl%myproc) then
         ic1a = ic1a+1
         nicas_blk%c1a_to_c1(ic1a) = ic1
         nicas_blk%c1u_to_c1a(ic1u) = ic1a
      end if
   end if
end do

! Conversion between subsets Sc0 et Sc1, halo A
do ic1a=1,nicas_blk%nc1a
   ic1 = nicas_blk%c1a_to_c1(ic1a)
   ic0 = nicas_blk%c1_to_c0(ic1)
   ic0a = geom%c0_to_c0a(ic0)
   nicas_blk%c1a_to_c0a(ic1a) = ic0a
end do

! Fields on subset Sc1
do ic1a=1,nicas_blk%nc1a
   ic0a = nicas_blk%c1a_to_c0a(ic1a)
   nicas_blk%lon_c1a(ic1a) = geom%lon_c0a(ic0a)
   nicas_blk%lat_c1a(ic1a) = geom%lat_c0a(ic0a)
   do il1=1,nicas_blk%nl1
      il0 = nicas_blk%l1_to_l0(il1)
      nicas_blk%gmask_c1a(ic1a,il1) = geom%gmask_c0a(ic0a,il0)
   end do
end do

! Setup subset Sc1 communication, local to universe
call nicas_blk%com_AU%setup(mpl,'com_AU',nicas_blk%nc1a,nicas_blk%nc1u,nicas_blk%nc1,nicas_blk%c1a_to_c1,nicas_blk%c1u_to_c1)

end subroutine nicas_blk_compute_mpi_a

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
integer :: ic1a,ic0a,il0,il1,ic1,ic1u,ic2,is,isu,ic0,iproc
integer :: nc1a_gmask(nicas_blk%nl1),nc1_gmask(nicas_blk%nl1),c1ul1_to_s(nicas_blk%nc1u,nicas_blk%nl1)
integer,allocatable :: c2_to_c1(:)
real(kind_real) :: rhs_c1a(nicas_blk%nc1a)
character(len=1024),parameter :: subr = 'nicas_blk_compute_sampling_c2'

! Allocation
allocate(nicas_blk%nc2(nicas_blk%nl1))
allocate(nicas_blk%nc2u(nicas_blk%nl1))
allocate(nicas_blk%gmask_c2u(nicas_blk%nc1u,nicas_blk%nl1))

! Initialization
do il1=1,nicas_blk%nl1
   il0 = nicas_blk%l1_to_l0(il1)
   nc1a_gmask(il1) = count(nicas_blk%gmask_c1a(:,il1))
end do
call mpl%f_comm%allreduce(nc1a_gmask,nc1_gmask,fckit_mpi_sum())
is = 0
c1ul1_to_s = mpl%msv%vali

! Vertically dependent horizontal subsampling
if ((trim(nicas_blk%subsamp)=='h').or.(trim(nicas_blk%subsamp)=='vh').or.(trim(nicas_blk%subsamp)=='hvh')) then
   write(mpl%info,'(a10,a)') '','Compute vertically dependent horizontal subsampling: '
   if (nicas_blk%verbosity) call mpl%flush
end if
do il1=1,nicas_blk%nl1
   il0 = nicas_blk%l1_to_l0(il1)
   if (nicas_blk%vlev(il0).and.((trim(nicas_blk%subsamp)=='h').or.(trim(nicas_blk%subsamp)=='vh').or. &
 & (trim(nicas_blk%subsamp)=='hvh'))) then
      write(mpl%info,'(a13,a,i3,a)') '','Level ',il1,':'
      if (nicas_blk%verbosity) call mpl%flush

      ! Compute nc2
      nicas_blk%nc2(il1) = floor(2.0*geom%area(il0)*nam%resol**2/(sqrt(3.0)*nicas_blk%rhs_avg(il0)**2))
      write(mpl%info,'(a16,a,i8)') '','Estimated nc2 from horizontal support radius: ',nicas_blk%nc2(il1)
      if (nicas_blk%verbosity) call mpl%flush
      if (nicas_blk%nc2(il1)<3) call mpl%abort(subr,'nicas_blk%nc2 lower than 3')
      nicas_blk%nc2(il1) = min(nicas_blk%nc2(il1),nc1_gmask(il1))
      write(mpl%info,'(a16,a,i8)') '','Final nc2: ',nicas_blk%nc2(il1)
      if (nicas_blk%verbosity) call mpl%flush
   else
      ! No C2 subsampling
      nicas_blk%nc2(il1) = nc1_gmask(il1)
   end if

   ! Compute horizontal subset C2
   write(mpl%info,'(a16,a)') '','Compute horizontal subset C2: '
   if (nicas_blk%verbosity) call mpl%flush(.false.)

   ! Allocation
   allocate(c2_to_c1(nicas_blk%nc2(il1)))

   ! Initialization
   do ic1a=1,nicas_blk%nc1a
      ic0a = nicas_blk%c1a_to_c0a(ic1a)
      rhs_c1a(ic1a) = cmat_blk%rhs(ic0a,il0)
   end do

   ! Initialize sampling
   call initialize_sampling(mpl,rng,geom%area(il0),nicas_blk%nc1a,nicas_blk%lon_c1a,nicas_blk%lat_c1a, &
 & nicas_blk%gmask_c1a(:,il1),rhs_c1a,nicas_blk%c1a_to_c1,nam%ntry,nam%nrep,nicas_blk%nc2(il1),c2_to_c1, &
 & fast=nam%fast_sampling,verbosity=nicas_blk%verbosity)

   ! Fill subset Sc2 mask
   nicas_blk%gmask_c2u(:,il1) = .false.
   do ic2=1,nicas_blk%nc2(il1)
      is = is+1
      ic1 = c2_to_c1(ic2)
      ic0 = nicas_blk%c1_to_c0(ic1)
      iproc = geom%c0_to_proc(ic0)
      if (geom%myuniverse(iproc)) then
         ic1u = nicas_blk%c1_to_c1u(ic1)
         nicas_blk%gmask_c2u(ic1u,il1) = .true.
         c1ul1_to_s(ic1u,il1) = is
      end if
   end do

   ! Number of subset Sc2 points in universe
   nicas_blk%nc2u(il1) = count(nicas_blk%gmask_c2u(:,il1))

   ! Release memory
   deallocate(c2_to_c1)
end do

! Size
nicas_blk%ns = sum(nicas_blk%nc2)
nicas_blk%nsu = sum(nicas_blk%nc2u)

! Allocation
allocate(nicas_blk%su_to_c1u(nicas_blk%nsu))
allocate(nicas_blk%su_to_l1(nicas_blk%nsu))
allocate(nicas_blk%c1ul1_to_su(nicas_blk%nc1u,nicas_blk%nl1))
allocate(nicas_blk%su_to_s(nicas_blk%nsu))

! Conversions
isu = 0
do il1=1,nicas_blk%nl1
   do ic1u=1,nicas_blk%nc1u
      if (nicas_blk%gmask_c2u(ic1u,il1)) then
         isu = isu+1
         is = c1ul1_to_s(ic1u,il1)
         nicas_blk%su_to_c1u(isu) = ic1u
         nicas_blk%su_to_l1(isu) = il1
         nicas_blk%c1ul1_to_su(ic1u,il1) = isu
         nicas_blk%su_to_s(isu) = is
      end if
   end do
end do

end subroutine nicas_blk_compute_sampling_c2

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
integer :: il0i,ic0,ic1,ic1a,ic1u,ic1b_h,ic1b,il0,il1,isa,isb,i_s,isu,is,jc1u,jsu,nc1b_h,ifmt,iproc
integer,allocatable :: c1b_h_to_c1u(:),c1b_h_to_c1b(:),sa_to_c1a(:),sa_to_l1(:),sb_to_s(:)
real(kind_real),allocatable :: lon_c1b_h(:),lat_c1b_h(:)
logical :: lcheck_c1b_h(nicas_blk%nc1u),lcheck_c1b(nicas_blk%nc1u),inside
logical,allocatable :: gmask_c1b_h(:,:)
character(len=1024),parameter :: subr = 'nicas_blk_compute_mpi_ab'

! Allocation
allocate(nicas_blk%h(geom%nl0i))
allocate(nicas_blk%s(nicas_blk%nl1))
allocate(nicas_blk%lcheck_sa(nicas_blk%nsu))
allocate(nicas_blk%lcheck_sb(nicas_blk%nsu))

! Initialization
ifmt = 0
if (nicas_blk%verbosity) ifmt = 10

! Compute interpolation
do il0i=1,geom%nl0i
   write(nicas_blk%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
   call nicas_blk%h(il0i)%interp(mpl,rng,nam,geom,il0i,nicas_blk%nc1u,nicas_blk%lon_c1u,nicas_blk%lat_c1u, &
 & nicas_blk%gmask_c1u(:,il0i),geom%nc0a,geom%lon_c0a,geom%lat_c0a,geom%gmask_c0a(:,il0i),ifmt)
end do

! Define halo A
nicas_blk%lcheck_sa = .false.
do isu=1,nicas_blk%nsu
   ic1u = nicas_blk%su_to_c1u(isu)
   ic1 = nicas_blk%c1u_to_c1(ic1u)
   ic0 = nicas_blk%c1_to_c0(ic1)
   iproc = geom%c0_to_proc(ic0)
   if (iproc==mpl%myproc) nicas_blk%lcheck_sa(isu) = .true.
end do
nicas_blk%nsa = count(nicas_blk%lcheck_sa)

! Define halo B (after first horizontal interpolation)
lcheck_c1b_h = .false.
do il0i=1,geom%nl0i
   do i_s=1,nicas_blk%h(il0i)%n_s
      jc1u = nicas_blk%h(il0i)%col(i_s)
      lcheck_c1b_h(jc1u) = .true.
   end do
end do
nc1b_h = count(lcheck_c1b_h)

! Allocation
allocate(c1b_h_to_c1u(nc1b_h))
allocate(lon_c1b_h(nc1b_h))
allocate(lat_c1b_h(nc1b_h))
allocate(gmask_c1b_h(nc1b_h,nicas_blk%nl1))

! Conversion
ic1b_h = 0
do ic1u=1,nicas_blk%nc1u
   if (lcheck_c1b_h(ic1u)) then
      ic1b_h = ic1b_h+1
      c1b_h_to_c1u(ic1b_h) = ic1u
   end if
end do

! Fields
do ic1b_h=1,nc1b_h
   ic1u = c1b_h_to_c1u(ic1b_h)
   lon_c1b_h(ic1b_h) = nicas_blk%lon_c1u(ic1u)
   lat_c1b_h(ic1b_h) = nicas_blk%lat_c1u(ic1u)
   do il1=1,nicas_blk%nl1
      il0 = nicas_blk%l1_to_l0(il1)
      gmask_c1b_h(ic1b_h,il1) = nicas_blk%gmask_c1u(ic1u,il0)
   end do
end do

! Compute interpolation
do il1=1,nicas_blk%nl1
   write(nicas_blk%s(il1)%prefix,'(a,i3.3)') 's_',il1

   ! Compute interpolation
   il0 = nicas_blk%l1_to_l0(il1)
   call nicas_blk%s(il1)%interp(mpl,rng,nam,geom,il0,nicas_blk%nc1u,nicas_blk%lon_c1u,nicas_blk%lat_c1u, &
 & nicas_blk%gmask_c2u(:,il1),nc1b_h,lon_c1b_h,lat_c1b_h,gmask_c1b_h(:,il1),ifmt)
end do

! Define halo B (required for the second horizontal interpolation)
nicas_blk%lcheck_sb = nicas_blk%lcheck_sa
lcheck_c1b = lcheck_c1b_h
do il1=1,nicas_blk%nl1
   do i_s=1,nicas_blk%s(il1)%n_s
      jc1u = nicas_blk%s(il1)%col(i_s)
      jsu = nicas_blk%c1ul1_to_su(jc1u,il1)
      lcheck_c1b(jc1u) = .true.
      nicas_blk%lcheck_sb(jsu) = .true.
   end do
end do
nicas_blk%nc1b = count(lcheck_c1b)
nicas_blk%nsb = count(nicas_blk%lcheck_sb)

! Check halos consistency
do isu=1,nicas_blk%nsu
   if (nicas_blk%lcheck_sa(isu).and.(.not.nicas_blk%lcheck_sb(isu))) call mpl%abort(subr,'point in halo A but not in halo B')
end do

! Allocation
allocate(sa_to_c1a(nicas_blk%nsa))
allocate(sa_to_l1(nicas_blk%nsa))
allocate(nicas_blk%sa_to_su(nicas_blk%nsa))
allocate(nicas_blk%sa_to_s(nicas_blk%nsa))
allocate(nicas_blk%hash_sa(nicas_blk%nsa))
allocate(nicas_blk%c1b_to_c1u(nicas_blk%nc1b))
allocate(nicas_blk%c1u_to_c1b(nicas_blk%nc1u))
allocate(sb_to_s(nicas_blk%nsb))
allocate(nicas_blk%sb_to_su(nicas_blk%nsb))
allocate(nicas_blk%su_to_sb(nicas_blk%nsu))
allocate(c1b_h_to_c1b(nc1b_h))
allocate(nicas_blk%sb_to_c1b(nicas_blk%nsb))
allocate(nicas_blk%sb_to_l1(nicas_blk%nsb))
allocate(nicas_blk%c1bl1_to_sb(nicas_blk%nc1b,nicas_blk%nl1))
allocate(nicas_blk%vbot(nicas_blk%nc1b))
allocate(nicas_blk%vtop(nicas_blk%nc1b))

! Conversions
isa = 0
do isu=1,nicas_blk%nsu
   if (nicas_blk%lcheck_sa(isu)) then
      isa = isa+1
      is = nicas_blk%su_to_s(isu)
      ic1u = nicas_blk%su_to_c1u(isu)
      ic1a = nicas_blk%c1u_to_c1a(ic1u)
      il1 = nicas_blk%su_to_l1(isu)
      nicas_blk%sa_to_su(isa) = isu
      nicas_blk%sa_to_s(isa) = is
      sa_to_l1(isa) = il1
      sa_to_c1a(isa) = ic1a
   end if
end do
nicas_blk%c1u_to_c1b = mpl%msv%vali
ic1b = 0
do ic1u=1,nicas_blk%nc1u
   if (lcheck_c1b(ic1u)) then
      ic1b = ic1b+1
      nicas_blk%c1b_to_c1u(ic1b) = ic1u
      nicas_blk%c1u_to_c1b(ic1u) = ic1b
   end if
end do
nicas_blk%su_to_sb = mpl%msv%vali
isb = 0
do isu=1,nicas_blk%nsu
   if (nicas_blk%lcheck_sb(isu)) then
      isb = isb+1
      is = nicas_blk%su_to_s(isu)
      sb_to_s(isb) = is
      nicas_blk%sb_to_su(isb) = isu
      nicas_blk%su_to_sb(isu) = isb
   end if
end do
do ic1b_h=1,nc1b_h
   ic1u = c1b_h_to_c1u(ic1b_h)
   ic1b = nicas_blk%c1u_to_c1b(ic1u)
   c1b_h_to_c1b(ic1b_h) = ic1b
end do

! Local interpolation source
do il0i=1,geom%nl0i
   nicas_blk%h(il0i)%n_src = nicas_blk%nc1b
   do i_s=1,nicas_blk%h(il0i)%n_s
      nicas_blk%h(il0i)%col(i_s) = nicas_blk%c1u_to_c1b(nicas_blk%h(il0i)%col(i_s))
   end do
end do
do il1=1,nicas_blk%nl1
   nicas_blk%s(il1)%n_src = nicas_blk%nc1b
   nicas_blk%s(il1)%n_dst = nicas_blk%nc1b
   do i_s=1,nicas_blk%s(il1)%n_s
      nicas_blk%s(il1)%row(i_s) = c1b_h_to_c1b(nicas_blk%s(il1)%row(i_s))
      nicas_blk%s(il1)%col(i_s) = nicas_blk%c1u_to_c1b(nicas_blk%s(il1)%col(i_s))
   end do
end do

! Setup communications
call nicas_blk%com_AB%setup(mpl,'com_AB',nicas_blk%nsa,nicas_blk%nsb,nicas_blk%ns,nicas_blk%sa_to_s,sb_to_s)

! Conversions
nicas_blk%c1bl1_to_sb = mpl%msv%vali
do isb=1,nicas_blk%nsb
   isu = nicas_blk%sb_to_su(isb)
   ic1u = nicas_blk%su_to_c1u(isu)
   il1 = nicas_blk%su_to_l1(isu)
   ic1b = nicas_blk%c1u_to_c1b(ic1u)
   nicas_blk%sb_to_c1b(isb) = ic1b
   nicas_blk%sb_to_l1(isb) = il1
   nicas_blk%c1bl1_to_sb(ic1b,il1) = isb
end do

! Hash value
do isa=1,nicas_blk%nsa
   ic1a = sa_to_c1a(isa)
   il1 = sa_to_l1(isa)
   nicas_blk%hash_sa(isa) = lonlathash(nicas_blk%lon_c1a(ic1a),nicas_blk%lat_c1a(ic1a),il1)
end do

! Find bottom and top for each point of subset Sc1, halo B
if (nicas_blk%var2d) then
   do il0=1,geom%nl0
      if (nicas_blk%vlev(il0)) then
         nicas_blk%vbot = il0
         nicas_blk%vtop = il0
      end if
   end do
else
   ! Initialization
   nicas_blk%vbot = mpl%msv%vali
   nicas_blk%vtop = mpl%msv%vali

   ! Loop over points
   !$omp parallel do schedule(static) private(ic1b,ic1u,inside,il1,il0)
   do ic1b=1,nicas_blk%nc1b
      ic1u = nicas_blk%c1b_to_c1u(ic1b)
      inside = .false.
      nicas_blk%vtop(ic1b) = geom%nl0
      do il1=1,nicas_blk%nl1
         il0 = nicas_blk%l1_to_l0(il1)
         if (.not.inside.and.nicas_blk%gmask_c1u(ic1u,il0)) then
            ! Bottom level
            nicas_blk%vbot(ic1b) = il0
            inside = .true.
         end if
         if (inside.and.(.not.nicas_blk%gmask_c1u(ic1u,il0))) then
            ! Top level
            nicas_blk%vtop(ic1b) = il0
            inside = .false.
         end if
      end do
      if (mpl%msv%is(nicas_blk%vbot(ic1b))) call mpl%abort(subr,'bottom level not found')
      if (mpl%msv%is(nicas_blk%vtop(ic1b))) call mpl%abort(subr,'top level not found')
      if (nicas_blk%vbot(ic1b)>nicas_blk%vtop(ic1b)) call mpl%abort(subr,'non contiguous mask')
   end do
   !$omp end parallel do
end if

! Release memory
deallocate(c1b_h_to_c1u)
deallocate(lon_c1b_h)
deallocate(lat_c1b_h)
deallocate(gmask_c1b_h)
deallocate(sa_to_c1a)
deallocate(sa_to_l1)
deallocate(sb_to_s)
deallocate(c1b_h_to_c1b)

end subroutine nicas_blk_compute_mpi_ab

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_interp_v
! Purpose: compute vertical interpolation
!----------------------------------------------------------------------
subroutine nicas_blk_compute_interp_v(nicas_blk,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(geom_type),intent(in) :: geom               ! Geometry

! Local variables
integer :: ic1b,ic1u,jl0,il0,jl1,il0inf,il0sup

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
            ic1u = nicas_blk%c1b_to_c1u(ic1b)
            nicas_blk%v%Svec(nicas_blk%v%n_s,ic1b) = abs(nicas_blk%vunit_c1u(ic1u,il0sup)-nicas_blk%vunit_c1u(ic1u,il0)) &
 & /abs(nicas_blk%vunit_c1u(ic1u,il0sup)-nicas_blk%vunit_c1u(ic1u,il0inf))
         end do
         nicas_blk%v%n_s = nicas_blk%v%n_s+1
         nicas_blk%v%row(nicas_blk%v%n_s) = il0
         nicas_blk%v%col(nicas_blk%v%n_s) = il0sup
         do ic1b=1,nicas_blk%nc1b
            ic1u = nicas_blk%c1b_to_c1u(ic1b)
            nicas_blk%v%Svec(nicas_blk%v%n_s,ic1b) = abs(nicas_blk%vunit_c1u(ic1u,il0)-nicas_blk%vunit_c1u(ic1u,il0inf)) &
 & /abs(nicas_blk%vunit_c1u(ic1u,il0sup)-nicas_blk%vunit_c1u(ic1u,il0inf))
         end do
      end do
      il0inf = jl0
   end if
end do

! Conversion
nicas_blk%v%col = nicas_blk%l0_to_l1(nicas_blk%v%col)

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
integer :: n_s_max,ithread,isu,ic1u,jc1u,il1,il0,j,jsu,isb,ic1b,ic0a,ic1a,i_s,jc,kc,ksu,jbd,jl1,ic1bb,isbb
integer :: c_n_s(mpl%nthread)
integer,allocatable :: nn(:),nn_index(:),inec(:),c_ind(:,:)
real(kind_real),allocatable :: rh_c1a(:,:),rv_c1a(:,:)
real(kind_real),allocatable :: H11_c1a(:,:),H22_c1a(:,:),H33_c1a(:,:),H12_c1a(:,:),Hcoef_c1a(:,:),Hcoef_c1u(:,:)
real(kind_real),allocatable :: Hcoef(:,:)
real(kind_real),allocatable :: c_S(:,:),c_S_conv(:)
logical :: add_op,lcheck_c1bb(nicas_blk%nc1u)
type(linop_type) :: ctmp,c(mpl%nthread)

! Associate
associate(ib=>nicas_blk%ib)

! Set anisotropic parameter
nicas_blk%anisotropic = cmat_blk%anisotropic

! Allocation
call nicas_blk%tree%alloc(mpl,nicas_blk%nc1u)

! Initialization
call nicas_blk%tree%init(nicas_blk%lon_c1u,nicas_blk%lat_c1u)

! Find largest possible radius
call mpl%f_comm%allreduce(maxval(cmat_blk%rh,mask=mpl%msv%isnot(cmat_blk%rh)),nicas_blk%rhmax,fckit_mpi_max())
if ((.not.nicas_blk%smoother).and.nicas_blk%sqrt_rescaling) then
   ! Square-root rescaling of the largest possible radius
   if (nicas_blk%anisotropic) then
      nicas_blk%rhmax = nicas_blk%rhmax*sqrt_h
   else
      nicas_blk%rhmax = nicas_blk%rhmax*sqrt_r
   end if
end if

if ((nicas_blk%lsqrt==1).or.nicas_blk%smoother) then
   ! Copy
   nicas_blk%nc1bb = nicas_blk%nc1b
else
   ! Allocation
   allocate(nn(nicas_blk%nc1b))

   do ic1b=1,nicas_blk%nc1b
      ! Indices
      ic1u = nicas_blk%c1b_to_c1u(ic1b)

      ! Count nearest neighbors
      call nicas_blk%tree%count_nearest_neighbors(nicas_blk%lon_c1u(ic1u),nicas_blk%lat_c1u(ic1u),nicas_blk%rhmax,nn(ic1b))
   end do

   ! Initialization
   write(mpl%info,'(a10,a)') '','Define extended halo: '
   if (nicas_blk%verbosity) call mpl%flush(.false.)
   if (nicas_blk%verbosity) call mpl%prog_init(nicas_blk%nc1b)
   lcheck_c1bb = .false.

   do ic1b=1,nicas_blk%nc1b
      ! Indices
      ic1u = nicas_blk%c1b_to_c1u(ic1b)

      ! Allocation
      allocate(nn_index(nn(ic1b)))

      ! Find nearest neighbors
      call nicas_blk%tree%find_nearest_neighbors(nicas_blk%lon_c1u(ic1u),nicas_blk%lat_c1u(ic1u),nn(ic1b),nn_index)

      ! Fill mask
      do j=1,nn(ic1b)
         jc1u = nn_index(j)
         lcheck_c1bb(jc1u) = .true.
      end do

      ! Release memory
      deallocate(nn_index)

      ! Update
      if (nicas_blk%verbosity) call mpl%prog_print(ic1b)
   end do
   if (nicas_blk%verbosity) call mpl%prog_final

   ! Halo size
   nicas_blk%nc1bb = count(lcheck_c1bb)
   write(mpl%info,'(a10,a,i6,a,i6)') '','Halo sizes nc1b / nc1bb: ',nicas_blk%nc1b,' / ',nicas_blk%nc1bb
   if (nicas_blk%verbosity) call mpl%flush

   ! Release memory
   deallocate(nn)
end if

! Allocation
allocate(nicas_blk%c1bb_to_c1u(nicas_blk%nc1bb))
allocate(nicas_blk%c1u_to_c1bb(nicas_blk%nc1u))

if ((nicas_blk%lsqrt==1).or.nicas_blk%smoother) then
   ! Copy
   nicas_blk%c1bb_to_c1u = nicas_blk%c1b_to_c1u
   nicas_blk%c1u_to_c1bb = nicas_blk%c1u_to_c1b
   nicas_blk%nsbb = nicas_blk%nsb
else
   ! Conversions
   nicas_blk%c1u_to_c1bb = mpl%msv%vali
   ic1bb = 0
   do ic1u=1,nicas_blk%nc1u
      if (lcheck_c1bb(ic1u)) then
         ic1bb = ic1bb+1
         nicas_blk%c1bb_to_c1u(ic1bb) = ic1u
         nicas_blk%c1u_to_c1bb(ic1u) = ic1bb
      end if
   end do

   ! Count points in extended halo
   nicas_blk%nsbb = 0
   do isu=1,nicas_blk%nsu
      ic1u = nicas_blk%su_to_c1u(isu)
      if (lcheck_c1bb(ic1u)) nicas_blk%nsbb = nicas_blk%nsbb+1
   end do
   write(mpl%info,'(a10,a,i6,a,i6)') '','Halo sizes nsb / nsbb:   ',nicas_blk%nsb,' / ',nicas_blk%nsbb
   if (nicas_blk%verbosity) call mpl%flush
end if

! Allocation
allocate(nicas_blk%sbb_to_su(nicas_blk%nsbb))

if ((nicas_blk%lsqrt==1).or.nicas_blk%smoother) then
   ! Copy
   nicas_blk%sbb_to_su = nicas_blk%sb_to_su
else
   ! Conversions
   isbb = 0
   do isu=1,nicas_blk%nsu
      ic1u = nicas_blk%su_to_c1u(isu)
      if (lcheck_c1bb(ic1u)) then
         isbb = isbb+1
         nicas_blk%sbb_to_su(isbb) = isu
      end if
   end do
end if

! Compute horizontal and vertical parameters
write(mpl%info,'(a10,a)') '','Compute horizontal and vertical parameters'
if (nicas_blk%verbosity) call mpl%flush

! Allocation
allocate(rh_c1a(nicas_blk%nc1a,nicas_blk%nl1))
if (.not.nicas_blk%horizontal) allocate(rv_c1a(nicas_blk%nc1a,nicas_blk%nl1))
allocate(nicas_blk%rh_c1u(nicas_blk%nc1u,nicas_blk%nl1))
if (.not.nicas_blk%horizontal) allocate(nicas_blk%rv_c1u(nicas_blk%nc1u,nicas_blk%nl1))
if (nicas_blk%anisotropic) then
   allocate(H11_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(H22_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   if (.not.nicas_blk%horizontal) allocate(H33_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(H12_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(Hcoef_c1a(nicas_blk%nc1a,nicas_blk%nl1))
   allocate(nicas_blk%H11_c1u(nicas_blk%nc1u,nicas_blk%nl1))
   allocate(nicas_blk%H22_c1u(nicas_blk%nc1u,nicas_blk%nl1))
   if (.not.nicas_blk%horizontal) allocate(nicas_blk%H33_c1u(nicas_blk%nc1u,nicas_blk%nl1))
   allocate(nicas_blk%H12_c1u(nicas_blk%nc1u,nicas_blk%nl1))
   allocate(Hcoef_c1u(nicas_blk%nc1u,nicas_blk%nl1))
   allocate(nicas_blk%Hcoef(nicas_blk%nsbb))
end if

! Copy and rescale
write(mpl%info,'(a13,a)') '','Copy and rescale'
if (nicas_blk%verbosity) call mpl%flush
do il1=1,nicas_blk%nl1
   do ic1a=1,nicas_blk%nc1a
      ! Indices
      ic0a = nicas_blk%c1a_to_c0a(ic1a)
      il0 = nicas_blk%l1_to_l0(il1)

      if (geom%gmask_c0a(ic0a,il0)) then
         ! Copy
         rh_c1a(ic1a,il1) = cmat_blk%rh(ic0a,il0)
         if (.not.nicas_blk%horizontal) rv_c1a(ic1a,il1) = cmat_blk%rv(ic0a,il0)
         if (nicas_blk%anisotropic) then
            H11_c1a(ic1a,il1) = cmat_blk%H11(ic0a,il0)
            H22_c1a(ic1a,il1) = cmat_blk%H22(ic0a,il0)
            if (.not.nicas_blk%horizontal) H33_c1a(ic1a,il1) = cmat_blk%H33(ic0a,il0)
            H12_c1a(ic1a,il1) = cmat_blk%H12(ic0a,il0)
            Hcoef_c1a(ic1a,il1) = cmat_blk%Hcoef(ic0a,il0)
         end if

         if ((.not.nicas_blk%smoother).and.nicas_blk%sqrt_rescaling) then
            ! Square-root rescaling
            if (nicas_blk%anisotropic) then
               rh_c1a(ic1a,il1) = rh_c1a(ic1a,il1)*sqrt_h
               if (.not.nicas_blk%horizontal) rv_c1a(ic1a,il1) = rv_c1a(ic1a,il1)*sqrt_h
               H11_c1a(ic1a,il1) = H11_c1a(ic1a,il1)/sqrt_h**2
               H22_c1a(ic1a,il1) = H22_c1a(ic1a,il1)/sqrt_h**2
               if (.not.nicas_blk%horizontal) H33_c1a(ic1a,il1) = H33_c1a(ic1a,il1)/sqrt_h**2
               H12_c1a(ic1a,il1) = H12_c1a(ic1a,il1)/sqrt_h**2
            else
               rh_c1a(ic1a,il1) = rh_c1a(ic1a,il1)*sqrt_r
               if (.not.nicas_blk%horizontal) rv_c1a(ic1a,il1) = rv_c1a(ic1a,il1)*sqrt_r
            end if
         end if
      else
         ! Missing values
         rh_c1a(ic1a,il1) = mpl%msv%valr
         if (.not.nicas_blk%horizontal) rv_c1a(ic1a,il1) = mpl%msv%valr
         if (nicas_blk%anisotropic) then
            H11_c1a(ic1a,il1) = mpl%msv%valr
            H22_c1a(ic1a,il1) = mpl%msv%valr
            if (.not.nicas_blk%horizontal) H33_c1a(ic1a,il1) = mpl%msv%valr
            H12_c1a(ic1a,il1) = mpl%msv%valr
            Hcoef_c1a(ic1a,il1) = mpl%msv%valr
         end if
      end if
   end do
end do

! Communication
write(mpl%info,'(a13,a)') '','Communication'
if (nicas_blk%verbosity) call mpl%flush
call nicas_blk%com_AU%ext(mpl,nicas_blk%nl1,rh_c1a,nicas_blk%rh_c1u)
if (.not.nicas_blk%horizontal) call nicas_blk%com_AU%ext(mpl,nicas_blk%nl1,rv_c1a,nicas_blk%rv_c1u)
if (nicas_blk%anisotropic) then
   call nicas_blk%com_AU%ext(mpl,nicas_blk%nl1,H11_c1a,nicas_blk%H11_c1u)
   call nicas_blk%com_AU%ext(mpl,nicas_blk%nl1,H22_c1a,nicas_blk%H22_c1u)
   if (.not.nicas_blk%horizontal) call nicas_blk%com_AU%ext(mpl,nicas_blk%nl1,H33_c1a,nicas_blk%H33_c1u)
   call nicas_blk%com_AU%ext(mpl,nicas_blk%nl1,H12_c1a,nicas_blk%H12_c1u)
   call nicas_blk%com_AU%ext(mpl,nicas_blk%nl1,Hcoef_c1a,Hcoef_c1u)
end if

! Release memory
deallocate(rh_c1a)
if (.not.nicas_blk%horizontal) deallocate(rv_c1a)
if (nicas_blk%anisotropic) then
   deallocate(H11_c1a)
   deallocate(H22_c1a)
   if (.not.nicas_blk%horizontal) deallocate(H33_c1a)
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
deallocate(nicas_blk%rh_c1u)
if (.not.nicas_blk%horizontal) deallocate(nicas_blk%rv_c1u)
if (nicas_blk%anisotropic) then
   deallocate(nicas_blk%H11_c1u)
   deallocate(nicas_blk%H22_c1u)
   if (.not.nicas_blk%horizontal) deallocate(nicas_blk%H33_c1u)
   deallocate(nicas_blk%H12_c1u)
end if

if (nicas_blk%anisotropic) then
   ! Compute ball data
   write(mpl%info,'(a13,a)') '','Compute ball data'
   if (nicas_blk%verbosity) call mpl%flush
   !$omp parallel do schedule(static) private(isbb,isu,ic1u,il1,jbd,jc1u,jl1) firstprivate(Hcoef)
   do isbb=1,nicas_blk%nsbb
      ! Indices
      isu = nicas_blk%sbb_to_su(isbb)
      ic1u = nicas_blk%su_to_c1u(isu)
      il1 = nicas_blk%su_to_l1(isu)

      ! Allocation
      allocate(Hcoef(nicas_blk%nc1u,nicas_blk%nl1))

      ! Initialization
      Hcoef = mpl%msv%valr

      do jbd=1,nicas_blk%distnorm(isbb)%nbd
         ! Indices
         jc1u = nicas_blk%distnorm(isbb)%bd_to_c1u(jbd)
         jl1 = nicas_blk%distnorm(isbb)%bd_to_l1(jbd)
         Hcoef(jc1u,jl1) = sqrt(Hcoef_c1u(ic1u,il1)*Hcoef_c1u(jc1u,jl1))
      end do

      ! Pack data
      call nicas_blk%Hcoef(isbb)%pack(mpl,nicas_blk%nc1u,nicas_blk%nl1,Hcoef)

      ! Release memory
      deallocate(Hcoef)
   end do
   !$omp end parallel do

   ! Release memory
   deallocate(Hcoef_c1u)
end if

! Allocation
if (nicas_blk%smoother) allocate(nicas_blk%smoother_norm(nicas_blk%nsbb))

! Compute weights
call nicas_blk%compute_convol_weights(mpl,geom,ctmp)

! Release memory
do isbb=1,nicas_blk%nsbb
   call nicas_blk%distnorm(isbb)%dealloc
   if (nicas_blk%anisotropic) call nicas_blk%Hcoef(isbb)%dealloc
end do
deallocate(nicas_blk%distnorm)
if (nicas_blk%anisotropic) deallocate(nicas_blk%Hcoef)

if ((nicas_blk%lsqrt==1).or.nicas_blk%smoother) then
   ! Copy
   call nicas_blk%c%copy(ctmp)
else
   ! Compute convolution inverse mapping
   allocate(inec(nicas_blk%ns))
   inec = 0
   do i_s=1,ctmp%n_s
      isu = ctmp%col(i_s)
      inec(isu) = inec(isu)+1
   end do

   allocate(c_ind(maxval(inec),nicas_blk%nsu))
   allocate(c_S(maxval(inec),nicas_blk%nsu))
   c_ind = mpl%msv%vali
   c_S = mpl%msv%valr
   inec = 0
   do i_s=1,ctmp%n_s
      isu = ctmp%col(i_s)
      jsu = ctmp%row(i_s)
      inec(isu) = inec(isu)+1
      c_ind(inec(isu),isu) = jsu
      c_S(inec(isu),isu) = ctmp%S(i_s)
   end do

   ! Initialization
   write(mpl%info,'(a10,a)') '','Second pass:     '
   if (nicas_blk%verbosity) call mpl%flush(.false.)
   if (nicas_blk%verbosity) call mpl%prog_init(nicas_blk%nsb)
   n_s_max = 10*nint(real(geom%nc0u*geom%nl0,kind_real)/real(mpl%nthread*mpl%nproc,kind_real))
   c_n_s = 0
   do ithread=1,mpl%nthread
      c(ithread)%n_s = n_s_max
      call c(ithread)%alloc
      c(ithread)%row = mpl%msv%vali
      c(ithread)%col = mpl%msv%vali
      c(ithread)%S = mpl%msv%valr
   end do

   ! Apply convolution
   !$omp parallel do schedule(static) private(isb,isu,ithread,jc,jsu,kc,ksu,add_op) firstprivate(c_S_conv)
   do isb=1,nicas_blk%nsb
      ! Indices
      isu = nicas_blk%sb_to_su(isb)
      ithread = 1
!$    ithread = omp_get_thread_num()+1

      ! Allocation
      allocate(c_S_conv(nicas_blk%nsu))

      ! Initialization
      c_S_conv = 0.0

      ! Loop twice over points
      do jc=1,inec(isu)
         jsu = c_ind(jc,isu)
         do kc=1,inec(jsu)
            ksu = c_ind(kc,jsu)
            c_S_conv(ksu) = c_S_conv(ksu)+c_S(jc,isu)*c_S(kc,jsu)
         end do
      end do

      ! Store coefficient for convolution
      do jsu=1,nicas_blk%nsu
         add_op = .false.
         if (nicas_blk%mpicom==1) then
            add_op = (nicas_blk%lcheck_sb(jsu).and.(isu<=jsu)).or.(.not.nicas_blk%lcheck_sb(jsu))
         elseif (nicas_blk%mpicom==2) then
            add_op = nicas_blk%lcheck_sa(isu).and.((nicas_blk%lcheck_sa(jsu).and.(isu<=jsu)) &
 & .or.(.not.nicas_blk%lcheck_sa(jsu)))
         end if
         if (add_op) call c(ithread)%add_op(c_n_s(ithread),isu,jsu,c_S_conv(jsu))
      end do

      ! Release memory
      deallocate(c_S_conv)

      ! Update
      if (nicas_blk%verbosity) call mpl%prog_print(isb)
   end do
   !$omp end parallel do
   if (nicas_blk%verbosity) call mpl%prog_final

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
integer :: isu,ic1u,il0,jl1,np,np_new,j,k,ip,kc1u,jc1u,il1,dkl1,kl1,jp,isbb,djl1,jl0
integer,allocatable :: plist(:,:),plist_new(:,:)
real(kind_real) :: disttest
real(kind_real) :: dnb,dx,dy,dz,disthsq,distvsq,rhsq,rvsq,H11,H22,H33,H12
real(kind_real),allocatable :: distnorm(:,:),net_dnb(:,:,:,:)
logical :: add_to_front
logical,allocatable :: net_arc(:,:,:)
type(mesh_type) :: mesh

! Allocation
call mesh%alloc(nicas_blk%nc1u)

! Initialization
call mesh%init(mpl,rng,nicas_blk%lon_c1u,nicas_blk%lat_c1u)

! Allocation
allocate(net_arc(nicas_blk%nc1u,nicas_blk%nl1,mesh%maxcols))
allocate(net_dnb(nicas_blk%nc1u,nicas_blk%nl1,mesh%maxcols,-1:1))

! Find mesh neighbors
write(mpl%info,'(a10,a)') '','Find mesh neighbors: '
if (nicas_blk%verbosity) call mpl%flush(.false.)
if (nicas_blk%verbosity) call mpl%prog_init(nicas_blk%nc1u)
net_arc = .false.
do ic1u=1,nicas_blk%nc1u
   do j=1,mesh%rows(ic1u)%cols
      ! Index
      jc1u = mesh%rows(ic1u)%nodes(j)

      ! Check arc
      if (nam%mask_check) then
         do il1=1,nicas_blk%nl1
            il0 = nicas_blk%l1_to_l0(il1)
            call geom%check_arc(mpl,il0,nicas_blk%lon_c1u(ic1u),nicas_blk%lat_c1u(ic1u),nicas_blk%lon_c1u(jc1u), &
 & nicas_blk%lat_c1u(jc1u),net_arc(ic1u,il1,j))
         end do
      else
         net_arc(ic1u,:,j) = .true.
      end if
   end do

   ! Update
   if (nicas_blk%verbosity) call mpl%prog_print(ic1u)
end do
if (nicas_blk%verbosity) call mpl%prog_final

! Compute mesh edges distances
write(mpl%info,'(a10,a)') '','Compute mesh edges distances: '
if (nicas_blk%verbosity) call mpl%flush(.false.)
if (nicas_blk%verbosity) call mpl%prog_init(nicas_blk%nc1u)
net_dnb = 1.0
!$omp parallel do schedule(static) private(ic1u,j,jc1u,dnb,dx,dy,dz,il1,il0,djl1,jl1,jl0,H11,H22,H33), &
!$omp&                             private(H12,disthsq,distvsq,rhsq,rvsq)
do ic1u=1,nicas_blk%nc1u
   do j=1,mesh%rows(ic1u)%cols
      ! Indices
      jc1u = mesh%rows(ic1u)%nodes(j)

      if (nicas_blk%gmask_hor_c1u(jc1u)) then
         if (nicas_blk%anisotropic) then
            ! Compute longitude/latitude differences
            dx = nicas_blk%lon_c1u(jc1u)-nicas_blk%lon_c1u(ic1u)
            dy = nicas_blk%lat_c1u(jc1u)-nicas_blk%lat_c1u(ic1u)
            call lonlatmod(dx,dy)
            dx = dx*cos(0.5*(nicas_blk%lat_c1u(ic1u)+nicas_blk%lat_c1u(jc1u)))
         else
            ! Compute horizontal distance
            call sphere_dist(nicas_blk%lon_c1u(ic1u),nicas_blk%lat_c1u(ic1u),nicas_blk%lon_c1u(jc1u),nicas_blk%lat_c1u(jc1u),dnb)
         end if

         do il1=1,nicas_blk%nl1
            ! Index
            il0 = nicas_blk%l1_to_l0(il1)

            do djl1=-1,1
               ! Index
               jl1 = max(1,min(il1+djl1,nicas_blk%nl1))
               jl0 = nicas_blk%l1_to_l0(jl1)

               ! Check valid arc for both levels
               if (nicas_blk%gmask_c1u(ic1u,il0).and.nicas_blk%gmask_c1u(jc1u,jl0).and.net_arc(ic1u,il1,j) &
 & .and.net_arc(ic1u,jl1,j)) then
                  ! Squared support radii
                  if (nicas_blk%anisotropic) then
                     H11 = 0.5*(nicas_blk%H11_c1u(ic1u,il1)+nicas_blk%H11_c1u(jc1u,jl1))
                     H22 = 0.5*(nicas_blk%H22_c1u(ic1u,il1)+nicas_blk%H22_c1u(jc1u,jl1))
                     H12 = 0.5*(nicas_blk%H12_c1u(ic1u,il1)+nicas_blk%H12_c1u(jc1u,jl1))
                     net_dnb(ic1u,il1,j,djl1) = H11*dx**2+H22*dy**2+2.0*H12*dx*dy
                     if (nicas_blk%horizontal) then
                        if (il0/=jl0) net_dnb(ic1u,il1,j,djl1) = net_dnb(ic1u,il1,j,djl1)+0.5*huge(1.0)
                     else
                        dz = nicas_blk%vunit_c1u(ic1u,il0)-nicas_blk%vunit_c1u(jc1u,jl0)
                        H33 = 0.5*(nicas_blk%H33_c1u(ic1u,il1)+nicas_blk%H33_c1u(jc1u,jl1))
                        net_dnb(ic1u,il1,j,djl1) = net_dnb(ic1u,il1,j,djl1)+H33*dz**2
                     end if
                     net_dnb(ic1u,il1,j,djl1) = sqrt(net_dnb(ic1u,il1,j,djl1))*gc2gau
                  else
                     disthsq = dnb**2
                     rhsq = 0.5*(nicas_blk%rh_c1u(ic1u,il1)**2+nicas_blk%rh_c1u(jc1u,jl1)**2)
                     if (rhsq>0.0) then
                        net_dnb(ic1u,il1,j,djl1) = disthsq/rhsq
                     elseif (disthsq>0.0) then
                        net_dnb(ic1u,il1,j,djl1) = 0.5*huge_real
                     else
                        net_dnb(ic1u,il1,j,djl1) = 0.0
                     end if
                     if (nicas_blk%horizontal) then
                        if (il0/=jl0) net_dnb(ic1u,il1,j,djl1) = net_dnb(ic1u,il1,j,djl1)+0.5*huge_real
                     else
                        distvsq = (nicas_blk%vunit_c1u(ic1u,il0)-nicas_blk%vunit_c1u(jc1u,jl0))**2
                        rvsq = 0.5*(nicas_blk%rv_c1u(ic1u,il1)**2+nicas_blk%rv_c1u(jc1u,jl1)**2)
                        if (rvsq>0.0) then
                           net_dnb(ic1u,il1,j,djl1) = net_dnb(ic1u,il1,j,djl1)+distvsq/rvsq
                        elseif (disthsq>0.0) then
                           net_dnb(ic1u,il1,j,djl1) = net_dnb(ic1u,il1,j,djl1)+0.5*huge_real
                        end if
                     end if
                     net_dnb(ic1u,il1,j,djl1) = sqrt(net_dnb(ic1u,il1,j,djl1))
                  end if
               end if
            end do
         end do
      end if
   end do

   ! Update
   if (nicas_blk%verbosity) call mpl%prog_print(ic1u)
end do
!$omp end parallel do
if (nicas_blk%verbosity) call mpl%prog_final

! Compute distances
write(mpl%info,'(a10,a)') '','Compute distances: '
if (nicas_blk%verbosity) call mpl%flush(.false.)
if (nicas_blk%verbosity) call mpl%prog_init(nicas_blk%nsbb)
!$omp parallel do schedule(static) private(isbb,isu,ic1u,il1,np,np_new,jc1u,jl1,k,kc1u,dkl1,kl1,disttest,add_to_front,jp), &
!$omp&                             firstprivate(distnorm,plist,plist_new)
do isbb=1,nicas_blk%nsbb
   ! Indices
   isu = nicas_blk%sbb_to_su(isbb)
   ic1u = nicas_blk%su_to_c1u(isu)
   il1 = nicas_blk%su_to_l1(isu)

   ! Allocation
   allocate(distnorm(nicas_blk%nc1u,nicas_blk%nl1))
   allocate(plist(nicas_blk%nc1u*nicas_blk%nl1,2))
   allocate(plist_new(nicas_blk%nc1u*nicas_blk%nl1,2))

   ! Initialize the front
   np = 1
   plist(1,1) = ic1u
   plist(1,2) = il1
   distnorm = 1.0
   distnorm(ic1u,il1) = 0.0

   do while (np>0)
      ! Propagate the front
      np_new = 0

      do ip=1,np
         ! Indices of the central point
         jc1u = plist(ip,1)
         jl1 = plist(ip,2)

         ! Loop over neighbors
         do k=1,mesh%rows(jc1u)%cols
            kc1u = mesh%rows(jc1u)%nodes(k)
            do dkl1=-1,1
               kl1 = max(1,min(jl1+dkl1,nicas_blk%nl1))
               if (nicas_blk%gmask_c2u(kc1u,kl1)) then
                  disttest = distnorm(jc1u,jl1)+net_dnb(jc1u,jl1,k,dkl1)
                  if (inf(disttest,1.0_kind_real)) then
                     ! Point is inside the support
                     if (inf(disttest,distnorm(kc1u,kl1))) then
                        ! Update distance
                        distnorm(kc1u,kl1) = disttest

                        ! Check if the point should be added to the front (avoid duplicates)
                        add_to_front = .true.
                        do jp=1,np_new
                           if ((plist_new(jp,1)==kc1u).and.(plist_new(jp,2)==kl1)) then
                              add_to_front = .false.
                              exit
                           end if
                        end do

                        if (add_to_front) then
                           ! Add point to the front
                           np_new = np_new+1
                           plist_new(np_new,1) = kc1u
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
      do ic1u=1,nicas_blk%nc1u
         if (supeq(distnorm(ic1u,il1),1.0_kind_real)) distnorm(ic1u,il1) = mpl%msv%valr
      end do
   end do
   call nicas_blk%distnorm(isbb)%pack(mpl,nicas_blk%nc1u,nicas_blk%nl1,distnorm)

   ! Release memory
   deallocate(distnorm)
   deallocate(plist)
   deallocate(plist_new)

   ! Update
   if (nicas_blk%verbosity) call mpl%prog_print(isbb)
end do
!$omp end parallel do
if (nicas_blk%verbosity) call mpl%prog_final

! Release memory
deallocate(net_arc)
deallocate(net_dnb)

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
integer :: nnmax,isu,ic1u,jc1u,il1,il0,j,jsu,jl0,jl1,ic1bb,isbb
integer :: nn(nicas_blk%nc1bb)
integer,allocatable :: nn_index(:,:)
real(kind_real) :: nnavg,rr,disthsq,distvsq,rhsq,rvsq
real(kind_real) :: dx,dy,dz,H11,H22,H33,H12
real(kind_real),allocatable :: distnorm(:,:),nn_dist(:,:)
logical,allocatable :: valid_arc(:,:,:)

! Count nearest neighbors
write(mpl%info,'(a10,a)') '','Count nearest neighbors: '
if (nicas_blk%verbosity) call mpl%flush(.false.)
if (nicas_blk%verbosity) call mpl%prog_init(nicas_blk%nc1bb)
do ic1bb=1,nicas_blk%nc1bb
   ! Indices
   ic1u = nicas_blk%c1bb_to_c1u(ic1bb)

   ! Research radius
   rr = sqrt(0.5*(maxval(nicas_blk%rh_c1u(ic1u,:))**2+nicas_blk%rhmax**2))

   ! Count nearest neighbors
   call nicas_blk%tree%count_nearest_neighbors(nicas_blk%lon_c1u(ic1u),nicas_blk%lat_c1u(ic1u),rr,nn(ic1bb))

   ! Update
   if (nicas_blk%verbosity) call mpl%prog_print(ic1bb)
end do
if (nicas_blk%verbosity) call mpl%prog_final

! Print results
if (nicas_blk%nc1bb>0) then
   nnmax = maxval(nn)
   nnavg = real(sum(nn),kind_real)/real(nicas_blk%nc1bb,kind_real)
else
   nnmax = 0
   nnavg = 0
end if
write(mpl%info,'(a10,a,f8.1,a,f5.1,a,i6,a,f5.1,a)') '','Average / maximum number of neighbors to find for this task: ',nnavg, &
 & ' (', nnavg/real(nicas_blk%nc1u,kind_real)*100.0,'%)  /',nnmax,' (', &
 & real(nnmax,kind_real)/real(nicas_blk%nc1u,kind_real)*100.0,'%)'
if (nicas_blk%verbosity) call mpl%flush

! Allocation
allocate(nn_index(nnmax,nicas_blk%nc1bb))
allocate(nn_dist(nnmax,nicas_blk%nc1bb))
allocate(valid_arc(nnmax,nicas_blk%nc1bb,nicas_blk%nl1))

! Find nearest neighbors
write(mpl%info,'(a10,a)') '','Find nearest neighbors: '
if (nicas_blk%verbosity) call mpl%flush(.false.)
if (nicas_blk%verbosity) call mpl%prog_init(nicas_blk%nc1bb)
do ic1bb=1,nicas_blk%nc1bb
   ! Indices
   ic1u = nicas_blk%c1bb_to_c1u(ic1bb)

   ! Find nearest neighbors
   call nicas_blk%tree%find_nearest_neighbors(nicas_blk%lon_c1u(ic1u),nicas_blk%lat_c1u(ic1u),nn(ic1bb), &
 & nn_index(1:nn(ic1bb),ic1bb),nn_dist(1:nn(ic1bb),ic1bb))

   ! Check arc validity
   do j=1,nn(ic1bb)
      jc1u = nn_index(j,ic1bb)
      do il1=1,nicas_blk%nl1
         il0 = nicas_blk%l1_to_l0(il1)
         valid_arc(j,ic1bb,il1) = (nicas_blk%gmask_c1u(ic1u,il0).and.nicas_blk%gmask_c1u(jc1u,il0))
         if (nam%mask_check.and.valid_arc(j,ic1bb,il1)) call geom%check_arc(mpl,il0,nicas_blk%lon_c1u(ic1u), &
 & nicas_blk%lat_c1u(ic1u),nicas_blk%lon_c1u(jc1u),nicas_blk%lat_c1u(jc1u),valid_arc(j,ic1bb,il1))
      end do
   end do

   ! Update
   if (nicas_blk%verbosity) call mpl%prog_print(ic1bb)
end do
if (nicas_blk%verbosity) call mpl%prog_final

! Compute distances
write(mpl%info,'(a10,a)') '','Compute distances: '
if (nicas_blk%verbosity) call mpl%flush(.false.)
if (nicas_blk%verbosity) call mpl%prog_init(nicas_blk%nsbb)
!$omp parallel do schedule(static) private(isbb,isu,ic1u,ic1bb,il1,il0,j,jl1,jl0,jsu,jc1u,disthsq,distvsq,rhsq,rvsq), &
!$omp&                             private(dx,dy,dz,H11,H22,H33,H12) firstprivate(distnorm)
do isbb=1,nicas_blk%nsbb
   ! Indices
   isu = nicas_blk%sbb_to_su(isbb)
   ic1u = nicas_blk%su_to_c1u(isu)
   ic1bb = nicas_blk%c1u_to_c1bb(ic1u)
   il1 = nicas_blk%su_to_l1(isu)
   il0 = nicas_blk%l1_to_l0(il1)

   ! Allocation
   allocate(distnorm(nicas_blk%nc1u,nicas_blk%nl1))

   ! Initialization
   distnorm = mpl%msv%valr

   ! Loop on nearest neighbors
   do j=1,nn(ic1bb)
      jc1u = nn_index(j,ic1bb)
      do jl1=1,nicas_blk%nl1
         jl0 = nicas_blk%l1_to_l0(jl1)
         if (nicas_blk%gmask_c2u(jc1u,jl1)) then
            ! Check valid arc for both levels
            if (valid_arc(j,ic1bb,il1).and.valid_arc(j,ic1bb,jl1)) then
               ! Normalized distance
               if (nicas_blk%anisotropic) then
                  dx = nicas_blk%lon_c1u(jc1u)-nicas_blk%lon_c1u(ic1u)
                  dy = nicas_blk%lat_c1u(jc1u)-nicas_blk%lat_c1u(ic1u)
                  call lonlatmod(dx,dy)
                  dx = dx*cos(0.5*(nicas_blk%lat_c1u(ic1u)+nicas_blk%lat_c1u(jc1u)))
                  H11 = 0.5*(nicas_blk%H11_c1u(ic1u,il1)+nicas_blk%H11_c1u(jc1u,jl1))
                  H22 = 0.5*(nicas_blk%H22_c1u(ic1u,il1)+nicas_blk%H22_c1u(jc1u,jl1))
                  H12 = 0.5*(nicas_blk%H12_c1u(ic1u,il1)+nicas_blk%H12_c1u(jc1u,jl1))
                  distnorm(jc1u,jl1) = H11*dx**2+H22*dy**2+2.0*H12*dx*dy
                  if (nicas_blk%horizontal) then
                     if (il0/=jl0) distnorm(jc1u,jl1) = distnorm(jc1u,jl1)+0.5*huge_real
                  else
                     dz = nicas_blk%vunit_c1u(ic1u,il0)-nicas_blk%vunit_c1u(jc1u,jl0)
                     H33 = 0.5*(nicas_blk%H33_c1u(ic1u,il1)+nicas_blk%H33_c1u(jc1u,jl1))
                     distnorm(jc1u,jl1) = distnorm(jc1u,jl1)+H33*dz**2
                  end if
                  distnorm(jc1u,jl1) = sqrt(distnorm(jc1u,jl1))*gc2gau
               else
                  disthsq = nn_dist(j,ic1bb)**2
                  rhsq = 0.5*(nicas_blk%rh_c1u(ic1u,il1)**2+nicas_blk%rh_c1u(jc1u,jl1)**2)
                  if (rhsq>0.0) then
                     distnorm(jc1u,jl1) = disthsq/rhsq
                  elseif (disthsq>0.0) then
                     distnorm(jc1u,jl1) = 0.5*huge_real
                  else
                     distnorm(jc1u,jl1) = 0.0
                  end if
                  if (nicas_blk%horizontal) then
                     if (il0/=jl0) distnorm(jc1u,jl1) = distnorm(jc1u,jl1)+0.5*huge_real
                  else
                     distvsq = (nicas_blk%vunit_c1u(ic1u,il0)-nicas_blk%vunit_c1u(jc1u,jl0))**2
                     rvsq = 0.5*(nicas_blk%rv_c1u(ic1u,il1)**2+nicas_blk%rv_c1u(jc1u,jl1)**2)
                     if (rvsq>0.0) then
                        distnorm(jc1u,jl1) = distnorm(jc1u,jl1)+distvsq/rvsq
                     elseif (distvsq>0.0) then
                        distnorm(jc1u,jl1) = distnorm(jc1u,jl1)+0.5*huge_real
                     end if
                  end if
                  distnorm(jc1u,jl1) = sqrt(distnorm(jc1u,jl1))
               end if
            end if
         end if
      end do
   end do

   ! Pack data
   do il1=1,nicas_blk%nl1
      do ic1u=1,nicas_blk%nc1u
         if (supeq(distnorm(ic1u,il1),1.0_kind_real)) distnorm(ic1u,il1) = mpl%msv%valr
      end do
   end do
   call nicas_blk%distnorm(isbb)%pack(mpl,nicas_blk%nc1u,nicas_blk%nl1,distnorm)

   ! Release memory
   deallocate(distnorm)

   ! Update
   if (nicas_blk%verbosity) call mpl%prog_print(isbb)
end do
!$omp end parallel do
if (nicas_blk%verbosity) call mpl%prog_final

! Release memory
deallocate(nn_index)
deallocate(nn_dist)
deallocate(valid_arc)

end subroutine nicas_blk_compute_convol_distance

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_convol_weights
! Purpose: compute convolution weights
!----------------------------------------------------------------------
subroutine nicas_blk_compute_convol_weights(nicas_blk,mpl,geom,ctmp)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(geom_type),intent(in) :: geom               ! Geometry
type(linop_type),intent(inout) :: ctmp           ! Convolution operator

! Local variables
integer :: n_s_max,ithread,isu,jc1u,jl1,jbd,jsu,isbb
integer :: c_n_s(mpl%nthread),c_nor_n_s(mpl%nthread)
real(kind_real) :: S_test,S_test_tot(nicas_blk%nsbb,mpl%nthread)
logical :: add_op
type(linop_type) :: c(mpl%nthread),c_nor(mpl%nthread)

! Allocation
n_s_max = 10*nint(real(geom%nc0u*geom%nl0)/real(mpl%nthread*mpl%nproc))
do ithread=1,mpl%nthread
   c(ithread)%n_s = n_s_max
   call c(ithread)%alloc
   c_nor(ithread)%n_s = n_s_max
   call c_nor(ithread)%alloc
end do

! Initialization
write(mpl%info,'(a10,a)') '','Compute weights: '
if (nicas_blk%verbosity) call mpl%flush(.false.)
if (nicas_blk%verbosity) call mpl%prog_init(nicas_blk%nsbb)
c_n_s = 0
c_nor_n_s = 0
S_test_tot = 0.0

! Compute weights
!$omp parallel do schedule(static) private(isbb,isu,ithread,jbd,jc1u,jl1,jsu,S_test,add_op)
do isbb=1,nicas_blk%nsbb
   ! Indices
   isu = nicas_blk%sbb_to_su(isbb)
   ithread = 1
!$ ithread = omp_get_thread_num()+1

   ! Count convolution operations
   do jbd=1,nicas_blk%distnorm(isbb)%nbd
      ! Indices
      jc1u = nicas_blk%distnorm(isbb)%bd_to_c1u(jbd)
      jl1 = nicas_blk%distnorm(isbb)%bd_to_l1(jbd)
      jsu = nicas_blk%c1ul1_to_su(jc1u,jl1)

      ! Fit function
      if (nicas_blk%anisotropic) then
         S_test = nicas_blk%Hcoef(isbb)%val(jbd)*fit_func(mpl,nicas_blk%distnorm(isbb)%val(jbd))
      else
         S_test = fit_func(mpl,nicas_blk%distnorm(isbb)%val(jbd))
      end if

      if (sup(S_test,S_inf)) then
         ! Store coefficient for convolution
         if (nicas_blk%lsqrt==1) then
            add_op = .false.
            if (nicas_blk%mpicom==1) then
               add_op = (nicas_blk%lcheck_sb(jsu).and.(isu<=jsu)).or.(.not.nicas_blk%lcheck_sb(jsu))
            elseif (nicas_blk%mpicom==2) then
               add_op = nicas_blk%lcheck_sa(isu).and.((nicas_blk%lcheck_sa(jsu).and.(isu<=jsu)) &
 & .or.(.not.nicas_blk%lcheck_sa(jsu)))
            end if
         else
            add_op = .true.
         end if
         if (add_op) call c(ithread)%add_op(c_n_s(ithread),isu,jsu,S_test)
         if (nicas_blk%smoother) S_test_tot(isbb,ithread) = S_test_tot(isbb,ithread)+S_test

         if (.not.nicas_blk%smoother) then
            ! Store coefficient for normalization
            add_op = nicas_blk%lcheck_sb(isu).and.((nicas_blk%lcheck_sb(jsu).and.(isu<=jsu)).or.(.not.nicas_blk%lcheck_sb(jsu)))
            if (add_op) call c_nor(ithread)%add_op(c_nor_n_s(ithread),isu,jsu,S_test)
         end if
      end if
   end do

   ! Update
   if (nicas_blk%verbosity) call mpl%prog_print(isbb)
end do
!$omp end parallel do
if (nicas_blk%verbosity) call mpl%prog_final

! Gather data from threads
call ctmp%gather(mpl,c_n_s,c)
if (nicas_blk%smoother) then
   nicas_blk%smoother_norm = sum(S_test_tot,dim=2)
else
   call nicas_blk%c_nor%gather(mpl,c_nor_n_s,c_nor)
endif

! Release memory
do ithread=1,mpl%nthread
   call c(ithread)%dealloc
   if (.not.nicas_blk%smoother) call c_nor(ithread)%dealloc
end do

end subroutine nicas_blk_compute_convol_weights

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_mpi_c
! Purpose: compute NICAS MPI distribution, halo C
!----------------------------------------------------------------------
subroutine nicas_blk_compute_mpi_c(nicas_blk,mpl)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data

! Local variables
integer :: isa,isb,isc,i_s,is,isu,jsu
integer,allocatable :: sc_to_s(:),su_to_sc(:),sc_nor_to_s(:),su_to_sc_nor(:)
logical,allocatable :: lcheck_sc(:),lcheck_sc_nor(:)
character(len=1024),parameter :: subr = 'nicas_blk_compute_mpi_c'

! Allocation
allocate(lcheck_sc(nicas_blk%nsu))
if (.not.nicas_blk%smoother) allocate(lcheck_sc_nor(nicas_blk%nsu))

! Define halo C
lcheck_sc = nicas_blk%lcheck_sb
do i_s=1,nicas_blk%c%n_s
   isu = nicas_blk%c%row(i_s)
   jsu = nicas_blk%c%col(i_s)
   lcheck_sc(isu) = .true.
   lcheck_sc(jsu) = .true.
end do
nicas_blk%nsc = count(lcheck_sc)
if (.not.nicas_blk%smoother) then
   lcheck_sc_nor = nicas_blk%lcheck_sb
   do i_s=1,nicas_blk%c_nor%n_s
      isu = nicas_blk%c_nor%row(i_s)
      jsu = nicas_blk%c_nor%col(i_s)
      lcheck_sc_nor(isu) = .true.
      lcheck_sc_nor(jsu) = .true.
   end do
   nicas_blk%nsc_nor = count(lcheck_sc_nor)
end if

! Check halos consistency
do isu=1,nicas_blk%nsu
   if (nicas_blk%lcheck_sa(isu).and.(.not.lcheck_sc(isu))) call mpl%abort(subr,'point in halo A but not in halo C')
   if (nicas_blk%lcheck_sb(isu).and.(.not.lcheck_sc(isu))) call mpl%abort(subr,'point in halo B but not in halo C')
end do

! Allocation
allocate(sc_to_s(nicas_blk%nsc))
allocate(nicas_blk%sc_to_su(nicas_blk%nsc))
allocate(su_to_sc(nicas_blk%nsu))
allocate(nicas_blk%sa_to_sc(nicas_blk%nsa))
allocate(nicas_blk%sb_to_sc(nicas_blk%nsb))
if (.not.nicas_blk%smoother) then
   allocate(sc_nor_to_s(nicas_blk%nsc_nor))
   allocate(su_to_sc_nor(nicas_blk%nsu))
   allocate(nicas_blk%sa_to_sc_nor(nicas_blk%nsa))
   allocate(nicas_blk%sb_to_sc_nor(nicas_blk%nsb))
end if

! Global-local conversion for halo C
su_to_sc = mpl%msv%vali
isc = 0
do isu=1,nicas_blk%nsu
   if (lcheck_sc(isu)) then
      isc = isc+1
      is = nicas_blk%su_to_s(isu)
      sc_to_s(isc) = is
      nicas_blk%sc_to_su(isc) = isu
      su_to_sc(isu) = isc
   end if
end do

if (.not.nicas_blk%smoother) then
   ! Global-local conversion for halo C (normalization)
   su_to_sc_nor = mpl%msv%vali
   isc = 0
   do isu=1,nicas_blk%nsu
      if (lcheck_sc_nor(isu)) then
         isc = isc+1
         is = nicas_blk%su_to_s(isu)
         sc_nor_to_s(isc) = is
         su_to_sc_nor(isu) = isc
      end if
   end do
end if

! Halos A-C conversions
do isa=1,nicas_blk%nsa
   isu = nicas_blk%sa_to_su(isa)
   isc = su_to_sc(isu)
   nicas_blk%sa_to_sc(isa) = isc
   if (.not.nicas_blk%smoother) then
      isc = su_to_sc_nor(isu)
      nicas_blk%sa_to_sc_nor(isa) = isc
   end if
end do
do isb=1,nicas_blk%nsb
   isu = nicas_blk%sb_to_su(isb)
   isc = su_to_sc(isu)
   nicas_blk%sb_to_sc(isb) = isc
   if (.not.nicas_blk%smoother) then
      isc = su_to_sc_nor(isu)
      nicas_blk%sb_to_sc_nor(isb) = isc
   end if
end do

! Local convolutions source and destination
nicas_blk%c%n_src = nicas_blk%nsc
nicas_blk%c%n_dst = nicas_blk%nsc
do i_s=1,nicas_blk%c%n_s
   if (nicas_blk%smoother) then
      isu = nicas_blk%c%row(i_s)
      isb = nicas_blk%su_to_sb(isu)
      if (nicas_blk%smoother_norm(isb)>0.0) then
         nicas_blk%c%S(i_s) = nicas_blk%c%S(i_s)/nicas_blk%smoother_norm(isb)
      else
        if (nicas_blk%c%S(i_s)>0.0) call mpl%abort(subr,'error in smoother_norm')
      end if
   end if
   nicas_blk%c%row(i_s) = su_to_sc(nicas_blk%c%row(i_s))
   nicas_blk%c%col(i_s) = su_to_sc(nicas_blk%c%col(i_s))
end do
if (.not.nicas_blk%smoother) then
   nicas_blk%c_nor%n_src = nicas_blk%nsc
   nicas_blk%c_nor%n_dst = nicas_blk%nsc
   do i_s=1,nicas_blk%c_nor%n_s
      nicas_blk%c_nor%row(i_s) = su_to_sc_nor(nicas_blk%c_nor%row(i_s))
      nicas_blk%c_nor%col(i_s) = su_to_sc_nor(nicas_blk%c_nor%col(i_s))
   end do
end if

! Setup communications
call nicas_blk%com_AC%setup(mpl,'com_AC',nicas_blk%nsa,nicas_blk%nsc,nicas_blk%ns,nicas_blk%sa_to_s,sc_to_s)
if (.not.nicas_blk%smoother) call nicas_blk%com_AC_nor%setup(mpl,'com_AC_nor',nicas_blk%nsa,nicas_blk%nsc_nor, &
 & nicas_blk%ns,nicas_blk%sa_to_s,sc_nor_to_s)

! Release memory
deallocate(lcheck_sc)
deallocate(sc_to_s)
deallocate(su_to_sc)
if (.not.nicas_blk%smoother) then
   deallocate(lcheck_sc_nor)
   deallocate(sc_nor_to_s)
   deallocate(su_to_sc_nor)
end if

end subroutine nicas_blk_compute_mpi_c

!----------------------------------------------------------------------
! Subroutine: nicas_blk_compute_internal_normalization
! Purpose: compute internal normalization
!----------------------------------------------------------------------
subroutine nicas_blk_compute_internal_normalization(nicas_blk,mpl)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl              ! MPI data

! Local variables
integer :: i_s,isa,isc,jsc
integer,allocatable :: inec(:)
real(kind_real),allocatable :: c_S(:,:),norm_a(:)

! Compute convolution inverse mapping
allocate(inec(nicas_blk%nsc_nor))
inec = 0
do i_s=1,nicas_blk%c_nor%n_s
   isc = nicas_blk%c_nor%col(i_s)
   jsc = nicas_blk%c_nor%row(i_s)
   if (jsc/=isc) then
      inec(isc) = inec(isc)+1
      inec(jsc) = inec(jsc)+1
   end if
end do
allocate(c_S(maxval(inec),nicas_blk%nsc_nor))
c_S = mpl%msv%valr
inec = 0
do i_s=1,nicas_blk%c_nor%n_s
   isc = nicas_blk%c_nor%col(i_s)
   jsc = nicas_blk%c_nor%row(i_s)
   if (jsc/=isc) then
      inec(isc) = inec(isc)+1
      c_S(inec(isc),isc) = nicas_blk%c_nor%S(i_s)
      inec(jsc) = inec(jsc)+1
      c_S(inec(jsc),jsc) = nicas_blk%c_nor%S(i_s)
   end if
end do

! Allocation
allocate(norm_a(nicas_blk%nsa))
allocate(nicas_blk%inorm(nicas_blk%nsc))
allocate(nicas_blk%inorm_nor(nicas_blk%nsc_nor))

! Compute normalization weights
do isa=1,nicas_blk%nsa
   ! Index
   isc = nicas_blk%sa_to_sc_nor(isa)

   ! Sum of squared values
   norm_a(isa) = 1.0+sum(c_S(1:inec(isc),isc)**2)

   ! Normalization factor
   norm_a(isa) = 1.0/sqrt(norm_a(isa))
end do

! Halo extension from zone A to zone C
call nicas_blk%com_AC%ext(mpl,norm_a,nicas_blk%inorm)
call nicas_blk%com_AC_nor%ext(mpl,norm_a,nicas_blk%inorm_nor)

! Release memory
deallocate(inec)
deallocate(c_S)
deallocate(norm_a)

end subroutine nicas_blk_compute_internal_normalization

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
integer :: il0i,i_s,ic1b,jc1b,isc,jsb,jsc,is,ic0a,il0,il1,ih,jv,nlr,ilr,ic,isc_add
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
      inec(isc) = inec(isc)+1
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
      inec(isc) = inec(isc)+1
      c_ind(inec(isc),isc) = jsc
      c_S(inec(isc),isc) = nicas_blk%c_nor%S(i_s)
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
      if (nicas_blk%verbosity) call mpl%flush(.false.)
      if (nicas_blk%verbosity) call mpl%prog_init(geom%nc0a)

      !$omp parallel do schedule(static) private(ic0a,nlr,isc_add,S_add,ih,ic1b,jv,il1,ilr,ic,isc,jsc,is), &
      !$omp&                             firstprivate(isc_list,S_list,S_list_tmp)
      do ic0a=1,geom%nc0a
         ! Index
         if (geom%gmask_c0a(ic0a,il0)) then
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
               do jv=1,inev(il0)
                  il1 = v_col(jv,il0)
                  if ((nicas_blk%l1_to_l0(il1)>=nicas_blk%vbot(ic1b)).and.(nicas_blk%l1_to_l0(il1)<=nicas_blk%vtop(ic1b))) then
                     do is=1,ines(ic1b,il1)
                        isc_add = s_col(is,ic1b,il1)
                        S_add = h_S(ih,ic0a,il0i)*v_S(jv,ic1b,il0)*s_S(is,ic1b,il1)*nicas_blk%inorm_nor(isc_add)
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
            if (nicas_blk%verbosity) call mpl%prog_print(ic0a)

            ! Release memory
            deallocate(isc_list)
            deallocate(S_list)
            deallocate(S_list_tmp)
         end if
      end do
      !$omp end parallel do
      if (nicas_blk%verbosity) call mpl%prog_final
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
subroutine nicas_blk_compute_grids(nicas_blk,nam)

implicit none

! Passed variables
class(nicas_blk_type),intent(inout) :: nicas_blk ! NICAS data block
type(nam_type),intent(in) :: nam                 ! Namelist

! Local variables
integer :: isa,isb,isc,isu,ic1u,il1,il0

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
      isu = nicas_blk%sa_to_su(isa)
      ic1u = nicas_blk%su_to_c1u(isu)
      il1 = nicas_blk%su_to_l1(isu)
      il0 = nicas_blk%l1_to_l0(il1)
      nicas_blk%lon_sa(isa) = nicas_blk%lon_c1u(ic1u)*rad2deg
      nicas_blk%lat_sa(isa) = nicas_blk%lat_c1u(ic1u)*rad2deg
      nicas_blk%lev_sa(isa) = nam%levs(il0)
   end do
end if
if (nicas_blk%nsb>0) then
   do isb=1,nicas_blk%nsb
      isu = nicas_blk%sb_to_su(isb)
      ic1u = nicas_blk%su_to_c1u(isu)
      il1 = nicas_blk%su_to_l1(isu)
      il0 = nicas_blk%l1_to_l0(il1)
      nicas_blk%lon_sb(isb) = nicas_blk%lon_c1u(ic1u)*rad2deg
      nicas_blk%lat_sb(isb) = nicas_blk%lat_c1u(ic1u)*rad2deg
      nicas_blk%lev_sb(isb) = nam%levs(il0)
   end do
end if
if (nicas_blk%nsc>0) then
   do isc=1,nicas_blk%nsc
      isu = nicas_blk%sc_to_su(isc)
      ic1u = nicas_blk%su_to_c1u(isu)
      il1 = nicas_blk%su_to_l1(isu)
      il0 = nicas_blk%l1_to_l0(il1)
      nicas_blk%lon_sc(isc) = nicas_blk%lon_c1u(ic1u)*rad2deg
      nicas_blk%lat_sc(isc) = nicas_blk%lat_c1u(ic1u)*rad2deg
      nicas_blk%lev_sc(isc) = nam%levs(il0)
   end do
end if

end subroutine nicas_blk_compute_grids

!----------------------------------------------------------------------
! Subroutine: nicas_blk_apply
! Purpose: apply NICAS method
!----------------------------------------------------------------------
subroutine nicas_blk_apply(nicas_blk,mpl,geom,fld)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk            ! NICAS data block
type(mpl_type),intent(inout) :: mpl                      ! MPI data
type(geom_type),intent(in) :: geom                       ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0) ! Field

! Local variables
integer :: il0
real(kind_real) :: sums(geom%nl0),sume(geom%nl0)
real(kind_real) :: alpha_a(nicas_blk%nsa),alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

if (nicas_blk%smoother) then
   ! Save global sum for each level
   call mpl%f_comm%allreduce(sum(fld,dim=1,mask=geom%gmask_c0a),sums,fckit_mpi_sum())
else
   ! Normalization
   fld = fld*nicas_blk%norm
end if

! Adjoint interpolation
call nicas_blk%apply_interp_ad(mpl,geom,fld,alpha_b)

! Communication
if (nicas_blk%mpicom==1) then
   ! Initialization
   alpha_c = 0.0

   ! Copy zone B into zone C
   alpha_c(nicas_blk%sb_to_sc) = alpha_b
elseif (nicas_blk%mpicom==2) then
   ! Halo reduction from zone B to zone A
   call nicas_blk%com_AB%red(mpl,alpha_b,alpha_a)

   ! Initialization
   alpha_c = 0.0

   ! Copy zone A into zone C
   alpha_c(nicas_blk%sa_to_sc) = alpha_a
end if

! Internal normalization
if (.not.nicas_blk%smoother) alpha_c = alpha_c*nicas_blk%inorm

! Convolution
call nicas_blk%apply_convol(mpl,alpha_c)

! Internal normalization
if (.not.nicas_blk%smoother) alpha_c = alpha_c*nicas_blk%inorm

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(mpl,alpha_c,alpha_a)

! Halo extension from zone A to zone B
call nicas_blk%com_AB%ext(mpl,alpha_a,alpha_b)

! Interpolation
call nicas_blk%apply_interp(mpl,geom,alpha_b,fld)

if (nicas_blk%smoother) then
   ! Reset global sum for each level
   call mpl%f_comm%allreduce(sum(fld,dim=1,mask=geom%gmask_c0a),sume,fckit_mpi_sum())
   do il0=1,geom%nl0
      fld(:,il0) = fld(:,il0)*sums(il0)/sume(il0)
   end do
else
   ! Normalization
   fld = fld*nicas_blk%norm
end if

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

! Internal normalization
if (.not.nicas_blk%smoother) alpha_c = alpha_c*nicas_blk%inorm

! Halo reduction from zone C to zone A
call nicas_blk%com_AC%red(mpl,alpha_c,alpha_a)

! Halo extension from zone A to zone B
call nicas_blk%com_AB%ext(mpl,alpha_a,alpha_b)

! Interpolation
call nicas_blk%apply_interp(mpl,geom,alpha_b,fld)

! Normalization
if (.not.nicas_blk%smoother) fld = fld*nicas_blk%norm

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
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0),alpha_b(nicas_blk%nsb),alpha_c(nicas_blk%nsc)

! Initialization
fld_tmp = fld

! Normalization
if (.not.nicas_blk%smoother) fld_tmp = fld_tmp*nicas_blk%norm

! Adjoint interpolation
call nicas_blk%apply_interp_ad(mpl,geom,fld_tmp,alpha_b)

! Halo reduction from zone B to zone A
call nicas_blk%com_AB%red(mpl,alpha_b,alpha)

! Initialization
alpha_c = 0.0

! Copy zone A into zone C
alpha_c(nicas_blk%sa_to_sc) = alpha

! Internal normalization
if (.not.nicas_blk%smoother) alpha_c = alpha_c*nicas_blk%inorm

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

! Horizontal interpolation
call nicas_blk%apply_interp_h_ad(mpl,geom,fld,delta)

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

! Local variables
real(kind_real),allocatable :: alpha_tmp(:)

if (nicas_blk%smoother) then
   ! Allocation
   allocate(alpha_tmp(nicas_blk%nsc))

   ! Apply linear operator
   alpha_tmp = alpha
   alpha = 0.0
   call nicas_blk%c%apply_ad(mpl,alpha_tmp,alpha)

   ! Release memory
   deallocate(alpha_tmp)
else
   ! Apply linear operator, symmetric
   call nicas_blk%c%apply_sym(mpl,alpha)
end if

end subroutine nicas_blk_apply_convol

!----------------------------------------------------------------------
! Subroutine: nicas_blk_test_adjoint
! Purpose: test NICAS adjoint accuracy
!----------------------------------------------------------------------
subroutine nicas_blk_test_adjoint(nicas_blk,mpl,rng,geom)

implicit none

! Passed variables
class(nicas_blk_type),intent(in) :: nicas_blk ! NICAS data block
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(rng_type),intent(inout) :: rng           ! Random number generator
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
if (nicas_blk%verbosity) call mpl%flush

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
if (nicas_blk%verbosity) call mpl%flush

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
if (nicas_blk%verbosity) call mpl%flush

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
if (nicas_blk%verbosity) call mpl%flush

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
if (nicas_blk%verbosity) call mpl%flush

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
if (nicas_blk%verbosity) call mpl%flush

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
if (nicas_blk%verbosity) call mpl%flush

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
if (nicas_blk%verbosity) call mpl%flush

! Initialization
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld2_save)
fld1 = fld1_save
fld2 = fld2_save

! Adjoint test
if (nicas_blk%lsqrt==1) then
   call nicas_blk%apply_from_sqrt(mpl,geom,fld1)
   call nicas_blk%apply_from_sqrt(mpl,geom,fld2)
else
   call nicas_blk%apply(mpl,geom,fld1)
   call nicas_blk%apply(mpl,geom,fld2)
end if

! Print result
call mpl%dot_prod(fld1,fld2_save,sum1)
call mpl%dot_prod(fld2,fld1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','NICAS adjoint test:                       ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
if (nicas_blk%verbosity) call mpl%flush

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
if (nicas_blk%lsqrt==1) then
   call nicas_blk%apply_from_sqrt(mpl,geom,fld)
else
   call nicas_blk%apply(mpl,geom,fld)
end if

! Write field
filename = trim(nam%prefix)//'_dirac'
call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)
call io%fld_write(mpl,nam,geom,filename,'nicas_blk',fld,bpar%blockname(ib))

! Print results
write(mpl%info,'(a7,a)') '','Values at dirac points:'
if (nicas_blk%verbosity) call mpl%flush
do idir=1,geom%ndir
   if (geom%iprocdir(idir)==mpl%myproc) val = fld(geom%ic0adir(idir),geom%il0dir(idir))
   call mpl%f_comm%broadcast(val,geom%iprocdir(idir)-1)
   if (mpl%msv%isnot(val)) then
      write(mpl%info,'(a10,f6.1,a,f6.1,a,f10.7)') '',geom%londir(idir)*rad2deg,' / ',geom%latdir(idir)*rad2deg,': ',val
      if (nicas_blk%verbosity) call mpl%flush
   else
      write(mpl%info,'(a10,f6.1,a,f6.1,a)') '',geom%londir(idir)*rad2deg,' / ',geom%latdir(idir)*rad2deg,': missing value'
      if (nicas_blk%verbosity) call mpl%flush
   end if
end do
write(mpl%info,'(a7,a)') '','Min - max: '
if (nicas_blk%verbosity) call mpl%flush
do il0=1,geom%nl0
   call mpl%f_comm%allreduce(minval(fld(:,il0),mask=geom%gmask_c0a(:,il0)),valmin_tot,fckit_mpi_min())
   call mpl%f_comm%allreduce(maxval(fld(:,il0),mask=geom%gmask_c0a(:,il0)),valmax_tot,fckit_mpi_max())
   if (mpl%msv%isnot(valmin_tot).or.mpl%msv%isnot(valmax_tot)) then
      write(mpl%info,'(a10,a,i3,a,f10.7,a,f10.7)') '','Level ',nam%levs(il0),': ',valmin_tot,' - ',valmax_tot
      if (nicas_blk%verbosity) call mpl%flush
   else
      write(mpl%info,'(a10,a,i3,a)') '','Level ',nam%levs(il0),': missing values'
      if (nicas_blk%verbosity) call mpl%flush
   end if
end do

! End associate
end associate

end subroutine nicas_blk_test_dirac

end module type_nicas_blk
