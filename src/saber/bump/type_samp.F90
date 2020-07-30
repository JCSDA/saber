!----------------------------------------------------------------------
! Module: type_samp
! Purpose: sampling derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_samp

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_status
use netcdf
!$ use omp_lib
use tools_const, only: pi,req,reqkm,deg2rad,rad2deg
use tools_func, only: sphere_dist,fit_func
use tools_kinds, only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use tools_repro, only: eq,inf
use tools_samp, only: initialize_sampling
use type_com, only: com_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_tree, only: tree_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

! Sampling derived type
type samp_type
   ! Parameters
   character(len=1024) :: name                      ! Sampling name
   logical :: sc2                                   ! Subset Sc2 flag
   logical :: sc3                                   ! Subset Sc3 flag

   ! Subset Sc0
   logical,allocatable :: smask_c0u(:,:)            ! Mask on subset Sc0, universe
   logical,allocatable :: smask_hor_c0u(:)          ! Union of horizontal masks on subset Sc0, universe
   logical,allocatable :: smask_c0a(:,:)            ! Mask on subset Sc0, halo A
   logical,allocatable :: smask_hor_c0a(:)          ! Union of horizontal masks on subset Sc0, halo A
   integer,allocatable :: nc0_smask(:)              ! Horizontal mask size on subset Sc0
   integer :: nc0c                                  ! Number of points in subset Sc0, halo C

   ! Subset Sc1
   integer,allocatable :: c1_to_c0(:)               ! Subset Sc1, global, to subset Sc0, global
   integer,allocatable :: c1_to_proc(:)             ! Subset Sc1, global to processor
   integer :: nc1u                                  ! Number of points in subset Sc1, universe
   integer,allocatable :: c1u_to_c1(:)              ! Subset Sc1, universe, to subset Sc1, global
   integer,allocatable :: c1_to_c1u(:)              ! Subset Sc1, global, to subset Sc1, universe
   integer,allocatable :: c1u_to_c1a(:)             ! Subset Sc1, universe, to subset Sc1, halo A
   integer,allocatable :: c1u_to_c0u(:)             ! Subset Sc1, universe, to subset Sc0, universe
   real(kind_real),allocatable :: lon_c1u(:)        ! Longitudes on subset Sc1, universe
   real(kind_real),allocatable :: lat_c1u(:)        ! Latitudes on subset Sc1, universe
   logical,allocatable :: smask_c1u(:,:)            ! Mask on subset Sc1, universe
   logical,allocatable :: smask_hor_c1u(:)          ! Union of horizontal masks on subset Sc1, universe
   integer :: nc1a                                  ! Number of points in subset Sc1, halo A
   integer,allocatable :: c1a_to_c1(:)              ! Subset Sc1, halo A, to subset Sc1, global
   integer,allocatable :: c1a_to_c1u(:)             ! Subset Sc1, halo A, to subset Sc1, universe
   integer,allocatable :: c1a_to_c0a(:)             ! Subset Sc1, halo A, to subset Sc0, halo A
   logical,allocatable :: c1al0_check(:,:)          ! Mask boundaries checking activation
   real(kind_real),allocatable :: lon_c1a(:)        ! Longitudes on subset Sc1, halo A
   real(kind_real),allocatable :: lat_c1a(:)        ! Latitudes on subset Sc1, halo A
   logical,allocatable :: smask_c1a(:,:)            ! Mask on subset Sc1, halo A
   logical,allocatable :: smask_hor_c1a(:)          ! Union of horizontal masks on subset Sc1, halo A
   integer,allocatable :: c1ac3_to_c0u(:,:)         ! Subsets Sc1 and Sc3, halo A, to subset Sc0, universe
   logical,allocatable :: smask_c1ac3(:,:,:)        ! Mask on subset Sc1 and Sc3, halo A
   logical,allocatable :: smask_c1dc3(:,:,:)        ! Mask on subset Sc1 and Sc3, halo D
   integer,allocatable :: c1a_to_c0c(:)             ! Subset Sc1, halo A, to subset Sc0, halo C
   integer,allocatable :: c1ac3_to_c0c(:,:)         ! Subsets Sc1 and Sc3, halo A, to subset Sc0, halo C
   integer :: nc1d                                  ! Number of points in subset Sc1, halo D
   integer,allocatable :: c1d_to_c1u(:)             ! Subset Sc1, halo D to universe
   integer :: nc1e                                  ! Number of points in subset Sc1, halo E
   integer,allocatable :: c1e_to_c1u(:)             ! Subset Sc1, halo E to universe

   ! Subset Sc2
   integer,allocatable :: c2_to_c1(:)               ! Subset Sc2, global, to subset Sc1, global
   integer,allocatable :: c2_to_proc(:)             ! Subset Sc2, global to processor
   integer :: nc2u                                  ! Number of points in subset Sc2, universe
   integer,allocatable :: c2u_to_c2(:)              ! Subset Sc2, universe, to subset Sc2, global
   integer,allocatable :: c2_to_c2u(:)              ! Subset Sc2, global, to subset Sc2, universe
   integer,allocatable :: c2u_to_c1u(:)             ! Subset Sc2, universe to subset Sc1, universe
   integer,allocatable :: c2u_to_c0u(:)             ! Subset Sc2, universe to subset Sc0, universe
   real(kind_real),allocatable :: lon_c2u(:)        ! Longitudes on subset Sc2, universe
   real(kind_real),allocatable :: lat_c2u(:)        ! Latitudes on subset Sc2, universe
   logical,allocatable :: smask_c2u(:,:)            ! Mask on subset Sc2, universe
   logical,allocatable :: smask_hor_c2u(:)          ! Union of horizontal masks on subset Sc2, universe
   integer,allocatable :: proc_to_nc2a(:)           ! Processor to subset Sc2 size, halo A
   integer,allocatable :: proc_to_c2_offset(:)      ! Processor to offset on subset Sc2
   integer :: nc2a                                  ! Number of points in subset Sc2, halo A
   integer,allocatable :: c2a_to_c2(:)              ! Subset Sc2, halo A, to subset Sc2, global
   integer,allocatable :: c2a_to_c2u(:)             ! Subset Sc2, halo A, to subset Sc2, universe
   integer,allocatable :: c2a_to_c1a(:)             ! Subset Sc2, halo A, to subset Sc1, halo A
   integer,allocatable :: c2a_to_c0a(:)             ! Subset Sc2, halo A, to subset Sc0, halo A
   real(kind_real),allocatable :: lon_c2a(:)        ! Longitudes on subset Sc2, halo A
   real(kind_real),allocatable :: lat_c2a(:)        ! Latitudes on subset Sc2, halo A
   logical,allocatable :: smask_c2a(:,:)            ! Mask on subset Sc2, halo A
   logical,allocatable :: smask_hor_c2a(:)          ! Union of horizontal masks on subset Sc2, halo A
   integer :: nc2b                                  ! Number of points in subset Sc2, halo B
   integer,allocatable :: c2b_to_c2u(:)             ! Subset Sc2, halo B to universe

   ! Local data
   logical,allocatable :: vbal_mask(:,:)            ! Vertical balance mask
   logical,allocatable :: local_mask(:,:)           ! Local mask
   integer,allocatable :: nn_c2a_index(:,:)         ! Nearest diagnostic neighbors from diagnostic points
   real(kind_real),allocatable :: nn_c2a_dist(:,:)  ! Nearest diagnostic neighbors distance from diagnostic points

   ! Forced points
   integer,allocatable :: ldwv_to_proc(:)           ! Local diagnostics profiles to task
   integer,allocatable :: ldwv_to_c0a(:)            ! Local diagnostics profiles to subset Sc0, halo A

   ! Sampling mesh
   type(mesh_type) :: mesh                          ! Sampling mesh

   ! Advection
   integer,allocatable :: adv_nn(:)                 ! Number of nearest neighbors inside search radius
   integer :: adv_nnmax                             ! Maximum number of nearest neighbors inside search radius
   integer,allocatable :: adv_nn_index(:,:)         ! Index of nearest neighbors inside search radius
   real(kind_real),allocatable :: adv_lon(:,:,:)    ! Interpolated advected longitude
   real(kind_real),allocatable :: adv_lat(:,:,:)    ! Interpolated advected latitude

   ! Interpolations
   type(linop_type),allocatable :: h(:)             ! Horizontal interpolation from Sc2 to Sc0 (local)
   type(linop_type),allocatable :: d(:,:)           ! Advection interpolation

   ! Communications
   type(com_type) :: com_AB                         ! Communication between halos A and B
   type(com_type) :: com_AC                         ! Communication between halos A and C
   type(com_type) :: com_AD                         ! Communication between halos A and D (diagnostic)
   type(com_type) :: com_AE                         ! Communication between halos A and E (vertical balance)
   type(com_type) :: com_AU                         ! Communication between halo A and universe on subset Sc2
contains
   procedure :: samp_alloc_mask
   procedure :: samp_alloc_other
   generic :: alloc => samp_alloc_mask,samp_alloc_other
   procedure :: partial_dealloc => samp_partial_dealloc
   procedure :: dealloc => samp_dealloc
   procedure :: read => samp_read
   procedure :: write => samp_write
   procedure :: write_grids => samp_write_grids
   procedure :: samp_setup_1
   procedure :: samp_setup_2
   generic :: setup => samp_setup_1,samp_setup_2
   procedure :: compute_mask => samp_compute_mask
   procedure :: compute_c1 => samp_compute_c1
   procedure :: compute_mpi_c1a => samp_compute_mpi_c1a
   procedure :: compute_c3 => samp_compute_c3
   procedure :: check_mask => samp_check_mask
   procedure :: compute_c2 => samp_compute_c2
   procedure :: compute_mpi_c2a => samp_compute_mpi_c2a
   procedure :: compute_mpi_c2b => samp_compute_mpi_c2b
   procedure :: compute_mesh_c2 => samp_compute_mesh_c2
   procedure :: compute_mpi_c => samp_compute_mpi_c
   procedure :: compute_mpi_d => samp_compute_mpi_d
   procedure :: compute_mpi_e => samp_compute_mpi_e
   procedure :: diag_filter => samp_diag_filter
   procedure :: diag_fill => samp_diag_fill
end type samp_type

private
public :: samp_type

contains

!----------------------------------------------------------------------
! Subroutine: samp_alloc_mask
! Purpose: allocation for mask
!----------------------------------------------------------------------
subroutine samp_alloc_mask(samp,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(geom_type),intent(in) :: geom     ! Geometry

! Allocation
allocate(samp%smask_c0u(geom%nc0u,geom%nl0))
allocate(samp%smask_c0a(geom%nc0a,geom%nl0))
allocate(samp%smask_hor_c0u(geom%nc0u))
allocate(samp%smask_hor_c0a(geom%nc0a))
allocate(samp%nc0_smask(0:geom%nl0))

end subroutine samp_alloc_mask

!----------------------------------------------------------------------
! Subroutine: samp_alloc_other
! Purpose: allocation for other variables
!----------------------------------------------------------------------
subroutine samp_alloc_other(samp,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Initialization
samp%sc2 = (trim(samp%name)=='vbal').or.(trim(samp%name)=='lct') &
 & .or.((trim(samp%name)=='hdiag').and.(nam%local_diag.or.nam%adv_diag))
samp%sc3 = (trim(samp%name)=='hdiag').or.(trim(samp%name)=='lct')

! Allocation
allocate(samp%c1_to_c0(nam%nc1))
allocate(samp%c1_to_proc(nam%nc1))
allocate(samp%c1_to_c1u(nam%nc1))
if (samp%sc2) then
   allocate(samp%c2_to_c1(nam%nc2))
   allocate(samp%c2_to_proc(nam%nc2))
   allocate(samp%c2_to_c2u(nam%nc2))
end if
if (nam%adv_diag) then
   allocate(samp%adv_lon(geom%nc0a,geom%nl0,nam%nts))
   allocate(samp%adv_lat(geom%nc0a,geom%nl0,nam%nts))
end if

end subroutine samp_alloc_other

!----------------------------------------------------------------------
! Subroutine: samp_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine samp_partial_dealloc(samp)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling

! Release memory
if (allocated(samp%smask_c0u)) deallocate(samp%smask_c0u)
if (allocated(samp%smask_hor_c0u)) deallocate(samp%smask_hor_c0u)
if (allocated(samp%smask_c0a)) deallocate(samp%smask_c0a)
if (allocated(samp%smask_hor_c0a)) deallocate(samp%smask_hor_c0a)
if (allocated(samp%nc0_smask)) deallocate(samp%nc0_smask)
if (allocated(samp%c1_to_c0)) deallocate(samp%c1_to_c0)
if (allocated(samp%c1_to_proc)) deallocate(samp%c1_to_proc)
if (allocated(samp%c1u_to_c1)) deallocate(samp%c1u_to_c1)
if (allocated(samp%c1_to_c1u)) deallocate(samp%c1_to_c1u)
if (allocated(samp%c1u_to_c1a)) deallocate(samp%c1u_to_c1a)
if (allocated(samp%c1u_to_c0u)) deallocate(samp%c1u_to_c0u)
if (allocated(samp%lon_c1u)) deallocate(samp%lon_c1u)
if (allocated(samp%lat_c1u)) deallocate(samp%lat_c1u)
if (allocated(samp%smask_hor_c1u)) deallocate(samp%smask_hor_c1u)
if (allocated(samp%c1a_to_c1)) deallocate(samp%c1a_to_c1)
if (allocated(samp%c1a_to_c1u)) deallocate(samp%c1a_to_c1u)
if (allocated(samp%c1al0_check)) deallocate(samp%c1al0_check)
if (allocated(samp%lon_c1a)) deallocate(samp%lon_c1a)
if (allocated(samp%lat_c1a)) deallocate(samp%lat_c1a)
if (allocated(samp%smask_hor_c1a)) deallocate(samp%smask_hor_c1a)
if (allocated(samp%c2_to_c1)) deallocate(samp%c2_to_c1)
if (allocated(samp%c2_to_proc)) deallocate(samp%c2_to_proc)
if (allocated(samp%c2_to_c2u)) deallocate(samp%c2_to_c2u)
if (allocated(samp%c2u_to_c1u)) deallocate(samp%c2u_to_c1u)
if (allocated(samp%c2u_to_c0u)) deallocate(samp%c2u_to_c0u)
if (allocated(samp%lon_c2u)) deallocate(samp%lon_c2u)
if (allocated(samp%lat_c2u)) deallocate(samp%lat_c2u)
if (allocated(samp%smask_c2u)) deallocate(samp%smask_c2u)
if (allocated(samp%smask_hor_c2u)) deallocate(samp%smask_hor_c2u)
if (allocated(samp%lon_c2a)) deallocate(samp%lon_c2a)
if (allocated(samp%lat_c2a)) deallocate(samp%lat_c2a)
if (allocated(samp%smask_hor_c2a)) deallocate(samp%smask_hor_c2a)
if (allocated(samp%c2b_to_c2u)) deallocate(samp%c2b_to_c2u)
if (allocated(samp%ldwv_to_proc)) deallocate(samp%ldwv_to_proc)

end subroutine samp_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: samp_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine samp_dealloc(samp)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling

! Local variables
integer :: il0,its

! Release memory
call samp%partial_dealloc
if (allocated(samp%c1a_to_c0c)) deallocate(samp%c1a_to_c0c)
if (allocated(samp%c1ac3_to_c0c)) deallocate(samp%c1ac3_to_c0c)
if (allocated(samp%c1a_to_c0a)) deallocate(samp%c1a_to_c0a)
if (allocated(samp%smask_c1a)) deallocate(samp%smask_c1a)
if (allocated(samp%c1ac3_to_c0u)) deallocate(samp%c1ac3_to_c0u)
if (allocated(samp%smask_c1ac3)) deallocate(samp%smask_c1ac3)
if (allocated(samp%smask_c1dc3)) deallocate(samp%smask_c1dc3)
if (allocated(samp%smask_c1u)) deallocate(samp%smask_c1u)
if (allocated(samp%c1d_to_c1u)) deallocate(samp%c1d_to_c1u)
if (allocated(samp%c1e_to_c1u)) deallocate(samp%c1e_to_c1u)
if (allocated(samp%c2u_to_c2)) deallocate(samp%c2u_to_c2)
if (allocated(samp%proc_to_nc2a)) deallocate(samp%proc_to_nc2a)
if (allocated(samp%proc_to_c2_offset)) deallocate(samp%proc_to_c2_offset)
if (allocated(samp%c2a_to_c2)) deallocate(samp%c2a_to_c2)
if (allocated(samp%c2a_to_c2u)) deallocate(samp%c2a_to_c2u)
if (allocated(samp%c2a_to_c1a)) deallocate(samp%c2a_to_c1a)
if (allocated(samp%c2a_to_c0a)) deallocate(samp%c2a_to_c0a)
if (allocated(samp%smask_c2a)) deallocate(samp%smask_c2a)
if (allocated(samp%vbal_mask)) deallocate(samp%vbal_mask)
if (allocated(samp%local_mask)) deallocate(samp%local_mask)
if (allocated(samp%nn_c2a_index)) deallocate(samp%nn_c2a_index)
if (allocated(samp%nn_c2a_dist)) deallocate(samp%nn_c2a_dist)
if (allocated(samp%ldwv_to_c0a)) deallocate(samp%ldwv_to_c0a)
call samp%mesh%dealloc
if (allocated(samp%adv_lon)) deallocate(samp%adv_lon)
if (allocated(samp%adv_lat)) deallocate(samp%adv_lat)
if (allocated(samp%h)) then
   do il0=1,size(samp%h)
      call samp%h(il0)%dealloc
   end do
   deallocate(samp%h)
end if
if (allocated(samp%d)) then
   do its=1,size(samp%d,2)
      do il0=1,size(samp%d,1)
         call samp%d(il0,its)%dealloc
      end do
   end do
   deallocate(samp%d)
end if
call samp%com_AB%dealloc
call samp%com_AC%dealloc
call samp%com_AD%dealloc
call samp%com_AE%dealloc
call samp%com_AU%dealloc

end subroutine samp_dealloc

!----------------------------------------------------------------------
! Subroutine: samp_read
! Purpose: read
!----------------------------------------------------------------------
subroutine samp_read(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: il0,ic1a,ic1u,jc3,grid_hash
integer :: ncid,grpid,c1a_to_c0a_id,smask_c1u_id,c1ac3_to_c0u_id,smask_c1ac3_id,c2u_to_c1u_id,c2u_to_c0u_id
integer,allocatable :: smask_c1uint(:,:),smask_c1ac3int(:,:,:)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'samp_read'

! Open file
write(filename,'(a,a,i6.6,a,i6.6)') trim(nam%prefix),'_sampling_',mpl%nproc,'-',mpl%myproc
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))

! Check grid hash
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'grid_hash',grid_hash))
if (grid_hash/=geom%grid_hash) call mpl%abort(subr,'wrong grid hash')

! Get group
call mpl%ncerr(subr,nf90_inq_grp_ncid(ncid,samp%name,grpid))

! Get or check dimensions
call mpl%nc_dim_check(subr,ncid,'nl0',geom%nl0)
samp%nc1a = mpl%nc_dim_inquire(subr,grpid,'nc1a')
samp%nc1u = mpl%nc_dim_inquire(subr,grpid,'nc1u')
if (samp%sc2) samp%nc2u = mpl%nc_dim_inquire(subr,grpid,'nc2u')
if (samp%sc3) call mpl%nc_dim_check(subr,ncid,'nc3',nam%nc3)

! Get variables
call mpl%ncerr(subr,nf90_inq_varid(grpid,'c1a_to_c0a',c1a_to_c0a_id))
call mpl%ncerr(subr,nf90_inq_varid(grpid,'smask_c1u',smask_c1u_id))
if (samp%sc2) then
   call mpl%ncerr(subr,nf90_inq_varid(grpid,'c2u_to_c1u',c2u_to_c1u_id))
   call mpl%ncerr(subr,nf90_inq_varid(grpid,'c2u_to_c0u',c2u_to_c0u_id))
end if
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_inq_varid(grpid,'c1ac3_to_c0u',c1ac3_to_c0u_id))
   call mpl%ncerr(subr,nf90_inq_varid(grpid,'smask_c1ac3',smask_c1ac3_id))
end if

! Allocation
allocate(samp%c1a_to_c0a(samp%nc1a))
allocate(smask_c1uint(samp%nc1u,geom%nl0))
allocate(samp%smask_c1u(samp%nc1u,geom%nl0))
if (samp%sc2) then
   allocate(samp%c2u_to_c1u(samp%nc2u))
   allocate(samp%c2u_to_c0u(samp%nc2u))
end if
if (samp%sc3) then
   allocate(samp%c1ac3_to_c0u(samp%nc1a,nam%nc3))
   allocate(smask_c1ac3int(samp%nc1a,nam%nc3,geom%nl0))
   allocate(samp%smask_c1ac3(samp%nc1a,nam%nc3,geom%nl0))
end if

! Read variables
call mpl%ncerr(subr,nf90_get_var(grpid,c1a_to_c0a_id,samp%c1a_to_c0a))
call mpl%ncerr(subr,nf90_get_var(grpid,smask_c1u_id,smask_c1uint))
do il0=1,geom%nl0
   do ic1u=1,samp%nc1u
      if (smask_c1uint(ic1u,il0)==0) then
         samp%smask_c1u(ic1u,il0) = .false.
      else if (smask_c1uint(ic1u,il0)==1) then
         samp%smask_c1u(ic1u,il0) = .true.
      else
         print*, samp%nc1u,smask_c1uint(ic1u,il0)
         call mpl%abort(subr,'wrong smask_c1u')
      end if
   end do
end do
if (samp%sc2) then
   call mpl%ncerr(subr,nf90_get_var(grpid,c2u_to_c1u_id,samp%c2u_to_c1u))
   call mpl%ncerr(subr,nf90_get_var(grpid,c2u_to_c0u_id,samp%c2u_to_c0u))
end if
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_get_var(grpid,c1ac3_to_c0u_id,samp%c1ac3_to_c0u))
   call mpl%ncerr(subr,nf90_get_var(grpid,smask_c1ac3_id,smask_c1ac3int))
   do il0=1,geom%nl0
      do jc3=1,nam%nc3
         do ic1a=1,samp%nc1a
            if (smask_c1ac3int(ic1a,jc3,il0)==0) then
               samp%smask_c1ac3(ic1a,jc3,il0) = .false.
            else if (smask_c1ac3int(ic1a,jc3,il0)==1) then
               samp%smask_c1ac3(ic1a,jc3,il0) = .true.
            else
               call mpl%abort(subr,'wrong smask_c1ac3')
            end if
         end do
      end do
   end do
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! Release memory
deallocate(smask_c1uint)
if (samp%sc3) deallocate(smask_c1ac3int)

! Other arrays

! Allocation
allocate(samp%smask_c1a(samp%nc1a,geom%nl0))

! Initialization
do il0=1,geom%nl0
   do ic1a=1,samp%nc1a
      samp%smask_c1a(ic1a,il0) = samp%smask_c1ac3(ic1a,1,il0)
   end do
end do

end subroutine samp_read

!----------------------------------------------------------------------
! Subroutine: samp_write
! Purpose: write
!----------------------------------------------------------------------
subroutine samp_write(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(in) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl ! MPI data
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry

! Local variables
integer :: il0,ic1a,ic1u,jc3
integer :: ncid,grpid,nl0_id,nc1a_id,nc1u_id,nc2u_id,nc3_id
integer :: c1a_to_c0a_id,smask_c1u_id,c1ac3_to_c0u_id,smask_c1ac3_id
integer :: c2u_to_c1u_id,c2u_to_c0u_id
integer,allocatable :: smask_c1uint(:,:),smask_c1ac3int(:,:,:)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'samp_write'

write(mpl%info,'(a7,a)') '','Write sampling'
call mpl%flush

! Allocation
allocate(smask_c1uint(samp%nc1u,geom%nl0))
if (samp%sc3) allocate(smask_c1ac3int(samp%nc1a,nam%nc3,geom%nl0))

! Define file
write(filename,'(a,a,i6.6,a,i6.6)') trim(nam%prefix),'_sampling_',mpl%nproc,'-',mpl%myproc
ncid = mpl%nc_file_create_or_open(subr,trim(nam%datadir)//'/'//trim(filename)//'.nc')

! Write grid hash
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'grid_hash',geom%grid_hash))

! Write namelist parameters
call nam%write(mpl,ncid)

! Define group
grpid = mpl%nc_group_define_or_get(subr,ncid,samp%name)

! Define dimensions
nl0_id = mpl%nc_dim_define_or_get(subr,ncid,'nl0',geom%nl0)
nc1a_id = mpl%nc_dim_define_or_get(subr,grpid,'nc1a',samp%nc1a)
nc1u_id = mpl%nc_dim_define_or_get(subr,grpid,'nc1u',samp%nc1u)
if (samp%sc2) nc2u_id = mpl%nc_dim_define_or_get(subr,grpid,'nc2u',samp%nc2u)
if (samp%sc3) nc3_id = mpl%nc_dim_define_or_get(subr,ncid,'nc3',nam%nc3)

! Define variables
c1a_to_c0a_id = mpl%nc_var_define_or_get(subr,grpid,'c1a_to_c0a',nf90_int,(/nc1a_id/))
smask_c1u_id = mpl%nc_var_define_or_get(subr,grpid,'smask_c1u',nf90_int,(/nc1u_id,nl0_id/))
if (samp%sc2) then
   c2u_to_c1u_id = mpl%nc_var_define_or_get(subr,grpid,'c2u_to_c1u',nf90_int,(/nc2u_id/))
   c2u_to_c0u_id = mpl%nc_var_define_or_get(subr,grpid,'c2u_to_c10',nf90_int,(/nc2u_id/))
end if
if (samp%sc3) then
   c1ac3_to_c0u_id = mpl%nc_var_define_or_get(subr,grpid,'c1ac3_to_c0u',nf90_int,(/nc1a_id,nc3_id/))
   smask_c1ac3_id = mpl%nc_var_define_or_get(subr,grpid,'smask_c1ac3',nf90_int,(/nc1a_id,nc3_id,nl0_id/))
end if

! Convert data
do il0=1,geom%nl0
   do ic1u=1,samp%nc1u
      if (samp%smask_c1u(ic1u,il0)) then
         smask_c1uint(ic1u,il0) = 1
      else
         smask_c1uint(ic1u,il0) = 0
      end if
   end do
   if (samp%sc3) then
      do jc3=1,nam%nc3
         do ic1a=1,samp%nc1a
            if (samp%smask_c1ac3(ic1a,jc3,il0)) then
               smask_c1ac3int(ic1a,jc3,il0) = 1
            else
               smask_c1ac3int(ic1a,jc3,il0) = 0
            end if
         end do
      end do
   end if
end do

! Write variables
call mpl%ncerr(subr,nf90_put_var(grpid,c1a_to_c0a_id,samp%c1a_to_c0a))
call mpl%ncerr(subr,nf90_put_var(grpid,smask_c1u_id,smask_c1uint))
if (samp%sc2) then
   call mpl%ncerr(subr,nf90_put_var(grpid,c2u_to_c1u_id,samp%c2u_to_c1u))
   call mpl%ncerr(subr,nf90_put_var(grpid,c2u_to_c0u_id,samp%c2u_to_c0u))
end if
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_put_var(grpid,c1ac3_to_c0u_id,samp%c1ac3_to_c0u))
   call mpl%ncerr(subr,nf90_put_var(grpid,smask_c1ac3_id,smask_c1ac3int))
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! Release memory
deallocate(smask_c1uint)
if (samp%sc3) deallocate(smask_c1ac3int)

end subroutine samp_write

!----------------------------------------------------------------------
! Subroutine: samp_write_grids
! Purpose: write
!----------------------------------------------------------------------
subroutine samp_write_grids(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(in) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl ! MPI data
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry

! Local variables
integer :: il0,jc3,ic1a,ic2u,ic2a,ic2b,jc1u,jc1d,jc1e,ic0a,ic0u,jc0u,j,nc1max,nc1max_tot
integer :: ncid,grpid,nc0a_id,nl0_id,nc1a_id,nc3_id,nc2a_id,nc2b_id,nc1max_id
integer :: lon_c0a_id,lat_c0a_id,gmask_c0a_id
integer :: lon_id,lat_id,lon_local_id,lat_local_id,lon_vbal_id,lat_vbal_id
integer :: igmask_c0a(geom%nc0a,geom%nl0)
real(kind_real),allocatable :: lon(:,:,:),lat(:,:,:)
real(kind_real),allocatable :: lon_local(:,:,:),lat_local(:,:,:)
real(kind_real),allocatable :: lon_vbal(:,:,:),lat_vbal(:,:,:)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'samp_write_grids'

write(mpl%info,'(a7,a)') '','Write sampling grids'
call mpl%flush

! Define file
write(filename,'(a,a,i6.6,a,i6.6)') trim(nam%prefix),'_sampling_grids_',mpl%nproc,'-',mpl%myproc
ncid = mpl%nc_file_create_or_open(subr,trim(nam%datadir)//'/'//trim(filename)//'.nc')

! Write namelist parameters
call nam%write(mpl,ncid)

! Define group
grpid = mpl%nc_group_define_or_get(subr,ncid,samp%name)

! Define dimensions
nc0a_id = mpl%nc_dim_define_or_get(subr,ncid,'nc0a',geom%nc0a)
nl0_id = mpl%nc_dim_define_or_get(subr,ncid,'nl0',geom%nl0)
if (samp%sc3) then
   nc3_id = mpl%nc_dim_define_or_get(subr,ncid,'nc3',nam%nc3)
   nc1a_id = mpl%nc_dim_define_or_get(subr,grpid,'nc1a',samp%nc1a)
end if
if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   nc1max = 0
   do ic2a=1,samp%nc2a
      nc1max = max(count(samp%local_mask(:,ic2a)),nc1max)
   end do
   call mpl%f_comm%allreduce(nc1max,nc1max_tot,fckit_mpi_sum())
   nc1max_tot = nc1max_tot+1
   nc1max_id = mpl%nc_dim_define_or_get(subr,grpid,'nc1max',nc1max_tot)
   nc2a_id = mpl%nc_dim_define_or_get(subr,grpid,'nc2a',samp%nc2a)
end if
if (trim(samp%name)=='vbal') then
   nc1max = 0
   do ic2b=1,samp%nc2b
      nc1max = max(count(samp%vbal_mask(:,ic2b)),nc1max)
   end do
   call mpl%f_comm%allreduce(nc1max,nc1max_tot,fckit_mpi_sum())
   nc1max_tot = nc1max_tot+1
   nc1max_id = mpl%nc_dim_define_or_get(subr,grpid,'nc1max',nc1max_tot)
   nc2b_id = mpl%nc_dim_define_or_get(subr,grpid,'nc2b',samp%nc2b)
end if

! Define variables
lon_c0a_id = mpl%nc_var_define_or_get(subr,ncid,'lon_c0a',nc_kind_real,(/nc0a_id/))
lat_c0a_id = mpl%nc_var_define_or_get(subr,ncid,'lat_c0a',nc_kind_real,(/nc0a_id/))
gmask_c0a_id = mpl%nc_var_define_or_get(subr,ncid,'gmask_c0a',nf90_int,(/nc0a_id,nl0_id/))
if (samp%sc3) then
   lon_id = mpl%nc_var_define_or_get(subr,grpid,'lon',nc_kind_real,(/nc1a_id,nc3_id,nl0_id/))
   lat_id = mpl%nc_var_define_or_get(subr,grpid,'lat',nc_kind_real,(/nc1a_id,nc3_id,nl0_id/))
end if
if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   lon_local_id = mpl%nc_var_define_or_get(subr,grpid,'lon_local',nc_kind_real,(/nc1max_id,nc2a_id,nl0_id/))
   lat_local_id = mpl%nc_var_define_or_get(subr,grpid,'lat_local',nc_kind_real,(/nc1max_id,nc2a_id,nl0_id/))
end if
if (trim(samp%name)=='vbal') then
   lon_vbal_id = mpl%nc_var_define_or_get(subr,grpid,'lon_vbal',nc_kind_real,(/nc1max_id,nc2b_id,nl0_id/))
   lat_vbal_id = mpl%nc_var_define_or_get(subr,grpid,'lat_vbal',nc_kind_real,(/nc1max_id,nc2b_id,nl0_id/))
end if

! Convert data
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      if (geom%gmask_c0a(ic0a,il0)) then
         igmask_c0a(ic0a,il0) = 1
      else
         igmask_c0a(ic0a,il0) = 0
      end if
   end do
end do
if (samp%sc3) then
   ! Allocation
   allocate(lon(samp%nc1a,nam%nc3,geom%nl0))
   allocate(lat(samp%nc1a,nam%nc3,geom%nl0))

   ! Distant points
   lon = mpl%msv%valr
   lat = mpl%msv%valr
   do il0=1,geom%nl0
      do jc3=1,nam%nc3
         do ic1a=1,samp%nc1a
            if (samp%smask_c1ac3(ic1a,jc3,il0)) then
               ic0u = samp%c1ac3_to_c0u(ic1a,jc3)
               lon(ic1a,jc3,il0) = geom%lon_c0u(ic0u)*rad2deg
               lat(ic1a,jc3,il0) = geom%lat_c0u(ic0u)*rad2deg
            end if
         end do
      end do
   end do
end if
if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   ! Allocation
   allocate(lon_local(nc1max_tot,samp%nc2a,geom%nl0))
   allocate(lat_local(nc1max_tot,samp%nc2a,geom%nl0))

   ! Initialization
   lon_local = mpl%msv%valr
   lat_local = mpl%msv%valr

   ! Fill valid points
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         if (samp%smask_c2a(ic2a,il0)) then
            j = 1
            lon_local(j,ic2a,il0) = samp%lon_c2a(ic2a)*rad2deg
            lat_local(j,ic2a,il0) = samp%lat_c2a(ic2a)*rad2deg
            do jc1d=1,samp%nc1d
               jc1u = samp%c1d_to_c1u(jc1d)
               if (samp%local_mask(jc1u,ic2a).and.samp%smask_c1u(jc1u,il0)) then
                  j = j+1
                  jc0u = samp%c1u_to_c0u(jc1u)
                  lon_local(j,ic2a,il0) = geom%lon_c0u(jc0u)*rad2deg
                  lat_local(j,ic2a,il0) = geom%lat_c0u(jc0u)*rad2deg
               end if
            end do
         end if
      end do
   end do
end if
if (trim(samp%name)=='vbal') then
   ! Allocation
   allocate(lon_vbal(nc1max_tot,samp%nc2b,geom%nl0))
   allocate(lat_vbal(nc1max_tot,samp%nc2b,geom%nl0))

   ! Initialization
   lon_vbal = mpl%msv%valr
   lat_vbal = mpl%msv%valr

   ! Fill valid points
   do il0=1,geom%nl0
      do ic2b=1,samp%nc2b
         ic2u = samp%c2b_to_c2u(ic2b)
         if (samp%smask_c2u(ic2u,il0)) then
            j = 1
            lon_vbal(j,ic2b,il0) = samp%lon_c2u(ic2u)*rad2deg
            lat_vbal(j,ic2b,il0) = samp%lat_c2u(ic2u)*rad2deg
            do jc1e=1,samp%nc1e
               jc1u = samp%c1e_to_c1u(jc1e)
               if (samp%vbal_mask(jc1u,ic2b).and.samp%smask_c1u(jc1u,il0)) then
                  j = j+1
                  jc0u = samp%c1u_to_c0u(jc1u)
                  lon_vbal(j,ic2b,il0) = geom%lon_c0u(jc0u)*rad2deg
                  lat_vbal(j,ic2b,il0) = geom%lat_c0u(jc0u)*rad2deg
               end if
            end do
         end if
      end do
   end do
end if

! Write variables
call mpl%ncerr(subr,nf90_put_var(ncid,lon_c0a_id,geom%lon_c0a*rad2deg))
call mpl%ncerr(subr,nf90_put_var(ncid,lat_c0a_id,geom%lat_c0a*rad2deg))
call mpl%ncerr(subr,nf90_put_var(ncid,gmask_c0a_id,igmask_c0a))
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_put_var(grpid,lon_id,lon))
   call mpl%ncerr(subr,nf90_put_var(grpid,lat_id,lat))
end if
if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   call mpl%ncerr(subr,nf90_put_var(grpid,lon_local_id,lon_local))
   call mpl%ncerr(subr,nf90_put_var(grpid,lat_local_id,lat_local))
end if
if (trim(samp%name)=='vbal') then
   call mpl%ncerr(subr,nf90_put_var(grpid,lon_vbal_id,lon_vbal))
   call mpl%ncerr(subr,nf90_put_var(grpid,lat_vbal_id,lat_vbal))
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! Release memory
if (samp%sc3) then
   deallocate(lon)
   deallocate(lat)
end if
if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   deallocate(lon_local)
   deallocate(lat_local)
end if
if (trim(samp%name)=='vbal') then
   deallocate(lon_vbal)
   deallocate(lat_vbal)
end if

end subroutine samp_write_grids

!----------------------------------------------------------------------
! Subroutine: samp_setup_1
! Purpose: setup sampling, first step
!----------------------------------------------------------------------
subroutine samp_setup_1(samp,sname,mpl,rng,nam,geom,ens)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
character(len=*),intent(in) :: sname   ! Sampling name
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(ens_type),intent(in) :: ens       ! Ensemble

! Local variables
integer :: il0,jc3,ildwv,jldwv,ival,nc1_valid
real(kind_real),allocatable :: ldwv_to_lon(:),ldwv_to_lat(:)
logical :: valid
character(len=8) :: ivalformat
character(len=1024) :: color
character(len=1024),parameter :: subr = 'samp_compute_c1'

! Set sampling name
samp%name = sname

! Allocation
call samp%alloc(geom)

! Compute sampling mask
call samp%compute_mask(mpl,nam,geom,ens)

! Compute nearest neighbors for local diagnostics output
if (nam%nldwv>0) then
   write(mpl%info,'(a7,a)') '','Compute local diagnostics locations:'
   call mpl%flush

   ! Allocation
   allocate(samp%ldwv_to_proc(nam%nldwv))
   allocate(samp%ldwv_to_c0a(nam%nldwv))
   allocate(ldwv_to_lon(nam%nldwv))
   allocate(ldwv_to_lat(nam%nldwv))

   ! Initialization
   ldwv_to_lon = 0.0
   ldwv_to_lat = 0.0

   do ildwv=1,nam%nldwv
      ! Get index from lon/lat
      call geom%index_from_lonlat(mpl,nam%lon_ldwv(ildwv),nam%lat_ldwv(ildwv),0,samp%ldwv_to_proc(ildwv),samp%ldwv_to_c0a(ildwv), &
 & valid)
      if (valid) then
         if (samp%ldwv_to_proc(ildwv)==mpl%myproc) then
            ldwv_to_lon(ildwv) = geom%lon_c0a(samp%ldwv_to_c0a(ildwv))
            ldwv_to_lat(ildwv) = geom%lat_c0a(samp%ldwv_to_c0a(ildwv))
         end if
         call mpl%f_comm%broadcast(ldwv_to_lon(ildwv),samp%ldwv_to_proc(ildwv)-1)
         call mpl%f_comm%broadcast(ldwv_to_lat(ildwv),samp%ldwv_to_proc(ildwv)-1)
      else
         call mpl%abort(subr,'profile '//trim(nam%name_ldwv(ildwv))//' is masked or out of the domain')
      end if
   end do

   do ildwv=1,nam%nldwv
      ! Check redundancy
      do jldwv=1,ildwv-1
         if (eq(ldwv_to_lon(ildwv),ldwv_to_lon(jldwv)).and.eq(ldwv_to_lat(ildwv),ldwv_to_lat(jldwv))) &
 & call mpl%abort(subr,'profiles'//trim(nam%name_ldwv(ildwv))//' and '//trim(nam%name_ldwv(jldwv))// &
 & 'are located grid point')
      end do

      ! Print results
      write(mpl%info,'(a10,a,f6.1,a,f6.1)') '','Profile '//trim(nam%name_ldwv(ildwv))//' computed at lon/lat: ', &
 & ldwv_to_lon(ildwv)*rad2deg,' / ',ldwv_to_lat(ildwv)*rad2deg
      call mpl%flush
   end do

   ! Release memory
   deallocate(ldwv_to_lon)
   deallocate(ldwv_to_lat)
end if

! Check subsampling size
if (nam%nc1>samp%nc0_smask(0)) then
   ! Not enough points remaining in the sampling mask
   call mpl%warning(subr,'not enough points remaining in sampling mask, resetting nc1 to the largest possible value')
   nam%nc1 = samp%nc0_smask(0)
end if
if (nam%nc2>nam%nc1) then
   ! Subsampling should have less points
   call mpl%warning(subr,'subsampling should have less points, resetting nc2 to nc1')
   nam%nc2 = nam%nc1
end if

! Allocation
call samp%alloc(nam,geom)

if (nam%sam_read) then
   ! Read sampling
   write(mpl%info,'(a7,a)') '','Read sampling'
   call mpl%flush
   call samp%read(mpl,nam,geom)
else
   ! Compute sampling, subset Sc1
   write(mpl%info,'(a7,a,i5,a)') '','Compute sampling, subset Sc1 (nc1 = ',nam%nc1,')'
   call mpl%flush
   call samp%compute_c1(mpl,rng,nam,geom)

   ! Compute MPI distribution, halo A, subset Sc1
   write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo A, subset Sc1'
   call mpl%flush
   call samp%compute_mpi_c1a(mpl,nam,geom)

   if (samp%sc3) then
      ! Compute sampling, subset Sc3
      write(mpl%info,'(a7,a,i5,a)') '','Compute sampling, subset Sc3 (nc3 = ',nam%nc3,')'
      call mpl%flush
      call samp%compute_c3(mpl,rng,nam,geom)

      ! Check sampling mask
      call samp%check_mask(mpl,nam,geom)
   end if

   if (samp%sc2) then
      ! Compute sampling, subset Sc2
      write(mpl%info,'(a7,a,i5,a)') '','Compute sampling, subset Sc2 (nc2 = ',nam%nc2,')'
      call mpl%flush
      call samp%compute_c2(mpl,rng,nam,geom)
   end if
end if

if (samp%sc2) then
   ! Compute MPI distribution, halo A, subset Sc2
   write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo A, subset Sc2'
   call mpl%flush
   call samp%compute_mpi_c2a(mpl,nam,geom)

   ! Compute MPI distribution, halo B, subset Sc2
   write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo B, subset Sc2'
   call mpl%flush
   call samp%compute_mpi_c2b(mpl,rng,nam,geom)

   if (nam%adv_diag) then
      ! Compute mesh sampling, subset Sc2
      write(mpl%info,'(a7,a)') '','Compute sampling mesh, subset Sc2'
      call mpl%flush
      call samp%compute_mesh_c2(mpl,rng)
   end if
else
   ! Compute MPI distribution, halo A, subset Sc2
   samp%nc2a = 0
end if

if (samp%sc3) then
   ! Print results
   write(mpl%info,'(a7,a)') '','Sampling efficiency (%):'
   call mpl%flush
   do il0=1,geom%nl0
      write(mpl%info,'(a10,a,i3,a3)') '','Level ',nam%levs(il0),' ~>'
      call mpl%flush(.false.)
      do jc3=1,nam%nc3
         call mpl%f_comm%allreduce(count(samp%smask_c1ac3(:,jc3,il0)),nc1_valid,fckit_mpi_sum())
         ival = int(100.0*real(nc1_valid,kind_real)/real(nam%nc1,kind_real))
         if (ival==100) then
            ivalformat = '(a,i3,a)'
         else
            ivalformat = '(a,i2,a)'
         end if
         if (nc1_valid>=nam%nc1/2) then
            ! Sucessful sampling
            color = mpl%green
         else
            ! Insufficient sampling
            color = mpl%peach
         end if
         if (jc3==1) color = ' '//trim(color)
         write(mpl%info,ivalformat) trim(color),ival,trim(mpl%black)
         call mpl%flush(.false.)
         if (jc3<nam%nc3) then
            write(mpl%info,'(a)') '-'
            call mpl%flush(.false.)
         end if
      end do
      write(mpl%info,'(a)') ''
      call mpl%flush
   end do
end if

end subroutine samp_setup_1

!----------------------------------------------------------------------
! Subroutine: samp_setup_2
! Purpose: setup sampling, second step
!----------------------------------------------------------------------
subroutine samp_setup_2(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

if ((trim(samp%name)=='hdiag').or.(trim(samp%name)=='lct')) then
   ! Compute MPI distribution, halo C
   write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo C'
   call mpl%flush
   call samp%compute_mpi_c(mpl,rng,nam,geom)
end if

if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   ! Compute MPI distribution, halos D
   write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo D'
   call mpl%flush
   call samp%compute_mpi_d(mpl,nam,geom)
end if

if (trim(samp%name)=='vbal') then
   ! Compute MPI distribution, halo E
   write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo E'
   call mpl%flush
   call samp%compute_mpi_e(mpl,nam)
end if

! Write sampling data
if (nam%sam_write) then
   write(mpl%info,'(a7,a)') '','Write sampling data'
   call mpl%flush
   call samp%write(mpl,nam,geom)
   if (nam%sam_write_grids) call samp%write_grids(mpl,nam,geom)
end if

! Release memory (partial)
call samp%partial_dealloc

end subroutine samp_setup_2

!----------------------------------------------------------------------
! Subroutine: samp_compute_mask
! Purpose: compute mask
!----------------------------------------------------------------------
subroutine samp_compute_mask(samp,mpl,nam,geom,ens)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry
type(ens_type),intent(in) :: ens       ! Ensemble

! Local variables
integer :: nsmask,nsmask_tot,ic0a,il0,ildwv,iv,its,ie,ncontig,ncontigmax,latmin,latmax
integer :: nc0_smask(0:geom%nl0)
real(kind_real) :: dist
real(kind_real),allocatable :: m2(:,:,:,:)
logical :: valid
character(len=1024),parameter :: subr = 'samp_compute_mask'

! Count extra masked points in sampling
nsmask = count(geom%gmask_c0a.and..not.geom%smask_c0a)
call mpl%f_comm%allreduce(nsmask,nsmask_tot,fckit_mpi_sum())

if ((nsmask_tot>0).or.(trim(nam%mask_type)/='none').or.(nam%ncontig_th>0)) then
   ! Compute sampling mask
   write(mpl%info,'(a7,a)') '','Compute sampling mask'
   call mpl%flush

   ! Copy geometry mask
   samp%smask_c0a = geom%gmask_c0a
   if (allocated(geom%smask_c0a)) samp%smask_c0a = samp%smask_c0a.and.geom%smask_c0a

   ! Mask restriction
   if (nam%mask_type(1:3)=='lat') then
      ! Latitude band
      read(nam%mask_type(4:6),'(i3)') latmin
      read(nam%mask_type(7:9),'(i3)') latmax
      write(mpl%info,'(a10,a,i3,a,i3)') '','Latitude band between ',latmin,' and ',latmax
      call mpl%flush
      if (latmin>=latmax) call mpl%abort(subr,'latmin should be lower than latmax')
      do ic0a=1,geom%nc0a
         valid = (geom%lat_c0a(ic0a)>=real(latmin,kind_real)*deg2rad).and.(geom%lat_c0a(ic0a)<=real(latmax,kind_real)*deg2rad)
         do il0=1,geom%nl0
            samp%smask_c0a(ic0a,il0) = samp%smask_c0a(ic0a,il0).and.valid
         end do
      end do
   elseif (trim(nam%mask_type)=='ldwv') then
      ! Disk around vertical diagnostic points
      write(mpl%info,'(a10,a,e10.3,a)') '','Disk of ',1.1*nam%local_rad*reqkm,' km around vertical diagonstic points'
      call mpl%flush
      do ic0a=1,geom%nc0a
         valid = .false.
         do ildwv=1,nam%nldwv
            call sphere_dist(nam%lon_ldwv(ildwv),nam%lat_ldwv(ildwv),geom%lon_c0a(ic0a),geom%lat_c0a(ic0a),dist)
            valid = valid.or.(dist<1.1*nam%local_rad)
         end do
         do il0=1,geom%nl0
            samp%smask_c0a(ic0a,il0) = samp%smask_c0a(ic0a,il0).and.valid
         end do
      end do
   elseif (trim(nam%mask_type)=='stddev') then
      ! Standard-deviation threshold

      ! Allocation
      allocate(m2(geom%nc0a,geom%nl0,nam%nv,nam%nts))

      ! Compute variance and fourth-order moment
      m2 = 0.0
      do ie=1,ens%ne
         m2 = m2+ens%mem(ie)%fld**2
      end do
      m2 = m2/real(ens%ne-ens%nsub,kind_real)

      ! Check standard-deviation value
      do iv=1,nam%nv
         write(mpl%info,'(a10,a,e10.3,a)') '','Threshold ',nam%mask_th(iv),' used as a '//trim(nam%mask_lu(iv)) &
 & //' bound for standard-deviation'
         call mpl%flush
         do its=1,nam%nts
            if (trim(nam%mask_lu(iv))=='lower') then
               samp%smask_c0a = samp%smask_c0a.and.(m2(:,:,iv,its)>nam%mask_th(iv)**2)
            elseif (trim(nam%mask_lu(iv))=='upper') then
               samp%smask_c0a = samp%smask_c0a.and.(m2(:,:,iv,its)<nam%mask_th(iv)**2)
            end if
         end do
      end do

      ! Release memory
      deallocate(m2)
   else
      if (.not.allocated(geom%smask_c0a)) call mpl%abort(subr,'mask_type not recognized')
   end if

   ! Check vertically contiguous points
   if (nam%ncontig_th>0) then
      write(mpl%info,'(a10,a,i3,a)') '','Mask restricted with at least ',min(nam%ncontig_th,geom%nl0), &
 &  ' vertically contiguous points'
      call mpl%flush
      do ic0a=1,geom%nc0a
         ncontig = 0
         ncontigmax = 0
         do il0=1,geom%nl0
            if (samp%smask_c0a(ic0a,il0)) then
               ncontig = ncontig+1
            else
               ncontig = 0
            end if
            if (ncontig>ncontigmax) ncontigmax = ncontig
         end do
         samp%smask_c0a(ic0a,:) = samp%smask_c0a(ic0a,:).and.(ncontigmax>=min(nam%ncontig_th,geom%nl0))
      end do
   end if
else
   ! Copy geometry mask
   write(mpl%info,'(a7,a)') '','Copy geometry mask'
   call mpl%flush
   samp%smask_c0a = geom%gmask_c0a
end if

! Commnication
call geom%com_AU%ext(mpl,geom%nl0,samp%smask_c0a,samp%smask_c0u)

! Related masks
samp%smask_hor_c0a = any(samp%smask_c0a,dim=2)
samp%smask_hor_c0u = any(samp%smask_c0u,dim=2)
nc0_smask(0) = count(samp%smask_hor_c0a)
nc0_smask(1:geom%nl0) = count(samp%smask_c0a,dim=1)
call mpl%f_comm%allreduce(nc0_smask,samp%nc0_smask,fckit_mpi_sum())

! Check mask size
if (samp%nc0_smask(0)==0) call mpl%abort(subr,'no more points in the sampling mask')

! Print results
write(mpl%info,'(a7,a)') '','Sampling valid points (% of domain mask):'
call mpl%flush
do il0=1,geom%nl0
   write(mpl%info,'(a10,a,i3,a,f5.1,a)') '','Level ',nam%levs(il0),' ~> ',100.0*real(samp%nc0_smask(il0),kind_real) &
 & /real(geom%nc0_gmask(il0),kind_real),'%'
   call mpl%flush
end do

end subroutine samp_compute_mask

!----------------------------------------------------------------------
! Subroutine: samp_compute_c1
! Purpose: compute sampling, subset Sc1
!----------------------------------------------------------------------
subroutine samp_compute_c1(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic0,ic0a,ic0u,ic1,il0i,ildwv,ic1u,iproc
real(kind_real) :: rh_c0a(geom%nc0a)
logical :: smask_hor_c0a(geom%nc0a)
character(len=1024),parameter :: subr = 'samp_compute_c1'

! Select draw type
select case (trim(nam%draw_type))
case ('random_uniform')
   ! Random draw
    rh_c0a = 1.0
case ('random_coast')
   ! More points around coasts
   if (all(geom%gmask_c0a)) call mpl%abort(subr,'random_coast is not relevant if there is no coast')
   rh_c0a = 0.0
   do il0i=1,geom%nl0i
      do ic0a=1,geom%nc0a
         if (geom%gmask_c0a(ic0a,il0i)) then
            rh_c0a(ic0a) = rh_c0a(ic0a)+exp(-geom%mdist_c0a(ic0a,il0i)/nam%Lcoast)
         else
            rh_c0a(ic0a) = rh_c0a(ic0a)+1.0
         end if
      end do
   end do
   rh_c0a = nam%rcoast+(1.0-nam%rcoast)*(1.0-rh_c0a/real(geom%nl0i,kind_real))
end select

! Initialize sampling mask
smask_hor_c0a = samp%smask_hor_c0a

! Initialization
samp%c1_to_c0 = mpl%msv%vali

! Insert local diagnostic points and update sampling mask
do ildwv=1,nam%nldwv
   iproc = samp%ldwv_to_proc(ildwv)
   if (iproc==mpl%myproc) then
      ic0a = samp%ldwv_to_c0a(ildwv)
      ic0 = geom%c0a_to_c0(ic0a)
      samp%c1_to_c0(ildwv) = ic0
      smask_hor_c0a(ic0a) = .false.
   end if
   call mpl%f_comm%broadcast(samp%c1_to_c0(ildwv),iproc-1)
end do

! Initialize sampling
write(mpl%info,'(a7,a)') '','Compute horizontal subset C1: '
call mpl%flush(.false.)
call initialize_sampling(mpl,rng,maxval(geom%area),geom%nc0a,geom%lon_c0a,geom%lat_c0a,smask_hor_c0a,rh_c0a,geom%c0a_to_c0, &
 & nam%ntry,nam%nrep,nam%nc1-nam%nldwv,samp%c1_to_c0(nam%nldwv+1:nam%nc1),n_uni=geom%nc0u,uni_to_loc=geom%c0u_to_c0a, &
 & tree_uni=geom%tree_c0u)

! Count Sc1 point in universe
samp%nc1u = 0
do ic1=1,nam%nc1
   ic0 = samp%c1_to_c0(ic1)
   iproc = geom%c0_to_proc(ic0)
   samp%c1_to_proc(ic1) = iproc
   if (geom%myuniverse(iproc)) samp%nc1u = samp%nc1u+1
end do

! Allocation
allocate(samp%c1u_to_c1(samp%nc1u))
allocate(samp%c1u_to_c0u(samp%nc1u))
allocate(samp%lon_c1u(samp%nc1u))
allocate(samp%lat_c1u(samp%nc1u))
allocate(samp%smask_c1u(samp%nc1u,geom%nl0))
allocate(samp%smask_hor_c1u(samp%nc1u))

! Conversions
samp%c1_to_c1u = mpl%msv%vali
ic1u = 0
do ic1=1,nam%nc1
   ic0 = samp%c1_to_c0(ic1)
   iproc = geom%c0_to_proc(ic0)
   if (geom%myuniverse(iproc)) then
      ic1u = ic1u+1
      ic0u = geom%c0_to_c0u(ic0)
      samp%c1u_to_c1(ic1u) = ic1
      samp%c1_to_c1u(ic1) = ic1u
      samp%c1u_to_c0u(ic1u) = ic0u
      samp%lon_c1u(ic1u) = geom%lon_c0u(ic0u)
      samp%lat_c1u(ic1u) = geom%lat_c0u(ic0u)
      samp%smask_c1u(ic1u,:) = samp%smask_c0u(ic0u,:)
      samp%smask_hor_c1u(ic1u) = samp%smask_hor_c0u(ic0u)
   end if
end do

end subroutine samp_compute_c1

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_c1a
! Purpose: compute MPI distribution, halo A, subset Sc1
!----------------------------------------------------------------------
subroutine samp_compute_mpi_c1a(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic0a,ic0u,ic1a,ic1u,ic1,iproc

! Define universe and halo A on subset Sc1
samp%nc1a = 0
do ic1=1,nam%nc1
   iproc = samp%c1_to_proc(ic1)
   if (iproc==mpl%myproc) samp%nc1a = samp%nc1a+1
end do

! Allocation
allocate(samp%c1a_to_c1(samp%nc1a))
allocate(samp%c1a_to_c1u(samp%nc1a))
allocate(samp%c1u_to_c1a(samp%nc1u))
allocate(samp%c1a_to_c0a(samp%nc1a))
allocate(samp%lon_c1a(samp%nc1a))
allocate(samp%lat_c1a(samp%nc1a))
allocate(samp%smask_c1a(samp%nc1a,geom%nl0))
allocate(samp%smask_hor_c1a(samp%nc1a))
allocate(samp%c1al0_check(samp%nc1a,geom%nl0))

! Conversions
samp%c1u_to_c1a = mpl%msv%vali
ic1a = 0
do ic1=1,nam%nc1
   iproc = samp%c1_to_proc(ic1)
   if (iproc==mpl%myproc) then
      ic1a = ic1a+1
      ic1u = samp%c1_to_c1u(ic1)
      ic0u = samp%c1u_to_c0u(ic1u)
      ic0a = geom%c0u_to_c0a(ic0u)
      samp%c1a_to_c1(ic1a) = ic1
      samp%c1a_to_c1u(ic1a) = ic1u
      samp%c1u_to_c1a(ic1u) = ic1a
      samp%c1a_to_c0a(ic1a) = ic0a
      samp%lon_c1a(ic1a) = geom%lon_c0a(ic0a)
      samp%lat_c1a(ic1a) = geom%lat_c0a(ic0a)
      samp%smask_c1a(ic1a,:) = samp%smask_c0a(ic0a,:)
      samp%smask_hor_c1a(ic1a) = samp%smask_hor_c0a(ic0a)
      samp%c1al0_check(ic1a,:) = (nam%mask_check.and.samp%smask_c0a(ic0a,:))
   end if
end do

! Print results
write(mpl%info,'(a10,a,i8)') '','nc1 = ',nam%nc1
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1a = ',samp%nc1a
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1u = ',samp%nc1u
call mpl%flush

end subroutine samp_compute_mpi_c1a

!----------------------------------------------------------------------
! Subroutine: samp_compute_c3
! Purpose: compute sampling, subset Sc3
!----------------------------------------------------------------------
subroutine samp_compute_c3(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: jc3,ic1a,ir,irtmp,jc0a,jc0u,jc0,icinf,icsup,ictest,nn_index(nam%nc3),il0,iproc
real(kind_real) :: d,nn_dist(nam%nc3)
logical :: found,proc_to_done(mpl%nproc)
character(len=1024),parameter :: subr = 'samp_compute_c3'

! Allocation
allocate(samp%c1ac3_to_c0u(samp%nc1a,nam%nc3))
allocate(samp%smask_c1ac3(samp%nc1a,nam%nc3,geom%nl0))

! Initialization
samp%c1ac3_to_c0u = mpl%msv%vali

if (trim(samp%name)=='hdiag') then
   ! First class
   samp%c1ac3_to_c0u(:,1) = samp%c1u_to_c0u(samp%c1a_to_c1u)

   ! Resynchronize random number generator
   call rng%resync(mpl)

   if (nam%nc3>1) then
      ! Initialization
      write(mpl%info,'(a7,a)') '','Compute HDIAG pairs: '
      call mpl%flush(.false.)
      call mpl%prog_init(nam%nc3*samp%nc1a)
      do ic1a=1,samp%nc1a
         mpl%done((ic1a-1)*nam%nc3+1) = .true.
      end do
      call mpl%f_comm%allgather(all(mpl%done),proc_to_done)
      ir = 0

      ! Sample classes of positive separation
      do while ((.not.all(proc_to_done)).and.(ir<=nam%irmax))
         ! Define a random geographical point
         call geom%rand_point(mpl,rng,0,iproc,jc0a,irtmp)
         ir = ir+irtmp

         if (geom%myuniverse(iproc)) then
            ! Indices
            jc0 = geom%proc_to_c0_offset(iproc)+jc0a
            jc0u = geom%c0_to_c0u(jc0)

            ! Fill classes
            !$omp parallel do schedule(static) private(ic1a,d,jc3,icinf,icsup,found,ictest)
            do ic1a=1,samp%nc1a
               ! Compute the distance
               call sphere_dist(samp%lon_c1a(ic1a),samp%lat_c1a(ic1a),geom%lon_c0u(jc0u),geom%lat_c0u(jc0u),d)

               ! Find the class (dichotomy method)
               if ((d>0.0).and.(d<(real(nam%nc3,kind_real)-0.5)*nam%dc)) then
                  jc3 = 1
                  icinf = 1
                  icsup = nam%nc3
                  found = .false.
                  do while (.not.found)
                     ! New value
                     ictest = (icsup+icinf)/2

                     ! Update
                     if (d<(real(ictest-1,kind_real)-0.5)*nam%dc) icsup = ictest
                     if (d>(real(ictest-1,kind_real)-0.5)*nam%dc) icinf = ictest

                     ! Exit test
                     if (icsup==icinf+1) then
                        if (abs(real(icinf-1,kind_real)*nam%dc-d)<abs(real(icsup-1,kind_real)*nam%dc-d)) then
                           jc3 = icinf
                        else
                           jc3 = icsup
                        end if

                        ! Check class
                        if (d<max((real(jc3-1,kind_real)-0.5)*nam%dc,0.0_kind_real)) call mpl%abort(subr,'jc3 is too high')
                        if (d>(real(jc3,kind_real)-0.5)*nam%dc) call mpl%abort(subr,'jc3 is too low')
                        found = .true.
                     end if
                  end do

                  ! Find if this class has not been aready filled
                  if ((jc3/=1).and.(mpl%msv%is(samp%c1ac3_to_c0u(ic1a,jc3)))) then
                     samp%c1ac3_to_c0u(ic1a,jc3) = jc0u
                     mpl%done((ic1a-1)*nam%nc3+jc3) = .true.   
                  end if
                end if
            end do
            !$omp end parallel do

            ! Update
            call mpl%prog_print
         end if
         call mpl%f_comm%allgather(all(mpl%done),proc_to_done)
      end do
      call mpl%prog_final
   end if

   ! Desynchronize random number generator
   call rng%desync(mpl)
elseif (trim(samp%name)=='lct') then
   ! Initialization
   write(mpl%info,'(a7,a)') '','Compute LCT neighborhood: '
   call mpl%flush(.false.)
   call mpl%prog_init(samp%nc1a)

   do ic1a=1,samp%nc1a
      ! Find neighbors
      call geom%tree_c0u%find_nearest_neighbors(samp%lon_c1a(ic1a),samp%lat_c1a(ic1a),nam%nc3,nn_index,nn_dist)

      ! Copy neighbor index
      do jc3=1,nam%nc3
         samp%c1ac3_to_c0u(ic1a,jc3) = nn_index(jc3)
      end do

      ! Update
      call mpl%prog_print(ic1a)
   end do
   call mpl%prog_final

   ! Define whether this point should be checked for boundaries
   do il0=1,geom%nl0
      do ic1a=1,samp%nc1a
         if (samp%c1al0_check(ic1a,il0)) samp%c1al0_check(ic1a,il0) = any(.not.samp%smask_c0u(samp%c1ac3_to_c0u(ic1a,:),il0))
      end do
   end do
endif

end subroutine samp_compute_c3

!----------------------------------------------------------------------
! Subroutine: samp_check_mask
! Purpose: check sampling mask
!----------------------------------------------------------------------
subroutine samp_check_mask(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: il0,jc3,ic1a,jc0u
logical :: valid

! Check sampling mask
write(mpl%info,'(a7,a)') '','Check sampling mask: '
call mpl%flush(.false.)

! Second point
call mpl%prog_init(samp%nc1a)
do ic1a=1,samp%nc1a
   do il0=1,geom%nl0
      do jc3=1,nam%nc3
         ! Index
         jc0u = samp%c1ac3_to_c0u(ic1a,jc3)

         ! Check point index
         valid = mpl%msv%isnot(jc0u)

         ! Check sampling mask
         if (valid) valid = samp%smask_c1a(ic1a,il0).and.samp%smask_c0u(jc0u,il0)

         ! Check mask boundaries
         if (valid.and.nam%mask_check) call geom%check_arc(mpl,il0,samp%lon_c1a(ic1a),samp%lat_c1a(ic1a),geom%lon_c0u(jc0u), &
 & geom%lat_c0u(jc0u),valid)

         ! Copy validity
         samp%smask_c1ac3(ic1a,jc3,il0) = valid
      end do
   end do

   ! Update
   call mpl%prog_print(ic1a)
end do
call mpl%prog_final

end subroutine samp_check_mask

!----------------------------------------------------------------------
! Subroutine: samp_compute_c2
! Purpose: compute sampling, subset Sc2
!----------------------------------------------------------------------
subroutine samp_compute_c2(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic0u,ic1a,ic1u,ic1,ic2u,ic2,iproc,ildwv
real(kind_real) :: rh_c1a(samp%nc1a)
logical :: smask_hor_c1a(samp%nc1a)

! Initialization
smask_hor_c1a = samp%smask_hor_c1a
rh_c1a = 1.0

! Insert local diagnostic points
do ildwv=1,nam%nldwv
   ic1 = ildwv
   iproc = samp%c1_to_proc(ildwv)
   if (iproc==mpl%myproc) then
      ic1u = samp%c1_to_c1u(ic1)
      ic1a = samp%c1u_to_c1a(ic1u)
      smask_hor_c1a(ic1a) = .false.
   end if
   samp%c2_to_c1(ic1) = ic1
end do

! Initialize subsampling
write(mpl%info,'(a7,a)') '','Compute horizontal subset C2: '
call mpl%flush(.false.)
call initialize_sampling(mpl,rng,maxval(geom%area),samp%nc1a,samp%lon_c1a,samp%lat_c1a,smask_hor_c1a,rh_c1a,samp%c1a_to_c1, &
 & nam%ntry,nam%nrep,nam%nc2-nam%nldwv,samp%c2_to_c1(nam%nldwv+1:nam%nc2))

! Count Sc2 point in universe
samp%nc2u = 0
do ic2=1,nam%nc2
   ic1 = samp%c2_to_c1(ic2)
   iproc = samp%c1_to_proc(ic1)
   if (geom%myuniverse(iproc)) samp%nc2u = samp%nc2u+1
end do

! Allocation
allocate(samp%c2u_to_c2(samp%nc2u))
allocate(samp%c2u_to_c1u(samp%nc2u))
allocate(samp%c2u_to_c0u(samp%nc2u))
allocate(samp%lon_c2u(samp%nc2u))
allocate(samp%lat_c2u(samp%nc2u))
allocate(samp%smask_c2u(samp%nc2u,geom%nl0))
allocate(samp%smask_hor_c2u(samp%nc2u))

! Conversions
samp%c2_to_c2u = mpl%msv%vali
ic2u = 0
do ic2=1,nam%nc2
   ic1 = samp%c2_to_c1(ic2)
   iproc = samp%c1_to_proc(ic1)
   samp%c2_to_proc(ic2) = iproc
   if (geom%myuniverse(iproc)) then
      ic2u = ic2u+1
      ic1u = samp%c1_to_c1u(ic1)
      ic0u = samp%c1u_to_c0u(ic1u)
      samp%c2_to_c2u(ic2) = ic2u
      samp%c2u_to_c2(ic2u) = ic2
      samp%c2u_to_c1u(ic2u) = ic1u
      samp%c2u_to_c0u(ic2u) = ic0u
      samp%lon_c2u(ic2u) = geom%lon_c0u(ic0u)
      samp%lat_c2u(ic2u) = geom%lat_c0u(ic0u)
      samp%smask_c2u(ic2u,:) = samp%smask_c0u(ic0u,:)
      samp%smask_hor_c2u(ic2u) = samp%smask_hor_c0u(ic0u)
   end if
end do

end subroutine samp_compute_c2

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_c2a
! Purpose: compute sampling MPI distribution, halo A, subset Sc2
!----------------------------------------------------------------------
subroutine samp_compute_mpi_c2a(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic0a,ic0u,ic1a,ic1u,ic2a,ic2u,ic2,iproc
type(tree_type) :: tree

! Allocation
allocate(samp%proc_to_nc2a(mpl%nproc))
allocate(samp%proc_to_c2_offset(mpl%nproc))

! Define universe and halo A on subset Sc1
samp%proc_to_nc2a = 0
do ic2=1,nam%nc2
   iproc = samp%c2_to_proc(ic2)
   samp%proc_to_nc2a(iproc) = samp%proc_to_nc2a(iproc)+1
end do
samp%nc2a = samp%proc_to_nc2a(mpl%myproc)

! Offset
samp%proc_to_c2_offset(1) = 0
do iproc=2,mpl%nproc
   samp%proc_to_c2_offset(iproc) = samp%proc_to_c2_offset(iproc-1)+samp%proc_to_nc2a(iproc-1)
end do

! Allocation
allocate(samp%c2a_to_c2(samp%nc2a))
allocate(samp%c2a_to_c2u(samp%nc2a))
allocate(samp%c2a_to_c1a(samp%nc2a))
allocate(samp%c2a_to_c0a(samp%nc2a))
allocate(samp%lon_c2a(samp%nc2a))
allocate(samp%lat_c2a(samp%nc2a))
allocate(samp%smask_c2a(samp%nc2a,geom%nl0))
allocate(samp%smask_hor_c2a(samp%nc2a))

! Conversions
ic2a = 0
do ic2=1,nam%nc2
   iproc = samp%c2_to_proc(ic2)
   if (iproc==mpl%myproc) then
      ic2a = ic2a+1
      ic2u = samp%c2_to_c2u(ic2)
      ic1u = samp%c2u_to_c1u(ic2u)
      ic0u = samp%c2u_to_c0u(ic2u)
      ic1a = samp%c1u_to_c1a(ic1u)
      ic0a = geom%c0u_to_c0a(ic0u)
      samp%c2a_to_c2(ic2a) = ic2
      samp%c2a_to_c2u(ic2a) = ic2u
      samp%c2a_to_c1a(ic2a) = ic1a
      samp%c2a_to_c0a(ic2a) = ic0a
      samp%lon_c2a(ic2a) = geom%lon_c0a(ic0a)
      samp%lat_c2a(ic2a) = geom%lat_c0a(ic0a)
      samp%smask_c2a(ic2a,:) = samp%smask_c0a(ic0a,:)
      samp%smask_hor_c2a(ic2a) = samp%smask_hor_c0a(ic0a)
   end if
end do

! Find nearest neighbors

! Allocation
allocate(samp%nn_c2a_index(samp%nc2u,samp%nc2a))
allocate(samp%nn_c2a_dist(samp%nc2u,samp%nc2a))
call tree%alloc(mpl,samp%nc2u)

! Initialization
call tree%init(samp%lon_c2u,samp%lat_c2u)

! Find nearest neighbors
do ic2a=1,samp%nc2a
   ic0a = samp%c2a_to_c0a(ic2a)
   call tree%find_nearest_neighbors(samp%lon_c2a(ic2a),samp%lat_c2a(ic2a),samp%nc2u,samp%nn_c2a_index(:,ic2a), &
 & samp%nn_c2a_dist(:,ic2a))
end do

! Release memory
call tree%dealloc

! Setup subset Sc2 communication, local to universe
call samp%com_AU%setup(mpl,'com_AU',samp%nc2a,samp%nc2u,nam%nc2,samp%c2a_to_c2,samp%c2u_to_c2)

! Print results
write(mpl%info,'(a10,a,i8)') '','nc2 = ',nam%nc2
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc2a = ',samp%nc2a
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc2u = ',samp%nc2u
call mpl%flush

end subroutine samp_compute_mpi_c2a

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_c2b
! Purpose: compute sampling MPI distribution, halo B
!----------------------------------------------------------------------
subroutine samp_compute_mpi_c2b(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic2b,ic2,ic2u,jc2u,i_s,il0i,iproc
integer,allocatable :: c2b_to_c2(:),c2u_to_c2b(:)
logical :: lcheck_c2b(samp%nc2u)

! Allocation
allocate(samp%h(geom%nl0i))

! Compute interpolation
do il0i=1,geom%nl0i
   write(samp%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
   call samp%h(il0i)%interp(mpl,rng,nam,geom,il0i,samp%nc2u,samp%lon_c2u,samp%lat_c2u,samp%smask_c2u(:,il0i),geom%nc0a, &
 & geom%lon_c0a,geom%lat_c0a,geom%gmask_c0a(:,il0i),7)
end do

! Define halo B
lcheck_c2b = .false.
do ic2u=1,samp%nc2u
   ic2 = samp%c2u_to_c2(ic2u)
   iproc = samp%c2_to_proc(ic2)
   if (iproc==mpl%myproc) lcheck_c2b(ic2u) = .true.
end do
do il0i=1,geom%nl0i
   do i_s=1,samp%h(il0i)%n_s
      jc2u = samp%h(il0i)%col(i_s)
      lcheck_c2b(jc2u) = .true.
   end do
end do
samp%nc2b = count(lcheck_c2b)

! Allocation
allocate(c2b_to_c2(samp%nc2b))
allocate(samp%c2b_to_c2u(samp%nc2b))
allocate(c2u_to_c2b(samp%nc2u))

! Global-local conversion for halo B
c2u_to_c2b = mpl%msv%vali
ic2b = 0
do ic2u=1,samp%nc2u
   if (lcheck_c2b(ic2u)) then
      ic2b = ic2b+1
      ic2 = samp%c2u_to_c2(ic2u)
      c2b_to_c2(ic2b) = ic2
      samp%c2b_to_c2u(ic2b) = ic2
      c2u_to_c2b(ic2u) = ic2b
   end if
end do

do il0i=1,geom%nl0i
   ! Local interpolation source
   samp%h(il0i)%n_src = samp%nc2b
   do i_s=1,samp%h(il0i)%n_s
      samp%h(il0i)%col(i_s) = c2u_to_c2b(samp%h(il0i)%col(i_s))
   end do
end do

! Setup communications
call samp%com_AB%setup(mpl,'com_AB',samp%nc2a,samp%nc2b,nam%nc2,samp%c2a_to_c2,c2b_to_c2)

! Release memory
deallocate(c2b_to_c2)

! Print results
write(mpl%info,'(a10,a,i8)') '','nc2b = ',samp%nc2b
do il0i=1,geom%nl0i
   write(mpl%info,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',samp%h(il0i)%n_s
   call mpl%flush
end do

end subroutine samp_compute_mpi_c2b

!----------------------------------------------------------------------
! Subroutine: samp_compute_mesh_c2
! Purpose: compute sampling mesh, subset Sc2
!----------------------------------------------------------------------
subroutine samp_compute_mesh_c2(samp,mpl,rng)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator

! Alocation
call samp%mesh%alloc(samp%nc2u)

! Initialization
call samp%mesh%init(mpl,rng,samp%lon_c2u,samp%lat_c2u)

! Compute triangles list
write(mpl%info,'(a7,a)') '','Compute triangles list '
call mpl%flush
call samp%mesh%trlist(mpl)

! Find boundary nodes
write(mpl%info,'(a7,a)') '','Find boundary nodes'
call mpl%flush
call samp%mesh%bnodes(mpl)

end subroutine samp_compute_mesh_c2

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_c
! Purpose: compute sampling MPI distribution, halo C
!----------------------------------------------------------------------
subroutine samp_compute_mpi_c(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: jc3,ic0,ic0a,ic0c,ic0u,jc0u,jc0c,ic1a,its,il0,i_s,iproc
integer :: c0u_to_c0c(geom%nc0u)
integer,allocatable :: c0c_to_c0(:)
real(kind_real),allocatable :: lon_c1a(:),lat_c1a(:)
logical :: lcheck_c0c(geom%nc0u)

if ((trim(samp%name)=='hdiag').and.nam%adv_diag) then
   write(mpl%info,'(a7,a)') '','Compute advection interpolation'
   call mpl%flush

   ! Allocation
   allocate(lon_c1a(samp%nc1a))
   allocate(lat_c1a(samp%nc1a))
   allocate(samp%d(geom%nl0,nam%nts))

   ! Compute interpolation
   do its=1,nam%nts
      do il0=1,geom%nl0
         do ic1a=1,samp%nc1a
            ic0a = samp%c1a_to_c0a(ic1a)
            lon_c1a(ic1a) = samp%adv_lon(ic0a,il0,its)
            lat_c1a(ic1a) = samp%adv_lat(ic0a,il0,its)
         end do
         write(samp%d(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
         call samp%d(il0,its)%interp(mpl,rng,nam,geom,il0,geom%nc0u,geom%lon_c0u,geom%lat_c0u,geom%gmask_c0u(:,il0),samp%nc1a, &
 & lon_c1a,lat_c1a,samp%smask_c1a(:,il0),10)
      end do
   end do

   ! Release memory
   deallocate(lon_c1a)
   deallocate(lat_c1a)
end if

! Define halo C
lcheck_c0c = .false.
do ic0u=1,geom%nc0u
   ic0 = geom%c0u_to_c0(ic0u)
   iproc = geom%c0_to_proc(ic0)
   if (iproc==mpl%myproc) lcheck_c0c(ic0u) = .true.
end do
do jc3=1,nam%nc3
   do ic1a=1,samp%nc1a
      if (any(samp%smask_c1ac3(ic1a,jc3,:))) then
         ic0u = samp%c1ac3_to_c0u(ic1a,jc3)
         lcheck_c0c(ic0u) = .true.
      end if
   end do
end do
if ((trim(samp%name)=='hdiag').and.nam%adv_diag) then
   do its=1,nam%nts
      do il0=1,geom%nl0
         do i_s=1,samp%d(il0,its)%n_s
            jc0u = samp%d(il0,its)%col(i_s)
            lcheck_c0c(jc0u) = .true.
         end do
      end do
   end do
end if
samp%nc0c = count(lcheck_c0c)

! Allocation
allocate(c0c_to_c0(samp%nc0c))
allocate(samp%c1a_to_c0c(samp%nc1a))
allocate(samp%c1ac3_to_c0c(samp%nc1a,nam%nc3))

! Initialization
c0u_to_c0c = mpl%msv%vali
samp%c1ac3_to_c0c = mpl%msv%vali

! Conversions
ic0c = 0
do ic0u=1,geom%nc0u
   if (lcheck_c0c(ic0u)) then
      ic0c = ic0c+1
      ic0 = geom%c0u_to_c0(ic0u)
      c0c_to_c0(ic0c) = ic0
      c0u_to_c0c(ic0u) = ic0c
   end if
end do
do ic1a=1,samp%nc1a
   ic0a = samp%c1a_to_c0a(ic1a)
   ic0u = geom%c0a_to_c0u(ic0a)
   ic0c = c0u_to_c0c(ic0u)
   samp%c1a_to_c0c(ic1a) = ic0c
   do jc3=1,nam%nc3
      if (any(samp%smask_c1ac3(ic1a,jc3,:))) then
         jc0u = samp%c1ac3_to_c0u(ic1a,jc3)
         jc0c = c0u_to_c0c(jc0u)
         samp%c1ac3_to_c0c(ic1a,jc3) = jc0c
      end if
   end do
end do

if ((trim(samp%name)=='hdiag').and.nam%adv_diag) then
   ! Local interpolation source
   do its=1,nam%nts
      do il0=1,geom%nl0
         samp%d(il0,its)%n_src = samp%nc0c
         do i_s=1,samp%d(il0,its)%n_s
            samp%d(il0,its)%col(i_s) = c0u_to_c0c(samp%d(il0,its)%col(i_s))
         end do
      end do
   end do
end if

! Setup communications
call samp%com_AC%setup(mpl,'com_AC',geom%nc0a,samp%nc0c,geom%nc0,geom%c0a_to_c0,c0c_to_c0)

! Release memory
deallocate(c0c_to_c0)

! Print results
write(mpl%info,'(a7,a,i6)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0c =      ',samp%nc0c
call mpl%flush

end subroutine samp_compute_mpi_c

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_d
! Purpose: compute sampling MPI distribution, halo D
!----------------------------------------------------------------------
subroutine samp_compute_mpi_d(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic2a,nn,i,ic1d,ic1,ic1a,ic1u,jc1u,npack,ipack,il0,jc3
integer,allocatable :: nn_index(:),c1d_to_c1(:)
logical :: lcheck_c1d(samp%nc1u)
logical,allocatable :: sbuf(:,:),rbuf(:,:)
type(tree_type) :: tree

! Allocation
allocate(samp%local_mask(samp%nc1u,samp%nc2a))
call tree%alloc(mpl,samp%nc1u,mask=samp%smask_hor_c1u)

! Initialization
samp%local_mask = .false.
lcheck_c1d = .false.
do ic1a=1,samp%nc1a
   ic1u = samp%c1a_to_c1u(ic1a)
   lcheck_c1d(ic1u) = .true.
end do
call tree%init(samp%lon_c1u,samp%lat_c1u)

! Define masks
do ic2a=1,samp%nc2a
   ! Count nearest neighbors
   call tree%count_nearest_neighbors(samp%lon_c2a(ic2a),samp%lat_c2a(ic2a),nam%local_rad,nn)
   nn = max(nn,1)

   ! Allocation
   allocate(nn_index(nn))

   ! Find nearest neighbors
   call tree%find_nearest_neighbors(samp%lon_c2a(ic2a),samp%lat_c2a(ic2a),nn,nn_index)

   ! Update masks
   do i=1,nn
      jc1u = nn_index(i)
      samp%local_mask(jc1u,ic2a) = .true.
      lcheck_c1d(jc1u) = .true.
   end do

   ! Release memory
   deallocate(nn_index)
end do
samp%nc1d = count(lcheck_c1d)

! Release memory
call tree%dealloc

! Allocation
allocate(c1d_to_c1(samp%nc1d))
allocate(samp%c1d_to_c1u(samp%nc1d))

! Halo D
ic1d = 0
do ic1u=1,samp%nc1u
   if (lcheck_c1d(ic1u)) then
      ic1d = ic1d+1
      ic1 = samp%c1u_to_c1(ic1u)
      c1d_to_c1(ic1d) = ic1
      samp%c1d_to_c1u(ic1d) = ic1u
   end if
end do

! Setup communications
call samp%com_AD%setup(mpl,'com_AD',samp%nc1a,samp%nc1d,nam%nc1,samp%c1a_to_c1,c1d_to_c1)

! Allocation
npack = nam%nc3*geom%nl0
allocate(sbuf(samp%nc1a,npack))
allocate(rbuf(samp%nc1d,npack))
allocate(samp%smask_c1dc3(samp%nc1d,nam%nc3,geom%nl0))

! Pack
ipack = 0
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      ipack = ipack+1
      sbuf(:,ipack) = samp%smask_c1ac3(:,jc3,il0)
   end do
end do

! Communication
call samp%com_AD%ext(mpl,npack,sbuf,rbuf)

! Unpack
ipack = 0
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      ipack = ipack+1
      samp%smask_c1dc3(:,jc3,il0) = rbuf(:,ipack)
   end do
end do

! Release memory
deallocate(c1d_to_c1)
deallocate(sbuf)
deallocate(rbuf)

! Print results
write(mpl%info,'(a7,a,i6)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1d =      ',samp%nc1d
call mpl%flush

end subroutine samp_compute_mpi_d

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_e
! Purpose: compute sampling MPI distribution, halo E
!----------------------------------------------------------------------
subroutine samp_compute_mpi_e(samp,mpl,nam)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist

! Local variables
integer :: ic2b,ic2u,nn,i,ic1,ic1e,ic1u,ic1a,jc1u
integer,allocatable :: nn_index(:),c1e_to_c1(:)
logical :: lcheck_c1e(samp%nc1u)
character(len=1024),parameter :: subr = 'samp_compute_mpi_e'
type(tree_type) :: tree

! Allocation
allocate(samp%vbal_mask(samp%nc1u,samp%nc2b))
call tree%alloc(mpl,samp%nc1u,mask=samp%smask_hor_c1u)

! Initialization
samp%vbal_mask = .false.
lcheck_c1e = .false.
do ic1a=1,samp%nc1a
   ic1u = samp%c1a_to_c1u(ic1a)
   lcheck_c1e(ic1u) = .true.
end do
call tree%init(samp%lon_c1u,samp%lat_c1u)

! Halo E
do ic2b=1,samp%nc2b
   ! Indices
   ic2u = samp%c2b_to_c2u(ic2b)
   ic1u = samp%c2u_to_c1u(ic2u)

   ! Origin point
   samp%vbal_mask(ic1u,ic2b) = .true.
   lcheck_c1e(ic1u) = .true.

   if (nam%vbal_rad>0.0) then
      ! Count nearest neighbors
      call tree%count_nearest_neighbors(samp%lon_c2u(ic2u),samp%lat_c2u(ic2u),nam%vbal_rad,nn)

      ! Allocation
      allocate(nn_index(nn))

      ! Find nearest neighbors
      call tree%find_nearest_neighbors(samp%lon_c2u(ic2u),samp%lat_c2u(ic2u),nn,nn_index)

      ! Update masks
      do i=1,nn
         jc1u = nn_index(i)
         samp%vbal_mask(jc1u,ic2b) = .true.
         lcheck_c1e(jc1u) = .true.
      end do

      ! Release memory
      deallocate(nn_index)
   elseif (nam%vbal_dlat>0.0) then
      ! Update masks
      do jc1u=1,samp%nc1u
         if (abs(samp%lat_c2u(ic2u)-samp%lat_c1u(jc1u))<nam%vbal_dlat) then
            samp%vbal_mask(jc1u,ic2b) = .true.
            lcheck_c1e(jc1u) = .true.
         end if
      end do
   else
      call mpl%abort(subr,'vbal_rad or vbal_dlat should be positive')
   end if
end do
samp%nc1e = count(lcheck_c1e)

! Release memory
call tree%dealloc

! Halo E
allocate(c1e_to_c1(samp%nc1e))
allocate(samp%c1e_to_c1u(samp%nc1e))
ic1e = 0
do ic1u=1,samp%nc1u
   if (lcheck_c1e(ic1u)) then
      ic1e = ic1e+1
      ic1 = samp%c1u_to_c1(ic1u)
      c1e_to_c1(ic1e) = ic1
      samp%c1e_to_c1u(ic1e) = ic1u
   end if
end do

! Setup communications
call samp%com_AE%setup(mpl,'com_AE',samp%nc1a,samp%nc1e,nam%nc1,samp%c1a_to_c1,c1e_to_c1)

! Release memory
deallocate(c1e_to_c1)

! Print results
write(mpl%info,'(a7,a,i6)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1e =      ',samp%nc1e
call mpl%flush

end subroutine samp_compute_mpi_e

!----------------------------------------------------------------------
! Subroutine: samp_diag_filter
! Purpose: filter diagnostics
!----------------------------------------------------------------------
subroutine samp_diag_filter(samp,mpl,nam,filter_type,rflt,diag_c2a,val_c2a)

implicit none

! Passed variables
class(samp_type),intent(in) :: samp                       ! Sampling
type(mpl_type),intent(inout) :: mpl                       ! MPI data
type(nam_type),intent(in) :: nam                          ! Namelist
character(len=*),intent(in) :: filter_type                ! Filter type
real(kind_real),intent(in) :: rflt                        ! Filter support radius
real(kind_real),intent(inout) :: diag_c2a(samp%nc2a)      ! Filtered diagnostic
real(kind_real),intent(in),optional :: val_c2a(samp%nc2a) ! Useful value for filtering

! Local variables
integer :: ic2,ic2a,nc2f,ic2f,ic2u,jc2u,nc2eff,ic2eff,kc2u,kc2f
integer :: c2u_to_c2f(samp%nc2u)
integer,allocatable :: c2f_to_c2(:),order(:)
real(kind_real) :: distnorm,norm,wgt
real(kind_real),allocatable :: diag_c2f(:),diag_eff(:),diag_eff_dist(:)
real(kind_real),allocatable :: val_c2f(:),val_eff(:)
logical :: lcheck_c2f(samp%nc2u)
character(len=1024),parameter :: subr = 'samp_diag_filter'
type(com_type) :: com_AF

if (rflt>0.0) then
   ! Define halo F
   lcheck_c2f = .false.
   do ic2a=1,samp%nc2a
      ic2u = samp%c2a_to_c2u(ic2a)
      lcheck_c2f(ic2u) = .true.
      jc2u = 1
      do while (inf(samp%nn_c2a_dist(jc2u,ic2a),rflt))
         kc2u = samp%nn_c2a_index(jc2u,ic2a)
         lcheck_c2f(kc2u) = .true.
         jc2u = jc2u+1
         if (jc2u>samp%nc2u) exit
      end do
   end do
   nc2f = count(lcheck_c2f)

   ! Allocation
   allocate(c2f_to_c2(nc2f))

   ! Global-local conversion for halo F
   c2u_to_c2f = mpl%msv%vali
   ic2f = 0
   do ic2u=1,samp%nc2u
      if (lcheck_c2f(ic2u)) then
         ic2f = ic2f+1
         ic2 = samp%c2u_to_c2(ic2u)
         c2f_to_c2(ic2f) = ic2
         c2u_to_c2f(ic2u) = ic2f
      end if
   end do

   ! Setup communications
   call com_AF%setup(mpl,'com_AF',samp%nc2a,nc2f,nam%nc2,samp%c2a_to_c2,c2f_to_c2)

   ! Allocation
   allocate(diag_c2f(nc2f))
   if (present(val_c2a)) allocate(val_c2f(nc2f))

   ! Communication
   call com_AF%ext(mpl,diag_c2a,diag_c2f)
   if (present(val_c2a)) call com_AF%ext(mpl,val_c2a,val_c2f)

   !$omp parallel do schedule(static) private(ic2a,nc2eff,ic2eff,jc2u,kc2u,kc2f,distnorm,norm,wgt), &
   !$omp&                             firstprivate(diag_eff,diag_eff_dist,val_eff,order)
   do ic2a=1,samp%nc2a
      ! Count involved points
      nc2eff = 0
      jc2u = 1
      do while (inf(samp%nn_c2a_dist(jc2u,ic2a),rflt))
         ! Check the point validity
         kc2u = samp%nn_c2a_index(jc2u,ic2a)
         kc2f = c2u_to_c2f(kc2u)
         if (mpl%msv%isnot(diag_c2f(kc2f))) nc2eff = nc2eff+1
         jc2u = jc2u+1
         if (jc2u>samp%nc2u) exit
      end do

      ! Allocation
      allocate(diag_eff(nc2eff))
      allocate(diag_eff_dist(nc2eff))
      if (present(val_c2a)) allocate(val_eff(nc2eff))

      ! Build diag_eff of valid points
      ic2eff = 0
      jc2u = 1
      do while (inf(samp%nn_c2a_dist(jc2u,ic2a),rflt))
         ! Check the point validity
         kc2u = samp%nn_c2a_index(jc2u,ic2a)
         kc2f = c2u_to_c2f(kc2u)
         if (mpl%msv%isnot(diag_c2f(kc2f))) then
            ic2eff = ic2eff+1
            diag_eff(ic2eff) = diag_c2f(kc2f)
            diag_eff_dist(ic2eff) = samp%nn_c2a_dist(jc2u,ic2a)
            if (present(val_c2a)) val_eff(ic2eff) = val_c2f(kc2f)
         end if
         jc2u = jc2u+1
         if (jc2u>samp%nc2u) exit
      end do

      ! Apply filter
      if (nc2eff>0) then
         select case (trim(filter_type))
         case ('average')
            ! Compute average
            diag_c2a(ic2a) = sum(diag_eff)/real(nc2eff,kind_real)
         case ('gc99')
            ! Gaspari-Cohn (1999) kernel
            diag_c2a(ic2a) = 0.0
            norm = 0.0
            do ic2eff=1,nc2eff
               distnorm = diag_eff_dist(ic2eff)/rflt
               wgt = fit_func(mpl,distnorm)
               diag_c2a(ic2a) = diag_c2a(ic2a)+wgt*diag_eff(ic2eff)
               norm = norm+wgt
            end do
            if (norm>0.0) diag_c2a(ic2a) = diag_c2a(ic2a)/norm
         case ('median')
            ! Compute median
            allocate(order(nc2eff))
            if (present(val_c2a)) then
               ! Use external value
               call qsort(nc2eff,val_eff,order)
               diag_eff = diag_eff(order)
            else
               ! Use diagnostic value
               call qsort(nc2eff,diag_eff,order)
            end if
            if (mod(nc2eff,2)==0) then
               diag_c2a(ic2a) = 0.5*(diag_eff(nc2eff/2)+diag_eff(nc2eff/2+1))
            else
               diag_c2a(ic2a) = diag_eff((nc2eff+1)/2)
            end if
            deallocate(order)
         case default
            ! Wrong filter
            call mpl%abort(subr,'wrong filter type')
         end select
      else
         diag_c2a(ic2a) = mpl%msv%valr
      end if

      ! Release memory
      deallocate(diag_eff)
      deallocate(diag_eff_dist)
      if (present(val_c2a)) deallocate(val_eff)
   end do
   !$omp end parallel do

   ! Release memory
   deallocate(c2f_to_c2)
   call com_AF%dealloc
   deallocate(diag_c2f)
   if (present(val_c2a)) deallocate(val_c2f)
end if

end subroutine samp_diag_filter

!----------------------------------------------------------------------
! Subroutine: samp_diag_fill
! Purpose: fill diagnostics missing values
!----------------------------------------------------------------------
subroutine samp_diag_fill(samp,mpl,diag_c2a)

implicit none

! Passed variables
class(samp_type),intent(in) :: samp                  ! Sampling
type(mpl_type),intent(inout) :: mpl                  ! MPI data
real(kind_real),intent(inout) :: diag_c2a(samp%nc2a) ! Filtered diagnostic

! Local variables
integer :: nmsr,nmsr_tot,ic2a,jc2u,kc2u
real(kind_real),allocatable :: diag_c2u(:)

! Count missing points
if (samp%nc2a>0) then
   nmsr = count(mpl%msv%is(diag_c2a))
else
   nmsr = 0
end if
call mpl%f_comm%allreduce(nmsr,nmsr_tot,fckit_mpi_sum())

if (nmsr_tot>0) then
   ! Allocation
   allocate(diag_c2u(samp%nc2u))

   ! Communication
   call samp%com_AU%ext(mpl,diag_c2a,diag_c2u)

   ! Fill points
   do ic2a=1,samp%nc2a
      jc2u = 1
      do while (mpl%msv%is(diag_c2a(ic2a)))
         kc2u = samp%nn_c2a_index(jc2u,ic2a)
         if (mpl%msv%isnot(diag_c2u(kc2u))) diag_c2a(ic2a) = diag_c2u(kc2u)
         jc2u = jc2u+1
         if (jc2u>samp%nc2u) exit
      end do
   end do

   ! Release memory
   deallocate(diag_c2u)
end if

end subroutine samp_diag_fill

end module type_samp
