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
   ! Sampling
   character(len=1024) :: name                      ! Sampling name
   logical :: new_sampling                          ! New sampling flag
   logical,allocatable :: smask_c0u(:,:)            ! Mask on subset Sc0, universe
   logical,allocatable :: smask_c0a(:,:)            ! Mask on subset Sc0, halo A
   logical,allocatable :: smask_hor_c0u(:)          ! Union of horizontal masks on subset Sc0, universe
   logical,allocatable :: smask_hor_c0a(:)          ! Union of horizontal masks on subset Sc0, halo A
   integer,allocatable :: nc0_smask(:)              ! Horizontal mask size on subset Sc0
   integer,allocatable :: c1_to_c0(:)               ! First sampling index
   logical,allocatable :: c1al0_check(:,:)          ! Mask boundaries checking activation
   logical,allocatable :: smask_c1(:,:)             ! Log for the first sampling index
   integer,allocatable :: c1ac3_to_c0(:,:)          ! Second horizontal sampling index, halo A
   logical,allocatable :: smask_c1ac3(:,:,:)        ! Log for the second horizontal sampling index, halo A
   logical,allocatable :: smask_c1dc3(:,:,:)        ! Log for the second horizontal sampling index, halo D
   integer,allocatable :: c2_to_c1(:)               ! Subgrid to diagnostic points
   integer,allocatable :: c2_to_c0(:)               ! Subgrid to grid
   real(kind_real),allocatable :: lon_c2(:)         ! Longitudes on subset Sc2
   real(kind_real),allocatable :: lat_c2(:)         ! Latitudes on subset Sc2
   logical,allocatable :: mask_c2(:,:)              ! Mask on subset Sc2

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

   ! MPI distribution
   integer :: nc0c                                  ! Number of points in subset Sc0, halo C
   integer :: nc1a                                  ! Number of points in subset Sc1, halo A
   integer :: nc1d                                  ! Number of points in subset Sc1, halo D
   integer :: nc1e                                  ! Number of points in subset Sc1, halo E
   logical :: sc2                                   ! Subset Sc2 flag
   logical :: sc3                                   ! Subset Sc3 flag
   integer :: nc2a                                  ! Number of points in subset Sc2, halo A
   integer :: nc2b                                  ! Number of points in subset Sc2, halo B
   integer :: nc2f                                  ! Number of points in subset Sc2, halo F
   logical,allocatable :: lcheck_c0a(:)             ! Detection of halo A on subset Sc0
   logical,allocatable :: lcheck_c1a(:)             ! Detection of halo A on subset Sc1
   logical,allocatable :: lcheck_c2a(:)             ! Detection of halo A on subset Sc2
   integer,allocatable :: c0_to_c0c(:)              ! Subset Sc0, global to halo C
   integer,allocatable :: c1a_to_c1(:)              ! Subset Sc1, halo A to global
   integer,allocatable :: c1_to_c1a(:)              ! Subset Sc1, global to halo A
   integer,allocatable :: c1_to_proc(:)             ! Subset Sc1, global to processor
   integer,allocatable :: c1d_to_c1(:)              ! Subset Sc1, halo D to global
   integer,allocatable :: c1_to_c1d(:)              ! Subset Sc1, global to halo D
   integer,allocatable :: c1e_to_c1(:)              ! Subset Sc1, halo E to global
   integer,allocatable :: c1_to_c1e(:)              ! Subset Sc1, global to halo E
   integer,allocatable :: c2a_to_c2(:)              ! Subset Sc2, halo A to global
   integer,allocatable :: c2_to_c2a(:)              ! Subset Sc2, global to halo A
   logical,allocatable :: mask_c2a(:,:)             ! Mask on subset Sc2, halo A
   integer,allocatable :: c2a_to_c1a(:)             ! Subgrid to diagnostic points, halo A
   integer,allocatable :: c2b_to_c2(:)              ! Subset Sc2, halo B to global
   integer,allocatable :: c2_to_c2f(:)              ! Subset Sc2, global to halo F
   integer,allocatable :: c2_to_proc(:)             ! Subset Sc2, global to processor
   type(com_type) :: com_AB                         ! Communication between halos A and B
   type(com_type) :: com_AC                         ! Communication between halos A and C
   type(com_type) :: com_AD                         ! Communication between halos A and D (diagnostic)
   type(com_type) :: com_AE                         ! Communication between halos A and E (vertical balance)
   type(com_type) :: com_AF                         ! Communication between halos A and F (filtering)
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
   procedure :: compute_zs => samp_compute_zs
   procedure :: compute_mpi_c1a => samp_compute_mpi_c1a
   procedure :: compute_ps => samp_compute_ps
   procedure :: check_mask => samp_check_mask
   procedure :: compute_c2 => samp_compute_c2
   procedure :: compute_mesh_c2 => samp_compute_mesh_c2
   procedure :: compute_mpi_c2a => samp_compute_mpi_c2a
   procedure :: compute_mpi_b => samp_compute_mpi_b
   procedure :: compute_mpi_c => samp_compute_mpi_c
   procedure :: compute_mpi_d => samp_compute_mpi_d
   procedure :: compute_mpi_e => samp_compute_mpi_e
   procedure :: compute_mpi_f => samp_compute_mpi_f
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
allocate(samp%smask_c0(geom%nc0,geom%nl0))
allocate(samp%smask_c0a(geom%nc0a,geom%nl0))
allocate(samp%smask_hor_c0(geom%nc0))
allocate(samp%smask_hor_c0a(geom%nc0a))
allocate(samp%nc0_smask(0:geom%nl0))

end subroutine samp_alloc_mask

!----------------------------------------------------------------------
! Subroutine: samp_alloc_other
! Purpose: allocation for other variables
!------------------------------------------------g---------------------
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
allocate(samp%smask_c1(nam%nc1,geom%nl0))
if (samp%sc2) then
   allocate(samp%c2_to_c1(nam%nc2))
   allocate(samp%c2_to_c0(nam%nc2))
   allocate(samp%lon_c2(nam%nc2))
   allocate(samp%lat_c2(nam%nc2))
   allocate(samp%smask_c2(nam%nc2,geom%nl0))
   call samp%mesh%alloc(nam%nc2)
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
if (allocated(samp%smask_c0)) deallocate(samp%smask_c0)
if (allocated(samp%smask_c0a)) deallocate(samp%smask_c0a)
if (allocated(samp%smask_hor_c0)) deallocate(samp%smask_hor_c0)
if (allocated(samp%smask_hor_c0a)) deallocate(samp%smask_hor_c0a)
if (allocated(samp%nc0_smask)) deallocate(samp%nc0_smask)
if (allocated(samp%c1al0_check)) deallocate(samp%c1al0_check)
if (allocated(samp%adv_nn)) deallocate(samp%adv_nn)
if (allocated(samp%adv_nn_index)) deallocate(samp%adv_nn_index)
if (allocated(samp%lcheck_c1a)) deallocate(samp%lcheck_c1a)
if (allocated(samp%lcheck_c2a)) deallocate(samp%lcheck_c2a)

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
if (allocated(samp%c1_to_c0)) deallocate(samp%c1_to_c0)
if (allocated(samp%smask_c1)) deallocate(samp%smask_c1)
if (allocated(samp%c1ac3_to_c0)) deallocate(samp%c1ac3_to_c0)
if (allocated(samp%smask_c1ac3)) deallocate(samp%smask_c1ac3)
if (allocated(samp%c2_to_c1)) deallocate(samp%c2_to_c1)
if (allocated(samp%c2_to_c0)) deallocate(samp%c2_to_c0)
if (allocated(samp%smask_c2)) deallocate(samp%smask_c2)
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
if (allocated(samp%lcheck_c0a)) deallocate(samp%lcheck_c0a)
if (allocated(samp%c0_to_c0c)) deallocate(samp%c0_to_c0c)
if (allocated(samp%c1a_to_c1)) deallocate(samp%c1a_to_c1)
if (allocated(samp%c1_to_c1a)) deallocate(samp%c1_to_c1a)
if (allocated(samp%c1_to_proc)) deallocate(samp%c1_to_proc)
if (allocated(samp%c1d_to_c1)) deallocate(samp%c1d_to_c1)
if (allocated(samp%c1e_to_c1)) deallocate(samp%c1e_to_c1)
if (allocated(samp%c2a_to_c2)) deallocate(samp%c2a_to_c2)
if (allocated(samp%c2_to_c2a)) deallocate(samp%c2_to_c2a)
if (allocated(samp%smask_c2a)) deallocate(samp%smask_c2a)
if (allocated(samp%c2a_to_c1a)) deallocate(samp%c2a_to_c1a)
if (allocated(samp%c2b_to_c2)) deallocate(samp%c2b_to_c2)
if (allocated(samp%c2_to_c2f)) deallocate(samp%c2_to_c2f)
if (allocated(samp%c2_to_proc)) deallocate(samp%c2_to_proc)
call samp%com_AB%dealloc
call samp%com_AC%dealloc
call samp%com_AD%dealloc
call samp%com_AE%dealloc
call samp%com_AF%dealloc

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
integer :: il0,ic1,ic1a,jc3
integer :: info,ncid,nl0_id,nc1_id,nc1a_id,nc2_id,nc3_id
integer :: c1_to_c0_id,smask_c1_id,c1ac3_to_c0_id,smask_c1ac3_id
integer :: c2_to_c1_id,c2_to_c0_id
integer :: new_sampling,new_sampling_tot
integer,allocatable :: smask_c1int(:,:),smask_c1ac3int(:,:,:)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'samp_read'

! Initialization
new_sampling = 0

! Allocation
allocate(smask_c1int(nam%nc1,geom%nl0))
if (samp%sc3) allocate(smask_c1ac3int(samp%nc1a,nam%nc3,geom%nl0))

! Open file
write(filename,'(a,a,i4.4,a,i4.4)') trim(nam%prefix),'_sampling_',mpl%nproc,'-',mpl%myproc
info = nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid)
if (info/=nf90_noerr) then
   call mpl%warning(subr,'cannot find sampling to read, recomputing sampling')
   new_sampling = 1
end if
call mpl%f_comm%allreduce(new_sampling,new_sampling_tot,fckit_mpi_sum())
samp%new_sampling = (new_sampling_tot>0)
if (samp%new_sampling) return

! Check dimensions
nl0_id = mpl%ncdimcheck(subr,ncid,'nl0',geom%nl0,.false.)
nc1_id = mpl%ncdimcheck(subr,ncid,'nc1',nam%nc1,.false.)
nc1a_id = mpl%ncdimcheck(subr,ncid,'nc1a',samp%nc1a,.false.)
if (samp%sc2) then
   nc2_id = mpl%ncdimcheck(subr,ncid,'nc2',nam%nc2,.false.)
else
   nc2_id = 0
end if
if (samp%sc3) then
   nc3_id = mpl%ncdimcheck(subr,ncid,'nc3',nam%nc3,.false.)
else
   nc3_id = 0
end if
if (mpl%msv%is(nl0_id).or.mpl%msv%is(nc1a_id).or.mpl%msv%is(nc2_id).or.mpl%msv%is(nc3_id)) then
   call mpl%warning(subr,'wrong dimension when reading sampling, recomputing sampling')
   call mpl%ncerr(subr,nf90_close(ncid))
   new_sampling = 1
end if
call mpl%f_comm%allreduce(new_sampling,new_sampling_tot,fckit_mpi_sum())
samp%new_sampling = (new_sampling_tot>0)
if (samp%new_sampling) return

! Compute halo A
call samp%compute_mpi_c1a(mpl,nam,geom)

! Read sampling
write(mpl%info,'(a7,a)') '','Read sampling'
call mpl%flush

! Get arrays ID
call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(samp%name)//'_c1_to_c0',c1_to_c0_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(samp%name)//'_smask_c1',smask_c1_id))
if (samp%sc2) then
   call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(samp%name)//'_c2_to_c1',c2_to_c1_id))
   call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(samp%name)//'_c2_to_c0',c2_to_c0_id))
end if
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(samp%name)//'_c1ac3_to_c0',c1ac3_to_c0_id))
   call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(samp%name)//'_smask_c1ac3',smask_c1ac3_id))
end if

! Read arrays
call mpl%ncerr(subr,nf90_get_var(ncid,c1_to_c0_id,samp%c1_to_c0))
call mpl%ncerr(subr,nf90_get_var(ncid,smask_c1_id,smask_c1int))
do il0=1,geom%nl0
   do ic1=1,nam%nc1
      if (smask_c1int(ic1,il0)==0) then
         samp%smask_c1(ic1,il0) = .false.
      else if (smask_c1int(ic1,il0)==1) then
         samp%smask_c1(ic1,il0) = .true.
      else
         call mpl%abort(subr,'wrong smask_c1')
      end if
   end do
end do
if (samp%sc2) then
   call mpl%ncerr(subr,nf90_get_var(ncid,c2_to_c1_id,samp%c2_to_c1))
   call mpl%ncerr(subr,nf90_get_var(ncid,c2_to_c0_id,samp%c2_to_c0))
end if
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_get_var(ncid,c1ac3_to_c0_id,samp%c1ac3_to_c0))
   call mpl%ncerr(subr,nf90_get_var(ncid,smask_c1ac3_id,smask_c1ac3int))
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
deallocate(smask_c1int)
if (samp%sc3) deallocate(smask_c1ac3int)

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
integer :: il0,ic1,ic1a,jc3
integer :: ncid,nl0_id,nc1_id,nc1a_id,nc2_id,nc3_id
integer :: c1_to_c0_id,smask_c1_id,c1ac3_to_c0_id,smask_c1ac3_id
integer :: c2_to_c1_id,c2_to_c0_id
integer,allocatable :: smask_c1int(:,:),smask_c1ac3int(:,:,:)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'samp_write'

! Allocation
allocate(smask_c1int(nam%nc1,geom%nl0))
if (samp%sc3) allocate(smask_c1ac3int(samp%nc1a,nam%nc3,geom%nl0))

! Create file
write(mpl%info,'(a7,a)') '','Write sampling'
call mpl%flush
write(filename,'(a,a,i4.4,a,i4.4)') trim(nam%prefix),'_sampling_',mpl%nproc,'-',mpl%myproc
call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call nam%write(mpl,ncid)

! Define dimensions
call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1',nam%nc1,nc1_id))
call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1a',samp%nc1a,nc1a_id))
if (samp%sc2) call mpl%ncerr(subr,nf90_def_dim(ncid,'nc2',nam%nc2,nc2_id))
if (samp%sc3) call mpl%ncerr(subr,nf90_def_dim(ncid,'nc3',nam%nc3,nc3_id))

! Define variables
call mpl%ncerr(subr,nf90_def_var(ncid,trim(samp%name)//'_c1_to_c0',nf90_int,(/nc1_id/),c1_to_c0_id))
call mpl%ncerr(subr,nf90_put_att(ncid,c1_to_c0_id,'_FillValue',mpl%msv%vali))
call mpl%ncerr(subr,nf90_def_var(ncid,trim(samp%name)//'_smask_c1',nf90_int,(/nc1_id,nl0_id/),smask_c1_id))
call mpl%ncerr(subr,nf90_put_att(ncid,smask_c1_id,'_FillValue',mpl%msv%vali))
if (samp%sc2) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(samp%name)//'_c2_to_c1',nf90_int,(/nc2_id/),c2_to_c1_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,c2_to_c1_id,'_FillValue',mpl%msv%vali))
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(samp%name)//'_c2_to_c0',nf90_int,(/nc2_id/),c2_to_c0_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,c2_to_c0_id,'_FillValue',mpl%msv%vali))
end if
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(samp%name)//'_c1ac3_to_c0',nf90_int,(/nc1a_id,nc3_id/),c1ac3_to_c0_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,c1ac3_to_c0_id,'_FillValue',mpl%msv%vali))
   call mpl%ncerr(subr,nf90_def_var(ncid,trim(samp%name)//'_smask_c1ac3',nf90_int,(/nc1a_id,nc3_id,nl0_id/),smask_c1ac3_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,smask_c1ac3_id,'_FillValue',mpl%msv%vali))
end if

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Convert data
do il0=1,geom%nl0
   do ic1=1,nam%nc1
      if (samp%smask_c1(ic1,il0)) then
         smask_c1int(ic1,il0) = 1
      else
         smask_c1int(ic1,il0) = 0
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
call mpl%ncerr(subr,nf90_put_var(ncid,c1_to_c0_id,samp%c1_to_c0))
call mpl%ncerr(subr,nf90_put_var(ncid,smask_c1_id,smask_c1int))
if (samp%sc2) then
   call mpl%ncerr(subr,nf90_put_var(ncid,c2_to_c1_id,samp%c2_to_c1))
   call mpl%ncerr(subr,nf90_put_var(ncid,c2_to_c0_id,samp%c2_to_c0))
end if
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_put_var(ncid,c1ac3_to_c0_id,samp%c1ac3_to_c0))
   call mpl%ncerr(subr,nf90_put_var(ncid,smask_c1ac3_id,smask_c1ac3int))
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! Release memory
deallocate(smask_c1int)
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
integer :: il0,ic1,jc3,ic1a,ic2,ic2a,ic2b,jc1,jc1d,jc1e,ic0,jc0,j,nc1max_local,nc1max_vbal
integer :: ncid,nl0_id,nc1a_id,nc3_id,nc2a_id,nc2b_id,nc1max_local_id,nc1max_vbal_id
integer :: lon_id,lat_id,lon_local_id,lat_local_id,lon_vbal_id,lat_vbal_id
integer :: lon_ori_id,lat_ori_id,lon_local_ori_id,lat_local_ori_id,lon_vbal_ori_id,lat_vbal_ori_id
real(kind_real),allocatable :: lon_ori(:),lat_ori(:),lon(:,:,:),lat(:,:,:)
real(kind_real),allocatable :: lon_local_ori(:),lat_local_ori(:),lon_local(:,:),lat_local(:,:)
real(kind_real),allocatable :: lon_vbal_ori(:),lat_vbal_ori(:),lon_vbal(:,:),lat_vbal(:,:)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'samp_write_grids'

! Create files
write(mpl%info,'(a7,a)') '','Write sampling grids'
call mpl%flush
write(filename,'(a,a,i4.4,a,i4.4)') trim(nam%prefix),'_sampling_grids_',mpl%nproc,'-',mpl%myproc
call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Write namelist parameters
call nam%write(mpl,ncid)

! Define dimensions
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc3',nam%nc3,nc3_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1a',samp%nc1a,nc1a_id))
end if
if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   nc1max_local = 0
   do ic2a=1,samp%nc2a
      nc1max_local = max(count(samp%local_mask(:,ic2a)),nc1max_local)
   end do
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1max_local',nc1max_local,nc1max_local_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc2a',samp%nc2a,nc2a_id))
end if
if (trim(samp%name)=='vbal') then
   nc1max_vbal = 0
   do ic2b=1,samp%nc2b
      nc1max_vbal = max(count(samp%vbal_mask(:,ic2b)),nc1max_vbal)
   end do
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc1max_vbal',nc1max_vbal,nc1max_vbal_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nc2b',samp%nc2b,nc2b_id))
end if

! Define variables
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_ori',nc_kind_real,(/nc1a_id/),lon_ori_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_ori_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_ori',nc_kind_real,(/nc1a_id/),lat_ori_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_ori_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon',nc_kind_real,(/nc1a_id,nc3_id,nl0_id/),lon_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat',nc_kind_real,(/nc1a_id,nc3_id,nl0_id/),lat_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',mpl%msv%valr))
end if
if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_local_ori',nc_kind_real,(/nc2a_id/),lon_local_ori_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_local_ori_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_local_ori',nc_kind_real,(/nc2a_id/),lat_local_ori_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_local_ori_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_local',nc_kind_real,(/nc1max_local_id,nc2a_id/),lon_local_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_local_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_local',nc_kind_real,(/nc1max_local_id,nc2a_id/),lat_local_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_local_id,'_FillValue',mpl%msv%valr))
end if
if (trim(samp%name)=='vbal') then
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_vbal_ori',nc_kind_real,(/nc2b_id/),lon_vbal_ori_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_vbal_ori_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_vbal_ori',nc_kind_real,(/nc2b_id/),lat_vbal_ori_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_vbal_ori_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lon_vbal',nc_kind_real,(/nc1max_vbal_id,nc2b_id/),lon_vbal_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lon_vbal_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'lat_vbal',nc_kind_real,(/nc1max_vbal_id,nc2b_id/),lat_vbal_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,lat_vbal_id,'_FillValue',mpl%msv%valr))
end if

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Convert data
if (samp%sc3) then
   ! Allocation
   allocate(lon_ori(samp%nc1a))
   allocate(lat_ori(samp%nc1a))
   allocate(lon(samp%nc1a,nam%nc3,geom%nl0))
   allocate(lat(samp%nc1a,nam%nc3,geom%nl0))

   ! Origin points
   do ic1a=1,samp%nc1a
      ic1 = samp%c1a_to_c1(ic1a)
      ic0 = samp%c1_to_c0(ic1)
      lon_ori(ic1a) = geom%lon_c0(ic0)*rad2deg
      lat_ori(ic1a) = geom%lat_c0(ic0)*rad2deg
   end do

   ! Distant points
   lon = mpl%msv%valr
   lat = mpl%msv%valr
   do il0=1,geom%nl0
      do jc3=1,nam%nc3
         do ic1a=1,samp%nc1a
            if (samp%smask_c1ac3(ic1a,jc3,il0)) then
               ic0 = samp%c1ac3_to_c0(ic1a,jc3)
               lon(ic1a,jc3,il0) = geom%lon_c0(ic0)*rad2deg
               lat(ic1a,jc3,il0) = geom%lat_c0(ic0)*rad2deg
            end if
         end do
      end do
   end do
end if
if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   ! Allocation
   allocate(lon_local_ori(samp%nc2a))
   allocate(lat_local_ori(samp%nc2a))
   allocate(lon_local(nc1max_local,samp%nc2a))
   allocate(lat_local(nc1max_local,samp%nc2a))

   ! Origin points
   do ic2a=1,samp%nc2a
      ic2 = samp%c2a_to_c2(ic2a)
      ic0 = samp%c2_to_c0(ic2)
      lon_local_ori(ic2a) = geom%lon_c0(ic0)*rad2deg
      lat_local_ori(ic2a) = geom%lat_c0(ic0)*rad2deg
   end do

   ! Distant points
   lon_local = mpl%msv%valr
   lat_local = mpl%msv%valr
   do ic2a=1,samp%nc2a
      j = 0
      do jc1d=1,samp%nc1d
         jc1 = samp%c1d_to_c1(jc1d)
         if (samp%local_mask(jc1,ic2a)) then
            j = j+1
            jc0 = samp%c1_to_c0(jc1)
            lon_local(j,ic2a) = geom%lon_c0(jc0)*rad2deg
            lat_local(j,ic2a) = geom%lat_c0(jc0)*rad2deg
         end if
      end do
   end do
end if
if (trim(samp%name)=='vbal') then
   ! Allocation
   allocate(lon_vbal_ori(samp%nc2b))
   allocate(lat_vbal_ori(samp%nc2b))
   allocate(lon_vbal(nc1max_vbal,samp%nc2b))
   allocate(lat_vbal(nc1max_vbal,samp%nc2b))

   ! Origin points
   do ic2b=1,samp%nc2b
      ic2 = samp%c2b_to_c2(ic2b)
      ic0 = samp%c2_to_c0(ic2)
      lon_vbal_ori(ic2b) = geom%lon_c0(ic0)*rad2deg
      lat_vbal_ori(ic2b) = geom%lat_c0(ic0)*rad2deg
   end do

   ! Distant points
   lon_vbal = mpl%msv%valr
   lat_vbal = mpl%msv%valr
   do ic2b=1,samp%nc2b
      j = 0
      do jc1e=1,samp%nc1e
         jc1 = samp%c1e_to_c1(jc1e)
         if (samp%vbal_mask(jc1,ic2b)) then
            j = j+1
            jc0 = samp%c1_to_c0(jc1)
            lon_vbal(j,ic2b) = geom%lon_c0(jc0)*rad2deg
            lat_vbal(j,ic2b) = geom%lat_c0(jc0)*rad2deg
         end if
      end do
   end do
end if

! Write variables
if (samp%sc3) then
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_ori_id,lon_ori))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_ori_id,lat_ori))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,lon))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,lat))
end if
if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_local_ori_id,lon_local_ori))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_local_ori_id,lat_local_ori))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_local_id,lon_local))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_local_id,lat_local))
end if
if (trim(samp%name)=='vbal') then
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_vbal_ori_id,lon_vbal_ori))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_vbal_ori_id,lat_vbal_ori))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_vbal_id,lon_vbal))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_vbal_id,lat_vbal))
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! Release memory
if (samp%sc3) then
   deallocate(lon_ori)
   deallocate(lat_ori)
   deallocate(lon)
   deallocate(lat)
end if
if ((trim(samp%name)=='hdiag').and.nam%local_diag) then
   deallocate(lon_local_ori)
   deallocate(lat_local_ori)
   deallocate(lon_local)
   deallocate(lat_local)
end if
if (trim(samp%name)=='vbal') then
   deallocate(lon_vbal_ori)
   deallocate(lat_vbal_ori)
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
integer :: ic0,jc0,il0,jc3,ildwv,jldwv,ival,nn_index(geom%nc0),nc1_valid
real(kind_real),allocatable :: ldwv_to_lon(:),ldwv_to_lat(:)
character(len=8) :: ivalformat
character(len=1024) :: color
character(len=1024),parameter :: subr = 'samp_compute_c1'

! Set sampling name
samp%name = trim(sname)

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
      call geom%index_from_lonlat(mpl,nam%lon_ldwv(ildwv),nam%lat_ldwv(ildwv),0,samp%ldwv_to_c0a(ildwv),samp%ldwv_to_proc(ildwv))
      if (mpl%msv%is(samp%ldwv_to_proc(ildwv))) then
         ldwv_to_lon(ildwv) = mpl%msv%valr
         ldwv_to_lat(ildwv) = mpl%msv%valr
      else
         if (samp%ldwv_to_proc(ildwv)==mpl%myproc) then
            ldwv_to_lon(ildwv) = geom%lon_c0a(samp%ldwv_to_c0a(ildwv))
            ldwv_to_lat(ildwv) = geom%lat_c0a(samp%ldwv_to_c0a(ildwv))
         end if
         call mpl%f_comm%broadcast(ldwv_to_lon(ildwv),samp%ldwv_to_proc(ildwv)-1)
         call mpl%f_comm%broadcast(ldwv_to_lat(ildwv),samp%ldwv_to_proc(ildwv)-1)
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

! Read or compute sampling data
samp%new_sampling = .true.
if (nam%sam_read) then
   call samp%read(mpl,nam,geom)
   if (samp%new_sampling) nam%sam_write = .true.
end if

if (samp%new_sampling) then
   ! Compute zero separation sampling
   write(mpl%info,'(a7,a,i5,a)') '','Compute zero separation sampling (nc1 = ',nam%nc1,')'
   call mpl%flush
   call samp%compute_zs(mpl,rng,nam,geom)

   ! Compute MPI distribution, halo A, subset Sc1
   write(mpl%info,'(a7,a)') '','Compute MPI distribution, halos A, subset Sc1'
   call mpl%flush
   call samp%compute_mpi_c1a(mpl,nam,geom)

   if (samp%sc3) then
      ! Compute positive separation sampling
      write(mpl%info,'(a7,a,i5,a)') '','Compute positive separation sampling (nc3 = ',nam%nc3,')'
      call mpl%flush
      call samp%compute_ps(mpl,rng,nam,geom)
   end if

   ! Check sampling mask
   call samp%check_mask(mpl,nam,geom)

   if (samp%sc2) then
      ! Compute sampling, subset Sc2
      write(mpl%info,'(a7,a,i5,a)') '','Compute sampling, subset Sc2 (nc2 = ',nam%nc2,')'
      call mpl%flush
      call samp%compute_c2(mpl,rng,nam,geom)
   end if
end if

if (samp%sc2) then
   ! Compute mesh sampling, subset Sc2
   write(mpl%info,'(a7,a)') '','Compute sampling mesh, subset Sc2'
   call mpl%flush
   call samp%compute_mesh_c2(mpl,rng,nam,geom)

   ! Compute MPI distribution, halo A, subset Sc2
   write(mpl%info,'(a7,a)') '','Compute MPI distribution, halos A, subset Sc2'
   call mpl%flush
   call samp%compute_mpi_c2a(mpl,nam,geom)

   ! Compute MPI distribution, halo B, subset Sc2
   write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo B, subset Sc2'
   call mpl%flush
   call samp%compute_mpi_b(mpl,rng,nam,geom)
else
   ! Compute MPI distribution, halo A, subset Sc2
   samp%nc2a = 0
end if

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
   call samp%compute_mpi_e(mpl,nam,geom)
end if

if ((((trim(samp%name)=='hdiag').and.(nam%local_diag.or.nam%adv_diag)).or.(trim(samp%name)=='lct')).and.(nam%diag_rhflt>0.0)) then
   ! Compute MPI distribution, halo F
   write(mpl%info,'(a7,a)') '','Compute MPI distribution, halo F'
   call mpl%flush
   call samp%compute_mpi_f(mpl,nam)
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
integer :: nsmask,nsmask_tot,ic0u,il0,ildwv,iv,its,ie,ncontig,ncontigmax,latmin,latmax
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

   ! Commnication
   call geom%com_AU%ext(mpl,geom%nl0,samp%smask_c0a,samp%smask_c0u)

   ! Other masks
   samp%smask_hor_c0a = any(samp%smask_c0a,dim=2)
   samp%smask_hor_c0u = any(samp%smask_c0u,dim=2)
   nc0_smask(0) = count(samp%smask_hor_c0a)
   nc0_smask(1:geom%nl0) = count(samp%smask_c0a,dim=1)
   call mpl%f_comm%allreduce(nc0_smask,samp%nc0_smask,fckit_mpi_sum())

   ! Check mask size
   if (samp%nc0_smask(0)==0) call mpl%abort(subr,'no more points in the sampling mask')
else
   ! Copy geometry mask
   write(mpl%info,'(a7,a)') '','Copy geometry mask'
   call mpl%flush
   samp%smask_c0u = geom%gmask_c0u
   samp%smask_c0a = geom%gmask_c0a
   samp%smask_hor_c0u = geom%gmask_hor_c0u
   samp%smask_hor_c0a = geom%gmask_hor_c0a
   samp%nc0_smask = geom%nc0_gmask
end if

! Print results
write(mpl%info,'(a7,a)') '','Sampling valid points (% of domain mask):'
call mpl%flush
do il0=1,geom%nl0
   write(mpl%info,'(a10,a,i3,a,f5.1,a)') '','Level ',nam%levs(il0),' ~> ',100.0*real(samp%nc0_smask(il0),kind_real) &
 & /real(geom%nc0_gmask(il0),kind_real),'%'
   call mpl%flush
end do

end subroutine samp_compute_mask

! TODO: HERE

!----------------------------------------------------------------------
! Subroutine: samp_compute_zs
! Purpose: compute zero-separation sampling
!----------------------------------------------------------------------
subroutine samp_compute_zs(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic0,ic0a,ic1,il0i,ildwv
real(kind_real) :: rh_c0a(geom%nc0a)
logical :: mask_hor_c0a(geom%nc0a)
character(len=1024),parameter :: subr = 'samp_compute_zs'

! Initialization
samp%c1_to_c0 = mpl%msv%vali

! Compute subset
if (samp%nc0_smask(0)>nam%nc1) then
   write(mpl%info,'(a7,a)') '','Compute horizontal subset C1: '
   call mpl%flush(.false.)
   select case (trim(nam%draw_type))
   case ('random_uniform')
      ! Random draw
      rh_c0a = 1.0
   case ('random_coast')
      ! More points around coasts
      if (all(geom%gmask_c0)) call mpl%abort(subr,'random_coast is not relevant if there is no coast')
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
   mask_hor_c0a = samp%smask_hor_c0a
   do ildwv=1,nam%nldwv
      ic0 = samp%ldwv_to_c0(ildwv)
      samp%c1_to_c0(ildwv) = ic0
      if (mpl%myproc==geom%c0_to_proc(ic0)) then
         ic0a = geom%c0_to_c0a(ic0)
         mask_hor_c0a(ic0a) = .false.
      end if
   end do

   ! Initialize sampling 
   call initialize_sampling(mpl,rng,maxval(geom%area),geom%nc0a,geom%lon_c0a,geom%lat_c0a,mask_hor_c0a,rh_c0a,geom%c0a_to_c0, &
 & nam%ntry,nam%nrep,nam%nc1-nam%nldwv,samp%c1_to_c0(nam%nldwv+1:nam%nc1),n_uni=geom%nc0,uni_to_proc=geom%c0_to_proc, &
 & uni_to_loc=geom%c0_to_c0a,tree_uni=geom%tree)
elseif (samp%nc0_smask(0)==nam%nc1) then
   ! Keep all remaining points
   ic1 = 0
   do ic0=1,geom%nc0
      if (samp%smask_hor_c0(ic0)) then
         ic1 = ic1+1
         samp%c1_to_c0(ic1) = ic0
      end if
   end do
else
   ! Not enough points remaining in the sampling mask
   call mpl%abort(subr,'not enough points remaining in sampling mask')
end if

end subroutine samp_compute_zs

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_c1a
! Purpose: compute sampling MPI distribution, halo A, subset Sc1
!----------------------------------------------------------------------
subroutine samp_compute_mpi_c1a(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic0a,ic0,ic1a,ic1,il0

! Allocation
allocate(samp%lcheck_c0a(geom%nc0))
allocate(samp%lcheck_c1a(nam%nc1))

! Define halo A
samp%lcheck_c0a = .false.
samp%lcheck_c1a = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   if (geom%c0_to_proc(ic0)==mpl%myproc) samp%lcheck_c0a(ic0) = .true.
end do
do ic1=1,nam%nc1
   ic0 = samp%c1_to_c0(ic1)
   if (geom%c0_to_proc(ic0)==mpl%myproc) samp%lcheck_c1a(ic1) = .true.
end do
samp%nc1a = count(samp%lcheck_c1a)

! Allocation
allocate(samp%c1a_to_c1(samp%nc1a))
allocate(samp%c1_to_c1a(nam%nc1))
allocate(samp%c1_to_proc(nam%nc1))
allocate(samp%c1al0_check(samp%nc1a,geom%nl0))

! Global-local conversion for halo A
ic1a = 0
do ic1=1,nam%nc1
   if (samp%lcheck_c1a(ic1)) then
      ic1a = ic1a+1
      samp%c1a_to_c1(ic1a) = ic1
      samp%c1_to_c1a(ic1) = ic1a
   end if
end do
call mpl%glb_to_loc_index(samp%nc1a,samp%c1a_to_c1,nam%nc1,samp%c1_to_c1a)

! Subset Sc1 to processor
samp%c1_to_proc = geom%c0_to_proc(samp%c1_to_c0)

! Define whether this point should be checked for boundaries
do il0=1,geom%nl0
   do ic1a=1,samp%nc1a
      ic1 = samp%c1a_to_c1(ic1a)
      ic0 = samp%c1_to_c0(ic1)
      ic0a = geom%c0_to_c0a(ic0)
      samp%c1al0_check(ic1a,il0) = (nam%mask_check.and.samp%smask_c0a(ic0a,il0))
   end do
end do

! Print results
write(mpl%info,'(a10,a,i8)') '','nc1 = ',nam%nc1
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1a = ',samp%nc1a
call mpl%flush

end subroutine samp_compute_mpi_c1a

!----------------------------------------------------------------------
! Subroutine: samp_compute_ps
! Purpose: compute positive separation sampling
!----------------------------------------------------------------------
subroutine samp_compute_ps(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: jc3,ic1a,ic1,ir,irtmp,ic0,jc0,icinf,icsup,ictest,nn_index(nam%nc3),il0
real(kind_real) :: d,nn_dist(nam%nc3)
logical :: valid,found
character(len=1024),parameter :: subr = 'samp_compute_ps'

! Allocation
allocate(samp%c1ac3_to_c0(samp%nc1a,nam%nc3))
allocate(samp%smask_c1ac3(samp%nc1a,nam%nc3,geom%nl0))

! Initialization
samp%c1ac3_to_c0 = mpl%msv%vali

if (trim(samp%name)=='hdiag') then
   ! First class
   samp%c1ac3_to_c0(:,1) = samp%c1_to_c0(samp%c1a_to_c1)

   ! Synchronize random number generator
   call rng%sync(mpl)

   if (nam%nc3>1) then
      ! Initialization
      write(mpl%info,'(a7,a)') '','Compute positive separation sampling: '
      call mpl%flush(.false.)
      call mpl%prog_init(nam%nc3*samp%nc1a)
      do ic1a=1,samp%nc1a
         mpl%done((ic1a-1)*nam%nc3+1) = .true.
      end do
      ir = 0

      ! Sample classes of positive separation
      do while ((.not.all(mpl%done)).and.(ir<=nam%irmax))
         ! Define a random geographical point
         call geom%rand_point(mpl,rng,0,iproc,jc0a,irtmp)
         ir = ir+irtmp
         valid = samp%smask_hor_c0(jc0) ! TODO : HERE

         if (valid) then
            ! Fill classes
            !$omp parallel do schedule(static) private(ic1a,ic1,ic0,d,jc3,icinf,icsup,found,ictest)
            do ic1a=1,samp%nc1a
               ! Index
               ic1 = samp%c1a_to_c1(ic1a)

               ! Check if there is a valid first point
               if (mpl%msv%isnot(samp%c1_to_c0(ic1))) then
                  ! Compute the distance
                  ic0 = samp%c1_to_c0(ic1)
                  call sphere_dist(geom%lon_c0(ic0),geom%lat_c0(ic0),geom%lon_c0(jc0),geom%lat_c0(jc0),d)

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
                     if ((jc3/=1).and.(mpl%msv%is(samp%c1ac3_to_c0(ic1a,jc3)))) then
                        samp%c1ac3_to_c0(ic1a,jc3) = jc0
                        mpl%done((ic1a-1)*nam%nc3+jc3) = .true.
                     end if
                  end if
               end if
            end do
            !$omp end parallel do

            ! Update
            call mpl%prog_print
         end if
      end do
      call mpl%prog_final
   end if

   ! Reset seeds
   call rng%reseed(mpl)
elseif (trim(samp%name)=='lct') then
   ! Initialization
   write(mpl%info,'(a7,a)') '','Compute LCT sampling: '
   call mpl%flush(.false.)
   call mpl%prog_init(samp%nc1a)

   do ic1a=1,samp%nc1a
      ! Index
      ic1 = samp%c1a_to_c1(ic1a)

      ! Check location validity
      if (mpl%msv%isnot(samp%c1_to_c0(ic1))) then
         ! Find neighbors
         call geom%tree%find_nearest_neighbors(geom%lon_c0(samp%c1_to_c0(ic1)),geom%lat_c0(samp%c1_to_c0(ic1)), &
       & nam%nc3,nn_index,nn_dist)

         ! Copy neighbor index
         do jc3=1,nam%nc3
            jc0 = nn_index(jc3)
            samp%c1ac3_to_c0(ic1a,jc3) = nn_index(jc3)
         end do
      end if

      ! Update
      call mpl%prog_print(ic1a)
   end do
   call mpl%prog_final

   ! Define whether this point should be checked for boundaries
   do il0=1,geom%nl0
      do ic1a=1,samp%nc1a
         if (samp%c1al0_check(ic1a,il0)) samp%c1al0_check(ic1a,il0) = any(.not.samp%smask_c0(samp%c1ac3_to_c0(ic1a,:),il0))
      end do
   end do
endif

end subroutine samp_compute_ps

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
integer :: il0,jc3,ic1a,ic1,ic0,jc0
logical :: valid

! First point
do il0=1,geom%nl0
   samp%smask_c1(:,il0) = samp%smask_c0(samp%c1_to_c0,il0)
end do

if (samp%sc3) then
   ! Check sampling mask
   write(mpl%info,'(a7,a)') '','Check sampling mask: '
   call mpl%flush(.false.)

   ! Second point
   call mpl%prog_init(samp%nc1a)
   do ic1a=1,samp%nc1a
      ! Indices
      ic1 = samp%c1a_to_c1(ic1a)
      ic0 = samp%c1_to_c0(ic1)

      do il0=1,geom%nl0
         do jc3=1,nam%nc3
            ! Index
            jc0 = samp%c1ac3_to_c0(ic1a,jc3)

            ! Check point index
            valid = mpl%msv%isnot(ic0).and.mpl%msv%isnot(jc0)

            ! Check sampling mask
            if (valid) valid = samp%smask_c0(ic0,il0).and.samp%smask_c0(jc0,il0)

            ! Check mask boundaries
            if (valid.and.nam%mask_check) call geom%check_arc(mpl,il0,geom%lon_c0(ic0),geom%lat_c0(ic0),geom%lon_c0(jc0), &
          & geom%lat_c0(jc0),valid)

            ! Copy validity
            samp%smask_c1ac3(ic1a,jc3,il0) = valid
         end do
      end do

      ! Update
      call mpl%prog_print(ic1a)
   end do
   call mpl%prog_final
end if

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
integer :: ic1a,ic1,ildwv,ic0
real(kind_real),allocatable :: lon_c1a(:),lat_c1a(:)
real(kind_real),allocatable :: rh_c1a(:)
logical,allocatable :: mask_c1a(:)

! Allocation
allocate(lon_c1a(samp%nc1a))
allocate(lat_c1a(samp%nc1a))
allocate(mask_c1a(samp%nc1a))
allocate(rh_c1a(samp%nc1a))

! Define subsampling
write(mpl%info,'(a7,a)') '','Define subsampling:'
call mpl%flush(.false.)
lon_c1a = geom%lon_c0(samp%c1_to_c0(samp%c1a_to_c1))
lat_c1a = geom%lat_c0(samp%c1_to_c0(samp%c1a_to_c1))
mask_c1a = any(samp%smask_c1(samp%c1a_to_c1,:),dim=2)
rh_c1a = 1.0
do ildwv=1,nam%nldwv
   ic1 = ildwv
   samp%c2_to_c1(ic1) = ic1
   ic0 = samp%c1_to_c0(ic1)
   if (mpl%myproc==geom%c0_to_proc(ic0)) then
      ic1a = samp%c1_to_c1a(ic1)
      mask_c1a(ic1a) = .false.
   end if
end do
call initialize_sampling(mpl,rng,maxval(geom%area),samp%nc1a,lon_c1a,lat_c1a,mask_c1a,rh_c1a,samp%c1a_to_c1, &
 & nam%ntry,nam%nrep,nam%nc2-nam%nldwv,samp%c2_to_c1(nam%nldwv+1:nam%nc2))
samp%c2_to_c0 = samp%c1_to_c0(samp%c2_to_c1)

! Release memory
deallocate(lon_c1a)
deallocate(lat_c1a)
deallocate(mask_c1a)
deallocate(rh_c1a)

end subroutine samp_compute_c2

!----------------------------------------------------------------------
! Subroutine: samp_compute_mesh_c2
! Purpose: compute sampling mesh, subset Sc2
!----------------------------------------------------------------------
subroutine samp_compute_mesh_c2(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(inout) :: nam    ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: il0,ic1,ic2

! Initialization
samp%lon_c2 = geom%lon_c0(samp%c2_to_c0)
samp%lat_c2 = geom%lat_c0(samp%c2_to_c0)
do il0=1,geom%nl0
   do ic2=1,nam%nc2
      ic1 = samp%c2_to_c1(ic2)
      samp%smask_c2(ic2,il0) = samp%smask_c1(ic1,il0)
   end do
end do

if (nam%adv_diag) then
   ! Initialization
   call samp%mesh%init(mpl,rng,samp%lon_c2,samp%lat_c2,.true.)

   ! Compute triangles list
   write(mpl%info,'(a7,a)') '','Compute triangles list '
   call mpl%flush
   call samp%mesh%trlist(mpl)

   ! Find boundary nodes
   write(mpl%info,'(a7,a)') '','Find boundary nodes'
   call mpl%flush
   call samp%mesh%bnodes(mpl)
end if

end subroutine samp_compute_mesh_c2

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
integer :: ic0,ic1a,ic1,ic2a,ic2,il0
type(tree_type) :: tree

! Allocation
allocate(samp%lcheck_c2a(nam%nc2))

! Define halo A
samp%lcheck_c2a = .false.
do ic2=1,nam%nc2
   ic0 = samp%c2_to_c0(ic2)
   if (geom%c0_to_proc(ic0)==mpl%myproc) samp%lcheck_c2a(ic2) = .true.
end do
samp%nc2a = count(samp%lcheck_c2a)

! Allocation
allocate(samp%c2_to_proc(nam%nc2))
allocate(samp%c2a_to_c2(samp%nc2a))
allocate(samp%c2_to_c2a(nam%nc2))
allocate(samp%c2a_to_c1a(samp%nc2a))

! Global-local conversion for halo A
ic2a = 0
do ic2=1,nam%nc2
   if (samp%lcheck_c2a(ic2)) then
      ic2a = ic2a+1
      samp%c2a_to_c2(ic2a) = ic2
   end if
end do
call mpl%glb_to_loc_index(samp%nc2a,samp%c2a_to_c2,nam%nc2,samp%c2_to_c2a)

! Local conversion on halo A
do ic2a=1,samp%nc2a
   ic2 = samp%c2a_to_c2(ic2a)
   ic1 = samp%c2_to_c1(ic2)
   ic1a = samp%c1_to_c1a(ic1)
   samp%c2a_to_c1a(ic2a) = ic1a
end do

! MPI splitting on subset Sc2
do ic2=1,nam%nc2
   ic0 = samp%c2_to_c0(ic2)
   samp%c2_to_proc(ic2) = geom%c0_to_proc(ic0)
end do

! Mask on subset Sc2, halo A
allocate(samp%smask_c2a(samp%nc2a,geom%nl0))
do il0=1,geom%nl0
   do ic2a=1,samp%nc2a
      ic2 = samp%c2a_to_c2(ic2a)
      samp%smask_c2a(ic2a,il0) = samp%smask_c2(ic2,il0)
   end do
end do

! Find nearest neighbors

! Allocation
allocate(samp%nn_c2a_index(nam%nc2,samp%nc2a))
allocate(samp%nn_c2a_dist(nam%nc2,samp%nc2a))
call tree%alloc(mpl,nam%nc2)

! Initialization
call tree%init(samp%lon_c2,samp%lat_c2)

! Find nearest neighbors
do ic2a=1,samp%nc2a
   ic2 = samp%c2a_to_c2(ic2a)
   ic0 = samp%c2_to_c0(ic2)
   call tree%find_nearest_neighbors(geom%lon_c0(ic0),geom%lat_c0(ic0),nam%nc2,samp%nn_c2a_index(:,ic2a),samp%nn_c2a_dist(:,ic2a))
end do

! Release memory
call tree%dealloc

! Print results
write(mpl%info,'(a10,a,i8)') '','nc2 = ',nam%nc2
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc2a = ',samp%nc2a
call mpl%flush

end subroutine samp_compute_mpi_c2a

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_b
! Purpose: compute sampling MPI distribution, halo B
!----------------------------------------------------------------------
subroutine samp_compute_mpi_b(samp,mpl,rng,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic2a,ic2b,ic2,jc2,i_s,il0i
integer,allocatable :: c2_to_c2b(:),c2a_to_c2b(:)
logical :: lcheck_c2b(nam%nc2)

! Allocation
allocate(samp%h(geom%nl0i))

! Compute interpolation
do il0i=1,geom%nl0i
   write(samp%h(il0i)%prefix,'(a,i3.3)') 'h_',il0i
   call samp%h(il0i)%interp(mpl,rng,nam,geom,il0i,nam%nc2,samp%lon_c2,samp%lat_c2,samp%smask_c2(:,il0i),geom%nc0a, &
 & geom%lon_c0a,geom%lat_c0a,geom%gmask_c0a(:,il0i),7)
end do

! Define halo B
lcheck_c2b = samp%lcheck_c2a
do il0i=1,geom%nl0i
   do i_s=1,samp%h(il0i)%n_s
      jc2 = samp%h(il0i)%col(i_s)
      lcheck_c2b(jc2) = .true.
   end do
end do
samp%nc2b = count(lcheck_c2b)

! Allocation
allocate(samp%c2b_to_c2(samp%nc2b))
allocate(c2_to_c2b(nam%nc2))
allocate(c2a_to_c2b(samp%nc2a))

! Global-local conversion for halo B
ic2b = 0
c2_to_c2b = mpl%msv%vali
do ic2=1,nam%nc2
   if (lcheck_c2b(ic2)) then
      ic2b = ic2b+1
      samp%c2b_to_c2(ic2b) = ic2
      c2_to_c2b(ic2) = ic2b
   end if
end do

! Halos A-B conversion
do ic2a=1,samp%nc2a
   ic2 = samp%c2a_to_c2(ic2a)
   ic2b = c2_to_c2b(ic2)
   c2a_to_c2b(ic2a) = ic2b
end do

do il0i=1,geom%nl0i
   ! Local interpolation source
   samp%h(il0i)%n_src = samp%nc2b
   do i_s=1,samp%h(il0i)%n_s
      samp%h(il0i)%col(i_s) = c2_to_c2b(samp%h(il0i)%col(i_s))
   end do
end do

! Setup communications
call samp%com_AB%setup(mpl,'com_AB',nam%nc2,samp%nc2a,samp%nc2b,samp%nc2a,samp%c2b_to_c2,c2a_to_c2b,samp%c2_to_proc,samp%c2_to_c2a)

! Release memory
deallocate(c2_to_c2b)
deallocate(c2a_to_c2b)

! Print results
write(mpl%info,'(a10,a,i8)') '','nc2b = ',samp%nc2b
do il0i=1,geom%nl0i
   write(mpl%info,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',samp%h(il0i)%n_s
   call mpl%flush
end do

end subroutine samp_compute_mpi_b

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
integer :: jc3,ic0,ic0a,ic0c,ic1,ic1a,its,il0,i_s
integer,allocatable :: c0c_to_c0(:),c0a_to_c0c(:)
real(kind_real),allocatable :: lon_c1a(:),lat_c1a(:)
logical :: lcheck_c0c(geom%nc0)
logical,allocatable :: mask_c1a(:)

if ((trim(samp%name)=='hdiag').and.nam%adv_diag) then
   write(mpl%info,'(a7,a)') '','Compute advection interpolation'
   call mpl%flush

   ! Allocation
   samp%nc1a = count(geom%c0_to_proc(samp%c1_to_c0)==mpl%myproc)
   allocate(lon_c1a(samp%nc1a))
   allocate(lat_c1a(samp%nc1a))
   allocate(mask_c1a(samp%nc1a))
   allocate(samp%d(geom%nl0,nam%nts))

   ! Compute interpolation
   do its=1,nam%nts
      do il0=1,geom%nl0
         do ic1a=1,samp%nc1a
            ic1 = samp%c1a_to_c1(ic1a)
            ic0 = samp%c1_to_c0(ic1)
            ic0a = geom%c0_to_c0a(ic0)
            lon_c1a(ic1a) = samp%adv_lon(ic0a,il0,its)
            lat_c1a(ic1a) = samp%adv_lat(ic0a,il0,its)
            mask_c1a(ic1a) = samp%smask_c1(ic1,il0)
         end do
         write(samp%d(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
         call samp%d(il0,its)%interp(mpl,rng,nam,geom,il0,geom%nc0,geom%lon_c0,geom%lat_c0,geom%gmask_c0(:,il0),samp%nc1a, &
       & lon_c1a,lat_c1a,mask_c1a,10)
      end do
   end do

   ! Release memory
   deallocate(lon_c1a)
   deallocate(lat_c1a)
   deallocate(mask_c1a)
end if

! Define halo C
lcheck_c0c = samp%lcheck_c0a
do jc3=1,nam%nc3
   do ic1a=1,samp%nc1a
      if (any(samp%smask_c1ac3(ic1a,jc3,:))) then
         ic0 = samp%c1ac3_to_c0(ic1a,jc3)
         lcheck_c0c(ic0) = .true.
      end if
   end do
end do
if ((trim(samp%name)=='hdiag').and.nam%adv_diag) then
   do its=1,nam%nts
      do il0=1,geom%nl0
         do i_s=1,samp%d(il0,its)%n_s
            ic0 = samp%d(il0,its)%col(i_s)
            ic1a = samp%d(il0,its)%row(i_s)
            ic1 = samp%c1a_to_c1(ic1a)
            if (samp%lcheck_c1a(ic1)) lcheck_c0c(ic0) = .true.
         end do
      end do
   end do
end if
samp%nc0c = count(lcheck_c0c)

! Global-local conversion for halo C
allocate(c0c_to_c0(samp%nc0c))
allocate(samp%c0_to_c0c(geom%nc0))
samp%c0_to_c0c = mpl%msv%vali
ic0c = 0
do ic0=1,geom%nc0
   if (lcheck_c0c(ic0)) then
      ic0c = ic0c+1
      c0c_to_c0(ic0c) = ic0
      samp%c0_to_c0c(ic0) = ic0c
   end if
end do

! Halos A-C conversion
allocate(c0a_to_c0c(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0c = samp%c0_to_c0c(ic0)
   c0a_to_c0c(ic0a) = ic0c
end do

if ((trim(samp%name)=='hdiag').and.nam%adv_diag) then
   ! Local interpolation source
   do its=1,nam%nts
      do il0=1,geom%nl0
         samp%d(il0,its)%n_src = samp%nc0c
         do i_s=1,samp%d(il0,its)%n_s
            samp%d(il0,its)%col(i_s) = samp%c0_to_c0c(samp%d(il0,its)%col(i_s))
         end do
      end do
   end do
end if

! Setup communications
call samp%com_AC%setup(mpl,'com_AC',geom%nc0,geom%nc0a,samp%nc0c,geom%nc0a,c0c_to_c0,c0a_to_c0c,geom%c0_to_proc,geom%c0_to_c0a)

! Release memory
deallocate(c0c_to_c0)
deallocate(c0a_to_c0c)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
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
integer :: ic2a,ic2,ic0,nn,i,ic1d,ic1,ic1a,jc1,npack,ipack,il0,jc3
integer :: c1_to_c1d(nam%nc1)
integer,allocatable :: nn_index(:),c1a_to_c1d(:),c1_to_proc(:)
real(kind_real) :: lon_c1(nam%nc1),lat_c1(nam%nc1)
logical :: mask_c1(nam%nc1),lcheck_c1d(nam%nc1)
logical,allocatable :: sbuf(:,:),rbuf(:,:)
type(tree_type) :: tree

! Allocation
allocate(samp%local_mask(nam%nc1,samp%nc2a))

! Initialization
mask_c1 = any(samp%smask_c1,dim=2)
lon_c1 = geom%lon_c0(samp%c1_to_c0)
lat_c1 = geom%lat_c0(samp%c1_to_c0)
samp%local_mask = .false.
lcheck_c1d = samp%lcheck_c1a

! Allocation
call tree%alloc(mpl,nam%nc1,mask=mask_c1)

! Initialization
call tree%init(lon_c1,lat_c1)

! Halo D
do ic2a=1,samp%nc2a
   ! Indices
   ic2 = samp%c2a_to_c2(ic2a)
   ic1 = samp%c2_to_c1(ic2)
   ic0 = samp%c2_to_c0(ic2)

   ! Count nearest neighbors
   call tree%count_nearest_neighbors(geom%lon_c0(ic0),geom%lat_c0(ic0),nam%local_rad,nn)

   ! Allocation
   allocate(nn_index(nn))

   ! Find nearest neighbors
   call tree%find_nearest_neighbors(geom%lon_c0(ic0),geom%lat_c0(ic0),nn,nn_index)

   ! Update masks
   samp%local_mask(ic1,ic2a) = .true.
   do i=1,nn
      jc1 = nn_index(i)
      samp%local_mask(jc1,ic2a) = .true.
      lcheck_c1d(jc1) = .true.
   end do

   ! Release memory
   deallocate(nn_index)
end do
samp%nc1d = count(lcheck_c1d)

! Release memory
call tree%dealloc

! Global <-> local conversions for fields

! Halo D
allocate(samp%c1d_to_c1(samp%nc1d))
c1_to_c1d = mpl%msv%vali
ic1d = 0
do ic1=1,nam%nc1
   if (lcheck_c1d(ic1)) then
      ic1d = ic1d+1
      samp%c1d_to_c1(ic1d) = ic1
      c1_to_c1d(ic1) = ic1d
   end if
end do

! Inter-halo conversions
allocate(c1a_to_c1d(samp%nc1a))
do ic1a=1,samp%nc1a
   ic1 = samp%c1a_to_c1(ic1a)
   ic1d = c1_to_c1d(ic1)
   c1a_to_c1d(ic1a) = ic1d
end do

! MPI splitting
allocate(c1_to_proc(nam%nc1))
do ic1=1,nam%nc1
   ic0 = samp%c1_to_c0(ic1)
   c1_to_proc(ic1) = geom%c0_to_proc(ic0)
end do

! Setup communications
call samp%com_AD%setup(mpl,'com_AD',nam%nc1,samp%nc1a,samp%nc1d,samp%nc1a,samp%c1d_to_c1,c1a_to_c1d,c1_to_proc,samp%c1_to_c1a)

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
deallocate(c1a_to_c1d)
deallocate(c1_to_proc)
deallocate(sbuf)
deallocate(rbuf)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1d =      ',samp%nc1d
call mpl%flush

end subroutine samp_compute_mpi_d

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_e
! Purpose: compute sampling MPI distribution, halo E
!----------------------------------------------------------------------
subroutine samp_compute_mpi_e(samp,mpl,nam,geom)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ic2b,ic2,ic0,nn,i,ic1e,ic1,ic1a,jc1,jc0
integer :: c1_to_c1e(nam%nc1)
integer,allocatable :: nn_index(:),c1a_to_c1e(:),c1_to_proc(:)
real(kind_real) :: lon_c1(nam%nc1),lat_c1(nam%nc1)
logical :: mask_c1(nam%nc1),lcheck_c1e(nam%nc1)
type(tree_type) :: tree
character(len=1024),parameter :: subr = 'samp_compute_mpi_e'

! Allocation
allocate(samp%vbal_mask(nam%nc1,samp%nc2b))

! Initialization
mask_c1 = any(samp%smask_c1,dim=2)
lon_c1 = geom%lon_c0(samp%c1_to_c0)
lat_c1 = geom%lat_c0(samp%c1_to_c0)
samp%vbal_mask = .false.
lcheck_c1e = samp%lcheck_c1a

! Allocation
call tree%alloc(mpl,nam%nc1,mask=mask_c1)

! Initialization
call tree%init(lon_c1,lat_c1)

! Halo E
do ic2b=1,samp%nc2b
   ! Indices
   ic2 = samp%c2b_to_c2(ic2b)
   ic1 = samp%c2_to_c1(ic2)
   ic0 = samp%c2_to_c0(ic2)

   ! Origin point
   samp%vbal_mask(ic1,ic2b) = .true.
   lcheck_c1e(ic1) = .true.

   if (nam%vbal_rad>0.0) then
      ! Count nearest neighbors
      call tree%count_nearest_neighbors(geom%lon_c0(ic0),geom%lat_c0(ic0),nam%vbal_rad,nn)

      ! Allocation
      allocate(nn_index(nn))

      ! Find nearest neighbors
      call tree%find_nearest_neighbors(geom%lon_c0(ic0),geom%lat_c0(ic0),nn,nn_index)

      ! Update masks
      do i=1,nn
         jc1 = nn_index(i)
         samp%vbal_mask(jc1,ic2b) = .true.
         lcheck_c1e(jc1) = .true.
      end do

      ! Release memory
      deallocate(nn_index)
   elseif (nam%vbal_dlat>0.0) then
      ! Update masks
      do jc1=1,nam%nc1
         jc0 = samp%c1_to_c0(jc1)
         if (abs(geom%lat_c0(ic0)-geom%lat_c0(jc0))<nam%vbal_dlat) then
            samp%vbal_mask(jc1,ic2b) = .true.
            lcheck_c1e(jc1) = .true.
         end if
      end do
   else
      call mpl%abort(subr,'vbal_rad or vbal_dlat should be positive')
   end if
end do
samp%nc1e = count(lcheck_c1e)

! Release memory
call tree%dealloc

! Global <-> local conversions for fields

! Halo E
allocate(samp%c1e_to_c1(samp%nc1e))
c1_to_c1e = mpl%msv%vali
ic1e = 0
do ic1=1,nam%nc1
   if (lcheck_c1e(ic1)) then
      ic1e = ic1e+1
      samp%c1e_to_c1(ic1e) = ic1
      c1_to_c1e(ic1) = ic1e
   end if
end do

! Inter-halo conversions
allocate(c1a_to_c1e(samp%nc1a))
do ic1a=1,samp%nc1a
   ic1 = samp%c1a_to_c1(ic1a)
   ic1e = c1_to_c1e(ic1)
   c1a_to_c1e(ic1a) = ic1e
end do

! MPI splitting
allocate(c1_to_proc(nam%nc1))
do ic1=1,nam%nc1
   ic0 = samp%c1_to_c0(ic1)
   c1_to_proc(ic1) = geom%c0_to_proc(ic0)
end do

! Setup communications
call samp%com_AE%setup(mpl,'com_AE',nam%nc1,samp%nc1a,samp%nc1e,samp%nc1a,samp%c1e_to_c1,c1a_to_c1e,c1_to_proc,samp%c1_to_c1a)

! Release memory
deallocate(c1a_to_c1e)
deallocate(c1_to_proc)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc1e =      ',samp%nc1e
call mpl%flush

end subroutine samp_compute_mpi_e

!----------------------------------------------------------------------
! Subroutine: samp_compute_mpi_f
! Purpose: compute sampling MPI distribution, halo F
!----------------------------------------------------------------------
subroutine samp_compute_mpi_f(samp,mpl,nam)

implicit none

! Passed variables
class(samp_type),intent(inout) :: samp ! Sampling
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist

! Local variables
integer :: ic2a,ic2f,ic2,jc2,kc2
integer,allocatable :: c2f_to_c2(:),c2a_to_c2f(:)
logical :: lcheck_c2f(nam%nc2)

! Define halo F
lcheck_c2f = samp%lcheck_c2a
do ic2a=1,samp%nc2a
   jc2 = 1
   do while (samp%nn_c2a_dist(jc2,ic2a)<nam%diag_rhflt)
      kc2 = samp%nn_c2a_index(jc2,ic2a)
      lcheck_c2f(kc2) = .true.
      jc2 = jc2+1
      if (jc2>nam%nc2) exit
   end do
end do
samp%nc2f = count(lcheck_c2f)

! Global-local conversion for halo F
allocate(c2f_to_c2(samp%nc2f))
allocate(samp%c2_to_c2f(nam%nc2))
samp%c2_to_c2f = mpl%msv%vali
ic2f = 0
do ic2=1,nam%nc2
   if (lcheck_c2f(ic2)) then
      ic2f = ic2f+1
      c2f_to_c2(ic2f) = ic2
      samp%c2_to_c2f(ic2) = ic2f
   end if
end do

! Halos A-F conversion
allocate(c2a_to_c2f(samp%nc2a))
do ic2a=1,samp%nc2a
   ic2 = samp%c2a_to_c2(ic2a)
   ic2f = samp%c2_to_c2f(ic2)
   c2a_to_c2f(ic2a) = ic2f
end do

! Setup communications
call samp%com_AF%setup(mpl,'com_AF',nam%nc2,samp%nc2a,samp%nc2f,samp%nc2a,c2f_to_c2,c2a_to_c2f,samp%c2_to_proc,samp%c2_to_c2a)

! Release memory
deallocate(c2f_to_c2)
deallocate(c2a_to_c2f)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc2f =       ',samp%nc2f
call mpl%flush

end subroutine samp_compute_mpi_f

!----------------------------------------------------------------------
! Subroutine: samp_diag_filter
! Purpose: filter diagnostics
!----------------------------------------------------------------------
subroutine samp_diag_filter(samp,mpl,nam,filter_type,rflt,diag,val)

implicit none

! Passed variables
class(samp_type),intent(in) :: samp                   ! Sampling
type(mpl_type),intent(inout) :: mpl                   ! MPI data
type(nam_type),intent(in) :: nam                      ! Namelist
character(len=*),intent(in) :: filter_type            ! Filter type
real(kind_real),intent(in) :: rflt                    ! Filter support radius
real(kind_real),intent(inout) :: diag(samp%nc2a)      ! Filtered diagnostic
real(kind_real),intent(in),optional :: val(samp%nc2a) ! Useful value for filtering

! Local variables
integer :: ic2a,jc2,nc2eff,kc2,kc2_glb
integer,allocatable :: order(:)
real(kind_real) :: distnorm,norm,wgt
real(kind_real),allocatable :: diag_glb(:),diag_eff(:),diag_eff_dist(:)
real(kind_real),allocatable :: val_glb(:),val_eff(:)
logical :: nam_rad
character(len=1024),parameter :: subr = 'samp_diag_filter'

! Check radius
nam_rad = eq(rflt,nam%diag_rhflt)

if (rflt>0.0) then
   if (nam_rad) then
      ! Allocation
      allocate(diag_glb(samp%nc2f))
      if (present(val)) allocate(val_glb(samp%nc2f))

      ! Communication
      call samp%com_AF%ext(mpl,diag,diag_glb)
      if (present(val)) call samp%com_AF%ext(mpl,diag,diag_glb)
   else
      ! Allocation
      allocate(diag_glb(nam%nc2))
      if (present(val)) allocate(val_glb(nam%nc2))

      ! Local to global
      call mpl%loc_to_glb(samp%nc2a,diag,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.true.,diag_glb)
      if (present(val)) call mpl%loc_to_glb(samp%nc2a,val,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.true.,val_glb)
   end if

   !$omp parallel do schedule(static) private(ic2a,nc2eff,jc2,kc2,kc2_glb,distnorm,norm,wgt), &
   !$omp&                             firstprivate(diag_eff,diag_eff_dist,val_eff,order)
   do ic2a=1,samp%nc2a
      ! Allocation
      allocate(diag_eff(nam%nc2))
      allocate(diag_eff_dist(nam%nc2))
      if (present(val)) allocate(val_eff(nam%nc2))

      ! Build diag_eff of valid points
      nc2eff = 0
      jc2 = 1
      do while (inf(samp%nn_c2a_dist(jc2,ic2a),rflt))
         ! Check the point validity
         kc2 = samp%nn_c2a_index(jc2,ic2a)
         if (nam_rad) then
            kc2_glb = samp%c2_to_c2f(kc2)
         else
            kc2_glb = kc2
         end if
         if (mpl%msv%isnot(diag_glb(kc2_glb))) then
            nc2eff = nc2eff+1
            diag_eff(nc2eff) = diag_glb(kc2_glb)
            diag_eff_dist(nc2eff) = samp%nn_c2a_dist(jc2,ic2a)
            if (present(val)) val_eff(nc2eff) = val_glb(kc2_glb)
         end if
         jc2 = jc2+1
         if (jc2>nam%nc2) exit
      end do

      ! Apply filter
      if (nc2eff>0) then
         select case (trim(filter_type))
         case ('average')
            ! Compute average
            diag(ic2a) = sum(diag_eff(1:nc2eff))/real(nc2eff,kind_real)
         case ('gc99')
            ! Gaspari-Cohn (1999) kernel
            diag(ic2a) = 0.0
            norm = 0.0
            do jc2=1,nc2eff
               distnorm = diag_eff_dist(jc2)/rflt
               wgt = fit_func(mpl,distnorm)
               diag(ic2a) = diag(ic2a)+wgt*diag_eff(jc2)
               norm = norm+wgt
            end do
            if (norm>0.0) diag(ic2a) = diag(ic2a)/norm
         case ('median')
            ! Compute median
            allocate(order(nc2eff))
            if (present(val)) then
               ! Use external value
               call qsort(nc2eff,val_eff(1:nc2eff),order)
               diag_eff = diag_eff(order)
            else
               ! Use diagnostic value
               call qsort(nc2eff,diag_eff(1:nc2eff),order)
            end if
            if (mod(nc2eff,2)==0) then
               diag(ic2a) = 0.5*(diag_eff(nc2eff/2)+diag_eff(nc2eff/2+1))
            else
               diag(ic2a) = diag_eff((nc2eff+1)/2)
            end if
            deallocate(order)
         case default
            ! Wrong filter
            call mpl%abort(subr,'wrong filter type')
         end select
      else
         diag(ic2a) = mpl%msv%valr
      end if

      ! Release memory
      deallocate(diag_eff)
      deallocate(diag_eff_dist)
      if (present(val)) deallocate(val_eff)
   end do
   !$omp end parallel do

   ! Release memory
   deallocate(diag_glb)
   if (present(val)) deallocate(val_glb)
end if

end subroutine samp_diag_filter

!----------------------------------------------------------------------
! Subroutine: samp_diag_fill
! Purpose: fill diagnostics missing values
!----------------------------------------------------------------------
subroutine samp_diag_fill(samp,mpl,nam,diag)

implicit none

! Passed variables
class(samp_type),intent(in) :: samp              ! Sampling
type(mpl_type),intent(inout) :: mpl              ! MPI data
type(nam_type),intent(in) :: nam                 ! Namelist
real(kind_real),intent(inout) :: diag(samp%nc2a) ! Filtered diagnostic

! Local variables
integer :: nmsr,nmsr_tot,ic2a,ic2,jc2,kc2
real(kind_real),allocatable :: diag_glb(:)

! Count missing points
if (samp%nc2a>0) then
   nmsr = count(mpl%msv%is(diag))
else
   nmsr = 0
end if
call mpl%f_comm%allreduce(nmsr,nmsr_tot,fckit_mpi_sum())

if (nmsr_tot>0) then
   ! Allocation
   allocate(diag_glb(nam%nc2))

   ! Local to global
   call mpl%loc_to_glb(samp%nc2a,diag,nam%nc2,samp%c2_to_proc,samp%c2_to_c2a,.true.,diag_glb)

   do ic2a=1,samp%nc2a
      ic2 = samp%c2a_to_c2(ic2a)
      jc2 = 1
      do while (mpl%msv%is(diag(ic2a)))
         kc2 = samp%nn_c2a_index(jc2,ic2a)
         if (mpl%msv%isnot(diag_glb(kc2))) diag(ic2a) = diag_glb(kc2)
         jc2 = jc2+1
         if (jc2>nam%nc2) exit
      end do
   end do

   ! Release memory
   deallocate(diag_glb)
end if

end subroutine samp_diag_fill

end module type_samp
