!----------------------------------------------------------------------
! Module: type_geom
! Purpose: geometry derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_geom

use atlas_module
use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_min,fckit_mpi_max,fckit_mpi_status
#if SABER_USE_METIS
use metis_interface, only: metis_setdefaultoptions,metis_partgraphrecursive,metis_partgraphkway, &
                         & metis_noptions,metis_option_numbering,metis_ok
#endif
use netcdf
use tools_atlas, only: create_atlas_function_space,field_to_fld
use tools_const, only: pi,req,deg2rad,rad2deg,reqkm
use tools_func, only: lonlatmod,sphere_dist,lonlat2xyz,xyz2lonlat,vector_product,vector_triple_product
use tools_kinds, only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use tools_repro, only: inf,eq
use type_com, only: com_type
use type_tree, only: tree_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

! Geometry derived type
type geom_type
   ! Number of points and levels
   integer :: nmg                                 ! Number of model grid points
   integer :: nc0                                 ! Number of points in subset Sc0
   integer :: nl0                                 ! Number of levels in subset Sl0
   integer :: nl0i                                ! Number of independent levels in subset Sl0

   ! Basic geometry data on model grid, halo A
   real(kind_real),allocatable :: lon_mga(:)       ! Longitudes
   real(kind_real),allocatable :: lat_mga(:)       ! Latitudes
   real(kind_real),allocatable :: area_mga(:,:)    ! Domain area
   real(kind_real),allocatable :: vunit_mga(:,:)   ! Vertical unit
   logical,allocatable :: gmask_mga(:,:)           ! Geometry mask
   logical,allocatable :: smask_mga(:,:)           ! Sampling ma

   ! ATLAS function space
   type(atlas_functionspace) :: afunctionspace_mg  ! ATLAS function space of model grid

   ! Basic geometry data
   real(kind_real),allocatable :: lon(:)           ! Longitudes on subset Sc0, global
   real(kind_real),allocatable :: lon_c0a(:)       ! Longitudes on subset Sc0, halo A
   real(kind_real),allocatable :: lat(:)           ! Latitudes on subset Sc0, global
   real(kind_real),allocatable :: lat_c0a(:)       ! Latitudes on subset Sc0, halo A
   real(kind_real),allocatable :: area(:)          ! Domain area
   real(kind_real),allocatable :: vunit_c0(:,:)    ! Vertical unit on subset Sc0, global
   real(kind_real),allocatable :: vunit_c0a(:,:)   ! Vertical unit on subset Sc0, halo A
   real(kind_real),allocatable :: vunitavg(:)      ! Averaged vertical unit
   real(kind_real),allocatable :: disth(:)         ! Horizontal distance

   ! Masks
   logical,allocatable :: mask_c0(:,:)             ! Mask on subset Sc0, global
   logical,allocatable :: mask_c0a(:,:)            ! Mask on subset Sc0, halo A
   logical,allocatable :: mask_hor_c0(:)           ! Union of horizontal masks on subset Sc0, global
   logical,allocatable :: mask_hor_c0a(:)          ! Union of horizontal masks on subset Sc0, halo A
   logical,allocatable :: mask_ver_c0(:)           ! Union of vertical masks
   integer,allocatable :: nc0_mask(:)              ! Horizontal mask size on subset Sc0
   logical,allocatable :: smask_c0a(:,:)           ! Sampling mask on subset Sc0, halo A  
   real(kind_real),allocatable :: mdist(:,:)       ! Minimum distance to mask

   ! Mesh
   type(mesh_type) :: mesh                         ! Mesh

   ! Tree
   type(tree_type) :: tree                         ! Tree

   ! Boundary fields
   integer,allocatable :: nbnda(:)                 ! Number of boundary arcs
   real(kind_real),allocatable :: v1bnda(:,:,:)    ! Boundary arcs, first vector
   real(kind_real),allocatable :: v2bnda(:,:,:)    ! Boundary arcs, second vector
   real(kind_real),allocatable :: vabnda(:,:,:)    ! Boundary arcs, orthogonal vector

   ! Dirac information
   integer :: ndir                                 ! Number of valid Dirac points
   real(kind_real),allocatable :: londir(:)        ! Dirac longitude
   real(kind_real),allocatable :: latdir(:)        ! Dirac latitude
   integer,allocatable :: iprocdir(:)              ! Dirac processor
   integer,allocatable :: ic0adir(:)               ! Dirac gridpoint
   integer,allocatable :: il0dir(:)                ! Dirac level
   integer,allocatable :: ivdir(:)                 ! Dirac variable  
   integer,allocatable :: itsdir(:)                ! Dirac timeslot

   ! MPI distribution
   integer :: nmga                                 ! Halo A size for model grid
   integer :: nc0a                                 ! Halo A size for subset Sc0
   integer,allocatable :: proc_to_nmga(:)          ! Halo A size for each proc
   integer,allocatable :: c0_to_proc_init(:)       ! Subset Sc0 to local task, initial version
   integer,allocatable :: c0_to_proc(:)            ! Subset Sc0 to local task
   logical :: same_grid                            ! Same grid and distribution flag
   integer,allocatable :: c0_to_c0a(:)             ! Subset Sc0, global to halo A
   integer,allocatable :: c0a_to_c0(:)             ! Subset Sc0, halo A to global
   integer,allocatable :: proc_to_nc0a(:)          ! Halo A size for each proc
   integer,allocatable :: c0a_to_mga(:)            ! Subset Sc0 to model grid, halo A
   type(com_type) :: com_mg                        ! Communication between subset Sc0 and model grid
contains
   procedure :: partial_dealloc => geom_partial_dealloc
   procedure :: dealloc => geom_dealloc
   procedure :: setup => geom_setup
   procedure :: from_atlas => geom_from_atlas
   procedure :: remap => geom_remap
   procedure :: define_dirac => geom_define_dirac
   procedure :: check_arc => geom_check_arc
   procedure :: copy_c0a_to_mga => geom_copy_c0a_to_mga
   procedure :: geom_copy_mga_to_c0a_real
   procedure :: geom_copy_mga_to_c0a_logical
   generic :: copy_mga_to_c0a => geom_copy_mga_to_c0a_real,geom_copy_mga_to_c0a_logical
   procedure :: compute_deltas => geom_compute_deltas
end type geom_type

character(len=1024),parameter :: metis_partgraph = 'kway'

private
public :: geom_type

contains

!----------------------------------------------------------------------
! Subroutine: geom_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine geom_partial_dealloc(geom)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry

! Release memory
if (allocated(geom%lon)) deallocate(geom%lon)
if (allocated(geom%lon_c0a)) deallocate(geom%lon_c0a)
if (allocated(geom%lat)) deallocate(geom%lat)
if (allocated(geom%lat_c0a)) deallocate(geom%lat_c0a)
if (allocated(geom%area)) deallocate(geom%area)
if (allocated(geom%vunit_c0)) deallocate(geom%vunit_c0)
if (allocated(geom%vunit_c0a)) deallocate(geom%vunit_c0a)
if (allocated(geom%vunitavg)) deallocate(geom%vunitavg)
if (allocated(geom%disth)) deallocate(geom%disth)
if (allocated(geom%mask_c0)) deallocate(geom%mask_c0)
if (allocated(geom%mask_hor_c0)) deallocate(geom%mask_hor_c0)
if (allocated(geom%mask_hor_c0a)) deallocate(geom%mask_hor_c0a)
if (allocated(geom%mask_ver_c0)) deallocate(geom%mask_ver_c0)
if (allocated(geom%nc0_mask)) deallocate(geom%nc0_mask)
if (allocated(geom%smask_c0a)) deallocate(geom%smask_c0a)
if (allocated(geom%mdist)) deallocate(geom%mdist)
call geom%mesh%dealloc
call geom%tree%dealloc
if (allocated(geom%nbnda)) deallocate(geom%nbnda)
if (allocated(geom%v1bnda)) deallocate(geom%v1bnda)
if (allocated(geom%v2bnda)) deallocate(geom%v2bnda)
if (allocated(geom%vabnda)) deallocate(geom%vabnda)
if (allocated(geom%c0_to_proc_init)) deallocate(geom%c0_to_proc_init)
if (allocated(geom%c0_to_proc)) deallocate(geom%c0_to_proc)
if (allocated(geom%c0_to_c0a)) deallocate(geom%c0_to_c0a)
if (allocated(geom%c0a_to_c0)) deallocate(geom%c0a_to_c0)
if (allocated(geom%proc_to_nc0a)) deallocate(geom%proc_to_nc0a)

end subroutine geom_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: geom_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine geom_dealloc(geom)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry

! Release memory
call geom%partial_dealloc
call geom%afunctionspace_mg%final()
if (allocated(geom%mask_c0a)) deallocate(geom%mask_c0a)
if (allocated(geom%c0a_to_mga)) deallocate(geom%c0a_to_mga)
call geom%com_mg%dealloc

end subroutine geom_dealloc

!----------------------------------------------------------------------
! Subroutine: geom_setup
! Purpose: setup geometry
!----------------------------------------------------------------------
subroutine geom_setup(geom,mpl,rng,nam,afunctionspace,afieldset)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom                 ! Geometry
type(mpl_type),intent(inout) :: mpl                    ! MPI data
type(rng_type),intent(inout) :: rng                    ! Random number generator
type(nam_type),intent(in) :: nam                       ! Namelists
type(atlas_functionspace),intent(in) :: afunctionspace ! ATLAS function space
type(atlas_fieldset),intent(in),optional :: afieldset  ! ATLAS fieldset

! Local variables
integer :: ic0,jc0,kc0,i,j,k,ic0a,jc3,il0,il0i,offset,iproc,img,imga,iend,ibnda,nn_index(1),nc0own,ic0own
integer,allocatable :: proc_to_nmga(:),mg_to_mga(:),mga_to_mg(:),c0_to_mg(:),redundant(:),mg_to_c0(:),mga_to_c0(:),mg_to_proc(:)
integer,allocatable :: c0own_to_mga(:),order(:),order_inv(:),bnda_to_c0(:,:)
real(kind_real) :: lat_arc(2),lon_arc(2),xbnda(2),ybnda(2),zbnda(2)
real(kind_real),allocatable :: sbuf(:),rbuf(:),lon_mg(:),lat_mg(:),area_mg(:),vunit_mg(:,:),list(:)
logical :: same_mask,init,imask,jmask,kmask
logical,allocatable :: sbufl(:),rbufl(:),gmask_mg(:,:),mask_hor_mg(:),not_mask_c0(:)
character(len=1024),parameter :: subr = 'geom_setup'
type(fckit_mpi_status) :: status
type(tree_type) :: tree

! Number of levels
geom%nl0 = nam%nl

! Copy function space pointer
geom%afunctionspace_mg = atlas_functionspace(afunctionspace%c_ptr())
if (present(afieldset)) then
   call geom%from_atlas(mpl,afunctionspace,afieldset)
else
   call geom%from_atlas(mpl,afunctionspace)
end if

! No mask
if (nam%nomask) geom%gmask_mga = .true.

! Allocation
allocate(proc_to_nmga(mpl%nproc))

! Communication
if (mpl%main) then
   ! Receive data on rootproc
   do iproc=1,mpl%nproc
      if (iproc==mpl%rootproc) then
         ! Copy data
         proc_to_nmga(iproc) = geom%nmga
      else
         ! Receive data
         call mpl%f_comm%receive(proc_to_nmga(iproc),iproc-1,mpl%tag,status)
      end if
   end do
else
   ! Send data to rootproc
   call mpl%f_comm%send(geom%nmga,mpl%rootproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Allocation
allocate(sbuf(geom%nmga*(3+geom%nl0)))
allocate(sbufl(geom%nmga*geom%nl0))

! Prepare buffers
sbuf(0*geom%nmga+1:1*geom%nmga) = geom%lon_mga
sbuf(1*geom%nmga+1:2*geom%nmga) = geom%lat_mga
sbuf(2*geom%nmga+1:3*geom%nmga) = geom%area_mga(:,1)
do il0=1,geom%nl0
   sbuf((2+il0)*geom%nmga+1:(3+il0)*geom%nmga) = geom%vunit_mga(:,il0)
   sbufl((il0-1)*geom%nmga+1:il0*geom%nmga) = geom%gmask_mga(:,il0)
end do

! Communication of model grid points
if (mpl%main) then
   ! Global number of model grid points
   geom%nmg = sum(proc_to_nmga)

   ! Allocation
   allocate(lon_mg(geom%nmg))
   allocate(lat_mg(geom%nmg))
   allocate(area_mg(geom%nmg))
   allocate(vunit_mg(geom%nmg,geom%nl0))
   allocate(gmask_mg(geom%nmg,geom%nl0))

   ! Receive data on rootproc
   offset = 0
   do iproc=1,mpl%nproc
      ! Allocation
      allocate(rbuf(proc_to_nmga(iproc)*(3+geom%nl0)))
      allocate(rbufl(proc_to_nmga(iproc)*geom%nl0))

      if (iproc==mpl%rootproc) then
         ! Copy data
         rbuf = sbuf
         rbufl = sbufl
      else
         ! Receive data
         call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
         call mpl%f_comm%receive(rbufl,iproc-1,mpl%tag+1,status)
      end if

      ! Copy buffers
      lon_mg(offset+1:offset+proc_to_nmga(iproc)) = rbuf(0*proc_to_nmga(iproc)+1:1*proc_to_nmga(iproc))
      lat_mg(offset+1:offset+proc_to_nmga(iproc)) = rbuf(1*proc_to_nmga(iproc)+1:2*proc_to_nmga(iproc))
      area_mg(offset+1:offset+proc_to_nmga(iproc)) = rbuf(2*proc_to_nmga(iproc)+1:3*proc_to_nmga(iproc))
      do il0=1,geom%nl0
         vunit_mg(offset+1:offset+proc_to_nmga(iproc),il0) = rbuf((2+il0)*proc_to_nmga(iproc)+1:(3+il0)*proc_to_nmga(iproc))
         gmask_mg(offset+1:offset+proc_to_nmga(iproc),il0) = rbufl((il0-1)*proc_to_nmga(iproc)+1:il0*proc_to_nmga(iproc))
      end do

      ! Update
      offset = offset+proc_to_nmga(iproc)

      ! Release memory
      deallocate(rbuf)
      deallocate(rbufl)
   end do
else
   ! Send data to rootproc
   call mpl%f_comm%send(sbuf,mpl%rootproc-1,mpl%tag)
   call mpl%f_comm%send(sbufl,mpl%rootproc-1,mpl%tag+1)
end if
call mpl%update_tag(2)

if (mpl%main) then
   ! Allocation
   allocate(mg_to_c0(geom%nmg))
   allocate(redundant(geom%nmg))
   allocate(mask_hor_mg(geom%nmg))
   allocate(list(geom%nmg))
   allocate(order(geom%nmg))

   ! Initialization
   redundant = mpl%msv%vali

   ! Look for redundant points
   write(mpl%info,'(a7,a)') '','Look for redundant points in the model grid'
   call mpl%flush

   ! Define points order
   do img=1,geom%nmg
      list(img) = aint(abs(lon_mg(img)+pi)*1.0e6)+abs(lat_mg(img)+0.5*pi)*1.0e-1
   end do
   call qsort(geom%nmg,list,order)

   ! Look for redundant points
   do img=2,geom%nmg
      if (eq(list(img-1),list(img))) redundant(order(img)) = order(img-1)
   end do

   ! Check for successive redundant points
   do img=1,geom%nmg
      if (mpl%msv%isnot(redundant(img))) then
         do while (mpl%msv%isnot(redundant(redundant(img))))
            redundant(img) = redundant(redundant(img))
         end do
      end if
   end do

   ! Horizontal model grid mask
   mask_hor_mg = mpl%msv%is(redundant)

   ! Count subset Sc0 points
   geom%nc0 = count(mask_hor_mg)

   ! Allocation
   allocate(c0_to_mg(geom%nc0))

   ! Initialization
   mg_to_c0 = mpl%msv%vali

   ! Conversion
   ic0 = 0
   do img=1,geom%nmg
      if (mask_hor_mg(img)) then
         ic0 = ic0+1
         c0_to_mg(ic0) = img
         mg_to_c0(img) = ic0
      end if
   end do

   ! Deal with redundant points
   do img=1,geom%nmg
      if (mpl%msv%isnot(redundant(img))) mg_to_c0(img) = mg_to_c0(redundant(img))
   end do

   ! Release memory
   deallocate(list)
   deallocate(order)

   ! Deal with mask on redundant points
   do il0=1,geom%nl0
      do img=1,geom%nmg
         if (mpl%msv%isnot(redundant(img))) gmask_mg(img,il0) = gmask_mg(img,il0).or.gmask_mg(redundant(img),il0)
      end do
   end do

   ! Allocation
   allocate(geom%lon(geom%nc0))
   allocate(geom%lat(geom%nc0))
   allocate(geom%area(geom%nl0))
   allocate(geom%vunit_c0(geom%nc0,geom%nl0))
   allocate(geom%mask_c0(geom%nc0,geom%nl0))

   ! Remove redundant points
   geom%lon = lon_mg(c0_to_mg)
   geom%lat = lat_mg(c0_to_mg)
   do il0=1,geom%nl0
      geom%area(il0) = sum(area_mg(c0_to_mg),gmask_mg(c0_to_mg,il0))/req**2
      geom%vunit_c0(:,il0) = vunit_mg(c0_to_mg,il0)
      geom%mask_c0(:,il0) = gmask_mg(c0_to_mg,il0)
   end do

   ! Release memory
   deallocate(lon_mg)
   deallocate(lat_mg)
   deallocate(area_mg)
   deallocate(vunit_mg)
   deallocate(gmask_mg)
   deallocate(redundant)
   deallocate(mask_hor_mg)
else
   ! Send data to rootproc
   call mpl%f_comm%send(sbuf,mpl%rootproc-1,mpl%tag)
   call mpl%f_comm%send(sbufl,mpl%rootproc-1,mpl%tag+1)
end if
call mpl%update_tag(2)

! Broadcast
call mpl%f_comm%broadcast(geom%nmg,mpl%rootproc-1)
call mpl%f_comm%broadcast(geom%nc0,mpl%rootproc-1)

if (.not.mpl%main) then
   ! Allocation
   allocate(mg_to_c0(geom%nmg))
   allocate(c0_to_mg(geom%nc0))
   allocate(geom%lon(geom%nc0))
   allocate(geom%lat(geom%nc0))
   allocate(geom%area(geom%nl0))
   allocate(geom%vunit_c0(geom%nc0,geom%nl0))
   allocate(geom%mask_c0(geom%nc0,geom%nl0))
end if

! Broadcast
call mpl%f_comm%broadcast(mg_to_c0,mpl%rootproc-1)
call mpl%f_comm%broadcast(c0_to_mg,mpl%rootproc-1)
call mpl%f_comm%broadcast(geom%lon,mpl%rootproc-1)
call mpl%f_comm%broadcast(geom%lat,mpl%rootproc-1)
call mpl%f_comm%broadcast(geom%area,mpl%rootproc-1)
call mpl%f_comm%broadcast(geom%vunit_c0,mpl%rootproc-1)
call mpl%f_comm%broadcast(geom%mask_c0,mpl%rootproc-1)
call mpl%f_comm%broadcast(proc_to_nmga,mpl%rootproc-1)

! Allocation
allocate(mg_to_mga(geom%nmg))
allocate(mga_to_mg(geom%nmga))
allocate(mg_to_proc(geom%nmg))

! Global-local conversion for halo A on model grid
img = 0
do iproc=1,mpl%nproc
   do imga=1,proc_to_nmga(iproc)
      img = img+1
      mg_to_mga(img) = imga
      if (iproc==mpl%myproc) mga_to_mg(imga) = img
      mg_to_proc(img) = iproc
   end do
end do

! Allocation
allocate(geom%mask_hor_c0(geom%nc0))
allocate(geom%mask_ver_c0(geom%nl0))
allocate(geom%nc0_mask(0:geom%nl0))
allocate(geom%vunitavg(geom%nl0))
allocate(geom%c0_to_proc_init(geom%nc0))
allocate(geom%c0_to_proc(geom%nc0))

! Other fields
geom%mask_hor_c0 = any(geom%mask_c0,dim=2)
geom%mask_ver_c0 = any(geom%mask_c0,dim=1)
geom%nc0_mask(0) = count(geom%mask_hor_c0)
geom%nc0_mask(1:geom%nl0) = count(geom%mask_c0,dim=1)
do il0=1,geom%nl0
   if (geom%mask_ver_c0(il0)) then
      geom%vunitavg(il0) = sum(geom%vunit_c0(:,il0),geom%mask_c0(:,il0))/real(geom%nc0_mask(il0),kind_real)
   else
      geom%vunitavg(il0) = 0.0
   end if
end do

! Define halo A on subset Sc0

! Local distribution
do ic0=1,geom%nc0
   img = c0_to_mg(ic0)
   iproc = mg_to_proc(img)
   geom%c0_to_proc_init(ic0) = iproc
end do
geom%c0_to_proc = geom%c0_to_proc_init

if (nam%remap) then
   if (mpl%nproc>1) then
      ! Modify distribution
      write(mpl%info,'(a7,a)') '','Remap points to improve load balance'
      call mpl%flush
#if SABER_USE_METIS
      call geom%remap(mpl,rng,nam)
#else
      call geom%remap(mpl,nam)
#endif
   else
      ! Warning message
      call mpl%warning(subr,'remapping with only one MPI task is not relevant')
   end if
end if

! Check grid similarity
geom%same_grid = (all(geom%c0_to_proc==geom%c0_to_proc_init).and.(geom%nc0==geom%nmg))

! Allocation
allocate(geom%proc_to_nc0a(mpl%nproc))

! Count points on each halo A on subset Sc0
do iproc=1,mpl%nproc
   geom%proc_to_nc0a(iproc) = count(geom%c0_to_proc==iproc)
end do
geom%nc0a = geom%proc_to_nc0a(mpl%myproc)

! Allocation
allocate(geom%c0a_to_c0(geom%nc0a))
allocate(geom%c0_to_c0a(geom%nc0))

! Subset Sc0 conversions
ic0a = 0
do ic0=1,geom%nc0
   iproc = geom%c0_to_proc(ic0)
   if (iproc==mpl%myproc) then
      ic0a = ic0a+1
      geom%c0a_to_c0(ic0a) = ic0
   end if
end do
call mpl%glb_to_loc_index(geom%nc0a,geom%c0a_to_c0,geom%nc0,geom%c0_to_c0a)

! Allocation
allocate(mga_to_c0(geom%nmga))
allocate(geom%c0a_to_mga(geom%nc0a))

! Inter-halo conversions
do imga=1,geom%nmga
   img = mga_to_mg(imga)
   ic0 = mg_to_c0(img)
   mga_to_c0(imga) = ic0
end do
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   img = c0_to_mg(ic0)
   imga = mg_to_mga(img)
   geom%c0a_to_mga(ic0a) = imga
end do

! Release memory
deallocate(proc_to_nmga)
deallocate(mga_to_mg)
deallocate(mg_to_c0)

! Allocation
allocate(order(geom%nc0))
allocate(order_inv(geom%nc0))
allocate(list(geom%nc0))

! Define Sc0 points order
do ic0=1,geom%nc0
   list(ic0) = aint(abs(geom%lon(ic0)+pi)*1.0e6)+abs(geom%lat(ic0)+0.5*pi)*1.0e-1
end do
call qsort(geom%nc0,list,order)
do ic0=1,geom%nc0
   order_inv(order(ic0)) = ic0
end do

! Reorder Sc0 points
geom%c0_to_proc = geom%c0_to_proc(order)
geom%c0_to_c0a = geom%c0_to_c0a(order)
geom%c0a_to_c0 = order_inv(geom%c0a_to_c0)
geom%lon = geom%lon(order)
geom%lat = geom%lat(order)
do il0=1,geom%nl0
   geom%vunit_c0(:,il0) = geom%vunit_c0(order,il0)
   geom%mask_c0(:,il0) = geom%mask_c0(order,il0)
end do
geom%mask_hor_c0 = geom%mask_hor_c0(order)
mga_to_c0 = order_inv(mga_to_c0)
c0_to_mg = c0_to_mg(order)

! Go through Sc0 subset distribution, first pass
nc0own = 0
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   img = c0_to_mg(ic0)
   iproc = mg_to_proc(img)
   if (mpl%myproc==iproc) nc0own = nc0own+1
end do

! Allocation
allocate(c0own_to_mga(nc0own))

! Go through Sc0 subset distribution, second pass
ic0own = 0
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   img = c0_to_mg(ic0)
   iproc = mg_to_proc(img)
   if (mpl%myproc==iproc) then
      ic0own = ic0own+1
      imga = mg_to_mga(img)
      c0own_to_mga(ic0own) = imga
   end if
end do

! Setup redundant points communication
call geom%com_mg%setup(mpl,'com_mg',geom%nc0,geom%nc0a,geom%nmga,nc0own,mga_to_c0,c0own_to_mga,geom%c0_to_proc, &
 & geom%c0_to_c0a)

! Release memory
deallocate(order)
deallocate(order_inv)
deallocate(list)
deallocate(c0_to_mg)
deallocate(mg_to_proc)
deallocate(mg_to_mga)
deallocate(mga_to_c0)
deallocate(c0own_to_mga)

! Allocation
allocate(geom%lon_c0a(geom%nc0a))
allocate(geom%lat_c0a(geom%nc0a))
allocate(geom%vunit_c0a(geom%nc0a,geom%nl0))
allocate(geom%mask_c0a(geom%nc0a,geom%nl0))
allocate(geom%mask_hor_c0a(geom%nc0a))

! Define other fields
geom%lon_c0a = geom%lon(geom%c0a_to_c0)
geom%lat_c0a = geom%lat(geom%c0a_to_c0)
geom%vunit_c0a = geom%vunit_c0(geom%c0a_to_c0,:)
geom%mask_c0a = geom%mask_c0(geom%c0a_to_c0,:)
geom%mask_hor_c0a = geom%mask_hor_c0(geom%c0a_to_c0)

! Allocation
call geom%mesh%alloc(geom%nc0)

! Initialization
call geom%mesh%init(mpl,rng,geom%lon,geom%lat,.true.)

! Compute boundary nodes
call geom%mesh%bnodes(mpl,nam%adv_diag)

! Check whether the mask is the same for all levels
same_mask = .true.
do il0=2,geom%nl0
   same_mask = same_mask.and.(all((geom%mask_c0(:,il0).and.geom%mask_c0(:,1)) &
             & .or.(.not.geom%mask_c0(:,il0).and..not.geom%mask_c0(:,1))))
end do

! Define number of independent levels
if (same_mask) then
   geom%nl0i = 1
else
   geom%nl0i = geom%nl0
end if
write(mpl%info,'(a7,a,i3)') '','Number of independent levels: ',geom%nl0i
call mpl%flush

if ((trim(nam%draw_type)=='random_coast').or.(nam%adv_diag)) then
   ! Define minimum distance to mask
   allocate(geom%mdist(geom%nc0,geom%nl0i))
   geom%mdist = pi
   do il0i=1,geom%nl0i
      ! Check mask
      if (any(.not.geom%mask_c0(:,il0i))) then
         ! Allocation
         allocate(not_mask_c0(geom%nc0))
         not_mask_c0 = .not.geom%mask_c0(:,il0i)
         call tree%alloc(mpl,geom%nc0,mask=not_mask_c0)

         ! Initialization
         call tree%init(geom%lon,geom%lat)

         ! Find nearest neighbors
         do ic0=1,geom%nc0
            if (geom%mask_c0(ic0,il0i)) call tree%find_nearest_neighbors(geom%lon(ic0),geom%lat(ic0),1,nn_index, &
          & geom%mdist(ic0,il0i))
         end do

         ! Release memory
         deallocate(not_mask_c0)
         call tree%dealloc
      end if
   end do
end if

! Allocation
call geom%tree%alloc(mpl,geom%nc0)

! Initialization
call geom%tree%init(geom%lon,geom%lat)

! Horizontal distance
allocate(geom%disth(nam%nc3))
do jc3=1,nam%nc3
   geom%disth(jc3) = real(jc3-1,kind_real)*nam%dc
end do

! Define dirac points
if (nam%new_cortrack.or.nam%new_corstats.or.nam%check_dirac.and.(nam%ndir>0)) call geom%define_dirac(mpl,nam)

if (nam%mask_check) then
   ! Allocation
   allocate(geom%nbnda(0:geom%nl0))

   ! Count boundary arcs
   do il0=0,geom%nl0
      geom%nbnda(il0) = 0
      do i=1,geom%mesh%n
         ic0 = geom%mesh%order(i)
         if (il0==0) then
            imask = geom%mask_hor_c0(ic0)
         else
            imask = geom%mask_c0(ic0,il0)
         end if
         if (.not.imask) then
            iend = geom%mesh%lend(i)
            init = .true.
            do while ((iend/=geom%mesh%lend(i)).or.init)
               j = abs(geom%mesh%list(iend))
               k = abs(geom%mesh%list(geom%mesh%lptr(iend)))
               jc0 = geom%mesh%order(j)
               kc0 = geom%mesh%order(k)
               if (il0==0) then
                   jmask = geom%mask_hor_c0(jc0)
                   kmask = geom%mask_hor_c0(kc0)
               else
                   jmask = geom%mask_c0(jc0,il0)
                   kmask = geom%mask_c0(kc0,il0)
               end if
               if (.not.jmask.and.kmask) geom%nbnda(il0) = geom%nbnda(il0)+1
               iend = geom%mesh%lptr(iend)
               init = .false.
            end do
         end if
      end do
   end do

   ! Allocation
   allocate(geom%v1bnda(3,maxval(geom%nbnda),0:geom%nl0))
   allocate(geom%v2bnda(3,maxval(geom%nbnda),0:geom%nl0))
   allocate(geom%vabnda(3,maxval(geom%nbnda),0:geom%nl0))
   allocate(bnda_to_c0(2,maxval(geom%nbnda)))

   do il0=1,geom%nl0
      ! Define boundary arcs
      ibnda = 0
      do i=1,geom%mesh%n
         ic0 = geom%mesh%order(i)
         if (il0==0) then
            imask = geom%mask_hor_c0(ic0)
         else
            imask = geom%mask_c0(ic0,il0)
         end if
         if (.not.imask) then
            iend = geom%mesh%lend(i)
            init = .true.
            do while ((iend/=geom%mesh%lend(i)).or.init)
               j = abs(geom%mesh%list(iend))
               k = abs(geom%mesh%list(geom%mesh%lptr(iend)))
               jc0 = geom%mesh%order(j)
               kc0 = geom%mesh%order(k)
               if (il0==0) then
                   jmask = geom%mask_hor_c0(jc0)
                   kmask = geom%mask_hor_c0(kc0)
               else
                   jmask = geom%mask_c0(jc0,il0)
                   kmask = geom%mask_c0(kc0,il0)
               end if
               if (.not.jmask.and.kmask) then
                  ibnda = ibnda+1
                  bnda_to_c0(1,ibnda) = ic0
                  bnda_to_c0(2,ibnda) = jc0
               end if
               iend = geom%mesh%lptr(iend)
               init = .false.
            end do
         end if
      end do

      ! Compute boundary arcs coordinates
      do ibnda=1,geom%nbnda(il0)
         lon_arc = geom%lon(bnda_to_c0(:,ibnda))
         lat_arc = geom%lat(bnda_to_c0(:,ibnda))
         call lonlat2xyz(mpl,lon_arc(1),lat_arc(1),xbnda(1),ybnda(1),zbnda(1))
         call lonlat2xyz(mpl,lon_arc(2),lat_arc(2),xbnda(2),ybnda(2),zbnda(2))
         geom%v1bnda(:,ibnda,il0) = (/xbnda(1),ybnda(1),zbnda(1)/)
         geom%v2bnda(:,ibnda,il0) = (/xbnda(2),ybnda(2),zbnda(2)/)
         call vector_product(geom%v1bnda(:,ibnda,il0),geom%v2bnda(:,ibnda,il0),geom%vabnda(:,ibnda,il0))
      end do
   end do
end if

! Print summary
write(mpl%info,'(a7,a,i8)') '','Model grid size:         ',geom%nmg
call mpl%flush
write(mpl%info,'(a7,a,i8)') '','Subset Sc0 size:         ',geom%nc0
call mpl%flush
write(mpl%info,'(a7,a,i6,a,f6.2,a)') '','Number of redundant points:    ',(geom%nmg-geom%nc0), &
 & ' (',real(geom%nmg-geom%nc0,kind_real)/real(geom%nmg,kind_real)*100.0,'%)'
call mpl%flush
write(mpl%info,'(a7,a,f7.1,a,f7.1)') '','Min. / max. longitudes:',minval(geom%lon)*rad2deg,' / ',maxval(geom%lon)*rad2deg
call mpl%flush
write(mpl%info,'(a7,a,f7.1,a,f7.1)') '','Min. / max. latitudes: ',minval(geom%lat)*rad2deg,' / ',maxval(geom%lat)*rad2deg
call mpl%flush
write(mpl%info,'(a7,a,f5.1,a)') '','Domain area (% of Earth area):',100.0*maxval(geom%area)/(4.0*pi),'%'
call mpl%flush
write(mpl%info,'(a7,a)') '','Valid points (% of total domain):'
call mpl%flush
do il0=1,geom%nl0
   write(mpl%info,'(a10,a,i3,a,f5.1,a)') '','Level ',nam%levs(il0),' ~> ', &
 & 100.0*real(geom%nc0_mask(il0),kind_real)/real(geom%nc0,kind_real),'%'
   call mpl%flush
end do
write(mpl%info,'(a7,a)') '','Vertical unit:'
call mpl%flush
do il0=1,geom%nl0
   write(mpl%info,'(a10,a,i3,a,e10.3,a)') '','Level ',nam%levs(il0),' ~> ',geom%vunitavg(il0),' vert. unit'
   call mpl%flush
end do
write(mpl%info,'(a10,a)') '','Distribution summary:'
call mpl%flush
do iproc=1,mpl%nproc
   write(mpl%info,'(a13,a,i4,a,i8,a,i8)') '','Task ',iproc,': ',count(geom%c0_to_proc_init==iproc),' ~> ', &
 & count(geom%c0_to_proc==iproc)
   call mpl%flush
end do

! Deallocate memory
deallocate(geom%lon_mga)
deallocate(geom%lat_mga)
deallocate(geom%area_mga)
deallocate(geom%vunit_mga)
deallocate(geom%gmask_mga)
deallocate(geom%smask_mga)

end subroutine geom_setup

!----------------------------------------------------------------------
! Subroutine: geom_from_atlas
! Purpose: set geometry from ATLAS fieldset
!----------------------------------------------------------------------
subroutine geom_from_atlas(geom,mpl,afunctionspace,afieldset)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom                 ! Geometry
type(mpl_type),intent(inout) :: mpl                    ! MPI data
type(atlas_functionspace),intent(in) :: afunctionspace ! ATLAS function space
type(atlas_fieldset),intent(in),optional :: afieldset  ! ATLAS fieldset

! Local variables
integer :: il0,i,j,imga
real(kind_real) :: lonlat(2)
real(kind_real),pointer :: real_ptr(:,:)
character(len=1024),parameter :: subr = 'geom_from_atlas'
type(atlas_field) :: afield,afield_lonlat
type(atlas_functionspace_nodecolumns) :: afunctionspace_nc
type(atlas_functionspace_structuredcolumns) :: afunctionspace_sc
type(atlas_mesh_nodes) :: anodes
type(atlas_structuredgrid) :: asgrid

select case (afunctionspace%name())
case ('NodeColumns')
    ! Get node columns function space
    afunctionspace_nc = atlas_functionspace_nodecolumns(afunctionspace%c_ptr())

    ! Get number of nodes
    geom%nmga = afunctionspace_nc%nb_nodes()

    ! Allocation
    allocate(geom%lon_mga(geom%nmga))
    allocate(geom%lat_mga(geom%nmga))

    ! Get lon/lat field
    anodes = afunctionspace_nc%nodes()
    afield_lonlat = anodes%lonlat()
    call afield_lonlat%data(real_ptr)
    geom%lon_mga = real_ptr(1,:)*deg2rad
    geom%lat_mga = real_ptr(2,:)*deg2rad
case ('StructuredColumns')
    ! Get structured columns function space
    afunctionspace_sc = atlas_functionspace_structuredcolumns(afunctionspace%c_ptr())

    ! Get number of nodes
    geom%nmga = afunctionspace_sc%size_owned()

    ! Allocation
    allocate(geom%lon_mga(geom%nmga))
    allocate(geom%lat_mga(geom%nmga))

    ! Get lon/lat
    asgrid = afunctionspace_sc%grid()
    imga = 0
    do j=afunctionspace_sc%j_begin(),afunctionspace_sc%j_end()
       do i=afunctionspace_sc%i_begin(j),afunctionspace_sc%i_end(j)
          imga = imga+1
          lonlat = asgrid%lonlat(i,j)
          geom%lon_mga(imga) = lonlat(1)*deg2rad
          geom%lat_mga(imga) = lonlat(2)*deg2rad
       end do
    end do
case default
   call mpl%abort(subr,'wrong function space: '//afunctionspace%name())
end select

! Enforce proper bounds
do imga=1,geom%nmga
   call lonlatmod(geom%lon_mga(imga),geom%lat_mga(imga))
end do

! Allocation
allocate(geom%area_mga(geom%nmga,1))
allocate(geom%vunit_mga(geom%nmga,geom%nl0))
allocate(geom%gmask_mga(geom%nmga,geom%nl0))
allocate(geom%smask_mga(geom%nmga,geom%nl0))

! Default values
geom%area_mga = 1.0_kind_real
do il0=1,geom%nl0
   geom%vunit_mga(:,il0) = real(il0,kind_real)
end do
geom%gmask_mga = .true.
geom%smask_mga = .true.

if (present(afieldset)) then
   ! Get area
   if (afieldset%has_field('area')) then
      afield = afieldset%field('area')
      call field_to_fld(mpl,afield,geom%area_mga)
      call afield%final()
   end if

   ! Get vertical unit
   if (afieldset%has_field('vunit')) then
      afield = afieldset%field('vunit')
      call field_to_fld(mpl,afield,geom%vunit_mga)
      call afield%final()
   end if

   ! Get geometry mask
   if (afieldset%has_field('gmask')) then
      afield = afieldset%field('gmask')
      call field_to_fld(mpl,afield,geom%gmask_mga)
      call afield%final()
   end if

   ! Get sampling mask
   if (afieldset%has_field('smask')) then
      afield = afieldset%field('smask')
      call field_to_fld(mpl,afield,geom%smask_mga)
      call afield%final()
   end if
end if

end subroutine geom_from_atlas

!----------------------------------------------------------------------
! Subroutine: geom_remap
! Purpose: remap points to improve load balance
!----------------------------------------------------------------------
#if SABER_USE_METIS
subroutine geom_remap(geom,mpl,rng,nam)
#else
subroutine geom_remap(geom,mpl,nam)
#endif

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data
#if SABER_USE_METIS
type(rng_type),intent(inout) :: rng    ! Random number generator
#endif
type(nam_type),intent(in) :: nam       ! Namelist

! Local variables
integer :: ncid,nc0_id,lon_c0_id,lat_c0_id,c0_to_proc_init_id,c0_to_proc_id,info
integer :: max_work,ic0,iproc
integer :: c0_to_work(geom%nc0),proc_to_work(mpl%nproc)
#if SABER_USE_METIS
integer :: i,j,jc0,iadj,iend,ierr,im,jm,jproc,imax(2),ipart,ipart_max,iproc_max,nadj,objval,part_to_proc(mpl%nproc)
integer :: affinity(mpl%nproc,mpl%nproc),c0_to_part(geom%nc0)
integer,allocatable :: adjncy(:),adjwgt(:),full_to_m(:),options(:),part(:),vwgt(:),xadj(:)
logical :: init,affinity_mask(mpl%nproc,mpl%nproc)
type(mesh_type) :: mesh
#endif
real(kind_real) :: avg_work
character(len=4) :: nprocchar
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'geom_remap'

! Initialization
do ic0=1,geom%nc0
   c0_to_work(ic0) = count(geom%mask_c0(ic0,:))
end do
avg_work = real(sum(c0_to_work),kind_real)/real(mpl%nproc,kind_real)

! Compute initial load unbalance
do iproc=1,mpl%nproc
   proc_to_work(iproc) = sum(c0_to_work,mask=(geom%c0_to_proc==iproc).and.geom%mask_hor_c0)
end do
max_work = maxval(proc_to_work)
write(mpl%info,'(a10,a,f6.2,a)') '','Initial load unbalance: ',100.0*real(max_work,kind_real)/avg_work,'%'
call mpl%flush

if (mpl%main) then
#if SABER_USE_METIS
   ! Allocation
   call mesh%alloc(geom%nc0)
   allocate(full_to_m(geom%nc0))

   ! Initialization
   call mesh%init(mpl,rng,geom%lon,geom%lat,.false.)
   full_to_m = mpl%msv%vali

   ! Get boundary nodes
   call mesh%bnodes(mpl)

   ! Build adjency list
   nadj = 0
   im = 0
   do i=1,mesh%n
      ic0 = mesh%order(i)
      if (geom%mask_hor_c0(ic0)) then
         im = im+1
         full_to_m(i) = im
         iend = mesh%lend(i)
         init = .true.
         do while ((iend/=mesh%lend(i)).or.init)
            j = mesh%list(iend)
            jc0 = mesh%order(j)
            if ((j>0).and.geom%mask_hor_c0(jc0)) nadj = nadj+1
            iend = mesh%lptr(iend)
            init = .false.
         end do
      end if
   end do

   ! Allocation
   allocate(options(0:metis_noptions-1))
   allocate(xadj(geom%nc0_mask(0)+1))
   allocate(vwgt(geom%nc0_mask(0)))
   allocate(adjncy(nadj))
   allocate(adjwgt(nadj))
   allocate(part(geom%nc0_mask(0)))

   ! Build adjency list
   iadj = 1
   do i=1,mesh%n
      ic0 = mesh%order(i)
      if (geom%mask_hor_c0(ic0)) then
         im = full_to_m(i)
         xadj(im) = iadj
         vwgt(im) = c0_to_work(ic0)
         iproc = geom%c0_to_proc_init(ic0)
         iend = mesh%lend(i)
         init = .true.
         do while ((iend/=mesh%lend(i)).or.init)
            j = mesh%list(iend)
            jc0 = mesh%order(j)
            if ((j>0).and.geom%mask_hor_c0(jc0)) then
               jm = full_to_m(j)
               adjncy(iadj) = jm
               jproc = geom%c0_to_proc_init(jc0)
               if (iproc==jproc) then
                  adjwgt(iadj) = 100
               else
                  adjwgt(iadj) = 1
               end if
               iadj = iadj+1
            end if
            iend = mesh%lptr(iend)
            init = .false.
         end do
      end if
   end do
   xadj(geom%nc0_mask(0)+1) = iadj

   ! Set default METIS options
   ierr = metis_setdefaultoptions(options)
   if (ierr/=metis_ok) call mpl%abort(subr,'metis_setdefaultoptions failed')
   options(metis_option_numbering) = 1

   ! Partition graph using METIS
   if (trim(metis_partgraph)=='recursive') then
      ierr = metis_partgraphrecursive(geom%nc0_mask(0),1,xadj=xadj,adjncy=adjncy,vwgt=vwgt,adjwgt=adjwgt,nparts=mpl%nproc, &
           & options=options,objval=objval,part=part)
   elseif (trim(metis_partgraph)=='kway') then
      ierr = metis_partgraphkway(geom%nc0_mask(0),1,xadj=xadj,adjncy=adjncy,vwgt=vwgt,adjwgt=adjwgt,nparts=mpl%nproc, &
           & options=options,objval=objval,part=part)
   else
      call mpl%abort(subr,'wrong metis_partgraph')
   end if
   if (ierr/=metis_ok) call mpl%abort(subr,'metis_partgraph'//trim(metis_partgraph)//' failed')

   ! Get partitioning
   c0_to_part = mpl%msv%vali
   do i=1,mesh%n
      ic0 = mesh%order(i)
      if (geom%mask_hor_c0(ic0)) then
         im = full_to_m(i)
         c0_to_part(ic0) = part(im)
      end if
   end do

   ! Compute affinity
   affinity = 0
   do ic0=1,geom%nc0
      if (geom%mask_hor_c0(ic0)) then
         iproc = geom%c0_to_proc_init(ic0)
         ipart = c0_to_part(ic0)
         affinity(iproc,ipart) = affinity(iproc,ipart)+1
      end if
   end do

   ! Reassign partitioning
   affinity_mask = .true.
   do iproc=1,mpl%nproc
      imax = maxloc(affinity,mask=affinity_mask)
      iproc_max = imax(1)
      ipart_max = imax(2)
      part_to_proc(ipart_max) = iproc_max
      affinity_mask(iproc_max,:) = .false.
      affinity_mask(:,ipart_max) = .false.
   end do
   do ic0=1,geom%nc0
      if (geom%mask_hor_c0(ic0)) then
         ipart = c0_to_part(ic0)
         iproc = part_to_proc(ipart)
         geom%c0_to_proc(ic0) = iproc
      end if
   end do

   ! Release memory
   deallocate(options)
   deallocate(full_to_m)
   deallocate(xadj)
   deallocate(vwgt)
   deallocate(adjncy)
   deallocate(adjwgt)
   deallocate(part)
#else
   ! Warning message
   call mpl%warning(subr,'METIS has not been compiled, no remapping')
#endif
   ! Define file name
   write(nprocchar,'(i4.4)') mpl%nproc
   filename = trim(nam%prefix)//'_distribution_'//nprocchar

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

   ! Define dimension
   nc0_id = mpl%ncdimcheck(subr,ncid,'nc0',geom%nc0,.true.)

   ! Define variables
   info = nf90_def_var(ncid,'lon_c0',nc_kind_real,(/nc0_id/),lon_c0_id)
   if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_inq_varid(ncid,'lon_c0',lon_c0_id))
   info = nf90_def_var(ncid,'lat_c0',nc_kind_real,(/nc0_id/),lat_c0_id)
   if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_inq_varid(ncid,'lat_c0',lat_c0_id))
   info = nf90_def_var(ncid,'c0_to_proc_init',nc_kind_real,(/nc0_id/),c0_to_proc_init_id)
   if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_inq_varid(ncid,'c0_to_proc_init',c0_to_proc_init_id))
   info = nf90_def_var(ncid,'c0_to_proc',nc_kind_real,(/nc0_id/),c0_to_proc_id)
   if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_inq_varid(ncid,'c0_to_proc',c0_to_proc_id))

   ! End definition mode
   call mpl%ncerr(subr,nf90_enddef(ncid))

   ! Write variables
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c0_id,geom%lon*rad2deg))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c0_id,geom%lat*rad2deg))
   call mpl%ncerr(subr,nf90_put_var(ncid,c0_to_proc_init_id,geom%c0_to_proc_init))
   call mpl%ncerr(subr,nf90_put_var(ncid,c0_to_proc_id,geom%c0_to_proc))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))
end if

! Broadcast results
call mpl%f_comm%broadcast(geom%c0_to_proc,mpl%rootproc-1)

! Print results
do iproc=1,mpl%nproc
   proc_to_work(iproc) = sum(c0_to_work,mask=(geom%c0_to_proc==iproc).and.geom%mask_hor_c0)
end do
max_work = maxval(proc_to_work)
write(mpl%info,'(a10,a,f6.2,a)') '','Final load unbalance: ',100.0*real(maxval(proc_to_work),kind_real)/avg_work,'%'
call mpl%flush
write(mpl%info,'(a10,a,f6.2,a)') '','Final communications: ', &
 & 100.0*real(count(geom%c0_to_proc/=geom%c0_to_proc_init),kind_real)/real(geom%nc0,kind_real),'%'
call mpl%flush

end subroutine geom_remap

!----------------------------------------------------------------------
! Subroutine: geom_define_dirac
! Purpose: define dirac indices
!----------------------------------------------------------------------
subroutine geom_define_dirac(geom,mpl,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist

! Local variables
integer :: idir,il0,nn_index(1),ic0dir,il0dir
character(len=1024),parameter :: subr = 'geom_define_dirac'

! Allocation
allocate(geom%londir(nam%ndir))
allocate(geom%latdir(nam%ndir))
allocate(geom%iprocdir(nam%ndir))
allocate(geom%ic0adir(nam%ndir))
allocate(geom%il0dir(nam%ndir))
allocate(geom%ivdir(nam%ndir))
allocate(geom%itsdir(nam%ndir))

! Initialization
geom%ndir = 0
do idir=1,nam%ndir
   ! Find level
   il0dir = mpl%msv%vali
   do il0=1,geom%nl0
      if (nam%levs(il0)==nam%levdir(idir)) il0dir = il0
   end do
   if (mpl%msv%is(il0dir)) call mpl%abort(subr,'impossible to find the Dirac level')

   ! Find nearest neighbor
   call geom%tree%find_nearest_neighbors(nam%londir(idir),nam%latdir(idir),1,nn_index)
   ic0dir = nn_index(1)

   if (geom%mask_c0(ic0dir,il0dir)) then
      ! Add valid dirac point
      geom%ndir = geom%ndir+1
      geom%londir(geom%ndir) = nam%londir(idir)
      geom%latdir(geom%ndir) = nam%latdir(idir)
      geom%iprocdir(geom%ndir) = geom%c0_to_proc(ic0dir)
      geom%ic0adir(geom%ndir) = geom%c0_to_c0a(ic0dir)
      geom%il0dir(geom%ndir) = il0dir
      geom%ivdir(geom%ndir) = nam%ivdir(idir)
      geom%itsdir(geom%ndir) = nam%itsdir(idir)
   end if
end do

end subroutine geom_define_dirac

!----------------------------------------------------------------------
! Subroutine: geom_check_arc
! Purpose: check if an arc is crossing boundaries
!----------------------------------------------------------------------
subroutine geom_check_arc(geom,mpl,il0,lon_s,lat_s,lon_e,lat_e,valid)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(in) :: il0           ! Level
real(kind_real),intent(in) :: lon_s ! First point longitude
real(kind_real),intent(in) :: lat_s ! First point latitude
real(kind_real),intent(in) :: lon_e ! Second point longitude
real(kind_real),intent(in) :: lat_e ! Second point latitude
logical,intent(out) :: valid        ! True for valid arcs

! Local variables
integer :: ibnda
real(kind_real) :: x(2),y(2),z(2),v1(3),v2(3),va(3),vp(3),t(4)
logical :: cflag(4)

! Initialization
valid = .true.

! Transform to cartesian coordinates
call lonlat2xyz(mpl,lon_s,lat_s,x(1),y(1),z(1))
call lonlat2xyz(mpl,lon_e,lat_e,x(2),y(2),z(2))

! Compute arc orthogonal vector
v1 = (/x(1),y(1),z(1)/)
v2 = (/x(2),y(2),z(2)/)
call vector_product(v1,v2,va)

! Check if arc is crossing boundary arcs
do ibnda=1,geom%nbnda(il0)
   call vector_product(va,geom%vabnda(:,ibnda,il0),vp)
   v1 = (/x(1),y(1),z(1)/)
   call vector_triple_product(v1,va,vp,t(1),cflag(1))
   v1 = (/x(2),y(2),z(2)/)
   call vector_triple_product(v1,va,vp,t(2),cflag(2))
   call vector_triple_product(geom%v1bnda(:,ibnda,il0),geom%vabnda(:,ibnda,il0),vp,t(3),cflag(3))
   call vector_triple_product(geom%v2bnda(:,ibnda,il0),geom%vabnda(:,ibnda,il0),vp,t(4),cflag(4))
   t(1) = -t(1)
   t(3) = -t(3)
   if ((all(t>0.0).or.all(t<0.0)).and.all(cflag)) then
      valid = .false.
      exit
   end if
end do

end subroutine geom_check_arc

!----------------------------------------------------------------------
! Subroutine: geom_copy_c0a_to_mga
! Purpose: copy from subset Sc0 to model grid, halo A
!----------------------------------------------------------------------
subroutine geom_copy_c0a_to_mga(geom,mpl,fld_c0a,fld_mga)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom                        ! Geometry
type(mpl_type),intent(inout) :: mpl                        ! MPI data
real(kind_real),intent(in) :: fld_c0a(geom%nc0a,geom%nl0)  ! Field on subset Sc0, halo A
real(kind_real),intent(out) :: fld_mga(geom%nmga,geom%nl0) ! Field on model grid, halo A

! Local variables
integer :: ic0a,il0
real(kind_real) :: fld_c0a_masked(geom%nc0a,geom%nl0)

! Set masked values at missing value
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      if (geom%mask_c0a(ic0a,il0)) then
         fld_c0a_masked(ic0a,il0) = fld_c0a(ic0a,il0)
      else
         fld_c0a_masked(ic0a,il0) = mpl%msv%valr
      end if
   end do
end do

if (geom%same_grid) then
   ! Same grid
   fld_mga = fld_c0a
else
   ! Extend subset Sc0 to model grid
   call geom%com_mg%ext(mpl,geom%nl0,fld_c0a_masked,fld_mga)
end if

end subroutine geom_copy_c0a_to_mga

!----------------------------------------------------------------------
! Subroutine: geom_copy_mga_to_c0a_real
! Purpose: copy from model grid to subset Sc0, halo A, real
!----------------------------------------------------------------------
subroutine geom_copy_mga_to_c0a_real(geom,mpl,fld_mga,fld_c0a)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom                        ! Geometry
type(mpl_type),intent(inout) :: mpl                        ! MPI data
real(kind_real),intent(in) :: fld_mga(geom%nmga,geom%nl0)  ! Field on model grid, halo A
real(kind_real),intent(out) :: fld_c0a(geom%nc0a,geom%nl0) ! Field on subset Sc0, halo A

! Local variables
integer :: ic0a,il0

if (geom%same_grid) then
   ! Same grid
   fld_c0a = fld_mga
else
   ! Reduce model grid to subset Sc0
   call geom%com_mg%red(mpl,geom%nl0,fld_mga,fld_c0a,.true.)
end if

! Set masked values at missing value
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      if (.not.geom%mask_c0a(ic0a,il0)) fld_c0a(ic0a,il0) = mpl%msv%valr
   end do
end do

end subroutine geom_copy_mga_to_c0a_real

!----------------------------------------------------------------------
! Subroutine: geom_copy_mga_to_c0a_logical
! Purpose: copy from model grid to subset Sc0, halo A, logical
!----------------------------------------------------------------------
subroutine geom_copy_mga_to_c0a_logical(geom,mpl,fld_mga,fld_c0a)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom                ! Geometry
type(mpl_type),intent(inout) :: mpl                ! MPI data
logical,intent(in) :: fld_mga(geom%nmga,geom%nl0)  ! Field on model grid, halo A
logical,intent(out) :: fld_c0a(geom%nc0a,geom%nl0) ! Field on subset Sc0, halo A

! Local variables
integer :: ic0a,imga,il0
real(kind_real) :: fld_mga_real(geom%nmga,geom%nl0),fld_c0a_real(geom%nc0a,geom%nl0)

! Convert array into real values
do il0=1,geom%nl0
   do imga=1,geom%nmga
      if (fld_mga(imga,il0)) then
         fld_mga_real(imga,il0) = 1.0
      else
         fld_mga_real(imga,il0) = 0.0
      end if
   end do
end do

! Copy
call geom%copy_mga_to_c0a(mpl,fld_mga_real,fld_c0a_real)

! Convert array back to logical values
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      fld_c0a(ic0a,il0) = (fld_c0a_real(ic0a,il0)>0.5)
   end do
end do

end subroutine geom_copy_mga_to_c0a_logical

!----------------------------------------------------------------------
! Subroutine: geom_compute_deltas
! Purpose: compute deltas for LCT definition
!----------------------------------------------------------------------
subroutine geom_compute_deltas(geom,ic0,il0,jc0,jl0,dx,dy,dz)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
integer,intent(in) :: ic0           ! First horizontal index
integer,intent(in) :: il0           ! First vertical index
integer,intent(in) :: jc0           ! Second horizontal index
integer,intent(in) :: jl0           ! Second vertical index
real(kind_real),intent(out) :: dx   ! Longitude delta
real(kind_real),intent(out) :: dy   ! Latitude delta
real(kind_real),intent(out) :: dz   ! Altitude delta

! Compute deltas
dx = geom%lon(jc0)-geom%lon(ic0)
dy = geom%lat(jc0)-geom%lat(ic0)
call lonlatmod(dx,dy)
dx = dx*cos(geom%lat(ic0))
dz = real(geom%vunit_c0(ic0,jl0)-geom%vunit_c0(ic0,il0),kind_real)

end subroutine geom_compute_deltas

end module type_geom
