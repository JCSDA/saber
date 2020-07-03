!----------------------------------------------------------------------
! Module: type_geom
! Purpose: geometry derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_geom

use atlas_module, only: atlas_field,atlas_fieldset,atlas_functionspace,atlas_functionspace_nodecolumns, &
 & atlas_functionspace_pointcloud,atlas_functionspace_structuredcolumns,atlas_mesh_nodes,atlas_structuredgrid
use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_min,fckit_mpi_max
use tools_atlas, only: field_to_fld
use tools_const, only: pi,req,deg2rad,rad2deg,reqkm
use tools_func, only: fletcher32,lonlatmod,lonlathash,sphere_dist,lonlat2xyz,xyz2lonlat,vector_product,vector_triple_product
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
   ! Number of processors
   integer :: nproc                               ! Number of processors

   ! ATLAS function space
   type(atlas_functionspace) :: afunctionspace_mg ! ATLAS function space of model grid

   ! Geometry data on model grid, halo A
   integer :: nmga                                ! Halo A size for model grid
   real(kind_real),allocatable :: lon_mga(:)      ! Longitudes
   real(kind_real),allocatable :: lat_mga(:)      ! Latitudes
   real(kind_real),allocatable :: area_mga(:)     ! Areas
   real(kind_real),allocatable :: vunit_mga(:,:)  ! Vertical unit
   logical,allocatable :: gmask_mga(:,:)          ! Geometry mask
   logical,allocatable :: smask_mga(:,:)          ! Sampling mask
   real(kind_real),allocatable :: hash_mga(:)     ! Longitudes/latitudes hash
   logical :: area_provided                       ! Activated if areas are provided

   ! Link between model grid and subset Sc0 on halo A
   type(com_type) :: com_mg                       ! Communication between subset Sc0 and model grid

   ! Geometry data on model grid, global
   integer :: nmg                                 ! Number of model grid points
   integer,allocatable :: proc_to_mg_offset(:)    ! Processor to offset on model grid

   ! Universe
   logical,allocatable :: myuniverse(:)           ! MPI tasks in the universe of the local task

   ! Geometry data on subset Sc0, halo A
   integer,allocatable :: proc_to_nc0a(:)         ! Processor to halo A size for subset Sc0
   integer :: nc0a                                ! Halo A size for subset Sc0
   logical :: same_grid                           ! Same grid and distribution flag
   real(kind_real),allocatable :: lon_c0a(:)      ! Longitudes
   real(kind_real),allocatable :: lat_c0a(:)      ! Latitudes
   real(kind_real),allocatable :: hash_c0a(:)     ! Longitudes/latitudes hash
   real(kind_real),allocatable :: area_c0a(:)     ! Areas
   real(kind_real),allocatable :: vunit_c0a(:,:)  ! Vertical unit
   logical,allocatable :: gmask_c0a(:,:)          ! Geometry mask
   logical,allocatable :: smask_c0a(:,:)          ! Sampling mask
   logical,allocatable :: gmask_hor_c0a(:)        ! Union of horizontal geometry masks
   real(kind_real),allocatable :: mdist_c0a(:,:)  ! Minimum distance to mask
   integer :: grid_hash                           ! Local grid hash

   ! Geometry data on subset Sc0, universe
   integer,allocatable :: proc_to_nc0u(:)         ! Processor to universe size for subset Sc0
   integer :: nc0u                                ! Universe size for subset Sc0
   real(kind_real),allocatable :: lon_c0u(:)      ! Longitudes
   real(kind_real),allocatable :: lat_c0u(:)      ! Latitudes
   real(kind_real),allocatable :: vunit_c0u(:,:)  ! Vertical unit
   logical,allocatable :: gmask_c0u(:,:)          ! Geometry mask
   logical,allocatable :: gmask_hor_c0u(:)        ! Union of horizontal geometry masks
   real(kind_real),allocatable :: mdist_c0u(:,:)  ! Minimum distance to mask

   ! Geometry data on subset Sc0, global
   integer :: nc0                                 ! Number of subset Sc0 points
   integer,allocatable :: proc_to_c0_offset(:)    ! Processor to offset on subset Sc0

   ! Link between halo A and universe on subset Sc0
   integer,allocatable :: c0a_to_c0u(:)           ! Halo A to universe on subset Sc0
   integer,allocatable :: c0u_to_c0a(:)           ! Universe to halo A on subset Sc0
   type(com_type) :: com_AU                       ! Communication between halo A and universe on subset Sc0

   ! Link between halo A and global on subset Sc0
   integer,allocatable :: c0a_to_c0(:)            ! Subset Sc0, halo A to global

   ! Link between universe and global on subset Sc0
   integer,allocatable :: c0u_to_c0(:)            ! Subset Sc0, universe to global

   ! Number of levels
   integer :: nl0                                 ! Number of levels in subset Sl0
   integer :: nl0i                                ! Number of independent levels in subset Sl0

   ! Other fields
   integer,allocatable :: nc0_gmask(:)            ! Horizontal mask size on subset Sc0
   real(kind_real),allocatable :: area(:)         ! Global area
   real(kind_real),allocatable :: vunitavg(:)     ! Averaged vertical unit
   real(kind_real),allocatable :: disth(:)        ! Horizontal distance

   ! Mesh
   type(mesh_type) :: mesh                        ! Mesh

   ! Tree
   type(tree_type) :: tree                        ! Tree

   ! Boundary fields
   integer,allocatable :: nbnda(:)                ! Number of boundary arcs
   real(kind_real),allocatable :: v1bnda(:,:,:)   ! Boundary arcs, first vector
   real(kind_real),allocatable :: v2bnda(:,:,:)   ! Boundary arcs, second vector
   real(kind_real),allocatable :: vabnda(:,:,:)   ! Boundary arcs, orthogonal vector

   ! Dirac fields
   integer :: ndir                                ! Number of valid Dirac points
   real(kind_real),allocatable :: londir(:)       ! Dirac longitude
   real(kind_real),allocatable :: latdir(:)       ! Dirac latitude
   integer,allocatable :: iprocdir(:)             ! Dirac processor
   integer,allocatable :: ic0adir(:)              ! Dirac gridpoint
   integer,allocatable :: il0dir(:)               ! Dirac level
   integer,allocatable :: ivdir(:)                ! Dirac variable
   integer,allocatable :: itsdir(:)               ! Dirac timeslots
contains
   procedure :: partial_dealloc => geom_partial_dealloc
   procedure :: dealloc => geom_dealloc
   procedure :: setup => geom_setup
   procedure :: from_atlas => geom_from_atlas
   procedure :: define_universe => geom_define_universe
   procedure :: setup_c0 => geom_setup_c0
   procedure :: setup_independent_levels => geom_setup_independent_levels
   procedure :: setup_mask_distance => geom_setup_mask_distance
   procedure :: setup_mask_check => geom_setup_mask_check
   procedure :: index_from_lonlat => geom_index_from_lonlat
   procedure :: define_dirac => geom_define_dirac
   procedure :: check_arc => geom_check_arc
   procedure :: copy_c0a_to_mga => geom_copy_c0a_to_mga
   procedure :: geom_copy_mga_to_c0a_real
   procedure :: geom_copy_mga_to_c0a_logical
   generic :: copy_mga_to_c0a => geom_copy_mga_to_c0a_real,geom_copy_mga_to_c0a_logical
   procedure :: compute_deltas => geom_compute_deltas
   procedure :: rand_level => geom_rand_level
   procedure :: rand_point => geom_rand_point
   procedure :: mg_to_proc => geom_mg_to_proc
   procedure :: c0_to_c0a => geom_c0_to_c0a
   procedure :: c0_to_proc => geom_c0_to_proc
   procedure :: c0_to_c0u => geom_c0_to_c0u
end type geom_type

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
if (allocated(geom%lon_mga)) deallocate(geom%lon_mga)
if (allocated(geom%lat_mga)) deallocate(geom%lat_mga)
if (allocated(geom%area_mga)) deallocate(geom%area_mga)
if (allocated(geom%vunit_mga)) deallocate(geom%vunit_mga)
if (allocated(geom%gmask_mga)) deallocate(geom%gmask_mga)
if (allocated(geom%smask_mga)) deallocate(geom%smask_mga)
if (allocated(geom%hash_mga)) deallocate(geom%hash_mga)
if (allocated(geom%proc_to_mg_offset)) deallocate(geom%proc_to_mg_offset)
if (allocated(geom%myuniverse)) deallocate(geom%myuniverse)
if (allocated(geom%proc_to_nc0a)) deallocate(geom%proc_to_nc0a)
if (allocated(geom%lon_c0a)) deallocate(geom%lon_c0a)
if (allocated(geom%lat_c0a)) deallocate(geom%lat_c0a)
if (allocated(geom%hash_c0a)) deallocate(geom%hash_c0a)
if (allocated(geom%area_c0a)) deallocate(geom%area_c0a)
if (allocated(geom%vunit_c0a)) deallocate(geom%vunit_c0a)
if (allocated(geom%smask_c0a)) deallocate(geom%smask_c0a)
if (allocated(geom%gmask_hor_c0a)) deallocate(geom%gmask_hor_c0a)
if (allocated(geom%mdist_c0a)) deallocate(geom%mdist_c0a)
if (allocated(geom%proc_to_nc0u)) deallocate(geom%proc_to_nc0u)
if (allocated(geom%lon_c0u)) deallocate(geom%lon_c0u)
if (allocated(geom%lat_c0u)) deallocate(geom%lat_c0u)
if (allocated(geom%vunit_c0u)) deallocate(geom%vunit_c0u)
if (allocated(geom%gmask_c0u)) deallocate(geom%gmask_c0u)
if (allocated(geom%gmask_hor_c0u)) deallocate(geom%gmask_hor_c0u)
if (allocated(geom%mdist_c0u)) deallocate(geom%mdist_c0u)
if (allocated(geom%proc_to_c0_offset)) deallocate(geom%proc_to_c0_offset)
if (allocated(geom%c0a_to_c0u)) deallocate(geom%c0a_to_c0u)
if (allocated(geom%c0u_to_c0a)) deallocate(geom%c0u_to_c0a)
call geom%com_AU%dealloc
if (allocated(geom%c0a_to_c0)) deallocate(geom%c0a_to_c0)
if (allocated(geom%c0u_to_c0)) deallocate(geom%c0u_to_c0)
if (allocated(geom%nc0_gmask)) deallocate(geom%nc0_gmask)
if (allocated(geom%area)) deallocate(geom%area)
if (allocated(geom%vunitavg)) deallocate(geom%vunitavg)
if (allocated(geom%disth)) deallocate(geom%disth)
call geom%mesh%dealloc
call geom%tree%dealloc
if (allocated(geom%nbnda)) deallocate(geom%nbnda)
if (allocated(geom%v1bnda)) deallocate(geom%v1bnda)
if (allocated(geom%v2bnda)) deallocate(geom%v2bnda)
if (allocated(geom%vabnda)) deallocate(geom%vabnda)
if (allocated(geom%londir)) deallocate(geom%londir)
if (allocated(geom%latdir)) deallocate(geom%latdir)
if (allocated(geom%iprocdir)) deallocate(geom%iprocdir)
if (allocated(geom%ic0adir)) deallocate(geom%ic0adir)
if (allocated(geom%il0dir)) deallocate(geom%il0dir)
if (allocated(geom%ivdir)) deallocate(geom%ivdir)
if (allocated(geom%itsdir)) deallocate(geom%itsdir)

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
if (allocated(geom%gmask_c0a)) deallocate(geom%gmask_c0a)
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
integer :: jc3,il0,iproc
real(kind_real) :: lon_min,lon_max,lat_min,lat_max

! Number of processors
geom%nproc = mpl%nproc

! Number of levels
geom%nl0 = nam%nl

! Copy function space pointer
geom%afunctionspace_mg = atlas_functionspace(afunctionspace%c_ptr())
if (present(afieldset)) then
   call geom%from_atlas(mpl,afunctionspace,afieldset)
else
   call geom%from_atlas(mpl,afunctionspace)
end if

! No mask option
if (nam%nomask) then
   geom%gmask_mga = .true.
   geom%smask_mga = .true.
end if

! Define universe
call geom%define_universe(mpl,rng,nam)
if (nam%default_seed) call rng%reseed(mpl)

! Setup subset Sc0
call geom%setup_c0(mpl)

! Create geometry mesh
write(mpl%info,'(a7,a)') '','Create geometry mesh'
call mpl%flush

! Allocation
call geom%mesh%alloc(geom%nc0u)

! Initialization
call geom%mesh%init(mpl,rng,geom%lon_c0u,geom%lat_c0u,.true.)

! Compute boundary nodes
call geom%mesh%bnodes(mpl,nam%adv_diag)

! Define number of independent levels
call geom%setup_independent_levels(mpl)

! Define minimum distance to mask
if ((trim(nam%draw_type)=='random_coast').or.(nam%adv_diag)) call geom%setup_mask_distance(mpl)

! Create geometry KD-tree
write(mpl%info,'(a7,a)') '','Create geometry KD-tree'
call mpl%flush

! Allocation
call geom%tree%alloc(mpl,geom%nc0u)

! Initialization
call geom%tree%init(geom%lon_c0u,geom%lat_c0u)

! Setup horizontal distance
write(mpl%info,'(a7,a)') '','Setup horizontal distance'
call mpl%flush
allocate(geom%disth(nam%nc3))
do jc3=1,nam%nc3
   geom%disth(jc3) = real(jc3-1,kind_real)*nam%dc
end do

! Define dirac points
if (nam%new_cortrack.or.nam%new_corstats.or.nam%check_dirac.and.(nam%ndir>0)) call geom%define_dirac(mpl,nam)

! Setup mask check
if (nam%mask_check) call geom%setup_mask_check(mpl)

! Summary data
call mpl%f_comm%allreduce(minval(geom%lon_c0a),lon_min,fckit_mpi_min())
call mpl%f_comm%allreduce(maxval(geom%lon_c0a),lon_max,fckit_mpi_max())
call mpl%f_comm%allreduce(minval(geom%lat_c0a),lat_min,fckit_mpi_min())
call mpl%f_comm%allreduce(maxval(geom%lat_c0a),lat_max,fckit_mpi_max())

! Print summary
write(mpl%info,'(a7,a)') '','Geometry summary:'
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','Model grid size:         ',geom%nmg
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','Subset Sc0 size:         ',geom%nc0
call mpl%flush
write(mpl%info,'(a10,a,i6,a,f6.2,a)') '','Number of redundant points:    ',(geom%nmg-geom%nc0), &
 & ' (',real(geom%nmg-geom%nc0,kind_real)/real(geom%nmg,kind_real)*100.0,'%)'
call mpl%flush
write(mpl%info,'(a10,a,f7.1,a,f7.1)') '','Min. / max. longitudes:',lon_min*rad2deg,' / ',lon_max*rad2deg
call mpl%flush
write(mpl%info,'(a10,a,f7.1,a,f7.1)') '','Min. / max. latitudes: ',lat_min*rad2deg,' / ',lat_max*rad2deg
call mpl%flush
write(mpl%info,'(a10,a,f5.1,a)') '','Domain area (% of Earth area):',100.0*maxval(geom%area)/(4.0*pi),'%'
call mpl%flush
write(mpl%info,'(a10,a)') '','Valid points (% of total domain):'
call mpl%flush
do il0=1,geom%nl0
   write(mpl%info,'(a13,a,i3,a,f5.1,a)') '','Level ',nam%levs(il0),' ~> ', &
 & 100.0*real(geom%nc0_gmask(il0),kind_real)/real(geom%nc0,kind_real),'%'
   call mpl%flush
end do
write(mpl%info,'(a10,a)') '','Vertical unit:'
call mpl%flush
do il0=1,geom%nl0
   write(mpl%info,'(a13,a,i3,a,e10.3,a)') '','Level ',nam%levs(il0),' ~> ',geom%vunitavg(il0),' vert. unit'
   call mpl%flush
end do
write(mpl%info,'(a10,a)') '','Distribution (local / universe):'
call mpl%flush
do iproc=1,mpl%nproc
   write(mpl%info,'(a13,a,i6,a,i8,a,i8)') '','Task ',iproc,': ',geom%proc_to_nc0a(iproc),' / ',geom%proc_to_nc0u(iproc)
   call mpl%flush
end do

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
real(kind_real),allocatable :: area_mga(:,:)
real(kind_real),pointer :: real_ptr(:,:)
character(len=1024),parameter :: subr = 'geom_from_atlas'
type(atlas_field) :: afield,afield_lonlat
type(atlas_functionspace_nodecolumns) :: afunctionspace_nc
type(atlas_functionspace_pointcloud) :: afunctionspace_pc
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
case ('PointCloud')
    ! Get point cloud function space
    afunctionspace_pc = atlas_functionspace_pointcloud(afunctionspace%c_ptr())

    ! Get number of points
    geom%nmga = afunctionspace_pc%size()

    ! Allocation
    allocate(geom%lon_mga(geom%nmga))
    allocate(geom%lat_mga(geom%nmga))

    ! Get lon/lat field
    afield_lonlat = afunctionspace_pc%lonlat()
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
allocate(geom%hash_mga(geom%nmga))
allocate(geom%area_mga(geom%nmga))
allocate(geom%vunit_mga(geom%nmga,geom%nl0))
allocate(geom%gmask_mga(geom%nmga,geom%nl0))
allocate(geom%smask_mga(geom%nmga,geom%nl0))

! Compute hash
do imga=1,geom%nmga
   geom%hash_mga(imga) = lonlathash(geom%lon_mga(imga),geom%lat_mga(imga))
end do

! Default values
geom%area_provided = .false.
do il0=1,geom%nl0
   geom%vunit_mga(:,il0) = real(il0,kind_real)
end do
geom%gmask_mga = .true.
geom%smask_mga = .true.

if (present(afieldset)) then
   ! Get area
   if (afieldset%has_field('area')) then
      afield = afieldset%field('area')
      allocate(area_mga(geom%nmga,1))
      call field_to_fld(mpl,afield,area_mga)
      geom%area_mga = area_mga(:,1)/req**2
      deallocate(area_mga)
      call afield%final()
      geom%area_provided = .true.
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
! Subroutine: geom_define_universe
! Purpose: define universe
!----------------------------------------------------------------------
subroutine geom_define_universe(geom,mpl,rng,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(rng_type),intent(inout) :: rng    ! Random number generator
type(nam_type),intent(in) :: nam       ! Namelist

! Local variables
integer :: imga,nnr,inr,iproc,jproc,nb_min,nb_tot,ib,nn_index(1)
integer :: order(geom%nmga),redundant(geom%nmga),proc_to_nb(mpl%nproc),displs(mpl%nproc),proc_to_universe_size(mpl%nproc)
real(kind_real) :: list(geom%nmga),nn_dist(1)
real(kind_real),allocatable :: lon_nr(:),lat_nr(:),lon_b(:),lat_b(:)
logical :: myuniverse
type(mesh_type) :: mesh
type(tree_type) :: tree

write(mpl%info,'(a7,a)') '','Define universe'
call mpl%flush

! Allocation
allocate(geom%myuniverse(mpl%nproc))

if (mpl%nproc<=2) then
   ! Single- or dual-processors case
   geom%myuniverse = .true.
else
   ! Define points order
   list = geom%hash_mga
   call qsort(geom%nmga,list,order)

   ! Look for redundant points
   redundant = mpl%msv%vali
   do imga=2,geom%nmga
      if (eq(list(imga-1),list(imga))) redundant(order(imga)) = order(imga-1)
   end do

   ! Allocation
   nnr = count(mpl%msv%is(redundant))
   allocate(lon_nr(nnr))
   allocate(lat_nr(nnr))
   call mesh%alloc(nnr)
   call tree%alloc(mpl,nnr)

   ! Initialization
   inr = 0
   do imga=1,geom%nmga
      if (mpl%msv%is(redundant(imga))) then
         inr = inr+1
         lon_nr(inr) = geom%lon_mga(imga)
         lat_nr(inr) = geom%lat_mga(imga)
      end if
   end do
   call mesh%init(mpl,rng,lon_nr,lat_nr)
   call tree%init(lon_nr,lat_nr)
   geom%myuniverse = .false.
   geom%myuniverse(mpl%myproc) = .true.

   ! Compute boundary nodes
   call mesh%bnodes(mpl)

   ! Communication
   call mpl%f_comm%allgather(mesh%nb,proc_to_nb)
   nb_min = minval(proc_to_nb)

   if (nb_min>0) then
      ! Allocation
      nb_tot = sum(proc_to_nb)
      allocate(lon_b(nb_tot))
      allocate(lat_b(nb_tot))

      ! Communication
      displs(1) = 0
      do iproc=2,mpl%nproc
         displs(iproc) = displs(iproc-1)+proc_to_nb(iproc-1)
      end do
      call mpl%f_comm%allgather(mesh%lon(mesh%bnd),lon_b,mesh%nb,proc_to_nb,displs)
      call mpl%f_comm%allgather(mesh%lat(mesh%bnd),lat_b,mesh%nb,proc_to_nb,displs)

      ! Check distances
      do iproc=1,mpl%nproc
         ! Check nearest neighbor
         do ib=1,proc_to_nb(iproc)
            if (.not.geom%myuniverse(iproc)) then
               call tree%find_nearest_neighbors(lon_b(displs(iproc)+ib),lat_b(displs(iproc)+ib),1,nn_index,nn_dist)
               if (inf(nn_dist(1),nam%universe_rad)) geom%myuniverse(iproc) = .true.
            end if
         end do
      end do

      ! Comunication
      do iproc=1,mpl%nproc
         if (iproc==mpl%myproc) then
            do jproc=1,mpl%nproc
               ! Receive data
               if (jproc/=iproc) then
                  call mpl%f_comm%receive(myuniverse,jproc-1,mpl%tag)
                  geom%myuniverse(jproc) = geom%myuniverse(jproc).or.myuniverse
               end if
            end do
         else
            ! Send data
            call mpl%f_comm%send(geom%myuniverse(iproc),iproc-1,mpl%tag)
         end if
         call mpl%update_tag(1)
      end do

      ! Release memory
      call mesh%dealloc
      call tree%dealloc
      deallocate(lon_b)
      deallocate(lat_b) 
   else
      ! Impossible to disentangle the local domains
      geom%myuniverse = .true.
   endif
end if

! Share universe size
call mpl%f_comm%allgather(count(geom%myuniverse),proc_to_universe_size)

! Print results
write(mpl%info,'(a7,a)') '','Size of universes:'
call mpl%flush
do iproc=1,mpl%nproc
   write(mpl%info,'(a10,a,i6,a,f6.2,a)') '','Task ',iproc,': ',100.0*real(proc_to_universe_size(iproc),kind_real) &
    & /real(mpl%nproc,kind_real),'%'
   call mpl%flush
end do

end subroutine geom_define_universe

!----------------------------------------------------------------------
! Subroutine: geom_setup_c0
! Purpose: setup subset Sc0
!----------------------------------------------------------------------
subroutine geom_setup_c0(geom,mpl)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data

! Local variables
integer :: nmgu,iproc,img,jmg,imga,imgu,imgu_s,imgu_e,imgu_min,ic0a,jc0a,ic0u,ic0,diff_grid,diff_grid_tot,il0,nr,nra,ir,ira
integer :: proc_to_nmga(mpl%nproc),proc_to_nra(mpl%nproc)
integer :: mga_to_mg(geom%nmga),mga_to_mgu(geom%nmga),mga_to_c0(geom%nmga),nc0_gmask(0:geom%nl0)
integer,allocatable :: mgu_to_mg(:),mgu_to_proc(:),redundant(:),order(:)
integer,allocatable :: c0a_to_mg(:),ra_to_r(:),r_to_mg(:),r_to_mg_tot(:),r_to_c0(:),r_to_c0_tot(:)
real(kind_real) :: areasum(geom%nl0),vunitsum(geom%nl0),vunitsum_tot(geom%nl0),norm(geom%nl0),norm_tot(geom%nl0)
real(kind_real),allocatable :: hash_mgu(:),lonlat_c0a(:),hash_c0u(:)
logical :: diff,chain
character(len=1024),parameter :: subr = 'geom_setup_c0'
type(com_type) :: com_AU

write(mpl%info,'(a7,a)') '','Setup subset Sc0'
call mpl%flush

! Allocation
allocate(geom%proc_to_nc0a(mpl%nproc))

! Communication
call mpl%f_comm%allgather(geom%nmga,proc_to_nmga)

! Model grid sizes
geom%nmg = sum(proc_to_nmga)
nmgu = sum(proc_to_nmga,mask=geom%myuniverse)

! Allocation
allocate(geom%proc_to_mg_offset(mpl%nproc))
allocate(mgu_to_mg(nmgu))
allocate(hash_mgu(nmgu))
allocate(order(nmgu))
allocate(redundant(nmgu))

! Model grid offset for halo A
geom%proc_to_mg_offset(1) = 0
do iproc=2,mpl%nproc
   geom%proc_to_mg_offset(iproc) = geom%proc_to_mg_offset(iproc-1)+proc_to_nmga(iproc-1)
end do

! Model grid conversions
img = 0
imgu = 0
do iproc=1,mpl%nproc
   do imga=1,proc_to_nmga(iproc)
      img = img+1
      if (iproc==mpl%myproc) mga_to_mg(imga) = img
      if (geom%myuniverse(iproc)) then
         imgu = imgu+1
         if (iproc==mpl%myproc) mga_to_mgu(imga) = imgu
         mgu_to_mg(imgu) = img
      end if
   end do
end do

! Setup model grid communication, local to universe
call com_AU%setup(mpl,'com_AU',geom%nmga,nmgu,geom%nmg,mga_to_mg,mgu_to_mg)

! Extend hash on model grid, halo A to universe
call com_AU%ext(mpl,geom%hash_mga,hash_mgu)

! Look for redundant points
write(mpl%info,'(a10,a)') '','Look for redundant points in the model grid'
call mpl%flush

! Define points order
call qsort(nmgu,hash_mgu,order)

! Look for redundant points
redundant = mpl%msv%vali
imgu_s = 1
chain = .false.
do imgu_e=2,nmgu+1
   if (imgu_e==nmgu+1) then
      diff = .true.
   else
      diff = inf(hash_mgu(imgu_s),hash_mgu(imgu_e))
   end if
   if (diff) then
      ! Different hash value
      if (chain) then
         ! Allocation
         allocate(mgu_to_proc(imgu_s:imgu_e-1))

         ! Get the processor index
         do imgu=imgu_s,imgu_e-1
            img = mgu_to_mg(order(imgu))
            iproc = geom%mg_to_proc(img)
            mgu_to_proc(imgu) = iproc
         end do

         ! Find the the smaller processor index with the smaller model grid index
         imgu_min = imgu_s
         do imgu=imgu_s+1,imgu_e-1
            if (mgu_to_proc(imgu)<=mgu_to_proc(imgu_min)) then
               if (mgu_to_mg(order(imgu))<mgu_to_mg(order(imgu_min))) imgu_min = imgu
            end if
         end do

         ! All the others are considered redundant
         do imgu=imgu_s,imgu_e-1
            if (imgu/=imgu_min) redundant(order(imgu)) = mgu_to_mg(order(imgu_min))
         end do

         ! Chain done
         chain = .false.

         ! Release memory
         deallocate(mgu_to_proc)
      end if

      ! Update
      imgu_s = imgu_e
   else
      ! Same hash value, this is a chain
      chain = .true.
   endif
end do

! Define subset Sc0
write(mpl%info,'(a10,a)') '','Define subset Sc0'
call mpl%flush

! Count points for subset Sc0, halo A
geom%nc0a = 0
do imga=1,geom%nmga
   imgu = mga_to_mgu(imga)
   if (mpl%msv%is(redundant(imgu))) geom%nc0a = geom%nc0a+1
end do

! Check grid similarity
if (geom%nc0a==geom%nmga) then
   diff_grid = 0
else
   diff_grid = 1
end if
call mpl%f_comm%allreduce(diff_grid,diff_grid_tot,fckit_mpi_sum())
geom%same_grid = (diff_grid_tot==0)

! Communication
call mpl%f_comm%allgather(geom%nc0a,geom%proc_to_nc0a)

! Subset Sc0 global size
geom%nc0 = sum(geom%proc_to_nc0a)

! Allocation
allocate(geom%proc_to_c0_offset(mpl%nproc))
allocate(geom%c0a_to_c0(geom%nc0a))
allocate(c0a_to_mg(geom%nc0a))
allocate(geom%lon_c0a(geom%nc0a))
allocate(geom%lat_c0a(geom%nc0a))
allocate(geom%hash_c0a(geom%nc0a))
if (geom%area_provided) allocate(geom%area_c0a(geom%nc0a))
allocate(geom%vunit_c0a(geom%nc0a,geom%nl0))
allocate(geom%gmask_c0a(geom%nc0a,geom%nl0))
allocate(geom%smask_c0a(geom%nc0a,geom%nl0))
allocate(geom%gmask_hor_c0a(geom%nc0a))
allocate(geom%nc0_gmask(0:geom%nl0))
allocate(geom%area(geom%nl0))
allocate(geom%vunitavg(geom%nl0))
allocate(lonlat_c0a(2*geom%nc0a))

! Subset Sc0 offset for halo A
geom%proc_to_c0_offset(1) = 0
do iproc=2,mpl%nproc
   geom%proc_to_c0_offset(iproc) = geom%proc_to_c0_offset(iproc-1)+geom%proc_to_nc0a(iproc-1)
end do

! Communicate data from model grid to subset Sc0
write(mpl%info,'(a10,a)') '','Communicate data from model grid to subset Sc0'
call mpl%flush

! Conversions
ic0a = 0
nra = 0
do imga=1,geom%nmga
   imgu = mga_to_mgu(imga)
   if (mpl%msv%is(redundant(imgu))) then
      ! Valid Sc0 point
      ic0a = ic0a+1
      ic0 = geom%proc_to_c0_offset(mpl%myproc)+ic0a
      geom%c0a_to_c0(ic0a) = ic0
      img = geom%proc_to_mg_offset(mpl%myproc)+imga
      c0a_to_mg(ic0a) = img
   else
      ! Redundant point
      nra = nra+1
   end if
end do

! Communication
call mpl%f_comm%allgather(nra,proc_to_nra)

! Global number of redundant points
nr = sum(proc_to_nra)

! Allocation
allocate(ra_to_r(nra))
allocate(r_to_mg(nr))
allocate(r_to_mg_tot(nr))
allocate(r_to_c0(nr))
allocate(r_to_c0_tot(nr))

! Conversion
ir = 0
do iproc=1,mpl%nproc
   do ira=1,proc_to_nra(iproc)
      ir = ir+1
      if (iproc==mpl%myproc) ra_to_r(ira) = ir
   end do
end do

! Copy model grid index of redundant points
r_to_mg = 0
ira = 0
do imga=1,geom%nmga
   imgu = mga_to_mgu(imga)
   if (mpl%msv%isnot(redundant(imgu))) then
      ira = ira+1
      ir = ra_to_r(ira)
      r_to_mg(ir) = redundant(imgu)
   end if
end do

! Communication
call mpl%f_comm%allreduce(r_to_mg,r_to_mg_tot,fckit_mpi_sum())

! Find Sc0 point for each redundant point
r_to_c0 = 0
ic0a = 0
do ir=1,nr
   img = r_to_mg_tot(ir)
   iproc = geom%mg_to_proc(img)
   if (iproc==mpl%myproc) then
      jc0a = 1
      do while (jc0a<=geom%nc0a)
         jmg = c0a_to_mg(jc0a)
         if (img==jmg) then
            ic0 = geom%c0a_to_c0(jc0a)
            r_to_c0(ir) = ic0
         end if
         jc0a = jc0a+1
      end do
   end if
end do

! Communication
call mpl%f_comm%allreduce(r_to_c0,r_to_c0_tot,fckit_mpi_sum())

! Conversion
ic0a = 0
ira = 0
do imga=1,geom%nmga
   imgu = mga_to_mgu(imga)
   if (mpl%msv%is(redundant(imgu))) then
      ! Valid Sc0 point
      ic0a = ic0a+1
      ic0 = geom%c0a_to_c0(ic0a)
   else
      ! Redundant point
      ira = ira+1
      ir = ra_to_r(ira)
      ic0 = r_to_c0_tot(ir)
   end if
   mga_to_c0(imga) = ic0
end do

! Setup redundant points communication
call geom%com_mg%setup(mpl,'com_mg',geom%nc0a,geom%nmga,geom%nc0,geom%c0a_to_c0,mga_to_c0)

! Reduce fields from model grid to subset Sc0 on halo A
call geom%com_mg%red(mpl,geom%lon_mga,geom%lon_c0a,.true.)
call geom%com_mg%red(mpl,geom%lat_mga,geom%lat_c0a,.true.)
call geom%com_mg%red(mpl,geom%hash_mga,geom%hash_c0a,.true.)
if (geom%area_provided) call geom%com_mg%red(mpl,geom%area_mga,geom%area_c0a,.true.)
call geom%com_mg%red(mpl,geom%nl0,geom%vunit_mga,geom%vunit_c0a,.true.)
call geom%com_mg%red(mpl,geom%nl0,geom%gmask_mga,geom%gmask_c0a,.true.)
call geom%com_mg%red(mpl,geom%nl0,geom%smask_mga,geom%smask_c0a,.true.)

! Related fields
geom%gmask_hor_c0a = any(geom%gmask_c0a,dim=2)
nc0_gmask(0) = count(geom%gmask_hor_c0a)
nc0_gmask(1:geom%nl0) = count(geom%gmask_c0a,dim=1)
call mpl%f_comm%allreduce(nc0_gmask,geom%nc0_gmask,fckit_mpi_sum())
do il0=1,geom%nl0
   norm(il0) = count(geom%gmask_c0a(:,il0))
   if (norm(il0)>0) then
      if (geom%area_provided) areasum(il0) = sum(geom%area_c0a,mask=geom%gmask_c0a(:,il0))
      vunitsum(il0) = sum(geom%vunit_c0a(:,il0),mask=geom%gmask_c0a(:,il0))
   else
      vunitsum(il0) = 0.0
   end if
end do
call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
if (geom%area_provided) then
   call mpl%f_comm%allreduce(areasum,geom%area,fckit_mpi_sum())
else
   geom%area = 4.0*pi
end if
call mpl%f_comm%allreduce(vunitsum,vunitsum_tot,fckit_mpi_sum())
do il0=1,geom%nl0
   if (norm_tot(il0)>0) then
      geom%vunitavg(il0) = vunitsum_tot(il0)/real(norm_tot(il0),kind_real)
   else
      geom%vunitavg(il0) = mpl%msv%valr
   end if
end do

! Grid hash
lonlat_c0a(1:geom%nc0a) = geom%lon_c0a
lonlat_c0a(geom%nc0a+1:2*geom%nc0a) = geom%lat_c0a
geom%grid_hash = fletcher32(lonlat_c0a)

! Deallocate memory
deallocate(geom%lon_mga)
deallocate(geom%lat_mga)
deallocate(geom%area_mga)
deallocate(geom%vunit_mga)
deallocate(geom%gmask_mga)
deallocate(geom%smask_mga)
deallocate(mgu_to_mg)
deallocate(hash_mgu)
deallocate(order)
deallocate(redundant)
deallocate(lonlat_c0a)
deallocate(c0a_to_mg)
deallocate(ra_to_r)
deallocate(r_to_mg)
deallocate(r_to_mg_tot)
deallocate(r_to_c0)
deallocate(r_to_c0_tot)

! Communicate data from halo A to universe on subset Sc0
write(mpl%info,'(a10,a)') '','Communicate data from halo A to universe on subset Sc0'
call mpl%flush

! Subset Sc0 universe size
geom%nc0u = sum(geom%proc_to_nc0a,mask=geom%myuniverse)

! Allocation
allocate(geom%proc_to_nc0u(mpl%nproc))

! Communication
call mpl%f_comm%allgather(geom%nc0u,geom%proc_to_nc0u)

! Check universe size
if (mpl%nproc>1) then
   do iproc=1,mpl%nproc
      if (geom%proc_to_nc0u(iproc)==geom%proc_to_nc0a(iproc)) call mpl%abort(subr,'universe is not larger than halo A')
   end do
end if

! Allocation
allocate(geom%c0a_to_c0u(geom%nc0a))
allocate(geom%c0u_to_c0a(geom%nc0u))
allocate(geom%c0u_to_c0(geom%nc0u))
allocate(geom%lon_c0u(geom%nc0u))
allocate(geom%lat_c0u(geom%nc0u))
allocate(hash_c0u(geom%nc0u))
allocate(geom%vunit_c0u(geom%nc0u,geom%nl0))
allocate(geom%gmask_c0u(geom%nc0u,geom%nl0))
allocate(geom%gmask_hor_c0u(geom%nc0u))
allocate(order(geom%nc0u))

! Conversions
geom%c0u_to_c0a = mpl%msv%vali
ic0u = 0
do ic0=1,geom%nc0
   iproc = geom%c0_to_proc(ic0)
   if (geom%myuniverse(iproc)) then
      ic0u = ic0u+1
      if (iproc==mpl%myproc) then
         ic0a = geom%c0_to_c0a(ic0)
         geom%c0a_to_c0u(ic0a) = ic0u
         geom%c0u_to_c0a(ic0u) = ic0a
      end if
      geom%c0u_to_c0(ic0u) = ic0
   end if
end do

! Setup subset Sc0 communication, local to universe
call geom%com_AU%setup(mpl,'com_AU',geom%nc0a,geom%nc0u,geom%nc0,geom%c0a_to_c0,geom%c0u_to_c0)

! Extend fields from halo A to universe on subset Sc0
call geom%com_AU%ext(mpl,geom%lon_c0a,geom%lon_c0u)
call geom%com_AU%ext(mpl,geom%lat_c0a,geom%lat_c0u)
call geom%com_AU%ext(mpl,geom%hash_c0a,hash_c0u)
call geom%com_AU%ext(mpl,geom%nl0,geom%vunit_c0a,geom%vunit_c0u)
call geom%com_AU%ext(mpl,geom%nl0,geom%gmask_c0a,geom%gmask_c0u)
geom%gmask_hor_c0u = any(geom%gmask_c0u,dim=2)


! Check that Sc0 points in universe are not redundant
write(mpl%info,'(a7,a)') '','Check that Sc0 points in universe are not redundant'
call mpl%flush
call qsort(geom%nc0u,hash_c0u,order)
do ic0u=2,geom%nc0u
   if (eq(hash_c0u(ic0u),hash_c0u(ic0u-1))) call mpl%abort(subr,'redundant points in Sc0 point on universe, check the universe')
end do

! Release memory
deallocate(order)
deallocate(hash_c0u)

end subroutine geom_setup_c0

!----------------------------------------------------------------------
! Subroutine: geom_setup_independent_levels
! Purpose: setup independent levels
!----------------------------------------------------------------------
subroutine geom_setup_independent_levels(geom,mpl)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data

! Local variables
integer :: il0
integer :: diff_mask,diff_mask_tot
logical :: same_mask

write(mpl%info,'(a7,a)') '','Define number of independent levels'
call mpl%flush

! Check local mask similarity
same_mask = .true.
do il0=2,geom%nl0
   same_mask = same_mask.and.(all((geom%gmask_c0a(:,il0).and.geom%gmask_c0a(:,1)) &
             & .or.(.not.geom%gmask_c0a(:,il0).and..not.geom%gmask_c0a(:,1))))
end do
if (same_mask) then
   diff_mask = 0
else
   diff_mask = 1
end if

! Communication
call mpl%f_comm%allreduce(diff_mask,diff_mask_tot,fckit_mpi_sum())

! Global mask similarity
same_mask = (diff_mask_tot==0)

! Define number of independent levels
if (same_mask) then
   geom%nl0i = 1
else
   geom%nl0i = geom%nl0
end if
write(mpl%info,'(a7,a,i3)') '','Number of independent levels: ',geom%nl0i
call mpl%flush

end subroutine geom_setup_independent_levels

!----------------------------------------------------------------------
! Subroutine: geom_setup_mask_distance
! Purpose: setup minimum distance to mask
!----------------------------------------------------------------------
subroutine geom_setup_mask_distance(geom,mpl)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data

! Local variables
integer :: il0i,ic0u,ic0a,nn_index(1)
logical :: not_mask_c0u(geom%nc0u)
type(tree_type) :: tree

write(mpl%info,'(a7,a)') '','Define minimum distance to mask'
call mpl%flush

! Allocation
allocate(geom%mdist_c0u(geom%nc0u,geom%nl0i))
allocate(geom%mdist_c0a(geom%nc0a,geom%nl0i))

! Initialization
geom%mdist_c0u = pi

do il0i=1,geom%nl0i
   ! Check mask
   if (any(.not.geom%gmask_c0u(:,il0i))) then
      ! Allocation
      not_mask_c0u = .not.geom%gmask_c0u(:,il0i)
      call tree%alloc(mpl,geom%nc0u,mask=not_mask_c0u)

      ! Initialization
      call tree%init(geom%lon_c0u,geom%lat_c0u)

      ! Find nearest neighbors
      do ic0u=1,geom%nc0u
         if (geom%gmask_c0u(ic0u,il0i)) call tree%find_nearest_neighbors(geom%lon_c0u(ic0u),geom%lat_c0u(ic0u),1,nn_index, &
       & geom%mdist_c0u(ic0u,il0i))
      end do

      ! Release memory
      call tree%dealloc
   end if
end do

! Local field
do ic0a=1,geom%nc0a
   ic0u = geom%c0a_to_c0u(ic0a)
   geom%mdist_c0a(ic0a,:) = geom%mdist_c0u(ic0u,:)
end do

end subroutine geom_setup_mask_distance

!----------------------------------------------------------------------
! Subroutine: geom_setup_mask_check
! Purpose: setup mask checking tool
!----------------------------------------------------------------------
subroutine geom_setup_mask_check(geom,mpl)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data

! Local variables
integer :: il0,ibnda,nbndamax
integer,allocatable :: bnda_to_c0u(:,:)
real(kind_real) :: lon_arc(2),lat_arc(2),xbnda(2),ybnda(2),zbnda(2)

write(mpl%info,'(a7,a)') '','Setup mask check'
call mpl%flush

! Allocation
allocate(geom%nbnda(0:geom%nl0))

! Count boundary arcs
do il0=0,geom%nl0
   if (il0==0) then
      call geom%mesh%count_bnda(geom%gmask_hor_c0u,geom%nbnda(il0))
   else
      call geom%mesh%count_bnda(geom%gmask_c0u(:,il0),geom%nbnda(il0))
   end if
end do

! Allocation
nbndamax = maxval(geom%nbnda)
allocate(bnda_to_c0u(2,nbndamax))
allocate(geom%v1bnda(3,nbndamax,0:geom%nl0))
allocate(geom%v2bnda(3,nbndamax,0:geom%nl0))
allocate(geom%vabnda(3,nbndamax,0:geom%nl0))

do il0=1,geom%nl0
   ! Get boundary arcs
   if (il0==0) then
      call geom%mesh%get_bnda(geom%gmask_hor_c0u,geom%nbnda(il0),bnda_to_c0u)
   else
      call geom%mesh%get_bnda(geom%gmask_c0u(:,il0),geom%nbnda(il0),bnda_to_c0u)
   end if

   ! Compute boundary arcs coordinates
   do ibnda=1,geom%nbnda(il0)
      lon_arc = geom%lon_c0u(bnda_to_c0u(:,ibnda))
      lat_arc = geom%lat_c0u(bnda_to_c0u(:,ibnda))
      call lonlat2xyz(mpl,lon_arc(1),lat_arc(1),xbnda(1),ybnda(1),zbnda(1))
      call lonlat2xyz(mpl,lon_arc(2),lat_arc(2),xbnda(2),ybnda(2),zbnda(2))
      geom%v1bnda(:,ibnda,il0) = (/xbnda(1),ybnda(1),zbnda(1)/)
      geom%v2bnda(:,ibnda,il0) = (/xbnda(2),ybnda(2),zbnda(2)/)
      call vector_product(geom%v1bnda(:,ibnda,il0),geom%v2bnda(:,ibnda,il0),geom%vabnda(:,ibnda,il0))
   end do
end do

end subroutine geom_setup_mask_check

!----------------------------------------------------------------------
! Subroutine: geom_index_from_lonlat
! Purpose: get nearest neighbor index from longitude/latitude/level
!----------------------------------------------------------------------
subroutine geom_index_from_lonlat(geom,mpl,lon,lat,il0,iproc,ic0a,gmask)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl ! MPI data
real(kind_real),intent(in) :: lon   ! Longitude
real(kind_real),intent(in) :: lat   ! Latitude
integer,intent(in) :: il0           ! Level index
integer,intent(out) :: iproc        ! Task index
integer,intent(out) :: ic0a         ! Local index
logical,intent(out) :: gmask        ! Local mask

! Local variables
integer :: nn_index(1),ic0,nn_proc,proc_to_nn_proc(mpl%nproc),jproc,ic0u
real(kind_real) :: nn_dist(1),proc_to_nn_dist(mpl%nproc),distmin
logical :: valid
character(len=1024),parameter :: subr = 'geom_index_from_lonlat'

integer :: istart,iend,proc_to_time(mpl%nproc),iloc(1)

! Find nearest neighbor

! TODO: finaud! créer un mesh sur la halo A, extraire les point de bordure, créer un petit arbre avec cette bordure, et vérifier que la distance mini entre cette bordure et le point lon/lat est inférieure à la distance max entre deux points de la bordure.
call mpl%f_comm%barrier()
call system_clock(count=istart)
call geom%tree%find_nearest_neighbors(lon,lat,1,nn_index,nn_dist)
call system_clock(count=iend)
call mpl%f_comm%allgather(iend-istart,proc_to_time)
iloc = maxloc(proc_to_time)
if ((mpl%myproc==iloc(1)).and.(maxval(proc_to_time)>8)) print*, 'bad',iend-istart,nn_dist*reqkm

ic0 = geom%c0u_to_c0(nn_index(1))
nn_proc = geom%c0_to_proc(ic0)

! Communication
call mpl%f_comm%allgather(nn_dist(1),proc_to_nn_dist)
call mpl%f_comm%allgather(nn_proc,proc_to_nn_proc)

! The correct processor should handle its own nearest neighbor, with the minimum distance
distmin = minval(proc_to_nn_dist)
iproc = mpl%msv%vali
do jproc=1,mpl%nproc
   if ((proc_to_nn_proc(jproc)==jproc).and.eq(proc_to_nn_dist(jproc),distmin)) iproc = jproc
end do
if (mpl%msv%is(iproc)) call mpl%abort(subr,'cannot find root processor')

if (iproc==mpl%myproc) then
   ! Check whether the location is in the convex hull      
   call geom%mesh%inside(mpl,lon,lat,valid)

   ! Local index
   ic0u = nn_index(1)
   ic0a = geom%c0u_to_c0a(ic0u)

   ! Check mask
   if (il0==0) then
      gmask = geom%gmask_hor_c0a(ic0a)
   else
      gmask = geom%gmask_c0a(ic0a,il0)
   end if
end if
call mpl%f_comm%broadcast(valid,iproc-1)

if (valid) then
   ! Broadcast data
   call mpl%f_comm%broadcast(ic0a,iproc-1)
   call mpl%f_comm%broadcast(gmask,iproc-1)
else
   ! Missing values
   iproc = mpl%msv%vali
   ic0a = mpl%msv%vali
   gmask = .false.
end if

end subroutine geom_index_from_lonlat

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
integer :: idir,il0,il0dir,iprocdir,ic0adir
logical :: valid
character(len=1024),parameter :: subr = 'geom_define_dirac'

write(mpl%info,'(a7,a)') '','Define Dirac parameters'
call mpl%flush

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

   ! Index from lon/lat/level
   call geom%index_from_lonlat(mpl,nam%londir(idir),nam%latdir(idir),il0dir,iprocdir,ic0adir,valid)

   if (valid) then
      ! Add valid dirac point
      geom%ndir = geom%ndir+1
      geom%londir(geom%ndir) = nam%londir(idir)
      geom%latdir(geom%ndir) = nam%latdir(idir)
      geom%iprocdir(geom%ndir) = iprocdir
      geom%ic0adir(geom%ndir) = ic0adir
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
      if (geom%gmask_c0a(ic0a,il0)) then
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
      if (.not.geom%gmask_c0a(ic0a,il0)) fld_c0a(ic0a,il0) = mpl%msv%valr
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
subroutine geom_compute_deltas(geom,ic0u,il0,jc0u,jl0,dx,dy,dz)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
integer,intent(in) :: ic0u          ! First horizontal index, universe
integer,intent(in) :: il0           ! First vertical index
integer,intent(in) :: jc0u          ! Second horizontal index, universe
integer,intent(in) :: jl0           ! Second vertical index
real(kind_real),intent(out) :: dx   ! Longitude delta
real(kind_real),intent(out) :: dy   ! Latitude delta
real(kind_real),intent(out) :: dz   ! Altitude delta

! Compute deltas
dx = geom%lon_c0u(jc0u)-geom%lon_c0u(ic0u)
dy = geom%lat_c0u(jc0u)-geom%lat_c0u(ic0u)
call lonlatmod(dx,dy)
dx = dx*cos(geom%lat_c0u(ic0u))
dz = real(geom%vunit_c0u(ic0u,jl0)-geom%vunit_c0u(ic0u,il0),kind_real)

end subroutine geom_compute_deltas

!----------------------------------------------------------------------
! Subroutine: geom_rand_level
! Purpose: select random level
!----------------------------------------------------------------------
subroutine geom_rand_level(geom,mpl,rng,il0)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl ! MPI data
type(rng_type),intent(inout) :: rng ! Random number generator
integer,intent(out) :: il0           ! Level

! Resynchronize random number generator
call rng%resync(mpl)

! Generate random level
call rng%rand_integer(1,geom%nl0,il0)

! Desynchronize random number generator
call rng%desync(mpl)

end subroutine geom_rand_level

!----------------------------------------------------------------------
! Subroutine: geom_rand_point
! Purpose: select random valid point on the horizontal grid
!----------------------------------------------------------------------
subroutine geom_rand_point(geom,mpl,rng,il0,iproc,ic0a,nr)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl ! MPI data
type(rng_type),intent(inout) :: rng ! Random number generator
integer,intent(in) :: il0           ! Level
integer,intent(out) :: iproc        ! Processor
integer,intent(out) :: ic0a         ! Local index
integer,intent(out),optional :: nr  ! Number of random tries

! Local variables
integer :: lnr
real(kind_real) :: lon,lat
logical :: valid

! Initialization
valid = .false.
lnr = 0

! Loop
do while (.not.valid)
   ! Generate random lon/lat
   call rng%rand_real(-pi,pi,lon)
   call rng%rand_real(-1.0_kind_real,1.0_kind_real,lat)
   lat = asin(lat)

   ! Get index from lon/lat
   call geom%index_from_lonlat(mpl,lon,lat,il0,iproc,ic0a,valid)

   ! Update number of tries
   lnr = lnr+1
end do

! Set number of tries
if (present(nr)) nr = lnr

end subroutine geom_rand_point

!----------------------------------------------------------------------
! Function: geom_mg_to_proc
! Purpose: conversion from global to processor on model grid
!----------------------------------------------------------------------
function geom_mg_to_proc(geom,img)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
integer,intent(in) :: img           ! Global index

! Returned variable
integer :: geom_mg_to_proc

! Find processor
do geom_mg_to_proc=1,geom%nproc-1
   if ((geom%proc_to_mg_offset(geom_mg_to_proc)<img).and.(img<=geom%proc_to_mg_offset(geom_mg_to_proc+1))) return
end do

end function geom_mg_to_proc

!----------------------------------------------------------------------
! Function: geom_c0_to_c0a
! Purpose: conversion from global to halo A on subset Sc0
!----------------------------------------------------------------------
function geom_c0_to_c0a(geom,ic0)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
integer,intent(in) :: ic0           ! Global index

! Returned variable
integer :: geom_c0_to_c0a

! Local variable
integer :: iproc

! Find processor
iproc = geom%c0_to_proc(ic0)

! Get halo A index
geom_c0_to_c0a = ic0-geom%proc_to_c0_offset(iproc)

end function geom_c0_to_c0a

!----------------------------------------------------------------------
! Function: geom_c0_to_proc
! Purpose: conversion from global to processor on subset Sc0
!----------------------------------------------------------------------
function geom_c0_to_proc(geom,ic0)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
integer,intent(in) :: ic0           ! Global index

! Returned variable
integer :: geom_c0_to_proc

! Find processor
do geom_c0_to_proc=1,geom%nproc-1
   if ((geom%proc_to_c0_offset(geom_c0_to_proc)<ic0).and.(ic0<=geom%proc_to_c0_offset(geom_c0_to_proc+1))) return
end do

end function geom_c0_to_proc

!----------------------------------------------------------------------
! Function: geom_c0_to_c0u
! Purpose: conversion from global to universe on subset Sc0
!----------------------------------------------------------------------
function geom_c0_to_c0u(geom,ic0)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
integer,intent(in) :: ic0           ! Global index

! Returned variable
integer :: geom_c0_to_c0u

! Local variable
integer :: iproc,ic0a,offset,jproc

! Find processor
iproc = geom%c0_to_proc(ic0)

if (geom%myuniverse(iproc)) then
   ! Get halo A index
   ic0a = ic0-geom%proc_to_c0_offset(iproc)

   ! Compute universe offset
   offset = 0
   do jproc=1,iproc-1
      if (geom%myuniverse(jproc)) offset = offset+geom%proc_to_nc0a(jproc)
   end do

   ! Get universe index
   geom_c0_to_c0u = offset+ic0a
else
   ! Not in my universe
   geom_c0_to_c0u = 0
end if

end function geom_c0_to_c0u

end module type_geom
