!----------------------------------------------------------------------
! Module: type_geom
! Purpose: geometry derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_geom

use atlas_module
use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_min,fckit_mpi_max
use netcdf
use tools_atlas, only: field_to_fld
use tools_const, only: pi,req,deg2rad,rad2deg,reqkm
use tools_func, only: lonlatmod,lonlathash,sphere_dist,lonlat2xyz,xyz2lonlat,vector_product,vector_triple_product
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

   ! Geometry data on model grid, global
   integer :: nmg                                 ! Number of model grid points

   ! Universe
   logical,allocatable :: myuniverse(:)           ! MPI tasks in the universe of the local task

   ! Geometry data on model grid, universe
   integer :: nmgu                                ! Universe size for model grid

   ! Link between halo A and universe on model grid
   integer,allocatable :: mga_to_mgu(:)           ! Halo A to universe on model grid
   integer,allocatable :: mgu_to_mga(:)           ! Universe to halo A on model grid

   ! Geometry data on subset Sc0, universe
   integer :: nc0u                                ! Universe size for subset Sc0
   integer,allocatable :: c0u_to_proc(:)          ! Task
   real(kind_real),allocatable :: lon_c0u(:)      ! Longitudes
   real(kind_real),allocatable :: lat_c0u(:)      ! Latitudes
   real(kind_real),allocatable :: hash_c0u(:)     ! Longitudes/latitudes hash
   real(kind_real),allocatable :: area_c0u(:)     ! Areas
   real(kind_real),allocatable :: vunit_c0u(:,:)  ! Vertical unit
   logical,allocatable :: gmask_c0u(:,:)          ! Geometry mask
   logical,allocatable :: smask_c0u(:,:)          ! Sampling mask
   logical,allocatable :: gmask_hor_c0u(:)        ! Union of horizontal geometry masks
   real(kind_real),allocatable :: mdist_c0u(:,:)  ! Minimum distance to mask

   ! Link between model grid and subset Sc0 on universe
   integer,allocatable :: mgu_to_c0u(:)           ! Model grid to subset Sc0 on universe
   integer,allocatable :: c0u_to_mgu(:)           ! Subset Sc0 to model grid on universe
   logical :: same_grid                           ! Same grid and distribution flag

   ! Geometry data on subset Sc0, halo A
   integer :: nc0a                                ! Halo A size for subset Sc0
   real(kind_real),allocatable :: lon_c0a(:)      ! Longitudes
   real(kind_real),allocatable :: lat_c0a(:)      ! Latitudes
   real(kind_real),allocatable :: hash_c0a(:)     ! Longitudes/latitudes hash
   real(kind_real),allocatable :: area_c0a(:)     ! Areas
   real(kind_real),allocatable :: vunit_c0a(:,:)  ! Vertical unit
   logical,allocatable :: gmask_c0a(:,:)          ! Geometry mask
   logical,allocatable :: gmask_hor_c0a(:)        ! Union of horizontal geometry masks
   logical,allocatable :: smask_c0a(:,:)          ! Sampling mask
   real(kind_real),allocatable :: mdist_c0a(:,:)  ! Minimum distance to mask

   ! Link between halo A and universe on subset Sc0
   integer,allocatable :: c0a_to_c0u(:)           ! Halo A to universe on subset Sc0
   integer,allocatable :: c0u_to_c0a(:)           ! Universe to halo A on subset Sc0
   type(com_type) :: com_AU                       ! Communication between subset Sc0 and model grid

   ! Geometry data on subset Sc0, global
   integer :: nc0                                 ! Number of subset Sc0 points
   integer,allocatable :: c0_to_proc(:)           ! Task
   integer,allocatable :: proc_to_c0_offset(:)    ! Offset for a given task

   ! Link between halo A and global on subset Sc0
   integer,allocatable :: c0_to_c0a(:)            ! Subset Sc0, global to halo A
   integer,allocatable :: c0a_to_c0(:)            ! Subset Sc0, halo A to global

   ! Link between universe and global on subset Sc0 TODO: remove that
   integer,allocatable :: c0_to_c0u(:)            ! Subset Sc0, global to universe
   integer,allocatable :: c0u_to_c0(:)            ! Subset Sc0, universe to global

   ! Link between model grid and subset Sc0 on halo A
   type(com_type) :: com_mg                       ! Communication between subset Sc0 and model grid

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
   integer,allocatable :: itsdir(:)               ! Dirac timeslot
contains
   procedure :: partial_dealloc => geom_partial_dealloc
   procedure :: dealloc => geom_dealloc
   procedure :: setup => geom_setup
   procedure :: from_atlas => geom_from_atlas
   procedure :: define_universe => geom_define_universe
   procedure :: setup_universe => geom_setup_universe
   procedure :: setup_local => geom_setup_local
   procedure :: setup_independent_levels => geom_setup_independent_levels
   procedure :: setup_mask_distance => geom_setup_mask_distance
   procedure :: setup_mask_check => geom_setup_mask_check
   procedure :: index_from_lonlat => geom_index_from_lonlat
!   procedure :: remap => geom_remap
   procedure :: define_dirac => geom_define_dirac
   procedure :: check_arc => geom_check_arc
   procedure :: copy_c0a_to_mga => geom_copy_c0a_to_mga
   procedure :: geom_copy_mga_to_c0a_real
   procedure :: geom_copy_mga_to_c0a_logical
   generic :: copy_mga_to_c0a => geom_copy_mga_to_c0a_real,geom_copy_mga_to_c0a_logical
   procedure :: compute_deltas => geom_compute_deltas
   procedure :: rand_point => geom_rand_point
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
! TODO

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
call geom%define_universe(mpl,nam)

! Setup universe
call geom%setup_universe(mpl)

! Setup local distribution and communications
call geom%setup_local(mpl)

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

! Allocation
call geom%tree%alloc(mpl,geom%nc0u)

! Initialization
call geom%tree%init(geom%lon_c0u,geom%lat_c0u)

! Horizontal distance
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
write(mpl%info,'(a7,a,i8)') '','Model grid size:         ',geom%nmg
call mpl%flush
write(mpl%info,'(a7,a,i8)') '','Subset Sc0 size:         ',geom%nc0
call mpl%flush
write(mpl%info,'(a7,a,i6,a,f6.2,a)') '','Number of redundant points:    ',(geom%nmg-geom%nc0), &
 & ' (',real(geom%nmg-geom%nc0,kind_real)/real(geom%nmg,kind_real)*100.0,'%)'
call mpl%flush
write(mpl%info,'(a7,a,f7.1,a,f7.1)') '','Min. / max. longitudes:',lon_min*rad2deg,' / ',lon_max*rad2deg
call mpl%flush
write(mpl%info,'(a7,a,f7.1,a,f7.1)') '','Min. / max. latitudes: ',lat_min*rad2deg,' / ',lat_max*rad2deg
call mpl%flush
write(mpl%info,'(a7,a,f5.1,a)') '','Domain area (% of Earth area):',100.0*maxval(geom%area)/(4.0*pi),'%'
call mpl%flush
write(mpl%info,'(a7,a)') '','Valid points (% of total domain):'
call mpl%flush
do il0=1,geom%nl0
   write(mpl%info,'(a10,a,i3,a,f5.1,a)') '','Level ',nam%levs(il0),' ~> ', &
 & 100.0*real(geom%nc0_gmask(il0),kind_real)/real(geom%nc0,kind_real),'%'
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
   write(mpl%info,'(a13,a,i4,a,i8)') '','Task ',iproc,': ',count(geom%c0_to_proc==iproc)
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
allocate(geom%area_mga(geom%nmga))
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
      allocate(area_mga(geom%nmga,1))
      call field_to_fld(mpl,afield,area_mga)
      geom%area_mga = area_mga(:,1)
      deallocate(area_mga)
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
! Subroutine: geom_define_universe
! Purpose: define universe
!----------------------------------------------------------------------
subroutine geom_define_universe(geom,mpl,nam)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist

! Local variables
integer :: imga,iproc
real(kind_real) :: xbar,ybar,zbar,x,y,z,lonbar,latbar,dist
real(kind_real),allocatable :: proc_to_lonbar(:),proc_to_latbar(:)

! Allocation
allocate(proc_to_lonbar(mpl%nproc))
allocate(proc_to_latbar(mpl%nproc))
allocate(geom%myuniverse(mpl%nproc))

! Compute local barycenters coordinates
xbar = 0.0
ybar = 0.0
zbar = 0.0
do imga=1,geom%nmga
   call lonlat2xyz(mpl,geom%lon_mga(imga),geom%lat_mga(imga),x,y,z)
   xbar = xbar+x
   ybar = ybar+y
   zbar = zbar+z
end do
if (geom%nmga>0) then
   xbar = xbar/real(geom%nmga,kind_real)
   ybar = ybar/real(geom%nmga,kind_real)
   zbar = zbar/real(geom%nmga,kind_real)
   call xyz2lonlat(mpl,xbar,ybar,zbar,lonbar,latbar)
else
   lonbar = mpl%msv%vali
   latbar = mpl%msv%vali
end if
write(mpl%info,'(a7,a,f6.1,a,f6.1)') '','Local barycenter coordinates: ',lonbar*rad2deg,' / ',latbar*rad2deg
call mpl%flush

! Communication
call mpl%f_comm%allgather(lonbar,proc_to_lonbar)
call mpl%f_comm%allgather(latbar,proc_to_latbar)

! Compute distances between local barycenters to set universe
write(mpl%info,'(a7,a)') '','Tasks in my universe: '
call mpl%flush(.false.)
do iproc=1,mpl%nproc
   call sphere_dist(lonbar,latbar,proc_to_lonbar(iproc),proc_to_latbar(iproc),dist)
   geom%myuniverse(iproc) = (dist<0.4*pi) ! TODO: nam parameter
   if (geom%myuniverse(iproc)) then
      write(mpl%info,'(i5)') iproc
      call mpl%flush(.false.)
   end if
end do
write(mpl%info,'(a)') ''
call mpl%flush

end subroutine geom_define_universe

!----------------------------------------------------------------------
! Subroutine: geom_setup_universe
! Purpose: setup geometry on universe
!----------------------------------------------------------------------
subroutine geom_setup_universe(geom,mpl)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data

! Local variables
integer :: nmgu,iproc,img,imga,imgu,ic0u,diff_grid,diff_grid_tot,il0
integer :: proc_to_nmga(mpl%nproc)
integer,allocatable :: mga_to_mg(:),mgu_to_mg(:),mgu_to_proc(:),redundant(:),order(:)
real(kind_real),allocatable :: lon_mgu(:),lat_mgu(:),area_mgu(:),vunit_mgu(:,:),list(:)
logical,allocatable :: gmask_mgu(:,:),smask_mgu(:,:),mask_hor_mgu(:)
type(com_type) :: com_AU

! Communication
call mpl%f_comm%allgather(geom%nmga,proc_to_nmga)

! Model grid sizes
geom%nmg = sum(proc_to_nmga)
nmgu = sum(proc_to_nmga,mask=geom%myuniverse)

! Allocation
allocate(mga_to_mg(geom%nmga))
allocate(geom%mga_to_mgu(geom%nmga))
allocate(geom%mgu_to_mga(geom%nmga))
allocate(mgu_to_mg(nmgu))
allocate(mgu_to_proc(nmgu))
allocate(lon_mgu(nmgu))
allocate(lat_mgu(nmgu))
allocate(area_mgu(nmgu))
allocate(vunit_mgu(nmgu,geom%nl0))
allocate(gmask_mgu(nmgu,geom%nl0))
allocate(smask_mgu(nmgu,geom%nl0))

! Model grid conversions
img = 0
imgu = 0
do iproc=1,mpl%nproc
   do imga=1,proc_to_nmga(iproc)
      img = img+1
      if (iproc==mpl%myproc) mga_to_mg(imga) = img
      if (geom%myuniverse(iproc)) then
         imgu = imgu+1
         if (iproc==mpl%myproc) geom%mga_to_mgu(imga) = imgu
         geom%mgu_to_mga(imgu) = imga
         mgu_to_mg(imgu) = img
         mgu_to_proc(imgu) = iproc
      end if
   end do
end do

! Setup model grid communication, local to universe
call com_AU%setup(mpl,'com_AU',geom%nmga,geom%nmga,nmgu,geom%nmg,mga_to_mg,mga_to_mg,mgu_to_mg)

! Extend model grid, halo A to universe
call com_AU%ext(mpl,geom%lon_mga,lon_mgu)
call com_AU%ext(mpl,geom%lat_mga,lat_mgu)
call com_AU%ext(mpl,geom%area_mga,area_mgu)
call com_AU%ext(mpl,geom%nl0,geom%vunit_mga,vunit_mgu)
call com_AU%ext(mpl,geom%nl0,geom%gmask_mga,gmask_mgu)
call com_AU%ext(mpl,geom%nl0,geom%smask_mga,smask_mgu)

! Deallocate memory
deallocate(geom%lon_mga)
deallocate(geom%lat_mga)
deallocate(geom%area_mga)
deallocate(geom%vunit_mga)
deallocate(geom%gmask_mga)
deallocate(geom%smask_mga)

! Allocation
allocate(list(nmgu))
allocate(order(nmgu))
allocate(redundant(nmgu))
allocate(mask_hor_mgu(nmgu))

! Initialization
redundant = mpl%msv%vali

! Look for redundant points
write(mpl%info,'(a7,a)') '','Look for redundant points in the model grid'
call mpl%flush

! Define points order
do imgu=1,nmgu
   list(imgu) = lonlathash(lon_mgu(imgu),lat_mgu(imgu))
end do
call qsort(nmgu,list,order)

! Look for redundant points
do imgu=2,nmgu
   if (eq(list(imgu-1),list(imgu))) redundant(order(imgu)) = order(imgu-1)
end do

! Check for successive redundant points
do imgu=1,nmgu
   if (mpl%msv%isnot(redundant(imgu))) then
      do while (mpl%msv%isnot(redundant(redundant(imgu))))
         redundant(imgu) = redundant(redundant(imgu))
      end do
   end if
end do

! Horizontal model grid mask
mask_hor_mgu = mpl%msv%is(redundant)

! Count subset Sc0 points on universe
geom%nc0u = count(mask_hor_mgu)

! Check grid similarity
if (geom%nc0u==geom%nmgu) then
   diff_grid = 0
else
   diff_grid = 1
end if
call mpl%f_comm%allreduce(diff_grid,diff_grid_tot,fckit_mpi_sum())
geom%same_grid = (diff_grid_tot==0)

! Allocation
allocate(geom%mgu_to_c0u(geom%nmgu))
allocate(geom%c0u_to_mgu(geom%nc0u))

! Initialization
geom%mgu_to_c0u = mpl%msv%vali

! Conversion
ic0u = 0
do imgu=1,nmgu
   if (mask_hor_mgu(imgu)) then
      ic0u = ic0u+1
      geom%mgu_to_c0u(imgu) = ic0u
      geom%c0u_to_mgu(ic0u) = imgu
   end if
end do

! Deal with successive redundant points
do imgu=1,nmgu
   if (mpl%msv%isnot(redundant(imgu))) geom%mgu_to_c0u(imgu) = geom%mgu_to_c0u(redundant(imgu))
end do

! Deal with mask on redundant points
do il0=1,geom%nl0
   do imgu=1,nmgu
      if (mpl%msv%isnot(redundant(imgu))) gmask_mgu(imgu,il0) = gmask_mgu(imgu,il0).or.gmask_mgu(redundant(imgu),il0)
   end do
end do

! Allocation
allocate(geom%lon_c0u(geom%nc0u))
allocate(geom%lat_c0u(geom%nc0u))
allocate(geom%hash_c0u(geom%nc0u))
allocate(geom%area_c0u(geom%nl0))
allocate(geom%vunit_c0u(geom%nc0u,geom%nl0))
allocate(geom%gmask_c0u(geom%nc0u,geom%nl0))
allocate(geom%smask_c0u(geom%nc0u,geom%nl0))
allocate(geom%gmask_hor_c0u(geom%nc0u))
allocate(geom%c0u_to_proc(geom%nc0u))

! Remove redundant points
geom%lon_c0u = lon_mgu(geom%c0u_to_mgu)
geom%lat_c0u = lat_mgu(geom%c0u_to_mgu)
geom%area_c0u = area_mgu(geom%c0u_to_mgu)/req**2
geom%vunit_c0u = vunit_mgu(geom%c0u_to_mgu,:)
geom%gmask_c0u = gmask_mgu(geom%c0u_to_mgu,:)
geom%smask_c0u = smask_mgu(geom%c0u_to_mgu,:)
geom%gmask_hor_c0u = any(geom%gmask_c0u,dim=2)

! Hash function
do ic0u=1,geom%nc0u
   geom%hash_c0u(ic0u) = lonlathash(geom%lon_c0u(ic0u),geom%lat_c0u(ic0u))
end do

! Local distribution
do ic0u=1,geom%nc0u
   imgu = geom%c0u_to_mgu(ic0u)
   iproc = mgu_to_proc(imgu)
   geom%c0u_to_proc(ic0u) = iproc
end do

! Release memory
deallocate(mga_to_mg)
deallocate(mgu_to_mg)
deallocate(mgu_to_proc)
deallocate(lon_mgu)
deallocate(lat_mgu)
deallocate(area_mgu)
deallocate(vunit_mgu)
deallocate(gmask_mgu)
deallocate(smask_mgu)
deallocate(redundant)
deallocate(mask_hor_mgu)
deallocate(list)
deallocate(order)

end subroutine geom_setup_universe

!----------------------------------------------------------------------
! Subroutine: geom_setup_local
! Purpose: setup geometry on local task
!----------------------------------------------------------------------
subroutine geom_setup_local(geom,mpl)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl    ! MPI data

! Local variables
integer :: ic0,ic0a,ic0u,iproc,imga,imgu,il0
integer :: proc_to_nc0a(mpl%nproc),ic0a_arr(mpl%nproc),nc0_gmask(0:geom%nl0),norm(geom%nl0),norm_tot(geom%nl0)
integer,allocatable :: mga_to_c0(:)
real(kind_real) :: areasum(geom%nl0),vunitsum(geom%nl0),vunitsum_tot(geom%nl0)

! Number of points on subset Sc0, halo A
geom%nc0a = count(geom%c0u_to_proc==mpl%myproc)

! Communication
call mpl%f_comm%allgather(geom%nc0a,proc_to_nc0a)

! Subset Sc0 offset for halo A
geom%proc_to_c0_offset(1) = 0
do iproc=2,mpl%nproc
   geom%proc_to_c0_offset(iproc) = geom%proc_to_c0_offset(iproc-1)+proc_to_nc0a(iproc-1)
end do

! Subset Sc0 global size
geom%nc0 = sum(proc_to_nc0a)

! Allocation
allocate(geom%c0a_to_c0u(geom%nc0a))
allocate(geom%c0u_to_c0a(geom%nc0u))
allocate(geom%c0_to_c0u(geom%nc0))
allocate(geom%c0u_to_c0(geom%nc0u))
allocate(geom%c0a_to_c0(geom%nc0a))
allocate(geom%c0_to_c0a(geom%nc0))
allocate(geom%c0_to_proc(geom%nc0))
allocate(mga_to_c0(geom%nmga))

! Subset Sc0 conversions
geom%c0u_to_c0a = mpl%msv%vali
geom%c0_to_c0u = mpl%msv%vali
ic0a_arr = 0
do ic0u=1,geom%nc0u
   ! Get task
   iproc = geom%c0u_to_proc(ic0u)

   ! Update local index
   ic0a_arr(iproc) = ic0a_arr(iproc)+1

   ! Global index
   ic0 = geom%proc_to_c0_offset(iproc)+ic0a_arr(iproc)

   ! Conversions between local and universe
   if (iproc==mpl%myproc) geom%c0a_to_c0u(ic0a_arr(iproc)) = ic0u
   geom%c0u_to_c0a(ic0u) = ic0a_arr(iproc)

   ! Conversions between global and universe
   geom%c0_to_c0u(ic0) = ic0u
   geom%c0u_to_c0(ic0) = ic0

   ! Conversions between global and local
   if (iproc==mpl%myproc) geom%c0a_to_c0(ic0a_arr(iproc)) = ic0
   geom%c0_to_c0a(ic0) = ic0a_arr(iproc)
   geom%c0_to_proc(ic0) = iproc
end do

! Other conversions
do imga=1,geom%nmga
   imgu = geom%mga_to_mgu(imga)
   ic0u = geom%mgu_to_c0u(imgu)
   ic0a = geom%c0u_to_c0a(ic0u)
   ic0 = geom%c0a_to_c0(ic0a)
   mga_to_c0(imga) = ic0
end do

! Setup subset Sc0 communication, local to universe
call geom%com_AU%setup(mpl,'com_AU',geom%nc0a,geom%nc0a,geom%nc0u,geom%nc0,geom%c0a_to_c0,geom%c0a_to_c0,geom%c0u_to_c0)

! Setup redundant points communication
call geom%com_mg%setup(mpl,'com_mg',geom%nc0a,geom%nc0a,geom%nmga,geom%nc0,geom%c0a_to_c0,geom%c0a_to_c0,mga_to_c0)

! Allocation
allocate(geom%lon_c0a(geom%nc0a))
allocate(geom%lat_c0a(geom%nc0a))
allocate(geom%hash_c0a(geom%nc0a))
allocate(geom%area_c0a(geom%nc0a))
allocate(geom%vunit_c0a(geom%nc0a,geom%nl0))
allocate(geom%gmask_c0a(geom%nc0a,geom%nl0))
allocate(geom%smask_c0a(geom%nc0a,geom%nl0))
allocate(geom%gmask_hor_c0a(geom%nc0a))
allocate(geom%nc0_gmask(0:geom%nl0))
allocate(geom%area(geom%nl0))
allocate(geom%vunitavg(geom%nl0))

! Define fields on halo A
geom%lon_c0a = geom%lon_c0u(geom%c0a_to_c0u)
geom%lat_c0a = geom%lat_c0u(geom%c0a_to_c0u)
geom%hash_c0a = geom%hash_c0u(geom%c0a_to_c0u)
geom%area_c0a = geom%area_c0a(geom%c0a_to_c0u)
geom%vunit_c0a = geom%vunit_c0u(geom%c0a_to_c0u,:)
geom%gmask_c0a = geom%gmask_c0u(geom%c0a_to_c0u,:)
geom%smask_c0a = geom%smask_c0u(geom%c0a_to_c0u,:)
geom%gmask_hor_c0a = geom%gmask_hor_c0u(geom%c0a_to_c0u)
nc0_gmask(0) = count(geom%gmask_hor_c0a)
nc0_gmask(1:geom%nl0) = count(geom%gmask_c0a,dim=1)
call mpl%f_comm%allreduce(nc0_gmask,geom%nc0_gmask,fckit_mpi_sum())
do il0=1,geom%nl0
   norm(il0) = count(geom%gmask_c0a(:,il0))
   if (norm(il0)>0) then
      areasum(il0) = sum(geom%area_c0a,mask=geom%gmask_c0a(:,il0))
      vunitsum(il0) = sum(geom%vunit_c0a(:,il0),mask=geom%gmask_c0a(:,il0))
   else
      vunitsum(il0) = 0.0
   end if
end do
call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
call mpl%f_comm%allreduce(areasum,geom%area,fckit_mpi_sum())
call mpl%f_comm%allreduce(vunitsum,vunitsum_tot,fckit_mpi_sum())
do il0=1,geom%nl0
   if (norm_tot(il0)>0) then
      geom%vunitavg(il0) = vunitsum_tot(il0)/real(norm_tot(il0),kind_real)
   else
      geom%vunitavg(il0) = mpl%msv%valr
   end if
end do

! Release memory
deallocate(mga_to_c0)

end subroutine geom_setup_local

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
integer :: il0,i,j,k,ic0u,jc0u,kc0u,iend,ibnda
integer,allocatable :: bnda_to_c0u(:,:)
real(kind_real) :: lon_arc(2),lat_arc(2),xbnda(2),ybnda(2),zbnda(2)
logical :: imask,jmask,kmask,init

! Allocation
allocate(geom%nbnda(0:geom%nl0))

! Count boundary arcs
do il0=0,geom%nl0
   geom%nbnda(il0) = 0
   do i=1,geom%mesh%n
      ic0u = geom%mesh%order(i)
      if (il0==0) then
         imask = geom%gmask_hor_c0u(ic0u)
      else
         imask = geom%gmask_c0u(ic0u,il0)
      end if
      if (.not.imask) then
         iend = geom%mesh%lend(i)
         init = .true.
         do while ((iend/=geom%mesh%lend(i)).or.init)
            j = abs(geom%mesh%list(iend))
            k = abs(geom%mesh%list(geom%mesh%lptr(iend)))
            jc0u = geom%mesh%order(j)
            kc0u = geom%mesh%order(k)
            if (il0==0) then
               jmask = geom%gmask_hor_c0u(jc0u)
               kmask = geom%gmask_hor_c0u(kc0u)
            else
               jmask = geom%gmask_c0u(jc0u,il0)
               kmask = geom%gmask_c0u(kc0u,il0)
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
allocate(bnda_to_c0u(2,maxval(geom%nbnda)))

do il0=1,geom%nl0
   ! Define boundary arcs
   ibnda = 0
   do i=1,geom%mesh%n
      ic0u = geom%mesh%order(i)
      if (il0==0) then
         imask = geom%gmask_hor_c0u(ic0u)
      else
         imask = geom%gmask_c0u(ic0u,il0)
      end if
      if (.not.imask) then
         iend = geom%mesh%lend(i)
         init = .true.
         do while ((iend/=geom%mesh%lend(i)).or.init)
            j = abs(geom%mesh%list(iend))
            k = abs(geom%mesh%list(geom%mesh%lptr(iend)))
            jc0u = geom%mesh%order(j)
            kc0u = geom%mesh%order(k)
            if (il0==0) then
                jmask = geom%gmask_hor_c0u(jc0u)
                kmask = geom%gmask_hor_c0u(kc0u)
            else
                jmask = geom%gmask_c0u(jc0u,il0)
                kmask = geom%gmask_c0u(kc0u,il0)
            end if
            if (.not.jmask.and.kmask) then
               ibnda = ibnda+1
               bnda_to_c0u(1,ibnda) = ic0u
               bnda_to_c0u(2,ibnda) = jc0u
            end if
            iend = geom%mesh%lptr(iend)
            init = .false.
         end do
      end if
   end do

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
integer :: nn_index(1),nn_proc,proc_to_nn_proc(mpl%nproc),jproc
real(kind_real) :: nn_dist(1),proc_to_nn_dist(mpl%nproc),distmin
logical :: valid

! Find nearest neighbor
call geom%tree%find_nearest_neighbors(lon,lat,1,nn_index,nn_dist)
nn_proc = geom%c0u_to_proc(nn_index(1))

! Communication
call mpl%f_comm%allgather(nn_dist(1),proc_to_nn_dist)
call mpl%f_comm%allgather(nn_proc,proc_to_nn_proc)

! The correct task should handle its own nearest neighbor, with the minimum distance
distmin = minval(proc_to_nn_dist)
do jproc=1,mpl%nproc
   if ((proc_to_nn_proc(jproc)==jproc).and.eq(proc_to_nn_dist(jproc),distmin)) iproc = jproc
end do

if (iproc==mpl%myproc) then
   ! Check whether the location is in the convex hull      
   call geom%mesh%inside(mpl,lon,lat,valid)

   ! Local index
   ic0a = nn_index(1)

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
! Subroutine: geom_rand_point
! Purpose: select random point on the grid
!----------------------------------------------------------------------
subroutine geom_rand_point(geom,mpl,rng,il0,iproc,ic0a,nr)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom ! Geometry
type(mpl_type),intent(inout) :: mpl ! MPI data
type(rng_type),intent(inout) :: rng ! Random number generator
integer,intent(in) :: il0           ! Level
integer,intent(out) :: iproc        ! Task
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

end module type_geom
