!----------------------------------------------------------------------
! Module: type_mesh
! Purpose: mesh derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_mesh


use atlas_module, only: atlas_unstructuredgrid,atlas_meshgenerator,atlas_mesh,atlas_build_edges, &
 & atlas_build_node_to_edge_connectivity,atlas_mesh_nodes,atlas_connectivity
use iso_fortran_env, only: int64
!$ use omp_lib
use tools_const, only: pi,req,rad2deg,reqkm
use tools_func, only: lonlathash,lonlatmod,sphere_dist,lonlat2xyz,xyz2lonlat,vector_product,det
use tools_kinds, only: kind_real
use tools_qsort, only: qsort
use tools_repro, only: repro
use type_mpl, only: mpl_type
use type_rng, only: rng_type
use type_tree, only: tree_type

implicit none

logical,parameter :: shuffle = .true. ! Shuffle mesh order (more efficient to compute the Delaunay triangulation)
integer,parameter :: nnmax = 20       ! Maximum number of nearest neighbors in the triangle search

! Connectivity table row derived type
type row_type
   integer :: cols                         ! Nomber of columns
   integer,allocatable :: nodes(:)         ! Nodes indices
   real(kind_real),allocatable :: dists(:) ! Nodes indices
end type row_type

! Mesh derived type
type mesh_type
   integer :: n                            ! Number of points
   integer,allocatable :: order(:)         ! Order of shuffled points
   integer,allocatable :: order_inv(:)     ! Inverse order of shuffled points
   real(kind_real),allocatable :: lon(:)   ! Points longitudes
   real(kind_real),allocatable :: lat(:)   ! Points latitudes
   real(kind_real),allocatable :: xyz(:,:) ! Cartesian coordinates
   type(atlas_mesh) :: amesh               ! ATLAS mesh
   type(row_type),allocatable :: rows(:)   ! Connectivity table rows
   integer :: maxcols                      ! Maximum number of columns
contains
   procedure :: alloc => mesh_alloc
   procedure :: init => mesh_init
   procedure :: dealloc => mesh_dealloc
   procedure :: copy => mesh_copy
   procedure :: barycentric => mesh_barycentric
   procedure :: count_bnda => mesh_count_bnda
   procedure :: get_bnda => mesh_get_bnda
end type mesh_type

private
public :: mesh_type

contains

!----------------------------------------------------------------------
! Subroutine: mesh_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine mesh_alloc(mesh,n)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh ! Mesh
integer,intent(in) :: n                ! Mesh size

! Allocation
mesh%n = n

! Allocation
allocate(mesh%order(mesh%n))
allocate(mesh%order_inv(mesh%n))
allocate(mesh%lon(mesh%n))
allocate(mesh%lat(mesh%n))
allocate(mesh%xyz(3,mesh%n))
allocate(mesh%rows(mesh%n))

end subroutine mesh_alloc

!----------------------------------------------------------------------
! Subroutine: mesh_init
! Purpose: intialization
!----------------------------------------------------------------------
subroutine mesh_init(mesh,mpl,rng,lon,lat)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh            ! Mesh
type(mpl_type),intent(inout) :: mpl               ! MPI data
type(rng_type),intent(inout) :: rng               ! Random number generator
real(kind_real),intent(in) :: lon(mesh%n)         ! Longitudes
real(kind_real),intent(in) :: lat(mesh%n)         ! Latitudes

! Local variables
integer :: i,j,k,cols,maxedge,jj
integer,allocatable :: jtab(:),nodes(:,:),order(:)
integer,pointer :: row(:)
integer(int64) :: hash_int64
real(kind_real) :: lonlat(2,mesh%n),hash,pert(2),rvec(3),costheta,sintheta,rvecxv(3)
real(kind_real),allocatable :: list(:),v(:,:)
character(len=1024),parameter :: subr = 'mesh_init'
type(atlas_unstructuredgrid) :: agrid
type(atlas_meshgenerator) :: ameshgenerator
type(atlas_mesh_nodes) :: anodes
type(atlas_connectivity) :: aconnectivity

! Points order
do i=1,mesh%n
   mesh%order(i) = i
end do

! Allocation
allocate(list(mesh%n))

! Reorder points
do i=1,mesh%n
   list(i) = lonlathash(lon(i),lat(i))
end do
call qsort(mesh%n,list,mesh%order)

! Release memory
deallocate(list)

if (shuffle) then
   ! Allocation
   allocate(jtab(mesh%n))

   ! Shuffle order (more efficient to compute the Delaunay triangulation)
   call rng%resync(mpl)
   call rng%rand_integer(1,mesh%n,jtab)
   call rng%desync(mpl)
   do i=mesh%n,2,-1
      k = mesh%order(jtab(i))
      mesh%order(jtab(i)) = mesh%order(i)
      mesh%order(i) = k
   end do

   ! Release memory
   deallocate(jtab)
end if

! Restrictive inverse order
mesh%order_inv = mpl%msv%vali
do i=1,mesh%n
   mesh%order_inv(mesh%order(i)) = i
end do

! Store coordinates
mesh%lon = lon(mesh%order)
mesh%lat = lat(mesh%order)

! Transform to cartesian coordinates
do i=1,mesh%n
   call lonlat2xyz(mpl,mesh%lon(i),mesh%lat(i),mesh%xyz(1,i),mesh%xyz(2,i),mesh%xyz(3,i))
end do

! Unstructured grid coordinates
do i=1,mesh%n
   ! Copy coordinates
   lonlat(1,i) = mesh%lon(i)
   lonlat(2,i) = mesh%lat(i)

   if (repro) then
      ! Compute hash
      hash = lonlathash(mesh%lon(i),mesh%lat(i))

      ! Generate pseudo-random numbers between -1e-12 and 1e-12 from the hash
      hash_int64 = int(hash*1.0e6,int64)
      call rng%lcg(pert(1),hash_int64)
      call rng%lcg(pert(2),hash_int64)
      pert = 1.0e-12*(2.0*pert-1.0)

      ! Apply perturbation for reproducibility with square cells
      lonlat(:,i) = lonlat(:,i)*(1.0+pert)

      ! Apply bounds
      call lonlatmod(lonlat(1,i),lonlat(2,i))
   end if

   ! Convert to degrees
   lonlat(:,i) = lonlat(:,i)*rad2deg
end do

! Create unstructured grid
agrid = atlas_unstructuredgrid(lonlat)

! Mesh generator
ameshgenerator = atlas_meshgenerator('delaunay')

! Create mesh from grid
mesh%amesh = ameshgenerator%generate(agrid)

! Build connectivity
call atlas_build_edges(mesh%amesh)
call atlas_build_node_to_edge_connectivity(mesh%amesh)

! Get connectivity
anodes = mesh%amesh%nodes()
aconnectivity = anodes%edge_connectivity()

! Check number of rows
if (aconnectivity%rows()/=mesh%n) call mpl%abort(subr,'number of rows should be equal to the number of points')

! Allocation
mesh%maxcols = aconnectivity%maxcols()
maxedge = 0
do i=1,mesh%n
   call aconnectivity%row(i,row,cols)
   mesh%rows(i)%cols = cols
   maxedge = max(maxedge,maxval(row(1:cols)))
   allocate(mesh%rows(i)%nodes(mesh%rows(i)%cols))
   allocate(mesh%rows(i)%dists(mesh%rows(i)%cols))
end do
allocate(nodes(maxedge,2))

! Fill nodes for each edge
nodes = mpl%msv%vali
do i=1,mesh%n
   call aconnectivity%row(i,row,cols)
   do j=1,mesh%rows(i)%cols
      k = row(j)
      if (mpl%msv%is(nodes(k,1))) then
         ! First point
         nodes(k,1) = i
      elseif (mpl%msv%is(nodes(k,2))) then
         ! Second point
         nodes(k,2) = i
      else
         ! Error
         call mpl%abort(subr,'third node for the same edge')
      end if
   end do
end do

! Fill nodes index in connectivity list
do i=1,mesh%n
   call aconnectivity%row(i,row,cols)
   do j=1,mesh%rows(i)%cols
      k = row(j)
      if (i==nodes(k,1)) then
         ! Node is first point
         mesh%rows(i)%nodes(j) = nodes(k,2)
      elseif (i==nodes(k,2)) then
         ! Node is second point
         mesh%rows(i)%nodes(j) = nodes(k,1)
      else
         ! Error
         call mpl%abort(subr,'no node for the this edge')
      end if
   end do
end do

! Release memory
deallocate(nodes)

! Sort nodes in counter-clockwise order
do i=1,mesh%n
   ! Allocation
   allocate(v(3,mesh%rows(i)%cols))
   allocate(list(mesh%rows(i)%cols))
   allocate(order(mesh%rows(i)%cols))

   ! Rotation vector in cartesian coordinates
   call lonlat2xyz(mpl,mesh%lon(i)-0.5*pi,0.0_kind_real,rvec(1),rvec(2),rvec(3))

   ! Rotation angle
   costheta = cos(0.5*pi-mesh%lat(i))
   sintheta = sin(0.5*pi-mesh%lat(i))


   ! Compute angle
   do j=1,mesh%rows(i)%cols
      jj = mesh%rows(i)%nodes(j)

      ! Rodrigues' rotation
      call vector_product(rvec,mesh%xyz(:,jj),rvecxv)
      v(:,j) = mesh%xyz(:,jj)*costheta+rvecxv*sintheta+rvec*sum(rvec*mesh%xyz(:,jj))*(1.0-costheta)

      ! Angle
      list(j) = atan2(v(2,j),v(1,j))
   end do

   ! Sort angles in counter-clockwise order
   call qsort(mesh%rows(i)%cols,list,order)
   mesh%rows(i)%nodes = mesh%rows(i)%nodes(order)

   ! Release memory
   deallocate(v)
   deallocate(list)
   deallocate(order)
end do

! Compute edge distances
do i=1,mesh%n
   do j=1,mesh%rows(i)%cols
      ! Neighbor node index
      k = mesh%rows(i)%nodes(j)

      ! Compute distance
      call sphere_dist(mesh%lon(i),mesh%lat(i),mesh%lon(k),mesh%lat(k),mesh%rows(i)%dists(j))
   end do
end do

end subroutine mesh_init

!----------------------------------------------------------------------
! Subroutine: mesh_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine mesh_dealloc(mesh)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh ! Mesh

! Local variables
integer :: i

! Release memory
if (allocated(mesh%order)) deallocate(mesh%order)
if (allocated(mesh%order_inv)) deallocate(mesh%order_inv)
if (allocated(mesh%lon)) deallocate(mesh%lon)
if (allocated(mesh%lat)) deallocate(mesh%lat)
if (allocated(mesh%xyz)) deallocate(mesh%xyz)
if (allocated(mesh%rows)) then
   do i=1,mesh%n
      if (allocated(mesh%rows(i)%nodes)) deallocate(mesh%rows(i)%nodes)
      if (allocated(mesh%rows(i)%dists)) deallocate(mesh%rows(i)%dists)
   end do
   deallocate(mesh%rows)
end if

end subroutine mesh_dealloc

!----------------------------------------------------------------------
! Subroutine: mesh_copy
! Purpose: copy
!----------------------------------------------------------------------
subroutine mesh_copy(mesh_out,mesh_in)

implicit none

! Passed variables
class(mesh_type),intent(inout) :: mesh_out ! Output mesh
type(mesh_type),intent(in) :: mesh_in      ! Input mesh

! Local variables
integer :: i

! Release memory
call mesh_out%dealloc

! Allocation
call mesh_out%alloc(mesh_in%n)
do i=1,mesh_out%n
   if (allocated(mesh_in%rows(i)%nodes)) allocate(mesh_out%rows(i)%nodes(mesh_in%rows(i)%cols))
   if (allocated(mesh_in%rows(i)%dists)) allocate(mesh_out%rows(i)%dists(mesh_in%rows(i)%cols))
end do

! Copy data
mesh_out%order = mesh_in%order
mesh_out%order_inv = mesh_in%order_inv
mesh_out%lon = mesh_in%lon
mesh_out%lat = mesh_in%lat
mesh_out%xyz = mesh_in%xyz
do i=1,mesh_out%n
   mesh_out%rows(i)%cols = mesh_in%rows(i)%cols
   if (allocated(mesh_out%rows(i)%nodes)) mesh_out%rows(i)%nodes = mesh_in%rows(i)%nodes
   if (allocated(mesh_out%rows(i)%dists)) mesh_out%rows(i)%dists = mesh_in%rows(i)%dists
end do

end subroutine mesh_copy

!----------------------------------------------------------------------
! Subroutine: mesh_barycentric
! Purpose: compute barycentric coordinates (not optimal but sure to converge)
!----------------------------------------------------------------------
subroutine mesh_barycentric(mesh,mpl,lon,lat,tree,b,ib)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh ! Mesh
type(mpl_type),intent(inout) :: mpl ! MPI data
real(kind_real),intent(in) :: lon   ! Longitude
real(kind_real),intent(in) :: lat   ! Latitude
type(tree_type),intent(in) :: tree  ! Tree
real(kind_real),intent(out) :: b(3) ! Barycentric weights
integer,intent(out) :: ib(3)        ! Barycentric indices

! Local variables
integer :: i,j,jj,kk,nn,l
integer,allocatable :: nn_index(:)
real(kind_real) :: xyz(3)
real(kind_real),allocatable :: nn_dist(:)
logical :: found,cflag(3),valid
character(len=1024),parameter :: subr = 'mesh_barycentric'

! Initialization
found = .false.
nn = 1
ib = mpl%msv%vali

! Transform to cartesian coordinates
call lonlat2xyz(mpl,lon,lat,xyz(1),xyz(2),xyz(3))

do while ((.not.found).and.(nn<=nnmax))
   ! Allocation
   if (nn>mesh%n) call mpl%abort(subr,'cannot find a triangle containing this point')
   allocate(nn_index(nn))
   allocate(nn_dist(nn))

   ! Find the next nearest neighbor
   call tree%find_nearest_neighbors(lon,lat,nn,nn_index,nn_dist)
   i = mesh%order_inv(nn_index(nn))

   ! Check distance
   if (nn_dist(nn)>0.0) then
      ! Loop over neigbors
      do j=1,mesh%rows(i)%cols
         ! Find the two other vertices of the triangle
         jj = mesh%rows(i)%nodes(j)
         if (j<mesh%rows(i)%cols) then
            kk = mesh%rows(i)%nodes(j+1)
         else
            kk = mesh%rows(i)%nodes(1)
         end if

         ! Compute weights
         call det(mesh%xyz(:,jj),mesh%xyz(:,kk),xyz,b(1),cflag(1))
         call det(mesh%xyz(:,kk),mesh%xyz(:,i),xyz,b(2),cflag(2))
         call det(mesh%xyz(:,i),mesh%xyz(:,jj),xyz,b(3),cflag(3))

         ! Check if the points are colinear
         if (count(cflag)>=2) then
            ! At least two non-flat triangles, check if weights are positive
            valid = .true.
            do l=1,3
               if (cflag(l)) valid = valid.and.(b(l)>0.0_kind_real)
            end do
            if (valid) then
               ib(1) = i
               ib(2) = jj
               ib(3) = kk
               b = max(b,0.0_kind_real)
               b = b/sum(b)
               found = .true.
               exit
            end if
         elseif (count(cflag)==0) then
            ! Only flat triangles, pick the two closest nodes
            nn = 3
            deallocate(nn_index)
            deallocate(nn_dist)
            allocate(nn_index(nn))
            allocate(nn_dist(nn))
            call tree%find_nearest_neighbors(lon,lat,nn,nn_index,nn_dist)
            ib(1) = mesh%order_inv(nn_index(1))
            ib(2) = mesh%order_inv(nn_index(2))
            ib(3) = mesh%order_inv(nn_index(3))
            b(1) = nn_dist(2)
            b(2) = nn_dist(1)
            b(3) = 0.0
            b = b/sum(b)
            found = .true.
            exit
         end if
      end do
   else
      ib(1) = i
      ib(2) = i
      ib(3) = i
      b(1) = 1.0
      b(2) = 0.0
      b(3) = 0.0
   end if

   ! Release memory
   deallocate(nn_index)
   deallocate(nn_dist)

   ! Update number of neighbors
   nn = nn+1
end do

! Transform indices
if (mpl%msv%isallnot(ib)) ib = mesh%order(ib)

end subroutine mesh_barycentric

!----------------------------------------------------------------------
! Subroutine: mesh_count_bnda
! Purpose: count boundary arcs
!----------------------------------------------------------------------
subroutine mesh_count_bnda(mesh,gmask,nbnda)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh ! Mesh
logical,intent(in) :: gmask(mesh%n) ! Mask
integer,intent(out) :: nbnda        ! Number of boundary nodes

! Local variables
integer :: i,j,k,kk

! Initialiation
nbnda = 0

! Loop over points
do i=1,mesh%n
   if (.not.gmask(mesh%order(i))) then
      ! Loop over neigbors
      do j=1,mesh%rows(i)%cols
         k = mesh%rows(i)%nodes(j)
         if (j<mesh%rows(i)%cols) then
            kk = mesh%rows(i)%nodes(j+1)
         else
            kk = mesh%rows(i)%nodes(1)
         end if
         if (.not.gmask(mesh%order(k)).and.gmask(mesh%order(kk))) nbnda = nbnda+1
      end do
   end if
end do

end subroutine mesh_count_bnda

!----------------------------------------------------------------------
! Subroutine: mesh_get_bnda
! Purpose: get boundary arcs
!----------------------------------------------------------------------
subroutine mesh_get_bnda(mesh,gmask,nbnda,bnda_index)

implicit none

! Passed variables
class(mesh_type),intent(in) :: mesh        ! Mesh
logical,intent(in) :: gmask(mesh%n)        ! Mask
integer,intent(in) :: nbnda                ! Number of boundary nodes
integer,intent(out) :: bnda_index(2,nbnda) ! Boundary node index

! Local variables
integer :: i,j,k,kk,ibnda

! Initialiation
ibnda = 0
bnda_index = 0

! Loop over points
do i=1,mesh%n
   if (.not.gmask(mesh%order(i))) then
      ! Loop over neigbors
      do j=1,mesh%rows(i)%cols
         k = mesh%rows(i)%nodes(j)
         if (j<mesh%rows(i)%cols) then
            kk = mesh%rows(i)%nodes(j+1)
         else
            kk = mesh%rows(i)%nodes(1)
         end if
         if (.not.gmask(mesh%order(k)).and.gmask(mesh%order(kk))) then
            ibnda = ibnda+1
            bnda_index(1,ibnda) = mesh%order(i)
            bnda_index(2,ibnda) = mesh%order(k)
         end if
      end do
   end if
end do

end subroutine mesh_get_bnda

end module type_mesh
