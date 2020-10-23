!----------------------------------------------------------------------
! Module: type_model
!> Model routines
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_model

use atlas_module, only: atlas_field,atlas_fieldset,atlas_integer,atlas_real,atlas_functionspace
use netcdf
use tools_atlas, only: create_atlas_function_space
use tools_const, only: deg2rad,rad2deg,req,ps,pi
use tools_func, only: lonlatmod
use tools_kinds,only: kind_int,kind_real,nc_kind_real
use tools_qsort, only: qsort
use tools_repro, only: inf
use type_fieldset, only: fieldset_type
use type_tree, only: tree_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type,nvmax
use type_rng, only: rng_type

implicit none

real(kind_real),parameter :: qmin = 1.0e-6        ! Minimum specific humidity value in the log variable change
character(len=1024) :: zone = 'C+I'               ! Computation zone for AROME ('C', 'C+I' or 'C+I+E')

! Model derived type
type model_type
   ! Global dimensions
   integer :: nlon                             !< Longitude size
   integer :: nlat                             !< Latitude size
   integer :: ntile                            !< Number of tiles
   integer :: nmg                              !< Number of model grid points
   integer :: nlev                             !< Number of levels
   integer :: nl0                              !< Number of levels in subset Sl0

   ! Packing arrays
   integer,allocatable :: mg_to_lon(:)         !< Model grid to longitude index
   integer,allocatable :: mg_to_lat(:)         !< Model grid to latgitude index
   integer,allocatable :: mg_to_tile(:)        !< Model grid to tile index

   ! Coordinates
   real(kind_real),allocatable :: lon(:)       !< Longitude
   real(kind_real),allocatable :: lat(:)       !< Latitude
   real(kind_real),allocatable :: area(:)      !< Area
   real(kind_real),allocatable :: vunit(:,:)   !< Vertical unit
   logical,allocatable :: mask(:,:)            !< Mask

   ! Local distribution
   integer :: nmga                             !< Halo A size for model grid
   integer,allocatable :: mg_to_proc(:)        !< Model grid to local task
   integer,allocatable :: mg_to_mga(:)         !< Model grid, global to halo A
   integer,allocatable :: mga_to_mg(:)         !< Model grid, halo A to global

   ! ATLAS node columns
   type(atlas_functionspace) :: afunctionspace !< ATLAS function space

   ! Fieldset
   type(fieldset_type) :: fieldset             !< Fieldset

   ! Tiles distribution
   logical,allocatable :: tilepool(:,:)        !< Pool of task for each task
   integer :: mytile                           !< Tile handled by a given task
   integer,allocatable :: ioproc(:)            !< I/O task for each tile
   integer :: nmgt                             !< Number of model grid point on each tile
   integer,allocatable :: mga_to_mgt(:)        !< Model grid, halo A, to model grid on a tile
   integer,allocatable :: mgt_to_mg(:)         !< Model grid on a tile to model grid, global

   ! Ensembles
   type(fieldset_type),allocatable :: ens1(:)  !< Ensemble 1 members
   type(fieldset_type),allocatable :: ens2(:)  !< Ensemble 2 members

   ! Observations locations
   integer :: nobsa                            !< Number of observations, halo A
   real(kind_real),allocatable :: lonobs(:)    !< Observations longitudes, halo A
   real(kind_real),allocatable :: latobs(:)    !< Observations latitudes, halo A
contains
   ! Model specific procedures
   procedure :: aro_coord => model_aro_coord
   procedure :: aro_read => model_aro_read
   procedure :: arp_coord => model_arp_coord
   procedure :: arp_read => model_arp_read
   procedure :: fv3_coord => model_fv3_coord
   procedure :: fv3_read => model_fv3_read
   procedure :: gem_coord => model_gem_coord
   procedure :: gem_read => model_gem_read
   procedure :: geos_coord => model_geos_coord
   procedure :: geos_read => model_geos_read
   procedure :: gfs_coord => model_gfs_coord
   procedure :: gfs_read => model_gfs_read
   procedure :: ifs_coord => model_ifs_coord
   procedure :: ifs_read => model_ifs_read
   procedure :: mpas_coord => model_mpas_coord
   procedure :: mpas_read => model_mpas_read
   procedure :: nemo_coord => model_nemo_coord
   procedure :: nemo_read => model_nemo_read
   procedure :: norcpm_coord => model_norcpm_coord
   procedure :: norcpm_read => model_norcpm_read
   procedure :: qg_coord => model_qg_coord
   procedure :: qg_read => model_qg_read
   procedure :: res_coord => model_res_coord
   procedure :: res_read => model_res_read
   procedure :: wrf_coord => model_wrf_coord
   procedure :: wrf_read => model_wrf_read

   ! Generic procedures
   procedure :: alloc => model_alloc
   procedure :: dealloc => model_dealloc
   procedure :: setup => model_setup
   procedure :: read => model_read
   procedure :: read_member => model_read_member
   procedure :: load_ens => model_load_ens
   procedure :: generate_obs => model_generate_obs
end type model_type

private
public :: model_type

contains

! Include model interfaces
include 'model/model_aro.inc'
include 'model/model_arp.inc'
include 'model/model_fv3.inc'
include 'model/model_gem.inc'
include 'model/model_geos.inc'
include 'model/model_gfs.inc'
include 'model/model_ifs.inc'
include 'model/model_mpas.inc'
include 'model/model_nemo.inc'
include 'model/model_norcpm.inc'
include 'model/model_qg.inc'
include 'model/model_res.inc'
include 'model/model_wrf.inc'

!----------------------------------------------------------------------
! Subroutine: model_alloc
!> Allocation
!----------------------------------------------------------------------
subroutine model_alloc(model)

implicit none

! Passed variables
class(model_type),intent(inout) :: model !< Model

! Allocation
allocate(model%mg_to_lon(model%nmg))
allocate(model%mg_to_lat(model%nmg))
allocate(model%mg_to_tile(model%nmg))
allocate(model%lon(model%nmg))
allocate(model%lat(model%nmg))
allocate(model%area(model%nmg))
allocate(model%vunit(model%nmg,model%nl0))
allocate(model%mask(model%nmg,model%nl0))

end subroutine model_alloc

!----------------------------------------------------------------------
! Subroutine: model_dealloc
!> Release memory
!----------------------------------------------------------------------
subroutine model_dealloc(model)

implicit none

! Passed variables
class(model_type),intent(inout) :: model !< Model

! Local variables
integer :: ie

! Release memory
if (allocated(model%lon)) deallocate(model%lon)
if (allocated(model%lat)) deallocate(model%lat)
if (allocated(model%area)) deallocate(model%area)
if (allocated(model%mask)) deallocate(model%mask)
if (allocated(model%mg_to_proc)) deallocate(model%mg_to_proc)
if (allocated(model%mg_to_mga)) deallocate(model%mg_to_mga)
if (allocated(model%mga_to_mg)) deallocate(model%mga_to_mg)
if (allocated(model%tilepool)) deallocate(model%tilepool)
if (allocated(model%ioproc)) deallocate(model%ioproc)
if (allocated(model%mga_to_mgt)) deallocate(model%mga_to_mgt)
if (allocated(model%mgt_to_mg)) deallocate(model%mgt_to_mg)
if (allocated(model%ens1)) then
   do ie=1,size(model%ens1)
      call model%ens1(ie)%final()
   end do
   deallocate(model%ens1)
end if
if (allocated(model%ens2)) then
   do ie=1,size(model%ens2)
      call model%ens2(ie)%final()
   end do
   deallocate(model%ens2)
end if
if (allocated(model%lonobs)) deallocate(model%lonobs)
if (allocated(model%latobs)) deallocate(model%latobs)
call model%afunctionspace%final()
call model%fieldset%final()

end subroutine model_dealloc

!----------------------------------------------------------------------
! Subroutine: model_setup
!> Setup model
!----------------------------------------------------------------------
subroutine model_setup(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model !< Model
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(nam_type),intent(inout) :: nam      !< Namelist variables

! Local variables
integer :: img,info,iproc,imga,il0,nmga,ny,nres,iy,delta,ix,i,nv_save,ildw,itile,ilon,ilat,ilonsub,ilatsub,imgt
integer :: deltalon,deltalat,nprocpertile,nreslon,nreslat,nxmax
integer :: ncid,nmg_id,mg_to_proc_id,mg_to_mga_id,lon_id,lat_id
integer,allocatable :: nx(:),nlatsub(:),nlonsub(:,:),lonlattile_to_proc(:,:,:),imga_arr(:)
integer(kind_int),pointer :: gmask_ptr(:,:),smask_ptr(:,:)
real(kind_real) :: dlat,dlon
real(kind_real),allocatable :: lon_mga(:),lat_mga(:)
real(kind_real),allocatable :: lon_inf(:),lon_sup(:),lat_inf(:),lat_sup(:)
real(kind_real),pointer :: area_ptr(:),vunit_ptr(:,:),real_ptr(:,:)
logical :: maskval
logical :: found
character(len=1024) :: variables,filename
character(len=1024),dimension(nvmax) :: variables_save
character(len=1024),parameter :: subr = 'model_define'
type(atlas_field) :: afield,afield_area,afield_vunit,afield_gmask,afield_smask
type(fieldset_type) :: fieldset

! Number of levels
model%nl0 = nam%nl

! Select model
if (trim(nam%model)=='aro') call model%aro_coord(mpl,nam)
if (trim(nam%model)=='arp') call model%arp_coord(mpl,nam)
if (trim(nam%model)=='fv3') call model%fv3_coord(mpl,nam)
if (trim(nam%model)=='gem') call model%gem_coord(mpl,nam)
if (trim(nam%model)=='geos') call model%geos_coord(mpl,nam)
if (trim(nam%model)=='gfs') call model%gfs_coord(mpl,nam)
if (trim(nam%model)=='ifs') call model%ifs_coord(mpl,nam)
if (trim(nam%model)=='mpas') call model%mpas_coord(mpl,nam)
if (trim(nam%model)=='nemo') call model%nemo_coord(mpl,nam)
if (trim(nam%model)=='norcpm') call model%norcpm_coord(mpl,nam)
if (trim(nam%model)=='qg') call model%qg_coord(mpl,nam)
if (trim(nam%model)=='res') call model%res_coord(mpl,nam)
if (trim(nam%model)=='wrf') call model%wrf_coord(mpl,nam)

! Set longitude and latitude bounds
do img=1,model%nmg
   call lonlatmod(model%lon(img),model%lat(img))
end do

! Allocation
allocate(model%mg_to_proc(model%nmg))
allocate(model%mg_to_mga(model%nmg))

! Define distribution
if (mpl%nproc==1) then
   ! All points on a single processor
   model%mg_to_proc = 1
   do img=1,model%nmg
      model%mg_to_mga(img) = img
   end do
elseif (mpl%nproc>1) then
   ! Set file name
   write(filename,'(a,a,i6.6)') trim(nam%prefix),'_distribution_',mpl%nproc

   if (nam%check_no_point) then
      ! No points on the last task
      info = nf90_noerr+1
   else
      if (mpl%main) then
         ! Open file
         info = nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid)
      end if
      call mpl%f_comm%broadcast(info,mpl%rootproc-1)
   end if

   if (info==nf90_noerr) then
      ! Read local distribution
      write(mpl%info,'(a7,a,i6,a)') '','Read local distribution for: ',mpl%nproc,' MPI tasks'
      call mpl%flush

      if (mpl%main) then
         ! Get variables ID
         call mpl%ncerr(subr,nf90_inq_varid(ncid,'mg_to_proc',mg_to_proc_id))
         call mpl%ncerr(subr,nf90_inq_varid(ncid,'mg_to_mga',mg_to_mga_id))

         ! Read varaibles
         call mpl%ncerr(subr,nf90_get_var(ncid,mg_to_proc_id,model%mg_to_proc))
         call mpl%ncerr(subr,nf90_get_var(ncid,mg_to_mga_id,model%mg_to_mga))

         ! Close file
         call mpl%ncerr(subr,nf90_close(ncid))
      end if

      ! Broadcast distribution
      call mpl%f_comm%broadcast(model%mg_to_proc,mpl%rootproc-1)
      call mpl%f_comm%broadcast(model%mg_to_mga,mpl%rootproc-1)

      ! Check
      if (maxval(model%mg_to_proc)>mpl%nproc) call mpl%abort(subr,'wrong distribution')
   else
      ! Generate a distribution
      if (mpl%msv%isnot(model%ntile)) then
         ! Special case for model with tiles
         if (mod(mpl%nproc,model%ntile)/=0) call mpl%abort(subr, &
 & 'the number of MPI tasks should be a multiple of the number of tiles')
         nprocpertile = mpl%nproc/model%ntile
      else
         ! General case
         nprocpertile = mpl%nproc
      end if

      ! Define number of subtiles in x and y directions
      ny = nint(sqrt(real(nprocpertile,kind_real)))
      if (ny**2<nprocpertile) ny = ny+1
      allocate(nx(ny))
      nres = nprocpertile
      do iy=1,ny
         delta = nprocpertile/ny
         if (nres>(ny-iy+1)*delta) delta = delta+1
         nx(iy) = delta
         nres = nres-delta
      end do
      if (sum(nx)/=nprocpertile) call mpl%abort(subr,'wrong number of subtiles in define_distribution')

      if (mpl%msv%isnot(model%nlon).and.mpl%msv%isnot(model%nlat)) then
         ! Regular grid

         ! Allocation
         nxmax = maxval(nx)
         allocate(nlatsub(ny))
         allocate(nlonsub(nxmax,ny))
         allocate(lonlattile_to_proc(model%nlon,model%nlat,model%ntile))

         ! Split tile
         nlatsub = 0
         nlonsub = 0
         nreslat = model%nlat
         do iy=1,ny
            deltalat = model%nlat/ny
            if (nreslat>(ny-iy+1)*deltalat) deltalat = deltalat+1
            nlatsub(iy) = deltalat
            nreslat = nreslat-deltalat
            nreslon = model%nlon
            do ix=1,nx(iy)
               deltalon = model%nlon/nx(iy)
               if (nreslon>(nx(iy)-ix+1)*deltalon) deltalon = deltalon+1
               nlonsub(ix,iy) = deltalon
               nreslon = nreslon-deltalon
            end do
         end do

         ! Assign task
         iproc = 0
         do itile=1,model%ntile
            do iy=1,ny
               do ix=1,nx(iy)
                  iproc = iproc+1
                  do ilatsub=1,nlatsub(iy)
                     do ilonsub=1,nlonsub(ix,iy)
                        ilat = ilatsub
                        if (iy>1) ilat = ilat+sum(nlatsub(1:iy-1))
                        ilon = ilonsub
                        if (ix>1) ilon = ilon+sum(nlonsub(1:ix-1,iy))
                        lonlattile_to_proc(ilon,ilat,itile) = iproc
                     end do
                  end do
               end do
            end do
         end do

         ! Pack task
         do img=1,model%nmg
            itile = model%mg_to_tile(img)
            ilon = model%mg_to_lon(img)
            ilat = model%mg_to_lat(img)
            model%mg_to_proc(img) = lonlattile_to_proc(ilon,ilat,itile)
         end do

         ! Release memory
         deallocate(nlatsub)
         deallocate(nlonsub)
         deallocate(lonlattile_to_proc)
      else
         ! Unstructured grid
         if ((model%ntile>1).and.inf(sum(model%area),4.0*pi)) call mpl%abort(subr, &
 & 'unstructured grid requires a single tile and a global domain')

         ! Allocation
         allocate(lon_inf(nprocpertile))
         allocate(lon_sup(nprocpertile))
         allocate(lat_inf(nprocpertile))
         allocate(lat_sup(nprocpertile))

         ! Define distribution bounds using a regular splitting
         dlat = (maxval(model%lat)-minval(model%lat))/ny
         iproc = 0
         do iy=1,ny
            dlon = (maxval(model%lon)-minval(model%lon))/nx(iy)
            do ix=1,nx(iy)
               iproc = iproc+1
               lon_inf(iproc) = minval(model%lon)+real(ix-1,kind_real)*dlon
               lon_sup(iproc) = minval(model%lon)+real(ix,kind_real)*dlon
               lat_inf(iproc) = minval(model%lat)+real(iy-1,kind_real)*dlat
               lat_sup(iproc) = minval(model%lat)+real(iy,kind_real)*dlat
            end do
         end do

         ! Define distribution
         do img=1,model%nmg
            ! Processor index
            iproc = 1
            found = .false.
            do while (.not.found)
               if (iproc<=nprocpertile) then
                  if (((model%lon(img)>=lon_inf(iproc)).and.(model%lon(img)<=lon_sup(iproc)) &
 & .and.(model%lat(img)>=lat_inf(iproc)).and.(model%lat(img)<=lat_sup(iproc)))) then
                     model%mg_to_proc(img) = iproc
                     found = .true.
                  end if
               else
                  model%mg_to_proc(img) = nprocpertile
                  found = .true.
               end if
               iproc = iproc+1
            end do
         end do

         ! Release memory
         deallocate(lon_inf)
         deallocate(lon_sup)
         deallocate(lat_inf)
         deallocate(lat_sup)
      end if

      ! Allocation
      allocate(imga_arr(mpl%nproc))

      ! Local index
      imga_arr = 0
      do img=1,model%nmg
         iproc = model%mg_to_proc(img)
         imga_arr(iproc) = imga_arr(iproc)+1
         model%mg_to_mga(img) = imga_arr(iproc)
      end do

      ! Release memory
      deallocate(imga_arr)

      ! No points on the last task
      if (nam%check_no_point) then
         ! Count points on the penultimate task
         nmga = count(model%mg_to_proc==mpl%nproc-1)

         ! Move all points from the last to the penultimate task
         do img=1,model%nmg
            if (model%mg_to_proc(img)==mpl%nproc) then
               nmga = nmga+1
               model%mg_to_proc(img) = mpl%nproc-1
               model%mg_to_mga(img) = nmga
            end if
         end do
      end if

      ! Write distribution
      if (mpl%main) then
         ! Create file
         call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

         ! Write namelist parameters
         call nam%write(mpl,ncid)

         ! Define dimension
         call mpl%ncerr(subr,nf90_def_dim(ncid,'nmg',model%nmg,nmg_id))

         ! Define variables
         call mpl%ncerr(subr,nf90_def_var(ncid,'lon',nc_kind_real,(/nmg_id/),lon_id))
         call mpl%ncerr(subr,nf90_def_var(ncid,'lat',nc_kind_real,(/nmg_id/),lat_id))
         call mpl%ncerr(subr,nf90_def_var(ncid,'mg_to_proc',nf90_int,(/nmg_id/),mg_to_proc_id))
         call mpl%ncerr(subr,nf90_def_var(ncid,'mg_to_mga',nf90_int,(/nmg_id/),mg_to_mga_id))

         ! End definition mode
         call mpl%ncerr(subr,nf90_enddef(ncid))

         ! Write variables
         call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,model%lon*rad2deg))
         call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,model%lat*rad2deg))
         call mpl%ncerr(subr,nf90_put_var(ncid,mg_to_proc_id,model%mg_to_proc))
         call mpl%ncerr(subr,nf90_put_var(ncid,mg_to_mga_id,model%mg_to_mga))

         ! Close file
         call mpl%ncerr(subr,nf90_close(ncid))
      end if
   end if
end if

! Size of tiles
model%nmga = count(model%mg_to_proc==mpl%myproc)

! Allocation
allocate(model%mga_to_mg(model%nmga))
allocate(lon_mga(model%nmga))
allocate(lat_mga(model%nmga))

! Conversion
imga = 0
do img=1,model%nmg
   if (model%mg_to_proc(img)==mpl%myproc) then
      imga = imga+1
      model%mga_to_mg(imga) = img
      lon_mga(imga) = model%lon(img)
      lat_mga(imga) = model%lat(img)
   end if
end do

! Create ATLAS function space
call create_atlas_function_space(model%nmga,lon_mga,lat_mga,model%afunctionspace)

! Create fieldset
model%fieldset = atlas_fieldset()

! Add area
afield_area = model%afunctionspace%create_field(name='area',kind=atlas_real(kind_real),levels=0)
call afield_area%data(area_ptr)
call model%fieldset%add(afield_area)

! Add vertical unit
afield_vunit = model%afunctionspace%create_field(name='vunit',kind=atlas_real(kind_real),levels=model%nl0)
call afield_vunit%data(vunit_ptr)
call model%fieldset%add(afield_vunit)

! Add geometry mask
afield_gmask = model%afunctionspace%create_field(name='gmask',kind=atlas_integer(kind_int),levels=model%nl0)
call afield_gmask%data(gmask_ptr)
call model%fieldset%add(afield_gmask)

! Add sampling mask
afield_smask = model%afunctionspace%create_field(name='smask',kind=atlas_integer(kind_int),levels=model%nl0)
call afield_smask%data(smask_ptr)
call model%fieldset%add(afield_smask)

! Conversion
imga = 0
do img=1,model%nmg
   if (model%mg_to_proc(img)==mpl%myproc) then
      imga = imga+1
      model%mga_to_mg(imga) = img
      area_ptr(imga) = model%area(img)*req**2
      vunit_ptr(:,imga) = model%vunit(img,:)
      do il0=1,model%nl0
         if (model%mask(img,il0).and.(.not.nam%check_no_point_mask)) then
            gmask_ptr(il0,imga) = 1
         else
            gmask_ptr(il0,imga) = 0
         end if
      end do
   end if
end do

! Define sampling mask
smask_ptr = gmask_ptr
select case(trim(nam%mask_type))
case ('none','stddev')
   ! All points accepted in sampling
case ('ldwv')
   ! All points accepted in sampling
   do ildw=1,nam%nldwv
      if (nam%img_ldwv(ildw)>0) then
         ! Find lon/lat based on model grid index
         if (nam%img_ldwv(ildw)<=model%nmg) then
            nam%lon_ldwv(ildw) = model%lon(nam%img_ldwv(ildw))*rad2deg
            nam%lat_ldwv(ildw) = model%lat(nam%img_ldwv(ildw))*rad2deg
         else
            call mpl%abort(subr,'model grid index is positive but too large')
         end if
      end if
      write(mpl%info,'(a7,a,e15.8,a,e15.8)') '','Profile '//trim(nam%name_ldwv(ildw))//' required at lon/lat: ', &
 & nam%lon_ldwv(ildw),' / ',nam%lat_ldwv(ildw)
      call mpl%flush
      if (.not.any(model%mask(nam%img_ldwv(ildw),:))) call mpl%warning(subr,'profile '//trim(nam%name_ldwv(ildw))//' is not valid')
   end do
case default
   if (nam%mask_type(1:3)=='lat') then
      ! All points accepted in sampling
   else
      i = index(trim(nam%mask_type),'@')
      if (i>1) then
         ! Get variable and filename
         variables = nam%mask_type(1:i-1)
         filename = nam%mask_type(i+1:)
         write(mpl%info,'(a7,a)') '','Read sampling mask from variable '//trim(variables)//' in file '//trim(filename)
         call mpl%flush
         write(mpl%info,'(a7,a,e10.3,a)') '','Threshold ',nam%mask_th(1),' used as a '//trim(nam%mask_lu(1))//' bound'
         call mpl%flush

         ! Save namelist parameters
         nv_save = nam%nv
         variables_save = nam%variables

         ! Set namelist parameters
         nam%nv = 1
         nam%variables(1) = variables

         ! Read file
         call model%read(mpl,nam,filename,fieldset)

         ! Get data
         afield = fieldset%field(variables)
         call afield%data(real_ptr)

         ! Compute mask
         do il0=1,model%nl0
            do imga=1,model%nmga
               if (trim(nam%mask_lu(1))=='lower') then
                  maskval = (real_ptr(il0,imga)>nam%mask_th(1))
               elseif (trim(nam%mask_lu(1))=='upper') then
                  maskval = (real_ptr(il0,imga)<nam%mask_th(1))
               else
                  call mpl%abort(subr,'mask_lu not recognized')
               end if
               if (maskval) then
                  smask_ptr(il0,imga) = 1
               else
                  smask_ptr(il0,imga) = 0
               end if
            end do
         end do

         ! Reset namelist parameters
         nam%nv = nv_save
         nam%variables = variables_save

         ! Release pointers
         call afield%final()
         call fieldset%final()
      else
         call mpl%abort(subr,'mask_type should be formatted as VARIABLE@FILE to read a mask from file')
      end if
   end if
end select

if (model%ntile>1) then
   ! Allocation
   allocate(model%tilepool(mpl%nproc,model%ntile))
   allocate(model%ioproc(model%ntile))

   ! Tile/task bind
   model%tilepool = .false.
   do img=1,model%nmg
      itile = model%mg_to_tile(img)
      iproc = model%mg_to_proc(img)
      model%tilepool(iproc,itile) = .true.
   end do

   ! Check that a given task handles one tile only
   do iproc=1,mpl%nproc
      if (count(model%tilepool(iproc,:))>1) call mpl%abort(subr,'a given task should handle one tile only')
   end do

   do itile=1,model%ntile
      ! Assign tiles
      if (model%tilepool(mpl%myproc,itile)) model%mytile = itile

      ! Assign reading task
      do iproc=1,mpl%nproc
         if (model%tilepool(iproc,itile)) then
            model%ioproc(itile) = iproc
            exit
         end if
      end do
   end do

   ! Count points on tile
   model%nmgt = count(model%mg_to_tile==model%mytile)

   ! Allocation
   allocate(model%mga_to_mgt(model%nmgt))
   allocate(model%mgt_to_mg(model%nmgt))

   ! Conversion
   imgt = 0
   do img=1,model%nmg
      if (model%mg_to_tile(img)==model%mytile) then
         imgt = imgt+1
         model%mgt_to_mg(imgt) = img
         if (model%mg_to_proc(img)==mpl%myproc) then
            imga = model%mg_to_mga(img)
            model%mga_to_mgt(imga) = imgt
         end if
      end if
   end do
end if

! Release memory
deallocate(lon_mga)
deallocate(lat_mga)
call afield_area%final()
call afield_vunit%final()
call afield_gmask%final()
call afield_smask%final()

end subroutine model_setup

!----------------------------------------------------------------------
! Subroutine: model_read
!> Read member field
!----------------------------------------------------------------------
subroutine model_read(model,mpl,nam,filename,fieldset)

implicit none

! Passed variables
class(model_type),intent(inout) :: model      !< Model
type(mpl_type),intent(inout) :: mpl           !< MPI data
type(nam_type),intent(in) :: nam              !< Namelist
character(len=*),intent(in) :: filename       !< File name
type(fieldset_type),intent(inout) :: fieldset !< Fieldset

! Local variables
integer :: iv
real(kind_real) :: fld_mga(model%nmga,model%nl0,nam%nv)
real(kind_real),pointer :: real_ptr(:,:)
type(atlas_field) :: afield

! Select model
if (trim(nam%model)=='aro') call model%aro_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='arp') call model%arp_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='fv3') call model%fv3_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='gem') call model%gem_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='geos') call model%geos_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='gfs') call model%gfs_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='ifs') call model%ifs_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='mpas') call model%mpas_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='nemo') call model%nemo_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='norcpm') call model%norcpm_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='qg') call model%qg_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='res') call model%res_read(mpl,nam,filename,fld_mga)
if (trim(nam%model)=='wrf') call model%wrf_read(mpl,nam,filename,fld_mga)

! Add data into fieldset
do iv=1,nam%nv
   ! Create field
   afield = model%afunctionspace%create_field(name=nam%variables(iv),kind=atlas_real(kind_real),levels=model%nl0)

   ! Add field
   call fieldset%add(afield)

   ! Copy data
   call afield%data(real_ptr)
   real_ptr = transpose(fld_mga(:,:,iv))
end do

end subroutine model_read

!----------------------------------------------------------------------
! Subroutine: model_read_member
!> Read member field
!----------------------------------------------------------------------
subroutine model_read_member(model,mpl,nam,filename,ie,fieldset)

implicit none

! Passed variables
class(model_type),intent(inout) :: model     !< Model
type(mpl_type),intent(inout) :: mpl          !< MPI data
type(nam_type),intent(in) :: nam             !< Namelist
character(len=*),intent(in) :: filename      !< File name
integer,intent(in) :: ie                     !< Ensemble member index
type(fieldset_type),intent(out) :: fieldset  !< Fieldset

! Local variables
character(len=1024) :: fullname

! Create fieldset
fieldset = atlas_fieldset()

! Define filename
write(fullname,'(a,i6.6)') trim(filename)//'_',ie

! Read file
call model%read(mpl,nam,fullname,fieldset)

end subroutine model_read_member

!----------------------------------------------------------------------
! Subroutine: model_load_ens
!> Load ensemble data
!----------------------------------------------------------------------
subroutine model_load_ens(model,mpl,nam,filename)

implicit none

! Passed variables
class(model_type),intent(inout) :: model !< Model
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(nam_type),intent(in) :: nam         !< Namelist
character(len=*),intent(in) :: filename  !< Filename ('ens1' or 'ens2')

! Local variables
integer :: ne,ie,nsub,isub,ie_sub
character(len=1024),parameter :: subr = 'model_load_ens'

! Initialization
ne = 0
nsub = 0

select case (trim(filename))
case ('ens1')
   ! Initialization
   ne = nam%ens1_ne
   nsub = nam%ens1_nsub

   ! Allocation
   if (ne>0) allocate(model%ens1(ne))
case ('ens2')
   ! Initialization
   ne = nam%ens2_ne
   nsub = nam%ens2_nsub

   ! Allocation
   if (ne>0) allocate(model%ens2(ne))
case default
   call mpl%abort(subr,'wrong filename in model_load_ens')
end select

! Loop over sub-ensembles
do isub=1,nsub
   write(mpl%info,'(a7,a)') '','Reading ensemble member:'
   call mpl%flush(.false.)

   ! Loop over members for a given sub-ensemble
   do ie_sub=1,ne/nsub
      write(mpl%info,'(i6)') ie_sub
      call mpl%flush(.false.)

      ! Read member
      ie = ie_sub+(isub-1)*ne/nsub
      select case (trim(filename))
      case ('ens1')
         call model%read_member(mpl,nam,filename,ie_sub,model%ens1(ie))
      case ('ens2')
         call model%read_member(mpl,nam,filename,ie_sub,model%ens2(ie))
      end select
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
end do

end subroutine model_load_ens

!----------------------------------------------------------------------
! Subroutine: model_generate_obs
!> Generate observations locations
!----------------------------------------------------------------------
subroutine model_generate_obs(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model !< Model
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(nam_type),intent(in) :: nam         !< Namelist

! Local variables
integer :: nres,delta,iproc,iobs,proc_to_nobsa(mpl%nproc),nn_index(10),img,inb
real(kind_real) :: nn_dist(10),dist_obs
character(len=1024),parameter :: subr = 'model_generate_obs'
type(rng_type) :: rng
type(tree_type) :: tree

! Initialize random number generator
call rng%init(mpl,nam)

! Set number of observations for each processor
nres = nam%nobs
do iproc=1,mpl%nproc
   delta = nam%nobs/mpl%nproc
   if (nres>(mpl%nproc-iproc+1)*delta) delta = delta+1
   proc_to_nobsa(iproc) = delta
   nres = nres-delta
end do

! No observations on the last task
if (nam%check_no_obs) then
   if (mpl%nproc==1) call mpl%abort(subr,'at least 2 MPI tasks required for test_no_obs')
   proc_to_nobsa(mpl%nproc-1) = proc_to_nobsa(mpl%nproc-1)+proc_to_nobsa(mpl%nproc)
   proc_to_nobsa(mpl%nproc) = 0
end if

! Allocation
model%nobsa = proc_to_nobsa(mpl%myproc)
allocate(model%lonobs(model%nobsa))
allocate(model%latobs(model%nobsa))
call tree%alloc(mpl,model%nmg)

! Initialization
call tree%init(model%lon,model%lat)

! Generate random observation network
iobs = 1
do while (iobs<=model%nobsa)
   ! Generate observation
   call rng%rand_real(-pi,pi,model%lonobs(iobs))
   call rng%rand_real(-0.5*pi,0.5*pi,model%latobs(iobs))

   ! Check distance with nearest neighbor
   call tree%find_nearest_neighbors(model%lonobs(iobs),model%latobs(iobs),1,nn_index(1:1),nn_dist(1:1))
   img = nn_index(1)
   dist_obs = nn_dist(1)
   call tree%find_nearest_neighbors(model%lon(img),model%lat(img),10,nn_index,nn_dist)
   do inb=1,10
      if (nn_dist(inb)>0.0) then
         if (dist_obs<nn_dist(inb)) iobs=iobs+1
         exit
      end if
   end do
end do

! Convert to degrees
model%lonobs = model%lonobs*rad2deg
model%latobs = model%latobs*rad2deg

! Release memory
call tree%dealloc

end subroutine model_generate_obs

end module type_model
