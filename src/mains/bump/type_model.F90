!----------------------------------------------------------------------
! Module: type_model
! Purpose: model routines
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_model

use netcdf
use tools_const, only: deg2rad,rad2deg,req,ps,pi
use tools_func, only: lonlatmod
use tools_kinds,only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use type_tree, only: tree_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type,nvmax
use type_rng, only: rng_type

implicit none

character(len=1024) :: zone = 'C+I'               ! Computation zone for AROME ('C', 'C+I' or 'C+I+E')
integer,parameter :: ntile = 6                    ! Number of tiles for FV3

! Member field derived type
type member_field_type
   real(kind_real),allocatable :: fld(:,:,:,:)    ! Ensemble perturbation
end type member_field_type

! Model derived type
type model_type
   ! Global dimensions
   integer :: nlon                                ! Longitude size
   integer :: nlat                                ! Latitude size
   integer :: nmg                                 ! Number of model grid points
   integer :: nlev                                ! Number of levels
   integer :: nl0                                 ! Number of levels in subset Sl0

   ! Packing arrays
   integer,allocatable :: mg_to_lon(:)            ! Model grid to longitude index
   integer,allocatable :: mg_to_lat(:)            ! Model grid to latgitude index
   integer,allocatable :: mg_to_tile(:)           ! Model grid to tile index

   ! Coordinates
   real(kind_real),allocatable :: lon(:)          ! Longitude
   real(kind_real),allocatable :: lat(:)          ! Latitude
   real(kind_real),allocatable :: area(:)         ! Area
   real(kind_real),allocatable :: vunit(:,:)      ! Vertical unit
   logical,allocatable :: mask(:,:)               ! Mask

   ! Local distribution
   integer :: nmga                                ! Halo A size for model grid
   integer,allocatable :: mg_to_proc(:)           ! Model grid to local task
   integer,allocatable :: mg_to_mga(:)            ! Model grid, global to halo A
   integer,allocatable :: mga_to_mg(:)            ! Model grid, halo A to global

   ! Local coordinates
   real(kind_real),allocatable :: lon_mga(:)      ! Longitude on model grid, halo A
   real(kind_real),allocatable :: lat_mga(:)      ! Latitude on model grid, halo A
   real(kind_real),allocatable :: area_mga(:)     ! Area on model grid, halo A
   real(kind_real),allocatable :: vunit_mga(:,:)  ! Vertical unit on model grid, halo A
   logical,allocatable :: mask_mga(:,:)           ! Mask on model grid, halo A
   logical,allocatable :: smask_mga(:,:)          ! Sampling mask on model grid, halo A

   ! Ensembles
   integer :: ens1_ne                             ! Ensemble 1 size
   integer :: ens1_nsub                           ! Ensemble 1 sub-ensembles number
   integer :: ens2_ne                             ! Ensemble 2 size
   integer :: ens2_nsub                           ! Ensemble 2 sub-ensembles number
   type(member_field_type),allocatable :: ens1(:) ! Ensemble 1 members
   type(member_field_type),allocatable :: ens2(:) ! Ensemble 2 members

   ! Observations locations
   integer :: nobsa                               ! Number of observations, halo A 
   real(kind_real),allocatable :: lonobs(:)       ! Observations longitudes, halo A
   real(kind_real),allocatable :: latobs(:)       ! Observations latitudes, halo A
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
include "model/model_aro.inc"
include "model/model_arp.inc"
include "model/model_fv3.inc"
include "model/model_gem.inc"
include "model/model_geos.inc"
include "model/model_gfs.inc"
include "model/model_ifs.inc"
include "model/model_mpas.inc"
include "model/model_nemo.inc"
include "model/model_qg.inc"
include "model/model_res.inc"
include "model/model_wrf.inc"

!----------------------------------------------------------------------
! Subroutine: model_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine model_alloc(model)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model

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
! Purpose: release memory
!----------------------------------------------------------------------
subroutine model_dealloc(model)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model

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
if (allocated(model%lon_mga)) deallocate(model%lon_mga)
if (allocated(model%lat_mga)) deallocate(model%lat_mga)
if (allocated(model%area_mga)) deallocate(model%area_mga)
if (allocated(model%mask_mga)) deallocate(model%mask_mga)
if (allocated(model%smask_mga)) deallocate(model%smask_mga)
if (allocated(model%ens1)) then
   do ie=1,size(model%ens1)
      if (allocated(model%ens1(ie)%fld)) deallocate(model%ens1(ie)%fld)
   end do
   deallocate(model%ens1)
end if
if (allocated(model%ens2)) then
   do ie=1,size(model%ens2)
      if (allocated(model%ens2(ie)%fld)) deallocate(model%ens2(ie)%fld)
   end do
   deallocate(model%ens2)
end if
if (allocated(model%lonobs)) deallocate(model%lonobs)
if (allocated(model%latobs)) deallocate(model%latobs)

end subroutine model_dealloc

!----------------------------------------------------------------------
! Subroutine: model_setup
! Purpose: setup model
!----------------------------------------------------------------------
subroutine model_setup(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(inout) :: nam      ! Namelist variables

! Local variables
integer :: iv,img,info,iproc,imga,nmga,ny,nres,iy,delta,ix,i,nv_save,ildw
integer :: ncid,nmg_id,mg_to_proc_id,mg_to_mga_id,lon_id,lat_id
integer :: nn_index(1)
integer,allocatable :: nx(:),imga_arr(:)
real(kind_real) :: dlat,dlon
real(kind_real),allocatable :: lon_center(:),lat_center(:),fld(:,:,:)
character(len=4) :: nprocchar
character(len=1024) :: varname,filename
character(len=1024),dimension(nvmax) :: varname_save,addvar2d_save
character(len=1024),parameter :: subr = 'model_define_distribution'
type(tree_type) :: tree

! Number of levels
model%nl0 = nam%nl
do iv=1,nam%nv
   if (trim(nam%addvar2d(iv))/='') model%nl0 = nam%nl+1
end do

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
   if (mpl%main) then
      ! Open file
      write(nprocchar,'(i4.4)') mpl%nproc
      filename = trim(nam%prefix)//'_distribution_'//nprocchar//'.nc'
      info = nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid)
   end if
   call mpl%f_comm%broadcast(info,mpl%rootproc-1)

   ! No points on the last task
   if (nam%check_no_point) then
      info = nf90_noerr+1
      if (mpl%main) call mpl%ncerr(subr,nf90_close(ncid))
   end if

   if (info==nf90_noerr) then
      ! Read local distribution
      write(mpl%info,'(a7,a,i4,a)') '','Read local distribution for: ',mpl%nproc,' MPI tasks'
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

      ! Allocation
      allocate(lon_center(mpl%nproc))
      allocate(lat_center(mpl%nproc))
      allocate(imga_arr(mpl%nproc))

      ! Define distribution centers using a regular splitting
      ny = nint(sqrt(real(mpl%nproc,kind_real)))
      if (ny**2<mpl%nproc) ny = ny+1
      allocate(nx(ny))
      nres = mpl%nproc
      do iy=1,ny
         delta = mpl%nproc/ny
         if (nres>(ny-iy+1)*delta) delta = delta+1
         nx(iy) = delta
         nres = nres-delta
      end do
      if (sum(nx)/=mpl%nproc) call mpl%abort(subr,'wrong number of tiles in define_distribution')
      dlat = (maxval(model%lat)-minval(model%lat))/ny
      iproc = 0
      do iy=1,ny
         dlon = (maxval(model%lon)-minval(model%lon))/nx(iy)
         do ix=1,nx(iy)
            iproc = iproc+1
            lat_center(iproc) = minval(model%lat)+(real(iy,kind_real)-0.5)*dlat
            lon_center(iproc) = minval(model%lon)+(real(ix,kind_real)-0.5)*dlon
         end do
      end do

      if (mpl%main) then
         ! Allocation
         call tree%alloc(mpl,mpl%nproc)

         ! Initialization
         call tree%init(lon_center,lat_center)

         ! Local processor
         do img=1,model%nmg
            call tree%find_nearest_neighbors(model%lon(img),model%lat(img),1,nn_index)
            model%mg_to_proc(img) = nn_index(1)
         end do

         ! Local index
         imga_arr = 0
         do img=1,model%nmg
            iproc = model%mg_to_proc(img)
            imga_arr(iproc) = imga_arr(iproc)+1
            model%mg_to_mga(img) = imga_arr(iproc)
         end do
      end if

      ! Broadcast distribution
      call mpl%f_comm%broadcast(model%mg_to_proc,mpl%rootproc-1)
      call mpl%f_comm%broadcast(model%mg_to_mga,mpl%rootproc-1)

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
         call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))

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

      ! Release memory
      deallocate(lon_center)
      deallocate(lat_center)
      deallocate(imga_arr)
   end if
end if

! Size of tiles
model%nmga = count(model%mg_to_proc==mpl%myproc)

! Allocation
allocate(model%mga_to_mg(model%nmga))
allocate(model%lon_mga(model%nmga))
allocate(model%lat_mga(model%nmga))
allocate(model%area_mga(model%nmga))
allocate(model%vunit_mga(model%nmga,model%nl0))
allocate(model%mask_mga(model%nmga,model%nl0))
allocate(model%smask_mga(model%nmga,model%nl0))

! Conversion
imga = 0
do img=1,model%nmg
   if (model%mg_to_proc(img)==mpl%myproc) then
      imga = imga+1
      model%mga_to_mg(imga) = img
      model%lon_mga(imga) = model%lon(img)
      model%lat_mga(imga) = model%lat(img)
      model%area_mga(imga) = model%area(img)
      model%vunit_mga(imga,:) = model%vunit(img,:)
      model%mask_mga(imga,:) = model%mask(img,:)
   end if
end do

! All points of the last task are masked
if (nam%check_no_point_mask) then
   if (mpl%myproc==mpl%nproc) model%mask_mga = .false.
end if

! Define sampling mask
model%smask_mga = model%mask_mga
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
         varname = nam%mask_type(1:i-1)
         filename = trim(nam%mask_type(i+1:))
         write(mpl%info,'(a7,a)') '','Read sampling mask from variable '//trim(varname)//' in file '//trim(filename)
         call mpl%flush
         write(mpl%info,'(a7,a,e10.3,a)') '','Threshold ',nam%mask_th(1),' used as a '//trim(nam%mask_lu(1))//' bound'
         call mpl%flush
   
         ! Save namelist parameters
         nv_save = nam%nv
         varname_save = nam%varname
         addvar2d_save = nam%addvar2d
   
         ! Set namelist parameters
         nam%nv = 1
         nam%varname(1) = varname
         if (nam%nl<model%nl0) nam%addvar2d(1) = varname
   
         ! Read file
         allocate(fld(model%nmga,model%nl0,nam%nv))
         call model%read(mpl,nam,filename,1,fld)

         ! Compute mask
         if (trim(nam%mask_lu(1))=='lower') then
            model%smask_mga = model%smask_mga.and.(fld(:,:,1)>nam%mask_th(1))
         elseif (trim(nam%mask_lu(1))=='upper') then
            model%smask_mga = model%smask_mga.and.(fld(:,:,1)<nam%mask_th(1))
         else
            call mpl%abort(subr,'mask_lu not recognized')
         end if

         ! Reset namelist parameters
         nam%nv = nv_save
         nam%varname = varname_save
         nam%addvar2d = addvar2d_save
      else
         call mpl%abort(subr,'mask_type should be formatted as VARIABLE@FILE to read a mask from file')
      end if
   end if
end select

end subroutine model_setup

!----------------------------------------------------------------------
! Subroutine: model_read
! Purpose: read member field
!----------------------------------------------------------------------
subroutine model_read(model,mpl,nam,filename,its,fld)

implicit none

! Passed variables
class(model_type),intent(inout) :: model                        ! Model
type(mpl_type),intent(inout) :: mpl                             ! MPI data
type(nam_type),intent(in) :: nam                                ! Namelist
character(len=*),intent(in) :: filename                         ! File name
integer,intent(in) :: its                                       ! Timeslot index
real(kind_real),intent(out) :: fld(model%nmga,model%nl0,nam%nv) ! Field

! Initialization
fld = mpl%msv%valr

! Select model
if (trim(nam%model)=='aro') call model%aro_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='arp') call model%arp_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='fv3') call model%fv3_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='gem') call model%gem_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='geos') call model%geos_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='gfs') call model%gfs_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='ifs') call model%ifs_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='mpas') call model%mpas_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='nemo') call model%nemo_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='qg') call model%qg_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='res') call model%res_read(mpl,nam,filename,its,fld)
if (trim(nam%model)=='wrf') call model%wrf_read(mpl,nam,filename,its,fld)

end subroutine model_read

!----------------------------------------------------------------------
! Subroutine: model_read_member
! Purpose: read member field
!----------------------------------------------------------------------
subroutine model_read_member(model,mpl,nam,filename,ie,fld)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl                                     ! MPI data
type(nam_type),intent(in) :: nam                                        ! Namelist
character(len=*),intent(in) :: filename                                 ! File name
integer,intent(in) :: ie                                                ! Ensemble member index
real(kind_real),intent(out) :: fld(model%nmga,model%nl0,nam%nv,nam%nts) ! Field

! Local variables
integer :: its
character(len=1024) :: fullname

! Initialization
fld = mpl%msv%valr

do its=1,nam%nts
   ! Define filename
   write(fullname,'(a,a,i2.2,a,i4.4,a)') trim(filename),'_',nam%timeslot(its),'_',ie,'.nc'

   ! Read file
   call model%read(mpl,nam,fullname,its,fld(:,:,:,its))
end do

end subroutine model_read_member

!----------------------------------------------------------------------
! Subroutine: model_load_ens
! Purpose: load ensemble data
!----------------------------------------------------------------------
subroutine model_load_ens(model,mpl,nam,filename)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist
character(len=*),intent(in) :: filename  ! Filename ('ens1' or 'ens2')

! Local variables
integer :: ne,ie,nsub,isub,ie_sub
real(kind_real),allocatable :: mean(:,:,:,:,:)
character(len=1024),parameter :: subr = 'model_load_ens'

! Initialization
ne = 0
nsub = 0

select case (trim(filename))
case ('ens1')
   ! Initialization
   ne = nam%ens1_ne
   nsub = nam%ens1_nsub
   model%ens1_ne = nam%ens1_ne
   model%ens1_nsub = nam%ens1_nsub

   if (ne>0) then
      ! Allocation
      allocate(model%ens1(ne))
      do ie=1,ne
         allocate(model%ens1(ie)%fld(model%nmga,model%nl0,nam%nv,nam%nts))
      end do

      ! Initialization
      do ie=1,ne
         model%ens1(ie)%fld = mpl%msv%valr
      end do
   end if
case ('ens2')
   ! Initialization
   ne = nam%ens2_ne
   nsub = nam%ens2_nsub
   model%ens2_ne = nam%ens2_ne
   model%ens2_nsub = nam%ens2_nsub

   if (ne>0) then
      allocate(model%ens2(ne))
      do ie=1,ne
         allocate(model%ens2(ie)%fld(model%nmga,model%nl0,nam%nv,nam%nts))
      end do

      ! Initialization
      do ie=1,ne
         model%ens2(ie)%fld = mpl%msv%valr
      end do
   end if
case default
   call mpl%abort(subr,'wrong filename in model_load_ens')
end select

! Allocation
if (ne>0) allocate(mean(model%nmga,model%nl0,nam%nv,nam%nts,nsub))

! Loop over sub-ensembles
do isub=1,nsub
   write(mpl%info,'(a7,a)') '','Reading ensemble member:'
   call mpl%flush(.false.)

   ! Loop over members for a given sub-ensemble
   do ie_sub=1,ne/nsub
      write(mpl%info,'(i4)') ie_sub
      call mpl%flush(.false.)

      ! Read member
      ie = ie_sub+(isub-1)*ne/nsub
      select case (trim(filename))
      case ('ens1')
         call model%read_member(mpl,nam,filename,ie_sub,model%ens1(ie)%fld)
      case ('ens2')
         call model%read_member(mpl,nam,filename,ie_sub,model%ens2(ie)%fld)
      end select
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
end do

end subroutine model_load_ens

!----------------------------------------------------------------------
! Subroutine: model_generate_obs
! Purpose: generate observations locations
!----------------------------------------------------------------------
subroutine model_generate_obs(model,mpl,rng,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(rng_type),intent(inout) :: rng      ! Random number generator
type(nam_type),intent(in) :: nam         ! Namelist

! Local variables
integer :: nres,delta,iproc,iobs,proc_to_nobsa(mpl%nproc),nn_index(10),img,inb
real(kind_real) :: nn_dist(10),dist_obs
logical :: valid
character(len=1024),parameter :: subr = 'model_generate_obs'
type(tree_type) :: tree

! Check observation number
if (nam%nobs<1) call mpl%abort(subr,'nobs should be positive for offline observation operator')

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

! Release memory
call tree%dealloc

end subroutine model_generate_obs

end module type_model
