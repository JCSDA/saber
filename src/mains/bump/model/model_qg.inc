!----------------------------------------------------------------------
! Subroutine: model_qg_coord
! Purpose: get QG coordinates
!----------------------------------------------------------------------
subroutine model_qg_coord(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist

! Local variables
integer :: img,ilon,ilat,il0
integer :: ncid,nx_id,ny_id,nz_id,lon_id,lat_id,z_id,area_id,lmask_id
real(kind_real),allocatable :: lon(:,:),lat(:,:),z(:),area(:,:),lmask(:,:,:)
character(len=1024),parameter :: subr = 'model_qg_coord'

! Open file and get dimensions
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nx',nx_id))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'ny',ny_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nx_id,len=model%nlon))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,ny_id,len=model%nlat))
model%nmg = (model%nlon+1)*model%nlat
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nz',nz_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nz_id,len=model%nlev))

! Allocation
call model%alloc
allocate(lon(model%nlon+1,model%nlat))
allocate(lat(model%nlon+1,model%nlat))
allocate(z(model%nlev))
allocate(area(model%nlon+1,model%nlat))
allocate(lmask(model%nlon+1,model%nlat,model%nlev))

! Read data and close file
call mpl%ncerr(subr,nf90_inq_varid(ncid,'lon',lon_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'lat',lat_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'z',z_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'area',area_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'lmask',lmask_id))
call mpl%ncerr(subr,nf90_get_var(ncid,lon_id,lon(1:model%nlon,:)))
call mpl%ncerr(subr,nf90_get_var(ncid,lat_id,lat(1:model%nlon,:)))
call mpl%ncerr(subr,nf90_get_var(ncid,z_id,z))
call mpl%ncerr(subr,nf90_get_var(ncid,area_id,area(1:model%nlon,:)))
call mpl%ncerr(subr,nf90_get_var(ncid,lmask_id,lmask(1:model%nlon,:,:)))
call mpl%ncerr(subr,nf90_close(ncid))

! Add redundant longitude for tests
lon(model%nlon+1,:) = lon(1,:)
lat(model%nlon+1,:) = lat(1,:)
area(model%nlon+1,:) = area(1,:)
lmask(model%nlon+1,:,:) = lmask(1,:,:)

! Convert to radian
lon = lon*deg2rad
lat = lat*deg2rad

! Model grid
img = 0
do ilon=1,model%nlon+1
   do ilat=1,model%nlat
      img = img+1
      model%mg_to_lon(img) = ilon
      model%mg_to_lat(img) = ilat
      model%lon(img) = lon(ilon,ilat)
      model%lat(img) = lat(ilon,ilat)
      model%area(img) = area(ilon,ilat)/req**2
      do il0=1,model%nl0
        model%mask(img,il0) = (lmask(ilon,ilat,nam%levs(il0))>0.5)
      end do
   end do
end do

! Vertical unit
do il0=1,model%nl0
   model%vunit(1:model%nmg,il0) = z(nam%levs(il0))
end do

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(z)
deallocate(area)
deallocate(lmask)

end subroutine model_qg_coord

!----------------------------------------------------------------------
! Subroutine: model_qg_read
! Purpose: read QG field
!----------------------------------------------------------------------
subroutine model_qg_read(model,mpl,nam,filename,its,fld)

implicit none

! Passed variables
class(model_type),intent(inout) :: model                        ! Model
type(mpl_type),intent(inout) :: mpl                             ! MPI data
type(nam_type),intent(in) :: nam                                ! Namelist
character(len=*),intent(in) :: filename                         ! File name
integer,intent(in) :: its                                       ! Timeslot index
real(kind_real),intent(out) :: fld(model%nmga,model%nl0,nam%nv) ! Field

! Local variables
integer :: iv,il0,img,ilon,ilat
integer :: ncid,fld_id
real(kind_real) :: fld_mg(model%nmg,model%nl0)
real(kind_real),allocatable :: fld_tmp(:,:,:)
character(len=1024),parameter :: subr = 'model_qg_read'

if (mpl%main) then
   ! Allocation
   allocate(fld_tmp(model%nlon+1,model%nlat,model%nl0))

   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))
end if

do iv=1,nam%nv
   if (mpl%main) then
      ! 3d variable

      ! Get variable id
      call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(nam%varname(iv)),fld_id))

      ! Read data
      do il0=1,nam%nl
         call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(1:model%nlon,:,il0),(/1,1,nam%levs(il0)/), &
       & (/model%nlon,model%nlat,1/)))
         fld_tmp(model%nlon+1,:,il0) = fld_tmp(1,:,il0)
      end do

      if (trim(nam%addvar2d(iv))/='') then
         ! 2d variable

         ! Get id
         call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

         ! Read data
         call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(1:model%nlon,:,model%nl0)))
         fld_tmp(model%nlon+1,:,model%nl0) = fld_tmp(1,:,model%nl0)
      end if

      ! Pack
      do il0=1,model%nl0
         do img=1,model%nmg
            ilon = model%mg_to_lon(img)
            ilat = model%mg_to_lat(img)
            fld_mg(img,il0) = fld_tmp(ilon,ilat,il0)
         end do
      end do
   end if
   call mpl%glb_to_loc(model%nl0,model%nmg,model%mg_to_proc,model%mg_to_mga,fld_mg,model%nmga,fld(:,:,iv))
end do

if (mpl%main) then
   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))

   ! Release memory
   deallocate(fld_tmp)
end if

end subroutine model_qg_read
