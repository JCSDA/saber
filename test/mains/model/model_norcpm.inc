!----------------------------------------------------------------------
! Subroutine: model_norcpm_coord
!> Get NorCPM coordinates
!----------------------------------------------------------------------
subroutine model_norcpm_coord(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model !< Model
type(mpl_type),intent(inout) :: mpl      !< MPI data
type(nam_type),intent(in) :: nam         !< Namelist

! Local variables
integer :: img,il0,ilon,ilat
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,area_id,gmask_id,dp_id
integer,allocatable :: gmask(:,:)
real(kind=8),allocatable :: lon(:,:),lat(:,:),area(:,:),dp(:,:,:)
character(len=1024),parameter :: subr = 'model_norcpm_coord'

! Open file and get dimensions
model%ntile = 1
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'x',nlon_id))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'y',nlat_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=model%nlon))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=model%nlat))
model%nmg = model%nlon*model%nlat
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'kk2',nlev_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=model%nlev))
model%nlev = model%nlev/2

! Allocation
call model%alloc
allocate(lon(model%nlon,model%nlat))
allocate(lat(model%nlon,model%nlat))
allocate(area(model%nlon,model%nlat))
allocate(gmask(model%nlon,model%nlat))
allocate(dp(model%nlon,model%nlat,model%nlev))

! Read data and close file
call mpl%ncerr(subr,nf90_inq_varid(ncid,'plon',lon_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'plat',lat_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'parea',area_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'pmask',gmask_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'dp',dp_id))
call mpl%ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call mpl%ncerr(subr,nf90_get_var(ncid,area_id,area))
call mpl%ncerr(subr,nf90_get_var(ncid,gmask_id,gmask))
call mpl%ncerr(subr,nf90_get_var(ncid,dp_id,dp,(/1,1,1/),(/1,1,model%nlev/)))
call mpl%ncerr(subr,nf90_close(ncid))

! Model grid
img = 0
do ilon=1,model%nlon
   do ilat=1,model%nlat
      img = img+1
      model%mg_to_tile(img) = 1
      model%mg_to_lon(img) = ilon
      model%mg_to_lat(img) = ilat
      model%lon(img) = real(lon(ilon,ilat),kind_real)*deg2rad
      model%lat(img) = real(lat(ilon,ilat),kind_real)*deg2rad
      model%area(img) = real(area(ilon,ilat),kind_real)/req**2
      model%mask(img,:) = (gmask(ilon,ilat)>0)
   end do
end do

! Vertical unit
do img=1,model%nmg
   ilon = model%mg_to_lon(img)
   ilat = model%mg_to_lat(img)
      
   do il0=1,nam%nl
      if (nam%levs(il0)==1) then
         model%vunit(img,il0) = ps+0.5*dp(ilon,ilat,1)
      else
         model%vunit(img,il0) = ps+sum(dp(ilon,ilat,1:nam%levs(il0)-1))+0.5*dp(ilon,ilat,nam%levs(il0))
      end if
   end do
end do

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(area)
deallocate(gmask)
deallocate(dp)

end subroutine model_norcpm_coord

!----------------------------------------------------------------------
! Subroutine: model_norcpm_read
!> Read NorCPM field
!----------------------------------------------------------------------
subroutine model_norcpm_read(model,mpl,nam,filename,fld)

implicit none

! Passed variables
class(model_type),intent(inout) :: model                        !< Model
type(mpl_type),intent(inout) :: mpl                             !< MPI data
type(nam_type),intent(in) :: nam                                !< Namelist
character(len=*),intent(in) :: filename                         !< File name
real(kind_real),intent(out) :: fld(model%nmga,model%nl0,nam%nv) !< Field

! Local variables
integer :: iv,il0,img,ilon,ilat,ndims
integer :: ncid,fld_id
real(kind_real) :: fld_mg(model%nmg,model%nl0)
real(kind_real),allocatable :: fld_tmp(:,:,:)
character(len=1024),parameter :: subr = 'model_norcpm_read'

if (mpl%main) then
   ! Allocation
   allocate(fld_tmp(model%nlon,model%nlat,model%nl0))

   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))
end if

do iv=1,nam%nv
   if (mpl%main) then
      ! Get variable id
      call mpl%ncerr(subr,nf90_inq_varid(ncid,nam%variables(iv),fld_id))

      ! Check field size
      call mpl%ncerr(subr,nf90_inquire_variable(ncid,fld_id,ndims=ndims))

      ! Read data
      select case (ndims)
      case (3)
         ! 2D data
         fld_tmp = mpl%msv%valr
         if (trim(nam%lev2d)=='first') then
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,1),(/1,1,1/), &
       & (/model%nlon,model%nlat,1/)))
         elseif (trim(nam%lev2d)=='last') then
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,model%nl0),(/1,1,1/), &
       & (/model%nlon,model%nlat,1/)))
         end if
      case (4) 
         ! 3D data
         do il0=1,nam%nl
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,il0),(/1,1,nam%levs(il0),1/), &
       & (/model%nlon,model%nlat,1,1/)))
         end do
      case default
         call mpl%abort(subr,'wrong number of dimensions for variable '//trim(nam%variables(iv)))
      end select

      ! Pack
      do il0=1,model%nl0
         do img=1,model%nmg
            ilon = model%mg_to_lon(img)
            ilat = model%mg_to_lat(img)
            fld_mg(img,il0) = fld_tmp(ilon,ilat,il0)
         end do
      end do
   end if
   call mpl%glb_to_loc(model%nl0,model%nmga,model%nmg,model%mga_to_mg,fld_mg,fld(:,:,iv))
end do

if (mpl%main) then
   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))

   ! Release memory
   deallocate(fld_tmp)
end if

end subroutine model_norcpm_read
