!----------------------------------------------------------------------
! Subroutine: model_wrf_coord
! Purpose: get WRF coordinates
!----------------------------------------------------------------------
subroutine model_wrf_coord(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist

! Local variables
integer :: img,ilon,ilat
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,pres_id
real(kind_real) :: dx,dy
real(kind_real),allocatable :: lon(:,:),lat(:,:),pres(:)
character(len=1024),parameter :: subr = 'model_wrf_coord'

! Open file and get dimensions
model%ntile = 1
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'west_east',nlon_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=model%nlon))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'south_north',nlat_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=model%nlat))
model%nmg = model%nlon*model%nlat
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'bottom_top',nlev_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=model%nlev))

! Allocation
call model%alloc
allocate(lon(model%nlon,model%nlat))
allocate(lat(model%nlon,model%nlat))
allocate(pres(model%nlev))

! Read data and close file
call mpl%ncerr(subr,nf90_inq_varid(ncid,'XLONG',lon_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'XLAT',lat_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'PB',pres_id))
call mpl%ncerr(subr,nf90_get_var(ncid,lon_id,lon,(/1,1,1/),(/model%nlon,model%nlat,1/)))
call mpl%ncerr(subr,nf90_get_var(ncid,lat_id,lat,(/1,1,1/),(/model%nlon,model%nlat,1/)))
call mpl%ncerr(subr,nf90_get_var(ncid,pres_id,pres))
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'DX',dx))
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'DY',dy))
call mpl%ncerr(subr,nf90_close(ncid))

! Convert to radian
lon = lon*deg2rad
lat = lat*deg2rad

! Model grid
img = 0
do ilon=1,model%nlon
   do ilat=1,model%nlat
      img = img+1
      model%mg_to_tile(img) = 1
      model%mg_to_lon(img) = ilon
      model%mg_to_lat(img) = ilat
      model%lon(img) = lon(ilon,ilat)
      model%lat(img) = lat(ilon,ilat)
   end do
end do
model%area = dx*dy/req**2
model%mask = .true.

! Vertical unit
do img=1,model%nmg
   if (nam%logpres) then
      model%vunit(img,1:nam%nl) = log(pres(nam%levs(1:nam%nl)))
      if (model%nl0>nam%nl) model%vunit(img,model%nl0) = log(ps)
   else
      model%vunit(img,:) = real(nam%levs(1:model%nl0),kind_real)
   end if
end do

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(pres)

end subroutine model_wrf_coord

!----------------------------------------------------------------------
! Subroutine: model_wrf_read
! Purpose: read WRF field
!----------------------------------------------------------------------
subroutine model_wrf_read(model,mpl,nam,filename,fld)

implicit none

! Passed variables
class(model_type),intent(inout) :: model                        ! Model
type(mpl_type),intent(inout) :: mpl                             ! MPI data
type(nam_type),intent(in) :: nam                                ! Namelist
character(len=*),intent(in) :: filename                         ! File name
real(kind_real),intent(out) :: fld(model%nmga,model%nl0,nam%nv) ! Field

! Local variables
integer :: iv,il0,img,ilon,ilat,ndims
integer :: ncid,fld_id
real(kind_real) :: fld_tmp2,fld_mg(model%nmg,model%nl0)
real(kind_real),allocatable :: fld_tmp(:,:,:)
character(len=1024),parameter :: subr = 'model_wrf_read'

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
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,1),(/1,1,1/),(/model%nlon,model%nlat,1/)))
         elseif (trim(nam%lev2d)=='last') then
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,model%nl0),(/1,1,1/),(/model%nlon,model%nlat,1/)))
         end if
      case (4) 
         ! 3D data
         do il0=1,nam%nl
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,il0),(/1,1,nam%levs(il0),1/),(/model%nlon,model%nlat,1,1/)))
            select case (trim(nam%variables(iv)))
            case ('U')
               do ilat=1,model%nlat
                  do ilon=1,model%nlon
                     call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/ilon+1,ilat,nam%levs(il0),1/)))
                     fld_tmp(ilon,ilat,il0) = 0.5*(fld_tmp(ilon,ilat,il0)+fld_tmp2)
                  end do
               end do
            case ('V')
               do ilat=1,model%nlat
                  do ilon=1,model%nlon
                     call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/ilon,ilat+1,nam%levs(il0),1/)))
                     fld_tmp(ilon,ilat,il0) = 0.5*(fld_tmp(ilon,ilat,il0)+fld_tmp2)
                  end do
               end do
            end select
         end do
      case default
         call mpl%abort(subr,'wrong number of dimensions for variable '//trim(nam%variables(iv)))
      end select


      ! Read data
      do il0=1,nam%nl
         call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,il0),(/1,1,nam%levs(il0),1/),(/model%nlon,model%nlat,1,1/)))
         select case (trim(nam%variables(iv)))
         case ('U')
            do ilat=1,model%nlat
               do ilon=1,model%nlon
                  call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/ilon+1,ilat,nam%levs(il0),1/)))
                  fld_tmp(ilon,ilat,il0) = 0.5*(fld_tmp(ilon,ilat,il0)+fld_tmp2)
               end do
            end do
         case ('V')
            do ilat=1,model%nlat
               do ilon=1,model%nlon
                  call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/ilon,ilat+1,nam%levs(il0),1/)))
                  fld_tmp(ilon,ilat,il0) = 0.5*(fld_tmp(ilon,ilat,il0)+fld_tmp2)
               end do
            end do
         end select
      end do

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

end subroutine model_wrf_read
