!----------------------------------------------------------------------
! Subroutine: model_nemo_coord
! Purpose: get NEMO coordinates
!----------------------------------------------------------------------
subroutine model_nemo_coord(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist

! Local variables
integer :: il0,img,ilat,ilon
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,tmask_id,e1t_id,e2t_id,e3t_id
integer(kind=1),allocatable :: tmask(:,:,:)
real(kind_real),allocatable :: lon(:,:),lat(:,:),e1t(:,:),e2t(:,:),e3t(:,:,:)
character(len=1024),parameter :: subr = 'model_nemo_coord'

! Open file and get dimensions
model%ntile = 1
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'x',nlon_id))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'y',nlat_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=model%nlon))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=model%nlat))
model%nmg = model%nlon*model%nlat
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'z',nlev_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=model%nlev))

! Allocation
call model%alloc
allocate(lon(model%nlon,model%nlat))
allocate(lat(model%nlon,model%nlat))
allocate(tmask(model%nlon,model%nlat,model%nl0))
allocate(e1t(model%nlon,model%nlat))
allocate(e2t(model%nlon,model%nlat))
allocate(e3t(model%nlon,model%nlat,model%nlev))

! Read data and close file
call mpl%ncerr(subr,nf90_inq_varid(ncid,'nav_lon',lon_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'nav_lat',lat_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'tmask',tmask_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'e1t',e1t_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'e2t',e2t_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'e3t',e3t_id))
call mpl%ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_get_var(ncid,lat_id,lat))
do il0=1,model%nl0
   call mpl%ncerr(subr,nf90_get_var(ncid,tmask_id,tmask(:,:,il0),(/1,1,nam%levs(il0),1/),(/model%nlon,model%nlat,1,1/)))
end do
call mpl%ncerr(subr,nf90_get_var(ncid,e1t_id,e1t,(/1,1,1/),(/model%nlon,model%nlat,1/)))
call mpl%ncerr(subr,nf90_get_var(ncid,e2t_id,e2t,(/1,1,1/),(/model%nlon,model%nlat,1/)))
do il0=1,model%nlev
   call mpl%ncerr(subr,nf90_get_var(ncid,e3t_id,e3t(:,:,il0),(/1,1,il0,1/),(/model%nlon,model%nlat,1,1/)))
end do
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
      model%area(img) = e1t(ilon,ilat)*e2t(ilon,ilat)/req**2
      do il0=1,model%nl0
         model%mask(img,il0) = (tmask(ilon,ilat,il0)>0)
      end do
   end do
end do

! Vertical unit
do img=1,model%nmg
   if (nam%logpres) then
      ilon = model%mg_to_lon(img)
      ilat = model%mg_to_lat(img)
      do il0=1,nam%nl
         if (nam%levs(il0)==1) then
            model%vunit(img,il0) = -0.5*e3t(ilon,ilat,1)
         else
            model%vunit(img,il0) = -sum(e3t(ilon,ilat,1:nam%levs(il0)-1))-0.5*e3t(ilon,ilat,nam%levs(il0))
         end if
      end do
      if (model%nl0>nam%nl) model%vunit(img,model%nl0) = 0.0
   else
      model%vunit(img,:) = real(nam%levs(1:model%nl0),kind_real)
   end if
end do

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(tmask)

end subroutine model_nemo_coord

!----------------------------------------------------------------------
! Subroutine: model_nemo_read
! Purpose: read NEMO field
!----------------------------------------------------------------------
subroutine model_nemo_read(model,mpl,nam,filename,its,fld)

implicit none

! Passed variables
class(model_type),intent(inout) :: model                        ! Model
type(mpl_type),intent(inout) :: mpl                             ! MPI data
type(nam_type),intent(in) :: nam                                ! Namelist
character(len=*),intent(in) :: filename                         ! File name
integer,intent(in) :: its                                       ! Timeslot index
real(kind_real),intent(out) :: fld(model%nmga,model%nl0,nam%nv) ! Field

! Local variables
integer :: iv,il0,img,ilon,ilat,timeslot,ndims
integer :: ncid,fld_id
real(kind_real) :: fld_tmp2,fld_mg(model%nmg,model%nl0)
real(kind_real),allocatable :: fld_tmp(:,:,:)
character(len=1024),parameter :: subr = 'model_nemo_read'

if (mpl%main) then
   ! Allocation
   allocate(fld_tmp(model%nlon,model%nlat,model%nl0))

   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))

   ! Timeslot index
   read(nam%timeslots(its),'(i2)') timeslot
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
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,1),(/1,1,timeslot/), &
          & (/model%nlon,model%nlat,1/)))
         elseif (trim(nam%lev2d)=='last') then
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,model%nl0),(/1,1,timeslot/), &
          & (/model%nlon,model%nlat,1/)))
         end if
      case (4) 
         ! 3D data
         do il0=1,nam%nl
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,il0),(/1,1,nam%levs(il0),timeslot/), &
          & (/model%nlon,model%nlat,1,1/)))
            select case (trim(nam%variables(iv)))
            case ('un')
               do ilat=1,model%nlat
                  do ilon=1,model%nlon
                     if (ilon==1) then
                        call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/model%nlon,ilat,nam%levs(il0),timeslot/)))
                     else
                        call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/ilon-1,ilat,nam%levs(il0),timeslot/)))
                     end if
                     fld_tmp(ilon,ilat,il0) = 0.5*(fld_tmp(ilon,ilat,il0)+fld_tmp2)
                  end do
               end do
            case ('vn')
               do ilat=1,model%nlat
                  do ilon=1,model%nlon
                     if (ilat==1) then
                        call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/ilon,model%nlat,nam%levs(il0),timeslot/)))
                     else
                        call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp2,(/ilon,ilat-1,nam%levs(il0),timeslot/)))
                     end if
                     fld_tmp(ilon,ilat,il0) = 0.5*(fld_tmp(ilon,ilat,il0)+fld_tmp2)
                  end do
               end do
            end select
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

end subroutine model_nemo_read
