!----------------------------------------------------------------------
! Subroutine: model_aro_coord
! Purpose: load AROME coordinates
!----------------------------------------------------------------------
subroutine model_aro_coord(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist

! Local variables
integer :: img,ilon,ilat
integer :: ncid,nlon_id,nlat_id,nlev_id,pp_id,lon_id,lat_id,cmask_id,a_id,b_id
real(kind_real) :: dx,dy
real(kind_real),allocatable :: lon(:,:),lat(:,:),cmask(:,:),a(:),b(:)
character(len=1024),parameter :: subr = 'model_aro_coord'

! Open file and get dimensions
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'X',nlon_id))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'Y',nlat_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=model%nlon))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=model%nlat))
model%nmg = model%nlon*model%nlat
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'Z',nlev_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=model%nlev))

! Allocation
call model%alloc
allocate(lon(model%nlon,model%nlat))
allocate(lat(model%nlon,model%nlat))
allocate(cmask(model%nlon,model%nlat))
allocate(a(model%nlev+1))
allocate(b(model%nlev+1))

! Read data and close file
call mpl%ncerr(subr,nf90_inq_varid(ncid,'longitude',lon_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'latitude',lat_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'cmask',cmask_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'hybrid_coef_A',a_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'hybrid_coef_B',b_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'Projection_parameters',pp_id))
call mpl%ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call mpl%ncerr(subr,nf90_get_var(ncid,cmask_id,cmask))
call mpl%ncerr(subr,nf90_get_var(ncid,a_id,a))
call mpl%ncerr(subr,nf90_get_var(ncid,b_id,b))
call mpl%ncerr(subr,nf90_get_att(ncid,pp_id,'x_resolution',dx))
call mpl%ncerr(subr,nf90_get_att(ncid,pp_id,'y_resolution',dy))
call mpl%ncerr(subr,nf90_close(ncid))

! Convert to radian
lon = lon*deg2rad
lat = lat*deg2rad

! Model grid
img = 0
do ilon=1,model%nlon
   do ilat=1,model%nlat
      img = img+1
      model%mg_to_lon(img) = ilon
      model%mg_to_lat(img) = ilat
      model%lon(img) = lon(ilon,ilat)
      model%lat(img) = lat(ilon,ilat)
      select case (trim(zone))
      case ('C')
         model%mask(img,:) = (cmask(ilon,ilat)>0.75)
      case ('C+I')
         model%mask(img,:) = (cmask(ilon,ilat)>0.25)
      case ('C+I+E')
         model%mask(img,:) = .true.
      case default
         call mpl%abort(subr,'wrong AROME zone')
      end select
   end do
end do
model%area = dx*dy/req**2

! Vertical unit
do img=1,model%nmg
   if (nam%logpres) then
      model%vunit(img,1:nam%nl) = log(0.5*(a(nam%levs(1:nam%nl))+a(nam%levs(1:nam%nl)+1)) &
                                & +0.5*(b(nam%levs(1:nam%nl))+b(nam%levs(1:nam%nl)+1))*ps)
      if (model%nl0>nam%nl) model%vunit(img,model%nl0) = log(ps)
   else
      model%vunit(img,:) = real(nam%levs(1:model%nl0),kind_real)
   end if
end do

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(cmask)
deallocate(a)
deallocate(b)

end subroutine model_aro_coord

!----------------------------------------------------------------------
! Subroutine: model_aro_read
! Purpose: read AROME field
!----------------------------------------------------------------------
subroutine model_aro_read(model,mpl,nam,filename,its,fld)

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
character(len=3) :: ilchar
character(len=1024),parameter :: subr = 'model_aro_read'

if (mpl%main) then
   ! Allocation
   allocate(fld_tmp(model%nlon,model%nlat,model%nl0))

   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))
end if

do iv=1,nam%nv
   if (mpl%main) then
      ! 3d variable
      do il0=1,nam%nl
         ! Get id
         write(ilchar,'(i3.3)') nam%levs(il0)
         call mpl%ncerr(subr,nf90_inq_varid(ncid,'S'//ilchar//trim(nam%varname(iv)),fld_id))

         ! Read data
         call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,il0)))
      end do

      if (trim(nam%addvar2d(iv))/='') then
         ! 2d variable

         ! Get id
         call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

         ! Read data
         call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,model%nl0)))

         ! Variable change for surface pressure
         if (trim(nam%addvar2d(iv))=='SURFPRESSION') fld_tmp(:,:,model%nl0) = exp(fld_tmp(:,:,model%nl0))
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

end subroutine model_aro_read
