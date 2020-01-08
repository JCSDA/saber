!----------------------------------------------------------------------
! Subroutine: model_fv3_coord
! Purpose: get FV3 coordinates
!----------------------------------------------------------------------
subroutine model_fv3_coord(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist

! Local variables
integer :: img,ilon,ilat,itile
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,phis_id,a_id,b_id
real(kind_real) :: z,h,pss
real(kind_real),allocatable :: lon(:,:,:),lat(:,:,:),phis(:,:,:),a(:),b(:)
real(kind_real),allocatable :: phis_mg(:)
character(len=1024),parameter :: subr = 'model_fv3_coord'

! Open file and get dimensions
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'xaxis_1',nlon_id))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'yaxis_2',nlat_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=model%nlon))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=model%nlat))
model%nmg = model%nlon*model%nlat*ntile
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'xaxis_1_ab',nlev_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=model%nlev))

! Allocation
call model%alloc
allocate(lon(model%nlon,model%nlat,ntile))
allocate(lat(model%nlon,model%nlat,ntile))
allocate(phis(model%nlon,model%nlat,ntile))
allocate(a(model%nlev))
allocate(b(model%nlev))
allocate(phis_mg(model%nmg))

! Read data and close file
call mpl%ncerr(subr,nf90_inq_varid(ncid,'grid_lon',lon_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'grid_lat',lat_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'phis',phis_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'ak',a_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'bk',b_id))
call mpl%ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call mpl%ncerr(subr,nf90_get_var(ncid,phis_id,phis))
call mpl%ncerr(subr,nf90_get_var(ncid,a_id,a))
call mpl%ncerr(subr,nf90_get_var(ncid,b_id,b))
call mpl%ncerr(subr,nf90_close(ncid))

! Convert to radian
lon = lon*deg2rad
lat = lat*deg2rad

! Model grid
img = 0
do itile=1,ntile
   do ilon=1,model%nlon
      do ilat=1,model%nlat
         img = img+1
         model%mg_to_lon(img) = ilon
         model%mg_to_lat(img) = ilat
         model%mg_to_tile(img) = itile
         model%lon(img) = lon(ilon,ilat,itile)
         model%lat(img) = lat(ilon,ilat,itile)
         call lonlatmod(model%lon(img),model%lat(img))
         phis_mg(img) = phis(ilon,ilat,itile)
      end do
   end do
end do
model%area = 4.0*pi/real(model%nmg,kind_real)
model%mask = .true.

! Vertical unit
do img=1,model%nmg
   if (nam%logpres) then
      if (.true.) then
         ! Get surface pressure from surface geopotential (https://github.com/jgoppert/pdas)
         z = phis_mg(img)/9.81
         h = z*req/(z+req)*1.0e-3
         pss = ((1.0+(-6.5/288.15)*h)**(34.163195/6.5))*ps
      else
         ! Flat surface
         pss = ps
      end if
      model%vunit(img,1:nam%nl) = log(0.5*(a(nam%levs(1:nam%nl))+a(nam%levs(1:nam%nl)+1)) &
                                  & +0.5*(b(nam%levs(1:nam%nl))+b(nam%levs(1:nam%nl)+1))*pss)
      if (model%nl0>nam%nl) model%vunit(img,model%nl0) = log(pss)
   else
      model%vunit(img,:) = real(nam%levs(1:model%nl0),kind_real)
   end if
end do

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(phis)
deallocate(a)
deallocate(b)
deallocate(phis_mg)

end subroutine model_fv3_coord

!----------------------------------------------------------------------
! Subroutine: model_fv3_read
! Purpose: read FV3 field
!----------------------------------------------------------------------
subroutine model_fv3_read(model,mpl,nam,filename,its,fld)

implicit none

! Passed variables
class(model_type),intent(inout) :: model                        ! Model
type(mpl_type),intent(inout) :: mpl                             ! MPI data
type(nam_type),intent(in) :: nam                                ! Namelist
character(len=*),intent(in) :: filename                         ! File name
integer,intent(in) :: its                                       ! Timeslot index
real(kind_real),intent(out) :: fld(model%nmga,model%nl0,nam%nv) ! Field

! Local variables
integer :: iv,il0,img,ilon,ilat,itile
integer :: ncid,fld_id
real(kind_real) :: fld_mg(model%nmg,model%nl0)
real(kind_real),allocatable :: fld_tmp(:,:,:,:)
character(len=1024),parameter :: subr = 'model_fv3_read'

if (mpl%main) then
   ! Allocation
   allocate(fld_tmp(model%nlon,model%nlat,model%nl0,ntile))

   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))
end if

do iv=1,nam%nv
   if (mpl%main) then
      ! 3d variable

      ! Get variable id
      call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(nam%varname(iv)),fld_id))

      ! Read data
      do itile=1,ntile
         do il0=1,nam%nl
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,il0,itile),(/1,1,nam%levs(il0),itile/), &
          & (/model%nlon,model%nlat,1,1/)))
         end do
      end do

      if (trim(nam%addvar2d(iv))/='') then
         ! 2d variable

         ! Get id
         call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

         ! Read data
         do itile=1,ntile
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,model%nl0,itile),(/1,1,itile/),(/model%nlon,model%nlat,1/)))
         end do
      end if

      ! Pack
      do il0=1,model%nl0
         do img=1,model%nmg
            ilon = model%mg_to_lon(img)
            ilat = model%mg_to_lat(img)
            itile = model%mg_to_tile(img)
            fld_mg(img,il0) = fld_tmp(ilon,ilat,il0,itile)
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

end subroutine model_fv3_read
