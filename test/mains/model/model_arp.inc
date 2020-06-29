!----------------------------------------------------------------------
! Subroutine: model_arp_coord
! Purpose: get ARPEGE coordinates
!----------------------------------------------------------------------
subroutine model_arp_coord(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist

! Local variables
integer :: img,ilon,ilat
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,a_id,b_id
real(kind_real),allocatable :: lon(:,:),lat(:,:),a(:),b(:)
character(len=1024),parameter :: subr = 'model_arp_coord'

! Open file and get dimensions
model%ntile = 1
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'longitude',nlon_id))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'latitude',nlat_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=model%nlon))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=model%nlat))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'Z',nlev_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=model%nlev))

! Allocation
call model%alloc
allocate(lon(model%nlon,model%nlat))
allocate(lat(model%nlon,model%nlat))
allocate(a(model%nlev+1))
allocate(b(model%nlev+1))

! Read data and close file
call mpl%ncerr(subr,nf90_inq_varid(ncid,'longitude',lon_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'latitude',lat_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'hybrid_coef_A',a_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'hybrid_coef_B',b_id))
call mpl%ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call mpl%ncerr(subr,nf90_get_var(ncid,a_id,a))
call mpl%ncerr(subr,nf90_get_var(ncid,b_id,b))
call mpl%ncerr(subr,nf90_close(ncid))

! Grid size
model%nmg = count(lon>-1000.0)

! Convert to radian
lon = lon*deg2rad
lat = lat*deg2rad

! Model grid
img = 0
do ilon=1,model%nlon
   do ilat=1,model%nlat
      if (lon(ilon,ilat)>-1000.0) then
         img = img+1
         model%mg_to_tile(img) = 1
         model%mg_to_lon(img) = ilon
         model%mg_to_lat(img) = ilat
         model%lon(img) = lon(ilon,ilat)
         model%lat(img) = lat(ilon,ilat)
      end if
   end do
end do
model%area = 4.0*pi/real(model%nmg,kind_real)
model%mask = .true.

! Vertical unit
do img=1,model%nmg
   if (nam%logpres) then
      model%vunit(img,1:nam%nl) = log(0.5*(a(nam%levs(1:nam%nl))+a(nam%levs(1:nam%nl)+1)) &
                                & +0.5*(b(nam%levs(1:nam%nl))+b(nam%levs(1:nam%nl)+1))*ps)
      if (model%nl0>nam%nl) model%vunit(img,model%nl0) = log(ps)
   else
      model%vunit(img,:) = real(nam%levs(1:model%nl0))
   end if
end do

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(a)
deallocate(b)

end subroutine model_arp_coord

!----------------------------------------------------------------------
! Subroutine: model_arp_read
! Purpose: read ARPEGE field
!----------------------------------------------------------------------
subroutine model_arp_read(model,mpl,nam,filename,fld)

implicit none

! Passed variables
class(model_type),intent(inout) :: model                        ! Model
type(mpl_type),intent(inout) :: mpl                             ! MPI data
type(nam_type),intent(in) :: nam                                ! Namelist
character(len=*),intent(in) :: filename                         ! File name
real(kind_real),intent(out) :: fld(model%nmga,model%nl0,nam%nv) ! Field

! Local variables
integer :: iv,il0,img,ilon,ilat,info
integer :: ncid,fld_id
real(kind_real) :: fld_mg(model%nmg,model%nl0)
real(kind_real),allocatable :: fld_tmp(:,:,:)
character(len=3) :: ilchar
character(len=1024),parameter :: subr = 'model_arp_read'

if (mpl%main) then
   ! Allocation
   allocate(fld_tmp(model%nlon,model%nlat,model%nl0))

   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))
end if

do iv=1,nam%nv
   if (mpl%main) then
      do il0=1,nam%nl
         ! Get id
         write(ilchar,'(i3.3)') nam%levs(il0)
         info = nf90_inq_varid(ncid,'S'//ilchar//trim(nam%variables(iv)),fld_id)

         ! Read data
         if (info==nf90_noerr) then
            ! 3D data
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,il0)))
         else
            ! Check if the variable exists as 2D data
            call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(nam%variables(iv)),fld_id))

            ! 2D data
            if (((trim(nam%lev2d)=='first').and.(il0==1)).or.((trim(nam%lev2d)=='last').and.(il0==model%nl0))) then
               fld_tmp = mpl%msv%valr
               call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp(:,:,il0)))
            end if
         end if
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
   call mpl%glb_to_loc(model%nl0,model%nmg,model%mg_to_proc,model%mg_to_mga,fld_mg,model%nmga,fld(:,:,iv))
end do

if (mpl%main) then
   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))

   ! Release memory
   deallocate(fld_tmp)
end if

end subroutine model_arp_read
