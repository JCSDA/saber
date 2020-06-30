!----------------------------------------------------------------------
! Subroutine: model_mpas_coord
! Purpose: get MPAS coordinates
!----------------------------------------------------------------------
subroutine model_mpas_coord(model,mpl,nam)

implicit none

! Passed variables
class(model_type),intent(inout) :: model ! Model
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist

! Local variables
integer :: img
integer :: ncid,nmg_id,nlev_id,lon_id,lat_id,pres_id
real(kind=4),allocatable :: lon(:),lat(:),pres(:)
character(len=1024),parameter :: subr = 'model_mpas_coord'

! Open file and get dimensions
model%ntile = 1
model%nlon = mpl%msv%vali
model%nlat = mpl%msv%vali
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_share,ncid))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nCells',nmg_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nmg_id,len=model%nmg))
call mpl%ncerr(subr,nf90_inq_dimid(ncid,'nVertLevels',nlev_id))
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=model%nlev))

! Allocation
call model%alloc
allocate(lon(model%nmg))
allocate(lat(model%nmg))
allocate(pres(model%nlev))

! Read data and close file
call mpl%ncerr(subr,nf90_inq_varid(ncid,'lonCell',lon_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'latCell',lat_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,'pressure_base',pres_id))
call mpl%ncerr(subr,nf90_get_var(ncid,lon_id,lon))
call mpl%ncerr(subr,nf90_get_var(ncid,lat_id,lat))
call mpl%ncerr(subr,nf90_get_var(ncid,pres_id,pres))
call mpl%ncerr(subr,nf90_close(ncid))

! Model grid
model%mg_to_tile = 1
model%lon = lon
model%lat = lat
model%area = 4.0*pi/real(model%nmg,kind_real)
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

end subroutine model_mpas_coord

!----------------------------------------------------------------------
! Subroutine: model_mpas_read
! Purpose: read MPAS field
!----------------------------------------------------------------------
subroutine model_mpas_read(model,mpl,nam,filename,its,fld)

implicit none

! Passed variables
class(model_type),intent(inout) :: model                        ! Model
type(mpl_type),intent(inout) :: mpl                             ! MPI data
type(nam_type),intent(in) :: nam                                ! Namelist
character(len=*),intent(in) :: filename                         ! File name
integer,intent(in) :: its                                       ! Timeslot index
real(kind_real),intent(out) :: fld(model%nmga,model%nl0,nam%nv) ! Field

! Local variables
integer :: iv,il0,timeslot,ndims
integer :: ncid,fld_id
real(kind_real) :: fld_mg(model%nmg,model%nl0)
character(len=1024),parameter :: subr = 'model_mpas_read'

if (mpl%main) then
   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))

   ! Timeslot index
   read(nam%timeslots(its),'(i2)') timeslot
end if

do iv=1,nam%nv
   if (mpl%main) then
      ! Get variable id
      call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(nam%variables(iv)),fld_id))

      ! Check field size
      call mpl%ncerr(subr,nf90_inquire_variable(ncid,fld_id,ndims=ndims))

      ! Read data
      select case (ndims)
      case (2)
         ! 2D data
         fld_mg = mpl%msv%valr
         if (trim(nam%lev2d)=='first') then
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_mg(:,1),(/1,timeslot/),(/model%nmg,1/)))
         elseif (trim(nam%lev2d)=='last') then
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_mg(:,model%nl0),(/1,timeslot/),(/model%nmg,1/)))
         end if
      case (3) 
         ! 3D data
         do il0=1,nam%nl
            call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_mg(:,il0),(/nam%levs(il0),1,timeslot/),(/1,model%nmg,1/)))
         end do
      case default
         call mpl%abort(subr,'wrong number of dimensions for variable '//trim(nam%variables(iv)))
      end select
   end if
   call mpl%glb_to_loc(model%nl0,model%nmg,model%mg_to_proc,model%mg_to_mga,fld_mg,model%nmga,fld(:,:,iv))
end do

! Close file
if (mpl%main) call mpl%ncerr(subr,nf90_close(ncid))

end subroutine model_mpas_read
