!----------------------------------------------------------------------
! Module: type_io
! Purpose: I/O derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_io

use fckit_mpi_module, only: fckit_mpi_comm
use netcdf
use tools_const, only: pi,deg2rad,rad2deg,reqkm
use tools_kinds, only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use type_com, only: com_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

! I/O derived type
type io_type
   ! Field output
   integer :: nc0io                              ! I/O chunk size on subset Sc0
   integer :: ic0_s                              ! Starting index for the I/O chunk on subset Sc0
   integer,allocatable :: c0io_to_c0(:)          ! Subset Sc0, I/O chunk to global
   integer,allocatable :: procio_to_proc(:)      ! I/O task to main communicator task
   type(com_type) :: com_AIO                     ! Communication between halo A and I/O chunk

   ! Field regridding and output
   integer :: nlon                               ! Number of output grid longitudes
   integer :: nlat                               ! Number of output grid latitudes
   real(kind_real),allocatable :: lon(:)         ! Output grid longitudes
   real(kind_real),allocatable :: lat(:)         ! Output grid latitudes
   integer,allocatable :: og_to_lon(:)           ! Output grid to longitude index
   integer,allocatable :: og_to_lat(:)           ! Output grid to latitude index
   integer :: nog                                ! Output grid size
   integer :: noga                               ! Output grid size, halo A
   integer :: nc0b                               ! Subset Sc0 size, halo B
   integer,allocatable :: og_to_proc(:)          ! Output grid to processor
   integer,allocatable :: proc_to_noga(:)        ! Number of points in output grid, halo A, for each processor
   integer,allocatable :: oga_to_og(:)           ! Output grid, halo A to global
   integer,allocatable :: og_to_oga(:)           ! Output grid, global to halo A
   integer,allocatable :: c0b_to_c0(:)           ! Subset Sc0, halo B to global
   integer,allocatable :: c0_to_c0b(:)           ! Subset Sc0, global to halo B
   integer,allocatable :: c0a_to_c0b(:)          ! Subset Sc0, halo A to halo B
   logical,allocatable :: mask(:,:)              ! Mask on output grid
   type(linop_type) :: og                        ! Subset Sc0 to grid interpolation
   type(com_type) :: com_AB                      ! Communication between halos A and B
   integer :: nlonio                             ! I/O chunk size on output grid
   integer :: ilon_s                             ! Starting index for the I/O chunk on output grid
   integer,allocatable :: ogio_to_og(:)          ! Output grid, I/O chunk to global
   integer,allocatable :: procio_to_proc_grid(:) ! I/O task to main communicator task
   type(com_type) :: com_AIO_grid                ! Communication between halo A and I/O chunk
contains
   procedure :: dealloc => io_dealloc
   procedure :: fld_read => io_fld_read
   procedure :: fld_write => io_fld_write
   procedure :: fld_write_grid => io_fld_write_grid
   procedure :: init => io_init
   procedure :: init_grid => io_init_grid
end type io_type

private
public :: io_type

contains

!----------------------------------------------------------------------
! Subroutine: io_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine io_dealloc(io)

implicit none

! Passed variables
class(io_type),intent(inout) :: io ! I/O

! Release memory
if (allocated(io%c0io_to_c0)) deallocate(io%c0io_to_c0)
if (allocated(io%procio_to_proc)) deallocate(io%procio_to_proc)
call io%com_AIO%dealloc
if (allocated(io%lon)) deallocate(io%lon)
if (allocated(io%lat)) deallocate(io%lat)
if (allocated(io%og_to_lon)) deallocate(io%og_to_lon)
if (allocated(io%og_to_lat)) deallocate(io%og_to_lat)
if (allocated(io%og_to_proc)) deallocate(io%og_to_proc)
if (allocated(io%proc_to_noga)) deallocate(io%proc_to_noga)
if (allocated(io%oga_to_og)) deallocate(io%oga_to_og)
if (allocated(io%og_to_oga)) deallocate(io%og_to_oga)
if (allocated(io%c0b_to_c0)) deallocate(io%c0b_to_c0)
if (allocated(io%c0_to_c0b)) deallocate(io%c0_to_c0b)
if (allocated(io%c0a_to_c0b)) deallocate(io%c0a_to_c0b)
call io%og%dealloc
call io%com_AB%dealloc
if (allocated(io%mask)) deallocate(io%mask)
if (allocated(io%ogio_to_og)) deallocate(io%ogio_to_og)
if (allocated(io%procio_to_proc_grid)) deallocate(io%procio_to_proc_grid)
call io%com_AIO_grid%dealloc

end subroutine io_dealloc

!----------------------------------------------------------------------
! Subroutine: io_fld_read
! Purpose: write field
!----------------------------------------------------------------------
subroutine io_fld_read(io,mpl,nam,geom,filename,varname,fld)

implicit none

! Passed variables
class(io_type),intent(in) :: io                        ! I/O
type(mpl_type),intent(inout) :: mpl                    ! MPI data
type(nam_type),intent(in) :: nam                       ! Namelist
type(geom_type),intent(in) :: geom                     ! Geometry
character(len=*),intent(in) :: filename                ! File name
character(len=*),intent(in) :: varname                 ! Variable name
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) ! Field

! Local variables
integer :: ncid,fld_id
real(kind_real),allocatable :: fld_c0io(:,:)
character(len=1024),parameter :: subr = 'io_fld_read'

! Allocation
allocate(fld_c0io(io%nc0io,geom%nl0))

if (any(io%procio_to_proc==mpl%myproc).and.(io%nc0io>0)) then
   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))

   ! Get variable id
   call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(varname),fld_id))

   ! Get data
   call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_c0io,(/io%ic0_s,1/),(/io%nc0io,geom%nl0/)))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))
end if

! Communication
call io%com_AIO%red(mpl,geom%nl0,fld_c0io,fld)

! Release memory
deallocate(fld_c0io)

end subroutine io_fld_read

!----------------------------------------------------------------------
! Subroutine: io_fld_write
! Purpose: write field
!----------------------------------------------------------------------
subroutine io_fld_write(io,mpl,nam,geom,filename,varname,fld)

implicit none

! Passed variables
class(io_type),intent(in) :: io                       ! I/O
type(mpl_type),intent(inout) :: mpl                   ! MPI data
type(nam_type),intent(in) :: nam                      ! Namelist
type(geom_type),intent(in) :: geom                    ! Geometry
character(len=*),intent(in) :: filename               ! File name
character(len=*),intent(in) :: varname                ! Variable name
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) ! Field

! Local variables
integer :: ic0a,il0,info,info_coord,color
integer :: ncid,nc0_id,nl0_id,fld_id,lon_id,lat_id
real(kind_real) :: fld_c0a(geom%nc0a,geom%nl0)
real(kind_real),allocatable :: fld_c0io(:,:),lon_c0io(:),lat_c0io(:)
character(len=1024) :: cname
character(len=1024),parameter :: subr = 'io_fld_write'
type(fckit_mpi_comm) :: f_comm

! Apply mask
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      if (geom%mask_c0a(ic0a,il0)) then
         fld_c0a(ic0a,il0) = fld(ic0a,il0)
      else
         fld_c0a(ic0a,il0) = mpl%msv%valr
      end if
   end do
end do

! Allocation
allocate(fld_c0io(io%nc0io,geom%nl0))
allocate(lon_c0io(io%nc0io))
allocate(lat_c0io(io%nc0io))

! Communication
call io%com_AIO%ext(mpl,geom%nl0,fld_c0a,fld_c0io)
call io%com_AIO%ext(mpl,geom%lon_c0a,lon_c0io)
call io%com_AIO%ext(mpl,geom%lat_c0a,lat_c0io)

! Create communicator for parallel NetCDF I/O
if (any(io%procio_to_proc==mpl%myproc).and.(io%nc0io>0)) then
   color = 1
   cname = trim(mpl%f_comm%name())//'_'//trim(nam%prefix)//'_io'
else
   color = 0
   cname = trim(mpl%f_comm%name())//'_'//trim(nam%prefix)//'_no_io'
endif
f_comm = mpl%f_comm%split(color,cname)

! Parallel I/O
if (any(io%procio_to_proc==mpl%myproc).and.(io%nc0io>0)) then
   ! Check if the file exists
   info = nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',ior(nf90_noclobber,ior(nf90_netcdf4,nf90_mpiio)),ncid, &
        & comm=f_comm%communicator(),info=f_comm%info_null())
   if (info/=nf90_noerr) then
      ! Open file
      call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc', &
    & ior(nf90_write,ior(nf90_netcdf4,nf90_mpiio)),ncid,comm=f_comm%communicator(),info=f_comm%info_null()))

      ! Enter definition mode
      call mpl%ncerr(subr,nf90_redef(ncid))
   end if

   ! Create dimension
   nc0_id = mpl%ncdimcheck(subr,ncid,'nc0',geom%nc0,.true.)
   nl0_id = mpl%ncdimcheck(subr,ncid,'nl0',geom%nl0,.true.,.true.)

   ! Define coordinates if necessary
   info_coord = nf90_inq_varid(ncid,'lon',lon_id)
   if (info_coord/=nf90_noerr) then
      call mpl%ncerr(subr,nf90_def_var(ncid,'lon',nc_kind_real,(/nc0_id/),lon_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',mpl%msv%valr))
      call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'unit','degrees_north'))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lat',nc_kind_real,(/nc0_id/),lat_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',mpl%msv%valr))
      call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'unit','degrees_east'))
   end if

   ! Define variable if necessary
   info = nf90_inq_varid(ncid,trim(varname),fld_id)
   if (info/=nf90_noerr) then
      call mpl%ncerr(subr,nf90_def_var(ncid,trim(varname),nc_kind_real,(/nc0_id,nl0_id/),fld_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',mpl%msv%valr))
   end if

   ! End definition mode
   call mpl%ncerr(subr,nf90_enddef(ncid))

   ! Write coordinates if necessary
   if (info_coord/=nf90_noerr) then
      ! Convert to degrees
      lon_c0io = lon_c0io*rad2deg
      lat_c0io = lat_c0io*rad2deg

      ! Write data
      call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,lon_c0io,(/io%ic0_s/),(/io%nc0io/)))
      call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,lat_c0io,(/io%ic0_s/),(/io%nc0io/)))
   end if

   ! Write variable
   call mpl%ncerr(subr,nf90_put_var(ncid,fld_id,fld_c0io,(/io%ic0_s,1/),(/io%nc0io,geom%nl0/)))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))
end if

! Write global attributes
if (io%procio_to_proc(1)==mpl%myproc) then
   ! Reopen file with one task only
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_write,ncid))

   ! Write namelist parameters
   call nam%write(mpl,ncid)

   ! Define attribute
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',mpl%msv%valr))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))
end if

! Release memory
deallocate(fld_c0io)
deallocate(lon_c0io)
deallocate(lat_c0io)

call f_comm%delete()

! Gridded field output
if (nam%grid_output) call io%fld_write_grid(mpl,nam,geom,trim(filename)//'_gridded',varname,fld_c0a)

! Wait for everybody
call mpl%f_comm%barrier()

end subroutine io_fld_write

!----------------------------------------------------------------------
! Subroutine: io_fld_write_grid
! Purpose: interpolate and write gridded field
!----------------------------------------------------------------------
subroutine io_fld_write_grid(io,mpl,nam,geom,filename,varname,fld)

implicit none

! Passed variables
class(io_type),intent(in) :: io                       ! I/O
type(mpl_type),intent(inout) :: mpl                   ! MPI data
type(nam_type),intent(in) :: nam                      ! Namelist
type(geom_type),intent(in) :: geom                    ! Geometry
character(len=*),intent(in) :: filename               ! File name
character(len=*),intent(in) :: varname                ! Variable name
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) ! Field

! Local variables
integer :: il0,info,info_coord,ilonio,ilat,iogio,ioga,color
integer :: ncid,nlon_id,nlat_id,nlev_id,fld_id,lon_id,lat_id,lev_id
real(kind_real) :: fld_c0b(io%nc0b,geom%nl0)
real(kind_real) :: fld_oga(io%noga,geom%nl0)
real(kind_real),allocatable :: fld_ogio(:,:),fld_loniolat(:,:,:),lon(:),lat(:)
character(len=1024) :: cname
character(len=1024),parameter :: subr = 'io_fld_write_grid'
type(fckit_mpi_comm) :: f_comm

! Halo extension and interpolation
call io%com_AB%ext(mpl,geom%nl0,fld,fld_c0b)
do il0=1,geom%nl0
   call io%og%apply(mpl,fld_c0b(:,il0),fld_oga(:,il0),mssrc=.true.)
end do

! Apply mask
do il0=1,geom%nl0
   do ioga=1,io%noga
      if (.not.io%mask(ioga,il0)) fld_oga(ioga,il0) = mpl%msv%valr
   end do
end do

! Allocation
allocate(fld_ogio(io%nlonio*io%nlat,geom%nl0))
allocate(fld_loniolat(io%nlonio,io%nlat,geom%nl0))

! Communication
call io%com_AIO_grid%ext(mpl,geom%nl0,fld_oga,fld_ogio)

! Reshape
do ilonio=1,io%nlonio
   do ilat=1,io%nlat
      iogio = (ilat-1)*io%nlonio+ilonio
      fld_loniolat(ilonio,ilat,:) = fld_ogio(iogio,:)
   end do
end do

! Create communicator for parallel NetCDF I/O
if (any(io%procio_to_proc_grid==mpl%myproc).and.(io%nlonio>0)) then
   color = 1
   cname = trim(mpl%f_comm%name())//'_'//trim(nam%prefix)//'_io_grid'
else
   color = 0
   cname = trim(mpl%f_comm%name())//'_'//trim(nam%prefix)//'_no_io_grid'
endif
f_comm = mpl%f_comm%split(color,cname)

! Parallel I/O
if (any(io%procio_to_proc_grid==mpl%myproc).and.(io%nlonio>0)) then
   ! Check if the file exists
   info = nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',ior(nf90_noclobber,ior(nf90_netcdf4,nf90_mpiio)),ncid, &
        & comm=f_comm%communicator(),info=f_comm%info_null())
   if (info/=nf90_noerr) then
      ! Open file
      call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc', &
    & ior(nf90_write,ior(nf90_netcdf4,nf90_mpiio)),ncid,comm=f_comm%communicator(),info=f_comm%info_null()))

      ! Enter definition mode
      call mpl%ncerr(subr,nf90_redef(ncid))
   end if

   ! Define dimensions
   nlon_id = mpl%ncdimcheck(subr,ncid,'nlon',io%nlon,.true.)
   nlat_id = mpl%ncdimcheck(subr,ncid,'nlat',io%nlat,.true.)
   nlev_id = mpl%ncdimcheck(subr,ncid,'nlev',geom%nl0,.true.,.true.)

   ! Define coordinates if necessary
   info_coord = nf90_inq_varid(ncid,'lon',lon_id)
   if (info_coord/=nf90_noerr) then
      call mpl%ncerr(subr,nf90_def_var(ncid,'lon',nc_kind_real,(/nlon_id/),lon_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',mpl%msv%valr))
      call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'unit','degrees_north'))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lat',nc_kind_real,(/nlat_id/),lat_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',mpl%msv%valr))
      call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'unit','degrees_east'))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lev',nc_kind_real,(/nlev_id/),lev_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,lev_id,'_FillValue',mpl%msv%valr))
      call mpl%ncerr(subr,nf90_put_att(ncid,lev_id,'unit','layer'))
   end if

   ! Define variable if necessary
   info = nf90_inq_varid(ncid,trim(varname),fld_id)
   if (info/=nf90_noerr) then
      call mpl%ncerr(subr,nf90_def_var(ncid,trim(varname),nc_kind_real,(/nlon_id,nlat_id,nlev_id/),fld_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',mpl%msv%valr))
   end if

   ! End definition mode
   call mpl%ncerr(subr,nf90_enddef(ncid))

   ! Write coordinates if necessary
   if (info_coord/=nf90_noerr) then
      ! Allocation
      allocate(lon(io%nlon))
      allocate(lat(io%nlat))

      ! Convert to degrees
      lon = io%lon*rad2deg
      lat = io%lat*rad2deg

      ! Write data
      call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,lon))
      call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,lat))
      do il0=1,geom%nl0
         call mpl%ncerr(subr,nf90_put_var(ncid,lev_id,real(il0,kind_real),(/il0/)))
      end do

      ! Release memory
      deallocate(lon)
      deallocate(lat)
   end if

   ! Write variable
   call mpl%ncerr(subr,nf90_put_var(ncid,fld_id,fld_loniolat,(/io%ilon_s,1,1/),(/io%nlonio,io%nlat,geom%nl0/)))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))
end if

! Release memory
deallocate(fld_ogio)
deallocate(fld_loniolat)
call f_comm%delete()

end subroutine io_fld_write_grid

!----------------------------------------------------------------------
! Subroutine: io_init
! Purpose: initialize fields output
!----------------------------------------------------------------------
subroutine io_init(io,mpl,rng,nam,geom)

implicit none

! Passed variables
class(io_type),intent(inout) :: io  ! I/O
type(mpl_type),intent(inout) :: mpl ! MPI data
type(rng_type),intent(inout) :: rng ! Random number generator
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry

! Local variables
integer :: nres,iprocio,delta,ic0io,ic0,ic0_s,ic0_e,nc0own,ic0own,iproc,jproc
integer,allocatable :: procio_to_nc0io(:),list_proc(:),order_proc(:),order_c0(:),c0own_to_c0io(:)
real(kind_real),allocatable :: list_c0(:)
logical,allocatable :: proc_isio(:)

! Allocation
allocate(procio_to_nc0io(nam%nprocio))
allocate(list_proc(mpl%nproc))
allocate(order_proc(mpl%nproc))
allocate(proc_isio(mpl%nproc))
allocate(io%procio_to_proc(nam%nprocio))
allocate(list_c0(geom%nc0))
allocate(order_c0(geom%nc0))

! Initialization
nres = geom%nc0
proc_isio = .false.
io%nc0io = 0

! Use hash order to order points
list_c0 = geom%hash_c0
call qsort(geom%nc0,list_c0,order_c0)

! Define Sc0 subset splitting among I/O processors
do iprocio=1,nam%nprocio
   ! Define contiguous chunk
   delta = geom%nc0/nam%nprocio
   if (nres>(nam%nprocio-iprocio+1)*delta) delta = delta+1
   procio_to_nc0io(iprocio) = delta
   nres = nres-delta
   if (iprocio==1) then
      ic0_s = 1
   else
      ic0_s = sum(procio_to_nc0io(1:iprocio-1))+1
   end if
   ic0_e = sum(procio_to_nc0io(1:iprocio))

   ! Order processors given their implication in this chunk
   do iproc=1,mpl%nproc
      list_proc(iproc) = count(geom%c0_to_proc(order_c0(ic0_s:ic0_e))==iproc)
   end do
   call qsort(mpl%nproc,list_proc,order_proc)

   ! Select I/O processor for this chunk
   do iproc=mpl%nproc,1,-1
      jproc = order_proc(iproc)
      if (.not.proc_isio(jproc)) then
         proc_isio(jproc) = .true.
         io%procio_to_proc(iprocio) = jproc
         if (mpl%myproc==jproc) then
            io%nc0io = procio_to_nc0io(iprocio)
            io%ic0_s = ic0_s
         end if
         exit
      end if
   end do
end do

! Allocation
allocate(io%c0io_to_c0(io%nc0io))

! Initialization
iprocio = 1
iproc = io%procio_to_proc(iprocio)
ic0io = 0
nc0own = 0

! Go through Sc0 subset distribution, first pass
do ic0=1,geom%nc0
   ic0io = ic0io+1
   if (ic0io>procio_to_nc0io(iprocio)) then
      iprocio = iprocio+1
      iproc = io%procio_to_proc(iprocio)
      ic0io = 1
   end if
   if (mpl%myproc==iproc) then
      io%c0io_to_c0(ic0io) = ic0
      if (geom%c0_to_proc(order_c0(ic0))==iproc) nc0own = nc0own+1
   end if
end do

! Allocation
allocate(c0own_to_c0io(nc0own))

! Initialization
iprocio = 1
iproc = io%procio_to_proc(iprocio)
ic0io = 0
ic0own = 0

! Go through Sc0 subset distribution, second pass
do ic0=1,geom%nc0
   ic0io = ic0io+1
   if (ic0io>procio_to_nc0io(iprocio)) then
      iprocio = iprocio+1
      iproc = io%procio_to_proc(iprocio)
      ic0io = 1
   end if
   if (mpl%myproc==iproc) then
      if (geom%c0_to_proc(order_c0(ic0))==iproc) then
         ic0own = ic0own+1
         c0own_to_c0io(ic0own) = ic0io
      end if
   end if
end do

! Print results
write(mpl%info,'(a7,a)') '','I/O splitting:'
call mpl%flush
do iprocio=1,nam%nprocio
   write(mpl%info,'(a10,a,i4,a,i8,a,i4)') '','Chunk ',iprocio,' ~> ',procio_to_nc0io(iprocio),' grid-points handled by task ', &
 & io%procio_to_proc(iprocio)
   call mpl%flush
end do

! Setup Sc0 subset I/O communication
call io%com_AIO%setup(mpl,'com_AIO',geom%nc0,geom%nc0a,io%nc0io,nc0own,io%c0io_to_c0,c0own_to_c0io,geom%c0_to_proc(order_c0), &
 & geom%c0_to_c0a(order_c0))

! Release memory
deallocate(list_proc)
deallocate(order_proc)
deallocate(proc_isio)
deallocate(procio_to_nc0io)
deallocate(list_c0)
deallocate(order_c0)
deallocate(c0own_to_c0io)

if (nam%grid_output) call io%init_grid(mpl,rng,nam,geom)

end subroutine io_init

!----------------------------------------------------------------------
! Subroutine: io_init_grid
! Purpose: initialize fields gridding
!----------------------------------------------------------------------
subroutine io_init_grid(io,mpl,rng,nam,geom)

implicit none

! Passed variables
class(io_type),intent(inout) :: io  ! I/O
type(mpl_type),intent(inout) :: mpl ! MPI data
type(rng_type),intent(inout) :: rng ! Random number generator
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry

! Local variables
integer :: ilon,ilat,i_s,iog,ic0,ic0a,ic0b,ioga,il0,nn_index(1)
integer :: nres,iprocio,delta,ilonio,ilon_s,ilon_e,nogown,iogown,iproc,jproc,iogio
integer,allocatable :: procio_to_nlonio(:),list(:),order(:),ogown_to_ogio(:)
real(kind_real) :: dlon,dlat
real(kind_real),allocatable :: lon_og(:),lat_og(:),lon_oga(:),lat_oga(:)
logical :: mask_c0(geom%nc0)
logical,allocatable :: mask_oga(:),lcheck_c0b(:),proc_isio(:)
character(len=1024),parameter :: subr = 'io_init_grid'

! Grid size
io%nlat = nint(pi/nam%grid_resol)
io%nlon = 2*io%nlat
dlon = 2.0*pi/real(io%nlon,kind_real)
dlat = pi/real(io%nlat,kind_real)

! Print results
write(mpl%info,'(a7,a)') '','Output grid:'
call mpl%flush
write(mpl%info,'(a10,a,f7.2,a,f5.2,a)') '','Effective resolution: ',0.5*(dlon+dlat)*reqkm,' km (', &
 & 0.5*(dlon+dlat)*rad2deg,' deg.)'
call mpl%flush
write(mpl%info,'(a10,a,i4,a,i4)') '',      'Size (nlon x nlat):   ',io%nlon,' x ',io%nlat
call mpl%flush

! Allocation
allocate(io%lon(io%nlon))
allocate(io%lat(io%nlat))
io%nog = io%nlon*io%nlat
allocate(io%og_to_lon(io%nog))
allocate(io%og_to_lat(io%nog))
allocate(lon_og(io%nog))
allocate(lat_og(io%nog))
allocate(io%og_to_proc(io%nog))
allocate(io%proc_to_noga(mpl%nproc))
allocate(io%og_to_oga(io%nog))

! Define lon/lat
do ilat=1,io%nlat
   do ilon=1,io%nlon
      io%lon(ilon) = (-pi+dlon/2)+real(ilon-1,kind_real)*dlon
      io%lat(ilat) = (-pi/2+dlat/2)+real(ilat-1,kind_real)*dlat
   end do
end do

! Setup packed output grid
do ilat=1,io%nlat
   do ilon=1,io%nlon
      iog = (ilat-1)*io%nlon+ilon
      lon_og(iog) = io%lon(ilon)
      lat_og(iog) = io%lat(ilat)
      io%og_to_lon(iog) = ilon
      io%og_to_lat(iog) = ilat
      call geom%tree%find_nearest_neighbors(lon_og(iog),lat_og(iog),1,nn_index)
      io%og_to_proc(iog) = geom%c0_to_proc(nn_index(1))
   end do
end do
if (mpl%msv%isany(io%og_to_proc)) call mpl%abort(subr,'some output grid points do not have a processorw')
do iproc=1,mpl%nproc
   io%proc_to_noga(iproc) = count(io%og_to_proc==iproc)
end do
io%noga = io%proc_to_noga(mpl%myproc)

! Allocation
allocate(io%oga_to_og(io%noga))
allocate(lon_oga(io%noga))
allocate(lat_oga(io%noga))
allocate(mask_oga(io%noga))

! Global-local conversions for halo A
ioga = 0
do iog=1,io%nog
   iproc = io%og_to_proc(iog)
   if (iproc==mpl%myproc) then
      ioga = ioga+1
      io%oga_to_og(ioga) = iog
   end if
end do
call mpl%glb_to_loc_index(io%noga,io%oga_to_og,io%nog,io%og_to_oga)

! Compute interpolation
lon_oga = lon_og(io%oga_to_og)
lat_oga = lat_og(io%oga_to_og)
do ioga=1,io%noga
   ! Indices
   iog = io%oga_to_og(ioga)
   ilon = io%og_to_lon(iog)
   ilat = io%og_to_lat(iog)

   ! Check that the interpolation point is inside the domain
   call geom%mesh%inside(mpl,io%lon(ilon),io%lat(ilat),mask_oga(ioga))

   ! Check poles
   if (abs(io%lat(ilat))>maxval(abs(geom%lat_c0))) mask_oga(ioga) = .false.
end do
mask_c0 = .true.
write(io%og%prefix,'(a,i3.3)') 'og'
call io%og%interp(mpl,rng,nam,geom,0,geom%nc0,geom%lon_c0,geom%lat_c0,mask_c0,io%noga,lon_oga,lat_oga,mask_oga,10)

! Allocation
allocate(lcheck_c0b(geom%nc0))

! Define halo B
lcheck_c0b = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   lcheck_c0b(ic0) = .true.
end do
do i_s=1,io%og%n_s
   ic0 = io%og%col(i_s)
   lcheck_c0b(ic0) = .true.
end do
io%nc0b = count(lcheck_c0b)

! Allocation
allocate(io%c0b_to_c0(io%nc0b))
allocate(io%c0_to_c0b(geom%nc0))
allocate(io%c0a_to_c0b(geom%nc0a))

! Global-local conversions for halo B
ic0b = 0
do ic0=1,geom%nc0
   if (lcheck_c0b(ic0)) then
      ic0b = ic0b+1
      io%c0b_to_c0(ic0b) = ic0
      io%c0_to_c0b(ic0) = ic0b
   end if
end do

! Halos A-B conversion
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0b = io%c0_to_c0b(ic0)
   io%c0a_to_c0b(ic0a) = ic0b
end do

! Local interpolation source
io%og%n_src = io%nc0b
do i_s=1,io%og%n_s
   io%og%col(i_s) = io%c0_to_c0b(io%og%col(i_s))
end do

! Setup communications
call io%com_AB%setup(mpl,'com_AB',geom%nc0,geom%nc0a,io%nc0b,geom%nc0a,io%c0b_to_c0,io%c0a_to_c0b,geom%c0_to_proc,geom%c0_to_c0a)

! Compute output grid mask
allocate(io%mask(io%noga,geom%nl0))
io%mask = .true.
do i_s=1,io%og%n_s
   ic0b = io%og%col(i_s)
   ioga = io%og%row(i_s)
   ic0 = io%c0b_to_c0(ic0b)
   do il0=1,geom%nl0
      if (.not.geom%mask_c0(ic0,il0)) io%mask(ioga,il0) = .false.
   end do
end do

! I/O splitting

! Allocation
allocate(procio_to_nlonio(nam%nprocio))
allocate(list(mpl%nproc))
allocate(order(mpl%nproc))
allocate(proc_isio(mpl%nproc))
allocate(io%procio_to_proc_grid(nam%nprocio))

! Initialization
nres = io%nlon
proc_isio = .false.
io%nlonio = 0

! Define output grid splitting among I/O processors
do iprocio=1,nam%nprocio
   ! Define contiguous chunk
   delta = io%nlon/nam%nprocio
   if (nres>(nam%nprocio-iprocio+1)*delta) delta = delta+1
   procio_to_nlonio(iprocio) = delta
   nres = nres-delta
   if (iprocio==1) then
      ilon_s = 1
   else
      ilon_s = sum(procio_to_nlonio(1:iprocio-1))+1
   end if
   ilon_e = sum(procio_to_nlonio(1:iprocio))

   ! Order processors given their implication in this chunk
   list = 0
   do ilat=1,io%nlat
      do ilon=ilon_s,ilon_e
         iog = (ilat-1)*io%nlon+ilon
         iproc = io%og_to_proc(iog)
         list(iproc) = list(iproc)+1
      end do
   end do
   call qsort(mpl%nproc,list,order)

   ! Select I/O processor for this chunk
   do iproc=mpl%nproc,1,-1
      jproc = order(iproc)
      if (.not.proc_isio(jproc)) then
         proc_isio(jproc) = .true.
         io%procio_to_proc_grid(iprocio) = jproc
         if (mpl%myproc==jproc) then
            io%nlonio = procio_to_nlonio(iprocio)
            io%ilon_s = ilon_s
         end if
         exit
      end if
   end do
end do

! Allocation
allocate(io%ogio_to_og(io%nlonio*io%nlat))

! Initialization
iprocio = 1
iproc = io%procio_to_proc_grid(iprocio)
ilonio = 0
nogown = 0

! Go through output grid distribution, first pass
do ilon=1,io%nlon
   ilonio = ilonio+1
   if (ilonio>procio_to_nlonio(iprocio)) then
      iprocio = iprocio+1
      iproc = io%procio_to_proc_grid(iprocio)
      ilonio = 1
   end if
   if (mpl%myproc==iproc) then
      do ilat=1,io%nlat
         iogio = (ilat-1)*io%nlonio+ilonio
         iog = (ilat-1)*io%nlon+ilon
         io%ogio_to_og(iogio) = iog
         if (io%og_to_proc(iog)==iproc) nogown = nogown+1
      end do
   end if
end do

! Allocation
allocate(ogown_to_ogio(nogown))

! Initialization
iprocio = 1
iproc = io%procio_to_proc_grid(iprocio)
ilonio = 0
iogown = 0

! Go through output grid distribution, second pass
do ilon=1,io%nlon
   ilonio = ilonio+1
   if (ilonio>procio_to_nlonio(iprocio)) then
      iprocio = iprocio+1
      iproc = io%procio_to_proc_grid(iprocio)
      ilonio = 1
   end if
   if (mpl%myproc==iproc) then
      do ilat=1,io%nlat
         iogio = (ilat-1)*io%nlonio+ilonio
         iog = (ilat-1)*io%nlon+ilon
         if (io%og_to_proc(iog)==iproc) then
            iogown = iogown+1
            ogown_to_ogio(iogown) = iogio
         end if
      end do
   end if
end do

! Print results
write(mpl%info,'(a7,a)') '','I/O splitting:'
call mpl%flush
do iprocio=1,nam%nprocio
   write(mpl%info,'(a10,a,i4,a,i8,a,i4)') '','Chunk ',iprocio,' ~> ',procio_to_nlonio(iprocio),' longitudes handled by task ', &
 & io%procio_to_proc_grid(iprocio)
   call mpl%flush
end do

! Setup Sc0 subset I/O communication
call io%com_AIO_grid%setup(mpl,'com_AIO_grid',io%nog,io%noga,io%nlonio*io%nlat,nogown,io%ogio_to_og,ogown_to_ogio,io%og_to_proc, &
 & io%og_to_oga)

! Release memory
deallocate(procio_to_nlonio)
deallocate(list)
deallocate(order)
deallocate(proc_isio)
deallocate(ogown_to_ogio)

! Print results
write(mpl%info,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0 =    ',geom%nc0
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0a =   ',geom%nc0a
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nc0b =   ',io%nc0b
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','nog =    ',io%nog
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','noga =   ',io%noga
call mpl%flush
write(mpl%info,'(a10,a,i8)') '','og%n_s = ',io%og%n_s
call mpl%flush

! Release memory
deallocate(lon_og)
deallocate(lat_og)
deallocate(lon_oga)
deallocate(lat_oga)
deallocate(mask_oga)
deallocate(lcheck_c0b)

end subroutine io_init_grid

end module type_io
