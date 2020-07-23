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
   integer,allocatable :: c0u_to_c0b(:)          ! Subset Sc0, universe to halo B
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
   procedure :: init => io_init
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
if (allocated(io%c0u_to_c0b)) deallocate(io%c0u_to_c0b)
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
subroutine io_fld_read(io,mpl,nam,geom,filename,variable,fld,groupname)

implicit none

! Passed variables
class(io_type),intent(in) :: io                        ! I/O
type(mpl_type),intent(inout) :: mpl                    ! MPI data
type(nam_type),intent(in) :: nam                       ! Namelist
type(geom_type),intent(in) :: geom                     ! Geometry
character(len=*),intent(in) :: filename                ! File name
character(len=*),intent(in) :: variable                ! Variable name
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) ! Field
character(len=*),intent(in),optional :: groupname      ! Group name

! Local variables
integer :: ncid,grpid,fld_id
real(kind_real),allocatable :: fld_c0io(:,:)
character(len=1024),parameter :: subr = 'io_fld_read'

! Allocation
allocate(fld_c0io(io%nc0io,geom%nl0))

if (any(io%procio_to_proc==mpl%myproc).and.(io%nc0io>0)) then
   ! Open file
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))

   ! Get group
   if (present(groupname)) then
      call mpl%ncerr(subr,nf90_inq_grp_ncid(ncid,trim(groupname),grpid))
   else
      grpid = ncid
   end if

   ! Check dimension
   call mpl%nc_dim_check(subr,ncid,'nc0',geom%nc0)
   call mpl%nc_dim_check(subr,ncid,'nl0',geom%nl0)

   ! Get variable
   call mpl%ncerr(subr,nf90_inq_varid(grpid,trim(variable),fld_id))

   ! Get data
   call mpl%ncerr(subr,nf90_get_var(grpid,fld_id,fld_c0io,(/io%ic0_s,1/),(/io%nc0io,geom%nl0/)))

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
subroutine io_fld_write(io,mpl,nam,geom,filename,variable,fld,groupname,subgroupname)

implicit none

! Passed variables
class(io_type),intent(in) :: io                       ! I/O
type(mpl_type),intent(inout) :: mpl                   ! MPI data
type(nam_type),intent(in) :: nam                      ! Namelist
type(geom_type),intent(in) :: geom                    ! Geometry
character(len=*),intent(in) :: filename               ! File name
character(len=*),intent(in) :: variable               ! Variable name
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) ! Field
character(len=*),intent(in),optional :: groupname     ! Group name
character(len=*),intent(in),optional :: subgroupname  ! Subgroup name

! Local variables
integer :: ic0a,il0,info,color
integer :: ncid,grpid,subgrpid,nc0_id,nl0_id,fld_id,lon_id,lat_id
real(kind_real) :: fld_c0a(geom%nc0a,geom%nl0)
real(kind_real),allocatable :: fld_c0io(:,:),lon_c0io(:),lat_c0io(:)
character(len=1024) :: cname
character(len=1024),parameter :: subr = 'io_fld_write'
type(fckit_mpi_comm) :: f_comm

! Apply mask
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      if (geom%gmask_c0a(ic0a,il0)) then
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
   ! Define file
   ncid = mpl%nc_file_create_or_open(subr,trim(nam%datadir)//'/'//trim(filename)//'.nc',f_comm)

   ! Define group
   if (present(groupname)) then
      grpid = mpl%nc_group_define_or_get(subr,ncid,trim(groupname))

      ! Define sub-group
      if (present(subgroupname)) then
         subgrpid = mpl%nc_group_define_or_get(subr,grpid,trim(subgroupname))
      else
         subgrpid = grpid
      end if
   else
      subgrpid = ncid
   end if

   ! Define dimensions
   nc0_id = mpl%nc_dim_define_or_get(subr,ncid,'nc0',geom%nc0)
   nl0_id = mpl%nc_dim_define_or_get(subr,ncid,'nl0',geom%nl0)

   ! Define coordinates
   lon_id = mpl%nc_var_define_or_get(subr,ncid,'lon',nc_kind_real,(/nc0_id/),'degrees_north')
   lat_id = mpl%nc_var_define_or_get(subr,ncid,'lat',nc_kind_real,(/nc0_id/),'degrees_east')

   ! Define variable
   fld_id = mpl%nc_var_define_or_get(subr,subgrpid,trim(variable),nc_kind_real,(/nc0_id,nl0_id/))

   ! Convert to degrees
   lon_c0io = lon_c0io*rad2deg
   lat_c0io = lat_c0io*rad2deg

   ! Write coordinates
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,lon_c0io,(/io%ic0_s/),(/io%nc0io/)))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,lat_c0io,(/io%ic0_s/),(/io%nc0io/)))

   ! Write variable
   call mpl%ncerr(subr,nf90_put_var(subgrpid,fld_id,fld_c0io,(/io%ic0_s,1/),(/io%nc0io,geom%nl0/)))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))
end if

! Write global attributes
if (io%procio_to_proc(1)==mpl%myproc) then
   ! Reopen file with one task only
   call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_write,ncid))

   ! Write namelist parameters
   call nam%write(mpl,ncid)

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))
end if

! Release memory
deallocate(fld_c0io)
deallocate(lon_c0io)
deallocate(lat_c0io)

! Delete communicator
call f_comm%delete()

! Wait for everybody
call mpl%f_comm%barrier()

end subroutine io_fld_write

!----------------------------------------------------------------------
! Subroutine: io_init
! Purpose: initialize fields output
!----------------------------------------------------------------------
subroutine io_init(io,mpl,nam,geom)

implicit none

! Passed variables
class(io_type),intent(inout) :: io  ! I/O
type(mpl_type),intent(inout) :: mpl ! MPI data
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry

! Local variables
integer :: nres,iprocio,delta,ic0io,ic0,ic0_s,ic0_e,jc0,nc0own,ic0own,iproc,jproc
integer,allocatable :: procio_to_nc0io(:),list_proc(:),order_proc(:),order_c0(:),c0own_to_c0(:)
real(kind_real),allocatable :: hash_c0(:)
logical,allocatable :: proc_isio(:)

! Allocation
allocate(procio_to_nc0io(nam%nprocio))
allocate(list_proc(mpl%nproc))
allocate(order_proc(mpl%nproc))
allocate(proc_isio(mpl%nproc))
allocate(io%procio_to_proc(nam%nprocio))
if (nam%repro) then
   allocate(hash_c0(geom%nc0))
   allocate(order_c0(geom%nc0))
end if

! Initialization
nres = geom%nc0
proc_isio = .false.
io%nc0io = 0

if (nam%repro) then
   ! Communication
   call mpl%loc_to_glb(geom%nc0a,geom%nc0,geom%c0a_to_c0,geom%hash_c0a,hash_c0,.true.)

   ! Use hash order to order points
   call qsort(geom%nc0,hash_c0,order_c0)
end if

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
   list_proc = 0
   do iproc=1,mpl%nproc
      do ic0=ic0_s,ic0_e
         if (nam%repro) then
            jc0 = order_c0(ic0)
         else
            jc0 = ic0
         end if
         if (geom%c0_to_proc(jc0)==iproc) list_proc(iproc) = list_proc(iproc)+1
      end do
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
   if (nam%repro) then
      jc0 = order_c0(ic0)
   else
      jc0 = ic0
   end if
   if (iproc==mpl%myproc) then
      io%c0io_to_c0(ic0io) = jc0
      if (geom%c0_to_proc(jc0)==iproc) nc0own = nc0own+1
   end if
end do

! Allocation
allocate(c0own_to_c0(nc0own))

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
   if (nam%repro) then
      jc0 = order_c0(ic0)
   else
      jc0 = ic0
   end if
   if (iproc==mpl%myproc) then
      if (geom%c0_to_proc(jc0)==iproc) then
         ic0own = ic0own+1
         c0own_to_c0(ic0own) = jc0
      end if
   end if
end do

! Print results
write(mpl%info,'(a7,a)') '','I/O splitting:'
call mpl%flush
do iprocio=1,nam%nprocio
   write(mpl%info,'(a10,a,i6,a,i8,a,i6)') '','Chunk ',iprocio,' ~> ',procio_to_nc0io(iprocio),' grid-points handled by task ', &
 & io%procio_to_proc(iprocio)
   call mpl%flush
end do

! Setup Sc0 subset I/O communication
call io%com_AIO%setup(mpl,'com_AIO',geom%nc0a,io%nc0io,geom%nc0,geom%c0a_to_c0,io%c0io_to_c0,c0own_to_c0)

! Release memory
deallocate(list_proc)
deallocate(order_proc)
deallocate(proc_isio)
deallocate(procio_to_nc0io)
if (nam%repro) then
   deallocate(hash_c0)
   deallocate(order_c0)
end if
deallocate(c0own_to_c0)

end subroutine io_init

end module type_io
