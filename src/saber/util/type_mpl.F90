!----------------------------------------------------------------------
! Module: type_mpl
! Purpose: MPI parameters derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_mpl

use iso_fortran_env, only : output_unit
use fckit_mpi_module, only: fckit_mpi_comm,fckit_mpi_sum,fckit_mpi_status
use fckit_log_module, only: fckit_log
use netcdf
!$ use omp_lib
use tools_kinds, only: kind_real,nc_kind_real
use type_msv, only: msv_type

implicit none

integer,parameter :: lunit_min=10   ! Minimum unit number
integer,parameter :: lunit_max=1000 ! Maximum unit number
integer,parameter :: ddis = 5       ! Progression display step

type mpl_type
   ! MPI communicator
   type(fckit_mpi_comm) :: f_comm   ! MPI data
   integer :: nproc                 ! Number of MPI tasks
   integer :: myproc                ! MPI task index
   integer :: rootproc              ! Main task index
   logical :: main                  ! Main task logical
   integer :: tag                   ! MPI tag

   ! Number of OpenMP threads
   integer :: nthread               ! Number of OpenMP threads

   ! Missing values
   type(msv_type) :: msv            ! Missing values

   ! Display parameters
   character(len=1024) :: verbosity ! Verbosity level
   character(len=1024) :: info      ! Info buffer
   character(len=1024) :: trace     ! Trace buffer
   character(len=1024) :: stats     ! Stats buffer
   character(len=1024) :: test      ! Test buffer
   integer :: lunit                 ! Listing unit

   ! Display colors
   character(len=1024) :: black     ! Black color code
   character(len=1024) :: green     ! Green color code
   character(len=1024) :: peach     ! Peach color code
   character(len=1024) :: aqua      ! Aqua color code
   character(len=1024) :: purple    ! Purple color code
   character(len=1024) :: err       ! Error color code
   character(len=1024) :: wng       ! Warning color code

   ! Progression print
   integer :: nprog                 ! Progression array size
   integer :: progint               ! Progression integer
   logical,allocatable :: done(:)   ! Progression array
contains
   procedure :: newunit => mpl_newunit
   procedure :: init => mpl_init
   procedure :: final => mpl_final
   procedure :: flush => mpl_flush
   procedure :: abort => mpl_abort
   procedure :: warning => mpl_warning
   procedure :: update_tag => mpl_update_tag
   procedure :: broadcast => mpl_broadcast_string_1d
   procedure :: mpl_dot_prod_1d
   procedure :: mpl_dot_prod_2d
   procedure :: mpl_dot_prod_3d
   procedure :: mpl_dot_prod_4d
   generic :: dot_prod => mpl_dot_prod_1d,mpl_dot_prod_2d,mpl_dot_prod_3d,mpl_dot_prod_4d
   procedure :: glb_to_loc_index => mpl_glb_to_loc_index
   procedure :: mpl_glb_to_loc_integer_1d
   procedure :: mpl_glb_to_loc_integer_2d
   procedure :: mpl_glb_to_loc_real_1d
   procedure :: mpl_glb_to_loc_real_2d
   procedure :: mpl_glb_to_loc_logical_1d
   procedure :: mpl_glb_to_loc_logical_2d
   generic :: glb_to_loc => mpl_glb_to_loc_integer_1d,mpl_glb_to_loc_integer_2d,mpl_glb_to_loc_real_1d,mpl_glb_to_loc_real_2d, &
                          & mpl_glb_to_loc_logical_1d,mpl_glb_to_loc_logical_2d
   procedure :: mpl_loc_to_glb_integer_1d
   procedure :: mpl_loc_to_glb_integer_2d
   procedure :: mpl_loc_to_glb_real_1d
   procedure :: mpl_loc_to_glb_real_2d
   procedure :: mpl_loc_to_glb_logical_1d
   procedure :: mpl_loc_to_glb_logical_2d
   generic :: loc_to_glb => mpl_loc_to_glb_integer_1d,mpl_loc_to_glb_integer_2d,mpl_loc_to_glb_real_1d,mpl_loc_to_glb_real_2d, &
                          & mpl_loc_to_glb_logical_1d,mpl_loc_to_glb_logical_2d
   procedure :: prog_init => mpl_prog_init
   procedure :: prog_print => mpl_prog_print
   procedure :: prog_final => mpl_prog_final
   procedure :: nc_file_create_or_open => mpl_nc_file_create_or_open
   procedure :: nc_group_define_or_get => mpl_nc_group_define_or_get
   procedure :: nc_dim_define_or_get => mpl_nc_dim_define_or_get
   procedure :: nc_dim_inquire => mpl_nc_dim_inquire
   procedure :: nc_dim_check => mpl_nc_dim_check
   procedure :: nc_var_define_or_get => mpl_nc_var_define_or_get
   procedure :: ncerr => mpl_ncerr
   procedure :: mpl_write_integer
   procedure :: mpl_write_integer_array
   procedure :: mpl_write_real
   procedure :: mpl_write_real_array
   procedure :: mpl_write_logical
   procedure :: mpl_write_logical_array
   procedure :: mpl_write_string
   procedure :: mpl_write_string_array
   generic :: write => mpl_write_integer,mpl_write_integer_array,mpl_write_real,mpl_write_real_array, &
            & mpl_write_logical,mpl_write_logical_array,mpl_write_string,mpl_write_string_array
end type mpl_type

private
public :: mpl_type

contains

!----------------------------------------------------------------------
! Subroutine: mpl_newunit
! Purpose: find a free unit
!----------------------------------------------------------------------
subroutine mpl_newunit(mpl,lunit)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(out) :: lunit         ! New unit

! Local variables
integer :: lun
logical :: lopened
character(len=1024),parameter :: subr = 'mpl_newunit'

! Loop over possible units
do lun=lunit_min,lunit_max
   inquire(unit=lun,opened=lopened)
   if (.not.lopened) then
      lunit=lun
      exit
   end if
end do

! Check
if (lopened) call mpl%abort(subr,'cannot find a free unit')

end subroutine mpl_newunit

!----------------------------------------------------------------------
! Subroutine: mpl_init
! Purpose: initialize MPL object
!----------------------------------------------------------------------
subroutine mpl_init(mpl,f_comm)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl             ! MPI data
type(fckit_mpi_comm),intent(in),target :: f_comm ! FCKIT MPI communicator wrapper

! Get BUMP own FCKIT MPI communicator wrapper
mpl%f_comm = f_comm

! Get MPI size
mpl%nproc = mpl%f_comm%size()

! Get MPI rank
mpl%myproc = mpl%f_comm%rank()+1

! Define main task
mpl%rootproc = 1
mpl%main = (mpl%myproc==mpl%rootproc)

! Time-based tag
if (mpl%main) then
   call system_clock(count=mpl%tag)
end if
call mpl%f_comm%broadcast(mpl%tag,mpl%rootproc-1)
call mpl%update_tag(0)

! Set max number of OpenMP threads
mpl%nthread = 1
!$ mpl%nthread = omp_get_max_threads()
!$ call omp_set_num_threads(mpl%nthread)

! Set log at 'no message'
mpl%info = 'no_message'
mpl%trace = 'no_message'
mpl%stats = 'no_message'
mpl%test = 'no_message'

! Set default listing
mpl%black = ' '
mpl%green = ' '
mpl%peach = ' '
mpl%aqua = ' '
mpl%purple = ' '
mpl%err = ' '
mpl%wng = ' '
mpl%verbosity = 'all'
mpl%lunit = mpl%msv%vali

end subroutine mpl_init

!----------------------------------------------------------------------
! Subroutine: mpl_final
! Purpose: finalize MPI
!----------------------------------------------------------------------
subroutine mpl_final(mpl)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data

! Release memory
if (allocated(mpl%done)) deallocate(mpl%done)

end subroutine mpl_final

!----------------------------------------------------------------------
! Subroutine: mpl_flush
! Purpose: flush listings
!----------------------------------------------------------------------
subroutine mpl_flush(mpl,advance_flag)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl        ! MPI data
logical,intent(in),optional :: advance_flag ! Advance flag

! Local variables
logical :: ladvance_flag

if ((trim(mpl%verbosity)=='all').or.((trim(mpl%verbosity)=='main').and.mpl%main)) then
   ! Set advance flag
   ladvance_flag = .true.
   if (present(advance_flag)) ladvance_flag = advance_flag

   ! Check info message
   if (trim(mpl%info)/='no_message') then
      ! Write message
      if (ladvance_flag) then
         if (mpl%msv%is(mpl%lunit)) then
            call fckit_log%info(mpl%info)
         else
            write(mpl%lunit,'(a)') trim(mpl%info)
            call flush(mpl%lunit)
         end if
      else
         if (mpl%msv%is(mpl%lunit)) then
            call fckit_log%info(mpl%info,newl=.false.)
         else
            write(mpl%lunit,'(a)',advance='no') trim(mpl%info)
            call flush(mpl%lunit)
         end if
      end if

      ! Set at 'no message'
      mpl%info = 'no_message'
   end if

   ! Check trace message
   if (trim(mpl%trace)/='no_message') then
      ! Write message
      if (ladvance_flag) then
         if (mpl%msv%is(mpl%lunit)) then
            call fckit_log%info('OOPS_TRACE: '//trim(mpl%trace))
         else
            write(mpl%lunit,'(a)') 'OOPS_TRACE: '//trim(mpl%trace)
            call flush(mpl%lunit)
         end if
      else
         if (mpl%msv%is(mpl%lunit)) then
            call fckit_log%info('OOPS_TRACE: '//trim(mpl%trace),newl=.false.)
         else
            write(mpl%lunit,'(a)',advance='no') 'OOPS_TRACE: '//trim(mpl%trace)
            call flush(mpl%lunit)
         end if
      end if

      ! Set at 'no message'
      mpl%trace = 'no_message'
   end if

   ! Check stats message
   if (trim(mpl%stats)/='no_message') then
      ! Write message
      if (ladvance_flag) then
         if (mpl%msv%is(mpl%lunit)) then
            call fckit_log%info('OOPS_STATS: '//trim(mpl%stats))
         else
            write(mpl%lunit,'(a)') 'OOPS_STATS: '//trim(mpl%stats)
            call flush(mpl%lunit)
         end if
      else
         if (mpl%msv%is(mpl%lunit)) then
            call fckit_log%info('OOPS_STATS: '//trim(mpl%stats),newl=.false.)
         else
            write(mpl%lunit,'(a)',advance='no') 'OOPS_STATS: '//trim(mpl%stats)
            call flush(mpl%lunit)
         end if
      end if

      ! Set at 'no message'
      mpl%stats = 'no_message'
   end if

   ! Check test message
   if (trim(mpl%test)/='no_message') then
      ! Write message
      if (ladvance_flag) then
         if (mpl%msv%is(mpl%lunit)) then
            call fckit_log%info('Test     : '//trim(mpl%test))
         else
            write(mpl%lunit,'(a)') 'Test     : '//trim(mpl%test)
            call flush(mpl%lunit)
         end if
      else
         if (mpl%msv%is(mpl%lunit)) then
            call fckit_log%info('Test     : '//trim(mpl%test),newl=.false.)
         else
            write(mpl%lunit,'(a)',advance='no') 'Test     : '//trim(mpl%test)
            call flush(mpl%lunit)
         end if
      end if

      ! Set at 'no message'
      mpl%test = 'no_message'
   end if
end if

end subroutine mpl_flush

!----------------------------------------------------------------------
! Subroutine: mpl_abort
! Purpose: clean MPI abort
!----------------------------------------------------------------------
subroutine mpl_abort(mpl,subr,message)

implicit none

! Passed variable
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: subr    ! Calling subroutine
character(len=*),intent(in) :: message ! Message

! Write standard output message
write(output_unit,'(a,i4.4,a)') '!!! ABORT in '//trim(subr)//' on task #',mpl%myproc,': '//trim(message)
call flush(output_unit)

! Abort with proper MPI communicator
call mpl%f_comm%abort(1)

end subroutine mpl_abort

!----------------------------------------------------------------------
! Subroutine: mpl_warning
! Purpose: print warning message
!----------------------------------------------------------------------
subroutine mpl_warning(mpl,subr,message)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: subr    ! Calling subroutine
character(len=*),intent(in) :: message ! Message

! Print warning message
write(mpl%info,'(a)') trim(mpl%wng)//'!!! Warning in '//trim(subr)//': '//trim(message)//trim(mpl%black)
call mpl%flush

end subroutine mpl_warning


!----------------------------------------------------------------------
! Subroutine: mpl_update_tag
! Purpose: update MPI tag
!----------------------------------------------------------------------
subroutine mpl_update_tag(mpl,add)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: add              ! Tag update incrememnt

! Update tag
mpl%tag = mpl%tag+add

! Apply bounds (between 1 and 10000)
mpl%tag = mod(mpl%tag,10000)
mpl%tag = max(mpl%tag,1)

end subroutine mpl_update_tag

!----------------------------------------------------------------------
! Subroutine: mpl_broadcast_string_1d
! Purpose: broadcast 1d string array
!----------------------------------------------------------------------
subroutine mpl_broadcast_string_1d(mpl,var,root)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl                  ! MPI data
character(len=*),dimension(:),intent(inout) :: var ! Logical array, 1d
integer,intent(in) :: root                         ! Root task

! Local variable
integer :: i

! Broadcast one string at a time
do i=1,size(var)
   call mpl%f_comm%broadcast(var(i),root)
end do

end subroutine mpl_broadcast_string_1d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_1d
! Purpose: global dot product over local fields, 1d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_1d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl     ! MPI data
real(kind_real),intent(in) :: fld1(:) ! Field 1
real(kind_real),intent(in) :: fld2(:) ! Field 2
real(kind_real),intent(out) :: dp     ! Global dot product

! Local variable
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = 0.0
dp_loc(1) = sum(fld1*fld2,mask=mpl%msv%isnot(fld1).and.mpl%msv%isnot(fld2))

! Allreduce
dp_out = 0.0
call mpl%f_comm%allreduce(dp_loc,dp_out,fckit_mpi_sum())
dp = dp_out(1)

! Broadcast
call mpl%f_comm%broadcast(dp,mpl%rootproc-1)

end subroutine mpl_dot_prod_1d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_2d
! Purpose: global dot product over local fields, 2d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_2d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl       ! MPI data
real(kind_real),intent(in) :: fld1(:,:) ! Field 1
real(kind_real),intent(in) :: fld2(:,:) ! Field 2
real(kind_real),intent(out) :: dp       ! Global dot product

! Local variable
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = 0.0
dp_loc(1) = sum(fld1*fld2,mask=mpl%msv%isnot(fld1).and.mpl%msv%isnot(fld2))

! Allreduce
dp_out = 0.0
call mpl%f_comm%allreduce(dp_loc,dp_out,fckit_mpi_sum())
dp = dp_out(1)

! Broadcast
call mpl%f_comm%broadcast(dp,mpl%rootproc-1)

end subroutine mpl_dot_prod_2d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_3d
! Purpose: global dot product over local fields, 3d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_3d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl         ! MPI data
real(kind_real),intent(in) :: fld1(:,:,:) ! Field 1
real(kind_real),intent(in) :: fld2(:,:,:) ! Field 2
real(kind_real),intent(out) :: dp         ! Global dot product

! Local variable
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = 0.0
dp_loc(1) = sum(fld1*fld2,mask=mpl%msv%isnot(fld1).and.mpl%msv%isnot(fld2))

! Allreduce
dp_out = 0.0
call mpl%f_comm%allreduce(dp_loc,dp_out,fckit_mpi_sum())
dp = dp_out(1)

! Broadcast
call mpl%f_comm%broadcast(dp,mpl%rootproc-1)

end subroutine mpl_dot_prod_3d

!----------------------------------------------------------------------
! Subroutine: mpl_dot_prod_4d
! Purpose: global dot product over local fields, 4d
!----------------------------------------------------------------------
subroutine mpl_dot_prod_4d(mpl,fld1,fld2,dp)

implicit none

! Passed variables
class(mpl_type),intent(in) :: mpl           ! MPI data
real(kind_real),intent(in) :: fld1(:,:,:,:) ! Field 1
real(kind_real),intent(in) :: fld2(:,:,:,:) ! Field 2
real(kind_real),intent(out) :: dp           ! Global dot product

! Local variable
real(kind_real) :: dp_loc(1),dp_out(1)

! Product and sum
dp_loc(1) = 0.0
dp_loc(1) = sum(fld1*fld2,mask=mpl%msv%isnot(fld1).and.mpl%msv%isnot(fld2))

! Allreduce
dp_out = 0.0
call mpl%f_comm%allreduce(dp_loc,dp_out,fckit_mpi_sum())
dp = dp_out(1)

! Broadcast
call mpl%f_comm%broadcast(dp,mpl%rootproc-1)

end subroutine mpl_dot_prod_4d

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_index
! Purpose: communicate global index to local index
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_index(mpl,n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc,rootproc,pool)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl           ! MPI data
integer,intent(in) :: n_loc                    ! Local dimension
integer,intent(in) :: loc_to_glb(n_loc)        ! Local to global index
integer,intent(in) :: n_glb                    ! Global dimension
integer,intent(out) :: glb_to_loc(:)           ! Global to local index
integer,intent(out) :: glb_to_proc(:)          ! Global to processor
integer,intent(in),optional :: rootproc        ! Root task
logical,intent(in),optional :: pool(mpl%nproc) ! Tasks pool

! Local variables
integer :: iproc,i_loc,n_loc_tmp,lrootproc
integer,allocatable :: loc_to_glb_tmp(:)
logical :: lpool(mpl%nproc)
character(len=1024),parameter :: subr = 'mpl_glb_to_loc_index'
type(fckit_mpi_status) :: status

! Get local rootproc and pool
lrootproc = mpl%rootproc
if (present(rootproc)) lrootproc = rootproc
lpool = .true.
if (present(pool)) lpool = pool

! Check global array size
if (mpl%myproc==lrootproc) then
   if (.not.lpool(lrootproc)) call mpl%abort(subr,'root task should be in the tasks pool')
end if

if (mpl%myproc==lrootproc) then
   ! Check global array size
   if (size(glb_to_loc)/=n_glb) call mpl%abort(subr,'wrong dimension for the glb_to_loc')
   if (size(glb_to_proc)/=n_glb) call mpl%abort(subr,'wrong dimension for the glb_to_proc')

   ! Initialization
   glb_to_loc = mpl%msv%vali
   glb_to_proc = mpl%msv%vali

   do iproc=1,mpl%nproc
      if (lpool(iproc)) then
         if (iproc==lrootproc) then
            ! Copy dimension
            n_loc_tmp = n_loc
         else
            ! Receive dimension on rootproc
            call mpl%f_comm%receive(n_loc_tmp,iproc-1,mpl%tag,status)
         end if

         ! Allocation
         allocate(loc_to_glb_tmp(n_loc_tmp))

         if (iproc==lrootproc) then
            ! Copy data
            loc_to_glb_tmp = loc_to_glb
         else
            ! Receive data on rootproc
            call mpl%f_comm%receive(loc_to_glb_tmp,iproc-1,mpl%tag+1,status)
         end if

         ! Fill glb_to_loc and glb_to_proc if required
         do i_loc=1,n_loc_tmp
            glb_to_loc(loc_to_glb_tmp(i_loc)) = i_loc
            glb_to_proc(loc_to_glb_tmp(i_loc)) = iproc
         end do

         ! Release memory
         deallocate(loc_to_glb_tmp)
      end if
   end do
else
   if (lpool(mpl%myproc)) then
      ! Send dimensions to rootproc
      call mpl%f_comm%send(n_loc,lrootproc-1,mpl%tag)

      ! Send data to rootproc
      call mpl%f_comm%send(loc_to_glb,lrootproc-1,mpl%tag+1)
   end if
end if
call mpl%update_tag(2)

end subroutine mpl_glb_to_loc_index

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_integer_1d
! Purpose: global to local, 1d array
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_integer_1d(mpl,n_loc,n_glb,loc_to_glb,glb,loc)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl    ! MPI data
integer,intent(in) :: n_loc             ! Local array size
integer,intent(in) :: n_glb             ! Global array size
integer,intent(in) :: loc_to_glb(n_loc) ! Local to global
integer,intent(in) :: glb(:)            ! Global array
integer,intent(out) :: loc(n_loc)       ! Local array

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp
integer,allocatable :: glb_to_loc(:),glb_to_proc(:),sbuf(:)
character(len=1024),parameter :: subr = 'mpl_glb_to_loc_integer_1d'
type(fckit_mpi_status) :: status

! Check global array size
if (mpl%main) then
   if (size(glb)/=n_glb) call mpl%abort(subr,'wrong dimension for the global array in mpl_glb_to_loc_integer_1d')
end if

! Allocation
if (mpl%main) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(sbuf(n_loc_tmp))

      ! Prepare buffers
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            sbuf(i_loc) = glb(i_glb)
         end if
      end do

      if (iproc==mpl%rootproc) then
         ! Copy data
         loc = sbuf
      else
         ! Send data to iproc
         call mpl%f_comm%send(sbuf,iproc-1,mpl%tag)
      end if

      ! Release memory
      deallocate(sbuf)
   end do
else
   ! Receive data from rootproc
   call mpl%f_comm%receive(loc,mpl%rootproc-1,mpl%tag,status)
end if
call mpl%update_tag(1)

! Release memory
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_glb_to_loc_integer_1d

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_integer_2d
! Purpose: global to local, 2d array
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_integer_2d(mpl,nl,n_loc,n_glb,loc_to_glb,glb,loc,rootproc,pool)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl           ! MPI data
integer,intent(in) :: nl                       ! Number of levels
integer,intent(in) :: n_loc                    ! Local array size
integer,intent(in) :: n_glb                    ! Global array size
integer,intent(in) :: loc_to_glb(n_loc)        ! Local to global
integer,intent(in) :: glb(:,:)                 ! Global array
integer,intent(out) :: loc(n_loc,nl)           ! Local array
integer,intent(in),optional :: rootproc        ! Root task
logical,intent(in),optional :: pool(mpl%nproc) ! Tasks pool

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp,il,lrootproc
integer,allocatable :: glb_to_loc(:),glb_to_proc(:),sbuf(:),rbuf(:)
logical :: lpool(mpl%nproc)
character(len=1024),parameter :: subr = 'mpl_glb_to_loc_integer_2d'
type(fckit_mpi_status) :: status

! Get local rootproc and pool
lrootproc = mpl%rootproc
if (present(rootproc)) lrootproc = rootproc
lpool = .true.
if (present(pool)) lpool = pool

! Check global array size
if (mpl%myproc==lrootproc) then
   if (.not.lpool(lrootproc)) call mpl%abort(subr,'root task should be in the tasks pool')
   if (size(glb,1)/=n_glb) call mpl%abort(subr,'wrong first dimension for the global array in mpl_glb_to_loc_integer_2d')
   if (size(glb,2)/=nl) call mpl%abort(subr,'wrong second dimension for the global array in mpl_glb_to_loc_integer_2d')
end if

! Allocation
if (mpl%myproc==lrootproc) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

if (lpool(mpl%myproc)) then
   ! Allocation
   allocate(rbuf(n_loc*nl))
end if

if (mpl%myproc==lrootproc) then
   do iproc=1,mpl%nproc
      if (lpool(iproc)) then
         ! Allocation 
         n_loc_tmp = count(glb_to_proc==iproc)
         allocate(sbuf(n_loc_tmp*nl))

         ! Prepare buffers
         do i_glb=1,n_glb
            jproc = glb_to_proc(i_glb)
            if (iproc==jproc) then
               i_loc = glb_to_loc(i_glb)
               do il=1,nl
                  sbuf((il-1)*n_loc_tmp+i_loc) = glb(i_glb,il)
               end do
            end if
         end do

         if (iproc==lrootproc) then
            ! Copy data
            rbuf = sbuf
         else
            ! Send data to iproc
            call mpl%f_comm%send(sbuf,iproc-1,mpl%tag)
         end if

         ! Release memory
         deallocate(sbuf)
      end if
   end do
else
   if (lpool(mpl%myproc)) then
      ! Receive data from rootproc
      call mpl%f_comm%receive(rbuf,lrootproc-1,mpl%tag,status)
   end if
end if
call mpl%update_tag(1)

if (lpool(mpl%myproc)) then
   ! Unpack buffer
   do il=1,nl
      do i_loc=1,n_loc
         loc(i_loc,il) = rbuf((il-1)*n_loc+i_loc)
      end do
   end do
end if

! Release memory
deallocate(rbuf)
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_glb_to_loc_integer_2d

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_real_1d
! Purpose: global to local, 1d array
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_real_1d(mpl,n_loc,n_glb,loc_to_glb,glb,loc)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl      ! MPI data
integer,intent(in) :: n_loc               ! Local array size
integer,intent(in) :: n_glb               ! Global array size
integer,intent(in) :: loc_to_glb(n_loc)   ! Local to global
real(kind_real),intent(in) :: glb(:)      ! Global array
real(kind_real),intent(out) :: loc(n_loc) ! Local array

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp
integer,allocatable :: glb_to_loc(:),glb_to_proc(:)
real(kind_real),allocatable :: sbuf(:)
character(len=1024),parameter :: subr = 'mpl_glb_to_loc_real_1d'
type(fckit_mpi_status) :: status

! Check global array size
if (mpl%main) then
   if (size(glb)/=n_glb) call mpl%abort(subr,'wrong dimension for the global array in mpl_glb_to_loc_real_1d')
end if

! Allocation
if (mpl%main) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(sbuf(n_loc_tmp))

      ! Prepare buffers
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            sbuf(i_loc) = glb(i_glb)
         end if
      end do

      if (iproc==mpl%rootproc) then
         ! Copy data
         loc = sbuf
      else
         ! Send data to iproc
         call mpl%f_comm%send(sbuf,iproc-1,mpl%tag)
      end if

      ! Release memory
      deallocate(sbuf)
   end do
else
   ! Receive data from rootproc
   call mpl%f_comm%receive(loc,mpl%rootproc-1,mpl%tag,status)
end if
call mpl%update_tag(1)

! Release memory
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_glb_to_loc_real_1d

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_real_2d
! Purpose: global to local, 2d array
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_real_2d(mpl,nl,n_loc,n_glb,loc_to_glb,glb,loc,rootproc,pool)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl           ! MPI data
integer,intent(in) :: nl                       ! Number of levels
integer,intent(in) :: n_loc                    ! Local array size
integer,intent(in) :: n_glb                    ! Global array size
integer,intent(in) :: loc_to_glb(n_loc)        ! Local to global
real(kind_real),intent(in) :: glb(:,:)         ! Global array
real(kind_real),intent(out) :: loc(n_loc,nl)   ! Local array
integer,intent(in),optional :: rootproc        ! Root task
logical,intent(in),optional :: pool(mpl%nproc) ! Tasks pool

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp,il,lrootproc
integer,allocatable :: glb_to_loc(:),glb_to_proc(:)
real(kind_real),allocatable :: sbuf(:),rbuf(:)
logical :: lpool(mpl%nproc)
character(len=1024),parameter :: subr = 'mpl_glb_to_loc_real_2d'
type(fckit_mpi_status) :: status

! Get local rootproc and pool
lrootproc = mpl%rootproc
if (present(rootproc)) lrootproc = rootproc
lpool = .true.
if (present(pool)) lpool = pool

! Check global array size
if (mpl%myproc==lrootproc) then
   if (.not.lpool(lrootproc)) call mpl%abort(subr,'root task should be in the tasks pool')
   if (size(glb,1)/=n_glb) call mpl%abort(subr,'wrong first dimension for the global array in mpl_glb_to_loc_real_2d')
   if (size(glb,2)/=nl) call mpl%abort(subr,'wrong second dimension for the global array in mpl_glb_to_loc_real_2d')
end if

! Allocation
if (mpl%myproc==lrootproc) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc,rootproc,pool)

if (lpool(mpl%myproc)) then
   ! Allocation
   allocate(rbuf(n_loc*nl))
end if

if (mpl%myproc==lrootproc) then
   do iproc=1,mpl%nproc
      if (lpool(iproc)) then
         ! Allocation 
         n_loc_tmp = count(glb_to_proc==iproc)
         allocate(sbuf(n_loc_tmp*nl))

         ! Prepare buffers
         do i_glb=1,n_glb
            jproc = glb_to_proc(i_glb)
            if (iproc==jproc) then
               i_loc = glb_to_loc(i_glb)
               do il=1,nl
                  sbuf((il-1)*n_loc_tmp+i_loc) = glb(i_glb,il)
               end do
            end if
         end do

         if (iproc==lrootproc) then
            ! Copy data
            rbuf = sbuf
         else
            ! Send data to iproc
            call mpl%f_comm%send(sbuf,iproc-1,mpl%tag)
         end if

         ! Release memory
         deallocate(sbuf)
      end if
   end do
else
   if (lpool(mpl%myproc)) then
      ! Receive data from rootproc
      call mpl%f_comm%receive(rbuf,lrootproc-1,mpl%tag,status)
   end if
end if
call mpl%update_tag(1)

if (lpool(mpl%myproc)) then
   ! Unpack buffer
   do il=1,nl
      do i_loc=1,n_loc
         loc(i_loc,il) = rbuf((il-1)*n_loc+i_loc)
      end do
   end do
end if

! Release memory
deallocate(rbuf)
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_glb_to_loc_real_2d

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_logical_1d
! Purpose: global to local, 1d array
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_logical_1d(mpl,n_loc,n_glb,loc_to_glb,glb,loc)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl     ! MPI data
integer,intent(in) :: n_loc             ! Local array size
integer,intent(in) :: n_glb             ! Global array size
integer,intent(in) :: loc_to_glb(n_loc) ! Local to global
logical,intent(in) :: glb(:)            ! Global array
logical,intent(out) :: loc(n_loc)       ! Local array

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp
integer,allocatable :: glb_to_loc(:),glb_to_proc(:)
logical,allocatable :: sbuf(:)
character(len=1024),parameter :: subr = 'mpl_glb_to_loc_logical_1d'
type(fckit_mpi_status) :: status

! Check global array size
if (mpl%main) then
   if (size(glb)/=n_glb) call mpl%abort(subr,'wrong dimension for the global array in mpl_glb_to_loc_logical_1d')
end if

! Allocation
if (mpl%main) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(sbuf(n_loc_tmp))

      ! Prepare buffers
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            sbuf(i_loc) = glb(i_glb)
         end if
      end do

      if (iproc==mpl%rootproc) then
         ! Copy data
         loc = sbuf
      else
         ! Send data to iproc
         call mpl%f_comm%send(sbuf,iproc-1,mpl%tag)
      end if

      ! Release memory
      deallocate(sbuf)
   end do
else
   ! Receive data from rootproc
   call mpl%f_comm%receive(loc,mpl%rootproc-1,mpl%tag,status)
end if
call mpl%update_tag(1)

! Release memory
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_glb_to_loc_logical_1d

!----------------------------------------------------------------------
! Subroutine: mpl_glb_to_loc_logical_2d
! Purpose: global to local, 2d array
!----------------------------------------------------------------------
subroutine mpl_glb_to_loc_logical_2d(mpl,nl,n_loc,n_glb,loc_to_glb,glb,loc,rootproc,pool)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl           ! MPI data
integer,intent(in) :: nl                       ! Number of levels
integer,intent(in) :: n_loc                    ! Local array size
integer,intent(in) :: n_glb                    ! Global array size
integer,intent(in) :: loc_to_glb(n_loc)        ! Local to global
logical,intent(in) :: glb(:,:)                 ! Global array
logical,intent(out) :: loc(n_loc,nl)           ! Local array
integer,intent(in),optional :: rootproc        ! Root task
logical,intent(in),optional :: pool(mpl%nproc) ! Tasks pool

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp,il,lrootproc
integer,allocatable :: glb_to_loc(:),glb_to_proc(:)
logical,allocatable :: sbuf(:),rbuf(:)
logical :: lpool(mpl%nproc)
character(len=1024),parameter :: subr = 'mpl_glb_to_loc_logical_2d'
type(fckit_mpi_status) :: status

! Get local rootproc and pool
lrootproc = mpl%rootproc
if (present(rootproc)) lrootproc = rootproc
lpool = .true.
if (present(pool)) lpool = pool

! Check global array size
if (mpl%myproc==lrootproc) then
   if (.not.lpool(lrootproc)) call mpl%abort(subr,'root task should be in the tasks pool')
   if (size(glb,1)/=n_glb) call mpl%abort(subr,'wrong first dimension for the global array in mpl_glb_to_loc_logical_2d')
   if (size(glb,2)/=nl) call mpl%abort(subr,'wrong second dimension for the global array in mpl_glb_to_loc_logical_2d')
end if

! Allocation
if (mpl%myproc==lrootproc) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

if (lpool(mpl%myproc)) then
   ! Allocation
   allocate(rbuf(n_loc*nl))
end if

if (mpl%myproc==lrootproc) then
   do iproc=1,mpl%nproc
      if (lpool(iproc)) then
         ! Allocation 
         n_loc_tmp = count(glb_to_proc==iproc)
         allocate(sbuf(n_loc_tmp*nl))

         ! Prepare buffers
         do i_glb=1,n_glb
            jproc = glb_to_proc(i_glb)
            if (iproc==jproc) then
               i_loc = glb_to_loc(i_glb)
               do il=1,nl
                  sbuf((il-1)*n_loc_tmp+i_loc) = glb(i_glb,il)
               end do
            end if
         end do

         if (iproc==lrootproc) then
            ! Copy data
            rbuf = sbuf
         else
            ! Send data to iproc
            call mpl%f_comm%send(sbuf,iproc-1,mpl%tag)
         end if

         ! Release memory
         deallocate(sbuf)
      end if
   end do
else
   if (lpool(mpl%myproc)) then
      ! Receive data from rootproc
      call mpl%f_comm%receive(rbuf,lrootproc-1,mpl%tag,status)
   end if
end if
call mpl%update_tag(1)

if (lpool(mpl%myproc)) then
   ! Unpack buffer
   do il=1,nl
      do i_loc=1,n_loc
         loc(i_loc,il) = rbuf((il-1)*n_loc+i_loc)
      end do
   end do
end if

! Release memory
deallocate(rbuf)
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_glb_to_loc_logical_2d

!----------------------------------------------------------------------
! Subroutine: mpl_loc_to_glb_integer_1d
! Purpose: local to global, 1d array
!----------------------------------------------------------------------
subroutine mpl_loc_to_glb_integer_1d(mpl,n_loc,n_glb,loc_to_glb,loc,glb,bcast)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl    ! MPI data
integer,intent(in) :: n_loc             ! Local array size
integer,intent(in) :: n_glb             ! Global array size
integer,intent(in) :: loc_to_glb(n_loc) ! Local to global
integer,intent(in) :: loc(n_loc)        ! Local array
integer,intent(out) :: glb(:)           ! Global array
logical,intent(in),optional :: bcast    ! Broadcast option

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp
integer,allocatable :: glb_to_loc(:),glb_to_proc(:),rbuf(:)
logical :: lbcast
character(len=1024),parameter :: subr = 'mpl_loc_to_glb_integer_1d'
type(fckit_mpi_status) :: status

! Get local bcast
lbcast = .false.
if (present(bcast)) lbcast = bcast

! Check global array size
if (mpl%main.or.lbcast) then
   if (size(glb)/=n_glb) call mpl%abort(subr,'wrong dimension for the global array in mpl_loc_to_glb_integer_1d')
end if

! Allocation
if (mpl%main) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(rbuf(n_loc_tmp))

      if (iproc==mpl%rootproc) then
          ! Copy data
          rbuf = loc
      else
          ! Receive data from iproc
          call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
      end if

      ! Add data to glb
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            glb(i_glb) = rbuf(i_loc)
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Send data to rootproc
   call mpl%f_comm%send(loc,mpl%rootproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Broadcast
if (lbcast) call mpl%f_comm%broadcast(glb,mpl%rootproc-1)

! Release memory
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_loc_to_glb_integer_1d

!----------------------------------------------------------------------
! Subroutine: mpl_loc_to_glb_integer_2d
! Purpose: local to global, 2d array
!----------------------------------------------------------------------
subroutine mpl_loc_to_glb_integer_2d(mpl,nl,n_loc,n_glb,loc_to_glb,loc,glb,bcast)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl    ! MPI data
integer,intent(in) :: nl                ! Number of levels
integer,intent(in) :: n_loc             ! Local array size
integer,intent(in) :: n_glb             ! Global array size
integer,intent(in) :: loc_to_glb(n_loc) ! Local to global
integer,intent(in) :: loc(n_loc,nl)     ! Local array
integer,intent(out) :: glb(:,:)         ! Global array
logical,intent(in),optional :: bcast    ! Broadcast option

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp,il
integer,allocatable :: glb_to_loc(:),glb_to_proc(:),rbuf(:),sbuf(:)
logical :: lbcast
character(len=1024),parameter :: subr = 'mpl_loc_to_glb_integer_2d'
type(fckit_mpi_status) :: status

! Get local bcast
lbcast = .false.
if (present(bcast)) lbcast = bcast

! Check global array size
if (mpl%main.or.lbcast) then
   if (size(glb,1)/=n_glb) call mpl%abort(subr,'wrong first dimension for the global array in mpl_loc_to_glb_integer_2d')
   if (size(glb,2)/=nl) call mpl%abort(subr,'wrong second dimension for the global array in mpl_loc_to_glb_integer_1d')
end if

! Allocation
if (mpl%main) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

! Allocation
allocate(sbuf(n_loc*nl))

! Prepare buffer
do il=1,nl
   do i_loc=1,n_loc
      sbuf((il-1)*n_loc+i_loc) = loc(i_loc,il)
   end do
end do

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(rbuf(n_loc_tmp*nl))

      if (iproc==mpl%rootproc) then
          ! Copy data
          rbuf = sbuf
      else
          ! Receive data from iproc
          call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
      end if

      ! Add data to glb
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            do il=1,nl
               glb(i_glb,il) = rbuf((il-1)*n_loc_tmp+i_loc)
            end do
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Send data to rootproc
   call mpl%f_comm%send(sbuf,mpl%rootproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Broadcast
if (lbcast) call mpl%f_comm%broadcast(glb,mpl%rootproc-1)

! Release memory
deallocate(sbuf)
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_loc_to_glb_integer_2d

!----------------------------------------------------------------------
! Subroutine: mpl_loc_to_glb_real_1d
! Purpose: local to global, 1d array
!----------------------------------------------------------------------
subroutine mpl_loc_to_glb_real_1d(mpl,n_loc,n_glb,loc_to_glb,loc,glb,bcast)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl     ! MPI data
integer,intent(in) :: n_loc              ! Local array size
integer,intent(in) :: n_glb              ! Global array size
integer,intent(in) :: loc_to_glb(n_loc)  ! Local to global
real(kind_real),intent(in) :: loc(n_loc) ! Local array
real(kind_real),intent(out) :: glb(:)    ! Global array
logical,intent(in),optional :: bcast     ! Broadcast option

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp
integer,allocatable :: glb_to_loc(:),glb_to_proc(:)
real(kind_real),allocatable :: rbuf(:)
logical :: lbcast
character(len=1024),parameter :: subr = 'mpl_loc_to_glb_real_1d'
type(fckit_mpi_status) :: status

! Get local bcast
lbcast = .false.
if (present(bcast)) lbcast = bcast

! Check global array size
if (mpl%main.or.lbcast) then
   if (size(glb)/=n_glb) call mpl%abort(subr,'wrong dimension for the global array in mpl_loc_to_glb_real_1d')
end if

! Allocation
if (mpl%main) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(rbuf(n_loc_tmp))

      if (iproc==mpl%rootproc) then
          ! Copy data
          rbuf = loc
      else
          ! Receive data from iproc
          call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
      end if

      ! Add data to glb
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            glb(i_glb) = rbuf(i_loc)
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Send data to rootproc
   call mpl%f_comm%send(loc,mpl%rootproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Broadcast
if (lbcast) call mpl%f_comm%broadcast(glb,mpl%rootproc-1)

! Release memory
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_loc_to_glb_real_1d

!----------------------------------------------------------------------
! Subroutine: mpl_loc_to_glb_real_2d
! Purpose: local to global, 2d array
!----------------------------------------------------------------------
subroutine mpl_loc_to_glb_real_2d(mpl,nl,n_loc,n_glb,loc_to_glb,loc,glb,bcast)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl        ! MPI data
integer,intent(in) :: nl                    ! Number of levels
integer,intent(in) :: n_loc                 ! Local array size
integer,intent(in) :: n_glb                 ! Global array size
integer,intent(in) :: loc_to_glb(n_loc)     ! Local to global
real(kind_real),intent(in) :: loc(n_loc,nl) ! Local array
real(kind_real),intent(out) :: glb(:,:)     ! Global array
logical,intent(in),optional :: bcast        ! Broadcast option

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp,il
integer,allocatable :: glb_to_loc(:),glb_to_proc(:)
real(kind_real),allocatable :: rbuf(:),sbuf(:)
logical :: lbcast
character(len=1024),parameter :: subr = 'mpl_loc_to_glb_real_2d'
type(fckit_mpi_status) :: status

! Get local bcast
lbcast = .false.
if (present(bcast)) lbcast = bcast

! Check global array size
if (mpl%main.or.lbcast) then
   if (size(glb,1)/=n_glb) call mpl%abort(subr,'wrong first dimension for the global array in mpl_loc_to_glb_real_2d')
   if (size(glb,2)/=nl) call mpl%abort(subr,'wrong second dimension for the global array in mpl_loc_to_glb_real_2d')
end if

! Allocation
if (mpl%main) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

! Allocation
allocate(sbuf(n_loc*nl))

! Prepare buffer
do il=1,nl
   do i_loc=1,n_loc
      sbuf((il-1)*n_loc+i_loc) = loc(i_loc,il)
   end do
end do

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(rbuf(n_loc_tmp*nl))

      if (iproc==mpl%rootproc) then
          ! Copy data
          rbuf = sbuf
      else
          ! Receive data from iproc
          call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
      end if

      ! Add data to glb
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            do il=1,nl
               glb(i_glb,il) = rbuf((il-1)*n_loc_tmp+i_loc)
            end do
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Send data to rootproc
   call mpl%f_comm%send(sbuf,mpl%rootproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Broadcast
if (lbcast) call mpl%f_comm%broadcast(glb,mpl%rootproc-1)

! Release memory
deallocate(sbuf)
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_loc_to_glb_real_2d

!----------------------------------------------------------------------
! Subroutine: mpl_loc_to_glb_logical_1d
! Purpose: local to global, 1d array
!----------------------------------------------------------------------
subroutine mpl_loc_to_glb_logical_1d(mpl,n_loc,n_glb,loc_to_glb,loc,glb,bcast)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl     ! MPI data
integer,intent(in) :: n_loc             ! Local array size
integer,intent(in) :: n_glb             ! Global array size
integer,intent(in) :: loc_to_glb(n_loc) ! Local to global
logical,intent(in) :: loc(n_loc)        ! Local array
logical,intent(out) :: glb(:)           ! Global array
logical,intent(in),optional :: bcast    ! Broadcast option

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp
integer,allocatable :: glb_to_loc(:),glb_to_proc(:)
logical :: lbcast
logical,allocatable :: rbuf(:)
character(len=1024),parameter :: subr = 'mpl_loc_to_glb_logical_1d'
type(fckit_mpi_status) :: status

! Get local bcast
lbcast = .false.
if (present(bcast)) lbcast = bcast

! Check global array size
if (mpl%main.or.lbcast) then
   if (size(glb)/=n_glb) call mpl%abort(subr,'wrong dimension for the global array in mpl_loc_to_glb_logical_1d')
end if

! Allocation
if (mpl%main) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(rbuf(n_loc_tmp))

      if (iproc==mpl%rootproc) then
          ! Copy data
          rbuf = loc
      else
          ! Receive data from iproc
          call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
      end if

      ! Add data to glb
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            glb(i_glb) = rbuf(i_loc)
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Send data to rootproc
   call mpl%f_comm%send(loc,mpl%rootproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Broadcast
if (lbcast) call mpl%f_comm%broadcast(glb,mpl%rootproc-1)

! Release memory
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_loc_to_glb_logical_1d

!----------------------------------------------------------------------
! Subroutine: mpl_loc_to_glb_logical_2d
! Purpose: local to global for a logical, 2d array
!----------------------------------------------------------------------
subroutine mpl_loc_to_glb_logical_2d(mpl,nl,n_loc,n_glb,loc_to_glb,loc,glb,bcast)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl     ! MPI data
integer,intent(in) :: nl                ! Number of levels
integer,intent(in) :: n_loc             ! Local array size
integer,intent(in) :: n_glb             ! Global array size
integer,intent(in) :: loc_to_glb(n_loc) ! Local to global
logical,intent(in) :: loc(n_loc,nl)     ! Local array
logical,intent(out) :: glb(:,:)         ! Global array
logical,intent(in),optional :: bcast    ! Broadcast option

! Local variables
integer :: iproc,jproc,i_glb,i_loc,n_loc_tmp,il
integer,allocatable :: glb_to_loc(:),glb_to_proc(:)
logical :: lbcast
logical,allocatable :: rbuf(:),sbuf(:)
character(len=1024),parameter :: subr = 'mpl_loc_to_glb_logical_2d'
type(fckit_mpi_status) :: status

! Get local bcast
lbcast = .false.
if (present(bcast)) lbcast = bcast

! Check global array size
if (mpl%main.or.lbcast) then
   if (size(glb,1)/=n_glb) call mpl%abort(subr,'wrong first dimension for the global array in mpl_loc_to_glb_logical_2d')
   if (size(glb,2)/=nl) call mpl%abort(subr,'wrong second dimension for the global array in mpl_loc_to_glb_logical_2d')
end if

! Allocation
if (mpl%main) then
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
else
   allocate(glb_to_loc(0))
   allocate(glb_to_proc(0))
end if

! Get global index and processor
call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)

! Allocation
allocate(sbuf(n_loc*nl))

! Prepare buffer
do il=1,nl
   do i_loc=1,n_loc
      sbuf((il-1)*n_loc+i_loc) = loc(i_loc,il)
   end do
end do

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Allocation
      n_loc_tmp = count(glb_to_proc==iproc)
      allocate(rbuf(n_loc_tmp*nl))

      if (iproc==mpl%rootproc) then
          ! Copy data
          rbuf = sbuf
      else
          ! Receive data from iproc
          call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
      end if

      ! Add data to glb
      do i_glb=1,n_glb
         jproc = glb_to_proc(i_glb)
         if (iproc==jproc) then
            i_loc = glb_to_loc(i_glb)
            do il=1,nl
               glb(i_glb,il) = rbuf((il-1)*n_loc_tmp+i_loc)
            end do
         end if
      end do

      ! Release memory
      deallocate(rbuf)
   end do
else
   ! Send data to rootproc
   call mpl%f_comm%send(sbuf,mpl%rootproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Broadcast
if (lbcast) call mpl%f_comm%broadcast(glb,mpl%rootproc-1)

! Release memory
deallocate(sbuf)
deallocate(glb_to_loc)
deallocate(glb_to_proc)

end subroutine mpl_loc_to_glb_logical_2d

!----------------------------------------------------------------------
! Subroutine: mpl_prog_init
! Purpose: initialize progression display
!----------------------------------------------------------------------
subroutine mpl_prog_init(mpl,nprog)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(in) :: nprog          ! Array size

! Local variables
integer :: ithread

! Print message
ithread = 0
!$ ithread = omp_get_thread_num()
if (ithread==0) then
   write(mpl%info,'(a)') ' 0%'
   call mpl%flush(.false.)
end if

! Allocation
allocate(mpl%done(nprog))

! Initialization
mpl%nprog = nprog
mpl%progint = ddis
mpl%done = .false.

end subroutine mpl_prog_init

!----------------------------------------------------------------------
! Subroutine: mpl_prog_print
! Purpose: print progression display
!----------------------------------------------------------------------
subroutine mpl_prog_print(mpl,i)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(in),optional :: i     ! Index

! Local variables
integer :: ithread
real(kind_real) :: prog

! Update progression array
if (present(i)) mpl%done(i) = .true.

! Print message
prog = 100.0*real(count(mpl%done),kind_real)/real(mpl%nprog,kind_real)
ithread = 0
!$ ithread = omp_get_thread_num()
do while ((int(prog)>mpl%progint).and.(ithread==0))
   if (mpl%progint<100) then
      if (mpl%progint<10) then
         write(mpl%info,'(i2,a)') mpl%progint,'% '
      else
         write(mpl%info,'(i3,a)') mpl%progint,'% '
      end if
      call mpl%flush(.false.)
   end if
   mpl%progint = mpl%progint+ddis
end do

end subroutine mpl_prog_print

!----------------------------------------------------------------------
! Subroutine: mpl_prog_final
! Purpose: finalize progression display
!----------------------------------------------------------------------
subroutine mpl_prog_final(mpl,advance_flag)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl        ! MPI data
logical,intent(in),optional :: advance_flag ! Advance flag

! Local variables
integer :: ithread
logical :: ladvance_flag

! Set advance flag
ladvance_flag = .true.
if (present(advance_flag)) ladvance_flag = advance_flag

! Print message
ithread = 0
!$ ithread = omp_get_thread_num()
do while ((mpl%progint<=100).and.(ithread==0))
   if (mpl%progint<100) then
      if (mpl%progint<10) then
         write(mpl%info,'(i2,a)') mpl%progint,'% '
      else
         write(mpl%info,'(i3,a)') mpl%progint,'% '
      end if
      call mpl%flush(.false.)
   else
      write(mpl%info,'(a)') ' 100%'
      call mpl%flush(ladvance_flag)
   end if
   mpl%progint = mpl%progint+ddis
end do

! Release memory
deallocate(mpl%done)

end subroutine mpl_prog_final

!----------------------------------------------------------------------
! Function: mpl_nc_file_create_or_open
! Purpose: create or open NetCDF file
!----------------------------------------------------------------------
function mpl_nc_file_create_or_open(mpl,subr,filename,f_comm) result (ncid)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl               ! MPI data
character(len=*),intent(in) :: subr                ! Calling subroutine
character(len=*),intent(in) :: filename            ! File name
type(fckit_mpi_comm),intent(in),optional :: f_comm ! Communicator

! Returned variable
integer :: ncid                                    ! NetCDF file ID

! Local variables
integer :: info

! Create file
if (present(f_comm)) then
   info = nf90_create(filename,ior(nf90_noclobber,ior(nf90_netcdf4,nf90_mpiio)),ncid, &
 & comm=f_comm%communicator(),info=f_comm%info_null())
else
   info = nf90_create(filename,ior(nf90_noclobber,nf90_netcdf4),ncid)
end if

if (info==nf90_noerr) then
   ! Set real missing value
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',mpl%msv%valr))
else
   ! Open file
   if (present(f_comm)) then
      call mpl%ncerr(subr,nf90_open(filename,nf90_write,ncid, &
 & comm=f_comm%communicator(),info=f_comm%info_null()))
   else
      call mpl%ncerr(subr,nf90_open(filename,nf90_write,ncid))
   end if
end if

end function mpl_nc_file_create_or_open

!----------------------------------------------------------------------
! Function: mpl_nc_group_define_or_get
! Purpose: define or get group
!----------------------------------------------------------------------
function mpl_nc_group_define_or_get(mpl,subr,ncid,grpname) result (grpid)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: subr    ! Calling subroutine
integer,intent(in) :: ncid             ! NetCDF file ID
character(len=*),intent(in) :: grpname ! Group name

! Returned variable
integer :: grpid                       ! NetCDF group ID

! Local variables
integer :: info

! Define group
info = nf90_def_grp(ncid,grpname,grpid)

! Get group
if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_inq_grp_ncid(ncid,grpname,grpid))

end function mpl_nc_group_define_or_get

!----------------------------------------------------------------------
! Function: mpl_nc_dim_define_or_get
! Purpose: define or get (and check) NetCDF dimension
!----------------------------------------------------------------------
function mpl_nc_dim_define_or_get(mpl,subr,ncid,dimname,dimsize) result (dimid)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: subr    ! Calling subroutine
integer,intent(in) :: ncid             ! NetCDF file ID
character(len=*),intent(in) :: dimname ! Dimension name
integer,intent(in) :: dimsize          ! Dimension size

! Returned variable
integer :: dimid                       ! NetCDF dimension ID

! Local variables
integer :: info,dimsize_test

! Define dimension
info = nf90_def_dim(ncid,dimname,dimsize,dimid)

if (info/=nf90_noerr) then
   ! Get dimension
   call mpl%ncerr(subr,nf90_inq_dimid(ncid,dimname,dimid))

   ! Get dimension size
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,dimid,len=dimsize_test))

   ! Check dimension size
   if (dimsize_test/=dimsize) call mpl%abort(subr,'dimension '//trim(dimname)//' has a different size in file')
end if

end function mpl_nc_dim_define_or_get

!----------------------------------------------------------------------
! Function: mpl_nc_dim_inquire
! Purpose: inquire NetCDF file dimension size
!----------------------------------------------------------------------
function mpl_nc_dim_inquire(mpl,subr,ncid,dimname) result (dimsize)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: subr    ! Calling subroutine
integer,intent(in) :: ncid             ! NetCDF file ID
character(len=*),intent(in) :: dimname ! Dimension name

! Returned variable
integer :: dimsize                     ! NetCDF dimension size

! Local variables
integer :: info,dimid

! Get dimension ID
info = nf90_inq_dimid(ncid,dimname,dimid)
if (info==nf90_noerr) then
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,dimid,len=dimsize))
else
   dimsize = 0
end if

end function mpl_nc_dim_inquire

!----------------------------------------------------------------------
! Subroutine: mpl_nc_dim_check
! Purpose: check if NetCDF file dimension exists and has the right size
!----------------------------------------------------------------------
subroutine mpl_nc_dim_check(mpl,subr,ncid,dimname,dimsize)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
character(len=*),intent(in) :: subr    ! Calling subroutine
integer,intent(in) :: ncid             ! NetCDF file ID
character(len=*),intent(in) :: dimname ! Dimension name
integer,intent(in) :: dimsize          ! Expected dimension size

! Local variables
integer :: dimid,dimsize_test

! Get dimension
call mpl%ncerr(subr,nf90_inq_dimid(ncid,dimname,dimid))

! Get dimension size
call mpl%ncerr(subr,nf90_inquire_dimension(ncid,dimid,len=dimsize_test))

! Check dimension size
if (dimsize_test/=dimsize) call mpl%abort(subr,'dimension '//trim(dimname)//' has a different size in file')

end subroutine mpl_nc_dim_check

!----------------------------------------------------------------------
! Function: mpl_nc_var_define_or_get
! Purpose: define or get NetCDF variable
!----------------------------------------------------------------------
function mpl_nc_var_define_or_get(mpl,subr,ncid,varname,varkind,varshape,unitname) result (varid)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl             ! MPI data
character(len=*),intent(in) :: subr              ! Calling subroutine
integer,intent(in) :: ncid                       ! NetCDF file ID
character(len=*),intent(in) :: varname           ! Variable name
integer,intent(in) :: varkind                    ! Variable kind
integer,intent(in) :: varshape(:)                ! Variable shape
character(len=*),intent(in),optional :: unitname ! Unit name

! Returned variable
integer :: varid                       ! NetCDF variable ID

! Local variables
integer :: info

! Define variable
info = nf90_def_var(ncid,varname,varkind,varshape,varid)

if (info==nf90_noerr) then
   ! Set missing value attribute
   if (varkind==nf90_int) then
      call mpl%ncerr(subr,nf90_put_att(ncid,varid,'_FillValue',mpl%msv%vali))
   elseif (varkind==nc_kind_real) then
      call mpl%ncerr(subr,nf90_put_att(ncid,varid,'_FillValue',mpl%msv%valr))
   else
      call mpl%abort(subr,'wrong variable kind')
   end if

   ! Set unit
   if (present(unitname)) call mpl%ncerr(subr,nf90_put_att(ncid,varid,'unit',unitname))
else
   ! Get variable
   call mpl%ncerr(subr,nf90_inq_varid(ncid,varname,varid))
end if

end function mpl_nc_var_define_or_get

!----------------------------------------------------------------------
! Subroutine: mpl_ncerr
! Purpose: handle NetCDF error
!----------------------------------------------------------------------
subroutine mpl_ncerr(mpl,subr,info)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl ! MPI data
character(len=*),intent(in) :: subr  ! Calling subroutine
integer,intent(in) :: info           ! Info index

! Check status
if (info/=nf90_noerr) call mpl%abort(subr,nf90_strerror(info))

end subroutine mpl_ncerr

!----------------------------------------------------------------------
! Subroutine: mpl_write_integer
! Purpose: write integer into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_integer(mpl,ncid,prefix,variables,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: prefix  ! Prefix
character(len=*),intent(in) :: variables ! Variable name
integer,intent(in) :: var              ! Integer

! Local variables
integer :: delta
character(len=1024) :: for
character(len=1024),parameter :: subr = 'mpl_write_integer'

if (mpl%msv%is(ncid)) then
   ! Write integer into a log file
   delta = 2
   if (var<0) delta = delta+1
   if (abs(var)>0) then
      write(for,'(a,i4.4,a,i4.4,a,i4.4,a)') '(a10,a',len_trim(variables),',a',25-len_trim(variables),',a,i', &
       & floor(log(abs(real(var,kind_real)))/log(10.0))+delta,',a)'
   else
      write(for,'(a,i4.4,a,i4.4,a,i4.4,a)') '(a10,a',len_trim(variables),',a',25-len_trim(variables),',a,i',delta,',a)'
   end if
   write(mpl%info,for) '',trim(variables),'',':',var
   call mpl%flush
else
   ! Write integer into a NetCDF file
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(prefix)//'_'//trim(variables),var))
end if

end subroutine mpl_write_integer

!----------------------------------------------------------------------
! Subroutine: mpl_write_integer_array
! Purpose: write integer array into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_integer_array(mpl,ncid,prefix,variables,n,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: prefix  ! Prefix
character(len=*),intent(in) :: variables ! Variable name
integer,intent(in) :: n                ! Integer array size
integer,intent(in) :: var(n)           ! Integer array

! Local variables
integer :: i,delta
character(len=1024) :: for,fullstr
character(len=1024),parameter :: subr = 'mpl_write_integer_array'

if (mpl%msv%is(ncid)) then
   ! Write integer array into a log file
   write(for,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(variables),',a',25-len_trim(variables),',a)'
   write(mpl%info,for) '',trim(variables),'',':'
   call mpl%flush(.false.)
   do i=1,n
      delta = 2
      if (var(i)<0) delta = delta+1
      if (abs(var(i))>0) then
         write(for,'(a,i4.4,a)') '(i',floor(log(abs(real(var(i),kind_real)))/log(10.0))+delta,',a)'
      else
         write(for,'(a,i4.4,a)') '(i',delta,',a)'
      end if
      write(mpl%info,for) var(i),','
      call mpl%flush(.false.)
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
else
   ! Write integer array into a NetCDF file
   if (n>0) then
      write(fullstr,'(i3.3)') var(1)
      do i=2,n
         write(fullstr,'(a,i3.3)') trim(fullstr(1:1024-4))//':',var(i)
      end do
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(prefix)//'_'//trim(variables),fullstr))
   end if
end if

end subroutine mpl_write_integer_array

!----------------------------------------------------------------------
! Subroutine: mpl_write_real
! Purpose: write real into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_real(mpl,ncid,prefix,variables,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: prefix  ! Prefix
character(len=*),intent(in) :: variables ! Variable name
real(kind_real),intent(in) :: var      ! Real

! Local variables
integer :: delta
character(len=1024) :: for
character(len=1024),parameter :: subr = 'mpl_write_real'

if (mpl%msv%is(ncid)) then
   ! Write real into a log file
   delta = 10
   if (var<0.0) delta = delta+1
   write(for,'(a,i4.4,a,i4.4,a,i2,a)') '(a10,a',len_trim(variables),',a',25-len_trim(variables),',a,e',delta,'.3,a)'
   write(mpl%info,for) '',trim(variables),'',':',var
   call mpl%flush
else
   ! Write real into a NetCDF file
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(prefix)//'_'//trim(variables),var))
end if

end subroutine mpl_write_real

!----------------------------------------------------------------------
! Subroutine: mpl_write_real_array
! Purpose: write real array into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_real_array(mpl,ncid,prefix,variables,n,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: prefix  ! Prefix
character(len=*),intent(in) :: variables ! Variable name
integer,intent(in) :: n                ! Real array size
real(kind_real),intent(in) :: var(n)   ! Real array

! Local variables
integer :: i,delta
character(len=1024) :: for,fullstr
character(len=1024),parameter :: subr = 'mpl_write_real_array'

if (mpl%msv%is(ncid)) then
   ! Write real array into a log file
   write(for,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(variables),',a',25-len_trim(variables),',a)'
   write(mpl%info,for) '',trim(variables),'',':'
   call mpl%flush(.false.)
   do i=1,n
      delta = 10
      if (var(i)<0.0) delta = delta+1
      write(for,'(a,i2,a)') '(e',delta,'.3,a)'
      write(mpl%info,for) var(i),','
      call mpl%flush(.false.)
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
else
   ! Write real array into a NetCDF file
   if (n>0) then
      write(fullstr,'(e10.3)') var(1)
      do i=2,n
         write(fullstr,'(a,e10.3)') trim(fullstr(1:1024-11))//':',var(i)
      end do
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(prefix)//'_'//trim(variables),fullstr))
   end if
end if

end subroutine mpl_write_real_array

!----------------------------------------------------------------------
! Subroutine: mpl_write_logical
! Purpose: write logical into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_logical(mpl,ncid,prefix,variables,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: prefix  ! Prefix
character(len=*),intent(in) :: variables ! Variable name
logical,intent(in) :: var              ! Logical

! Local variables
character(len=1024) :: for
character(len=1024),parameter :: subr = 'mpl_write_logical'

if (mpl%msv%is(ncid)) then
   ! Write logical into a log file
   write(for,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(variables),',a',25-len_trim(variables),',a,l2,a)'
   write(mpl%info,for) '',trim(variables),'',':',var
   call mpl%flush
else
   ! Write logical into a NetCDF file
   if (var) then
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(prefix)//'_'//trim(variables),'.true.'))
   else
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(prefix)//'_'//trim(variables),'.false.'))
   end if
end if

end subroutine mpl_write_logical

!----------------------------------------------------------------------
! Subroutine: mpl_write_logical_array
! Purpose: write logical array into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_logical_array(mpl,ncid,prefix,variables,n,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: prefix  ! Prefix
character(len=*),intent(in) :: variables ! Variable name
integer,intent(in) :: n                ! Real array size
logical,intent(in) :: var(n)           ! Logical array

! Local variables
integer :: i
character(len=1024) :: for,fullstr
character(len=1024),parameter :: subr = 'mpl_write_logical_array'

if (mpl%msv%is(ncid)) then
   ! Write logical array into a log file
   write(for,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(variables),',a',25-len_trim(variables),',a)'
   write(mpl%info,for) '',trim(variables),'',':'
   call mpl%flush(.false.)
   do i=1,n
      write(for,'(a)') '(l2,a)'
      write(mpl%info,for) var(i),','
      call mpl%flush(.false.)
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
else
   ! Write real array into a NetCDF file
   if (n>0) then
      if (var(1)) then
         fullstr = '.true.'
      else
         fullstr = '.false.'
      end if
      do i=2,n
         if (var(i)) then
            fullstr = trim(fullstr(1:1024-7))//':'//'.true.'
         else
            fullstr = trim(fullstr(1:1024-8))//':'//'.false.'
         end if
      end do
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(prefix)//'_'//trim(variables),fullstr))
   end if
end if

end subroutine mpl_write_logical_array

!----------------------------------------------------------------------
! Subroutine: mpl_write_string
! Purpose: write string into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_string(mpl,ncid,prefix,variables,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: prefix  ! Prefix
character(len=*),intent(in) :: variables ! Variable name
character(len=*),intent(in) :: var     ! String

! Local variables
character(len=1024) :: str
character(len=1024),parameter :: subr = 'mpl_write_string'

if (mpl%msv%is(ncid)) then
   ! Write string into a log file
   if (len_trim(var)>0) then
      write(str,'(a,i4.4,a,i4.4,a,i4.4,a)') '(a10,a',len_trim(variables),',a',25-len_trim(variables),',a,a1,a',len_trim(var),')'
      write(mpl%info,str) '',trim(variables),'',':','',trim(var)
   else
      write(str,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(variables),',a',25-len_trim(variables),',a)'
      write(mpl%info,str) '',trim(variables),'',':'
   end if
   call mpl%flush
else
   ! Write string into a NetCDF file
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(prefix)//'_'//trim(variables),var))
end if

end subroutine mpl_write_string

!----------------------------------------------------------------------
! Subroutine: mpl_write_string_array
! Purpose: write string array into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mpl_write_string_array(mpl,ncid,prefix,variables,n,var)

implicit none

! Passed variables
class(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid             ! NetCDF file id
character(len=*),intent(in) :: prefix  ! Prefix
character(len=*),intent(in) :: variables ! Variable name
integer,intent(in) :: n                ! String array size
character(len=*),intent(in) :: var(n)  ! String array

! Local variables
integer :: i
character(len=1024) :: for,fullstr
character(len=1024),parameter :: subr = 'mpl_write_string_array'

if (mpl%msv%is(ncid)) then
   ! Write string array into a log file
   write(for,'(a,i4.4,a,i4.4,a)') '(a10,a',len_trim(variables),',a',25-len_trim(variables),',a)'
   write(mpl%info,for) '',trim(variables),'',':'
   call mpl%flush(.false.)
   do i=1,n
      if (len_trim(var(i))>0) then
         write(for,'(a,i4.4,a)') '(a1,a',len_trim(var(i)),',a)'
         write(mpl%info,for) '',trim(var(i)),','
      else
         write(mpl%info,'(a1,a)') '',','
      end if
      call mpl%flush(.false.)
   end do
   write(mpl%info,'(a)') ''
   call mpl%flush
else
   ! Write string array into a NetCDF file
   if (n>0) then
      fullstr = trim(var(1))
      do i=2,n
         fullstr = trim(fullstr)//':'//trim(var(i))
      end do
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(prefix)//'_'//trim(variables),fullstr))
   end if
end if

end subroutine mpl_write_string_array

end module type_mpl
