!----------------------------------------------------------------------
! Module: type_com
! Purpose: communications derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_com

use fckit_mpi_module, only: fckit_mpi_status
use netcdf
!$ use omp_lib
use tools_kinds, only: kind_real,nc_kind_real
use type_mpl, only: mpl_type
use fckit_mpi_module, only: fckit_mpi_status

implicit none

! Communication derived type
type com_type
   ! Setup data
   integer,allocatable :: ext_to_proc(:) ! Extended index to processor
   integer,allocatable :: ext_to_red(:)  ! Extended index to reduced index

   ! Communication data
   character(len=1024) :: prefix         ! Communication prefix
   integer :: nred                       ! Reduction size
   integer :: next                       ! Extension size
   integer :: nown                       ! Own data size
   integer,allocatable :: own_to_ext(:)  ! Own data to extension conversion
   integer,allocatable :: own_to_red(:)  ! Own data to reduction conversion
   integer :: nhalo                      ! Halo buffer size
   integer :: nexcl                      ! Exclusive interior buffer size
   integer,allocatable :: jhalocounts(:) ! Halo counts
   integer,allocatable :: jexclcounts(:) ! Exclusive interior counts
   integer,allocatable :: jhalodispls(:) ! Halo displacement
   integer,allocatable :: jexcldispls(:) ! Exclusive interior displacement
   integer,allocatable :: halo(:)        ! Halo buffer
   integer,allocatable :: excl(:)        ! Exclusive interior buffer
contains
   procedure :: dealloc => com_dealloc
   procedure :: read => com_read
   procedure :: write => com_write
   procedure :: receive => com_receive
   procedure :: send => com_send
   procedure :: setup => com_setup
   procedure :: com_ext_1d
   procedure :: com_ext_2d
   generic :: ext => com_ext_1d,com_ext_2d
   procedure :: com_red_1d
   procedure :: com_red_2d
   generic :: red => com_red_1d,com_red_2d
end type com_type

integer,parameter :: com_ntag = 2 ! Number of communication steps to send/receive

private
public :: com_type,com_ntag

contains

!----------------------------------------------------------------------
! Subroutine: com_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine com_dealloc(com)

implicit none

! Passed variables
class(com_type),intent(inout) :: com ! Communication data

! Release memory
if (allocated(com%ext_to_proc)) deallocate(com%ext_to_proc)
if (allocated(com%ext_to_red)) deallocate(com%ext_to_red)
if (allocated(com%own_to_ext)) deallocate(com%own_to_ext)
if (allocated(com%own_to_red)) deallocate(com%own_to_red)
if (allocated(com%jhalocounts)) deallocate(com%jhalocounts)
if (allocated(com%jexclcounts)) deallocate(com%jexclcounts)
if (allocated(com%jhalodispls)) deallocate(com%jhalodispls)
if (allocated(com%jexcldispls)) deallocate(com%jexcldispls)
if (allocated(com%halo)) deallocate(com%halo)
if (allocated(com%excl)) deallocate(com%excl)

end subroutine com_dealloc

!----------------------------------------------------------------------
! Subroutine: com_read
! Purpose: read communications from a NetCDF file
!----------------------------------------------------------------------
subroutine com_read(com,mpl,ncid)

implicit none

! Passed variables
class(com_type),intent(inout) :: com  ! Communication data
type(mpl_type),intent(inout) :: mpl   ! MPI data
integer,intent(in) :: ncid            ! NetCDF file id

! Local variables
integer :: info
integer :: nred_id,next_id,nown_id,own_to_ext_id,own_to_red_id,nhalo_id,nexcl_id
integer :: jhalocounts_id,jexclcounts_id,jhalodispls_id,jexcldispls_id,halo_id,excl_id
character(len=1024),parameter :: subr = 'com_read'

! Get dimensions
info = nf90_inq_dimid(ncid,trim(com%prefix)//'_nred',nred_id)
if (info==nf90_noerr) then
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nred_id,len=com%nred))
else
   com%nred = 0
end if
info = nf90_inq_dimid(ncid,trim(com%prefix)//'_next',next_id)
if (info==nf90_noerr) then
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,next_id,len=com%next))
else
   com%next = 0
end if
info = nf90_inq_dimid(ncid,trim(com%prefix)//'_nown',nown_id)
if (info==nf90_noerr) then
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nown_id,len=com%nown))
else
   com%nown = 0
end if
info = nf90_inq_dimid(ncid,trim(com%prefix)//'_nhalo',nhalo_id)
if (info==nf90_noerr) then
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nhalo_id,len=com%nhalo))
else
   com%nhalo = 0
end if
info = nf90_inq_dimid(ncid,trim(com%prefix)//'_nexcl',nexcl_id)
if (info==nf90_noerr) then
   call mpl%ncerr(subr,nf90_inquire_dimension(ncid,nexcl_id,len=com%nexcl))
else
   com%nexcl = 0
end if

! Allocation
if (com%nown>0) allocate(com%own_to_ext(com%nown))
if (com%nown>0) allocate(com%own_to_red(com%nown))
allocate(com%jhalocounts(mpl%nproc))
allocate(com%jexclcounts(mpl%nproc))
allocate(com%jhalodispls(mpl%nproc))
allocate(com%jexcldispls(mpl%nproc))
if (com%nhalo>0) allocate(com%halo(com%nhalo))
if (com%nexcl>0) allocate(com%excl(com%nexcl))

! Get variables id
if (com%nown>0) call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(com%prefix)//'_own_to_ext',own_to_ext_id))
if (com%nown>0) call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(com%prefix)//'_own_to_red',own_to_red_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(com%prefix)//'_jhalocounts',jhalocounts_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(com%prefix)//'_jexclcounts',jexclcounts_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(com%prefix)//'_jhalodispls',jhalodispls_id))
call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(com%prefix)//'_jexcldispls',jexcldispls_id))
if (com%nhalo>0) call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(com%prefix)//'_halo',halo_id))
if (com%nexcl>0) call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(com%prefix)//'_excl',excl_id))

! Get variable
if (com%nown>0) call mpl%ncerr(subr,nf90_get_var(ncid,own_to_ext_id,com%own_to_ext))
if (com%nown>0) call mpl%ncerr(subr,nf90_get_var(ncid,own_to_red_id,com%own_to_red))
call mpl%ncerr(subr,nf90_get_var(ncid,jhalocounts_id,com%jhalocounts))
call mpl%ncerr(subr,nf90_get_var(ncid,jexclcounts_id,com%jexclcounts))
call mpl%ncerr(subr,nf90_get_var(ncid,jhalodispls_id,com%jhalodispls))
call mpl%ncerr(subr,nf90_get_var(ncid,jexcldispls_id,com%jexcldispls))
if (com%nhalo>0) call mpl%ncerr(subr,nf90_get_var(ncid,halo_id,com%halo))
if (com%nexcl>0) call mpl%ncerr(subr,nf90_get_var(ncid,excl_id,com%excl))

end subroutine com_read

!----------------------------------------------------------------------
! Subroutine: com_write
! Purpose: write communications to a NetCDF file
!----------------------------------------------------------------------
subroutine com_write(com,mpl,ncid)

implicit none

! Passed variables
class(com_type),intent(in) :: com   ! Communication data
type(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(in) :: ncid          ! NetCDF file id

! Local variables
integer :: info
integer :: nproc_id,nred_id,next_id,nown_id,own_to_ext_id,own_to_red_id,nhalo_id,nexcl_id
integer :: jhalocounts_id,jexclcounts_id,jhalodispls_id,jexcldispls_id,halo_id,excl_id
character(len=1024),parameter :: subr = 'com_write'

! Start definition mode
call mpl%ncerr(subr,nf90_redef(ncid))

! Define dimensions
info = nf90_inq_dimid(ncid,'nproc',nproc_id)
if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_def_dim(ncid,'nproc',mpl%nproc,nproc_id))
if (com%nred>0) call mpl%ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nred',com%nred,nred_id))
if (com%next>0) call mpl%ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_next',com%next,next_id))
if (com%nown>0) call mpl%ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nown',com%nown,nown_id))
if (com%nhalo>0) call mpl%ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nhalo',com%nhalo,nhalo_id))
if (com%nexcl>0) call mpl%ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nexcl',com%nexcl,nexcl_id))

! Define variables
if (com%nown>0) call mpl%ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_own_to_ext',nf90_int,(/nown_id/),own_to_ext_id))
if (com%nown>0) call mpl%ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_own_to_red',nf90_int,(/nown_id/),own_to_red_id))
call mpl%ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jhalocounts',nf90_int,(/nproc_id/),jhalocounts_id))
call mpl%ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jexclcounts',nf90_int,(/nproc_id/),jexclcounts_id))
call mpl%ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jhalodispls',nf90_int,(/nproc_id/),jhalodispls_id))
call mpl%ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jexcldispls',nf90_int,(/nproc_id/),jexcldispls_id))
if (com%nhalo>0) call mpl%ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_halo',nf90_int,(/nhalo_id/),halo_id))
if (com%nexcl>0) call mpl%ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_excl',nf90_int,(/nexcl_id/),excl_id))

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Put variables
if (com%nown>0) call mpl%ncerr(subr,nf90_put_var(ncid,own_to_ext_id,com%own_to_ext))
if (com%nown>0) call mpl%ncerr(subr,nf90_put_var(ncid,own_to_red_id,com%own_to_red))
call mpl%ncerr(subr,nf90_put_var(ncid,jhalocounts_id,com%jhalocounts))
call mpl%ncerr(subr,nf90_put_var(ncid,jexclcounts_id,com%jexclcounts))
call mpl%ncerr(subr,nf90_put_var(ncid,jhalodispls_id,com%jhalodispls))
call mpl%ncerr(subr,nf90_put_var(ncid,jexcldispls_id,com%jexcldispls))
if (com%nhalo>0) call mpl%ncerr(subr,nf90_put_var(ncid,halo_id,com%halo))
if (com%nexcl>0) call mpl%ncerr(subr,nf90_put_var(ncid,excl_id,com%excl))

end subroutine com_write

!----------------------------------------------------------------------
! Subroutine: com_receive
! Purpose: receive
!----------------------------------------------------------------------
subroutine com_receive(com,mpl,iproc,tag_start)

implicit none

! Passed variables
class(com_type),intent(inout) :: com ! Communication data
type(mpl_type),intent(inout) :: mpl  ! MPI data
integer,intent(in) :: iproc          ! Source task index
integer,intent(in) :: tag_start      ! MPI tag

! Local variables
integer :: tag,n_dim,n_int,offset_int
integer,allocatable :: rbuf_dim(:),rbuf_int(:)
type(fckit_mpi_status) :: status

! Initialization
tag = tag_start
n_dim = 0
n_int = 0
offset_int = 0

! Define buffer size
n_dim = n_dim+5

! Allocation
allocate(rbuf_dim(n_dim))

! Receive buffer
call mpl%f_comm%receive(rbuf_dim,iproc-1,tag,status)
tag = tag+1

! Copy data
com%nred = rbuf_dim(1)
com%next = rbuf_dim(2)
com%nown = rbuf_dim(3)
com%nhalo = rbuf_dim(4)
com%nexcl = rbuf_dim(5)

! Release memory
deallocate(rbuf_dim)

! Allocation
if (com%nown>0) allocate(com%own_to_ext(com%nown))
if (com%nown>0) allocate(com%own_to_red(com%nown))
allocate(com%jhalocounts(mpl%nproc))
allocate(com%jexclcounts(mpl%nproc))
allocate(com%jhalodispls(mpl%nproc))
allocate(com%jexcldispls(mpl%nproc))
if (com%nhalo>0) allocate(com%halo(com%nhalo))
if (com%nexcl>0) allocate(com%excl(com%nexcl))

! Define buffer sizes
n_int = n_int+2*com%nown+4*mpl%nproc+com%nhalo+com%nexcl

! Allocation
allocate(rbuf_int(n_int))

! Receive buffers
call mpl%f_comm%receive(rbuf_int,iproc-1,tag,status)
tag = tag+1

! Copy data
if (com%nown>0) com%own_to_ext = rbuf_int(offset_int+1:offset_int+com%nown)
offset_int = offset_int+com%nown
if (com%nown>0) com%own_to_red = rbuf_int(offset_int+1:offset_int+com%nown)
offset_int = offset_int+com%nown
com%jhalocounts = rbuf_int(offset_int+1:offset_int+mpl%nproc)
offset_int = offset_int+mpl%nproc
com%jexclcounts = rbuf_int(offset_int+1:offset_int+mpl%nproc)
offset_int = offset_int+mpl%nproc
com%jhalodispls = rbuf_int(offset_int+1:offset_int+mpl%nproc)
offset_int = offset_int+mpl%nproc
com%jexcldispls = rbuf_int(offset_int+1:offset_int+mpl%nproc)
offset_int = offset_int+mpl%nproc
if (com%nhalo>0) com%halo = rbuf_int(offset_int+1:offset_int+com%nhalo)
offset_int = offset_int+com%nhalo
if (com%nexcl>0) com%excl = rbuf_int(offset_int+1:offset_int+com%nexcl)
offset_int = offset_int+com%nexcl

! Release memory
deallocate(rbuf_int)

end subroutine com_receive

!----------------------------------------------------------------------
! Subroutine: com_send
! Purpose: send
!----------------------------------------------------------------------
subroutine com_send(com,mpl,iproc,tag_start)

implicit none

! Passed variables
class(com_type),intent(in) :: com   ! Communication data
type(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(in) :: iproc         ! Destination task index
integer,intent(in) :: tag_start     ! MPI tag

! Local variables
integer :: tag,n_dim,n_int,offset_int
integer,allocatable :: sbuf_dim(:),sbuf_int(:)

! Initialization
tag = tag_start
n_dim = 0
n_int = 0
offset_int = 0

! Define buffer size
n_dim = n_dim+5

! Allocation
allocate(sbuf_dim(n_dim))

! Copy data
sbuf_dim(1) = com%nred
sbuf_dim(2) = com%next
sbuf_dim(3) = com%nown
sbuf_dim(4) = com%nhalo
sbuf_dim(5) = com%nexcl

! Send buffer
call mpl%f_comm%send(sbuf_dim,iproc-1,tag)
tag = tag+1

! Release memory
deallocate(sbuf_dim)

! Define buffer sizes
n_int = n_int+2*com%nown+4*mpl%nproc+com%nhalo+com%nexcl

! Allocation
allocate(sbuf_int(n_int))

! Copy data
if (com%nown>0) sbuf_int(offset_int+1:offset_int+com%nown) = com%own_to_ext
offset_int = offset_int+com%nown
if (com%nown>0) sbuf_int(offset_int+1:offset_int+com%nown) = com%own_to_red
offset_int = offset_int+com%nown
sbuf_int(offset_int+1:offset_int+mpl%nproc) = com%jhalocounts
offset_int = offset_int+mpl%nproc
sbuf_int(offset_int+1:offset_int+mpl%nproc) = com%jexclcounts
offset_int = offset_int+mpl%nproc
sbuf_int(offset_int+1:offset_int+mpl%nproc) = com%jhalodispls
offset_int = offset_int+mpl%nproc
sbuf_int(offset_int+1:offset_int+mpl%nproc) = com%jexcldispls
offset_int = offset_int+mpl%nproc
if (com%nhalo>0) sbuf_int(offset_int+1:offset_int+com%nhalo) = com%halo
offset_int = offset_int+com%nhalo
if (com%nexcl>0) sbuf_int(offset_int+1:offset_int+com%nexcl) = com%excl
offset_int = offset_int+com%nexcl

! Send buffers
call mpl%f_comm%send(sbuf_int,iproc-1,tag)
tag = tag+1

! Release memory
deallocate(sbuf_int)

end subroutine com_send

!----------------------------------------------------------------------
! Subroutine: com_setup
! Purpose: setup communications
!----------------------------------------------------------------------
subroutine com_setup(com_out,mpl,prefix,nglb,nred,next,nown,ext_to_glb,own_to_ext,glb_to_proc,glb_to_red)

implicit none

! Passed variables
class(com_type),intent(inout) :: com_out ! Communication data
type(mpl_type),intent(inout) :: mpl      ! MPI data
character(len=*),intent(in) :: prefix    ! Prefix
integer,intent(in) :: nglb               ! Global size
integer,intent(in) :: nred               ! Reduced halo size
integer,intent(in) :: next               ! Extended halo size
integer,intent(in) :: nown               ! Own data size
integer,intent(in) :: ext_to_glb(next)   ! Extended halo to global
integer,intent(in) :: own_to_ext(nown)   ! Own data to extended halo
integer,intent(in) :: glb_to_proc(nglb)  ! Global to processor
integer,intent(in) :: glb_to_red(nglb)   ! Global to reduced halo

! Local variables
integer :: iproc,jproc,iext,iglb,ired,icount,nred_tmp,next_tmp,nown_tmp
integer,allocatable :: ext_to_glb_tmp(:),own_to_ext_tmp(:),own_to_red(:),own_to_red_tmp(:)
type(com_type) :: com_in(mpl%nproc)
type(fckit_mpi_status) :: status

! Allocation
allocate(own_to_red(nown))

! Initialization
own_to_red = glb_to_red(ext_to_glb(own_to_ext))

if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%rootproc) then
         ! Copy dimensions
         nred_tmp = nred
         next_tmp = next
         nown_tmp = nown
      else
         ! Receive dimensions on rootproc
         call mpl%f_comm%receive(nred_tmp,iproc-1,mpl%tag,status)
         call mpl%f_comm%receive(next_tmp,iproc-1,mpl%tag+1,status)
         call mpl%f_comm%receive(nown_tmp,iproc-1,mpl%tag+2,status)
      end if

      ! Allocation
      allocate(ext_to_glb_tmp(next_tmp))
      allocate(own_to_ext_tmp(nown_tmp))
      allocate(own_to_red_tmp(nown_tmp))

      ! Communicate data
      if (iproc==mpl%rootproc) then
         ! Copy data
         ext_to_glb_tmp = ext_to_glb
         own_to_ext_tmp = own_to_ext
         own_to_red_tmp = own_to_red
      else
         ! Receive data on rootproc
         call mpl%f_comm%receive(ext_to_glb_tmp,iproc-1,mpl%tag+3,status)
         call mpl%f_comm%receive(own_to_ext_tmp,iproc-1,mpl%tag+4,status)
         call mpl%f_comm%receive(own_to_red_tmp,iproc-1,mpl%tag+5,status)
      end if

      ! Allocation
      com_in(iproc)%nred = nred_tmp
      com_in(iproc)%next = next_tmp
      com_in(iproc)%nown = nown_tmp
      allocate(com_in(iproc)%ext_to_proc(com_in(iproc)%next))
      allocate(com_in(iproc)%ext_to_red(com_in(iproc)%next))
      allocate(com_in(iproc)%own_to_ext(com_in(iproc)%nown))
      allocate(com_in(iproc)%own_to_red(com_in(iproc)%nown))

      ! Communication parameters
      do iext=1,next_tmp
         iglb = ext_to_glb_tmp(iext)
         com_in(iproc)%ext_to_proc(iext) = glb_to_proc(iglb)
         ired = glb_to_red(iglb)
         com_in(iproc)%ext_to_red(iext) = ired
      end do
      com_in(iproc)%own_to_ext = own_to_ext_tmp
      com_in(iproc)%own_to_red = own_to_red_tmp

      ! Release memory
      deallocate(ext_to_glb_tmp)
      deallocate(own_to_ext_tmp)
      deallocate(own_to_red_tmp)
   end do
else
   ! Send dimensions to rootproc
   call mpl%f_comm%send(nred,mpl%rootproc-1,mpl%tag)
   call mpl%f_comm%send(next,mpl%rootproc-1,mpl%tag+1)
   call mpl%f_comm%send(nown,mpl%rootproc-1,mpl%tag+2)

   ! Send data to rootproc
   call mpl%f_comm%send(ext_to_glb,mpl%rootproc-1,mpl%tag+3)
   call mpl%f_comm%send(own_to_ext,mpl%rootproc-1,mpl%tag+4)
   call mpl%f_comm%send(own_to_red,mpl%rootproc-1,mpl%tag+5)
end if
call mpl%update_tag(6)

! Release memory
deallocate(own_to_red)

if (mpl%main) then
   ! Allocation
   do iproc=1,mpl%nproc
      allocate(com_in(iproc)%jhalocounts(mpl%nproc))
      allocate(com_in(iproc)%jexclcounts(mpl%nproc))
      allocate(com_in(iproc)%jhalodispls(mpl%nproc))
      allocate(com_in(iproc)%jexcldispls(mpl%nproc))
   end do

   ! Initialization
   do iproc=1,mpl%nproc
      com_in(iproc)%jhalocounts = 0
      com_in(iproc)%jexclcounts = 0
      com_in(iproc)%jhalodispls = 0
      com_in(iproc)%jexcldispls = 0
   end do

   ! Compute counts
   do iproc=1,mpl%nproc
      do iext=1,com_in(iproc)%next
         jproc = com_in(iproc)%ext_to_proc(iext)
         if (jproc/=iproc) then
            ! Count of points sent from IPROC to JPROC
            com_in(iproc)%jhalocounts(jproc) = com_in(iproc)%jhalocounts(jproc)+1

            ! Count of points received on JPROC from IPROC
            com_in(jproc)%jexclcounts(iproc) = com_in(jproc)%jexclcounts(iproc)+1
         end if
      end do
   end do

   ! Compute displacement
   do iproc=1,mpl%nproc
      com_in(iproc)%jhalodispls(1) = 0
      com_in(iproc)%jexcldispls(1) = 0
      do jproc=2,mpl%nproc
         com_in(iproc)%jhalodispls(jproc) = com_in(iproc)%jhalodispls(jproc-1)+com_in(iproc)%jhalocounts(jproc-1)
         com_in(iproc)%jexcldispls(jproc) = com_in(iproc)%jexcldispls(jproc-1)+com_in(iproc)%jexclcounts(jproc-1)
      end do
   end do

   ! Allocation
   do iproc=1,mpl%nproc
      com_in(iproc)%nhalo = sum(com_in(iproc)%jhalocounts)
      com_in(iproc)%nexcl = sum(com_in(iproc)%jexclcounts)
      allocate(com_in(iproc)%halo(com_in(iproc)%nhalo))
      allocate(com_in(iproc)%excl(com_in(iproc)%nexcl))
   end do

   ! Fill halo array
   do iproc=1,mpl%nproc
      com_in(iproc)%jhalocounts = 0
   end do
   do iproc=1,mpl%nproc
      do iext=1,com_in(iproc)%next
         ! Check for halo points
         jproc = com_in(iproc)%ext_to_proc(iext)
         if (jproc/=iproc) then
            ! Count of points sent from IPROC to JPROC
            com_in(iproc)%jhalocounts(jproc) = com_in(iproc)%jhalocounts(jproc)+1

            ! Local index of points sent from IPROC to JPROC
            com_in(iproc)%halo(com_in(iproc)%jhalodispls(jproc)+com_in(iproc)%jhalocounts(jproc)) = iext
         end if
      end do
   end do

   ! Fill excl array
   do jproc=1,mpl%nproc
      ! Loop over processors sending data to JPROC
      do iproc=1,mpl%nproc
         do icount=1,com_in(iproc)%jhalocounts(jproc)
            ! Local index of points received on JPROC from IPROC
            com_in(jproc)%excl(com_in(jproc)%jexcldispls(iproc)+icount) = &
          & com_in(iproc)%ext_to_red(com_in(iproc)%halo(com_in(iproc)%jhalodispls(jproc)+icount))
         end do
      end do
   end do

   ! Communicate dimensions
   do iproc=1,mpl%nproc
      if (iproc==mpl%rootproc) then
         ! Copy dimensions
         com_out%nred = com_in(iproc)%nred
         com_out%next = com_in(iproc)%next
         com_out%nown = com_in(iproc)%nown
         com_out%nhalo = com_in(iproc)%nhalo
         com_out%nexcl = com_in(iproc)%nexcl
      else
         ! Send dimensions to iproc
         call mpl%f_comm%send(com_in(iproc)%nred,iproc-1,mpl%tag)
         call mpl%f_comm%send(com_in(iproc)%next,iproc-1,mpl%tag+1)
         call mpl%f_comm%send(com_in(iproc)%nown,iproc-1,mpl%tag+2)
         call mpl%f_comm%send(com_in(iproc)%nhalo,iproc-1,mpl%tag+3)
         call mpl%f_comm%send(com_in(iproc)%nexcl,iproc-1,mpl%tag+4)
      end if
   end do
else
   ! Receive dimensions from rootproc
   call mpl%f_comm%receive(com_out%nred,mpl%rootproc-1,mpl%tag,status)
   call mpl%f_comm%receive(com_out%next,mpl%rootproc-1,mpl%tag+1,status)
   call mpl%f_comm%receive(com_out%nown,mpl%rootproc-1,mpl%tag+2,status)
   call mpl%f_comm%receive(com_out%nhalo,mpl%rootproc-1,mpl%tag+3,status)
   call mpl%f_comm%receive(com_out%nexcl,mpl%rootproc-1,mpl%tag+4,status)
end if
call mpl%update_tag(4)

! Allocation
allocate(com_out%own_to_ext(com_out%nown))
allocate(com_out%own_to_red(com_out%nown))
allocate(com_out%jhalocounts(mpl%nproc))
allocate(com_out%jexclcounts(mpl%nproc))
allocate(com_out%jhalodispls(mpl%nproc))
allocate(com_out%jexcldispls(mpl%nproc))
allocate(com_out%halo(com_out%nhalo))
allocate(com_out%excl(com_out%nexcl))

! Communicate data
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%rootproc) then
         ! Copy data
         if (com_out%nown>0)com_out%own_to_ext = com_in(iproc)%own_to_ext
         if (com_out%nown>0)com_out%own_to_red = com_in(iproc)%own_to_red
         com_out%jhalocounts = com_in(iproc)%jhalocounts
         com_out%jexclcounts = com_in(iproc)%jexclcounts
         com_out%jhalodispls = com_in(iproc)%jhalodispls
         com_out%jexcldispls = com_in(iproc)%jexcldispls
         com_out%halo = com_in(iproc)%halo
         com_out%excl = com_in(iproc)%excl
      else
         ! Send dimensions to iproc
         if (com_in(iproc)%nown>0) call mpl%f_comm%send(com_in(iproc)%own_to_ext,iproc-1,mpl%tag)
         if (com_in(iproc)%nown>0) call mpl%f_comm%send(com_in(iproc)%own_to_red,iproc-1,mpl%tag+1)
         call mpl%f_comm%send(com_in(iproc)%jhalocounts,iproc-1,mpl%tag+2)
         call mpl%f_comm%send(com_in(iproc)%jexclcounts,iproc-1,mpl%tag+3)
         call mpl%f_comm%send(com_in(iproc)%jhalodispls,iproc-1,mpl%tag+4)
         call mpl%f_comm%send(com_in(iproc)%jexcldispls,iproc-1,mpl%tag+5)
         if (com_in(iproc)%nhalo>0) call mpl%f_comm%send(com_in(iproc)%halo,iproc-1,mpl%tag+6)
         if (com_in(iproc)%nexcl>0) call mpl%f_comm%send(com_in(iproc)%excl,iproc-1,mpl%tag+7)
      end if
   end do
else
   ! Receive dimensions from rootproc
   if (com_out%nown>0) call mpl%f_comm%receive(com_out%own_to_ext,mpl%rootproc-1,mpl%tag,status)
   if (com_out%nown>0) call mpl%f_comm%receive(com_out%own_to_red,mpl%rootproc-1,mpl%tag+1,status)
   call mpl%f_comm%receive(com_out%jhalocounts,mpl%rootproc-1,mpl%tag+2,status)
   call mpl%f_comm%receive(com_out%jexclcounts,mpl%rootproc-1,mpl%tag+3,status)
   call mpl%f_comm%receive(com_out%jhalodispls,mpl%rootproc-1,mpl%tag+4,status)
   call mpl%f_comm%receive(com_out%jexcldispls,mpl%rootproc-1,mpl%tag+5,status)
   if (com_out%nhalo>0) call mpl%f_comm%receive(com_out%halo,mpl%rootproc-1,mpl%tag+6,status)
   if (com_out%nexcl>0) call mpl%f_comm%receive(com_out%excl,mpl%rootproc-1,mpl%tag+7,status)
end if
call mpl%update_tag(8)

! Set prefix
com_out%prefix = trim(prefix)

end subroutine com_setup

!----------------------------------------------------------------------
! Subroutine: com_ext_1d
! Purpose: communicate field to halo (extension), 1d
!----------------------------------------------------------------------
subroutine com_ext_1d(com,mpl,vec_red,vec_ext)

implicit none

! Passed variables
class(com_type),intent(in) :: com                ! Communication data
type(mpl_type),intent(inout) :: mpl              ! MPI data
real(kind_real),intent(in) :: vec_red(com%nred)  ! Reduced vector
real(kind_real),intent(out) :: vec_ext(com%next) ! Extended vector

! Local variables
integer :: iexcl,iown,ihalo
real(kind_real),allocatable :: sbuf(:),rbuf(:)

! Allocation
allocate(sbuf(com%nexcl))
allocate(rbuf(com%nhalo))

! Prepare buffers to send
!$omp parallel do schedule(static) private(iexcl)
do iexcl=1,com%nexcl
   sbuf(iexcl) = vec_red(com%excl(iexcl))
end do
!$omp end parallel do

! Communication
call mpl%f_comm%alltoall(sbuf,com%jexclcounts,com%jexcldispls,rbuf,com%jhalocounts,com%jhalodispls)

! Copy interior
vec_ext = 0.0
!$omp parallel do schedule(static) private(iown)
do iown=1,com%nown
   vec_ext(com%own_to_ext(iown)) = vec_red(com%own_to_red(iown))
end do
!$omp end parallel do

! Copy halo
!$omp parallel do schedule(static) private(ihalo)
do ihalo=1,com%nhalo
   vec_ext(com%halo(ihalo)) = rbuf(ihalo)
end do
!$omp end parallel do

! Release memory
deallocate(sbuf)
deallocate(rbuf)

end subroutine com_ext_1d

!----------------------------------------------------------------------
! Subroutine: com_ext_2d
! Purpose: communicate field to halo (extension), 2d
!----------------------------------------------------------------------
subroutine com_ext_2d(com,mpl,nl,vec_red,vec_ext)

implicit none

! Passed variables
class(com_type),intent(in) :: com                   ! Communication data
type(mpl_type),intent(inout) :: mpl                 ! MPI data
integer,intent(in) :: nl                            ! Number of levels
real(kind_real),intent(in) :: vec_red(com%nred,nl)  ! Reduced vector
real(kind_real),intent(out) :: vec_ext(com%next,nl) ! Extended vector

! Local variables
integer :: il,iexcl,iown,ihalo
integer :: jexclcounts(mpl%nproc),jexcldispls(mpl%nproc),jhalocounts(mpl%nproc),jhalodispls(mpl%nproc)
real(kind_real),allocatable :: sbuf(:),rbuf(:)

! Allocation
allocate(sbuf(com%nexcl*nl))
allocate(rbuf(com%nhalo*nl))

! Prepare buffers to send
!$omp parallel do schedule(static) private(il,iexcl)
do il=1,nl
   do iexcl=1,com%nexcl
      sbuf((iexcl-1)*nl+il) = vec_red(com%excl(iexcl),il)
   end do
end do
!$omp end parallel do

! Communication
jexclcounts = com%jexclcounts*nl
jexcldispls = com%jexcldispls*nl
jhalocounts = com%jhalocounts*nl
jhalodispls = com%jhalodispls*nl
call mpl%f_comm%alltoall(sbuf,jexclcounts,jexcldispls,rbuf,jhalocounts,jhalodispls)

! Copy interior
vec_ext = 0.0
!$omp parallel do schedule(static) private(il,iown)
do il=1,nl
   do iown=1,com%nown
      vec_ext(com%own_to_ext(iown),il) = vec_red(com%own_to_red(iown),il)
   end do
end do
!$omp end parallel do

! Copy halo
!$omp parallel do schedule(static) private(il,ihalo)
do il=1,nl
   do ihalo=1,com%nhalo
      vec_ext(com%halo(ihalo),il) = rbuf((ihalo-1)*nl+il)
   end do
end do
!$omp end parallel do

! Release memory
deallocate(sbuf)
deallocate(rbuf)

end subroutine com_ext_2d

!----------------------------------------------------------------------
! Subroutine: com_red_1d
! Purpose: communicate vector from halo (reduction)
!----------------------------------------------------------------------
subroutine com_red_1d(com,mpl,vec_ext,vec_red)

implicit none

! Passed variables
class(com_type),intent(in) :: com                ! Communication data
type(mpl_type),intent(inout) :: mpl              ! MPI data
real(kind_real),intent(in) :: vec_ext(com%next)  ! Extended vector
real(kind_real),intent(out) :: vec_red(com%nred) ! Reduced vector

! Local variables
integer :: ihalo,iown,ired,iexcl,ithread
real(kind_real),allocatable :: sbuf(:),rbuf(:),vec_red_arr(:,:)

! Allocation
allocate(sbuf(com%nhalo))
allocate(rbuf(com%nexcl))
allocate(vec_red_arr(com%nred,mpl%nthread))

! Prepare buffers to send
!$omp parallel do schedule(static) private(ihalo)
do ihalo=1,com%nhalo
   sbuf(ihalo) = vec_ext(com%halo(ihalo))
end do
!$omp end parallel do

! Communication
call mpl%f_comm%alltoall(sbuf,com%jhalocounts,com%jhalodispls,rbuf,com%jexclcounts,com%jexcldispls)

! Copy interior
vec_red = 0.0
!$omp parallel do schedule(static) private(iown)
do iown=1,com%nown
   vec_red(com%own_to_red(iown)) = vec_ext(com%own_to_ext(iown))
end do
!$omp end parallel do

! Copy halo
vec_red_arr = 0.0
!$omp parallel do schedule(static) private(iexcl,ithread)
do iexcl=1,com%nexcl
   ithread = 1
!$ ithread = omp_get_thread_num()+1
   vec_red_arr(com%excl(iexcl),ithread) = vec_red_arr(com%excl(iexcl),ithread)+rbuf(iexcl)
end do
!$omp end parallel do

! Sum over threads
do ithread=1,mpl%nthread
   do ired=1,com%nred
      vec_red(ired) = vec_red(ired)+vec_red_arr(ired,ithread)
   end do
end do

! Release memory
deallocate(sbuf)
deallocate(rbuf)
deallocate(vec_red_arr)

end subroutine com_red_1d

!----------------------------------------------------------------------
! Subroutine: com_red_2d
! Purpose: communicate vector from halo (reduction)
!----------------------------------------------------------------------
subroutine com_red_2d(com,mpl,nl,vec_ext,vec_red)

implicit none

! Passed variables
class(com_type),intent(in) :: com                   ! Communication data
type(mpl_type),intent(inout) :: mpl                 ! MPI data
integer,intent(in) :: nl                            ! Number of levels
real(kind_real),intent(in) :: vec_ext(com%next,nl)  ! Extended vector
real(kind_real),intent(out) :: vec_red(com%nred,nl) ! Reduced vector

! Local variables
integer :: il,ihalo,iown,ired,iexcl,ithread
real(kind_real),allocatable :: sbuf(:),rbuf(:),vec_red_arr(:,:,:)

! Allocation
allocate(sbuf(com%nhalo*nl))
allocate(rbuf(com%nexcl*nl))
allocate(vec_red_arr(com%nred,nl,mpl%nthread))

! Prepare buffers to send
!$omp parallel do schedule(static) private(il,ihalo)
do il=1,nl
   do ihalo=1,com%nhalo
      sbuf((ihalo-1)*nl+il) = vec_ext(com%halo(ihalo),il)
   end do
end do
!$omp end parallel do

! Communication
call mpl%f_comm%alltoall(sbuf,com%jhalocounts*nl,com%jhalodispls*nl,rbuf,com%jexclcounts*nl,com%jexcldispls*nl)

! Copy interior
vec_red = 0.0
!$omp parallel do schedule(static) private(il,iown)
do il=1,nl
   do iown=1,com%nown
      vec_red(com%own_to_red(iown),il) = vec_ext(com%own_to_ext(iown),il)
   end do
end do
!$omp end parallel do

! Copy halo
vec_red_arr = 0.0
!$omp parallel do schedule(static) private(il,iexcl,ithread)
do il=1,nl
   do iexcl=1,com%nexcl
      ithread = 1
!$    ithread = omp_get_thread_num()+1
      vec_red_arr(com%excl(iexcl),il,ithread) = vec_red_arr(com%excl(iexcl),il,ithread)+rbuf((iexcl-1)*nl+il)
   end do
end do
!$omp end parallel do

! Sum over threads
do ithread=1,mpl%nthread
   do il=1,nl
      do ired=1,com%nred
         vec_red(ired,il) = vec_red(ired,il)+vec_red_arr(ired,il,ithread)
      end do
   end do
end do

! Release memory
deallocate(sbuf)
deallocate(rbuf)
deallocate(vec_red_arr)

end subroutine com_red_2d

end module type_com
