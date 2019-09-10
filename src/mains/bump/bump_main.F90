!----------------------------------------------------------------------
! subroutine: bump_main
! Purpose: call to the BUMP library
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
subroutine bump_main(n1,arg1,n2,arg2) bind (c,name='bump_main_f90')

use iso_c_binding
use iso_fortran_env, only : output_unit
use tools_const, only: rad2deg,req
use tools_kinds,only: kind_real
use type_bump, only: bump_type
use type_model, only: model_type
use type_mpl, only: mpl_type
use type_rng, only: rng_type
use type_timer, only: timer_type

implicit none

! Passed variables
integer,intent(in) :: n1
character,intent(in) :: arg1(n1)
integer,intent(in) :: n2
character,intent(in) :: arg2(n2)

! Local variables
integer :: i,iproc,ie,ifileunit
character(len=1024) :: namelname,logdir,filename
type(bump_type) :: bump
type(model_type) :: model
type(mpl_type) :: mpl
type(rng_type) :: rng
type(timer_type) :: timer

! Copy namelname and logdir
namelname = ''
do i=1,n1
   namelname(i:i) = arg1(i)
end do
logdir = ''
do i=1,n2
   logdir(i:i) = arg2(i)
end do

! Set missing values
call mpl%msv%init(-999,-999.0_kind_real)

! Initialize MPL
call mpl%init

! Initialize timer
call timer%start(mpl)

! Initialize, read and broadcast namelist
call bump%nam%init
call bump%nam%read(mpl,namelname)
call bump%nam%bcast(mpl)

! Define info unit and open file
do iproc=1,mpl%nproc
   if ((trim(bump%nam%verbosity)=='all').or.((trim(bump%nam%verbosity)=='main').and.(iproc==mpl%ioproc))) then
      if (iproc==mpl%myproc) then
         ! Find a free unit
         call mpl%newunit(mpl%lunit)

         ! Open listing file
         write(filename,'(a,i4.4,a)') trim(bump%nam%prefix)//'.',mpl%myproc-1,'.out'
         inquire(file=filename,number=ifileunit)
         if (ifileunit<0) then
            open(unit=mpl%lunit,file=trim(logdir)//'/'//trim(filename),action='write',status='replace')
         else
            close(ifileunit)
            open(unit=mpl%lunit,file=trim(logdir)//'/'//trim(filename),action='write',status='replace')
         end if
      end if
      call mpl%f_comm%barrier
   end if
end do

! Header
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- You are running the BUMP main program -------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Author: Benjamin Menetrier ------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT -----'
call mpl%flush

! Initialize random number generator
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Initialize random number generator'
call mpl%flush
call rng%init(mpl,bump%nam)

! Model setup
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Setup model'
call mpl%flush
call model%setup(mpl,rng,bump%nam)
if (bump%nam%default_seed) call rng%reseed(mpl)

! Load ensembles
if (bump%nam%ens1_ne>0) then
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Load ensemble 1'
   call mpl%flush
   call model%load_ens(mpl,bump%nam,'ens1')
else
   model%ens1_ne = 0
   model%ens1_nsub = 0
end if
if (bump%nam%ens2_ne>0) then
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Load ensemble 2'
   call mpl%flush
   call model%load_ens(mpl,bump%nam,'ens2')
else
   model%ens2_ne = 0
   model%ens2_nsub = 0
end if

if (bump%nam%new_obsop) then
   ! Generate observations locations
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Generate observations locations'
   call mpl%flush
   call model%generate_obs(mpl,rng,bump%nam)
   if (bump%nam%default_seed) call rng%reseed(mpl)
else
   model%nobsa = 0
   allocate(model%lonobs(model%nobsa))
   allocate(model%latobs(model%nobsa))
end if

! BUMP setup
call bump%setup_online(model%nmga,model%nl0,bump%nam%nv,bump%nam%nts, &
                     & model%lon_mga*rad2deg,model%lat_mga*rad2deg,model%area_mga*req**2,model%vunit_mga,model%mask_mga, &
                     & smask=model%smask_mga, &
                     & ens1_ne=model%ens1_ne,ens1_nsub=model%ens1_nsub,ens2_ne=model%ens2_ne,ens2_nsub=model%ens2_nsub, &
                     & nobs=model%nobsa,lonobs=model%lonobs*rad2deg,latobs=model%latobs*rad2deg, &
                     & lunit=mpl%lunit,msvali=mpl%msv%vali,msvalr=mpl%msv%valr)

! Transfer members
if (bump%nam%ens1_ne>0) then
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Transfer members of ensemble 1'
   call mpl%flush
   do ie=1,model%ens1_ne
      write(mpl%info,'(a7,a,i4,a,i4)') '','Member ',ie,' of ',model%ens1_ne
      call mpl%flush
      call bump%add_member(model%ens1(ie)%fld,ie,1)
      deallocate(model%ens1(ie)%fld)
   end do
end if
if (bump%nam%ens2_ne>0) then
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Transfer members of ensemble 2'
   call mpl%flush
   do ie=1,model%ens2_ne
      write(mpl%info,'(a7,a,i4,a,i4)') '','Member ',ie,' of ',model%ens2_ne
      call mpl%flush
      call bump%add_member(model%ens2(ie)%fld,ie,2)
      deallocate(model%ens2(ie)%fld)
   end do
end if

! Run drivers
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Run drivers'
call mpl%flush
call bump%run_drivers

! Execution stats
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Execution stats'
call timer%display(mpl)
call mpl%flush

if ((trim(bump%nam%verbosity)=='all').or.((trim(bump%nam%verbosity)=='main').and.mpl%main)) then
   ! Close listings
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Close listings'
   call mpl%flush
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   close(unit=mpl%lunit)
end if

! Finalize MPL
call mpl%final

! Release memory
call bump%dealloc
call model%dealloc

end subroutine bump_main
