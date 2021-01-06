!----------------------------------------------------------------------
! subroutine: bump_main
!> Call to the BUMP library
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
subroutine bump_main(n1,arg1,n2,arg2) bind (c,name='bump_main_f90')

use fckit_configuration_module, only: fckit_configuration,fckit_yamlconfiguration
use fckit_mpi_module, only: fckit_mpi_comm
use fckit_pathname_module, only : fckit_pathname
use iso_c_binding
use iso_fortran_env, only : output_unit
use type_bump, only: bump_type
use type_model, only: model_type
use type_mpl, only: mpl_type
use type_timer, only: timer_type

implicit none

! Passed variables
integer(c_int),intent(in) :: n1          !< First argument size
character(c_char),intent(in) :: arg1(n1) !< First argument
integer(c_int),intent(in) :: n2          !< Second argument size
character(c_char),intent(in) :: arg2(n2) !< Second argument

! Local variables
integer :: i,ppos,iproc,ie,ifileunit
character(len=1024) :: inputfile,logdir,ext,filename
type(fckit_configuration) :: conf

type(bump_type) :: bump
type(fckit_mpi_comm) :: f_comm
type(model_type) :: model
type(mpl_type) :: mpl
type(timer_type) :: timer

! Initialize MPI
f_comm = fckit_mpi_comm()

! Copy inputfile and logdir
inputfile = ''
do i=1,n1
   inputfile(i:i) = arg1(i)
end do
logdir = ''
do i=1,n2
   logdir(i:i) = arg2(i)
end do

! Set missing values
mpl%msv%vali = -999
mpl%msv%valr = -999.0

! Initialize MPL
call mpl%init(f_comm)

! Initialize timer
call timer%start(mpl)

! Initialize namelist
call bump%nam%init(mpl%nproc)

! Find whether input file is a namelist (xxx.nam) or a yaml (xxxx.yaml) and read it
ppos = scan(inputfile,".",back=.true.)
ext = inputfile(ppos+1:)
select case (trim(ext))
case ('nam')
   ! Namelist
   call bump%nam%read(mpl,inputfile)
case ('yaml')
   ! YAML file
   if (mpl%main) then
      ! Set fckit configuration from input file
      conf = fckit_yamlconfiguration(fckit_pathname(inputfile))

      ! Convert fckit configuration to namelist
      call bump%nam%from_conf(f_comm,conf)
   end if
case default
   ! Wrong extension
   write(output_unit,'(a)') 'Error: input file has a wrong extension (should be .nam or .yaml)'
   call flush(output_unit)
   error stop 3
end select

! Broadcast namelist
call bump%nam%bcast(mpl)

! Define info unit and open file
do iproc=1,mpl%nproc
   if ((trim(bump%nam%verbosity)=='all').or.((trim(bump%nam%verbosity)=='main').and.(iproc==mpl%rootproc))) then
      if (iproc==mpl%myproc) then
         ! Find a free unit
         call mpl%newunit(mpl%lunit)

         ! Open listing file
         write(filename,'(a,i6.6,a)') trim(bump%nam%prefix)//'.',mpl%myproc-1,'.out'
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

! Model setup
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Setup model'
call mpl%flush
call model%setup(mpl,bump%nam)

! Load ensembles
if (bump%nam%ens1_ne>0) then
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Load ensemble 1'
   call mpl%flush
   call model%load_ens(mpl,bump%nam,'ens1')
end if
if (bump%nam%ens2_ne>0) then
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Load ensemble 2'
   call mpl%flush
   call model%load_ens(mpl,bump%nam,'ens2')
end if

! BUMP setup
call bump%setup(f_comm,model%afunctionspace,model%fieldset,lunit=mpl%lunit,msvali=mpl%msv%vali,msvalr=mpl%msv%valr)

! Add members
if (bump%nam%ens1_ne>0) then
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Add members of ensemble 1'
   call mpl%flush
   do ie=1,bump%nam%ens1_ne
      write(mpl%info,'(a7,a,i4,a,i4)') '','Member ',ie,' / ',bump%nam%ens1_ne
      call mpl%flush
      call bump%add_member(model%ens1(ie),ie,1)
   end do
end if
if (bump%nam%ens2_ne>0) then
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Add members of ensemble 2'
   call mpl%flush
   do ie=1,bump%nam%ens2_ne
      write(mpl%info,'(a7,a,i4,a,i4)') '','Member ',ie,' / ',bump%nam%ens2_ne
      call mpl%flush
      call bump%add_member(model%ens2(ie),ie,2)
   end do
end if

! Test set_parameters interfaces
if (bump%nam%check_set_param_cor.or.bump%nam%check_set_param_hyb.or.bump%nam%check_set_param_lct) then
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Test set_parameters interfaces'
   call bump%mpl%flush()
   call bump%test_set_parameter
   if (bump%nam%default_seed) call bump%rng%reseed(mpl)
end if

! Run drivers
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Run drivers'
call mpl%flush
call bump%run_drivers

! Test get_parameter interfaces
if (bump%nam%check_get_param_cor.or.bump%nam%check_get_param_hyb.or.bump%nam%check_get_param_Dloc &
 & .or.bump%nam%check_get_param_lct) then
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Test get_parameter interfaces'
   call bump%mpl%flush()
   call bump%test_get_parameter
   if (bump%nam%default_seed) call bump%rng%reseed(mpl)
end if

! Release memory (partial)
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Release memory (partial)'
call mpl%flush
call bump%partial_dealloc

! Test interfaces
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Test apply interfaces'
call mpl%flush
call bump%test_apply_interfaces

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
