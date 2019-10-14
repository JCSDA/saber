!----------------------------------------------------------------------
! subroutine: bump_main
! Purpose: call to the BUMP library
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
subroutine bump_main(n1,arg1,n2,arg2) bind (c,name='bump_main_f90')

use fckit_mpi_module, only: fckit_mpi_comm
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
integer(c_int),intent(in) :: n1
character(c_char),intent(in) :: arg1(n1)
integer(c_int),intent(in) :: n2
character(c_char),intent(in) :: arg2(n2)

! Local variables
integer :: i,narg,info,ppos,iproc,ie,ifileunit,iv,its
real(kind_real),allocatable :: fld_c0(:,:),fld_c0a(:,:),fld_mga(:,:,:,:)
character(len=1024) :: inputfile,logdir,ext,filename
type(bump_type) :: bump
type(fckit_mpi_comm) :: f_comm
type(model_type) :: model
type(mpl_type) :: mpl
type(rng_type) :: rng
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
call mpl%msv%init(-999,-999.0_kind_real)

! Initialize MPL
call mpl%init(f_comm)

! Initialize timer
call timer%start(mpl)

! Initialize namelist
call bump%nam%init

! Find whether input file is a namelist (xxx.nam) or a yaml (xxxx.yaml) and read it
ppos = scan(trim(inputfile),".",back=.true.)
ext = trim(inputfile(ppos+1:))
select case (trim(ext))
case ('nam')
   ! Namelist
   call bump%nam%read(mpl,inputfile)
case ('yaml')
   ! yaml
   call bump%nam%read_yaml(mpl,inputfile)
case default
   ! Wrong extension
   write(output_unit,'(a)') 'Error: input file has a wrong extension (should be .nam or .yaml)'
   call flush(output_unit)
   error stop 3
end select

! Broadcast namelist
call bump%nam%bcast(mpl)
if (mpl%main) call bump%nam%write(mpl)

! Define info unit and open file
do iproc=1,mpl%nproc
   if ((trim(bump%nam%verbosity)=='all').or.((trim(bump%nam%verbosity)=='main').and.(iproc==mpl%rootproc))) then
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
call model%setup(mpl,bump%nam)

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

! Proper dimensions
model%lon_mga = model%lon_mga*rad2deg
model%lat_mga = model%lat_mga*rad2deg
model%area_mga = model%area_mga*req**2

! BUMP setup
call bump%setup_online(f_comm,model%nmga,model%nl0,bump%nam%nv,bump%nam%nts, &
                     & model%lon_mga,model%lat_mga,model%area_mga,model%vunit_mga,model%mask_mga, &
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

! Test set_parameters interfaces
if (bump%nam%check_set_param_cor.or.bump%nam%check_set_param_hyb.or.bump%nam%check_set_param_lct) then
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(bump%mpl%info,'(a)') '--- Set parameters into BUMP'
   call bump%mpl%flush()

   ! Allocation
   if (mpl%main) allocate(fld_c0(bump%geom%nc0,bump%geom%nl0))
   allocate(fld_c0a(bump%geom%nc0a,bump%geom%nl0))
   allocate(fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts))

   ! Random initialization
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         if (mpl%main) call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_c0)
         call bump%mpl%glb_to_loc(bump%geom%nl0,bump%geom%nc0,bump%geom%c0_to_proc,bump%geom%c0_to_c0a, &
       & fld_c0,bump%geom%nc0a,fld_c0a)
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a,fld_mga(:,:,iv,its))
      end do
   end do

   ! Set parameter
   if (bump%nam%check_set_param_cor) then
      call bump%set_parameter('var',fld_mga)
      call bump%set_parameter('cor_rh',fld_mga*req)
      call bump%set_parameter('cor_rv',fld_mga)
      call bump%set_parameter('cor_rv_rfac',fld_mga)
      call bump%set_parameter('cor_rv_coef',fld_mga)
   elseif (bump%nam%check_set_param_hyb) then
      call bump%set_parameter('loc_coef',fld_mga)
      call bump%set_parameter('loc_rh',fld_mga*req)
      call bump%set_parameter('loc_rv',fld_mga)
      call bump%set_parameter('hyb_coef',fld_mga)
   elseif (bump%nam%check_set_param_lct) then
      call bump%set_parameter('D11',fld_mga*req**2)
      call bump%set_parameter('D22',fld_mga*req**2)
      call bump%set_parameter('D33',fld_mga)
      call bump%set_parameter('D12',fld_mga)
      call bump%set_parameter('Dcoef',fld_mga)
   end if

   ! Release memory
   deallocate(fld_mga)
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
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(bump%mpl%info,'(a)') '--- Get parameters from BUMP'
   call bump%mpl%flush()

   ! Allocation
   allocate(fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts))

   ! Get parameter
   if (bump%nam%check_get_param_cor) then
      call bump%get_parameter('var',fld_mga)
      call bump%get_parameter('cor_rh',fld_mga)
      call bump%get_parameter('cor_rv',fld_mga)
      call bump%get_parameter('cor_rv_rfac',fld_mga)
      call bump%get_parameter('cor_rv_coef',fld_mga)
   elseif (bump%nam%check_get_param_hyb) then
      call bump%get_parameter('loc_coef',fld_mga)
      call bump%get_parameter('loc_rh',fld_mga)
      call bump%get_parameter('loc_rv',fld_mga)
      call bump%get_parameter('hyb_coef',fld_mga)
   elseif (bump%nam%check_get_param_Dloc) then
      call bump%get_parameter('loc_D11',fld_mga)
      call bump%get_parameter('loc_D22',fld_mga)
      call bump%get_parameter('loc_D33',fld_mga)
      call bump%get_parameter('loc_D12',fld_mga)
      call bump%get_parameter('loc_Dcoef',fld_mga)
      call bump%get_parameter('loc_DLh',fld_mga)
   elseif (bump%nam%check_get_param_lct) then
      call bump%get_parameter('D11_1',fld_mga)
      call bump%get_parameter('D22_1',fld_mga)
      call bump%get_parameter('D33_1',fld_mga)
      call bump%get_parameter('D12_1',fld_mga)
      call bump%get_parameter('Dcoef_1',fld_mga)
      call bump%get_parameter('DLh_1',fld_mga)
      call bump%get_parameter('D11_2',fld_mga)
      call bump%get_parameter('D22_2',fld_mga)
      call bump%get_parameter('D33_2',fld_mga)
      call bump%get_parameter('D12_2',fld_mga)
      call bump%get_parameter('Dcoef_2',fld_mga)
      call bump%get_parameter('DLh_2',fld_mga)
   end if

   ! Release memory
   deallocate(fld_mga)
end if

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
