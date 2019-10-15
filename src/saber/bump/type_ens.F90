!----------------------------------------------------------------------
! Module: type_ens
! Purpose: ensemble derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_ens

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_max
use netcdf
use tools_const, only: deg2rad,rad2deg
use tools_func, only: sphere_dist
use tools_kinds, only: kind_real,nc_kind_real
use type_geom, only: geom_type
use type_io, only: io_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

! Member field derived type
type member_field_type
   real(kind_real),allocatable :: fld(:,:,:,:)    ! Ensemble perturbation
end type member_field_type

! Ensemble derived type
type ens_type
   ! Attributes
   integer :: ne                                  ! Ensemble size
   integer :: nsub                                ! Number of sub-ensembles
   logical :: allocated                           ! Allocation flag

   ! Data
   type(member_field_type),allocatable :: mem(:)  ! Members
   type(member_field_type),allocatable :: mean(:) ! Ensemble mean
contains
   procedure :: set_att => ens_set_att
   procedure :: alloc => ens_alloc
   procedure :: dealloc => ens_dealloc
   procedure :: copy => ens_copy
   procedure :: remove_mean => ens_remove_mean
   procedure :: apply_bens => ens_apply_bens
   procedure :: cortrack => ens_cortrack
   procedure :: corstats => ens_corstats
end type ens_type

private
public :: ens_type

contains

!----------------------------------------------------------------------
! Subroutine: ens_set_att
! Purpose: set attributes
!----------------------------------------------------------------------
subroutine ens_set_att(ens,ne,nsub)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble
integer,intent(in) :: ne             ! Ensemble size
integer,intent(in) :: nsub           ! Number of sub-ensembles

! Copy attributes
ens%ne = ne
ens%nsub = nsub
ens%allocated = .false.

end subroutine ens_set_att

!----------------------------------------------------------------------
! Subroutine: ens_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine ens_alloc(ens,nam,geom,ne,nsub)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
integer,intent(in) :: ne             ! Ensemble size
integer,intent(in) :: nsub           ! Number of sub-ensembles

! Local variables
integer :: isub

! Copy attributes
call ens%set_att(ne,nsub)

! Allocation
if (ne>0) then
   allocate(ens%mem(ne))
   allocate(ens%mean(nsub))
   do isub=1,nsub
      allocate(ens%mean(isub)%fld(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   end do
   ens%allocated = .true.
end if

end subroutine ens_alloc

!----------------------------------------------------------------------
! Subroutine: ens_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine ens_dealloc(ens)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble

! Local variables
integer :: ie,isub

! Release memory
if (allocated(ens%mem)) then
   do ie=1,ens%ne
      if (allocated(ens%mem(ie)%fld)) deallocate(ens%mem(ie)%fld)
   end do
   deallocate(ens%mem)
end if
if (allocated(ens%mean)) then
   do isub=1,ens%nsub
      if (allocated(ens%mean(isub)%fld)) deallocate(ens%mean(isub)%fld)
   end do
   deallocate(ens%mean)
end if
ens%allocated = .false.

end subroutine ens_dealloc

!----------------------------------------------------------------------
! Subroutine: ens_copy
! Purpose: copy
!----------------------------------------------------------------------
subroutine ens_copy(ens_out,ens_in)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens_out ! Output ensemble
type(ens_type),intent(in) :: ens_in      ! Input ensemble

! Local variables
integer :: ie,isub

! Copy data
if (allocated(ens_in%mem)) then
   do ie=1,ens_in%ne
      if (.not.allocated(ens_out%mem(ie)%fld)) allocate(ens_out%mem(ie)%fld(size(ens_in%mem(ie)%fld,1), &
                                                                          & size(ens_in%mem(ie)%fld,2), &
                                                                          & size(ens_in%mem(ie)%fld,3), &
                                                                          & size(ens_in%mem(ie)%fld,4)))
      ens_out%mem(ie)%fld = ens_in%mem(ie)%fld
   end do
end if
if (allocated(ens_in%mean)) then
   do isub=1,ens_in%nsub
      ens_out%mean(isub)%fld = ens_in%mean(isub)%fld
   end do
end if

end subroutine ens_copy

!----------------------------------------------------------------------
! Subroutine: ens_remove_mean
! Purpose: remove ensemble mean
!----------------------------------------------------------------------
subroutine ens_remove_mean(ens)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble

! Local variables
integer :: isub,ie_sub,ie

if (ens%allocated) then
   ! Loop over sub-ensembles
   do isub=1,ens%nsub
      ! Compute mean
      ens%mean(isub)%fld = 0.0
      do ie_sub=1,ens%ne/ens%nsub
         ie = ie_sub+(isub-1)*ens%ne/ens%nsub
         ens%mean(isub)%fld = ens%mean(isub)%fld+ens%mem(ie)%fld
      end do
      ens%mean(isub)%fld = ens%mean(isub)%fld/(ens%ne/ens%nsub)

      ! Remove mean
      do ie_sub=1,ens%ne/ens%nsub
         ie = ie_sub+(isub-1)*ens%ne/ens%nsub
         ens%mem(ie)%fld = ens%mem(ie)%fld-ens%mean(isub)%fld
      end do
   end do
end if

end subroutine ens_remove_mean

!----------------------------------------------------------------------
! Subroutine: ens_apply_bens
! Purpose: apply raw ensemble covariance
!----------------------------------------------------------------------
subroutine ens_apply_bens(ens,mpl,nam,geom,fld)

implicit none

! Passed variables
class(ens_type),intent(in) :: ens                                       ! Ensemble
type(mpl_type),intent(inout) :: mpl                                     ! MPI data
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Field

! Local variable
integer :: ie,ic0a,il0,iv,its
real(kind_real) :: alpha,norm
real(kind_real) :: fld_copy(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: pert(geom%nc0a,geom%nl0,nam%nv,nam%nts)

! Initialization
fld_copy = fld

! Apply localized ensemble covariance formula
fld = 0.0
norm = sqrt(real(nam%ens1_ne-1,kind_real))
do ie=1,nam%ens1_ne
   ! Compute perturbation
   !$omp parallel do schedule(static) private(its,iv,il0,ic0a)
   do its=1,nam%nts
      do iv=1,nam%nv
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask_c0a(ic0a,il0)) then
                  pert(ic0a,il0,iv,its) = ens%mem(ie)%fld(ic0a,il0,iv,its)/norm
               else
                  pert(ic0a,il0,iv,its) = mpl%msv%valr
               end if
            end do
         end do
      end do
   end do
   !$omp end parallel do

   ! Dot product
   call mpl%dot_prod(pert,fld_copy,alpha)

   ! Schur product
   !$omp parallel do schedule(static) private(its,iv,il0,ic0a)
   do its=1,nam%nts
      do iv=1,nam%nv
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%mask_c0a(ic0a,il0)) fld(ic0a,il0,iv,its) = fld(ic0a,il0,iv,its)+alpha*pert(ic0a,il0,iv,its)
            end do
         end do
      end do
   end do
   !$omp end parallel do
end do

end subroutine ens_apply_bens

!----------------------------------------------------------------------
! Subroutine: ens_cortrack
! Purpose: correlation tracker
!----------------------------------------------------------------------
subroutine ens_cortrack(ens,mpl,nam,geom,io)

implicit none

! Passed variables
class(ens_type),intent(in) :: ens   ! Ensemble
type(mpl_type),intent(inout) :: mpl ! MPI data
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry
type(io_type),intent(in) :: io      ! I/O

! Local variable
integer :: ic0a,ic0,il0,ie,its,iproc(1),ind(2)
real(kind_real) :: proc_to_val(mpl%nproc),val,var_dirac
real(kind_real) :: var(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: dirac(geom%nc0a,geom%nl0,nam%nv,nam%nts),cor(geom%nc0a,geom%nl0,nam%nv,nam%nts)
character(len=2) :: timeslotchar
character(len=1024) :: filename

! File name
filename = trim(nam%prefix)//'_cortrack'
call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)

! Compute variance
write(mpl%info,'(a7,a)') '','Compute variance'
call mpl%flush
var = 0.0
do ie=1,ens%ne
   var = var+ens%mem(ie)%fld**2
end do
var = var/real(ens%ne-ens%nsub,kind_real)

! Generate dirac field
dirac = 0.0
if (geom%iprocdir(1)==mpl%myproc) dirac(geom%ic0adir(1),geom%il0dir(1),1,geom%itsdir(1)) = 1.0

! Apply raw ensemble covariance
write(mpl%info,'(a10,a)') '','Apply raw ensemble covariance'
call mpl%flush
cor = dirac
call ens%apply_bens(mpl,nam,geom,cor)

! Normalize correlation
write(mpl%info,'(a10,a)') '','Normalize correlation'
call mpl%flush
call mpl%f_comm%allreduce(sum(dirac*var),var_dirac,fckit_mpi_sum())
cor = cor/sqrt(var*var_dirac)

! Correlation maximum displacement
write(mpl%info,'(a10,a)') '','Correlation maximum displacement'
call mpl%flush
do its=1,nam%nts
   val = maxval(cor(:,:,1,its))
   call mpl%f_comm%allgather(val,proc_to_val)
   iproc = maxloc(proc_to_val)
   if (mpl%myproc==iproc(1)) then
      ind = maxloc(cor(:,:,1,its))
      ic0a = ind(1)
      il0 = ind(2)
      ic0 = geom%c0a_to_c0(ic0a)
   end if
   call mpl%f_comm%broadcast(ic0,iproc(1)-1)
   call mpl%f_comm%broadcast(il0,iproc(1)-1)
   write(mpl%info,'(a13,a,i2,a,f6.1,a,f6.1,a,i3,a,f6.2)') '','Timeslot ',nam%timeslot(its),' ~> lon / lat / lev / val: ', &
 & geom%lon(ic0)*rad2deg,' / ',geom%lat(ic0)*rad2deg,' / ',nam%levs(il0),' / ',proc_to_val(iproc(1))
   call mpl%flush
end do

! Write correlation
write(mpl%info,'(a10,a)') '','Write correlation'
call mpl%flush
do its=1,nam%nts
   write(timeslotchar,'(i2.2)') nam%timeslot(its)
   call io%fld_write(mpl,nam,geom,filename,'cor_'//timeslotchar,cor(:,:,:,its))
end do

! Correlation tracker
write(mpl%info,'(a7,a)') '','Correlation tracker'
call mpl%flush
write(timeslotchar,'(i2.2)') nam%timeslot(1)
call io%fld_write(mpl,nam,geom,filename,'tracker_'//timeslotchar//'_0',cor(:,:,:,1))
call io%fld_write(mpl,nam,geom,filename,'tracker_'//timeslotchar//'_1',cor(:,:,:,2))
do its=2,nam%nts-1
   ! Correlation maximum displacement
   ic0a = mpl%msv%vali
   val = maxval(cor(:,:,1,its))
   call mpl%f_comm%allgather(val,proc_to_val)
   iproc = maxloc(proc_to_val)
   if (mpl%myproc==iproc(1)) then
      ind = maxloc(cor(:,:,1,its))
      ic0a = ind(1)
      il0 = ind(2)
      ic0 = geom%c0a_to_c0(ic0a)
   end if
   call mpl%f_comm%broadcast(ic0,iproc(1)-1)
   call mpl%f_comm%broadcast(il0,iproc(1)-1)
   write(mpl%info,'(a10,a,i2,a,f6.1,a,f6.1,a,i3,a,f6.2)') '','Timeslot ',nam%timeslot(its),' ~> lon / lat / lev / val: ', &
 & geom%lon(ic0)*rad2deg,' / ',geom%lat(ic0)*rad2deg,' / ',nam%levs(il0),' / ',proc_to_val(iproc(1))
   call mpl%flush

   ! Generate dirac field
   dirac = 0.0
   if (iproc(1)==mpl%myproc) dirac(ic0a,il0,1,its) = 1.0

   ! Apply raw ensemble covariance
   write(mpl%info,'(a13,a)') '','Apply raw ensemble covariance'
   call mpl%flush
   cor = dirac
   call ens%apply_bens(mpl,nam,geom,cor)

   ! Normalize correlation tracker
   write(mpl%info,'(a13,a)') '','Normalize correlation tracker'
   call mpl%flush
   call mpl%f_comm%allreduce(sum(dirac*var),var_dirac,fckit_mpi_sum())
   cor = cor/sqrt(var*var_dirac)

   ! Write correlation tracker
   write(mpl%info,'(a13,a)') '','Write correlation tracker'
   call mpl%flush
   write(timeslotchar,'(i2.2)') nam%timeslot(its)
   call io%fld_write(mpl,nam,geom,filename,'tracker_'//timeslotchar//'_0',cor(:,:,:,its))
   call io%fld_write(mpl,nam,geom,filename,'tracker_'//timeslotchar//'_1',cor(:,:,:,its+1))
end do

end subroutine ens_cortrack

!----------------------------------------------------------------------
! Subroutine: ens_corstats
! Purpose: correlation tracker
!----------------------------------------------------------------------
subroutine ens_corstats(ens,mpl,rng,nam,geom)

implicit none

! Passed variables
class(ens_type),intent(in) :: ens   ! Ensemble
type(mpl_type),intent(inout) :: mpl ! MPI data
type(rng_type),intent(inout) :: rng ! Random number generator
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry

! Local variable
integer,parameter :: ntest = 1000
integer :: ie,il0,itest,ic0dir,iprocdir,ic0adir,its,iv
integer :: ncid,ntest_id,nl0_id,nv_id,nts_id,cor_max_id,cor_max_avg_id,cor_max_std_id
real(kind_real) :: var_dirac,var(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: dirac(geom%nc0a,geom%nl0,nam%nv,nam%nts),cor(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: cor_max(ntest,geom%nl0,nam%nv,nam%nts),cor_max_avg(geom%nl0,nam%nv,nam%nts),cor_max_std(geom%nl0,nam%nv,nam%nts)
character(len=1024) :: filename
character(len=1024) :: subr = 'ens_corstats'

! Compute variance
write(mpl%info,'(a7,a)') '','Compute variance'
call mpl%flush
var = 0.0
do ie=1,ens%ne
   var = var+ens%mem(ie)%fld**2
end do
var = var/real(ens%ne-ens%nsub,kind_real)

! Compute correlation maximum statistics
write(mpl%info,'(a7,a)') '','Compute correlation maximum statistics'
call mpl%flush
do il0=1,geom%nl0
   write(mpl%info,'(a10,a,i3,a)') '','Level ',il0,': '
   call mpl%flush(.false.)

   ! Initialization
   itest = 1
   call mpl%prog_init(ntest)

   do while (itest<=ntest)
      ! Generate random dirac point
      if (mpl%main) call rng%rand_integer(1,geom%nc0,ic0dir)
      call mpl%f_comm%broadcast(ic0dir,mpl%rootproc-1)

      if (geom%mask_c0(ic0dir,il0)) then
         ! Get processor and local index
         iprocdir = geom%c0_to_proc(ic0dir)
         if (iprocdir==mpl%myproc) ic0adir = geom%c0_to_c0a(ic0dir)

         do iv=1,nam%nv
            ! Generate dirac field
            dirac = 0.0
            if (iprocdir==mpl%myproc) dirac(ic0adir,il0,iv,nam%nts) = 1.0

            ! Apply raw ensemble covariance
            cor = dirac
            call ens%apply_bens(mpl,nam,geom,cor)

            ! Normalize correlation
            call mpl%f_comm%allreduce(sum(dirac*var),var_dirac,fckit_mpi_sum())
            cor = cor/sqrt(var*var_dirac)

            ! Save correlation maximum
            do its=1,nam%nts
               call mpl%f_comm%allreduce(maxval(cor(:,il0,iv,its)),cor_max(itest,il0,iv,its),fckit_mpi_max())
            end do
         end do

         ! Update
         call mpl%prog_print(itest)
         itest = itest+1
      end if
   end do
   call mpl%prog_final

   ! Compute average and standard-deviation
   do its=1,nam%nts
      do iv=1,nam%nv
         cor_max_avg(il0,iv,its) = sum(cor_max(:,il0,iv,its))/real(ntest,kind_real)
         cor_max_std(il0,iv,its) = sqrt(sum((cor_max(:,il0,iv,its)-cor_max_avg(il0,iv,its))**2)/real(ntest-1,kind_real))
      end do
   end do
end do

if (mpl%main) then
   ! Create file
   filename = trim(nam%prefix)//'_corstats'
   call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

   ! Write namelist parameters
   call nam%write(mpl,ncid)

   ! Define dimensions
   ntest_id = mpl%ncdimcheck(subr,ncid,'ntest',ntest,.true.)
   nl0_id = mpl%ncdimcheck(subr,ncid,'nl0',geom%nl0,.true.)
   nv_id = mpl%ncdimcheck(subr,ncid,'nv',nam%nv,.true.)
   nts_id = mpl%ncdimcheck(subr,ncid,'nts',nam%nts,.true.)

   ! Define variables
   call mpl%ncerr(subr,nf90_def_var(ncid,'cor_max',nc_kind_real,(/ntest_id,nl0_id,nv_id,nts_id/),cor_max_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,cor_max_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'cor_max_avg',nc_kind_real,(/nl0_id,nv_id,nts_id/),cor_max_avg_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,cor_max_avg_id,'_FillValue',mpl%msv%valr))
   call mpl%ncerr(subr,nf90_def_var(ncid,'cor_max_std',nc_kind_real,(/nl0_id,nv_id,nts_id/),cor_max_std_id))
   call mpl%ncerr(subr,nf90_put_att(ncid,cor_max_std_id,'_FillValue',mpl%msv%valr))


   ! End definition mode
   call mpl%ncerr(subr,nf90_enddef(ncid))

   ! Write variables
   call mpl%ncerr(subr,nf90_put_var(ncid,cor_max_id,cor_max))
   call mpl%ncerr(subr,nf90_put_var(ncid,cor_max_avg_id,cor_max_avg))
   call mpl%ncerr(subr,nf90_put_var(ncid,cor_max_std_id,cor_max_std))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))
end if

end subroutine ens_corstats

end module type_ens
