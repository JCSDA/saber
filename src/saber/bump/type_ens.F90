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
use tools_const, only: deg2rad,rad2deg,req
use tools_func, only: sphere_dist,lonlat2xyz,xyz2lonlat
use tools_kinds, only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use type_geom, only: geom_type
use type_io, only: io_type
use type_linop, only: linop_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

! Member field derived type
type member_field_type
   real(kind_real),allocatable :: fld(:,:,:)      ! Ensemble perturbation
end type member_field_type

! Ensemble derived type
type ens_type
   ! Attributes
   integer :: ne                                  ! Ensemble size
   integer :: nsub                                ! Number of sub-ensembles
   logical :: allocated = .false.                 ! Allocation flag

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
   procedure :: apply_bens_dirac => ens_apply_bens_dirac
   procedure :: normality => ens_normality
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
      allocate(ens%mean(isub)%fld(geom%nc0a,geom%nl0,nam%nv))
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
 & size(ens_in%mem(ie)%fld,2),size(ens_in%mem(ie)%fld,3)))
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
class(ens_type),intent(in) :: ens                               ! Ensemble
type(mpl_type),intent(inout) :: mpl                             ! MPI data
type(nam_type),intent(in) :: nam                                ! Namelist
type(geom_type),intent(in) :: geom                              ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv) ! Field

! Local variable
integer :: ie,ic0a,il0,iv
real(kind_real) :: alpha,norm
real(kind_real) :: fld_copy(geom%nc0a,geom%nl0,nam%nv)
real(kind_real) :: pert(geom%nc0a,geom%nl0,nam%nv)

! Initialization
fld_copy = fld

! Apply ensemble covariance formula
fld = 0.0
norm = 1.0/real(ens%ne-1,kind_real)
do ie=1,ens%ne
   ! Set perturbation
   !$omp parallel do schedule(static) private(iv,il0,ic0a)
   do iv=1,nam%nv
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%gmask_c0a(ic0a,il0)) then
               pert(ic0a,il0,iv) = ens%mem(ie)%fld(ic0a,il0,iv)
            else
               pert(ic0a,il0,iv) = mpl%msv%valr
            end if
          end do
      end do
   end do
   !$omp end parallel do

   ! Dot product
   call mpl%dot_prod(pert,fld_copy,alpha)

   ! Schur product
   !$omp parallel do schedule(static) private(iv,il0,ic0a)
   do iv=1,nam%nv
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%gmask_c0a(ic0a,il0)) fld(ic0a,il0,iv) = fld(ic0a,il0,iv)+alpha*pert(ic0a,il0,iv)*norm
         end do
      end do
   end do
   !$omp end parallel do
end do

end subroutine ens_apply_bens

!----------------------------------------------------------------------
! Subroutine: ens_apply_bens_dirac
! Purpose: apply raw ensemble covariance to a Dirac (faster formulation)
!----------------------------------------------------------------------
subroutine ens_apply_bens_dirac(ens,mpl,nam,geom,iprocdir,ic0adir,il0dir,ivdir,fld)

implicit none

! Passed variables
class(ens_type),intent(in) :: ens                             ! Ensemble
type(mpl_type),intent(inout) :: mpl                           ! MPI data
type(nam_type),intent(in) :: nam                              ! Namelist
type(geom_type),intent(in) :: geom                            ! Geometry
integer,intent(in) :: iprocdir                                ! Processor index for dirac function
integer,intent(in) :: ic0adir                                 ! Subset Sc0, halo A index for dirac function
integer,intent(in) :: il0dir                                  ! Subset Sl0 index for dirac function
integer,intent(in) :: ivdir                                   ! Variable index for dirac function
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv) ! Field

! Local variable
integer :: ie,ic0a,il0,iv
real(kind_real) :: alpha(ens%ne),norm

! Apply ensemble covariance formula for a Dirac function
norm = 1.0/real(ens%ne-1,kind_real)
if (mpl%myproc==iprocdir) then
   do ie=1,ens%ne
      alpha(ie) = ens%mem(ie)%fld(ic0adir,il0dir,ivdir)
   end do
end if
call mpl%f_comm%broadcast(alpha,iprocdir-1)
fld = 0.0
do ie=1,ens%ne
   !$omp parallel do schedule(static) private(iv,il0,ic0a)
   do iv=1,nam%nv
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%gmask_c0a(ic0a,il0)) then
               fld(ic0a,il0,iv) = fld(ic0a,il0,iv)+alpha(ie)*ens%mem(ie)%fld(ic0a,il0,iv)*norm
            else
               fld(ic0a,il0,iv) = mpl%msv%valr
            end if
         end do
      end do
   end do
   !$omp end parallel do
end do

end subroutine ens_apply_bens_dirac

!----------------------------------------------------------------------
! Subroutine: ens_normality
! Purpose: perform some normality diagnostics
!----------------------------------------------------------------------
subroutine ens_normality(ens,mpl,nam,geom,io)

implicit none

! Passed variables
class(ens_type),intent(in) :: ens   ! Ensemble
type(mpl_type),intent(inout) :: mpl ! MPI data
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry
type(io_type),intent(in) :: io      ! I/O

! Local variables
integer :: ncid,nloc_id,ne_id,nem1_id,ic0a_id,il0_id,iv_id,order_id,ens_norm_id,ens_step_id
integer :: iv,il0,ic0a,ie,nloc,iloc,nglb
integer,allocatable :: ic0a_loc(:),il0_loc(:),iv_loc(:),order(:,:)
real(kind_real) :: norm_m2,norm_m4,norm
real(kind_real) :: m2(geom%nc0a,geom%nl0,nam%nv),m4(geom%nc0a,geom%nl0,nam%nv),kurt(geom%nc0a,geom%nl0,nam%nv)
real(kind_real),allocatable :: ens_loc(:),ens_norm(:,:),ens_step(:,:)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'ens_normality'

! Set file name
filename = trim(nam%prefix)//'_umf'

! Write vertical unit
call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)

! Initialization
norm_m2 = 1.0/real(ens%ne-1,kind_real)
norm_m4 = 1.0/real(ens%ne,kind_real)

! Compute variance and kurtosis
write(mpl%info,'(a7,a)') '','Compute variance and kurtosis'
call mpl%flush
m2 = 0.0
m4 = 0.0
do ie=1,ens%ne
   m2 = m2+ens%mem(ie)%fld**2
   m4 = m4+ens%mem(ie)%fld**4
end do
m2 = m2*norm_m2
m4 = m4*norm_m4
kurt = mpl%msv%valr
do iv=1,nam%nv
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         if (m2(ic0a,il0,iv)>0.0) kurt(ic0a,il0,iv) = m4(ic0a,il0,iv)/m2(ic0a,il0,iv)**2
      end do
   end do
end do
call io%fld_write(mpl,nam,geom,filename,'m2',m2)
call io%fld_write(mpl,nam,geom,filename,'m4',m4)
call io%fld_write(mpl,nam,geom,filename,'kurt',kurt)

! Allocation
nloc = count(mpl%msv%isnot(kurt).and.(kurt>nam%gen_kurt_th))
allocate(ic0a_loc(nloc))
allocate(il0_loc(nloc))
allocate(iv_loc(nloc))
allocate(order(ens%ne,nloc))
allocate(ens_loc(ens%ne))
allocate(ens_norm(ens%ne,nloc))
allocate(ens_step(ens%ne-1,nloc))
call mpl%f_comm%allreduce(nloc,nglb,fckit_mpi_sum())

! Save ensemble
write(mpl%info,'(a7,a,i6,a,i6,a)') '','Save ensemble for ',nloc,' points (',nglb,' total)'
call mpl%flush
iloc = 0
do iv=1,nam%nv
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         if (mpl%msv%isnot(kurt(ic0a,il0,iv)).and.(kurt(ic0a,il0,iv)>nam%gen_kurt_th)) then
            ! Update index
            iloc = iloc+1

            ! Copy data
            ic0a_loc(iloc) = ic0a
            il0_loc(iloc) = il0
            iv_loc(iloc) = iv
            do ie=1,ens%ne
               ens_loc(ie) = ens%mem(ie)%fld(ic0a,il0,iv)
            end do

            ! Sort ensemble
            call qsort(ens%ne,ens_loc,order(:,iloc))

            ! Normalize ensemble
            norm = 1.0/(maxval(ens_loc)-minval(ens_loc))
            ens_norm(:,iloc) = (ens_loc-minval(ens_loc))*norm

            ! Compute ensemble steps
            do ie=1,ens%ne-1
               ens_step(ie,iloc) = ens_norm(ie+1,iloc)-ens_norm(ie,iloc)
            end do
         end if
      end do
   end do
end do

! Write ensemble
write(mpl%info,'(a7,a)') '','Write ensemble'
call mpl%flush

! Define file
write(filename,'(a,a,i6.6,a,i6.6)') trim(nam%prefix),'_normality_',mpl%nproc,'-',mpl%myproc
ncid = mpl%nc_file_create_or_open(subr,trim(nam%datadir)//'/'//trim(filename)//'.nc')

! Define dimensions
nloc_id = mpl%nc_dim_define_or_get(subr,ncid,'nloc',nloc)
ne_id = mpl%nc_dim_define_or_get(subr,ncid,'ne',ens%ne)
nem1_id = mpl%nc_dim_define_or_get(subr,ncid,'nem1',ens%ne-1)

if (nloc>0) then
   ! Define variables
   ic0a_id = mpl%nc_var_define_or_get(subr,ncid,'ic0a',nf90_int,(/nloc_id/))
   il0_id = mpl%nc_var_define_or_get(subr,ncid,'il0',nf90_int,(/nloc_id/))
   iv_id = mpl%nc_var_define_or_get(subr,ncid,'iv',nf90_int,(/nloc_id/))
   order_id = mpl%nc_var_define_or_get(subr,ncid,'order',nf90_int,(/ne_id,nloc_id/))
   ens_norm_id = mpl%nc_var_define_or_get(subr,ncid,'ens_norm',nc_kind_real,(/ne_id,nloc_id/))
   ens_step_id = mpl%nc_var_define_or_get(subr,ncid,'ens_step',nc_kind_real,(/nem1_id,nloc_id/))
end if

if (nloc>0) then
   ! Write variables
   call mpl%ncerr(subr,nf90_put_var(ncid,ic0a_id,ic0a_loc))
   call mpl%ncerr(subr,nf90_put_var(ncid,il0_id,il0_loc))
   call mpl%ncerr(subr,nf90_put_var(ncid,iv_id,iv_loc))
   call mpl%ncerr(subr,nf90_put_var(ncid,order_id,order))
   call mpl%ncerr(subr,nf90_put_var(ncid,ens_norm_id,ens_norm))
   call mpl%ncerr(subr,nf90_put_var(ncid,ens_step_id,ens_step))
end if

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

! Release memory
deallocate(ic0a_loc)
deallocate(il0_loc)
deallocate(iv_loc)
deallocate(order)
deallocate(ens_loc)
deallocate(ens_norm)
deallocate(ens_step)

end subroutine ens_normality

end module type_ens
