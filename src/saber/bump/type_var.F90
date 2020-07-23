!----------------------------------------------------------------------
! Module: type_var
! Purpose: variance derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_var

use fckit_mpi_module, only: fckit_mpi_sum
!$ use omp_lib
use tools_const, only: reqkm
use tools_kinds, only: kind_real,huge_real
use type_cmat_blk, only: cmat_blk_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_nicas_blk, only: nicas_blk_type
use type_rng, only: rng_type

implicit none

! Variance derived type
type var_type
   integer :: ne                                  ! Ensemble size
   real(kind_real),allocatable :: m2(:,:,:,:)     ! Variance
   real(kind_real),allocatable :: m4(:,:,:,:)     ! Fourth-order centered moment
   real(kind_real),allocatable :: m2flt(:,:,:,:)  ! Filtered variance
   real(kind_real),allocatable :: m2sqrt(:,:,:,:) ! Variance square-root
contains
   procedure :: alloc => var_alloc
   procedure :: partial_dealloc => var_partial_dealloc
   procedure :: dealloc => var_dealloc
   procedure :: read => var_read
   procedure :: write => var_write
   procedure :: run_var => var_run_var
   procedure :: filter => var_filter
   procedure :: apply_sqrt => var_apply_sqrt
   procedure :: apply_sqrt_inv => var_apply_sqrt_inv
end type var_type

private
public :: var_type

contains

!----------------------------------------------------------------------
! Subroutine: var_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine var_alloc(var,nam,geom)

implicit none

! Passed variables
class(var_type),intent(inout) :: var ! Variance
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry

! Allocation
allocate(var%m2(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(var%m4(geom%nc0a,geom%nl0,nam%nv,nam%nts))
if (nam%var_filter) allocate(var%m2flt(geom%nc0a,geom%nl0,nam%nv,nam%nts))
allocate(var%m2sqrt(geom%nc0a,geom%nl0,nam%nv,nam%nts))

end subroutine var_alloc

!----------------------------------------------------------------------
! Subroutine: var_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine var_partial_dealloc(var)

implicit none

! Passed variables
class(var_type),intent(inout) :: var ! Variance

! Release memory
if (allocated(var%m2)) deallocate(var%m2)
if (allocated(var%m4)) deallocate(var%m4)
if (allocated(var%m2flt)) deallocate(var%m2flt)

end subroutine var_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: var_dealloc
! Purpose: release memory (full)
!----------------------------------------------------------------------
subroutine var_dealloc(var)

implicit none

! Passed variables
class(var_type),intent(inout) :: var ! Variance

! Release memory
call var%partial_dealloc
if (allocated(var%m2sqrt)) deallocate(var%m2sqrt)

end subroutine var_dealloc

!----------------------------------------------------------------------
! Subroutine: var_read
! Purpose: read
!----------------------------------------------------------------------
subroutine var_read(var,mpl,nam,geom,io)

implicit none

! Passed variables
class(var_type),intent(inout) :: var ! Variance
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(io_type),intent(in) :: io       ! I/O

! Local variables
integer :: iv,its
character(len=1024) :: filename,grpname

! Read variance
write(mpl%info,'(a7,a)') '','Read variance'
call mpl%flush

! Allocation
call var%alloc(nam,geom)

! Set filename
filename = trim(nam%prefix)//'_var'

! Read raw variance, fourth-order moment, filtered variance and standard-deviation
do its=1,nam%nts
   do iv=1,nam%nv
      grpname = trim(nam%variables(iv))//'_'//trim(nam%timeslots(its))
      call io%fld_read(mpl,nam,geom,filename,'m2',var%m2(:,:,iv,its),trim(grpname))
      call io%fld_read(mpl,nam,geom,filename,'m4',var%m4(:,:,iv,its),trim(grpname))
      if (nam%var_filter) call io%fld_read(mpl,nam,geom,filename,'m2flt',var%m2flt(:,:,iv,its),trim(grpname))
      call io%fld_read(mpl,nam,geom,filename,'m2sqrt',var%m2sqrt(:,:,iv,its),trim(grpname))
   end do
end do

end subroutine var_read

!----------------------------------------------------------------------
! Subroutine: var_write
! Purpose: write
!----------------------------------------------------------------------
subroutine var_write(var,mpl,nam,geom,io)

implicit none

! Passed variables
class(var_type),intent(inout) :: var ! Variance
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(io_type),intent(in) :: io       ! I/O

! Local variables
integer :: iv,its
character(len=1024) :: filename,grpname

! Write variance
write(mpl%info,'(a7,a)') '','Write variance'
call mpl%flush

! Set filename
filename = trim(nam%prefix)//'_var'

! Write vertical unit
call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)

! Write raw variance, fourth-order moment, filtered variance and standard-deviation
do its=1,nam%nts
   do iv=1,nam%nv
      grpname = trim(nam%variables(iv))//'_'//trim(nam%timeslots(its))
      call io%fld_write(mpl,nam,geom,filename,'m2',var%m2(:,:,iv,its),trim(grpname))
      call io%fld_write(mpl,nam,geom,filename,'m4',var%m4(:,:,iv,its),trim(grpname))
      if (nam%var_filter) call io%fld_write(mpl,nam,geom,filename,'m2flt',var%m2flt(:,:,iv,its),trim(grpname))
      call io%fld_write(mpl,nam,geom,filename,'m2sqrt',var%m2sqrt(:,:,iv,its),trim(grpname))
   end do
end do

end subroutine var_write

!----------------------------------------------------------------------
! Subroutine: var_run_var
! Purpose: compute variance
!----------------------------------------------------------------------
subroutine var_run_var(var,mpl,rng,nam,geom,ens,io)

implicit none

! Passed variables
class(var_type),intent(inout) :: var ! Variance
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(rng_type),intent(inout) :: rng  ! Random number generator
type(nam_type),intent(inout) :: nam  ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(ens_type), intent(in) :: ens    ! Ensemble
type(io_type),intent(in) :: io       ! I/O

! Local variables
integer :: ie,ic0a,il0
real(kind_real) :: norm_m2,norm_m4

! Allocation
call var%alloc(nam,geom)

! Initialization
norm_m2 = 1.0/real(ens%ne-1,kind_real)
norm_m4 = 1.0/real(ens%ne,kind_real)
var%ne = ens%ne
var%m2 = 0.0
var%m4 = 0.0

! Compute variance
write(mpl%info,'(a7,a)') '','Compute variance'
call mpl%flush
do ie=1,ens%ne
   var%m2 = var%m2+ens%mem(ie)%fld**2
   var%m4 = var%m4+ens%mem(ie)%fld**4
end do
var%m2 = var%m2*norm_m2
var%m4 = var%m4*norm_m4

! Apply mask
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      if (.not.geom%gmask_c0a(ic0a,il0)) then
         var%m2(ic0a,il0,:,:) = mpl%msv%valr
         var%m4(ic0a,il0,:,:) = mpl%msv%valr
      end if
   end do
end do

if (nam%var_filter) then
   ! Filter variance
   call var%filter(mpl,rng,nam,geom)

   ! Take square-root
   var%m2sqrt = sqrt(var%m2flt)
else
   ! Take square-root
   var%m2sqrt = sqrt(var%m2)
end if

! Write variance
call var%write(mpl,nam,geom,io)

end subroutine var_run_var

!----------------------------------------------------------------------
! Subroutine: var_filter
! Purpose: filter variance
!----------------------------------------------------------------------
subroutine var_filter(var,mpl,rng,nam,geom)

implicit none

! Passed variables
class(var_type),intent(inout) :: var ! Variance
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(rng_type),intent(inout) :: rng  ! Random number generator
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry

! Local variables
integer :: n,its,iv,il0,iter
real(kind_real) :: P9,P20,P21,diff,diff_abs_min
real(kind_real) :: m2sq(geom%nl0),m2sq_tot(geom%nl0),m4(geom%nl0),m4_tot(geom%nl0),m2sqasy(geom%nl0)
real(kind_real) :: rhflt(geom%nl0),drhflt(geom%nl0),m2prod(geom%nl0),m2prod_tot(geom%nl0)
real(kind_real) :: m2_ini(geom%nc0a,geom%nl0),m2(geom%nc0a,geom%nl0)
logical :: dichotomy(geom%nl0),convergence(geom%nl0)
type(nicas_blk_type) :: nicas_blk

write(mpl%info,'(a7,a)') '','Filter variance'
call mpl%flush

! Ensemble size-dependent coefficients
n = var%ne
P9 = -real(n,kind_real)/real((n-2)*(n-3),kind_real)
P20 = real((n-1)*(n**2-3*n+3),kind_real)/real(n*(n-2)*(n-3),kind_real)
P21 = real(n-1,kind_real)/real(n+1,kind_real)

do its=1,nam%nts
   write(mpl%info,'(a10,a,a)') '','Timeslot ',trim(nam%timeslots(its))
   call mpl%flush

   do iv=1,nam%nv
      write(mpl%info,'(a13,a,a)') '','Variable ',trim(nam%variables(iv))
      call mpl%flush

      ! Global sum
      do il0=1,geom%nl0
         m2sq(il0) = sum(var%m2(:,il0,iv,its)**2,mask=geom%gmask_c0a(:,il0))
         m4(il0) = sum(var%m4(:,il0,iv,its),mask=geom%gmask_c0a(:,il0))
      end do
      call mpl%f_comm%allreduce(m2sq,m2sq_tot,fckit_mpi_sum())
      call mpl%f_comm%allreduce(m4,m4_tot,fckit_mpi_sum())

      ! Asymptotic statistics
      if (nam%gau_approx) then
         ! Gaussian approximation
         m2sqasy = P21*m2sq_tot
      else
         ! General case
         m2sqasy = P20*m2sq_tot+P9*m4_tot
      end if

      ! Dichotomy initialization
      m2_ini = var%m2(:,:,iv,its)
      convergence = .true.
      dichotomy = .false.
      rhflt = nam%var_rhflt
      drhflt = rhflt
      diff_abs_min = huge_real

      do iter=1,nam%var_niter
         ! Copy initial value
         m2 = m2_ini

         ! Set smoother parameters
         call nicas_blk%compute_parameters(mpl,rng,nam,geom,rhflt)

         ! Apply smoother
         call nicas_blk%apply(mpl,geom,m2)

         ! Global product
         do il0=1,geom%nl0
            m2prod(il0) = sum(m2(:,il0)*m2_ini(:,il0),mask=geom%gmask_c0a(:,il0))
         end do
         call mpl%f_comm%allreduce(m2prod,m2prod_tot,fckit_mpi_sum())

         ! Print results
         write(mpl%info,'(a16,a,i2)') '','Iteration ',iter
         call mpl%flush
         do il0=1,geom%nl0
            if (m2sqasy(il0)>0.0) then
               write(mpl%info,'(a19,a,i3,a,f10.2,a,e12.5)') '','Level ',il0,': rhflt = ',rhflt(il0)*reqkm,' km, rel. diff. = ', &
 & (m2prod_tot(il0)-m2sqasy(il0))/m2sqasy(il0)
               call mpl%flush
            end if
         end do

         ! Update support radius
         do il0=1,geom%nl0
            diff = m2prod_tot(il0)-m2sqasy(il0)
            if (diff>0.0) then
               ! Increase filtering support radius
               if (dichotomy(il0)) then
                  drhflt(il0) = 0.5*drhflt(il0)
                  rhflt(il0) = rhflt(il0)+drhflt(il0)
               else
                  convergence(il0) = .false.
                  rhflt(il0) = rhflt(il0)+drhflt(il0)
                  drhflt(il0) = 2.0*drhflt(il0)
               end if
            else
               ! Convergence
               convergence(il0) = .true.

               ! Change dichotomy status
               if (.not.dichotomy(il0)) then
                  dichotomy(il0) = .true.
                  drhflt(il0) = 0.5*drhflt(il0)
               end if

               ! Decrease filtering support radius
               drhflt(il0) = 0.5*drhflt(il0)
               rhflt(il0) = rhflt(il0)-drhflt(il0)
            end if
            if ((abs(diff)<diff_abs_min).or.(nam%var_niter==1)) then
               ! Copy best result
               diff_abs_min = abs(diff)
               var%m2flt(:,:,iv,its) = m2
            end if
         end do

         ! Release memory
         call nicas_blk%dealloc
      end do
   end do
end do

end subroutine var_filter

!----------------------------------------------------------------------
! Subroutine: var_apply_sqrt
! Purpose: apply square-root variance
!----------------------------------------------------------------------
subroutine var_apply_sqrt(var,nam,geom,fld)

implicit none

! Passed variables
class(var_type),intent(in) :: var                                       ! Variance
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Source/destination vector

! Local variables
integer :: ic0a,il0

! Apply variance
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      if (geom%gmask_c0a(ic0a,il0)) fld(ic0a,il0,:,:) = fld(ic0a,il0,:,:)*var%m2sqrt(ic0a,il0,:,:)
   end do
end do

end subroutine var_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: var_apply_sqrt_inv
! Purpose: apply square-root variance inverse
!----------------------------------------------------------------------
subroutine var_apply_sqrt_inv(var,nam,geom,fld)

implicit none

! Passed variables
class(var_type),intent(in) :: var                                       ! Variance
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Source/destination vector

! Local variables
integer :: ic0a,il0

! Apply inverse variance
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      if (geom%gmask_c0a(ic0a,il0)) fld(ic0a,il0,:,:) = fld(ic0a,il0,:,:)/var%m2sqrt(ic0a,il0,:,:)
   end do
end do

end subroutine var_apply_sqrt_inv

end module type_var
