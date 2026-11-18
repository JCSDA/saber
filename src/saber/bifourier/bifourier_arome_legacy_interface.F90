!----------------------------------------------------------------------
! Module: bifourier_arome_legacy_interface
!> Bifourier AROME legacy Fortran interface
! Author: Benjamin Menetrier
! Copyright 2025 Meteorologisk Institutt
!----------------------------------------------------------------------
module bifourier_arome_legacy_interface

use fckit_configuration_module, only: fckit_configuration
use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double
use bifourier_arome_legacy_mod, only: bifourier_arome_legacy_read_balance, bifourier_arome_legacy_write_balance, &
 & bifourier_arome_legacy_read_covariance, bifourier_arome_legacy_write_covariance

implicit none

private

contains

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_read_balance_c(c_conf,c_attr,nsmaxp1,nflev,kspec2g,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu, &
 & fact1) bind(c,name="bifourier_arome_legacy_read_balance_f90")

implicit none

! Passed variables
type(c_ptr),intent(in),value :: c_conf
type(c_ptr),intent(in),value :: c_attr
integer(c_int),intent(in) :: nsmaxp1
integer(c_int),intent(in) :: nflev
integer(c_int),intent(in) :: kspec2g
real(c_double),intent(inout) :: sdivpb(nsmaxp1*nflev*nflev)
real(c_double),intent(inout) :: stpspb(nsmaxp1*nflev*(nflev+1))
real(c_double),intent(inout) :: stpsdivu(nsmaxp1*nflev*(nflev+1))
real(c_double),intent(inout) :: sqpb(nsmaxp1*nflev*nflev)
real(c_double),intent(inout) :: sqdivu(nsmaxp1*nflev*nflev)
real(c_double),intent(inout) :: sqtpsu(nsmaxp1*(nflev+1)*nflev)
real(c_double),intent(inout) :: fact1(kspec2g)

! Local variables
type(fckit_configuration) :: f_conf
type(fckit_configuration) :: f_attr

! Interface
f_conf = fckit_configuration(c_conf)
f_attr = fckit_configuration(c_attr)

! Call Fortran
call bifourier_arome_legacy_read_balance(f_conf,f_attr,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu,fact1)

! Release memory
call f_conf%final()
call f_attr%final()

end subroutine bifourier_arome_legacy_read_balance_c

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_write_balance_c(c_conf,c_attr,nsmaxp1,nflev,kspec2g,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu, &
 & fact1) bind(c,name="bifourier_arome_legacy_write_balance_f90")

implicit none

! Passed variables
type(c_ptr),intent(in),value :: c_conf
type(c_ptr),intent(in),value :: c_attr
integer(c_int),intent(in) :: nsmaxp1
integer(c_int),intent(in) :: nflev
integer(c_int),intent(in) :: kspec2g
real(c_double),intent(in) :: sdivpb(nsmaxp1*nflev*nflev)
real(c_double),intent(in) :: stpspb(nsmaxp1*nflev*(nflev+1))
real(c_double),intent(in) :: stpsdivu(nsmaxp1*nflev*(nflev+1))
real(c_double),intent(in) :: sqpb(nsmaxp1*nflev*nflev)
real(c_double),intent(in) :: sqdivu(nsmaxp1*nflev*nflev)
real(c_double),intent(in) :: sqtpsu(nsmaxp1*(nflev+1)*nflev)
real(c_double),intent(in) :: fact1(kspec2g)

! Local variables
type(fckit_configuration) :: f_conf
type(fckit_configuration) :: f_attr

! Interface
f_conf = fckit_configuration(c_conf)
f_attr = fckit_configuration(c_attr)

! Call Fortran
call bifourier_arome_legacy_write_balance(f_conf,f_attr,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu,fact1)

! Release memory
call f_conf%final()
call f_attr%final()

end subroutine bifourier_arome_legacy_write_balance_c

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_read_covariance_c(c_conf,c_attr,nsmaxp1,nflev,vorcov,divucov,tpsucov,qucov) &
 & bind(c,name="bifourier_arome_legacy_read_covariance_f90")

implicit none

! Passed variables
type(c_ptr),intent(in),value :: c_conf
type(c_ptr),intent(in),value :: c_attr
integer(c_int),intent(in) :: nsmaxp1
integer(c_int),intent(in) :: nflev
real(c_double),intent(inout) :: vorcov(nsmaxp1*nflev*nflev)
real(c_double),intent(inout) :: divucov(nsmaxp1*nflev*nflev)
real(c_double),intent(inout) :: tpsucov(nsmaxp1*(nflev+1)*(nflev+1))
real(c_double),intent(inout) :: qucov(nsmaxp1*nflev*nflev)

! Local variables
type(fckit_configuration) :: f_conf
type(fckit_configuration) :: f_attr

! Interface
f_conf = fckit_configuration(c_conf)
f_attr = fckit_configuration(c_attr)

! Call Fortran
call bifourier_arome_legacy_read_covariance(f_conf,f_attr,vorcov,divucov,tpsucov,qucov)

! Release memory
call f_conf%final()
call f_attr%final()

end subroutine bifourier_arome_legacy_read_covariance_c

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_write_covariance_c(c_conf,c_attr,nsmaxp1,nflev,vorcov,divucov,tpsucov,qucov) &
 & bind(c,name="bifourier_arome_legacy_write_covariance_f90")

implicit none

! Passed variables
type(c_ptr),intent(in),value :: c_conf
type(c_ptr),intent(in),value :: c_attr
integer(c_int),intent(in) :: nsmaxp1
integer(c_int),intent(in) :: nflev
real(c_double),intent(in) :: vorcov(nsmaxp1*nflev*nflev)
real(c_double),intent(in) :: divucov(nsmaxp1*nflev*nflev)
real(c_double),intent(in) :: tpsucov(nsmaxp1*(nflev+1)*(nflev+1))
real(c_double),intent(in) :: qucov(nsmaxp1*nflev*nflev)

! Local variables
type(fckit_configuration) :: f_conf
type(fckit_configuration) :: f_attr

! Interface
f_conf = fckit_configuration(c_conf)
f_attr = fckit_configuration(c_attr)

! Call Fortran
call bifourier_arome_legacy_write_covariance(f_conf,f_attr,vorcov,divucov,tpsucov,qucov)

! Release memory
call f_conf%final()
call f_attr%final()

end subroutine bifourier_arome_legacy_write_covariance_c

!----------------------------------------------------------------------

end module bifourier_arome_legacy_interface
