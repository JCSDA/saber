!----------------------------------------------------------------------
! Module: bifourier_arome_legacy_interface
!> Bifourier AROME legacy Fortran interface
! Author: Benjamin Menetrier
! Copyright 2025 Meteorologisk Institutt
!----------------------------------------------------------------------
module bifourier_arome_legacy_interface

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding, only: c_ptr, c_int, c_double
use bifourier_arome_legacy_mod, only: bifourier_arome_legacy_vortopb, bifourier_arome_legacy_balance, &
 & bifourier_arome_legacy_covariance

implicit none

private

contains

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_vortopb_c(c_conf,nial,fact1) &
 & bind(c,name='bifourier_arome_legacy_vortopb_f90')
implicit none

! Passed variables
type(c_ptr),intent(in),value :: c_conf
integer(c_int),intent(in) :: nial
real(c_double),intent(inout) :: fact1(nial)

! Local variables
type(fckit_configuration) :: f_conf

! Interface
f_conf = fckit_configuration(c_conf)

! Call Fortran
call bifourier_arome_legacy_vortopb(f_conf,nial,fact1)

! Release memory
call f_conf%final()

end subroutine bifourier_arome_legacy_vortopb_c

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_balance_c(c_conf,nwglb,nflev,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu) &
 & bind(c,name='bifourier_arome_legacy_balance_f90')
implicit none

! Passed variables
type(c_ptr),intent(in),value :: c_conf
integer(c_int),intent(in) :: nwglb
integer(c_int),intent(in) :: nflev
real(c_double),intent(inout) :: sdivpb(nwglb*nflev*nflev)
real(c_double),intent(inout) :: stpspb(nwglb*nflev*(nflev+1))
real(c_double),intent(inout) :: stpsdivu(nwglb*nflev*(nflev+1))
real(c_double),intent(inout) :: sqpb(nwglb*nflev*nflev)
real(c_double),intent(inout) :: sqdivu(nwglb*nflev*nflev)
real(c_double),intent(inout) :: sqtpsu(nwglb*(nflev+1)*nflev)

! Local variables
type(fckit_configuration) :: f_conf

! Interface
f_conf = fckit_configuration(c_conf)

! Call Fortran
call bifourier_arome_legacy_balance(f_conf,nwglb,nflev,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu)

! Release memory
call f_conf%final()

end subroutine bifourier_arome_legacy_balance_c

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_covariance_c(c_conf,nwglb,nflev,vorcov,divucov,tpsucov,qucov) &
 & bind(c,name='bifourier_arome_legacy_covariance_f90')

implicit none

! Passed variables
type(c_ptr),intent(in),value :: c_conf
integer(c_int),intent(in) :: nwglb
integer(c_int),intent(in) :: nflev
real(c_double),intent(inout) :: vorcov(nwglb*nflev*nflev)
real(c_double),intent(inout) :: divucov(nwglb*nflev*nflev)
real(c_double),intent(inout) :: tpsucov(nwglb*(nflev+1)*(nflev+1))
real(c_double),intent(inout) :: qucov(nwglb*nflev*nflev)

! Local variables
type(fckit_configuration) :: f_conf

! Interface
f_conf = fckit_configuration(c_conf)

! Call Fortran
call bifourier_arome_legacy_covariance(f_conf,nwglb,nflev,vorcov,divucov,tpsucov,qucov)

! Release memory
call f_conf%final()

end subroutine bifourier_arome_legacy_covariance_c

!----------------------------------------------------------------------

end module bifourier_arome_legacy_interface
