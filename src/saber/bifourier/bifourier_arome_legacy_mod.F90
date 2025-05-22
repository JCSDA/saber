!----------------------------------------------------------------------
! Module: bifourier_arome_legacy_mod
!> Bifourier AROME legacy Fortran module
! Author: Benjamin Menetrier
! Source: readjbbal.F90 and readjbdat96.F90
! Copyright 2025 Meteorologisk Institutt
!----------------------------------------------------------------------
module bifourier_arome_legacy_mod

use atlas_module, only: atlas_fieldset,atlas_field,atlas_real
use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_int, kind_real

implicit none

private
public :: bifourier_arome_legacy_vortopb,bifourier_arome_legacy_balance,bifourier_arome_legacy_covariance

contains

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_vortopb(conf,nial,fact1)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
integer(kind_int),intent(in) :: nial
real(kind_real),intent(inout) :: fact1(nial)

! Local variables
integer(kind_int),parameter :: iultmp = 10
integer(kind_int) :: ichkwd,idate,idim1,idim2,ilendef,inbmat,inbset,iorig,ipar1,ipar2,isetdist,itime,itypdi1,itypdi2, &
 & itypmat,iweight,jj,jk,jn,idgl,idgux,idlon,idlux,ksmax,kmsmax,kflevg
real(kind_real) :: zlat0,zlat1,zlat2,zlon0,zlon1,zlon2
character(len=10) :: clid
character(len=70) :: clcom
character(len=1024) :: cdfile
character(len=:),allocatable :: str

! Get filename from configuration
call conf%get_or_die("input file",str)
cdfile = str

! Open file
open(iultmp,file=cdfile,form='unformatted',convert='big_endian')

! Read and check clid
read(iultmp) clid
write(*,'(a,a)') 'Info     : - GSA ID: ',clid
if (clid/='ALADIN98') call abor1_ftn('bad id in gsa file')

! Read description
read(iultmp) clcom
write(*,'(a,a)') 'Info     : - Description : ',clcom

! Read center and date
read(iultmp) iorig,idate,itime,inbset
write(*,'(a,i3)') 'Info     : - Center: ',iorig
write(*,'(a,i8,a,i6)') 'Info     : - Date/time: ',idate,' / ',itime

! Read gsa set 0: model geometry definition
write(*,'(a)') 'Info     : - Reading gsa set 0: model geometry definition'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (itypmat/=0) call abor1_ftn('no model geometry description')
if ((idim1/=1).or.(idim2/=13).or.(ipar1/=50).or.(ipar2/=0).or.(itypdi1/=0).or.(itypdi2/=0)) &
 & call abor1_ftn('nonexpected parameters for model geometry description')
read(iultmp) zlon1,zlat1,zlon2,zlat2,zlon0,zlat0,idgl,idlon,idgux,idlux,ksmax,kmsmax,kflevg,ichkwd
if (ichkwd/=3141592) call abor1_ftn('bad gsa control word')
write(*,'(a,f12.8,a,f12.8)') 'Info     : - File geometry : zlat1 =',zlat1,' / zlat2 =',zlat2
write(*,'(a,i5,a,i5)') 'Info     :                   ksmax =',ksmax,' / kmsmax =',kmsmax

! Read gsa set 1: header
write(*,'(a)') 'Info     : - Reading gsa set 1: fact1'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
read(iultmp)
if (itypmat/=4) call abor1_ftn('no horizontal balance in gsa set 1')

! Check size
if (idim2 /= nial) call abor1_ftn('inconsistent number of wavenumbers in fact1 file')

! Read gsa set 1: fact1
read(iultmp) (fact1(jj),jj=1,idim2),ichkwd
if (ichkwd/=3141592) call abor1_ftn('bad gsa control word')

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_vortopb

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_balance(conf,nwglb,nflev,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
integer(kind_int),intent(in) :: nwglb
integer(kind_int),intent(in) :: nflev
real(kind_real),intent(inout) :: sdivpb(nwglb*nflev*nflev)
real(kind_real),intent(inout) :: stpspb(nwglb*nflev*(nflev+1))
real(kind_real),intent(inout) :: stpsdivu(nwglb*nflev*(nflev+1))
real(kind_real),intent(inout) :: sqpb(nwglb*nflev*nflev)
real(kind_real),intent(inout) :: sqdivu(nwglb*nflev*nflev)
real(kind_real),intent(inout) :: sqtpsu(nwglb*(nflev+1)*nflev)

! Local variables
integer(kind_int),parameter :: iultmp = 10
integer(kind_int) :: ichkwd,idate,idim1,idim2,ilendef,inbmat,inbset,iorig,ipar1,ipar2,isetdist,itime,itypdi1,itypdi2, &
 & itypmat,iweight,jj,jk,jn,idgl,idgux,idlon,idlux,ksmax,kmsmax,kflevg
real(kind_real) :: zlat0,zlat1,zlat2,zlon0,zlon1,zlon2
real(kind_real),allocatable :: fact1(:)
character(len=10) :: clid
character(len=70) :: clcom
character(len=1024) :: cdfile
character(len=:),allocatable :: str

! Get filename from configuration
call conf%get_or_die("input file",str)
cdfile = str

! Open file
open(iultmp,file=cdfile,form='unformatted',convert='big_endian')

! Read and check clid
read(iultmp) clid
write(*,'(a,a)') 'Info     : - GSA ID: ',clid
if (clid/='ALADIN98') call abor1_ftn('bad id in gsa file')

! Read description
read(iultmp) clcom
write(*,'(a,a)') 'Info     : - Description : ',clcom

! Read center and date
read(iultmp) iorig,idate,itime,inbset
write(*,'(a,i3)') 'Info     : - Center: ',iorig
write(*,'(a,i8,a,i6)') 'Info     : - Date/time: ',idate,' / ',itime

! Read gsa set 0: model geometry definition
write(*,'(a)') 'Info     : - Reading gsa set 0: model geometry definition'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (itypmat/=0) call abor1_ftn('no model geometry description')
if ((idim1/=1).or.(idim2/=13).or.(ipar1/=50).or.(ipar2/=0).or.(itypdi1/=0).or.(itypdi2/=0)) &
 & call abor1_ftn('nonexpected parameters for model geometry description')
read(iultmp) zlon1,zlat1,zlon2,zlat2,zlon0,zlat0,idgl,idlon,idgux,idlux,ksmax,kmsmax,kflevg,ichkwd
if (ichkwd/=3141592) call abor1_ftn('bad gsa control word')
write(*,'(a,i5,a,i3)') 'Info     : - File geometry : nsmax =',ksmax,' / nflev =',kflevg

! Read gsa set 1: header
write(*,'(a)') 'Info     : - Reading gsa set 1: fact1'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
read(iultmp)
if (itypmat/=4) call abor1_ftn('no horizontal balance in gsa set 1')

! Allocation
allocate(fact1(idim2))

! Read gsa set 1: fact1
read(iultmp) (fact1(jj),jj=1,idim2),ichkwd
if (ichkwd/=3141592) call abor1_ftn('bad gsa control word')

! Release memory
deallocate(fact1)

! Read gsa set 2: header
write(*,'(a)') 'Info     : - Reading gsa set 2: sdivpb'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (itypmat/=5) call abor1_ftn('not vert balance in gsa set 2')
if ((idim1/=kflevg).or.(idim2/=kflevg)) call abor1_ftn('bad vertical resolution in gsa set 2')
if ((ipar1/=11).or.(ipar2/=15)) call abor1_ftn('not pb->divb operator in gsa set 2')

! Check sizes
if (kflevg /= nflev) call abor1_ftn('inconsistent number of levels in balance file')
if (ksmax+1 /= nwglb) call abor1_ftn('inconsistent number of total wavenumbers in balance file')

! Read gsa set 2: sdivpb
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((sdivpb((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg),jj=1,kflevg),ichkwd
  if (ichkwd /= 3141592) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 3: header
write(*,'(a)') 'Info     : - Reading gsa set 3: stpspb'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1/=13).or.(ipar2/=15)) call abor1_ftn('no pb->tpsb operator in gsa set 3')

! Read gsa set 3: stpspb
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((stpspb((jn-1)*kflevg*(kflevg+1)+(jk-1)*(kflevg+1)+jj),jk=1,kflevg),jj=1,kflevg+1),ichkwd
  if(ichkwd/=3141592) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 4: header
write(*,'(a)') 'Info     : - Reading gsa set 4: stpsdivu'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1/=13).or.(ipar2/=12)) call abor1_ftn('not divu->tpsb operator in gsa set 4')

! Read gsa set 4: stpsdivu
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((stpsdivu((jn-1)*kflevg*(kflevg+1)+(jk-1)*(kflevg+1)+jj),jk=1,kflevg),jj=1,kflevg+1),ichkwd
  if(ichkwd /= 3141592) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 5: header
write(*,'(a)') 'Info     : - Reading gsa set 5: sqpb'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1/=16).or.(ipar2/=15)) call abor1_ftn('no pb->qb operator in gsa set 5')

! Read gsa set 5: sqpb
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((sqpb((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg),jj=1,kflevg),ichkwd
  if (ichkwd/=3141592) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 6: header
write(*,'(a)') 'Info     : - Reading gsa set 6: sqdivu'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1/=16).or.(ipar2/=12)) call abor1_ftn('no divu->qb operator in gsa set 6')

! Read gsa set 6: sqdivu
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((sqdivu((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg),jj=1,kflevg),ichkwd
  if (ichkwd/=3141592) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 7: header
write(*,'(a)') 'Info     : - Reading gsa set 7: sqtpsu'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1/=16).or.(ipar2/=14)) call abor1_ftn('no tpsu->qb operator in gsa set 7')

! Read gsa set 7: sqtpsu
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((sqtpsu((jn-1)*(kflevg+1)*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg+1),jj=1,kflevg),ichkwd
  if(ichkwd /= 3141592) call abor1_ftn('bad gsa control word')
end do

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_balance

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_covariance(conf,nwglb,nflev,vorcov,divucov,tpsucov,qucov)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
integer(kind_int),intent(in) :: nwglb
integer(kind_int),intent(in) :: nflev
real(kind_real),intent(inout) :: vorcov(nwglb*nflev*nflev)
real(kind_real),intent(inout) :: divucov(nwglb*nflev*nflev)
real(kind_real),intent(inout) :: tpsucov(nwglb*(nflev+1)*(nflev+1))
real(kind_real),intent(inout) :: qucov(nwglb*nflev*nflev)

! Local variables
integer(kind_int),parameter :: iultmp = 10
integer(kind_int) :: ichkwd,idate,idim1,idim2,ilendef,inbmat,inbset,iorig,ipar1,ipar2,isetdist,itime,itypdi1,itypdi2,&
 & itypmat,iweight,jj,jk,jn,idgl,idgux,idlon,idlux,ksmax,kmsmax,kflevg
real(kind_real) :: zlat0,zlat1,zlat2,zlon0,zlon1,zlon2,zdummy
real(kind_real) :: zpdat(nflev+1)
character(len=10) :: clid
character(len=70) :: clcom
character(len=1024) :: cdfile
character(len=:),allocatable :: str

! Get filename from configuration
call conf%get_or_die("input file",str)
cdfile = str

! Open file
open(iultmp,file=cdfile,form='unformatted',convert='big_endian')

! Read and check clid
read(iultmp) clid
write(*,'(a,a)') 'Info     : - GSA ID: ',clid
if (clid/='ALADIN98') call abor1_ftn('bad id in gsa file')

! Read description
read(iultmp) clcom
write(*,'(a,a)') 'Info     : - Description : ',clcom

! Read center and date
read(iultmp) iorig,idate,itime,inbset
write(*,'(a,i3)') 'Info     : - Center: ',iorig
write(*,'(a,i8,a,i6)') 'Info     : - Date/time: ',idate,' / ',itime

! Read gsa set 0: model geometry definition
write(*,'(a)') 'Info     : - Reading gsa set 0: model geometry definition'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (itypmat/=0) call abor1_ftn('no model geometry description')
if ((idim1/=1).or.(idim2/=13).or.(ipar1/=50).or.(ipar2/=0).or.(itypdi1/=0).or.(itypdi2/=0)) &
 & call abor1_ftn('nonexpected parameters for model geometry description')
read(iultmp) zlon1,zlat1,zlon2,zlat2,zlon0,zlat0,idgl,idlon,idgux,idlux,ksmax,kmsmax,kflevg,ichkwd
if (ichkwd/=3141592) call abor1_ftn('bad gsa control word')
write(*,'(a,i5,a,i3)') 'Info     : - File geometry : nsmax =',ksmax,' / nflev =',kflevg

! Check sizes
if (kflevg /= nflev) call abor1_ftn('inconsistent number of levels in covariance file')
if (ksmax+1 /= nwglb) call abor1_ftn('inconsistent number of total wavenumbers in covariance file')

! Read gsa set 1: header
write(*,'(a)') 'Info     : - Reading gsa set 1: vorCov'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((idim1/=idim2).or.(ipar1/=ipar2).or.(itypdi1/=itypdi2)) call abor1_ftn('nonsymmetric matrix')
if (idim1<=0) call abor1_ftn('bad matrix dimensions')
if (idim1/=kflevg) call abor1_ftn('code/data dim mismatch')
if (itypdi1/=1) call abor1_ftn('matrix not on pressure levels')
if (ipar1/=4) call abor1_ftn('not vorticity in gsa set 1')

! Read gsa set 1: vorCov
do jn=1,ksmax+1
  read(iultmp) zdummy
  read(iultmp) ((vorcov((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jj=1,kflevg),jk=1,kflevg),ichkwd
  if (ichkwd /= 3141592) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 2: header
write(*,'(a)') 'Info     : - Reading gsa set 2: divuCov'
read(iultmp)
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (ipar1/=12) call abor1_ftn('not unbal div in gsa set 2')

! Read gsa set 2: divuCov
do jn=1,ksmax+1
  read(iultmp) zdummy
  read(iultmp) ((divuCov((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jj=1,kflevg),jk=1,kflevg),ichkwd
  if (ichkwd /= 3141592) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 3: header
write(*,'(a)') 'Info     : - Reading gsa set 3: tPsuCov'
read(iultmp)
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp) (zpdat(jj),jj=1,idim1)
read(iultmp)
if (ipar1/=14) call abor1_ftn('not unbal t,lnps in gsa set 3')
if (idim1/=kflevg+1) call abor1_ftn('code/data dim mismatch')

! Read gsa set 3: tPsuCov
do jn=1,ksmax+1
  read(iultmp) zdummy
  read(iultmp) ((tPsuCov((jn-1)*(kflevg+1)*(kflevg+1)+(jk-1)*(kflevg+1)+jj),jj=1,kflevg+1),jk=1,kflevg+1),ichkwd
  if(ichkwd/=3141592) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 4: header
write(*,'(a)') 'Info     : - Reading gsa set 4: quCov'
read(iultmp)
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (ipar1/=17) call abor1_ftn('not unbal q in gsa set 4')

! Read gsa set 4: quCov
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((quCov((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jj=1,kflevg),jk=1,kflevg),ichkwd
  if(ichkwd /= 3141592) call abor1_ftn('bad gsa control word')
end do

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_covariance

!----------------------------------------------------------------------

end module bifourier_arome_legacy_mod
