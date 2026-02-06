!----------------------------------------------------------------------
! Module: bifourier_arome_legacy_mod
!> Bifourier AROME legacy Fortran module
! Author: Benjamin Menetrier
! Source: readjbbal.F90, ewgsabal.F90, readjbdat96.F90 and ewgsacov.F90
! Copyright 2025 Meteorologisk Institutt
!----------------------------------------------------------------------
module bifourier_arome_legacy_mod

use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_int, kind_real

implicit none

integer(kind_int),parameter :: ichkwd = 3141592

private
public :: bifourier_arome_legacy_read_balance, bifourier_arome_legacy_write_balance, &
 & bifourier_arome_legacy_read_covariance, bifourier_arome_legacy_write_covariance

contains

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_read_balance(conf,nwglb,nflev,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu,nial,fact1)

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
integer(kind_int),intent(in) :: nial
real(kind_real),intent(inout) :: fact1(nial)

! Local variables
integer(kind_int),parameter :: iultmp = 10
integer(kind_int) :: itestwd,idate,idim1,idim2,ilendef,inbmat,inbset,iorig,ipar1,ipar2,isetdist,itime,itypdi1,itypdi2, &
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
if (clid /= 'ALADIN98') call abor1_ftn('bad id in gsa file')

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
if (itypmat /= 0) call abor1_ftn('no model geometry description')
if ((idim1 /= 1).or.(idim2 /= 13).or.(ipar1 /= 50).or.(ipar2 /= 0).or.(itypdi1 /= 0).or.(itypdi2 /= 0)) &
 & call abor1_ftn('nonexpected parameters for model geometry description')
read(iultmp) zlon1,zlat1,zlon2,zlat2,zlon0,zlat0,idgl,idlon,idgux,idlux,ksmax,kmsmax,kflevg,itestwd
if (itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
write(*,'(a,i5,a,i3)') 'Info     : - File geometry : nsmax =',ksmax,' / nflev =',kflevg

! Read gsa set 1: header
write(*,'(a)') 'Info     : - Reading gsa set 1: fact1'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (itypmat /= 4) call abor1_ftn('no horizontal balance in gsa set 1')

! Check size
if (idim2 /= nial) call abor1_ftn('inconsistent number of wavenumbers in fact1 file')

! Read gsa set 1: fact1
read(iultmp)
read(iultmp) (fact1(jj),jj=1,idim2),itestwd
if (itestwd /= ichkwd) call abor1_ftn('bad gsa control word')

! Read gsa set 2: header
write(*,'(a)') 'Info     : - Reading gsa set 2: sdivpb'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (itypmat /= 5) call abor1_ftn('not vert balance in gsa set 2')
if ((idim1 /= kflevg).or.(idim2/=kflevg)) call abor1_ftn('bad vertical resolution in gsa set 2')
if ((ipar1 /= 11).or.(ipar2 /= 15)) call abor1_ftn('not pb->divb operator in gsa set 2')

! Check sizes
if (kflevg /= nflev) call abor1_ftn('inconsistent number of levels in balance file')
if (ksmax+1 /= nwglb) call abor1_ftn('inconsistent number of total wavenumbers in balance file')

! Read gsa set 2: sdivpb
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((sdivpb((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg),jj=1,kflevg),itestwd
  if (itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 3: header
write(*,'(a)') 'Info     : - Reading gsa set 3: stpspb'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1 /= 13).or.(ipar2 /= 15)) call abor1_ftn('no pb->tpsb operator in gsa set 3')

! Read gsa set 3: stpspb
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((stpspb((jn-1)*kflevg*(kflevg+1)+(jk-1)*(kflevg+1)+jj),jk=1,kflevg),jj=1,kflevg+1),itestwd
  if(itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 4: header
write(*,'(a)') 'Info     : - Reading gsa set 4: stpsdivu'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1 /= 13).or.(ipar2 /= 12)) call abor1_ftn('not divu->tpsb operator in gsa set 4')

! Read gsa set 4: stpsdivu
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((stpsdivu((jn-1)*kflevg*(kflevg+1)+(jk-1)*(kflevg+1)+jj),jk=1,kflevg),jj=1,kflevg+1),itestwd
  if(itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 5: header
write(*,'(a)') 'Info     : - Reading gsa set 5: sqpb'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1 /= 16).or.(ipar2 /= 15)) call abor1_ftn('no pb->qb operator in gsa set 5')

! Read gsa set 5: sqpb
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((sqpb((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg),jj=1,kflevg),itestwd
  if (itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 6: header
write(*,'(a)') 'Info     : - Reading gsa set 6: sqdivu'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1 /= 16).or.(ipar2 /= 12)) call abor1_ftn('no divu->qb operator in gsa set 6')

! Read gsa set 6: sqdivu
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((sqdivu((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg),jj=1,kflevg),itestwd
  if (itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 7: header
write(*,'(a)') 'Info     : - Reading gsa set 7: sqtpsu'
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1 /= 16).or.(ipar2 /= 14)) call abor1_ftn('no tpsu->qb operator in gsa set 7')

! Read gsa set 7: sqtpsu
do jn=1,ksmax+1
  read(iultmp)
  read(iultmp) ((sqtpsu((jn-1)*(kflevg+1)*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg+1),jj=1,kflevg),itestwd
  if(itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_read_balance

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_write_balance(conf,nwglb,nflev,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu,nial,fact1)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
integer(kind_int),intent(in) :: nwglb
integer(kind_int),intent(in) :: nflev
real(kind_real),intent(in) :: sdivpb(nwglb*nflev*nflev)
real(kind_real),intent(in) :: stpspb(nwglb*nflev*(nflev+1))
real(kind_real),intent(in) :: stpsdivu(nwglb*nflev*(nflev+1))
real(kind_real),intent(in) :: sqpb(nwglb*nflev*nflev)
real(kind_real),intent(in) :: sqdivu(nwglb*nflev*nflev)
real(kind_real),intent(in) :: sqtpsu(nwglb*(nflev+1)*nflev)
integer(kind_int),intent(in) :: nial
real(kind_real),intent(in) :: fact1(nial)

! Local variables
integer(kind_int),parameter :: iultmp = 10
integer(kind_int) :: idate,itime,iweight,jj,jk,jn,idgl,idgux,idlon,idlux,ksmax,kmsmax,kflevg,kspec2g
real(kind_real) :: zlat0,zlat1,zlat2,zlon0,zlon1,zlon2
real(kind_real) :: zpres(nflev+1)
character(len=10) :: clid
character(len=70) :: clcom
character(len=1024) :: cdfile
character(len=:),allocatable :: str

! Prepare zpres
do jj=1,nflev+1
  zpres(jj) = real(jj,kind=kind_real)
end do

! Get filename from configuration
call conf%get_or_die("output file",str)
cdfile = str

! Set relevant parameters
ksmax = nwglb-1
kflevg = nflev
kspec2g = nial

! Set dummy parameters
idate = 0
itime = 0
iweight = 0
zlon1 = 0.0
zlon2 = 0.0
zlat1 = 0.0
zlat2 = 0.0
zlon0 = 0.0
zlat0 = 0.0
idgl = 0
idlon = 0
idgux = 0
idlux = 0
kmsmax = 0

! Open file
open(iultmp,file=cdfile,form='unformatted',convert='big_endian')

! Write clid
clid = 'ALADIN98'
write(iultmp) clid

! Write description
clcom = ' Balanced statistcs for a LAM, after L. Berre 1998'
write(iultmp) clcom


! Write center and date
write(iultmp) 85,idate,itime,8

! Write gsa set 0: model geometry definition
write(*,'(a)') 'Info     : - Writing gsa set 0: model geometry definition'
write(iultmp) 1,iweight,0,1,0
write(iultmp) 1,13,50,0,0,0
write(iultmp)
write(iultmp)
write(iultmp) zlon1,zlat1,zlon2,zlat2,zlon0,zlat0,idgl,idlon,idgux,idlux,ksmax,kmsmax,kflevg,ichkwd

! Write gsa set 1: header
write(*,'(a)') 'Info     : - Writing gsa set 1: fact1'
write(iultmp) 1,iweight,4,4,1
write(iultmp) 1,kspec2g,15,4,0,3
write(iultmp)
write(iultmp)

! Write gsa set 1: fact1
write(iultmp) real(1,kind=kind_real)
write(iultmp) (fact1(jj),jj=1,kspec2g),ichkwd
if (ichkwd/=ichkwd) call abor1_ftn('bad gsa control word')

! Write gsa set 2: header
write(*,'(a)') 'Info     : - Writing gsa set 2: sdivpb'
write(iultmp) ksmax+1,iweight,5,2,1
write(iultmp) nflev,nflev,11,15,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 2: sdivpb
do jn=1,ksmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((sdivpb((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg),jj=1,kflevg),ichkwd
end do

! Write gsa set 3: header
write(*,'(a)') 'Info     : - Writing gsa set 3: stpspb'
write(iultmp) ksmax+1,iweight,5,2,1
write(iultmp) nflev+1,nflev,13,15,1,1
write(iultmp) (zpres(jj),jj=1,nflev+1)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 3: stpspb
do jn=1,ksmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((stpspb((jn-1)*kflevg*(kflevg+1)+(jk-1)*(kflevg+1)+jj),jk=1,kflevg),jj=1,kflevg+1),ichkwd
end do

! Write gsa set 4: header
write(*,'(a)') 'Info     : - Writing gsa set 4: stpsdivu'
write(iultmp) ksmax+1,iweight,5,2,1
write(iultmp) nflev+1,nflev,13,12,1,1
write(iultmp) (zpres(jj),jj=1,nflev+1)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 4: stpsdivu
do jn=1,ksmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((stpsdivu((jn-1)*kflevg*(kflevg+1)+(jk-1)*(kflevg+1)+jj),jk=1,kflevg),jj=1,kflevg+1),ichkwd
end do

! Write gsa set 5: header
write(*,'(a)') 'Info     : - Writing gsa set 5: sqpb'
write(iultmp) ksmax+1,iweight,5,2,1
write(iultmp) nflev,nflev,16,15,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 5: sqpb
do jn=1,ksmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((sqpb((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg),jj=1,kflevg),ichkwd
end do

! Write gsa set 6: header
write(*,'(a)') 'Info     : - Writing gsa set 6: sqdivu'
write(iultmp) ksmax+1,iweight,5,2,1
write(iultmp) nflev,nflev,16,12,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 6: sqdivu
do jn=1,ksmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((sqdivu((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg),jj=1,kflevg),ichkwd
end do

! Write gsa set 7: header
write(*,'(a)') 'Info     : - Writing gsa set 7: sqtpsu'
write(iultmp) ksmax+1,iweight,5,2,1
write(iultmp) nflev,nflev+1,16,14,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev+1)

! Write gsa set 7: sqtpsu
do jn=1,ksmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp)((sqtpsu((jn-1)*(kflevg+1)*kflevg+(jk-1)*kflevg+jj),jk=1,kflevg+1),jj=1,kflevg),ichkwd
end do

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_write_balance

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_read_covariance(conf,nwglb,nflev,vorcov,divucov,tpsucov,qucov)

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
integer(kind_int) :: itestwd,idate,idim1,idim2,ilendef,inbmat,inbset,iorig,ipar1,ipar2,isetdist,itime,itypdi1,itypdi2,&
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
if (clid /= 'ALADIN98') call abor1_ftn('bad id in gsa file')

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
if (itypmat /= 0) call abor1_ftn('no model geometry description')
if ((idim1 /= 1).or.(idim2 /= 13).or.(ipar1 /= 50).or.(ipar2 /= 0).or.(itypdi1 /= 0).or.(itypdi2 /= 0)) &
 & call abor1_ftn('nonexpected parameters for model geometry description')
read(iultmp) zlon1,zlat1,zlon2,zlat2,zlon0,zlat0,idgl,idlon,idgux,idlux,ksmax,kmsmax,kflevg,itestwd
if (itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
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
if ((idim1 /= idim2).or.(ipar1 /= ipar2).or.(itypdi1 /= itypdi2)) call abor1_ftn('nonsymmetric matrix')
if (idim1 <= 0) call abor1_ftn('bad matrix dimensions')
if (idim1 /= kflevg) call abor1_ftn('code/data dim mismatch')
if (itypdi1 /= 1) call abor1_ftn('matrix not on pressure levels')
if (ipar1 /= 4) call abor1_ftn('not vorticity in gsa set 1')

! Read gsa set 1: vorCov
do jn=1,ksmax+1
  read(iultmp) zdummy
  read(iultmp) ((vorcov((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jj=1,kflevg),jk=1,kflevg),itestwd
  if (itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 2: header
write(*,'(a)') 'Info     : - Reading gsa set 2: divuCov'
read(iultmp)
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (ipar1 /= 12) call abor1_ftn('not unbal div in gsa set 2')

! Read gsa set 2: divuCov
do jn=1,ksmax+1
  read(iultmp) zdummy
  read(iultmp) ((divuCov((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jj=1,kflevg),jk=1,kflevg),itestwd
  if (itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 3: header
write(*,'(a)') 'Info     : - Reading gsa set 3: tPsuCov'
read(iultmp)
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp) (zpdat(jj),jj=1,idim1)
read(iultmp)
if (ipar1 /= 14) call abor1_ftn('not unbal t,lnps in gsa set 3')
if (idim1 /= kflevg+1) call abor1_ftn('code/data dim mismatch')

! Read gsa set 3: tPsuCov
do jn=1,ksmax+1
  read(iultmp) zdummy
  read(iultmp) ((tPsuCov((jn-1)*(kflevg+1)*(kflevg+1)+(jk-1)*(kflevg+1)+jj),jj=1,kflevg+1),jk=1,kflevg+1),itestwd
  if(itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Read gsa set 4: header
write(*,'(a)') 'Info     : - Reading gsa set 4: quCov'
read(iultmp)
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (ipar1 /= 17) call abor1_ftn('not unbal q in gsa set 4')

! Read gsa set 4: quCov
do jn=1,ksmax+1
  read(iultmp) zdummy
  read(iultmp) ((quCov((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jj=1,kflevg),jk=1,kflevg),itestwd
  if(itestwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_read_covariance

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_write_covariance(conf,nwglb,nflev,vorcov,divucov,tpsucov,qucov)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
integer(kind_int),intent(in) :: nwglb
integer(kind_int),intent(in) :: nflev
real(kind_real),intent(in) :: vorcov(nwglb*nflev*nflev)
real(kind_real),intent(in) :: divucov(nwglb*nflev*nflev)
real(kind_real),intent(in) :: tpsucov(nwglb*(nflev+1)*(nflev+1))
real(kind_real),intent(in) :: qucov(nwglb*nflev*nflev)

! Local variables
integer(kind_int),parameter :: iultmp = 10
integer(kind_int) :: idate,itime,iweight,jj,jk,jn,idgl,idgux,idlon,idlux,ksmax,kmsmax,kflevg
real(kind_real) :: zlat0,zlat1,zlat2,zlon0,zlon1,zlon2
real(kind_real) :: zpres(nflev+1)
character(len=10) :: clid
character(len=70) :: clcom
character(len=1024) :: cdfile
character(len=:),allocatable :: str

! Prepare zpres
do jj=1,nflev+1
  zpres(jj) = real(jj,kind=kind_real)
end do

! Get filename from configuration
call conf%get_or_die("output file",str)
cdfile = str

! Set relevant parameters
ksmax = nwglb-1
kflevg = nflev

! Set dummy parameters
idate = 0
itime = 0
iweight = 0
zlon1 = 0.0
zlon2 = 0.0
zlat1 = 0.0
zlat2 = 0.0
zlon0 = 0.0
zlat0 = 0.0
idgl = 0
idlon = 0
idgux = 0
idlux = 0
kmsmax = 0

! Open file
open(iultmp,file=cdfile,form='unformatted',convert='big_endian')

! Write clid
clid = 'ALADIN98'
write(iultmp) clid

! Write description
clcom = ' Balanced statistcs for a LAM, after L. Berre 1998'
write(iultmp) clcom

! Write center and date
write(iultmp) 85,idate,itime,8

! Write gsa set 0: model geometry definition
write(*,'(a)') 'Info     : - Writing gsa set 0: model geometry definition'
write(iultmp) 1,iweight,0,1,0
write(iultmp) 1,13,50,0,0,0
write(iultmp)
write(iultmp)
write(iultmp) zlon1,zlat1,zlon2,zlat2,zlon0,zlat0,idgl,idlon,idgux,idlux,ksmax,kmsmax,kflevg,ichkwd

! Write gsa set 1: header
write(*,'(a)') 'Info     : - Writing gsa set 1: vorCov'
write(iultmp) ksmax+1,45,1,2,1
write(iultmp) nflev,nflev,4,4,1,1
write(iultmp)
write(iultmp)

! Write gsa set 1: vorCov
do jn=1,ksmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((vorcov((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jj=1,kflevg),jk=1,kflevg),ichkwd
end do

! Write gsa set 2: header
write(*,'(a)') 'Info     : - Writing gsa set 2: divuCov'
write(iultmp) ksmax+1,45,1,2,1
write(iultmp) nflev,nflev,12,12,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 2: divuCov
do jn=1,ksmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((divuCov((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jj=1,kflevg),jk=1,kflevg),ichkwd
end do

! Write gsa set 3: header
write(*,'(a)') 'Info     : - Writing gsa set 3: tPsuCov'
write(iultmp) ksmax+1,45,1,2,1
write(iultmp) nflev+1,nflev+1,14,14,1,1
write(iultmp) (zpres(jj),jj=1,nflev+1)
write(iultmp) (zpres(jj),jj=1,nflev+1)

! Write gsa set 3: tPsuCov
do jn=1,ksmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((tPsuCov((jn-1)*(kflevg+1)*(kflevg+1)+(jk-1)*(kflevg+1)+jj),jj=1,kflevg+1),jk=1,kflevg+1),ichkwd
end do

! Write gsa set 4: header
write(*,'(a)') 'Info     : - Writing gsa set 4: quCov'
write(iultmp) ksmax+1,45,1,2,1
write(iultmp) nflev,nflev,17,17,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 4: quCov
do jn=1,ksmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((quCov((jn-1)*kflevg*kflevg+(jk-1)*kflevg+jj),jj=1,kflevg),jk=1,kflevg),ichkwd
  if(ichkwd /= ichkwd) call abor1_ftn('bad gsa control word')
end do

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_write_covariance

!----------------------------------------------------------------------

end module bifourier_arome_legacy_mod
