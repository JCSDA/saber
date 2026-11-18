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

subroutine bifourier_arome_legacy_read_balance(conf,attr,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu,fact1)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
type(fckit_configuration),intent(in) :: attr
real(kind_real),intent(inout) :: sdivpb(:)
real(kind_real),intent(inout) :: stpspb(:)
real(kind_real),intent(inout) :: stpsdivu(:)
real(kind_real),intent(inout) :: sqpb(:)
real(kind_real),intent(inout) :: sqdivu(:)
real(kind_real),intent(inout) :: sqtpsu(:)
real(kind_real),intent(inout) :: fact1(:)

! Local variables
integer(kind_int),parameter :: iultmp = 10
integer(kind_int) :: itestwd,idate,idim1,idim2,ilendef,inbmat,inbset,iorig,ipar1,ipar2,isetdist,itime,itypdi1,itypdi2, &
 & itypmat,iweight,jj,jk,jn,ndgl,ndlon,ndgux,ndlux,nsmax,nmsmax,nflev,kspec2g,nsmax_file,nflev_file,kspec2g_file
real(kind_real) :: elat0,elat1,elat2,elon0,elon1,elon2
character(len=10) :: clid
character(len=70) :: clcom
character(len=1024) :: cdfile
character(len=:),allocatable :: str

! Get filename from configuration
call conf%get_or_die("input file",str)
cdfile = str

! Get attributes
call attr%get_or_die("nsmax",nsmax)
call attr%get_or_die("nflev",nflev)
call attr%get_or_die("kspec2g",kspec2g)

! Open file
open(iultmp,file=cdfile,form="unformatted",convert="big_endian")

! Read and check clid
read(iultmp) clid
write(*,"(a,a)") "Info     : - GSA ID: ",clid
if (clid /= "ALADIN98") call abor1_ftn("bad id in gsa file")

! Read description
read(iultmp) clcom
write(*,"(a,a)") "Info     : - Description : ",clcom

! Read center and date
read(iultmp) iorig,idate,itime,inbset
write(*,"(a,i3)") "Info     : - Center: ",iorig
write(*,"(a,i8,a,i6)") "Info     : - Date/time: ",idate," / ",itime

! Read gsa set 0: model geometry definition
write(*,"(a)") "Info     : - Reading gsa set 0: model geometry definition"
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (itypmat /= 0) call abor1_ftn("no model geometry description")
if ((idim1 /= 1).or.(idim2 /= 13).or.(ipar1 /= 50).or.(ipar2 /= 0).or.(itypdi1 /= 0).or.(itypdi2 /= 0)) then
  call abor1_ftn("nonexpected parameters for model geometry description")
end if
read(iultmp) elon1,elat1,elon2,elat2,elon0,elat0,ndgl,ndlon,ndgux,ndlux,nsmax_file,nmsmax,nflev_file,itestwd
if (itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
write(*,"(a,i5,a,i3)") "Info     : - File geometry : nsmax =",nsmax_file," / nflev =",nflev_file

! Check sizes
if (nsmax_file /= nsmax_file) call abor1_ftn("inconsistent number of total wavenumbers in balance file")
if (nflev /= nflev_file) call abor1_ftn("inconsistent number of levels in balance file")

! Read gsa set 1: header
write(*,"(a)") "Info     : - Reading gsa set 1: fact1"
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,kspec2g_file,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (itypmat /= 4) call abor1_ftn("no horizontal balance in gsa set 1")

! Check size
if (kspec2g /= kspec2g_file) call abor1_ftn("inconsistent number of wavenumbers in fact1 file")

! Read gsa set 1: fact1
read(iultmp)
read(iultmp) (fact1(jj),jj=1,kspec2g),itestwd
if (itestwd /= ichkwd) call abor1_ftn("bad gsa control word")

! Read gsa set 2: header
write(*,"(a)") "Info     : - Reading gsa set 2: sdivpb"
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (itypmat /= 5) call abor1_ftn("not vert balance in gsa set 2")
if ((idim1 /= nflev).or.(idim2/=nflev)) call abor1_ftn("bad vertical resolution in gsa set 2")
if ((ipar1 /= 11).or.(ipar2 /= 15)) call abor1_ftn("not pb->divb operator in gsa set 2")

! Read gsa set 2: sdivpb
do jn=1,nsmax+1
  read(iultmp)
  read(iultmp) ((sdivpb((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jk=1,nflev),jj=1,nflev),itestwd
  if (itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Read gsa set 3: header
write(*,"(a)") "Info     : - Reading gsa set 3: stpspb"
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1 /= 13).or.(ipar2 /= 15)) call abor1_ftn("no pb->tpsb operator in gsa set 3")

! Read gsa set 3: stpspb
do jn=1,nsmax+1
  read(iultmp)
  read(iultmp) ((stpspb((jn-1)*nflev*(nflev+1)+(jk-1)*(nflev+1)+jj),jk=1,nflev),jj=1,nflev+1),itestwd
  if(itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Read gsa set 4: header
write(*,"(a)") "Info     : - Reading gsa set 4: stpsdivu"
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1 /= 13).or.(ipar2 /= 12)) call abor1_ftn("not divu->tpsb operator in gsa set 4")

! Read gsa set 4: stpsdivu
do jn=1,nsmax+1
  read(iultmp)
  read(iultmp) ((stpsdivu((jn-1)*nflev*(nflev+1)+(jk-1)*(nflev+1)+jj),jk=1,nflev),jj=1,nflev+1),itestwd
  if(itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Read gsa set 5: header
write(*,"(a)") "Info     : - Reading gsa set 5: sqpb"
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1 /= 16).or.(ipar2 /= 15)) call abor1_ftn("no pb->qb operator in gsa set 5")

! Read gsa set 5: sqpb
do jn=1,nsmax+1
  read(iultmp)
  read(iultmp) ((sqpb((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jk=1,nflev),jj=1,nflev),itestwd
  if (itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Read gsa set 6: header
write(*,"(a)") "Info     : - Reading gsa set 6: sqdivu"
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1 /= 16).or.(ipar2 /= 12)) call abor1_ftn("no divu->qb operator in gsa set 6")

! Read gsa set 6: sqdivu
do jn=1,nsmax+1
  read(iultmp)
  read(iultmp) ((sqdivu((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jk=1,nflev),jj=1,nflev),itestwd
  if (itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Read gsa set 7: header
write(*,"(a)") "Info     : - Reading gsa set 7: sqtpsu"
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((ipar1 /= 16).or.(ipar2 /= 14)) call abor1_ftn("no tpsu->qb operator in gsa set 7")

! Read gsa set 7: sqtpsu
do jn=1,nsmax+1
  read(iultmp)
  read(iultmp) ((sqtpsu((jn-1)*(nflev+1)*nflev+(jk-1)*nflev+jj),jk=1,nflev+1),jj=1,nflev),itestwd
  if(itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_read_balance

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_write_balance(conf,attr,sdivpb,stpspb,stpsdivu,sqpb,sqdivu,sqtpsu,fact1)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
type(fckit_configuration),intent(in) :: attr
real(kind_real),intent(in) :: sdivpb(:)
real(kind_real),intent(in) :: stpspb(:)
real(kind_real),intent(in) :: stpsdivu(:)
real(kind_real),intent(in) :: sqpb(:)
real(kind_real),intent(in) :: sqdivu(:)
real(kind_real),intent(in) :: sqtpsu(:)
real(kind_real),intent(in) :: fact1(:)

! Local variables
integer(kind_int),parameter :: iultmp = 10
integer(kind_int) :: iorig,idate,itime,iweight,jj,jk,jn,ndgl,ndlon,ndgux,ndlux,nsmax,nmsmax,nflev,kspec2g
real(kind_real) :: elat0,elat1,elat2,elon0,elon1,elon2
real(kind_real),allocatable :: zpres(:)
character(len=10) :: clid
character(len=70) :: clcom
character(len=1024) :: cdfile
character(len=:),allocatable :: str

! Get filename from configuration
call conf%get_or_die("output file",str)
cdfile = str

! Get attributes
call attr%get_or_die("clid",str)
clid = str
call attr%get_or_die("clcom",str)
clcom = str
call attr%get_or_die("iorig",iorig)
idate = 0
itime = 0
iweight = 0
call attr%get_or_die("elon0",elon0)
call attr%get_or_die("elat0",elat0)
call attr%get_or_die("elon1",elon1)
call attr%get_or_die("elat1",elat1)
call attr%get_or_die("elon2",elon2)
call attr%get_or_die("elat2",elat2)
call attr%get_or_die("ndgl",ndgl)
call attr%get_or_die("ndlon",ndlon)
call attr%get_or_die("ndgux",ndgux)
call attr%get_or_die("ndlux",ndlux)
call attr%get_or_die("nsmax",nsmax)
call attr%get_or_die("nmsmax",nmsmax)
call attr%get_or_die("nflev",nflev)
call attr%get_or_die("kspec2g",kspec2g)

! Allocation
allocate(zpres(nflev+1))

! Prepare zpres
do jj=1,nflev+1
  zpres(jj) = real(jj,kind=kind_real)
end do

! Open file
open(iultmp,file=cdfile,form="unformatted",convert="big_endian")

! Write clid
write(iultmp) clid

! Write description
write(iultmp) clcom

! Write center and date
write(iultmp) iorig,idate,itime,8

! Write gsa set 0: model geometry definition
write(*,"(a)") "Info     : - Writing gsa set 0: model geometry definition"
write(iultmp) 1,iweight,0,1,0
write(iultmp) 1,13,50,0,0,0
write(iultmp)
write(iultmp)
write(iultmp) elon1,elat1,elon2,elat2,elon0,elat0,ndgl,ndlon,ndgux,ndlux,nsmax,nmsmax,nflev,ichkwd

! Write gsa set 1: header
write(*,"(a)") "Info     : - Writing gsa set 1: fact1"
write(iultmp) 1,iweight,4,4,1
write(iultmp) 1,kspec2g,15,4,0,3
write(iultmp)
write(iultmp)

! Write gsa set 1: fact1
write(iultmp) real(1,kind=kind_real)
write(iultmp) (fact1(jj),jj=1,kspec2g),ichkwd
if (ichkwd/=ichkwd) call abor1_ftn("bad gsa control word")

! Write gsa set 2: header
write(*,"(a)") "Info     : - Writing gsa set 2: sdivpb"
write(iultmp) nsmax+1,iweight,5,2,1
write(iultmp) nflev,nflev,11,15,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 2: sdivpb
do jn=1,nsmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((sdivpb((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jk=1,nflev),jj=1,nflev),ichkwd
end do

! Write gsa set 3: header
write(*,"(a)") "Info     : - Writing gsa set 3: stpspb"
write(iultmp) nsmax+1,iweight,5,2,1
write(iultmp) nflev+1,nflev,13,15,1,1
write(iultmp) (zpres(jj),jj=1,nflev+1)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 3: stpspb
do jn=1,nsmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((stpspb((jn-1)*nflev*(nflev+1)+(jk-1)*(nflev+1)+jj),jk=1,nflev),jj=1,nflev+1),ichkwd
end do

! Write gsa set 4: header
write(*,"(a)") "Info     : - Writing gsa set 4: stpsdivu"
write(iultmp) nsmax+1,iweight,5,2,1
write(iultmp) nflev+1,nflev,13,12,1,1
write(iultmp) (zpres(jj),jj=1,nflev+1)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 4: stpsdivu
do jn=1,nsmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((stpsdivu((jn-1)*nflev*(nflev+1)+(jk-1)*(nflev+1)+jj),jk=1,nflev),jj=1,nflev+1),ichkwd
end do

! Write gsa set 5: header
write(*,"(a)") "Info     : - Writing gsa set 5: sqpb"
write(iultmp) nsmax+1,iweight,5,2,1
write(iultmp) nflev,nflev,16,15,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 5: sqpb
do jn=1,nsmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((sqpb((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jk=1,nflev),jj=1,nflev),ichkwd
end do

! Write gsa set 6: header
write(*,"(a)") "Info     : - Writing gsa set 6: sqdivu"
write(iultmp) nsmax+1,iweight,5,2,1
write(iultmp) nflev,nflev,16,12,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 6: sqdivu
do jn=1,nsmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((sqdivu((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jk=1,nflev),jj=1,nflev),ichkwd
end do

! Write gsa set 7: header
write(*,"(a)") "Info     : - Writing gsa set 7: sqtpsu"
write(iultmp) nsmax+1,iweight,5,2,1
write(iultmp) nflev,nflev+1,16,14,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev+1)

! Write gsa set 7: sqtpsu
do jn=1,nsmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp)((sqtpsu((jn-1)*(nflev+1)*nflev+(jk-1)*nflev+jj),jk=1,nflev+1),jj=1,nflev),ichkwd
end do

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_write_balance

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_read_covariance(conf,attr,vorcov,divucov,tpsucov,qucov)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
type(fckit_configuration),intent(in) :: attr
real(kind_real),intent(inout) :: vorcov(:)
real(kind_real),intent(inout) :: divucov(:)
real(kind_real),intent(inout) :: tpsucov(:)
real(kind_real),intent(inout) :: qucov(:)

! Local variables
integer(kind_int),parameter :: iultmp = 10
integer(kind_int) :: itestwd,idate,idim1,idim2,ilendef,inbmat,inbset,iorig,ipar1,ipar2,isetdist,itime,itypdi1,itypdi2, &
 & itypmat,iweight,jj,jk,jn,ndgl,ndlon,ndgux,ndlux,nsmax,nmsmax,nflev,nsmax_file,nflev_file
real(kind_real) :: elat0,elat1,elat2,elon0,elon1,elon2,zdummy
real(kind_real),allocatable :: zpdat(:)
character(len=10) :: clid
character(len=70) :: clcom
character(len=1024) :: cdfile
character(len=:),allocatable :: str

! Get filename from configuration
call conf%get_or_die("input file",str)
cdfile = str

! Get attributes
call attr%get_or_die("nsmax",nsmax)
call attr%get_or_die("nflev",nflev)

! Allocation
allocate(zpdat(nflev+1))

! Open file
open(iultmp,file=cdfile,form="unformatted",convert="big_endian")

! Read and check clid
read(iultmp) clid
write(*,"(a,a)") "Info     : - GSA ID: ",clid
if (clid /= "ALADIN98") call abor1_ftn("bad id in gsa file")

! Read description
read(iultmp) clcom
write(*,"(a,a)") "Info     : - Description : ",clcom

! Read center and date
read(iultmp) iorig,idate,itime,inbset
write(*,"(a,i3)") "Info     : - Center: ",iorig
write(*,"(a,i8,a,i6)") "Info     : - Date/time: ",idate," / ",itime

! Read gsa set 0: model geometry definition
write(*,"(a)") "Info     : - Reading gsa set 0: model geometry definition"
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (itypmat /= 0) call abor1_ftn("no model geometry description")
if ((idim1 /= 1).or.(idim2 /= 13).or.(ipar1 /= 50).or.(ipar2 /= 0).or.(itypdi1 /= 0).or.(itypdi2 /= 0)) then
  call abor1_ftn("nonexpected parameters for model geometry description")
end if
read(iultmp) elon1,elat1,elon2,elat2,elon0,elat0,ndgl,ndlon,ndgux,ndlux,nsmax_file,nmsmax,nflev_file,itestwd
if (itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
write(*,"(a,i5,a,i3)") "Info     : - File geometry : nsmax =",nsmax_file," / nflev =",nflev_file

! Check sizes
if (nsmax_file /= nsmax_file) call abor1_ftn("inconsistent number of total wavenumbers in covariance file")
if (nflev /= nflev_file) call abor1_ftn("inconsistent number of levels in covariance file")

! Read gsa set 1: header
write(*,"(a)") "Info     : - Reading gsa set 1: vorCov"
read(iultmp) inbmat,iweight,itypmat,isetdist,ilendef
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if ((idim1 /= idim2).or.(ipar1 /= ipar2).or.(itypdi1 /= itypdi2)) call abor1_ftn("nonsymmetric matrix")
if (idim1 <= 0) call abor1_ftn("bad matrix dimensions")
if (idim1 /= nflev) call abor1_ftn("code/data dim mismatch")
if (itypdi1 /= 1) call abor1_ftn("matrix not on pressure levels")
if (ipar1 /= 4) call abor1_ftn("not vorticity in gsa set 1")

! Read gsa set 1: vorCov
do jn=1,nsmax+1
  read(iultmp) zdummy
  read(iultmp) ((vorcov((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jj=1,nflev),jk=1,nflev),itestwd
  if (itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Read gsa set 2: header
write(*,"(a)") "Info     : - Reading gsa set 2: divuCov"
read(iultmp)
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (ipar1 /= 12) call abor1_ftn("not unbal div in gsa set 2")

! Read gsa set 2: divuCov
do jn=1,nsmax+1
  read(iultmp) zdummy
  read(iultmp) ((divuCov((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jj=1,nflev),jk=1,nflev),itestwd
  if (itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Read gsa set 3: header
write(*,"(a)") "Info     : - Reading gsa set 3: tPsuCov"
read(iultmp)
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp) (zpdat(jj),jj=1,idim1)
read(iultmp)
if (ipar1 /= 14) call abor1_ftn("not unbal t,lnps in gsa set 3")
if (idim1 /= nflev+1) call abor1_ftn("code/data dim mismatch")

! Read gsa set 3: tPsuCov
do jn=1,nsmax+1
  read(iultmp) zdummy
  read(iultmp) ((tPsuCov((jn-1)*(nflev+1)*(nflev+1)+(jk-1)*(nflev+1)+jj),jj=1,nflev+1),jk=1,nflev+1),itestwd
  if(itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Read gsa set 4: header
write(*,"(a)") "Info     : - Reading gsa set 4: quCov"
read(iultmp)
read(iultmp) idim1,idim2,ipar1,ipar2,itypdi1,itypdi2
read(iultmp)
read(iultmp)
if (ipar1 /= 17) call abor1_ftn("not unbal q in gsa set 4")

! Read gsa set 4: quCov
do jn=1,nsmax+1
  read(iultmp) zdummy
  read(iultmp) ((quCov((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jj=1,nflev),jk=1,nflev),itestwd
  if(itestwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_read_covariance

!----------------------------------------------------------------------

subroutine bifourier_arome_legacy_write_covariance(conf,attr,vorcov,divucov,tpsucov,qucov)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
type(fckit_configuration),intent(in) :: attr
real(kind_real),intent(in) :: vorcov(:)
real(kind_real),intent(in) :: divucov(:)
real(kind_real),intent(in) :: tpsucov(:)
real(kind_real),intent(in) :: qucov(:)

! Local variables
integer(kind_int),parameter :: iultmp = 10
integer(kind_int) :: iorig,idate,itime,iweight,jj,jk,jn,ndgl,ndlon,ndgux,ndlux,nsmax,nmsmax,nflev
real(kind_real) :: elat0,elat1,elat2,elon0,elon1,elon2
real(kind_real),allocatable :: zpres(:)
character(len=10) :: clid
character(len=70) :: clcom
character(len=1024) :: cdfile
character(len=:),allocatable :: str

! Get filename from configuration
call conf%get_or_die("output file",str)
cdfile = str

! Get attributes
call attr%get_or_die("clid",str)
clid = str
call attr%get_or_die("clcom",str)
clcom = str
call attr%get_or_die("iorig",iorig)
idate = 0
itime = 0
iweight = 0
call attr%get_or_die("elon0",elon0)
call attr%get_or_die("elat0",elat0)
call attr%get_or_die("elon1",elon1)
call attr%get_or_die("elat1",elat1)
call attr%get_or_die("elon2",elon2)
call attr%get_or_die("elat2",elat2)
call attr%get_or_die("ndgl",ndgl)
call attr%get_or_die("ndlon",ndlon)
call attr%get_or_die("ndgux",ndgux)
call attr%get_or_die("ndlux",ndlux)
call attr%get_or_die("nsmax",nsmax)
call attr%get_or_die("nmsmax",nmsmax)
call attr%get_or_die("nflev",nflev)

! Allocation
allocate(zpres(nflev+1))

! Prepare zpres
do jj=1,nflev+1
  zpres(jj) = real(jj,kind=kind_real)
end do

! Open file
open(iultmp,file=cdfile,form="unformatted",convert="big_endian")

! Write clid
write(iultmp) clid

! Write description
write(iultmp) clcom

! Write center and date
write(iultmp) iorig,idate,itime,8

! Write gsa set 0: model geometry definition
write(*,"(a)") "Info     : - Writing gsa set 0: model geometry definition"
write(iultmp) 1,iweight,0,1,0
write(iultmp) 1,13,50,0,0,0
write(iultmp)
write(iultmp)
write(iultmp) elon1,elat1,elon2,elat2,elon0,elat0,ndgl,ndlon,ndgux,ndlux,nsmax,nmsmax,nflev,ichkwd

! Write gsa set 1: header
write(*,"(a)") "Info     : - Writing gsa set 1: vorCov"
write(iultmp) nsmax+1,45,1,2,1
write(iultmp) nflev,nflev,4,4,1,1
write(iultmp)
write(iultmp)

! Write gsa set 1: vorCov
do jn=1,nsmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((vorcov((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jj=1,nflev),jk=1,nflev),ichkwd
end do

! Write gsa set 2: header
write(*,"(a)") "Info     : - Writing gsa set 2: divuCov"
write(iultmp) nsmax+1,45,1,2,1
write(iultmp) nflev,nflev,12,12,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 2: divuCov
do jn=1,nsmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((divuCov((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jj=1,nflev),jk=1,nflev),ichkwd
end do

! Write gsa set 3: header
write(*,"(a)") "Info     : - Writing gsa set 3: tPsuCov"
write(iultmp) nsmax+1,45,1,2,1
write(iultmp) nflev+1,nflev+1,14,14,1,1
write(iultmp) (zpres(jj),jj=1,nflev+1)
write(iultmp) (zpres(jj),jj=1,nflev+1)

! Write gsa set 3: tPsuCov
do jn=1,nsmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((tPsuCov((jn-1)*(nflev+1)*(nflev+1)+(jk-1)*(nflev+1)+jj),jj=1,nflev+1),jk=1,nflev+1),ichkwd
end do

! Write gsa set 4: header
write(*,"(a)") "Info     : - Writing gsa set 4: quCov"
write(iultmp) nsmax+1,45,1,2,1
write(iultmp) nflev,nflev,17,17,1,1
write(iultmp) (zpres(jj),jj=1,nflev)
write(iultmp) (zpres(jj),jj=1,nflev)

! Write gsa set 4: quCov
do jn=1,nsmax+1
  write(iultmp) real(jn-1,kind=kind_real)
  write(iultmp) ((quCov((jn-1)*nflev*nflev+(jk-1)*nflev+jj),jj=1,nflev),jk=1,nflev),ichkwd
  if(ichkwd /= ichkwd) call abor1_ftn("bad gsa control word")
end do

! Close file
close(iultmp)

end subroutine bifourier_arome_legacy_write_covariance

!----------------------------------------------------------------------

end module bifourier_arome_legacy_mod
