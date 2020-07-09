!----------------------------------------------------------------------
! Module: type_obsop
! Purpose: observation operator data derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_obsop

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_min,fckit_mpi_max,fckit_mpi_status
use netcdf
use tools_const, only: pi,deg2rad,rad2deg,reqkm
use tools_func, only: lonlatmod,sphere_dist
use tools_kinds, only: kind_real,nc_kind_real,huge_real
use tools_repro, only: rth
use type_com, only: com_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

! Observation operator data derived type
type obsop_type
   ! Observations
   integer :: nobs                          ! Number of observations
   real(kind_real),allocatable :: lonobs(:) ! Observations longitudes
   real(kind_real),allocatable :: latobs(:) ! Observations latitudes
   integer,allocatable :: obsa_to_obs(:)    ! Local to global observation

   ! Required data to apply an observation operator

   ! Number of points
   integer :: nc0b                          ! Halo B size

   ! Number of observations
   integer :: nobsa                         ! Local number of observations

   ! Interpolation data
   type(linop_type) :: h                    ! Interpolation data

   ! Communication data
   type(com_type) :: com                    ! Communication data
contains
   procedure :: partial_dealloc => obsop_partial_dealloc
   procedure :: dealloc => obsop_dealloc
   procedure :: read => obsop_read
   procedure :: write => obsop_write
   procedure :: from => obsop_from
   procedure :: run_obsop => obsop_run_obsop
   procedure :: run_obsop_tests => obsop_run_obsop_tests
   procedure :: apply => obsop_apply
   procedure :: apply_ad => obsop_apply_ad
   procedure :: test_adjoint => obsop_test_adjoint
   procedure :: test_accuracy => obsop_test_accuracy
end type obsop_type

private
public :: obsop_type

contains

!----------------------------------------------------------------------
! Subroutine: obsop_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine obsop_partial_dealloc(obsop)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data

! Release memory
if (allocated(obsop%lonobs)) deallocate(obsop%lonobs)
if (allocated(obsop%latobs)) deallocate(obsop%latobs)

end subroutine obsop_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: obsop_dealloc
! Purpose: release memory (full)
!----------------------------------------------------------------------
subroutine obsop_dealloc(obsop)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data

! Release memory
call obsop%partial_dealloc
if (allocated(obsop%obsa_to_obs)) deallocate(obsop%obsa_to_obs)
call obsop%h%dealloc
call obsop%com%dealloc

end subroutine obsop_dealloc

!----------------------------------------------------------------------
! Subroutine: obsop_read
! Purpose: read observations locations
!----------------------------------------------------------------------
subroutine obsop_read(obsop,mpl,nam,geom)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist
type(geom_type),intent(in) :: geom       ! Geometry

! Local variables
integer :: ncid,grid_hash
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'obsop_read'

! Create file
write(filename,'(a,a,i6.6,a,i6.6)') trim(nam%prefix),'_obs_',mpl%nproc,'-',mpl%myproc
call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))

! Check grid hash
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'grid_hash',grid_hash))
if (grid_hash/=geom%grid_hash) call mpl%abort(subr,'wrong grid hash')

! Get attributes
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'nc0b',obsop%nc0b))
call mpl%ncerr(subr,nf90_get_att(ncid,nf90_global,'nobsa',obsop%nobsa))

! Read interpolation
obsop%h%prefix = 'o'
call obsop%h%read(mpl,ncid)

! Read communication
obsop%com%prefix = 'com'
call obsop%com%read(mpl,ncid)

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

end subroutine obsop_read

!----------------------------------------------------------------------
! Subroutine: obsop_write
! Purpose: write observations locations
!----------------------------------------------------------------------
subroutine obsop_write(obsop,mpl,nam,geom)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist
type(geom_type),intent(in) :: geom       ! Geometry

! Local variables
integer :: ncid
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'obsop_write'

! Create file
write(filename,'(a,a,i6.6,a,i6.6)') trim(nam%prefix),'_obs_',mpl%nproc,'-',mpl%myproc
call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Write grid hash
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'grid_hash',geom%grid_hash))

! Write namelist parameters
call nam%write(mpl,ncid)

! Write attributes
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'nc0b',obsop%nc0b))
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'nobsa',obsop%nobsa))

! End definition mode
call mpl%ncerr(subr,nf90_enddef(ncid))

! Write interpolation
call obsop%h%write(mpl,ncid)

! Write communication
call obsop%com%write(mpl,ncid)

! Close file
call mpl%ncerr(subr,nf90_close(ncid))

end subroutine obsop_write

!----------------------------------------------------------------------
! Subroutine: obsop_from
! Purpose: copy observation operator data
!----------------------------------------------------------------------
subroutine obsop_from(obsop,nobsa,lonobs,latobs)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop    ! Observation operator data
integer,intent(in) :: nobsa                 ! Number of observations
real(kind_real),intent(in) :: lonobs(nobsa) ! Observations longitudes (in degrees)
real(kind_real),intent(in) :: latobs(nobsa) ! Observations latitudes (in degrees)

! Local variables
integer :: iobsa

! Get size
obsop%nobsa = nobsa

! Release memory
call obsop%dealloc

! Allocation
allocate(obsop%lonobs(obsop%nobsa))
allocate(obsop%latobs(obsop%nobsa))

if (obsop%nobsa>0) then
   ! Copy
   obsop%lonobs = lonobs*deg2rad
   obsop%latobs = latobs*deg2rad

   ! Enforce correct bounds
   do iobsa=1,obsop%nobsa
      call lonlatmod(obsop%lonobs(iobsa),obsop%latobs(iobsa))
   end do
end if

end subroutine obsop_from

!----------------------------------------------------------------------
! Subroutine: obsop_run_obsop
! Purpose: observation operator driver
!----------------------------------------------------------------------
subroutine obsop_run_obsop(obsop,mpl,rng,nam,geom)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(rng_type),intent(inout) :: rng      ! Random number generator
type(nam_type),intent(in) :: nam         ! Namelist
type(geom_type),intent(in) :: geom       ! Geometry

! Local variables
integer :: iobsa,iproc,i_s,ic0,ic0u,jc0u,ic0b,ic0a,nobsa_eff
integer :: nobs_eff,nn_index(1),proc_to_nobsa(mpl%nproc),proc_to_nobsa_eff(mpl%nproc)
integer :: c0u_to_c0b(geom%nc0u)
integer,allocatable :: c0b_to_c0(:)
real(kind_real) :: nn_dist(1),N_max,C_max
logical :: maskobsa(obsop%nobsa),lcheck_nc0b(geom%nc0u)
character(len=1024),parameter :: subr = 'obsop_run_obsop'

! Check that universe is global
if (any(.not.geom%myuniverse)) call mpl%abort(subr,'universe should be global for obsop')

! Check whether observations are inside the mesh
if (obsop%nobsa>0) then
   do iobsa=1,obsop%nobsa
      call geom%mesh_c0u%inside(mpl,obsop%lonobs(iobsa),obsop%latobs(iobsa),maskobsa(iobsa))
      if (.not.maskobsa(iobsa)) then
         ! Check for very close points
         call geom%tree_c0u%find_nearest_neighbors(obsop%lonobs(iobsa),obsop%latobs(iobsa),1,nn_index,nn_dist)
         if (nn_dist(1)<rth) maskobsa(iobsa) = .true.
      end if
   end do
   nobsa_eff = count(maskobsa)
else
   nobsa_eff = 0
end if

! Get global number of observations
call mpl%f_comm%allgather(obsop%nobsa,proc_to_nobsa)
call mpl%f_comm%allgather(nobsa_eff,proc_to_nobsa_eff)
obsop%nobs = sum(proc_to_nobsa)
nobs_eff = sum(proc_to_nobsa_eff)

! Print input
write(mpl%info,'(a7,a)') '','Number of observations / valid observations per MPI task:'
call mpl%flush
do iproc=1,mpl%nproc
   write(mpl%info,'(a10,a,i3,a,i8,a,i8)') '','Task ',iproc,': ',proc_to_nobsa(iproc),' / ',proc_to_nobsa_eff(iproc)
   call mpl%flush
end do
write(mpl%info,'(a10,a,i8,a,i8)') '','Total   : ',obsop%nobs,' / ',nobs_eff
call mpl%flush

! Compute interpolation
obsop%h%prefix = 'o'
write(mpl%info,'(a7,a)') '','Single level:'
call mpl%flush
call obsop%h%interp(mpl,rng,nam,geom,0,geom%nc0u,geom%lon_c0u,geom%lat_c0u,geom%gmask_hor_c0u,obsop%nobsa,obsop%lonobs, &
 & obsop%latobs,maskobsa,10)

! Define halo B
lcheck_nc0b = .false.
do ic0a=1,geom%nc0a
   ic0u = geom%c0a_to_c0u(ic0a)
   lcheck_nc0b(ic0u) = .true.
end do
do iobsa=1,obsop%nobsa
   do i_s=1,obsop%h%n_s
      jc0u = obsop%h%col(i_s)
      lcheck_nc0b(jc0u) = .true.
   end do
end do
obsop%nc0b = count(lcheck_nc0b)

! Allocation
allocate(c0b_to_c0(obsop%nc0b))

! Global-local conversion for halo B
c0u_to_c0b = mpl%msv%vali
ic0b = 0
do ic0u=1,geom%nc0u
   if (lcheck_nc0b(ic0u)) then
      ic0b = ic0b+1
      ic0 = geom%c0u_to_c0(ic0u)
      c0b_to_c0(ic0b) = ic0
      c0u_to_c0b(ic0u) = ic0b
   end if
end do

! Local interpolation source
obsop%h%n_src = obsop%nc0b
do i_s=1,obsop%h%n_s
   obsop%h%col(i_s) = c0u_to_c0b(obsop%h%col(i_s))
end do

! Setup communications
call obsop%com%setup(mpl,'com',geom%nc0a,obsop%nc0b,geom%nc0,geom%c0a_to_c0,c0b_to_c0)

! Compute scores, only if there observations present globally
if ( nobs_eff > 0 ) then
  call mpl%f_comm%allreduce(real(obsop%com%nhalo,kind_real),C_max,fckit_mpi_max())
  C_max = C_max/(3.0*real(nobs_eff,kind_real)/real(mpl%nproc,kind_real))
  N_max = real(maxval(proc_to_nobsa_eff),kind_real)/(real(nobs_eff,kind_real)/real(mpl%nproc,kind_real))

  ! Print results
  write(mpl%info,'(a7,a,f5.1,a)') '','Observation repartition imbalance: ',100.0*real(maxval(proc_to_nobsa_eff) &
  & -minval(proc_to_nobsa_eff),kind_real)/(real(sum(proc_to_nobsa_eff),kind_real)/real(mpl%nproc,kind_real)),' %'
  call mpl%flush
  write(mpl%info,'(a7,a,i8,a,i8,a,i8)') '','Number of grid points / halo size / number of received values: ', &
  & obsop%com%nred,' / ',obsop%com%next,' / ',obsop%com%nhalo
  call mpl%flush
  write(mpl%info,'(a7,a,f10.2,a,f10.2)') '','Scores (N_max / C_max):',N_max,' / ',C_max
  call mpl%flush
end if

! Write observation operator
if (nam%write_obsop) call obsop%write(mpl,nam,geom)

! Release memory
deallocate(c0b_to_c0)

end subroutine obsop_run_obsop

!----------------------------------------------------------------------
! Subroutine: obsop_run_obsop_tests
! Purpose: observation operator tests driver
!----------------------------------------------------------------------
subroutine obsop_run_obsop_tests(obsop,mpl,nam,rng,geom)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist
type(rng_type),intent(inout) :: rng      ! Random number generator
type(geom_type),intent(in) :: geom       ! Geometry

if (nam%check_adjoints) then
   ! Test adjoints
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Test observation operator adjoint'
   call mpl%flush
   call obsop%test_adjoint(mpl,rng,geom)
end if

if (nam%check_obsop.and.allocated(obsop%obsa_to_obs)) then
   ! Test precision
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Test observation operator precision'
   call mpl%flush
   call obsop%test_accuracy(mpl,geom)
end if

end subroutine obsop_run_obsop_tests

!----------------------------------------------------------------------
! Subroutine: obsop_apply
! Purpose: observation operator interpolation
!----------------------------------------------------------------------
subroutine obsop_apply(obsop,mpl,geom,fld,obs)

implicit none

! Passed variables
class(obsop_type),intent(in) :: obsop                    ! Observation operator data
type(mpl_type),intent(inout) :: mpl                      ! MPI data
type(geom_type),intent(in) :: geom                       ! Geometry
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0)    ! Field
real(kind_real),intent(out) :: obs(obsop%nobsa,geom%nl0) ! Observations columns

! Local variables
integer :: il0
real(kind_real) :: fld_ext(obsop%nc0b,geom%nl0)

! Halo extension
call obsop%com%ext(mpl,geom%nl0,fld,fld_ext)

if (obsop%nobsa>0) then
   ! Horizontal interpolation
   !$omp parallel do schedule(static) private(il0) shared(geom,obsop,mpl,fld_ext,obs)
   do il0=1,geom%nl0
      call obsop%h%apply(mpl,fld_ext(:,il0),obs(:,il0))
   end do
   !$omp end parallel do
end if

end subroutine obsop_apply

!----------------------------------------------------------------------
! Subroutine: obsop_apply_ad
! Purpose: observation operator interpolation adjoint
!----------------------------------------------------------------------
subroutine obsop_apply_ad(obsop,mpl,geom,obs,fld)

implicit none

! Passed variables
class(obsop_type),intent(in) :: obsop                   ! Observation operator data
type(mpl_type),intent(inout) :: mpl                     ! MPI data
type(geom_type),intent(in) :: geom                      ! Geometry
real(kind_real),intent(in) :: obs(obsop%nobsa,geom%nl0) ! Observations columns
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0)  ! Field

! Local variables
integer :: il0
real(kind_real) :: fld_ext(obsop%nc0b,geom%nl0)

if (obsop%nobsa>0) then
   ! Horizontal interpolation
   !$omp parallel do schedule(static) private(il0) shared(geom,obsop,mpl,obs,fld_ext)
   do il0=1,geom%nl0
      call obsop%h%apply_ad(mpl,obs(:,il0),fld_ext(:,il0))
   end do
   !$omp end parallel do
else
   ! No observation on this task
   fld_ext = 0.0
end if

! Halo reduction
call obsop%com%red(mpl,geom%nl0,fld_ext,fld)

end subroutine obsop_apply_ad

!----------------------------------------------------------------------
! Subroutine: obsop_test_adjoint
! Purpose: test observation operator adjoints accuracy
!----------------------------------------------------------------------
subroutine obsop_test_adjoint(obsop,mpl,rng,geom)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(rng_type),intent(inout) :: rng      ! Random number generator
type(geom_type),intent(in) :: geom       ! Geometry

! Local variables
real(kind_real) :: sum1,sum2_loc,sum2
real(kind_real) :: fld(geom%nc0a,geom%nl0),fld_save(geom%nc0a,geom%nl0)
real(kind_real) :: yobs(obsop%nobsa,geom%nl0),yobs_save(obsop%nobsa,geom%nl0)

! Generate random fields
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld_save)
if (obsop%nobsa>0) call rng%rand_real(0.0_kind_real,1.0_kind_real,yobs_save)

! Apply direct and adjoint obsservation operators
call obsop%apply(mpl,geom,fld_save,yobs)
call obsop%apply_ad(mpl,geom,yobs_save,fld)

! Compute adjoint test
call mpl%dot_prod(fld,fld_save,sum1)
if (obsop%nobsa>0) then
   sum2_loc = sum(yobs*yobs_save,mask=mpl%msv%isnot(yobs))
else
   sum2_loc = 0.0
end if
call mpl%f_comm%allreduce(sum2_loc,sum2,fckit_mpi_sum())

! Print results
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Observation operator adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush

end subroutine obsop_test_adjoint

!----------------------------------------------------------------------
! Subroutine: obsop_test_accuracy
! Purpose: test observation operator accuracy
!----------------------------------------------------------------------
subroutine obsop_test_accuracy(obsop,mpl,geom)

implicit none

! Passed variables
class(obsop_type),intent(inout) :: obsop ! Observation operator data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(geom_type),intent(in) :: geom       ! Geometry

! Local variables
integer :: ic0a,iobsa
integer :: iprocmax(1),iobsamax(1)
real(kind_real) :: lonmax,latmax,ylonmax,ylatmax
real(kind_real) :: norm,distmin,distmax,distsum
real(kind_real) :: norm_tot,distmin_tot,proc_to_distmax(mpl%nproc),distsum_tot
real(kind_real) :: lon(geom%nc0a,geom%nl0),lat(geom%nc0a,geom%nl0)
real(kind_real) :: ylon(obsop%nobsa,geom%nl0),ylat(obsop%nobsa,geom%nl0)
real(kind_real) :: dist(obsop%nobsa)
character(len=1024),parameter :: subr = 'obsop_test_accuracy'

! Initialization
do ic0a=1,geom%nc0a
   lon(ic0a,:) = geom%lon_c0a(ic0a)
   lat(ic0a,:) = geom%lat_c0a(ic0a)
end do

! Apply obsop
call obsop%apply(mpl,geom,lon,ylon)
call obsop%apply(mpl,geom,lat,ylat)

if (obsop%nobsa>0) then
   ! Remove points close to the longitude discontinuity and to the poles
   dist = mpl%msv%valr
   do iobsa=1,obsop%nobsa
      if (mpl%msv%isnot(ylon(iobsa,1)).and.mpl%msv%isnot(ylat(iobsa,1))) then
         if ((abs(obsop%lonobs(iobsa))<0.8*pi).and.(abs(obsop%latobs(iobsa))<0.4*pi)) then
            call sphere_dist(ylon(iobsa,1),ylat(iobsa,1),obsop%lonobs(iobsa),obsop%latobs(iobsa),dist(iobsa))
            dist(iobsa) = dist(iobsa)*reqkm
         end if
      end if
   end do
   norm = real(count(mpl%msv%isnot(dist)),kind_real)
   if (norm>0) then
      distmin = minval(dist,mask=mpl%msv%isnot(dist))
      distmax = maxval(dist,mask=mpl%msv%isnot(dist))
      distsum = sum(dist,mask=mpl%msv%isnot(dist))
   else
      distmin = huge_real
      distmax = 0.0
      distsum = 0.0
   end if
else
   ! No observation on this task
   norm = 0
   distmin = huge_real
   distmax = 0.0
   distsum = 0.0
end if

! Gather results
call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
call mpl%f_comm%allreduce(distmin,distmin_tot,fckit_mpi_min())
call mpl%f_comm%allgather(distmax,proc_to_distmax)
call mpl%f_comm%allreduce(distsum,distsum_tot,fckit_mpi_sum())

! Maximum error detail
iprocmax = maxloc(proc_to_distmax)
if (iprocmax(1)==mpl%myproc) then
   iobsamax = maxloc(dist,mask=mpl%msv%isnot(dist))
   lonmax = obsop%lonobs(iobsamax(1))
   latmax = obsop%latobs(iobsamax(1))
   ylonmax = ylon(iobsamax(1),1)
   ylatmax = ylat(iobsamax(1),1)
end if

! Broadcast results
call mpl%f_comm%broadcast(lonmax,iprocmax(1)-1)
call mpl%f_comm%broadcast(latmax,iprocmax(1)-1)
call mpl%f_comm%broadcast(ylonmax,iprocmax(1)-1)
call mpl%f_comm%broadcast(ylatmax,iprocmax(1)-1)

! Print results
if (norm_tot>0.0) then
   write(mpl%info,'(a7,a,f10.2,a,f10.2,a,f10.2,a)') '','Interpolation error (min/mean/max): ',distmin_tot*reqkm, &
 & ' km / ',distsum_tot*reqkm/norm_tot,' km / ',maxval(proc_to_distmax),' km'
   call mpl%flush
   write(mpl%info,'(a7,a)') '','Max. interpolation error location (lon/lat): '
   call mpl%flush
   write(mpl%info,'(a10,a14,f10.2,a,f10.2,a)') '','Observation:  ',lonmax*rad2deg,' deg. / ' ,latmax*rad2deg,' deg.'
   call mpl%flush
   write(mpl%info,'(a10,a14,f10.2,a,f10.2,a)') '','Interpolation:',ylonmax*rad2deg,' deg. / ' ,ylatmax*rad2deg,' deg.'
   call mpl%flush
else
   call mpl%abort(subr,'all observations are out of the test windows')
end if

end subroutine obsop_test_accuracy

end module type_obsop
