!----------------------------------------------------------------------
! Module: type_bump
! Purpose: BUMP derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_bump

use atlas_module, only: atlas_field,atlas_fieldset,atlas_integer,atlas_real,atlas_functionspace
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm,fckit_mpi_sum,fckit_mpi_min,fckit_mpi_max
use tools_atlas, only: create_atlas_fieldset,create_atlas_function_space,atlas_to_fld,fld_to_atlas
use tools_const, only: req,deg2rad
use tools_func, only: sphere_dist,lct_r2d
use tools_kinds,only: kind_int,kind_real
use tools_repro,only: repro
use type_bpar, only: bpar_type
use type_cmat, only: cmat_type
use type_cv, only: cv_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_hdiag, only: hdiag_type
use type_io, only: io_type
use type_lct, only: lct_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_nicas, only: nicas_type
use type_obsop, only: obsop_type
use type_rng, only: rng_type
use type_var, only: var_type
use type_vbal, only: vbal_type

implicit none

! BUMP derived type
type bump_type
   type(bpar_type) :: bpar
   type(cmat_type) :: cmat
   type(ens_type) :: ens1
   type(ens_type) :: ens1u
   type(ens_type) :: ens2
   type(geom_type) :: geom
   type(hdiag_type) :: hdiag
   type(io_type) :: io
   type(lct_type) :: lct
   type(mpl_type) :: mpl
   type(nam_type) :: nam
   type(nicas_type) :: nicas
   type(obsop_type) :: obsop
   type(rng_type) :: rng
   type(var_type) :: var
   type(vbal_type) :: vbal
   real(kind_real),allocatable :: fld_uv(:,:,:,:)
contains
   procedure :: create => bump_create
   procedure :: setup => bump_setup
   procedure :: setup_online => bump_setup_online_deprecated
   procedure :: run_drivers => bump_run_drivers
   procedure :: add_member => bump_add_member
   procedure :: remove_member => bump_remove_member
   procedure :: apply_vbal => bump_apply_vbal
   procedure :: apply_vbal_inv => bump_apply_vbal_inv
   procedure :: apply_vbal_ad => bump_apply_vbal_ad
   procedure :: apply_vbal_inv_ad => bump_apply_vbal_inv_ad
   procedure :: apply_stddev => bump_apply_stddev
   procedure :: apply_stddev_inv => bump_apply_stddev_inv
   procedure :: bump_apply_nicas
   procedure :: bump_apply_nicas_deprecated
   generic :: apply_nicas => bump_apply_nicas,bump_apply_nicas_deprecated
   procedure :: get_cv_size => bump_get_cv_size
   procedure :: bump_apply_nicas_sqrt
   procedure :: bump_apply_nicas_sqrt_deprecated
   generic :: apply_nicas_sqrt => bump_apply_nicas_sqrt,bump_apply_nicas_sqrt_deprecated
   procedure :: apply_nicas_sqrt_ad => bump_apply_nicas_sqrt_ad
   procedure :: randomize => bump_randomize
   procedure :: bump_apply_obsop
   procedure :: bump_apply_obsop_deprecated
   generic :: apply_obsop => bump_apply_obsop,bump_apply_obsop_deprecated
   procedure :: bump_apply_obsop_ad
   procedure :: bump_apply_obsop_ad_deprecated
   generic :: apply_obsop_ad => bump_apply_obsop_ad,bump_apply_obsop_ad_deprecated
   procedure :: get_parameter => bump_get_parameter
   procedure :: copy_to_field => bump_copy_to_field
   procedure :: test_get_parameter => bump_test_get_parameter
   procedure :: bump_set_parameter
   procedure :: bump_set_parameter_deprecated
   generic :: set_parameter => bump_set_parameter,bump_set_parameter_deprecated
   procedure :: copy_from_field => bump_copy_from_field
   procedure :: test_set_parameter => bump_test_set_parameter
   procedure :: test_apply_interfaces => bump_test_apply_interfaces
   procedure :: partial_dealloc => bump_partial_dealloc
   procedure :: dealloc => bump_dealloc
   final :: dummy
end type bump_type

integer,parameter :: dmsvali = -999           ! Default missing value for integers
real(kind_real),parameter :: dmsvalr = -999.0 ! Default missing value for reals
logical,parameter :: print_member = .false.   ! Print info when adding member to BUMP
logical,parameter :: write_member = .false.   ! Write file when adding member to BUMP

private
public :: bump_type
public :: bump_registry

! BUMP registry
#define LISTED_TYPE bump_type

! Linked list interface - defines registry_t type
#include "saber/util/linkedList_i.f"

! Global registry
type(registry_t) :: bump_registry

contains

!----------------------------------------------------------------------
! Linked list implementation
!----------------------------------------------------------------------
#include "saber/util/linkedList_c.f"

!----------------------------------------------------------------------
! Subroutine: bump_create
! Purpose: create
!----------------------------------------------------------------------
subroutine bump_create(bump,comm,afunctionspace,afieldset,conf,grid)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                 ! BUMP
type(fckit_mpi_comm),intent(in) :: comm                ! FCKIT MPI communicator wrapper
type(atlas_functionspace),intent(in) :: afunctionspace ! ATLAS function space
type(atlas_fieldset),intent(in) :: afieldset           ! ATLAS fieldset  (containing geometry features: area, vunit, gmask, smask, wind)
type(fckit_configuration),intent(in) :: conf           ! FCKIT configuration
type(fckit_configuration),intent(in) :: grid           ! FCKIT grid configuration

! Local variables
integer :: lmsvali, llunit
real(kind_real) :: lmsvalr

! Initialize namelist
call bump%nam%init(comm%size())

! Read configuration
call bump%nam%from_conf(conf)

! Read grid configuration
call bump%nam%from_conf(grid)

! Set missing values
lmsvali = dmsvali
lmsvalr = dmsvalr
if (conf%has('msvali')) call conf%get_or_die('msvali',lmsvali)
if (conf%has('msvalr')) call conf%get_or_die('msvalr',lmsvalr)

! Set log unit
llunit = lmsvali
if (conf%has('lunit')) call conf%get_or_die('lunit',llunit)

! Setup BUMP
call bump%setup(comm,afunctionspace,afieldset,lunit=llunit,msvali=lmsvali,msvalr=lmsvalr)

end subroutine bump_create

!----------------------------------------------------------------------
! Subroutine: bump_setup
! Purpose: setup
!----------------------------------------------------------------------
subroutine bump_setup(bump,f_comm,afunctionspace,afieldset,nobs,lonobs,latobs,lunit,msvali,msvalr)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                 ! BUMP
type(fckit_mpi_comm),intent(in) :: f_comm              ! FCKIT MPI communicator wrapper
type(atlas_functionspace),intent(in) :: afunctionspace ! ATLAS functionspace
type(atlas_fieldset),intent(in),optional :: afieldset  ! ATLAS fieldset (containing geometry features: area, vunit, gmask, smask, wind)
integer,intent(in),optional :: nobs                    ! Number of observations
real(kind_real),intent(in),optional :: lonobs(:)       ! Observations longitude (in degrees)
real(kind_real),intent(in),optional :: latobs(:)       ! Observations latitude (in degrees)
integer,intent(in),optional :: lunit                   ! Listing unit
integer,intent(in),optional :: msvali                  ! Missing value for integers
real(kind_real),intent(in),optional :: msvalr          ! Missing value for reals

! Local variables
integer :: iv,its
real(kind_real),allocatable :: fld_uv(:,:)
real(kind_real),pointer :: real_ptr(:,:)
character(len=1024) :: fieldname
character(len=1024),parameter :: subr = 'bump_setup'
type(atlas_field) :: afield

! Initialize MPL
call bump%mpl%init(f_comm)

! Set missing values
bump%mpl%msv%vali = dmsvali
bump%mpl%msv%valr = dmsvalr
if (present(msvali)) bump%mpl%msv%vali = msvali
if (present(msvalr)) bump%mpl%msv%valr = msvalr

! Set reproducibility parameter
repro = bump%nam%repro

! Initialize listing
bump%mpl%lunit = bump%mpl%msv%vali
if (present(lunit)) bump%mpl%lunit = lunit
bump%mpl%verbosity = bump%nam%verbosity
if (bump%nam%colorlog) then
   bump%mpl%black = char(27)//'[0;0m'
   bump%mpl%green = char(27)//'[0;32m'
   bump%mpl%peach = char(27)//'[1;91m'
   bump%mpl%aqua = char(27)//'[1;36m'
   bump%mpl%purple = char(27)//'[1;35m'
   bump%mpl%err = char(27)//'[0;37;41;1m'
   bump%mpl%wng = char(27)//'[0;37;42;1m'
else
   bump%mpl%black = ' '
   bump%mpl%green = ' '
   bump%mpl%peach = ' '
   bump%mpl%aqua = ' '
   bump%mpl%purple = ' '
   bump%mpl%err = ' '
   bump%mpl%wng = ' '
end if

! Header
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- You are running the BUMP library ------------------------------'
call bump%mpl%flush

! Check namelist parameters
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Check namelist parameters'
call bump%mpl%flush
call bump%nam%check(bump%mpl)
call bump%nam%write(bump%mpl)

! Write parallel setup
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a,i3,a,i2,a)') '--- Parallelization with ',bump%mpl%nproc,' MPI tasks and ', &
 & bump%mpl%nthread,' OpenMP threads'
call bump%mpl%flush

! Initialize random number generator
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Initialize random number generator'
call bump%mpl%flush
call bump%rng%init(bump%mpl,bump%nam)

! Initialize allocation flags
bump%cmat%allocated = .false.
bump%lct%allocated = .false.
bump%nicas%allocated = .false.

! Initialize geometry
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Initialize geometry'
call bump%mpl%flush
if (present(afieldset)) then
   call bump%geom%setup(bump%mpl,bump%rng,bump%nam,afunctionspace,afieldset)
else
   call bump%geom%setup(bump%mpl,bump%rng,bump%nam,afunctionspace)
end if
if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

! Initialize fields output
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Initialize fields output'
call bump%mpl%flush
call bump%io%init(bump%mpl,bump%nam,bump%geom)
if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

! Initialize block parameters
write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
call bump%mpl%flush
write(bump%mpl%info,'(a)') '--- Initialize block parameters'
call bump%mpl%flush
call bump%bpar%alloc(bump%nam,bump%geom)
call bump%bpar%init(bump%mpl,bump%nam,bump%geom)

if (bump%nam%ens1_ne>0) then
   ! Initialize ensemble 1
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Initialize ensemble 1'
   call bump%mpl%flush
   call bump%ens1%alloc(bump%nam,bump%geom,bump%nam%ens1_ne,bump%nam%ens1_nsub)
else
   call bump%ens1%set_att(bump%nam%ens1_ne,bump%nam%ens1_nsub)
end if

if (bump%nam%ens2_ne>0) then
   ! Initialize ensemble 2
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Initialize ensemble 2'
   call bump%mpl%flush
   call bump%ens2%alloc(bump%nam,bump%geom,bump%nam%ens2_ne,bump%nam%ens2_nsub)
else
   call bump%ens2%set_att(bump%nam%ens2_ne,bump%nam%ens2_nsub)
end if

if (bump%nam%new_cortrack.or.(trim(bump%nam%adv_type)=='wind').or.(trim(bump%nam%adv_type)=='windmax')) then
   ! Initialize wind fields
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Initialize wind fields'
   call bump%mpl%flush

   ! Check that afieldset is present
   if (.not.present(afieldset)) call bump%mpl%abort(subr,'afieldset required to initialize wind fields')

   ! Allocation
   allocate(fld_uv(bump%geom%nmga,bump%geom%nl0))
   allocate(bump%fld_uv(bump%geom%nc0a,bump%geom%nl0,2,bump%nam%nts))

   ! Get field from ATLAS fieldset
   do its=1,bump%nam%nts
      do iv=1,2
         ! Get field
         fieldname = trim(bump%nam%wind_variables(iv))//'_'//trim(bump%nam%timeslots(its))
         afield = afieldset%field(trim(fieldname))

         ! Get data
         call afield%data(real_ptr)

         ! Copy to BUMP
         fld_uv = transpose(real_ptr)
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_uv,bump%fld_uv(:,:,iv,its))
      end do
   end do

   ! Normalization
   bump%fld_uv = bump%fld_uv/req

   ! Release memory
   deallocate(fld_uv)
end if

if (present(nobs)) then
   ! Check arguments consistency
   if ((.not.present(lonobs)).or.(.not.present(latobs))) call bump%mpl%abort(subr,'lonobs and latobs are missing')

   ! Check sizes consistency
   if (size(lonobs)/=nobs) call bump%mpl%abort(subr,'wrong size for lonobs')
   if (size(latobs)/=nobs) call bump%mpl%abort(subr,'wrong size for latobs')

   ! Initialize observations locations
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Initialize observations locations'
   call bump%mpl%flush
   call bump%obsop%from(nobs,lonobs,latobs)
end if

end subroutine bump_setup

!----------------------------------------------------------------------
! Subroutine: bump_setup_online_deprecated
! Purpose: online setup (deprecated)
!----------------------------------------------------------------------
subroutine bump_setup_online_deprecated(bump,f_comm,nmga,nl0,nv,nts,lon,lat,area,vunit,gmask,smask,ens1_ne,ens1_nsub, &
 & ens2_ne,ens2_nsub,nobs,lonobs,latobs,lunit,msvali,msvalr)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump            ! BUMP
type(fckit_mpi_comm),intent(in) :: f_comm         ! FCKIT MPI communicator wrapper
integer,intent(in) :: nmga                        ! Halo A size
integer,intent(in) :: nl0                         ! Number of levels in subset Sl0
integer,intent(in) :: nv                          ! Number of variables
integer,intent(in) :: nts                         ! Number of time slots
real(kind_real),intent(in) :: lon(nmga)           ! Longitude (in degrees: -180 to 180)
real(kind_real),intent(in) :: lat(nmga)           ! Latitude (in degrees: -90 to 90)
real(kind_real),intent(in) :: area(nmga)          ! Area (in m^2)
real(kind_real),intent(in) :: vunit(nmga,nl0)     ! Vertical unit
logical,intent(in) :: gmask(nmga,nl0)             ! Geometry mask
logical,intent(in),optional :: smask(nmga,nl0)    ! Sampling mask
integer,intent(in),optional :: ens1_ne            ! Ensemble 1 size
integer,intent(in),optional :: ens1_nsub          ! Ensemble 1 number of sub-ensembles
integer,intent(in),optional :: ens2_ne            ! Ensemble 2 size
integer,intent(in),optional :: ens2_nsub          ! Ensemble 2 size of sub-ensembles
integer,intent(in),optional :: nobs               ! Number of observations
real(kind_real),intent(in),optional :: lonobs(:)  ! Observations longitude (in degrees: -180 to 180)
real(kind_real),intent(in),optional :: latobs(:)  ! Observations latitude (in degrees: -90 to 90)
integer,intent(in),optional :: lunit              ! Listing unit
integer,intent(in),optional :: msvali             ! Missing value for integers
real(kind_real),intent(in),optional :: msvalr     ! Missing value for reals

! Local variables
integer :: lnobs,lmsvali,llunit,imga,il0,iv,its
integer(kind_int),pointer :: int_ptr_2(:,:)
real(kind_real) :: lmsvalr
real(kind_real),allocatable :: llonobs(:),llatobs(:)
real(kind_real),pointer :: real_ptr_1(:),real_ptr_2(:,:)
character(len=1024),parameter :: subr = 'bump_setup_online_deprecated'
type(atlas_field) :: afield
type(atlas_fieldset) :: afieldset
type(atlas_functionspace) :: afunctionspace

! Force optional parameters
lnobs = 0
lmsvali = dmsvali
lmsvalr = dmsvalr
if (present(nobs)) lnobs = nobs
allocate(llonobs(lnobs))
allocate(llatobs(lnobs))
if (present(lonobs).and.present(latobs)) then
  llonobs = lonobs
  llatobs = latobs
end if
if (present(msvali)) lmsvali = msvali
if (present(msvalr)) lmsvalr = msvalr
llunit = lmsvali
if (present(lunit)) llunit = lunit

! Set namelist parameters
bump%nam%nl = nl0
bump%nam%nv = nv
do iv=1,bump%nam%nv
   write(bump%nam%variables(iv),'(a,i2.2)') 'var_',iv
end do
bump%nam%nts = nts
do its=1,bump%nam%nts
   write(bump%nam%timeslots(its),'(a,i2.2)') 'ts_',its
end do
bump%nam%lev2d = 'first'
if (present(ens1_ne)) bump%nam%ens1_ne = ens1_ne
if (present(ens1_nsub)) bump%nam%ens1_nsub = ens1_nsub
if (present(ens2_ne)) bump%nam%ens2_ne = ens2_ne
if (present(ens2_nsub)) bump%nam%ens2_nsub = ens2_nsub

! Create ATLAS function space
call create_atlas_function_space(nmga,lon*deg2rad,lat*deg2rad,afunctionspace)

! Create ATLAS fieldset with empty fields
call create_atlas_fieldset(afunctionspace,bump%nam%nl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts), &
 & afieldset)

! Set geometry

! Set area
afield = afunctionspace%create_field(name='area',kind=atlas_real(kind_real),levels=0)
call afield%data(real_ptr_1)
real_ptr_1 = area
call afieldset%add(afield)
call afield%final()

! Set vertical unit
afield = afunctionspace%create_field(name='vunit',kind=atlas_real(kind_real),levels=bump%nam%nl)
call afield%data(real_ptr_2)
real_ptr_2 = transpose(vunit)
call afieldset%add(afield)
call afield%final()

! Set geometry mask
afield = afunctionspace%create_field(name='gmask',kind=atlas_integer(kind_int),levels=bump%nam%nl)
call afield%data(int_ptr_2)
do il0=1,bump%nam%nl
   do imga=1,nmga
      if (gmask(imga,il0)) then
         int_ptr_2(il0,imga) = 1
      else
         int_ptr_2(il0,imga) = 0
      end if
   end do
end do
call afieldset%add(afield)
call afield%final()

if (present(smask)) then
   ! Set sampling mask
   afield = afunctionspace%create_field(name='smask',kind=atlas_integer(kind_int),levels=bump%nam%nl)
   call afield%data(int_ptr_2)
   do il0=1,bump%nam%nl
      do imga=1,nmga
         if (smask(imga,il0)) then
            int_ptr_2(il0,imga) = 1
         else
            int_ptr_2(il0,imga) = 0
         end if
      end do
   end do
   call afieldset%add(afield)
   call afield%final()
end if

! BUMP setup
call bump%setup(f_comm,afunctionspace,afieldset=afieldset,nobs=lnobs,lonobs=llonobs,latobs=llatobs, &
 & lunit=llunit,msvali=lmsvali,msvalr=lmsvalr)

! Deprecation warning
call bump%mpl%warning(subr,'this interface is deprecated, consider using the ATLAS-based interface')

end subroutine bump_setup_online_deprecated

!----------------------------------------------------------------------
! Subroutine: bump_run_drivers
! Purpose: run drivers
!----------------------------------------------------------------------
subroutine bump_run_drivers(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP

if (bump%nam%ens1_ne>0) then
   ! Finalize ensemble 1
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Finalize ensemble 1'
   call bump%mpl%flush
   call bump%ens1%remove_mean
end if

if (bump%nam%ens2_ne>0) then
   ! Finalize ensemble 2
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Finalize ensemble 2'
   call bump%mpl%flush
   call bump%ens2%remove_mean
end if

if (bump%nam%new_normality) then
   ! Run normality tests
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run normality tests'
   call bump%mpl%flush
   call bump%ens1%normality(bump%mpl,bump%nam,bump%geom,bump%io)
end if

if (bump%nam%new_cortrack) then
   ! Run correlation tracker
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run correlation tracker'
   call bump%mpl%flush
   if (allocated(bump%fld_uv)) then
      call bump%ens1%cortrack(bump%mpl,bump%rng,bump%nam,bump%geom,bump%io,bump%fld_uv)
   else
      call bump%ens1%cortrack(bump%mpl,bump%rng,bump%nam,bump%geom,bump%io)
   end if
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
end if

if (bump%nam%new_corstats) then
   ! Run correlation statistics
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run correlation statistics'
   call bump%mpl%flush
   call bump%ens1%corstats(bump%mpl,bump%rng,bump%nam,bump%geom)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
end if

if (bump%nam%new_vbal) then
   ! Run vertical balance driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run vertical balance driver'
   call bump%mpl%flush
   call bump%vbal%run_vbal(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%ens1,bump%ens1u)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
elseif (bump%nam%load_vbal) then
   ! Read vertical balance
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Read vertical balance'
   call bump%mpl%flush
   call bump%vbal%read(bump%mpl,bump%nam,bump%geom,bump%bpar)
end if

if (bump%nam%new_vbal.or.bump%nam%load_vbal) then
   ! Run vertical balance tests driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run vertical balance tests driver'
   call bump%mpl%flush
   call bump%vbal%run_vbal_tests(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
end if

if (bump%nam%new_var) then
   ! Run variance driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run variance driver'
   call bump%mpl%flush
   call bump%var%run_var(bump%mpl,bump%rng,bump%nam,bump%geom,bump%ens1,bump%io)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
elseif (bump%nam%load_var) then
   ! Read variance
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Read variance'
   call bump%mpl%flush
   call bump%var%read(bump%mpl,bump%nam,bump%geom,bump%io)
end if

if (bump%nam%new_hdiag) then
   ! Run HDIAG driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run HDIAG driver'
   call bump%mpl%flush
   if ((trim(bump%nam%method)=='hyb-rnd').or.(trim(bump%nam%method)=='dual-ens')) then
      if (allocated(bump%fld_uv)) then
         call bump%hdiag%run_hdiag(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1,ens2=bump%ens2, &
 & fld_uv=bump%fld_uv)
      else
         call bump%hdiag%run_hdiag(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1,ens2=bump%ens2)
      end if
   else
      if (allocated(bump%fld_uv)) then
         call bump%hdiag%run_hdiag(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1,fld_uv=bump%fld_uv)
      else
         call bump%hdiag%run_hdiag(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1)
      end if
   end if
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Copy HDIAG into C matrix
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Copy HDIAG into C matrix'
   call bump%mpl%flush
   call bump%cmat%from_hdiag(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%hdiag)

   ! Release memory
   call bump%hdiag%dealloc
end if

if (bump%nam%new_lct) then
   ! Run LCT driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run LCT driver'
   call bump%mpl%flush
   call bump%lct%run_lct(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)

   ! Copy LCT into C matrix
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Copy LCT into C matrix'
   call bump%mpl%flush
   call bump%cmat%from_lct(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%lct)

   ! Release memory (partial)
   call bump%lct%partial_dealloc
end if

if (bump%nam%load_cmat) then
   ! Read C matrix
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Read C matrix'
   call bump%mpl%flush
   call bump%cmat%read(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%io)
else
   if (bump%nam%forced_radii) then
      ! Copy namelist support radii into C matrix
      write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
      call bump%mpl%flush
      write(bump%mpl%info,'(a)') '--- Copy namelist support radii into C matrix'
      call bump%mpl%flush
      call bump%cmat%from_nam(bump%mpl,bump%nam,bump%geom,bump%bpar)
   end if
end if

if (bump%cmat%allocated.or.bump%nam%new_nicas) then
   ! Get C matrix from BUMP interface
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Get C matrix from BUMP interface'
   call bump%mpl%flush
   call bump%cmat%from_bump(bump%mpl,bump%nam,bump%geom,bump%bpar)

   ! Setup C matrix sampling
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Setup C matrix sampling'
   call bump%mpl%flush
   call bump%cmat%setup_sampling(bump%mpl,bump%nam,bump%geom,bump%bpar)

   if (bump%nam%write_cmat) then
      ! Write C matrix
      write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
      call bump%mpl%flush
      write(bump%mpl%info,'(a)') '--- Write C matrix'
      call bump%mpl%flush
      call bump%cmat%write(bump%mpl,bump%nam,bump%geom,bump%bpar,bump%io)
   end if
end if

if (bump%nam%new_nicas) then
   ! Run NICAS driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run NICAS driver'
   call bump%mpl%flush
   call bump%nicas%run_nicas(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%cmat)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
elseif (bump%nam%load_nicas) then
   ! Read NICAS parameters
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Read NICAS parameters'
   call bump%mpl%flush
   call bump%nicas%read(bump%mpl,bump%nam,bump%geom,bump%bpar)
end if

! Release memory (partial)
call bump%cmat%partial_dealloc

if (bump%nam%new_nicas.or.bump%nam%load_nicas) then
   ! Run NICAS tests driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run NICAS tests driver'
   call bump%mpl%flush
   call bump%nicas%run_nicas_tests(bump%mpl,bump%rng,bump%nam,bump%geom,bump%bpar,bump%io,bump%ens1)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
end if

if (bump%nam%new_obsop) then
   ! Run observation operator driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run observation operator driver'
   call bump%mpl%flush
   call bump%obsop%run_obsop(bump%mpl,bump%rng,bump%nam,bump%geom)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
elseif (bump%nam%load_obsop) then
   ! Read observation operator
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Read observation operator'
   call bump%mpl%flush
   call bump%obsop%read(bump%mpl,bump%nam,bump%geom)
end if

if (bump%nam%new_obsop.or.bump%nam%load_obsop) then
   ! Run observation operator tests driver
   write(bump%mpl%info,'(a)') '-------------------------------------------------------------------'
   call bump%mpl%flush
   write(bump%mpl%info,'(a)') '--- Run observation operator tests driver'
   call bump%mpl%flush
   call bump%obsop%run_obsop_tests(bump%mpl,bump%nam,bump%rng,bump%geom)
   if (bump%nam%default_seed) call bump%rng%reseed(bump%mpl)
end if

end subroutine bump_run_drivers

!----------------------------------------------------------------------
! Subroutine: bump_add_member
! Purpose: add member into bump%ens[1,2]
!----------------------------------------------------------------------
subroutine bump_add_member(bump,afieldset,ie,iens)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset
integer,intent(in) :: ie                        ! Member index
integer,intent(in) :: iens                      ! Ensemble number

! Local variables
integer :: its,iv,nnonzero,nzero,nmask,nnonzero_tot,nzero_tot,nmask_tot
real(kind_real) :: norm,norm_tot,fld_c0a(bump%geom%nc0a,bump%geom%nl0)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
character(len=1024) :: filename,variables
character(len=1024),parameter :: subr = 'bump_add_member'

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

! Check ensemble number
if ((iens/=1).and.(iens/=2)) call bump%mpl%abort(subr,'wrong ensemble number')

! Allocate member
if (iens==1) then
   if (.not.allocated(bump%ens1%mem(ie)%fld)) &
 & allocate(bump%ens1%mem(ie)%fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts))
elseif (iens==2) then
   if (.not.allocated(bump%ens2%mem(ie)%fld)) &
 & allocate(bump%ens2%mem(ie)%fld(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts))
end if

! Add member
write(bump%mpl%info,'(a7,a,i3,a,i1)') '','Member ',ie,' added to ensemble ',iens
call bump%mpl%flush
do its=1,bump%nam%nts
   do iv=1,bump%nam%nv
      ! Model grid to subset Sc0
      call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a)

      ! Copy to ensemble structure
      if (iens==1) then
         bump%ens1%mem(ie)%fld(:,:,iv,its) = fld_c0a
      elseif (iens==2) then
         bump%ens2%mem(ie)%fld(:,:,iv,its) = fld_c0a
      end if

      if (print_member) then
         ! Print norm
         norm = sum(fld_c0a**2,mask=bump%geom%gmask_c0a)
         call bump%mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
         write(bump%mpl%info,'(a10,a,i2,a,i2,a,e9.2)') '','Local norm for variable ',iv,' and timeslot ',its,': ',norm
         call bump%mpl%flush
         write(bump%mpl%info,'(a10,a,i2,a,i2,a,e9.2)') '','Global norm for variable ',iv,' and timeslot ',its,': ',norm_tot
         call bump%mpl%flush
         if (bump%geom%nc0a>0) then
            nnonzero = count((abs(fld_c0a)>0.0).and.bump%geom%gmask_c0a)
            nzero = count((.not.(abs(fld_c0a)>0.0)).and.bump%geom%gmask_c0a)
            nmask = count(.not.bump%geom%gmask_c0a)
         else
            nnonzero = 0
            nzero = 0
            nmask = bump%geom%nc0a
         end if
         call bump%mpl%f_comm%allreduce(nnonzero,nnonzero_tot,fckit_mpi_sum())
         call bump%mpl%f_comm%allreduce(nzero,nzero_tot,fckit_mpi_sum())
         call bump%mpl%f_comm%allreduce(nmask,nmask_tot,fckit_mpi_sum())
         write(bump%mpl%info,'(a10,a,i8,a,i8,a,i8,a,i8)') '','Local total / non-zero / zero / masked points: ', &
 & bump%geom%nc0a*bump%geom%nl0,' / ',nnonzero,' / ',nzero,' / ',nmask
         call bump%mpl%flush
         write(bump%mpl%info,'(a10,a,i8,a,i8,a,i8,a,i8)') '','Global total / non-zero / zero / masked points: ', &
 & bump%geom%nc0*bump%geom%nl0,' / ',nnonzero_tot,' / ',nzero_tot,' / ',nmask_tot
         call bump%mpl%flush
      end if

      if (write_member) then
         ! Write member
         write(filename,'(a,a,i6.6,a,i6.6)') trim(bump%nam%prefix),'_member_',iens,'-',ie
         variables = trim(bump%nam%variables(iv))//'_'//trim(bump%nam%timeslots(its))
         call bump%io%fld_write(bump%mpl,bump%nam,bump%geom,filename,variables,fld_c0a)
      end if
   end do
end do

end subroutine bump_add_member

!----------------------------------------------------------------------
! Subroutine: bump_remove_member
! Purpose: remove member into bump%ens[1,2]
!----------------------------------------------------------------------
subroutine bump_remove_member(bump,afieldset,ie,iens)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset
integer,intent(in) :: ie                        ! Member index
integer,intent(in) :: iens                      ! Ensemble number

! Local variables
integer :: its,iv,isub
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
character(len=1024),parameter :: subr = 'bump_remove_member'

! Check ensemble number
if ((iens/=1).and.(iens/=2)) call bump%mpl%abort(subr,'wrong ensemble number')

! Remove member
write(bump%mpl%info,'(a7,a,i3,a,i1)') '','Member ',ie,' removed from ensemble ',iens
call bump%mpl%flush
do its=1,bump%nam%nts
   do iv=1,bump%nam%nv
      ! Copy from ensemble structure and add mean
      if (iens==1) then
         isub = (ie-1)*bump%ens1%nsub/bump%ens1%ne+1
         fld_c0a = bump%ens1%mem(ie)%fld(:,:,iv,its)+bump%ens1%mean(isub)%fld(:,:,iv,its)
      elseif (iens==2) then
         isub = (ie-1)*bump%ens2%nsub/bump%ens2%ne+1
         fld_c0a = bump%ens2%mem(ie)%fld(:,:,iv,its)+bump%ens2%mean(isub)%fld(:,:,iv,its)
      end if

      ! Model grid to subset Sc0
      call bump%geom%copy_mga_to_c0a(bump%mpl,fld_c0a,fld_mga(:,:,iv,its))
   end do
end do

! Release memory
if (iens==1) then
   deallocate(bump%ens1%mem(ie)%fld)
elseif (iens==2) then
   deallocate(bump%ens2%mem(ie)%fld)
end if

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_remove_member

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal
! Purpose: vertical balance application
!----------------------------------------------------------------------
subroutine bump_apply_vbal(bump,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

do its=1,bump%nam%nts
   if (bump%geom%same_grid) then
       ! Apply vertical balance
      call bump%vbal%apply(bump%nam,bump%geom,bump%bpar,fld_mga(:,:,:,its))
   else
      ! Model grid to subset Sc0
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv))
      end do

      ! Apply vertical balance
      call bump%vbal%apply(bump%nam,bump%geom,bump%bpar,fld_c0a)

      ! Subset Sc0 to model grid
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv),fld_mga(:,:,iv,its))
      end do
   end if
end do

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_apply_vbal

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_inv
! Purpose: vertical balance application, inverse
!----------------------------------------------------------------------
subroutine bump_apply_vbal_inv(bump,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

do its=1,bump%nam%nts
   if (bump%geom%same_grid) then
      ! Apply vertical balance, inverse
      call bump%vbal%apply_inv(bump%nam,bump%geom,bump%bpar,fld_mga(:,:,:,its))
   else
      ! Model grid to subset Sc0
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv))
      end do

      ! Apply vertical balance, inverse
      call bump%vbal%apply_inv(bump%nam,bump%geom,bump%bpar,fld_c0a)

      ! Subset Sc0 to model grid
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv),fld_mga(:,:,iv,its))
      end do
   end if
end do

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_apply_vbal_inv

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_ad
! Purpose: vertical balance application, adjoint
!----------------------------------------------------------------------
subroutine bump_apply_vbal_ad(bump,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

do its=1,bump%nam%nts
   if (bump%geom%same_grid) then
      ! Apply vertical balance, adjoint
      call bump%vbal%apply_ad(bump%nam,bump%geom,bump%bpar,fld_mga(:,:,:,its))
   else
      ! Model grid to subset Sc0
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv))
      end do

      ! Apply vertical balance, adjoint
      call bump%vbal%apply_ad(bump%nam,bump%geom,bump%bpar,fld_c0a)

      ! Subset Sc0 to model grid
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv),fld_mga(:,:,iv,its))
      end do
   end if
end do

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_apply_vbal_ad

!----------------------------------------------------------------------
! Subroutine: bump_apply_vbal_inv_ad
! Purpose: vertical balance application, inverse adjoint
!----------------------------------------------------------------------
subroutine bump_apply_vbal_inv_ad(bump,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

do its=1,bump%nam%nts
   if (bump%geom%same_grid) then
      ! Apply vertical balance, inverse adjoint
      call bump%vbal%apply_inv_ad(bump%nam,bump%geom,bump%bpar,fld_mga(:,:,:,its))
   else
      ! Model grid to subset Sc0
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv))
      end do

      ! Apply vertical balance, inverse adjoint
      call bump%vbal%apply_inv_ad(bump%nam,bump%geom,bump%bpar,fld_c0a)

      ! Subset Sc0 to model grid
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv),fld_mga(:,:,iv,its))
      end do
   end if
end do

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_apply_vbal_inv_ad

!----------------------------------------------------------------------
! Subroutine: bump_apply_stddev
! Purpose: standard-deviation application
!----------------------------------------------------------------------
subroutine bump_apply_stddev(bump,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

 if (bump%geom%same_grid) then
   ! Apply standard-deviation
   call bump%var%apply_sqrt(bump%nam,bump%geom,fld_mga)
else
   ! Model grid to subset Sc0
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv,its))
      end do
   end do

   ! Apply standard-deviation
   call bump%var%apply_sqrt(bump%nam,bump%geom,fld_c0a)

   ! Subset Sc0 to model grid
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv,its),fld_mga(:,:,iv,its))
      end do
   end do
end if

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_apply_stddev

!----------------------------------------------------------------------
! Subroutine: bump_apply_stddev_inv
! Purpose: standard-deviation application, inverse
!----------------------------------------------------------------------
subroutine bump_apply_stddev_inv(bump,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

 if (bump%geom%same_grid) then
   ! Apply standard-deviation inverse
   call bump%var%apply_sqrt_inv(bump%nam,bump%geom,fld_mga)
else
   ! Model grid to subset Sc0
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv,its))
      end do
   end do

   ! Apply standard-deviation
   call bump%var%apply_sqrt_inv(bump%nam,bump%geom,fld_c0a)

   ! Subset Sc0 to model grid
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv,its),fld_mga(:,:,iv,its))
      end do
   end do
end if

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_apply_stddev_inv

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas
! Purpose: NICAS application
!----------------------------------------------------------------------
subroutine bump_apply_nicas(bump,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

if (bump%geom%same_grid) then
   ! Apply NICAS
   if (bump%nam%lsqrt) then
      call bump%nicas%apply_from_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_mga)
   else
      call bump%nicas%apply(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_mga)
   end if
else
   ! Model grid to subset Sc0
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv,its))
      end do
   end do

   ! Apply NICAS
   if (bump%nam%lsqrt) then
      call bump%nicas%apply_from_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_c0a)
   else
      call bump%nicas%apply(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_c0a)
   end if

   ! Subset Sc0 to model grid
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv,its),fld_mga(:,:,iv,its))
      end do
   end do
end if

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_apply_nicas

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_deprecated
! Purpose: NICAS application (deprecated)
!----------------------------------------------------------------------
subroutine bump_apply_nicas_deprecated(bump,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                          ! BUMP
real(kind_real),intent(inout) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
character(len=1024),parameter :: subr = 'bump_apply_nicas_deprecated'

! Deprecation warning
call bump%mpl%warning(subr,'this interface is deprecated, consider using the ATLAS-based interface')

if (bump%geom%same_grid) then
   ! Apply NICAS
   if (bump%nam%lsqrt) then
      call bump%nicas%apply_from_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_mga)
   else
      call bump%nicas%apply(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_mga)
   end if
else
   ! Model grid to subset Sc0
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv,its))
      end do
   end do

   ! Apply NICAS
   if (bump%nam%lsqrt) then
      call bump%nicas%apply_from_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_c0a)
   else
      call bump%nicas%apply(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_c0a)
   end if

   ! Subset Sc0 to model grid
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv,its),fld_mga(:,:,iv,its))
      end do
   end do
end if

end subroutine bump_apply_nicas_deprecated

!----------------------------------------------------------------------
! Subroutine: bump_get_cv_size
! Purpose: get control variable size
!----------------------------------------------------------------------
subroutine bump_get_cv_size(bump,n)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP
integer,intent(out) :: n               ! Control variable size

! Local variables
type(cv_type) :: cv

! Allocate control variable
call bump%nicas%alloc_cv(bump%mpl,bump%bpar,cv,getsizeonly=.true.)

! Copy size
n = cv%n

end subroutine bump_get_cv_size

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_sqrt
! Purpose: NICAS square-root application
!----------------------------------------------------------------------
subroutine bump_apply_nicas_sqrt(bump,pcv,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
real(kind_real),intent(in) :: pcv(:)            ! Packed control variable
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
character(len=1024),parameter :: subr = 'bump_apply_nicas_sqrt'
type(cv_type) :: cv

! Allocation
call bump%nicas%alloc_cv(bump%mpl,bump%bpar,cv)

! Check dimension
if (size(pcv)==cv%n) then
   ! Unpack control variable
   call cv%unpack(pcv)
else
   call bump%mpl%abort(subr,'wrong control variable size in bump_apply_nicas_sqrt')
end if

if (bump%geom%same_grid) then
   ! Apply NICAS square-root
   call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld_mga)
else
   ! Apply NICAS square-root
   call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld_c0a)

   ! Subset Sc0 to model grid
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv,its),fld_mga(:,:,iv,its))
      end do
   end do
end if

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_apply_nicas_sqrt

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_sqrt_deprecated
! Purpose: NICAS square-root application (deprecated)
!----------------------------------------------------------------------
subroutine bump_apply_nicas_sqrt_deprecated(bump,pcv,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                          ! BUMP
real(kind_real),intent(in) :: pcv(:)                                                            ! Packed control variable
real(kind_real),intent(inout) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
character(len=1024),parameter :: subr = 'bump_apply_nicas_sqrt_deprecated'
type(cv_type) :: cv

! Deprecation warning
call bump%mpl%warning(subr,'this interface is deprecated, consider using the ATLAS-based interface')

! Allocation
call bump%nicas%alloc_cv(bump%mpl,bump%bpar,cv)

! Check dimension
if (size(pcv)==cv%n) then
   ! Unpack control variable
   call cv%unpack(pcv)
else
   call bump%mpl%abort(subr,'wrong control variable size in bump_apply_nicas_sqrt')
end if

if (bump%geom%same_grid) then
   ! Apply NICAS square-root
   call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld_mga)
else
   ! Apply NICAS square-root
   call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld_c0a)

   ! Subset Sc0 to model grid
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv,its),fld_mga(:,:,iv,its))
      end do
   end do
end if

end subroutine bump_apply_nicas_sqrt_deprecated

!----------------------------------------------------------------------
! Subroutine: bump_apply_nicas_sqrt_ad
! Purpose: NICAS square-root adjoint application
!----------------------------------------------------------------------
subroutine bump_apply_nicas_sqrt_ad(bump,afieldset,pcv)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset
real(kind_real),intent(inout) :: pcv(:)         ! Packed control variable

! Local variables
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
character(len=1024),parameter :: subr = 'bump_apply_nicas_sqrt_ad'
type(cv_type) :: cv

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

if (bump%geom%same_grid) then
   ! Apply NICAS square-root adjoint
   call bump%nicas%apply_sqrt_ad(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_mga,cv)
else
   ! Model grid to subset Sc0
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,iv,its),fld_c0a(:,:,iv,its))
      end do
   end do

   ! Apply NICAS square-root adjoint
   call bump%nicas%apply_sqrt_ad(bump%mpl,bump%nam,bump%geom,bump%bpar,fld_c0a,cv)
end if

! Check dimension
if (size(pcv)==cv%n) then
   ! Pack control variable
   call cv%pack(pcv)
else
   call bump%mpl%abort(subr,'wrong control variable size in bump_apply_nicas_sqrt_ad')
end if

end subroutine bump_apply_nicas_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: bump_randomize
! Purpose: NICAS randomization
!----------------------------------------------------------------------
subroutine bump_randomize(bump,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variable
integer :: its,iv
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)
type(cv_type) :: cv

! Generate random control vector
call bump%nicas%random_cv(bump%mpl,bump%rng,bump%bpar,cv)

if (bump%geom%same_grid) then
   ! Apply NICAS square-root
   call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld_mga)
else
   ! Apply NICAS square-root
   call bump%nicas%apply_sqrt(bump%mpl,bump%nam,bump%geom,bump%bpar,cv,fld_c0a)

   ! Subset Sc0 to model grid
   do its=1,bump%nam%nts
      do iv=1,bump%nam%nv
         call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a(:,:,iv,its),fld_mga(:,:,iv,its))
      end do
   end do
end if

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_randomize

!----------------------------------------------------------------------
! Subroutine: bump_apply_obsop
! Purpose: observation operator application
!----------------------------------------------------------------------
subroutine bump_apply_obsop(bump,afieldset,obs)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                             ! BUMP
type(atlas_fieldset),intent(inout) :: afieldset                    ! ATLAS fieldset
real(kind_real),intent(out) :: obs(bump%obsop%nobsa,bump%geom%nl0) ! Observations columns

! Local variables
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,1,1)
character(len=1024),parameter :: subr = 'bump_apply_obsop'

! Test dimensions
if (bump%nam%nv>1) call bump%mpl%abort(subr,'only one variable to call bump_apply_obsop')
if (bump%nam%nts>1) call bump%mpl%abort(subr,'only one timeslot to call bump_apply_obsop')

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

if (bump%geom%same_grid) then
   ! Apply observation operator
   call bump%obsop%apply(bump%mpl,bump%geom,fld_mga(:,:,1,1),obs)
else
   ! Model grid to subset Sc0
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga(:,:,1,1),fld_c0a)

   ! Apply observation operator
   call bump%obsop%apply(bump%mpl,bump%geom,fld_c0a,obs)
end if

end subroutine bump_apply_obsop

!----------------------------------------------------------------------
! Subroutine: bump_apply_obsop_deprecated
! Purpose: observation operator application (deprecated)
!----------------------------------------------------------------------
subroutine bump_apply_obsop_deprecated(bump,fld_mga,obs)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                              ! BUMP
real(kind_real),intent(in) :: fld_mga(bump%geom%nmga,bump%geom%nl0) ! Field
real(kind_real),intent(out) :: obs(bump%obsop%nobsa,bump%geom%nl0)  ! Observations columns

! Local variables
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0)
character(len=1024),parameter :: subr = 'bump_apply_obsop_deprecated'

! Deprecation warning
call bump%mpl%warning(subr,'this interface is deprecated, consider using the ATLAS-based interface')

! Test dimensions
if (bump%nam%nv>1) call bump%mpl%abort(subr,'only one variable to call bump_apply_obsop')
if (bump%nam%nts>1) call bump%mpl%abort(subr,'only one timeslot to call bump_apply_obsop')

if (bump%geom%same_grid) then
   ! Apply observation operator
   call bump%obsop%apply(bump%mpl,bump%geom,fld_mga,obs)
else
   ! Model grid to subset Sc0
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,fld_c0a)

   ! Apply observation operator
   call bump%obsop%apply(bump%mpl,bump%geom,fld_c0a,obs)
end if

end subroutine bump_apply_obsop_deprecated

!----------------------------------------------------------------------
! Subroutine: bump_apply_obsop_ad
! Purpose: observation operator adjoint application
!----------------------------------------------------------------------
subroutine bump_apply_obsop_ad(bump,obs,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                            ! BUMP
real(kind_real),intent(in) :: obs(bump%obsop%nobsa,bump%geom%nl0) ! Observations columns
type(atlas_fieldset),intent(inout) :: afieldset                   ! ATLAS fieldset

! Local variables
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0)
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,1,1)
character(len=1024),parameter :: subr = 'bump_apply_obsop_ad'

! Test dimensions
if (bump%nam%nv>1) call bump%mpl%abort(subr,'only one variable to call bump_apply_obsop_ad')
if (bump%nam%nts>1) call bump%mpl%abort(subr,'only one timeslot to call bump_apply_obsop_ad')

if (bump%geom%same_grid) then
   ! Apply observation operator adjoint
   call bump%obsop%apply_ad(bump%mpl,bump%geom,obs,fld_mga(:,:,1,1))
else
   ! Apply observation operator adjoint
   call bump%obsop%apply_ad(bump%mpl,bump%geom,obs,fld_c0a)

   ! Subset Sc0 to model grid
   call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a,fld_mga(:,:,1,1))
end if

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_apply_obsop_ad

!----------------------------------------------------------------------
! Subroutine: bump_apply_obsop_ad_deprecated
! Purpose: observation operator adjoint application (deprecated)
!----------------------------------------------------------------------
subroutine bump_apply_obsop_ad_deprecated(bump,obs,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                               ! BUMP
real(kind_real),intent(in) :: obs(bump%obsop%nobsa,bump%geom%nl0)    ! Observations columns
real(kind_real),intent(out) :: fld_mga(bump%geom%nmga,bump%geom%nl0) ! Field

! Local variables
real(kind_real) :: fld_c0a(bump%geom%nc0a,bump%geom%nl0)
character(len=1024),parameter :: subr = 'bump_apply_obsop_ad_deprecated'

! Deprecation warning
call bump%mpl%warning(subr,'this interface is deprecated, consider using the ATLAS-based interface')

! Test dimensions
if (bump%nam%nv>1) call bump%mpl%abort(subr,'only one variable to call bump_apply_obsop_ad')
if (bump%nam%nts>1) call bump%mpl%abort(subr,'only one timeslot to call bump_apply_obsop_ad')

if (bump%geom%same_grid) then
   ! Apply observation operator adjoint
   call bump%obsop%apply_ad(bump%mpl,bump%geom,obs,fld_mga)
else
   ! Apply observation operator adjoint
   call bump%obsop%apply_ad(bump%mpl,bump%geom,obs,fld_c0a)

   ! Subset Sc0 to model grid
   call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a,fld_mga)
end if

end subroutine bump_apply_obsop_ad_deprecated

!----------------------------------------------------------------------
! Subroutine: bump_get_parameter
! Purpose: get a parameter
!----------------------------------------------------------------------
subroutine bump_get_parameter(bump,param,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
character(len=*),intent(in) :: param            ! Parameter
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variables
integer :: ib,iv,jv,its,jts
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)

write(bump%mpl%info,'(a7,a,a)') '','Get ',trim(param)
call bump%mpl%flush

select case (trim(param))
case ('cor_rh','cor_rv','loc_coef','loc_rh','loc_rv','hyb_coef','loc_D11','loc_D22','loc_D33', &
 & 'loc_D12','loc_Dcoef','loc_DLh')
   select case (trim(bump%nam%strategy))
   case ('specific_univariate','specific_multivariate')
      do ib=1,bump%bpar%nb
         ! Get indices
         iv = bump%bpar%b_to_v1(ib)
         jv = bump%bpar%b_to_v2(ib)
         its = bump%bpar%b_to_ts1(ib)
         jts = bump%bpar%b_to_ts2(ib)

         ! Copy to field
         if ((iv==jv).and.(its==jts)) call bump%copy_to_field(param,ib,fld_mga(:,:,iv,its))
      end do
   case ('common','common_univariate','common_weighted')
      ! Set common index
      ib = bump%bpar%nbe

      do its=1,bump%nam%nts
         do iv=1,bump%nam%nv
            ! Copy to field
            call bump%copy_to_field(param,ib,fld_mga(:,:,iv,its))
         end do
      end do
   end select
case default
   do ib=1,bump%bpar%nb
      ! Get indices
      iv = bump%bpar%b_to_v1(ib)
      jv = bump%bpar%b_to_v2(ib)
      its = bump%bpar%b_to_ts1(ib)
      jts = bump%bpar%b_to_ts2(ib)

      ! Copy to field
      if ((iv==jv).and.(its==jts)) call bump%copy_to_field(param,ib,fld_mga(:,:,iv,its))
   end do
end select

! Field to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

end subroutine bump_get_parameter

!----------------------------------------------------------------------
! Subroutine: bump_copy_to_field
! Purpose: copy to field
!----------------------------------------------------------------------
subroutine bump_copy_to_field(bump,param,ib,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                               ! BUMP
character(len=*),intent(in) :: param                                 ! Parameter
integer,intent(in) :: ib                                             ! Block index
real(kind_real),intent(out) :: fld_mga(bump%geom%nmga,bump%geom%nl0) ! Field

! Local variables
integer :: iscales,ie,imga,il0,iv,its
real(kind_real) :: tmp
character(len=1024),parameter :: subr = 'bump_copy_to_field'

! Check allocation / parameter existence
select case (trim(param))
case ('stddev')
   if (.not.allocated(bump%var%m2sqrt)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
case ('cor_rh','cor_rv','loc_coef','loc_rh','loc_rv','hyb_coef','loc_D11','loc_D22','loc_D33','loc_D12', &
 & 'loc_Dcoef','loc_DLh')
   if (.not.allocated(bump%cmat%blk)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
case default
   select case (param(1:4))
   case ('D11_','D22_','D33_','D12_','Dcoe','DLh_')
      if (.not.allocated(bump%lct%blk)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   case default
      if (param(1:6)=='ens1u_') then
         if (.not.allocated(bump%ens1u%mem)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
      else
         call bump%mpl%abort(subr,'parameter '//trim(param)//' not yet implemented in get_parameter')
      end if
   end select
end select

! Select parameter from cmat
select case (trim(param))
case ('stddev')
   iv = bump%bpar%b_to_v1(ib)
   its = bump%bpar%b_to_ts1(ib)
   if (.not.allocated(bump%var%m2sqrt)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%var%m2sqrt(:,:,iv,its),fld_mga)
case ('cor_rh')
   if (.not.allocated(bump%cmat%blk(ib)%rh)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rh,fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnot(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req
      end do
   end do
case ('cor_rv')
   if (.not.allocated(bump%cmat%blk(ib)%rv)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv,fld_mga)
case ('loc_coef')
   if (.not.allocated(bump%cmat%blk(ib)%coef_ens)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%coef_ens,fld_mga)
case ('loc_rh')
   if (.not.allocated(bump%cmat%blk(ib)%rh)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rh,fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnot(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req
      end do
   end do
case ('loc_rv')
   if (.not.allocated(bump%cmat%blk(ib)%rv)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv,fld_mga)
case ('hyb_coef')
   if (.not.allocated(bump%cmat%blk(ib)%coef_sta)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%coef_sta,fld_mga)
case ('loc_D11','loc_D22')
   if (.not.allocated(bump%cmat%blk(ib)%rh)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rh,fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnot(fld_mga(imga,il0))) then
            tmp = fld_mga(imga,il0)*req
            call lct_r2d(tmp,fld_mga(imga,il0))
         end if
      end do
   end do
case ('loc_D33')
   if (.not.allocated(bump%cmat%blk(ib)%rv)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rv,fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnot(fld_mga(imga,il0))) then
            tmp = fld_mga(imga,il0)
            call lct_r2d(tmp,fld_mga(imga,il0))
         end if
      end do
   end do
case ('loc_D12')
   if (.not.allocated(bump%cmat%blk(ib)%rh)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   fld_mga = 0.0
case ('loc_Dcoef')
   if (.not.allocated(bump%cmat%blk(ib)%coef_ens)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%coef_ens,fld_mga)
case ('loc_DLh')
   if (.not.allocated(bump%cmat%blk(ib)%rh)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%cmat%blk(ib)%rh,fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnot(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req
      end do
   end do
end select

! Select parameter from lct
select case (param(1:min(len(param),4)))
case ('D11_')
   if (.not.allocated(bump%lct%blk(ib)%D11)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   read(param(5:5),'(i1)') iscales
   if (iscales>size(bump%lct%blk(ib)%D11,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D11(:,:,iscales),fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnot(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req**2
      end do
   end do
case ('D22_')
   if (.not.allocated(bump%lct%blk(ib)%D22)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   read(param(5:5),'(i1)') iscales
   if (iscales>size(bump%lct%blk(ib)%D22,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D22(:,:,iscales),fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnot(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req**2
      end do
   end do
case ('D33_')
   if (.not.allocated(bump%lct%blk(ib)%D33)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   read(param(5:5),'(i1)') iscales
   if (iscales>size(bump%lct%blk(ib)%D33,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D33(:,:,iscales),fld_mga)
case ('D12_')
   if (.not.allocated(bump%lct%blk(ib)%D12)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   read(param(5:5),'(i1)') iscales
   if (iscales>size(bump%lct%blk(ib)%D12,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%D12(:,:,iscales),fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnot(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req**2
      end do
   end do
case ('Dcoe')
   if (.not.allocated(bump%lct%blk(ib)%Dcoef)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   read(param(7:7),'(i1)') iscales
   if (iscales>size(bump%lct%blk(ib)%Dcoef,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%Dcoef(:,:,iscales),fld_mga)
case ('DLh_')
   if (.not.allocated(bump%lct%blk(ib)%DLh)) call bump%mpl%abort(subr,trim(param)//' is not allocated in bump%copy_to_field')
   read(param(5:5),'(i1)') iscales
   if (iscales>size(bump%lct%blk(ib)%DLh,3)) call bump%mpl%abort(subr,trim(param)//' has fewer scales in bump%copy_to_field')
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%lct%blk(ib)%DLh(:,:,iscales),fld_mga)
   do il0=1,bump%geom%nl0
      do imga=1,bump%geom%nmga
         if (bump%mpl%msv%isnot(fld_mga(imga,il0))) fld_mga(imga,il0) = fld_mga(imga,il0)*req
      end do
   end do
end select

! Select parameter from ens1u
if (param(1:min(6,len(param)))=='ens1u_') then
   read(param(7:12),'(i6.6)') ie
   if (ie>size(bump%ens1u%mem)) call bump%mpl%abort(subr,trim(param)//' has fewer members in bump%copy_to_field')
   iv = bump%bpar%b_to_v1(ib)
   its = bump%bpar%b_to_ts1(ib)
   call bump%geom%copy_c0a_to_mga(bump%mpl,bump%ens1u%mem(ie)%fld(:,:,iv,its),fld_mga)
end if

end subroutine bump_copy_to_field

!----------------------------------------------------------------------
! Subroutine: bump_test_get_parameter
! Purpose: test get_parameter
!----------------------------------------------------------------------
subroutine bump_test_get_parameter(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP

! Local variables
type(atlas_fieldset) :: afieldset

! Create ATLAS fieldset with empty fields
call create_atlas_fieldset(bump%geom%afunctionspace_mg,bump%geom%nl0,bump%nam%variables(1:bump%nam%nv), &
 & bump%nam%timeslots(1:bump%nam%nts),afieldset)

! Get parameter
if (bump%nam%check_get_param_stddev) then
   call bump%get_parameter('stddev',afieldset)
elseif (bump%nam%check_get_param_cor) then
   call bump%get_parameter('cor_rh',afieldset)
   call bump%get_parameter('cor_rv',afieldset)
elseif (bump%nam%check_get_param_hyb) then
   call bump%get_parameter('loc_coef',afieldset)
   call bump%get_parameter('loc_rh',afieldset)
   call bump%get_parameter('loc_rv',afieldset)
   call bump%get_parameter('hyb_coef',afieldset)
elseif (bump%nam%check_get_param_Dloc) then
   call bump%get_parameter('loc_D11',afieldset)
   call bump%get_parameter('loc_D22',afieldset)
   call bump%get_parameter('loc_D33',afieldset)
   call bump%get_parameter('loc_D12',afieldset)
   call bump%get_parameter('loc_Dcoef',afieldset)
   call bump%get_parameter('loc_DLh',afieldset)
elseif (bump%nam%check_get_param_lct) then
   call bump%get_parameter('D11_1',afieldset)
   call bump%get_parameter('D22_1',afieldset)
   call bump%get_parameter('D33_1',afieldset)
   call bump%get_parameter('D12_1',afieldset)
   call bump%get_parameter('Dcoef_1',afieldset)
   call bump%get_parameter('DLh_1',afieldset)
   call bump%get_parameter('D11_2',afieldset)
   call bump%get_parameter('D22_2',afieldset)
   call bump%get_parameter('D33_2',afieldset)
   call bump%get_parameter('D12_2',afieldset)
   call bump%get_parameter('Dcoef_2',afieldset)
   call bump%get_parameter('DLh_2',afieldset)
end if

! Release memory
call afieldset%final()

end subroutine bump_test_get_parameter

!----------------------------------------------------------------------
! Subroutine: bump_set_parameter
! Purpose: set a parameter
!----------------------------------------------------------------------
subroutine bump_set_parameter(bump,param,afieldset)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump          ! BUMP
character(len=*),intent(in) :: param            ! Parameter
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset

! Local variables
integer :: ib,iv,jv,its,jts
real(kind_real) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts)

! ATLAS fieldset to field
call atlas_to_fld(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),afieldset,fld_mga,bump%nam%lev2d)

write(bump%mpl%info,'(a7,a,a)') '','Set ',trim(param)
call bump%mpl%flush

select case (trim(param))
case ('cor_rh','cor_rv','loc_coef','loc_rh','loc_rv','hyb_coef')
   select case (trim(bump%nam%strategy))
   case ('specific_univariate','specific_multivariate')
      do ib=1,bump%bpar%nb
         ! Get indices
         iv = bump%bpar%b_to_v1(ib)
         jv = bump%bpar%b_to_v2(ib)
         its = bump%bpar%b_to_ts1(ib)
         jts = bump%bpar%b_to_ts2(ib)

         ! Copy to field
         if ((iv==jv).and.(its==jts)) call bump%copy_from_field(param,ib,fld_mga(:,:,iv,its))
      end do
   case ('common','common_univariate','common_weighted')
      ! Set common index
      ib = bump%bpar%nbe

      do its=1,bump%nam%nts
         do iv=1,bump%nam%nv
            ! Copy to field
            call bump%copy_from_field(param,ib,fld_mga(:,:,iv,its))
         end do
      end do
   end select
case default
   do ib=1,bump%bpar%nb
      ! Get indices
      iv = bump%bpar%b_to_v1(ib)
      jv = bump%bpar%b_to_v2(ib)
      its = bump%bpar%b_to_ts1(ib)
      jts = bump%bpar%b_to_ts2(ib)

      ! Copy to field
      if ((iv==jv).and.(its==jts)) call bump%copy_from_field(param,ib,fld_mga(:,:,iv,its))
   end do
end select

end subroutine bump_set_parameter

!----------------------------------------------------------------------
! Subroutine: bump_set_parameter_deprecated
! Purpose: set a parameter (deprecated)
!----------------------------------------------------------------------
subroutine bump_set_parameter_deprecated(bump,param,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                                                       ! BUMP
character(len=*),intent(in) :: param                                                         ! Parameter
real(kind_real),intent(in) :: fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts) ! Field

! Local variables
integer :: ib,iv,jv,its,jts
character(len=1024),parameter :: subr = 'bump_set_parameter_deprecated'

! Deprecation warning
call bump%mpl%warning(subr,'this interface is deprecated, consider using the ATLAS-based interface')

write(bump%mpl%info,'(a7,a,a)') '','Set ',trim(param)
call bump%mpl%flush

select case (trim(param))
case ('cor_rh','cor_rv','loc_coef','loc_rh','loc_rv','hyb_coef')
   select case (trim(bump%nam%strategy))
   case ('specific_univariate','specific_multivariate')
      do ib=1,bump%bpar%nb
         ! Get indices
         iv = bump%bpar%b_to_v1(ib)
         jv = bump%bpar%b_to_v2(ib)
         its = bump%bpar%b_to_ts1(ib)
         jts = bump%bpar%b_to_ts2(ib)

         ! Copy to field
         if ((iv==jv).and.(its==jts)) call bump%copy_from_field(param,ib,fld_mga(:,:,iv,its))
      end do
   case ('common','common_univariate','common_weighted')
      ! Set common index
      ib = bump%bpar%nbe

      do its=1,bump%nam%nts
         do iv=1,bump%nam%nv
            ! Copy to field
            call bump%copy_from_field(param,ib,fld_mga(:,:,iv,its))
         end do
      end do
   end select
case default
   do ib=1,bump%bpar%nb
      ! Get indices
      iv = bump%bpar%b_to_v1(ib)
      jv = bump%bpar%b_to_v2(ib)
      its = bump%bpar%b_to_ts1(ib)
      jts = bump%bpar%b_to_ts2(ib)

      ! Copy to field
      if ((iv==jv).and.(its==jts)) call bump%copy_from_field(param,ib,fld_mga(:,:,iv,its))
   end do
end select

end subroutine bump_set_parameter_deprecated

!----------------------------------------------------------------------
! Subroutine: bump_copy_from_field
! Purpose: copy from field
!----------------------------------------------------------------------
subroutine bump_copy_from_field(bump,param,ib,fld_mga)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump                              ! BUMP
character(len=*),intent(in) :: param                                ! Parameter
integer,intent(in) :: ib                                            ! Block index
real(kind_real),intent(in) :: fld_mga(bump%geom%nmga,bump%geom%nl0) ! Field

! Local variables
integer :: ic0a,il0
character(len=1024),parameter :: subr = 'bump_copy_from_field'

! Check allocation / parameter existence
select case (trim(param))
case ('cor_rh','cor_rv','loc_coef','loc_rh','loc_rv','hyb_coef','D11','D22','D33','D12','Dcoef')
   if (.not.allocated(bump%cmat%blk)) allocate(bump%cmat%blk(bump%bpar%nbe))
case default
   call bump%mpl%abort(subr,'parameter '//trim(param)//' not yet implemented in set_parameter')
end select

! Select parameter from cmat
select case (trim(param))
case ('cor_rh')
   if (.not.allocated(bump%cmat%blk(ib)%bump_rh)) allocate(bump%cmat%blk(ib)%bump_rh(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_rh)
   do il0=1,bump%geom%nl0
      do ic0a=1,bump%geom%nc0a
         if (bump%mpl%msv%isnot(bump%cmat%blk(ib)%bump_rh(ic0a,il0))) &
 & bump%cmat%blk(ib)%bump_rh(ic0a,il0) = bump%cmat%blk(ib)%bump_rh(ic0a,il0)/req
      end do
   end do
case ('cor_rv')
   if (.not.allocated(bump%cmat%blk(ib)%bump_rv)) allocate(bump%cmat%blk(ib)%bump_rv(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_rv)
case ('loc_coef')
   if (.not.allocated(bump%cmat%blk(ib)%bump_coef_ens)) allocate(bump%cmat%blk(ib)%bump_coef_ens(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_coef_ens)
case ('loc_rh')
   if (.not.allocated(bump%cmat%blk(ib)%bump_rh)) allocate(bump%cmat%blk(ib)%bump_rh(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_rh)
   do il0=1,bump%geom%nl0
      do ic0a=1,bump%geom%nc0a
         if (bump%mpl%msv%isnot(bump%cmat%blk(ib)%bump_rh(ic0a,il0))) &
 & bump%cmat%blk(ib)%bump_rh(ic0a,il0) = bump%cmat%blk(ib)%bump_rh(ic0a,il0)/req
      end do
   end do
case ('loc_rv')
   if (.not.allocated(bump%cmat%blk(ib)%bump_rv)) allocate(bump%cmat%blk(ib)%bump_rv(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_rv)
case ('hyb_coef')
   if (.not.allocated(bump%cmat%blk(ib)%bump_coef_sta)) allocate(bump%cmat%blk(ib)%bump_coef_sta(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_coef_sta)
case ('D11')
   if (.not.allocated(bump%cmat%blk(ib)%bump_D11)) allocate(bump%cmat%blk(ib)%bump_D11(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_D11)
   do il0=1,bump%geom%nl0
      do ic0a=1,bump%geom%nc0a
         if (bump%mpl%msv%isnot(bump%cmat%blk(ib)%bump_D11(ic0a,il0))) &
 & bump%cmat%blk(ib)%bump_D11(ic0a,il0) = bump%cmat%blk(ib)%bump_D11(ic0a,il0)/req**2
      end do
   end do
case ('D22')
   if (.not.allocated(bump%cmat%blk(ib)%bump_D22)) allocate(bump%cmat%blk(ib)%bump_D22(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_D22)
   do il0=1,bump%geom%nl0
      do ic0a=1,bump%geom%nc0a
         if (bump%mpl%msv%isnot(bump%cmat%blk(ib)%bump_D22(ic0a,il0))) &
 & bump%cmat%blk(ib)%bump_D22(ic0a,il0) = bump%cmat%blk(ib)%bump_D22(ic0a,il0)/req**2
      end do
   end do
case ('D33')
   if (.not.allocated(bump%cmat%blk(ib)%bump_D33)) allocate(bump%cmat%blk(ib)%bump_D33(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_D33)
case ('D12')
   if (.not.allocated(bump%cmat%blk(ib)%bump_D12)) allocate(bump%cmat%blk(ib)%bump_D12(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_D12)
   do il0=1,bump%geom%nl0
      do ic0a=1,bump%geom%nc0a
         if (bump%mpl%msv%isnot(bump%cmat%blk(ib)%bump_D12(ic0a,il0))) &
 & bump%cmat%blk(ib)%bump_D12(ic0a,il0) = bump%cmat%blk(ib)%bump_D12(ic0a,il0)/req**2
      end do
   end do
case ('Dcoef')
   if (.not.allocated(bump%cmat%blk(ib)%bump_Dcoef)) allocate(bump%cmat%blk(ib)%bump_Dcoef(bump%geom%nc0a,bump%geom%nl0))
   call bump%geom%copy_mga_to_c0a(bump%mpl,fld_mga,bump%cmat%blk(ib)%bump_Dcoef)
case default
   call bump%mpl%abort(subr,'parameter '//trim(param)//' not yet implemented in set_parameter')
end select

end subroutine bump_copy_from_field

!----------------------------------------------------------------------
! Subroutine: bump_test_set_parameter
! Purpose: test set_parameter
!----------------------------------------------------------------------
subroutine bump_test_set_parameter(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP

! Local variables
integer :: iv,its,ic0a,ic0
real(kind_real) :: hash_min,hash_max,hash_spread_inv
real(kind_real),allocatable :: fld_c0a(:,:),fld_mga(:,:,:,:)
type(atlas_fieldset) :: afieldset,afieldset_req,afieldset_reqsq,afieldset_vert,afieldset_vertsq

! Allocation
allocate(fld_c0a(bump%geom%nc0a,bump%geom%nl0))
allocate(fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts))

! Initialization
call bump%mpl%f_comm%allreduce(minval(bump%geom%hash_c0a),hash_min,fckit_mpi_min())
call bump%mpl%f_comm%allreduce(maxval(bump%geom%hash_c0a),hash_max,fckit_mpi_max())
hash_spread_inv = 1.0/(hash_max-hash_min)
do its=1,bump%nam%nts
   do iv=1,bump%nam%nv
      do ic0a=1,bump%geom%nc0a
         ic0 = bump%geom%c0a_to_c0(ic0a)
         fld_c0a(ic0a,:) = max(min(1.0e-6_kind_real,(bump%geom%hash_c0a(ic0a)-hash_min)*hash_spread_inv),1.0-1.0e-6_kind_real)
      end do
      call bump%geom%copy_c0a_to_mga(bump%mpl,fld_c0a,fld_mga(:,:,iv,its))
   end do
end do

! Create ATLAS fieldset with empty fields
call create_atlas_fieldset(bump%geom%afunctionspace_mg,bump%geom%nl0,bump%nam%variables(1:bump%nam%nv), &
 & bump%nam%timeslots(1:bump%nam%nts),afieldset)
call create_atlas_fieldset(bump%geom%afunctionspace_mg,bump%geom%nl0,bump%nam%variables(1:bump%nam%nv), &
 & bump%nam%timeslots(1:bump%nam%nts),afieldset_req)
call create_atlas_fieldset(bump%geom%afunctionspace_mg,bump%geom%nl0,bump%nam%variables(1:bump%nam%nv), &
 & bump%nam%timeslots(1:bump%nam%nts),afieldset_reqsq)
call create_atlas_fieldset(bump%geom%afunctionspace_mg,bump%geom%nl0,bump%nam%variables(1:bump%nam%nv), &
 & bump%nam%timeslots(1:bump%nam%nts),afieldset_vert)
call create_atlas_fieldset(bump%geom%afunctionspace_mg,bump%geom%nl0,bump%nam%variables(1:bump%nam%nv), &
 & bump%nam%timeslots(1:bump%nam%nts),afieldset_vertsq)

! Convert to ATLAS fieldset
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga*req, &
 & afieldset_req,bump%nam%lev2d)
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga*req**2, &
 & afieldset_reqsq,bump%nam%lev2d)
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts), &
 & (1.0+fld_mga)*(maxval(bump%geom%vunitavg)-minval(bump%geom%vunitavg)),afieldset_vert,bump%nam%lev2d)
call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts), &
 & ((1.0+fld_mga)*(maxval(bump%geom%vunitavg)-minval(bump%geom%vunitavg))),afieldset_vertsq,bump%nam%lev2d)

! Set parameter
if (bump%nam%check_set_param_cor) then
   call bump%set_parameter('cor_rh',afieldset_req)
   call bump%set_parameter('cor_rv',afieldset_vert)
elseif (bump%nam%check_set_param_hyb) then
   call bump%set_parameter('loc_coef',afieldset)
   call bump%set_parameter('loc_rh',afieldset_req)
   call bump%set_parameter('loc_rv',afieldset_vert)
   call bump%set_parameter('hyb_coef',afieldset)
elseif (bump%nam%check_set_param_lct) then
   call bump%set_parameter('D11',afieldset_reqsq)
   call bump%set_parameter('D22',afieldset_reqsq)
   call bump%set_parameter('D33',afieldset_vertsq)
   call bump%set_parameter('D12',afieldset)
   call bump%set_parameter('Dcoef',afieldset)
end if

! Release memory
deallocate(fld_c0a)
deallocate(fld_mga)
call afieldset%final()
call afieldset_req%final()
call afieldset_reqsq%final()

end subroutine bump_test_set_parameter

!----------------------------------------------------------------------
! Subroutine: bump_test_apply_interfaces
! Purpose: test BUMP apply interfaces
!----------------------------------------------------------------------
subroutine bump_test_apply_interfaces(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP

! Local variables
integer :: n,nv_save,nts_save
real(kind_real),allocatable :: fld_mga(:,:,:,:),pcv(:),obs(:,:)
type(atlas_fieldset) :: afieldset

! Test apply_vbal
if (bump%nam%check_apply_vbal) then
   write(bump%mpl%info,'(a7,a)') '','Test apply_vbal'
   call bump%mpl%flush

   ! Allocation
   allocate(fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts))

   ! Initialization
   call bump%rng%rand_real(0.0_kind_real,1.0_kind_real,fld_mga)

   ! Create ATLAS fieldset with empty fields
   call create_atlas_fieldset(bump%geom%afunctionspace_mg,bump%geom%nl0,bump%nam%variables(1:bump%nam%nv), &
 & bump%nam%timeslots(1:bump%nam%nts),afieldset)

   ! Convert to ATLAS fieldset
   call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

   ! Calls
   call bump%apply_vbal(afieldset)
   call bump%apply_vbal_inv(afieldset)
   call bump%apply_vbal_ad(afieldset)
   call bump%apply_vbal_inv_ad(afieldset)

   ! Release memory
   deallocate(fld_mga)
   call afieldset%final()
end if

! Test apply_stddev
if (bump%nam%check_apply_stddev) then
   write(bump%mpl%info,'(a7,a)') '','Test apply_stddev'
   call bump%mpl%flush

   ! Allocation
   allocate(fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts))

   ! Initialization
   call bump%rng%rand_real(0.0_kind_real,1.0_kind_real,fld_mga)

   ! Create ATLAS fieldset with empty fields
   call create_atlas_fieldset(bump%geom%afunctionspace_mg,bump%geom%nl0,bump%nam%variables,bump%nam%timeslots,afieldset)

   ! Convert to ATLAS fieldset
   call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

   ! Calls
   call bump%apply_stddev(afieldset)
   call bump%apply_stddev_inv(afieldset)

   ! Release memory
   deallocate(fld_mga)
   call afieldset%final()
end if

! Test apply_nicas
if (bump%nam%check_apply_nicas) then
   write(bump%mpl%info,'(a7,a)') '','Test apply_nicas'
   call bump%mpl%flush

   ! Get control variable size
   call bump%get_cv_size(n)

   ! Allocation
   allocate(fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts))
   allocate(pcv(n))

   ! Initialization
   call bump%rng%rand_real(0.0_kind_real,1.0_kind_real,fld_mga)

   ! Create ATLAS fieldset with empty fields
   call create_atlas_fieldset(bump%geom%afunctionspace_mg,bump%geom%nl0,bump%nam%variables(1:bump%nam%nv), &
 & bump%nam%timeslots(1:bump%nam%nts),afieldset)

   ! Convert to ATLAS fieldset
   call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

   ! Calls
   call bump%apply_nicas(afieldset)
   call bump%apply_nicas_sqrt(pcv,afieldset)
   call bump%apply_nicas_sqrt_ad(afieldset,pcv)
   call bump%randomize(afieldset)

   ! Release memory
   deallocate(fld_mga)
   deallocate(pcv)
   call afieldset%final()
end if

! Test apply_obsop
if (bump%nam%check_apply_obsop) then
   write(bump%mpl%info,'(a7,a)') '','Test apply_obsop'
   call bump%mpl%flush

   ! Save namelist parameters
   nv_save = bump%nam%nv
   nts_save = bump%nam%nts

   ! Set namelist parameters
   bump%nam%nv = 1
   bump%nam%nts = 1

   ! Allocation
   allocate(fld_mga(bump%geom%nmga,bump%geom%nl0,bump%nam%nv,bump%nam%nts))
   allocate(obs(bump%obsop%nobsa,bump%geom%nl0))

   ! Initialization
   call bump%rng%rand_real(0.0_kind_real,1.0_kind_real,fld_mga)

   ! Create ATLAS fieldset with empty fields
   call create_atlas_fieldset(bump%geom%afunctionspace_mg,bump%geom%nl0,bump%nam%variables(1:bump%nam%nv), &
 & bump%nam%timeslots(1:bump%nam%nts),afieldset)

   ! Convert to ATLAS fieldset
   call fld_to_atlas(bump%mpl,bump%nam%variables(1:bump%nam%nv),bump%nam%timeslots(1:bump%nam%nts),fld_mga,afieldset,bump%nam%lev2d)

   ! Calls
   call bump%apply_obsop(afieldset,obs)
   call bump%apply_obsop_ad(obs,afieldset)

   ! Reset namelist parameters
   bump%nam%nv = nv_save
   bump%nam%nts = nts_save

   ! Release memory
   deallocate(fld_mga)
   deallocate(obs)
   call afieldset%final()
end if

end subroutine bump_test_apply_interfaces

!----------------------------------------------------------------------
! Subroutine: bump_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine bump_partial_dealloc(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP

! Release memory
call bump%cmat%partial_dealloc
call bump%ens1%dealloc
call bump%ens1u%dealloc
call bump%ens2%dealloc
call bump%geom%partial_dealloc
call bump%hdiag%dealloc
call bump%io%dealloc
call bump%lct%partial_dealloc
call bump%nicas%partial_dealloc
call bump%obsop%partial_dealloc
call bump%var%partial_dealloc
call bump%vbal%partial_dealloc
if (allocated(bump%fld_uv)) deallocate(bump%fld_uv)

end subroutine bump_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: bump_dealloc
! Purpose: release memory (full)
!----------------------------------------------------------------------
subroutine bump_dealloc(bump)

implicit none

! Passed variables
class(bump_type),intent(inout) :: bump ! BUMP

! Release memory
call bump%bpar%dealloc
call bump%cmat%dealloc
call bump%ens1%dealloc
call bump%ens1u%dealloc
call bump%ens2%dealloc
call bump%geom%dealloc
call bump%hdiag%dealloc
call bump%io%dealloc
call bump%lct%dealloc
call bump%nicas%dealloc
call bump%obsop%dealloc
call bump%var%dealloc
call bump%vbal%dealloc
if (allocated(bump%fld_uv)) deallocate(bump%fld_uv)

end subroutine bump_dealloc

!----------------------------------------------------------------------
! Subroutine: dummy
! Purpose: dummy finalization
!----------------------------------------------------------------------
subroutine dummy(bump)

implicit none

! Passed variables
type(bump_type),intent(inout) :: bump ! BUMP

end subroutine dummy

end module type_bump
