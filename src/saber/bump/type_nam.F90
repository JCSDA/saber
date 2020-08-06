!----------------------------------------------------------------------
! Module: type_nam
! Purpose: namelist derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_nam

use fckit_configuration_module, only: fckit_configuration,fckit_yamlconfiguration
use fckit_pathname_module, only : fckit_pathname
use iso_c_binding
use tools_const, only: pi,req,deg2rad,rad2deg
use tools_kinds,only: kind_real
use type_mpl, only: mpl_type

implicit none

integer,parameter :: nvmax = 20      ! Maximum number of variables
integer,parameter :: ntsmax = 99     ! Maximum number of time slots
integer,parameter :: nlmax = 200     ! Maximum number of levels
integer,parameter :: nc3max = 1000   ! Maximum number of classes
integer,parameter :: nscalesmax = 5  ! Maximum number of variables
integer,parameter :: ndirmax = 300   ! Maximum number of diracs
integer,parameter :: nldwvmax = 99   ! Maximum number of local diagnostic profiles
integer,parameter :: nprociomax = 20 ! Maximum number of I/O tasks

type nam_type
   ! general_param
   character(len=1024) :: datadir                       ! Data directory
   character(len=1024) :: prefix                        ! Files prefix
   character(len=1024) :: model                         ! Model name ('aro', 'arp', 'fv3', 'gem', 'geos', 'gfs', 'ifs', 'mpas', 'nemo', 'qg, 'res' or 'wrf')
   character(len=1024) :: verbosity                     ! Verbosity level ('all', 'main' or 'none')
   logical :: colorlog                                  ! Add colors to the log (for display on terminal)
   logical :: default_seed                              ! Default seed for random numbers
   logical :: repro                                     ! Inter-compilers reproducibility
   integer :: nprocio                                   ! Number of IO processors
   logical :: remap                                     ! Remap points to improve load balance

   ! driver_param
   character(len=1024) :: method                        ! Localization/hybridization to compute ('cor', 'loc', 'hyb-avg', 'hyb-rnd' or 'dual-ens')
   character(len=1024) :: strategy                      ! Localization strategy ('diag_all', 'common', 'common_univariate', 'common_weighted', 'specific_univariate' or 'specific_multivariate')
   logical :: new_normality                             ! New normality test
   logical :: new_cortrack                              ! New correlation tracker
   logical :: new_corstats                              ! New correlation statistics
   logical :: new_vbal                                  ! Compute new vertical balance operator
   logical :: load_vbal                                 ! Load existing vertical balance operator
   logical :: write_vbal                                ! Write vertical balance operator
   logical :: new_var                                   ! Compute new variance
   logical :: load_var                                  ! Load existing variance
   logical :: write_var                                 ! Write variance
   logical :: new_mom                                   ! Compute new sample moments
   logical :: load_mom                                  ! Load sample moments
   logical :: write_mom                                 ! Write sample moments
   logical :: new_hdiag                                 ! Compute new HDIAG diagnostics
   logical :: write_hdiag                               ! Write HDIAG diagnostics
   logical :: new_lct                                   ! Compute new LCT
   logical :: write_lct                                 ! Write LCT
   logical :: load_cmat                                 ! Load existing C matrix
   logical :: write_cmat                                ! Write existing C matrix
   logical :: new_nicas                                 ! Compute new NICAS parameters
   logical :: load_nicas                                ! Load existing NICAS parameters
   logical :: write_nicas                               ! Write NICAS parameters
   logical :: new_obsop                                 ! Compute new observation operator
   logical :: load_obsop                                ! Load existing observation operator
   logical :: write_obsop                               ! Write observation operator
   logical :: check_vbal                                ! Test vertical balance inverse and adjoint
   logical :: check_adjoints                            ! Test NICAS adjoints
   logical :: check_dirac                               ! Test NICAS application on diracs
   logical :: check_randomization                       ! Test NICAS randomization
   logical :: check_consistency                         ! Test HDIAG-NICAS consistency
   logical :: check_optimality                          ! Test HDIAG optimality
   logical :: check_obsop                               ! Test observation operator
   logical :: check_no_obs                              ! Test observation operator with no observation on the last MPI task
   logical :: check_no_point                            ! Test BUMP with no grid point on the last MPI task
   logical :: check_no_point_mask                       ! Test BUMP with all grid points masked on the last MPI task
   logical :: check_no_point_nicas                      ! Test NICAS with no subgrid point on the last MPI task
   logical :: check_set_param_cor                       ! Test set_parameter interface for correlation
   logical :: check_set_param_hyb                       ! Test set_parameter interface for hybrid case
   logical :: check_set_param_lct                       ! Test set_parameter interface for LCT
   logical :: check_get_param_stddev                    ! Test get_parameter interface for standard-deviation
   logical :: check_get_param_cor                       ! Test get_parameter interface for correlation
   logical :: check_get_param_hyb                       ! Test get_parameter interface for hybrid case
   logical :: check_get_param_Dloc                      ! Test get_parameter interface for anisotropic localization
   logical :: check_get_param_lct                       ! Test get_parameter interface for LCT
   logical :: check_apply_vbal                          ! Test apply_vbal interfaces
   logical :: check_apply_stddev                        ! Test apply_stddev interfaces
   logical :: check_apply_nicas                         ! Test apply_nicas interfaces
   logical :: check_apply_obsop                         ! Test apply_obsop interfaces

   ! model_param
   integer :: nl                                        ! Number of levels
   integer :: levs(nlmax)                               ! Levels
   character(len=1024) :: lev2d                         ! Level for 2D variables ('first' or 'last')
   logical :: logpres                                   ! Use pressure logarithm as vertical coordinate (model level if .false.)
   integer :: nv                                        ! Number of variables
   character(len=1024),dimension(nvmax) :: variables    ! Variables names
   integer :: nts                                       ! Number of time slots
   character(len=1024),dimension(ntsmax) :: timeslots   ! Timeslots
   real(kind_real) :: dts                               ! Timeslots width [in s]
   logical :: nomask                                    ! Do not use geometry mask
   character(len=1024) :: wind_filename                 ! Wind field file name
   character(len=1024) :: wind_variables(2)               ! Wind field variables names (u and v)

   ! ens1_param
   integer :: ens1_ne                                   ! Ensemble 1 size
   integer :: ens1_nsub                                 ! Ensemble 1 sub-ensembles number

   ! ens2_param
   integer :: ens2_ne                                   ! Ensemble 2 size
   integer :: ens2_nsub                                 ! Ensemble 2 sub-ensembles number

   ! sampling_param
   logical :: sam_write                                 ! Write sampling
   logical :: sam_write_grids                           ! Write sampling grids
   logical :: sam_read                                  ! Read sampling
   character(len=1024) :: mask_type                     ! Mask restriction type
   character(len=1024),dimension(nvmax) :: mask_lu      ! Mask threshold side ("lower" if mask_th is the lower bound, "upper" if mask_th is the upper bound)
   real(kind_real),dimension(nvmax) :: mask_th          ! Mask threshold
   integer :: ncontig_th                                ! Threshold on vertically contiguous points for sampling mask (0 to skip the test)
   logical :: mask_check                                ! Check that sampling couples and interpolations do not cross mask boundaries
   character(len=1024) :: draw_type                     ! Sampling draw type ('random_uniform','random_coast' or 'icosahedron')
   real(kind_real) :: Lcoast                            ! Length-scale to increase sampling density along coasts [in meters]
   real(kind_real) :: rcoast                            ! Minimum value to increase sampling density along coasts
   integer :: nc1                                       ! Number of sampling points
   integer :: nc2                                       ! Number of diagnostic points
   integer :: ntry                                      ! Number of tries to get the most separated point for the zero-separation sampling
   integer :: nrep                                      ! Number of replacement to improve homogeneity of the zero-separation sampling
   integer :: nc3                                       ! Number of classes
   real(kind_real) :: dc                                ! Class size (for sam_type='hor'), should be larger than the typical grid cell size [in meters]
   integer :: nl0r                                      ! Reduced number of levels for diagnostics
   integer :: irmax                                     ! Maximum number of random number draws

   ! diag_param
   integer :: ne                                        ! Ensemble size
   real(kind_real) :: gen_kurt_th                       ! Threshold on generalized kurtosis (3.0 = Gaussian distribution)
   logical :: gau_approx                                ! Gaussian approximation for asymptotic quantities
   integer :: avg_nbins                                 ! Number of bins for averaged statistics histograms
   logical :: vbal_block(nvmax*(nvmax-1)/2)             ! Activation of vertical balance (ordered line by line in the lower triangular formulation)
   real(kind_real) :: vbal_rad                          ! Vertical balance diagnostic radius [in meters]
   real(kind_real) :: vbal_dlat                         ! Vertical balance diagnostic latitude band half-width [in degrees]
   logical :: vbal_diag_auto(nvmax*(nvmax-1)/2)         ! Diagonal auto-covariance for the inversion
   logical :: vbal_diag_reg(nvmax*(nvmax-1)/2)          ! Diagonal regression
   logical :: var_filter                                ! Filter variances
   integer :: var_niter                                 ! Number of iteration for the variances filtering
   real(kind_real) :: var_rhflt                         ! Variances initial filtering support radius [in meters]
   logical :: local_diag                                ! Activate local diagnostics
   real(kind_real) :: local_rad                         ! Local diagnostics calculation radius [in meters]
   logical :: adv_diag                                  ! Activate advection diagnostic
   character(len=1024) :: adv_type                      ! Advection diagnostic type ('max', 'wind' or 'windmax')
   real(kind_real) :: adv_rad                           ! Advection diagnostic calculation radius [in meters]
   integer :: adv_niter                                 ! Number of iteration for the advection filtering
   real(kind_real) :: adv_rhflt                         ! Advection initial filtering support radius [in meters]
   real(kind_real) :: adv_valid                         ! Required proportion of valid points for filtering convergence

   ! fit_param
   character(len=1024) :: minim_algo                    ! Minimization algorithm ('none', 'fast' or 'hooke')
   real(kind_real) :: diag_rhflt                        ! Horizontal filtering suport radius [in meters]
   real(kind_real) :: diag_rvflt                        ! Vertical filtering support radius
   real(kind_real) :: smoothness_penalty                ! Smoothness penalty weight (default 0.01)
   integer :: fit_dl0                                   ! Number of levels between interpolation levels
   integer :: lct_nscales                               ! Number of LCT scales
   real(kind_real) :: lct_scale_ratio                   ! Factor between diffusion scales
   real(kind_real) :: lct_cor_min                       ! Minimum relevant correlation for LCT first guess
   logical :: lct_diag(nscalesmax)                      ! Diagnostic of diagonal LCT components only
   real(kind_real) :: lct_qc_th                         ! LCT quality control threshold
   real(kind_real) :: lct_qc_max                        ! LCT quality control maximum
   logical :: lct_write_cor                             ! Write full correlations

   ! nicas_param
   logical :: nonunit_diag                              ! Non-unit diagonal for the NICAS application
   logical :: lsqrt                                     ! Square-root formulation
   real(kind_real) :: resol                             ! Resolution
   integer :: nc1max                                    ! Maximum size of the Sc1 subset
   logical :: fast_sampling                             ! Fast sampling flag
   character(len=1024) :: subsamp                       ! Subsampling structure ('h', 'hv', 'vh' or 'hvh')
   logical :: network                                   ! Network-base convolution calculation (distance-based if false)
   integer :: mpicom                                    ! Number of communication steps
   integer :: adv_mode                                  ! Advection mode (1: direct, -1: direct and inverse)
   logical :: forced_radii                              ! Force specific support radii
   real(kind_real) :: rh                                ! Forced horizontal support radius [in meters]
   real(kind_real) :: rv                                ! Forced vertical support radius
   logical :: pos_def_test                              ! Positive-definiteness test
   logical :: write_grids                               ! Write NICAS grids

   ! dirac_param
   integer :: ndir                                      ! Number of Diracs
   real(kind_real) :: londir(ndirmax)                   ! Diracs longitudes [in degrees]
   real(kind_real) :: latdir(ndirmax)                   ! Diracs latitudes [in degrees]
   integer :: levdir(ndirmax)                           ! Diracs level
   integer :: ivdir(ndirmax)                            ! Diracs variable indices
   integer :: itsdir(ndirmax)                           ! Diracs timeslots indices

   ! obsop_param
   integer :: nobs                                      ! Number of observations

   ! output_param
   integer :: nldwv                                     ! Number of local diagnostics profiles to write (for local_diag = .true.)
   integer :: img_ldwv(nldwvmax)                        ! Index on model grid of the local diagnostics profiles to write
   real(kind_real) :: lon_ldwv(nldwvmax)                ! Longitudes of the local diagnostics profiles to write [in degrees]
   real(kind_real) :: lat_ldwv(nldwvmax)                ! Latitudes of the local diagnostics profiles to write [in degrees]
   character(len=1024),dimension(nldwvmax) :: name_ldwv ! Name of the local diagnostics profiles to write
   logical :: grid_output                               ! Write regridded fields
   real(kind_real) :: grid_resol                        ! Regridded fields resolution [in meters]
contains
   procedure :: init => nam_init
   procedure :: read => nam_read
   procedure :: read_yaml => nam_read_yaml
   procedure :: bcast => nam_bcast
   procedure :: from_conf => nam_from_conf
   procedure :: check => nam_check
   procedure :: write => nam_write
end type nam_type

private
public :: nvmax,ntsmax,nlmax,nc3max,nscalesmax,ndirmax,nldwvmax
public :: nam_type

contains

!----------------------------------------------------------------------
! Subroutine: nam_init
! Purpose: intialize
!----------------------------------------------------------------------
subroutine nam_init(nam,nproc)

implicit none

! Passed variable
class(nam_type),intent(out) :: nam ! Namelist
integer,intent(in) :: nproc        ! Number of MPI task

! Local variable
integer :: il,iv,ildwv

! general_param default
nam%datadir = '.'
nam%prefix = ''
nam%model = 'online'
nam%verbosity = 'all'
nam%colorlog = .false.
nam%default_seed = .true.
nam%repro = .true.
nam%nprocio = min(nproc,nprociomax)
nam%remap = .false.

! driver_param default
nam%method = ''
nam%strategy = ''
nam%new_normality = .false.
nam%new_cortrack = .false.
nam%new_corstats = .false.
nam%new_vbal = .false.
nam%load_vbal = .false.
nam%write_vbal = .true.
nam%new_var = .false.
nam%load_var = .false.
nam%write_var = .true.
nam%new_mom = .true.
nam%load_mom = .false.
nam%write_mom = .false.
nam%new_hdiag = .false.
nam%write_hdiag = .true.
nam%new_lct = .false.
nam%write_lct = .true.
nam%load_cmat = .false.
nam%write_cmat = .true.
nam%new_nicas = .false.
nam%load_nicas = .false.
nam%write_nicas = .true.
nam%new_obsop = .false.
nam%load_obsop = .false.
nam%write_obsop = .true.
nam%check_vbal = .false.
nam%check_adjoints = .false.
nam%check_dirac = .false.
nam%check_randomization = .false.
nam%check_consistency = .false.
nam%check_optimality = .false.
nam%check_obsop = .false.
nam%check_no_obs = .false.
nam%check_no_point = .false.
nam%check_no_point_mask = .false.
nam%check_no_point_nicas = .false.
nam%check_set_param_cor = .false.
nam%check_set_param_hyb = .false.
nam%check_set_param_lct = .false.
nam%check_get_param_stddev = .false.
nam%check_get_param_cor = .false.
nam%check_get_param_hyb = .false.
nam%check_get_param_Dloc = .false.
nam%check_get_param_lct = .false.
nam%check_apply_vbal = .false.
nam%check_apply_stddev = .false.
nam%check_apply_nicas = .false.
nam%check_apply_obsop = .false.

! model_param default
nam%nl = 0
do il=1,nlmax
   nam%levs(il) = il
end do
nam%lev2d = 'first'
nam%logpres = .false.
nam%nv = 0
do iv=1,nvmax
   nam%variables(iv) = ''
end do
nam%nts = 0
nam%timeslots = ''
nam%dts = 3600.0
nam%nomask = .false.
nam%wind_filename = ''
nam%wind_variables = (/'',''/)

! ens1_param default
nam%ens1_ne = 0
nam%ens1_nsub = 1

! ens2_param default
nam%ens2_ne = 0
nam%ens2_nsub = 1

! sampling_param default
nam%sam_write = .false.
nam%sam_write_grids = .false.
nam%sam_read = .false.
nam%mask_type = 'none'
do iv=1,nvmax
   nam%mask_lu(iv) = 'lower'
   nam%mask_th(iv) = 0.0
end do
nam%ncontig_th = 0
nam%mask_check = .false.
nam%draw_type = 'random_uniform'
nam%Lcoast = 0.0
nam%rcoast = 0.0
nam%nc1 = 0
nam%nc2 = 0
nam%ntry = 0
nam%nrep = 0
nam%nc3 = 0
nam%dc = 0.0
nam%nl0r = 0
nam%irmax = 10000

! diag_param default
nam%ne = 0
nam%gen_kurt_th = huge(1.0)
nam%gau_approx = .false.
nam%avg_nbins = 0
do iv=1,nvmax*(nvmax-1)/2
   nam%vbal_block(iv) = .false.
end do
nam%vbal_rad = 0.0
nam%vbal_dlat = 0.0
do iv=1,nvmax*(nvmax-1)/2
   nam%vbal_diag_auto(iv) = .false.
end do
do iv=1,nvmax*(nvmax-1)/2
   nam%vbal_diag_reg(iv) = .false.
end do
nam%var_filter = .false.
nam%var_niter = 0
nam%var_rhflt = 0.0
nam%local_diag = .false.
nam%local_rad = 0.0
nam%adv_diag = .false.
nam%adv_type = ''
nam%adv_rad = 0.0
nam%adv_niter = 0
nam%adv_rhflt = 0.0
nam%adv_valid = 0.99

! fit_param default
nam%minim_algo = 'hooke'
nam%diag_rhflt = 0.0
nam%diag_rvflt = 0.0
nam%smoothness_penalty = 0.01
nam%fit_dl0 = 1
nam%lct_nscales = 0
nam%lct_scale_ratio = 10.0
nam%lct_cor_min = 0.5
nam%lct_diag = .false.
nam%lct_qc_th = 0.0
nam%lct_qc_max = 1.0
nam%lct_write_cor = .false.

! nicas_param default
nam%nonunit_diag = .false.
nam%lsqrt = .true.
nam%resol = 0.0
nam%nc1max = 15000
nam%fast_sampling = .false.
nam%subsamp = 'hv'
nam%network = .false.
nam%mpicom = 0
nam%adv_mode = 0
nam%forced_radii = .false.
nam%rh = 0.0
nam%rv = 0.0
nam%pos_def_test = .false.
nam%write_grids = .false.

! dirac_param default
nam%ndir = 0
nam%londir = 0.0
nam%latdir = 0.0
nam%levdir = 0
nam%ivdir = 0
nam%itsdir = 0

! obsop_param default
nam%nobs = 0

! output_param default
nam%nldwv = 0
nam%img_ldwv = 0
nam%lon_ldwv = 0.0
nam%lat_ldwv = 0.0
do ildwv=1,nldwvmax
   nam%name_ldwv(ildwv) = ''
end do
nam%grid_output = .false.
nam%grid_resol = 0.0

end subroutine nam_init

!----------------------------------------------------------------------
! Subroutine: nam_read
! Purpose: read
!----------------------------------------------------------------------
subroutine nam_read(nam,mpl,namelname)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam     ! Namelist
type(mpl_type),intent(inout) :: mpl      ! MPI data
character(len=*),intent(in) :: namelname ! Namelist name

! Local variables
integer :: iv
character(len=1024),parameter :: subr = 'nam_read'

! Namelist variables
character(len=1024) :: datadir
character(len=1024) :: prefix
character(len=1024) :: model
character(len=1024) :: verbosity
logical :: colorlog
logical :: default_seed
logical :: repro
integer :: nprocio
logical :: remap
character(len=1024) :: method
character(len=1024) :: strategy
logical :: new_normality
logical :: new_cortrack
logical :: new_corstats
logical :: new_vbal
logical :: load_vbal
logical :: write_vbal
logical :: new_var
logical :: load_var
logical :: write_var
logical :: new_mom
logical :: load_mom
logical :: write_mom
logical :: new_hdiag
logical :: write_hdiag
logical :: new_lct
logical :: write_lct
logical :: load_cmat
logical :: write_cmat
logical :: new_nicas
logical :: load_nicas
logical :: write_nicas
logical :: new_obsop
logical :: load_obsop
logical :: write_obsop
logical :: check_vbal
logical :: check_adjoints
logical :: check_dirac
logical :: check_randomization
logical :: check_consistency
logical :: check_optimality
logical :: check_obsop
logical :: check_no_obs
logical :: check_no_point
logical :: check_no_point_mask
logical :: check_no_point_nicas
logical :: check_set_param_cor
logical :: check_set_param_hyb
logical :: check_set_param_lct
logical :: check_get_param_stddev
logical :: check_get_param_cor
logical :: check_get_param_hyb
logical :: check_get_param_Dloc
logical :: check_get_param_lct
logical :: check_apply_vbal
logical :: check_apply_stddev
logical :: check_apply_nicas
logical :: check_apply_obsop
integer :: nl
integer :: levs(nlmax)
character(len=1024) :: lev2d
logical :: logpres
integer :: nv
character(len=1024),dimension(nvmax) :: variables
integer :: nts
character(len=1024),dimension(ntsmax) :: timeslots
real(kind_real) :: dts
logical :: nomask
character(len=1024) :: wind_filename
character(len=1024) :: wind_variables(2)
integer :: ens1_ne
integer :: ens1_nsub
integer :: ens2_ne
integer :: ens2_nsub
logical :: sam_write
logical :: sam_write_grids
logical :: sam_read
character(len=1024) :: mask_type
character(len=1024),dimension(nvmax) :: mask_lu
real(kind_real),dimension(nvmax) :: mask_th
integer :: ncontig_th
logical :: mask_check
character(len=1024) :: draw_type
real(kind_real) :: Lcoast
real(kind_real) :: rcoast
integer :: nc1
integer :: nc2
integer :: ntry
integer :: nrep
integer :: nc3
real(kind_real) :: dc
integer :: nl0r
integer :: irmax
integer :: ne
real(kind_real) :: gen_kurt_th
logical :: gau_approx
integer :: avg_nbins
logical :: vbal_block(nvmax*(nvmax-1)/2)
real(kind_real) :: vbal_rad
real(kind_real) :: vbal_dlat
logical :: vbal_diag_auto(nvmax*(nvmax-1)/2)
logical :: vbal_diag_reg(nvmax*(nvmax-1)/2)
logical :: var_filter
integer :: var_niter
real(kind_real) :: var_rhflt
logical :: local_diag
real(kind_real) :: local_rad
logical :: adv_diag
character(len=1024) :: adv_type
real(kind_real) :: adv_rad
integer :: adv_niter
real(kind_real) :: adv_rhflt
real(kind_real) :: adv_valid
character(len=1024) :: minim_algo
real(kind_real) :: diag_rhflt
real(kind_real) :: diag_rvflt
real(kind_real) :: smoothness_penalty
integer :: fit_dl0
integer :: lct_nscales
real(kind_real) :: lct_scale_ratio
real(kind_real) :: lct_cor_min
logical :: lct_diag(nscalesmax)
real(kind_real) :: lct_qc_th
real(kind_real) :: lct_qc_max
logical :: lct_write_cor
logical :: nonunit_diag
logical :: lsqrt
real(kind_real) :: resol
integer :: nc1max
logical :: fast_sampling
character(len=1024) :: subsamp
logical :: network
integer :: mpicom
integer :: adv_mode
logical :: forced_radii
real(kind_real) :: rh
real(kind_real) :: rv
logical :: pos_def_test
logical :: write_grids
integer :: ndir
real(kind_real) :: londir(ndirmax)
real(kind_real) :: latdir(ndirmax)
integer :: levdir(ndirmax)
integer :: ivdir(ndirmax)
integer :: itsdir(ndirmax)
integer :: nobs
integer :: nldwv
integer :: img_ldwv(nldwvmax)
real(kind_real) :: lon_ldwv(nldwvmax)
real(kind_real) :: lat_ldwv(nldwvmax)
character(len=1024),dimension(nldwvmax) :: name_ldwv
logical :: grid_output
real(kind_real) :: grid_resol

! Local variables
integer :: il,ildwv,lunit

! Namelist blocks
namelist/general_param/datadir, &
                     & prefix, &
                     & model, &
                     & verbosity, &
                     & colorlog, &
                     & default_seed, &
                     & repro, &
                     & nprocio, &
                     & remap
namelist/driver_param/method, &
                    & strategy, &
                    & new_normality, &
                    & new_cortrack, &
                    & new_corstats, &
                    & new_vbal, &
                    & load_vbal, &
                    & write_vbal, &
                    & new_var, &
                    & load_var, &
                    & write_var, &
                    & new_mom, &
                    & load_mom, &
                    & write_mom, &
                    & new_hdiag, &
                    & write_hdiag, &
                    & new_lct, &
                    & write_lct, &
                    & load_cmat, &
                    & write_cmat, &
                    & new_nicas, &
                    & load_nicas, &
                    & write_nicas, &
                    & new_obsop, &
                    & load_obsop, &
                    & write_obsop, &
                    & check_vbal, &
                    & check_adjoints, &
                    & check_dirac, &
                    & check_randomization, &
                    & check_consistency, &
                    & check_optimality, &
                    & check_obsop, &
                    & check_no_obs, &
                    & check_no_point, &
                    & check_no_point_mask, &
                    & check_no_point_nicas, &
                    & check_set_param_cor, &
                    & check_set_param_hyb, &
                    & check_set_param_lct, &
                    & check_get_param_stddev, &
                    & check_get_param_cor, &
                    & check_get_param_hyb, &
                    & check_get_param_Dloc, &
                    & check_get_param_lct, &
                    & check_apply_vbal, &
                    & check_apply_stddev, &
                    & check_apply_nicas, &
                    & check_apply_obsop
namelist/model_param/nl, &
                   & levs, &
                   & lev2d, &
                   & logpres, &
                   & nv, &
                   & variables, &
                   & nts, &
                   & timeslots, &
                   & dts, &
                   & nomask, &
                   & wind_filename, &
                   & wind_variables
namelist/ens1_param/ens1_ne, &
                  & ens1_nsub
namelist/ens2_param/ens2_ne, &
                  & ens2_nsub
namelist/sampling_param/sam_write, &
                      & sam_write_grids, &
                      & sam_read, &
                      & mask_type, &
                      & mask_lu, &
                      & mask_th, &
                      & ncontig_th, &
                      & mask_check, &
                      & draw_type, &
                      & Lcoast, &
                      & rcoast, &
                      & nc1, &
                      & nc2, &
                      & ntry, &
                      & nrep, &
                      & nc3, &
                      & dc, &
                      & nl0r, &
                      & irmax
namelist/diag_param/ne, &
                  & gen_kurt_th, &
                  & gau_approx, &
                  & avg_nbins, &
                  & vbal_block, &
                  & vbal_rad, &
                  & vbal_dlat, &
                  & vbal_diag_auto, &
                  & vbal_diag_reg, &
                  & var_filter, &
                  & var_niter, &
                  & var_rhflt, &
                  & local_diag, &
                  & local_rad, &
                  & adv_diag, &
                  & adv_type, &
                  & adv_rad, &
                  & adv_niter, &
                  & adv_rhflt, &
                  & adv_valid
namelist/fit_param/minim_algo, &
                 & diag_rhflt, &
                 & diag_rvflt, &
                 & smoothness_penalty, &
                 & fit_dl0, &
                 & lct_nscales, &
                 & lct_scale_ratio, &
                 & lct_cor_min, &
                 & lct_diag, &
                 & lct_qc_th, &
                 & lct_qc_max, &
                 & lct_write_cor
namelist/nicas_param/nonunit_diag, &
                   & lsqrt, &
                   & resol, &
                   & nc1max, &
                   & fast_sampling, &
                   & subsamp, &
                   & network, &
                   & mpicom, &
                   & adv_mode, &
                   & forced_radii, &
                   & rh, &
                   & rv, &
                   & pos_def_test, &
                   & write_grids, &
                   & ndir, &
                   & londir, &
                   & latdir, &
                   & levdir, &
                   & ivdir, &
                   & itsdir
namelist/obsop_param/nobs
namelist/output_param/nldwv, &
                    & img_ldwv, &
                    & lon_ldwv, &
                    & lat_ldwv, &
                    & name_ldwv, &
                    & diag_rhflt, &
                    & grid_output, &
                    & grid_resol

if (mpl%main) then
   ! general_param default
   datadir = '.'
   prefix = ''
   model = 'online'
   verbosity = 'all'
   colorlog = .false.
   default_seed = .true.
   repro = .true.
   nprocio = min(mpl%nproc,nprociomax)
   remap = .false.

   ! driver_param default
   method = ''
   strategy = ''
   new_normality = .false.
   new_cortrack = .false.
   new_corstats = .false.
   new_vbal = .false.
   load_vbal = .false.
   write_vbal = .true.
   new_var = .false.
   load_var = .false.
   write_var = .true.
   new_mom = .true.
   load_mom = .false.
   write_mom = .false.
   new_hdiag = .false.
   write_hdiag = .true.
   new_lct = .false.
   write_lct = .true.
   load_cmat = .false.
   write_cmat = .true.
   new_nicas = .false.
   load_nicas = .false.
   write_nicas = .true.
   new_obsop = .false.
   load_obsop = .false.
   write_obsop = .true.
   check_vbal = .false.
   check_adjoints = .false.
   check_dirac = .false.
   check_randomization = .false.
   check_consistency = .false.
   check_optimality = .false.
   check_obsop = .false.
   check_no_obs = .false.
   check_no_point = .false.
   check_no_point_mask = .false.
   check_no_point_nicas = .false.
   check_set_param_cor = .false.
   check_set_param_hyb = .false.
   check_set_param_lct = .false.
   check_get_param_stddev = .false.
   check_get_param_cor = .false.
   check_get_param_hyb = .false.
   check_get_param_Dloc = .false.
   check_get_param_lct = .false.
   check_apply_vbal = .false.
   check_apply_stddev = .false.
   check_apply_nicas = .false.
   check_apply_obsop = .false.

   ! model_param default
   nl = 0
   do il=1,nlmax
      levs(il) = il
   end do
   lev2d = 'first'
   logpres = .false.
   nv = 0
   do iv=1,nvmax
      variables(iv) = ''
   end do
   nts = 0
   timeslots = ''
   dts = 3600.0
   nomask = .false.
   wind_filename = ''
   wind_variables = (/'',''/)

   ! ens1_param default
   ens1_ne = 0
   ens1_nsub = 1

   ! ens2_param default
   ens2_ne = 0
   ens2_nsub = 1

   ! sampling_param default
   sam_write = .false.
   sam_write_grids = .false.
   sam_read = .false.
   mask_type = 'none'
   do iv=1,nvmax
      mask_lu(iv) = 'lower'
      mask_th(iv) = 0.0
   end do
   ncontig_th = 0
   mask_check = .false.
   draw_type = 'random_uniform'
   Lcoast = 0.0
   rcoast = 0.0
   nc1 = 0
   nc2 = 0
   ntry = 0
   nrep = 0
   nc3 = 0
   dc = 0.0
   nl0r = 0
   irmax = 10000

   ! diag_param default
   ne = 0
   gen_kurt_th = huge(1.0)
   gau_approx = .false.
   avg_nbins = 0
   do iv=1,nvmax*(nvmax-1)/2
      vbal_block(iv) = .false.
   end do
   vbal_rad = 0.0
   vbal_dlat = 0.0
   do iv=1,nvmax*(nvmax-1)/2
      vbal_diag_auto(iv) = .true.
   end do
   do iv=1,nvmax*(nvmax-1)/2
      vbal_diag_reg(iv) = .true.
   end do
   var_filter = .false.
   var_niter = 0
   var_rhflt = 0.0
   local_diag = .false.
   local_rad = 0.0
   adv_diag = .false.
   adv_type = ''
   adv_rad = 0.0
   adv_niter = 0
   adv_rhflt = 0.0
   adv_valid = 0.99

   ! fit_param default
   minim_algo = 'hooke'
   diag_rhflt = 0.0
   diag_rvflt = 0.0
   smoothness_penalty = 0.01
   fit_dl0 = 1
   lct_nscales = 0
   lct_scale_ratio = 10.0
   lct_cor_min = 0.5
   lct_diag = .false.
   lct_qc_th = 0.0
   lct_qc_max = 1.0
   lct_write_cor = .false.

   ! nicas_param default
   nonunit_diag = .false.
   lsqrt = .true.
   resol = 0.0
   nc1max = 15000
   fast_sampling = .false.
   subsamp = 'hv'
   network = .false.
   mpicom = 0
   adv_mode = 0
   forced_radii = .false.
   rh = 0.0
   rv = 0.0
   pos_def_test = .false.
   write_grids = .false.

   ! dirac_param default
   ndir = 0
   londir = 0.0
   latdir = 0.0
   levdir = 0
   ivdir = 0
   itsdir = 0

   ! obsop_param default
   nobs = 0

   ! output_param default
   nldwv = 0
   img_ldwv = 0
   lon_ldwv = 0.0
   lat_ldwv = 0.0
   do ildwv=1,nldwvmax
      name_ldwv(ildwv) = ''
   end do
   grid_output = .false.
   grid_resol = 0.0

   ! Open namelist
   call mpl%newunit(lunit)
   open(unit=lunit,file=trim(namelname),status='old',action='read')

   ! general_param
   read(lunit,nml=general_param)
   nam%datadir = datadir
   nam%prefix = prefix
   nam%model = model
   nam%verbosity = verbosity
   nam%colorlog = colorlog
   nam%default_seed = default_seed
   nam%repro = repro
   nam%nprocio = nprocio
   nam%remap = remap

   ! driver_param
   read(lunit,nml=driver_param)
   nam%method = method
   nam%strategy = strategy
   nam%new_normality = new_normality
   nam%new_cortrack = new_cortrack
   nam%new_corstats = new_corstats
   nam%new_vbal = new_vbal
   nam%load_vbal = load_vbal
   nam%write_vbal = write_vbal
   nam%new_var = new_var
   nam%load_var = load_var
   nam%write_var = write_var
   nam%new_mom = new_mom
   nam%load_mom = load_mom
   nam%write_mom = write_mom
   nam%new_hdiag = new_hdiag
   nam%write_hdiag = write_hdiag
   nam%new_lct = new_lct
   nam%write_lct = write_lct
   nam%load_cmat = load_cmat
   nam%write_cmat = write_cmat
   nam%new_nicas = new_nicas
   nam%load_nicas = load_nicas
   nam%write_nicas = write_nicas
   nam%new_obsop = new_obsop
   nam%load_obsop = load_obsop
   nam%write_obsop = write_obsop
   nam%check_vbal = check_vbal
   nam%check_adjoints = check_adjoints
   nam%check_dirac = check_dirac
   nam%check_randomization = check_randomization
   nam%check_consistency = check_consistency
   nam%check_optimality = check_optimality
   nam%check_obsop = check_obsop
   nam%check_no_obs = check_no_obs
   nam%check_no_point = check_no_point
   nam%check_no_point_mask = check_no_point_mask
   nam%check_no_point_nicas = check_no_point_nicas
   nam%check_set_param_cor = check_set_param_cor
   nam%check_set_param_hyb = check_set_param_hyb
   nam%check_set_param_lct = check_set_param_lct
   nam%check_get_param_stddev = check_get_param_stddev
   nam%check_get_param_cor = check_get_param_cor
   nam%check_get_param_hyb = check_get_param_hyb
   nam%check_get_param_Dloc = check_get_param_Dloc
   nam%check_get_param_lct = check_get_param_lct
   nam%check_apply_vbal = check_apply_vbal
   nam%check_apply_stddev = check_apply_stddev
   nam%check_apply_nicas = check_apply_nicas
   nam%check_apply_obsop = check_apply_obsop

   ! model_param
   read(lunit,nml=model_param)
   if (nl>nlmax) call mpl%abort(subr,'nl is too large')
   if (nv>nvmax) call mpl%abort(subr,'nv is too large')
   if (nts>ntsmax) call mpl%abort(subr,'nts is too large')
   nam%nl = nl
   if (nl>0) nam%levs(1:nl) = levs(1:nl)
   nam%lev2d = lev2d
   nam%logpres = logpres
   nam%nv = nv
   if (nv>0) nam%variables(1:nv) = variables(1:nv)
   nam%nts = nts
   if (nts>0) nam%timeslots(1:nts) = timeslots(1:nts)
   nam%dts = dts
   nam%nomask = nomask
   nam%wind_filename = wind_filename
   nam%wind_variables = wind_variables

   ! ens1_param
   read(lunit,nml=ens1_param)
   nam%ens1_ne = ens1_ne
   nam%ens1_nsub = ens1_nsub

   ! ens2_param
   read(lunit,nml=ens2_param)
   nam%ens2_ne = ens2_ne
   nam%ens2_nsub = ens2_nsub

   ! sampling_param
   read(lunit,nml=sampling_param)
   if (nc3>nc3max) call mpl%abort(subr,'nc3 is too large')
   nam%sam_write = sam_write
   nam%sam_write_grids = sam_write_grids
   nam%sam_read = sam_read
   nam%mask_type = mask_type
   if (nv>0) nam%mask_lu(1:nam%nv) = mask_lu(1:nam%nv)
   if (nv>0) nam%mask_th(1:nam%nv) = mask_th(1:nam%nv)
   nam%ncontig_th = ncontig_th
   nam%mask_check = mask_check
   nam%draw_type = draw_type
   nam%Lcoast = Lcoast
   nam%rcoast = rcoast
   nam%nc1 = nc1
   nam%nc2 = nc2
   nam%ntry = ntry
   nam%nrep = nrep
   nam%nc3 = nc3
   nam%dc = dc
   nam%nl0r = nl0r
   nam%irmax = irmax

   ! diag_param
   read(lunit,nml=diag_param)
   nam%ne = ne
   nam%gen_kurt_th = gen_kurt_th
   nam%gau_approx = gau_approx
   nam%avg_nbins = avg_nbins
   if (nv>1) nam%vbal_block(1:nam%nv*(nam%nv-1)/2) = vbal_block(1:nam%nv*(nam%nv-1)/2)
   nam%vbal_rad = vbal_rad
   nam%vbal_dlat = vbal_dlat
   if (nv>1) nam%vbal_diag_auto(1:nam%nv*(nam%nv-1)/2) = vbal_diag_auto(1:nam%nv*(nam%nv-1)/2)
   if (nv>1) nam%vbal_diag_reg(1:nam%nv*(nam%nv-1)/2) = vbal_diag_reg(1:nam%nv*(nam%nv-1)/2)
   nam%var_filter = var_filter
   nam%var_niter = var_niter
   nam%var_rhflt = var_rhflt
   nam%local_diag = local_diag
   nam%local_rad = local_rad
   nam%adv_diag = adv_diag
   nam%adv_type = adv_type
   nam%adv_rad = adv_rad
   nam%adv_niter = adv_niter
   nam%adv_rhflt = adv_rhflt
   nam%adv_valid = adv_valid

   ! fit_param
   read(lunit,nml=fit_param)
   if (lct_nscales>nscalesmax) call mpl%abort(subr,'lct_nscales is too large')
   nam%minim_algo = minim_algo
   nam%diag_rhflt = diag_rhflt
   nam%diag_rvflt = diag_rvflt
   nam%smoothness_penalty = smoothness_penalty
   nam%fit_dl0 = fit_dl0
   nam%lct_nscales = lct_nscales
   nam%lct_scale_ratio = lct_scale_ratio
   nam%lct_cor_min = lct_cor_min
   if (lct_nscales>0) nam%lct_diag(1:lct_nscales) = lct_diag(1:lct_nscales)
   nam%lct_qc_th = lct_qc_th
   nam%lct_qc_max = lct_qc_max
   nam%lct_write_cor = lct_write_cor

   ! nicas_param
   read(lunit,nml=nicas_param)
   nam%nonunit_diag = nonunit_diag
   nam%lsqrt = lsqrt
   nam%resol = resol
   nam%nc1max = nc1max
   nam%fast_sampling = fast_sampling
   nam%subsamp = subsamp
   nam%network = network
   nam%mpicom = mpicom
   nam%adv_mode = adv_mode
   nam%forced_radii = forced_radii
   nam%rh = rh
   nam%rv = rv
   nam%pos_def_test = pos_def_test
   nam%write_grids = write_grids

   ! dirac_param
   if (ndir>ndirmax) call mpl%abort(subr,'ndir is too large')
   nam%ndir = ndir
   if (ndir>0) nam%londir(1:ndir) = londir(1:ndir)
   if (ndir>0) nam%latdir(1:ndir) = latdir(1:ndir)
   if (ndir>0) nam%levdir(1:ndir) = levdir(1:ndir)
   if (ndir>0) nam%ivdir(1:ndir) = ivdir(1:ndir)
   if (ndir>0) nam%itsdir(1:ndir) = itsdir(1:ndir)

   ! obsop_param
   read(lunit,nml=obsop_param)
   nam%nobs = nobs

   ! output_param
   read(lunit,nml=output_param)
   nam%nldwv = nldwv
   if (nldwv>0) then
      nam%img_ldwv(1:nldwv) = img_ldwv(1:nldwv)
      nam%lon_ldwv(1:nldwv) = lon_ldwv(1:nldwv)
      nam%lat_ldwv(1:nldwv) = lat_ldwv(1:nldwv)
      nam%name_ldwv(1:nldwv) = name_ldwv(1:nldwv)
   end if
   nam%grid_output = grid_output
   nam%grid_resol = grid_resol

   ! Close namelist
   close(unit=lunit)
end if

end subroutine nam_read

!----------------------------------------------------------------------
! Subroutine: nam_read_yaml
! Purpose: read YAML file
!----------------------------------------------------------------------
subroutine nam_read_yaml(nam,mpl,yamlname)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam       ! Namelist
type(mpl_type),intent(inout) :: mpl        ! MPI data
character(len=*),intent(inout) :: yamlname ! YAML name

! Local variables
type(fckit_configuration) :: conf

if (mpl%main) then
   ! Set fckit configuration from yamlname
   conf = fckit_yamlconfiguration(fckit_pathname(trim(yamlname)))

   ! Convert fckit configuration to namelist
   call nam%from_conf(conf)
end if

end subroutine nam_read_yaml

!----------------------------------------------------------------------
! Subroutine: nam_bcast
! Purpose: broadcast
!----------------------------------------------------------------------
subroutine nam_bcast(nam,mpl)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam ! Namelist
type(mpl_type),intent(inout) :: mpl  ! MPI data

! general_param
call mpl%f_comm%broadcast(nam%datadir,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%prefix,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%model,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%verbosity,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%colorlog,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%default_seed,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%repro,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nprocio,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%remap,mpl%rootproc-1)

! driver_param
call mpl%f_comm%broadcast(nam%method,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%strategy,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%new_normality,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%new_cortrack,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%new_corstats,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%new_vbal,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%load_vbal,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%write_vbal,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%new_var,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%load_var,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%write_var,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%new_mom,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%load_mom,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%write_mom,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%new_hdiag,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%write_hdiag,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%new_lct,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%write_lct,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%load_cmat,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%write_cmat,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%new_nicas,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%load_nicas,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%write_nicas,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%new_obsop,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%load_obsop,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%write_obsop,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_vbal,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_adjoints,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_dirac,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_randomization,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_consistency,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_optimality,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_obsop,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_no_obs,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_no_point,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_no_point_mask,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_no_point_nicas,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_set_param_cor,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_set_param_hyb,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_set_param_lct,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_get_param_stddev,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_get_param_cor,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_get_param_hyb,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_get_param_Dloc,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_get_param_lct,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_apply_vbal,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_apply_stddev,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_apply_nicas,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%check_apply_obsop,mpl%rootproc-1)

! model_param
call mpl%f_comm%broadcast(nam%nl,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%levs,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lev2d,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%logpres,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nv,mpl%rootproc-1)
call mpl%broadcast(nam%variables,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nts,mpl%rootproc-1)
call mpl%broadcast(nam%timeslots,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%dts,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nomask,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%wind_filename,mpl%rootproc-1)
call mpl%broadcast(nam%wind_variables,mpl%rootproc-1)

! ens1_param
call mpl%f_comm%broadcast(nam%ens1_ne,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%ens1_nsub,mpl%rootproc-1)

! ens2_param
call mpl%f_comm%broadcast(nam%ens2_ne,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%ens2_nsub,mpl%rootproc-1)

! sampling_param
call mpl%f_comm%broadcast(nam%sam_write,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%sam_write_grids,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%sam_read,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%mask_type,mpl%rootproc-1)
call mpl%broadcast(nam%mask_lu,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%mask_th,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%ncontig_th,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%mask_check,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%draw_type,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%Lcoast,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%rcoast,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nc1,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nc2,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%ntry,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nrep,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nc3,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%dc,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nl0r,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%irmax,mpl%rootproc-1)

! diag_param
call mpl%f_comm%broadcast(nam%ne,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%gen_kurt_th,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%gau_approx,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%avg_nbins,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%vbal_block,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%vbal_rad,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%vbal_dlat,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%vbal_diag_auto,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%vbal_diag_reg,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%var_filter,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%var_niter,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%var_rhflt,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%local_diag,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%local_rad,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%adv_diag,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%adv_type,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%adv_rad,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%adv_niter,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%adv_rhflt,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%adv_valid,mpl%rootproc-1)

! fit_param
call mpl%f_comm%broadcast(nam%minim_algo,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%diag_rhflt,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%diag_rvflt,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%smoothness_penalty,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%fit_dl0,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lct_nscales,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lct_scale_ratio,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lct_cor_min,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lct_diag,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lct_qc_th,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lct_qc_max,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lct_write_cor,mpl%rootproc-1)

! nicas_param
call mpl%f_comm%broadcast(nam%nonunit_diag,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lsqrt,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%resol,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nc1max,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%fast_sampling,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%subsamp,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%network,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%mpicom,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%adv_mode,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%forced_radii,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%rh,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%rv,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%pos_def_test,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%write_grids,mpl%rootproc-1)

! dirac_param
call mpl%f_comm%broadcast(nam%ndir,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%londir,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%latdir,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%levdir,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%ivdir,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%itsdir,mpl%rootproc-1)

! obsop_param
call mpl%f_comm%broadcast(nam%nobs,mpl%rootproc-1)

! output_param
call mpl%f_comm%broadcast(nam%nldwv,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%img_ldwv,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lon_ldwv,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lat_ldwv,mpl%rootproc-1)
call mpl%broadcast(nam%name_ldwv,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%grid_output,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%grid_resol,mpl%rootproc-1)

end subroutine nam_bcast

!----------------------------------------------------------------------
! Subroutine: nam_from_conf
! Purpose: intialize from configuration
!----------------------------------------------------------------------
subroutine nam_from_conf(nam,conf)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam         ! Namelist
type(fckit_configuration),intent(in) :: conf ! Configuration

! Local variables
integer,allocatable :: integer_array(:)
real(kind_real),allocatable :: real_array(:)
logical,allocatable :: logical_array(:)
character(len=:),allocatable :: str
character(len=:),allocatable :: str_array(:)

! general_param
if (conf%has("datadir")) then
   call conf%get_or_die("datadir",str)
   nam%datadir = str
end if
if (conf%has("prefix")) then
   call conf%get_or_die("prefix",str)
   nam%prefix = str
end if
if (conf%has("model")) then
   call conf%get_or_die("model",str)
   nam%model = str
end if
if (conf%has("verbosity")) then
   call conf%get_or_die("verbosity",str)
   nam%verbosity = str
end if
if (conf%has("colorlog")) call conf%get_or_die("colorlog",nam%colorlog)
if (conf%has("default_seed")) call conf%get_or_die("default_seed",nam%default_seed)
if (conf%has("repro")) call conf%get_or_die("repro",nam%repro)
if (conf%has("nprocio")) call conf%get_or_die("nprocio",nam%nprocio)
if (conf%has("remap")) call conf%get_or_die("remap",nam%remap)

! driver_param
if (conf%has("method")) then
   call conf%get_or_die("method",str)
   nam%method = str
end if
if (conf%has("strategy")) then
   call conf%get_or_die("strategy",str)
   nam%strategy = str
end if
if (conf%has("new_normality")) call conf%get_or_die("new_normality",nam%new_normality)
if (conf%has("new_cortrack")) call conf%get_or_die("new_cortrack",nam%new_cortrack)
if (conf%has("new_corstats")) call conf%get_or_die("new_corstats",nam%new_corstats)
if (conf%has("new_vbal")) call conf%get_or_die("new_vbal",nam%new_vbal)
if (conf%has("load_vbal")) call conf%get_or_die("load_vbal",nam%load_vbal)
if (conf%has("write_vbal")) call conf%get_or_die("write_vbal",nam%write_vbal)
if (conf%has("new_var")) call conf%get_or_die("new_var",nam%new_var)
if (conf%has("load_var")) call conf%get_or_die("load_var",nam%load_var)
if (conf%has("write_var")) call conf%get_or_die("write_var",nam%write_var)
if (conf%has("new_mom")) call conf%get_or_die("new_mom",nam%new_mom)
if (conf%has("load_mom")) call conf%get_or_die("load_mom",nam%load_mom)
if (conf%has("write_mom")) call conf%get_or_die("write_mom",nam%write_mom)
if (conf%has("new_hdiag")) call conf%get_or_die("new_hdiag",nam%new_hdiag)
if (conf%has("write_hdiag")) call conf%get_or_die("write_hdiag",nam%write_hdiag)
if (conf%has("new_lct")) call conf%get_or_die("new_lct",nam%new_lct)
if (conf%has("write_lct")) call conf%get_or_die("write_lct",nam%write_lct)
if (conf%has("load_cmat")) call conf%get_or_die("load_cmat",nam%load_cmat)
if (conf%has("write_cmat")) call conf%get_or_die("write_cmat",nam%write_cmat)
if (conf%has("new_nicas")) call conf%get_or_die("new_nicas",nam%new_nicas)
if (conf%has("load_nicas")) call conf%get_or_die("load_nicas",nam%load_nicas)
if (conf%has("write_nicas")) call conf%get_or_die("write_nicas",nam%write_nicas)
if (conf%has("new_obsop")) call conf%get_or_die("new_obsop",nam%new_obsop)
if (conf%has("load_obsop")) call conf%get_or_die("load_obsop",nam%load_obsop)
if (conf%has("write_obsop")) call conf%get_or_die("write_obsop",nam%write_obsop)
if (conf%has("check_vbal")) call conf%get_or_die("check_vbal",nam%check_vbal)
if (conf%has("check_adjoints")) call conf%get_or_die("check_adjoints",nam%check_adjoints)
if (conf%has("check_dirac")) call conf%get_or_die("check_dirac",nam%check_dirac)
if (conf%has("check_randomization")) call conf%get_or_die("check_randomization",nam%check_randomization)
if (conf%has("check_consistency")) call conf%get_or_die("check_consistency",nam%check_consistency)
if (conf%has("check_optimality")) call conf%get_or_die("check_optimality",nam%check_optimality)
if (conf%has("check_obsop")) call conf%get_or_die("check_obsop",nam%check_obsop)
if (conf%has("check_no_obs")) call conf%get_or_die("check_no_obs",nam%check_no_obs)
if (conf%has("check_no_point")) call conf%get_or_die("check_no_point",nam%check_no_point)
if (conf%has("check_no_point_mask")) call conf%get_or_die("check_no_point_mask",nam%check_no_point_mask)
if (conf%has("check_no_point_nicas")) call conf%get_or_die("check_no_point_nicas",nam%check_no_point_nicas)
if (conf%has("check_set_param_cor")) call conf%get_or_die("check_set_param_cor",nam%check_set_param_cor)
if (conf%has("check_set_param_hyb")) call conf%get_or_die("check_set_param_hyb",nam%check_set_param_hyb)
if (conf%has("check_set_param_lct")) call conf%get_or_die("check_set_param_lct",nam%check_set_param_lct)
if (conf%has("check_get_param_stddev")) call conf%get_or_die("check_get_param_stddev",nam%check_get_param_stddev)
if (conf%has("check_get_param_cor")) call conf%get_or_die("check_get_param_cor",nam%check_get_param_cor)
if (conf%has("check_get_param_hyb")) call conf%get_or_die("check_get_param_hyb",nam%check_get_param_hyb)
if (conf%has("check_get_param_Dloc")) call conf%get_or_die("check_get_param_Dloc",nam%check_get_param_Dloc)
if (conf%has("check_get_param_lct")) call conf%get_or_die("check_get_param_lct",nam%check_get_param_lct)
if (conf%has("check_apply_vbal")) call conf%get_or_die("check_apply_vbal",nam%check_apply_vbal)
if (conf%has("check_apply_stddev")) call conf%get_or_die("check_apply_stddev",nam%check_apply_stddev)
if (conf%has("check_apply_nicas")) call conf%get_or_die("check_apply_nicas",nam%check_apply_nicas)
if (conf%has("check_apply_obsop")) call conf%get_or_die("check_apply_obsop",nam%check_apply_obsop)

! model_param
if (conf%has("nl")) call conf%get_or_die("nl",nam%nl)
if (conf%has("levs")) then
   call conf%get_or_die("levs",integer_array)
   nam%levs(1:nam%nl) = integer_array(1:nam%nl)
end if
if (conf%has("lev2d")) then
   call conf%get_or_die("lev2d",str)
   nam%lev2d = str
end if
if (conf%has("logpres")) call conf%get_or_die("logpres",nam%logpres)
if (conf%has("nv")) call conf%get_or_die("nv",nam%nv)
if (conf%has("variables")) then
   call conf%get_or_die("variables",str_array)
   nam%variables(1:nam%nv) = str_array(1:nam%nv)
end if
if (conf%has("nts")) call conf%get_or_die("nts",nam%nts)
if (conf%has("timeslots")) then
   call conf%get_or_die("timeslots",str_array)
   nam%timeslots(1:nam%nts) = str_array(1:nam%nts)
end if
if (conf%has("dts")) call conf%get_or_die("dts",nam%dts)
if (conf%has("nomask")) call conf%get_or_die("nomask",nam%nomask)
if (conf%has("wind_filename")) then
   call conf%get_or_die("wind_filename",str)
   nam%wind_filename = str
end if
if (conf%has("wind_variables")) then
   call conf%get_or_die("wind_variables",str_array)
   nam%wind_variables(1:2) = str_array(1:2)
end if

! ens1_param
if (conf%has("ens1_ne")) call conf%get_or_die("ens1_ne",nam%ens1_ne)
if (conf%has("ens1_nsub")) call conf%get_or_die("ens1_nsub",nam%ens1_nsub)

! ens2_param
if (conf%has("ens2_ne")) call conf%get_or_die("ens2_ne",nam%ens2_ne)
if (conf%has("ens2_nsub")) call conf%get_or_die("ens2_nsub",nam%ens2_nsub)

! sampling_param
if (conf%has("sam_read")) call conf%get_or_die("sam_read",nam%sam_read)
if (conf%has("sam_write")) call conf%get_or_die("sam_write",nam%sam_write)
if (conf%has("sam_write_grids")) call conf%get_or_die("sam_write_grids",nam%sam_write_grids)
if (conf%has("mask_type")) then
   call conf%get_or_die("mask_type",str)
   nam%mask_type = str
end if
if (conf%has("mask_lu")) then
   call conf%get_or_die("mask_lu",str_array)
   nam%mask_lu(1:nam%nv) = str_array(1:nam%nv)
end if
if (conf%has("mask_th")) then
   call conf%get_or_die("mask_th",real_array)
   nam%mask_th(1:nam%nv) = real_array(1:nam%nv)
end if
if (conf%has("ncontig_th")) call conf%get_or_die("ncontig_th",nam%ncontig_th)
if (conf%has("mask_check")) call conf%get_or_die("mask_check",nam%mask_check)
if (conf%has("draw_type")) then
   call conf%get_or_die("draw_type",str)
   nam%draw_type = str
end if
if (conf%has("Lcoast")) call conf%get_or_die("Lcoast",nam%Lcoast)
if (conf%has("rcoast")) call conf%get_or_die("rcoast",nam%rcoast)
if (conf%has("nc1")) call conf%get_or_die("nc1",nam%nc1)
if (conf%has("nc2")) call conf%get_or_die("nc2",nam%nc2)
if (conf%has("ntry")) call conf%get_or_die("ntry",nam%ntry)
if (conf%has("nrep")) call conf%get_or_die("nrep",nam%nrep)
if (conf%has("nc3")) call conf%get_or_die("nc3",nam%nc3)
if (conf%has("dc")) call conf%get_or_die("dc",nam%dc)
if (conf%has("nl0r")) call conf%get_or_die("nl0r",nam%nl0r)
if (conf%has("irmax")) call conf%get_or_die("irmax",nam%irmax)

! diag_param
if (conf%has("ne")) call conf%get_or_die("ne",nam%ne)
if (conf%has("gen_kurt_th")) call conf%get_or_die("gen_kurt_th",nam%gen_kurt_th)
if (conf%has("gau_approx")) call conf%get_or_die("gau_approx",nam%gau_approx)
if (conf%has("avg_nbins")) call conf%get_or_die("avg_nbins",nam%avg_nbins)
if (conf%has("vbal_block")) then
   call conf%get_or_die("vbal_block",logical_array)
   nam%vbal_block(1:nam%nv*(nam%nv-1)/2) = logical_array(1:nam%nv*(nam%nv-1)/2)
end if
if (conf%has("vbal_rad")) call conf%get_or_die("vbal_rad",nam%vbal_rad)
if (conf%has("vbal_dlat")) call conf%get_or_die("vbal_dlat",nam%vbal_dlat)
if (conf%has("vbal_diag_auto")) then
   call conf%get_or_die("vbal_diag_auto",logical_array)
   nam%vbal_diag_auto(1:nam%nv*(nam%nv-1)/2) = logical_array(1:nam%nv*(nam%nv-1)/2)
end if
if (conf%has("vbal_diag_reg")) then
   call conf%get_or_die("vbal_diag_reg",logical_array)
   nam%vbal_diag_reg(1:nam%nv*(nam%nv-1)/2) = logical_array(1:nam%nv*(nam%nv-1)/2)
end if
if (conf%has("var_filter")) call conf%get_or_die("var_filter",nam%var_filter)
if (conf%has("var_niter")) call conf%get_or_die("var_niter",nam%var_niter)
if (conf%has("var_rhflt")) call conf%get_or_die("var_rhflt",nam%var_rhflt)
if (conf%has("local_diag")) call conf%get_or_die("local_diag",nam%local_diag)
if (conf%has("local_rad")) call conf%get_or_die("local_rad",nam%local_rad)
if (conf%has("adv_diag")) call conf%get_or_die("adv_diag",nam%adv_diag)
if (conf%has("adv_type")) then
   call conf%get_or_die("adv_type",str)
   nam%adv_type = str
end if
if (conf%has("adv_rad")) call conf%get_or_die("adv_rad",nam%adv_rad)
if (conf%has("adv_niter")) call conf%get_or_die("adv_niter",nam%adv_niter)
if (conf%has("adv_rhflt")) call conf%get_or_die("adv_rhflt",nam%adv_rhflt)
if (conf%has("adv_valid")) call conf%get_or_die("adv_valid",nam%adv_valid)

! fit_param
if (conf%has("minim_algo")) then
   call conf%get_or_die("minim_algo",str)
   nam%minim_algo = str
end if
if (conf%has("diag_rhflt")) call conf%get_or_die("diag_rhflt",nam%diag_rhflt)
if (conf%has("diag_rvflt")) call conf%get_or_die("diag_rvflt",nam%diag_rvflt)
if (conf%has("smoothness_penalty")) call conf%get_or_die("smoothness_penalty",nam%smoothness_penalty)
if (conf%has("fit_dl0")) call conf%get_or_die("fit_dl0",nam%fit_dl0)
if (conf%has("lct_nscales")) call conf%get_or_die("lct_nscales",nam%lct_nscales)
if (conf%has("lct_scale_ratio")) call conf%get_or_die("lct_scale_ratio",nam%lct_scale_ratio)
if (conf%has("lct_cor_min")) call conf%get_or_die("lct_cor_min",nam%lct_cor_min)
if (conf%has("lct_diag")) then
   call conf%get_or_die("lct_diag",logical_array)
   nam%lct_diag(1:nam%lct_nscales) = logical_array(1:nam%lct_nscales)
end if
if (conf%has("lct_qc_th")) call conf%get_or_die("lct_qc_th",nam%lct_qc_th)
if (conf%has("lct_qc_max")) call conf%get_or_die("lct_qc_max",nam%lct_qc_max)
if (conf%has("lct_write_cor")) call conf%get_or_die("lct_write_cor",nam%lct_write_cor)

! nicas_param
if (conf%has("nonunit_diag")) call conf%get_or_die("nonunit_diag",nam%nonunit_diag)
if (conf%has("lsqrt")) call conf%get_or_die("lsqrt",nam%lsqrt)
if (conf%has("resol")) call conf%get_or_die("resol",nam%resol)
if (conf%has("nc1max")) call conf%get_or_die("nc1max",nam%nc1max)
if (conf%has("fast_sampling")) call conf%get_or_die("fast_sampling",nam%fast_sampling)
if (conf%has("subsamp")) then
   call conf%get_or_die("subsamp",str)
   nam%subsamp = str
end if
if (conf%has("network")) call conf%get_or_die("network",nam%network)
if (conf%has("mpicom")) call conf%get_or_die("mpicom",nam%mpicom)
if (conf%has("adv_mode")) call conf%get_or_die("adv_mode",nam%adv_mode)
if (conf%has("forced_radii")) call conf%get_or_die("forced_radii",nam%forced_radii)
if (conf%has("rh")) call conf%get_or_die("rh",nam%rh)
if (conf%has("rv")) call conf%get_or_die("rv",nam%rv)
if (conf%has("pos_def_test")) call conf%get_or_die("pos_def_test",nam%pos_def_test)
if (conf%has("write_grids")) call conf%get_or_die("write_grids",nam%write_grids)

! dirac_param
if (conf%has("ndir")) call conf%get_or_die("ndir",nam%ndir)
if (conf%has("londir")) then
   call conf%get_or_die("londir",real_array)
   nam%londir(1:nam%ndir) = real_array(1:nam%ndir)
end if
if (conf%has("latdir")) then
   call conf%get_or_die("latdir",real_array)
   nam%latdir(1:nam%ndir) = real_array(1:nam%ndir)
end if
if (conf%has("levdir")) then
   call conf%get_or_die("levdir",integer_array)
   nam%levdir(1:nam%ndir) = integer_array(1:nam%ndir)
end if
if (conf%has("ivdir")) then
   call conf%get_or_die("ivdir",integer_array)
   nam%ivdir(1:nam%ndir) = integer_array(1:nam%ndir)
end if
if (conf%has("itsdir")) then
   call conf%get_or_die("itsdir",integer_array)
   nam%itsdir(1:nam%ndir) = integer_array(1:nam%ndir)
end if

! obsop_param
if (conf%has("nobs")) call conf%get_or_die("nobs",nam%nobs)

! output_param
if (conf%has("nldwv")) call conf%get_or_die("nldwv",nam%nldwv)
if (conf%has("img_ldwv")) then
   call conf%get_or_die("img_ldwv",integer_array)
   nam%img_ldwv(1:nam%nldwv) = integer_array(1:nam%nldwv)
end if
if (conf%has("lon_ldwv")) then
   call conf%get_or_die("lon_ldwv",real_array)
   nam%lon_ldwv(1:nam%nldwv) = real_array(1:nam%nldwv)
end if
if (conf%has("lat_ldwv")) then
   call conf%get_or_die("lat_ldwv",real_array)
   nam%lat_ldwv(1:nam%nldwv) = real_array(1:nam%nldwv)
end if
if (conf%has("name_ldwv")) then
   call conf%get_or_die("name_ldwv",str_array)
   nam%name_ldwv(1:nam%nldwv) = str_array(1:nam%nldwv)
end if
if (conf%has("grid_output")) call conf%get_or_die("grid_output",nam%grid_output)
if (conf%has("grid_resol")) call conf%get_or_die("grid_resol",nam%grid_resol)

end subroutine nam_from_conf

!----------------------------------------------------------------------
! Subroutine: nam_check
! Purpose: check namelist parameters
!----------------------------------------------------------------------
subroutine nam_check(nam,mpl)

implicit none

! Passed variable
class(nam_type),intent(inout) :: nam ! Namelist
type(mpl_type),intent(inout) :: mpl  ! MPI data

! Local variables
integer :: iv,its,il,idir,ildwv
character(len=2) :: ivchar,itschar,ildwvchar
character(len=1024),parameter :: subr = 'nam_check'

! Check maximum sizes
if (nam%nl>nlmax) call mpl%abort(subr,'nl is too large')
if (nam%nv>nvmax) call mpl%abort(subr,'nv is too large')
if (nam%nts>ntsmax) call mpl%abort(subr,'nts is too large')
if (nam%nc3>nc3max) call mpl%abort(subr,'nc3 is too large')
if (nam%lct_nscales>nscalesmax) call mpl%abort(subr,'lct_nscales is too large')
if (nam%ndir>ndirmax) call mpl%abort(subr,'ndir is too large')
if (nam%nldwv>nldwvmax) call mpl%abort(subr,'nldwv is too large')

! Namelist parameters normalization (meters to radians and degrees to radians)
nam%Lcoast = nam%Lcoast/req
nam%dc = nam%dc/req
nam%vbal_rad = nam%vbal_rad/req
nam%vbal_dlat = nam%vbal_dlat*deg2rad
nam%var_rhflt = nam%var_rhflt/req
nam%local_rad = nam%local_rad/req
nam%adv_rad = nam%adv_rad/req
nam%adv_rhflt = nam%adv_rhflt/req
nam%diag_rhflt = nam%diag_rhflt/req
nam%rh = nam%rh/req
if (nam%ndir>0) nam%londir(1:nam%ndir) = nam%londir(1:nam%ndir)*deg2rad
if (nam%ndir>0) nam%latdir(1:nam%ndir) = nam%latdir(1:nam%ndir)*deg2rad
if (nam%nldwv>0) nam%lon_ldwv(1:nam%nldwv) = nam%lon_ldwv(1:nam%nldwv)*deg2rad
if (nam%nldwv>0) nam%lat_ldwv(1:nam%nldwv) = nam%lat_ldwv(1:nam%nldwv)*deg2rad
nam%grid_resol = nam%grid_resol/req

! Check general_param
if (trim(nam%datadir)=='') call mpl%abort(subr,'datadir not specified')
if (trim(nam%prefix)=='') call mpl%abort(subr,'prefix not specified')
select case (trim(nam%model))
case ('aro','arp','fv3','gem','geos','gfs','ifs','mpas','nemo','online','qg','res','wrf')
case default
   call mpl%abort(subr,'wrong model')
end select
select case (trim(nam%verbosity))
case ('all','main','none')
case default
   call mpl%abort(subr,'wrong verbosity level')
end select
if (nam%nprocio<1) call mpl%abort(subr,'number of I/O tasks should be positive')
if (nam%nprocio>mpl%nproc) then
   call mpl%warning(subr,'number of I/O tasks should be smaller than the total number of tasks, resetting nprocio')
   nam%nprocio = mpl%nproc
end if

! Check driver_param
if (nam%new_hdiag.or.nam%check_optimality) then
   select case (trim(nam%method))
   case ('cor','loc','hyb-avg','hyb-rnd','dual-ens')
   case default
      call mpl%abort(subr,'wrong method')
   end select
end if
if (nam%new_lct) then
   if (trim(nam%method)/='cor') call mpl%abort(subr,'cor method required for new_lct')
end if
if (nam%new_hdiag.or.nam%new_lct.or.nam%load_cmat.or.nam%new_nicas.or.nam%load_nicas) then
   select case (trim(nam%strategy))
   case ('diag_all','common','common_univariate','common_weighted','specific_univariate','specific_multivariate')
   case default
      call mpl%abort(subr,'wrong strategy')
   end select
end if
if (nam%new_vbal.and.nam%load_vbal) call mpl%abort(subr,'new_vbal and load_vbal are exclusive')
if (nam%new_var.and.nam%load_var) call mpl%abort(subr,'new_var and load_var are exclusive')
if (nam%new_mom.and.nam%load_mom) call mpl%abort(subr,'new_mom and load_mom are exclusive')
if (nam%new_hdiag.and.nam%new_lct) call mpl%abort(subr,'new_hdiag and new_lct are exclusive')
if ((nam%new_hdiag.or.nam%new_lct).and.nam%load_cmat) call mpl%abort(subr,'new_hdiag or new_lct and load_cmat are exclusive')
if (nam%new_nicas.and.nam%load_nicas) call mpl%abort(subr,'new_nicas and load_nicas are exclusive')
if (nam%new_obsop.and.nam%load_obsop) call mpl%abort(subr,'new_obsop and load_obsop are exclusive')
if (nam%check_vbal.and..not.(nam%new_vbal.or.nam%load_vbal)) call mpl%abort(subr,'new_vbal or load_vbal required for check_vbal')
if ((nam%new_hdiag.or.nam%new_lct).and.(.not.(nam%new_mom.or.nam%load_mom))) &
 & call mpl%abort(subr,'new_mom or load_mom required for new_hdiag and new_lct')
if (nam%check_adjoints.and..not.(nam%new_vbal.or.nam%load_vbal.or.nam%new_nicas.or.nam%load_nicas.or.nam%new_obsop &
 & .or.nam%load_obsop)) call mpl%abort(subr,'new or load for vbal, nicas or obsop required for check_adjoints')
if (nam%check_dirac.and..not.(nam%new_vbal.or.nam%load_vbal.or.nam%new_nicas.or.nam%load_nicas)) &
 & call mpl%abort(subr,'new or load for vbal or nicas required for check_dirac')
if (nam%check_randomization) then
   if (trim(nam%method)/='cor') call mpl%abort(subr,'cor method required for check_randomization')
   if (.not.nam%new_nicas) call mpl%abort(subr,'new_nicas required for check_randomization')
end if
if (nam%check_consistency) then
   if (trim(nam%method)/='cor') call mpl%abort(subr,'cor method required for check_consistency')
   if (.not.nam%new_nicas) call mpl%abort(subr,'new_nicas required for check_consistency')
end if
if (nam%check_optimality) then
   if (trim(nam%method)/='cor') call mpl%abort(subr,'cor method required for check_optimality')
   if (.not.nam%new_nicas) call mpl%abort(subr,'new_nicas required for check_optimality')
   if (.not.nam%write_hdiag) call mpl%abort(subr,'write_hdiag required for check_optimality')
end if
if (nam%check_obsop.and..not.(nam%new_obsop.or.nam%load_obsop)) &
 & call mpl%abort(subr,'new or load obsop required for check_obsop')
if (nam%check_no_obs.and..not.(nam%new_obsop.or.nam%load_obsop)) &
 & call mpl%abort(subr,'new or load_obsop required for check_no_obs')
if (nam%check_no_obs.and.(mpl%nproc<2)) call mpl%abort(subr,'at least 2 MPI tasks required for check_no_obs')
if (nam%check_no_point.and.nam%check_no_point_mask) &
 & call mpl%abort(subr,'check_no_point and check_no_point_mask are exclusive')
if (nam%check_no_point.and.(mpl%nproc<2)) call mpl%abort(subr,'at least 2 MPI tasks required for check_no_point')
if (nam%check_no_point_mask.and.(mpl%nproc<2)) call mpl%abort(subr,'at least 2 MPI tasks required for check_no_point_mask')
if (nam%check_no_point_nicas.and..not.(nam%new_nicas.or.nam%load_nicas)) &
 & call mpl%abort(subr,'new_nicas or load_nicas required for check_no_point_nicas')
if (nam%check_no_point_nicas.and.(mpl%nproc<2)) call mpl%abort(subr,'at least 2 MPI tasks required for check_no_point_nicas')
if ((nam%check_set_param_cor.or.nam%check_set_param_cor.or.nam%check_set_param_cor).and..not.nam%new_nicas) &
 & call mpl%abort(subr,'new_nicas required for check_set_param_[...]')
if (nam%check_get_param_stddev.and..not.(nam%new_var)) &
 & call mpl%abort(subr,'new_var required for check_get_param_stddev')
if (nam%check_get_param_cor.and..not.(nam%new_hdiag.and.(trim(nam%method)=='cor'))) &
 & call mpl%abort(subr,'new_hdiag and cor method required for check_get_param_cor')
if (nam%check_get_param_hyb.and..not.(nam%new_hdiag.and.(trim(nam%method)=='hyb-avg'))) &
 & call mpl%abort(subr,'new_hdiag and hyb-avg method required for check_get_param_hyb')
if (nam%check_get_param_Dloc.and..not.(nam%new_hdiag.and.(trim(nam%method)=='loc'))) &
 & call mpl%abort(subr,'new_hdiag and loc method required for check_get_param_Dloc')
if (nam%check_get_param_lct.and..not.(nam%new_lct.and.(nam%lct_nscales==2))) &
 & call mpl%abort(subr,'new_lct and lct_nscales = 2 required for check_get_param_lct')
if (nam%check_apply_vbal.and..not.(nam%new_vbal.or.nam%load_vbal)) &
 & call mpl%abort(subr,'new_vbal or load_vbal required for check_apply_vbal')
if (nam%check_apply_stddev.and..not.(nam%new_var.or.nam%load_var)) &
 & call mpl%abort(subr,'new_var or load_var required for check_apply_stddev')
if (nam%check_apply_nicas.and..not.(nam%new_nicas.or.nam%load_nicas)) &
 & call mpl%abort(subr,'new_nicas or load_nicas required for check_apply_nicas')
if (nam%check_apply_obsop.and..not.(nam%new_obsop.or.nam%load_obsop)) &
 & call mpl%abort(subr,'new_obsop or load_obsop required for check_apply_obsop')

! Check model_param
if (nam%nl<=0) call mpl%abort(subr,'nl should be positive')
do il=1,nam%nl
   if (nam%levs(il)<=0) call mpl%abort(subr,'levs should be positive')
   if (count(nam%levs(1:nam%nl)==nam%levs(il))>1) call mpl%abort(subr,'redundant levels')
end do
if ((trim(nam%lev2d)/='first').and.(trim(nam%lev2d)/='last')) call mpl%abort(subr,'wrong lev2d value')
if (nam%new_vbal.or.nam%load_vbal.or.nam%new_var.or.nam%load_var.or.nam%new_hdiag.or.nam%new_lct.or.nam%load_cmat &
 & .or.nam%new_nicas.or.nam%load_nicas) then
   if (nam%nv<=0) call mpl%abort(subr,'nv should be positive')
   do iv=1,nam%nv
      write(ivchar,'(i2.2)') iv
      if (trim(nam%variables(iv))=='') call mpl%abort(subr,'variables not specified for variable '//ivchar)
   end do
   if (nam%nts<=0) call mpl%abort(subr,'nts should be positive')
   do its=1,nam%nts
      write(itschar,'(i2.2)') its
      if (trim(nam%timeslots(its))=='') call mpl%abort(subr,'timeslots not specified for '//itschar)
   end do
   if (.not.(nam%dts>0.0)) call mpl%abort(subr,'dts should be positive')
end if

! Check ens1_param
if (nam%new_normality.or.nam%new_cortrack.or.nam%new_corstats.or.nam%new_vbal.or.nam%new_var.or.nam%new_hdiag.or.nam%new_lct &
 & .or.nam%check_randomization.or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%ens1_nsub<1) call mpl%abort(subr,'ens1_nsub should be positive')
   if (mod(nam%ens1_ne,nam%ens1_nsub)/=0) call mpl%abort(subr,'ens1_nsub should be a divider of ens1_ne')
   if (nam%ens1_ne/nam%ens1_nsub<=3) call mpl%abort(subr,'ens1_ne/ens1_nsub should be larger than 3')
end if

! Check ens2_param
if (nam%new_hdiag.and.((trim(nam%method)=='hyb-rnd').or.(trim(nam%method)=='dual-ens'))) then
   if (nam%ens2_nsub<1) call mpl%abort(subr,'ens2_nsub should be non-negative')
   if (mod(nam%ens2_ne,nam%ens2_nsub)/=0) call mpl%abort(subr,'ens2_nsub should be a divider of ens2_ne')
   if (nam%ens2_ne/nam%ens2_nsub<=3) call mpl%abort(subr,'ens2_ne/ens2_nsub should be larger than 3')
end if

! Check sampling_param
if (nam%new_vbal.or.nam%new_hdiag.or.nam%new_lct.or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%sam_write.and.nam%sam_read) call mpl%abort(subr,'sam_write and sam_read are both true')
   if (nam%sam_write_grids.and.(.not.nam%sam_write)) call mpl%abort(subr,'sam_write required for sam_write_grids')
   select case (trim(nam%draw_type))
   case ('random_uniform')
   case ('random_coast')
      if (.not.(nam%Lcoast>0.0)) call mpl%abort (subr,'Lcoast should be positive')
      if (.not.(nam%rcoast>0.0)) call mpl%abort (subr,'rcoast should be positive')
   case default
      call mpl%abort(subr,'wrong draw_type')
   end select
   select case (trim(nam%mask_type))
   case ('ldwv')
      if (nam%nldwv<=0) call mpl%abort(subr,'nldwv should not be negative for mask_type = ldwv')
   case ('stddev')
      do iv=1,nam%nv
         select case (trim(nam%mask_lu(iv)))
         case ('lower','upper')
         case default
            call mpl%abort(subr,'wrong mask_lu')
         end select
      end do
   end select
   if (nam%nc1<3) call mpl%abort(subr,'nc1 should be larger than 2')
   if (nam%new_vbal.or.(nam%new_hdiag.and.(nam%local_diag.or.nam%adv_diag))) then
      if (nam%nc2<3) call mpl%abort(subr,'nc2 should be larger than 2')
   else
      if (nam%nc2<0) then
          call mpl%warning(subr,'nc2 should be set non-negative, resetting nc2 to zero')
          nam%nc2 = 0
      end if
   end if
   if (nam%new_lct) then
      if (nam%nc2/=nam%nc1) then
         call mpl%warning(subr,'nc2 should be equal to nc1 for new_lct, resetting nc2 to nc1')
         nam%nc2 = nam%nc1
      end if
   end if
end if
if (nam%new_vbal.or.nam%new_hdiag.or.nam%new_lct.or.nam%new_nicas) then
   if (nam%ntry<=0) call mpl%abort(subr,'ntry should be positive')
   if (nam%nrep<0) call mpl%abort(subr,'nrep should be non-negative')
end if
if (nam%new_hdiag.or.nam%new_lct.or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%nc3<=0) call mpl%abort(subr,'nc3 should be positive')
end if
if (nam%new_hdiag.or.nam%check_consistency.or.nam%check_optimality) then
   if (.not.(nam%dc>0.0)) call mpl%abort(subr,'dc should be positive')
end if
if (nam%new_hdiag.or.nam%new_lct.or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%nl0r<1) call mpl%abort (subr,'nl0r should be positive')
end if
if (nam%new_hdiag.or.nam%check_consistency.or.nam%check_optimality) then
   if (nam%irmax<1) call mpl%abort (subr,'irmax should be positive')
end if

! Check diag_param
if (nam%new_vbal) then
   if (nam%nv<2) call mpl%abort(subr,'at least two variables required to diagnose vertical balance')
   if (.not.(any(nam%vbal_block(1:nam%nv*(nam%nv-1)/2)))) &
 & call mpl%abort(subr,'no block selected for the vertical balance diagnostics')
   if ((.not.(nam%vbal_rad>0.0)).and.(.not.(nam%vbal_dlat>0.0))) call mpl%abort(subr,'vbal_rad or vbal_dlat should be positive')
end if
if (nam%new_var) then
   if (nam%var_filter) then
      if (nam%var_niter<=0) call mpl%abort(subr,'var_niter should be positive')
      if (.not.(nam%var_rhflt>0.0)) call mpl%abort(subr,'var_rhflt should be positive')
   end if
end if
if (nam%new_hdiag.or.nam%check_consistency.or.nam%check_optimality) then
   select case (trim(nam%method))
   case ('loc','hyb-avg','hyb-rnd','dual-ens')
      if (nam%ne<=3) call mpl%abort(subr,'ne should be larger than 3')
   end select
   if (.not.(nam%gen_kurt_th>0.0)) call mpl%abort(subr,'gen_kurt_th should be positive')
   if (nam%local_diag) then
      if (.not.(nam%local_rad>0.0)) call mpl%abort(subr,'local_rad should be positive')
   end if
   if (nam%adv_diag) then
      if (.not.(nam%adv_rad>0.0)) call mpl%abort(subr,'adv_rad should be positive')
      if (nam%adv_niter<=0) call mpl%abort(subr,'adv_niter should be positive')
      if (.not.(nam%adv_rhflt>0.0)) call mpl%abort(subr,'adv_rhflt should be positive')
      if (nam%adv_valid<0.0) call mpl%abort(subr,'adv_valid should be non-negative')
      if (nam%adv_valid>1.0) call mpl%abort(subr,'adv_valid should be not be higher than 1.0')
   end if
end if

! Check fit_param
if (nam%new_hdiag.or.nam%new_lct.or.nam%check_consistency.or.nam%check_optimality) then
   select case (trim(nam%minim_algo))
   case ('none','fast','hooke','praxis')
   case default
      call mpl%abort(subr,'wrong minim_algo')
   end select
   if (nam%new_lct.and.((trim(nam%minim_algo)=='none').or.(trim(nam%minim_algo)=='fast'))) &
 & call mpl%abort(subr,'wrong minim_algo for LCT')
   if (nam%diag_rhflt<0.0) call mpl%abort(subr,'diag_rhflt should be non-negative')
   if (nam%diag_rvflt<0.0) call mpl%abort(subr,'diag_rvflt should be non-negative')
   if (nam%smoothness_penalty<0.0) call mpl%abort(subr,'smoothness_penalty should be non-negative')
   if (nam%fit_dl0<=0) call mpl%abort(subr,'fit_dl0 should be postive')
end if
if (nam%new_lct) then
   if (nam%lct_nscales<=0) call mpl%abort(subr,'lct_nscales should be postive')
   if (.not.(nam%lct_scale_ratio>0.0)) call mpl%abort(subr,'lct_scale_ratio should be postive')
   if (nam%lct_cor_min<0) call mpl%abort(subr,'lct_cor_min should be non-negative')
   if ((nam%lct_qc_th<-1.0).or.(nam%lct_qc_th>1.0)) call mpl%abort(subr,'lct_qc_th should be between -1 and 1')
   if ((nam%lct_qc_max<-1.0).or.(nam%lct_qc_max>1.0)) call mpl%abort(subr,'lct_qc_max should be between -1 and 1')
end if

! Check ensemble sizes
if (nam%new_hdiag) then
   if (trim(nam%method)/='cor') then
      if (nam%ne>nam%ens1_ne) call mpl%warning(subr,'ensemble size larger than ens1_ne (might enhance sampling noise)')
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd','dual-ens')
         if (nam%ne>nam%ens2_ne) call mpl%warning(subr,'ensemble size larger than ens2_ne (might enhance sampling noise)')
      end select
   end if
end if

! Check nicas_param
if (nam%new_nicas.or.nam%load_nicas) then
   if (nam%nonunit_diag) then
      if (nam%method=='cor') call mpl%abort(subr,'nonunit_diag is inconsistent with correlation, use variance operator instead')
   end if
   if (nam%lsqrt) then
      if (nam%mpicom==1) call mpl%abort(subr,'mpicom should be 2 for square-root application')
   end if
   if (trim(nam%method)=='specific_multivariate') then
      if (.not.nam%lsqrt) call mpl%abort(subr,'square-root formulation required for specific multivariate strategy')
   end if
   if (nam%check_randomization) then
      if (.not.nam%lsqrt) call mpl%abort(subr,'lsqrt required for check_randomization')
      if (.not.nam%forced_radii) call mpl%abort(subr,'forced_radii required for check_randomization')
   end if
   if (nam%check_consistency) then
      if (.not.nam%lsqrt) call mpl%abort(subr,'lsqrt required for check_consistency')
      if (.not.nam%forced_radii) call mpl%abort(subr,'forced_radii required for check_consistency')
   end if
   if (nam%check_optimality) then
      if (.not.nam%lsqrt) call mpl%abort(subr,'lsqrt required for check_optimality')
      if (.not.nam%forced_radii) call mpl%abort(subr,'forced_radii required for check_optimality')
   end if
   if (nam%new_nicas) then
      if (.not.(nam%resol>0.0)) call mpl%abort(subr,'resol should be positive')
      if (nam%nc1max<=0) call mpl%abort(subr,'nc1max should be positive')
   end if
   if (nam%new_nicas.or.nam%load_nicas) then
      if ((nam%mpicom/=1).and.(nam%mpicom/=2)) call mpl%abort(subr,'mpicom should be 1 or 2')
   end if
   if (nam%forced_radii) then
      if (nam%new_hdiag) then
         select case (trim(nam%method))
         case ('hyb-avg','hyb-rnd')
         case default
           call mpl%abort(subr,'new_diag forbidden for forced_radii (except for static hybridization)')
         end select
      end if
      if (nam%new_lct.or.nam%load_cmat) call mpl%abort(subr,'new_lct and load_cmat forbidden for forced_radii')
      if (nam%rh<0.0) call mpl%abort(subr,'rh should be non-negative')
      if (nam%rv<0.0) call mpl%abort(subr,'rv should be non-negative')
   end if
   if (abs(nam%adv_mode)>1) call mpl%abort(subr,'nam%adv_mode should be -1, 0 or 1')
   select case (trim(nam%subsamp))
   case ('h','hv','vh','hvh')
   case default
      call mpl%abort(subr,'wrong subsampling structure for NICAS')
   end select
end if
if (nam%write_grids.and.(.not.nam%new_nicas)) call mpl%abort(subr,'new_nicas required for write_grids')

! Check dirac_param
if (nam%new_cortrack.or.nam%check_dirac) then
   if (nam%ndir<1) call mpl%abort(subr,'ndir should be positive')
   do idir=1,nam%ndir
      if ((nam%londir(idir)<-pi).or.(nam%londir(idir)>pi)) call mpl%abort(subr,'londir should lie between -180 and 180')
      if ((nam%latdir(idir)<-0.5*pi).or.(nam%latdir(idir)>0.5*pi)) call mpl%abort(subr,'latdir should lie between -90 and 90')
      if (.not.any(nam%levs(1:nam%nl)==nam%levdir(idir))) call mpl%abort(subr,'wrong level for a Dirac')
      if ((nam%ivdir(idir)<1).or.(nam%ivdir(idir)>nam%nv)) call mpl%abort(subr,'wrong variable for a Dirac')
      if ((nam%itsdir(idir)<1).or.(nam%itsdir(idir)>nam%nts)) call mpl%abort(subr,'wrong timeslots for a Dirac')
   end do
end if

! Check output_param
if (nam%new_hdiag) then
   if (nam%nldwv<0) call mpl%abort(subr,'nldwv should be non-negative')
   if (nam%nldwv>0) then
      if (.not.nam%local_diag) call mpl%abort(subr,'local_diag required for nldwv>0')
         if (.not.all(nam%img_ldwv(1:nam%nldwv)>0)) then
         if (any(nam%lon_ldwv(1:nam%nldwv)<-pi).or.any(nam%lon_ldwv(1:nam%nldwv)>pi)) &
       & call mpl%abort(subr,'lon_ldwv should lie between -180 and 180')
         if (any(nam%lat_ldwv(1:nam%nldwv)<-0.5*pi).or.any(nam%lat_ldwv(1:nam%nldwv)>0.5*pi)) &
       & call mpl%abort(subr,'lat_ldwv should lie between -90 and 90')
         do ildwv=1,nam%nldwv
            write(ildwvchar,'(i2.2)') ildwv
            if (trim(nam%name_ldwv(ildwv))=='') call mpl%abort(subr,'name_ldwv not specified for profile '//ildwvchar)
         end do
      end if
   end if
end if
if (nam%new_hdiag.or.nam%new_nicas.or.nam%check_adjoints.or.nam%check_dirac.or.nam%check_randomization.or.nam%new_lct) then
   if (nam%grid_output) then
      if (.not.(nam%grid_resol>0.0)) call mpl%abort(subr,'grid_resol should be positive')
   end if
end if

end subroutine nam_check

!----------------------------------------------------------------------
! Subroutine: nam_write
! Purpose: write namelist parameters into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine nam_write(nam,mpl,ncid)

implicit none

! Passed variable
class(nam_type),intent(in) :: nam   ! Namelist
type(mpl_type),intent(inout) :: mpl ! MPI data
integer,intent(in),optional :: ncid ! NetCDF file ID

! Local variables
integer :: lncid
real(kind_real),allocatable :: londir(:),latdir(:),lon_ldwv(:),lat_ldwv(:)

! Set ncid
lncid = mpl%msv%vali
if (present(ncid)) lncid = ncid

! general_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','General parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','datadir',nam%datadir)
call mpl%write(lncid,'nam','prefix',nam%prefix)
call mpl%write(lncid,'nam','model',nam%model)
call mpl%write(lncid,'nam','verbosity',nam%verbosity)
call mpl%write(lncid,'nam','colorlog',nam%colorlog)
call mpl%write(lncid,'nam','default_seed',nam%default_seed)
call mpl%write(lncid,'nam','repro',nam%repro)
call mpl%write(lncid,'nam','nprocio',nam%nprocio)
call mpl%write(lncid,'nam','remap',nam%remap)

! driver_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','Driver parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','method',nam%method)
call mpl%write(lncid,'nam','strategy',nam%strategy)
call mpl%write(lncid,'nam','new_normality',nam%new_normality)
call mpl%write(lncid,'nam','new_cortrack',nam%new_cortrack)
call mpl%write(lncid,'nam','new_corstats',nam%new_corstats)
call mpl%write(lncid,'nam','new_vbal',nam%new_vbal)
call mpl%write(lncid,'nam','load_vbal',nam%load_vbal)
call mpl%write(lncid,'nam','write_vbal',nam%write_vbal)
call mpl%write(lncid,'nam','new_var',nam%new_var)
call mpl%write(lncid,'nam','load_var',nam%load_var)
call mpl%write(lncid,'nam','write_var',nam%write_var)
call mpl%write(lncid,'nam','new_mom',nam%new_mom)
call mpl%write(lncid,'nam','load_mom',nam%load_mom)
call mpl%write(lncid,'nam','write_mom',nam%write_mom)
call mpl%write(lncid,'nam','new_hdiag',nam%new_hdiag)
call mpl%write(lncid,'nam','write_hdiag',nam%write_hdiag)
call mpl%write(lncid,'nam','new_lct',nam%new_lct)
call mpl%write(lncid,'nam','write_lct',nam%write_lct)
call mpl%write(lncid,'nam','load_cmat',nam%load_cmat)
call mpl%write(lncid,'nam','write_cmat',nam%write_cmat)
call mpl%write(lncid,'nam','new_nicas',nam%new_nicas)
call mpl%write(lncid,'nam','load_nicas',nam%load_nicas)
call mpl%write(lncid,'nam','write_nicas',nam%write_nicas)
call mpl%write(lncid,'nam','new_obsop',nam%new_obsop)
call mpl%write(lncid,'nam','load_obsop',nam%load_obsop)
call mpl%write(lncid,'nam','write_obsop',nam%write_obsop)
call mpl%write(lncid,'nam','check_vbal',nam%check_vbal)
call mpl%write(lncid,'nam','check_adjoints',nam%check_adjoints)
call mpl%write(lncid,'nam','check_dirac',nam%check_dirac)
call mpl%write(lncid,'nam','check_randomization',nam%check_randomization)
call mpl%write(lncid,'nam','check_consistency',nam%check_consistency)
call mpl%write(lncid,'nam','check_optimality',nam%check_optimality)
call mpl%write(lncid,'nam','check_obsop',nam%check_obsop)
call mpl%write(lncid,'nam','check_no_obs',nam%check_no_obs)
call mpl%write(lncid,'nam','check_no_point',nam%check_no_point)
call mpl%write(lncid,'nam','check_no_point_mask',nam%check_no_point_mask)
call mpl%write(lncid,'nam','check_no_point_nicas',nam%check_no_point_nicas)
call mpl%write(lncid,'nam','check_set_param_cor',nam%check_set_param_cor)
call mpl%write(lncid,'nam','check_set_param_hyb',nam%check_set_param_hyb)
call mpl%write(lncid,'nam','check_set_param_lct',nam%check_set_param_lct)
call mpl%write(lncid,'nam','check_get_param_stddev',nam%check_get_param_stddev)
call mpl%write(lncid,'nam','check_get_param_cor',nam%check_get_param_cor)
call mpl%write(lncid,'nam','check_get_param_hyb',nam%check_get_param_hyb)
call mpl%write(lncid,'nam','check_get_param_Dloc',nam%check_get_param_Dloc)
call mpl%write(lncid,'nam','check_get_param_lct',nam%check_get_param_lct)
call mpl%write(lncid,'nam','check_apply_vbal',nam%check_apply_vbal)
call mpl%write(lncid,'nam','check_apply_stddev',nam%check_apply_stddev)
call mpl%write(lncid,'nam','check_apply_nicas',nam%check_apply_nicas)
call mpl%write(lncid,'nam','check_apply_obsop',nam%check_apply_obsop)

! model_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','Model parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','nl',nam%nl)
call mpl%write(lncid,'nam','levs',nam%nl,nam%levs(1:nam%nl))
call mpl%write(lncid,'nam','lev2d',nam%lev2d)
call mpl%write(lncid,'nam','logpres',nam%logpres)
call mpl%write(lncid,'nam','nv',nam%nv)
call mpl%write(lncid,'nam','variables',nam%nv,nam%variables(1:nam%nv))
call mpl%write(lncid,'nam','nts',nam%nts)
call mpl%write(lncid,'nam','timeslots',nam%nts,nam%timeslots(1:nam%nts))
call mpl%write(lncid,'nam','dts',nam%dts)
call mpl%write(lncid,'nam','nomask',nam%nomask)
call mpl%write(lncid,'nam','wind_filename',nam%wind_filename)
call mpl%write(lncid,'nam','wind_variables',2,nam%wind_variables(1:2))

! ens1_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','Ensemble 1 parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','ens1_ne',nam%ens1_ne)
call mpl%write(lncid,'nam','ens1_nsub',nam%ens1_nsub)

! ens2_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','Ensemble 2 parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','ens2_ne',nam%ens2_ne)
call mpl%write(lncid,'nam','ens2_nsub',nam%ens2_nsub)

! sampling_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','Sampling parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','sam_write',nam%sam_write)
call mpl%write(lncid,'nam','sam_write_grids',nam%sam_write_grids)
call mpl%write(lncid,'nam','sam_read',nam%sam_read)
call mpl%write(lncid,'nam','mask_type',nam%mask_type)
call mpl%write(lncid,'nam','mask_lu',nam%nv,nam%mask_lu(1:nam%nv))
call mpl%write(lncid,'nam','mask_th',nam%nv,nam%mask_th(1:nam%nv))
call mpl%write(lncid,'nam','ncontig_th',nam%ncontig_th)
call mpl%write(lncid,'nam','mask_check',nam%mask_check)
call mpl%write(lncid,'nam','draw_type',nam%draw_type)
call mpl%write(lncid,'nam','Lcoast',nam%Lcoast*req)
call mpl%write(lncid,'nam','rcoast',nam%rcoast)
call mpl%write(lncid,'nam','nc1',nam%nc1)
call mpl%write(lncid,'nam','nc2',nam%nc2)
call mpl%write(lncid,'nam','ntry',nam%ntry)
call mpl%write(lncid,'nam','nrep',nam%nrep)
call mpl%write(lncid,'nam','nc3',nam%nc3)
call mpl%write(lncid,'nam','dc',nam%dc*req)
call mpl%write(lncid,'nam','nl0r',nam%nl0r)
call mpl%write(lncid,'nam','irmax',nam%irmax)

! diag_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','Diagnostics parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','ne',nam%ne)
call mpl%write(lncid,'nam','gen_kurt_th',nam%gen_kurt_th)
call mpl%write(lncid,'nam','gau_approx',nam%gau_approx)
call mpl%write(lncid,'nam','avg_nbins',nam%avg_nbins)
call mpl%write(lncid,'nam','vbal_block',nam%nv*(nam%nv-1)/2,nam%vbal_block(1:nam%nv*(nam%nv-1)/2))
call mpl%write(lncid,'nam','vbal_rad',nam%vbal_rad)
call mpl%write(lncid,'nam','vbal_dlat',nam%vbal_dlat*rad2deg)
call mpl%write(lncid,'nam','vbal_diag_auto',nam%nv*(nam%nv-1)/2,nam%vbal_diag_auto(1:nam%nv*(nam%nv-1)/2))
call mpl%write(lncid,'nam','vbal_diag_reg',nam%nv*(nam%nv-1)/2,nam%vbal_diag_reg(1:nam%nv*(nam%nv-1)/2))
call mpl%write(lncid,'nam','var_filter',nam%var_filter)
call mpl%write(lncid,'nam','var_niter',nam%var_niter)
call mpl%write(lncid,'nam','var_rhflt',nam%var_rhflt*req)
call mpl%write(lncid,'nam','local_diag',nam%local_diag)
call mpl%write(lncid,'nam','local_rad',nam%local_rad*req)
call mpl%write(lncid,'nam','adv_diag',nam%adv_diag)
call mpl%write(lncid,'nam','adv_type',nam%adv_type)
call mpl%write(lncid,'nam','adv_rad',nam%adv_rad*req)
call mpl%write(lncid,'nam','adv_niter',nam%adv_niter)
call mpl%write(lncid,'nam','adv_rhflt',nam%adv_rhflt*req)
call mpl%write(lncid,'nam','adv_valid',nam%adv_valid)

! fit_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','Fit parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','minim_algo',nam%minim_algo)
call mpl%write(lncid,'nam','diag_rhflt',nam%diag_rhflt*req)
call mpl%write(lncid,'nam','diag_rvflt',nam%diag_rvflt)
call mpl%write(lncid,'nam','smoothness_penalty',nam%smoothness_penalty)
call mpl%write(lncid,'nam','fit_dl0',nam%fit_dl0)
call mpl%write(lncid,'nam','lct_nscales',nam%lct_nscales)
call mpl%write(lncid,'nam','lct_scale_ratio',nam%lct_scale_ratio)
call mpl%write(lncid,'nam','lct_cor_min',nam%lct_cor_min)
call mpl%write(lncid,'nam','lct_diag',nam%lct_nscales,nam%lct_diag(1:nam%lct_nscales))
call mpl%write(lncid,'nam','lct_qc_th',nam%lct_qc_th)
call mpl%write(lncid,'nam','lct_qc_max',nam%lct_qc_max)
call mpl%write(lncid,'nam','lct_write_cor',nam%lct_write_cor)

! nicas_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','NICAS parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','nonunit_diag',nam%nonunit_diag)
call mpl%write(lncid,'nam','lsqrt',nam%lsqrt)
call mpl%write(lncid,'nam','resol',nam%resol)
call mpl%write(lncid,'nam','nc1max',nam%nc1max)
call mpl%write(lncid,'nam','fast_sampling',nam%fast_sampling)
call mpl%write(lncid,'nam','subsamp',nam%subsamp)
call mpl%write(lncid,'nam','network',nam%network)
call mpl%write(lncid,'nam','mpicom',nam%mpicom)
call mpl%write(lncid,'nam','adv_mode',nam%adv_mode)
call mpl%write(lncid,'nam','forced_radii',nam%forced_radii)
call mpl%write(lncid,'nam','rh',nam%rh)
call mpl%write(lncid,'nam','rv',nam%rv)
call mpl%write(lncid,'nam','pos_def_test',nam%pos_def_test)
call mpl%write(lncid,'nam','write_grids',nam%write_grids)

! dirac_param
call mpl%write(lncid,'nam','ndir',nam%ndir)
allocate(londir(nam%ndir))
allocate(latdir(nam%ndir))
if (nam%ndir>0) then
   londir = nam%londir(1:nam%ndir)*rad2deg
   latdir = nam%latdir(1:nam%ndir)*rad2deg
end if
call mpl%write(lncid,'nam','londir',nam%ndir,londir)
call mpl%write(lncid,'nam','latdir',nam%ndir,latdir)
call mpl%write(lncid,'nam','levdir',nam%ndir,nam%levdir(1:nam%ndir))
call mpl%write(lncid,'nam','ivdir',nam%ndir,nam%ivdir(1:nam%ndir))
call mpl%write(lncid,'nam','itsdir',nam%ndir,nam%itsdir(1:nam%ndir))

! obsop_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','Observation operator parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','nobs',nam%nobs)

! output_param
if (mpl%msv%is(lncid)) then
   write(mpl%info,'(a7,a)') '','Output parameters'
   call mpl%flush
end if
call mpl%write(lncid,'nam','nldwv',nam%nldwv)
call mpl%write(lncid,'nam','img_ldwv',nam%nldwv,nam%img_ldwv(1:nam%nldwv))
allocate(lon_ldwv(nam%nldwv))
allocate(lat_ldwv(nam%nldwv))
if (nam%nldwv>0) then
   lon_ldwv = nam%lon_ldwv(1:nam%nldwv)*rad2deg
   lat_ldwv = nam%lat_ldwv(1:nam%nldwv)*rad2deg
end if
call mpl%write(lncid,'nam','lon_ldwv',nam%nldwv,lon_ldwv)
call mpl%write(lncid,'nam','lat_ldwv',nam%nldwv,lat_ldwv)
call mpl%write(lncid,'nam','name_ldwv',nam%nldwv,nam%name_ldwv(1:nam%nldwv))
call mpl%write(lncid,'nam','grid_output',nam%grid_output)
call mpl%write(lncid,'nam','grid_resol',nam%grid_resol*req)

! Release memory
deallocate(londir)
deallocate(latdir)
deallocate(lon_ldwv)
deallocate(lat_ldwv)

end subroutine nam_write

end module type_nam
