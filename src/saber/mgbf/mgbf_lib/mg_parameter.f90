module mg_parameter
!$$$  submodule documentation block
!                .      .    .                                       .
! module:   mg_parameter
!   prgmmr: rancic           org: NCEP/EMC            date: 2022
!
! abstract:  Set resolution, grid and decomposition (offset version)
!
! module history log:
!   2023-04-19  lei     - object-oriented coding
!   2024-01-11  rancic  - optimization for ensemble localization
!   2024-02-20  yokota  - refactoring to apply for GSI
!
! Subroutines Included:
!   init_mg_parameter -
!   def_maxgen -
!   def_ngens -
!
! Functions Included:
!
! remarks:
!   ixm(1)=nxm, jym(1)=nym
!   If mod(nxm,2)=0 then mod(im0,2)=0
!   If mod(nxm,2)>0 then mod(im0,8)=0 (for 4 generations)
!   (This will keep the right boundary of all decompmisitions
!   at same physical location)
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use mgbf_kinds, only: i_kind,r_kind
use mgbf_utils, only: mgbf_abort
use jp_pietc, only: u1
use phint1, only: make_ssgrid, zsigtossig
use mpi

implicit none
public
integer(i_kind),parameter :: lm_max=200
type::  mg_parameter_type
  !-----------------------------------------------------------------------
  !***
  logical:: l_for_localization=.false.  !used for localizaiton while multiple variates need additional treeatment
  logical:: l_mgbf_inhomogeneous=.false.  !used inhomogeneous mgbf
  !*** Namelist parameters
  !***
  real(r_kind):: mg_ampl01,mg_ampl02,mg_ampl03
  real(r_kind):: mg_weig1,mg_weig2,mg_weig3,mg_weig4
  ! avoid a global version of it to avoid memory usage
  integer(i_kind):: mgbf_proc   !1-2: 3D filter                  (1: radial, 2: line)
  !3-5: 2D filter for static B     (3: radial, 4: line, 5: isotropic line)
  !6-8: 2D filter for localization (6: radial, 7: line, 8: isotropic line)
  logical:: mgbf_line
  integer(i_kind):: nxPE,nyPE,im_filt,jm_filt
  logical:: lquart,lhelm

  !***
  !*** Number of generations
  !***
  integer(i_kind):: gm
  integer(i_kind):: gm_max   !clt should be removed?

  !***
  !*** Horizontal resolution
  !***

  !
  ! Original number of data on GSI analysis grid
  !
  integer(i_kind):: nA_max0
  integer(i_kind):: mA_max0

  !
  ! Global number of data on Analysis grid
  !
  integer(i_kind):: nm0
  integer(i_kind):: mm0

  !
  ! Number of PEs on Analysis grid
  !
  integer(i_kind):: nxm
  integer(i_kind):: nym

  !
  ! Number of data on local Analysis grid
  !
  integer(i_kind):: nm
  integer(i_kind):: mm

  !
  ! Number of data on global Filter grid
  !
  integer(i_kind):: im00
  integer(i_kind):: jm00

  !
  ! Number of data on local  Filter grid
  !
  integer(i_kind):: im
  integer(i_kind):: jm

  !
  ! Initial index on local  Filter grid
  !
  integer(i_kind):: i0
  integer(i_kind):: j0
  !
  ! Initial index on local analysis grid
  !
  integer(i_kind):: n0
  integer(i_kind):: m0

  !
  ! Halo on local Filter grid
  !
  integer(i_kind):: ib
  integer(i_kind):: jb

  !
  ! Halo on local Analysis grid
  !
  integer(i_kind):: nb
  integer(i_kind):: mb

  integer(i_kind):: hx,hy,hz
  integer(i_kind):: p
  integer(i_kind):: nh,nfil
  real(r_kind):: pasp01,pasp02,pasp03
  real(r_kind):: pee2,rmom2_1,rmom2_2,rmom2_3,rmom2_4

  integer, allocatable, dimension(:):: maxpe_fgen
  integer, allocatable, dimension(:):: ixm,jym,nxy
  integer, allocatable, dimension(:):: im0,jm0
  integer, allocatable, dimension(:):: Fimax,Fjmax
  integer, allocatable, dimension(:):: FimaxL,FjmaxL
  real(r_kind), allocatable, dimension(:):: zofis  ! index of s(fitering grids) in analysis grids (its index is its coor)
  real(r_kind), allocatable, dimension(:):: isofz  ! index of z of analysis grids in the filtering grids

  integer(i_kind):: npes_filt
  integer(i_kind):: maxpe_filt

  integer(i_kind):: imL,jmL
  integer(i_kind):: imH,jmH
  integer(i_kind):: lm_a          ! number of vertical layers in analysis fields
  integer(i_kind):: lm            ! number of vertical layers in filter grids
  real(r_kind):: coef_normalization(lm_max)=1 !normalizaton coefficients
  real(r_kind):: coef_normalization_const=-9999.0 ! constant, if set, this contant will be
  ! assigned to all elements of coef_normalization

  integer(i_kind):: km2           ! number of 2d variables for filtering
  integer(i_kind):: km3           ! number of 3d variables for filtering
  integer(i_kind):: n_ens         ! number of ensemble members
  integer(i_kind):: km_a          ! total number of horizontal levels for analysis
  integer(i_kind):: km_all        ! total number of k levels of ensemble for filtering
  integer(i_kind):: km_a_all      ! total number of k levels of ensemble
  integer(i_kind):: km2_all       ! total number of k horizontal levels of ensemble for filtering
  integer(i_kind):: km3_all       ! total number of k vertical levels of ensemble
  logical :: l_loc=.false.               ! logical flag for localization
  logical :: l_filt_g1=.true.            ! logical flag for filtering of generation one
  logical :: l_lin_vertical=.true.       ! logical flag for linear interpolation in vertcial
  logical :: l_lin_horizontal=.true.     ! logical flag for linear interpolation in horizontal
  logical :: l_quad_horizontal=.false.    ! logical flag for quadratic interpolation in horizontal
  logical :: l_new_map            ! logical flag for new mapping between analysis and filter grid
  logical :: l_vertical_filter    ! logical flag for vertical filtering
  logical :: l_vert_stretched_filtgrid=.false.  ! true : filtering grids are stretched in tems of analysis grid unit
  logical :: l_anal_sub_of_filt   ! true : analysis grids and filtering grids are the same excpet for later has boundary points
  integer(i_kind):: km            ! number of vertically stacked all variables (km=km2+lm*km3)
  integer(i_kind):: km_4
  integer(i_kind):: km_16
  integer(i_kind):: km_64

  real(r_kind):: lengthx,lengthy,xa0,ya0,xf0,yf0
  real(r_kind):: dxf,dyf,dxa,dya
  real(r_kind),allocatable,dimension (:,:):: dxfm,dyfm  ! actual filtering grid intervals in meters
  ! The control filtering grid intervals corresponding to the constant horizontal aspect tensor.
  real(r_kind):: dxfmctrl=35000,dyfmctrl=35000
  real(r_kind):: dx_a2f_ratio=1,dy_a2f_ratio=1  !ratio between analsysis grids to filtering grids in x and y
  !it will be derived from other namelist parameters
  logical :: l_constant_aspt2 =.true. ! using constant horizontal aspect tensor : ampl02

  integer(i_kind):: npadx         ! x padding on analysis grid
  integer(i_kind):: mpady         ! y padding on analysis grid

  integer(i_kind):: ipadx         ! x padding on filter decomposition
  integer(i_kind):: jpady         ! y padding on filter deocmposition

  !
  ! Just for standalone test
  !
  logical:: ldelta

  !from mg_mppstuff.f90
  character(len=5):: c_mype
  integer(i_kind):: mype
  ! dTYPE is the MPI datatype used by the halo exchanges in mg_bocos.
  integer(i_kind):: npes,dTYPE,ierr,ierror
  ! mpi_comm_comp is the component communicator, borrowed from the covariance
  ! object: this type never duplicates or frees it.  The group and filter handles
  ! below ARE owned here, created in init_mg_MPI and released in mg_finalize.
  integer(i_kind):: mpi_comm_comp=MPI_COMM_NULL
  integer(i_kind):: filter_comm=MPI_COMM_NULL
  integer(i_kind):: component_group=MPI_GROUP_NULL
  integer(i_kind):: filter_group=MPI_GROUP_NULL
  integer(i_kind):: mype_gr,npes_gr
  integer(i_kind):: my_hgen
  integer(i_kind):: mype_hgen
  logical:: l_hgen
  integer(i_kind):: nx,my
  !from mg_domain.f90
  logical,dimension(2):: Flwest,Fleast,Flnorth,Flsouth
  integer(i_kind),dimension(2):: Fitarg_n,Fitarg_e,Fitarg_s,Fitarg_w
  integer(i_kind),dimension(2):: Fitarg_sw,Fitarg_se,Fitarg_ne,Fitarg_nw
  logical,dimension(2):: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne
  integer(i_kind),dimension(2):: Fitarg_up
  integer(i_kind):: itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw
  integer(i_kind):: itarg_wA,itarg_eA,itarg_sA,itarg_nA
  logical:: lwestA,leastA,lsouthA,lnorthA
  integer(i_kind):: ix,jy
  integer(i_kind),dimension(2):: mype_filt
  !from mg_domain_loc.f90
  integer(i_kind):: nsq21,nsq32,nsq43
  logical,dimension(4):: Flsouth_loc,Flnorth_loc,Flwest_loc,Fleast_loc
  integer(i_kind),dimension(4):: Fitarg_s_loc,Fitarg_n_loc,Fitarg_w_loc,Fitarg_e_loc
  integer(i_kind),dimension(4):: Fitargup_loc12
  integer(i_kind),dimension(4):: Fitargup_loc23
  integer(i_kind),dimension(4):: Fitargup_loc34
  integer(i_kind):: itargdn_sw_loc21,itargdn_se_loc21,itargdn_nw_loc21,itargdn_ne_loc21
  integer(i_kind):: itargdn_sw_loc32,itargdn_se_loc32,itargdn_nw_loc32,itargdn_ne_loc32
  integer(i_kind):: itargdn_sw_loc43,itargdn_se_loc43,itargdn_nw_loc43,itargdn_ne_loc43
  logical:: lsendup_sw_loc,lsendup_se_loc,lsendup_nw_loc,lsendup_ne_loc
  logical:: l_mg_weig_readin=.false.
  !clt for use of resolution-varied vertical filtering grids
  real(r_kind), allocatable,dimension(:):: aspect_vert_profile_angrid ! should be of size (lm)
  real(r_kind), allocatable,dimension(:):: aspect_vert_profile_filtgrid ! should be of size (lm)

  contains
    procedure :: init_mg_parameter
    !from mg_mppstuff.f90
    procedure :: init_mg_MPI
    procedure :: finishMPI
    procedure :: barrierMPI
    !from mg_domain.f90
    procedure :: init_mg_domain
    procedure :: init_domain
    procedure :: init_topology_2d
    procedure :: real_itarg
    !from mg_domain_loc.f90
    procedure :: init_domain_loc
    procedure :: sidesend_loc
    procedure :: targup_loc
    procedure :: targdn21_loc
    procedure :: targdn32_loc
    procedure :: targdn43_loc
    !from jp_pbfil.f90
    generic :: cholaspect => cholaspect1,cholaspect2,cholaspect3,cholaspect4
    procedure,nopass :: cholaspect1,cholaspect2,cholaspect3,cholaspect4
    generic :: getlinesum => getlinesum1,getlinesum1d,getlinesum2,getlinesum3
    procedure :: getlinesum1,getlinesum1d,getlinesum2,getlinesum3
    generic :: rbeta => rbeta1,rbeta3d_1,rbeta2,rbeta3,rbeta4,vrbeta1,vrbeta2,vrbeta3,vrbeta4
    procedure:: rbeta1,rbeta3d_1,rbeta2,rbeta3,rbeta4,vrbeta1,vrbeta2,vrbeta3,vrbeta4
    generic :: rbetaT => rbeta1t,rbeta3d_1t,rbeta2t,rbeta3t,rbeta4t,vrbeta1t,vrbeta2t,vrbeta3t,vrbeta4t
    procedure:: rbeta1t,rbeta3d_1t,rbeta2t,rbeta3t,rbeta4t,vrbeta1t,vrbeta2t,vrbeta3t,vrbeta4t
end type  mg_parameter_type

interface
  !from mg_mppstuff.f90
  module subroutine init_mg_MPI(this)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
  end subroutine init_mg_MPI
  module subroutine finishMPI(this)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
  end subroutine finishMPI
  module subroutine barrierMPI(this)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
  end subroutine barrierMPI
  !from mg_domain.f90
  module subroutine init_mg_domain(this)
    implicit none
    class(mg_parameter_type), intent(inout) :: this
  end subroutine init_mg_domain
  module subroutine init_domain(this)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
  end subroutine init_domain
  module subroutine init_topology_2d(this)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
  end subroutine init_topology_2d
  module subroutine real_itarg (this,itarg)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
    integer(i_kind), intent(inout):: itarg
  end subroutine real_itarg
  !from mg_domain_loc.f90
  module subroutine init_domain_loc(this)
    implicit none
    class(mg_parameter_type), intent(inout) :: this
  end subroutine init_domain_loc
  module subroutine sidesend_loc(this)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
  end subroutine sidesend_loc
  module subroutine targup_loc(this)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
  end subroutine targup_loc
  module subroutine targdn21_loc(this)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
  end subroutine targdn21_loc
  module subroutine targdn32_loc(this)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
  end subroutine targdn32_loc
  module subroutine targdn43_loc(this)
    implicit none
    class(mg_parameter_type), intent(inout), target :: this
  end subroutine targdn43_loc
  !from jp_pbfil.f90
  module subroutine cholaspect1(lx,mx, el)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    integer,                      intent(in   ):: lx,mx
    real(dp),dimension(1,1,lx:mx),intent(inout):: el
  end subroutine cholaspect1
  module subroutine cholaspect2(lx,mx, ly,my, el)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    integer,                            intent(in   ):: lx,mx, ly,my
    real(dp),dimension(2,2,lx:mx,ly:my),intent(inout):: el
    real(dp),dimension(2,2):: tel
  end subroutine cholaspect2
  module subroutine cholaspect3(lx,mx, ly,my, lz,mz, el)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    integer,                                  intent(in   ):: lx,mx, ly,my, lz,mz
    real(dp),dimension(3,3,lx:mx,ly:my,lz:mz),intent(inout):: el
    real(dp),dimension(3,3):: tel
  end subroutine cholaspect3
  module subroutine cholaspect4(lx,mx, ly,my, lz,mz, lw,mw,el)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    integer,                                        intent(in   ):: lx,mx, ly,my, lz,mz, lw,mw
    real(dp),dimension(4,4,lx:mx,ly:my,lz:mz,lw:mw),intent(inout):: el
    real(dp),dimension(4,4):: tel
  end subroutine cholaspect4
  module subroutine getlinesum1(this,hx,lx,mx, el, ss)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                      intent(in   ):: hx,Lx,mx
    real(dp),dimension(1,1,Lx:Mx),intent(in   ):: el
    real(dp),dimension(    lx:mx),intent(  out):: ss
  end subroutine getlinesum1
  module subroutine getlinesum1d(this,hx,lx,mx, el, ss)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                      intent(in   ):: hx,Lx,mx
    real(dp),dimension(Lx:Mx),intent(in   ):: el
    real(dp),dimension(    lx:mx),intent(  out):: ss
  end subroutine getlinesum1d
  module subroutine getlinesum2(this,hx,lx,mx, hy,ly,my, el, ss)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                            intent(in   ):: hx,Lx,mx, hy,ly,my
    real(dp),dimension(2,2,Lx:Mx,Ly:My),intent(in   ):: el
    real(dp),dimension(    lx:mx,ly:my),intent(  out):: ss
  end subroutine getlinesum2
  module subroutine getlinesum3(this,hx,lx,mx, hy,ly,my, hz,lz,mz, el, ss)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                                  intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz
    real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: el
    real(dp),dimension(    lx:mx,ly:my,lz:mz),intent(  out):: ss
  end subroutine getlinesum3
  module subroutine getlinesum4(this,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el, ss)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                                        intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw
    real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
    real(dp),dimension(    lx:mx,ly:my,lz:mz,Lw:Mw),intent(  out):: ss
  end subroutine getlinesum4
  module subroutine rbeta1(this,hx,lx,mx, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                  intent(in   ):: hx,Lx,mx
    real(dp),dimension(Lx:Mx),intent(in   ):: el
    real(dp),dimension(Lx:Mx),intent(in   ):: ss
    real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
  end subroutine rbeta1
  module subroutine rbeta3d_1(this,nz,hx,lx,mx, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                  intent(in   ):: nz, hx,Lx,mx
    real(dp),dimension(nz,Lx:Mx),intent(in   ):: el
    real(dp),dimension(nz,Lx:Mx),intent(in   ):: ss
    real(dp),dimension(nz,lx-hx:mx+hx),intent(inout):: a
  end subroutine rbeta3d_1
  module subroutine rbeta2(this,hx,lx,mx, hy,ly,my, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                            intent(in   ):: hx,Lx,mx, hy,ly,my
    real(dp),dimension(2,2,Lx:Mx,Ly:My),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My),intent(in   ):: ss
    real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
  end subroutine rbeta2
  module subroutine rbeta3(this,hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss,a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                                  intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz
    real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: ss
    real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz),intent(inout):: a
  end subroutine rbeta3
  module subroutine rbeta4(this,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el,ss,a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                                        intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw
    real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: ss
    real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz,lw-hw:mw+hw),intent(inout):: a
  end subroutine rbeta4
  module subroutine rbeta1T(this,hx,lx,mx, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                      intent(in   ):: hx,Lx,mx
    real(dp),dimension(1,1,Lx:Mx),intent(in   ):: el
    real(dp),dimension(    Lx:Mx),intent(in   ):: ss
    real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
  end subroutine rbeta1T
  module subroutine rbeta3d_1T(this,nz,hx,lx,mx, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                      intent(in   ):: nz, hx,Lx,mx
    real(dp),dimension(nz,Lx:Mx),intent(in   ):: el
    real(dp),dimension(nz, Lx:Mx),intent(in   ):: ss
    real(dp),dimension(nz,lx-hx:mx+hx),intent(inout):: a
  end subroutine rbeta3d_1T
  module subroutine rbeta2T(this,hx,lx,mx, hy,ly,my, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                            intent(in   ):: hx,Lx,mx, hy,ly,my
    real(dp),dimension(2,2,Lx:Mx,Ly:My),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My),intent(in   ):: ss
    real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
  end subroutine rbeta2T
  module subroutine rbeta3T(this,hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                                  intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz
    real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: ss
    real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz),intent(inout):: a
  end subroutine rbeta3T
  module subroutine rbeta4T(this,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                                        intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw
    real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: ss
    real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz,lw-hw:mw+hw),intent(inout):: a
  end subroutine rbeta4T
  module subroutine vrbeta1(this,nv,hx,lx,mx, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                      intent(in   ):: nv,hx,Lx,mx
    real(dp),dimension(1,1,Lx:Mx),intent(in   ):: el
    real(dp),dimension(    Lx:Mx),intent(in   ):: ss
    real(dp),dimension(nv,lx-hx:mx+hx),intent(inout):: a
  end subroutine vrbeta1
  module subroutine vrbeta2(this,nv,hx,lx,mx, hy,ly,my, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                            intent(in   ):: nv, hx,Lx,mx, hy,ly,my
    real(dp),dimension(2,2,Lx:Mx,Ly:My),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My),intent(in   ):: ss
    real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
  end subroutine vrbeta2
  module subroutine vrbeta3(this,nv, hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss,a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                                  intent(in   ):: nv, hx,Lx,mx, hy,ly,my, hz,lz,mz
    real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: ss
    real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz),intent(inout):: a
  end subroutine vrbeta3
  module subroutine vrbeta4(this,nv,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el,ss,a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                                        intent(in   ):: nv, hx,Lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw
    real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: ss
    real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz,lw-hw:mw+hw),intent(inout):: a
  end subroutine vrbeta4
  module subroutine vrbeta1T(this,nv, hx,lx,mx, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                      intent(in   ):: nv,hx,Lx,mx
    real(dp),dimension(1,1,Lx:Mx),intent(in   ):: el
    real(dp),dimension(    Lx:Mx),intent(in   ):: ss
    real(dp),dimension(nv,lx-hx:mx+hx),intent(inout):: a
  end subroutine vrbeta1T
  module subroutine vrbeta2T(this,nv,hx,lx,mx, hy,ly,my, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                            intent(in   ):: nv, hx,Lx,mx, hy,ly,my
    real(dp),dimension(2,2,Lx:Mx,Ly:My),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My),intent(in   ):: ss
    real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
  end subroutine vrbeta2T
  module subroutine vrbeta3T(this,nv,hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                                  intent(in   ):: nv, hx,Lx,mx, hy,ly,my, hz,lz,mz
    real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: ss
    real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz),intent(inout):: a
  end subroutine vrbeta3T
  module subroutine vrbeta4T(this,nv,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el,ss, a)
    use mgbf_kinds, only: dp=>r_kind
    implicit none
    class(mg_parameter_type), intent(inout) :: this
    integer,                                        intent(in   ):: nv, hx,Lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw
    real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
    real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: ss
    real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz,lw-hw:mw+hw),intent(inout):: a
  end subroutine vrbeta4T
end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine init_mg_parameter(this,inputfilename)
!**********************************************************************!
!                                                                      !
! Initialize ....                                                      !
!                                                                      !
!**********************************************************************!
implicit none
class(mg_parameter_type), intent(inout), target :: this
integer(i_kind):: g
character(*), intent(in) :: inputfilename

! Namelist parameters as local variable
real(r_kind):: mg_ampl01,mg_ampl02,mg_ampl03
real(r_kind):: mg_weig1,mg_weig2,mg_weig3,mg_weig4
integer(i_kind):: mgbf_proc
logical:: mgbf_line
integer(i_kind):: nxPE,nyPE,im_filt,jm_filt
logical, save :: lquart=.false.,lhelm=.false. !clt what should be the default
logical, save :: ldelta=.false.
logical, save :: l_for_localization=.false.
logical, save :: l_mgbf_inhomogeneous=.false.

integer(i_kind):: lm_a          ! number of vertical layers in analysis fields
integer(i_kind):: lm            ! number of vertical layers in filter grids
real(r_kind), save :: coef_normalization(lm_max)=1 !normalizaton coefficients
real(r_kind), save :: coef_normalization_const=-9999.0 ! constant, if set, this contant will be
! The control filtering grid intervals corresponding to the constant horizontal aspect tensor.
real(r_kind), save :: dxfmctrl=35000,dyfmctrl=35000
logical, save :: l_constant_aspt2 =.true. ! using constant horizontal aspect tensor : ampl02
character(len=256), save :: file_coef_normalization="XXXX"
character(len=256), save :: dir_coef_normalization="XXXX"
integer(i_kind):: km2           ! number of 2d variables for filtering
integer(i_kind):: km3           ! number of 3d variables for filtering
integer(i_kind), save :: n_ens=1         ! number of ensemble members
logical, save :: l_loc=.false.
logical, save :: l_filt_g1=.false.            ! logical flag for filtering of generation one
logical, save :: l_lin_vertical=.false.       ! logical flag for linear interpolation in vertcial
logical, save :: l_lin_horizontal=.false.     ! logical flag for linear interpolation in horizontal
logical, save :: l_quad_horizontal=.false.    ! logical flag for quadratic interpolation in horizontal
logical, save :: l_new_map=.false.            ! logical flag for new mapping between analysis and filter grid
logical, save :: l_vertical_filter=.true.    ! logical flag for vertical filtering
logical, save ::  l_anal_sub_of_filt=.false.
logical, save ::  l_vert_stretched_filtgrid=.false.
integer(i_kind), save :: gm_max=4   !clt by defaul

! Global number of data on Analysis grid
integer(i_kind):: nm0
integer(i_kind):: mm0

integer(i_kind):: hx,hy,hz
integer(i_kind):: p
logical, save :: l_mg_weig_readin=.false.
integer(i_kind), parameter       :: nf=20! refinement factor for z grid,used in make_ssgrid
integer(i_kind) :: myunit,i,item,mype,ierr
character(len=4) :: str_rank
integer :: n_sample_levelsx4normalization
logical :: l_exist
character(len=256) :: abort_msg
namelist /parameters_mgbeta/ mg_ampl01,mg_ampl02,mg_ampl03            &
    ,mg_weig1,mg_weig2,mg_weig3,mg_weig4      &
    ,hx,hy,hz,p                               &
    ,mgbf_line,mgbf_proc                      &
    ,lm_a,lm,coef_normalization               &
    ,coef_normalization_const &
    ,dir_coef_normalization  &
    ,file_coef_normalization  &
    , dxfmctrl,dyfmctrl       &
    , l_constant_aspt2        &
    ,km2,km3                                  &
    ,n_ens                                    &
    ,l_loc                                    &
    ,l_filt_g1                                &
    ,l_lin_vertical                           &
    ,l_lin_horizontal                         &
    ,l_quad_horizontal                        &
    ,l_new_map                                &
    ,l_vertical_filter                        &
    ,l_anal_sub_of_filt                       &
    ,l_vert_stretched_filtgrid                     &
    ,l_for_localization,ldelta,lquart,lhelm   &
    , l_mgbf_inhomogeneous                    &
    ,gm_max                                   &
    ,nm0,mm0                                  &
    ,nxPE,nyPE,im_filt,jm_filt ,              &
    l_mg_weig_readin

open(unit=10,file=trim(inputfilename),status="old",action="read")
read(10,nml=parameters_mgbeta)
close(unit=10)
!
allocate(this%zofis(lm))
allocate(this%isofz(lm_a))
this%l_vert_stretched_filtgrid=l_vert_stretched_filtgrid

if(lm_a /= lm ) then
  if (this%mype == 0) write(6,*) "MGBF l_vert_stretched_filtgrid ", this%l_vert_stretched_filtgrid
  call convert_vert_varied_aspt
  !in which the mg_ampl01 will be re-defined
end if
!-----------------------------------------------------------------
!for safety, copy all namelist loc vars to them of this object
this%mg_ampl01=mg_ampl01
this%mg_ampl02=mg_ampl02
this%mg_ampl03=mg_ampl03
this%mg_weig1=mg_weig1
this%mg_weig2=mg_weig2
this%mg_weig3=mg_weig3
this%mg_weig4=mg_weig4
this%hx=hx
this%hy=hy
this%hz=hz
this%p =p
this%mgbf_line=mgbf_line
this%mgbf_proc=mgbf_proc
this%lm_a=lm_a
this%lm=lm
if (coef_normalization_const >0 ) then  ! constant, if set, this contant will be

  if(trim(file_coef_normalization)=="XXXX" .and. trim(dir_coef_normalization)=="XXXX" ) then
    l_exist=.false.
    coef_normalization=coef_normalization_const
  else
    if (trim(dir_coef_normalization) /= "XXXX") then
      ! The profile is per subdomain, so it must be selected by the rank within this
      ! component.  Using the world rank would make a second component read the first
      ! component's profiles and filter with the wrong normalization.
      mype=this%mype
      write(str_rank, "(I4.4)") mype
      file_coef_normalization=trim(dir_coef_normalization)//"/profile_subdomain_"//str_rank//".txt"
    end if
    if (this%mype == 0) write(6,*) "MGBF normalization file is ", trim(file_coef_normalization)
    inquire(file=trim(file_coef_normalization),exist=l_exist)
    if(l_exist) then
      open(newunit=myunit,file=trim(file_coef_normalization),status="old",action="read")
      if(trim(dir_coef_normalization) /= "XXXX") then
        ! to use file slike profiles_out/profile_subdomain_0475.txt
        read(myunit,*)
        read(myunit,*)i,n_sample_levelsx4normalization
        read(myunit,*)
        do i=1,n_sample_levelsx4normalization
          read(myunit,*)
        end do
        read(myunit,*)
        do i=1,lm_a
          read(myunit,*)item, coef_normalization(i) !notice, the data in the file is reversed already
        end do
        close (myunit)
      else
        if (this%mype == 0) write(6,*) "the normalization profile file is ", trim(file_coef_normalization)
        !clt in the ../covairance/mgbf_covariance_mod.f90 the fldset is reversed in the vertical direction
        open(newunit=myunit,file=trim(file_coef_normalization),status="old",action="read")
        read(myunit,*)(coef_normalization(i),i=lm_a,1,-1)
        close (myunit)
      end if
      coef_normalization(1:lm_a)=coef_normalization(1:lm_a)*coef_normalization_const  !re-calc
    else

      call mgbf_abort("normalization profile file does not exist: " // &
          trim(file_coef_normalization))
    end if
  end if
else
  coef_normalization=1.0


end if




this%coef_normalization=coef_normalization
this%dxfmctrl=dxfmctrl; this%dyfmctrl=dyfmctrl
this%l_constant_aspt2 = l_constant_aspt2
this%km2=km2
this%km3=km3
this%n_ens=n_ens
this%l_loc=l_loc
this%l_filt_g1=l_filt_g1
this%l_lin_vertical=l_lin_vertical
this%l_lin_horizontal=l_lin_horizontal
this%l_quad_horizontal=l_quad_horizontal
this%l_new_map=l_new_map
this%l_vertical_filter=l_vertical_filter
this%l_anal_sub_of_filt=l_anal_sub_of_filt
this%l_for_localization=l_for_localization
this%l_mgbf_inhomogeneous = l_mgbf_inhomogeneous
this%ldelta=ldelta
this%lquart=lquart
this%lhelm=lhelm
this%nm0=nm0
this%mm0=mm0
this%nxPE=nxPE
this%nyPE=nyPE
this%im_filt=im_filt
this%jm_filt=jm_filt
this%nxm = nxPE
this%nym = nyPE

this%im = im_filt
this%jm = jm_filt

!-----------------------------------------------------------------
!
!
! For 168 PES
!
!
! For 256 PES
!
!
! For 336 PES
!
!
! For 448 PES
!
!
! For 512 PES
!
!
! For 704 PES
!
!
! For 768 PES
!
!
! For 924 PES
!
!
! For 1056 PES
!
!
! For 1408 PES
!
!
! For 1848 PES
!
!
! For 2464 PES
!

!
! Define total number of horizontal levels in the case of ensemble
!

this%km_a = this%km2+this%lm_a*this%km3
this%km   = this%km2+this%lm  *this%km3

this%km_a_all = this%km_a * this%n_ens
this%km_all   = this%km   * this%n_ens

this%km2_all = this%km2 * this%n_ens
this%km3_all = this%km3 * this%n_ens

this%km_4  = this%km/4
this%km_16 = this%km/16
this%km_64 = this%km/64
this%l_mg_weig_readin=l_mg_weig_readin

!
! Define maximum number of generations 'gm'
!

call def_maxgen(this%nxm,this%nym,this%gm)

! Restrict to gm_max

if(this%gm>gm_max) then
  this%gm=gm_max
end if
if(this%nxm*this%nym<=1) then
  this%gm=gm_max
end if

!***
!***     Analysis grid
!***

!
! Number of grid intervals on GSI grid for the reduced RTMA domain
! before padding
!

!
! Number of grid points on the analysis grid after padding
!
this%nm = this%nm0/this%nxm
this%mm = this%mm0/this%nym
this%dx_a2f_ratio=this%nm/this%im_filt
this%dy_a2f_ratio=this%mm/this%jm_filt
if(this%l_anal_sub_of_filt ) then
  if(this%im_filt/=this%nm.or.this%jm_filt/=this%mm) then
    call mgbf_abort("l_anal_sub_of_filter requires compatible analysis and filtering grids")
  end if
  if(l_lin_horizontal.or.l_quad_horizontal) then
    call mgbf_abort("l_anal_sub_of_filter currently requires lsqr")
  end if
end if

!***
!***     Filter grid
!***


!
! For 168 PES
!
!
! For 256 PES
!
!
!
! For 336 PES
!
!
! For 448 PES
!
!
! For 512 PES
!
!
! For 704 PES
!
!
! For 768 PES
!
!
! For 924 PES
!
!
! For 1056 PES
!
!
! For 1408 PES
!
!
! For 1848 PES
!
!
! For 2464 PES
!

this%im00 = this%nxm*this%im
this%jm00 = this%nym*this%jm

this%n0 = 1
this%m0 = 1

this%i0 = 1
this%j0 = 1

!
! Make sure that nm0 and mm0 and divisibvle with nxm and nym
!
if(this%nm*this%nxm /= this%nm0 ) then
  write(abort_msg,'("nm0 is not divisible by nxm: nm=",I0,", nxm=",I0,", nm0=",I0)') &
      this%nm, this%nxm, this%nm0
  call mgbf_abort(abort_msg)
end if

if(this%mm*this%nym /= this%mm0 ) then
  write(abort_msg,'("mm0 is not divisible by nym: mm=",I0,", nym=",I0,", mm0=",I0)') &
      this%mm, this%nym, this%mm0
  call mgbf_abort(abort_msg)
end if

!
! Set number of processors at higher generations
!

allocate(this%ixm(this%gm))
allocate(this%jym(this%gm))
allocate(this%nxy(this%gm))
allocate(this%maxpe_fgen(0:this%gm))
allocate(this%im0(this%gm))
allocate(this%jm0(this%gm))
allocate(this%Fimax(this%gm))
allocate(this%Fjmax(this%gm))
allocate(this%FimaxL(this%gm))
allocate(this%FjmaxL(this%gm))

call def_ngens(this%ixm,this%gm,this%nxm)
call def_ngens(this%jym,this%gm,this%nym)

!$omp parallel do private(g) schedule(static)
do g=1,this%gm
  this%nxy(g)=this%ixm(g)*this%jym(g)
end do
!$omp end parallel do

this%maxpe_fgen(0)= 0
do g=1,this%gm
  this%maxpe_fgen(g)=this%maxpe_fgen(g-1)+this%nxy(g)
end do

this%maxpe_filt=this%maxpe_fgen(this%gm)
this%npes_filt=this%maxpe_filt-this%nxy(1)

this%im0(1)=this%im00
do g=2,this%gm
  this%im0(g)=this%im0(g-1)/2
end do

this%jm0(1)=this%jm00
do g=2,this%gm
  this%jm0(g)=this%jm0(g-1)/2
end do

!$omp parallel do private(g) schedule(static)
do g=1,this%gm
  this%Fimax(g)=this%im0(g)-this%im*(this%ixm(g)-1)
  this%Fjmax(g)=this%jm0(g)-this%jm*(this%jym(g)-1)
end do
!$omp end parallel do

!$omp parallel do private(g) schedule(static)
do g=1,this%gm
  this%FimaxL(g)=this%Fimax(g)/2
  this%FjmaxL(g)=this%Fjmax(g)/2
end do
!$omp end parallel do

!***
!*** Filter related parameters
!**
this%lengthx = 1.*this%nm     ! arbitrary chosen scale of the domain
this%lengthy = 1.*this%mm     ! arbitrary chosen scale of the domain

this%ib=6
this%jb=6

this%dxa =this%lengthx/this%nm
this%dxf = this%lengthx/this%im
allocate(this%dxfm(this%im,this%jm))
this%nb = 2*this%dxf/this%dxa

this%dya = this%lengthy/this%mm
this%dyf = this%lengthy/this%jm
allocate(this%dyfm(this%im,this%jm))
this%mb = 2*this%dyf/this%dya

this%xa0 = this%dxa*0.5
this%ya0 = this%dya*0.5

this%xf0 = this%dxf*0.5
this%yf0 = this%dyf*0.5

this%imL=this%im/2
this%jmL=this%jm/2

this%imH=this%im0(this%gm)
this%jmH=this%jm0(this%gm)
this%pasp01 = mg_ampl01
this%pasp02 = mg_ampl02
this%pasp03 = mg_ampl03

this%nh= max(hx,hy,hz)
this%nfil = this%nh + 2

this%pee2=this%p*2
this%rmom2_1=u1/sqrt(this%pee2+3)
this%rmom2_2=u1/sqrt(this%pee2+4)
this%rmom2_3=u1/sqrt(this%pee2+5)
this%rmom2_4=u1/sqrt(this%pee2+6)
contains

subroutine convert_vert_varied_aspt

integer(i_kind) :: lm_tmp,iz,is,mype,ierr
real(r_kind):: sstop,dss
real (r_kind),allocatable,dimension(:):: sigofz
real (r_kind),allocatable,dimension(:):: sigofis
integer(i_kind):: user_mpi_real
real (r_kind) :: mg_ampl01_org

if( .not. allocated(this%aspect_vert_profile_angrid )) then
  allocate(this%aspect_vert_profile_angrid(lm_a),this%aspect_vert_profile_filtgrid(lm))
end if
allocate(sigofz(lm_a),sigofis(lm))
! Component rank, not world rank: this routine runs once per intstate, and the
! broadcast below is collective over the component communicator only.
mype=this%mype
if(this%l_vert_stretched_filtgrid) then
  if(mype==0) then
    open(newunit=myunit,file="mgbf_vert_aspt_profile.txt",status="old",iostat=ierr)
    if (ierr /= 0) call mgbf_abort("cannot open mgbf_vert_aspt_profile.txt")
    read(myunit,*)lm_tmp
    if(lm_tmp/=lm_a) then
      call mgbf_abort("lm_a differs from the size in mgbf_vert_aspt_profile.txt")
    end if
    do i=1,lm_a
      read(myunit,*)this%aspect_vert_profile_angrid(i)
    end do
    close(myunit)
  end if
  call MPI_Type_match_size(MPI_TYPECLASS_REAL, kind(this%aspect_vert_profile_angrid(1)), user_mpi_real, ierr)
  if (ierr /= MPI_SUCCESS) then
    write(abort_msg,'("no matching MPI type for real kind =",I0)') &
        kind(this%aspect_vert_profile_angrid(1))
    call mgbf_abort(abort_msg)
  end if
  call MPI_Bcast(this%aspect_vert_profile_angrid, lm_a, user_mpi_real, 0, this%mpi_comm_comp, ierr)
  if (ierr /= MPI_SUCCESS) then
    call mgbf_abort("MPI_Bcast of the vertical aspect profile failed")
  end if


  ! calibrate sigscale to make sigofz go to sigbottom at z=0:
  sigofz=sqrt(this%aspect_vert_profile_angrid)
  if(mype==0) then
    do iz=lm_a,1,-1
      write(6,*) iz, sigofz(iz)
    end do
  end if
else
  sigofz=sqrt(mg_ampl01)

end if

! Make the new grid whose resolution of the correlation scale sigofz
! is uniform throughout.
! isofz is the s-index coordinate of each of the original z-grid points.
! zofis is the z-index coordinate of each of the new s-grid points.
call make_ssgrid(lm_a-1,nf,lm-1,sigofz, sstop,dss,this%isofz,this%zofis)

! Use the new s-grid locations zofis, and the original profile of
! correlation scales sigofz, to interpolate, smoothly and positively,
! these scales sig to each of the new s-grid points:
call zsigtossig(lm_a-1,nf,lm-1,this%zofis,sigofz,sigofis)
mg_ampl01_org=mg_ampl01
mg_ampl01=(sum(sigofis**2)/size(sigofis))
if(.not.this%l_vert_stretched_filtgrid) then !the former could be only true when the latter is in effect
  if (this%mype == 0) write(6,*) " suggested and actual/original ampl01 is ", &
      mg_ampl01, " ", mg_ampl01_org
  mg_ampl01=mg_ampl01_org
end if
if (this%mype == 0) write(6,*) " the original and final ampl01 is ", &
    mg_ampl01_org, " ", mg_ampl01

do is=1,lm
  if (this%mype == 0) write(6,*) is, this%zofis(is), (sigofis(is))**2
end do
!clttodo need to access this from mgbf lib too
this%zofis=this%zofis(lm:1:-1)




deallocate(sigofz,sigofis)
end subroutine convert_vert_varied_aspt



!----------------------------------------------------------------------
end subroutine init_mg_parameter

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine def_maxgen &
!**********************************************************************
!                                                                     !
!  Given number of PEs in x and y direction decides what is the       !
!  maximum number of generations that a multigrid scheme can support  !
!                                                                     !
!                                                    M. Rancic 2020   !
!**********************************************************************
    (nxm,nym,gm)
!----------------------------------------------------------------------
implicit none
integer, intent(in):: nxm,nym
integer, intent(out):: gm
integer:: npx,npy,gx,gy

npx = nxm;  gx=1
Do
  npx = (npx + 1)/2
  gx = gx + 1
  if(npx == 1) exit
end do

npy = nym;  gy=1
Do
  npy = (npy + 1)/2
  gy = gy + 1
  if(npy == 1) exit
end do

gm = Min(gx,gy)


!----------------------------------------------------------------------
end subroutine def_maxgen

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine def_ngens &
!*********************************************************************!
!                                                                     !
!  Given number of generations, find number of PEs is s direction     !
!                                                                     !
!                                                    M. Rancic 2020   !
!*********************************************************************!
    (nsm,gm,nsm0)
!----------------------------------------------------------------------
implicit none
integer, intent(in):: gm,nsm0
integer, dimension(gm), intent(out):: nsm
integer:: g
!----------------------------------------------------------------------

nsm(1)=nsm0
Do g=2,gm
  nsm(g) = (nsm(g-1) + 1)/2
end do

!----------------------------------------------------------------------
end subroutine def_ngens

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module mg_parameter
