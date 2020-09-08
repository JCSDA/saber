! (C) Copyright 2019-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> \brief Bump Interpolation module
!!
!!
!! \date Jan, 2020: Created by M. Miesch (JCSDA/UCAR) based on previously exising type_obsop.F90
!! written by Benjamin Menetrier (JCSDA/UCAR, CERFACS, METEO-FRANCE, IRIT)
!!
module bump_interpolation_mod

use atlas_module
use atlas_metadata_module, only: atlas_Metadata
use fckit_configuration_module, only : fckit_configuration
use fckit_log_module, only: fckit_log
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum,fckit_mpi_min,&
                            fckit_mpi_max,fckit_mpi_status
use netcdf
use tools_atlas
use tools_const, only: pi,deg2rad,rad2deg
use tools_func, only: lonlatmod, sphere_dist
use tools_kinds, only: kind_real
use tools_qsort, only: qsort
use tools_repro, only: rth
use type_bump, only : bump_type
use type_com, only: com_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none
private

public :: bump_interpolator
public :: bump_interpolator_registry

integer, parameter :: max_string = 1024

! ------------------------------------------------------------------------------
!!
!! bint = shorthand for bump interpolator
!!
type bump_interpolator
private
  type(bump_type), public :: bump

  !> Grid information
  type(geom_type) :: outgeom    !< output grid - unstructured
  type(atlas_functionspace) :: in_funcspace  !< atlas functionspace for input grid
  type(atlas_functionspace) :: out_funcspace !< atlas functionspace for output grid

  !> Number of points
  integer :: nc0b                !< Halo B size
  integer, public :: nlev        !< number of levels

  integer, public :: nout         !< global number of output grid points
  integer, public :: nout_local   !< local number of output grid points

  !> Interpolation data (operator)
  type(linop_type) :: h

  !> Communication data
  type(com_type) :: com

contains
  private

  procedure, private :: driver => bint_driver

  procedure, public :: init => bint_init

  procedure, public :: apply => bint_apply
  procedure, public :: apply_ad => bint_apply_ad

  procedure :: deallocate_outgrid => bint_deallocate_outgrid
  procedure, public :: delete => bint_delete

  ! low-level apply methods
  procedure, private :: apply_interp, apply_interp_ad
  final :: dummy

end type bump_interpolator

! ------------------------------------------------------------------------------
!> Registry for bump_interpolator objects

#define LISTED_TYPE bump_interpolator

!> Linked list interface - defines registry_t type
#include "saber/util/linkedList_i.f"

!> Global registry
type(registry_t) :: bump_interpolator_registry

! ------------------------------------------------------------------------------

contains
!-------------------------------------------------------------------------------
!> Linked list implementation
#include "saber/util/linkedList_c.f"

! ------------------------------------------------------------------------------
!> Initialize interpolation object
!!
!! The input and output fields are atlas_FieldSet objects that are assumed
!! to be created from atlas functionspaces.  So, they have the grid and
!! mesh information built in.
!!
!! \param[in] config = configuration
!!
!! \param[in] in_afs = This is the input grid, rendered as an atlas
!! functionspace, so it includes information about the mesh (parallel
!! connectivity) as well.
!!
!! \param[in] out_afs = This is the output grid, also represented as
!! an altas_functionspace.
!!
!! \param[in] masks = This contains metadata needed for the interpolation.
!!  It is rendered as an atlas FieldSet with the following named fields.
!!  Each of these named fields is optional; if omitted default values will
!!  be provided
!! * area  : cell area
!! * vunit : vertical unit
!! * gmask : geometry mask
!! * smask : sampling mask
!!
subroutine bint_init(self, config, comm, in_funcspace, out_funcspace, masks)
  use fckit_configuration_module, only : fckit_configuration
  use iso_c_binding, only : c_char
  class(bump_interpolator), intent(inout) :: self
  class(fckit_configuration), intent(in)  :: config
  type(fckit_mpi_comm)      , intent(in)  :: comm
  class(atlas_functionspace), intent(in)  :: in_funcspace
  class(atlas_functionspace), intent(in)  :: out_funcspace
  type(atlas_fieldset)      , intent(in), optional :: masks

  ! local variables
  type(fckit_configuration) :: bump_config
  integer :: msvali, j
  real(kind_real) :: msvalr
  integer, allocatable :: levels(:)
  character(len=max_string) :: msg
  character(len=max_string) :: myname = "saber::interpolation::bump_interpolation_mod::bint_init "

  !--------------------------------------------------------------------------------
  ! set bump namelist parameters.

  call config%get_or_die("bump",bump_config)

  ! bump defaults
  call self%bump%nam%init(comm%size())  ! set up defaults

  ! modified defaults for bump interpolator
  self%bump%nam%default_seed = .true.
  self%bump%nam%prefix       = 'bump_interpolator'
  self%bump%nam%datadir      = 'testdata'

  ! optionally override defaults with config
  call self%bump%nam%from_conf(bump_config)

  ! optional missing values for integers and reals
  msvali = -999
  msvalr = -999.0
  If (config%has("missingvalue_int")) &
       call config%get_or_die("missingvalue_int",msvali)
  If (config%has("missingvalue_real")) &
       call config%get_or_die("missingvalue_real",msvalr)

  ! save these for future use
  self%in_funcspace = in_funcspace
  self%out_funcspace = out_funcspace

  !--------------------------------------------------------------------------------
  ! determine the number of vertical levels

  self%nlev = 1
  If (config%has("nlevels")) call config%get_or_die("nlevels",self%nlev)

  self%bump%nam%nl = self%nlev

  allocate(levels(1:self%nlev))
  if (config%has("levels")) then
     call config%get_or_die("levels",levels)
  else
     do j = 1, self%nlev
        levels(j) = j
     enddo
  endif
  self%bump%nam%levs(1:self%nlev) = levels
  deallocate(levels)

  !--------------------------------------------------------------------------------
  ! Initialize BUMP
  ! ---------------

  ! clear any pre-existing data
  call self%delete

  ! basic bump setup
  if (present(masks)) then
     call self%bump%setup(self%bump%mpl%f_comm, in_funcspace, masks, &
          msvali=msvali, msvalr=msvalr)
  else
     call self%bump%setup(self%bump%mpl%f_comm, in_funcspace, &
          msvali=msvali, msvalr=msvalr)
  endif

  !--------------------------------------------------------------------------------
  ! initialize output grid

  call fckit_log%info('-------------------------------------------------------------------')
  call fckit_log%info('--- Initialize output grid')

  self%outgeom%nl0 = self%nlev

  ! unpack output grid from output functionspace
  call self%outgeom%from_atlas(self%bump%mpl, out_funcspace)
  self%nout_local = size(self%outgeom%lon_mga)

  !--------------------------------------------------------------------------------
  ! Run basic BUMP drivers
  ! The only thing here that may be needed is to possibly run the nicas driver
  ! after we get things to work, try commenting this out and see if it still works
  call self%bump%run_drivers

  ! Run interpolation driver
  call fckit_log%info('-------------------------------------------------------------------')
  call fckit_log%info('--- Run bump_interpolation driver')
  call self%driver(self%bump%mpl,self%bump%rng,self%bump%nam,self%bump%geom)
  if (self%bump%nam%default_seed) call self%bump%rng%reseed(self%bump%mpl)

  !--------------------------------------------------------------------------------
  ! more initializations and checks that the setup is correct

  if (self%nout < 1) then
     write(msg,'(a)') trim(myname) // ": " // "Output grid not properly defined"
     call abor1_ftn(msg)
  endif

  if (self%nlev /= self%bump%geom%nl0) then
  write(msg,'(a)') trim(myname) // ": " // &
                             "number of levels does not match bump geom"
  call abor1_ftn(msg)
  endif

  if (self%nlev /= self%outgeom%nl0) then
     write(msg,'(a)') trim(myname) // ": " // &
          "ERROR: Output grid has different number of levels than input grid!"
     call abor1_ftn(msg)
  endif

  !--------------------------------------------------------------------------------

end subroutine bint_init

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
!> Initialize bump to perform interpolation.
!!
!----------------------------------------------------------------------
subroutine bint_driver(self,mpl,rng,nam,geom)
  class(bump_interpolator), intent(inout) :: self
  type(mpl_type),intent(inout) :: mpl      !< MPI data
  type(rng_type),intent(inout) :: rng      !< Random number generator
  type(nam_type),intent(in) :: nam         !< Namelist
  type(geom_type),intent(in) :: geom       !< Geometry

  ! Local variables
  integer :: iouta,iproc,i_s,ic0,ic0u,jc0u,ic0b,ic0a,nouta_eff
  integer :: nout_eff,nn_index(1),proc_to_nouta(mpl%nproc),proc_to_nouta_eff(mpl%nproc)
  integer :: c0u_to_c0b(geom%nc0u)
  integer,allocatable :: c0b_to_c0(:)
  real(kind_real) :: nn_dist(1),N_max,C_max
  logical :: maskouta(self%nout_local),lcheck_nc0b(geom%nc0)
  character(len=1024),parameter :: subr = 'bint_driver'

  ! Check that universe is global
  if (any(.not.geom%myuniverse)) call mpl%abort(subr,'universe should be global for interpolation')

  ! Check whether output grid points are inside the mesh
  if (self%nout_local > 0) then
     do iouta=1,self%nout_local
        call geom%mesh_c0u%inside(mpl,self%outgeom%lon_mga(iouta),self%outgeom%lat_mga(iouta),maskouta(iouta))
        if (.not.maskouta(iouta)) then
           ! Check for very close points
           call geom%tree_c0u%find_nearest_neighbors(self%outgeom%lon_mga(iouta),self%outgeom%lat_mga(iouta), &
 & 1,nn_index,nn_dist)
           if (nn_dist(1)<rth) maskouta(iouta) = .true.
        end if
     end do
     nouta_eff = count(maskouta)
  else
     nouta_eff = 0
  endif

  ! Get global number of output grid points
  call mpl%f_comm%allgather(self%nout_local,proc_to_nouta)
  call mpl%f_comm%allgather(nouta_eff,proc_to_nouta_eff)
  self%nout = sum(proc_to_nouta)
  nout_eff = sum(proc_to_nouta_eff)

  ! Print input
  write(mpl%info,'(a7,a)') '','Number of points in output grid / valid points per MPI task:'
  call mpl%flush
  do iproc=1,mpl%nproc
     write(mpl%info,'(a10,a,i3,a,i8,a,i8)') '','Task ',iproc,': ', &
           proc_to_nouta(iproc),' / ',proc_to_nouta_eff(iproc)
     call mpl%flush
  end do
  write(mpl%info,'(a10,a,i8,a,i8)') '','Total   : ',self%nout,' / ',nout_eff
  call mpl%flush

  ! Compute interpolation
  self%h%prefix = 'o'
  write(mpl%info,'(a7,a)') '','Single level:'
  call mpl%flush
  call self%h%interp(mpl,rng,nam,geom,0,geom%nc0u,geom%lon_c0u,geom%lat_c0u,&
                     geom%gmask_hor_c0u,self%nout_local,self%outgeom%lon_mga,&
                     self%outgeom%lat_mga,maskouta,10)

  ! Define halo B
  lcheck_nc0b = .false.
  do ic0a=1,geom%nc0a
     ic0u = geom%c0a_to_c0u(ic0a)
     lcheck_nc0b(ic0u) = .true.
  end do
  do iouta=1,self%nout_local
     do i_s=1,self%h%n_s
        jc0u = self%h%col(i_s)
        lcheck_nc0b(jc0u) = .true.
     end do
  end do
  self%nc0b = count(lcheck_nc0b)

  ! Allocation
  allocate(c0b_to_c0(self%nc0b))

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
  self%h%n_src = self%nc0b
  do i_s=1,self%h%n_s
     self%h%col(i_s) = c0u_to_c0b(self%h%col(i_s))
  end do

  ! Setup communications
  call self%com%setup(mpl,'com',geom%nc0a,self%nc0b,geom%nc0,geom%c0a_to_c0,c0b_to_c0)

  ! Compute scores
  if (nout_eff > 0) then
     call mpl%f_comm%allreduce(real(self%com%nhalo,kind_real),C_max,fckit_mpi_max())
     C_max = C_max/(3.0*real(nout_eff,kind_real)/real(mpl%nproc,kind_real))
     N_max = real(maxval(proc_to_nouta_eff),kind_real)/(real(nout_eff,kind_real)/real(mpl%nproc,kind_real))

     ! Print results
     write(mpl%info,'(a7,a,f5.1,a)') '','Output grid repartition imbalance: ',100.0*real(maxval(proc_to_nouta_eff) &
 & -minval(proc_to_nouta_eff),kind_real)/(real(sum(proc_to_nouta_eff),kind_real)/real(mpl%nproc,kind_real)),' %'
     call mpl%flush
     write(mpl%info,'(a7,a,i8,a,i8,a,i8)') '','Number of grid points / halo size / number of received values: ', &
 & self%com%nred,' / ',self%com%next,' / ',self%com%nhalo
     call mpl%flush
     write(mpl%info,'(a7,a,f10.2,a,f10.2)') '','Scores (N_max / C_max):',N_max,' / ',C_max
     call mpl%flush
  endif

  ! Release memory
  deallocate(c0b_to_c0)

end subroutine bint_driver

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
!> Apply interpolation
!!
!! \param[in] instate = input fields represented as an atlas fieldset,
!! created from a functionspace
!!
!! \param[out] outstate = output fields represented as an atlas fieldset.
!! If the fields that constitudte the fieldset are not already allocated by
!! the caller, then they will be created and allocated by this method.
!! So, the user can optionally pass this routine an empty output fieldset.
!!
subroutine bint_apply(self, infields, outfields)
  class(bump_interpolator), intent(inout) :: self
  type(atlas_fieldset),intent(in)         :: infields
  type(atlas_fieldset),intent(inout)      :: outfields

  ! local variables
  type(atlas_field) :: infield, outfield
  real(kind_real), allocatable :: infld_mga(:,:), infld_c0a(:,:)
  real(kind_real), allocatable :: outfld(:,:)
  integer :: ifield
  character(len=max_string) :: fieldname

  ! allocate bump arrays
  allocate(infld_mga(self%bump%geom%nmga,self%nlev))
  allocate(outfld(self%nout_local,self%nlev))

  if (.not.self%bump%geom%same_grid) then
     allocate(infld_c0a(self%bump%geom%nc0a, self%nlev))
  endif

  do ifield = 1, infields%size()

     infield = infields%field(ifield)
     fieldname = infield%name()

     !--------------------------------------------
     ! allocate output field if necessary

     if (.not. outfields%has_field(fieldname)) then
        outfield = self%out_funcspace%create_field(name=fieldname, &
             kind=atlas_real(kind_real),levels=self%nlev)
     else
        outfield = outfields%field(name=fieldname)
     endif

     !--------------------------------------------
     ! compute interpolation

     ! atlas field to bump fld
     call field_to_fld(self%bump%mpl, infield, infld_mga)

     if (self%bump%geom%same_grid) then
        call self%apply_interp(infld_mga,outfld)
     else
        ! Model grid to subset Sc0
        infld_c0a = 0.0_kind_real
        call self%bump%geom%copy_mga_to_c0a(self%bump%mpl,infld_mga,infld_c0a)

        ! Apply interpolation
        call self%apply_interp(infld_c0a,outfld)
     end if

     ! bump fld to atlas field
     call fld_to_field(self%bump%mpl, outfld, outfield)

     !--------------------------------------------
     ! Add output field to output fields
     if (.not. outfields%has_field(fieldname)) then
        call outfields%add(outfield)
     endif

     ! release pointers
     call infield%final()
     call outfield%final()

  enddo

  ! clean up
  if (allocated(infld_mga)) deallocate(infld_mga)
  if (allocated(infld_c0a)) deallocate(infld_c0a)
  if (allocated(outfld)) deallocate(outfld)

end subroutine bint_apply

!----------------------------------------------------------------------
!> Subroutine: apply_interp
!! Purpose: low-level routine to apply the interpolation to a single
!! field on a single level
!! \param[in]  infield: input field
!! \param[out] outfield: output field
!!
subroutine apply_interp(self,infield,outfield)
  class(bump_interpolator), intent(inout) :: self
  real(kind_real),intent(in)   :: infield(self%bump%geom%nc0a,self%nlev)
  real(kind_real),intent(out)  :: outfield(self%nout_local,self%nlev)

  integer :: ilev
  real(kind_real), allocatable :: infield_ext(:,:)

  allocate(infield_ext(self%nc0b,self%nlev))

  ! Halo extension
  call self%com%ext(self%bump%mpl, self%nlev, infield ,infield_ext)

  if (self%nout_local > 0) then
     ! Horizontal interpolation
     !$omp parallel do schedule(static) private(ilev)
     do ilev = 1, self%nlev
        call self%h%apply(self%bump%mpl,infield_ext(:,ilev),outfield(:,ilev))
     end do
     !$omp end parallel do
  end if

  deallocate(infield_ext)

end subroutine apply_interp

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!> Apply interpolator operator adjoint
!!
!! \param[in] fields_outgrid = These are the fields on the second grid,
!!  i.e. the target grid of the original interpolation.  For the adjoint,
!!  these fields are treated as an input.
!!
!! \param[inout] fields_ingrid = These are the fields defined on the
!! first grid, i.e. the source grid of the original interpolation.  For
!! the adjoint, these are treated as an output.  The caller can optionally
!! pass this arguement as an empty fieldset and the routine will create
!! and allocate each component of the fieldset.  Or, if the field components
!! of the fieldset are already allocated by the caller, then this routine
!! will merely replace the field values with the result of the computation.
!!
!----------------------------------------------------------------------
subroutine bint_apply_ad(self, fields_outgrid, fields_ingrid)
  class(bump_interpolator), intent(inout) :: self
  type(atlas_fieldset), intent(in)    :: fields_outgrid
  type(atlas_fieldset), intent(inout) :: fields_ingrid

  ! Local variables
  type(atlas_field) :: field_ingrid, field_outgrid
  real(kind_real), allocatable :: fld_ingrid_mga(:,:), fld_ingrid_c0a(:,:)
  real(kind_real), allocatable :: fld_outgrid(:,:)
  integer :: ifield
  character(len=max_string) :: fieldname

  ! allocate bump arrays
  allocate(fld_ingrid_mga(self%bump%geom%nmga,self%nlev))
  allocate(fld_outgrid(self%nout_local,self%nlev))

  if (.not.self%bump%geom%same_grid) then
     allocate(fld_ingrid_c0a(self%bump%geom%nc0a, self%nlev))
  endif

  do ifield = 1, fields_outgrid%size()

     field_outgrid = fields_outgrid%field(ifield)
     fieldname = field_outgrid%name()

     !--------------------------------------------
     ! allocate output field if necessary

     if (.not. fields_ingrid%has_field(fieldname)) then
        field_ingrid = self%in_funcspace%create_field(name=fieldname, &
             kind=atlas_real(kind_real),levels=self%nlev)
     else
        field_ingrid = fields_ingrid%field(name=fieldname)
     endif

     !--------------------------------------------
     ! compute interpolation adjoint

     ! atlas field to bump fld
     call field_to_fld(self%bump%mpl, field_outgrid, fld_outgrid)

     if (self%bump%geom%same_grid) then
        ! Apply observation operator adjoint
        call self%apply_interp_ad(fld_outgrid,fld_ingrid_mga)
     else
        ! Apply observation operator adjoint
        call self%apply_interp_ad(fld_outgrid,fld_ingrid_c0a)

        ! Subset Sc0 to model grid
        call self%bump%geom%copy_c0a_to_mga(self%bump%mpl,fld_ingrid_c0a,fld_ingrid_mga)
     end if

     ! bump fld to atlas field
     call fld_to_field(self%bump%mpl, fld_ingrid_mga, field_ingrid)

     !--------------------------------------------

     ! add field to result
     if (.not. fields_ingrid%has_field(fieldname)) then
        call fields_ingrid%add(field_ingrid)
     endif

     ! release pointers
     call field_ingrid%final()
     call field_outgrid%final()

  enddo

  ! clean up
  if (allocated(fld_ingrid_mga)) deallocate(fld_ingrid_mga)
  if (allocated(fld_ingrid_c0a)) deallocate(fld_ingrid_c0a)
  if (allocated(fld_outgrid)) deallocate(fld_outgrid)

end subroutine bint_apply_ad

!----------------------------------------------------------------------
!> Subroutine: apply_interp_ad
!! Purpose: low-level routine to apply the adjoint of the interpolation operator
!! to a single field on a single level
!!
!! \param[in]  mpl bump MPI data
!! \param[in]  geom bump geometry data
!! \param[in]  fld_outgrid field on output grid
!! \param[out] fld_ingrid field on input grid
!!

subroutine apply_interp_ad(self,fld_outgrid,fld_ingrid)
  class(bump_interpolator), intent(inout) :: self
  real(kind_real),intent(in)   :: fld_outgrid(self%nout_local,self%nlev)
  real(kind_real),intent(out)  :: fld_ingrid(self%bump%geom%nc0a,self%nlev)

  integer :: ilev
  real(kind_real) :: fld_ingrid_ext(self%nc0b,self%nlev)

  if (self%nout_local > 0) then
     ! Horizontal interpolation
     !$omp parallel do schedule(static) private(ilev)
     do ilev = 1, self%nlev
        call self%h%apply_ad(self%bump%mpl,fld_outgrid(:,ilev),fld_ingrid_ext(:,ilev))
     end do
     !$omp end parallel do
  else
     ! No observation on this task
     fld_ingrid_ext = 0.0
  end if

  ! Halo reduction
  call self%com%red(self%bump%mpl,self%nlev,fld_ingrid_ext,fld_ingrid)

end subroutine apply_interp_ad

!----------------------------------------------------------------------
!> Release memory (partial) by deallocating output grid
!----------------------------------------------------------------------
subroutine bint_deallocate_outgrid(self)
  class(bump_interpolator), intent(inout) :: self

  ! Release memory
  call self%outgeom%dealloc()

end subroutine bint_deallocate_outgrid

! ------------------------------------------------------------------------------
!> Release all memory
! ------------------------------------------------------------------------------

subroutine bint_delete(self)
  class(bump_interpolator), intent(inout) :: self

  call self%bump%dealloc()

  call self%deallocate_outgrid()
  call self%h%dealloc()
  call self%com%dealloc()

end subroutine bint_delete

!----------------------------------------------------------------------
! Subroutine: dummy
! Purpose: dummy finalization
!----------------------------------------------------------------------
subroutine dummy(bump)

implicit none

! Passed variables
type(bump_interpolator),intent(inout) :: bump ! BUMP

end subroutine dummy

end module bump_interpolation_mod
