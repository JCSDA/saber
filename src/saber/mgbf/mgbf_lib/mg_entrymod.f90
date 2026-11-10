submodule(mg_intstate) mg_entrymod
!$$$  submodule documentation block
!                .      .    .                                       .
! module:   mg_entrymod
!   prgmmr: rancic           org: NCEP/EMC            date: 2020
!
! abstract:  Initialize and finialize multigrid Beta filter
!            for modeling of background error covariance
!
! module history log:
!   2023-04-19  lei     - object-oriented coding
!   2024-01-11  rancic  - optimization for ensemble localization
!   2024-02-20  yokota  - refactoring to apply for GSI
!
! Subroutines Included:
!   mg_initialize -
!   mg_finalize -
!
! Functions Included:
!
! remarks:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use mgbf_kinds, only: r_kind,i_kind
use mgbf_utils, only: mgbf_abort

implicit none

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine mg_initialize(this,component_comm,n_owned_anl,anl_lonlat1d,inputfilename,obj_parameter)
use mpi, only: MPI_Comm_rank, MPI_Comm_size, MPI_SUCCESS
implicit none
!**********************************************************************!
!                                                                      !
!   Initialization subroutine                                          !
!                                                     M. Rancic (2020) !
!***********************************************************************
class (mg_intstate_type), intent(inout) :: this
integer(i_kind),intent(in):: component_comm
integer(i_kind),optional,intent(in):: n_owned_anl
real(r_kind),optional,intent(in):: anl_lonlat1d(:,:)
character(len=*),optional,intent(in) :: inputfilename

class(mg_parameter_type),optional,intent(in):: obj_parameter
character(len=256) :: abort_msg
integer(i_kind) :: ierr

!---------------------------------------------------------------------------
!
!               Firs set of subroutines is called only once and serves to
!               initialte the MGBF run
!
!---------------------------------------------------------------------------

!****
!**** Adopt the component communicator supplied by the caller.  It is borrowed:
!**** this object never duplicates it and never frees it.  Rank and size have to
!**** be known before the namelist is read, because the normalization filename
!**** and the vertical-profile broadcast are both component-scoped.
!****
this%mpi_comm_comp = component_comm
call MPI_Comm_rank(this%mpi_comm_comp, this%mype, ierr)
if (ierr /= MPI_SUCCESS) then
  call mgbf_abort("MPI_Comm_rank failed on the supplied component communicator")
end if
call MPI_Comm_size(this%mpi_comm_comp, this%npes, ierr)
if (ierr /= MPI_SUCCESS) then
  call mgbf_abort("MPI_Comm_size failed on the supplied component communicator")
end if

!****
!**** Initialize run multigrid Beta filter parameters
!****
if (present(inputfilename)) then
  call this%init_mg_parameter(inputfilename)
else if (present(obj_parameter)) then
  ! Intrinsic derived-type assignment would copy the MPI handles held by
  ! mg_parameter_type, giving two objects the same communicator and group and
  ! turning cleanup into a double free.  Configuration-only copying is not
  ! implemented, so the interface is kept but the path is refused.
  call mgbf_abort("initialization from obj_parameter is not implemented")
end if

if (present(anl_lonlat1d)) then
  if (.not. present(n_owned_anl)) then
    call mgbf_abort("anl_lonlat1d supplied without n_owned_anl")
  end if
  if (size(anl_lonlat1d,2) /= 2 .or. size(anl_lonlat1d,1) <  n_owned_anl) then
    write(abort_msg,'("anl_lonlat1d has wrong shape: (",I0,",",I0,"), n_owned_anl=",I0)') &
        size(anl_lonlat1d,1), size(anl_lonlat1d,2), n_owned_anl
    call mgbf_abort(abort_msg)
  end if

end if

!****
!**** Initialize the MPI topology.  This runs for every component size: a
!**** one-rank component needs its decomposition metadata established on the
!**** same path as a multi-rank one.  Only the reduced filter communicator is
!**** created conditionally, inside init_mg_MPI.
!****
call this%init_mg_MPI

!***
!*** Initialize integration domain
!***
call this%init_mg_domain
if(this%l_loc) then
  call this%init_domain_loc
end if

!---------------------------------------------------------------------------
!
!               All others are function of km2,km3,km,nm,mm,im,jm
!               and needs to be called separately for each application
!
!---------------------------------------------------------------------------
!***
!*** Define km and WORKA array based on input from mg_parameters and
!*** depending on specific application
!***

!***
!*** Allocate variables, define weights, prepare mapping
!*** between analysis and filter grid
!***

call this%allocate_mg_intstate

call this%def_offset_coef
if(present(n_owned_anl).and.present(anl_lonlat1d)) then
  call this%def_mg_weights(n_owned_anl=n_owned_anl,lonlat1d_anl=anl_lonlat1d)
else
  call this%def_mg_weights
end if

if(this%mgbf_line) then
  call this%init_mg_line
end if

call this%lsqr_mg_coef

call this%lwq_vertical_coef(this%lm_a,this%lm,this%cvf1,this%cvf2,this%cvf3,this%cvf4,this%lref)

!***
!*** Just for testing of standalone version. In GSI WORKA will be given
!*** through a separate subroutine
!***



!-----------------------------------------------------------------------
end subroutine mg_initialize

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine mg_finalize(this)
!**********************************************************************!
!                                                                      !
!   Finalize multigrid Beta Function                                   !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi, only: MPI_Comm_free, MPI_Group_free, MPI_COMM_NULL, &
    MPI_GROUP_NULL, MPI_SUCCESS
implicit none
class (mg_intstate_type), intent(inout) :: this

real(r_kind), allocatable, dimension(:,:):: PA, VA
integer(i_kind):: n,m,L
integer:: nm,mm,lm
integer(i_kind):: ierr
!-----------------------------------------------------------------------

if(this%ldelta) then
  !
  ! Horizontal cross-section
  !
  nm=this%nm
  mm=this%mm
  lm=this%lm
end if

! Release the MPI resources this object created from the borrowed component
! communicator.  MPI_Comm_create returns MPI_COMM_NULL on ranks outside
! filter_group, so every handle is checked before it is freed, and reset
! afterwards so a second finalize cannot double free.
!
! The component communicator itself is NOT freed here: it is borrowed from the
! covariance object, which owns it and frees it once every intstate is done.
!
! There is deliberately no barrier before this cleanup.  It bought nothing --
! freeing is already correctly ordered by the collective semantics of
! MPI_Comm_free -- and it added a place for an unwinding error path to deadlock.
if (this%filter_comm /= MPI_COMM_NULL) then
  call MPI_Comm_free(this%filter_comm, ierr)
  if (ierr /= MPI_SUCCESS) call mgbf_abort("MPI_Comm_free of the filter communicator failed")
  this%filter_comm = MPI_COMM_NULL
end if

if (this%filter_group /= MPI_GROUP_NULL) then
  call MPI_Group_free(this%filter_group, ierr)
  if (ierr /= MPI_SUCCESS) call mgbf_abort("MPI_Group_free of the filter group failed")
  this%filter_group = MPI_GROUP_NULL
end if

if (this%component_group /= MPI_GROUP_NULL) then
  call MPI_Group_free(this%component_group, ierr)
  if (ierr /= MPI_SUCCESS) call mgbf_abort("MPI_Group_free of the component group failed")
  this%component_group = MPI_GROUP_NULL
end if

this%mpi_comm_comp = MPI_COMM_NULL

call this%deallocate_mg_intstate

!-----------------------------------------------------------------------
end subroutine mg_finalize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end submodule mg_entrymod
