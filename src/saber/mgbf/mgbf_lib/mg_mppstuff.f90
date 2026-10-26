submodule(mg_parameter) mg_mppstuff
!$$$  submodule documentation block
!                .      .    .                                       .
! module:   mg_mppstuff
!   prgmmr: rancic           org: NCEP/EMC            date: 2020
!
! abstract:  Everything related to mpi communication
!
! module history log:
!   2023-04-19  lei     - object-oriented coding
!   2024-01-11  rancic  - optimization for ensemble localization
!   2024-02-20  yokota  - refactoring to apply for GSI
!
! Subroutines Included:
!   init_mg_MPI -
!   barrierMPI -
!   finishMPI -
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

use mgbf_kinds, only: i_kind
use mgbf_utils, only: mgbf_abort
implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine init_mg_MPI(this)
!***********************************************************************
!                                                                      !
!     Initialize mpi                                                   !
!     Create group for filter grid                                     !
!                                                                      !
!***********************************************************************
use mpi, only: MPI_BARRIER, MPI_COMM_CREATE, MPI_COMM_GROUP, &
    MPI_COMM_NULL, MPI_COMM_RANK, MPI_COMM_SIZE, &
    MPI_DOUBLE, MPI_GROUP_INCL, MPI_GROUP_NULL, &
    MPI_SUCCESS

implicit none
class(mg_parameter_type), intent(inout), target :: this
integer(i_kind):: g,m
integer(i_kind), dimension(this%npes_filt):: out_ranks
integer(i_kind):: nf
integer(i_kind)::ierr
include "type_parameter_locpointer.inc"
include "type_parameter_point2this.inc"
!-----------------------------------------------------------------------

!***
!***  No communicator is created here.  mype and npes were established from the
!***  borrowed component communicator in mg_initialize, before the namelist was
!***  read; the component is defined entirely by what SABER supplied.
!***
dTYPE = MPI_DOUBLE

!***
!*** Analysis grid
!***

nx = mod(mype,nxm)+1
my = (mype/nxm)+1

!***
!***  Define PEs that handle high generations
!***

mype_hgen=-1
my_hgen=-1

if( mype < maxpe_filt-nxy(1)) then
  mype_hgen=mype+nxy(1)
end if
do g=1,gm
  if(maxpe_fgen(g-1)<= mype_hgen .and. mype_hgen< maxpe_fgen(g)) then
    my_hgen=g
  end if
end do
l_hgen = mype_hgen >-1

!***
!***  Chars
!***
write(c_mype,1000) mype
1000 format(i5.5)

!-----------------------------------------------------------------------
!
call MPI_BARRIER(mpi_comm_comp,ierr)
!
!-----------------------------------------------------------------------
!***
!***  Define group communicator for higher generations
!***
!
!  Associate a group with communicator this@mpi_comm_comp
!
component_group = MPI_GROUP_NULL
filter_group    = MPI_GROUP_NULL
filter_comm     = MPI_COMM_NULL

call MPI_COMM_GROUP(mpi_comm_comp,component_group,ierr)
if (ierr /= MPI_SUCCESS) then
  call mgbf_abort("MPI_Comm_group failed on the component communicator")
end if
!
!  Create a new group out of exising group
!
!$omp parallel do private(nf) schedule(static)
do nf = 1,npes_filt
  out_ranks(nf)=nf-1
end do
!$omp end parallel do

call MPI_GROUP_INCL(component_group,npes_filt,out_ranks,filter_group,ierr)
if (ierr /= MPI_SUCCESS) then
  call mgbf_abort("MPI_Group_incl failed building the filter group")
end if
!
!  Now create a new communicator associated with new group.  This is collective
!  over the component communicator and returns MPI_COMM_NULL on the ranks that
!  are not members of filter_group.
!
call MPI_COMM_CREATE(mpi_comm_comp, filter_group, filter_comm, ierr)
if (ierr /= MPI_SUCCESS) then
  call mgbf_abort("MPI_Comm_create failed building the filter communicator")
end if

if( mype < npes_filt) then

  call MPI_COMM_RANK(filter_comm,mype_gr,ierr)
  call MPI_COMM_SIZE(filter_comm,npes_gr,ierr)

else

  mype_gr= -1
  npes_gr= npes_filt

end if

!-----------------------------------------------------------------------
!
call MPI_BARRIER(mpi_comm_comp,ierr)
!
!-----------------------------------------------------------------------
end subroutine init_mg_MPI

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine barrierMPI(this)
!***********************************************************************
!                                                                      !
!                                                                      !
!***********************************************************************
use mpi, only: MPI_BARRIER

implicit none
class(mg_parameter_type), intent(inout), target :: this
integer(i_kind):: ierr
include "type_parameter_locpointer.inc"
include "type_parameter_point2this.inc"
!-----------------------------------------------------------------------

call MPI_BARRIER(mpi_comm_comp,ierr)

!-----------------------------------------------------------------------
end subroutine barrierMPI

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine finishMPI(this)
!***********************************************************************
!                                                                      !
!     Finalize MPI                                                     !
!                                                                      !
!***********************************************************************

implicit none
class(mg_parameter_type), intent(inout), target :: this
!
! don't need mpi_finalize if mgbf is a lib to be called from outside
if (this%mype == 0) write(6,*) "mggbf%finishmpi should not be invoked if mgbf is used as a lib"
!
!
!-----------------------------------------------------------------------
end subroutine finishMPI

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end submodule mg_mppstuff
