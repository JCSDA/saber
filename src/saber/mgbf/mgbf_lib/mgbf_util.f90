module mgbf_utils
use, intrinsic :: iso_fortran_env, only: error_unit
! MPI_COMM_WORLD appears in this module only, and only for fatal handling.
use mpi, only: MPI_Abort, MPI_COMM_WORLD, MPI_Comm_rank, MPI_Initialized, &
    MPI_SUCCESS
use mgbf_kinds,only: r_kind, i_kind
implicit none
private
public :: contains_nonzero, mgbf_abort, mgbf_world_rank

interface contains_nonzero
  module procedure contains_nonzero_real
  module procedure contains_nonzero_integer
end interface

! Cached world rank; negative until MPI has been initialized.
integer, save :: cached_rank = -1

contains

! Rank of this task in MPI_COMM_WORLD, or 0 before MPI_Init.
!
! This is deliberately NOT the component rank.  Its only purpose is to say which
! process in the job hit a fatal error, so that a failure raised deep in a
! numerical kernel -- which holds no MGBF object and therefore no component
! communicator -- can still be attributed to a specific process.  Anything that
! selects a root for normal logging, builds a filename, or drives decomposition
! uses the component rank (this%mype) instead.
integer function mgbf_world_rank() result(rank)
integer :: ierr
logical :: is_initialized

if (cached_rank < 0) then
  call MPI_Initialized(is_initialized, ierr)
  if (ierr /= MPI_SUCCESS .or. .not. is_initialized) then
    rank = 0
    return
  end if
  call MPI_Comm_rank(MPI_COMM_WORLD, cached_rank, ierr)
  if (ierr /= MPI_SUCCESS) then
    cached_rank = -1
    rank = 0
    return
  end if
end if
rank = cached_rank
end function mgbf_world_rank

! Fatal error exit.  The message is written unconditionally and tagged with the
! world rank, because the conditions that reach here are detected per subdomain
! and are rarely hit on rank 0.  This procedure never returns: every call site
! relies on that, so fall back to a hard stop when MPI cannot abort.
!
! Aborting MPI_COMM_WORLD rather than the failing component's communicator is
! deliberate.  MGBF may be running on one of several disjoint components; if we
! terminated only our own, the surviving components would hang at their next
! collective instead of failing.  A rank-local failure ends the whole job.
subroutine mgbf_abort(message)
character(len=*), intent(in) :: message
integer :: ierr
logical :: is_initialized
character(len=16) :: rank_str

write(rank_str, "(I0)") mgbf_world_rank()
write(error_unit, "(A)") "MGBF fatal error [rank " // trim(rank_str) // "]: " // trim(message)
flush(error_unit)

call MPI_Initialized(is_initialized, ierr)
if (ierr == MPI_SUCCESS .and. is_initialized) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
error stop 1  ! reached only when MPI is unavailable or the abort itself failed
end subroutine mgbf_abort

logical function contains_nonzero_real(array) result(has_nonzero)
real(r_kind), intent(in) :: array(:, :, :)  ! Declare 3D array
integer(i_kind) :: i, j, k
real(r_kind), parameter :: zero_tol = epsilon(1.0_r_kind)

has_nonzero = .false.
do k = 1, size(array, 3)  ! Loop over third dimension
  do j = 1, size(array, 2)  ! Loop over second dimension
    do i = 1, size(array, 1)  ! Loop over first dimension
      if (abs(array(i, j, k)) > zero_tol) then
        has_nonzero = .true.
        return
      end if
    end do
  end do
end do
end function contains_nonzero_real

logical function contains_nonzero_integer(array) result(has_nonzero)
integer(i_kind), intent(in) :: array(:, :, :)  ! Declare 3D array
integer(i_kind) :: i, j, k

has_nonzero = .false.
do k = 1, size(array, 3)  ! Loop over third dimension
  do j = 1, size(array, 2)  ! Loop over second dimension
    do i = 1, size(array, 1)  ! Loop over first dimension
      if (array(i, j, k) /= 0) then
        has_nonzero = .true.
        return
      end if
    end do
  end do
end do
end function contains_nonzero_integer
end module mgbf_utils
