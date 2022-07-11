! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!>  Fortran interface to Logger

module logger_mod

use, intrinsic :: iso_c_binding
use fckit_c_interop_module, only : c_str_right_trim

implicit none

public :: oops_log

type :: oops_log_type
contains
  procedure, nopass, public :: info    !< Log to info channel
  procedure, nopass, public :: error   !< Log to error channel
  procedure, nopass, public :: warning !< Log to warning channel
  procedure, nopass, public :: debug   !< Log to debug channel
  procedure, nopass, public :: trace   !< Log to trace channel
  procedure, nopass, public :: stats   !< Log to stats channel
  procedure, nopass, public :: test    !< Log to test channel
  procedure, nopass, public :: timer   !< Log to timer channel
end type oops_log_type

type(oops_log_type) :: oops_log

interface
   subroutine c_log_info(msg,newl,flush) bind(c, name='log_info_f')
      use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_info
!-------------------------------------------------------------------------------
subroutine c_log_error(msg,newl,flush) bind(c, name='log_error_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_error
!-------------------------------------------------------------------------------
subroutine c_log_warning(msg,newl,flush) bind(c, name='log_warning_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_warning
!-------------------------------------------------------------------------------
subroutine c_log_debug(msg,newl,flush) bind(c, name='log_debug_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_debug
!-------------------------------------------------------------------------------
subroutine c_log_trace(msg,newl,flush) bind(c, name='log_trace_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_trace
!-------------------------------------------------------------------------------
subroutine c_log_stats(msg,newl,flush) bind(c, name='log_stats_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_stats
!-------------------------------------------------------------------------------
subroutine c_log_test(msg,newl,flush) bind(c, name='log_test_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_test
!-------------------------------------------------------------------------------
subroutine c_log_timer(msg,newl,flush) bind(c, name='log_timer_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_timer
end interface

private

contains

subroutine info(msg,newl,flush)

implicit none

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_info(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine info

!-------------------------------------------------------------------------------

subroutine error(msg,newl,flush)

implicit none

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_error(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine error

!-------------------------------------------------------------------------------

subroutine warning(msg,newl,flush)

implicit none

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_warning(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine warning

!-------------------------------------------------------------------------------

subroutine debug(msg,newl,flush)

implicit none

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_debug(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine debug

!-------------------------------------------------------------------------------

subroutine trace(msg,newl,flush)

implicit none

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_trace(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine trace

!-------------------------------------------------------------------------------

subroutine stats(msg,newl,flush)

implicit none

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_stats(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine stats

!-------------------------------------------------------------------------------

subroutine test(msg,newl,flush)

implicit none

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_test(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine test

!-------------------------------------------------------------------------------

subroutine timer(msg,newl,flush)

implicit none

! Passed variables
character(kind=c_char,len=*), intent(in) :: msg !< Message to be logged
logical, intent(in), optional :: newl           !< Add newline character (```\n```) after message. Default ```.true.```
logical, intent(in), optional :: flush          !< Flush channel after message. Default ```.true.```

! Local variables
integer(c_int32_t) :: opt_newl, opt_flush

! Local flags
opt_newl  = 1
if (present(newl)) then
  if (.not.newl) opt_newl  = 0
end if
opt_flush = 1
if (present(flush)) then
  if(.not.flush) opt_flush = 0
end if

! Call C++
call c_log_timer(c_str_right_trim(msg),opt_newl,opt_flush)

end subroutine timer

end module Logger_mod
