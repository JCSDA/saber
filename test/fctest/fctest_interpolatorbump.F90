! (C) Copyright 2022 UCAR
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! @author Benjamin Menetrier

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

TESTSUITE(fctest_interpolatorbump)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  use atlas_module
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  use atlas_module
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_interpolatorbump )
  use atlas_module, only: atlas_structuredgrid,atlas_functionspace,atlas_functionspace_structuredcolumns
  use fckit_mpi_module, only: fckit_mpi_comm,fckit_mpi_sum
  use interpolatorbump_mod, only: bump_interpolator
  use tools_kinds, only: kind_real
  use tools_repro, only: repro_th
  use type_fieldset, only: fieldset_type

  implicit none

  ! Local variables
  integer,parameter :: nl0 = 5
  integer :: nmga_in,nmga_out
  real(kind_real) :: dp_in,dp_out
  real(kind_real),allocatable :: array_in_1(:,:,:),array_in_2(:,:,:),array_out_1(:,:,:),array_out_2(:,:,:)
  logical,allocatable :: gmask_in(:,:),gmask_out(:,:)
  logical,dimension(2),parameter :: var2d = (/.false.,.false./)
  character(len=4),dimension(2),parameter :: variables = (/'var1','var2'/)
  character(len=5),parameter :: lev2d = 'first'
  type(fckit_mpi_comm) :: f_comm
  type(atlas_structuredgrid) :: grid_in,grid_out
  type(atlas_functionspace) :: fspace_in,fspace_out
  type(atlas_functionspace_structuredcolumns) :: fspace_in_sc,fspace_out_sc
  type(bump_interpolator) :: bumpinterp
  type(fieldset_type) :: fset_in_1,fset_in_2,fset_out_1,fset_out_2

  ! Initialize communicator
  f_comm = fckit_mpi_comm("world")

  ! Create input and output grids
  grid_in = atlas_structuredgrid("F15")
  grid_out = atlas_structuredgrid("F10")

  ! Create input and output function spaces
  fspace_in = atlas_functionspace_structuredcolumns(grid_in)
  fspace_out = atlas_functionspace_structuredcolumns(grid_out)

  ! Initialize interpolator
  call bumpinterp%init(f_comm,fspace_in,fspace_out,nl0)

  ! Create input fieldset
  fspace_in_sc = atlas_functionspace_structuredcolumns(fspace_in%c_ptr())
  nmga_in = fspace_in_sc%size_owned()
  allocate(gmask_in(nmga_in,nl0))
  gmask_in = .true.
  call fset_in_1%init(bumpinterp%bump%mpl,fspace_in,gmask_in,variables,lev2d,var2d)
  call fset_in_2%init(bumpinterp%bump%mpl,fspace_in,gmask_in,variables,lev2d,var2d)

  ! Initialize input fieldsets
  allocate(array_in_1(nmga_in,nl0,2))
  allocate(array_in_2(nmga_in,nl0,2))
  call bumpinterp%bump%rng%rand_gau(array_in_1)
  array_in_2 = array_in_1
  call fset_in_1%from_array(bumpinterp%bump%mpl,array_in_1)
  call fset_in_2%from_array(bumpinterp%bump%mpl,array_in_2)

  ! Create output fieldset
  fspace_out_sc = atlas_functionspace_structuredcolumns(fspace_out%c_ptr())
  nmga_out = fspace_out_sc%size_owned()
  allocate(gmask_out(nmga_out,nl0))
  gmask_out = .true.
  call fset_out_1%init(bumpinterp%bump%mpl,fspace_out,gmask_out,variables,lev2d,var2d)
  call fset_out_2%init(bumpinterp%bump%mpl,fspace_out,gmask_out,variables,lev2d,var2d)

  ! Initialize output fieldset
  allocate(array_out_1(nmga_out,nl0,2))
  allocate(array_out_2(nmga_out,nl0,2))
  call bumpinterp%bump%rng%rand_gau(array_out_1)
  array_out_2 = array_out_1
  call fset_out_1%from_array(bumpinterp%bump%mpl,array_out_1)
  call fset_out_2%from_array(bumpinterp%bump%mpl,array_out_2)

  ! Apply interpolator
  call bumpinterp%apply(fset_in_1,fset_out_1)

  ! Apply adjoint interpolator
  call bumpinterp%apply_ad(fset_out_2,fset_in_2)

  ! Adjoint test
  call fset_in_2%to_array(bumpinterp%bump%mpl,array_in_2)
  call fset_out_1%to_array(bumpinterp%bump%mpl,array_out_1)
  dp_in = sum(array_in_1*array_in_2)
  dp_out = sum(array_out_1*array_out_2)
  call bumpinterp%bump%mpl%f_comm%allreduce(dp_in,fckit_mpi_sum())
  call bumpinterp%bump%mpl%f_comm%allreduce(dp_out,fckit_mpi_sum())
  FCTEST_CHECK_CLOSE(dp_in,dp_out,repro_th)

  ! Delete interpolator
  call bumpinterp%delete()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
