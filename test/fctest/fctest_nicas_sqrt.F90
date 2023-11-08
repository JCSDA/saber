! (C) Copyright 2022 UCAR
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! @author Benjamin Menetrier

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

TESTSUITE(fctest_nicas_sqrt)

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
  use atlas_module, only: atlas_real,atlas_fieldset,atlas_field, &
 & atlas_structuredgrid,atlas_functionspace,atlas_functionspace_structuredcolumns
  use fckit_configuration_module, only: fckit_configuration
  use fckit_mpi_module, only: fckit_mpi_comm,fckit_mpi_sum
  use tools_kinds, only: kind_real
  use tools_repro, only: repro_th
  use type_bump, only: bump_type
  use type_fieldset, only: fieldset_type

  implicit none

  ! Local variables
  integer,parameter :: nl0 = 5
  integer :: n,nmga_out
  integer,dimension(2),parameter :: var2d_int = (/0,0/)
  real(kind_real) :: dp_in,dp_out
  real(kind_real),pointer :: ptr_1(:),ptr_2(:)
  real(kind_real),allocatable :: array_out_1(:,:,:),array_out_2(:,:,:)
  logical,allocatable :: gmask_out(:,:)
  logical,dimension(2),parameter :: var2d = (/.false.,.false./)
  character(len=4),dimension(2),parameter :: variables = (/'var1','var2'/)
  character(len=5),parameter :: lev2d = 'first'
  type(fckit_mpi_comm) :: f_comm
  type(atlas_field) :: cv_1,cv_2
  type(atlas_structuredgrid) :: grid_out
  type(atlas_functionspace) :: fspace_out
  type(atlas_functionspace_structuredcolumns) :: fspace_out_sc
  type(bump_type) :: bump
  type(fieldset_type) :: fset,fset_out_1,fset_out_2
  type(fckit_configuration) :: conf
  type(fckit_configuration) :: rh(1),rv(1)

  ! Initialize communicator
  f_comm = fckit_mpi_comm('world')

  ! Create output grid
  grid_out = atlas_structuredgrid('F10')

  ! Create output function space
  fspace_out = atlas_functionspace_structuredcolumns(grid_out)

  ! Create empty fieldset
  fset = atlas_fieldset()

  ! Create configuration
  conf = fckit_configuration()
  call conf%set('drivers.multivariate strategy','crossed')
  call conf%set('drivers.compute nicas',.true.)
  call conf%set('model.variables',variables)
  call conf%set('model.nl0',nl0)
  call conf%set('model.lev2d',lev2d)
  call conf%set('model.var2d',var2d_int)
  call conf%set('nicas.resolution',4.0_kind_real)
  call conf%set('nicas.explicit length-scales',.true.)
  rh(1) = fckit_configuration()
  call rh(1)%set('variables',(/'var1','var2'/))
  call rh(1)%set('value',1000.0e3_kind_real)
  call conf%set('nicas.horizontal length-scale',rh)
  rv(1) = fckit_configuration()
  call rv(1)%set('variables',(/'var1','var2'/))
  call rv(1)%set('value',3.0_kind_real)
  call conf%set('nicas.vertical length-scale',rv)

  ! Create BUMP
  call bump%create(f_comm,fspace_out,fset,conf)

  ! Run drivers
  call bump%run_drivers()

  ! Release memory (partial)
  call bump%partial_dealloc()

  ! Get control variable size
  call bump%get_cv_size(n)

  ! Initialize control variables
  cv_1 = atlas_field("cv1", atlas_real(kind_real), (/n/))
  cv_2 = atlas_field("cv2", atlas_real(kind_real), (/n/))
  call cv_1%data(ptr_1)
  call cv_2%data(ptr_2)
  call bump%rng%rand_gau(ptr_1)
  ptr_2 = ptr_1

  ! Create output fieldset
  fspace_out_sc = atlas_functionspace_structuredcolumns(fspace_out%c_ptr())
  nmga_out = fspace_out_sc%size_owned()
  allocate(gmask_out(nmga_out,nl0))
  gmask_out = .true.
  call fset_out_1%init(bump%mpl,fspace_out,gmask_out,variables,lev2d,var2d)
  call fset_out_2%init(bump%mpl,fspace_out,gmask_out,variables,lev2d,var2d)

  ! Initialize output fieldset
  allocate(array_out_1(nmga_out,nl0,2))
  allocate(array_out_2(nmga_out,nl0,2))
  call bump%rng%rand_gau(array_out_1)
  array_out_2 = array_out_1
  call fset_out_1%from_array(bump%mpl,array_out_1)
  call fset_out_2%from_array(bump%mpl,array_out_2)

  ! Apply NICAS square-root
  call bump%apply_nicas_sqrt(cv_1,fset_out_1,0)

  ! Apply adjoint NICAS square-root
  call bump%apply_nicas_sqrt_ad(fset_out_2,cv_2,0)

  ! Adjoint test
  call fset_out_1%to_array(bump%mpl,array_out_1)
  dp_in = sum(ptr_1*ptr_2)
  dp_out = sum(array_out_1*array_out_2)
  call bump%mpl%f_comm%allreduce(dp_in,fckit_mpi_sum())
  call bump%mpl%f_comm%allreduce(dp_out,fckit_mpi_sum())
  FCTEST_CHECK_CLOSE(dp_in,dp_out,repro_th)

  ! Release memory
  call bump%dealloc()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
