!----------------------------------------------------------------------
! Module: type_nicas
! Purpose: NICAS data derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_nicas

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_min
use netcdf
use tools_const, only: rad2deg,reqkm,pi
use tools_func, only: sphere_dist,cholesky,fit_diag
use tools_kinds, only: kind_real,nc_kind_real
use tools_qsort, only: qsort
use type_bpar, only: bpar_type
use type_cmat, only: cmat_type
use type_com, only: com_type,com_ntag
use type_cv, only: cv_type
use type_diag, only: diag_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_hdiag, only: hdiag_type
use type_io, only: io_type
use type_linop, only: linop_type,linop_ntag
use type_nicas_blk, only: nicas_blk_type,nicas_blk_ntag
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

integer,parameter :: nfac_rnd = 9 ! Number of ensemble size factors for randomization
integer,parameter :: nfac_opt = 4 ! Number of length-scale factors for optimization
integer,parameter :: ntest = 50   ! Number of tests


! NICAS derived type
type nicas_type
   character(len=1024) :: prefix              ! Prefix
   type(nicas_blk_type),allocatable :: blk(:) ! NICAS data blocks
   logical :: allocated                       ! Allocation flag
contains
   procedure :: alloc => nicas_alloc
   procedure :: partial_dealloc => nicas_partial_dealloc
   procedure :: dealloc => nicas_dealloc
   procedure :: read => nicas_read
   procedure :: write => nicas_write
   procedure :: run_nicas => nicas_run_nicas
   procedure :: run_nicas_tests => nicas_run_nicas_tests
   procedure :: alloc_cv => nicas_alloc_cv
   procedure :: random_cv => nicas_random_cv
   procedure :: apply => nicas_apply
   procedure :: apply_from_sqrt => nicas_apply_from_sqrt
   procedure :: apply_sqrt => nicas_apply_sqrt
   procedure :: apply_sqrt_ad => nicas_apply_sqrt_ad
   procedure :: randomize => nicas_randomize
   procedure :: apply_bens => nicas_apply_bens
   procedure :: test_adjoint => nicas_test_adjoint
   procedure :: test_dirac => nicas_test_dirac
   procedure :: test_randomization => nicas_test_randomization
   procedure :: test_consistency => nicas_test_consistency
   procedure :: test_optimality => nicas_test_optimality
end type nicas_type

private
public :: nicas_type

contains

!----------------------------------------------------------------------
! Subroutine: nicas_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine nicas_alloc(nicas,mpl,nam,bpar,prefix)

! Passed variables
class(nicas_type),intent(inout) :: nicas ! NICAS data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist
type(bpar_type),intent(in) :: bpar       ! Block parameters
character(len=*),intent(in) :: prefix    ! Prefix

! Local variable
integer :: ib

! Copy prefix
nicas%prefix = prefix

! Allocation
allocate(nicas%blk(bpar%nbe))

do ib=1,bpar%nbe
   ! Set block index
   nicas%blk(ib)%ib = ib

   ! Set number of communication steps
   nicas%blk(ib)%mpicom = nam%mpicom

   ! Set subsampling structure
   nicas%blk(ib)%subsamp = trim(nam%subsamp)

   ! Verbosity flag
   nicas%blk(ib)%verbosity = .true.

   ! Smoother flag
   nicas%blk(ib)%smoother = .false.

   ! Set name
   if (nam%lsqrt) then
      write(nicas%blk(ib)%name,'(a,i1,a,i4.4,a,i4.4,a,a)') trim(prefix)//'-',nam%mpicom,'-sqrt_',mpl%nproc,'-', &
    & mpl%myproc,'_',trim(bpar%blockname(ib))
   else
      write(nicas%blk(ib)%name,'(a,i1,a,i4.4,a,i4.4,a,a)') trim(prefix)//'-',nam%mpicom,'_',mpl%nproc,'-',mpl%myproc, &
    & '_',trim(bpar%blockname(ib))
   end if
end do

! Update allocation flag
nicas%allocated = .true.

end subroutine nicas_alloc

!----------------------------------------------------------------------
! Subroutine: nicas_partial_dealloc
! Purpose: release memory (partial)
!----------------------------------------------------------------------
subroutine nicas_partial_dealloc(nicas)

! Passed variables
class(nicas_type),intent(inout) :: nicas ! NICAS data

! Local variable
integer :: ib

! Release memory
if (allocated(nicas%blk)) then
   do ib=1,size(nicas%blk)
      call nicas%blk(ib)%partial_dealloc
   end do
end if

end subroutine nicas_partial_dealloc

!----------------------------------------------------------------------
! Subroutine: nicas_dealloc
! Purpose: release memory (full)
!----------------------------------------------------------------------
subroutine nicas_dealloc(nicas)

! Passed variables
class(nicas_type),intent(inout) :: nicas ! NICAS data

! Local variable
integer :: ib

! Release memory
if (allocated(nicas%blk)) then
   do ib=1,size(nicas%blk)
      call nicas%blk(ib)%dealloc
   end do
   deallocate(nicas%blk)
end if

! Update allocation flag
nicas%allocated = .false.

end subroutine nicas_dealloc

!----------------------------------------------------------------------
! Subroutine: nicas_read
! Purpose: read
!----------------------------------------------------------------------
subroutine nicas_read(nicas,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(inout) :: nicas ! NICAS data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(nam_type),intent(in) :: nam         ! Namelist
type(geom_type),intent(in) :: geom       ! Geometry
type(bpar_type),intent(in) :: bpar       ! Block parameters

! Local variables
integer :: ib,iproc,iprocio
type(nicas_type) :: nicas_tmp

! Allocation
call nicas%alloc(mpl,nam,bpar,'nicas')

! Write NICAS blocks
do ib=1,bpar%nbe
   if (bpar%B_block(ib)) then
      do iproc=1,mpl%nproc
         ! Reading task
         iprocio = mod(iproc,nam%nprocio)
         if (iprocio==0) iprocio = nam%nprocio

         if (mpl%myproc==iprocio) then
            if (iproc==iprocio) then
               write(mpl%info,'(a7,a,i4,a)') '','NICAS data of task ',iproc,' read by its own task'
               call mpl%flush

               ! Read data
               call nicas%blk(ib)%read(mpl,nam,geom,bpar)
            else
               write(mpl%info,'(a7,a,i4,a,i4,a)') '','NICAS data of task ',iproc,' send from task ',iprocio,' after reading'
               call mpl%flush

               ! Allocation
               call nicas_tmp%alloc(mpl,nam,bpar,'nicas')

               ! Set name
               if (nam%lsqrt) then
                  write(nicas_tmp%blk(ib)%name,'(a,i1,a,i4.4,a,i4.4,a,a)') trim(nicas_tmp%prefix)//'-',nam%mpicom,'-sqrt_', &
                & mpl%nproc,'-',iproc,'_',trim(bpar%blockname(ib))
               else
                  write(nicas_tmp%blk(ib)%name,'(a,i1,a,i4.4,a,i4.4,a,a)') trim(nicas_tmp%prefix)//'-',nam%mpicom,'_', &
                & mpl%nproc,'-',iproc,'_',trim(bpar%blockname(ib))
               end if

               ! Read data
               call nicas_tmp%blk(ib)%read(mpl,nam,geom,bpar)

               ! Send data to iproc
               call nicas_tmp%blk(ib)%send(mpl,nam,geom,bpar,iproc,mpl%tag)

               ! Release memory
               call nicas_tmp%dealloc
            end if
         elseif (mpl%myproc==iproc) then
            write(mpl%info,'(a7,a,i4,a,i4,a)') '','NICAS data of task ',iproc,' received from task ',iprocio,' after reading'
            call mpl%flush

            ! Receive data from iprocio
            call nicas%blk(ib)%receive(mpl,nam,geom,bpar,iprocio,mpl%tag)
         end if

         ! Update tag
         call mpl%update_tag(nicas_blk_ntag+4*com_ntag+(2+geom%nl0i+geom%nl0+(nam%nts-1)*geom%nl0)*linop_ntag)
      end do
   end if
end do

end subroutine nicas_read

!----------------------------------------------------------------------
! Subroutine: nicas_write
! Purpose: write
!----------------------------------------------------------------------
subroutine nicas_write(nicas,mpl,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas ! NICAS data
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(bpar_type),intent(in) :: bpar    ! Block parameters

! Local variables
integer :: ib,iproc,iprocio
type(nicas_type) :: nicas_tmp

! Write NICAS blocks
do ib=1,bpar%nbe
   if (bpar%B_block(ib)) then
      do iproc=1,mpl%nproc
         ! Writing task
         iprocio = mod(iproc,nam%nprocio)
         if (iprocio==0) iprocio = nam%nprocio
         if (mpl%myproc==iprocio) then
            if (iproc==iprocio) then
               write(mpl%info,'(a7,a,i4,a)') '','NICAS data of task ',iproc,' written by its own task'
               call mpl%flush

               ! Write data
               call nicas%blk(ib)%write(mpl,nam,geom,bpar)

               ! Write grids
               if (nam%write_grids) call nicas%blk(ib)%write_grids(mpl,nam,bpar)
            else
               write(mpl%info,'(a7,a,i4,a,i4,a)') '','NICAS data of task ',iproc,' received from task ',iprocio,' to be written'
               call mpl%flush

               ! Allocation
               call nicas_tmp%alloc(mpl,nam,bpar,'nicas')

               ! Set name
               if (nam%lsqrt) then
                  write(nicas_tmp%blk(ib)%name,'(a,i1,a,i4.4,a,i4.4,a,a)') trim(nicas_tmp%prefix)//'-',nam%mpicom,'-sqrt_', &
                & mpl%nproc,'-',iproc,'_',trim(bpar%blockname(ib))
               else
                  write(nicas_tmp%blk(ib)%name,'(a,i1,a,i4.4,a,i4.4,a,a)') trim(nicas_tmp%prefix)//'-',nam%mpicom,'_', &
                & mpl%nproc,'-',iproc,'_',trim(bpar%blockname(ib))
               end if

               ! Receive data from iproc
               call nicas_tmp%blk(ib)%receive(mpl,nam,geom,bpar,iproc,mpl%tag)

               ! Write data
               call nicas_tmp%blk(ib)%write(mpl,nam,geom,bpar)

               ! Write grids
               if (nam%write_grids) call nicas_tmp%blk(ib)%write_grids(mpl,nam,bpar)

               ! Release memory
               call nicas_tmp%dealloc
            end if
         elseif (mpl%myproc==iproc) then
            write(mpl%info,'(a7,a,i4,a,i4,a)') '','NICAS data of task ',iproc,' send to task ',iprocio,' for writing'
            call mpl%flush

            ! Send data to iprocio
            call nicas%blk(ib)%send(mpl,nam,geom,bpar,iprocio,mpl%tag)
         end if

         ! Update tag
         call mpl%update_tag(nicas_blk_ntag+4*com_ntag+(2+geom%nl0i+geom%nl0+(nam%nts-1)*geom%nl0)*linop_ntag)
      end do
   end if
end do

end subroutine nicas_write

!----------------------------------------------------------------------
! Subroutine: nicas_run_nicas
! Purpose: NICAS driver
!----------------------------------------------------------------------
subroutine nicas_run_nicas(nicas,mpl,rng,nam,geom,bpar,cmat)

implicit none

! Passed variables
class(nicas_type),intent(inout) :: nicas  ! NICAS data
type(mpl_type),intent(inout) :: mpl       ! MPI data
type(rng_type),intent(inout) :: rng       ! Random number generator
type(nam_type),intent(inout) :: nam       ! Namelist
type(geom_type),intent(inout) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar        ! Block parameters
type(cmat_type),intent(in) :: cmat        ! C matrix data

! Local variables
integer :: ib

! Allocation
call nicas%alloc(mpl,nam,bpar,'nicas')

! Compute NICAS parameters
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Compute NICAS parameters'
call mpl%flush

do ib=1,bpar%nbe
   if (bpar%nicas_block(ib).or.((ib==bpar%nbe).and.nam%adv_diag)) then
      write(mpl%info,'(a)') '-------------------------------------------------------------------'
      call mpl%flush
      write(mpl%info,'(a)') '--- Block: '//trim(bpar%blockname(ib))
      call mpl%flush
   end if

   ! NICAS parameters
   if (bpar%nicas_block(ib)) call nicas%blk(ib)%compute_parameters(mpl,rng,nam,geom,cmat%blk(ib))

   ! Advection
   if ((ib==bpar%nbe).and.nam%adv_diag) call nicas%blk(ib)%compute_adv(mpl,rng,nam,geom,cmat%blk(ib))

   ! Coefficient
   if (bpar%B_block(ib)) then
      ! Copy weights
      nicas%blk(ib)%wgt = cmat%blk(ib)%wgt
      if (bpar%nicas_block(ib)) then
         allocate(nicas%blk(ib)%coef_ens(geom%nc0a,geom%nl0))
         nicas%blk(ib)%coef_ens = cmat%blk(ib)%coef_ens
      end if
   end if
end do

if (nam%write_nicas) then
   ! Write NICAS parameters
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Write NICAS parameters'
   call mpl%flush
   call nicas%write(mpl,nam,geom,bpar)
end if

end subroutine nicas_run_nicas

!----------------------------------------------------------------------
! Subroutine: nicas_run_nicas_tests
! Purpose: NICAS tests driver
!----------------------------------------------------------------------
subroutine nicas_run_nicas_tests(nicas,mpl,rng,nam,geom,bpar,io,ens)

implicit none

! Passed variables
class(nicas_type),intent(inout) :: nicas  ! NICAS data
type(mpl_type),intent(inout) :: mpl       ! MPI data
type(rng_type),intent(inout) :: rng       ! Random number generator
type(nam_type),intent(inout) :: nam       ! Namelist
type(geom_type),intent(inout) :: geom     ! Geometry
type(bpar_type),intent(in) :: bpar        ! Block parameters
type(io_type),intent(in) :: io            ! I/O
type(ens_type),intent(in) :: ens          ! Ensemble

! Local variables
integer :: ib

if (nam%check_adjoints) then
   ! Test adjoint
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Test NICAS adjoint'
   call mpl%flush

   do ib=1,bpar%nbe
      if (bpar%nicas_block(ib)) then
         write(mpl%info,'(a)') '-------------------------------------------------------------------'
         call mpl%flush
         write(mpl%info,'(a)') '--- Block: '//trim(bpar%blockname(ib))
         call mpl%flush
         call nicas%blk(ib)%test_adjoint(mpl,rng,nam,geom)
      end if
   end do

   ! Test NICAS adjoint
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Test NICAS adjoint'
   call mpl%flush
   call nicas%test_adjoint(mpl,rng,nam,geom,bpar,ens)
end if

if (nam%check_dirac) then
   ! Apply NICAS to diracs
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Apply NICAS to diracs'
   call mpl%flush

   do ib=1,bpar%nbe
      if (bpar%nicas_block(ib)) then
         write(mpl%info,'(a)') '-------------------------------------------------------------------'
         call mpl%flush
         write(mpl%info,'(a)') '--- Block: '//trim(bpar%blockname(ib))
         call mpl%flush
         call nicas%blk(ib)%test_dirac(mpl,nam,geom,bpar,io)
      end if
   end do

   ! Apply NICAS to diracs
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Apply NICAS to diracs'
   call mpl%flush
   call nicas%test_dirac(mpl,nam,geom,bpar,io,ens)
end if

if (nam%check_randomization) then
   ! Test NICAS randomization
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Test NICAS randomization'
   call mpl%flush
   call nicas%test_randomization(mpl,rng,nam,geom,bpar)
end if

if (nam%check_consistency) then
   ! Test HDIAG-NICAS consistency
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Test HDIAG-NICAS consistency'
   call mpl%flush
   call nicas%test_consistency(mpl,rng,nam,geom,bpar,io)
end if

if (nam%check_optimality) then
   ! Test HDIAG optimality
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Test HDIAG optimality'
   call mpl%flush
   call nicas%test_optimality(mpl,rng,nam,geom,bpar,io)
end if

end subroutine nicas_run_nicas_tests

!----------------------------------------------------------------------
! Subroutine: nicas_alloc_cv
! Purpose: allocation
!----------------------------------------------------------------------
subroutine nicas_alloc_cv(nicas,mpl,bpar,cv,getsizeonly)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas      ! NICAS data
type(mpl_type),intent(inout) :: mpl        ! MPI data
type(bpar_type),intent(in) :: bpar         ! Block parameters
type(cv_type),intent(inout) :: cv          ! Control vector
logical,intent(in),optional :: getsizeonly ! Flag to get the control variable size only (no allocation)

! Local variables
integer :: ib
logical :: lgetsizeonly

! Check flag existence
lgetsizeonly = .false.
if (present(getsizeonly)) lgetsizeonly = getsizeonly

! Allocation
allocate(cv%blk(bpar%nbe))

! Initialization
cv%n = 0
cv%nbe = bpar%nbe

do ib=1,bpar%nbe
   if (mpl%msv%isnot(bpar%cv_block(ib))) then
      ! Copy block size
      cv%blk(ib)%n = nicas%blk(bpar%cv_block(ib))%nsa

      ! Update total size
      cv%n = cv%n+nicas%blk(bpar%cv_block(ib))%nsa

      if (.not.lgetsizeonly) then
         ! Allocation
         allocate(cv%blk(ib)%alpha(cv%blk(ib)%n))

         ! Initialization
         cv%blk(ib)%alpha = mpl%msv%valr
      end if
   else
      ! Set zero size
      cv%blk(ib)%n = 0
   end if
end do

end subroutine nicas_alloc_cv

!----------------------------------------------------------------------
! Subroutine: nicas_random_cv
! Purpose: generate a random control vector
!----------------------------------------------------------------------
subroutine nicas_random_cv(nicas,mpl,rng,bpar,cv)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas ! NICAS data
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(rng_type),intent(inout) :: rng   ! Random number generator
type(bpar_type),intent(in) :: bpar    ! Block parameters
type(cv_type),intent(out) :: cv       ! Control vector

! Local variables
integer :: ib,jb,isa,is
integer,allocatable :: order(:)
real(kind_real),allocatable :: hash_s(:),alpha(:)

! Allocation
call nicas%alloc_cv(mpl,bpar,cv)

! Resynchronize random number generator
call rng%resync(mpl)

! Random initialization
do ib=1,bpar%nbe
   ! CV block
   jb = bpar%cv_block(ib)

   if (mpl%msv%isnot(jb)) then
      ! Allocation
      allocate(hash_s(nicas%blk(jb)%ns))
      allocate(order(nicas%blk(jb)%ns))
      allocate(alpha(nicas%blk(jb)%ns))

      ! Random vector
      call rng%rand_gau(alpha)

      ! Get global hash
      call mpl%loc_to_glb(nicas%blk(jb)%nsa,nicas%blk(jb)%ns,nicas%blk(jb)%sa_to_s,nicas%blk(jb)%hash_sa,hash_s,.true.)

      ! Reorder random vector
      call qsort(nicas%blk(jb)%ns,hash_s,order)
      alpha = alpha(order)

      ! Copy local section     
      do isa=1,nicas%blk(jb)%nsa
         is = nicas%blk(jb)%sa_to_s(isa)
         cv%blk(ib)%alpha(isa) = alpha(is)
      end do

      ! Release memory
      deallocate(hash_s)
      deallocate(order)
      deallocate(alpha)
   end if
end do

! Desynchronize random number generator
call rng%desync(mpl)

end subroutine nicas_random_cv

!----------------------------------------------------------------------
! Subroutine: nicas_apply
! Purpose: apply NICAS
!----------------------------------------------------------------------
subroutine nicas_apply(nicas,mpl,nam,geom,bpar,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                   ! NICAS data
type(mpl_type),intent(inout) :: mpl                                     ! MPI data
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
type(bpar_type),intent(in) :: bpar                                      ! Block parameters
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Field

! Local variable
integer :: ib,its,iv,jv,il0,ic0a
real(kind_real) :: prod,prod_tot
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:)
real(kind_real),allocatable :: fld_save(:,:,:,:)
character(len=1024),parameter :: subr = 'nicas_apply'

if (nam%pos_def_test) then
   ! Save field for positive-definiteness test
   allocate(fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   fld_save = fld
end if

! Adjoint advection
if (nam%adv_mode==1) call nicas%blk(bpar%nbe)%apply_adv_ad(mpl,nam,geom,fld)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Sum product over variables and timeslots
   fld_3d = 0.0
   !$omp parallel do schedule(static) private(il0,ic0a,its,iv)
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         do its=1,nam%nts
            do iv=1,nam%nv
               fld_3d(ic0a,il0) = fld_3d(ic0a,il0)+fld(ic0a,il0,iv,its)
            end do
         end do
      end do
   end do
   !$omp end parallel do

   if (nam%nonunit_diag) then
      ! Apply common ensemble coefficient square-root
      !$omp parallel do schedule(static) private(il0,ic0a)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%gmask_c0a(ic0a,il0)) fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
         end do
      end do
      !$omp end parallel do
   end if

   ! Apply common NICAS
   call nicas%blk(bpar%nbe)%apply(mpl,geom,fld_3d)

   if (nam%nonunit_diag) then
      ! Apply common ensemble coefficient square-root
      !$omp parallel do schedule(static) private(il0,ic0a)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%gmask_c0a(ic0a,il0)) fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
         end do
      end do
      !$omp end parallel do
   end if

   ! Build final vector
   !$omp parallel do schedule(static) private(il0,ic0a,its,iv)
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         do its=1,nam%nts
            do iv=1,nam%nv
               fld(ic0a,il0,iv,its) = fld_3d(ic0a,il0)
            end do
         end do
      end do
   end do
   !$omp end parallel do

   ! Release memory
   deallocate(fld_3d)
case ('common_univariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   do iv=1,nam%nv
      if (nam%nonunit_diag) then
         ! Apply common ensemble coefficient square-root
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%gmask_c0a(ic0a,il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do
      end if

      ! Apply common NICAS
      call nicas%blk(bpar%nbe)%apply(mpl,geom,fld_4d(:,:,iv))

      if (nam%nonunit_diag) then
         ! Apply common ensemble coefficient square-root
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%gmask_c0a(ic0a,il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do

   ! Release memory
   deallocate(fld_4d)
case ('common_weighted')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   allocate(fld_4d_tmp(geom%nc0a,geom%nl0,nam%nv))
   allocate(wgt(nam%nv,nam%nv))
   allocate(wgt_diag(nam%nv))

   ! Copy weights
   wgt = 0.0
   wgt_diag = 0.0
   do ib=1,bpar%nb
      if (bpar%B_block(ib)) then
         ! Variable indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         wgt(iv,jv) = nicas%blk(ib)%wgt
         if (iv==jv) wgt_diag(iv) = wgt(iv,iv)
      end if
   end do

   ! Normalize weights
   do iv=1,nam%nv
      do jv=1,nam%nv
         wgt(iv,jv) = wgt(iv,jv)/sqrt(wgt_diag(iv)*wgt_diag(jv))
      end do
   end do

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld(:,:,:,its)
   end do

   do iv=1,nam%nv
      if (nam%nonunit_diag) then
         ! Apply common ensemble coefficient square-root
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%gmask_c0a(ic0a,il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do
      end if

      ! Apply common NICAS
      call nicas%blk(bpar%nbe)%apply(mpl,geom,fld_4d(:,:,iv))
      if (nam%nonunit_diag) then
         ! Apply common ensemble coefficient square-root
         !$omp parallel do schedule(static) private(il0,ic0a)
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%gmask_c0a(ic0a,il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
            end do
         end do
         !$omp end parallel do
      end if
   end do

   ! Apply weights
   fld_4d_tmp = 0.0
   do iv=1,nam%nv
      do jv=1,nam%nv
         fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)+wgt(iv,jv)*fld_4d(:,:,jv)
      end do
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d_tmp
   end do

   ! Release memory
   deallocate(fld_4d)
   deallocate(fld_4d_tmp)
   deallocate(wgt)
   deallocate(wgt_diag)
case ('specific_univariate')
   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)
         its = bpar%b_to_ts1(ib)

         if (nam%nonunit_diag) then
            ! Apply common ensemble coefficient square-root
            !$omp parallel do schedule(static) private(il0,ic0a)
            do il0=1,geom%nl0
               do ic0a=1,geom%nc0a
                  if (geom%gmask_c0a(ic0a,il0)) fld(ic0a,il0,iv,its) = fld(ic0a,il0,iv,its) &
                                                                    & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
               end do
            end do
            !$omp end parallel do
         end if

         ! Apply specific NICAS
         call nicas%blk(ib)%apply(mpl,geom,fld(:,:,iv,its))

         if (nam%nonunit_diag) then
            ! Apply common ensemble coefficient square-root
            !$omp parallel do schedule(static) private(il0,ic0a)
            do il0=1,geom%nl0
               do ic0a=1,geom%nc0a
                  if (geom%gmask_c0a(ic0a,il0)) fld(ic0a,il0,iv,its) = fld(ic0a,il0,iv,its) &
                                                                    & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
               end do
            end do
            !$omp end parallel do
         end if
      end if
   end do
case ('specific_multivariate')
   call mpl%abort(subr,'specific multivariate strategy should not be called from apply_NICAS (lsqrt required)')
end select

if (nam%pos_def_test) then
   ! Positive-definiteness test
   prod = sum(fld_save*fld)
   call mpl%f_comm%allreduce(prod,prod_tot,fckit_mpi_sum())
   if (prod_tot<0.0) call mpl%abort(subr,'negative result in nicas_apply')

   ! Release memory
   deallocate(fld_save)
end if

! Advection
if (nam%adv_mode==1) call nicas%blk(bpar%nbe)%apply_adv(mpl,nam,geom,fld)

end subroutine nicas_apply

!----------------------------------------------------------------------
! Subroutine: nicas_apply_from_sqrt
! Purpose: apply NICAS from square-root
!----------------------------------------------------------------------
subroutine nicas_apply_from_sqrt(nicas,mpl,nam,geom,bpar,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                   ! NICAS data
type(mpl_type),intent(inout) :: mpl                                     ! MPI data
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
type(bpar_type),intent(in) :: bpar                                      ! Block parameters
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Field

! Local variable
real(kind_real) :: prod,prod_tot
real(kind_real),allocatable :: fld_save(:,:,:,:)
character(len=1024),parameter :: subr = 'nicas_apply_from_sqrt'
type(cv_type) :: cv

if (nam%pos_def_test) then
   ! Save field for positivity test
   allocate(fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   fld_save = fld
end if

! Apply square-root adjoint
call nicas%apply_sqrt_ad(mpl,nam,geom,bpar,fld,cv)

! Apply square-root
call nicas%apply_sqrt(mpl,nam,geom,bpar,cv,fld)

if (nam%pos_def_test) then
   ! Positivity test
   prod = sum(fld_save*fld)
   call mpl%f_comm%allreduce(prod,prod_tot,fckit_mpi_sum())
   if (prod_tot<0.0) call mpl%abort(subr,'negative result in nicas_apply')

   ! Release memory
   deallocate(fld_save)
end if

end subroutine nicas_apply_from_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_apply_sqrt
! Purpose: apply NICAS square-root
!----------------------------------------------------------------------
subroutine nicas_apply_sqrt(nicas,mpl,nam,geom,bpar,cv,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                 ! NICAS data
type(mpl_type),intent(inout) :: mpl                                   ! MPI data
type(nam_type),intent(in) :: nam                                      ! Namelist
type(geom_type),intent(in) :: geom                                    ! Geometry
type(bpar_type),intent(in) :: bpar                                    ! Block parameters
type(cv_type),intent(in) :: cv                                        ! Control variable
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Field

! Local variable
integer :: ib,its,iv,jv,ic0a,il0,ierr
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:),wgt_u(:,:)
character(len=1024),parameter :: subr = 'nicas_apply_sqrt'

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Apply common NICAS
   call nicas%blk(bpar%nbe)%apply_sqrt(mpl,geom,cv%blk(bpar%nbe)%alpha,fld_3d)

   if (nam%nonunit_diag) then
      ! Apply common ensemble coefficient square-root
      !$omp parallel do schedule(static) private(il0,ic0a)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%gmask_c0a(ic0a,il0)) fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
         end do
      end do
     !$omp end parallel do
   end if

   ! Build final vector
   do its=1,nam%nts
      do iv=1,nam%nv
         fld(:,:,iv,its) = fld_3d
      end do
   end do

   ! Release memory
   deallocate(fld_3d)
case ('common_univariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   do ib=1,bpar%nb
      if (mpl%msv%isnot(bpar%cv_block(ib))) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply specific NICAS
         call nicas%blk(bpar%nbe)%apply_sqrt(mpl,geom,cv%blk(ib)%alpha,fld_4d(:,:,iv))

         if (nam%nonunit_diag) then
            ! Apply common ensemble coefficient square-root
            !$omp parallel do schedule(static) private(il0,ic0a)
            do il0=1,geom%nl0
               do ic0a=1,geom%nc0a
                  if (geom%gmask_c0a(ic0a,il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                   & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
               end do
            end do
            !$omp end parallel do
         end if
      end if
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d
   end do

   ! Release memory
   deallocate(fld_4d)
case ('common_weighted')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   allocate(fld_4d_tmp(geom%nc0a,geom%nl0,nam%nv))
   allocate(wgt(nam%nv,nam%nv))
   allocate(wgt_diag(nam%nv))
   allocate(wgt_u(nam%nv,nam%nv))

   ! Copy weights
   wgt = 0.0
   wgt_diag = 0.0
   do ib=1,bpar%nb
      if (bpar%B_block(ib)) then
         ! Variable indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         wgt(iv,jv) = nicas%blk(ib)%wgt
         if (iv==jv) wgt_diag(iv) = wgt(iv,iv)
      end if
   end do

   ! Normalize weights
   do iv=1,nam%nv
      do jv=1,nam%nv
         wgt(iv,jv) = wgt(iv,jv)/sqrt(wgt_diag(iv)*wgt_diag(jv))
      end do
   end do

   ! Cholesky decomposition
   call cholesky(mpl,nam%nv,wgt,wgt_u,ierr)
   if (ierr/=0) call mpl%abort(subr,'matrix is not positive semi-definite in Cholesky decomposition')

   do ib=1,bpar%nb
      if (mpl%msv%isnot(bpar%cv_block(ib))) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         ! Apply specific NICAS
         call nicas%blk(bpar%nbe)%apply_sqrt(mpl,geom,cv%blk(ib)%alpha,fld_4d(:,:,iv))

         if (nam%nonunit_diag) then
            ! Apply common ensemble coefficient square-root
            !$omp parallel do schedule(static) private(il0,ic0a)
            do il0=1,geom%nl0
               do ic0a=1,geom%nc0a
                  if (geom%gmask_c0a(ic0a,il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                   & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
               end do
            end do
            !$omp end parallel do
         end if
      end if
   end do

   ! Apply weights
   fld_4d_tmp = 0.0
   do iv=1,nam%nv
      do jv=1,iv
         fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)+wgt_u(iv,jv)*fld_4d(:,:,jv)
      end do
   end do

   ! Build final vector
   do its=1,nam%nts
      fld(:,:,:,its) = fld_4d_tmp
   end do

   ! Release memory
   deallocate(fld_4d)
   deallocate(fld_4d_tmp)
   deallocate(wgt)
   deallocate(wgt_diag)
   deallocate(wgt_u)
case ('specific_univariate')
   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)
         its = bpar%b_to_ts1(ib)

         ! Apply specific NICAS
         call nicas%blk(ib)%apply_sqrt(mpl,geom,cv%blk(ib)%alpha,fld(:,:,iv,its))

         if (nam%nonunit_diag) then
            ! Apply specific ensemble coefficient square-root
            !$omp parallel do schedule(static) private(il0,ic0a)
            do il0=1,geom%nl0
               do ic0a=1,geom%nc0a
                  if (geom%gmask_c0a(ic0a,il0)) fld(ic0a,il0,iv,its) = fld(ic0a,il0,iv,its) &
                                                                       & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
               end do
            end do
            !$omp end parallel do
         end if
      end if
   end do
case ('specific_multivariate')
   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)
         its = bpar%b_to_ts1(ib)

         ! Apply specific NICAS
         call nicas%blk(ib)%apply_sqrt(mpl,geom,cv%blk(1)%alpha,fld(:,:,iv,its))

         if (nam%nonunit_diag) then
            ! Apply specific ensemble coefficient square-root
            !$omp parallel do schedule(static) private(il0,ic0a)
            do il0=1,geom%nl0
               do ic0a=1,geom%nc0a
                  if (geom%gmask_c0a(ic0a,il0)) fld(ic0a,il0,iv,its) = fld(ic0a,il0,iv,its) &
                                                                    & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
               end do
            end do
            !$omp end parallel do
         end if
      end if
   end do
end select

! Advection
if (nam%adv_mode==1) call nicas%blk(bpar%nbe)%apply_adv(mpl,nam,geom,fld)

end subroutine nicas_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: nicas_apply_sqrt_ad
! Purpose: apply NICAS square-root, adjoint
!----------------------------------------------------------------------
subroutine nicas_apply_sqrt_ad(nicas,mpl,nam,geom,bpar,fld,cv)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                ! NICAS data
type(mpl_type),intent(inout) :: mpl                                  ! MPI data
type(nam_type),intent(in) :: nam                                     ! Namelist
type(geom_type),intent(in) :: geom                                   ! Geometry
type(bpar_type),intent(in) :: bpar                                   ! Block parameters
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Field
type(cv_type),intent(out) :: cv                                      ! Control variable

! Local variable
integer :: ib,its,iv,jv,ic0a,il0,ierr
real(kind_real),allocatable :: fld_3d(:,:),fld_4d(:,:,:),fld_4d_tmp(:,:,:),fld_5d(:,:,:,:)
real(kind_real),allocatable :: wgt(:,:),wgt_diag(:),wgt_u(:,:)
type(cv_type) :: cv_tmp
character(len=1024),parameter :: subr = 'nicas_apply_sqrt_ad'

! Allocation
allocate(fld_5d(geom%nc0a,geom%nl0,nam%nv,nam%nts))
call nicas%alloc_cv(mpl,bpar,cv)

! Copy
fld_5d = fld

! Adjoint advection
if (nam%adv_mode==1) call nicas%blk(bpar%nbe)%apply_adv_ad(mpl,nam,geom,fld_5d)

select case (nam%strategy)
case ('common')
   ! Allocation
   allocate(fld_3d(geom%nc0a,geom%nl0))

   ! Sum product over variables and timeslots
   fld_3d = 0.0
   do its=1,nam%nts
      do iv=1,nam%nv
         fld_3d = fld_3d+fld_5d(:,:,iv,its)
      end do
   end do

   if (nam%nonunit_diag) then
      ! Apply common ensemble coefficient square-root
      !$omp parallel do schedule(static) private(il0,ic0a)
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%gmask_c0a(ic0a,il0)) fld_3d(ic0a,il0) = fld_3d(ic0a,il0)*sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
         end do
      end do
      !$omp end parallel do
   end if

   ! Apply common NICAS
   call nicas%blk(bpar%nbe)%apply_sqrt_ad(mpl,geom,fld_3d,cv%blk(bpar%nbe)%alpha)

   ! Release memory
   deallocate(fld_3d)
case ('common_univariate')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld_5d(:,:,:,its)
   end do

   do ib=1,bpar%nb
      if (mpl%msv%isnot(bpar%cv_block(ib))) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         if (nam%nonunit_diag) then
            ! Apply common ensemble coefficient square-root
            !$omp parallel do schedule(static) private(il0,ic0a)
            do il0=1,geom%nl0
               do ic0a=1,geom%nc0a
                  if (geom%gmask_c0a(ic0a,il0)) fld_4d(ic0a,il0,iv) = fld_4d(ic0a,il0,iv) &
                                                                   & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
               end do
            end do
            !$omp end parallel do
         end if

         ! Apply specific NICAS
         call nicas%blk(bpar%nbe)%apply_sqrt_ad(mpl,geom,fld_4d(:,:,iv),cv%blk(ib)%alpha)
      end if
   end do

   ! Release memory
   deallocate(fld_4d)
case ('common_weighted')
   ! Allocation
   allocate(fld_4d(geom%nc0a,geom%nl0,nam%nv))
   allocate(fld_4d_tmp(geom%nc0a,geom%nl0,nam%nv))
   allocate(wgt(nam%nv,nam%nv))
   allocate(wgt_diag(nam%nv))
   allocate(wgt_u(nam%nv,nam%nv))

   ! Copy weights
   wgt = 0.0
   wgt_diag = 0.0
   do ib=1,bpar%nb
      if (bpar%B_block(ib)) then
         ! Variable indices
         iv = bpar%b_to_v1(ib)
         jv = bpar%b_to_v2(ib)
         wgt(iv,jv) = nicas%blk(ib)%wgt
         if (iv==jv) wgt_diag(iv) = wgt(iv,iv)
      end if
   end do

   ! Normalize weights
   do iv=1,nam%nv
      do jv=1,nam%nv
         wgt(iv,jv) = wgt(iv,jv)/sqrt(wgt_diag(iv)*wgt_diag(jv))
      end do
   end do

   ! Cholesky decomposition
   call cholesky(mpl,nam%nv,wgt,wgt_u,ierr)
   if (ierr/=0) call mpl%abort(subr,'matrix is not positive semi-definite in Cholesky decomposition')

   ! Sum product over timeslots
   fld_4d = 0.0
   do its=1,nam%nts
      fld_4d = fld_4d+fld_5d(:,:,:,its)
   end do

   ! Apply weights
   fld_4d_tmp = 0.0
   do iv=1,nam%nv
      do jv=iv,nam%nv
         fld_4d_tmp(:,:,iv) = fld_4d_tmp(:,:,iv)+wgt_u(jv,iv)*fld_4d(:,:,jv)
      end do
   end do

   do ib=1,bpar%nb
      if (mpl%msv%isnot(bpar%cv_block(ib))) then
         ! Variable index
         iv = bpar%b_to_v1(ib)

         if (nam%nonunit_diag) then
            ! Apply common ensemble coefficient square-root
            !$omp parallel do schedule(static) private(il0,ic0a)
            do il0=1,geom%nl0
               do ic0a=1,geom%nc0a
                  if (geom%gmask_c0a(ic0a,il0)) fld_4d_tmp(ic0a,il0,iv) = fld_4d_tmp(ic0a,il0,iv) &
                                                                       & *sqrt(nicas%blk(bpar%nbe)%coef_ens(ic0a,il0))
               end do
            end do
            !$omp end parallel do
         end if

         ! Apply specific NICAS
         call nicas%blk(bpar%nbe)%apply_sqrt_ad(mpl,geom,fld_4d_tmp(:,:,iv),cv%blk(ib)%alpha)
      end if
   end do

   ! Release memory
   deallocate(fld_4d)
   deallocate(fld_4d_tmp)
   deallocate(wgt)
   deallocate(wgt_diag)
   deallocate(wgt_u)
case ('specific_univariate')
   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Indices
         iv = bpar%b_to_v1(ib)
         its = bpar%b_to_ts1(ib)

         if (nam%nonunit_diag) then
            ! Apply common ensemble coefficient square-root
            !$omp parallel do schedule(static) private(il0,ic0a)
            do il0=1,geom%nl0
               do ic0a=1,geom%nc0a
                  if (geom%gmask_c0a(ic0a,il0)) fld_5d(ic0a,il0,iv,its) = fld_5d(ic0a,il0,iv,its) &
                                                                       & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
               end do
            end do
            !$omp end parallel do
         end if

         ! Apply specific NICAS
         call nicas%blk(ib)%apply_sqrt_ad(mpl,geom,fld_5d(:,:,iv,its),cv%blk(ib)%alpha)
      end if
   end do
case ('specific_multivariate')
   ! Allocation
   call nicas%alloc_cv(mpl,bpar,cv_tmp)

   ! Initialization
   cv%blk(1)%alpha = 0.0

   do ib=1,bpar%nb
      if (bpar%nicas_block(ib)) then
         ! Variable index
         iv = bpar%b_to_v1(ib)
         its = bpar%b_to_ts1(ib)

         if (nam%nonunit_diag) then
            ! Apply common ensemble coefficient square-root
            !$omp parallel do schedule(static) private(il0,ic0a)
            do il0=1,geom%nl0
               do ic0a=1,geom%nc0a
                  if (geom%gmask_c0a(ic0a,il0)) fld_5d(ic0a,il0,iv,its) = fld_5d(ic0a,il0,iv,its) &
                                                                       & *sqrt(nicas%blk(ib)%coef_ens(ic0a,il0))
               end do
            end do
            !$omp end parallel do
         end if

         ! Apply specific NICAS
         call nicas%blk(ib)%apply_sqrt_ad(mpl,geom,fld_5d(:,:,iv,its),cv_tmp%blk(1)%alpha)

         ! Sum control variable
         cv%blk(1)%alpha = cv%blk(1)%alpha+cv_tmp%blk(1)%alpha
      end if
   end do
end select

! Release memory
deallocate(fld_5d)

end subroutine nicas_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: nicas_randomize
! Purpose: randomize NICAS from square-root
!----------------------------------------------------------------------
subroutine nicas_randomize(nicas,mpl,rng,nam,geom,bpar,ne,ens)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas ! NICAS data
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(rng_type),intent(inout) :: rng   ! Random number generator
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(bpar_type),intent(in) :: bpar    ! Blocal parameters
integer,intent(in) :: ne              ! Number of members
type(ens_type),intent(out) :: ens     ! Ensemble

! Local variable
integer :: ie,ic0a,il0,its,iv
real(kind_real) :: std(geom%nc0a,geom%nl0,nam%nv,nam%nts)
type(cv_type) :: cv_ens(ne)

! Allocation
call ens%alloc(nam,geom,ne,1)

do ie=1,ne
   ! Generate random control vector
   call nicas%random_cv(mpl,rng,bpar,cv_ens(ie))

   ! Allocate member
   if (.not.allocated(ens%mem(ie)%fld)) allocate(ens%mem(ie)%fld(geom%nc0a,geom%nl0,nam%nv,nam%nts))

   ! Apply square-root
   call nicas%apply_sqrt(mpl,nam,geom,bpar,cv_ens(ie),ens%mem(ie)%fld)
end do

! Remove mean
call ens%remove_mean

! Compute standard deviation
!$omp parallel do schedule(static) private(its,iv,il0,ic0a,ie)
do its=1,nam%nts
   do iv=1,nam%nv
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if (geom%gmask_c0a(ic0a,il0)) then
               std(ic0a,il0,iv,its) = 0.0
               do ie=1,ne
                  std(ic0a,il0,iv,its) = std(ic0a,il0,iv,its)+ens%mem(ie)%fld(ic0a,il0,iv,its)**2
               end do
               std(ic0a,il0,iv,its) = sqrt(std(ic0a,il0,iv,its)/real(ne-1,kind_real))
            end if
         end do
      end do
   end do
end do
!$omp end parallel do

! Normalize perturbations
!$omp parallel do schedule(static) private(ie,its,iv,il0,ic0a)
do ie=1,ne
   do its=1,nam%nts
      do iv=1,nam%nv
         do il0=1,geom%nl0
            do ic0a=1,geom%nc0a
               if (geom%gmask_c0a(ic0a,il0)) ens%mem(ie)%fld(ic0a,il0,iv,its) = ens%mem(ie)%fld(ic0a,il0,iv,its) &
                                                                             & /std(ic0a,il0,iv,its)
            end do
         end do
      end do
   end do
end do
!$omp end parallel do

end subroutine nicas_randomize

!----------------------------------------------------------------------
! Subroutine: nicas_apply_bens
! Purpose: apply localized ensemble covariance
!----------------------------------------------------------------------
subroutine nicas_apply_bens(nicas,mpl,nam,geom,bpar,ens,fld)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas                                   ! NICAS data
type(mpl_type),intent(inout) :: mpl                                     ! MPI data
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
type(bpar_type),intent(in) :: bpar                                      ! Blocal parameters
type(ens_type),intent(in) :: ens                                        ! Ensemble
real(kind_real),intent(inout) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts) ! Field

! Local variable
integer :: ie
real(kind_real) :: fld_copy(geom%nc0a,geom%nl0,nam%nv,nam%nts),fld_tmp(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: pert(geom%nc0a,geom%nl0,nam%nv,nam%nts)

! Copy field
fld_copy = fld

! Adjoint advection
if (nam%adv_mode==-1) call nicas%blk(bpar%nbe)%apply_adv_ad(mpl,nam,geom,fld_copy)

! Apply localized ensemble covariance formula
fld = 0.0
do ie=1,ens%ne
   ! Compute perturbation
   pert = ens%mem(ie)%fld

   ! Inverse advection
   if (nam%adv_mode==-1) call nicas%blk(bpar%nbe)%apply_adv_inv(mpl,nam,geom,pert)

   ! Schur product
   fld_tmp = pert*fld_copy

   ! Apply NICAS
   if (nam%lsqrt) then
      call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld_tmp)
   else
      call nicas%apply(mpl,nam,geom,bpar,fld_tmp)
   end if

   ! Schur product
   fld = fld+fld_tmp*pert

   ! Normalization
   fld = fld/real(ens%ne-1,kind_real)
end do

! Advection
if (nam%adv_mode==-1) call nicas%blk(bpar%nbe)%apply_adv(mpl,nam,geom,fld)

end subroutine nicas_apply_bens

!----------------------------------------------------------------------
! Subroutine: nicas_test_adjoint
! Purpose: test NICAS adjoint
!----------------------------------------------------------------------
subroutine nicas_test_adjoint(nicas,mpl,rng,nam,geom,bpar,ens)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas     ! NICAS data
type(mpl_type),intent(inout) :: mpl       ! MPI data
type(rng_type),intent(inout) :: rng       ! Random number generator
type(nam_type),intent(in) :: nam          ! Namelist
type(geom_type),intent(in) :: geom        ! Geometry
type(bpar_type),intent(in) :: bpar        ! Block parameters
type(ens_type),intent(in) :: ens          ! Ensemble

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real) :: fld1_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts),fld1_save(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real) :: fld2_loc(geom%nc0a,geom%nl0,nam%nv,nam%nts),fld2_save(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real),allocatable :: fld1_adv(:,:,:,:),fld2_adv(:,:,:,:),fld1_bens(:,:,:,:),fld2_bens(:,:,:,:)

! Allocation
if (abs(nam%adv_mode)==1) then
   allocate(fld1_adv(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   allocate(fld2_adv(geom%nc0a,geom%nl0,nam%nv,nam%nts))
end if
if (ens%allocated) then
   allocate(fld1_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   allocate(fld2_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))
end if

! Generate random field
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld1_save)
call rng%rand_real(0.0_kind_real,1.0_kind_real,fld2_save)

! Adjoint test
fld1_loc = fld1_save
fld2_loc = fld2_save
if (nam%lsqrt) then
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld1_loc)
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld2_loc)
else
   call nicas%apply(mpl,nam,geom,bpar,fld1_loc)
   call nicas%apply(mpl,nam,geom,bpar,fld2_loc)
end if
if (abs(nam%adv_mode)==1) then
   fld1_adv = fld1_save
   fld2_adv = fld2_save
   call nicas%blk(bpar%nbe)%apply_adv(mpl,nam,geom,fld1_adv)
   call nicas%blk(bpar%nbe)%apply_adv_ad(mpl,nam,geom,fld2_adv)
end if
if (ens%allocated) then
   fld1_bens = fld1_save
   fld2_bens = fld2_save
   call nicas%apply_bens(mpl,nam,geom,bpar,ens,fld1_bens)
   call nicas%apply_bens(mpl,nam,geom,bpar,ens,fld2_bens)
end if

! Print result
call mpl%dot_prod(fld1_loc,fld2_save,sum1)
call mpl%dot_prod(fld2_loc,fld1_save,sum2)
write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','NICAS adjoint test:                       ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
call mpl%flush
if (abs(nam%adv_mode)==1) then
   call mpl%dot_prod(fld1_adv,fld2_save,sum1)
   call mpl%dot_prod(fld2_adv,fld1_save,sum2)
   write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Advection adjoint test:                   ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
   call mpl%flush
end if
if (ens%allocated) then
   call mpl%dot_prod(fld1_bens,fld2_save,sum1)
   call mpl%dot_prod(fld2_bens,fld1_save,sum2)
   write(mpl%info,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Ensemble B adjoint test:                  ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)
   call mpl%flush
end if

! Release memory
if (abs(nam%adv_mode)==1) then
   deallocate(fld1_adv)
   deallocate(fld2_adv)
end if
if (ens%allocated) then
   deallocate(fld1_bens)
   deallocate(fld2_bens)
end if

end subroutine nicas_test_adjoint

!----------------------------------------------------------------------
! Subroutine: nicas_test_dirac
! Purpose: apply NICAS to diracs
!----------------------------------------------------------------------
subroutine nicas_test_dirac(nicas,mpl,nam,geom,bpar,io,ens)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas     ! NICAS data
type(mpl_type),intent(inout) :: mpl       ! MPI data
type(nam_type),intent(in) :: nam          ! Namelist
type(geom_type),intent(in) :: geom        ! Geometry
type(bpar_type),intent(in) :: bpar        ! Block parameters
type(io_type),intent(in) :: io            ! I/O
type(ens_type),intent(in) :: ens          ! Ensemble

! Local variables
integer :: idir,iv,its
real(kind_real) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts)
real(kind_real),allocatable :: fld_bens(:,:,:,:)
character(len=2) :: itschar
character(len=1024) :: filename

! Generate dirac field
fld = 0.0
do idir=1,geom%ndir
   if (geom%iprocdir(idir)==mpl%myproc) fld(geom%ic0adir(idir),geom%il0dir(idir),geom%ivdir(idir),geom%itsdir(idir)) = 1.0
end do

! Allocation and initialization
if (ens%allocated.and.(trim(nam%method)/='cor')) then
   allocate(fld_bens(geom%nc0a,geom%nl0,nam%nv,nam%nts))
   fld_bens = fld
end if

! Apply NICAS to dirac
if (nam%lsqrt) then
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld)
else
   call nicas%apply(mpl,nam,geom,bpar,fld)
end if

! Apply localized ensemble covariance
if (ens%allocated.and.(trim(nam%method)/='cor')) call nicas%apply_bens(mpl,nam,geom,bpar,ens,fld_bens)

! Write field
filename = trim(nam%prefix)//'_dirac'
call io%fld_write(mpl,nam,geom,filename,'vunit',geom%vunit_c0a)
do its=1,nam%nts
   write(itschar,'(i2.2)') its
   do iv=1,nam%nv
      call io%fld_write(mpl,nam,geom,filename,'nicas_'//trim(nam%varname(iv))//'_'//itschar,fld(:,:,iv,its))
      if (ens%allocated.and.(trim(nam%method)/='cor')) call io%fld_write(mpl,nam,geom,filename, &
       & 'Bens_'//trim(nam%varname(iv))//'_'//itschar,fld_bens(:,:,iv,its))
   end do
end do

! Release memory
if (ens%allocated.and.(trim(nam%method)/='cor')) deallocate(fld_bens)

end subroutine nicas_test_dirac

!----------------------------------------------------------------------
! Subroutine: nicas_test_randomization
! Purpose: test NICAS randomization method with respect to theoretical error statistics
!----------------------------------------------------------------------
subroutine nicas_test_randomization(nicas,mpl,rng,nam,geom,bpar)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas ! NICAS data
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(rng_type),intent(inout) :: rng   ! Random number generator
type(nam_type),intent(inout) :: nam   ! Namelist variables
type(geom_type),intent(in) :: geom    ! Geometry
type(bpar_type),intent(in) :: bpar    ! Block parameters

! Local variables
integer :: ifac,itest,nefac(nfac_rnd),ens1_ne
integer :: ncid,ntest_id,nfac_id,mse_id,mse_th_id
real(kind_real) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts),mse(ntest,nfac_rnd),mse_th(ntest,nfac_rnd),mse_avg,mse_th_avg
real(kind_real),allocatable :: fld_ref(:,:,:,:,:),fld_save(:,:,:,:,:)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'nicas_test_randomization'
type(ens_type) :: ens

! Allocation
allocate(fld_ref(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest))
allocate(fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest))

! Define test vectors
write(mpl%info,'(a4,a)') '','Define test vectors'
call mpl%flush
call define_test_vectors(mpl,rng,nam,geom,ntest,fld_save)
if (nam%default_seed) call rng%reseed(mpl)

! Apply NICAS to test vectors
write(mpl%info,'(a4,a)') '','Apply NICAS to test vectors: '
call mpl%flush(.false.)
call mpl%prog_init(ntest)
fld_ref = fld_save
do itest=1,ntest
   ! Apply vector
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld_ref(:,:,:,:,itest))

   ! Update
   call mpl%prog_print(itest)
end do
call mpl%prog_final

! Save namelist variables
ens1_ne = nam%ens1_ne

write(mpl%info,'(a4,a)') '','Test randomization for various ensemble sizes:'
call mpl%flush
do ifac=1,nfac_rnd
   ! Ensemble size
   nefac(ifac) = max(int(real(ifac,kind_real)/real(nfac_rnd,kind_real)*real(ens1_ne,kind_real)),3)
   nam%ens1_ne = nefac(ifac)
   write(mpl%info,'(a7,a,i4,a)') '','Ensemble sizes: ',nefac(ifac),' members'
   call mpl%flush

   ! Randomize ensemble
   write(mpl%info,'(a10,a)') '','Randomization'
   call mpl%flush
   call nicas%randomize(mpl,rng,nam,geom,bpar,nefac(ifac),ens)

   ! Test randomized ensemble
   write(mpl%info,'(a10,a)') '','Apply NICAS to test vectors: '
   call mpl%flush(.false.)
   call mpl%prog_init(ntest)
   do itest=1,ntest
      ! Test NICAS
      fld = fld_save(:,:,:,:,itest)
      call ens%apply_bens(mpl,nam,geom,fld)

      ! RMSE
      fld = fld-fld_ref(:,:,:,:,itest)
      call mpl%dot_prod(fld,fld,mse(itest,ifac))
      call mpl%dot_prod(fld_ref(:,:,:,:,itest),fld_ref(:,:,:,:,itest),mse_th(itest,ifac))
      mse_th(itest,ifac) = 1.0/real(nam%ens1_ne-1,kind_real)*(mse_th(itest,ifac)+real(geom%nc0*geom%nl0*nam%nv*nam%nts,kind_real))

      ! Update
      call mpl%prog_print(itest)
   end do
   call mpl%prog_final

   ! Print scores
   mse_avg = sum(mse(:,ifac))/real(ntest,kind_real)
   mse_th_avg = sum(mse_th(:,ifac))/real(ntest,kind_real)
   write(mpl%info,'(a10,a,e15.8,a,e15.8,a,f5.3)') '','MSE (exp. / th. / ratio): ',mse_avg,' / ',mse_th_avg,' / ',mse_avg/mse_th_avg
   call mpl%flush

   ! Release memory
   call ens%dealloc
end do

! Reset namelist variables
nam%ens1_ne = ens1_ne

if (mpl%main) then
   ! Create file
   filename = trim(nam%prefix)//'_randomization'
   call mpl%ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

   ! Write namelist parameters
   call nam%write(mpl,ncid)

   ! Define dimensions
   call mpl%ncerr(subr,nf90_def_dim(ncid,'ntest',ntest,ntest_id))
   call mpl%ncerr(subr,nf90_def_dim(ncid,'nfac',nfac_rnd,nfac_id))

   ! Define variables
   call mpl%ncerr(subr,nf90_def_var(ncid,'mse',nc_kind_real,(/ntest_id,nfac_id/),mse_id))
   call mpl%ncerr(subr,nf90_def_var(ncid,'mse_th',nc_kind_real,(/ntest_id,nfac_id/),mse_th_id))

   ! End definition mode
   call mpl%ncerr(subr,nf90_enddef(ncid))

   ! Write variables
   call mpl%ncerr(subr,nf90_put_var(ncid,mse_id,mse))
   call mpl%ncerr(subr,nf90_put_var(ncid,mse_th_id,mse_th))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))
end if

! Release memory
deallocate(fld_ref)
deallocate(fld_save)

end subroutine nicas_test_randomization

!----------------------------------------------------------------------
! Subroutine: nicas_test_consistency
! Purpose: test HDIAG-NICAS consistency with a randomization method
!----------------------------------------------------------------------
subroutine nicas_test_consistency(nicas,mpl,rng,nam,geom,bpar,io)

implicit none

! Passed variables
class(nicas_type),intent(inout) :: nicas ! NICAS data
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(rng_type),intent(inout) :: rng      ! Random number generator
type(nam_type),intent(inout) :: nam      ! Namelist variables
type(geom_type),intent(inout) :: geom    ! Geometry
type(bpar_type),intent(in) :: bpar       ! Block parameters
type(io_type),intent(in) :: io           ! I/O

! Local variables
integer,parameter :: nrad = 5
integer :: irad,ib,il0
real(kind_real) :: rh,rv,rad(nrad),rh_diag(nrad),rv_diag(nrad),rh_norm,rv_norm
logical :: write_nicas
character(len=1024) :: prefix
type(cmat_type) :: cmat
type(ens_type) :: ens
type(hdiag_type) :: hdiag
type(nicas_type) :: nicas_test

! Initialization
cmat%allocated = .false.

! Save namelist parameters
prefix = nam%prefix
rh = nam%rh
rv = nam%rv
write_nicas = nam%write_nicas

do irad=1,nrad
   ! Set radius factor
   rad(irad) = real(irad,kind_real)/real(nrad,kind_real)

   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a,f4.2)') '--- Radii factor: ',rad(irad)
   call mpl%flush

   ! Copy namelist support radii into C matrix
   nam%rh = rh*rad(irad)
   nam%rv = rv*rad(irad)
   call cmat%from_nam(mpl,nam,geom,bpar)

   ! Setup C matrix sampling
   call cmat%setup_sampling(mpl,nam,geom,bpar)

   ! Run NICAS driver
   nam%write_nicas = .false.
   call nicas_test%run_nicas(mpl,rng,nam,geom,bpar,cmat)
   if (nam%default_seed) call rng%reseed(mpl)

   ! Randomize ensemble
   call nicas_test%randomize(mpl,rng,nam,geom,bpar,nam%ens1_ne,ens)
   if (nam%default_seed) call rng%reseed(mpl)

   ! Run HDIAG driver
   call hdiag%run_hdiag(mpl,rng,nam,geom,bpar,io,ens)
   if (nam%default_seed) call rng%reseed(mpl)

   ! Save result
   rh_diag(irad) = 0.0
   rv_diag(irad) = 0.0
   rh_norm = 0.0
   rv_norm = 0.0
   do ib=1,bpar%nbe
      if (bpar%nicas_block(ib)) then
         do il0=1,geom%nl0
            if (hdiag%cor_1%blk(0,ib)%fit_rh(il0)>0.0) then
               rh_diag(irad) = rh_diag(irad)+hdiag%cor_1%blk(0,ib)%fit_rh(il0)
               rh_norm = rh_norm+1.0
            end if
            if (hdiag%cor_1%blk(0,ib)%fit_rv(il0)>0.0) then
               rv_diag(irad) = rv_diag(irad)+hdiag%cor_1%blk(0,ib)%fit_rv(il0)
               rv_norm = rv_norm+1.0
            end if
         end do
      end if
   end do
   if (rh_norm>0.0) rh_diag(irad) = rh_diag(irad)/rh_norm
   if (rv_norm>0.0) rv_diag(irad) = rv_diag(irad)/rv_norm

   ! Write
   if (nam%write_hdiag) then
      write(nam%prefix,'(a,a,i2.2)') trim(prefix),'_',int(rad(irad)*10.0)
      call hdiag%cor_1%write(mpl,nam,geom,bpar,io,hdiag%samp)
      nam%prefix = prefix
   end if

   ! Release memory
   call cmat%dealloc
   call nicas_test%dealloc
   call ens%dealloc
   call hdiag%dealloc
end do

! Print factors
do irad=1,nrad
   write(mpl%info,'(a7,a,f4.2,a)') '','Radii factor: ',rad(irad),':'
   call mpl%flush
   if (rh_diag(irad)>0.0) then
      write(mpl%info,'(a10,a,f10.2,a,f10.2,a,f5.3,a)') '','Diagnostic for rh: ',rh*rad(irad)*reqkm,' ~> ',rh_diag(irad)*reqkm, &
    & ' (',rh*rad(irad)/rh_diag(irad),')'
      call mpl%flush
   end if
   if (rv_diag(irad)>0.0) then
      write(mpl%info,'(a10,a,f10.2,a,f10.2,a,f5.3,a)') '','Diagnostic for rv: ',rv*rad(irad),' ~> ',rv_diag(irad), &
    & ' (',rv*rad(irad)/rv_diag(irad),')'
      call mpl%flush
   end if
end do

! Reset namelist parameters
nam%prefix = prefix
nam%rh = rh
nam%rv = rv
nam%write_nicas = write_nicas

end subroutine nicas_test_consistency

!----------------------------------------------------------------------
! Subroutine: nicas_test_optimality
! Purpose: test HDIAG localization optimality with a randomization method
!----------------------------------------------------------------------
subroutine nicas_test_optimality(nicas,mpl,rng,nam,geom,bpar,io)

implicit none

! Passed variables
class(nicas_type),intent(in) :: nicas ! NICAS data
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(rng_type),intent(inout) :: rng   ! Random number generator
type(nam_type),intent(inout) :: nam   ! Namelist variables
type(geom_type),intent(in) :: geom    ! Geometry
type(bpar_type),intent(in) :: bpar    ! Block parameters
type(io_type),intent(in) :: io        ! I/O

! Local variables
integer :: ib,ifac,itest,il0
real(kind_real) :: fac(-nfac_opt:nfac_opt),mse(ntest,-nfac_opt:nfac_opt),mse_sum,mse_max,rh_sum,rh_tot,rv_sum,rv_tot
real(kind_real) :: fld_ref(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest),fld_save(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest)
real(kind_real) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts)
character(len=1024) :: method
type(cmat_type) :: cmat
type(diag_type) :: loc_opt
type(hdiag_type) :: hdiag
type(ens_type) :: ens_loc,ens_test
type(nicas_type) :: nicas_test

! Define test vectors
write(mpl%info,'(a4,a)') '','Define test vectors'
call mpl%flush
call define_test_vectors(mpl,rng,nam,geom,ntest,fld_save)
if (nam%default_seed) call rng%reseed(mpl)

! Apply NICAS to test vectors
write(mpl%info,'(a4,a)') '','Apply NICAS to test vectors'
call mpl%flush
fld_ref = fld_save
do itest=1,ntest
   call nicas%apply_from_sqrt(mpl,nam,geom,bpar,fld_ref(:,:,:,:,itest))
end do

! Randomize ensemble to compute localization
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Randomize ensemble to compute localization'
call mpl%flush
call nicas%randomize(mpl,rng,nam,geom,bpar,nam%ens1_ne,ens_loc)
if (nam%default_seed) call rng%reseed(mpl)

! Save namelist variables
method = nam%method

! Set namelist variables
nam%method = 'loc'

! Call HDIAG driver
call hdiag%run_hdiag(mpl,rng,nam,geom,bpar,io,ens_loc)
if (nam%default_seed) call rng%reseed(mpl)

! Randomize ensemble to test localization
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Randomize ensemble to test localization'
call mpl%flush
call nicas%randomize(mpl,rng,nam,geom,bpar,nam%ne,ens_test)
if (nam%default_seed) call rng%reseed(mpl)

! Allocation
call loc_opt%alloc(mpl,nam,geom,bpar,hdiag%samp,'loc_opt')

! Initialization
mse_max = huge(1.0)

do ifac=-nfac_opt,nfac_opt
   ! Copy HDIAG into C matrix
   call cmat%from_hdiag(mpl,nam,geom,bpar,hdiag)

   ! Setup C matrix sampling
   call cmat%setup_sampling(mpl,nam,geom,bpar)

   ! Multiplication factor
   fac(ifac) = 1.0+real(ifac,kind_real)/real(nfac_opt+1,kind_real)

   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a,f4.2,a)') '--- Generate NICAS with a multiplicative factor ',fac(ifac),' to length-scales'
   call mpl%flush

   ! Allocation
   call nicas_test%alloc(mpl,nam,bpar,'nicas_test')

   do ib=1,bpar%nbe
      if (bpar%nicas_block(ib)) then
         ! Length-scales multiplication
         cmat%blk(ib)%rhs = fac(ifac)*cmat%blk(ib)%rhs
         cmat%blk(ib)%rvs = fac(ifac)*cmat%blk(ib)%rvs
         cmat%blk(ib)%rh = fac(ifac)*cmat%blk(ib)%rh
         cmat%blk(ib)%rv = fac(ifac)*cmat%blk(ib)%rv

         ! Compute NICAS parameters
         call nicas_test%blk(ib)%compute_parameters(mpl,rng,nam,geom,cmat%blk(ib))
      end if

      if (bpar%B_block(ib)) then
         ! Copy weights
         nicas_test%blk(ib)%wgt = cmat%blk(ib)%wgt
         if (bpar%nicas_block(ib)) then
            allocate(nicas_test%blk(ib)%coef_ens(geom%nc0a,geom%nl0))
            nicas_test%blk(ib)%coef_ens = cmat%blk(ib)%coef_ens
         end if
      end if
   end do

   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a)') '--- Apply ensemble B to test vectors'
   call mpl%flush

   do itest=1,ntest
      ! Test NICAS
      fld = fld_save(:,:,:,:,itest)
      call nicas_test%apply_bens(mpl,nam,geom,bpar,ens_test,fld)

      ! RMSE
      mse_sum = sum((fld-fld_ref(:,:,:,:,itest))**2,mask=mpl%msv%isnot(fld_ref(:,:,:,:,itest)))
      call mpl%f_comm%allreduce(mse_sum,mse(itest,ifac),fckit_mpi_sum())
   end do

   ! Test score
   if (sum(mse(:,ifac))<mse_max) then
      mse_max = sum(mse(:,ifac))
      do ib=1,bpar%nbe
         if (bpar%nicas_block(ib)) then
            do il0=1,geom%nl0
               rh_sum = sum(cmat%blk(ib)%rh(:,il0),mask=mpl%msv%isnot(cmat%blk(ib)%rh(:,il0)))
               call mpl%f_comm%allreduce(rh_sum,rh_tot,fckit_mpi_sum())
               loc_opt%blk(0,ib)%fit_rh = rh_tot/real(geom%nc0_gmask(il0),kind_real)
               if ((nam%nl>1).and.(nam%rv>0.0)) then
                  rv_sum = sum(cmat%blk(ib)%rv(:,il0),mask=mpl%msv%isnot(cmat%blk(ib)%rv(:,il0)))
                  call mpl%f_comm%allreduce(rv_sum,rv_tot,fckit_mpi_sum())
                  loc_opt%blk(0,ib)%fit_rv = rv_tot/real(geom%nc0_gmask(il0),kind_real)
               end if
            end do
         end if
      end do
   end if

   ! Print scores
   write(mpl%info,'(a)') '-------------------------------------------------------------------'
   call mpl%flush
   write(mpl%info,'(a,f4.2,a,e15.8)') '--- Optimality results for a factor ',fac(ifac),', MSE: ', &
 & sum(mse(:,ifac))/real(ntest,kind_real)
   call mpl%flush

   ! Release memory
   call cmat%dealloc
   call nicas_test%dealloc
end do

! Print scores summary
write(mpl%info,'(a)') '-------------------------------------------------------------------'
call mpl%flush
write(mpl%info,'(a)') '--- Optimality results summary'
call mpl%flush
do ifac=-nfac_opt,nfac_opt
   write(mpl%info,'(a7,a,f4.2,a,e15.8)') '','Factor ',fac(ifac),', MSE: ',sum(mse(:,ifac))/real(ntest,kind_real)
   call mpl%flush
end do

! Best fit
do ib=1,bpar%nbe
   if (bpar%diag_block(ib)) then
      call fit_diag(mpl,nam%nc3,bpar%nl0r(ib),geom%nl0,bpar%l0rl0b_to_l0(:,:,ib),geom%disth,loc_opt%blk(0,ib)%distv, &
    & loc_opt%blk(0,ib)%coef_ens,loc_opt%blk(0,ib)%fit_rh,loc_opt%blk(0,ib)%fit_rv,loc_opt%blk(0,ib)%fit)
      call loc_opt%blk(0,ib)%write(mpl,nam,geom,bpar,trim(nam%prefix)//'_diag')
   end if
end do

! Reset namelist variables
nam%method = method

! Release memory
call ens_loc%dealloc
call ens_test%dealloc
call hdiag%dealloc

end subroutine nicas_test_optimality

!----------------------------------------------------------------------
! Subroutine: define_test_vectors
! Purpose: define test vectors
!----------------------------------------------------------------------
subroutine define_test_vectors(mpl,rng,nam,geom,ntest,fld)

! Passed variables
type(mpl_type),intent(inout) :: mpl                                         ! MPI data
type(rng_type),intent(inout) :: rng                                         ! Random number generator
type(nam_type),intent(in) :: nam                                            ! Namelist
type(geom_type),intent(in) :: geom                                          ! Geometry
integer,intent(in) :: ntest                                                 ! Number of vectors
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0,nam%nv,nam%nts,ntest) ! Field

! Local variables
integer :: itest
integer :: il0dir,iprocdir,ic0adir

do itest=1,ntest
   ! Define random dirac location
   call geom%rand_level(mpl,rng,il0dir)
   call geom%rand_point(mpl,rng,il0dir,iprocdir,ic0adir)

   ! Define test vector
   fld(:,:,:,:,itest) = 0.0
   if (iprocdir==mpl%myproc) fld(ic0adir,il0dir,1,1,itest) = 1.0
end do

end subroutine define_test_vectors

end module type_nicas
