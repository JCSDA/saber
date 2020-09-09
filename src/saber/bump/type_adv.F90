!----------------------------------------------------------------------
! Module: type_adv
! Purpose: advection derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_adv

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_status
use netcdf
!$ use omp_lib
use tools_const, only: pi,req,reqkm,rad2deg,deg2rad
use tools_func, only: lonlatmod,lonlathash,sphere_dist,reduce_arc,lonlat2xyz,xyz2lonlat,vector_product
use tools_kinds, only: kind_real,nc_kind_real,huge_real
use tools_qsort, only: qsort
use tools_repro, only: inf,sup,eq
use type_bpar, only: bpar_type
use type_com, only: com_type
use type_ens, only: ens_type
use type_geom, only: geom_type
use type_io, only: io_type
use type_linop, only: linop_type
use type_mesh, only: mesh_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_samp, only: samp_type
use type_rng, only: rng_type

implicit none

! Advection derived type
type adv_type
   integer :: niter                                   ! Number of stored iterations
   real(kind_real),allocatable :: lon_c2a(:)          ! Origin longitude
   real(kind_real),allocatable :: lat_c2a(:)          ! Origin latitude
   real(kind_real),allocatable :: x_c2a(:)            ! Origin x coordinate
   real(kind_real),allocatable :: y_c2a(:)            ! Origin y coordinate
   real(kind_real),allocatable :: z_c2a(:)            ! Origin z coordinate
   real(kind_real),allocatable :: lon_c2a_raw(:,:,:)  ! Raw advected longitude
   real(kind_real),allocatable :: lat_c2a_raw(:,:,:)  ! Raw advected latitude
   real(kind_real),allocatable :: dist_c2a_raw(:,:,:) ! Raw advection distance
   real(kind_real),allocatable :: dist_raw(:,:)       ! Averaged raw advection distance
   real(kind_real),allocatable :: valid_raw(:,:)      ! Averaged raw advection validity
   real(kind_real),allocatable :: lon_c2a_flt(:,:,:)  ! Filtered advected longitude
   real(kind_real),allocatable :: lat_c2a_flt(:,:,:)  ! Filtered advected latitude
   real(kind_real),allocatable :: dist_c2a_flt(:,:,:) ! Filtered advection distance
   real(kind_real),allocatable :: dist_flt(:,:)       ! Averaged filtered advection distance
   real(kind_real),allocatable :: valid_flt(:,:)      ! Averaged filtered advection validity
   real(kind_real),allocatable :: rhflt(:,:)          ! Optimal filtering support radius
   real(kind_real),allocatable :: cor_loc(:,:,:,:)    ! Local correlation
   real(kind_real),allocatable :: cor_adv(:,:,:,:)    ! Advected correlation
   real(kind_real),allocatable :: score_loc(:)        ! Local correlation score
   real(kind_real),allocatable :: score_adv(:)        ! Advected correlation score
contains
   procedure :: alloc => adv_alloc
   procedure :: init => adv_init
   procedure :: dealloc => adv_dealloc
   procedure :: compute => adv_compute
   procedure :: compute_max => adv_compute_max
   procedure :: compute_wind => adv_compute_wind
   procedure :: filter => adv_filter
   procedure :: interp => adv_interp
   procedure :: test => adv_test
   procedure :: write => adv_write
end type adv_type

integer,parameter :: nt = 6 ! Number of substeps for wind advection

private
public :: adv_type

contains

!----------------------------------------------------------------------
! Subroutine: adv_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine adv_alloc(adv,nam,geom,samp)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv ! Advection
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(samp_type),intent(in) :: samp   ! Sampling

! Allocation
allocate(adv%lon_c2a(samp%nc2a))
allocate(adv%lat_c2a(samp%nc2a))
allocate(adv%x_c2a(samp%nc2a))
allocate(adv%y_c2a(samp%nc2a))
allocate(adv%z_c2a(samp%nc2a))
allocate(adv%lon_c2a_raw(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%lat_c2a_raw(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%dist_c2a_raw(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%dist_raw(geom%nl0,2:nam%nts))
allocate(adv%valid_raw(geom%nl0,2:nam%nts))
allocate(adv%lon_c2a_flt(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%lat_c2a_flt(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%dist_c2a_flt(samp%nc2a,geom%nl0,2:nam%nts))
allocate(adv%dist_flt(geom%nl0,2:nam%nts))
allocate(adv%valid_flt(geom%nl0,2:nam%nts))
allocate(adv%rhflt(geom%nl0,2:nam%nts))
allocate(adv%cor_loc(geom%nc0a,geom%nl0,nam%nv,2:nam%nts))
allocate(adv%cor_adv(geom%nc0a,geom%nl0,nam%nv,2:nam%nts))
allocate(adv%score_loc(2:nam%nts))
allocate(adv%score_adv(2:nam%nts))

end subroutine adv_alloc

!----------------------------------------------------------------------
! Subroutine: adv_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine adv_dealloc(adv)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv ! Advection

! Release memory
if (allocated(adv%lon_c2a)) deallocate(adv%lon_c2a)
if (allocated(adv%lat_c2a)) deallocate(adv%lat_c2a)
if (allocated(adv%x_c2a)) deallocate(adv%x_c2a)
if (allocated(adv%y_c2a)) deallocate(adv%y_c2a)
if (allocated(adv%z_c2a)) deallocate(adv%z_c2a)
if (allocated(adv%lon_c2a_raw)) deallocate(adv%lon_c2a_raw)
if (allocated(adv%lat_c2a_raw)) deallocate(adv%lat_c2a_raw)
if (allocated(adv%dist_c2a_raw)) deallocate(adv%dist_c2a_raw)
if (allocated(adv%dist_raw)) deallocate(adv%dist_raw)
if (allocated(adv%valid_raw)) deallocate(adv%valid_raw)
if (allocated(adv%lon_c2a_flt)) deallocate(adv%lon_c2a_flt)
if (allocated(adv%lat_c2a_flt)) deallocate(adv%lat_c2a_flt)
if (allocated(adv%dist_c2a_flt)) deallocate(adv%dist_c2a_flt)
if (allocated(adv%dist_flt)) deallocate(adv%dist_flt)
if (allocated(adv%valid_flt)) deallocate(adv%valid_flt)
if (allocated(adv%rhflt)) deallocate(adv%rhflt)
if (allocated(adv%cor_loc)) deallocate(adv%cor_loc)
if (allocated(adv%cor_adv)) deallocate(adv%cor_adv)
if (allocated(adv%score_loc)) deallocate(adv%score_loc)
if (allocated(adv%score_adv)) deallocate(adv%score_adv)

end subroutine adv_dealloc

!----------------------------------------------------------------------
! Subroutine: adv_compute
! Purpose: compute advection
!----------------------------------------------------------------------
subroutine adv_compute(adv,mpl,rng,nam,geom,bpar,samp,io,ens,fld_uv)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv                                        ! Advection
type(mpl_type),intent(inout) :: mpl                                         ! MPI data
type(rng_type),intent(inout) :: rng                                         ! Random number generator
type(nam_type),intent(in) :: nam                                            ! Namelist
type(geom_type),intent(in) :: geom                                          ! Geometry
type(bpar_type),intent(in) :: bpar                                          ! Block parameters
type(samp_type),intent(inout) :: samp                                       ! Sampling
type(io_type),intent(in) :: io                                              ! I/O
type(ens_type), intent(in) :: ens                                           ! Ensemble
real(kind_real),intent(in),optional :: fld_uv(geom%nc0a,geom%nl0,2,nam%nts) ! Wind field

! Local variables
type(adv_type) :: adv_wind
character(len=1024),parameter :: subr = 'adv_compute'

! Allocation
call adv%alloc(nam,geom,samp)

! Allocation
call adv%init(mpl,samp)

select case (trim(nam%adv_type))
case ('max')
   ! Compute raw advection from maximum displacement, no guess
   call adv%compute_max(mpl,nam,geom,samp,ens)
case ('wind')
   ! Check wind field
   if (.not.present(fld_uv)) call mpl%abort(subr,'wind field absent')

   ! Compute raw advection from wind
   call adv%compute_wind(mpl,rng,nam,geom,samp,fld_uv)
case ('windmax')
   ! Check wind field
   if (.not.present(fld_uv)) call mpl%abort(subr,'wind field absent')

   ! Allocation
   call adv_wind%alloc(nam,geom,samp)

   ! Initialization
   call adv_wind%init(mpl,samp)

   ! Compute raw advection from wind
   call adv_wind%compute_wind(mpl,rng,nam,geom,samp,fld_uv)

   ! Compute raw advection from maximum displacement, using filtered wind advection as a guess
   call adv%compute_max(mpl,nam,geom,samp,ens,adv_wind)

   ! Release memory
   call adv_wind%dealloc
case default
   call mpl%abort(subr,'wrong advection diagnostic type')
end select

! Interpolate advection
call adv%interp(mpl,nam,geom,samp)

! Test advection efficiency
call adv%test(mpl,rng,nam,geom,samp,ens)

! Write
if (nam%write_hdiag) call adv%write(mpl,nam,geom,bpar,io,samp)

! Release memory
call adv%dealloc

end subroutine adv_compute

!----------------------------------------------------------------------
! Subroutine: adv_init
! Purpose: initialize advection computation
!----------------------------------------------------------------------
subroutine adv_init(adv,mpl,samp)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv ! Advection
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(samp_type),intent(in) :: samp   ! Sampling

! Local variables
integer :: ic2a

! Initialization
adv%lon_c2a = samp%lon_c2a
adv%lat_c2a = samp%lat_c2a
do ic2a=1,samp%nc2a
   call lonlat2xyz(mpl,adv%lon_c2a(ic2a),adv%lat_c2a(ic2a),adv%x_c2a(ic2a),adv%y_c2a(ic2a),adv%z_c2a(ic2a))
end do

end subroutine adv_init

!----------------------------------------------------------------------
! Subroutine: adv_compute_max
! Purpose: compute raw advection from correlation maximum
!----------------------------------------------------------------------
subroutine adv_compute_max(adv,mpl,nam,geom,samp,ens,adv_wind)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv           ! Advection
type(mpl_type),intent(inout) :: mpl            ! MPI data
type(nam_type),intent(in) :: nam               ! Namelist
type(geom_type),intent(in) :: geom             ! Geometry
type(samp_type),intent(inout) :: samp          ! Sampling
type(ens_type),intent(in) :: ens               ! Ensemble
type(adv_type),intent(in),optional :: adv_wind ! Wind advection

! Local variables
integer :: ic0u,ic0,ic2,ic2a,ic2u,jn,jc0u,il0,il0i,isub,iv,its,ie,ie_sub,ic0a,nc0d,ic0d,jc0d,jnmax,nnmax,iproc
integer :: nvalid,nmask,nvalid_tot,nmask_tot
integer :: nn(samp%nc2a,geom%nl0),ic0u_rac(1)
integer,allocatable :: jc0u_ra(:,:,:),c0u_to_c0d(:),c0d_to_c0(:)
real(kind_real) :: search_rad,m2m2,fld_1,fld_2,cov
real(kind_real) :: dist_sum,cor_avg_max,norm,norm_tot
real(kind_real) :: lon_rac(samp%nc2a,geom%nl0),lat_rac(samp%nc2a,geom%nl0)
real(kind_real) :: lon_c2u(samp%nc2u),lat_c2u(samp%nc2u)
real(kind_real),allocatable :: fld_ext_1(:,:),fld_ext_2(:,:)
real(kind_real),allocatable :: m2_1(:,:,:,:,:),m2_2(:,:,:,:,:)
real(kind_real),allocatable :: m11(:,:,:,:,:)
real(kind_real),allocatable :: cor(:,:),cor_avg(:)
logical :: lcheck_c0d(geom%nc0u),valid_c2u(samp%nc2u),valid_c2a(samp%nc2a)
logical,allocatable :: mask_nn(:,:,:),mask_nn_tmp(:)
type(com_type) :: com_AD
type(mesh_type) :: mesh

write(mpl%info,'(a7,a)') '','Compute raw advection using correlation maximum'
call mpl%flush

do its=2,nam%nts
   write(mpl%info,'(a10,a,i2)') '','Timeslot ',its
   call mpl%flush

   ! Find research area center
   write(mpl%info,'(a13,a)') '','Find research area center'
   call mpl%flush
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         if (samp%smask_c2a(ic2a,il0)) then
            ! Research area center
            if (present(adv_wind)) then
               ! Use wind-based estimation
               lon_rac(ic2a,il0) = adv_wind%lon_c2a_flt(ic2a,il0,its)
               lat_rac(ic2a,il0) = adv_wind%lat_c2a_flt(ic2a,il0,its)
            else
               if (its==2) then
                  ! Set at origin point
                  lon_rac(ic2a,il0) = adv%lon_c2a(ic2a)
                  lat_rac(ic2a,il0) = adv%lat_c2a(ic2a)
               else
                  ! Use previous timeslot filtered estimation
                  lon_rac(ic2a,il0) = adv%lon_c2a_flt(ic2a,il0,its-1)
                  lat_rac(ic2a,il0) = adv%lat_c2a_flt(ic2a,il0,its-1)
               end if
            end if
         end if
      end do
   end do

   ! Count nearest neighbors
   write(mpl%info,'(a13,a)') '','Count nearest neighbors:'
   call mpl%flush(.false.)
   call mpl%prog_init(samp%nc2a*geom%nl0)
   nn = mpl%msv%vali
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         if (samp%smask_c2a(ic2a,il0)) then
            ! Index
            call geom%tree_c0u%find_nearest_neighbors(lon_rac(ic2a,il0),lat_rac(ic2a,il0),1,ic0u_rac)

            ! Define search radius
            il0i = min(il0,geom%nl0i)
            search_rad = min(nam%adv_rad,min(geom%mdist_c0u(ic0u_rac(1),il0i),geom%mesh_c0u%bdist(ic0u_rac(1))))

            ! Count nearest neighbors
            call geom%tree_c0u%count_nearest_neighbors(lon_rac(ic2a,il0),lat_rac(ic2a,il0),search_rad,nn(ic2a,il0))

            ! Update
            call mpl%prog_print((il0-1)*samp%nc2a+ic2a)
         end if
      end do
   end do
   call mpl%prog_final

   ! Get maximum
   nnmax = maxval(nn)

   ! Allocation
   allocate(jc0u_ra(nnmax,samp%nc2a,geom%nl0))
   allocate(mask_nn(nnmax,samp%nc2a,geom%nl0))
   allocate(m2_1(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))
   allocate(m2_2(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))
   allocate(m11(nnmax,samp%nc2a,geom%nl0,nam%nv,ens%nsub))

   ! Find nearest neighbors
   write(mpl%info,'(a13,a)') '','Find nearest neighbors:'
   call mpl%flush(.false.)
   call mpl%prog_init(samp%nc2a*geom%nl0)
   lcheck_c0d = .false.
   do ic0a=1,geom%nc0a
      ic0u = geom%c0a_to_c0u(ic0a)
      lcheck_c0d(ic0u) = .true.
   end do
   jc0u_ra = mpl%msv%vali
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         if (samp%smask_c2a(ic2a,il0)) then
            ! Find nearest neighbors
            call geom%tree_c0u%find_nearest_neighbors(lon_rac(ic2a,il0),lat_rac(ic2a,il0),nn(ic2a,il0), &
 & jc0u_ra(1:nn(ic2a,il0),ic2a,il0))

            ! Check points
            do jn=1,nn(ic2a,il0)
               jc0u = jc0u_ra(jn,ic2a,il0)
               lcheck_c0d(jc0u) = .true.
            end do

            ! Update
            call mpl%prog_print((il0-1)*samp%nc2a+ic2a)
         end if
      end do
   end do
   call mpl%prog_final

   ! Halo D size
   nc0d = count(lcheck_c0d)

   ! Allocation
   allocate(c0u_to_c0d(geom%nc0u))
   allocate(c0d_to_c0(nc0d))
   allocate(fld_ext_1(nc0d,geom%nl0))
   allocate(fld_ext_2(nc0d,geom%nl0))

   ! Conversions
   c0u_to_c0d = mpl%msv%vali
   ic0d = 0
   do ic0u=1,geom%nc0u
      if (lcheck_c0d(ic0u)) then
         ic0d = ic0d+1
         ic0 = geom%c0u_to_c0(ic0u)
         c0u_to_c0d(ic0u) = ic0d
         c0d_to_c0(ic0d) = ic0
      end if
   end do

   ! Setup communications
   call com_AD%setup(mpl,'com_AD',geom%nc0a,nc0d,geom%nc0,geom%c0a_to_c0,c0d_to_c0)

   ! Mask
   mask_nn = .false.
   do il0=1,geom%nl0
      do ic2a=1,samp%nc2a
         if (samp%smask_c2a(ic2a,il0)) then
            do jn=1,nn(ic2a,il0)
               jc0u = jc0u_ra(jn,ic2a,il0)
               mask_nn(jn,ic2a,il0) = geom%gmask_c0u(jc0u,il0)
             end do
         end if
      end do
   end do

   ! Initialization
   m11 = 0.0
   m2_1 = 0.0
   m2_2 = 0.0

   ! Compute moments
   write(mpl%info,'(a13,a)') '','Compute moments'
   call mpl%flush
   do isub=1,ens%nsub
      if (ens%nsub==1) then
         write(mpl%info,'(a16,a)') '','Full ensemble, member:'
         call mpl%flush(.false.)
      else
         write(mpl%info,'(a16,a,i6,a)') '','Sub-ensemble ',isub,', member:'
         call mpl%flush(.false.)
      end if

      ! Compute covariances
      do ie_sub=1,ens%ne/ens%nsub
         write(mpl%info,'(i6)') ie_sub
         call mpl%flush(.false.)

         ! Full ensemble index
         ie = ie_sub+(isub-1)*ens%ne/ens%nsub

         do iv=1,nam%nv
            ! Halo extension
            call com_AD%ext(mpl,geom%nl0,ens%mem(ie)%fld(:,:,iv,1),fld_ext_1)
            call com_AD%ext(mpl,geom%nl0,ens%mem(ie)%fld(:,:,iv,its),fld_ext_2)

            !$omp parallel do schedule(static) private(il0,ic2a,jn,ic2,ic0a,ic0u,jc0u,ic0d,jc0d,fld_1,fld_2)
            do il0=1,geom%nl0
               do ic2a=1,samp%nc2a
                  if (samp%smask_c2a(ic2a,il0)) then
                     do jn=1,nn(ic2a,il0)
                        if (mask_nn(jn,ic2a,il0)) then
                           ! Indices
                           ic0a = samp%c2a_to_c0a(ic2a)
                           ic0u = geom%c0a_to_c0u(ic0a)
                           jc0u = jc0u_ra(jn,ic2a,il0)

                           ! Halo D indices
                           ic0d = c0u_to_c0d(ic0u)
                           jc0d = c0u_to_c0d(jc0u)

                           ! Copy points
                           fld_1 = fld_ext_1(ic0d,il0)
                           fld_2 = fld_ext_2(jc0d,il0)

                           ! Covariance
                           m11(jn,ic2a,il0,iv,isub) = m11(jn,ic2a,il0,iv,isub)+fld_1*fld_2

                           ! Variances
                           m2_1(jn,ic2a,il0,iv,isub) = m2_1(jn,ic2a,il0,iv,isub)+fld_1**2
                           m2_2(jn,ic2a,il0,iv,isub) = m2_2(jn,ic2a,il0,iv,isub)+fld_2**2
                        end if
                     end do
                  end if
               end do
            end do
            !$omp end parallel do
         end do
      end do
      write(mpl%info,'(a)') ''
      call mpl%flush
   end do

   ! Find correlation propagation
   write(mpl%info,'(a13,a)') '','Find correlation propagation'
   call mpl%flush

   do il0=1,geom%nl0
      write(mpl%info,'(a16,a,i3)') '','Level ',nam%levs(il0)
      call mpl%flush

      ! Initialization
      adv%lon_c2a_raw(:,il0,its) = adv%lon_c2a
      adv%lat_c2a_raw(:,il0,its) = adv%lat_c2a
      adv%dist_c2a_raw(:,il0,its) = mpl%msv%valr

      !$omp parallel do schedule(static) private(ic2a,jn,iv,cov,m2m2,jnmax,cor_avg_max,jc0u,ic2,ic0u,ic0a), &
      !$omp&                             firstprivate(cor,cor_avg,mask_nn_tmp)
      do ic2a=1,samp%nc2a
         if (samp%smask_c2a(ic2a,il0)) then
            if (nn(ic2a,il0)>0) then
               ! Allocation
               allocate(cor(nn(ic2a,il0),nam%nv))
               allocate(cor_avg(nn(ic2a,il0)))
               allocate(mask_nn_tmp(nn(ic2a,il0)))

               ! Initialization
               mask_nn_tmp = .false.
               jnmax = mpl%msv%vali
               cor_avg_max = 0.0

               do jn=1,nn(ic2a,il0)
                  if (mask_nn(jn,ic2a,il0)) then
                     ! Compute covariance and correlation
                     do iv=1,nam%nv
                        ! Covariance
                        cov = sum(m11(jn,ic2a,il0,iv,:))/real(ens%nsub,kind_real)

                        ! Correlation
                        m2m2 = sum(m2_1(jn,ic2a,il0,iv,:))*sum(m2_2(jn,ic2a,il0,iv,:))/real(ens%nsub**2,kind_real)
                        if (m2m2>0.0) then
                           cor(jn,iv) = cov/sqrt(m2m2)
                        else
                           cor(jn,iv) = mpl%msv%valr
                        end if
                     end do

                     ! Average correlation over variables
                     if (mpl%msv%isanynot(cor(jn,:))) then
                        cor_avg(jn) = sum(cor(jn,:),mask=mpl%msv%isnot(cor(jn,:)))/real(count(mpl%msv%isnot(cor(jn,:))),kind_real)
                     else
                        cor_avg(jn) = mpl%msv%valr
                     end if

                     ! Set mask
                     mask_nn_tmp(jn) = (cor_avg(jn)>0.0).and.mpl%msv%isnot(cor_avg(jn))
                  end if

                  ! Locate the maximum correlation
                  if (mask_nn_tmp(jn)) then
                     if (sup(cor_avg(jn),cor_avg_max)) then
                        jnmax = jn
                        cor_avg_max = cor_avg(jn)
                     end if
                  end if
               end do

               if (mpl%msv%isnot(jnmax)) then
                  ! Indices
                  jc0u = jc0u_ra(jnmax,ic2a,il0)

                  ! Save advection
                  adv%lon_c2a_raw(ic2a,il0,its) = geom%lon_c0u(jc0u)
                  adv%lat_c2a_raw(ic2a,il0,its) = geom%lat_c0u(jc0u)

                  ! Compute distance
                  call sphere_dist(adv%lon_c2a(ic2a),adv%lat_c2a(ic2a),adv%lon_c2a_raw(ic2a,il0,its), &
 & adv%lat_c2a_raw(ic2a,il0,its),adv%dist_c2a_raw(ic2a,il0,its))
               end if

               ! Release memory
               deallocate(cor)
               deallocate(cor_avg)
               deallocate(mask_nn_tmp)
            else
               ! Save advection
               adv%lon_c2a_raw(ic2a,il0,its) = lon_rac(ic2a,il0)
               adv%lat_c2a_raw(ic2a,il0,its) = lat_rac(ic2a,il0)

               ! Compute distance
               call sphere_dist(adv%lon_c2a(ic2a),adv%lat_c2a(ic2a),adv%lon_c2a_raw(ic2a,il0,its), &
 & adv%lat_c2a_raw(ic2a,il0,its),adv%dist_c2a_raw(ic2a,il0,its))
            end if
         end if
      end do
      !$omp end parallel do

      ! Check raw mesh
      call samp%com_AU%ext(mpl,adv%lon_c2a_raw(:,il0,its),lon_c2u)
      call samp%com_AU%ext(mpl,adv%lat_c2a_raw(:,il0,its),lat_c2u)
      call mesh%copy(samp%mesh)
      call mesh%store(mpl,lon_c2u,lat_c2u)
      call mesh%check(mpl,valid_c2u)
      do ic2u=1,samp%nc2u
         ic0u = samp%c2u_to_c0u(ic2u)
         ic0 = geom%c0u_to_c0(ic0u)
         iproc = geom%c0_to_proc(ic0)
         if (iproc==mpl%myproc) then
            valid_c2u(ic2u) = valid_c2u(ic2u).and.samp%mesh%valid(ic2u).and.samp%smask_c2u(ic2u,il0)
         else
            valid_c2u(ic2u) = .true.
         end if
      end do
      call samp%com_AU%red(mpl,valid_c2u,valid_c2a,op='and')
      nvalid = count(valid_c2a)
      nmask = count(samp%mesh%valid(samp%c2a_to_c2u).and.samp%smask_c2a(:,il0))
      call mpl%f_comm%allreduce(nvalid,nvalid_tot,fckit_mpi_sum())
      call mpl%f_comm%allreduce(nmask,nmask_tot,fckit_mpi_sum())
      adv%valid_raw(il0,its) = real(nvalid_tot,kind_real)/real(nmask_tot,kind_real)

      ! Average distance
      dist_sum = sum(adv%dist_c2a_raw(:,il0,its),mask=mpl%msv%isnot(adv%dist_c2a_raw(:,il0,its)))
      norm = real(count(mpl%msv%isnot(adv%dist_c2a_raw(:,il0,its))),kind_real)
      call mpl%f_comm%allreduce(dist_sum,adv%dist_raw(il0,its),fckit_mpi_sum())
      call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
      adv%dist_raw(il0,its) = adv%dist_raw(il0,its)/norm_tot

      ! Filter advection
      write(mpl%info,'(a16,a)') '','Filter advection'
      call mpl%flush
      call adv%filter(mpl,nam,geom,samp,il0,its)
   end do

   ! Release memory
   deallocate(jc0u_ra)
   deallocate(mask_nn)
   deallocate(m2_1)
   deallocate(m2_2)
   deallocate(m11)
   deallocate(c0d_to_c0)
   call com_AD%dealloc
   deallocate(fld_ext_1)
   deallocate(fld_ext_2)
end do

end subroutine adv_compute_max

!----------------------------------------------------------------------
! Subroutine: adv_compute_wind
! Purpose: compute raw advection from wind field
!----------------------------------------------------------------------
subroutine adv_compute_wind(adv,mpl,rng,nam,geom,samp,fld_uv)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv                               ! Advection
type(mpl_type),intent(inout) :: mpl                                ! MPI data
type(rng_type),intent(inout) :: rng                                ! Random number generator
type(nam_type),intent(in) :: nam                                   ! Namelist
type(geom_type),intent(in) :: geom                                 ! Geometry
type(samp_type),intent(inout) :: samp                              ! Sampling
real(kind_real),intent(in) :: fld_uv(geom%nc0a,geom%nl0,2,nam%nts) ! Wind field

! Local variables
integer :: ic0u,ic0a,ic0w,ic0own,ic0,i_s,ic2a,ic2u,il0,its,it,nc0w,nc0own,iproc
integer :: nvalid,nmask,nvalid_tot,nmask_tot
integer :: c0u_to_c0w(geom%nc0u)
integer,allocatable :: c0w_to_c0(:),c0own_to_c0(:)
real(kind_real) :: dist_sum,norm,norm_tot
real(kind_real) :: dtl,t,um(samp%nc2a),vm(samp%nc2a),up(samp%nc2a),vp(samp%nc2a)
real(kind_real) :: uxm,uym,uzm,uxp,uyp,uzp,ux,uy,uz,x,y,z
real(kind_real) :: lon_c2u(samp%nc2u),lat_c2u(samp%nc2u)
real(kind_real),allocatable :: fld_uv_ext(:,:,:)
logical :: lcheck_c0w(geom%nc0u),valid_c2u(samp%nc2u),valid_c2a(samp%nc2a)
type(com_type) :: com_AW
type(linop_type) :: h
type(mesh_type) :: mesh

write(mpl%info,'(a7,a)') '','Compute raw advection using wind field'
call mpl%flush

! Initialization
dtl = nam%dts/nt

do its=2,nam%nts
   write(mpl%info,'(a10,a,i2)') '','Timeslot ',its
   call mpl%flush

   do il0=1,geom%nl0
      write(mpl%info,'(a13,a,i3)') '','Level ',nam%levs(il0)
      call mpl%flush

      ! Initialization
      if (its==2) then
         adv%lon_c2a_raw(:,il0,its) = adv%lon_c2a
         adv%lat_c2a_raw(:,il0,its) = adv%lat_c2a
      else
         adv%lon_c2a_raw(:,il0,its) = adv%lon_c2a_flt(:,il0,its-1)
         adv%lat_c2a_raw(:,il0,its) = adv%lat_c2a_flt(:,il0,its-1)
      end if
      adv%dist_c2a_raw(:,il0,its) = mpl%msv%valr

      write(mpl%info,'(a16,a)') '','Propagate locations: '
      call mpl%flush(.false.)

      call mpl%prog_init(nt)
      do it=1,nt
         ! Define internal time
         t = real(it-1,kind_real)/real(nt,kind_real)

         ! Compute interpolation
         call h%interp(mpl,rng,nam,geom,il0,geom%nc0u,geom%lon_c0u,geom%lat_c0u,geom%gmask_c0u(:,il0), &
 & samp%nc2a,adv%lon_c2a_raw(:,il0,its),adv%lat_c2a_raw(:,il0,its),samp%smask_c2a(:,il0),0)

         ! Define halo W
         lcheck_c0w = .false.
         do i_s=1,h%n_s
            ic0u = h%col(i_s)
            lcheck_c0w(ic0u) = .true.
         end do
         nc0w = count(lcheck_c0w)

         ! Count owned points
         nc0own = 0
         do ic0a=1,geom%nc0a
            ic0u = geom%c0a_to_c0u(ic0a)
            if (lcheck_c0w(ic0u)) nc0own = nc0own+1
         end do

         ! Allocation
         allocate(c0w_to_c0(nc0w))
         allocate(c0own_to_c0(nc0own))
         allocate(fld_uv_ext(nc0w,2,2))

         ! Global-local conversion for halo W
         ic0w = 0
         ic0own = 0
         c0u_to_c0w = mpl%msv%vali
         do ic0u=1,geom%nc0u
            if (lcheck_c0w(ic0u)) then
               ic0w = ic0w+1
               ic0 = geom%c0u_to_c0(ic0u)
               c0w_to_c0(ic0w) = ic0
               c0u_to_c0w(ic0u) = ic0w
            end if
         end do

         ! Owned points conversion
         ic0own = 0
         do ic0a=1,geom%nc0a
            ic0u = geom%c0a_to_c0u(ic0a)
            if (lcheck_c0w(ic0u)) then
               ic0own = ic0own+1
               ic0 = geom%c0u_to_c0(ic0u)
               c0own_to_c0(ic0own) = ic0
            end if
         end do

         ! Local interpolation source
         h%n_src = nc0w
         do i_s=1,h%n_s
            h%col(i_s) = c0u_to_c0w(h%col(i_s))
         end do

         ! Setup communication
         call com_AW%setup(mpl,'com_AW',geom%nc0a,nc0w,geom%nc0,geom%c0a_to_c0,c0w_to_c0,c0own_to_c0)

         ! Communication
         call com_AW%ext(mpl,fld_uv(:,il0,1,its-1),fld_uv_ext(:,1,1))
         call com_AW%ext(mpl,fld_uv(:,il0,2,its-1),fld_uv_ext(:,2,1))
         call com_AW%ext(mpl,fld_uv(:,il0,1,its),fld_uv_ext(:,1,2))
         call com_AW%ext(mpl,fld_uv(:,il0,2,its),fld_uv_ext(:,2,2))

         ! Interpolate wind value at diagnostic points
         call h%apply(mpl,fld_uv_ext(:,1,1),um)
         call h%apply(mpl,fld_uv_ext(:,2,1),vm)
         call h%apply(mpl,fld_uv_ext(:,1,2),up)
         call h%apply(mpl,fld_uv_ext(:,2,2),vp)

         ! Release memory
         call h%dealloc
         deallocate(c0w_to_c0)
         deallocate(c0own_to_c0)
         call com_AW%dealloc
         deallocate(fld_uv_ext)

         do ic2a=1,samp%nc2a
            if (samp%smask_c2a(ic2a,il0)) then
               ! Transform wind to cartesian coordinates
               uxm = -sin(adv%lon_c2a_raw(ic2a,il0,its))*um(ic2a)- &
 & cos(adv%lon_c2a_raw(ic2a,il0,its))*sin(adv%lat_c2a_raw(ic2a,il0,its))*vm(ic2a)
               uym = cos(adv%lon_c2a_raw(ic2a,il0,its))*um(ic2a)- &
 & sin(adv%lon_c2a_raw(ic2a,il0,its))*sin(adv%lat_c2a_raw(ic2a,il0,its))*vm(ic2a)
               uzm = cos(adv%lat_c2a_raw(ic2a,il0,its))*vm(ic2a)
               uxp = -sin(adv%lon_c2a_raw(ic2a,il0,its))*up(ic2a)- &
 & cos(adv%lon_c2a_raw(ic2a,il0,its))*sin(adv%lat_c2a_raw(ic2a,il0,its))*vp(ic2a)
               uyp = cos(adv%lon_c2a_raw(ic2a,il0,its))*up(ic2a)- &
 & sin(adv%lon_c2a_raw(ic2a,il0,its))*sin(adv%lat_c2a_raw(ic2a,il0,its))*vp(ic2a)
               uzp = cos(adv%lat_c2a_raw(ic2a,il0,its))*vp(ic2a)

               ! Define wind in cartesian coordinates
               ux = (1.0-t)*uxm+t*uxp
               uy = (1.0-t)*uym+t*uyp
               uz = (1.0-t)*uzm+t*uzp

               ! Transform location to cartesian coordinates
               call lonlat2xyz(mpl,adv%lon_c2a_raw(ic2a,il0,its),adv%lat_c2a_raw(ic2a,il0,its),x,y,z)

               ! Propagate location
               x = x+ux*dtl
               y = y+uy*dtl
               z = z+uz*dtl

               ! Back to spherical coordinates
               call xyz2lonlat(mpl,x,y,z,adv%lon_c2a_raw(ic2a,il0,its),adv%lat_c2a_raw(ic2a,il0,its))
            end if
         end do

         ! Update
         call mpl%prog_print(it)
      end do
      call mpl%prog_final

      do ic2a=1,samp%nc2a
         if (samp%smask_c2a(ic2a,il0)) then
            ! Compute distance
            call sphere_dist(adv%lon_c2a(ic2a),adv%lat_c2a(ic2a),adv%lon_c2a_raw(ic2a,il0,its), &
 & adv%lat_c2a_raw(ic2a,il0,its),adv%dist_c2a_raw(ic2a,il0,its))
         end if
      end do

      ! Check raw mesh
      call samp%com_AU%ext(mpl,adv%lon_c2a_raw(:,il0,its),lon_c2u)
      call samp%com_AU%ext(mpl,adv%lat_c2a_raw(:,il0,its),lat_c2u)
      call mesh%copy(samp%mesh)
      call mesh%store(mpl,lon_c2u,lat_c2u)
      call mesh%check(mpl,valid_c2u)
      do ic2u=1,samp%nc2u
         ic0u = samp%c2u_to_c0u(ic2u)
         ic0 = geom%c0u_to_c0(ic0u)
         iproc = geom%c0_to_proc(ic0)
         if (iproc==mpl%myproc) then
            valid_c2u(ic2u) = valid_c2u(ic2u).and.samp%mesh%valid(ic2u).and.samp%smask_c2u(ic2u,il0)
         else
            valid_c2u(ic2u) = .true.
         end if
      end do
      call samp%com_AU%red(mpl,valid_c2u,valid_c2a,op='and')
      nvalid = count(valid_c2a)
      nmask = count(samp%mesh%valid(samp%c2a_to_c2u).and.samp%smask_c2a(:,il0))
      call mpl%f_comm%allreduce(nvalid,nvalid_tot,fckit_mpi_sum())
      call mpl%f_comm%allreduce(nmask,nmask_tot,fckit_mpi_sum())
      adv%valid_raw(il0,its) = real(nvalid_tot,kind_real)/real(nmask_tot,kind_real)

      ! Average distance
      dist_sum = sum(adv%dist_c2a_raw(:,il0,its),mask=mpl%msv%isnot(adv%dist_c2a_raw(:,il0,its)))
      norm = real(count(mpl%msv%isnot(adv%dist_c2a_raw(:,il0,its))),kind_real)
      call mpl%f_comm%allreduce(dist_sum,adv%dist_raw(il0,its),fckit_mpi_sum())
      call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
      adv%dist_raw(il0,its) = adv%dist_raw(il0,its)/norm_tot

      ! Filter advection
      write(mpl%info,'(a16,a)') '','Filter advection'
      call mpl%flush
      call adv%filter(mpl,nam,geom,samp,il0,its)
   end do
end do

end subroutine adv_compute_wind

!----------------------------------------------------------------------
! Subroutine: adv_filter
! Purpose: filter advection
!----------------------------------------------------------------------
subroutine adv_filter(adv,mpl,nam,geom,samp,il0,its)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv ! Advection
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
type(samp_type),intent(in) :: samp   ! Sampling
integer,intent(in) :: il0            ! Level index
integer,intent(in) :: its            ! Timeslot index

! Local variables
integer :: ic0,ic0a,ic0u,ic2a,ic2u,il0i,iter,iproc
integer :: nvalid,nmask,nvalid_tot,nmask_tot
real(kind_real) :: dist_sum,norm,norm_tot,valid_flt,dist_flt,rhflt,drhflt
real(kind_real) :: lon_c2a(samp%nc2a),lat_c2a(samp%nc2a),dist_c2a(samp%nc2a)
real(kind_real) :: dx_ini(samp%nc2a),dy_ini(samp%nc2a),dz_ini(samp%nc2a)
real(kind_real) :: dx(samp%nc2a),dy(samp%nc2a),dz(samp%nc2a),dd(samp%nc2a)
real(kind_real) :: lon_c2u(samp%nc2u),lat_c2u(samp%nc2u)
logical :: dichotomy,convergence,valid_c2u(samp%nc2u),valid_c2a(samp%nc2a)
character(len=1024),parameter :: subr = 'adv_filter'
type(mesh_type) :: mesh

if ((nam%adv_niter>0).and.(adv%valid_raw(il0,its)<nam%adv_valid)) then
   ! Convert to cartesian coordinates
   do ic2a=1,samp%nc2a
      if (samp%smask_c2a(ic2a,il0)) then
         call lonlat2xyz(mpl,adv%lon_c2a_raw(ic2a,il0,its),adv%lat_c2a_raw(ic2a,il0,its),dx_ini(ic2a),dy_ini(ic2a), &
 & dz_ini(ic2a))
         if (mpl%msv%isnot(dx_ini(ic2a))) dx_ini(ic2a) = dx_ini(ic2a)-adv%x_c2a(ic2a)
         if (mpl%msv%isnot(dy_ini(ic2a))) dy_ini(ic2a) = dy_ini(ic2a)-adv%y_c2a(ic2a)
         if (mpl%msv%isnot(dz_ini(ic2a))) dz_ini(ic2a) = dz_ini(ic2a)-adv%z_c2a(ic2a)
      else
         dx_ini(ic2a) = mpl%msv%valr
         dy_ini(ic2a) = mpl%msv%valr
         dz_ini(ic2a) = mpl%msv%valr
      end if
   end do

   ! Dichotomy initialization
   convergence = .false.
   dichotomy = .false.
   rhflt = nam%adv_rhflt
   drhflt = rhflt
   adv%rhflt(il0,its) = huge_real

   do iter=1,nam%adv_niter
      ! Copy advection
      dx = dx_ini
      dy = dy_ini
      dz = dz_ini

      ! Median filter to remove extreme values
      dd = dx**2+dy**2+dz**2
      call samp%diag_filter(mpl,nam,'median',rhflt,dx,dd)
      call samp%diag_filter(mpl,nam,'median',rhflt,dy,dd)
      call samp%diag_filter(mpl,nam,'median',rhflt,dz,dd)

      ! Average filter to smooth values
      call samp%diag_filter(mpl,nam,'gc99',rhflt,dx)
      call samp%diag_filter(mpl,nam,'gc99',rhflt,dy)
      call samp%diag_filter(mpl,nam,'gc99',rhflt,dz)

      ! Back to spherical coordinates
      dx = adv%x_c2a+dx
      dy = adv%y_c2a+dy
      dz = adv%z_c2a+dz
      do ic2a=1,samp%nc2a
         if (samp%smask_c2a(ic2a,il0)) then
            call xyz2lonlat(mpl,dx(ic2a),dy(ic2a),dz(ic2a),lon_c2a(ic2a),lat_c2a(ic2a))
         else
            lon_c2a(ic2a) = adv%lon_c2a(ic2a)
            lat_c2a(ic2a) = adv%lat_c2a(ic2a)
         end if
      end do

      ! Reduce distance with respect to the boundary
      il0i = min(il0,geom%nl0i)
      do ic2a=1,samp%nc2a
         ic0a = samp%c2a_to_c0a(ic2a)
         ic0u = geom%c0a_to_c0u(ic0a)
         call reduce_arc(adv%lon_c2a(ic2a),adv%lat_c2a(ic2a),lon_c2a(ic2a),lat_c2a(ic2a), &
 & min(geom%mdist_c0u(ic0u,il0i),geom%mesh_c0u%bdist(ic0u)),dist_c2a(ic2a))
      end do

      ! Check filtered mesh
      call samp%com_AU%ext(mpl,lon_c2a,lon_c2u)
      call samp%com_AU%ext(mpl,lat_c2a,lat_c2u)
      call mesh%copy(samp%mesh)
      call mesh%store(mpl,lon_c2u,lat_c2u)
      call mesh%check(mpl,valid_c2u)
      do ic2u=1,samp%nc2u
         ic0u = samp%c2u_to_c0u(ic2u)
         ic0 = geom%c0u_to_c0(ic0u)
         iproc = geom%c0_to_proc(ic0)
         if (iproc==mpl%myproc) then
            valid_c2u(ic2u) = valid_c2u(ic2u).and.samp%mesh%valid(ic2u).and.samp%smask_c2u(ic2u,il0)
         else
            valid_c2u(ic2u) = .true.
         end if
      end do
      call samp%com_AU%red(mpl,valid_c2u,valid_c2a,op='and')
      nvalid = count(valid_c2a)
      nmask = count(samp%mesh%valid(samp%c2a_to_c2u).and.samp%smask_c2a(:,il0))
      call mpl%f_comm%allreduce(nvalid,nvalid_tot,fckit_mpi_sum())
      call mpl%f_comm%allreduce(nmask,nmask_tot,fckit_mpi_sum())
      valid_flt = real(nvalid_tot,kind_real)/real(nmask_tot,kind_real)

      ! Print result
      write(mpl%info,'(a19,a,i2,a,f10.2,a,f7.3,a)') '','Iteration ',iter,': rhflt = ', &
 & rhflt*reqkm,' km, validity = ',100.0*valid_flt,'%'
      call mpl%flush

      ! Update support radius
      if (inf(valid_flt,nam%adv_valid)) then
         ! Increase filtering support radius
         if (dichotomy) then
            drhflt = 0.5*drhflt
            if (iter<nam%adv_niter) rhflt = rhflt+drhflt
         else
            convergence = .false.
            if (iter<nam%adv_niter) rhflt = rhflt+drhflt
            drhflt = 2.0*drhflt
         end if
      else
         if (inf(rhflt,adv%rhflt(il0,its))) then
            ! Average distance
            dist_sum = sum(dist_c2a,mask=mpl%msv%isnot(dist_c2a))
            norm = real(count(mpl%msv%isnot(dist_c2a)),kind_real)
            call mpl%f_comm%allreduce(dist_sum,dist_flt,fckit_mpi_sum())
            call mpl%f_comm%allreduce(norm,norm_tot,fckit_mpi_sum())
            dist_flt = dist_flt/norm_tot

            ! Save values
            adv%lon_c2a_flt(:,il0,its) = lon_c2a
            adv%lat_c2a_flt(:,il0,its) = lat_c2a
            adv%dist_c2a_flt(:,il0,its) = dist_c2a
            adv%valid_flt(il0,its) = valid_flt
            adv%dist_flt(il0,its) = dist_flt
            adv%rhflt(il0,its) = rhflt
         end if

         ! Convergence
         convergence = .true.

         ! Change dichotomy status
         if (.not.dichotomy) then
            dichotomy = .true.
            drhflt = 0.5*drhflt
         end if

         ! Decrease filtering support radius
         drhflt = 0.5*drhflt
         if (iter<nam%adv_niter) rhflt = rhflt-drhflt
      end if
   end do

   ! Check convergence
   if (.not.convergence) call mpl%abort(subr,'convergence failed in adv%compute')
else
   ! Advection filtering not required
   adv%lon_c2a_flt(:,il0,its) = adv%lon_c2a_raw(:,il0,its)
   adv%lat_c2a_flt(:,il0,its) = adv%lat_c2a_raw(:,il0,its)
   adv%dist_c2a_flt(:,il0,its) = adv%dist_c2a_raw(:,il0,its)
   adv%valid_flt(il0,its) = adv%valid_raw(il0,its)
   adv%dist_flt(il0,its) = adv%dist_raw(il0,its)
   adv%rhflt(il0,its) = 0.0
end if

if (nam%adv_niter>0) then
   ! Print results
   write(mpl%info,'(a19,a,f10.2,a)') '','Optimal filtering support radius: ',adv%rhflt(il0,its)*reqkm,' km'
   call mpl%flush
   write(mpl%info,'(a19,a,f5.1,a,f5.1,a)') '','Valid points:             ',100.0*adv%valid_raw(il0,its),'%  ~> ', &
 & 100.0*adv%valid_flt(il0,its),'%'
   call mpl%flush
    write(mpl%info,'(a19,a,f10.1,a,f10.1,a)') '','Advection distance: ',adv%dist_raw(il0,its)*reqkm,' km ~> ', &
 & adv%dist_flt(il0,its)*reqkm,' km'
   call mpl%flush
end if

end subroutine adv_filter

!----------------------------------------------------------------------
! Subroutine: adv_interp
! Purpose: interpolate advection
!----------------------------------------------------------------------
subroutine adv_interp(adv,mpl,nam,geom,samp)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv  ! Advection
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(samp_type),intent(inout) :: samp ! Sampling

! Local variables
integer :: ic0a,ic0u,ic2a,il0,il0i,its
real(kind_real) :: reduced_dist
real(kind_real) :: dx(samp%nc2a),dy(samp%nc2a),dz(samp%nc2a)
real(kind_real) :: dx_c2b(samp%nc2b),dy_c2b(samp%nc2b),dz_c2b(samp%nc2b)
real(kind_real) :: dx_c0a(geom%nc0a),dy_c0a(geom%nc0a),dz_c0a(geom%nc0a)
real(kind_real) :: x_ori_c0a(geom%nc0a),y_ori_c0a(geom%nc0a),z_ori_c0a(geom%nc0a)

write(mpl%info,'(a7,a)') '','Interpolate advection'
call mpl%flush

! Initialization
do ic0a=1,geom%nc0a
   call lonlat2xyz(mpl,geom%lon_c0a(ic0a),geom%lat_c0a(ic0a),x_ori_c0a(ic0a),y_ori_c0a(ic0a),z_ori_c0a(ic0a))
end do

! Advection for timeslot 1
do il0=1,geom%nl0
   do ic0a=1,geom%nc0a
      samp%adv_lon(ic0a,il0,1) = geom%lon_c0a(ic0a)
      samp%adv_lat(ic0a,il0,1) = geom%lat_c0a(ic0a)
   end do
end do

do its=2,nam%nts
   do il0=1,geom%nl0
      ! Convert to cartesian coordinates
      do ic2a=1,samp%nc2a
         if (samp%smask_c2a(ic2a,il0)) then
            call lonlat2xyz(mpl,adv%lon_c2a_flt(ic2a,il0,its),adv%lat_c2a_flt(ic2a,il0,its),dx(ic2a),dy(ic2a),dz(ic2a))
            if (mpl%msv%isnot(dx(ic2a))) dx(ic2a) = dx(ic2a)-adv%x_c2a(ic2a)
            if (mpl%msv%isnot(dy(ic2a))) dy(ic2a) = dy(ic2a)-adv%y_c2a(ic2a)
            if (mpl%msv%isnot(dz(ic2a))) dz(ic2a) = dz(ic2a)-adv%z_c2a(ic2a)
         else
            dx(ic2a) = 0.0
            dy(ic2a) = 0.0
            dz(ic2a) = 0.0
         end if
      end do

      ! Advection interpolation
      il0i = min(il0,geom%nl0i)
      call samp%com_AB%ext(mpl,dx,dx_c2b)
      call samp%com_AB%ext(mpl,dy,dy_c2b)
      call samp%com_AB%ext(mpl,dz,dz_c2b)
      call samp%h(il0i)%apply(mpl,dx_c2b,dx_c0a)
      call samp%h(il0i)%apply(mpl,dy_c2b,dy_c0a)
      call samp%h(il0i)%apply(mpl,dz_c2b,dz_c0a)

      ! Back to spherical coordinates
      dx_c0a = x_ori_c0a+dx_c0a
      dy_c0a = y_ori_c0a+dy_c0a
      dz_c0a = z_ori_c0a+dz_c0a
      do ic0a=1,geom%nc0a
         call xyz2lonlat(mpl,dx_c0a(ic0a),dy_c0a(ic0a),dz_c0a(ic0a),samp%adv_lon(ic0a,il0,its),samp%adv_lat(ic0a,il0,its))
      end do

      ! Reduce distance with respect to the boundary
      il0i = min(il0,geom%nl0i)
      do ic0a=1,geom%nc0a
         if (geom%gmask_c0a(ic0a,il0)) then
            ic0u = geom%c0a_to_c0u(ic0a)
            call reduce_arc(geom%lon_c0a(ic0a),geom%lat_c0a(ic0a),samp%adv_lon(ic0a,il0,its),samp%adv_lat(ic0a,il0,its), &
 & min(geom%mdist_c0a(ic0a,il0i),geom%mesh_c0u%bdist(ic0u)),reduced_dist)
         end if
      end do

      ! Deal with poles issues (back to origin point latitude by 10%)
      do ic0a=1,geom%nc0a
         if (eq(abs(samp%adv_lat(ic0a,il0,its)),0.5*pi)) samp%adv_lat(ic0a,il0,its) = geom%lat_c0a(ic0a) &
 & +0.9*(samp%adv_lat(ic0a,il0,its)-geom%lat_c0a(ic0a))
      end do
   end do
end do

end subroutine adv_interp

!----------------------------------------------------------------------
! Subroutine: adv_test
! Purpose: test advection efficiency
!----------------------------------------------------------------------
subroutine adv_test(adv,mpl,rng,nam,geom,samp,ens)

implicit none

! Passed variables
class(adv_type),intent(inout) :: adv  ! Advection
type(mpl_type),intent(inout) :: mpl   ! MPI data
type(rng_type),intent(inout) :: rng   ! Random number generator
type(nam_type),intent(in) :: nam      ! Namelist
type(geom_type),intent(in) :: geom    ! Geometry
type(samp_type),intent(inout) :: samp ! Sampling
type(ens_type), intent(in) :: ens     ! Ensemble

! Local variables
integer :: its,il0,ic0,ic0u,ic0a,i_s,jc0u,nc0d,ie,iv,ic0d
integer :: c0u_to_c0d(geom%nc0u)
integer,allocatable :: c0d_to_c0(:)
real(kind_real) :: score_loc,score_adv,norm_loc,norm_adv,norm_loc_tot,norm_adv_tot
real(kind_real) :: m2_1(geom%nc0a,geom%nl0),m2_2_loc(geom%nc0a,geom%nl0),m11_loc(geom%nc0a,geom%nl0)
real(kind_real) :: m2_2_adv(geom%nc0a,geom%nl0),m11_adv(geom%nc0a,geom%nl0)
real(kind_real) :: fld(geom%nc0a,geom%nl0)
real(kind_real),allocatable :: fld_d(:,:)
logical :: lcheck_c0d(geom%nc0u)
type(com_type) :: com_ADinv
type(linop_type) :: dinv(geom%nl0,2:nam%nts)

write(mpl%info,'(a7,a)') '','Test advection'
call mpl%flush

! Compute interpolation
do its=2,nam%nts
   do il0=1,geom%nl0
      write(dinv(il0,its)%prefix,'(a,i3.3,a,i2.2)') 'd_',il0,'_',its
      call dinv(il0,its)%interp(mpl,rng,nam,geom,il0,geom%nc0u,geom%lon_c0u,geom%lat_c0u,geom%gmask_c0u(:,il0),geom%nc0a, &
 & samp%adv_lon(:,il0,its),samp%adv_lat(:,il0,its),geom%gmask_c0a(:,il0),10)
   end do
end do

! Define halo D
lcheck_c0d = .false.
do ic0a=1,geom%nc0a
   ic0u = geom%c0a_to_c0u(ic0a)
   lcheck_c0d(ic0u) = .true.
end do
do its=2,nam%nts
   do il0=1,geom%nl0
      do i_s=1,dinv(il0,its)%n_s
         jc0u = dinv(il0,its)%col(i_s)
         lcheck_c0d(jc0u) = .true.
      end do
   end do
end do
nc0d = count(lcheck_c0d)

! Allocation
allocate(c0d_to_c0(nc0d))

! Global-local conversions for halo D
c0u_to_c0d = mpl%msv%vali
ic0d = 0
do ic0u=1,geom%nc0u
   if (lcheck_c0d(ic0u)) then
      ic0d = ic0d+1
      ic0 = geom%c0u_to_c0(ic0u)
      c0d_to_c0(ic0d) = ic0
      c0u_to_c0d(ic0u) = ic0d
   end if
end do

! Local interpolation source
do its=2,nam%nts
   do il0=1,geom%nl0
      dinv(il0,its)%n_src = nc0d
      do i_s=1,dinv(il0,its)%n_s
         dinv(il0,its)%col(i_s) = c0u_to_c0d(dinv(il0,its)%col(i_s))
      end do
   end do
end do

! Setup communications
call com_ADinv%setup(mpl,'com_ADinv',geom%nc0a,nc0d,geom%nc0,geom%c0a_to_c0,c0d_to_c0)

! Allocation
allocate(fld_d(nc0d,geom%nl0))

! Compute correlation at zero separation
do its=2,nam%nts
   do iv=1,nam%nv
      ! Initialization
      m2_1 = 0.0
      m2_2_loc = 0.0
      m2_2_adv = 0.0
      m11_loc = 0.0
      m11_adv = 0.0

      do ie=1,ens%ne
         ! Halo extension from zone A to zone D
         call com_ADinv%ext(mpl,geom%nl0,ens%mem(ie)%fld(:,:,iv,its),fld_d)

         ! Interpolation
         !$omp parallel do schedule(static) private(il0)
         do il0=1,geom%nl0
            call dinv(il0,its)%apply(mpl,fld_d(:,il0),fld(:,il0),msdst=.false.)
         end do
         !$omp end parallel do

         ! Compute variance/covariance
         m2_1 = m2_1+ens%mem(ie)%fld(:,:,iv,1)**2
         m2_2_loc = m2_2_loc+ens%mem(ie)%fld(:,:,iv,its)**2
         m2_2_adv = m2_2_adv+fld**2
         m11_loc = m11_loc+ens%mem(ie)%fld(:,:,iv,1)*ens%mem(ie)%fld(:,:,iv,its)
         m11_adv = m11_adv+ens%mem(ie)%fld(:,:,iv,1)*fld
      end do

      ! Compute correlation
      do il0=1,geom%nl0
         do ic0a=1,geom%nc0a
            if ((m2_1(ic0a,il0)>0.0).and.(m2_2_loc(ic0a,il0)>0.0)) then
               adv%cor_loc(ic0a,il0,iv,its) = m11_loc(ic0a,il0)/sqrt(m2_1(ic0a,il0)*m2_2_loc(ic0a,il0))
            else
               adv%cor_loc(ic0a,il0,iv,its) = mpl%msv%valr
            end if
            if ((m2_1(ic0a,il0)>0.0).and.(m2_2_adv(ic0a,il0)>0.0)) then
               adv%cor_adv(ic0a,il0,iv,its) = m11_adv(ic0a,il0)/sqrt(m2_1(ic0a,il0)*m2_2_adv(ic0a,il0))
            else
               adv%cor_adv(ic0a,il0,iv,its) = mpl%msv%valr
            end if
         end do
      end do
   end do
end do

! Print results
do its=2,nam%nts
   score_loc = sum(adv%cor_loc(:,:,:,its),mask=mpl%msv%isnot(adv%cor_loc(:,:,:,its)))
   score_adv = sum(adv%cor_adv(:,:,:,its),mask=mpl%msv%isnot(adv%cor_adv(:,:,:,its)))
   norm_loc = real(count(mpl%msv%isnot(adv%cor_loc(:,:,:,its))),kind_real)
   norm_adv = real(count(mpl%msv%isnot(adv%cor_adv(:,:,:,its))),kind_real)
   call mpl%f_comm%allreduce(score_loc,adv%score_loc(its),fckit_mpi_sum())
   call mpl%f_comm%allreduce(score_adv,adv%score_adv(its),fckit_mpi_sum())
   call mpl%f_comm%allreduce(norm_loc,norm_loc_tot,fckit_mpi_sum())
   call mpl%f_comm%allreduce(norm_adv,norm_adv_tot,fckit_mpi_sum())
   adv%score_loc(its) = adv%score_loc(its)/norm_loc_tot
   adv%score_adv(its) = adv%score_adv(its)/norm_adv_tot
   write(mpl%info,'(a10,a,i2,a,f6.3,a,f6.3,a,f6.3)') '','Timeslot ',its,': ',adv%score_loc(its),' ~> ',adv%score_adv(its)
   call mpl%flush
end do

end subroutine adv_test

!----------------------------------------------------------------------
! Subroutine: adv_write
! Purpose: write advection
!----------------------------------------------------------------------
subroutine adv_write(adv,mpl,nam,geom,bpar,io,samp)

implicit none

! Passed variables
class(adv_type),intent(in) :: adv   ! Advection
type(mpl_type),intent(inout) :: mpl ! MPI data
type(nam_type),intent(in) :: nam    ! Namelist
type(geom_type),intent(in) :: geom  ! Geometry
type(bpar_type),intent(in) :: bpar  ! Block parameters
type(io_type),intent(in) :: io      ! I/O
type(samp_type),intent(in) :: samp  ! Sampling

! Local variables
integer :: ncid,nc2_id,nl0_id,nts_id,na_id
integer :: lon_c2_id,lat_c2_id,lon_c2_raw_id,lat_c2_raw_id,dist_c2_raw_id,valid_raw_id,dist_raw_id
integer :: lon_c2_flt_id,lat_c2_flt_id,dist_c2_flt_id,valid_flt_id,dist_flt_id,rhflt_id
integer :: score_loc_id,score_adv_id
integer :: iproc,its,il0,ic2a,ic2,i,ib,iv,jv,jts
integer,allocatable :: order(:)
real(kind_real),allocatable :: sbuf(:),rbuf(:),lon_c2(:,:),lat_c2(:,:),hash_c2(:)
real(kind_real),allocatable :: lon_c2_raw(:,:,:),lat_c2_raw(:,:,:),dist_c2_raw(:,:,:)
real(kind_real),allocatable :: lon_c2_flt(:,:,:),lat_c2_flt(:,:,:),dist_c2_flt(:,:,:)
character(len=1024) :: filename
character(len=1024),parameter :: subr = 'adv_write'
type(fckit_mpi_status) :: status

write(mpl%info,'(a7,a)') '','Write advection'
call mpl%flush

! Allocation
allocate(sbuf(samp%nc2a*geom%nl0*(2+(nam%nts-1)*6)))

! Initialization
sbuf = mpl%msv%valr

! Prepare buffer
i = 1
do il0=1,geom%nl0
   do ic2a=1,samp%nc2a
      ! Apply mask
      if (samp%smask_c2a(ic2a,il0)) sbuf(i) = adv%lon_c2a(ic2a)
      i = i+1
      if (samp%smask_c2a(ic2a,il0)) sbuf(i) = adv%lat_c2a(ic2a)
      i = i+1
      do its=2,nam%nts
         if (samp%smask_c2a(ic2a,il0)) sbuf(i) = adv%lon_c2a_raw(ic2a,il0,its)
         i = i+1
         if (samp%smask_c2a(ic2a,il0)) sbuf(i) = adv%lat_c2a_raw(ic2a,il0,its)
         i = i+1
         if (samp%smask_c2a(ic2a,il0)) sbuf(i) = adv%dist_c2a_raw(ic2a,il0,its)
         i = i+1
         if (samp%smask_c2a(ic2a,il0)) sbuf(i) = adv%lon_c2a_flt(ic2a,il0,its)
         i = i+1
         if (samp%smask_c2a(ic2a,il0)) sbuf(i) = adv%lat_c2a_flt(ic2a,il0,its)
         i = i+1
         if (samp%smask_c2a(ic2a,il0)) sbuf(i) = adv%dist_c2a_flt(ic2a,il0,its)
         i = i+1
      end do
   end do
end do

if (mpl%main) then
   ! Allocation
   allocate(lon_c2(nam%nc2,geom%nl0))
   allocate(lat_c2(nam%nc2,geom%nl0))
   allocate(lon_c2_raw(nam%nc2,geom%nl0,nam%nts-1))
   allocate(lat_c2_raw(nam%nc2,geom%nl0,nam%nts-1))
   allocate(dist_c2_raw(nam%nc2,geom%nl0,nam%nts-1))
   allocate(lon_c2_flt(nam%nc2,geom%nl0,nam%nts-1))
   allocate(lat_c2_flt(nam%nc2,geom%nl0,nam%nts-1))
   allocate(dist_c2_flt(nam%nc2,geom%nl0,nam%nts-1))
   allocate(hash_c2(nam%nc2))
   allocate(order(nam%nc2))

   do iproc=1,mpl%nproc
      ! Allocation
      allocate(rbuf(samp%proc_to_nc2a(iproc)*geom%nl0*(2+(nam%nts-1)*6)))

      if (iproc==mpl%rootproc) then
         ! Copy buffer
         rbuf = sbuf
      else
         ! Receive data on rootproc
         call mpl%f_comm%receive(rbuf,iproc-1,mpl%tag,status)
      end if

      ! Copy data
      i = 1
      do il0=1,geom%nl0
         do ic2a=1,samp%proc_to_nc2a(iproc)
            ic2 = samp%proc_to_c2_offset(iproc)+ic2a
            lon_c2(ic2,il0) = rbuf(i)
            i = i+1
            lat_c2(ic2,il0) = rbuf(i)
            i = i+1
            do its=2,nam%nts
               lon_c2_raw(ic2,il0,its-1) = rbuf(i)
               i = i+1
               lat_c2_raw(ic2,il0,its-1) = rbuf(i)
               i = i+1
               dist_c2_raw(ic2,il0,its-1) = rbuf(i)
               i = i+1
               lon_c2_flt(ic2,il0,its-1) = rbuf(i)
               i = i+1
               lat_c2_flt(ic2,il0,its-1) = rbuf(i)
               i = i+1
               dist_c2_flt(ic2,il0,its-1) = rbuf(i)
               i = i+1
            end do
         end do
      end do

      ! Release memory
      deallocate(rbuf)
   end do

   do il0=1,geom%nl0
      ! Get order from hash
      do ic2=1,nam%nc2
         hash_c2(ic2) = lonlathash(lon_c2(ic2,il0),lat_c2(ic2,il0))
      end do
      call qsort(nam%nc2,hash_c2,order)

      ! Reorder data
      lon_c2(:,il0) = lon_c2(order,il0)
      lat_c2(:,il0) = lat_c2(order,il0)
      lon_c2_raw(:,il0,:) = lon_c2_raw(order,il0,:)
      lat_c2_raw(:,il0,:) = lat_c2_raw(order,il0,:)
      dist_c2_raw(:,il0,:) = dist_c2_raw(order,il0,:)
      lon_c2_flt(:,il0,:) = lon_c2_flt(order,il0,:)
      lat_c2_flt(:,il0,:) = lat_c2_flt(order,il0,:)
      dist_c2_flt(:,il0,:) = dist_c2_flt(order,il0,:)
   end do

   ! Rescale data
   lon_c2 = lon_c2*rad2deg
   lat_c2 = lat_c2*rad2deg
   do its=1,nam%nts-1
      do il0=1,geom%nl0
         do ic2=1,nam%nc2
            if (mpl%msv%isnot(lon_c2_raw(ic2,il0,its))) lon_c2_raw(ic2,il0,its) = lon_c2_raw(ic2,il0,its)*rad2deg
            if (mpl%msv%isnot(lat_c2_raw(ic2,il0,its))) lat_c2_raw(ic2,il0,its) = lat_c2_raw(ic2,il0,its)*rad2deg
            if (mpl%msv%isnot(dist_c2_raw(ic2,il0,its))) dist_c2_raw(ic2,il0,its) = dist_c2_raw(ic2,il0,its)*req
            if (mpl%msv%isnot(lon_c2_flt(ic2,il0,its))) lon_c2_flt(ic2,il0,its) = lon_c2_flt(ic2,il0,its)*rad2deg
            if (mpl%msv%isnot(lat_c2_flt(ic2,il0,its))) lat_c2_flt(ic2,il0,its) = lat_c2_flt(ic2,il0,its)*rad2deg
            if (mpl%msv%isnot(dist_c2_flt(ic2,il0,its))) dist_c2_flt(ic2,il0,its) = dist_c2_flt(ic2,il0,its)*req
         end do
      end do
   end do

   ! Release memory
   deallocate(hash_c2)
   deallocate(order)
else
   ! Send data to rootproc
   call mpl%f_comm%send(sbuf,mpl%rootproc-1,mpl%tag)
end if
call mpl%update_tag(1)

! Release memory
deallocate(sbuf)

if (mpl%main) then
   ! Define file
   filename = trim(nam%prefix)//'_adv_diag'
   ncid = mpl%nc_file_create_or_open(subr,trim(nam%datadir)//'/'//trim(filename)//'.nc')

   ! Write namelist parameters
   call nam%write(mpl,ncid)

   ! Define dimensions
   nc2_id = mpl%nc_dim_define_or_get(subr,ncid,'nc2',nam%nc2)
   nl0_id = mpl%nc_dim_define_or_get(subr,ncid,'nl0',geom%nl0)
   nts_id = mpl%nc_dim_define_or_get(subr,ncid,'nts',nam%nts-1)
   na_id = mpl%nc_dim_define_or_get(subr,ncid,'na',samp%mesh%na)

   ! Define variables
   lon_c2_id = mpl%nc_var_define_or_get(subr,ncid,'lon_c2',nc_kind_real,(/nc2_id,nl0_id/))
   lat_c2_id = mpl%nc_var_define_or_get(subr,ncid,'lat_c2',nc_kind_real,(/nc2_id,nl0_id/))
   lon_c2_raw_id = mpl%nc_var_define_or_get(subr,ncid,'lon_c2_raw',nc_kind_real,(/nc2_id,nl0_id,nts_id/))
   lat_c2_raw_id = mpl%nc_var_define_or_get(subr,ncid,'lat_c2_raw',nc_kind_real,(/nc2_id,nl0_id,nts_id/))
   dist_c2_raw_id = mpl%nc_var_define_or_get(subr,ncid,'dist_c2_raw',nc_kind_real,(/nc2_id,nl0_id,nts_id/))
   valid_raw_id = mpl%nc_var_define_or_get(subr,ncid,'valid_raw',nc_kind_real,(/nl0_id,nts_id/))
   dist_raw_id = mpl%nc_var_define_or_get(subr,ncid,'dist_raw',nc_kind_real,(/nl0_id,nts_id/))
   lon_c2_flt_id = mpl%nc_var_define_or_get(subr,ncid,'lon_c2_flt',nc_kind_real,(/nc2_id,nl0_id,nts_id/))
   lat_c2_flt_id = mpl%nc_var_define_or_get(subr,ncid,'lat_c2_flt',nc_kind_real,(/nc2_id,nl0_id,nts_id/))
   dist_c2_flt_id = mpl%nc_var_define_or_get(subr,ncid,'dist_c2_flt',nc_kind_real,(/nc2_id,nl0_id,nts_id/))
   valid_flt_id = mpl%nc_var_define_or_get(subr,ncid,'valid_flt',nc_kind_real,(/nl0_id,nts_id/))
   dist_flt_id = mpl%nc_var_define_or_get(subr,ncid,'dist_flt',nc_kind_real,(/nl0_id,nts_id/))
   rhflt_id = mpl%nc_var_define_or_get(subr,ncid,'rhflt',nc_kind_real,(/nl0_id,nts_id/))
   score_loc_id = mpl%nc_var_define_or_get(subr,ncid,'score_loc',nc_kind_real,(/nts_id/))
   score_adv_id = mpl%nc_var_define_or_get(subr,ncid,'score_adv',nc_kind_real,(/nts_id/))

   ! Write data
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c2_id,lon_c2))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c2_id,lat_c2))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c2_raw_id,lon_c2_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c2_raw_id,lat_c2_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_c2_raw_id,dist_c2_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,valid_raw_id,adv%valid_raw))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_raw_id,adv%dist_raw*req))
   call mpl%ncerr(subr,nf90_put_var(ncid,lon_c2_flt_id,lon_c2_flt))
   call mpl%ncerr(subr,nf90_put_var(ncid,lat_c2_flt_id,lat_c2_flt))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_c2_flt_id,dist_c2_flt))
   call mpl%ncerr(subr,nf90_put_var(ncid,valid_flt_id,adv%valid_flt))
   call mpl%ncerr(subr,nf90_put_var(ncid,dist_flt_id,adv%dist_flt*req))
   call mpl%ncerr(subr,nf90_put_var(ncid,rhflt_id,adv%rhflt*req))
   call mpl%ncerr(subr,nf90_put_var(ncid,score_loc_id,adv%score_loc))
   call mpl%ncerr(subr,nf90_put_var(ncid,score_adv_id,adv%score_adv))

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))

   ! Release memory
   deallocate(lon_c2)
   deallocate(lat_c2)
   deallocate(lon_c2_raw)
   deallocate(lat_c2_raw)
   deallocate(dist_c2_raw)
   deallocate(lon_c2_flt)
   deallocate(lat_c2_flt)
   deallocate(dist_c2_flt)
end if

! Set file name
filename = trim(nam%prefix)//'_adv_test'

! Write test
do ib=1,bpar%nb
   iv = bpar%b_to_v1(ib)
   jv = bpar%b_to_v2(ib)
   its = bpar%b_to_ts1(ib)
   jts = bpar%b_to_ts2(ib)
   if ((iv==jv).and.(its>=2).and.(its==jts)) then
      call io%fld_write(mpl,nam,geom,filename,'cor_loc',adv%cor_loc(:,:,iv,its),bpar%blockname(ib))
      call io%fld_write(mpl,nam,geom,filename,'cor_adv',adv%cor_adv(:,:,iv,its),bpar%blockname(ib))
   end if
end do

end subroutine adv_write

end module type_adv
