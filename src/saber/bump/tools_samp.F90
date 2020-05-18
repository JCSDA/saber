!----------------------------------------------------------------------
! Module: tools_samp
! Purpose: sampling functions
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_samp

use atlas_module
use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_status
use tools_const, only: pi,deg2rad
use tools_func, only: lonlathash,lonlatmod,sphere_dist
use tools_kinds, only: kind_real,kind_int
use tools_qsort, only: qsort
use tools_repro, only: repro,inf,sup,eq
use type_mpl, only: mpl_type
use type_rng, only: rng_type
use type_tree, only: tree_type

implicit none

private
public :: initialize_sampling

contains

!----------------------------------------------------------------------
! Subroutine: initialize_sampling
! Purpose: intialize sampling
!----------------------------------------------------------------------
subroutine initialize_sampling(mpl,rng,area,n_loc,lon_loc,lat_loc,mask_loc,rh_loc,loc_to_glb,ntry,nrep,ns2_glb,sam2_glb, &
 & fast,verbosity,n_uni,uni_to_proc,uni_to_loc,tree_uni)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl               ! MPI data
type(rng_type),intent(inout) :: rng               ! Random number generator
real(kind_real),intent(in) :: area                ! Global domain area
integer,intent(in) :: n_loc                       ! Number of points (local)
real(kind_real),intent(in) :: lon_loc(n_loc)      ! Longitudes (local)
real(kind_real),intent(in) :: lat_loc(n_loc)      ! Latitudes (local)
logical,intent(in) :: mask_loc(n_loc)             ! Mask (local)
real(kind_real),intent(in) :: rh_loc(n_loc)       ! Horizontal support radius (local)
integer,intent(in) :: loc_to_glb(n_loc)           ! Local to global index
integer,intent(in) :: ntry                        ! Number of tries
integer,intent(in) :: nrep                        ! Number of replacements
integer,intent(in) :: ns2_glb                     ! Number of samplings points (global)
integer,intent(out) :: sam2_glb(ns2_glb)          ! Horizontal sampling index (global)
logical,intent(in),optional :: fast               ! Fast sampling flag
logical,intent(in),optional :: verbosity          ! Verbosity flag
integer,intent(in),optional :: n_uni              ! Universe size
integer,intent(in),optional :: uni_to_proc(:)     ! Universe to processor
integer,intent(in),optional :: uni_to_loc(:)      ! Universe to local index
type(tree_type),intent(in),optional :: tree_uni   ! Universe KD-tree

! Local variables
integer :: n_glb,n_loc_eff,n_glb_eff,i_glb,i_loc,i_loc_eff,n,ix,iy,is2_glb,iproc,js,irep,irmax,itry,is1_glb,ir,ns1_loc,is1_loc,nfac
integer :: ns1_loc_tmp,irval,irvalmin,irvalmax,is2_glb_min,nrep_eff,nn_index(2),ns1_glb,is1_glb_eff,ns1_glb_eff,ns1_glb_val,offset
integer,allocatable :: glb_to_loc(:),glb_to_proc(:)
integer,allocatable :: s1_loc_to_glb(:),s1_loc_to_glb_tmp(:),sam1_loc(:),sam1_loc_tmp(:)
integer,allocatable :: sam1_glb(:),sam1_glb_eff(:),to_valid(:),sam2_glb_tmp(:),order(:)
real(kind_real) :: lonlat(2),d,distmax,distmin,nn_dist(2),cdf_norm,rr
real(kind_real),allocatable :: hash_glb(:),hash_loc(:)
real(kind_real),allocatable :: lon1_loc(:),lat1_loc(:),dist1_loc(:),rh1_loc(:)
real(kind_real),allocatable :: lon1_loc_tmp(:),lat1_loc_tmp(:),dist1_loc_tmp(:),rh1_loc_tmp(:)
real(kind_real),allocatable :: lon1_glb(:),lat1_glb(:),dist1_glb(:),rh1_glb(:)
real(kind_real),allocatable :: lon1_glb_eff(:),lat1_glb_eff(:),rh1_glb_eff(:)
real(kind_real),allocatable :: list(:),cdf(:)
real(kind_real),allocatable :: lon_rep(:),lat_rep(:),dist(:)
logical :: lfull_grid,lfast,lverbosity,update,retry
logical,allocatable :: mask_glb(:),lmask(:),smask(:),rmask(:)
character(len=6) :: gridid
character(len=1024),parameter :: subr = 'initialize_sampling'
type(fckit_mpi_status) :: status
type(tree_type) :: tree
type(atlas_structuredgrid) :: agrid

! Local flags
lfull_grid = present(n_uni).and.present(uni_to_proc).and.present(uni_to_loc).and.present(tree_uni)
lfast = .false.
lverbosity = .true.
if (present(fast)) lfast = fast
if (present(verbosity)) lverbosity = verbosity

! Global size
call mpl%f_comm%allreduce(n_loc,n_glb,fckit_mpi_sum())

! Number of effective points
n_loc_eff = count(mask_loc)
call mpl%f_comm%allreduce(n_loc_eff,n_glb_eff,fckit_mpi_sum())

! Check mask size
if (n_glb_eff==0) then
   call mpl%abort(subr,'empty mask in initialize sampling')
elseif (n_glb_eff<ns2_glb) then
   call mpl%abort(subr,'ns2_glb greater than n_glb_eff in initialize_sampling')
elseif (n_glb_eff==ns2_glb) then
   write(mpl%info,'(a)') ' all points are used'
   if (lverbosity) call mpl%flush

   ! Allocation
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
   allocate(hash_loc(n_loc))
   if (mpl%main) then
      allocate(hash_glb(n_glb))
      allocate(mask_glb(n_glb))
      allocate(list(ns2_glb))
      allocate(order(ns2_glb))
   else
      allocate(hash_glb(0))
      allocate(mask_glb(0))
   end if

   ! Compute hash
   do i_loc=1,n_loc
      hash_loc(i_loc) = lonlathash(lon_loc(i_loc),lat_loc(i_loc))
   end do

   ! Communication
   call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)
   call mpl%loc_to_glb(n_loc,hash_loc,n_glb,glb_to_proc,glb_to_loc,.false.,hash_glb)
   call mpl%loc_to_glb(n_loc,mask_loc,n_glb,glb_to_proc,glb_to_loc,.false.,mask_glb)

   if (mpl%main) then
      ! Use all valid points
      is2_glb = 0
      do i_glb=1,n_glb
         if (mask_glb(i_glb)) then
            is2_glb = is2_glb+1
            sam2_glb(is2_glb) = i_glb
            list(is2_glb) = hash_glb(i_glb)
         end if
      end do

      ! Define points order
      call qsort(ns2_glb,list,order)

      ! Reorder sampling
      sam2_glb = sam2_glb(order)
   end if

   ! Release memory
   deallocate(glb_to_loc)
   deallocate(glb_to_proc)
   deallocate(hash_loc)
   deallocate(hash_glb)
   deallocate(mask_glb)
   if (mpl%main) then
      deallocate(list)
      deallocate(order)
   end if
else
   ! First subsampling (local, using ATLAS octahedral grid)
   nfac = 1
   retry = .true.
   do while (retry)
      ! Update nfac
      nfac = 2*nfac

      if (lfull_grid) then
         ! Number of required points
         ns1_glb = nfac*int(real(ns2_glb,kind_real)*4.0*pi/area)

         ! Octahedral grid
         n = int(-4.5+sqrt(20.25+0.25*real(ns1_glb,kind_real)))+1
         write(gridid,'(a,i5.5)') 'O',n
         agrid = atlas_structuredgrid(gridid)
         ns1_glb = 4*n**2+36*n

         ! Allocation
         allocate(lon1_loc(ns1_glb))
         allocate(lat1_loc(ns1_glb))
         allocate(rh1_loc(ns1_glb))
         allocate(sam1_loc(ns1_glb))
         allocate(dist1_loc(ns1_glb))
         allocate(s1_loc_to_glb(ns1_glb))
         allocate(lmask(n_loc))

         ! Initialization
         lmask = mask_loc
         if (lverbosity) call mpl%prog_init(ns1_glb)

         ! Loop over octahedral grid points
         is1_loc = 0
         is1_glb = 0
         do iy=1,int(agrid%ny(),kind_int)
            do ix=1,int(agrid%nx(iy),kind_int)
               ! Global index
               is1_glb = is1_glb+1

               ! Get longitude/latitude
               lonlat = agrid%lonlat(ix,iy)*deg2rad

               ! Find nearest neighbor in universe
               call tree_uni%find_nearest_neighbors(lonlat(1),lonlat(2),1,nn_index(1:1),nn_dist(1:1))

               if (uni_to_proc(nn_index(1))==mpl%myproc) then
                  ! Get local index
                  i_loc = uni_to_loc(nn_index(1))

                  ! Keep valid points
                  if (lmask(i_loc)) then
                     is1_loc = is1_loc+1
                     lon1_loc(is1_loc) = lon_loc(i_loc)
                     lat1_loc(is1_loc) = lat_loc(i_loc)
                     rh1_loc(is1_loc) = rh_loc(i_loc)
                     sam1_loc(is1_loc) = loc_to_glb(i_loc)
                     dist1_loc(is1_loc) = nn_dist(1)
                     s1_loc_to_glb(is1_loc) = is1_glb
                     lmask(i_loc) = .false.
                  end if
               end if

               ! Update
               if (lverbosity) call mpl%prog_print(is1_glb)
            end do
         end do
         if (lverbosity) call mpl%prog_final(.false.)

         ! Number of local valid points
         ns1_loc = is1_loc

         ! Release memory
         deallocate(lmask)
      else
         ! All point are used
         ns1_glb = n_glb_eff
         ns1_loc = n_loc_eff

         ! Allocation
         allocate(lon1_loc(ns1_loc))
         allocate(lat1_loc(ns1_loc))
         allocate(rh1_loc(ns1_loc))
         allocate(sam1_loc(ns1_loc))

         ! Copy
         i_loc_eff = 0
         do i_loc=1,n_loc
            if (mask_loc(i_loc)) then
               i_loc_eff = i_loc_eff+1
               lon1_loc(i_loc_eff) = lon_loc(i_loc)
               lat1_loc(i_loc_eff) = lat_loc(i_loc)
               rh1_loc(i_loc_eff) = rh_loc(i_loc)
               sam1_loc(i_loc_eff) = loc_to_glb(i_loc)
            end if
         end do
      end if

      if (mpl%main) then
         if (lfull_grid) then
            ! Continue printing
            write(mpl%info,'(a)') ' => '
            if (lverbosity) call mpl%flush(.false.)
         end if

         ! Allocation
         allocate(lon1_glb(ns1_glb))
         allocate(lat1_glb(ns1_glb))
         allocate(rh1_glb(ns1_glb))
         allocate(sam1_glb(ns1_glb))
         if (lfull_grid) allocate(dist1_glb(ns1_glb))

         ! Initialization
         lon1_glb = mpl%msv%valr
         lat1_glb = mpl%msv%valr
         rh1_glb = mpl%msv%valr
         sam1_glb = mpl%msv%vali
         if (lfull_grid) then
            dist1_glb = mpl%msv%valr
         else
            offset = 0
         end if

         ! Receive and copy data on rootproc
         do iproc=1,mpl%nproc
            if (iproc==mpl%rootproc) then
               ! Copy dimension
               ns1_loc_tmp = ns1_loc
            else
               ! Receive dimension
               call mpl%f_comm%receive(ns1_loc_tmp,iproc-1,mpl%tag,status)
            end if

            ! Allocation
            allocate(lon1_loc_tmp(ns1_loc_tmp))
            allocate(lat1_loc_tmp(ns1_loc_tmp))
            allocate(rh1_loc_tmp(ns1_loc_tmp))
            allocate(sam1_loc_tmp(ns1_loc_tmp))
            if (lfull_grid) then
               allocate(dist1_loc_tmp(ns1_loc_tmp))
               allocate(s1_loc_to_glb_tmp(ns1_loc_tmp))
            end if

            if (iproc==mpl%rootproc) then
               ! Copy data
               lon1_loc_tmp = lon1_loc(1:ns1_loc)
               lat1_loc_tmp = lat1_loc(1:ns1_loc)
               rh1_loc_tmp = rh1_loc(1:ns1_loc)
               sam1_loc_tmp = sam1_loc(1:ns1_loc)
               if (lfull_grid) then
                  dist1_loc_tmp = dist1_loc(1:ns1_loc)
                  s1_loc_to_glb_tmp = s1_loc_to_glb(1:ns1_loc)
               end if
            else
               ! Receive data
               call mpl%f_comm%receive(lon1_loc_tmp,iproc-1,mpl%tag+1,status)
               call mpl%f_comm%receive(lat1_loc_tmp,iproc-1,mpl%tag+2,status)
               call mpl%f_comm%receive(rh1_loc_tmp,iproc-1,mpl%tag+3,status)
               call mpl%f_comm%receive(sam1_loc_tmp,iproc-1,mpl%tag+4,status)
               if (lfull_grid) then
                  call mpl%f_comm%receive(dist1_loc_tmp,iproc-1,mpl%tag+5,status)
                  call mpl%f_comm%receive(s1_loc_to_glb_tmp,iproc-1,mpl%tag+6,status)
               end if
            end if

            ! Update global array
            if (lfull_grid) then
               do is1_loc=1,ns1_loc_tmp
                  ! Global index
                  is1_glb = s1_loc_to_glb_tmp(is1_loc)

                  ! Check the global array status
                  update = mpl%msv%is(dist1_glb(is1_glb))
                  if (.not.update) update = (dist1_loc_tmp(is1_loc)<dist1_glb(is1_glb))
                  if (update) then
                     lon1_glb(is1_glb) = lon1_loc_tmp(is1_loc)
                     lat1_glb(is1_glb) = lat1_loc_tmp(is1_loc)
                     rh1_glb(is1_glb) = rh1_loc_tmp(is1_loc)
                     sam1_glb(is1_glb) = sam1_loc_tmp(is1_loc)
                     dist1_glb(is1_glb) = dist1_loc_tmp(is1_loc)
                  end if
               end do
            else
               ! Copy with offset
               lon1_glb(offset+1:offset+ns1_loc_tmp) = lon1_loc_tmp
               lat1_glb(offset+1:offset+ns1_loc_tmp) = lat1_loc_tmp
               rh1_glb(offset+1:offset+ns1_loc_tmp) = rh1_loc_tmp
               sam1_glb(offset+1:offset+ns1_loc_tmp) = sam1_loc_tmp
               offset = offset+ns1_loc_tmp
            end if

            ! Release memory
            deallocate(lon1_loc_tmp)
            deallocate(lat1_loc_tmp)
            deallocate(rh1_loc_tmp)
            deallocate(sam1_loc_tmp)
            if (lfull_grid) then
               deallocate(dist1_loc_tmp)
               deallocate(s1_loc_to_glb_tmp)
            end if
         end do

         ! Number of global valid points
         if (lfull_grid) then
            ns1_glb_eff = count(mpl%msv%isnot(dist1_glb))
         else
            ns1_glb_eff = ns1_glb
         end if
         retry = (ns1_glb_eff<ns2_glb)
      else
         ! Send data to rootproc
         call mpl%f_comm%send(ns1_loc,mpl%rootproc-1,mpl%tag)
         call mpl%f_comm%send(lon1_loc(1:ns1_loc),mpl%rootproc-1,mpl%tag+1)
         call mpl%f_comm%send(lat1_loc(1:ns1_loc),mpl%rootproc-1,mpl%tag+2)
         call mpl%f_comm%send(rh1_loc(1:ns1_loc),mpl%rootproc-1,mpl%tag+3)
         call mpl%f_comm%send(sam1_loc(1:ns1_loc),mpl%rootproc-1,mpl%tag+4)
         if (lfull_grid) then
            call mpl%f_comm%send(dist1_loc(1:ns1_loc),mpl%rootproc-1,mpl%tag+5)
            call mpl%f_comm%send(s1_loc_to_glb(1:ns1_loc),mpl%rootproc-1,mpl%tag+6)
         end if

         ! Stop printing
         write(mpl%info,'(a)') ''
         if (lverbosity) call mpl%flush
      end if
      if (lfull_grid) then
         call mpl%update_tag(7)
      else
         call mpl%update_tag(5)
      end if

      ! Broadcast
      call mpl%f_comm%broadcast(retry,mpl%rootproc-1)

      ! Check situation
      if ((.not.lfull_grid).and.retry) call mpl%abort(subr,'retry activated for limited grid')
      if ((nfac>=8).and.retry) call mpl%abort(subr,'retry activated for nfac>=8')

      ! Release memory
      deallocate(lon1_loc)
      deallocate(lat1_loc)
      deallocate(rh1_loc)
      deallocate(sam1_loc)
      if (lfull_grid) then
         deallocate(dist1_loc)
         deallocate(s1_loc_to_glb)
      end if
      if (mpl%main.and.retry) then
         deallocate(lon1_glb)
         deallocate(lat1_glb)
         deallocate(rh1_glb)
         deallocate(sam1_glb)
         if (lfull_grid) deallocate(dist1_glb)
      end if
   end do

   if (mpl%main) then
      ! Allocation
      allocate(lon1_glb_eff(ns1_glb_eff))
      allocate(lat1_glb_eff(ns1_glb_eff))
      allocate(rh1_glb_eff(ns1_glb_eff))
      allocate(sam1_glb_eff(ns1_glb_eff))
      allocate(list(ns1_glb_eff))
      allocate(order(ns1_glb_eff))

      ! Copy valid data
      if (lfull_grid) then
         is1_glb_eff = 0
         do is1_glb=1,ns1_glb
            if (mpl%msv%isnot(dist1_glb(is1_glb))) then
               is1_glb_eff = is1_glb_eff+1
               lon1_glb_eff(is1_glb_eff) = lon1_glb(is1_glb)
               lat1_glb_eff(is1_glb_eff) = lat1_glb(is1_glb)
               rh1_glb_eff(is1_glb_eff) = rh1_glb(is1_glb)
               sam1_glb_eff(is1_glb_eff) = sam1_glb(is1_glb)
            end if
         end do
      else
         lon1_glb_eff = lon1_glb
         lat1_glb_eff = lat1_glb
         rh1_glb_eff = rh1_glb
         sam1_glb_eff = sam1_glb
      end if

      ! Define points order
      do is1_glb=1,ns1_glb_eff
         list(is1_glb) = lonlathash(lon1_glb_eff(is1_glb),lat1_glb_eff(is1_glb))
      end do
      call qsort(ns1_glb_eff,list,order)

      ! Reorder data
      lon1_glb_eff = lon1_glb_eff(order)
      lat1_glb_eff = lat1_glb_eff(order)
      rh1_glb_eff = rh1_glb_eff(order)
      sam1_glb_eff = sam1_glb_eff(order)

      ! Release memory
      deallocate(lon1_glb)
      deallocate(lat1_glb)
      deallocate(rh1_glb)
      deallocate(sam1_glb)
      if (lfull_grid) deallocate(dist1_glb)
      deallocate(list)
      deallocate(order)

      ! Allocation
      nrep_eff = min(nrep,ns1_glb_eff-ns2_glb)
      allocate(sam2_glb_tmp(ns2_glb+nrep_eff))
      allocate(lmask(ns1_glb_eff))
      allocate(smask(ns1_glb_eff))
      allocate(to_valid(ns1_glb_eff))

      ! Initialization
      sam2_glb_tmp = mpl%msv%vali
      lmask = .true.
      smask = .false.
      to_valid = mpl%msv%vali
      do is1_glb=1,ns1_glb_eff
         to_valid(is1_glb) = is1_glb
      end do
      ns1_glb_val = ns1_glb_eff
      if (lverbosity) call mpl%prog_init(ns2_glb+nrep_eff)

      if (lfast) then
         ! Define sampling with a cumulative distribution function

         ! Allocation
         allocate(cdf(ns1_glb_eff))

         ! Initialization
         cdf(1) = 0.0
         do is1_glb=2,ns1_glb_eff
            if (lmask(is1_glb)) then
               cdf(is1_glb) = cdf(is1_glb-1)+1.0/rh1_glb_eff(is1_glb)**2
            end if
         end do
         cdf_norm = 1.0/cdf(ns1_glb_eff)
         cdf = cdf*cdf_norm

         do is2_glb=1,ns2_glb+nrep_eff
            ! Generate random number
            call rng%rand_real(0.0_kind_real,1.0_kind_real,rr)

            ! Dichotomy to find the value
            irvalmin = 1
            irvalmax = ns1_glb_val
            do while (irvalmax-irvalmin>1)
               irval = (irvalmin+irvalmax)/2
               if ((sup(cdf(irvalmin),rr).and.sup(cdf(irval),rr)).or.(inf(cdf(irvalmin),rr).and.inf(cdf(irval),rr))) then
                  irvalmin = irval
               else
                  irvalmax = irval
               end if
            end do

            ! New sampling point
            ir = to_valid(irval)
            sam2_glb_tmp(is2_glb) = ir

            ! Shift valid points array
            if (irval<ns1_glb_val) then
               cdf(irval:ns1_glb_val-1) = cdf(irval+1:ns1_glb_val)
               to_valid(irval:ns1_glb_val-1) = to_valid(irval+1:ns1_glb_val)
            end if
            ns1_glb_val = ns1_glb_val-1

            ! Renormalize cdf
            cdf_norm = 1.0/cdf(ns1_glb_val)
            cdf(1:ns1_glb_val) = cdf(1:ns1_glb_val)*cdf_norm

            ! Update
            if (lverbosity) call mpl%prog_print(is2_glb)
         end do
         if (lverbosity) call mpl%prog_final(.false.)

         ! Release memory
         deallocate(cdf)
      else
         ! Define sampling with a KD-tree
         do is2_glb=1,ns2_glb+nrep_eff
            if (is2_glb>2) then
               ! Allocation
               call tree%alloc(mpl,ns1_glb_eff,mask=smask)

               ! Initialization
               call tree%init(lon1_glb_eff,lat1_glb_eff)
            end if

            ! Initialization
            distmax = 0.0
            irmax = 0
            irvalmax = 0
            itry = 1

            ! Find a new point
            do itry=1,ntry
               ! Generate a random index among valid points
               call rng%rand_integer(1,ns1_glb_val,irval)
               ir = to_valid(irval)

               ! Check point validity
               if (is2_glb==1) then
                  ! Accept point
                  irvalmax = irval
                  irmax = ir
               else
                  if (is2_glb==2) then
                     ! Compute distance
                     nn_index(1) = sam2_glb_tmp(1)
                     call sphere_dist(lon1_glb_eff(ir),lat1_glb_eff(ir),lon1_glb_eff(nn_index(1)),lat1_glb_eff(nn_index(1)), &
                   & nn_dist(1))
                  else
                     ! Find nearest neighbor distance
                     call tree%find_nearest_neighbors(lon1_glb_eff(ir),lat1_glb_eff(ir),1,nn_index(1:1),nn_dist(1:1))
                  end if
                  d = nn_dist(1)**2/(rh1_glb_eff(ir)**2+rh1_glb_eff(nn_index(1))**2)

                  ! Check distance
                  if (sup(d,distmax)) then
                     distmax = d
                     irvalmax = irval
                     irmax = ir
                  end if
               end if
            end do

            ! Delete tree
            if (is2_glb>2) call tree%dealloc

            ! Add point to sampling
            if (irmax>0) then
               ! New sampling point
               sam2_glb_tmp(is2_glb) = irmax
               lmask(irmax) = .false.
               smask(irmax) = .true.

               ! Shift valid points array
               if (irvalmax<ns1_glb_val) to_valid(irvalmax:ns1_glb_val-1) = to_valid(irvalmax+1:ns1_glb_val)
               ns1_glb_val = ns1_glb_val-1
            end if

            ! Update
            if (lverbosity) call mpl%prog_print(is2_glb)
         end do
         if (lverbosity) call mpl%prog_final(.false.)
      end if

      if (nrep_eff>0) then
         ! Continue printing
         write(mpl%info,'(a)') ' => '
         if (lverbosity) call mpl%flush(.false.)

         ! Allocation
         allocate(rmask(ns2_glb+nrep_eff))
         allocate(lon_rep(ns2_glb+nrep_eff))
         allocate(lat_rep(ns2_glb+nrep_eff))
         allocate(dist(ns2_glb+nrep_eff))

         ! Initialization
         rmask = .true.
         do is2_glb=1,ns2_glb+nrep_eff
            lon_rep(is2_glb) = lon1_glb_eff(sam2_glb_tmp(is2_glb))
            lat_rep(is2_glb) = lat1_glb_eff(sam2_glb_tmp(is2_glb))
         end do
         dist = mpl%msv%valr
         if (lverbosity) call mpl%prog_init(nrep_eff)

         ! Remove closest points
         do irep=1,nrep_eff
            ! Allocation
            call tree%alloc(mpl,ns2_glb+nrep_eff,mask=rmask)

            ! Initialization
            call tree%init(lon_rep,lat_rep)

            ! Get minimum distance
            do is2_glb=1,ns2_glb+nrep_eff
               if (rmask(is2_glb)) then
                  ! Find nearest neighbor distance
                  call tree%find_nearest_neighbors(lon1_glb_eff(sam2_glb_tmp(is2_glb)),lat1_glb_eff(sam2_glb_tmp(is2_glb)), &
                & 2,nn_index,nn_dist)
                  if (nn_index(1)==is2_glb) then
                     dist(is2_glb) = nn_dist(2)
                  elseif (nn_index(2)==is2_glb) then
                     dist(is2_glb) = nn_dist(1)
                  else
                     call mpl%abort(subr,'wrong index in replacement')
                  end if
                  dist(is2_glb) = dist(is2_glb)**2/(rh1_glb_eff(sam2_glb_tmp(nn_index(1)))**2 &
                                & +rh1_glb_eff(sam2_glb_tmp(nn_index(2)))**2)
               end if
            end do

            ! Delete tree
            call tree%dealloc

            ! Remove worst point
            distmin = huge(1.0)
            is2_glb_min = mpl%msv%vali
            do is2_glb=1,ns2_glb+nrep_eff
               if (rmask(is2_glb)) then
                  if (inf(dist(is2_glb),distmin)) then
                     is2_glb_min = is2_glb
                     distmin = dist(is2_glb)
                  end if
               end if
            end do
            rmask(is2_glb_min) = .false.

             ! Update
            if (lverbosity) call mpl%prog_print(irep)
         end do
         if (lverbosity) call mpl%prog_final

         ! Copy sam2_glb
         js = 0
         do is2_glb=1,ns2_glb+nrep_eff
            if (rmask(is2_glb)) then
               js = js+1
               sam2_glb(js) = sam2_glb_tmp(is2_glb)
            end if
         end do

         ! Release memory
         deallocate(rmask)
         deallocate(lon_rep)
         deallocate(lat_rep)
         deallocate(dist)
      else
         ! Stop printing
         write(mpl%info,'(a)') ''
         if (lverbosity) call mpl%flush

         ! Copy sam2_glb
         sam2_glb = sam2_glb_tmp
      end if

      ! Apply first sampling step
      sam2_glb = sam1_glb_eff(sam2_glb)

      ! Release memory
      deallocate(lon1_glb_eff)
      deallocate(lat1_glb_eff)
      deallocate(rh1_glb_eff)
      deallocate(sam1_glb_eff)
      deallocate(sam2_glb_tmp)
      deallocate(lmask)
      deallocate(smask)
      deallocate(to_valid)
   end if
end if

! Broadcast
call mpl%f_comm%broadcast(sam2_glb,mpl%rootproc-1)

end subroutine initialize_sampling

end module tools_samp
