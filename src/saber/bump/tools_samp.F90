!----------------------------------------------------------------------
! Module: tools_samp
! Purpose: sampling functions
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_samp

use fckit_mpi_module, only: fckit_mpi_sum,fckit_mpi_status
use tools_const, only: pi
use tools_func, only: lonlatmod,sphere_dist
use tools_kinds, only: kind_real
use tools_qsort, only: qsort
use tools_repro, only: repro,inf,sup
use type_mpl, only: mpl_type
use type_rng, only: rng_type
use type_tree, only: tree_type

implicit none

logical,parameter :: nn_stats = .false. ! Compute and print subsampling statistics

private
public :: initialize_sampling

contains

!----------------------------------------------------------------------
! Subroutine: initialize_sampling
! Purpose: intialize sampling
!----------------------------------------------------------------------
subroutine initialize_sampling(mpl,rng,n_loc,lon_loc,lat_loc,mask_loc,rh_loc,loc_to_glb,ntry,nrep,ns2_glb,sam2_glb,fast,verbosity)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl          ! MPI data
type(rng_type),intent(inout) :: rng          ! Random number generator
integer,intent(in) :: n_loc                  ! Number of points (local)
real(kind_real),intent(in) :: lon_loc(n_loc) ! Longitudes (local)
real(kind_real),intent(in) :: lat_loc(n_loc) ! Latitudes (local)
logical,intent(in) :: mask_loc(n_loc)        ! Mask (local)
real(kind_real),intent(in) :: rh_loc(n_loc)  ! Horizontal support radius (local)
integer,intent(in) :: loc_to_glb(n_loc)      ! Local to global index
integer,intent(in) :: ntry                   ! Number of tries
integer,intent(in) :: nrep                   ! Number of replacements
integer,intent(in) :: ns2_glb                ! Number of samplings points (global)
integer,intent(out) :: sam2_glb(ns2_glb)     ! Horizontal sampling index (global)
logical,intent(in),optional :: fast          ! Fast sampling flag
logical,intent(in),optional :: verbosity     ! Verbosity flag

! Local variables
integer :: n_glb,n_glb_eff,i_glb,is2_glb,iproc,js,irep,irmax,itry,is1_glb,ir
integer :: irval,irvalmin,irvalmax,is2_glb_min,nrep_eff,nn_index(2),ns1_glb,is1_glb_eff,ns1_glb_eff,ns1_glb_val
integer,allocatable :: glb_to_loc(:),glb_to_proc(:),sam1_glb(:),sam1_glb_tmp(:),to_valid(:),sam2_glb_tmp(:),order(:)
real(kind_real) :: d,distmax,distmin,nn_dist(2),cdf_norm,rr
real(kind_real),allocatable :: lon1_glb(:),lat1_glb(:),dist1_glb(:),rh1_glb(:)
real(kind_real),allocatable :: lon1_glb_tmp(:),lat1_glb_tmp(:),dist1_glb_tmp(:),rh1_glb_tmp(:)
real(kind_real),allocatable :: list(:),cdf(:)
real(kind_real),allocatable :: lon_rep(:),lat_rep(:),dist(:)
logical :: lfast,lverbosity
logical,allocatable :: mask_glb(:),lmask(:),smask(:),rmask(:)
character(len=1024),parameter :: subr = 'initialize_sampling'
type(fckit_mpi_status) :: status
type(tree_type) :: tree

! Local verbosity flag
lverbosity = .true.
if (present(verbosity)) lverbosity = verbosity

! Number of effective points
call mpl%f_comm%allreduce(count(mask_loc),n_glb_eff,fckit_mpi_sum())

! Check mask size
if (n_glb_eff==0) then
   call mpl%abort(subr,'empty mask in initialize sampling')
elseif (n_glb_eff<ns2_glb) then
   call mpl%abort(subr,'ns2_glb greater than n_glb_eff in initialize_sampling')
elseif (n_glb_eff==ns2_glb) then
   write(mpl%info,'(a)') ' all points are used'
   if (lverbosity) call mpl%flush

   ! Global size
   call mpl%f_comm%allreduce(n_loc,n_glb,fckit_mpi_sum())

   ! Allocation
   allocate(glb_to_loc(n_glb))
   allocate(glb_to_proc(n_glb))
   allocate(mask_glb(n_glb))

   ! Communication
   call mpl%glb_to_loc_index(n_loc,loc_to_glb,n_glb,glb_to_loc,glb_to_proc)
   call mpl%loc_to_glb(n_loc,mask_loc,n_glb,glb_to_proc,glb_to_loc,.false.,mask_glb)

   if (mpl%main) then
      ! Use all valid points
      is2_glb = 0
      do i_glb=1,n_glb
         if (mask_glb(i_glb)) then
            is2_glb = is2_glb+1
            sam2_glb(is2_glb) = i_glb
         end if
      end do
   end if

   ! Release memory
   deallocate(glb_to_loc)
   deallocate(glb_to_proc)
   deallocate(mask_glb)
else
   ! Define global size
   ns1_glb = min(5*ns2_glb,n_glb_eff)

   ! Allocation
   call tree%alloc(mpl,n_loc)

   ! Initialization
   call tree%init(lon_loc,lat_loc)

   ! Allocation
   allocate(lon1_glb(ns1_glb))
   allocate(lat1_glb(ns1_glb))
   allocate(dist1_glb(ns1_glb))
   allocate(rh1_glb(ns1_glb))
   allocate(sam1_glb(ns1_glb))

   ! Synchronize seeds
   call rng%sync(mpl)

   ! Initialization
   if (lverbosity) call mpl%prog_init(ns1_glb)

   do is1_glb=1,ns1_glb
      ! Generate random coordinates
      call rng%rand_real(0.0_kind_real,1.0_kind_real,rr)
      lon1_glb(is1_glb) = 2.0*pi*rr
      call rng%rand_real(-1.0_kind_real,1.0_kind_real,rr)
      lat1_glb(is1_glb) = asin(rr)
      call lonlatmod(lon1_glb(is1_glb),lat1_glb(is1_glb))

      ! Find nearest neighbor
      call tree%find_nearest_neighbors(lon1_glb(is1_glb),lat1_glb(is1_glb),1,nn_index(1:1),dist1_glb(is1_glb:is1_glb))

      ! Copy rh and global index if the nearest neighbor is valid
      if (mask_loc(nn_index(1))) then
         rh1_glb(is1_glb) = rh_loc(nn_index(1))
         sam1_glb(is1_glb) = loc_to_glb(is1_glb)
      else
         rh1_glb(is1_glb) = mpl%msv%valr
         sam1_glb(is1_glb) = mpl%msv%vali
      end if

      ! Update
      if (lverbosity) call mpl%prog_print(is1_glb)
   end do
   if (lverbosity) call mpl%prog_final(.false.)

   ! Delete tree
   call tree%dealloc

   ! Desynchronize seeds
   call rng%desync(mpl)

   if (mpl%main) then
      ! Continue printing
      write(mpl%info,'(a)') ' => '
      if (lverbosity) call mpl%flush(.false.)

      ! Allocation
      allocate(lon1_glb_tmp(ns1_glb))
      allocate(lat1_glb_tmp(ns1_glb))
      allocate(dist1_glb_tmp(ns1_glb))
      allocate(rh1_glb_tmp(ns1_glb))
      allocate(sam1_glb_tmp(ns1_glb))

      ! Initialization
      lon1_glb_tmp = lon1_glb
      lat1_glb_tmp = lat1_glb
      dist1_glb_tmp = dist1_glb
      rh1_glb_tmp = rh1_glb
      sam1_glb_tmp = sam1_glb

      ! Receive and copy data on rootproc
      do iproc=1,mpl%nproc
         if (iproc/=mpl%rootproc) then
            ! Receive data
            call mpl%f_comm%receive(lon1_glb,iproc-1,mpl%tag,status)
            call mpl%f_comm%receive(lat1_glb,iproc-1,mpl%tag+1,status)
            call mpl%f_comm%receive(dist1_glb,iproc-1,mpl%tag+2,status)
            call mpl%f_comm%receive(rh1_glb,iproc-1,mpl%tag+3,status)
            call mpl%f_comm%receive(sam1_glb,iproc-1,mpl%tag+4,status)

            ! Copy data
            do is1_glb=1,ns1_glb
               ! Check if the nearest neighbor is vali
               if (mpl%msv%isnot(rh1_glb(is1_glb)).and.mpl%msv%isnot(sam1_glb(is1_glb))) then
                  ! Check if the distance is smaller
                  if (inf(dist1_glb_tmp(is1_glb),dist1_glb(is1_glb))) then
                     lon1_glb_tmp(is1_glb) = lon1_glb(is1_glb)
                     lat1_glb_tmp(is1_glb) = lat1_glb(is1_glb)
                     dist1_glb_tmp(is1_glb) = dist1_glb(is1_glb)
                     rh1_glb_tmp(is1_glb) = rh1_glb(is1_glb)
                     sam1_glb_tmp(is1_glb) = sam1_glb(is1_glb)
                  end if
               end if
            end do
         end if
      end do

      ! Remove unvalid points
      lon1_glb = mpl%msv%valr
      lat1_glb = mpl%msv%valr
      rh1_glb = mpl%msv%valr
      sam1_glb = mpl%msv%vali
      is1_glb_eff = 0
      do is1_glb=1,ns1_glb
         if (mpl%msv%isnot(rh1_glb_tmp(is1_glb)).and.mpl%msv%isnot(sam1_glb_tmp(is1_glb))) then
            is1_glb_eff = is1_glb_eff+1
            lon1_glb(is1_glb_eff) = lon1_glb_tmp(is1_glb)
            lat1_glb(is1_glb_eff) = lat1_glb_tmp(is1_glb)
            rh1_glb(is1_glb_eff) = rh1_glb_tmp(is1_glb)
            sam1_glb(is1_glb_eff) = sam1_glb_tmp(is1_glb)
         end if
      end do

      ! Allocation
      ns1_glb_eff = is1_glb_eff
      allocate(list(ns1_glb_eff))
      allocate(order(ns1_glb_eff))

      ! Define points order
      do is1_glb=1,ns1_glb_eff
         list(is1_glb) = aint(abs(lon1_glb(is1_glb)+pi)*1.0e6)+abs(lat1_glb(is1_glb)+0.5*pi)*1.0e-1
      end do
      call qsort(ns1_glb_eff,list,order)

      ! Reorder data
      lon1_glb(1:ns1_glb_eff) = lon1_glb(order)
      lat1_glb(1:ns1_glb_eff) = lat1_glb(order)
      rh1_glb(1:ns1_glb_eff) = rh1_glb(order)
      sam1_glb(1:ns1_glb_eff) = sam1_glb(order)

      ! Release memory
      deallocate(dist1_glb)
      deallocate(lon1_glb_tmp)
      deallocate(lat1_glb_tmp)
      deallocate(dist1_glb_tmp)
      deallocate(rh1_glb_tmp)
      deallocate(sam1_glb_tmp)
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
      lfast = .false.
      if (present(fast)) lfast = fast

      if (lfast) then
         ! Define sampling with a cumulative distribution function

         ! Allocation
         allocate(cdf(ns1_glb_eff))

         ! Initialization
         cdf(1) = 0.0
         do is1_glb=2,ns1_glb_eff
            if (lmask(is1_glb)) then
               cdf(is1_glb) = cdf(is1_glb-1)+1.0/rh1_glb(is1_glb)**2
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
               call tree%init(lon1_glb,lat1_glb)
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
                     call sphere_dist(lon1_glb(ir),lat1_glb(ir),lon1_glb(nn_index(1)),lat1_glb(nn_index(1)),nn_dist(1))
                  else
                     ! Find nearest neighbor distance
                     call tree%find_nearest_neighbors(lon1_glb(ir),lat1_glb(ir),1,nn_index(1:1),nn_dist(1:1))
                  end if
                  d = nn_dist(1)**2/(rh1_glb(ir)**2+rh1_glb(nn_index(1))**2)

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
            lon_rep(is2_glb) = lon1_glb(sam2_glb_tmp(is2_glb))
            lat_rep(is2_glb) = lat1_glb(sam2_glb_tmp(is2_glb))
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
                  call tree%find_nearest_neighbors(lon1_glb(sam2_glb_tmp(is2_glb)),lat1_glb(sam2_glb_tmp(is2_glb)), &
                & 2,nn_index,nn_dist)
                  if (nn_index(1)==is2_glb) then
                     dist(is2_glb) = nn_dist(2)
                  elseif (nn_index(2)==is2_glb) then
                     dist(is2_glb) = nn_dist(1)
                  else
                     call mpl%abort(subr,'wrong index in replacement')
                  end if
                  dist(is2_glb) = dist(is2_glb)**2/(rh1_glb(sam2_glb_tmp(nn_index(1)))**2+rh1_glb(sam2_glb_tmp(nn_index(2)))**2)
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
      sam2_glb = sam1_glb(sam2_glb)

      ! Release memory
      deallocate(lon1_glb)
      deallocate(lat1_glb)
      deallocate(rh1_glb)
      deallocate(sam1_glb)
      deallocate(sam2_glb_tmp)
      deallocate(lmask)
      deallocate(smask)
      deallocate(to_valid)
   else
      ! Send data to rootproc
      call mpl%f_comm%send(lon1_glb,mpl%rootproc-1,mpl%tag)
      call mpl%f_comm%send(lat1_glb,mpl%rootproc-1,mpl%tag+1)
      call mpl%f_comm%send(dist1_glb,mpl%rootproc-1,mpl%tag+2)
      call mpl%f_comm%send(rh1_glb,mpl%rootproc-1,mpl%tag+3)
      call mpl%f_comm%send(sam1_glb,mpl%rootproc-1,mpl%tag+4)

      ! Release memory
      deallocate(lon1_glb)
      deallocate(lat1_glb)
      deallocate(dist1_glb)
      deallocate(rh1_glb)
      deallocate(sam1_glb)

      ! Stop printing
      write(mpl%info,'(a)') ''
      if (lverbosity) call mpl%flush
   end if
   call mpl%update_tag(5)
end if

! Broadcast
call mpl%f_comm%broadcast(sam2_glb,mpl%rootproc-1)

end subroutine initialize_sampling

end module tools_samp
