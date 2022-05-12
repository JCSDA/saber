! (C) Copyright 2022 United States Government as represented by the Administrator of the National
!     Aeronautics and Space Administration
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module gsi_covariance_mod

! atlas
use atlas_module,                   only: atlas_fieldset, atlas_field

! fckit
use fckit_mpi_module,               only: fckit_mpi_comm
use fckit_configuration_module,     only: fckit_configuration

! oops
use kinds,                          only: kind_real
use random_mod

! saber
use gsi_grid_mod,                   only: gsi_grid

! gsibclim
use m_gsibclim,                     only: gsibclim_init
use m_gsibclim,                     only: gsibclim_cv_space
use m_gsibclim,                     only: gsibclim_sv_space
use m_gsibclim,                     only: gsibclim_befname
use m_gsibclim,                     only: gsibclim_final
use gsi_bundlemod,                  only: gsi_bundle
use gsi_bundlemod,                  only: gsi_bundlegetpointer
use control_vectors,                only: control_vector
use control_vectors,                only: cvars2d,cvars3d
use control_vectors,                only: allocate_cv
use control_vectors,                only: deallocate_cv
use control_vectors,                only: assignment(=)

use state_vectors,                  only: allocate_state
use state_vectors,                  only: svars2d,svars3d
use state_vectors,                  only: deallocate_state

implicit none
private
public gsi_covariance


! Fortran class header
type :: gsi_covariance
  type(gsi_grid) :: grid
  logical :: noGSI
  logical :: bypassGSIbe
  logical :: cv   ! cv=.true.; sv=.false.
  integer :: mp_comm_world
  integer :: rank
  integer :: lat2,lon2 ! these belog to gsi_grid
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: randomize
    procedure, public :: multiply
    procedure, public :: multiply_ad
end type gsi_covariance

character(len=*), parameter :: myname='gsi_covariance_mod'

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, comm, config)

! Arguments
class(gsi_covariance),     intent(inout) :: self
type(fckit_mpi_comm),      intent(in)    :: comm
type(fckit_configuration), intent(in)    :: config
!type(atlas_fieldset),      intent(in)    :: background   ! Uncomment once background available as Atlas fieldset
!type(atlas_fieldset),      intent(in)    :: first_guess

! Locals
character(len=*), parameter :: myname_=myname//'*create'
character(len=:), allocatable :: nml,bef
logical :: central
integer :: layout(2)

! Hold communicator
! -----------------
!self%mp_comm_world=comm%communicator()

! Create the grid
! ---------------
call self%grid%create(config, comm)
self%rank = comm%rank()

call config%get_or_die("debugging bypass gsi", self%noGSI)
if (.not. self%noGSI) then
  call config%get_or_die("saber central block", central)
  if (.not. central) then
     call abor1_ftn(myname_//": not ready to handle sqrt(B) case")
  endif
  call config%get_or_die("debugging deep bypass gsi B error", self%bypassGSIbe)

! Get required name of resources for GSI B error
! ----------------------------------------------
  call config%get_or_die("gsi berror namelist file",  nml)
  call config%get_or_die("gsi error covariance file", bef)

! Initialize GSI-Berror components
! --------------------------------
  layout=self%grid%layout
! layout=-1
  call gsibclim_init(self%cv,self%lat2,self%lon2,nmlfile=nml,befile=bef,layout=layout,comm=comm%communicator())
endif

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(gsi_covariance) :: self

! Locals

if (.not. self%noGSI) then
   call gsibclim_final(.false.)
endif

! Delete the grid
! ---------------
call self%grid%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine randomize(self, fields)

! Arguments
class(gsi_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

! Locals
type(atlas_field) :: afield
real(kind=kind_real), pointer :: psi(:,:), chi(:,:), t(:,:), q(:,:), qi(:,:), ql(:,:), o3(:,:)
real(kind=kind_real), pointer :: ps(:)

integer, parameter :: rseed = 3

! Get Atlas field
afield = fields%field('stream_function')
call afield%data(psi)

afield = fields%field('velocity_potential')
call afield%data(chi)

afield = fields%field('air_temperature')
call afield%data(t)

afield = fields%field('surface_pressure')
call afield%data(ps)

afield = fields%field('specific_humidity')
call afield%data(q)

afield = fields%field('cloud_liquid_ice')
call afield%data(qi)

afield = fields%field('cloud_liquid_water')
call afield%data(ql)

afield = fields%field('ozone_mass_mixing_ratio')
call afield%data(o3)


! Set fields to random numbers
call normal_distribution(psi, 0.0_kind_real, 1.0_kind_real, rseed)


end subroutine randomize

! --------------------------------------------------------------------------------------------------

subroutine multiply(self, fields)

! Arguments
class(gsi_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

! Locals
type(atlas_field) :: afield
real(kind=kind_real), pointer :: rank1(:)       =>NULL()
real(kind=kind_real), pointer :: rank2(:,:)     =>NULL()
real(kind=kind_real), pointer :: gsivar3d(:,:,:)=>NULL()
real(kind=kind_real), pointer :: gsivar2d(:,:)  =>NULL()
real(kind=kind_real), allocatable :: aux(:,:)
real(kind=kind_real), allocatable :: aux1(:)

type(control_vector) :: gsicv
type(gsi_bundle) :: gsisv(1)
integer :: isc,iec,jsc,jec,npz
integer :: iv,k,ier
integer,parameter :: hw=1
logical nosgi
integer, parameter :: rseed = 3

character(len=20), allocatable :: gvars2d(:),gvars3d(:)

! afield = fields%field('surface_pressure')
! call afield%data(rank1)
! rank1 = 0.0_kind_real
! rank1(int(size(rank1)/2)) = 1.0_kind_real
! return
if (self%noGSI) return

!   gsi-surface: k=1
!   quench-surface: k=1
!   fv3-surface: k=km
isc=self%grid%isc 
iec=self%grid%iec 
jsc=self%grid%jsc
jec=self%grid%jec
npz=self%grid%npz

! Allocate control vector as defined in GSI
! -----------------------------------------
if (self%cv) then
   allocate(gvars2d(size(cvars2d)),gvars3d(size(cvars3d)))
   gvars2d=cvars2d; gvars3d=cvars3d
   call allocate_cv(gsicv)
else
   allocate(gvars2d(size(svars2d)),gvars3d(size(svars3d)))
   gvars2d=svars2d; gvars3d=svars3d
   call allocate_state(gsisv(1))
endif

! Convert GSI bundle fields to Atlas (fieldsets)
! ----------------------------------------------
do iv=1,size(gvars2d)
   if (self%cv) then
      call gsi_bundlegetpointer (gsicv%step(1),gvars2d(iv),gsivar2d,ier)
   else
      call gsi_bundlegetpointer (     gsisv(1),gvars2d(iv),gsivar2d,ier)
   endif
   if (ier/=0) cycle
   call get_rank1_(ier)
   if (ier/=0) cycle
   allocate(aux(size(gsivar2d,1),size(gsivar2d,2)))
   call addhalo_(rank1,aux)
   gsivar2d=aux
   deallocate(aux)
enddo
do iv=1,size(gvars3d)
   if (self%cv) then
      call gsi_bundlegetpointer (gsicv%step(1),gvars3d(iv),gsivar3d,ier)
   else
      call gsi_bundlegetpointer (     gsisv(1),gvars3d(iv),gsivar3d,ier)
   endif
   if (ier/=0) cycle
   call get_rank2_(ier)
   if (ier/=0) cycle
 allocate(aux(size(gsivar3d,1),size(gsivar3d,2)))
   if (self%grid%vflip) then
      do k=1,npz
         call addhalo_(rank2(k,:),aux)
         gsivar3d(:,:,npz-k+1)=aux
      enddo
   else
      do k=1,npz
         call addhalo_(rank2(k,:),aux)
         gsivar3d(:,:,k)=aux
      enddo
   endif
 deallocate(aux)
enddo

! Apply GSI B-error operator
! --------------------------
if (self%cv) then
   call gsibclim_cv_space(gsicv,internalcv=.false.,bypassbe=self%bypassGSIbe)
else
   call gsibclim_sv_space(gsisv,internalsv=.false.,bypassbe=self%bypassGSIbe)
endif

! Convert back to Atlas Fields
! ----------------------------
do iv=1,size(gvars2d)
   if (self%cv) then
      call gsi_bundlegetpointer (gsicv%step(1),gvars2d(iv),gsivar2d,ier)
   else
      call gsi_bundlegetpointer (gsisv(1),gvars2d(iv),gsivar2d,ier)
   endif
   if (ier/=0) cycle
   call get_rank1_(ier)
   if (ier/=0) cycle
   allocate(aux1(size(rank1)))
   call remhalo_(gsivar2d,aux1)
   rank1=aux1
   deallocate(aux1)
enddo
do iv=1,size(gvars3d)
   if (self%cv) then
      call gsi_bundlegetpointer (gsicv%step(1),gvars3d(iv),gsivar3d,ier)
   else
      call gsi_bundlegetpointer (gsisv(1),gvars3d(iv),gsivar3d,ier)
   endif
   if (ier/=0) cycle
   call get_rank2_(ier)
   if (ier/=0) cycle
   allocate(aux1(size(rank2,2)))
   if (self%grid%vflip) then
      do k=1,npz
         call remhalo_(gsivar3d(:,:,k),aux1)
         rank2(npz-k+1,:)=aux1
      enddo
   else
      do k=1,npz
         call remhalo_(gsivar3d(:,:,k),aux1)
         rank2(k,:)=aux1
      enddo
   endif
 deallocate(aux1)
enddo

! Release pointer
! ---------------
if (self%cv) then
   call deallocate_cv(gsicv)
else
   call deallocate_state(gsisv(1))
endif
deallocate(gvars2d,gvars3d)
call afield%final()

contains
   subroutine get_rank1_(ier)
   integer ier
   ier=-1
   if (trim(gvars2d(iv)) == 'ps') then
      if (.not.fields%has_field('surface_pressure')) return
      afield = fields%field('surface_pressure')
      call afield%data(rank1)
      ier=0
   endif
   if (trim(gvars2d(iv)) == 'sst') then
      if (.not.fields%has_field('skin_surface_temperature')) return
      afield = fields%field('skin_surface_temperature')
      call afield%data(rank1)
      ier=0
   endif
   end subroutine get_rank1_

   subroutine get_rank2_(ier)
   integer ier
   ier=-1
   if (trim(gvars3d(iv)) == 'u') then
      if (.not.fields%has_field('eastward_wind')) return
      afield = fields%field('eastward_wind')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 'v') then
      if (.not.fields%has_field('northward_wind')) return
      afield = fields%field('northward_wind')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 'sf') then
      if (.not.fields%has_field('stream_function')) return
      afield = fields%field('stream_function')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 'vp') then
      if (.not.fields%has_field('velocity_potential')) return
      afield = fields%field('velocity_potential')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 't' .or. trim(gvars3d(iv)) == 'tv') then  ! this needs attention
      if (.not.fields%has_field('air_temperature')) return
      afield = fields%field('air_temperature')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 'q') then ! rh
      if (.not.fields%has_field('specific_humidity')) return
      afield = fields%field('specific_humidity')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 'qi') then
      if (.not.fields%has_field('cloud_liquid_ice')) return
      afield = fields%field('cloud_liquid_ice')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 'ql') then
      if (.not.fields%has_field('cloud_liquid_water')) return
      afield = fields%field('cloud_liquid_water')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 'qr') then
      if (.not.fields%has_field('cloud_liquid_rain')) return
      afield = fields%field('cloud_liquid_rain')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 'qs') then
      if (.not.fields%has_field('cloud_liquid_snow')) return
      afield = fields%field('cloud_liquid_snow')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 'cw') then
      if (.not.fields%has_field('cloud_water')) return
      afield = fields%field('cloud_water')
      call afield%data(rank2)
      ier=0
   endif
   if (trim(gvars3d(iv)) == 'oz') then
      if (.not.fields%has_field('ozone_mass_mixing_ratio')) return
      afield = fields%field('ozone_mass_mixing_ratio')
      call afield%data(rank2)
      ier=0
   endif
   end subroutine get_rank2_

   subroutine addhalo_(rank,var)
   real(kind=kind_real),intent(in) :: rank(:)
   real(kind=kind_real),intent(inout):: var(:,:)
   integer ii,jj,iii,jjj,jnode
   jnode=1
   do jj=2,self%lat2-1
      do ii=2,self%lon2-1
         var(jj,ii) = rank(jnode)
         jnode = jnode + 1
      enddo
   enddo
   end subroutine addhalo_

   subroutine remhalo_(var,rank)
   real(kind=kind_real),intent(in) :: var(:,:)
   real(kind=kind_real),intent(out):: rank(:)
   integer ii,jj,iii,jjj,jnode
   jnode=1
   do jj=2,self%lat2-1
      do ii=2,self%lon2-1
         rank(jnode) = var(jj,ii)
         jnode = jnode + 1
      enddo
   enddo
   end subroutine remhalo_


end subroutine multiply

! --------------------------------------------------------------------------------------------------

subroutine multiply_ad(self, fields)

! Arguments
class(gsi_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

! This routine only needed when B = G^T G (sqrt-factored)

! To do list for this method
! 1. Convert fields (Atlas fieldsets) to GSI bundle
! 2. Call GSI covariance operator adjoint (sqrt version)
!        afield = fields%field('stream_function')
!        call afield%data(var3d)
!        var3d=0.0_kind_real

end subroutine multiply_ad

! --------------------------------------------------------------------------------------------------

end module gsi_covariance_mod
