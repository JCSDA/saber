!----------------------------------------------------------------------
! Module: tools_atlas
! Purpose: random numbers generator derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------

module tools_atlas

use atlas_module, only: atlas_field,atlas_fieldset,atlas_integer,atlas_real,atlas_functionspace,atlas_functionspace_nodecolumns, &
 & atlas_functionspace_pointcloud,atlas_functionspace_structuredcolumns
use tools_const, only: rad2deg
use tools_kinds, only: kind_int,kind_real
use type_mpl, only: mpl_type

implicit none

private
public :: field_to_fld,fld_to_field,create_atlas_function_space,create_atlas_fieldset,atlas_to_fld,fld_to_atlas

interface field_to_fld
  module procedure field_to_fld_real
  module procedure field_to_fld_logical
end interface

interface fld_to_field
  module procedure fld_to_field_real
end interface

contains

!----------------------------------------------------------------------
! Subroutine: field_to_fld_real
! Purpose: convert ATLAS field to field, real
!----------------------------------------------------------------------
subroutine field_to_fld_real(mpl,afield,fld,lev2d)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(atlas_field),intent(in) :: afield        ! ATLAS field
real(kind_real),intent(out) :: fld(:,:)       ! Field
character(len=*),intent(in),optional :: lev2d ! Level for 2D variables

! Local variables
integer :: nmga,nl0
real(kind_real),pointer :: ptr_1(:),ptr_2(:,:)
character(len=1024) :: llev2d
character(len=1024),parameter :: subr = 'field_to_fld_real'
type(atlas_functionspace) :: afunctionspace
type(atlas_functionspace_nodecolumns) :: afunctionspace_nc
type(atlas_functionspace_pointcloud) :: afunctionspace_pc
type(atlas_functionspace_structuredcolumns) :: afunctionspace_sc

! Local lev2d
llev2d = 'first'
if (present(lev2d)) llev2d = lev2d

! Check kind
if (afield%kind()/=atlas_real(kind_real)) call mpl%abort(subr,'wrong kind for field '//afield%name())

! Get generic FunctionSpace
afunctionspace = afield%functionspace()

select case (afunctionspace%name())
case ('NodeColumns')
   ! Get NodeColumns function space
   afunctionspace_nc = afield%functionspace()

   ! Get number of nodes
   nmga = afunctionspace_nc%nb_nodes()
case ('PointCloud')
   ! Get PointCloud function space
   afunctionspace_pc = afield%functionspace()

   ! Get number of points
   nmga = afunctionspace_pc%size()
case ('StructuredColumns')
   ! Get StructuredColumns function space
   afunctionspace_sc = afield%functionspace()

   ! Get number of nodes
   nmga = afunctionspace_sc%size_owned()
case default
   call mpl%abort(subr,'wrong function space for field '//afield%name()//': '//afunctionspace%name())
end select

! Check number of nodes
if (nmga/=size(fld,1)) call mpl%abort(subr,'wrong number of nodes for field '//afield%name())

! Get number of levels
! - afield%levels() is 0 for 2D ATLAS fields, positive for 3D fields
! - the size of the second dimension of fld is always positive
! - to ensure that sizes are compatible for copying data, we use the minimum between the two
nl0 = min(afield%levels(),size(fld,2))

! Initialization
fld = 0.0

! Copy data
! For the 2D case (afield%levels()==0), the field is copied:
! - at the first level of fld if (lev2d=='first')
! - at the last level of fld if (lev2d=='last')
! NB: an ATLAS field with 1 level only (afield%levels()==1) is considered as a 3D field, so lev2d does not apply
if (nl0==0) then
   if (size(fld,2)>0) then
      call afield%data(ptr_1)
      if (trim(llev2d)=='first') then
         fld(1:nmga,1) = ptr_1(1:nmga)
      elseif (trim(llev2d)=='last') then
         fld(1:nmga,size(fld,2)) = ptr_1(1:nmga)
      end if
   end if
else
   call afield%data(ptr_2)
   fld(1:nmga,1:nl0) = transpose(ptr_2(1:nl0,1:nmga))
end if

end subroutine field_to_fld_real

!----------------------------------------------------------------------
! Subroutine: field_to_fld_logical
! Purpose: convert ATLAS field to field, logical
!----------------------------------------------------------------------
subroutine field_to_fld_logical(mpl,afield,fld,lev2d)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl           ! MPI data
type(atlas_field),intent(in) :: afield        ! ATLAS field
logical,intent(out) :: fld(:,:)               ! Field
character(len=*),intent(in),optional :: lev2d ! Level for 2D variables

! Local variables
integer :: nmga,nl0,imga,il0
integer,allocatable :: fld_int(:,:)
integer,pointer :: ptr_1(:),ptr_2(:,:)
character(len=1024) :: llev2d
character(len=1024),parameter :: subr = 'field_to_fld_logical'
type(atlas_functionspace) :: afunctionspace
type(atlas_functionspace_nodecolumns) :: afunctionspace_nc
type(atlas_functionspace_pointcloud) :: afunctionspace_pc
type(atlas_functionspace_structuredcolumns) :: afunctionspace_sc

! Local lev2d
llev2d = 'first'
if (present(lev2d)) llev2d = lev2d

! Check kind
if (afield%kind()/=atlas_integer(kind_int)) call mpl%abort(subr,'wrong kind for field '//afield%name())

! Get generic FunctionSpace
afunctionspace = afield%functionspace()

select case (afunctionspace%name())
case ('NodeColumns')
   ! Get NodeColumns function space
   afunctionspace_nc = afield%functionspace()

   ! Get number of nodes
   nmga = afunctionspace_nc%nb_nodes()
case ('PointCloud')
   ! Get PointCloud function space
   afunctionspace_pc = afield%functionspace()

   ! Get number of points
   nmga = afunctionspace_pc%size()
case ('StructuredColumns')
   ! Get StructuredColumns function space
   afunctionspace_sc = afield%functionspace()

   ! Get number of nodes
   nmga = afunctionspace_sc%size_owned()
case default
   call mpl%abort(subr,'wrong function space for field '//afield%name()//': '//afunctionspace%name())
end select

! Check number of nodes
if (nmga/=size(fld,1)) call mpl%abort(subr,'wrong number of nodes for field '//afield%name())

! Get number of levels
! - afield%levels() is 0 for 2D ATLAS fields, positive for 3D fields
! - the size of the second dimension of fld is always positive
! - to ensure that sizes are compatible for copying data, we use the minimum between the two
nl0 = min(afield%levels(),size(fld,2))

! Allocation
allocate(fld_int(size(fld,1),size(fld,2)))

! Initialization
fld_int = 0

! Copy data
! For the 2D case (afield%levels()==0), the 2D field is copied:
! - at the first level of fld if (lev2d=='first')
! - at the last level of fld if (lev2d=='last')
! NB: an ATLAS field with 1 level only (afield%levels()==1) is considered as a 3D field, so lev2d does not apply
if (nl0==0) then
   if (size(fld,2)>0) then
      call afield%data(ptr_1)
      if (trim(llev2d)=='first') then
         fld_int(1:nmga,1) = ptr_1(1:nmga)
      elseif (trim(llev2d)=='last') then
         fld_int(1:nmga,size(fld,2)) = ptr_1(1:nmga)
      end if
   end if
else
   call afield%data(ptr_2)
   fld_int(1:nmga,1:nl0) = transpose(ptr_2(1:nl0,1:nmga))
end if

! Integer to logical
do il0=1,size(fld,2)
   do imga=1,size(fld,1)
      if (fld_int(imga,il0)==0) then
         fld(imga,il0) = .false.
      elseif (fld_int(imga,il0)==1) then
         fld(imga,il0) = .true.
      else
         call mpl%abort(subr,'wrong value in 0-1 integer field for field '//afield%name())
      end if
   end do
end do

! Release memory
deallocate(fld_int)

end subroutine field_to_fld_logical

!----------------------------------------------------------------------
! Subroutine: fld_to_field_real
! Purpose: convert field to ATLAS field, real
!----------------------------------------------------------------------
subroutine fld_to_field_real(mpl,fld,afield,lev2d)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl           ! MPI data
real(kind_real),intent(in) :: fld(:,:)        ! Field
type(atlas_field),intent(inout) :: afield     ! ATLAS field
character(len=*),intent(in),optional :: lev2d ! Level for 2D variables

! Local variables
integer :: nmga,nl0
real(kind_real),pointer :: ptr_1(:),ptr_2(:,:)
character(len=1024) :: llev2d
character(len=1024),parameter :: subr = 'fld_to_field_real'
type(atlas_functionspace) :: afunctionspace
type(atlas_functionspace_nodecolumns) :: afunctionspace_nc
type(atlas_functionspace_pointcloud) :: afunctionspace_pc
type(atlas_functionspace_structuredcolumns) :: afunctionspace_sc

! Local lev2d
llev2d = 'first'
if (present(lev2d)) llev2d = lev2d

! Check kind
if (afield%kind()/=atlas_real(kind_real)) call mpl%abort(subr,'wrong kind for field '//afield%name())

! Get generic FunctionSpace
afunctionspace = afield%functionspace()

select case (afunctionspace%name())
case ('NodeColumns')
   ! Get NodeColumns function space
   afunctionspace_nc = afield%functionspace()

   ! Get number of nodes
   nmga = afunctionspace_nc%nb_nodes()
case ('PointCloud')
   ! Get PointCloud function space
   afunctionspace_pc = afield%functionspace()

   ! Get number of points
   nmga = afunctionspace_pc%size()
case ('StructuredColumns')
   ! Get StructuredColumns function space
   afunctionspace_sc = afield%functionspace()

   ! Get number of nodes
   nmga = afunctionspace_sc%size_owned()
case default
   call mpl%abort(subr,'wrong function space for field '//afield%name()//': '//afunctionspace%name())
end select

! Get number of levels
! - afield%levels() is 0 for 2D ATLAS fields, positive for 3D fields
! - the size of the second dimension of fld is always positive
! - to ensure that sizes are compatible for copying data, we use the minimum between the two
nl0 = min(afield%levels(),size(fld,2))

! Check number of nodes
if (nmga/=size(fld,1)) call mpl%abort(subr,'wrong number of nodes for field '//afield%name())

! Copy data
! For the 2D case (afield%levels()==0), the field is copied:
! - at the first level of fld if (lev2d=='first')
! - at the last level of fld if (lev2d=='last')
! NB: an ATLAS field with 1 level only (afield%levels()==1) is considered as a 3D field, so lev2d does not apply
if (nl0==0) then
   if (size(fld,2)>0) then
      call afield%data(ptr_1)
      if (trim(llev2d)=='first') then
         ptr_1(1:nmga) = fld(1:nmga,1)
      elseif (trim(llev2d)=='last') then
         ptr_1(1:nmga) = fld(1:nmga,size(fld,2))
      end if
   end if
else
   call afield%data(ptr_2)
   ptr_2(1:nl0,1:nmga) = transpose(fld(1:nmga,1:nl0))
end if

end subroutine fld_to_field_real

!----------------------------------------------------------------------
! Subroutine: create_atlas_function_space
! Purpose: create ATLAS function space
!----------------------------------------------------------------------
subroutine create_atlas_function_space(nmga,lon_mga,lat_mga,afunctionspace)

implicit none

! Passed variables
integer,intent(in) :: nmga                              ! Number of nodes
real(kind_real),intent(in) :: lon_mga(nmga)             ! Longitudes
real(kind_real),intent(in) :: lat_mga(nmga)             ! Latitudes
type(atlas_functionspace),intent(out) :: afunctionspace ! ATLAS function space

! Local variables
integer :: imga
real(kind_real),pointer :: real_ptr(:,:)
type(atlas_field) :: afield

! Create lon/lat field
afield = atlas_field(name="lonlat",kind=atlas_real(kind_real),shape=(/2,nmga/))
call afield%data(real_ptr)
do imga=1,nmga
   real_ptr(1,imga) = lon_mga(imga)*rad2deg
   real_ptr(2,imga) = lat_mga(imga)*rad2deg
end do

! Create function space PointCloud
afunctionspace = atlas_functionspace_pointcloud(afield)

end subroutine create_atlas_function_space

!----------------------------------------------------------------------
! Subroutine: create_atlas_fieldset
! Purpose: create ATLAS fieldset with empty fields
!----------------------------------------------------------------------
subroutine create_atlas_fieldset(afunctionspace,nl,variables,afieldset)

implicit none

! Passed variables
type(atlas_functionspace),intent(in) :: afunctionspace ! ATLAS function space
integer,intent(in) :: nl                               ! Number of levels
character(len=*),intent(in) :: variables(:)            ! Variables names
type(atlas_fieldset),intent(out) :: afieldset          ! ATLAS fieldset

! Local variables
integer :: iv
type(atlas_field) :: afield

! Set ATLAS fieldset
afieldset = atlas_fieldset()

! Create fields
do iv=1,size(variables)
   ! Create field
   afield = afunctionspace%create_field(name=variables(iv),kind=atlas_real(kind_real),levels=nl)

   ! Add field
   call afieldset%add(afield)

   ! Release pointer
   call afield%final()
end do

end subroutine create_atlas_fieldset

!----------------------------------------------------------------------
! Subroutine: atlas_to_fld
! Purpose: convert ATLAS fieldset to field
!----------------------------------------------------------------------
subroutine atlas_to_fld(mpl,variables,afieldset,fld,lev2d)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl             ! MPI data
character(len=*),intent(in) :: variables(:)     ! Variables names
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset
real(kind_real),intent(out) :: fld(:,:,:)       ! Field
character(len=*),intent(in),optional :: lev2d   ! Level for 2D variables

! Local variables
integer :: iv
character(len=1024) :: fieldname
character(len=1024),parameter :: subr = 'atlas_to_fld'
type(atlas_field) :: afield

! Check number of variables
if (size(variables)/=size(fld,3)) call mpl%abort(subr,'inconsistency in number of variables')

! Loop over fields
do iv=1,size(variables)
   ! Get field
   afield = afieldset%field(variables(iv))

   ! Get field data
   if (present(lev2d)) then
      call field_to_fld(mpl,afield,fld(:,:,iv),lev2d)
   else
      call field_to_fld(mpl,afield,fld(:,:,iv))
   end if
end do

end subroutine atlas_to_fld

!----------------------------------------------------------------------
! Subroutine: fld_to_atlas
! Purpose: convert field to ATLAS fieldset
!----------------------------------------------------------------------
subroutine fld_to_atlas(mpl,variables,fld,afieldset,lev2d)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl             ! MPI data
character(len=*),intent(in) :: variables(:)     ! Variables names
real(kind_real),intent(in) :: fld(:,:,:)        ! Field
type(atlas_fieldset),intent(inout) :: afieldset ! ATLAS fieldset
character(len=*),intent(in),optional :: lev2d   ! Level for 2D variables

! Local variables
integer :: iv,its
character(len=1024) :: fieldname
character(len=1024),parameter :: subr = 'fld_to_atlas'
type(atlas_field) :: afield

! Check number of variables
if (size(variables)/=size(fld,3)) call mpl%abort(subr,'inconsistency in number of variables')

! Loop over fields
do iv=1,size(variables)
   ! Get or create field
   afield = afieldset%field(variables(iv))

   ! Get field data
   if (present(lev2d)) then
      call fld_to_field(mpl,fld(:,:,iv),afield,lev2d)
   else
      call fld_to_field(mpl,fld(:,:,iv),afield)
   end if
end do

end subroutine fld_to_atlas

end module tools_atlas
