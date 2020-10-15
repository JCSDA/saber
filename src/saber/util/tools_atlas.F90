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

interface field_to_array
  module procedure field_to_array_real
  module procedure field_to_array_logical
end interface

interface field_from_array
  module procedure field_from_array_real
end interface

private
public :: field_to_array,field_from_array,create_atlas_function_space

contains

!----------------------------------------------------------------------
! Subroutine: field_to_array_real
! Purpose: convert ATLAS field to field, real
!----------------------------------------------------------------------
subroutine field_to_array_real(mpl,afield,fld,lev2d)

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
character(len=1024),parameter :: subr = 'field_to_array_real'
type(atlas_functionspace) :: afunctionspace
type(atlas_functionspace_nodecolumns) :: afunctionspace_nc
type(atlas_functionspace_pointcloud) :: afunctionspace_pc
type(atlas_functionspace_structuredcolumns) :: afunctionspace_sc

! Local lev2d
llev2d = 'first'
if (present(lev2d)) llev2d = lev2d

! Check kind
if (afield%kind()/=atlas_real(kind_real)) call mpl%abort(subr,'wrong kind for field '//afield%name())

! Get generic functionspace
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

end subroutine field_to_array_real

!----------------------------------------------------------------------
! Subroutine: field_to_array_logical
! Purpose: convert ATLAS field to field, logical
!----------------------------------------------------------------------
subroutine field_to_array_logical(mpl,afield,fld,lev2d)

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
character(len=1024),parameter :: subr = 'field_to_array_logical'
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

end subroutine field_to_array_logical

!----------------------------------------------------------------------
! Subroutine: field_from_array_real
! Purpose: convert field to ATLAS field, real
!----------------------------------------------------------------------
subroutine field_from_array_real(mpl,fld,afield,lev2d)

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
character(len=1024),parameter :: subr = 'field_from_array_real'
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

end subroutine field_from_array_real

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

end module tools_atlas
