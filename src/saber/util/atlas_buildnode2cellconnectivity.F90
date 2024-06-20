module atlas_buildnode2cellconnectivity

implicit none

!-------------------------------------------------------------------------------

interface
   subroutine atlas__build_node_to_cell_connectivity(mesh) bind(C,name="atlas__build_node_to_cell_connectivity")
      use iso_c_binding, only: c_ptr
      type(c_ptr), value :: mesh
   end subroutine
end interface

contains

subroutine atlas_build_node_to_cell_connectivity(mesh)
   use atlas_mesh_module, only: atlas_mesh
   type(atlas_mesh), intent(inout) :: mesh
   call atlas__build_node_to_cell_connectivity(mesh%c_ptr())
end subroutine atlas_build_node_to_cell_connectivity

!-------------------------------------------------------------------------------

end module
