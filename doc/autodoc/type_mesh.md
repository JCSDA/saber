# Module type_mesh

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [mesh_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L71) | allocation |
| subroutine | [mesh_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L100) | intialization |
| subroutine | [mesh_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L161) | release memory |
| subroutine | [mesh_copy](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L196) | copy |
| subroutine | [mesh_store](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L252) | store mesh cartesian coordinates |
| subroutine | [mesh_trlist](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L284) | compute triangle list, arc list |
| subroutine | [mesh_bnodes](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L343) | find boundary nodes |
| subroutine | [mesh_find_bdist](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L417) | find shortest distance to boundary arcs |
| subroutine | [mesh_check](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L474) | check whether the mesh is made of counter-clockwise triangles |
| subroutine | [mesh_inside](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L634) | find whether a point is inside the mesh |
| subroutine | [mesh_barycentric](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L666) | compute barycentric coordinates |
