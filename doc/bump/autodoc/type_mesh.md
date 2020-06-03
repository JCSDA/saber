# Module type_mesh

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [mesh_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L78) | allocation |
| subroutine | [mesh_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L108) | intialization |
| subroutine | [mesh_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L220) | release memory |
| subroutine | [mesh_copy](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L257) | copy |
| subroutine | [mesh_store](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L316) | store mesh cartesian coordinates |
| subroutine | [mesh_trlist](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L348) | compute triangle list, arc list |
| subroutine | [mesh_bnodes](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L407) | find boundary nodes |
| subroutine | [mesh_find_bdist](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L482) | find shortest distance to boundary arcs |
| subroutine | [mesh_check](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L539) | check whether the mesh is made of counter-clockwise triangles |
| subroutine | [mesh_inside](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L611) | find whether a point is inside the mesh |
| subroutine | [mesh_barycentric](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L643) | compute barycentric coordinates |
| subroutine | [mesh_count_bnda](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L679) | count boundary arcs |
| subroutine | [mesh_get_bnda](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mesh.F90#L722) | get boundary arcs |
