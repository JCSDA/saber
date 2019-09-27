# Module type_mesh

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [mesh%] [alloc](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L71) | allocation |
| subroutine | [mesh%] [init](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L100) | intialization |
| subroutine | [mesh%] [dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L161) | release memory |
| subroutine | [mesh%] [copy](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L196) | copy |
| subroutine | [mesh%] [store](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L245) | store mesh cartesian coordinates |
| subroutine | [mesh%] [trlist](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L277) | compute triangle list, arc list |
| subroutine | [mesh%] [bnodes](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L336) | find boundary nodes |
| subroutine | [mesh%] [find_bdist](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L410) | find shortest distance to boundary arcs |
| subroutine | [mesh%] [check](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L467) | check whether the mesh is made of counter-clockwise triangles |
| subroutine | [mesh%] [inside](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L627) | find whether a point is inside the mesh |
| subroutine | [mesh%] [barycentric](https://github.com/JCSDA/saber/src/saber/bump/type_mesh.F90#L659) | compute barycentric coordinates |
