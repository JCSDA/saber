# Module tools_stripack

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [addnod](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L27) | add a node to a triangulation |
| function | [areas](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L297) | compute the area of a spherical triangle on the unit sphere |
| subroutine | [bdyadd](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L442) | add a boundary node to a triangulation |
| subroutine | [bnodes](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L584) | return the boundary nodes of a triangulation |
| subroutine | [circum](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L720) | return the circumcenter of a spherical triangle |
| subroutine | [covsph](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L807) | connect an exterior node to boundary nodes, covering the sphere |
| subroutine | [det](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L919) | compute 3D determinant |
| subroutine | [crlist](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L956) | return triangle circumcenters and other information |
| subroutine | [insert](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L1596) | insert K as a neighbor of N1 |
| function | [inside](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L1655) | determine if a point is inside a polygonal region |
| subroutine | [intadd](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L2018) | add an interior node to a triangulation |
| subroutine | [intrsc](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L2120) | find the intersection of two great circles |
| subroutine | [jrand](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L2228) | return a random integer between 1 and N |
| subroutine | [left](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L2294) | determin whether a node is to the left of a plane through the origin |
| subroutine | [lstptr](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L2357) | return the index of NB in the adjacency list |
| function | [nbcnt](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L2438) | return the number of neighbors of a node |
| function | [nearnd](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L2514) | return the nearest node to a given point |
| subroutine | [swap](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L2838) | replace the diagonal arc of a quadrilateral with the other diagonal |
| subroutine | [swptst](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L2955) | decide whether to replace a diagonal arc by the other |
| subroutine | [trfind](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L3056) | locate a point relative to a triangulation |
| subroutine | [trlist](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L3570) | convert a triangulation data structure to a triangle list |
| subroutine | [trmesh](https://github.com/JCSDA/saber/src/bump/tools_stripack.F90#L3870) | create a Delaunay triangulation on the unit sphere |
