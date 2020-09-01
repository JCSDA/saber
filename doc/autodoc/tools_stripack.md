# Module tools_stripack

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [addnod](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L27) | add a node to a triangulation |
| subroutine | [bdyadd](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L281) | add a boundary node to a triangulation |
| subroutine | [bnodes](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L423) | return the boundary nodes of a triangulation |
| subroutine | [covsph](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L559) | connect an exterior node to boundary nodes, covering the sphere |
| subroutine | [det](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L671) | compute 3D determinant |
| subroutine | [insert](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L707) | insert K as a neighbor of N1 |
| function | [inside](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L766) | determine if a point is inside a polygonal region |
| subroutine | [intadd](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L1129) | add an interior node to a triangulation |
| subroutine | [intrsc](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L1231) | find the intersection of two great circles |
| subroutine | [jrand](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L1339) | return a random integer between 1 and N |
| subroutine | [left](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L1405) | determin whether a node is to the left of a plane through the origin |
| subroutine | [lstptr](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L1468) | return the index of NB in the adjacency list |
| subroutine | [swap](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L1549) | replace the diagonal arc of a quadrilateral with the other diagonal |
| subroutine | [swptst](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L1666) | decide whether to replace a diagonal arc by the other |
| subroutine | [trfind](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L1767) | locate a point relative to a triangulation |
| subroutine | [trlist](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L2281) | convert a triangulation data structure to a triangle list |
| subroutine | [trmesh](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_stripack.F90#L2581) | create a Delaunay triangulation on the unit sphere |
| subroutine | [create_gaugrid](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/tools_stripack.F90#L46) | Create Gaussian grid |
| subroutine | [delete_gaugrid](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/tools_stripack.F90#L62) | Delete Gaussian grid |
| subroutine | [gaugrid_alloc_coord](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/tools_stripack.F90#L76) | allocate Gaussian grid coordinate |
| subroutine | [gaugrid_dealloc_coord](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/tools_stripack.F90#L93) | deallocate Gaussian grid coordinate |
| subroutine | [gaugrid_alloc_field](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/tools_stripack.F90#L109) | allocate Gaussian grid field |
| subroutine | [gaugrid_dealloc_field](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/tools_stripack.F90#L123) | deallocate Gaussian grid field |
| subroutine | [gaugrid_calc_glb_latlon](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/tools_stripack.F90#L136) | calculate global Gaussian latitudes and longitudes |
| subroutine | [gaugrid_fld3d_pointer](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/tools_stripack.F90#L170) | Set 3D field pointer |
| subroutine | [gaugrid_fld2d_pointer](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/tools_stripack.F90#L187) | Set 2D field pointer |
| subroutine | [equals](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/tools_stripack.F90#L204) | create new gaussian grid from other |
| subroutine | [dummy](https://github.com/JCSDA/saber/tree/develop/src/saber/interpolation/tools_stripack.F90#L665) | dummy finalization |
