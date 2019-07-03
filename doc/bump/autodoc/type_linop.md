# Module type_linop

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [linop%] [alloc](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L64) | allocation |
| subroutine | [linop%] [dealloc](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L94) | release memory |
| subroutine | [linop%] [copy](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L113) | copy |
| subroutine | [linop%] [read](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L150) | read |
| subroutine | [linop%] [write](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L206) | write |
| subroutine | [linop%] [reorder](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L264) | reorder linear operator |
| subroutine | [linop%] [apply](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L323) | apply linear operator |
| subroutine | [linop%] [apply_ad](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L413) | apply linear operator, adjoint |
| subroutine | [linop%] [apply_sym](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L469) | apply linear operator, symmetric |
| subroutine | [linop%] [add_op](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L536) | add operation |
| subroutine | [linop%] [gather](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L581) | gather data from OpenMP threads |
| subroutine | [linop%] [interp_from_lat_lon](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L615) | compute horizontal interpolation from source latitude/longitude |
| subroutine | [linop%] [interp_from_mesh_tree](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L691) | compute horizontal interpolation from source mesh and tree |
| subroutine | [linop%] [interp_grid](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L905) | compute horizontal grid interpolation |
| subroutine | [linop%] [check_mask](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L1050) | check mask boundaries for linear operators |
| subroutine | [linop%] [interp_missing](https://github.com/JCSDA/saber/src/bump/type_linop.F90#L1118) | deal with missing interpolation points |
