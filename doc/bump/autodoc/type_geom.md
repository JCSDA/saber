# Module type_geom

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [geom%] [alloc](https://github.com/JCSDA/saber/src/bump/type_geom.F90#L116) | allocation |
| subroutine | [geom%] [dealloc](https://github.com/JCSDA/saber/src/bump/type_geom.F90#L142) | release memory |
| subroutine | [geom%] [setup](https://github.com/JCSDA/saber/src/bump/type_geom.F90#L192) | setup geometry |
| subroutine | [geom%] [find_sc0](https://github.com/JCSDA/saber/src/bump/type_geom.F90#L552) | find subset Sc0 points |
| subroutine | [geom%] [define_dirac](https://github.com/JCSDA/saber/src/bump/type_geom.F90#L658) | define dirac indices |
| subroutine | [geom%] [reorder_points](https://github.com/JCSDA/saber/src/bump/type_geom.F90#L713) | reorder Sc0 points based on lon/lat |
| subroutine | [geom%] [check_arc](https://github.com/JCSDA/saber/src/bump/type_geom.F90#L767) | check if an arc is crossing boundaries |
| subroutine | [geom%] [copy_c0a_to_mga](https://github.com/JCSDA/saber/src/bump/type_geom.F90#L820) | copy from subset Sc0 to model grid, halo A |
| subroutine | [geom%] [copy_mga_to_c0a](https://github.com/JCSDA/saber/src/bump/type_geom.F90#L859) | copy from model grid to subset Sc0, halo A |
| subroutine | [geom%] [compute_deltas](https://github.com/JCSDA/saber/src/bump/type_geom.F90#L921) | compute deltas for LCT definition |
