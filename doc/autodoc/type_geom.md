# Module type_geom

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [geom_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L129) | release memory (partial) |
| subroutine | [geom_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L177) | release memory |
| subroutine | [geom_setup](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L197) | setup geometry |
| subroutine | [geom_from_atlas](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L855) | set geometry from ATLAS fieldset |
| subroutine | [geom_remap](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L992) | remap points to improve load balance |
| subroutine | [geom_define_dirac](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1241) | define dirac indices |
| subroutine | [geom_check_arc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1296) | check if an arc is crossing boundaries |
| subroutine | [geom_copy_c0a_to_mga](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1350) | copy from subset Sc0 to model grid, halo A |
| subroutine | [geom_copy_mga_to_c0a_real](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1389) | copy from model grid to subset Sc0, halo A, real |
| subroutine | [geom_copy_mga_to_c0a_logical](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1423) | copy from model grid to subset Sc0, halo A, logical |
| subroutine | [geom_compute_deltas](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1464) | compute deltas for LCT definition |
