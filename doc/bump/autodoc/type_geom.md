# Module type_geom

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [geom_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L154) | release memory (partial) |
| subroutine | [geom_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L216) | release memory |
| subroutine | [geom_setup](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L235) | setup geometry |
| subroutine | [geom_from_atlas](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L356) | set geometry from ATLAS fieldset |
| subroutine | [geom_define_universe](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L497) | define universe |
| subroutine | [geom_setup_c0](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L628) | setup subset Sc0 |
| subroutine | [geom_setup_independent_levels](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L947) | setup independent levels |
| subroutine | [geom_setup_mask_distance](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L993) | setup minimum distance to mask |
| subroutine | [geom_setup_mask_check](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1046) | setup mask checking tool |
| subroutine | [geom_index_from_lonlat](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1104) | get nearest neighbor index from longitude/latitude/level |
| subroutine | [geom_define_dirac](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1174) | define dirac indices |
| subroutine | [geom_check_arc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1229) | check if an arc is crossing boundaries |
| subroutine | [geom_copy_c0a_to_mga](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1283) | copy from subset Sc0 to model grid, halo A |
| subroutine | [geom_copy_mga_to_c0a_real](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1322) | copy from model grid to subset Sc0, halo A, real |
| subroutine | [geom_copy_mga_to_c0a_logical](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1356) | copy from model grid to subset Sc0, halo A, logical |
| subroutine | [geom_compute_deltas](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1397) | compute deltas for LCT definition |
| subroutine | [geom_rand_level](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1424) | select random level |
| subroutine | [geom_rand_point](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1449) | select random valid point on the horizontal grid |
| function | [geom_c0_to_c0a](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1501) | conversion from global to halo A on subset Sc0 |
| function | [geom_c0_to_proc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1527) | conversion from global to processor on subset Sc0 |
| function | [geom_c0_to_c0u](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_geom.F90#L1549) | conversion from global to universe on subset Sc0 |
