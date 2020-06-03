# Module type_bump

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [bump_create](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L122) | create |
| subroutine | [bump_setup](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L166) | setup |
| subroutine | [bump_setup_online_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L371) | online setup (deprecated) |
| subroutine | [bump_run_drivers](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L511) | run drivers |
| subroutine | [bump_add_member](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L771) | add member into bump%ens[1,2] |
| subroutine | [bump_remove_member](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L861) | remove member into bump%ens[1,2] |
| subroutine | [bump_apply_vbal](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L915) | vertical balance application |
| subroutine | [bump_apply_vbal_inv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L960) | vertical balance application, inverse |
| subroutine | [bump_apply_vbal_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1005) | vertical balance application, adjoint |
| subroutine | [bump_apply_vbal_inv_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1050) | vertical balance application, inverse adjoint |
| subroutine | [bump_apply_stddev](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1095) | standard-deviation application |
| subroutine | [bump_apply_stddev_inv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1142) | standard-deviation application, inverse |
| subroutine | [bump_apply_nicas](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1189) | NICAS application |
| subroutine | [bump_apply_nicas_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1244) | NICAS application (deprecated) |
| subroutine | [bump_get_cv_size](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1296) | get control variable size |
| subroutine | [bump_apply_nicas_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1319) | NICAS square-root application |
| subroutine | [bump_apply_nicas_sqrt_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1370) | NICAS square-root application (deprecated) |
| subroutine | [bump_apply_nicas_sqrt_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1420) | NICAS square-root adjoint application |
| subroutine | [bump_randomize](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1468) | NICAS randomization |
| subroutine | [bump_apply_obsop](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1509) | observation operator application |
| subroutine | [bump_apply_obsop_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1547) | observation operator application (deprecated) |
| subroutine | [bump_apply_obsop_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1584) | observation operator adjoint application |
| subroutine | [bump_apply_obsop_ad_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1622) | observation operator adjoint application (deprecated) |
| subroutine | [bump_get_parameter](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1659) | get a parameter |
| subroutine | [bump_copy_to_field](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1723) | copy to field |
| subroutine | [bump_test_get_parameter](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1900) | test get_parameter |
| subroutine | [bump_set_parameter](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1956) | set a parameter |
| subroutine | [bump_set_parameter_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2019) | set a parameter (deprecated) |
| subroutine | [bump_copy_from_field](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2082) | copy from field |
| subroutine | [bump_test_set_parameter](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2179) | test set_parameter |
| subroutine | [bump_test_apply_interfaces](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2263) | test BUMP apply interfaces |
| subroutine | [bump_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2411) | release memory (partial) |
| subroutine | [bump_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2439) | release memory (full) |
| subroutine | [dummy](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2468) | dummy finalization |
