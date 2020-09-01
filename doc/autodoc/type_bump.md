# Module type_bump

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [bump_create](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L122) | create |
| subroutine | [bump_setup](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L166) | setup |
| subroutine | [bump_setup_online_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L368) | online setup (deprecated) |
| subroutine | [bump_run_drivers](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L509) | run drivers |
| subroutine | [bump_add_member](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L769) | add member into bump%ens[1,2] |
| subroutine | [bump_remove_member](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L859) | remove member into bump%ens[1,2] |
| subroutine | [bump_apply_vbal](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L913) | vertical balance application |
| subroutine | [bump_apply_vbal_inv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L958) | vertical balance application, inverse |
| subroutine | [bump_apply_vbal_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1003) | vertical balance application, adjoint |
| subroutine | [bump_apply_vbal_inv_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1048) | vertical balance application, inverse adjoint |
| subroutine | [bump_apply_stddev](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1093) | standard-deviation application |
| subroutine | [bump_apply_stddev_inv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1140) | standard-deviation application, inverse |
| subroutine | [bump_apply_nicas](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1187) | NICAS application |
| subroutine | [bump_apply_nicas_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1242) | NICAS application (deprecated) |
| subroutine | [bump_get_cv_size](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1294) | get control variable size |
| subroutine | [bump_apply_nicas_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1317) | NICAS square-root application |
| subroutine | [bump_apply_nicas_sqrt_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1368) | NICAS square-root application (deprecated) |
| subroutine | [bump_apply_nicas_sqrt_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1418) | NICAS square-root adjoint application |
| subroutine | [bump_randomize](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1466) | NICAS randomization |
| subroutine | [bump_apply_obsop](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1507) | observation operator application |
| subroutine | [bump_apply_obsop_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1545) | observation operator application (deprecated) |
| subroutine | [bump_apply_obsop_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1582) | observation operator adjoint application |
| subroutine | [bump_apply_obsop_ad_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1620) | observation operator adjoint application (deprecated) |
| subroutine | [bump_get_parameter](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1657) | get a parameter |
| subroutine | [bump_copy_to_field](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1721) | copy to field |
| subroutine | [bump_test_get_parameter](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1898) | test get_parameter |
| subroutine | [bump_set_parameter](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L1954) | set a parameter |
| subroutine | [bump_set_parameter_deprecated](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2017) | set a parameter (deprecated) |
| subroutine | [bump_copy_from_field](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2080) | copy from field |
| subroutine | [bump_test_set_parameter](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2177) | test set_parameter |
| subroutine | [bump_test_apply_interfaces](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2258) | test BUMP apply interfaces |
| subroutine | [bump_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2406) | release memory (partial) |
| subroutine | [bump_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2434) | release memory (full) |
| subroutine | [dummy](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump.F90#L2463) | dummy finalization |
