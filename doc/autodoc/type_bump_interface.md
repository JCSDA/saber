# Module type_bump_interface

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [bump_create_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L22) | create |
| subroutine | [bump_run_drivers_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L61) | run drivers |
| subroutine | [bump_add_member_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L83) | add member into bump%ens[1,2] |
| subroutine | [bump_remove_member_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L110) | remove member into bump%ens[1,2] |
| subroutine | [bump_apply_vbal_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L137) | vertical balance application |
| subroutine | [bump_apply_vbal_inv_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L162) | vertical balance application, inverse |
| subroutine | [bump_apply_vbal_ad_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L187) | vertical balance application, adjoint |
| subroutine | [bump_apply_vbal_inv_ad_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L212) | vertical balance application, inverse adjoint |
| subroutine | [bump_apply_stddev_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L237) | standard-deviation application |
| subroutine | [bump_apply_stddev_inv_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L262) | standard-deviation application, inverse |
| subroutine | [bump_apply_nicas_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L287) | NICAS application |
| subroutine | [bump_get_cv_size_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L312) | get control variable size |
| subroutine | [bump_apply_nicas_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L335) | NICAS square-root application |
| subroutine | [bump_apply_nicas_sqrt_ad_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L361) | NICAS square-root adjoint application |
| subroutine | [bump_randomize_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L387) | NICAS randomization |
| subroutine | [bump_get_parameter_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L412) | get a parameter |
| subroutine | [bump_set_parameter_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L445) | set a parameter |
| subroutine | [bump_dealloc_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bump_interface.F90#L478) | deallocation |
