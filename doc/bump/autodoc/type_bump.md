# Module type_bump

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [bump%] [setup_online](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L78) | online setup |
| subroutine | [bump%] [run_drivers](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L303) | run drivers |
| subroutine | [bump%] [add_member](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L516) | add member into bump0,000000e+00ns[1,2] |
| subroutine | [bump%] [apply_vbal](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L564) | vertical balance application |
| subroutine | [bump%] [apply_vbal_inv](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L602) | vertical balance application, inverse |
| subroutine | [bump%] [apply_vbal_ad](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L640) | vertical balance application, adjoint |
| subroutine | [bump%] [apply_vbal_inv_ad](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L678) | vertical balance application, inverse adjoint |
| subroutine | [bump%] [apply_nicas](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L716) | NICAS application |
| subroutine | [bump%] [get_cv_size](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L764) | get control variable size |
| subroutine | [bump%] [apply_nicas_sqrt](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L787) | NICAS square-root application |
| subroutine | [bump%] [apply_nicas_sqrt_ad](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L834) | NICAS square-root adjoint application |
| subroutine | [bump%] [randomize](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L878) | NICAS randomization |
| subroutine | [bump%] [apply_obsop](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L915) | observation operator application |
| subroutine | [bump%] [apply_obsop_ad](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L944) | observation operator adjoint application |
| subroutine | [bump%] [get_parameter](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L973) | get a parameter |
| subroutine | [bump%] [copy_to_field](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L1029) | copy to field |
| subroutine | [bump%] [set_parameter](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L1180) | set a parameter |
| subroutine | [bump%] [copy_from_field](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L1236) | copy from field |
| subroutine | [bump%] [partial_dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L1337) | release memory (partial) |
| subroutine | [bump%] [dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_bump.F90#L1364) | release memory (full) |
