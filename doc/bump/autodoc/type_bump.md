# Module type_bump

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [bump%] [setup_online](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L79) | online setup |
| subroutine | [bump%] [run_drivers](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L277) | run drivers |
| subroutine | [bump%] [add_member](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L481) | add member into bump0,000000e+00ns[1,2] |
| subroutine | [bump%] [apply_vbal](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L529) | vertical balance application |
| subroutine | [bump%] [apply_vbal_inv](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L567) | vertical balance application, inverse |
| subroutine | [bump%] [apply_vbal_ad](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L605) | vertical balance application, adjoint |
| subroutine | [bump%] [apply_vbal_inv_ad](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L643) | vertical balance application, inverse adjoint |
| subroutine | [bump%] [apply_nicas](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L681) | NICAS application |
| subroutine | [bump%] [get_cv_size](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L729) | get control variable size |
| subroutine | [bump%] [apply_nicas_sqrt](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L752) | NICAS square-root application |
| subroutine | [bump%] [apply_nicas_sqrt_ad](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L799) | NICAS square-root adjoint application |
| subroutine | [bump%] [randomize](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L843) | NICAS randomization |
| subroutine | [bump%] [apply_obsop](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L880) | observation operator application |
| subroutine | [bump%] [apply_obsop_ad](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L909) | observation operator adjoint application |
| subroutine | [bump%] [get_parameter](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L938) | get a parameter |
| subroutine | [bump%] [copy_to_field](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L994) | copy to field |
| subroutine | [bump%] [set_parameter](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L1145) | set a parameter |
| subroutine | [bump%] [copy_from_field](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L1201) | copy from field |
| subroutine | [bump%] [crtm_neighbors_3d](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L1302) | find nearest neighbors for CRTM, 3D |
| subroutine | [bump%] [crtm_neighbors_2d](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L1340) | find nearest neighbors for CRTM, 2D |
| subroutine | [bump%] [partial_dealloc](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L1376) | release memory (partial) |
| subroutine | [bump%] [dealloc](https://github.com/JCSDA/saber/src/bump/type_bump.F90#L1403) | release memory (full) |
