# Module type_nicas_blk

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [balldata_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L265) | allocation |
| subroutine | [balldata_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L283) | release memory |
| subroutine | [balldata_pack](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L301) | pack data into balldata object |
| subroutine | [nicas_blk_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L340) | release memory (partial) |
| subroutine | [nicas_blk_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L424) | release memory (full) |
| subroutine | [nicas_blk_read](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L490) | read |
| subroutine | [nicas_blk_write](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L684) | write |
| subroutine | [nicas_blk_write_grids](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L826) | write NICAS grids |
| subroutine | [nicas_blk_receive](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L910) | receive |
| subroutine | [nicas_blk_send](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1171) | send |
| subroutine | [nicas_blk_compute_parameters](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1378) | compute NICAS parameters |
| subroutine | [nicas_blk_compute_parameters_smoother](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1511) | compute NICAS parameters for a smoother |
| subroutine | [nicas_blk_compute_sampling_c1](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1589) | compute NICAS sampling, subset Sc1 |
| subroutine | [nicas_blk_compute_sampling_v](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1739) | compute NICAS sampling, vertical dimension |
| subroutine | [nicas_blk_compute_mpi_a](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1835) | compute NICAS MPI distribution, halos A |
| subroutine | [nicas_blk_compute_sampling_c2](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1905) | compute NICAS sampling, subset Sc2 |
| subroutine | [nicas_blk_compute_mpi_ab](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L2033) | compute NICAS MPI distribution, halos A-B |
| subroutine | [nicas_blk_compute_interp_v](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L2297) | compute vertical interpolation |
| subroutine | [nicas_blk_compute_convol](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L2380) | compute convolution |
| subroutine | [nicas_blk_compute_convol_network](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L2804) | compute convolution with a network approach |
| subroutine | [nicas_blk_compute_convol_distance](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3039) | compute convolution with a distance approach |
| subroutine | [nicas_blk_compute_convol_weights](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3213) | compute convolution weights |
| subroutine | [nicas_blk_compute_mpi_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3320) | compute NICAS MPI distribution, halo C |
| subroutine | [nicas_blk_compute_internal_normalization](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3470) | compute internal normalization |
| subroutine | [nicas_blk_compute_normalization](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3540) | compute normalization |
| subroutine | [nicas_blk_compute_grids](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3789) | compute grids |
| subroutine | [nicas_blk_compute_adv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3858) | compute advection |
| subroutine | [nicas_blk_apply](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4001) | apply NICAS method |
| subroutine | [nicas_blk_apply_from_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4080) | apply NICAS method from its square-root formulation |
| subroutine | [nicas_blk_apply_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4105) | apply NICAS method square-root |
| subroutine | [nicas_blk_apply_sqrt_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4149) | apply NICAS method square-root adjoint |
| subroutine | [nicas_blk_apply_interp](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4196) | apply interpolation |
| subroutine | [nicas_blk_apply_interp_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4225) | apply interpolation adjoint |
| subroutine | [nicas_blk_apply_interp_h](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4254) | apply horizontal interpolation |
| subroutine | [nicas_blk_apply_interp_h_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4285) | apply horizontal interpolation adjoint |
| subroutine | [nicas_blk_apply_interp_v](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4315) | apply vertical interpolation |
| subroutine | [nicas_blk_apply_interp_v_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4364) | apply vertical interpolation adjoint |
| subroutine | [nicas_blk_apply_interp_s](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4407) | apply subsampling interpolation |
| subroutine | [nicas_blk_apply_interp_s_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4444) | apply subsampling interpolation adjoint |
| subroutine | [nicas_blk_apply_convol](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4478) | apply convolution |
| subroutine | [nicas_blk_apply_adv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4512) | apply advection |
| subroutine | [nicas_blk_apply_adv_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4547) | apply advection |
| subroutine | [nicas_blk_apply_adv_inv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4582) | apply inverse advection |
| subroutine | [nicas_blk_test_adjoint](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4617) | test NICAS adjoint accuracy |
| subroutine | [nicas_blk_test_dirac](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4854) | apply NICAS to diracs |
