# Module type_nicas_blk

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [balldata_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L255) | allocation |
| subroutine | [balldata_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L273) | release memory |
| subroutine | [balldata_pack](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L291) | pack data into balldata object |
| subroutine | [nicas_blk_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L330) | release memory (partial) |
| subroutine | [nicas_blk_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L406) | release memory (full) |
| subroutine | [nicas_blk_read](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L471) | read |
| subroutine | [nicas_blk_write](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L661) | write |
| subroutine | [nicas_blk_write_grids](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L799) | write NICAS grids |
| subroutine | [nicas_blk_receive](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L883) | receive |
| subroutine | [nicas_blk_send](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1141) | send |
| subroutine | [nicas_blk_compute_parameters](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1346) | compute NICAS parameters |
| subroutine | [nicas_blk_compute_parameters_horizontal_smoother](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1488) | compute NICAS parameters for a horizontal smoother |
| subroutine | [nicas_blk_compute_sampling_c1](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1563) | compute NICAS sampling, subset Sc1 |
| subroutine | [nicas_blk_compute_sampling_v](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1681) | compute NICAS sampling, vertical dimension |
| subroutine | [nicas_blk_compute_sampling_c2](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1817) | compute NICAS sampling, subset Sc2 |
| subroutine | [nicas_blk_compute_mpi_a](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L1954) | compute NICAS MPI distribution, halos A |
| subroutine | [nicas_blk_compute_mpi_ab](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L2004) | compute NICAS MPI distribution, halos A-B |
| subroutine | [nicas_blk_compute_interp_v](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L2229) | compute vertical interpolation |
| subroutine | [nicas_blk_compute_convol](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L2318) | compute convolution |
| subroutine | [nicas_blk_compute_convol_network](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L2754) | compute convolution with a network approach |
| subroutine | [nicas_blk_compute_convol_distance](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3056) | compute convolution with a distance approach |
| subroutine | [nicas_blk_compute_convol_weights](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3249) | compute convolution weights |
| subroutine | [nicas_blk_compute_mpi_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3361) | compute NICAS MPI distribution, halo C |
| subroutine | [nicas_blk_compute_internal_normalization](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3509) | compute internal normalization |
| subroutine | [nicas_blk_compute_normalization](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3583) | compute normalization |
| subroutine | [nicas_blk_compute_grids](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3839) | compute grids |
| subroutine | [nicas_blk_compute_adv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L3912) | compute advection |
| subroutine | [nicas_blk_apply](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4066) | apply NICAS method |
| subroutine | [nicas_blk_apply_from_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4145) | apply NICAS method from its square-root formulation |
| subroutine | [nicas_blk_apply_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4170) | apply NICAS method square-root |
| subroutine | [nicas_blk_apply_sqrt_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4214) | apply NICAS method square-root adjoint |
| subroutine | [nicas_blk_apply_interp](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4261) | apply interpolation |
| subroutine | [nicas_blk_apply_interp_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4290) | apply interpolation adjoint |
| subroutine | [nicas_blk_apply_interp_h](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4319) | apply horizontal interpolation |
| subroutine | [nicas_blk_apply_interp_h_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4350) | apply horizontal interpolation adjoint |
| subroutine | [nicas_blk_apply_interp_v](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4380) | apply vertical interpolation |
| subroutine | [nicas_blk_apply_interp_v_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4429) | apply vertical interpolation adjoint |
| subroutine | [nicas_blk_apply_interp_s](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4472) | apply subsampling interpolation |
| subroutine | [nicas_blk_apply_interp_s_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4509) | apply subsampling interpolation adjoint |
| subroutine | [nicas_blk_apply_convol](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4543) | apply convolution |
| subroutine | [nicas_blk_apply_adv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4577) | apply advection |
| subroutine | [nicas_blk_apply_adv_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4612) | apply advection |
| subroutine | [nicas_blk_apply_adv_inv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4647) | apply inverse advection |
| subroutine | [nicas_blk_test_adjoint](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4682) | test NICAS adjoint accuracy |
| subroutine | [nicas_blk_test_dirac](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas_blk.F90#L4919) | apply NICAS to diracs |
