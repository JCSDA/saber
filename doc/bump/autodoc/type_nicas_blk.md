# Module type_nicas_blk

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [balldata_alloc](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L227) | allocation |
| subroutine | [balldata_dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L245) | release memory |
| subroutine | [balldata_pack](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L263) | pack data into balldata object |
| subroutine | [nicas_blk%] [partial_dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L302) | release memory (partial) |
| subroutine | [nicas_blk%] [dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L391) | release memory (full) |
| subroutine | [nicas_blk%] [read](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L446) | read |
| subroutine | [nicas_blk%] [write](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L606) | write |
| subroutine | [nicas_blk%] [receive](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L725) | receive |
| subroutine | [nicas_blk%] [send](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L913) | send |
| subroutine | [nicas_blk%] [compute_parameters](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L1073) | compute NICAS parameters |
| subroutine | [nicas_blk%] [compute_sampling_c1](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L1194) | compute NICAS sampling, subset Sc1 |
| subroutine | [nicas_blk%] [compute_sampling_v](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L1309) | compute NICAS sampling, vertical dimension |
| subroutine | [nicas_blk%] [compute_sampling_c2](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L1436) | compute NICAS sampling, subset Sc2 |
| subroutine | [nicas_blk%] [write_sampling](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L1572) | write NICAS sampling |
| subroutine | [nicas_blk%] [compute_mpi_a](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L1637) | compute NICAS MPI distribution, halos A |
| subroutine | [nicas_blk%] [compute_mpi_ab](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L1687) | compute NICAS MPI distribution, halos A-B |
| subroutine | [nicas_blk%] [compute_interp_v](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L1909) | compute vertical interpolation |
| subroutine | [nicas_blk%] [compute_convol](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L1998) | compute convolution |
| subroutine | [nicas_blk%] [compute_convol_network](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L2508) | compute convolution with a network approach |
| subroutine | [nicas_blk%] [compute_convol_distance](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L2801) | compute convolution with a distance approach |
| subroutine | [nicas_blk%] [compute_convol_weights](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L2982) | compute convolution weights |
| subroutine | [nicas_blk%] [compute_mpi_c](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3092) | compute NICAS MPI distribution, halo C |
| subroutine | [nicas_blk%] [compute_normalization](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3215) | compute normalization |
| subroutine | [nicas_blk%] [compute_adv](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3471) | compute advection |
| subroutine | [nicas_blk%] [apply](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3621) | apply NICAS method |
| subroutine | [nicas_blk%] [apply_from_sqrt](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3674) | apply NICAS method from its square-root formulation |
| subroutine | [nicas_blk%] [apply_sqrt](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3699) | apply NICAS method square-root |
| subroutine | [nicas_blk%] [apply_sqrt_ad](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3737) | apply NICAS method square-root adjoint |
| subroutine | [nicas_blk%] [apply_interp](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3775) | apply interpolation |
| subroutine | [nicas_blk%] [apply_interp_ad](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3807) | apply interpolation adjoint |
| subroutine | [nicas_blk%] [apply_interp_h](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3840) | apply horizontal interpolation |
| subroutine | [nicas_blk%] [apply_interp_h_ad](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3871) | apply horizontal interpolation adjoint |
| subroutine | [nicas_blk%] [apply_interp_v](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3901) | apply vertical interpolation |
| subroutine | [nicas_blk%] [apply_interp_v_ad](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3950) | apply vertical interpolation adjoint |
| subroutine | [nicas_blk%] [apply_interp_s](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L3993) | apply subsampling interpolation |
| subroutine | [nicas_blk%] [apply_interp_s_ad](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L4030) | apply subsampling interpolation adjoint |
| subroutine | [nicas_blk%] [apply_convol](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L4064) | apply convolution |
| subroutine | [nicas_blk%] [apply_adv](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L4082) | apply advection |
| subroutine | [nicas_blk%] [apply_adv_ad](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L4117) | apply advection |
| subroutine | [nicas_blk%] [apply_adv_inv](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L4152) | apply inverse advection |
| subroutine | [nicas_blk%] [test_adjoint](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L4187) | test NICAS adjoint accuracy |
| subroutine | [nicas_blk%] [test_dirac](https://github.com/JCSDA/saber/src/saber/bump/type_nicas_blk.F90#L4424) | apply NICAS to diracs |
