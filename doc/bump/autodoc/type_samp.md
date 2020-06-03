# Module type_samp

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [samp_alloc_mask](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L168) | allocation for mask |
| subroutine | [samp_alloc_other](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L189) | allocation for other variables |
| subroutine | [samp_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L223) | release memory (partial) |
| subroutine | [samp_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L272) | release memory |
| subroutine | [samp_read](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L336) | read |
| subroutine | [samp_write](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L465) | write |
| subroutine | [samp_write_grids](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L572) | write |
| subroutine | [samp_setup_1](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L801) | setup sampling, first step |
| subroutine | [samp_setup_2](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L993) | setup sampling, second step |
| subroutine | [samp_compute_mask](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1042) | compute mask |
| subroutine | [samp_compute_c1](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1189) | compute sampling, subset Sc1 |
| subroutine | [samp_compute_mpi_c1a](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1294) | compute MPI distribution, halo A, subset Sc1 |
| subroutine | [samp_compute_c3](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1361) | compute sampling, subset Sc3 |
| subroutine | [samp_check_mask](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1502) | check sampling mask |
| subroutine | [samp_compute_c2](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1554) | compute sampling, subset Sc2 |
| subroutine | [samp_compute_mpi_c2a](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1637) | compute sampling MPI distribution, halo A, subset Sc2 |
| subroutine | [samp_compute_mpi_c2b](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1738) | compute sampling MPI distribution, halo B |
| subroutine | [samp_compute_mesh_c2](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1824) | compute sampling mesh, subset Sc2 |
| subroutine | [samp_compute_mpi_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1856) | compute sampling MPI distribution, halo C |
| subroutine | [samp_compute_mpi_d](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1985) | compute sampling MPI distribution, halo D |
| subroutine | [samp_compute_mpi_e](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L2104) | compute sampling MPI distribution, halo E |
| subroutine | [samp_diag_filter](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L2210) | filter diagnostics |
| subroutine | [samp_diag_fill](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L2375) | fill diagnostics missing values |
