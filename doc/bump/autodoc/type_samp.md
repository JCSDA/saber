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
| subroutine | [samp_setup_1](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L789) | setup sampling, first step |
| subroutine | [samp_setup_2](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L981) | setup sampling, second step |
| subroutine | [samp_compute_mask](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1030) | compute mask |
| subroutine | [samp_compute_c1](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1177) | compute sampling, subset Sc1 |
| subroutine | [samp_compute_mpi_c1a](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1282) | compute MPI distribution, halo A, subset Sc1 |
| subroutine | [samp_compute_c3](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1349) | compute sampling, subset Sc3 |
| subroutine | [samp_check_mask](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1490) | check sampling mask |
| subroutine | [samp_compute_c2](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1542) | compute sampling, subset Sc2 |
| subroutine | [samp_compute_mpi_c2a](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1625) | compute sampling MPI distribution, halo A, subset Sc2 |
| subroutine | [samp_compute_mpi_c2b](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1726) | compute sampling MPI distribution, halo B |
| subroutine | [samp_compute_mesh_c2](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1812) | compute sampling mesh, subset Sc2 |
| subroutine | [samp_compute_mpi_c](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1844) | compute sampling MPI distribution, halo C |
| subroutine | [samp_compute_mpi_d](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L1973) | compute sampling MPI distribution, halo D |
| subroutine | [samp_compute_mpi_e](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L2092) | compute sampling MPI distribution, halo E |
| subroutine | [samp_diag_filter](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L2198) | filter diagnostics |
| subroutine | [samp_diag_fill](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_samp.F90#L2363) | fill diagnostics missing values |
