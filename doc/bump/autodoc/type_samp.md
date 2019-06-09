# Module type_samp

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [samp%] [alloc_mask](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L133) | allocation for mask |
| subroutine | [samp%] [alloc_other](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L152) | allocation for other variables |
| subroutine | [samp%] [dealloc](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L186) | release memory |
| subroutine | [samp%] [read](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L241) | read |
| subroutine | [samp%] [write](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L376) | write |
| subroutine | [samp%] [setup_sampling](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L517) | setup sampling |
| subroutine | [samp%] [compute_mask](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L750) | compute mask |
| subroutine | [samp%] [compute_sampling_zs](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L886) | compute zero-separation sampling |
| subroutine | [samp%] [compute_sampling_ps](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1011) | compute positive separation sampling |
| subroutine | [samp%] [compute_sampling_lct](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1142) | compute LCT sampling |
| subroutine | [samp%] [check_mask](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1205) | check sampling mask |
| subroutine | [samp%] [compute_mpi_a](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1277) | compute sampling MPI distribution, halo A |
| subroutine | [samp%] [compute_mpi_ab](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1345) | compute sampling MPI distribution, halos A-B |
| subroutine | [samp%] [compute_mpi_c](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1592) | compute sampling MPI distribution, halo C |
| subroutine | [samp%] [compute_mpi_f](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1772) | compute sampling MPI distribution, halo F |
| subroutine | [samp%] [diag_filter](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1841) | filter diagnostics |
| subroutine | [samp%] [diag_fill](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1974) | fill diagnostics missing values |
