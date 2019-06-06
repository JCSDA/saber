# Module type_samp

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [samp%] [alloc_mask](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L132) | allocation for mask |
| subroutine | [samp%] [alloc_other](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L151) | allocation for other variables |
| subroutine | [samp%] [dealloc](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L184) | release memory |
| subroutine | [samp%] [read](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L238) | read |
| subroutine | [samp%] [write](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L364) | write |
| subroutine | [samp%] [setup_sampling](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L497) | setup sampling |
| subroutine | [samp%] [compute_mask](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L730) | compute mask |
| subroutine | [samp%] [compute_sampling_zs](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L865) | compute zero-separation sampling |
| subroutine | [samp%] [compute_sampling_ps](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L982) | compute positive separation sampling |
| subroutine | [samp%] [compute_sampling_lct](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1113) | compute LCT sampling |
| subroutine | [samp%] [check_mask](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1188) | check sampling mask |
| subroutine | [samp%] [compute_mpi_a](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1263) | compute sampling MPI distribution, halo A |
| subroutine | [samp%] [compute_mpi_ab](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1331) | compute sampling MPI distribution, halos A-B |
| subroutine | [samp%] [compute_mpi_c](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1578) | compute sampling MPI distribution, halo C |
| subroutine | [samp%] [compute_mpi_f](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1758) | compute sampling MPI distribution, halo F |
| subroutine | [samp%] [diag_filter](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1827) | filter diagnostics |
| subroutine | [samp%] [diag_fill](https://github.com/JCSDA/saber/src/bump/type_samp.F90#L1960) | fill diagnostics missing values |
