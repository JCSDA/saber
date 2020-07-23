# Module type_samp

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [samp%] [alloc_mask](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L129) | allocation for mask |
| subroutine | [samp%] [alloc_other](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L151) | allocation for other variables |
| subroutine | [samp%] [partial_dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L186) | release memory (partial) |
| subroutine | [samp%] [dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L211) | release memory |
| subroutine | [samp%] [read](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L276) | read |
| subroutine | [samp%] [write](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L388) | write |
| subroutine | [samp%] [compute_sampling_c1](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L506) | compute sampling, subset Sc1 |
| subroutine | [samp%] [compute_sampling_c2](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L649) | compute sampling, subset Sc2 |
| subroutine | [samp%] [compute_mask](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L738) | compute mask |
| subroutine | [samp%] [compute_sampling_zs](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L897) | compute zero-separation sampling |
| subroutine | [samp%] [compute_sampling_ps](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L988) | compute positive separation sampling |
| subroutine | [samp%] [compute_sampling_lct](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L1124) | compute LCT sampling |
| subroutine | [samp%] [check_mask](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L1190) | check sampling mask |
| subroutine | [samp%] [compute_mpi_a](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L1275) | compute sampling MPI distribution, halo A |
| subroutine | [samp%] [compute_mpi_ab](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L1340) | compute sampling MPI distribution, halos A-B |
| subroutine | [samp%] [compute_mpi_c](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L1546) | compute sampling MPI distribution, halo C |
| subroutine | [samp%] [compute_mpi_d](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L1676) | compute sampling MPI distribution, halo D |
| subroutine | [samp%] [compute_mpi_f](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L1790) | compute sampling MPI distribution, halo F |
| subroutine | [samp%] [diag_filter](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L1857) | filter diagnostics |
| subroutine | [samp%] [diag_fill](https://github.com/JCSDA/saber/src/saber/bump/type_samp.F90#L1990) | fill diagnostics missing values |
