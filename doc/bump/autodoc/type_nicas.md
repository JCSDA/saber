# Module type_nicas

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [nicas%] [alloc](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L70) | allocation |
| subroutine | [nicas%] [partial_dealloc](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L114) | release memory (partial) |
| subroutine | [nicas%] [dealloc](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L135) | release memory (full) |
| subroutine | [nicas%] [read](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L160) | read |
| subroutine | [nicas%] [write](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L315) | write |
| subroutine | [nicas%] [write_mpi_summary](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L432) | write MPI related data summary |
| subroutine | [nicas%] [run_nicas](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L515) | NICAS driver |
| subroutine | [nicas%] [run_nicas_tests](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L589) | NICAS tests driver |
| subroutine | [nicas%] [alloc_cv](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L712) | allocation |
| subroutine | [nicas%] [random_cv](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L765) | generate a random control vector |
| subroutine | [nicas%] [apply](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L824) | apply NICAS |
| subroutine | [nicas%] [apply_from_sqrt](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L1095) | apply NICAS from square-root |
| subroutine | [nicas%] [apply_sqrt](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L1141) | apply NICAS square-root |
| subroutine | [nicas%] [apply_sqrt_ad](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L1349) | apply NICAS square-root, adjoint |
| subroutine | [nicas%] [randomize](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L1580) | randomize NICAS from square-root |
| subroutine | [nicas%] [apply_bens](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L1649) | apply localized ensemble covariance |
| subroutine | [nicas%] [test_adjoint](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L1708) | test NICAS adjoint |
| subroutine | [nicas%] [test_dirac](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L1801) | apply NICAS to diracs |
| subroutine | [nicas%] [test_randomization](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L1865) | test NICAS randomization method with respect to theoretical error statistics |
| subroutine | [nicas%] [test_consistency](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L1986) | test HDIAG-NICAS consistency with a randomization method |
| subroutine | [nicas%] [test_optimality](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L2053) | test HDIAG localization optimality with a randomization method |
| subroutine | [define_test_vectors](https://github.com/JCSDA/saber/src/bump/type_nicas.F90#L2242) | define test vectors |
