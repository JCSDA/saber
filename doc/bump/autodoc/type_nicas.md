# Module type_nicas

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [nicas%] [alloc](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L71) | allocation |
| subroutine | [nicas%] [partial_dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L115) | release memory (partial) |
| subroutine | [nicas%] [dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L136) | release memory (full) |
| subroutine | [nicas%] [read](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L161) | read |
| subroutine | [nicas%] [write](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L239) | write |
| subroutine | [nicas%] [write_mpi_summary](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L313) | write MPI related data summary |
| subroutine | [nicas%] [run_nicas](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L396) | NICAS driver |
| subroutine | [nicas%] [run_nicas_tests](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L470) | NICAS tests driver |
| subroutine | [nicas%] [alloc_cv](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L570) | allocation |
| subroutine | [nicas%] [random_cv](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L623) | generate a random control vector |
| subroutine | [nicas%] [apply](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L682) | apply NICAS |
| subroutine | [nicas%] [apply_from_sqrt](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L953) | apply NICAS from square-root |
| subroutine | [nicas%] [apply_sqrt](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L999) | apply NICAS square-root |
| subroutine | [nicas%] [apply_sqrt_ad](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L1207) | apply NICAS square-root, adjoint |
| subroutine | [nicas%] [randomize](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L1438) | randomize NICAS from square-root |
| subroutine | [nicas%] [apply_bens](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L1507) | apply localized ensemble covariance |
| subroutine | [nicas%] [test_adjoint](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L1566) | test NICAS adjoint |
| subroutine | [nicas%] [test_dirac](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L1659) | apply NICAS to diracs |
| subroutine | [nicas%] [test_randomization](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L1723) | test NICAS randomization method with respect to theoretical error statistics |
| subroutine | [nicas%] [test_consistency](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L1854) | test HDIAG-NICAS consistency with a randomization method |
| subroutine | [nicas%] [test_optimality](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L1921) | test HDIAG localization optimality with a randomization method |
| subroutine | [define_test_vectors](https://github.com/JCSDA/saber/src/saber/bump/type_nicas.F90#L2110) | define test vectors |
