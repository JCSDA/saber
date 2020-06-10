# Module type_nicas

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [nicas_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L71) | allocation |
| subroutine | [nicas_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L124) | release memory (partial) |
| subroutine | [nicas_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L145) | release memory (full) |
| subroutine | [nicas_read](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L170) | read |
| subroutine | [nicas_write](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L248) | write |
| subroutine | [nicas_run_nicas](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L328) | NICAS driver |
| subroutine | [nicas_run_nicas_tests](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L393) | NICAS tests driver |
| subroutine | [nicas_alloc_cv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L493) | allocation |
| subroutine | [nicas_random_cv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L546) | generate a random control vector |
| subroutine | [nicas_apply](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L611) | apply NICAS |
| subroutine | [nicas_apply_from_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L882) | apply NICAS from square-root |
| subroutine | [nicas_apply_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L928) | apply NICAS square-root |
| subroutine | [nicas_apply_sqrt_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1138) | apply NICAS square-root, adjoint |
| subroutine | [nicas_randomize](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1371) | randomize NICAS from square-root |
| subroutine | [nicas_apply_bens](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1448) | apply localized ensemble covariance |
| subroutine | [nicas_test_adjoint](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1507) | test NICAS adjoint |
| subroutine | [nicas_test_dirac](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1600) | apply NICAS to diracs |
| subroutine | [nicas_test_randomization](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1663) | test NICAS randomization method with respect to theoretical error statistics |
| subroutine | [nicas_test_consistency](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1797) | test HDIAG-NICAS consistency with a randomization method |
| subroutine | [nicas_test_optimality](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1924) | test HDIAG localization optimality with a randomization method |
| subroutine | [define_test_vectors](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L2114) | define test vectors |
