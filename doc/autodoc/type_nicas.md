# Module type_nicas

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [nicas_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L70) | allocation |
| subroutine | [nicas_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L126) | release memory (partial) |
| subroutine | [nicas_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L147) | release memory (full) |
| subroutine | [nicas_read](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L172) | read |
| subroutine | [nicas_write](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L254) | write |
| subroutine | [nicas_run_nicas](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L334) | NICAS driver |
| subroutine | [nicas_run_nicas_tests](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L399) | NICAS tests driver |
| subroutine | [nicas_alloc_cv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L499) | allocation |
| subroutine | [nicas_random_cv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L552) | generate a random control vector |
| subroutine | [nicas_apply](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L607) | apply NICAS |
| subroutine | [nicas_apply_from_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L878) | apply NICAS from square-root |
| subroutine | [nicas_apply_sqrt](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L924) | apply NICAS square-root |
| subroutine | [nicas_apply_sqrt_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1134) | apply NICAS square-root, adjoint |
| subroutine | [nicas_randomize](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1367) | randomize NICAS from square-root |
| subroutine | [nicas_apply_bens](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1444) | apply localized ensemble covariance |
| subroutine | [nicas_test_adjoint](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1503) | test NICAS adjoint |
| subroutine | [nicas_test_dirac](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1596) | apply NICAS to diracs |
| subroutine | [nicas_test_randomization](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1659) | test NICAS randomization method with respect to theoretical error statistics |
| subroutine | [nicas_test_consistency](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1790) | test HDIAG-NICAS consistency with a randomization method |
| subroutine | [nicas_test_optimality](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L1926) | test HDIAG localization optimality with a randomization method |
| subroutine | [define_test_vectors](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nicas.F90#L2115) | define test vectors |
