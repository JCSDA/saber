# Module type_obsop

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [obsop%] [partial_dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L66) | release memory (partial) |
| subroutine | [obsop%] [dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L83) | release memory (full) |
| subroutine | [obsop%] [read](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L102) | read observations locations |
| subroutine | [obsop%] [write](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L141) | write observations locations |
| subroutine | [obsop%] [from](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L184) | copy observation operator data |
| subroutine | [obsop%] [run_obsop](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L216) | observation operator driver |
| subroutine | [obsop%] [run_obsop_tests](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L340) | observation operator tests driver |
| subroutine | [obsop%] [apply](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L372) | observation operator interpolation |
| subroutine | [obsop%] [apply_ad](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L405) | observation operator interpolation adjoint |
| subroutine | [obsop%] [test_adjoint](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L441) | test observation operator adjoints accuracy |
| subroutine | [obsop%] [test_accuracy](https://github.com/JCSDA/saber/src/saber/bump/type_obsop.F90#L484) | test observation operator accuracy |
