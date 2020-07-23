# Module type_linop

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [linop%] [alloc](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L64) | allocation |
| subroutine | [linop%] [dealloc](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L94) | release memory |
| subroutine | [linop%] [copy](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L113) | copy |
| subroutine | [linop%] [read](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L155) | read |
| subroutine | [linop%] [write](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L211) | write |
| subroutine | [linop%] [receive](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L269) | linop_receive |
| subroutine | [linop%] [send](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L364) | linop_send |
| subroutine | [linop%] [reorder](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L455) | reorder linear operator |
| subroutine | [linop%] [apply](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L511) | apply linear operator |
| subroutine | [linop%] [apply_ad](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L601) | apply linear operator, adjoint |
| subroutine | [linop%] [apply_sym](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L657) | apply linear operator, symmetric |
| subroutine | [linop%] [add_op](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L724) | add operation |
| subroutine | [linop%] [gather](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L769) | gather data from OpenMP threads |
| subroutine | [linop%] [interp](https://github.com/JCSDA/saber/src/saber/bump/type_linop.F90#L808) | compute horizontal interpolation |
