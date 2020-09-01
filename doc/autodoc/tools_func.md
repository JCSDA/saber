# Module tools_func

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [lonlatmod](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L32) | set latitude between -pi/2 and pi/2 and longitude between -pi and pi |
| subroutine | [sphere_dist](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L62) | compute the great-circle distance between two points |
| subroutine | [reduce_arc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L88) | reduce arc to a given distance |
| subroutine | [lonlat2xyz](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L125) | convert longitude/latitude to cartesian coordinates |
| subroutine | [xyz2lonlat](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L164) | convert longitude/latitude to cartesian coordinates |
| subroutine | [vector_product](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L201) | compute normalized vector product |
| subroutine | [vector_triple_product](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L228) | compute vector triple product |
| subroutine | [add](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L266) | check if value missing and add if not missing |
| subroutine | [divide](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L296) | check if value missing and divide if not missing |
| subroutine | [fit_diag](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L318) | compute diagnostic fit function |
| function | [gc99](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L497) | Gaspari and Cohn (1999) function, with the support radius as a parameter |
| function | [fit_func](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L520) | fit_function |
| subroutine | [fit_lct](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L547) | LCT fit |
| subroutine | [lct_d2h](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L618) | from D (Daley tensor) to H (local correlation tensor) |
| subroutine | [lct_h2r](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L660) | from H (local correlation tensor) to support radii |
| subroutine | [lct_r2d](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L711) | from support radius to Daley tensor diagonal element |
| subroutine | [check_cond](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L728) | check tensor conditioning |
| function | [matern](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L769) | compute the normalized diffusion function from eq. (55) of Mirouze and Weaver (2013), for the 3d case (d = 3) |
| subroutine | [cholesky](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L812) | compute cholesky decomposition |
| subroutine | [syminv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L865) | compute inverse of a symmetric matrix |
| subroutine | [histogram](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/tools_func.F90#L918) | compute bins and histogram from a list of values |
