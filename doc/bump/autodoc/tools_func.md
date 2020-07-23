# Module tools_func

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [lonlatmod](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L32) | set latitude between -pi/2 and pi/2 and longitude between -pi and pi |
| subroutine | [sphere_dist](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L62) | compute the great-circle distance between two points |
| subroutine | [reduce_arc](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L82) | reduce arc to a given distance |
| subroutine | [lonlat2xyz](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L119) | convert longitude/latitude to cartesian coordinates |
| subroutine | [xyz2lonlat](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L154) | convert longitude/latitude to cartesian coordinates |
| subroutine | [vector_product](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L183) | compute normalized vector product |
| subroutine | [vector_triple_product](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L210) | compute vector triple product |
| subroutine | [add](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L248) | check if value missing and add if not missing |
| subroutine | [divide](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L278) | check if value missing and divide if not missing |
| subroutine | [fit_diag](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L300) | compute diagnostic fit function |
| subroutine | [fit_diag_dble](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L432) | compute diagnostic fit function |
| function | [gc99](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L580) | Gaspari and Cohn (1999) function, with the support radius as a parameter |
| function | [cres](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L603) | reservoir code correlation function, with the support radius as a parameter |
| function | [poly4](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L624) | 4th order polynomial correlation function, with the support radius as a parameter |
| function | [fb07_gc99](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L645) | normalized Furrer and Bengtsson (2007) localization function, with the support radius of the related gc99 function as a parameter |
| function | [fb07_cres](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L669) | normalized Furrer and Bengtsson (2007) localization function, with the support radius of the related cres function as a parameter |
| function | [fit_func](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L693) | fit_function |
| subroutine | [fit_lct](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L741) | LCT fit |
| subroutine | [lct_d2h](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L812) | from D (Daley tensor) to H (local correlation tensor) |
| subroutine | [lct_h2r](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L854) | from H (local correlation tensor) to support radii |
| subroutine | [lct_r2d](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L905) | from support radius to Daley tensor diagonal element |
| subroutine | [check_cond](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L922) | check tensor conditioning |
| function | [matern](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L963) | compute the normalized diffusion function from eq. (55) of Mirouze and Weaver (2013), for the 3d case (d = 3) |
| subroutine | [cholesky](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L1006) | compute cholesky decomposition |
| subroutine | [syminv](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L1058) | compute inverse of a symmetric matrix |
| subroutine | [histogram](https://github.com/JCSDA/saber/src/saber/bump/tools_func.F90#L1110) | compute bins and histogram from a list of values |
