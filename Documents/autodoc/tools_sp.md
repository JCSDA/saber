# Module tools_sp

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [splat](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_sp.F90#L39) | compute latitude functions | <b>IDRT</b> :  Grid identifier - integer - in <br><b>JMAX</b> :  Number of latitudes - integer - in <br><b>SLAT</b> :  Sines of latitude  - real(kind=kind_real) - out<br><b>WLAT</b> :  Gaussian weights- - real(kind=kind_real) - out |
| subroutine | [lubksb](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_sp.F90#L210) | solves a system of linear equations, follows call to LUDCMP | <b>NP</b> :  ? - integer - in<br><b>N</b> :  ? - integer - in<br><b>A</b> :  ? - real(kind=kind_real) - in<br><b>B</b> :  ? - real(kind=kind_real) - inout<br><b>INDX</b> :  ? - integer - in |
| subroutine | [ludcmp](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_sp.F90#L255) | replaces an NxN matrix A with the LU decomposition | <b>N</b> :  ? - integer - in<br><b>NP</b> :  ? - integer - in<br><b>A</b> :  ? - real(kind=kind_real) - inout<br><b>INDX</b> :  ?  - integer - out<br><b>D</b> :  ? - real(kind=kind_real) - out |
