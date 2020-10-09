# Module tools_sp

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [splat](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_sp.F90#L39) | compute latitude functions | **IDRT** :  Grid identifier - integer - in <br>**JMAX** :  Number of latitudes - integer - in <br>**SLAT** :  Sines of latitude  - real(kind=kind_real) - out<br>**WLAT** :  Gaussian weights- - real(kind=kind_real) - out |
| subroutine | [lubksb](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_sp.F90#L210) | solves a system of linear equations, follows call to LUDCMP | **NP** :  ? - integer - in<br>**N** :  ? - integer - in<br>**A** :  ? - real(kind=kind_real) - in<br>**B** :  ? - real(kind=kind_real) - inout<br>**INDX** :  ? - integer - in |
| subroutine | [ludcmp](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_sp.F90#L255) | replaces an NxN matrix A with the LU decomposition | **N** :  ? - integer - in<br>**NP** :  ? - integer - in<br>**A** :  ? - real(kind=kind_real) - inout<br>**INDX** :  ?  - integer - out<br>**D** :  ? - real(kind=kind_real) - out |
