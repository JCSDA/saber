# Module tools_asa007

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [asa007_cholesky](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_asa007.F90#L41) | compute cholesky decomposition | <b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>n</b> :  Matrix rank - integer - in<br><b>nn</b> :  Half-matrix size (n*(n-1)/2) - integer - in<br><b>a</b> :  Matrix - real(kind_real) - in<br><b>u</b> :  Matrix square-root - real(kind_real) - out<br><b>ierr</b> :  Error status - integer - out |
| subroutine | [asa007_syminv](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_asa007.F90#L120) | compute inverse of a symmetric matrix | <b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>n</b> :  Matrix rank - integer - in<br><b>nn</b> :  Half-matrix size (n*(n-1)/2) - integer - in<br><b>a</b> :  Matrix - real(kind_real) - in<br><b>c</b> :  Matrix inverse - real(kind_real) - out<br><b>ierr</b> :  Error status - integer - out |
