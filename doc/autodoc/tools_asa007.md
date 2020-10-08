# Module tools_asa007

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [asa007_cholesky](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_asa007.F90#L41) | compute cholesky decomposition | **mpl** :  MPI data - type(mpl_type) - inout<br>**n** :  Matrix rank - integer - in<br>**nn** :  Half-matrix size (n*(n-1)/2) - integer - in<br>**a** :  Matrix - real(kind_real) - in<br>**u** :  Matrix square-root - real(kind_real) - out<br>**ierr** :  Error status - integer - out |
| subroutine | [asa007_syminv](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_asa007.F90#L120) | compute inverse of a symmetric matrix | **mpl** :  MPI data - type(mpl_type) - inout<br>**n** :  Matrix rank - integer - in<br>**nn** :  Half-matrix size (n*(n-1)/2) - integer - in<br>**a** :  Matrix - real(kind_real) - in<br>**c** :  Matrix inverse - real(kind_real) - out<br>**ierr** :  Error status - integer - out |
