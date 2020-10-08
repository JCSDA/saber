# Module tools_asa007

| Type | Name | Purpose | Arguments |     | Type | Intent |
| :--: | :--: | :------ | ----: | :-------- | :--: | :----: |
| subroutine | [asa007_cholesky](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_asa007.F90#L41) | compute cholesky decomposition | **mpl**<br>**n**<br>**nn**<br>**a**<br>**u**<br>**ierr** |  MPI data<br> Matrix rank<br> Half-matrix size (n*(n-1)/2)<br> Matrix<br> Matrix square-root<br> Error status | type(mpl_type)<br>integer<br>integer<br>real(kind_real)<br>real(kind_real)<br>integer | inout<br>in<br>in<br>in<br>out<br>out |
| subroutine | [asa007_syminv](https://github.com/JCSDA/saber/tree/develop/src/saber/external/tools_asa007.F90#L120) | compute inverse of a symmetric matrix | **mpl**<br>**n**<br>**nn**<br>**a**<br>**c**<br>**ierr** |  MPI data<br> Matrix rank<br> Half-matrix size (n*(n-1)/2)<br> Matrix<br> Matrix inverse<br> Error status | type(mpl_type)<br>integer<br>integer<br>real(kind_real)<br>real(kind_real)<br>integer | inout<br>in<br>in<br>in<br>out<br>out |
