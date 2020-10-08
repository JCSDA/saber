# Module type_bpar

| Type | Name | Purpose | Arguments |     | Type | Intent |
| :--: | :--: | :------ | --------: | :-- | :--: | :----: |
| subroutine | [bpar_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bpar.F90#L57) | allocation | **bpar**<br>**nam**<br>**geom** |  Block parameters<br> Namelist<br> Geometry | class(bpar_type)<br>type(nam_type)<br>type(geom_type) | inout<br>in<br>in |
| subroutine | [bpar_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bpar.F90#L100) | initialization | **bpar**<br>**mpl**<br>**nam**<br>**geom** |  Block parameters<br> MPI data<br> Namelist<br> Geometry | class(bpar_type)<br>type(mpl_type)<br>type(nam_type)<br>type(geom_type) | inout<br>inout<br>in<br>in |
| subroutine | [bpar_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bpar.F90#L291) | release memory | **bpar** |  Block parameters | class(bpar_type) | inout |
