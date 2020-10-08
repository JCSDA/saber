# Module type_bpar

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [bpar_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bpar.F90#L57) | allocation | **bpar** :  Block parameters - class(bpar_type) - inout<br>**nam** :  Namelist - type(nam_type) - in<br>**geom** :  Geometry - type(geom_type) - in |
| subroutine | [bpar_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bpar.F90#L100) | initialization | **bpar** :  Block parameters - class(bpar_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**nam** :  Namelist - type(nam_type) - in<br>**geom** :  Geometry - type(geom_type) - in |
| subroutine | [bpar_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bpar.F90#L291) | release memory | **bpar** :  Block parameters - class(bpar_type) - inout |
