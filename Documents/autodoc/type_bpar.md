# Module type_bpar

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [bpar_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bpar.F90#L57) | allocation | <b>bpar</b> :  Block parameters - class(bpar_type) - inout<br><b>nam</b> :  Namelist - type(nam_type) - in<br><b>geom</b> :  Geometry - type(geom_type) - in |
| subroutine | [bpar_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bpar.F90#L100) | initialization | <b>bpar</b> :  Block parameters - class(bpar_type) - inout<br><b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>nam</b> :  Namelist - type(nam_type) - in<br><b>geom</b> :  Geometry - type(geom_type) - in |
| subroutine | [bpar_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_bpar.F90#L291) | release memory | <b>bpar</b> :  Block parameters - class(bpar_type) - inout |
