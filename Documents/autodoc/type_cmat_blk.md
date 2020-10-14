# Module type_cmat_blk

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [cmat_blk_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L72) | allocation | <b>cmat_blk</b> :  C matrix data block - class(cmat_blk_type) - inout<br><b>nam</b> :  Namelist - type(nam_type) - in<br><b>geom</b> :  Geometry - type(geom_type) - in<br><b>bpar</b> :  Block parameters - type(bpar_type) - in |
| subroutine | [cmat_blk_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L111) | initialization | <b>cmat_blk</b> :  C matrix data block - class(cmat_blk_type) - inout<br><b>mpl</b> :  MPI data - type(mpl_type) - in<br><b>nam</b> :  Namelist - type(nam_type) - in<br><b>bpar</b> :  Block parameters - type(bpar_type) - in |
| subroutine | [cmat_blk_partial_bump_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L147) | release memory (partial) | <b>cmat_blk</b> :  C matrix data block - class(cmat_blk_type) - inout |
| subroutine | [cmat_blk_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L171) | release memory (partial) | <b>cmat_blk</b> :  C matrix data block - class(cmat_blk_type) - inout |
| subroutine | [cmat_blk_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L193) | release memory | <b>cmat_blk</b> :  C matrix data block - class(cmat_blk_type) - inout |
