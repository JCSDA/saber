# Module type_cmat_blk

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [cmat_blk_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L72) | allocation | **cmat_blk** :  C matrix data block - class(cmat_blk_type) - inout<br>**nam** :  Namelist - type(nam_type) - in<br>**geom** :  Geometry - type(geom_type) - in<br>**bpar** :  Block parameters - type(bpar_type) - in |
| subroutine | [cmat_blk_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L111) | initialization | **cmat_blk** :  C matrix data block - class(cmat_blk_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - in<br>**nam** :  Namelist - type(nam_type) - in<br>**bpar** :  Block parameters - type(bpar_type) - in |
| subroutine | [cmat_blk_partial_bump_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L147) | release memory (partial) | **cmat_blk** :  C matrix data block - class(cmat_blk_type) - inout |
| subroutine | [cmat_blk_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L171) | release memory (partial) | **cmat_blk** :  C matrix data block - class(cmat_blk_type) - inout |
| subroutine | [cmat_blk_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L193) | release memory | **cmat_blk** :  C matrix data block - class(cmat_blk_type) - inout |
