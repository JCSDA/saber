# Module type_cmat_blk

| Type | Name | Purpose | Arguments |     | Type | Intent |
| :--: | :--: | :------ | --------: | :-- | :--: | :----: |
| subroutine | [cmat_blk_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L72) | allocation | **cmat_blk**<br>**nam**<br>**geom**<br>**bpar** |  C matrix data block<br> Namelist<br> Geometry<br> Block parameters | class(cmat_blk_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type) | inout<br>in<br>in<br>in |
| subroutine | [cmat_blk_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L111) | initialization | **cmat_blk**<br>**mpl**<br>**nam**<br>**bpar** |  C matrix data block<br> MPI data<br> Namelist<br> Block parameters | class(cmat_blk_type)<br>type(mpl_type)<br>type(nam_type)<br>type(bpar_type) | inout<br>in<br>in<br>in |
| subroutine | [cmat_blk_partial_bump_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L147) | release memory (partial) | **cmat_blk** |  C matrix data block | class(cmat_blk_type) | inout |
| subroutine | [cmat_blk_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L171) | release memory (partial) | **cmat_blk** |  C matrix data block | class(cmat_blk_type) | inout |
| subroutine | [cmat_blk_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_cmat_blk.F90#L193) | release memory | **cmat_blk** |  C matrix data block | class(cmat_blk_type) | inout |
