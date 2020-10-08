# Module type_mom_blk

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [mom_blk_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mom_blk.F90#L52) | allocation | **mom_blk** :  Moments block - class(mom_blk_type) - inout<br>**nc1** :  Subsampling size - integer - in<br>**geom** :  Geometry - type(geom_type) - in<br>**bpar** :  Block parameters - type(bpar_type) - in<br>**ne** :  Ensemble size - integer - in<br>**nsub** :  Number of sub-ensembles - integer - in |
| subroutine | [mom_blk_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mom_blk.F90#L81) | release memory | **mom_blk** :  Moments block - class(mom_blk_type) - inout |
| subroutine | [mom_blk_ext](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mom_blk.F90#L105) | halo extension | **mom_blk_out** :  Extended moments block - class(mom_blk_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**geom** :  Geometry - type(geom_type) - in<br>**bpar** :  Block parameters - type(bpar_type) - in<br>**samp** :  Sampling - type(samp_type) - in<br>**mom_blk_in** :  Reduced moments block - type(mom_blk_type) - in |
