# Module type_mom_blk

| Type | Name | Purpose | Arguments |     | Type | Intent |
| :--: | :--: | :------ | --------: | :-- | :--: | :----: |
| subroutine | [mom_blk_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mom_blk.F90#L52) | allocation | **mom_blk**<br>**nc1**<br>**geom**<br>**bpar**<br>**ne**<br>**nsub** |  Moments block<br> Subsampling size<br> Geometry<br> Block parameters<br> Ensemble size<br> Number of sub-ensembles | class(mom_blk_type)<br>integer<br>type(geom_type)<br>type(bpar_type)<br>integer<br>integer | inout<br>in<br>in<br>in<br>in<br>in |
| subroutine | [mom_blk_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mom_blk.F90#L81) | release memory | **mom_blk** |  Moments block | class(mom_blk_type) | inout |
| subroutine | [mom_blk_ext](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mom_blk.F90#L105) | halo extension | **mom_blk_out**<br>**mpl**<br>**geom**<br>**bpar**<br>**samp**<br>**mom_blk_in** |  Extended moments block<br> MPI data<br> Geometry<br> Block parameters<br> Sampling<br> Reduced moments block | class(mom_blk_type)<br>type(mpl_type)<br>type(geom_type)<br>type(bpar_type)<br>type(samp_type)<br>type(mom_blk_type) | inout<br>inout<br>in<br>in<br>in<br>in |
