# Module type_mom_blk

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [mom_blk_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mom_blk.F90#L52) | allocation | <b>mom_blk</b> :  Moments block - class(mom_blk_type) - inout<br><b>nc1</b> :  Subsampling size - integer - in<br><b>geom</b> :  Geometry - type(geom_type) - in<br><b>bpar</b> :  Block parameters - type(bpar_type) - in<br><b>ne</b> :  Ensemble size - integer - in<br><b>nsub</b> :  Number of sub-ensembles - integer - in |
| subroutine | [mom_blk_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mom_blk.F90#L81) | release memory | <b>mom_blk</b> :  Moments block - class(mom_blk_type) - inout |
| subroutine | [mom_blk_ext](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_mom_blk.F90#L105) | halo extension | <b>mom_blk_out</b> :  Extended moments block - class(mom_blk_type) - inout<br><b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>geom</b> :  Geometry - type(geom_type) - in<br><b>bpar</b> :  Block parameters - type(bpar_type) - in<br><b>samp</b> :  Sampling - type(samp_type) - in<br><b>mom_blk_in</b> :  Reduced moments block - type(mom_blk_type) - in |
