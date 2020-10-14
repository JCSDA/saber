# Module type_gaugrid

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [create_gaugrid](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/type_gaugrid.F90#L58) | Create Gaussian grid | <b>self</b> :  Gaussian grid - class(gaussian_grid) - inout |
| subroutine | [delete_gaugrid](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/type_gaugrid.F90#L75) | Delete Gaussian grid | <b>self</b> :  Gaussian grid - class(gaussian_grid) - inout |
| subroutine | [gaugrid_alloc_coord](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/type_gaugrid.F90#L91) | allocate Gaussian grid coordinate | <b>self</b> :  Gaussian grid - class(gaussian_grid) - inout |
| subroutine | [gaugrid_dealloc_coord](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/type_gaugrid.F90#L110) | deallocate Gaussian grid coordinate | <b>self</b> :  Gaussian grid - class(gaussian_grid) - inout |
| subroutine | [gaugrid_alloc_field](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/type_gaugrid.F90#L128) | allocate Gaussian grid field | <b>self</b> :  Gaussian grid - class(gaussian_grid) - inout |
| subroutine | [gaugrid_dealloc_field](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/type_gaugrid.F90#L144) | deallocate Gaussian grid field | <b>self</b> :  Gaussian grid - class(gaussian_grid) - inout |
| subroutine | [gaugrid_calc_glb_latlon](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/type_gaugrid.F90#L159) | calculate global Gaussian latitudes and longitudes | <b>self</b> :  Gaussian grid - class(gaussian_grid) - inout |
| subroutine | [gaugrid_fld3d_pointer](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/type_gaugrid.F90#L198) | Set 3D field pointer | <b>self</b> :  Gaussian grid - class(gaussian_grid) - inout<br><b>iv</b> :  Variable index - integer - in<br><b>var</b> :  Variable name - character(len=*) - in<br><b>fldpointer</b> :  Field pointer - real(kind_real) - inout |
| subroutine | [gaugrid_fld2d_pointer](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/type_gaugrid.F90#L217) | Set 2D field pointer | <b>self</b> :  Gaussian grid - class(gaussian_grid) - inout<br><b>iv</b> :  Variable index - integer - in<br><b>var</b> :  Variable name - character(len=*) - in<br><b>fldpointer</b> :  Field pointer - real(kind_real) - inout |
| subroutine | [equals](https://github.com/JCSDA/saber/tree/develop/src/saber/gaugrid/type_gaugrid.F90#L234) | create new gaussian grid from other | <b>self</b> :  Gaussian grid - class(gaussian_grid) - inout<br><b>rhs</b> :  Other Gaussian grid - type (gaussian_grid) - in |
| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [dummy](https://github.com/JCSDA/saber/tree/develop/src/saber/interpolation/type_gaugrid.F90#L668) | dummy finalization | <b>bump</b> :  BUMP - type(bump_interpolator) - inout |
