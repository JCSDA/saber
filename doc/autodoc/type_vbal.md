# Module type_vbal

| Type | Name | Purpose | Arguments |     | Type | Intent |
| :--: | :--: | :------ | ----: | :-------- | :--: | :----: |
| subroutine | [vbal_alloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L69) | allocation | **vbal**<br>**nam**<br>**geom**<br>**bpar** |  Vertical balance<br> Namelist<br> Geometry<br> Block parameters | class(vbal_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type) | inout<br>in<br>in<br>in |
| subroutine | [vbal_partial_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L101) | release memory (partial) | **vbal** |  Vertical balance | class(vbal_type) | inout |
| subroutine | [vbal_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L126) | release memory (full) | **vbal** |  Vertical balance | class(vbal_type) | inout |
| subroutine | [vbal_read](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L162) | read | **vbal**<br>**mpl**<br>**nam**<br>**geom**<br>**bpar** |  Vertical balance<br> MPI data<br> Namelist<br> Geometry<br> Block parameters | class(vbal_type)<br>type(mpl_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type) | inout<br>inout<br>in<br>in<br>in |
| subroutine | [vbal_write](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L233) | write | **vbal**<br>**mpl**<br>**nam**<br>**geom**<br>**bpar** |  Vertical balance<br> MPI data<br> Namelist<br> Geometry<br> Block parameters | class(vbal_type)<br>type(mpl_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type) | inout<br>inout<br>in<br>in<br>in |
| subroutine | [vbal_run_vbal](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L316) | compute vertical balance | **vbal**<br>**mpl**<br>**rng**<br>**nam**<br>**geom**<br>**bpar**<br>**ens**<br>**ensu** |  Vertical balance<br> MPI data<br> Random number generator<br> Namelist<br> Geometry<br> Block parameters<br> Ensemble<br> Unbalanced ensemble | class(vbal_type)<br>type(mpl_type)<br>type(rng_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type)<br>type(ens_type)<br>type(ens_type) | inout<br>inout<br>inout<br>inout<br>in<br>in<br>in<br>inout |
| subroutine | [vbal_run_vbal_tests](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L432) | compute vertical balance tests | **vbal**<br>**mpl**<br>**rng**<br>**nam**<br>**geom**<br>**bpar**<br>**io** |  Vertical balance<br> MPI data<br> Random number generator<br> Namelist<br> Geometry<br> Block parameters<br> I/O | class(vbal_type)<br>type(mpl_type)<br>type(rng_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type)<br>type(io_type) | inout<br>inout<br>inout<br>inout<br>in<br>in<br>in |
| subroutine | [vbal_apply](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L464) | apply vertical balance | **vbal**<br>**nam**<br>**geom**<br>**bpar**<br>**fld(geom%nc0a,geom%nl0,nam%nv)** |  Vertical balance<br> Namelist<br> Geometry<br> Block parameters<br> Source/destination vector | class(vbal_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type)<br>real(kind_real) | in<br>in<br>in<br>in<br>inout |
| subroutine | [vbal_apply_inv](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L502) | apply inverse vertical balance | **vbal**<br>**nam**<br>**geom**<br>**bpar**<br>**fld(geom%nc0a,geom%nl0,nam%nv)** |  Vertical balance<br> Namelist<br> Geometry<br> Block parameters<br> Source/destination vector | class(vbal_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type)<br>real(kind_real) | in<br>in<br>in<br>in<br>inout |
| subroutine | [vbal_apply_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L540) | apply adjoint vertical balance | **vbal**<br>**nam**<br>**geom**<br>**bpar**<br>**fld(geom%nc0a,geom%nl0,nam%nv)** |  Vertical balance<br> Namelist<br> Geometry<br> Block parameters<br> Source/destination vector | class(vbal_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type)<br>real(kind_real) | in<br>in<br>in<br>in<br>inout |
| subroutine | [vbal_apply_inv_ad](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L578) | apply inverse adjoint vertical balance | **vbal**<br>**nam**<br>**geom**<br>**bpar**<br>**fld(geom%nc0a,geom%nl0,nam%nv)** |  Vertical balance<br> Namelist<br> Geometry<br> Block parameters<br> Source/destination vector | class(vbal_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type)<br>real(kind_real) | in<br>in<br>in<br>in<br>inout |
| subroutine | [vbal_test_inverse](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L617) | test vertical balance inverse | **vbal**<br>**mpl**<br>**rng**<br>**nam**<br>**geom**<br>**bpar** |  Vertical balance<br> MPI data<br> Random number generator<br> Namelist<br> Geometry<br> Block parameters | class(vbal_type)<br>type(mpl_type)<br>type(rng_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type) | in<br>inout<br>inout<br>in<br>in<br>in |
| subroutine | [vbal_test_adjoint](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L678) | test vertical balance adjoint | **vbal**<br>**mpl**<br>**rng**<br>**nam**<br>**geom**<br>**bpar** |  Vertical balance<br> MPI data<br> Random number generator<br> Namelist<br> Geometry<br> Block parameters | class(vbal_type)<br>type(mpl_type)<br>type(rng_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type) | in<br>inout<br>inout<br>in<br>in<br>in |
| subroutine | [vbal_test_dirac](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_vbal.F90#L749) | apply vertical balance to diracs | **vbal**<br>**mpl**<br>**nam**<br>**geom**<br>**bpar**<br>**io** |  Vertical balance<br> MPI data<br> Namelist<br> Geometry<br> Block parameters<br> I/O | class(vbal_type)<br>type(mpl_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type)<br>type(io_type) | in<br>inout<br>in<br>in<br>in<br>in |