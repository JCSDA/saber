# Module type_hdiag

| Type | Name | Purpose | Arguments |     | Type | Intent |
| :--: | :--: | :------ | ----: | :-------- | :--: | :----: |
| subroutine | [hdiag_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_hdiag.F90#L58) | release memory | **hdiag** |  Hybrid diagnostics | class(hdiag_type) | inout |
| subroutine | [hdiag_run_hdiag](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_hdiag.F90#L94) | HDIAG driver | **hdiag**<br>**mpl**<br>**rng**<br>**nam**<br>**geom**<br>**bpar**<br>**io**<br>**ens1**<br>**ens2** |  Hybrid diagnostics<br> MPI data<br> Random number generator<br> Namelist<br> Geometry<br> Block parameters<br> I/O<br> Ensemble 1<br> Ensemble 2 | class(hdiag_type)<br>type(mpl_type)<br>type(rng_type)<br>type(nam_type)<br>type(geom_type)<br>type(bpar_type)<br>type(io_type)<br>type(ens_type)<br>type(ens_type) | inout<br>inout<br>inout<br>inout<br>in<br>in<br>in<br>inout<br>in |
