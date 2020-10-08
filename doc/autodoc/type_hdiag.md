# Module type_hdiag

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [hdiag_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_hdiag.F90#L58) | release memory | **hdiag** :  Hybrid diagnostics - class(hdiag_type) - inout |
| subroutine | [hdiag_run_hdiag](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_hdiag.F90#L94) | HDIAG driver | **hdiag** :  Hybrid diagnostics - class(hdiag_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**rng** :  Random number generator - type(rng_type) - inout<br>**nam** :  Namelist - type(nam_type) - inout<br>**geom** :  Geometry - type(geom_type) - in<br>**bpar** :  Block parameters - type(bpar_type) - in<br>**io** :  I/O - type(io_type) - in<br>**ens1** :  Ensemble 1 - type(ens_type) - inout<br>**ens2** :  Ensemble 2 - type(ens_type) - in |
