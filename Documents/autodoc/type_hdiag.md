# Module type_hdiag

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [hdiag_dealloc](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_hdiag.F90#L58) | release memory | <b>hdiag</b> :  Hybrid diagnostics - class(hdiag_type) - inout |
| subroutine | [hdiag_run_hdiag](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_hdiag.F90#L94) | HDIAG driver | <b>hdiag</b> :  Hybrid diagnostics - class(hdiag_type) - inout<br><b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>rng</b> :  Random number generator - type(rng_type) - inout<br><b>nam</b> :  Namelist - type(nam_type) - inout<br><b>geom</b> :  Geometry - type(geom_type) - in<br><b>bpar</b> :  Block parameters - type(bpar_type) - in<br><b>io</b> :  I/O - type(io_type) - in<br><b>ens1</b> :  Ensemble 1 - type(ens_type) - inout<br><b>ens2</b> :  Ensemble 2 - type(ens_type) - in |
