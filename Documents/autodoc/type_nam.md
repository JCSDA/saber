# Module type_nam

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [nam_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L218) | intialize | <b>nam</b> :  Namelist - class(nam_type) - out<br><b>nproc</b> :  Number of MPI task - integer - in |
| subroutine | [nam_read](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L415) | read | <b>nam</b> :  Namelist - class(nam_type) - inout<br><b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>namelname</b> :  Namelist name - character(len=*) - in |
| subroutine | [nam_read_yaml](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1097) | read YAML file | <b>nam</b> :  Namelist - class(nam_type) - inout<br><b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>yamlname</b> :  YAML name - character(len=*) - inout |
| subroutine | [nam_bcast](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1122) | broadcast | <b>nam</b> :  Namelist - class(nam_type) - inout<br><b>mpl</b> :  MPI data - type(mpl_type) - inout |
| subroutine | [nam_from_conf](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1299) | intialize from configuration | <b>nam</b> :  Namelist - class(nam_type) - inout<br><b>conf</b> :  Configuration - type(fckit_configuration) - in |
| subroutine | [nam_check](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1570) | check namelist parameters | <b>nam</b> :  Namelist - class(nam_type) - inout<br><b>mpl</b> :  MPI data - type(mpl_type) - inout |
| subroutine | [nam_write](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1939) | write namelist parameters into a log file or into a NetCDF file | <b>nam</b> :  Namelist - class(nam_type) - in<br><b>mpl</b> :  MPI data - type(mpl_type) - inout<br><b>ncid</b> :  NetCDF file - integer - in |
| subroutine | [nam_io_key_value](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L2187) | get I/O value from key | <b>nam</b> :  Namelist - class(nam_type) - in<br><b>io_key</b> :  I/O key - character(len=*) - in<br><b>io_value</b> :  I/O value - character(len=1024) - out |
