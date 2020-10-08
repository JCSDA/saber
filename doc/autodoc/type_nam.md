# Module type_nam

| Type | Name | Purpose | Arguments          |
| :--: | :--: | :------ | :----------------- |
| subroutine | [nam_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L218) | intialize | **nam** :  Namelist - class(nam_type) - out<br>**nproc** :  Number of MPI task - integer - in |
| subroutine | [nam_read](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L415) | read | **nam** :  Namelist - class(nam_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**namelname** :  Namelist name - character(len=*) - in |
| subroutine | [nam_read_yaml](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1097) | read YAML file | **nam** :  Namelist - class(nam_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**yamlname** :  YAML name - character(len=*) - inout |
| subroutine | [nam_bcast](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1122) | broadcast | **nam** :  Namelist - class(nam_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout |
| subroutine | [nam_from_conf](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1299) | intialize from configuration | **nam** :  Namelist - class(nam_type) - inout<br>**conf** :  Configuration - type(fckit_configuration) - in |
| subroutine | [nam_check](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1570) | check namelist parameters | **nam** :  Namelist - class(nam_type) - inout<br>**mpl** :  MPI data - type(mpl_type) - inout |
| subroutine | [nam_write](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1939) | write namelist parameters into a log file or into a NetCDF file | **nam** :  Namelist - class(nam_type) - in<br>**mpl** :  MPI data - type(mpl_type) - inout<br>**ncid** :  NetCDF file - integer - in |
| subroutine | [nam_io_key_value](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L2187) | get I/O value from key | **nam** :  Namelist - class(nam_type) - in<br>**io_key** :  I/O key - character(len=*) - in<br>**io_value** :  I/O value - character(len=1024) - out |
