# Module type_nam

| Type | Name | Purpose | Arguments |     | Type | Intent |
| :--: | :--: | :------ | --------: | :-- | :--: | :----: |
| subroutine | [nam_init](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L218) | intialize | **nam**<br>**nproc** |  Namelist<br> Number of MPI task | class(nam_type)<br>integer | out<br>in |
| subroutine | [nam_read](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L415) | read | **nam**<br>**mpl**<br>**namelname** |  Namelist<br> MPI data<br> Namelist name | class(nam_type)<br>type(mpl_type)<br>character(len=*) | inout<br>inout<br>in |
| subroutine | [nam_read_yaml](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1097) | read YAML file | **nam**<br>**mpl**<br>**yamlname** |  Namelist<br> MPI data<br> YAML name | class(nam_type)<br>type(mpl_type)<br>character(len=*) | inout<br>inout<br>inout |
| subroutine | [nam_bcast](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1122) | broadcast | **nam**<br>**mpl** |  Namelist<br> MPI data | class(nam_type)<br>type(mpl_type) | inout<br>inout |
| subroutine | [nam_from_conf](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1299) | intialize from configuration | **nam**<br>**conf** |  Namelist<br> Configuration | class(nam_type)<br>type(fckit_configuration) | inout<br>in |
| subroutine | [nam_check](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1570) | check namelist parameters | **nam**<br>**mpl** |  Namelist<br> MPI data | class(nam_type)<br>type(mpl_type) | inout<br>inout |
| subroutine | [nam_write](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L1939) | write namelist parameters into a log file or into a NetCDF file | **nam**<br>**mpl**<br>**ncid** |  Namelist<br> MPI data<br> NetCDF file | class(nam_type)<br>type(mpl_type)<br>integer | in<br>inout<br>in |
| subroutine | [nam_io_key_value](https://github.com/JCSDA/saber/tree/develop/src/saber/bump/type_nam.F90#L2187) | get I/O value from key | **nam**<br>**io_key**<br>**io_value** |  Namelist<br> I/O key<br> I/O value | class(nam_type)<br>character(len=*)<br>character(len=1024) | in<br>in<br>out |
