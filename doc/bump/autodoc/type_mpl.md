# Module type_mpl

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [mpl%] [newunit](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L115) | find a free unit |
| subroutine | [mpl%] [init](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L146) | initialize MPL object |
| subroutine | [mpl%] [final](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L201) | finalize MPI |
| subroutine | [mpl%] [flush](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L220) | flush listings |
| subroutine | [mpl%] [abort](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L335) | clean MPI abort |
| subroutine | [mpl%] [warning](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L364) | print warning message |
| subroutine | [mpl%] [prog_init](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L383) | initialize progression display |
| subroutine | [mpl%] [prog_print](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L409) | print progression display |
| subroutine | [mpl%] [prog_final](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L446) | finalize progression display |
| subroutine | [mpl%] [ncdimcheck](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L474) | check if NetCDF file dimension exists and has the right size |
| subroutine | [mpl%] [ncerr](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L542) | handle NetCDF error |
| subroutine | [mpl%] [update_tag](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L560) | update MPL tag |
| subroutine | [mpl%] [bcast_string_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L581) | broadcast 1d string array |
| subroutine | [mpl%] [dot_prod_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L604) | global dot product over local fields, 1d |
| subroutine | [mpl%] [dot_prod_2d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L635) | global dot product over local fields, 2d |
| subroutine | [mpl%] [dot_prod_3d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L666) | global dot product over local fields, 3d |
| subroutine | [mpl%] [dot_prod_4d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L697) | global dot product over local fields, 4d |
| subroutine | [mpl%] [split](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L728) | split array over different MPI tasks |
| subroutine | [mpl%] [share_integer_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L760) | share integer array over different MPI tasks, 1d |
| subroutine | [mpl%] [share_integer_2d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L797) | share integer array over different MPI tasks, 2d |
| subroutine | [mpl%] [share_real_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L855) | share real array over different MPI tasks, 1d |
| subroutine | [mpl%] [share_real_4d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L892) | share real array over different MPI tasks, 4d |
| subroutine | [mpl%] [share_logical_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L960) | share logical array over different MPI tasks, 1d |
| subroutine | [mpl%] [share_logical_3d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1017) | share logical array over different MPI tasks, 3d |
| subroutine | [mpl%] [share_logical_4d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1091) | share logical array over different MPI tasks, 4d |
| subroutine | [mpl%] [glb_to_loc_index](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1170) | communicate global index to local index |
| subroutine | [mpl%] [glb_to_loc_real_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1240) | global to local, 1d array |
| subroutine | [mpl%] [glb_to_loc_real_2d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1302) | global to local, 2d array |
| subroutine | [mpl%] [loc_to_glb_real_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1381) | local to global, 1d array |
| subroutine | [mpl%] [loc_to_glb_real_2d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1447) | local to global, 2d array |
| subroutine | [mpl%] [loc_to_glb_logical_2d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1530) | local to global for a logical, 2d array |
| subroutine | [mpl%] [write_integer](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1613) | write integer into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_integer_array](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1651) | write integer array into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_real](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1703) | write real into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_real_array](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1736) | write real array into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_logical](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1784) | write logical into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_logical_array](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1818) | write logical array into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_string](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1872) | write string into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_string_array](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1907) | write string array into a log file or into a NetCDF file |
