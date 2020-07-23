# Module type_mpl

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [mpl%] [newunit](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L116) | find a free unit |
| subroutine | [mpl%] [init](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L147) | initialize MPL object |
| subroutine | [mpl%] [final](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L206) | finalize MPI |
| subroutine | [mpl%] [flush](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L222) | flush listings |
| subroutine | [mpl%] [abort](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L337) | clean MPI abort |
| subroutine | [mpl%] [warning](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L359) | print warning message |
| subroutine | [mpl%] [update_tag](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L379) | update MPI tag |
| subroutine | [mpl%] [broadcast_string_1d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L400) | broadcast 1d string array |
| subroutine | [mpl%] [dot_prod_1d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L423) | global dot product over local fields, 1d |
| subroutine | [mpl%] [dot_prod_2d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L454) | global dot product over local fields, 2d |
| subroutine | [mpl%] [dot_prod_3d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L485) | global dot product over local fields, 3d |
| subroutine | [mpl%] [dot_prod_4d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L516) | global dot product over local fields, 4d |
| subroutine | [mpl%] [split_loop](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L547) | split loop over different MPI tasks |
| subroutine | [mpl%] [share_integer_1d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L579) | share integer array over different MPI tasks, 1d |
| subroutine | [mpl%] [share_integer_2d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L616) | share integer array over different MPI tasks, 2d |
| subroutine | [mpl%] [share_real_1d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L674) | share real array over different MPI tasks, 1d |
| subroutine | [mpl%] [share_real_4d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L711) | share real array over different MPI tasks, 4d |
| subroutine | [mpl%] [share_logical_1d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L779) | share logical array over different MPI tasks, 1d |
| subroutine | [mpl%] [share_logical_3d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L836) | share logical array over different MPI tasks, 3d |
| subroutine | [mpl%] [share_logical_4d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L910) | share logical array over different MPI tasks, 4d |
| subroutine | [mpl%] [glb_to_loc_index](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L989) | communicate global index to local index |
| subroutine | [mpl%] [glb_to_loc_real_1d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1059) | global to local, 1d array |
| subroutine | [mpl%] [glb_to_loc_real_2d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1121) | global to local, 2d array |
| subroutine | [mpl%] [loc_to_glb_real_1d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1200) | local to global, 1d array |
| subroutine | [mpl%] [loc_to_glb_real_2d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1266) | local to global, 2d array |
| subroutine | [mpl%] [loc_to_glb_logical_1d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1349) | local to global, 1d array |
| subroutine | [mpl%] [loc_to_glb_logical_2d](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1415) | local to global for a logical, 2d array |
| subroutine | [mpl%] [prog_init](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1498) | initialize progression display |
| subroutine | [mpl%] [prog_print](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1531) | print progression display |
| subroutine | [mpl%] [prog_final](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1568) | finalize progression display |
| subroutine | [mpl%] [ncdimcheck](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1611) | check if NetCDF file dimension exists and has the right size |
| subroutine | [mpl%] [ncerr](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1679) | handle NetCDF error |
| subroutine | [mpl%] [write_integer](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1697) | write integer into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_integer_array](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1736) | write integer array into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_real](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1788) | write real into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_real_array](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1822) | write real array into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_logical](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1870) | write logical into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_logical_array](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1905) | write logical array into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_string](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1959) | write string into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_string_array](https://github.com/JCSDA/saber/src/saber/bump/type_mpl.F90#L1995) | write string array into a log file or into a NetCDF file |
