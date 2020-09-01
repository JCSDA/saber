# Module type_mpl

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [mpl_newunit](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L115) | find a free unit |
| subroutine | [mpl_init](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L146) | initialize MPL object |
| subroutine | [mpl_final](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L202) | finalize MPI |
| subroutine | [mpl_flush](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L218) | flush listings |
| subroutine | [mpl_abort](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L333) | clean MPI abort |
| subroutine | [mpl_warning](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L355) | print warning message |
| subroutine | [mpl_update_tag](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L375) | update MPI tag |
| subroutine | [mpl_broadcast_string_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L396) | broadcast 1d string array |
| subroutine | [mpl_dot_prod_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L419) | global dot product over local fields, 1d |
| subroutine | [mpl_dot_prod_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L450) | global dot product over local fields, 2d |
| subroutine | [mpl_dot_prod_3d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L481) | global dot product over local fields, 3d |
| subroutine | [mpl_dot_prod_4d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L512) | global dot product over local fields, 4d |
| subroutine | [mpl_split_loop](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L543) | split loop over different MPI tasks |
| subroutine | [mpl_share_integer_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L575) | share integer array over different MPI tasks, 1d |
| subroutine | [mpl_share_integer_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L612) | share integer array over different MPI tasks, 2d |
| subroutine | [mpl_share_real_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L670) | share real array over different MPI tasks, 1d |
| subroutine | [mpl_share_real_4d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L707) | share real array over different MPI tasks, 4d |
| subroutine | [mpl_share_logical_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L775) | share logical array over different MPI tasks, 1d |
| subroutine | [mpl_share_logical_3d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L832) | share logical array over different MPI tasks, 3d |
| subroutine | [mpl_share_logical_4d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L906) | share logical array over different MPI tasks, 4d |
| subroutine | [mpl_glb_to_loc_index](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L985) | communicate global index to local index |
| subroutine | [mpl_glb_to_loc_real_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1055) | global to local, 1d array |
| subroutine | [mpl_glb_to_loc_real_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1117) | global to local, 2d array |
| subroutine | [mpl_loc_to_glb_real_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1214) | local to global, 1d array |
| subroutine | [mpl_loc_to_glb_real_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1280) | local to global, 2d array |
| subroutine | [mpl_loc_to_glb_logical_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1363) | local to global, 1d array |
| subroutine | [mpl_loc_to_glb_logical_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1429) | local to global for a logical, 2d array |
| subroutine | [mpl_prog_init](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1512) | initialize progression display |
| subroutine | [mpl_prog_print](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1545) | print progression display |
| subroutine | [mpl_prog_final](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1582) | finalize progression display |
| subroutine | [mpl_ncdimcheck](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1625) | check if NetCDF file dimension exists and has the right size |
| subroutine | [mpl_ncerr](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1693) | handle NetCDF error |
| subroutine | [mpl_write_integer](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1711) | write integer into a log file or into a NetCDF file |
| subroutine | [mpl_write_integer_array](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1750) | write integer array into a log file or into a NetCDF file |
| subroutine | [mpl_write_real](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1802) | write real into a log file or into a NetCDF file |
| subroutine | [mpl_write_real_array](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1836) | write real array into a log file or into a NetCDF file |
| subroutine | [mpl_write_logical](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1884) | write logical into a log file or into a NetCDF file |
| subroutine | [mpl_write_logical_array](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1919) | write logical array into a log file or into a NetCDF file |
| subroutine | [mpl_write_string](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1973) | write string into a log file or into a NetCDF file |
| subroutine | [mpl_write_string_array](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2009) | write string array into a log file or into a NetCDF file |
