# Module type_mpl

| Type | Name | Purpose |
| :--: | :--: | :---------- |
| subroutine | [mpl_newunit](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L123) | find a free unit |
| subroutine | [mpl_init](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L154) | initialize MPL object |
| subroutine | [mpl_final](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L210) | finalize MPI |
| subroutine | [mpl_flush](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L226) | flush listings |
| subroutine | [mpl_abort](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L341) | clean MPI abort |
| subroutine | [mpl_warning](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L363) | print warning message |
| subroutine | [mpl_update_tag](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L383) | update MPI tag |
| subroutine | [mpl_broadcast_string_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L404) | broadcast 1d string array |
| subroutine | [mpl_dot_prod_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L427) | global dot product over local fields, 1d |
| subroutine | [mpl_dot_prod_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L458) | global dot product over local fields, 2d |
| subroutine | [mpl_dot_prod_3d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L489) | global dot product over local fields, 3d |
| subroutine | [mpl_dot_prod_4d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L520) | global dot product over local fields, 4d |
| subroutine | [mpl_split_loop](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L551) | split loop over different MPI tasks |
| subroutine | [mpl_share_integer_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L583) | share integer array over different MPI tasks, 1d |
| subroutine | [mpl_share_integer_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L620) | share integer array over different MPI tasks, 2d |
| subroutine | [mpl_share_real_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L678) | share real array over different MPI tasks, 1d |
| subroutine | [mpl_share_real_4d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L715) | share real array over different MPI tasks, 4d |
| subroutine | [mpl_share_logical_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L783) | share logical array over different MPI tasks, 1d |
| subroutine | [mpl_share_logical_3d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L840) | share logical array over different MPI tasks, 3d |
| subroutine | [mpl_share_logical_4d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L914) | share logical array over different MPI tasks, 4d |
| subroutine | [mpl_glb_to_loc_index](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L993) | communicate global index to local index |
| subroutine | [mpl_glb_to_loc_integer_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1064) | global to local, 1d array |
| subroutine | [mpl_glb_to_loc_integer_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1141) | global to local, 2d array |
| subroutine | [mpl_glb_to_loc_real_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1251) | global to local, 1d array |
| subroutine | [mpl_glb_to_loc_real_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1329) | global to local, 2d array |
| subroutine | [mpl_glb_to_loc_logical_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1440) | global to local, 1d array |
| subroutine | [mpl_glb_to_loc_logical_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1518) | global to local, 2d array |
| subroutine | [mpl_loc_to_glb_integer_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1629) | local to global, 1d array |
| subroutine | [mpl_loc_to_glb_integer_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1715) | local to global, 2d array |
| subroutine | [mpl_loc_to_glb_real_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1816) | local to global, 1d array |
| subroutine | [mpl_loc_to_glb_real_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L1903) | local to global, 2d array |
| subroutine | [mpl_loc_to_glb_logical_1d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2005) | local to global, 1d array |
| subroutine | [mpl_loc_to_glb_logical_2d](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2092) | local to global for a logical, 2d array |
| subroutine | [mpl_prog_init](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2194) | initialize progression display |
| subroutine | [mpl_prog_print](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2227) | print progression display |
| subroutine | [mpl_prog_final](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2264) | finalize progression display |
| subroutine | [mpl_ncdimcheck](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2307) | check if NetCDF file dimension exists and has the right size |
| subroutine | [mpl_ncerr](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2375) | handle NetCDF error |
| subroutine | [mpl_write_integer](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2393) | write integer into a log file or into a NetCDF file |
| subroutine | [mpl_write_integer_array](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2432) | write integer array into a log file or into a NetCDF file |
| subroutine | [mpl_write_real](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2484) | write real into a log file or into a NetCDF file |
| subroutine | [mpl_write_real_array](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2518) | write real array into a log file or into a NetCDF file |
| subroutine | [mpl_write_logical](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2566) | write logical into a log file or into a NetCDF file |
| subroutine | [mpl_write_logical_array](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2601) | write logical array into a log file or into a NetCDF file |
| subroutine | [mpl_write_string](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2655) | write string into a log file or into a NetCDF file |
| subroutine | [mpl_write_string_array](https://github.com/JCSDA/saber/tree/develop/src/saber/util/type_mpl.F90#L2691) | write string array into a log file or into a NetCDF file |
