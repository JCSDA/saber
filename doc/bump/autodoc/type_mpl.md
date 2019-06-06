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
| subroutine | [mpl%] [ncerr](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L534) | handle NetCDF error |
| subroutine | [mpl%] [update_tag](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L552) | update MPL tag |
| subroutine | [mpl%] [bcast_string_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L573) | broadcast 1d string array |
| subroutine | [mpl%] [dot_prod_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L596) | global dot product over local fields, 1d |
| subroutine | [mpl%] [dot_prod_2d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L627) | global dot product over local fields, 2d |
| subroutine | [mpl%] [dot_prod_3d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L658) | global dot product over local fields, 3d |
| subroutine | [mpl%] [dot_prod_4d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L689) | global dot product over local fields, 4d |
| subroutine | [mpl%] [split](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L720) | split array over different MPI tasks |
| subroutine | [mpl%] [share_integer_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L752) | share integer array over different MPI tasks, 1d |
| subroutine | [mpl%] [share_integer_2d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L789) | share integer array over different MPI tasks, 2d |
| subroutine | [mpl%] [share_real_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L847) | share real array over different MPI tasks, 1d |
| subroutine | [mpl%] [share_real_4d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L884) | share real array over different MPI tasks, 4d |
| subroutine | [mpl%] [share_logical_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L952) | share logical array over different MPI tasks, 1d |
| subroutine | [mpl%] [share_logical_3d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1009) | share logical array over different MPI tasks, 3d |
| subroutine | [mpl%] [share_logical_4d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1083) | share logical array over different MPI tasks, 4d |
| subroutine | [mpl%] [glb_to_loc_index](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1162) | communicate global index to local index |
| subroutine | [mpl%] [glb_to_loc_real_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1232) | global to local, 1d array |
| subroutine | [mpl%] [glb_to_loc_real_2d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1294) | global to local, 2d array |
| subroutine | [mpl%] [loc_to_glb_real_1d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1373) | local to global, 1d array |
| subroutine | [mpl%] [loc_to_glb_real_2d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1439) | local to global, 2d array |
| subroutine | [mpl%] [loc_to_glb_logical_2d](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1522) | local to global for a logical, 2d array |
| subroutine | [mpl%] [write_integer](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1605) | write integer into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_integer_array](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1643) | write integer array into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_real](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1695) | write real into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_real_array](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1728) | write real array into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_logical](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1776) | write logical into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_logical_array](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1810) | write logical array into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_string](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1864) | write string into a log file or into a NetCDF file |
| subroutine | [mpl%] [write_string_array](https://github.com/JCSDA/saber/src/bump/type_mpl.F90#L1899) | write string array into a log file or into a NetCDF file |
