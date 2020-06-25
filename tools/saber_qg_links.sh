#!/bin/bash
#----------------------------------------------------------------------
# Bash script: saber_qg_links
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# dirac_bump_cov
ln -sf ../qg_parameters_bump_cov/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc testdata/qg_dirac_bump_cov/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc
ln -sf ../qg_parameters_bump_cov/test_00_var.nc testdata/qg_dirac_bump_cov/test_00_var.nc

# dirac_bump_lct
ln -sf ../qg_parameters_bump_lct/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc testdata/qg_dirac_bump_lct/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc
ln -sf ../qg_parameters_bump_lct/test_00_var.nc testdata/qg_dirac_bump_lct/test_00_var.nc

# dirac_bump_hyb (hybrid localization is not very good, take localization alone instead)
ln -sf ../qg_parameters_bump_loc_3d/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc testdata/qg_dirac_bump_hyb/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# dirac_bump_loc_3d
ln -sf ../qg_parameters_bump_loc_3d/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc testdata/qg_dirac_bump_loc_3d/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# dirac_bump_loc_4d
ln -sf ../qg_parameters_bump_loc_4d/test_00_nicas-2-sqrt_0001-0001_common.nc testdata/qg_dirac_bump_loc_4d/test_00_nicas-2-sqrt_0001-0001_common.nc

# dirac_id_cov
ln -sf ../qg_parameters_bump_cov/test_00_var.nc testdata/qg_dirac_id_cov/test_00_var.nc

# 3densvar_bump
ln -sf ../qg_parameters_bump_loc_3d/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc testdata/qg_3densvar_bump/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# 3dvar_bump
ln -sf ../qg_parameters_bump_cov/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc testdata/qg_3dvar_bump/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc
ln -sf ../qg_parameters_bump_cov/test_00_var.nc testdata/qg_3dvar_bump/test_00_var.nc

# 3dvar_hybrid_bump (hybrid localization is not very good, take localization alone instead)
ln -sf ../qg_parameters_bump_loc_3d/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc testdata/qg_3dvar_hybrid_bump/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# 4densvar_bump
ln -sf ../qg_parameters_bump_loc_4d/test_00_nicas-2-sqrt_0001-0001_common.nc testdata/qg_4densvar_bump/test_00_nicas-2-sqrt_0001-0001_common.nc

# 4densvar_advect_bump
ln -sf ../qg_parameters_bump_loc_4d/test_00_nicas-2-sqrt_0001-0001_common.nc testdata/qg_4densvar_advect_bump/test_00_nicas-2-sqrt_0001-0001_common.nc

# 4dvar_drplanczos_bump
ln -sf ../qg_parameters_bump_cov/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc testdata/qg_4dvar_drplanczos_bump/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc
ln -sf ../qg_parameters_bump_cov/test_00_var.nc testdata/qg_4dvar_drplanczos_bump/test_00_var.nc

# 4dvar_drplanczos_hybrid_bump (hybrid localization is not very good, take localization alone instead)
ln -sf ../qg_parameters_bump_loc_3d/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc testdata/qg_4dvar_drplanczos_hybrid_bump/test_00_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# Test passed!
echo -e "PASSED"
