#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_qg_links
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# dirac_bump_cov
ln -sf ../qg_parameters_bump_cov/test_nicas_000001-000001.nc testdata/qg_dirac_bump_cov/test_nicas_000001-000001.nc
ln -sf ../qg_parameters_bump_cov/test_var.nc testdata/qg_dirac_bump_cov/test_var.nc

# dirac_bump_cov_no_ens
ln -sf ../qg_parameters_bump_cov_no_ens/test_nicas_000001-000001.nc testdata/qg_dirac_bump_cov_no_ens/test_nicas_000001-000001.nc
ln -sf ../qg_parameters_bump_cov_no_ens/test_var.nc testdata/qg_dirac_bump_cov_no_ens/test_var.nc

# dirac_bump_lct
ln -sf ../qg_parameters_bump_lct/test_nicas_000001-000001.nc testdata/qg_dirac_bump_lct/test_nicas_000001-000001.nc
ln -sf ../qg_parameters_bump_lct/test_var.nc testdata/qg_dirac_bump_lct/test_var.nc

# dirac_bump_hyb (hybrid localization is not very good, take localization alone instead)
ln -sf ../qg_parameters_bump_loc/test_nicas_000001-000001.nc testdata/qg_dirac_bump_hyb/test_nicas_000001-000001.nc

# dirac_bump_loc_3d
ln -sf ../qg_parameters_bump_loc/test_nicas_000001-000001.nc testdata/qg_dirac_bump_loc_3d/test_nicas_000001-000001.nc

# dirac_bump_loc_4d
ln -sf ../qg_parameters_bump_loc/test_nicas_000001-000001.nc testdata/qg_dirac_bump_loc_4d/test_nicas_000001-000001.nc

# dirac_id_cov
ln -sf ../qg_parameters_bump_cov/test_var.nc testdata/qg_dirac_id_cov/test_var.nc

# 3densvar_bump
ln -sf ../qg_parameters_bump_loc/test_nicas_000001-000001.nc testdata/qg_3densvar_bump/test_nicas_000001-000001.nc

# 3dvar_bump
ln -sf ../qg_parameters_bump_cov/test_nicas_000001-000001.nc testdata/qg_3dvar_bump/test_nicas_000001-000001.nc
ln -sf ../qg_parameters_bump_cov/test_var.nc testdata/qg_3dvar_bump/test_var.nc

# 3dvar_hybrid_bump (hybrid localization is not very good, take localization alone instead)
ln -sf ../qg_parameters_bump_loc/test_nicas_000001-000001.nc testdata/qg_3dvar_hybrid_bump/test_nicas_000001-000001.nc

# 4densvar_bump
ln -sf ../qg_parameters_bump_loc/test_nicas_000001-000001.nc testdata/qg_4densvar_bump/test_nicas_000001-000001.nc

# 4dvar_drplanczos_bump
ln -sf ../qg_parameters_bump_cov/test_nicas_000001-000001.nc testdata/qg_4dvar_drplanczos_bump/test_nicas_000001-000001.nc
ln -sf ../qg_parameters_bump_cov/test_var.nc testdata/qg_4dvar_drplanczos_bump/test_var.nc

# 4dvar_drplanczos_hybrid_bump (hybrid localization is not very good, take localization alone instead)
ln -sf ../qg_parameters_bump_loc/test_nicas_000001-000001.nc testdata/qg_4dvar_drplanczos_hybrid_bump/test_nicas_000001-000001.nc

# Test passed!
echo -e "PASSED"
