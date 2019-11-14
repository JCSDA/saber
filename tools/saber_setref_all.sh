#!/bin/sh
#----------------------------------------------------------------------
# Shell script: saber_setref_all
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
testref=$1
testdata=$2

saber_test="
bump_cortrack
bump_hdiag_cor_diag_all
bump_hdiag_grid
bump_hdiag_hyb-rnd_common
bump_hdiag_mask_stddev_upper
bump_hdiag-nicas_cor_specific_univariate
bump_hdiag-nicas_local_diag_var_filter
bump_hdiag-nicas_loc_common
bump_hdiag-nicas_mask_check
bump_hdiag-nicas_network
bump_lct-nicas_one_scale
bump_nicas_fast_sampling
bump_nicas_subsamp_hvh
bump_read_cmat
bump_read_mom
bump_read_nicas
bump_read_obsop
bump_read_sampling
bump_read_vbal
bump_write_cmat
bump_write_mom
bump_write_nicas
bump_write_obsop
bump_write_sampling
bump_write_vbal
bump_get_param_cor
bump_get_param_Dloc
bump_get_param_hyb
bump_get_param_lct
bump_hdiag_draw_type
bump_hdiag_hyb-avg_common
bump_hdiag_ldwv
bump_hdiag_loc_gau_approx
bump_hdiag_loc_histograms
bump_hdiag_mask_lat
bump_hdiag_mask_stddev_lower
bump_hdiag_mask_stddev_ncontig
bump_hdiag-nicas_cor_common_univariate
bump_hdiag-nicas_double_fit
bump_hdiag-nicas_fit_type
bump_hdiag-nicas_lhom
bump_hdiag-nicas_loc_common_univariate
bump_hdiag-nicas_loc_common_weighted
bump_hdiag-nicas_loc_specific_multivariate
bump_hdiag-nicas_loc_specific_univariate
bump_hdiag-nicas_nonunit_diag
bump_hdiag-nicas_nprocio
bump_hdiag-nicas_nrep
bump_hdiag-nicas_rvflt
bump_interface_nicas
bump_interface_obsop
bump_interface_vbal
bump_lct-nicas_diagonal
bump_lct-nicas_mask_check
bump_lct-nicas_network
bump_lct_qc
bump_lct_two_scales
bump_lct_write_cor
bump_nicas_mpicom_lsqrt_a
bump_nicas_mpicom_lsqrt_b
bump_nicas_mpicom_lsqrt_c
bump_nicas_pos_def_test
bump_nicas_subsamp_h
bump_nicas_subsamp_hv
bump_nicas_subsamp_vh
bump_nicas_write_grids
bump_obsop
bump_set_param_cor
bump_set_param_hyb
bump_set_param_lct
bump_hdiag-nicas_loc_adv
bump_hdiag-nicas_loc_adv_cor_tracker
bump_hdiag-nicas_loc_adv_inv
bump_hdiag_optimality
bump_nicas_consistency
bump_nicas_randomization"

for test in ${saber_test}; do
   ./saber_setref.sh $1 $2 ${test}
done
