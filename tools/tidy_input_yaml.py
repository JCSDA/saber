#!/usr/bin/env python3

import argparse
import datetime
import os
import yaml

# Correct date formatting
def correct_date(d):
    if isinstance(d, dict):
        # Loop over dictionary items
        for k,v in d.items():
            d[k] = correct_date(v)
    elif isinstance(d, list):
        # Loop over list items
        i = 0
        for v in d:
            d[i] = correct_date(v)
            i += 1

    if isinstance(d, datetime.datetime):
        # Replace with string
        d = d.strftime("%Y-%m-%dT%H:%M:%SZ")
    return d

# BUMP items order
bump_order = ["datadir",
              "prefix",
              "model",
              "verbosity",
              "colorlog",
              "testing",
              "default_seed",
              "repro",
              "rth",
              "parallel_io",
              "nprocio",
              "universe_rad",
              "use_cgal",
              "write_c0",
              "method",
              "strategy",
              "new_normality",
              "new_vbal_cov",
              "update_vbal_cov",
              "load_vbal_cov",
              "write_vbal_cov",
              "new_vbal",
              "load_vbal",
              "write_vbal",
              "new_var",
              "update_var",
              "load_var",
              "write_var",
              "new_mom",
              "update_mom",
              "load_mom",
              "write_mom",
              "new_hdiag",
              "write_hdiag",
              "new_nicas",
              "load_nicas_local",
              "load_nicas_global",
              "write_nicas_local",
              "write_nicas_global",
              "new_wind",
              "load_wind_local",
              "write_wind_local",
              "check_vbal",
              "check_adjoints",
              "check_normalization",
              "check_dirac",
              "check_randomization",
              "check_consistency",
              "check_optimality",
              "check_no_point_mpi",
              "check_no_point_mask",
              "check_set_param",
              "check_get_param",
              "check_apply_vbal",
              "check_apply_stddev",
              "check_apply_nicas",
              "fname_var",
              "fname_samp",
              "fname_vbal_cov",
              "fname_vbal",
              "fname_mom",
              "fname_nicas",
              "fname_wind",
              "nl0",
              "levs",
              "lev2d",
              "logpres",
              "nv",
              "variables",
              "variable_change",
              "nomask",
              "io_keys",
              "io_values",
              "qg_regional",
              "qg_urban",
              "ens1_ne",
              "ens1_nsub",
              "ens2_ne",
              "ens2_nsub",
              "load_samp_local",
              "load_samp_global",
              "write_samp_local",
              "write_samp_global",
              "write_samp_grids",
              "mask_type",
              "mask_lu",
              "mask_th",
              "ncontig_th",
              "mask_check",
              "diag_draw_type",
              "Lcoast",
              "rcoast",
              "nc1",
              "nc2",
              "nc3",
              "nc3",
              "dc",
              "nl0r",
              "irmax",
              "ne",
              "ne_lr",
              "gen_kurt_th",
              "gau_approx",
              "avg_nbins",
              "vbal_block",
              "vbal_rad",
              "vbal_dlat",
              "vbal_diag_auto",
              "vbal_diag_reg",
              "vbal_pseudo_inv",
              "vbal_pseudo_inv_mmax",
              "vbal_pseudo_inv_var_th",
              "forced_var",
              "stddev",
              "var_filter",
              "var_niter",
              "var_npass",
              "var_rhflt",
              "local_diag",
              "local_rad",
              "local_dlat",
              "diag_rhflt",
              "diag_rvflt",
              "fit_dl0",
              "fit_ncmp",
              "write_hdiag_detail",
              "resol",
              "nc1max",
              "nicas_draw_type",
              "network",
              "forced_radii",
              "ncmp",
              "rh",
              "rv",
              "loc_wgt",
              "min_lev",
              "max_lev",
              "pos_def_test",
              "write_nicas_grids",
              "ndir",
              "londir",
              "latdir",
              "levdir",
              "ivdir",
              "full_grid_smoother_nn",
              "nldwv",
              "img_ldwv",
              "lon_ldwv",
              "lat_ldwv",
              "name_ldwv",
              "wind_streamfunction",
              "wind_velocity_potential",
              "wind_zonal",
              "wind_meridional",
              "wind_nlon",
              "wind_nlat",
              "wind_nsg",
              "wind_inflation",
              "universe radius",
              "input number of components",
              "input",
              "ensemble",
              "lowres ensemble",
              "msvalr",
              "grids",
              "output number of components",
              "output",
              "operators application"]

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="Yaml file name")
args = parser.parse_args()

# Read yaml file
with open(args.filename, "r") as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Correct date formatting recursively
correct_date(config)

# Reorder items in BUMP sections
if "bump" in config:
    bump = config["bump"]
    bump_reordered = {}
    for item in bump_order:
        if item in bump:
            bump_reordered[item] = bump[item]
    config["bump"] = bump_reordered

if "background error" in config:
    if "saber blocks" in config["background error"]:
        blocks = config["background error"]["saber blocks"]
        blocks_reordered = []
        for block in blocks:
            if "bump" in block:
                bump = block["bump"]
                bump_reordered = {}
                for item in bump_order:
                    if item in bump:
                        bump_reordered[item] = bump[item]
                block["bump"] = bump_reordered
            blocks_reordered.append(block)
        config["background error"]["saber blocks"] = blocks_reordered

# Write reordered yaml file
with open(args.filename, "w") as file:
    output = yaml.dump(config, file, sort_keys=False)
