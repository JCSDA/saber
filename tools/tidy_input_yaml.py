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
              "colorlog",
              "testing",
              "default_seed",
              "repro",
              "rth",
              "parallel_io",
              "nprocio",
              "universe_rad",
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
              "io_keys",
              "io_values",
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
              "nc1",
              "nc2",
              "nc3",
              "nc4",
              "dc",
              "nl0r",
              "irmax",
              "samp_interp_type",
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
              "vbal_id",
              "vbal_id_coef",
              "forced_var",
              "stddev",
              "var_filter",
              "var_niter",
              "var_npass",
              "var_rhflt",
              "local_diag",
              "local_rad",
              "local_dlat",
              "optimality_nfac",
              "optimality_delta",
              "optimality_ntest",
              "diag_rhflt",
              "diag_rvflt",
              "fit_dl0",
              "fit_ncmp",
              "write_hdiag_detail",
              "resol",
              "nc1max",
              "nicas_draw_type",
              "forced_radii",
              "ncmp",
              "rh",
              "rv",
              "loc_wgt",
              "min_lev",
              "max_lev",
              "nicas_interp_type",
              "pos_def_test",
              "write_nicas_grids",
              "interp_test",
              "ndir",
              "londir",
              "latdir",
              "levdir",
              "ivdir",
              "nldwv",
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
              "ensemble",
              "lowres ensemble",
              "grids",
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
