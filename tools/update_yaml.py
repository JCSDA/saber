#!/usr/bin/env python3

import argparse
import datetime
import os
import yaml
import numpy as np

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

# Find bump subsections and their indentation
def find_bump(d, bumps):
    if isinstance(d, dict):
        # Loop over dictionary items
        for k,v in d.items():
            # Check for bump subsection
            if k == "bump":
                bumps.append(v)
            else:
                find_bump(v, bumps)
    elif isinstance(d, list):
        # Loop over list items
        for v in d:
            find_bump(v, bumps)

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="Yaml file name")
parser.add_argument('--variables', nargs='+', help='Variables', default=['var1','var2','var3','var4'])
args = parser.parse_args()
print("File: " + args.filename)

# Read yaml file
with open(args.filename, "r") as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Correct date formatting recursively
correct_date(config)

# Look for bump sections
bumps = []
find_bump(config, bumps)

# BUMP structure
kv = []

general = {}
general["name"] = "general"
general["keys"] = ["colorlog", "testing", "default_seed", "repro", "rth", "universe_rad"]
kv.append(general)

io = {}
io["name"] = "io"
io["keys"] = ["datadir", "prefix", "parallel_io", "nprocio","fname_var", "fname_samp", "fname_vbal_cov", "fname_vbal", "fname_mom", "fname_nicas", "fname_wind"]
kv.append(io)

drivers = {}
drivers["name"] = "drivers"
drivers["keys"] = ["method", "strategy", "new_normality", "load_samp_local", "load_samp_global", "write_samp_local", "write_samp_global", "write_samp_grids", "new_vbal_cov", "update_vbal_cov", "load_vbal_cov", "write_vbal_cov", "new_vbal", "load_vbal", "write_vbal", "new_var", "update_var", "new_mom", "update_mom", "load_mom", "write_mom", "new_hdiag", "write_hdiag", "write_hdiag_detail", "new_nicas", "load_nicas_local", "load_nicas_global", "write_nicas_local", "write_nicas_global", "write_nicas_grids", "new_wind", "load_wind_local", "write_wind_local", "check_vbal", "check_adjoints", "check_normalization", "check_dirac", "check_randomization", "check_consistency", "check_optimality", "check_set_param", "check_get_param", "check_apply_vbal", "check_apply_stddev", "check_apply_nicas"]
kv.append(drivers)

model = {}
model["name"] = "model"
model["keys"] = ["nl0", "levs", "lev2d", "variables"]
kv.append(model)

ensembleSizes = {}
ensembleSizes["name"] = "ensemble sizes"
ensembleSizes["keys"] = ["ens1_ne", "ens1_nsub", "ens2_ne", "ens2_nsub"]
kv.append(ensembleSizes)

mask = {}
mask["name"] = "mask"
mask["keys"] = ["mask_type", "mask_lu", "mask_th", "ncontig_th", "mask_check"]
kv.append(mask)

sampling = {}
sampling["name"] = "sampling"
sampling["keys"] = ["nc1", "nc2", "nc3", "nc4", "dc", "nl0r", "local_diag", "local_rad", "local_dlat", "irmax"]
kv.append(sampling)

localization = {}
localization["name"] = "localization"
localization["keys"] = ["ne", "ne_lr", "gau_approx", "gen_kurt_th", "avg_nbins"]
kv.append(localization)

verticalBalance = {}
verticalBalance["name"] = "vertical balance"
verticalBalance["keys"] = ["vbal_rad", "vbal_dlat", "vbal_pseudo_inv", "vbal_pseudo_inv_mmax", "vbal_pseudo_inv_var_th", "vbal_id"]
kv.append(verticalBalance)

variance = {}
variance["name"] = "variance"
variance["keys"] = ["forced_var", "var_filter", "var_niter", "var_npass"]
kv.append(variance)

optimalityTest = {}
optimalityTest["name"] = "optimality test"
optimalityTest["keys"] = ["optimality_nfac", "optimality_delta", "optimality_ntest"]
kv.append(optimalityTest)

fit = {}
fit["name"] = "fit"
fit["keys"] = ["diag_rhflt", "diag_rvflt", "fit_dl0", "fit_ncmp"]
kv.append(fit)

localProfiles = {}
localProfiles["name"] = "local profiles"
localProfiles["keys"] = ["nldwv", "lon_ldwv", "lat_ldwv", "name_ldwv"]
kv.append(localProfiles)

nicas = {}
nicas["name"] = "nicas"
nicas["keys"] = ["resol", "nc1max", "nicas_draw_type", "forced_radii", "pos_def_test", "interp_test"]
kv.append(nicas)

dirac = {}
dirac["name"] = "dirac"
dirac["keys"] = ["ndir", "londir", "latdir", "levdir", "ivdir"]
kv.append(dirac)

wind = {}
wind["name"] = "wind"
wind["keys"] = ["wind_streamfunction", "wind_velocity_potential", "wind_zonal", "wind_meridional", "wind_nlon", "wind_nlat", "wind_nsg", "wind_inflation"]
kv.append(wind)

other_sections = ["ensemble", "lowres ensemble", "operators application"]

# Upgrade bump sections
for i in range(len(bumps)):
    # Copy existing keys
    old_bump = bumps[i]
    new_bump = {}
    for j in range(len(kv)):
        section = {}
        for item in kv[j]["keys"]:
            if item in old_bump:
                section[item] = old_bump[item]
        if section:
            new_bump[kv[j]["name"]] = section

    # Copy existing keys in grids
    if "grids" in old_bump:
        old_grids = old_bump["grids"]
        new_grids = []
        for old_grid in old_grids:
            new_grid = {}
            for j in range(len(kv)):
                section = {}
                for item in kv[j]["keys"]:
                    if item in old_grid:
                        section[item] = old_grid[item]
                if section:
                    new_grid[kv[j]["name"]] = section

            # Append grid
            new_grids.append(new_grid)

        # Reset grid
        new_bump["grids"] = new_grids

    # Update io_keys/io_values
    if "io_keys" in old_bump:
        if not "io" in new_bump:
            new_bump["io"] = {}
        alias = []
        for i in range(len(old_bump["io_keys"])):
            item = {}
            item["in code"] = old_bump["io_keys"][i]
            item["in file"] = old_bump["io_values"][i]
            alias.append(item)
        new_bump["io"]["alias"] = alias

    # Udpate diag_draw_type
    if "diag_draw_type" in old_bump:
        if not "sampling" in new_bump:
            new_bump["sampling"] = {}
        new_bump["sampling"]["draw_type"] = old_bump["diag_draw_type"]

    # Udpate samp_interp_type
    if "samp_interp_type" in old_bump:
        if not "sampling" in new_bump:
            new_bump["sampling"] = {}
        new_bump["sampling"]["interp_type"] = old_bump["samp_interp_type"]

    # Udpate vbal
    if "vbal_block" in old_bump:
        vbal_block = old_bump["vbal_block"]
        if "vbal_diag_auto" in old_bump:
            vbal_diag_auto = old_bump["vbal_diag_auto"]
        else:
            vbal_diag_auto = np.full((len(vbal_block)), False)
        if "vbal_diag_reg" in old_bump:
            vbal_diag_reg = old_bump["vbal_diag_reg"]
        else:
            vbal_diag_reg = np.full((len(vbal_block)), False)
        if "vbal_id_coef" in old_bump:
            vbal_id_coef = old_bump["vbal_id_coef"]
        else:
            vbal_id_coef = np.ones((len(vbal_block)))
        ib = 0
        vbal = []
        for ii in range(1, 10):
            for jj in range(0, ii):
                if ib < len(vbal_block):
                    if vbal_block[ib]:
                        block = {}
                        block["balanced"] = args.variables[ii]
                        block["unbalanced"] = args.variables[jj]
                        if vbal_diag_auto[ib]:
                            block["diag_auto"] = True
                        if vbal_diag_reg[ib]:
                            block["diag_reg"] = True
                        if "vbal_id_coef" in old_bump:
                            block["id_coef"] = vbal_id_coef[ib]
                        vbal.append(block)
                ib += 1
        if not "vertical balance" in new_bump:
            new_bump["vertical balance"] = {}
        new_bump["vertical balance"]["vbal"] = vbal

    # Update stddev and var_rhflt
    for key in ["stddev", "var_rhflt"]:
        if key in old_bump:
            vec = []
            for item in old_bump[key]:
                block = {}
                block["variables"] = [item]
                if len(old_bump[key][item]) == 1:
                    block["value"] = old_bump[key][item][0]
                else:
                    block["profile"] = old_bump[key][item]
                vec.append(block)
            if not "variance" in new_bump:
                new_bump["variance"] = {}
            new_bump["variance"][key] = vec

    # Update rh, rv, min_lev and max_lev
    for key in ["rh", "rv", "min_lev", "max_lev"]:
        if key in old_bump:
            vec = []
            for item in old_bump[key]:
                block = {}
                block["variables"] = [item]
                if len(old_bump[key][item]) == 1:
                    block["value"] = old_bump[key][item][0]
                else:
                    block["profile"] = old_bump[key][item]
                vec.append(block)
            if not "nicas" in new_bump:
                new_bump["nicas"] = {}
            new_bump["nicas"][key] = vec

    # Update nicas_interp_type
    if "nicas_interp_type" in old_bump:
        vec = []
        for item in old_bump["nicas_interp_type"]:
            block = {}
            block["variables"] = [item]
            block["type"] = old_bump["nicas_interp_type"][item]
            vec.append(block)
        if not "nicas" in new_bump:
            new_bump["nicas"] = {}
        new_bump["nicas"]["interp_type"] = vec

    # Update loc_wgt
    key = "loc_wgt"
    if key in old_bump:
        vec = []
        for item in old_bump[key]:
            block = {}
            block["row variables"] = [item.split('-')[0]]
            block["column variables"] = [item.split('-')[1]]
            block["value"] = old_bump[key][item]
            vec.append(block)
        if not "nicas" in new_bump:
            new_bump["nicas"] = {}
        new_bump["nicas"][key] = vec

    # Copy other sections
    for other_section in other_sections:
        if other_section in old_bump:
            new_bump[other_section] = old_bump[other_section]

    # Reset bump
    bumps[i] = new_bump

# Transform bump section into text vectors
bumps_text = []
for i in range(len(bumps)):
    file_tmp = open('tmpfile', 'w')
    yaml.dump(bumps[i], file_tmp, sort_keys=False)
    file_tmp.close()
    file_in = open('tmpfile', 'r')
    text = []
    for line in file_in:
       text.append(line)
    file_in.close()
    os.remove('tmpfile')
    bumps_text.append(text)

# Rename file
os.rename(args.filename, args.filename + ".bak")

# Read and rewrite file, updating the bump sections only
file_in = open(args.filename + ".bak", 'r')
file_out = open(args.filename, 'w')
i = 0
blank = ' '
ind_target = -1
for line in file_in:
    ind_current = len(line)-len(line.lstrip(' '))
    if "bump:" in line and not line.startswith("#"):
        ind_target = ind_current
        text = bumps_text[i]
        file_out.writelines(line)
        for newline in text:
            file_out.writelines((ind_target+2)*blank + newline)
        i += 1
    else:
        if ind_current <= ind_target:
            ind_target = -1
        if ind_target == -1:
            file_out.writelines(line)
file_in.close()
file_out.close()

# Remove backup file
os.remove(args.filename + ".bak")
