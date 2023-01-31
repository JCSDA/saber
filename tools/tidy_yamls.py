#!/usr/bin/env python3

import argparse
import datetime
import os
import yaml
import copy

# BUMP items order
template = {
  # General
  "general": {
    "color log": "value",
    "testing": "value",
    "default seed": "value",
    "reproducibility operators": "value",
    "reproducibility threshold": "value",
    "universe length-scale": "value"
  },

  # I/O
  "io": {
    "data directory": "value",
    "files prefix": "value",
    "parallel netcdf": "value",
    "io tasks": "value",
    "alias": {
      "subtemplate": "vector",
      "items": {
        "in code": "value",
        "in file": "value"
      }
    },
    "overriding sampling file": "value",
    "overriding vertical covariance file": "value",
    "overriding vertical balance file": "value",
    "overriding moments file": "value",
    "overriding lowres moments file": "value",
    "overriding nicas file": "value",
    "overriding psichitouv file": "value"
  },

  # Drivers
  "drivers": {
    "compute covariance": "value",
    "compute lowres covariance": "value",
    "compute correlation": "value",
    "compute lowres correlation": "value",
    "compute localization": "value",
    "compute lowres localization": "value",
    "compute hybrid weights": "value",
    "hybrid source": "value",
    "multivariate strategy": "value",
    "iterative algorithm": "value",
    "compute normality": "value",
    "read local sampling": "value",
    "read global sampling": "value",
    "write local sampling": "value",
    "write global sampling": "value",
    "write sampling grids": "value",
    "compute vertical covariance": "value",
    "read vertical covariance": "value",
    "write vertical covariance": "value",
    "compute vertical balance": "value",
    "read vertical balance": "value",
    "write vertical balance": "value",
    "compute variance": "value",
    "compute moments": "value",
    "read moments": "value",
    "write moments": "value",
    "write diagnostics": "value",
    "write diagnostics detail": "value",
    "compute nicas": "value",
    "read local nicas": "value",
    "read global nicas": "value",
    "write local nicas": "value",
    "write global nicas": "value",
    "write nicas grids": "value",
    "compute psichitouv": "value",
    "read local psichitouv": "value",
    "write local psichitouv": "value",
    "vertical balance inverse test": "value",
    "adjoints test": "value",
    "normalization test": "value",
    "internal dirac test": "value",
    "randomization test": "value",
    "internal consistency test": "value",
    "localization optimality test": "value"
  },

  # Model
  "model": {
    "level for 2d variables": "value",
    "variables": "value",
    "groups": {
      "subtemplate": "vector",
      "items": {
        "group name": "value",
        "variables": "value",
      }
    },
    "do not cross mask boundaries": "value"
  },

  # Ensemble size
  "ensemble sizes": {
    "total ensemble size": "value",
    "sub-ensembles": "value",
    "total lowres ensemble size": "value",
    "lowres sub-ensembles": "value"
  },

  # Sampling
  "sampling": {
    "computation grid size": "value",
    "diagnostic grid size": "value",
    "distance classes": "value",
    "angular sectors": "value",
    "distance class width": "value",
    "reduced levels": "value",
    "local diagnostic": "value",
    "averaging length-scale": "value",
    "averaging latitude width": "value",
    "grid type": "value",
    "max number of draws": "value",
    "interpolation type": "value",
    "masks": {
      "subtemplate": "vector",
      "items": {
        "type": "value",
        "threshold": "value",
        "side": "value",
        "variable": "value"
      }
    },
    "contiguous levels threshold": "value"
  },

  # Diagnostics
  "diagnostics": {
    "target ensemble size": "value",
    "target lowres ensemble size": "value",
    "gaussian approximation": "value",
    "generalized kurtosis threshold": "value",
    "histogram bins": "value"
  },
   
  # Vertical balance 
  "vertical balance": {
    "vbal": {
      "subtemplate": "vector",
      "items": {
        "balanced variable": "value",
        "unbalanced variable": "value",
        "diagonal autocovariance": "value",
        "diagonal regression": "value",
        "identity block weight": "value"
      }
    },
    "pseudo inverse": "value",
    "dominant mode": "value",
    "variance threshold": "value",
    "identity blocks": "value"
  },

  # Variance
  "variance": {
    "explicit stddev": "value",
    "stddev": {
      "subtemplate": "vector",
      "items": {
        "variables": "value",
        "value": "value",
        "profile": "value",
      }
    },
    "objective filtering": "value",
    "filtering iterations": "value",
    "filtering passes": "value",
    "initial length-scale": {
      "subtemplate": "vector",
      "items": {
        "variables": "value",
        "value": "value",
        "profile": "value",
      }
    },
  },

  # Optimality test
  "optimality test": {
    "half number of factors": "value",
    "factors increment": "value",
    "test vectors": "value"
  },

  # Fit
  "fit": {
    "horizontal filtering length-scale": "value",
    "vertical filtering length-scale": "value",
    "vertical stride": "value",
    "number of components": "value"
  },

  # Local profiles
  "local profiles": {
    "subtemplate": "vector",
    "items": {
      "longitude": "value",
      "latitude": "value",
      "name": "value",
    }
  },

  # NICAS
  "nicas": {
    "resolution": "value",
    "max horizontal grid size": "value",
    "grid type": "value",
    "explicit length-scales": "value",
    "horizontal length-scale": {
      "subtemplate": "vector",
      "items": {
        "groups": "value",
        "value": "value",
        "profile": "value",
      }
    },
    "vertical length-scale": {
      "subtemplate": "vector",
      "items": {
        "groups": "value",
        "value": "value",
        "profile": "value",
      }
    },
    "common localization weights": {
      "subtemplate": "vector",
      "items": {
        "row variables": "value",
        "column variables": "value",
        "value": "value",
      }
    },
    "minimum level": {
      "subtemplate": "vector",
      "items": {
        "groups": "value",
        "value": "value",
      }
    },
    "maximum level": {
      "subtemplate": "vector",
      "items": {
        "groups": "value",
        "value": "value",
      }
    },
    "interpolation type": {
      "subtemplate": "vector",
      "items": {
        "groups": "value",
        "type": "value",
      }
    },
    "positive-definiteness test": "value",
    "horizontal interpolation test": "value"
  },

  # Psichitouv
  "psichitouv": {
    "stream function": "value",
    "velocity potential": "value",
    "eastward wind": "value",
    "northward wind": "value",
    "longitudes": "value",
    "latitudes": "value",
    "savitzky-golay half width": "value",
    "wind inflation": "value"
  },

  # Dirac
  "dirac": {
    "subtemplate": "vector",
    "items": {
      "longitude": "value",
      "latitude": "value",
      "level": "value",
      "variable": "value",
    }
  },

  # External parameter
  "msvalr": "value",
  "ensemble": "value",
  "lowres ensemble": "value",
  "grids": "value",
  "operators application": "value"
}


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

# Update a bump section
def reorder_bump_from_template(config, template):
    output = {} 
    if isinstance(template, dict):
        for k,v in template.items():
            if k in config:
                if isinstance(v, dict):
                    if "subtemplate" in v:
                        if v["subtemplate"] == "vector":
                            output[k] = []
                            for item in config[k]:
                                output[k].append(reorder_bump_from_template(item, v["items"]))
                    else:
                        output[k] = reorder_bump_from_template(config[k], v)
                elif v == "value":
                    output[k] = config[k]
        return output

# Find bump subsections and update them recursively
def update_bump(d, template):
    if isinstance(d, dict):
        # Loop over dictionary items
        reordered = {}
        for k,v in d.items():
            # Check for bump subsection
            if k == "bump":
                # Reorder data
                reordered[k] = reorder_bump_from_template(v, template)
            else:
                # Recursive call
                if isinstance(v, dict):
                    reordered[k] = update_bump(v, template)
                elif isinstance(v, list):
                    vlist = []
                    for vv in v:
                        if isinstance(vv, dict):
                            vlist.append(update_bump(vv, template))
                        else:
                            vlist.append(copy.deepcopy(vv))
                    reordered[k] = copy.deepcopy(vlist)
                else:
                    reordered[k] = copy.deepcopy(v)
        return reordered
    else:
        print("ERROR!")
        print(d)
        exit()

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

# Look for bump sections and update them
reordered_config = update_bump(config, template)

# Write reordered yaml file
with open(args.filename, "w") as file:
    output = yaml.dump(reordered_config, file, sort_keys=False)
