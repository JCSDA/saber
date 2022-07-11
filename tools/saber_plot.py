#!/usr/bin/env python3
# Examples:
# * bump test: python saber_plot.py ${build}/bin bump ${build}/saber/test/testdata bump_hdiag-nicas_cor_specific_univariate 1 1
# * QG test: python saber_plot.py ${build}/bin qg fields ${build}/saber/test/testdata/qg_parameters_bump_cov/stddev.an.2010-01-01T12\:00\:00Z.nc
# * QUENCH test: python saber_plot.py ${build}/bin quench ${build}/test/testdata/quench_bump_nicas/dirac.nc
"""! Plot management script"""

import argparse
import os
import sys

# Main parser and sub-parser
parser = argparse.ArgumentParser()
parser.add_argument("bindir", help="Binary directory")
subparser = parser.add_subparsers(help="help", required=True, dest="case")

# BUMP sub-parser
parser_bump = subparser.add_parser("bump", help="BUMP parser")
parser_bump.add_argument("testdata", help="Test data directory")
parser_bump.add_argument("test", help="Test name")
parser_bump.add_argument("mpi", help="Number of MPI tasks")
parser_bump.add_argument("omp", help="Number of OpenMP threads")
parser_bump.add_argument("--output", help="Output file path")

# QG sub-parser (from oops/tools/plot.py)
parser_qg = subparser.add_parser("qg", help="QG parser")

# QG diagnostics sub-parsers
subparser_diagnostic_qg = parser_qg.add_subparsers(help="QG diagnostic help", required=True, dest="diagnostic")
parser_qg_cost = subparser_diagnostic_qg.add_parser("cost", help="QG cost parser")
parser_qg_fields = subparser_diagnostic_qg.add_parser("fields", help="QG fields parser")
parser_qg_obs = subparser_diagnostic_qg.add_parser("obs", help="QG obs parser")

# QG cost arguments
parser_qg_cost.add_argument("filepath", type=str, help="File path")
parser_qg_cost.add_argument("--output", help="Output file path")

# QG fields arguments
parser_qg_fields.add_argument("filepath", help="File path")
parser_qg_fields.add_argument("basefilepath", nargs="?", help="Base file path", default=None)
parser_qg_fields.add_argument("--plotwind", dest="plotwind", action="store_true", help="Plot wind")
parser_qg_fields.add_argument("--gif", help="Gif pattern values, separated with commas, replacing %id% in the file path")
parser_qg_fields.add_argument("--output", help="Output file path")

# QG obs arguments
parser_qg_obs.add_argument("filepath", type=str, help="File path")
parser_qg_obs.add_argument("--output", help="Output file path")

# QUENCH sub-parser
parser_quench = subparser.add_parser("quench", help="QUENCH parser")
parser_quench.add_argument("filepath", help="File path")
parser_quench.add_argument("variable", help="Variable")
parser_quench.add_argument("netcdf", nargs="?", help="NetCDF format input", default=True)
parser_quench.add_argument("--output", help="Output file path")
parser_quench.add_argument("--levels_min", help="Minimum level")
parser_quench.add_argument("--levels_max", help="Maximum level")
parser_quench.add_argument("--levels_step", help="Step between levels")

# Parse and print arguments
args = parser.parse_args()
print("Parameters:")
for arg in vars(args):
    if not arg is None:
        print(" - " + arg + ": " + str(getattr(args, arg)))

if (args.case == "bump"):
   # Insert path
   sys.path.insert(1, os.path.join(args.bindir, "saber_plot"))

   # Available plots list
   plot_list=["avg","diag","diag_dirac","dirac","normality","randomization","sampling_grids","umf","var"]
   done_list=["normality","sampling_grids","diag","diag_dirac"]
   done = {}
   alias_list = {
      "avg": "avg",
      "diag": "diag",
      "diag_dirac": "contour_centered",
      "dirac": "contour_centered",
      "normality": "normality",
      "randomization": "randomization",
      "sampling_grids": "sampling_grids",
      "umf": "umf",
      "var": "contour_positive"
   }
   for plot in plot_list:
      done[plot] = False

   # BUMP tests
   if args.test.find("bump_")==0:
      # Create jpg figures
      for f in sorted(os.listdir(os.path.join(args.testdata, args.test))):
         if os.path.isfile(os.path.join(args.testdata, args.test, f)) and f.find("test_" + args.mpi + "-" + args.omp)==0:
            suffix = f.split("test_" + args.mpi + "-" + args.omp + "_")[1][:-(len(f.rsplit(".")[-1])+1)]
            for plot in plot_list:
               if suffix.find(plot)==0 and not done[plot]:
                  print("Calling " + plot + " in " + args.test + " (" + args.mpi + "-" + args.omp + ")")
                  plot_alias = "bump_" + alias_list[plot]
                  module = __import__(plot_alias)
                  func = getattr(module, plot_alias)
                  if plot in done_list:
                     func(args, plot)
                     done[plot] = True
                  else:
                     func(args, suffix)

      if (os.path.isdir(os.path.join(args.testdata, args.test, "fig"))):
         # Create HTML page
         message = "<html><head></head><body><h1>" + args.test + "</h1><ul>"
         for f in sorted(os.listdir(os.path.join(args.testdata, args.test, "fig"))):
            if f.find("test_" + args.mpi + "-" + args.omp)==0:
               short_name = f.replace("test_" + args.mpi + "-" + args.omp + "_", "").replace(".jpg", "")
               message = message + "<li><a href=\"#" + short_name + "\">" + short_name + "</a></li>"
         message = message + "</ul>"
         for f in sorted(os.listdir(os.path.join(args.testdata, args.test, "fig"))):
            if f.find("test_" + args.mpi + "-" + args.omp)==0:
               short_name = f.replace("test_" + args.mpi + "-" + args.omp + "_", "").replace(".jpg", "")
               cmd = "mogrify -trim " + os.path.join(args.testdata, args.test, "fig", f)
               os.system(cmd)
               message = message + "<h2 id=\"" + short_name + "\">" + short_name + "</h2><img src=\"" + f + "\" width=800px><br><a href=\"#top\">Back to top</a>"

         message = message + "</body></html>"
         f = open(os.path.join(args.testdata, args.test, "fig", "index_" + args.mpi + "-" + args.omp + ".html"), "w")
         f.write(message)
         f.close()

if (args.case == "qg"):
   # Insert path
   sys.path.insert(1, os.path.join(args.bindir, "oops_plot"))

   # Run plot function
   module = __import__("qg_" + args.diagnostic)
   func = getattr(module, "func")
   print("Run script")
   func(args)

if (args.case == "quench"):
   # Insert path
   sys.path.insert(1, os.path.join(args.bindir, "saber_plot"))

   # Run plot function
   module = __import__("quench")
   func = getattr(module, "func")
   print("Run script")
   func(args)
