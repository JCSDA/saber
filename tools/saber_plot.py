#!/usr/bin/env python3
# Examples:
# * python saber_plot.py ${build}/saber/test/testdata error_covariance_training_hdiag_1 1 1
"""! Plot management script"""

import argparse
import os
import sys

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("testdata", help="Test data directory")
parser.add_argument("test", help="Test name")
parser.add_argument("mpi", help="Number of MPI tasks")
parser.add_argument("omp", help="Number of OpenMP threads")
parser.add_argument("--variable", help="Variable")
parser.add_argument("--levels_min", help="Minimum level")
parser.add_argument("--levels_max", help="Maximum level")
parser.add_argument("--levels_step", help="Step between levels")

# Parse and print arguments
args = parser.parse_args()
print("Parameters:")
for arg in vars(args):
    if not arg is None:
        print(" - " + arg + ": " + str(getattr(args, arg)))

# Insert path
sys.path.insert(1, os.path.join(os.path.dirname(os.path.realpath(__file__)), "saber_plot"))

# Available plots list
plot_list=["avg1","diag","normality","randomization","sampling_grids"]
done_list=["normality","sampling_grids"]
done = {}
alias_list = {
   "avg1": "avg1",
   "diag": "diag",
   "normality": "normality",
   "randomization": "randomization",
   "sampling_grids": "sampling_grids",
}
for plot in plot_list:
   done[plot] = False

# Unavailable plots list
noplot_list=["mom","nicas","sampling","vbal","wind"]

# Create jpg figures
for f in sorted(os.listdir(os.path.join(args.testdata, args.test))):
   if f.endswith(".nc"):
      if os.path.isfile(os.path.join(args.testdata, args.test, f)) and f.find(args.mpi + "-" + args.omp)==0:
         suffix = f.split(args.mpi + "-" + args.omp + "_")[1][:-(len(f.rsplit(".")[-1])+1)]
         quenchFile = True
         for plot in plot_list:
            if suffix.find(plot)==0:
               quenchFile = False
               if not done[plot]:
                  print("Calling " + plot + " in " + args.test + " (" + args.mpi + "-" + args.omp + ")")
                  plot_alias = "bump_" + alias_list[plot]
                  module = __import__(plot_alias)
                  func = getattr(module, plot_alias)
                  if plot in done_list:
                     func(args, plot)
                     done[plot] = True
                  else:
                     func(args, suffix)
         for plot in noplot_list:
            if suffix.find(plot)==0:
               quenchFile = False
         if quenchFile:
            print("Calling quench plot in " + args.test + " (" + args.mpi + "-" + args.omp + ")")
            module = __import__("quench")
            func = getattr(module, "func")
            func(args, suffix)

# Create HTML page
message = "<html><head></head><body><h1>" + args.test + "</h1><ul>"
if os.path.exists(args.test):
   for f in sorted(os.listdir(args.test)):
      if f.find(args.mpi + "-" + args.omp)==0:
         short_name = f.replace(args.mpi + "-" + args.omp + "_", "").replace(".jpg", "")
         message = message + "<li><a href=\"#" + short_name + "\">" + short_name + "</a></li>"
   message = message + "</ul>"
   for f in sorted(os.listdir(args.test)):
      if f.find(args.mpi + "-" + args.omp)==0:
         short_name = f.replace(args.mpi + "-" + args.omp + "_", "").replace(".jpg", "")
         cmd = "mogrify -trim " + os.path.join(args.test, f)
         os.system(cmd)
         message = message + "<h2 id=\"" + short_name + "\">" + short_name + "</h2><img src=\"" + f + "\" width=800px><br><a href=\"#top\">Back to top</a>"

message = message + "</body></html>"
f = open(os.path.join(args.test, "index_" + args.mpi + "-" + args.omp + ".html"), "w")
f.write(message)
f.close()
