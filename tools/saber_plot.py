#!/usr/bin/env python3

import argparse
from os import listdir
from os.path import isfile, islink, join
import sys

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("testdata", help="Test data directory")
parser.add_argument("test", help="Test name")
parser.add_argument("mpi", help="Number of MPI tasks")
parser.add_argument("omp", help="Number of OpenMP threads")
args = parser.parse_args()

# Insert path
sys.path.insert(1, 'plot')

# Available plots list
plot_list=["avg","corstats","cortrack","diag","dirac","lct","local_diag_cor_gridded","local_diag_cor","local_diag_loc_gridded","local_diag_loc","normality","umf"]
done_list=["normality"]
done = {}
for plot in plot_list:
   done[plot] = False

# BUMP tests
if args.test.find("bump_")==0:
   for f in listdir(join(args.testdata, args.test)):
      if isfile(join(args.testdata, args.test, f)) and f.find("test_" + args.mpi + "-" + args.omp)==0:
         suffix = f.split("test_" + args.mpi + "-" + args.omp + "_")[1][:-(len(f.rsplit(".")[-1])+1)]
         for plot in plot_list:
            if suffix.find(plot)==0 and not done[plot]:
               print("Calling " + plot + " in " + args.test + " (" + args.mpi + "-" + args.omp + ")")
               module = __import__(plot)
               func = getattr(module, plot)
               if plot in done_list:
                  func(args.testdata, args.test, args.mpi, args.omp, plot)
                  done[plot] = True
               else:
                  func(args.testdata, args.test, args.mpi, args.omp, suffix)
