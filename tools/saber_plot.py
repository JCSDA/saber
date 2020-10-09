#!/usr/bin/env python3

import argparse
import os
import sys

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("bindir", help="Binary directory")
parser.add_argument("testdata", help="Test data directory")
parser.add_argument("test", help="Test name")
parser.add_argument("mpi", help="Number of MPI tasks")
parser.add_argument("omp", help="Number of OpenMP threads")
args = parser.parse_args()

# Insert path
sys.path.insert(1, os.path.join(args.bindir, "saber_plot"))

# Available plots list
plot_list=["avg","diag","dirac","lct","lct_cor","local_diag_cor","local_diag_loc","normality","randomization","sampling_grids","umf","var"]
done_list=["normality","sampling_grids","lct","lct_cor"]
done = {}
alias_list = {
   "avg": "avg",
   "diag": "diag",
   "dirac": "contour_centered",
   "lct": "contour_positive",
   "lct_cor": "lct_cor",
   "local_diag_cor": "contour_positive",
   "local_diag_loc": "contour_positive",
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
   # Make output directory
   testfig = args.testdata + "/" + args.test + "/fig"
   if not os.path.exists(testfig):
      os.mkdir(testfig)

   # Create jpg figures
   for f in sorted(os.listdir(os.path.join(args.testdata, args.test))):
      if os.path.isfile(os.path.join(args.testdata, args.test, f)) and f.find("test_" + args.mpi + "-" + args.omp)==0:
         suffix = f.split("test_" + args.mpi + "-" + args.omp + "_")[1][:-(len(f.rsplit(".")[-1])+1)]
         for plot in plot_list:
            if suffix.find(plot)==0 and not done[plot]:
               print("Calling " + plot + " in " + args.test + " (" + args.mpi + "-" + args.omp + ")")
               plot_alias = alias_list[plot]
               module = __import__(plot_alias)
               func = getattr(module, plot_alias)
               if plot in done_list:
                  func(args.testdata, args.test, args.mpi, args.omp, plot, testfig)
                  done[plot] = True
               else:
                  func(args.testdata, args.test, args.mpi, args.omp, suffix, testfig)

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
