#!/bin/ksh
#----------------------------------------------------------------------
# Korn shell script: architecture
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

if type "dot" > /dev/null ; then
   # Directories
   doc=$1
   mkdir -p ${doc}
   architecture=$1/architecture
   mkdir -p ${architecture}

   # Build diagrams
   filename=setup_online
   echo "       File: "${filename}
cat<<EOFDIAGRAM > ${architecture}/${filename}.txt
digraph callgraph {
  "External code" -> "Initialize external namelist parameters"
  "Initialize external namelist parameters" -> "Initialize namelist :: nam%init";
  "Initialize external namelist parameters" -> "Set external namelist parameters (e.g. from a JSON file)";
  "External code" -> "Online setup :: bump%setup_online";
  "Online setup :: bump%setup_online" -> "Initialize MPL :: mpl%init";
  "Online setup :: bump%setup_online" -> "If namelist reading required";
  "If namelist reading required" -> "Read namelist :: nam%read";
  "If namelist reading required" -> "Broadcast namelist :: nam%bcast";
  "Online setup :: bump%setup_online" ->  "Set internal namelist parameters :: nam%setup_internal";
  "Online setup :: bump%setup_online" -> "Initialize listing :: mpl%init_listing";
  "Online setup :: bump%setup_online" -> "Generic setup :: setup_generic";
  "Online setup :: bump%setup_online" -> "Get coordinates from arguments :: geom%setup";
  "Online setup :: bump%setup_online" -> "Initialize geometry :: geom%setup";
  "Online setup :: bump%setup_online" -> "If fields regridding required";
  "If fields regridding required" -> "Initialize fields regridding :: io%grid_init";
  "Online setup :: bump%setup_online" -> "Initialize block parameters :: bpar%alloc";
  "Online setup :: bump%setup_online" -> "Initialize ensemble 1 :: ens1%alloc";
  "Online setup :: bump%setup_online" -> "Initialize ensemble 2 :: ens2%alloc";
  "Online setup :: bump%setup_online" -> "If observation operator required";
  "If observation operator required" -> "Initialize observations locations :: obsop%from";
}
EOFDIAGRAM
   cat ${architecture}/${filename}.txt | dot -Nshape="box" -Grankdir="LR" -Tsvg > ${architecture}/${filename}.svg
   rm -f ${architecture}/${filename}.txt

   filename=run_online
   echo "       File: "${filename}
cat<<EOFDIAGRAM > ${architecture}/${filename}.txt
digraph callgraph {
  "External code" -> "If ensembles required";
  "If ensembles required" -> "Add ensemble members :: bump%add_member";
  "External code" -> "If setting parameters required";
  "If setting parameters required" -> "Set parameters :: bump%set_parameter";
  "External code" -> "Run drivers :: bump%run_drivers";
}
EOFDIAGRAM
   cat ${architecture}/${filename}.txt | dot -Nshape="box" -Grankdir="LR" -Tsvg > ${architecture}/${filename}.svg
   rm -f ${architecture}/${filename}.txt

   filename=get_parameter
   echo "       File: "${filename}
cat<<EOFDIAGRAM > ${architecture}/${filename}.txt
digraph callgraph {
  "External code" -> "Get parameters :: bump%get_parameter";
}
EOFDIAGRAM
   cat ${architecture}/${filename}.txt | dot -Nshape="box" -Grankdir="LR" -Tsvg > ${architecture}/${filename}.svg
   rm -f ${architecture}/${filename}.txt

   filename=apply_vbal
   echo "       File: "${filename}
cat<<EOFDIAGRAM > ${architecture}/${filename}.txt
digraph callgraph {
  "External code" -> "Apply vertical balance :: bump%apply_vbal";
  "External code" -> "Apply vertical balance inverse :: bump%apply_vbal_inv";
  "External code" -> "Apply vertical balance adjoint :: bump%apply_vbal_ad";
  "External code" -> "Apply vertical balance inverse adjoint :: bump%apply_vbal_inv_ad";
}
EOFDIAGRAM
   cat ${architecture}/${filename}.txt | dot -Nshape="box" -Grankdir="LR" -Tsvg > ${architecture}/${filename}.svg
   rm -f ${architecture}/${filename}.txt

   filename=apply_nicas
   echo "       File: "${filename}
cat<<EOFDIAGRAM > ${architecture}/${filename}.txt
digraph callgraph {
  "External code" -> "Get control variable size :: bump_get_cv_size";
  "External code" -> "Apply NICAS square-root :: bump%apply_nicas_sqrt";
  "External code" -> "Apply NICAS square-root adjoint :: bump%apply_nicas_sqrt_ad";
  "External code" -> "Apply NICAS :: bump%apply_nicas";
}
EOFDIAGRAM
   cat ${architecture}/${filename}.txt | dot -Nshape="box" -Grankdir="LR" -Tsvg > ${architecture}/${filename}.svg
   rm -f ${architecture}/${filename}.txt

   filename=apply_obsop
   echo "       File: "${filename}
cat<<EOFDIAGRAM > ${architecture}/${filename}.txt
digraph callgraph {
  "External code" -> "Apply observation operator :: bump%apply_obsop";
  "External code" -> "Apply observation operator adjoint :: bump%apply_obsop_ad";
}
EOFDIAGRAM
   cat ${architecture}/${filename}.txt | dot -Nshape="box" -Grankdir="LR" -Tsvg > ${architecture}/${filename}.svg
   rm -f ${architecture}/${filename}.txt

   filename=deallocation
   echo "       File: "${filename}
cat<<EOFDIAGRAM > ${architecture}/${filename}.txt
digraph callgraph {
  "External code" -> "Deallocate BUMP object :: bump%dealloc";
}
EOFDIAGRAM
   cat ${architecture}/${filename}.txt | dot -Nshape="box" -Grankdir="LR" -Tsvg > ${architecture}/${filename}.svg
   rm -f ${architecture}/${filename}.txt

   filename=drivers
   echo "       File: "${filename}
cat<<EOFDIAGRAM > ${architecture}/${filename}.txt
digraph callgraph {
  "Run drivers :: bump%run_drivers" -> "Finalize ensemble 1 :: ens1%remove_mean";
  "Run drivers :: bump%run_drivers" -> "Finalize ensemble 2 :: ens2%remove_mean";
  "Run drivers :: bump%run_drivers" -> "If new vertical balance required";
  "If new vertical balance required" -> "Run vertical balance driver :: vbal%run_vbal"
  "Run drivers :: bump%run_drivers" -> "If vertical balance reading required";
  "If vertical balance reading required" -> "Read vertical balance :: vbal%read"
  "Run drivers :: bump%run_drivers" -> "If vertical balance test required";
  "If vertical balance test required" -> "Run vertical balance tests driver :: vbal%run_vbal_tests"
  "Run drivers :: bump%run_drivers" -> "If new HDIAG required";
  "If new HDIAG required" -> "Run HDIAG driver :: hdiag%run_hdiag"
  "If new HDIAG required" -> "Copy HDIAG into C matrix :: cmat%from_hdiag"
  "Run drivers :: bump%run_drivers" -> "If new LCT required";
  "If new LCT required" -> "Run LCT driver :: lct%run_lct"
  "If new LCT required" -> "Copy LCT into C matrix :: cmat%from_lct"
  "Run drivers :: bump%run_drivers" -> "If C matrix loading required";
  "If C matrix loading required" -> "Read C matrix :: cmat%read";
  "Run drivers :: bump%run_drivers" -> "If C matrix setting from namelist required";
  "If C matrix setting from namelist required" -> "Copy namelist support radii into C matrix :: cmat%read";
  "Run drivers :: bump%run_drivers" -> "If parameters have been set from external code";
  "If parameters have been set from external code" -> "Get C matrix from OOPS :: cmat%from_oops"
  "If parameters have been set from external code" -> "Setup C matrix sampling :: cmat%from_oops"
  "If parameters have been set from external code" -> "Write C matrix :: cmat%write"
  "Run drivers :: bump%run_drivers" -> "If new NICAS required";
  "If new NICAS required" -> "Run NICAS driver :: nicas%run_nicas";
  "Run drivers :: bump%run_drivers" -> "If NICAS loading required";
  "If NICAS loading required" -> "Read NICAS parameters :: nicas%read";
  "Run drivers :: bump%run_drivers" -> "If NICAS test required";
  "If NICAS test required" -> "Run NICAS tests driver :: nicas%run_nicas_tests";
  "Run drivers :: bump%run_drivers" -> "If new observation operator required";
  "If new observation operator required" -> "Run observation operator driver :: obsop%run_obsop";
  "Run drivers :: bump%run_drivers" -> "If observation operator loading required";
  "If observation operator loading required" -> "Read observation operator parameters :: obsop%read";
  "Run drivers :: bump%run_drivers" -> "If observation operator test required";
  "If observation operator test required" -> "Run observation operator tests driver :: obsop%run_obsop_tests";
  "Run drivers :: bump%run_drivers" -> "Close listings";
}
EOFDIAGRAM
   cat ${architecture}/${filename}.txt | dot -Nshape="box" -Grankdir="LR" -Tsvg > ${architecture}/${filename}.svg
   rm -f ${architecture}/${filename}.txt
else
   echo "dot not found: no architecture diagrams"
fi
