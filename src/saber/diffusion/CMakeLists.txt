# (C) Copyright 2024-2024 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND diffusion_src_files_list
Diffusion.cc
Diffusion.h
DiffusionParameters.h
)

set( diffusion_src_files ${diffusion_src_files_list} PARENT_SCOPE )
