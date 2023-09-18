# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

add_subdirectory( lib )
list( TRANSFORM bump_lib_src_files PREPEND lib/ )
list( APPEND bump_src_files_list ${bump_lib_src_files} )

list(APPEND bump_src_files_list
BUMPParameters.h
NICAS.cc
NICAS.h
PsiChiToUV.cc
PsiChiToUV.h
StdDev.cc
StdDev.h
VerticalBalance.cc
VerticalBalance.h
)

set( bump_src_files ${bump_src_files_list} PARENT_SCOPE )
