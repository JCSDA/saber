# (C) Copyright 2023 Meteorologisk Institutt.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND fastlam_src_files_list
# FastLAM block
FastLAM.cc
FastLAM.h
FastLAMParametersBase.h
InterpElement.h
LayerBase.cc
LayerBase.h
LayerHalo.cc
LayerHalo.h
LayerRC.cc
LayerRC.h
)

set( fastlam_src_files ${fastlam_src_files_list} PARENT_SCOPE )
