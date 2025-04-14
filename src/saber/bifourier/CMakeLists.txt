# (C) Copyright 2025 Meteorologisk Institutt.
#

list(APPEND bifourier_src_files_list
# Biperiodization block
Biperiodization.h
Biperiodization.cc
BiperiodizationImpl.h
BiperiodizationImpl.cc
)

set( bifourier_src_files ${bifourier_src_files_list} PARENT_SCOPE )
