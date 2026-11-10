# EMC/NCEP/NOAA 2024
# Explicit ordering required: parent modules must compile before their submodules
# so Intel Fortran can find the correct .smod files.
#
include(GNUInstallDirs)

# Dependency tree:
#   mg_parameter.f90  ->  jp_pbfil.f90, mg_domain.f90, mg_domain_loc.f90,
#                         mg_mppstuff.f90
#   mg_intstate.f90   ->  mg_bocos.f90, mg_entrymod.f90, mg_filtering.f90,
#                         mg_generations.f90, mg_interpolate.f90, mg_transfer.f90
#
# CMake cannot infer SUBMODULE(parent) compile-time ordering within a single
# target on older generators/compilers. Split into two targets so that
# add_dependencies enforces a hard TARGET-level barrier: every file in
# mgbflib_parents (including mg_parameter.f90 and mg_intstate.f90) is fully
# compiled — and their .mod/.smod files written — before a single file in
# mgbflib starts.

set( mgbf_lib_parent_src_files
  # standalone utility modules (no submodule dependencies)
  mgbf_kind.f90
  mgbf_util.f90
  jp_pkind.f90
  jp_pkind2.f90
  jp_pietc.f90
  jp_pietc_s.f90
  jp_pmat.f90
  jp_pmat4.f90
  jp_pbfil2.f90
  jp_pbfil3.f90
  phint.f90
  phint1.f90
  mg_timers.f90
  mg_tools.f90
  # parent modules (must come before their submodules)
  mg_parameter.f90
  mg_intstate.f90
)

set( mgbf_lib_submodule_src_files
  # submodules of mg_parameter
  jp_pbfil.f90
  mg_domain.f90
  mg_domain_loc.f90
  mg_mppstuff.f90
  # submodules of mg_intstate
  mg_entrymod.f90
  mg_bocos.f90
  mg_filtering.f90
  mg_generations.f90
  mg_interpolate.f90
  mg_transfer.f90
)

set(mgbf_lib_module_dir ${CMAKE_CURRENT_BINARY_DIR}/modules)

# --- Parent modules target ---------------------------------------------------
# Compiles all plain modules first; Intel ifort writes .mod AND .smod here.
add_library(mgbflib_parents STATIC ${mgbf_lib_parent_src_files})
set_target_properties(mgbflib_parents PROPERTIES
    Fortran_MODULE_DIRECTORY ${mgbf_lib_module_dir})
target_include_directories(mgbflib_parents PUBLIC
    $<BUILD_INTERFACE:${mgbf_lib_module_dir}>)
target_link_libraries(mgbflib_parents PUBLIC MPI::MPI_Fortran)

# --- Submodule target --------------------------------------------------------
# add_dependencies creates a hard barrier: mgbflib_parents must be 100 %
# complete (all .mod/.smod written) before any file in mgbflib is touched.
add_library(mgbflib STATIC ${mgbf_lib_submodule_src_files})
add_dependencies(mgbflib mgbflib_parents)
set_target_properties(mgbflib PROPERTIES
    Fortran_MODULE_DIRECTORY ${mgbf_lib_module_dir})
target_include_directories(mgbflib PUBLIC
    $<BUILD_INTERFACE:${mgbf_lib_module_dir}>)
target_link_libraries(mgbflib PUBLIC mgbflib_parents MPI::MPI_Fortran)

install(TARGETS mgbflib mgbflib_parents
        EXPORT saber-targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
