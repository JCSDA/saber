# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND bump_src_files
  bump/bump.F90
  bump/type_model.F90
  bump/type_timer.F90
)

ecbuild_add_executable( TARGET  bump.x
                        SOURCES ${bump_src_files}
                        LIBS    saber
                       )
