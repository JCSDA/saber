# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

option(BUILD_BUMP_DOC "Build BUMP documentation" OFF)

if (BUILD_BUMP_DOC)
    add_subdirectory( bump )
endif (BUILD_BUMP_DOC)
