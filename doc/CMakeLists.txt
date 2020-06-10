# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

ecbuild_add_option( FEATURE DOCS
                    DEFAULT OFF
                    DESCRIPTION "Generate documentation" )

if(HAVE_DOCS)
#    add_subdirectory( bump )
endif (HAVE_DOCS)
