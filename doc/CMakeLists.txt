# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Run dirlist script
add_custom_target( saber_doc_overview ALL
    COMMAND ${CMAKE_BINARY_DIR}/bin/saber_doc_overview.sh ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating SABER directory list"
    VERBATIM )
