# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Find Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
    # Generate the mainpage.h file
    message(STATUS "Build SABER documentation with saber_doc_overview.sh")
    add_custom_target( saber_doc_overview ALL
        COMMAND ${CMAKE_BINARY_DIR}/bin/saber_doc_overview.sh
                ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Build SABER documentation with saber_doc_overview.sh"
        VERBATIM )

    # Run Doxygen
    message(STATUS "Build SABER documentation with Doxygen")
    include(GNUInstallDirs)

    # Set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # Request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    # Note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen_saber ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating SABER documentation with Doxygen"
        DEPENDS saber_doc_overview
                ${PROJECT_NAME}
        VERBATIM )
    set( HTML_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/html )
    install(DIRECTORY ${HTML_BUILD_DIR} DESTINATION ${CMAKE_INSTALL_FULL_DATADIR}/${PROJECT_NAME}/doc COMPONENT Documentation)
else()
    message(WARNING "Build SABER documentation failure to find Doxygen")
endif()
