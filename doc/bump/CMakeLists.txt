# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Run architecture script
add_custom_target( bump_architecture ALL
   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/architecture.ksh ${CMAKE_CURRENT_SOURCE_DIR}
   COMMENT "Generating BUMP architecture documentation"
   VERBATIM )

# Run autodoc script
add_custom_target( bump_autodoc ALL
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/autodoc.ksh ${CMAKE_CURRENT_SOURCE_DIR}
   COMMENT "Generating BUMP auto-documentation"
   VERBATIM )

# Run cloc_report script
add_custom_target( bump_cloc_report ALL
   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cloc_report.ksh ${CMAKE_CURRENT_SOURCE_DIR}
   COMMENT "Generating BUMP cloc report"
   VERBATIM )

# Run pdf script
add_custom_target( bump_pdf_documentation ALL
   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/pdf.ksh ${CMAKE_CURRENT_SOURCE_DIR}
   COMMENT "Generating BUMP pdf documentation"
   VERBATIM )
