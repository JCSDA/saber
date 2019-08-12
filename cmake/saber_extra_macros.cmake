# (C) Copyright 2019 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Extra macros to eliminate repetition

# macro to prepend a prefix with relative path
function( PREPEND var prefix )
  set( listVar "" )
  foreach( f ${ARGN} )
    list( APPEND listVar "${prefix}/${f}" )
  endforeach( f )
  set( ${var} "${listVar}" PARENT_SCOPE )
endfunction( PREPEND )

# macro to create a symlink from src to dst
function(CREATE_SYMLINK src dst)
    foreach (FILENAME ${ARGN})
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${FILENAME} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK)
