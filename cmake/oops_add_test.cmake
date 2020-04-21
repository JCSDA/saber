function(oops_add_test)
  set( options COMPARE )
  set( single_value_args TESTNAME YAMLNAME EXENAME MPI CTOL IDIF RUN_FILE REF_FILE )
  set( multi_value_args  DEPENDS TEST_DEPENDS)
  cmake_parse_arguments( _p "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  # Set default values

  if (NOT _p_MPI)
    set(_p_MPI 1)
  endif()

  if (NOT _p_OMP)
    set(_p_OMP 1)
  endif()

  if (NOT _p_CTOL)
    set( _p_CTOL 0.0)  # Max relative diff
  endif()

  if (NOT _p_IDIF)
    set( _p_IDIF 0)  # Max diff
  endif()

  ecbuild_add_test( TARGET test_${_p_TESTNAME}
                    TYPE SCRIPT
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_wrapper.sh
                    ARGS ${MPIEXEC}
                         ${MPIEXEC_NUMPROC_FLAG}
                         ${_p_MPI}
                         ${_p_OMP}
                         ${CMAKE_BINARY_DIR}/bin/${_p_EXENAME}
                         ${_p_YAMLNAME}
                         ${_p_RUN_FILE}
                         ${CMAKE_BINARY_DIR}/bin/compare.py
                         ${_p_REF_FILE}
                         ${_p_CTOL}
                         ${_p_IDIF}
                    DEPENDS ${_p_DEPENDS}
                    TEST_DEPENDS ${_p_TEST_DEPENDS})
endfunction(oops_add_test)
