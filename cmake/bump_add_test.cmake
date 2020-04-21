function(bump_add_test)
    set(options VALGRIND)
    set(single_value_args TESTNAME MPI OMP)
    cmake_parse_arguments(_p "${options}" "${single_value_args}" "${multi_value_args}" ${_FIRST_ARG} ${ARGN})

    # Set default values
    if (NOT _p_MPI)
        set(_p_MPI 1)
    endif()
    if (NOT _p_OMP)
        set(_p_OMP 1)
    endif()

    # Add tests
    if (_p_VALGRIND)
        ecbuild_add_test(TARGET       test_${_p_TESTNAME}_${_p_MPI}-${_p_OMP}_valgrind
                         TYPE SCRIPT
                         COMMAND      ${CMAKE_SOURCE_DIR}/saber/tools/bump_valgrind.sh
                         ARGS         ${CMAKE_BINARY_DIR}/bin/bump.x
                                      testinput/${_p_TESTNAME}_${mpi}-${omp}.yaml
                                      testoutput
                         DEPENDS      bump.x
                         TEST_DEPENDS get_saber_data
                                      get_saber_ref)
     else()   
        ecbuild_add_test(TARGET       test_${_p_TESTNAME}_${_p_MPI}-${_p_OMP}
                         TYPE SCRIPT
                         COMMAND      ${CMAKE_SOURCE_DIR}/saber/tools/bump_wrapper.sh
                         ARGS         ${MPIEXEC}
                                      ${MPIEXEC_NUMPROC_FLAG}
                                      ${_p_MPI}
                                      ${_p_OMP}
                                      ${CMAKE_BINARY_DIR}/bin/bump.x
                                      testinput/${_p_TESTNAME}_${mpi}-${omp}.yaml
                                      testoutput
                                      ${CMAKE_BINARY_DIR}/bin/bump_compare.sh
                                      ${_p_TESTNAME}
                         DEPENDS      bump.x
                         TEST_DEPENDS get_saber_data
                                      get_saber_ref)
    endif()
endfunction()
