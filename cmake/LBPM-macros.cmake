# Copy an example folder
MACRO( INSTALL_EXAMPLE EXAMPLE )
    INSTALL( DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE}" DESTINATION "${LBPM_INSTALL_DIR}/example" )
ENDMACRO()


# Create an example test
CONFIGURE_FILE( "${${PROJ}_SOURCE_DIR}/cmake/CompareOutput.cmake" "${${PROJ}_BUILD_DIR}/CompareOutput.cmake" COPYONLY )
MACRO( TEST_EXAMPLE EXAMPLE EXEFILE PROCS ${ARGN} )
    SET( EXAMPLE_DIR "${CMAKE_CURRENT_BINARY_DIR}/${EXAMPLE}" )
    # Copy the example directory
    ADD_CUSTOM_TARGET(
        ${EXAMPLE} ALL
        ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE}" "${EXAMPLE_DIR}"
        DEPENDS ${EXEFILE}
    )
    # Create a wrapper script to run the test and copy the output to ${EXAMPLE}.out
    SET( FILENAME "${EXAMPLE_DIR}/run-${EXAMPLE}" )
    FILE(WRITE  "${FILENAME}" "# This is a automatically generated file to run example--${EXAMPLE}\n" )
    FILE(APPEND "${FILENAME}" "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} \"${LBPM_INSTALL_DIR}/bin/${EXEFILE}\" ${ARGN} 2>&1 | tee ${EXAMPLE}.out\n\n" )
    # Create the test to run the example
    SET( TESTNAME example--${EXAMPLE} )
    EXECUTE_PROCESS(COMMAND chmod 755 "${FILENAME}")
    ADD_TEST( 
        NAME ${TESTNAME}
        WORKING_DIRECTORY "${EXAMPLE_DIR}"
        COMMAND "${FILENAME}"
    )
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${PROCS} )
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES RESOURCE_LOCK ${EXEFILE} )
    # Create a test that checks the output against the data in EXAMPLE/OutputAns.txt
    IF ( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE}/ExampleOutput.txt" )
        ADD_TEST( 
            NAME ${TESTNAME}-output
            WORKING_DIRECTORY "${EXAMPLE_DIR}"
            COMMAND ${CMAKE_COMMAND} -DTEST=${EXAMPLE}.out -DGOLD=ExampleOutput.txt -P "${${PROJ}_BUILD_DIR}/CompareOutput.cmake"
        )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS 1 )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES DEPENDS ${TESTNAME} )
    ENDIF()
ENDMACRO()

