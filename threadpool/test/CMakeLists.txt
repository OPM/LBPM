# Add thread pool tests
ADD_LBPM_TEST( test_atomic )
ADD_LBPM_TEST( test_atomic_list )
SET_TESTS_PROPERTIES ( test_atomic PROPERTIES FAIL_REGULAR_EXPRESSION ".*FAILED.*" PROCESSORS 64 )
ADD_LBPM_TEST_THREAD_MPI( test_thread_pool 1 4 )
ADD_LBPM_TEST_THREAD_MPI( test_thread_pool 2 4 )
ADD_LBPM_TEST_THREAD_MPI( test_thread_pool 4 4 )
SET_PROPERTY( TEST test_thread_pool_1procs_4threads APPEND PROPERTY RUN_SERIAL 1 )
IF ( USE_MPI )
    SET_PROPERTY( TEST test_thread_pool_2procs_4threads APPEND PROPERTY RUN_SERIAL 1 )
    SET_PROPERTY( TEST test_thread_pool_4procs_4threads APPEND PROPERTY RUN_SERIAL 1 )
ENDIF()




