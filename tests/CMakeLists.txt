# Copy files for the tests
INSTALL_LBPM_EXE( lb2_Color_wia_mpi )
INSTALL_LBPM_EXE( lbpm_color_simulator )
INSTALL_LBPM_EXE( lbpm_sphere_pp )
INSTALL_LBPM_EXE( lbpm_disc_pp )
INSTALL_LBPM_EXE( lbpm_BlobAnalysis )
INSTALL_LBPM_EXE( TestBubble )
INSTALL_LBPM_EXE( BasicSimulator )
INSTALL_LBPM_EXE( BlobAnalysis )
INSTALL_LBPM_EXE( BlobIdentify )
INSTALL_LBPM_EXE( BlobIdentifyParallel )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cylindertest ${CMAKE_CURRENT_BINARY_DIR}/cylindertest COPYONLY )

# Add the tests
ADD_LBPM_TEST( pmmc_cylinder )
ADD_LBPM_TEST( TestBubble )
ADD_LBPM_TEST( TestCylinderAreas )
ADD_LBPM_TEST( TestInterfaceSpeed )
ADD_LBPM_TEST( TestSphereCurvature )
ADD_LBPM_TEST( TestContactAngle )
ADD_LBPM_TEST_1_2_4( TestTwoPhase )
ADD_LBPM_TEST_PARALLEL( TestTwoPhase 8 )
ADD_LBPM_TEST_1_2_4( testCommunication )
ADD_LBPM_TEST_1_2_4( testUtilities )
ADD_LBPM_TEST( TestWriter )

# Sample test that will run with 1, 2, and 4 processors, failing with 4 or more procs
ADD_LBPM_TEST_1_2_4( hello_world )
SET_TESTS_PROPERTIES( hello_world        PROPERTIES ENVIRONMENT "MPICH_RDMA_ENABLED_CUDA=0")
IF ( USE_MPI ) 
    SET_TESTS_PROPERTIES( hello_world_2procs PROPERTIES ENVIRONMENT "MPICH_RDMA_ENABLED_CUDA=0")
    SET_TESTS_PROPERTIES( hello_world_4procs PROPERTIES ENVIRONMENT "MPICH_RDMA_ENABLED_CUDA=0")
ENDIF()


# Add CPU/GPU specific test
IF ( USE_CUDA )
    ADD_SUBDIRECTORY( gpu )
ELSE()
    ADD_SUBDIRECTORY( cpu )
ENDIF()


