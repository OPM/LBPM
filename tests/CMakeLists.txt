# Copy files for the tests
ADD_LBPM_EXECUTABLE( lbpm_permeability_simulator )
ADD_LBPM_EXECUTABLE( lbpm_color_simulator )
ADD_LBPM_EXECUTABLE( lbpm_sphere_pp )
ADD_LBPM_EXECUTABLE( lbpm_random_pp )
ADD_LBPM_EXECUTABLE( lbpm_segmented_pp )
ADD_LBPM_EXECUTABLE( lbpm_segmented_decomp )
ADD_LBPM_EXECUTABLE( lbpm_disc_pp )
ADD_LBPM_EXECUTABLE( lbpm_captube_pp )
ADD_LBPM_EXECUTABLE( lbpm_squaretube_pp )
ADD_LBPM_EXECUTABLE( lbpm_BlobAnalysis )
ADD_LBPM_EXECUTABLE( TestBubble )
ADD_LBPM_EXECUTABLE( ComponentLabel )
ADD_LBPM_EXECUTABLE( ColorToBinary )
ADD_LBPM_EXECUTABLE( BlobAnalysis )
ADD_LBPM_EXECUTABLE( BlobIdentify )
ADD_LBPM_EXECUTABLE( BlobIdentifyParallel )
#ADD_LBPM_EXECUTABLE( BlobAnalyzeParallel )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cylindertest ${CMAKE_CURRENT_BINARY_DIR}/cylindertest COPYONLY )

# Add the tests
ADD_LBPM_TEST( pmmc_cylinder )
ADD_LBPM_TEST( TestBubble )
ADD_LBPM_TEST( TestTorus )
ADD_LBPM_TEST( TestInterfaceSpeed )
ADD_LBPM_TEST( TestSphereCurvature )
ADD_LBPM_TEST_1_2_4( TestTwoPhase )
ADD_LBPM_TEST_1_2_4( TestBlobIdentify )
ADD_LBPM_TEST_PARALLEL( TestTwoPhase 8 )
ADD_LBPM_TEST_PARALLEL( TestBlobAnalyze 8 )
ADD_LBPM_TEST_PARALLEL( TestSegDist 8 )
ADD_LBPM_TEST_PARALLEL( TestCommD3Q19 8 )
ADD_LBPM_TEST_PARALLEL( TestMassConservationD3Q7 1 )
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


