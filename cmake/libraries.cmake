MACRO ( CONFIGURE_LINE_COVERAGE )
    SET ( COVERAGE_LIBS )
    IF ( ENABLE_GCOV )
        ADD_DEFINITIONS ( -fprofile-arcs -ftest-coverage )
        SET ( COVERAGE_LIBS -lgcov -fprofile-arcs )
    ENDIF ()
ENDMACRO ()


# Macro to find and configure the MPI libraries
MACRO ( CONFIGURE_MPI )
    # Determine if we want to use MPI
    CHECK_ENABLE_FLAG(USE_MPI 1 )
    IF ( USE_MPI )
        # Check if we specified the MPI directory
        IF ( MPI_DIRECTORY )
            # Check the provided MPI directory for include files and the mpi executable
            VERIFY_PATH ( ${MPI_DIRECTORY} )
            SET ( MPI_INCLUDE_PATH ${MPI_DIRECTORY}/include )
            VERIFY_PATH ( ${MPI_INCLUDE_PATH} )
            IF ( NOT EXISTS ${MPI_INCLUDE_PATH}/mpi.h )
                MESSAGE ( FATAL_ERROR "mpi.h not found in ${MPI_INCLUDE_PATH}/include" )
            ENDIF ()
            INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH} )
            SET ( MPI_INCLUDE ${MPI_INCLUDE_PATH} )
            IF ( MPIEXEC ) 
                # User specified the MPI command directly, use as is
            ELSEIF ( MPIEXEC_CMD )
                # User specified the name of the MPI executable
                SET ( MPIEXEC ${MPI_DIRECTORY}/bin/${MPIEXEC_CMD} )
                IF ( NOT EXISTS ${MPIEXEC} )
                    MESSAGE ( FATAL_ERROR "${MPIEXEC_CMD} not found in ${MPI_DIRECTORY}/bin" )
                ENDIF ()
            ELSE ()
                # Search for the MPI executable in the current directory
                FIND_PROGRAM ( MPIEXEC  NAMES mpiexec mpirun lamexec  PATHS ${MPI_DIRECTORY}/bin  NO_DEFAULT_PATH )
                IF ( NOT MPIEXEC )
                    MESSAGE ( FATAL_ERROR "Could not locate mpi executable" )
                ENDIF()
            ENDIF ()
            # Set MPI flags
            IF ( NOT MPIEXEC_NUMPROC_FLAG )
                SET( MPIEXEC_NUMPROC_FLAG "-np" )
            ENDIF()
        ELSEIF ( MPI_COMPILER )
            # The mpi compiler should take care of everything
        ELSE()
            # Perform the default search for MPI
            INCLUDE ( FindMPI )
            IF ( NOT MPI_FOUND )
                MESSAGE ( FATAL_ERROR "Did not find MPI" )
            ENDIF ()
            INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH} )
            SET ( MPI_INCLUDE ${MPI_INCLUDE_PATH} )
        ENDIF()
        # Check if we need to use MPI for serial tests
        CHECK_ENABLE_FLAG( USE_MPI_FOR_SERIAL_TESTS 0 )
        # Set the definitions
        ADD_DEFINITIONS ( "-D USE_MPI" )  
        MESSAGE ( "Using MPI" )
        MESSAGE ( "  MPIEXEC = ${MPIEXEC}" )
        MESSAGE ( "  MPIEXEC_NUMPROC_FLAG = ${MPIEXEC_NUMPROC_FLAG}" )
        MESSAGE ( "  MPI_INCLUDE = ${MPI_INCLUDE}" )
        MESSAGE ( "  MPI_LINK_FLAGS = ${MPI_LINK_FLAGS}" )
        MESSAGE ( "  MPI_LIBRARIES = ${MPI_LIBRARIES}" )
    ELSE()
        SET( USE_MPI_FOR_SERIAL_TESTS 0 )
        SET( MPIEXEC "" )
        SET( MPIEXEC_NUMPROC_FLAG "" )
        SET( MPI_INCLUDE "" )
        SET( MPI_LINK_FLAGS "" )
        SET( MPI_LIBRARIES "" )
        MESSAGE ( "Not using MPI, all parallel tests will be disabled" )
    ENDIF()
ENDMACRO ()


# Macro to configure system-specific libraries and flags
MACRO ( CONFIGURE_SYSTEM )
    # First check/set the compile mode
    IF( NOT CMAKE_BUILD_TYPE )
        MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE is not set")
    ENDIF()
    # Remove extra library links
    # Get the compiler
    SET_COMPILER ()
    # Add the static flag if necessary
    CHECK_ENABLE_FLAG( USE_EXT_STATIC 0 )
    IF ( USE_EXT_STATIC )
        SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")    # Add static flag
        SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-static")  # Add static flag
    ENDIF()
    # Add system dependent flags
    IF ( USING_MICROSOFT )
        # Windows specific system libraries
        SET( SYSTEM_PATHS "C:/Program Files (x86)/Microsoft SDKs/Windows/v7.0A/Lib/x64" 
                          "C:/Program Files (x86)/Microsoft Visual Studio 8/VC/PlatformSDK/Lib/AMD64" )
        FIND_LIBRARY ( PSAPI_LIB    NAMES Psapi    PATHS ${SYSTEM_PATHS}  NO_DEFAULT_PATH )
        FIND_LIBRARY ( DBGHELP_LIB  NAMES DbgHelp  PATHS ${SYSTEM_PATHS}  NO_DEFAULT_PATH )
        SET( SYSTEM_LIBS ${PSAPI_LIB} ${DBGHELP_LIB} )
        MESSAGE("System libs: ${SYSTEM_LIBS}")
    ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Linux" )
        # Linux specific system libraries
        CHECK_C_COMPILER_FLAG("-rdynamic" RESULT)
        IF(RESULT)
            SET( SYSTEM_LIBS "-lpthread -lz -ldl -rdynamic" )
        ELSE()
            SET( SYSTEM_LIBS "-lpthread -lz -ldl" )
        ENDIF()
        IF ( USING_GCC )
            SET( SYSTEM_LIBS ${SYSTEM_LIBS} "-lgfortran" )
            SET(CMAKE_C_FLAGS   " ${CMAKE_C_FLAGS} -fPIC" )
            SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -fPIC" )
        ENDIF()
    ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin" )
        # Max specific system libraries
        SET( SYSTEM_LIBS "-lz -ldl" )
    ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Generic" )
        # Generic system libraries
    ELSE()
        MESSAGE( FATAL_ERROR "OS not detected" )
    ENDIF()
    # Set the compile flags based on the build
    SET_COMPILE_FLAGS()
ENDMACRO ()


# Macro to configure AtomicModel-specific options
MACRO ( CONFIGURE_LBPM )
    # Set the maximum number of processors for the tests
    IF ( NOT TEST_MAX_PROCS )
        SET( TEST_MAX_PROCS 32 )
    ENDIF()
ENDMACRO ()
