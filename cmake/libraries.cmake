INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Find_TIMER.cmake" )
INCLUDE(CheckCCompilerFlag)
INCLUDE(CheckCXXCompilerFlag)


FUNCTION( CONFIGURE_LINE_COVERAGE )
    SET( COVERAGE_FLAGS )
    SET( COVERAGE_LIBS )
    IF ( ENABLE_GCOV )
        SET( COVERAGE_FLAGS -DUSE_GCOV )
        SET( CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage )
        CHECK_CXX_SOURCE_COMPILES( "int main() { return 0;}" profile-arcs )
        IF ( profile-arcs )
            SET( COVERAGE_FLAGS "${COVERAGE_FLAGS} -fprofile-arcs -ftest-coverage" )
            SET( COVERAGE_LIBS ${COVERAGE_LIBS} -fprofile-arcs )
        ENDIF()
        SET( CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS} -lgcov" )
        CHECK_CXX_SOURCE_COMPILES( "int main() { return 0;}" lgcov )
        IF ( lgcov )
            SET( COVERAGE_LIBS -lgcov ${COVERAGE_LIBS} )
        ENDIF()
        MESSAGE("Enabling coverage:")
        MESSAGE("   COVERAGE_FLAGS = ${COVERAGE_FLAGS}")
        MESSAGE("   COVERAGE_LIBS = ${COVERAGE_LIBS}")
        ADD_DEFINITIONS( ${COVERAGE_FLAGS} )
        SET( COVERAGE_FLAGS ${COVERAGE_FLAGS} PARENT_SCOPE )
        SET( COVERAGE_LIBS  ${COVERAGE_LIBS}  PARENT_SCOPE )
    ENDIF()
ENDFUNCTION()


# Macro to configure MIC
MACRO( CONFIGURE_MIC )
    CHECK_ENABLE_FLAG( USE_MIC 0 )
    IF ( USE_MIC )
        ADD_DEFINITIONS( "-D USE_MIC" )
    ENDIF()
ENDMACRO()


# Macro to find and configure the MPI libraries
MACRO( CONFIGURE_MPI )
    # Determine if we want to use MPI
    CHECK_ENABLE_FLAG(USE_MPI 1 )
    IF ( USE_MPI )
        # Check if we specified the MPI directory
        IF ( MPI_DIRECTORY )
            # Check the provided MPI directory for include files
            VERIFY_PATH( "${MPI_DIRECTORY}" )
            IF ( EXISTS "${MPI_DIRECTORY}/include/mpi.h" )
                SET( MPI_INCLUDE_PATH "${MPI_DIRECTORY}/include" )
            ELSEIF ( EXISTS "${MPI_DIRECTORY}/Inc/mpi.h" )
                SET( MPI_INCLUDE_PATH "${MPI_DIRECTORY}/Inc" )
            ELSE()
                MESSAGE( FATAL_ERROR "mpi.h not found in ${MPI_DIRECTORY}/include" )
            ENDIF ()
            INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH} )
            SET ( MPI_INCLUDE ${MPI_INCLUDE_PATH} )
            # Set MPI libraries
            IF ( ${CMAKE_SYSTEM_NAME} STREQUAL "Windows" )
                FIND_LIBRARY( MSMPI_LIB     NAMES msmpi     PATHS "${MPI_DIRECTORY}/Lib/x64"    NO_DEFAULT_PATH )
                FIND_LIBRARY( MSMPI_LIB     NAMES msmpi     PATHS "${MPI_DIRECTORY}/Lib/amd64"  NO_DEFAULT_PATH )
                FIND_LIBRARY( MSMPIFEC_LIB  NAMES msmpifec  PATHS "${MPI_DIRECTORY}/Lib/x64"    NO_DEFAULT_PATH )
                FIND_LIBRARY( MSMPIFEC_LIB  NAMES msmpifec  PATHS "${MPI_DIRECTORY}/Lib/amd64"  NO_DEFAULT_PATH )
                FIND_LIBRARY( MSMPIFMC_LIB  NAMES msmpifmc  PATHS "${MPI_DIRECTORY}/Lib/x64"    NO_DEFAULT_PATH )
                FIND_LIBRARY( MSMPIFMC_LIB  NAMES msmpifmc  PATHS "${MPI_DIRECTORY}/Lib/amd64"  NO_DEFAULT_PATH )
                SET( MPI_LIBRARIES ${MSMPI_LIB} ${MSMPIFEC_LIB} ${MSMPIFMC_LIB} )
            ENDIF()
            # Set the mpi executable
            IF ( MPIEXEC ) 
                # User specified the MPI command directly, use as is
            ELSEIF ( MPIEXEC_CMD )
                # User specified the name of the MPI executable
                SET ( MPIEXEC ${MPI_DIRECTORY}/bin/${MPIEXEC_CMD} )
                IF ( NOT EXISTS ${MPIEXEC} )
                    MESSAGE( FATAL_ERROR "${MPIEXEC_CMD} not found in ${MPI_DIRECTORY}/bin" )
                ENDIF ()
            ELSE ()
                # Search for the MPI executable in the current directory
                FIND_PROGRAM ( MPIEXEC  NAMES mpiexec mpirun lamexec  PATHS ${MPI_DIRECTORY}/bin  NO_DEFAULT_PATH )
                IF ( NOT MPIEXEC )
                    MESSAGE( FATAL_ERROR "Could not locate mpi executable" )
                ENDIF()
            ENDIF ()
            # Set MPI flags
            IF ( NOT MPIEXEC_NUMPROC_FLAG )
                SET( MPIEXEC_NUMPROC_FLAG "-np" )
            ENDIF()
        ELSEIF ( MPI_COMPILER )
            # The mpi compiler should take care of everything
            IF ( MPI_INCLUDE )
                INCLUDE_DIRECTORIES( ${MPI_INCLUDE} )
            ENDIF()
        ELSE()
            # Perform the default search for MPI
            INCLUDE ( FindMPI )
            IF ( NOT MPI_FOUND )
                MESSAGE( "  MPI_INCLUDE = ${MPI_INCLUDE}" )
                MESSAGE( "  MPI_LINK_FLAGS = ${MPI_LINK_FLAGS}" )
                MESSAGE( "  MPI_LIBRARIES = ${MPI_LIBRARIES}" )
                MESSAGE( FATAL_ERROR "Did not find MPI" )
            ENDIF ()
            INCLUDE_DIRECTORIES( "${MPI_INCLUDE_PATH}" )
            SET( MPI_INCLUDE "${MPI_INCLUDE_PATH}" )
        ENDIF()
        # Check if we need to use MPI for serial tests
        CHECK_ENABLE_FLAG( USE_MPI_FOR_SERIAL_TESTS 0 )
        # Set defaults if they have not been set
        IF ( NOT MPIEXEC )
            SET( MPIEXEC mpirun )
        ENDIF()
        IF ( NOT MPIEXEC_NUMPROC_FLAG )
            SET( MPIEXEC_NUMPROC_FLAG "-np" )
        ENDIF()
        # Set the definitions
        ADD_DEFINITIONS( "-DUSE_MPI" )  
        MESSAGE( "Using MPI" )
        MESSAGE( "  MPIEXEC = ${MPIEXEC}" )
        MESSAGE( "  MPIEXEC_NUMPROC_FLAG = ${MPIEXEC_NUMPROC_FLAG}" )
        MESSAGE( "  MPI_INCLUDE = ${MPI_INCLUDE}" )
        MESSAGE( "  MPI_LINK_FLAGS = ${MPI_LINK_FLAGS}" )
        MESSAGE( "  MPI_LIBRARIES = ${MPI_LIBRARIES}" )
    ELSE()
        SET( USE_MPI_FOR_SERIAL_TESTS 0 )
        SET( MPIEXEC "" )
        SET( MPIEXEC_NUMPROC_FLAG "" )
        SET( MPI_INCLUDE "" )
        SET( MPI_LINK_FLAGS "" )
        SET( MPI_LIBRARIES "" )
        MESSAGE( "Not using MPI, all parallel tests will be disabled" )
    ENDIF()
ENDMACRO()


# Macro to find and configure hdf5
MACRO ( CONFIGURE_HDF5 )
    CHECK_ENABLE_FLAG( USE_HDF5 0 )
    IF ( USE_HDF5 )
        # Check if we specified the silo directory
        IF ( HDF5_DIRECTORY )
            VERIFY_PATH ( ${HDF5_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${HDF5_DIRECTORY}/include )
            SET ( HDF5_INCLUDE ${HDF5_DIRECTORY}/include )
            FIND_LIBRARY ( HDF5_LIB    NAMES hdf5    PATHS ${HDF5_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HDF5_HL_LIB NAMES hdf5_hl PATHS ${HDF5_DIRECTORY}/lib  NO_DEFAULT_PATH )
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for hdf5 is not yet supported.  Use -D HDF5_DIRECTORY=" )
        ENDIF()
        SET ( HDF5_LIBS
            ${HDF5_HL_LIB}
            ${HDF5_LIB}
        )
        ADD_DEFINITIONS ( -DUSE_HDF5 )  
        MESSAGE( "Using hdf5" )
        MESSAGE( "   ${HDF5_LIB}" )
    ENDIF()
ENDMACRO()


# Macro to find and configure netcdf
MACRO( CONFIGURE_NETCDF )
    CHECK_ENABLE_FLAG( USE_NETCDF 0 )
    IF ( USE_NETCDF )
        SET( USE_HDF5 1 )
        CONFIGURE_HDF5()
        IF ( NETCDF_DIRECTORY )
            VERIFY_PATH ( ${NETCDF_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${NETCDF_DIRECTORY}/include )
            SET ( NETCDF_INCLUDE ${NETCDF_DIRECTORY}/include )
            FIND_LIBRARY( NETCDF_NETCDF_LIB    NAMES netcdf    PATHS ${NETCDF_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY( NETCDF_HDF5_LIB      NAMES hdf5      PATHS ${NETCDF_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY( NETCDF_HL_LIB        NAMES hl        PATHS ${NETCDF_DIRECTORY}/lib  NO_DEFAULT_PATH )
            IF ( NOT NETCDF_NETCDF_LIB )
                MESSAGE( FATAL_ERROR "Did not find library for netcdf" )
            ENDIF()
            SET ( NETCDF_LIBS ${NETCDF_NETCDF_LIB} )
            IF ( NETCDF_HDF5_LIB )
                SET ( NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_HDF5_LIB} )
            ENDIF()
            IF ( NETCDF_HL_LIB )
                SET ( NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_HL_LIB} )
            ENDIF()
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for netcdf is not yet supported.  Use -D NETCDF_DIRECTORY=" )
        ENDIF()
        SET( EXTERNAL_LIBS ${NETCDF_LIBS} ${HDF5_LIBS} ${EXTERNAL_LIBS} )
        ADD_DEFINITIONS ( -DUSE_NETCDF )
        MESSAGE( "Using netcdf" )
        MESSAGE( "   ${NETCDF_LIBS}" )
    ENDIF()
ENDMACRO()


# Macro to find and configure the silo libraries
MACRO ( CONFIGURE_SILO )
    # Determine if we want to use silo
    CHECK_ENABLE_FLAG( USE_EXT_SILO 0 )
    IF ( USE_SILO )
        SET( USE_HDF5 1 )
        CONFIGURE_HDF5()
        # Check if we specified the silo directory
        IF ( SILO_DIRECTORY )
            VERIFY_PATH ( ${SILO_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${SILO_DIRECTORY}/include )
            SET ( SILO_INCLUDE ${SILO_DIRECTORY}/include )
            FIND_LIBRARY ( SILO_LIB  NAMES siloh5  PATHS ${SILO_DIRECTORY}/lib  NO_DEFAULT_PATH )
        ELSE()
            MESSAGE( "Default search for silo is not yet supported")
            MESSAGE( "Use -D SILO_DIRECTORY=" FATAL_ERROR)
        ENDIF()
        SET ( SILO_LIBS
            ${SILO_LIB}
        )
        SET( EXTERNAL_LIBS ${SILO_LIBS} ${HDF5_LIBS} ${EXTERNAL_LIBS} )
        ADD_DEFINITIONS ( -DUSE_SILO )  
        MESSAGE( "Using silo" )
        MESSAGE( "   ${SILO_LIB}" )
    ENDIF ()
ENDMACRO()


# Macro to configure system-specific libraries and flags
MACRO( CONFIGURE_SYSTEM )
    # First check/set the compile mode
    IF( NOT CMAKE_BUILD_TYPE )
        MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE is not set")
    ENDIF()
    # Remove extra library links
    # Get the compiler
    SET_COMPILER_FLAGS()
    CHECK_ENABLE_FLAG( USE_STATIC 0 )
    # Add system dependent flags
    MESSAGE("System is: ${CMAKE_SYSTEM_NAME}")
    IF ( ${CMAKE_SYSTEM_NAME} STREQUAL "Windows" )
        # Windows specific system libraries
        SET( SYSTEM_PATHS "C:/Program Files (x86)/Microsoft SDKs/Windows/v7.0A/Lib/x64" 
                          "C:/Program Files (x86)/Microsoft Visual Studio 8/VC/PlatformSDK/Lib/AMD64" 
                          "C:/Program Files (x86)/Microsoft Visual Studio 12.0/Common7/Packages/Debugger/X64" )
        FIND_LIBRARY( PSAPI_LIB    NAMES Psapi    PATHS ${SYSTEM_PATHS}  NO_DEFAULT_PATH )
        FIND_LIBRARY( DBGHELP_LIB  NAMES DbgHelp  PATHS ${SYSTEM_PATHS}  NO_DEFAULT_PATH )
        FIND_LIBRARY( DBGHELP_LIB  NAMES DbgHelp )
        IF ( PSAPI_LIB ) 
            ADD_DEFINITIONS( -D PSAPI )
            SET( SYSTEM_LIBS ${PSAPI_LIB} )
        ENDIF()
        IF ( DBGHELP_LIB ) 
            ADD_DEFINITIONS( -D DBGHELP )
            SET( SYSTEM_LIBS ${DBGHELP_LIB} )
        ELSE()
            MESSAGE( WARNING "Did not find DbgHelp, stack trace will not be availible" )
        ENDIF()
        MESSAGE("System libs: ${SYSTEM_LIBS}")
    ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Linux" )
        # Linux specific system libraries
        SET( SYSTEM_LIBS -lz -lpthread -ldl )
        IF ( NOT USE_STATIC )
            # Try to add rdynamic so we have names in backtrace
            SET( CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS} -rdynamic" )
            CHECK_CXX_SOURCE_COMPILES( "int main() { return 0;}" rdynamic )
            IF ( rdynamic )
                SET( SYSTEM_LDFLAGS ${SYSTEM_LDFLAGS} -rdynamic )
            ENDIF()
        ENDIF()
        # Try to add -fPIC
        SET( CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS} -fPIC" )
        CHECK_CXX_SOURCE_COMPILES( "int main() { return 0;}" fPIC )
        IF ( fPIC )
            SET( SYSTEM_LDFLAGS ${SYSTEM_LDFLAGS} -fPIC )
            SET( SYSTEM_LDFLAGS ${SYSTEM_LDFLAGS} -fPIC )
        ENDIF()
        IF ( USING_GCC )
            SET( SYSTEM_LIBS ${SYSTEM_LIBS} -lgfortran )
            SET(CFLAGS_EXTRA   " ${CFLAGS_EXTRA} -fPIC" )
            SET(CXXFLAGS_EXTRA " ${CXXFLAGS_EXTRA} -fPIC" )
        ENDIF()
    ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin" )
        # Max specific system libraries
        SET( SYSTEM_LIBS -lz -lpthread -ldl )
    ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Generic" )
        # Generic system libraries
    ELSE()
        MESSAGE( FATAL_ERROR "OS not detected" )
    ENDIF()
    # Add the static flag if necessary
    IF ( USE_STATIC )
        SET_STATIC_FLAGS()
    ENDIF()
    # Print some flags
    MESSAGE( "LDLIBS = ${LDLIBS}" )
ENDMACRO ()


# Macro to configure LBPM specific options
MACRO ( CONFIGURE_LBPM )
    # Set the maximum number of processors for the tests
    IF ( NOT TEST_MAX_PROCS )
        SET( TEST_MAX_PROCS 32 )
    ENDIF()
    # Add the correct paths to rpath in case we build shared libraries
    IF ( USE_STATIC )
        SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE )
        SET( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
    ELSE()
        SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
        SET( CMAKE_BUILD_WITH_INSTALL_RPATH TRUE )
        SET( CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} "${TIMER_DIRECTORY}" "${LBPM_INSTALL_DIR}/lib" )
    ENDIF()
ENDMACRO ()
