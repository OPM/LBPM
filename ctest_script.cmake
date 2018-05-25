# ctest script for building, running, and submitting the test results 
# Usage:  ctest -s script,build
#   build = debug / optimized / valgrind
# Note: this test will use use the number of processors defined in the variable N_PROCS,
#   the enviornmental variable N_PROCS, or the number of processors availible (if not specified)

# Set platform specific variables
SITE_NAME( HOSTNAME )
IF( ${HOSTNAME} STREQUAL "lap0086227" )
    SET( COVERAGE_COMMAND /usr/bin/gcov )
    SET( VALGRIND_COMMAND /usr/bin/valgrind )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles" )
    SET( CC "mpicc" )
    SET( CXX "mpicxx" )
    SET( C_FLAGS "-DCBUB" )
    SET( CXX_FLAGS "-DCBUB" )
    SET( MPIEXEC "mpirun" )
ELSEIF( ${HOSTNAME} MATCHES "vayu" ) 
    SET( COVERAGE_COMMAND /usr/bin/gcov )
    SET( VALGRIND_COMMAND /usr/local/bin/valgrind )
    SET( CUDA_FLAGS "--use_fast_math -Xptxas=-v -arch=sm_20" )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles" )
    SET( CC "mpicc" )
    SET( CXX "mpicxx" )
    SET( C_FLAGS "-DCBUB" )
    SET( CXX_FLAGS "-DCBUB" )
    SET( MPIEXEC "mpirun" )
ELSEIF( ${HOSTNAME} MATCHES "titan.*" ) 
    SET( COVERAGE_COMMAND "" )
    SET( VALGRIND_COMMAND "" )
    SET( CUDA_FLAGS "-arch sm_35" )
    SET( CUDA_HOST_COMPILER "/usr/bin/g++" )
    SET( CTEST_CMAKE_GENERATOR "Unix Makefiles" )
    SET( CC "cc" )
    SET( CXX "CC" )
    SET( C_FLAGS "-DCBUB" )
    SET( CXX_FLAGS "-DCBUB" )
    SET( MPIEXEC "aprun" )
    SET( N_PROCS 16 )
ELSE()
    MESSAGE( FATAL_ERROR "Unknown host: ${HOSTNAME}" )
ENDIF()


# Get the source directory based on the current directory
SET( LBPM_DIR "${CMAKE_CURRENT_LIST_DIR}" )


# Check that we specified the build type to run
IF( NOT CTEST_SCRIPT_ARG )
    MESSAGE(FATAL_ERROR "No build specified: ctest -S /path/to/script,build (debug/optimized/valgrind")
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "debug" )
    SET( CTEST_BUILD_NAME "LBPM-WIA-debug" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( CTEST_COVERAGE_COMMAND ${COVERAGE_COMMAND} )
    SET( ENABLE_GCOV "true" )
    SET( USE_VALGRIND FALSE )
    SET( USE_CUDA FALSE )
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "debug-cuda" )
    SET( CTEST_BUILD_NAME "LBPM-WIA-debug-cuda" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( CTEST_COVERAGE_COMMAND ${COVERAGE_COMMAND} )
    SET( ENABLE_GCOV "true" )
    SET( USE_VALGRIND FALSE )
    SET( USE_CUDA TRUE )
ELSEIF( (${CTEST_SCRIPT_ARG} STREQUAL "optimized") OR (${CTEST_SCRIPT_ARG} STREQUAL "opt") )
    SET( CTEST_BUILD_NAME "LBPM-WIA-opt" )
    SET( CMAKE_BUILD_TYPE "Release" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND FALSE )
    SET( USE_CUDA FALSE )
ELSEIF( (${CTEST_SCRIPT_ARG} STREQUAL "optimized-cuda") OR (${CTEST_SCRIPT_ARG} STREQUAL "opt-cuda") )
    SET( CTEST_BUILD_NAME "LBPM-WIA-opt-cuda" )
    SET( CMAKE_BUILD_TYPE "Release" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND FALSE )
    SET( USE_CUDA TRUE )
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "valgrind" )
    SET( CTEST_BUILD_NAME "LBPM-WIA-valgrind" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND TRUE )
    SET( USE_CUDA FALSE )
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "valgrind-cuda" )
    SET( CTEST_BUILD_NAME "LBPM-WIA-valgrind-cuda" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND TRUE )
    SET( USE_CUDA TRUE )
ELSE()
    MESSAGE(FATAL_ERROR "Invalid build (${CTEST_SCRIPT_ARG}): ctest -S /path/to/script,build (debug/opt/valgrind")
ENDIF()
IF ( NOT CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
ENDIF()


# Set the number of processors
IF( NOT DEFINED N_PROCS )
    SET( N_PROCS $ENV{N_PROCS} )
ENDIF()
IF( NOT DEFINED N_PROCS )
    SET(N_PROCS 1)
    # Linux:
    SET(cpuinfo_file "/proc/cpuinfo")
    IF(EXISTS "${cpuinfo_file}")
        FILE(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
        list(LENGTH procs N_PROCS)
    ENDIF()
    # Mac:
    IF(APPLE)
        find_program(cmd_sys_pro "system_profiler")
        if(cmd_sys_pro)
            execute_process(COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info)
            STRING(REGEX REPLACE "^.*Total Number Of Cores: ([0-9]+).*$" "\\1" N_PROCS "${info}")
        ENDIF()
    ENDIF()
    # Windows:
    IF(WIN32)
        SET(N_PROCS "$ENV{NUMBER_OF_PROCESSORS}")
    ENDIF()
ENDIF()


# Set basic variables
SET( CTEST_PROJECT_NAME "LBPM-WIA" )
SET( CTEST_SOURCE_DIRECTORY "${LBPM_DIR}" )
SET( CTEST_BINARY_DIRECTORY "." )
SET( CTEST_DASHBOARD "Nightly" )
SET( CTEST_TEST_TIMEOUT 300 )
SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 500 )
SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 500 )
SET( CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 10000 )
SET( CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 10000 )
SET( NIGHTLY_START_TIME "18:00:00 EST" )
SET( CTEST_NIGHTLY_START_TIME "22:00:00 EST" )
SET( CTEST_COMMAND "\"${CTEST_EXECUTABLE_NAME}\" -D ${CTEST_DASHBOARD}" )
SET( CTEST_BUILD_COMMAND "make -i -j ${N_PROCS} install" )


# Set valgrind options
#SET (VALGRIND_COMMAND_OPTIONS "--tool=memcheck --leak-check=yes --track-fds=yes --num-callers=50 --show-reachable=yes --track-origins=yes --malloc-fill=0xff --free-fill=0xfe --suppressions=${LBPM_DIR}/ValgrindSuppresionFile" )
SET( VALGRIND_COMMAND_OPTIONS  "--tool=memcheck --leak-check=yes --track-fds=yes --num-callers=50 --show-reachable=yes --suppressions=${LBPM_DIR}/ValgrindSuppresionFile" )
IF ( USE_VALGRIND )
    SET( MEMORYCHECK_COMMAND ${VALGRIND_COMMAND} )
    SET( MEMORYCHECKCOMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECK_COMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECKCOMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS ${VALGRIND_COMMAND_OPTIONS} )
    SET( CTEST_MEMORYCHECKCOMMAND_OPTIONS  ${VALGRIND_COMMAND_OPTIONS} )
ENDIF()


# Clear the binary directory and create an initial cache
CTEST_EMPTY_BINARY_DIRECTORY (${CTEST_BINARY_DIRECTORY})
FILE(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "CTEST_TEST_CTEST:BOOL=1")

# Set the configure options
SET( CTEST_OPTIONS )
SET( CTEST_OPTIONS "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}" )
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_C_COMPILER:PATH=${CC};-DCMAKE_C_FLAGS='${C_FLAGS}';" )
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DCMAKE_CXX_COMPILER:PATH=${CXX};-DCMAKE_CXX_FLAGS='${CXX_FLAGS}'" )
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DMPI_COMPILER:BOOL=true;-DMPIEXEC=${MPIEXEC};-DUSE_EXT_MPI_FOR_SERIAL_TESTS:BOOL=true")
IF ( USE_CUDA )
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DUSE_CUDA:BOOL=true;-DCUDA_NVCC_FLAGS='${CUDA_FLAGS}';-DCUDA_HOST_COMPILER=${CUDA_HOST_COMPILER}" )
ELSE()
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DUSE_CUDA:BOOL=false" )
ENDIF()
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DLDLIBS:STRING=\"${LDLIBS}\"" )
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DENABLE_GCOV:BOOL=${ENABLE_GCOV}" )

# Configure and run the tests
SET( CTEST_SITE ${HOSTNAME} )
CTEST_START("${CTEST_DASHBOARD}")
CTEST_UPDATE()
CTEST_CONFIGURE(
    BUILD   ${CTEST_BINARY_DIRECTORY}
    SOURCE  ${CTEST_SOURCE_DIRECTORY}
    OPTIONS "${CTEST_OPTIONS}"
)
CTEST_BUILD()
IF ( USE_VALGRIND_MATLAB )
    CTEST_TEST( INCLUDE MATLAB--test_hello_world  PARALLEL_LEVEL ${N_PROCS} )
ELSEIF ( USE_VALGRIND )
    CTEST_MEMCHECK( EXCLUDE procs   PARALLEL_LEVEL ${N_PROCS} )
ELSE()
    # CTEST_TEST( EXCLUDE WEEKLY  PARALLEL_LEVEL ${N_PROCS} )
    CTEST_TEST( PARALLEL_LEVEL ${N_PROCS} )
ENDIF()
IF( CTEST_COVERAGE_COMMAND )
    CTEST_COVERAGE()
ENDIF()


# Submit the results to oblivion
SET( CTEST_DROP_METHOD "http" )
SET( CTEST_DROP_SITE "oblivion.engr.colostate.edu" )
SET( CTEST_DROP_LOCATION "/CDash/submit.php?project=LBPM-WIA" )
SET( CTEST_DROP_SITE_CDASH TRUE )
SET( DROP_SITE_CDASH TRUE )
CTEST_SUBMIT()


# Clean up
# exec_program("make distclean")


