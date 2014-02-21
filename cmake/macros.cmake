INCLUDE(CheckCSourceCompiles)


# Macro to print all variables
MACRO( PRINT_ALL_VARIABLES )
    GET_CMAKE_PROPERTY(_variableNames VARIABLES)
    FOREACH(_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    ENDFOREACH()
ENDMACRO()


# Add a package to the LBPM library
MACRO( ADD_LBPM_LIBRARY PACKAGE )
    #INCLUDE_DIRECTORIES ( ${LBPM_INSTALL_DIR}/include/${PACKAGE} )
    ADD_SUBDIRECTORY( ${PACKAGE} )
ENDMACRO()


# Add an LBPM executable
MACRO( ADD_LBPM_EXECUTABLE PACKAGE )
    ADD_SUBDIRECTORY( ${PACKAGE} )
ENDMACRO()


# Initialize a package
MACRO (BEGIN_PACKAGE_CONFIG PACKAGE)
    SET( HEADERS "" )
    SET( CXXSOURCES "" )
    SET( CSOURCES "" )
    SET( FSOURCES "" )
    SET( M4FSOURCES "" )
    SET( CUSOURCES "" )
    SET( SOURCES "" )
    SET( CURPACKAGE ${PACKAGE} )
ENDMACRO ()


# Find the source files
MACRO (FIND_FILES)
    # Find the C/C++ headers
    SET( T_HEADERS "" )
    FILE( GLOB T_HEADERS "*.h" "*.hh" "*.hpp" "*.I" )
    # Find the C sources
    SET( T_CSOURCES "" )
    FILE( GLOB T_CSOURCES "*.c" )
    # Find the C++ sources
    SET( T_CXXSOURCES "" )
    FILE( GLOB T_CXXSOURCES "*.cc" "*.cpp" "*.cxx" "*.C" )
    # Find the C++ sources
    SET( T_CUSOURCES "" )
    FILE( GLOB T_CUSOURCES "*.cu" )
    # Add all found files to the current lists
    SET( HEADERS ${HEADERS} ${T_HEADERS} )
    SET( CXXSOURCES ${CXXSOURCES} ${T_CXXSOURCES} )
    SET( CSOURCES ${CSOURCES} ${T_CSOURCES} )
    SET( CUSOURCES ${CUSOURCES} ${T_CUSOURCES} )
    SET( SOURCES ${SOURCES} ${T_CXXSOURCES} ${T_CSOURCES} )
ENDMACRO()


# Find the source files
MACRO (FIND_FILES_PATH IN_PATH)
    # Find the C/C++ headers
    SET( T_HEADERS "" )
    FILE( GLOB T_HEADERS "${IN_PATH}/*.h" "${IN_PATH}/*.hh" "${IN_PATH}/*.hpp" "${IN_PATH}/*.I" )
    # Find the C sources
    SET( T_CSOURCES "" )
    FILE( GLOB T_CSOURCES "${IN_PATH}/*.c" )
    # Find the C++ sources
    SET( T_CXXSOURCES "" )
    FILE( GLOB T_CXXSOURCES "${IN_PATH}/*.cc" "${IN_PATH}/*.cpp" "${IN_PATH}/*.cxx" "${IN_PATH}/*.C" )
    # Find the CUDA sources
    SET( T_CUSOURCES "" )
    FILE( GLOB T_CUSOURCES "${IN_PATH}/*.cu" )
    # Add all found files to the current lists
    SET( HEADERS ${HEADERS} ${T_HEADERS} )
    SET( CXXSOURCES ${CXXSOURCES} ${T_CXXSOURCES} )
    SET( CSOURCES ${CSOURCES} ${T_CSOURCES} )
    SET( CUSOURCES ${CUSOURCES} ${T_CUSOURCES} )
    SET( SOURCES ${SOURCES} ${T_CXXSOURCES} ${T_CSOURCES} ${T_CUSOURCES} )
ENDMACRO()


# Add a subdirectory
MACRO( ADD_PACKAGE_SUBDIRECTORY SUBDIR )
    SET( FULLSUBDIR ${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR} )
    FIND_FILES_PATH( ${SUBDIR} )
    FILE( GLOB HFILES RELATIVE ${FULLSUBDIR} ${SUBDIR}/*.h ${SUBDIR}/*.hh ${SUBDIR}/*.hpp  ${SUBDIR}/*.I )
    FOREACH( HFILE ${HFILES} )
        CONFIGURE_FILE( ${FULLSUBDIR}/${HFILE} ${LBPM_INSTALL_DIR}/include/${SUBDIR}/${HFILE} COPYONLY )
        INCLUDE_DIRECTORIES( ${FULLSUBDIR} )
    ENDFOREACH()
    ADD_SUBDIRECTORY( ${SUBDIR} )
ENDMACRO()


# Install a package
MACRO( INSTALL_LBPM_TARGET PACKAGE )
    # Find all files in the current directory
    FIND_FILES()
    # Copy the header files to the include path
    FILE( GLOB HFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.hh ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp ${CMAKE_CURRENT_SOURCE_DIR}/*.I )
    FOREACH( HFILE ${HFILES} )
        #CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/${HFILE} ${LBPM_INSTALL_DIR}/include/${CURPACKAGE}/${HFILE} COPYONLY )
        CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/${HFILE} ${LBPM_INSTALL_DIR}/include/${HFILE} COPYONLY )
    ENDFOREACH()
    # Configure the CUDA files
    IF ( CUSOURCES )
        CUDA_COMPILE( CUOBJS ${CUSOURCES} )
    ENDIF()
    # Add the library
    ADD_LIBRARY( ${PACKAGE} ${LIB_TYPE} ${SOURCES} ${CUOBJS} )
    SET( TEST_DEP_LIST ${PACKAGE} ${TEST_DEP_LIST} )
    TARGET_LINK_LIBRARIES( ${PACKAGE} ${COVERAGE_LIBS} ${SYSTEM_LIBS} ${LDLIBS} )
    TARGET_LINK_LIBRARIES( ${PACKAGE} ${LAPACK_LIBS} ${BLAS_LIBS} )
    IF ( USE_MPI )
        TARGET_LINK_LIBRARIES( ${PACKAGE} ${MPI_LIBRARIES} )
    ENDIF()
    TARGET_LINK_LIBRARIES( ${PACKAGE} ${COVERAGE_LIBS} ${SYSTEM_LIBS} ${LDLIBS} )
    # Install the package
    INSTALL( TARGETS ${PACKAGE} DESTINATION ${LBPM_INSTALL_DIR}/lib )
    INSTALL( FILES ${HFILES} DESTINATION ${LBPM_INSTALL_DIR}/include )
    # Clear the sources
    SET( HEADERS "" )
    SET( CSOURCES "" )
    SET( CXXSOURCES "" )
ENDMACRO()


# Macro to verify that a variable has been set
MACRO( VERIFY_VARIABLE VARIABLE_NAME )
    IF( NOT ${VARIABLE_NAME} )
        MESSAGE( FATAL_ERROR "PLease set: " ${VARIABLE_NAME} )
    ENDIF()
ENDMACRO()


# Macro to verify that a path has been set
MACRO( VERIFY_PATH PATH_NAME )
    IF( NOT EXISTS "${PATH_NAME}" )
        MESSAGE( FATAL_ERROR "Path does not exist: " ${PATH_NAME} )
    ENDIF()
ENDMACRO()


# Macro to identify the compiler
MACRO( SET_COMPILER )
    # SET the C/C++ compiler
    IF( CMAKE_COMPILE_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
        SET( USING_GCC TRUE )
        MESSAGE("Using gcc")
    ELSEIF( MSVC OR MSVC_IDE OR MSVC60 OR MSVC70 OR MSVC71 OR MSVC80 OR CMAKE_COMPILER_2005 OR MSVC90 OR MSVC10 )
        IF( NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows" )
            MESSAGE( FATAL_ERROR "Using microsoft compilers on non-windows system?" )
        ENDIF()
        SET( USING_MICROSOFT TRUE )
        MESSAGE("Using Microsoft")
    ELSEIF( (${CMAKE_C_COMPILER_ID} MATCHES "Intel") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Intel") ) 
        SET(USING_ICC TRUE)
        MESSAGE("Using icc")
    ELSEIF( ${CMAKE_C_COMPILER_ID} MATCHES "PGI")
        SET(USING_PGCC TRUE)
        MESSAGE("Using pgCC")
    ELSE()
        SET(USING_DEFAULT TRUE)
        MESSAGE("${CMAKE_C_COMPILER_ID}")
        MESSAGE("Unknown C/C++ compiler, default flags will be used")
    ENDIF()
ENDMACRO()


# Macro to set the proper warnings
MACRO ( SET_WARNINGS )
  IF ( USING_GCC )
    ## Add gcc specific compiler options
    ##    -Wno-reorder:  warning: "" will be initialized after "" when initialized here
    #SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -Wall ") 
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall ")
  ELSEIF ( USING_MICROSOFT )
    # Add Microsoft specifc compiler options
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} /D _SCL_SECURE_NO_WARNINGS /D _CRT_SECURE_NO_WARNINGS /D _ITERATOR_DEBUG_LEVEL=0" )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _SCL_SECURE_NO_WARNINGS /D _CRT_SECURE_NO_WARNINGS /D _ITERATOR_DEBUG_LEVEL=0" )
  ELSEIF ( USING_ICC )
    ## Add Intel specifc compiler options
    #SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -Wall" )
    #SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall" )
    ## Disable warnings that I think are irrelavent (may need to be revisited)
    #SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -wd383 -wd981" )
    #SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -wd383 -wd981" )
  ELSEIF ( USING_DEFAULT )
    ## Add default compiler options
    #SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -Wall")
    #SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall")
  ENDIF ()
ENDMACRO ()


# Macro to add user compile flags
MACRO( ADD_USER_FLAGS )
    SET(CMAKE_C_FLAGS   " ${CMAKE_C_FLAGS} ${CFLAGS} ${LDFLAGS}" )
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} ${LDFLAGS}" )
ENDMACRO()


# Macro to set the flags for debug mode
MACRO( SET_COMPILE_FLAGS )
    ADD_USER_FLAGS()
    IF ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
        IF ( NOT DISABLE_GXX_DEBUG )
            SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -D_GLIBCXX_DEBUG" )
            SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -D_GLIBCXX_DEBUG" )
            SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -D_GLIBCXX_DEBUG_PEDANTIC" )
            SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -D_GLIBCXX_DEBUG_PEDANTIC" )
        ENDIF ()
        IF ( USING_MICROSOFT )
            SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -DDEBUG /DEBUG /Od" )
            SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -DDEBUG /DEBUG /Od" )
            SET(CONFIGURATION Debug )
        ELSE()
            SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -DDEBUG -g -O0" )
            SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -DDEBUG -g -O0" )
        ENDIF()
    ELSEIF ( ${CMAKE_BUILD_TYPE} STREQUAL "Release" )
        IF ( USING_MICROSOFT )
            SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} /O3" )
            SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} /O3" )
            SET(CONFIGURATION Release )
        ELSE()
            SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -O3" )
            SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -O3" )
        ENDIF()
    ELSE()
        MESSAGE(FATAL_ERROR "Unknown build type: ${CMAKE_BUILD_TYPE}")
    ENDIF()
    SET_WARNINGS()
ENDMACRO()




# Macro to add the dependencies and libraries to an executable
MACRO( ADD_LBPM_EXE_DEP EXE )
    # Add the package dependencies
    IF( LBPM_TEST_LIB_EXISTS )
        ADD_DEPENDENCIES ( ${EXE} ${PACKAGE_TEST_LIB} )
        TARGET_LINK_LIBRARIES ( ${EXE} ${PACKAGE_TEST_LIB} )
    ENDIF()
    # Add the executable to the dependencies of check and build-test
    ADD_DEPENDENCIES( check ${EXE} )
    ADD_DEPENDENCIES( build-test ${EXE} )
    # Add the libraries
    TARGET_LINK_LIBRARIES( ${EXE} ${LBPM_LIBS} )
    # Add external libraries
    TARGET_LINK_LIBRARIES( ${EXE} ${EXTERNAL_LIBS} )
    IF ( USE_MPI )
        TARGET_LINK_LIBRARIES( ${EXE} ${MPI_LINK_FLAGS} ${MPI_LIBRARIES} )
    ENDIF()
    TARGET_LINK_LIBRARIES( ${EXE} ${LAPACK_LIBS} ${BLAS_LIBS} )
    TARGET_LINK_LIBRARIES( ${EXE} ${COVERAGE_LIBS} ${SYSTEM_LIBS} ${LDLIBS} )
ENDMACRO()


# Add a executable
MACRO( INSTALL_LBPM_EXE EXE )
    SET( SOURCES ${EXE}.cpp )
    ADD_EXECUTABLE( ${EXE} ${SOURCES} )
    ADD_LBPM_EXE_DEP( ${EXE} )
    INSTALL( TARGETS ${EXE} DESTINATION ${LBPM_INSTALL_DIR}/bin )
ENDMACRO()


# Add a provisional test
FUNCTION( ADD_LBPM_PROVISIONAL_TEST EXEFILE )
    # Check if we actually want to add the test
    SET( EXCLUDE_TESTS_FROM_ALL 0 )
    # Check if test has already been added
    GET_TARGET_PROPERTY(tmp ${EXEFILE} LOCATION)
    IF ( NOT tmp )
        # The target has not been added
        SET( CXXFILE ${EXEFILE}.cpp )
        SET( TESTS_SO_FAR ${TESTS_SO_FAR} ${EXEFILE} )
        IF ( NOT EXCLUDE_TESTS_FROM_ALL )
            ADD_EXECUTABLE( ${EXEFILE} ${CXXFILE} )
        ELSE()
            ADD_EXECUTABLE( ${EXEFILE} EXCLUDE_FROM_ALL ${CXXFILE} )
        ENDIF()
        ADD_LBPM_EXE_DEP( ${EXEFILE} )
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/$(Configuration)/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/$(OutDir)/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSE()
        # We are trying to add 2 different tests with the same name
        MESSAGE ( "Existing test: ${tmp}" )
        MESSAGE ( "New test:      ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}" )
        MESSAGE ( FATAL_ERROR "Trying to add 2 different tests with the same name" )
    ENDIF()
ENDFUNCTION()


# Macro to create the test name
MACRO( CREATE_TEST_NAME TEST ${ARGN} )
    SET( TESTNAME "${TEST}" )
    FOREACH( tmp ${ARGN} )
        SET( TESTNAME "${TESTNAME}--${tmp}")
    endforeach()
    # STRING(REGEX REPLACE "--" "-" TESTNAME ${TESTNAME} )
ENDMACRO()


# Add a executable as a test
FUNCTION( ADD_LBPM_TEST EXEFILE ${ARGN} )
    ADD_LBPM_PROVISIONAL_TEST ( ${EXEFILE} )
    CREATE_TEST_NAME( ${EXEFILE} ${ARGN} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    IF ( USE_EXT_MPI_FOR_SERIAL_TESTS )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${EXE} ${ARGN} )
    ELSE()
        ADD_TEST( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
    ENDIF()
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION ".*FAILED.*" PROCESSORS 1 )
ENDFUNCTION()


# Add a executable as a weekly test
FUNCTION( ADD_LBPM_WEEKLY_TEST EXEFILE PROCS ${ARGN} )
    ADD_LBPM_PROVISIONAL_TEST ( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    IF( ${PROCS} STREQUAL "1" )
        CREATE_TEST_NAME( "${EXEFILE}_WEEKLY" ${ARGN} )
        IF( USE_EXT_MPI_FOR_SERIAL_TESTS )
            ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${EXE} ${ARGN} )
        ELSE()
            ADD_TEST( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
        ENDIF()
    ELSEIF( USE_MPI AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_WEEKLY" ${ARGN} )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
    ENDIF()
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION ".*FAILED.*" PROCESSORS ${PROCS} )
ENDFUNCTION()


# Add a executable as a parallel test
FUNCTION( ADD_LBPM_TEST_PARALLEL EXEFILE PROCS ${ARGN} )
    ADD_LBPM_PROVISIONAL_TEST ( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    IF ( USE_MPI )
        CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs" ${ARGN} )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION ".*FAILED.*" PROCESSORS ${PROCS} )
    ENDIF()
ENDFUNCTION()


# Add a parallel test on 1, 2, and 4 processors
MACRO( ADD_LBPM_TEST_1_2_4 EXENAME ${ARGN} )
    ADD_LBPM_TEST ( ${EXENAME} ${ARGN} )
    ADD_LBPM_TEST_PARALLEL ( ${EXENAME} 2 ${ARGN} )
    ADD_LBPM_TEST_PARALLEL ( ${EXENAME} 4 ${ARGN} )
ENDMACRO()


# Add a parallel test that may use both MPI and threads
# This allows us to correctly compute the number of processors used by the test
MACRO( ADD_LBPM_TEST_THREAD_MPI EXEFILE PROCS THREADS ${ARGN} )
    ADD_LBPM_PROVISIONAL_TEST( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_${THREADS}threads" ${ARGN} )
    MATH( EXPR TOT_PROCS "${PROCS} * ${THREADS}" )
    IF ( ${TOT_PROCS} GREATER ${TEST_MAX_PROCS} )
        MESSAGE("Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS}")
    ELSEIF ( ( ${PROCS} STREQUAL "1" ) AND NOT USE_EXT_MPI_FOR_SERIAL_TESTS )
        ADD_TEST ( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
    ELSEIF ( USE_MPI )
        ADD_TEST ( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
    ENDIF()
    SET_TESTS_PROPERTIES ( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION ".*FAILED.*" PROCESSORS ${TOT_PROCS} )
ENDMACRO()


# Macro to check if a flag is enabled
MACRO ( CHECK_ENABLE_FLAG FLAG DEFAULT )
    IF( NOT DEFINED ${FLAG} )
        SET( ${FLAG} ${DEFAULT} )
    ELSEIF( ${FLAG}  STREQUAL "" )
        SET( ${FLAG} ${DEFAULT} )
    ELSEIF( ( ${${FLAG}} STREQUAL "false" ) OR ( ${${FLAG}} STREQUAL "0" ) OR ( ${${FLAG}} STREQUAL "OFF" ) )
        SET( ${FLAG} 0 )
    ELSEIF( ( ${${FLAG}} STREQUAL "true" ) OR ( ${${FLAG}} STREQUAL "1" ) OR ( ${${FLAG}} STREQUAL "ON" ) )
        SET( ${FLAG} 1 )
    ELSE()
        MESSAGE( "Bad value for ${FLAG} (${${FLAG}}); use true or false" )
    ENDIF ()
ENDMACRO()


# Macro to check if a compiler flag is valid
MACRO (CHECK_C_COMPILER_FLAG _FLAG _RESULT)
    SET(SAFE_CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS}")
    SET(CMAKE_REQUIRED_DEFINITIONS "${_FLAG}")
    CHECK_C_SOURCE_COMPILES("int main() { return 0;}" ${_RESULT}
        # Some compilers do not fail with a bad flag
        FAIL_REGEX "error: bad value (.*) for .* switch"       # GNU
        FAIL_REGEX "argument unused during compilation"        # clang
        FAIL_REGEX "is valid for .* but not for C"             # GNU
        FAIL_REGEX "unrecognized .*option"                     # GNU
        FAIL_REGEX "ignoring unknown option"                   # MSVC
        FAIL_REGEX "[Uu]nknown option"                         # HP
        FAIL_REGEX "[Ww]arning: [Oo]ption"                     # SunPro
        FAIL_REGEX "command option .* is not recognized"       # XL
        FAIL_REGEX "WARNING: unknown flag:"                    # Open64
        FAIL_REGEX " #10159: "                                 # ICC
    )
    SET(CMAKE_REQUIRED_DEFINITIONS "${SAFE_CMAKE_REQUIRED_DEFINITIONS}")
ENDMACRO(CHECK_C_COMPILER_FLAG)


# Macro to add a latex file to the build
MACRO (ADD_LATEX_DOCS FILE)
    GET_FILENAME_COMPONENT(LATEX_TARGET ${FILE} NAME_WE)
    ADD_CUSTOM_TARGET( 
        ${LATEX_TARGET}_pdf 
        ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${FILE} ${CMAKE_CURRENT_BINARY_DIR}/.
        COMMAND pdflatex -interaction=batchmode -draftmode ${FILE}
        #COMMAND bibtex -terse ${LATEX_TARGET}
        COMMAND pdflatex -interaction=batchmode ${FILE}
        SOURCES ${FILE}
    )
    ADD_CUSTOM_COMMAND( 
        TARGET ${LATEX_TARGET}_pdf 
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/${LATEX_TARGET}.pdf ${AMR_MHD_INSTALL_DIR}/doc/.
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    ADD_DEPENDENCIES( latex_docs ${LATEX_TARGET}_pdf )
ENDMACRO()


# add custom target distclean
# cleans and removes cmake generated files etc.
MACRO( ADD_DISTCLEAN )
    SET(DISTCLEANED
        cmake.depends
        cmake.check_depends
        CMakeCache.txt
        CMakeFiles
        CMakeTmp
        cmake.check_cache
        *.cmake
        compile.log
        Doxyfile
        Makefile
        core core.*
        DartConfiguration.tcl
        install_manifest.txt
        Testing
        include
        doc
        lib
        tests
        bin
        liblbpm-wia.a
        liblbpm-wia.so
        cpu
        gpu
        example
        common
    )
    ADD_CUSTOM_TARGET (distclean @echo cleaning for source distribution)
    IF (UNIX)
        ADD_CUSTOM_COMMAND(
            DEPENDS clean
            COMMENT "distribution clean"
            COMMAND rm
            ARGS    -Rf ${DISTCLEANED}
            TARGET  distclean
        )
    ELSE()
        SET( DISTCLEANED
            ${DISTCLEANED}
            *.vcxproj*
            ipch
            x64
        )
        FILE(WRITE  ${CMAKE_CURRENT_BINARY_DIR}/distclean.bat "del /s /q /f " )
        FOREACH (fileToDelete ${DISTCLEANED})
            FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/distclean.bat "${fileToDelete} " )
        ENDFOREACH ()
        FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/distclean.bat "\n" )
        FOREACH (fileToDelete ${DISTCLEANED})
            FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/distclean.bat "for /d %%x in (${fileToDelete}) do rd /s /q \"%%x\"\n" )
        ENDFOREACH ()
        ADD_CUSTOM_COMMAND(
            DEPENDS clean
            COMMENT "distribution clean"
            COMMAND distclean.bat & del /s/q/f distclean.bat
            TARGET  distclean
        )
    ENDIF()
ENDMACRO()
