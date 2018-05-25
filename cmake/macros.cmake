INCLUDE(CheckCCompilerFlag)
INCLUDE(CheckCSourceCompiles)
INCLUDE(CheckCXXCompilerFlag)
INCLUDE(CheckCXXSourceCompiles)
IF ( NOT TEST_FAIL_REGULAR_EXPRESSION )
    # Note: we cannot check for "handles are still allocated" due to PETSc.  See static variable
    #   Petsc_Reduction_keyval on line 234 of comb.c
    #SET( TEST_FAIL_REGULAR_EXPRESSION "(FAILED)|(leaked context IDs detected)|(handles are still allocated)" )
    SET( TEST_FAIL_REGULAR_EXPRESSION "(FAILED)" )
ENDIF()


# Check that the PROJ and ${PROJ}_INSTALL_DIR variables are set 
# These variables are used to generate the ADD_PROJ_TEST macros
IF ( NOT PROJ )
    MESSAGE(FATAL_ERROR "PROJ must be set before including macros.cmake")
ENDIF()
IF ( NOT ${PROJ}_INSTALL_DIR )
    MESSAGE(FATAL_ERROR "${PROJ}_INSTALL_DIR must be set before including macros.cmake")
ENDIF()
IF ( NOT CMAKE_BUILD_TYPE )
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE must be set before including macros.cmake")
ENDIF()


# Enable json
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )


# Check for link time optimization (LTO)
IF ( ${CMAKE_BUILD_TYPE} STREQUAL "Release" AND ${CMAKE_VERSION} VERSION_GREATER 3.9.4
    AND NOT DISABLE_LTO AND NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION )
    CMAKE_MINIMUM_REQUIRED(VERSION 3.9)
    INCLUDE( CheckIPOSupported )
    CHECK_IPO_SUPPORTED(RESULT supported OUTPUT error)
    IF( supported )
        MESSAGE(STATUS "IPO / LTO enabled")
        SET( LTO TRUE )
        SET( CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE )
    ELSE()
        SET( LTO FALSE )
        MESSAGE(STATUS "IPO / LTO not supported")
    ENDIF()
ELSEIF( NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION )
    SET( LTO FALSE )
    MESSAGE(STATUS "IPO / LTO disabled")
ENDIF()


# Add some default targets if they do not exist
IF ( NOT TARGET copy-${PROJ}-Data )
    ADD_CUSTOM_TARGET( copy-${PROJ}-Data ALL )
ENDIF()
IF ( NOT TARGET copy-${PROJ}-include )
    ADD_CUSTOM_TARGET ( copy-${PROJ}-include ALL )
ENDIF()


# Dummy use to prevent unused cmake variable warning
MACRO( NULL_USE VAR )
    IF ( "${${VAR}}" STREQUAL "dummy" )
        MESSAGE( FATAL_ERROR "NULL_USE fail" )
    ENDIF()
ENDMACRO()
NULL_USE( CMAKE_C_FLAGS )


# Macro to set a global variable
MACRO(GLOBAL_SET VARNAME)
  SET(${VARNAME} ${ARGN} CACHE INTERNAL "")
ENDMACRO()


# Macro to print all variables
MACRO( PRINT_ALL_VARIABLES )
    GET_CMAKE_PROPERTY(_variableNames VARIABLES)
    FOREACH(_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    ENDFOREACH()
ENDMACRO()


# CMake assert
MACRO(ASSERT test comment)
    IF (NOT ${test})
        MESSSAGE(FATAL_ERROR "Assertion failed: ${comment}")
    ENDIF(NOT ${test})
ENDMACRO(ASSERT)


# Set the maximum number of processors
IF ( NOT TEST_MAX_PROCS )
    INCLUDE(ProcessorCount)
    ProcessorCount( TEST_MAX_PROCS )
    IF ( ${TEST_MAX_PROCS} EQUAL 0 )
        SET( TEST_MAX_PROCS 16 )
    ENDIF()
ENDIF()


# Macro to convert a m4 file
# This command converts a file of the format "global_path/file.m4"
# and convertes it to file.F.  It also requires the path.  
MACRO( CONVERT_M4_FORTRAN IN LOCAL_PATH OUT_PATH )
    STRING(REGEX REPLACE ${LOCAL_PATH} "" OUT ${IN} )
    STRING(REGEX REPLACE "/" "" OUT ${OUT} )
    STRING(REGEX REPLACE "(.fm4)|(.m4)" ".F" OUT "${CMAKE_CURRENT_BINARY_DIR}/${OUT_PATH}/${OUT}" )
    IF ( NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${OUT_PATH}" )
        FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${OUT_PATH}" )    
    ENDIF()
    CONFIGURE_FILE ( ${IN} ${IN} COPYONLY )
    IF ("${CMAKE_GENERATOR}" STREQUAL "Xcode")
        STRING(REGEX REPLACE ".F" ".o" OUT2 "${OUT}" )
        STRING(REGEX REPLACE ";" " " COMPILE_CMD "${CMAKE_Fortran_COMPILER} -c ${OUT} ${CMAKE_Fortran_FLAGS} -o ${OUT2}")
        STRING(REGEX REPLACE "\\\\" "" COMPILE_CMD "${COMPILE_CMD}")
        MESSAGE("COMPILE_CMD = ${COMPILE_CMD}")
        SET( COMPILE_CMD ${COMPILE_CMD} )
        add_custom_command(
            OUTPUT ${OUT2}
            COMMAND m4 -I${LOCAL_PATH} -I${SAMRAI_FORTDIR} ${M4DIRS} ${IN} > ${OUT}
            COMMAND ${COMPILE_CMD}
            DEPENDS ${IN}
            )
        set_source_files_properties(${OUT2} PROPERTIES GENERATED true)
        SET( SOURCES ${SOURCES} "${OUT2}" )
     ELSE()
        add_custom_command(
            OUTPUT ${OUT}
            COMMAND m4 -I${LOCAL_PATH} -I${SAMRAI_FORTDIR} ${M4DIRS} ${M4_OPTIONS} ${IN} > ${OUT}
            DEPENDS ${IN}
            )
         set_source_files_properties(${OUT} PROPERTIES GENERATED true)
         SET( SOURCES ${SOURCES} "${OUT}" )
     ENDIF()
ENDMACRO()


# Add a package to the project's library
MACRO( ADD_${PROJ}_LIBRARY PACKAGE )
    SET( CURRENT_LIBRARY ${PACKAGE} )
    ADD_SUBDIRECTORY( ${PACKAGE} )
    SET( CURRENT_LIBRARY )
ENDMACRO()


# Add a project executable
MACRO( ADD_${PROJ}_EXECUTABLE EXEFILE )
    ADD_PROJ_PROVISIONAL_TEST( ${EXEFILE} )
    INSTALL( TARGETS ${EXEFILE} DESTINATION ${${PROJ}_INSTALL_DIR}/bin )
ENDMACRO()


# Initialize a package
MACRO (BEGIN_PACKAGE_CONFIG PACKAGE)
    SET( HEADERS "" )
    SET( CXXSOURCES "" )
    SET( CSOURCES "" )
    SET( FSOURCES "" )
    SET( M4FSOURCES "" )
    SET( CUDASOURCES "" )
    SET( SOURCES "" )
    SET( CURPACKAGE ${PACKAGE} )
ENDMACRO ()


# Find the source files
MACRO (FIND_FILES)
    # Find the C/C++ headers
    SET( T_HEADERS "" )
    FILE( GLOB T_HEADERS "*.h" "*.H" "*.hh" "*.hpp" "*.I" )
    # Find the CUDA sources
    SET( T_CUDASOURCES "" )
    FILE( GLOB T_CUDASOURCES "*.cu" )
    # Find the C sources
    SET( T_CSOURCES "" )
    FILE( GLOB T_CSOURCES "*.c" )
    # Find the C++ sources
    SET( T_CXXSOURCES "" )
    FILE( GLOB T_CXXSOURCES "*.cc" "*.cpp" "*.cxx" "*.C" )
    # Find the Fortran sources
    SET( T_FSOURCES "" )
    FILE( GLOB T_FSOURCES "*.f" "*.f90" "*.F" "*.F90" )
    # Find the m4 fortran source (and convert)
    SET( T_M4FSOURCES "" )
    FILE( GLOB T_M4FSOURCES "*.m4" "*.fm4" )
    FOREACH( m4file ${T_M4FSOURCES} )
        CONVERT_M4_FORTRAN( ${m4file} ${CMAKE_CURRENT_SOURCE_DIR} "" )
    ENDFOREACH()
    # Add all found files to the current lists
    SET( HEADERS ${HEADERS} ${T_HEADERS} )
    SET( CXXSOURCES ${CXXSOURCES} ${T_CXXSOURCES} )
    SET( CUDASOURCES ${CUDASOURCES} ${T_CUDASOURCES} )
    SET( CSOURCES ${CSOURCES} ${T_CSOURCES} )
    SET( FSOURCES ${FSOURCES} ${T_FSOURCES} )
    SET( M4FSOURCES ${M4FSOURCES} ${T_M4FSOURCES} )
    SET( SOURCES ${SOURCES} ${T_CXXSOURCES} ${T_CSOURCES} ${T_FSOURCES} ${T_M4FSOURCES} )
ENDMACRO()


# Find the source files
MACRO( FIND_FILES_PATH IN_PATH )
    # Find the C/C++ headers
    SET( T_HEADERS "" )
    FILE( GLOB T_HEADERS "${IN_PATH}/*.h" "${IN_PATH}/*.H" "${IN_PATH}/*.hh" "${IN_PATH}/*.hpp" "${IN_PATH}/*.I" )
    # Find the CUDA sources
    SET( T_CUDASOURCES "" )
    FILE( GLOB T_CUDASOURCES "${IN_PATH}/*.cu" )
    # Find the C sources
    SET( T_CSOURCES "" )
    FILE( GLOB T_CSOURCES "${IN_PATH}/*.c" )
    # Find the C++ sources
    SET( T_CXXSOURCES "" )
    FILE( GLOB T_CXXSOURCES "${IN_PATH}/*.cc" "${IN_PATH}/*.cpp" "${IN_PATH}/*.cxx" "${IN_PATH}/*.C" )
    # Find the Fortran sources
    SET( T_FSOURCES "" )
    FILE( GLOB T_FSOURCES "${IN_PATH}/*.f" "${IN_PATH}/*.f90" "${IN_PATH}/*.F" "${IN_PATH}/*.F90" )
    # Find the m4 fortran source (and convert)
    SET( T_M4FSOURCES "" )
    FILE( GLOB T_M4FSOURCES "${IN_PATH}/*.m4" "${IN_PATH}/*.fm4" )
    FOREACH( m4file ${T_M4FSOURCES} )
        CONVERT_M4_FORTRAN( ${m4file} ${CMAKE_CURRENT_SOURCE_DIR}/${IN_PATH} ${IN_PATH} )
    ENDFOREACH ()
    # Add all found files to the current lists
    SET( HEADERS ${HEADERS} ${T_HEADERS} )
    SET( CXXSOURCES ${CXXSOURCES} ${T_CXXSOURCES} )
    SET( CUDASOURCES ${CUDASOURCES} ${T_CUDASOURCES} )
    SET( CSOURCES ${CSOURCES} ${T_CSOURCES} )
    SET( FSOURCES ${FSOURCES} ${T_FSOURCES} )
    SET( SOURCES ${SOURCES} ${T_CXXSOURCES} ${T_CSOURCES} ${T_FSOURCES} ${CUDASOURCES} )
ENDMACRO()


# Add a subdirectory
MACRO( ADD_PACKAGE_SUBDIRECTORY SUBDIR )
    FIND_FILES_PATH( ${SUBDIR} )
ENDMACRO()


# Add an external directory
# Note: this will add the files to compile but will not copy the headers
MACRO( ADD_EXTERNAL_DIRECTORY SUBDIR )
    FIND_FILES_PATH( ${SUBDIR} )
ENDMACRO()


# Install a package
MACRO( INSTALL_${PROJ}_TARGET PACKAGE )
    # Find all files in the current directory
    FIND_FILES()
    # Create the copy target
    STRING(REGEX REPLACE "${${PROJ}_SOURCE_DIR}/" "" COPY_TARGET "copy-${PROJ}-${CMAKE_CURRENT_SOURCE_DIR}-include" )
    STRING(REGEX REPLACE "/" "-" COPY_TARGET ${COPY_TARGET} )
    IF( NOT TARGET ${COPY_TARGET} )
        ADD_CUSTOM_TARGET( ${COPY_TARGET} ALL )
        ADD_DEPENDENCIES( copy-${PROJ}-include ${COPY_TARGET} )
    ENDIF()
    # Copy the header files to the include path
    IF ( HEADERS )
        FILE( GLOB HFILES RELATIVE "${${PROJ}_SOURCE_DIR}" ${HEADERS} )
        FOREACH( HFILE ${HFILES} )
            SET( SRC_FILE "${${PROJ}_SOURCE_DIR}/${HFILE}" )
            SET( DST_FILE "${${PROJ}_INSTALL_DIR}/include/${${PROJ}_INC}/${HFILE}" )
            # Only copy the headers if the exisit in the project source directory
            IF ( EXISTS "${SRC_FILE}" )
                ADD_CUSTOM_COMMAND(TARGET ${COPY_TARGET} 
                    PRE_BUILD 
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SRC_FILE}" "${DST_FILE}"
                    DEPENDS "${SRC_FILE}"
                )
            ENDIF()
        ENDFOREACH()
    ENDIF()
    # Add the library and install the package
    IF ( NOT ONLY_BUILD_DOCS AND SOURCES )
        # Set RPATH variables
        IF ( NOT CMAKE_RPATH_VARIABLES_SET AND NOT USE_STATIC )
            SET(CMAKE_RPATH_VARIABLES_SET ON)
            SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
            SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
            SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} "${CMAKE_INSTALL_PREFIX}/lib")
            SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
            SET(MACOSX_RPATH 0)
            LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
        ENDIF()
        # Add the library to the project libs list
        SET( ${PROJ}_LIBS ${${PROJ}_LIBS} ${PACKAGE} CACHE INTERNAL "")
        LIST( REMOVE_DUPLICATES ${PROJ}_LIBS )
        SET( ${PROJ}_LIBS ${${PROJ}_LIBS} CACHE INTERNAL "")
        # Create the library
        IF ( ${PROJ}_LIB )
            # We are using a single project library
            ADD_LIBRARY( ${PACKAGE} OBJECT ${SOURCES} )
        ELSE()
            # We are creating individual libraries
            ADD_LIBRARY( ${PACKAGE} ${LIB_TYPE} ${SOURCES} ${CUBINS} )
            TARGET_LINK_EXTERNAL_LIBRARIES( ${PACKAGE} )
        ENDIF()
        # Add coverage flags to target
        IF ( NOT DISABLE_TARGET_COVERAGE )
            TARGET_COMPILE_DEFINITIONS( ${PACKAGE} PUBLIC ${COVERAGE_FLAGS} )
        ENDIF()
        # Add target dependencies
        IF ( TARGET write_repo_version )
            ADD_DEPENDENCIES( ${PACKAGE} write_repo_version )
        ENDIF()
        ADD_DEPENDENCIES ( ${PACKAGE} copy-${PROJ}-include )
        IF ( NOT ${PROJ}_LIB )
            INSTALL( TARGETS ${PACKAGE} DESTINATION ${${PROJ}_INSTALL_DIR}/lib )
        ENDIF()
    ELSE()
        ADD_CUSTOM_TARGET( ${PACKAGE} ALL )
    ENDIF()
    # Clear the sources
    SET( HEADERS "" )
    SET( CSOURCES "" )
    SET( CXXSOURCES "" )
ENDMACRO()


# Install the project library
MACRO( INSTALL_PROJ_LIB )
    SET( tmp_link_list )
    FOREACH ( tmp ${${PROJ}_LIBS} )
        SET( tmp_link_list ${tmp_link_list} $<TARGET_OBJECTS:${tmp}> )
    ENDFOREACH()
    ADD_LIBRARY( ${${PROJ}_LIB} ${LIB_TYPE} ${tmp_link_list} )
    TARGET_LINK_EXTERNAL_LIBRARIES( ${${PROJ}_LIB} LINK_PUBLIC )
    INSTALL( TARGETS ${${PROJ}_LIB} DESTINATION ${${PROJ}_INSTALL_DIR}/lib )
ENDMACRO()


# Macro to verify that a variable has been set
MACRO( VERIFY_VARIABLE VARIABLE_NAME )
    IF ( NOT ${VARIABLE_NAME} )
        MESSAGE( FATAL_ERROR "Please set: " ${VARIABLE_NAME} )
    ENDIF()
ENDMACRO()


# Macro to verify that a path has been set
MACRO( VERIFY_PATH PATH_NAME )
    IF ("${PATH_NAME}" STREQUAL "")
        MESSAGE ( FATAL_ERROR "Path is not set: ${PATH_NAME}" )
    ENDIF()
    IF ( NOT EXISTS "${PATH_NAME}" )
        MESSAGE( FATAL_ERROR "Path does not exist: ${PATH_NAME}" )
    ENDIF()
ENDMACRO()


# Macro to tell cmake to use static libraries
MACRO( SET_STATIC_FLAGS )
    # Remove extra library links
    SET(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
    SET(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
    # Add the static flag if necessary
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static") # Add static flag
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -static ") 
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static ")
ENDMACRO()


# Macro to identify the compiler
MACRO( IDENTIFY_COMPILER )
    # SET the C/C++ compiler
    IF ( CMAKE_C_COMPILER_WORKS OR CMAKE_CXX_COMPILER_WORKS )
        IF( USING_GCC OR CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR
            (${CMAKE_C_COMPILER_ID} MATCHES "GNU") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU") )
            SET( USING_GCC TRUE )
            ADD_DEFINITIONS( -DUSING_GCC )
            MESSAGE("Using gcc")
        ELSEIF( USING_MSVC OR MSVC OR MSVC_IDE OR MSVC60 OR MSVC70 OR MSVC71 OR MSVC80 OR CMAKE_COMPILER_2005 OR MSVC90 OR MSVC10 )
            IF( NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows" )
                MESSAGE( FATAL_ERROR "Using microsoft compilers on non-windows system?" )
            ENDIF()
            SET( USING_MSVC TRUE )
            ADD_DEFINITIONS( -DUSING_MSVC )
            MESSAGE("Using Microsoft")
        ELSEIF( USING_ICC OR (${CMAKE_C_COMPILER_ID} MATCHES "Intel") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Intel") ) 
            SET(USING_ICC TRUE)
            ADD_DEFINITIONS( -DUSING_ICC )
            MESSAGE("Using icc")
        ELSEIF( USING_PGCC OR (${CMAKE_C_COMPILER_ID} MATCHES "PGI") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "PGI") )
            SET(USING_PGCC TRUE)
            ADD_DEFINITIONS( -DUSING_PGCC )
            MESSAGE("Using pgCC")
        ELSEIF( USING_CRAY OR (${CMAKE_C_COMPILER_ID} MATCHES "CRAY") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "CRAY") OR
                              (${CMAKE_C_COMPILER_ID} MATCHES "Cray") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Cray") )
            SET(USING_CRAY TRUE)
            ADD_DEFINITIONS( -DUSING_CRAY )
            MESSAGE("Using Cray")
        ELSEIF( USING_CLANG OR (${CMAKE_C_COMPILER_ID} MATCHES "CLANG") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "CLANG") OR
                               (${CMAKE_C_COMPILER_ID} MATCHES "Clang") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang") )
            SET(USING_CLANG TRUE)
            ADD_DEFINITIONS( -DUSING_CLANG )
            MESSAGE("Using Clang")
        ELSEIF( USING_XL OR (${CMAKE_C_COMPILER_ID} MATCHES "XL") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "XL") )
            SET(USING_XL TRUE)
            ADD_DEFINITIONS( -DUSING_XL )
            MESSAGE("Using XL")
        ELSE()
            MESSAGE( "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
            MESSAGE( "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
            MESSAGE( "CMAKE_C_COMPILER_ID=${CMAKE_C_COMPILER_ID}")
            MESSAGE( "CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
            MESSAGE(FATAL_ERROR "Unknown C/C++ compiler")
        ENDIF()
    ENDIF()
    # SET the Fortran compiler
    IF ( CMAKE_Fortran_COMPILER_WORKS )
        IF( CMAKE_COMPILER_IS_GNUG77 OR (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU") )
            SET( USING_GFORTRAN TRUE )
            MESSAGE("Using gfortran")
            IF ( NOT USING_GCC )
                LIST( REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES gcc )
            ENDIF()
        ELSEIF ( (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel") ) 
            SET(USING_IFORT TRUE)
            MESSAGE("Using ifort")
        ELSEIF ( ${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
            SET(USING_PGF90 TRUE)
            MESSAGE("Using pgf90")
        ELSEIF ( (${CMAKE_Fortran_COMPILER_ID} MATCHES "CLANG") OR (${CMAKE_Fortran_COMPILER_ID} MATCHES "Clang") OR
                 (${CMAKE_Fortran_COMPILER_ID} MATCHES "FLANG") OR (${CMAKE_Fortran_COMPILER_ID} MATCHES "Flang") )
            SET(USING_FLANG TRUE)
            MESSAGE("Using flang")
        ELSE()
            MESSAGE( "CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
            MESSAGE( "CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}")
            MESSAGE(FATAL_ERROR "Unknown Fortran compiler (${CMAKE_Fortran_COMPILER_ID})")
        ENDIF()
    ENDIF()
ENDMACRO()


# Macro to set the proper warnings
MACRO( SET_WARNINGS )
  IF ( USING_GCC )
    # Add gcc specific compiler options
    # Note: adding -Wlogical-op causes a wierd linking error on Titan using the nvcc wrapper:
    #    /usr/bin/ld: cannot find gical-op: No such file or directory
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -Wall -Wextra") 
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Woverloaded-virtual -Wsign-compare")
  ELSEIF ( USING_MSVC )
    # Add Microsoft specifc compiler options
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} /D _SCL_SECURE_NO_WARNINGS /D _CRT_SECURE_NO_WARNINGS /D _ITERATOR_DEBUG_LEVEL=0 /wd4267" )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _SCL_SECURE_NO_WARNINGS /D _CRT_SECURE_NO_WARNINGS /D _ITERATOR_DEBUG_LEVEL=0 /wd4267" )
  ELSEIF ( USING_ICC )
    # Add Intel specifc compiler options
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -Wall" )
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall" )
  ELSEIF ( USING_CRAY )
    # Add default compiler options
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS}")
  ELSEIF ( USING_PGCC )
    # Add default compiler options
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -lpthread")
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -lpthread -Minform=inform -Mlist --display_error_number")
    # Suppress unreachable code warning, it causes non-useful warnings with some tests/templates
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} --diag_suppress 111,128,185")
  ELSEIF ( USING_CLANG )
    # Add default compiler options
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -Wall")
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall -Wno-missing-braces -Wmissing-field-initializers -ftemplate-depth=1024")
  ELSEIF ( USING_XL )
    # Add default compiler options
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -Wall")
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall -ftemplate-depth=512")
  ELSE ( )
    MESSAGE("Compiler specific features are not set for this compiler")
  ENDIF()
ENDMACRO()


# Macro to add user compile flags
MACRO( ADD_USER_FLAGS )
    STRING( STRIP "${CMAKE_C_FLAGS} ${CFLAGS} ${CFLAGS_EXTRA}" CMAKE_C_FLAGS )
    STRING( STRIP "${CMAKE_CXX_FLAGS} ${CXXFLAGS} ${CXXFLAGS_EXTRA}" CMAKE_CXX_FLAGS )
    STRING( STRIP "${CMAKE_Fortran_FLAGS} ${FFLAGS} ${FFLAGS_EXTRA}" CMAKE_Fortran_FLAGS )
    STRING( STRIP "${LDFLAGS} ${LDFLAGS_EXTRA}" LDFLAGS )
    STRING( STRIP "${LDLIBS} ${LDLIBS_EXTRA}" LDLIBS )
ENDMACRO()


# Macro to set the compile/link flags
MACRO( SET_COMPILER_FLAGS )
    # Initilaize the compiler
    IDENTIFY_COMPILER()
    # Set the default flags for each build type
    IF ( USING_MSVC )
        SET(CMAKE_C_FLAGS_DEBUG   "-D_DEBUG /DEBUG /Od /EHsc /MDd /Z7" )
        SET(CMAKE_CXX_FLAGS_DEBUG "-D_DEBUG /DEBUG /Od /EHsc /MDd /Z7" )
        SET(CMAKE_C_FLAGS_RELEASE   "/O2 /EHsc /MD" )
        SET(CMAKE_CXX_FLAGS_RELEASE "/O2 /EHsc /MD" )
    ELSE()
    ENDIF()
    # Set the behavior of GLIBCXX flags
    CHECK_ENABLE_FLAG( ENABLE_GXX_DEBUG 0 )
    IF ( ENABLE_GXX_DEBUG ) 
        # Enable GLIBCXX_DEBUG flags
        SET( CMAKE_C_FLAGS_DEBUG   " ${CMAKE_C_FLAGS_DEBUG}   -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" )
        SET( CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" )
        SET( DISABLE_GXX_DEBUG OFF )
    ELSEIF ( DISABLE_GXX_DEBUG ) 
        # Disable GLIBCXX_DEBUG flags
        SET( DISABLE_GXX_DEBUG OFF )
    ELSE()
        # Default
        SET( DISABLE_GXX_DEBUG ON )
    ENDIF()
    # Set debug definitions
    IF ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" AND NOT ("${CMAKE_CXX_FLAGS_DEBUG}" MATCHES "-D_DEBUG") )
        SET( CMAKE_C_FLAGS_DEBUG   " ${CMAKE_C_FLAGS_DEBUG}   -DDEBUG -D_DEBUG" )
        SET( CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -D_DEBUG" )        
    ENDIF()
    # Save the debug/release specific flags to the cache
    SET( CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS_DEBUG}"     CACHE STRING "Debug flags"   FORCE)
    SET( CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}"   CACHE STRING "Release flags" FORCE)
    SET( CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}"   CACHE STRING "Debug flags"   FORCE)
    SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "Release flags" FORCE)
    # Add the user flags
    ADD_USER_FLAGS()
    # Set the warnings to use
    SET_WARNINGS()
    # Test the compile flags
    IF ( CMAKE_C_COMPILER_WORKS )
        CHECK_C_COMPILER_FLAG( "${CMAKE_C_FLAGS}" CHECK_C_FLAGS )
    ENDIF()
    IF ( CMAKE_CXX_COMPILER_WORKS )
        CHECK_CXX_COMPILER_FLAG( "${CMAKE_CXX_FLAGS}" CHECK_CXX_FLAGS )
    ENDIF()
    IF ( ( NOT CHECK_C_FLAGS ) OR ( NOT CHECK_CXX_FLAGS ) )
        IF ( USING_CRAY )
            MESSAGE(WARNING "Invalid C/CXX flags detected:\n"
                "C flags: ${CMAKE_C_FLAGS}\n" "CXX flags: ${CMAKE_CXX_FLAGS}\n" )
        ENDIF()
    ENDIF()
ENDMACRO()


# Macro to copy data file at build time
MACRO( COPY_DATA_FILE SRC_FILE DST_FILE )
    STRING(REGEX REPLACE "${${PROJ}_SOURCE_DIR}/" "" COPY_TARGET "copy-${PROJ}-${CMAKE_CURRENT_SOURCE_DIR}" )
    STRING(REGEX REPLACE "-${${PROJ}_SOURCE_DIR}" "" COPY_TARGET "${COPY_TARGET}" )
    STRING(REGEX REPLACE "/" "-" COPY_TARGET ${COPY_TARGET} )
    IF ( NOT TARGET ${COPY_TARGET} )
        ADD_CUSTOM_TARGET( ${COPY_TARGET} ALL )
        ADD_DEPENDENCIES( copy-${PROJ}-Data ${COPY_TARGET} )
    ENDIF()
    ADD_CUSTOM_COMMAND( TARGET ${COPY_TARGET}
        PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SRC_FILE}" "${DST_FILE}"
        DEPENDS "${SRC_FILE}"
    )
ENDMACRO()


# Macro to copy a data or input file
FUNCTION( COPY_TEST_FILE FILENAME ${ARGN} )
    SET( SEARCH_PATHS "${CMAKE_CURRENT_SOURCE_DIR}" )
    SET( SEARCH_PATHS ${SEARCH_PATHS} "${CMAKE_CURRENT_SOURCE_DIR}/data" )
    SET( SEARCH_PATHS ${SEARCH_PATHS} "${CMAKE_CURRENT_SOURCE_DIR}/inputs" )
    FOREACH( tmp ${ARGN} )
        SET( SEARCH_PATHS ${SEARCH_PATHS} "${tmp}" )
        SET( SEARCH_PATHS ${SEARCH_PATHS} "${CMAKE_CURRENT_SOURCE_DIR}/${tmp}" )
    ENDFOREACH()
    SET( FILE_TO_COPY )
    FOREACH( tmp ${SEARCH_PATHS} )
        IF ( EXISTS "${tmp}/${FILENAME}" )
            SET( FILE_TO_COPY "${tmp}/${FILENAME}" )
        ENDIF()
    ENDFOREACH()
    IF ( FILE_TO_COPY )
        SET( DESTINATION_NAME "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}" )
        COPY_DATA_FILE( ${FILE_TO_COPY} ${DESTINATION_NAME} )
    ELSE()
        SET( MSG_STR "Cannot find file: ${FILENAME}, searched:\n" )
        FOREACH( tmp ${SEARCH_PATHS} )
            SET( MSG_STR "${MSG_STR}   ${tmp}\n" )
        ENDFOREACH()
        MESSAGE( WARNING "test ${MSG_STR}" )
    ENDIF()
ENDFUNCTION()


# Macro to copy a data file
FUNCTION( RENAME_TEST_FILE SRC_NAME DST_NAME ${ARGN} )
    SET( SEARCH_PATHS "${CMAKE_CURRENT_SOURCE_DIR}" )
    SET( SEARCH_PATHS ${SEARCH_PATHS} "${CMAKE_CURRENT_SOURCE_DIR}/data" )
    SET( SEARCH_PATHS ${SEARCH_PATHS} "${CMAKE_CURRENT_SOURCE_DIR}/inputs" )
    FOREACH( tmp ${ARGN} )
        SET( SEARCH_PATHS ${SEARCH_PATHS} "${tmp}" )
        SET( SEARCH_PATHS ${SEARCH_PATHS} "${CMAKE_CURRENT_SOURCE_DIR}/${tmp}" )
    ENDFOREACH()
    SET( FILE_TO_COPY )
    FOREACH( tmp ${SEARCH_PATHS} )
        IF ( EXISTS "${tmp}/${SRC}" )
            SET( FILE_TO_COPY "${tmp}/${SRC_NAME}" )
        ENDIF()
    ENDFOREACH()
    IF ( FILE_TO_COPY )
        SET( DESTINATION_NAME ${CMAKE_CURRENT_BINARY_DIR}/${DST_NAME} )
        COPY_DATA_FILE( ${FILE_TO_COPY} ${DESTINATION_NAME} )
    ELSE()
        SET( MSG_STR "Cannot find file: ${SRC_NAME}, searched:\n" )
        FOREACH( tmp ${SEARCH_PATHS} )
            SET( MSG_STR "${MSG_STR}   ${tmp}\n" )
        ENDFOREACH()
        MESSAGE( WARNING "test ${MSG_STR}" )
    ENDIF()
ENDFUNCTION()


# Macro to copy a data file
FUNCTION( COPY_EXAMPLE_DATA_FILE FILENAME )
    SET( FILE_TO_COPY  ${CMAKE_CURRENT_SOURCE_DIR}/data/${FILENAME} )
    SET( DESTINATION1 ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
    SET( DESTINATION2 ${EXAMPLE_INSTALL_DIR}/${FILENAME} )
    IF ( EXISTS ${FILE_TO_COPY} )
        COPY_DATA_FILE( ${FILE_TO_COPY} ${DESTINATION1} )
        COPY_DATA_FILE( ${FILE_TO_COPY} ${DESTINATION2} )
    ELSE()
        MESSAGE( WARNING "Cannot find file: " ${FILE_TO_COPY} )
    ENDIF()
ENDFUNCTION()


# Macro to copy a mesh file
MACRO( COPY_MESH_FILE MESHNAME )
    # Check the local data directory
    FILE( GLOB MESHPATH "${CMAKE_CURRENT_SOURCE_DIR}/data/${MESHNAME}" )
    # Search all data paths
    FOREACH ( path ${${PROJ}_DATA} )
        # Check the DATA directory
        IF ( NOT MESHPATH )
            FILE( GLOB MESHPATH "${path}/${MESHNAME}" )
        ENDIF()
        # Check the DATA/vvu directory
        IF ( NOT MESHPATH )
            FILE( GLOB MESHPATH "${path}/vvu/meshes/${MESHNAME}" )
        ENDIF()
        # Check the DATA/meshes directory
        IF ( NOT MESHPATH )
            FILE( GLOB_RECURSE MESHPATH "${path}/meshes/*/${MESHNAME}" )
        ENDIF()
        # Check the entire DATA directory
        IF ( NOT MESHPATH )
            FILE( GLOB_RECURSE MESHPATH "${path}/*/${MESHNAME}" )
        ENDIF()
    ENDFOREACH()
    # We have either found the mesh or failed
    IF ( NOT MESHPATH )
        MESSAGE ( WARNING "Cannot find mesh: " ${MESHNAME} )
    ELSE ()
        SET( MESHPATH2 )
        FOREACH( tmp ${MESHPATH} )
            SET( MESHPATH2 "${tmp}" )
        ENDFOREACH()
        STRING(REGEX REPLACE "//${MESHNAME}" "" MESHPATH "${MESHPATH2}" )
        STRING(REGEX REPLACE "${MESHNAME}" "" MESHPATH "${MESHPATH}" )
        COPY_DATA_FILE( "${MESHPATH}/${MESHNAME}" "${CMAKE_CURRENT_BINARY_DIR}/${MESHNAME}" )
    ENDIF()
ENDMACRO()


# Link the libraries to the given target
MACRO( TARGET_LINK_EXTERNAL_LIBRARIES TARGET_NAME )
    FOREACH ( tmp ${TPL_LIBS} )
        TARGET_LINK_LIBRARIES( ${TARGET_NAME} ${ARGN} ${tmp} )
    ENDFOREACH()
    FOREACH ( tmp ${EXTERNAL_LIBS} )
        TARGET_LINK_LIBRARIES( ${TARGET_NAME} ${ARGN} ${tmp} )
    ENDFOREACH()
    FOREACH ( tmp ${LAPACK_LIBS} )
        TARGET_LINK_LIBRARIES( ${TARGET_NAME} ${ARGN} ${tmp} )
    ENDFOREACH()
    FOREACH ( tmp ${BLAS_LIBS} )
        TARGET_LINK_LIBRARIES( ${TARGET_NAME} ${ARGN} ${tmp} )
    ENDFOREACH()
    FOREACH ( tmp ${BLAS_LAPACK_LIBS} )
        TARGET_LINK_LIBRARIES( ${TARGET_NAME} ${ARGN} ${tmp} )
    ENDFOREACH()
    FOREACH ( MPI_LIBRARIES )
        TARGET_LINK_LIBRARIES( ${EXE} ${ARGN} ${tmp} )
    ENDFOREACH()
    FOREACH ( tmp ${CMAKE_C_IMPLICIT_LINK_LIBRARIES}
        ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} )
        TARGET_LINK_LIBRARIES( ${TARGET_NAME} ${ARGN} ${tmp} )
    ENDFOREACH()
ENDMACRO()


# Choose the debug or optimized library based on the build type
FUNCTION( KEEP_BUILD_LIBRARIES VAR )
    IF ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
        SET( build_type debug )
    ELSE()
        SET( build_type optimized )
    ENDIF()
    SET( build ${build_type} )
    SET( LIBS )
    FOREACH ( tmp ${${VAR}} )
        IF ( ( ${tmp} STREQUAL debug ) OR ( ${tmp} STREQUAL optimized ) )
            SET( build ${tmp} )
        ELSEIF ( ${build} STREQUAL ${build_type} )
            SET( LIBS ${LIBS} ${tmp} )
        ENDIF()
    ENDFOREACH()
    SET( ${VAR} ${LIBS} PARENT_SCOPE )
ENDFUNCTION()


# Macro to add the dependencies and libraries to an executable
MACRO( ADD_PROJ_EXE_DEP EXE )
    # Add the package dependencies
    IF( ${PROJ}_TEST_LIB_EXISTS )
        ADD_DEPENDENCIES ( ${EXE} ${PACKAGE_TEST_LIB} )
        TARGET_LINK_LIBRARIES ( ${EXE} ${PACKAGE_TEST_LIB} )
    ENDIF()
    # Add the executable to the dependencies of check and build-test
    ADD_DEPENDENCIES( check ${EXE} )
    ADD_DEPENDENCIES( build-test ${EXE} )
    # Add the file copy targets to the dependency list
    IF ( TARGET copy-${PROJ}-Data )
        ADD_DEPENDENCIES( ${EXE} copy-${PROJ}-Data )
    ENDIF()
    # Add the project libraries
    IF ( ${PROJ}_LIB )
        TARGET_LINK_LIBRARIES( ${EXE} ${${PROJ}_LIB} )
    ELSE()
        TARGET_LINK_LIBRARIES( ${EXE} ${${PROJ}_LIBS} ${${PROJ}_LIBS} )
    ENDIF()
    TARGET_LINK_LIBRARIES( ${EXE} ${${PROJECT_NAME}_LIBRARIES} )
    # Add coverage flags to target
    IF ( NOT DISABLE_TARGET_COVERAGE )
        TARGET_COMPILE_DEFINITIONS( ${EXE} PUBLIC ${COVERAGE_FLAGS} )
    ENDIF()
    # Link to external libraries
    TARGET_LINK_LIBRARIES( ${EXE} ${LINK_LIBRARIES} )
    TARGET_LINK_EXTERNAL_LIBRARIES( ${EXE} )
    TARGET_LINK_LIBRARIES( ${EXE} ${COVERAGE_LIBS} ${LDLIBS} ${LDLIBS_EXTRA} )
    TARGET_LINK_LIBRARIES( ${EXE} ${SYSTEM_LIBS} ${SYSTEM_LDFLAGS} )
    # Set extra link flags
    IF ( USE_MPI OR USE_EXT_MPI OR HAVE_MPI )
        SET_TARGET_PROPERTIES( ${EXE} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS} ${LDFLAGS} ${LDFLAGS_EXTRA}" )
    ELSE()
        SET_TARGET_PROPERTIES( ${EXE} PROPERTIES LINK_FLAGS "${LDFLAGS} ${LDFLAGS_EXTRA}" )
    ENDIF()
ENDMACRO()


# Check if we want to keep the test
FUNCTION( KEEP_TEST RESULT )
    SET( ${RESULT} 1 PARENT_SCOPE )
    IF ( NOT DEFINED ${PACKAGE_NAME}_ENABLE_TESTS )
        SET( ${PACKAGE_NAME}_ENABLE_TESTS 1 )
    ENDIF()
    IF ( NOT ${PACKAGE_NAME}_ENABLE_TESTS )
        SET( ${RESULT} 0 PARENT_SCOPE )
    ENDIF()
    IF ( ONLY_BUILD_DOCS )
        SET( ${RESULT} 0 PARENT_SCOPE )
    ENDIF()
ENDFUNCTION()


# Add a provisional test
FUNCTION( ADD_PROJ_PROVISIONAL_TEST EXEFILE )
    # Check if we actually want to add the test
    KEEP_TEST( RESULT )
    IF ( NOT RESULT )
        RETURN()
    ENDIF()
    # Check if the target does not exist (may not be added yet or we are re-configuring)
    IF ( NOT TARGET ${EXEFILE} )
        GLOBAL_SET( ${EXEFILE}-BINDIR )
    ENDIF()
    # Check if test has already been added
    IF ( NOT ${EXEFILE}-BINDIR )
        GLOBAL_SET( ${EXEFILE}-BINDIR "${CMAKE_CURRENT_BINARY_DIR}" )
        # The target has not been added
        SET( CXXFILE ${EXEFILE} )
        SET( TESTS_SO_FAR ${TESTS_SO_FAR} ${EXEFILE} )
        # Check if we want to add the test to all
        IF ( NOT EXCLUDE_TESTS_FROM_ALL )
            ADD_EXECUTABLE( ${EXEFILE} ${CXXFILE} )
        ELSE()
            ADD_EXECUTABLE( ${EXEFILE} EXCLUDE_FROM_ALL ${CXXFILE} )
        ENDIF()
        ADD_PROJ_EXE_DEP( ${EXEFILE} )
        IF ( CURRENT_LIBRARY )
            IF ( NOT TARGET ${CURRENT_LIBRARY}-test )
                ADD_CUSTOM_TARGET( ${CURRENT_LIBRARY}-test )
            ENDIF()
            ADD_DEPENDENCIES( ${CURRENT_LIBRARY}-test ${EXEFILE} )
        ENDIF()
    ELSEIF ( NOT ${EXEFILE-BINDIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR} )
        # We are trying to add 2 different tests with the same name
        MESSAGE( "Existing test: ${EXEFILE-BINDIR}/${EXEFILE}" )
        MESSAGE( "New test:      ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}" )
        MESSAGE( FATAL_ERROR "Trying to add 2 different tests with the same name" )
    ENDIF()
    SET( LAST_TEST ${EXEFILE} PARENT_SCOPE )
ENDFUNCTION()
FUNCTION( ADD_${PROJ}_PROVISIONAL_TEST EXEFILE )
    ADD_PROJ_PROVISIONAL_TEST( ${EXEFILE} )
ENDFUNCTION()


# Macro to create the test name
MACRO( CREATE_TEST_NAME TEST ${ARGN} )
    IF ( PACKAGE )
        SET( TESTNAME "${PACKAGE}::${TEST}" )
    ELSE()
        SET( TESTNAME "${TEST}" )
    ENDIF()
    FOREACH( tmp ${ARGN})
        SET( TESTNAME "${TESTNAME}--${tmp}")
    endforeach()
    # STRING(REGEX REPLACE "--" "-" TESTNAME ${TESTNAME} )
    SET( LAST_TESTNAME ${TESTNAME} PARENT_SCOPE )
ENDMACRO()


# Function to add the resource locks to an executable
FUNCTION( ADD_RESOURCE_LOCK TESTNAME EXEFILE ${ARGN} )
    IF ( NOT DISABLE_RESOURCE_LOCK )
        IF ( NOT ARGN )
            SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES RESOURCE_LOCK ${EXEFILE} )
        ELSE()
            FOREACH( tmp ${ARGN} )
                SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES RESOURCE_LOCK ${tmp} )
            ENDFOREACH()
        ENDIF()
    ENDIF()
ENDFUNCTION()


# Add a executable as a test
FUNCTION( ADD_${PROJ}_TEST EXEFILE ${ARGN} )
    # Check if we actually want to add the test
    KEEP_TEST( RESULT )
    IF ( NOT RESULT )
        RETURN()
    ENDIF()
    # Add the provisional test
    ADD_PROJ_PROVISIONAL_TEST( ${EXEFILE} )
    CREATE_TEST_NAME( ${EXEFILE} ${ARGN} )
    IF ( USE_MPI_FOR_SERIAL_TESTS )
        ADD_TEST( NAME ${TESTNAME} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${LAST_TEST}> ${ARGN} )
        SET_PROPERTY( TEST ${TESTNAME} APPEND PROPERTY ENVIRONMENT OMPI_MCA_hwloc_base_binding_policy=none )
    ELSE()
        ADD_TEST( NAME ${TESTNAME} COMMAND $<TARGET_FILE:${LAST_TEST}> ${ARGN} )
    ENDIF()
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS 1 )
    ADD_RESOURCE_LOCK( ${TESTNAME} ${EXEFILE} ${ARGN} )
ENDFUNCTION()


# Add a executable as a weekly test
FUNCTION( ADD_${PROJ}_WEEKLY_TEST EXEFILE PROCS ${ARGN} )
    # Check if we actually want to add the test
    KEEP_TEST( RESULT )
    IF ( NOT RESULT )
        RETURN()
    ENDIF()
    # Add the provisional test
    ADD_PROJ_PROVISIONAL_TEST( ${EXEFILE} )
    IF( ${PROCS} STREQUAL "1" )
        CREATE_TEST_NAME( "${EXEFILE}_WEEKLY" ${ARGN} )
    ELSEIF( (USE_MPI OR USE_EXT_MPI) AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_WEEKLY" ${ARGN} )
    ENDIF()
    IF ( ${PROCS} GREATER ${TEST_MAX_PROCS} )
        MESSAGE("Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS})")
    ELSEIF( ${PROCS} STREQUAL "1" )
        CREATE_TEST_NAME( "${EXEFILE}_WEEKLY" ${ARGN} )
        IF ( USE_MPI_FOR_SERIAL_TESTS )
            ADD_TEST( NAME ${TESTNAME} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${LAST_TEST}> ${ARGN} )
            SET_PROPERTY( TEST ${TESTNAME} APPEND PROPERTY ENVIRONMENT OMPI_MCA_hwloc_base_binding_policy=none )
        ELSE()
            ADD_TEST( NAME ${TESTNAME} COMMAND $<TARGET_FILE:${LAST_TEST}> ${ARGN} )
        ENDIF()
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS 1 )
        ADD_RESOURCE_LOCK( ${TESTNAME} ${EXEFILE} ${ARGN} )
    ELSEIF( (USE_MPI OR USE_EXT_MPI) AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_WEEKLY" ${ARGN} )
        ADD_TEST( NAME ${TESTNAME} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} $<TARGET_FILE:${LAST_TEST}> ${ARGN} )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${PROCS} )
        SET_PROPERTY( TEST ${TESTNAME} APPEND PROPERTY ENVIRONMENT OMPI_MCA_hwloc_base_binding_policy=none )
        ADD_RESOURCE_LOCK( ${TESTNAME} ${EXEFILE} ${ARGN} )
    ENDIF()
ENDFUNCTION()


# Add a executable as a parallel test
FUNCTION( ADD_${PROJ}_TEST_PARALLEL EXEFILE PROCS ${ARGN} )
    # Check if we actually want to add the test
    KEEP_TEST( RESULT )
    IF ( NOT RESULT )
        RETURN()
    ENDIF()
    # Add the provisional test
    ADD_PROJ_PROVISIONAL_TEST( ${EXEFILE} )
    CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs" ${ARGN} )
    IF ( NOT ( USE_MPI OR USE_EXT_MPI ) )
        MESSAGE("Disabling test ${TESTNAME} (configured without MPI)")
    ELSEIF ( ${PROCS} GREATER ${TEST_MAX_PROCS} )
        MESSAGE("Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS})")
    ELSE()
        ADD_TEST( NAME ${TESTNAME} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} $<TARGET_FILE:${LAST_TEST}> ${ARGN} )
        SET_PROPERTY( TEST ${TESTNAME} APPEND PROPERTY ENVIRONMENT OMPI_MCA_hwloc_base_binding_policy=none )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${PROCS} )
        ADD_RESOURCE_LOCK( ${TESTNAME} ${EXEFILE} ${ARGN} )
    ENDIF()
ENDFUNCTION()


# Add a parallel test that may use both MPI and threads
# This allows us to correctly compute the number of processors used by the test
MACRO( ADD_${PROJ}_TEST_THREAD_MPI EXEFILE PROCS THREADS ${ARGN} )
    ADD_PROJ_PROVISIONAL_TEST( ${EXEFILE} )
    CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_${THREADS}threads" ${ARGN} )
    MATH( EXPR TOT_PROCS "${PROCS} * ${THREADS}" )
    IF ( ${TOT_PROCS} GREATER ${TEST_MAX_PROCS} )
        MESSAGE("Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS})")
    ELSEIF ( ( ${PROCS} STREQUAL "1" ) AND NOT USE_MPI_FOR_SERIAL_TESTS )
        ADD_TEST( NAME ${TESTNAME} COMMAND $<TARGET_FILE:${LAST_TEST}> ${ARGN} )
        SET_TESTS_PROPERTIES ( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${TOT_PROCS} )
        ADD_RESOURCE_LOCK( ${TESTNAME} ${EXEFILE} ${ARGN} )
    ELSEIF ( USE_MPI OR USE_EXT_MPI )
        ADD_TEST( NAME ${TESTNAME} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} $<TARGET_FILE:${LAST_TEST}> ${ARGN} )
        SET_PROPERTY( TEST ${TESTNAME} APPEND PROPERTY ENVIRONMENT OMPI_MCA_hwloc_base_binding_policy=none )
        SET_TESTS_PROPERTIES ( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${TOT_PROCS} )
        ADD_RESOURCE_LOCK( ${TESTNAME} ${EXEFILE} ${ARGN} )
    ENDIF()
ENDMACRO()


# Add a executable as an example
FUNCTION( ADD_${PROJ}_EXAMPLE EXEFILE PROCS ${ARGN} )
    # Add the file to the example doxygen file
    SET( VALUE 0 )
    FOREACH(_variableName ${EXAMPLE_LIST})
        IF ( "${_variableName}" STREQUAL "${EXEFILE}" )
            SET( VALUE 1 )
        ENDIF()
    ENDFOREACH()
    IF ( NOT ${VALUE} )
        FILE(APPEND ${EXAMPLE_INSTALL_DIR}/examples.h "* \\ref ${EXEFILE} \"${EXEFILE}\"\n" )
        SET( EXAMPLE_LIST ${EXAMPLE_LIST} ${EXEFILE} CACHE INTERNAL "example_list" FORCE )
    ENDIF()
    # Check if we actually want to add the test
    IF ( ONLY_BUILD_DOCS )
        RETURN()
    ENDIF()
    # Add the provisional test
    ADD_PROJ_PROVISIONAL_TEST( ${EXEFILE} )
    ADD_DEPENDENCIES( build-examples ${EXEFILE} )
    ADD_CUSTOM_COMMAND( TARGET ${EXEFILE} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${LAST_TEST}> "${EXAMPLE_INSTALL_DIR}/${EXEFILE}"
    )
    IF( ${PROCS} STREQUAL "1" AND (NOT USE_MPI_FOR_SERIAL_TESTS) )
        CREATE_TEST_NAME( "example--${EXEFILE}" ${ARGN} )
        ADD_TEST( NAME ${TESTNAME} COMMAND $<TARGET_FILE:${LAST_TEST}> ${ARGN} )
    ELSEIF ( USE_EXT_MPI AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "example--${EXEFILE}_${PROCS}procs" ${ARGN} )
        ADD_TEST( NAME ${TESTNAME} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} $<TARGET_FILE:${LAST_TEST}> ${ARGN} )
        SET_PROPERTY( TEST ${TESTNAME} APPEND PROPERTY ENVIRONMENT OMPI_MCA_hwloc_base_binding_policy=none )
    ENDIF()
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${PROCS} )
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES RESOURCE_LOCK ${EXEFILE} )
ENDFUNCTION()


# Begin configure for the examples for a package
MACRO( BEGIN_EXAMPLE_CONFIG PACKAGE )
    # Set example install dir
    SET( EXAMPLE_INSTALL_DIR ${${PROJ}_INSTALL_DIR}/examples/${PACKAGE} )
    # Create list of examples
    SET( EXAMPLE_LIST "dummy" CACHE INTERNAL "example_list" FORCE )
    # Create doxygen input file for examples
    SET( DOXYFILE_EXTRA_SOURCES ${DOXYFILE_EXTRA_SOURCES} ${EXAMPLE_INSTALL_DIR} CACHE INTERNAL "doxyfile_extra_sources") 
    FILE(WRITE  ${EXAMPLE_INSTALL_DIR}/examples.h "// Include file for doxygen providing the examples for ${PACKAGE}\n")
    FILE(APPEND ${EXAMPLE_INSTALL_DIR}/examples.h "/*! \\page Examples_${PACKAGE}\n" )
ENDMACRO()


# Install the examples
MACRO( INSTALL_${PROJ}_EXAMPLE PACKAGE )
    FILE(APPEND ${EXAMPLE_INSTALL_DIR}/examples.h "*/\n" )
    SET( EXAMPLE_INSTALL_DIR "" )
ENDMACRO()


# Macro to check if a flag is enabled
MACRO( CHECK_ENABLE_FLAG FLAG DEFAULT )
    IF( NOT DEFINED ${FLAG} )
        SET( ${FLAG} ${DEFAULT} )
    ELSEIF( ${FLAG}  STREQUAL "" )
        SET( ${FLAG} ${DEFAULT} )
    ELSEIF( ( ${${FLAG}} STREQUAL "FALSE" ) OR ( ${${FLAG}} STREQUAL "false" ) OR ( ${${FLAG}} STREQUAL "0" ) OR ( ${${FLAG}} STREQUAL "OFF" ) )
        SET( ${FLAG} 0 )
    ELSEIF( ( ${${FLAG}} STREQUAL "TRUE" ) OR ( ${${FLAG}} STREQUAL "true" ) OR ( ${${FLAG}} STREQUAL "1" ) OR ( ${${FLAG}} STREQUAL "ON" ) )
        SET( ${FLAG} 1 )
    ELSE()
        MESSAGE( "Bad value for ${FLAG} (${${FLAG}}); use true or false" )
    ENDIF()
ENDMACRO()


# Macro to add a latex file to the build
MACRO( ADD_LATEX_DOCS FILE )
    GET_FILENAME_COMPONENT(LATEX_TARGET ${FILE} NAME_WE)
    ADD_CUSTOM_TARGET( 
        ${LATEX_TARGET}_pdf 
        ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${FILE} ${CMAKE_CURRENT_BINARY_DIR}/.
        COMMAND pdflatex -interaction=batchmode -draftmode ${FILE} ";" echo ""
        COMMAND bibtex -terse ${LATEX_TARGET} ";" echo ""
        COMMAND pdflatex -interaction=batchmode ${FILE} ";" echo ""
        SOURCES ${FILE}
    )
    ADD_CUSTOM_COMMAND( 
        TARGET ${LATEX_TARGET}_pdf 
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/${LATEX_TARGET}.pdf ${${PROJ}_INSTALL_DIR}/doc/.
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    ADD_DEPENDENCIES( latex_docs ${LATEX_TARGET}_pdf )
ENDMACRO()


# Add a matlab mex file
FUNCTION( ADD_MATLAB_MEX SOURCE )
    STRING( REGEX REPLACE "[.]cpp" "" TARGET ${SOURCE} )
    STRING( REGEX REPLACE "[.]c"   "" TARGET ${TARGET} )
    MATLAB_ADD_MEX(
        NAME ${TARGET}
        SRC ${SOURCE}
    )
    TARGET_LINK_LIBRARIES( ${TARGET} ${MATLAB_TARGET} )
    ADD_PROJ_EXE_DEP( ${TARGET} )
    ADD_DEPENDENCIES( mex ${TARGET} )
    INSTALL( TARGETS ${TARGET} DESTINATION ${${PROJ}_INSTALL_DIR}/mex )
    ADD_DEPENDENCIES( mex ${TARGET} )
    SET( MEX_FILES2 ${MEX_FILES} "${${PROJ}_INSTALL_DIR}/mex/${TARGET}.${Matlab_MEX_EXTENSION}" )
    LIST( REMOVE_DUPLICATES MEX_FILES2 )
    SET( MEX_FILES ${MEX_FILES2} CACHE INTERNAL "" )
ENDFUNCTION()


# Add a matlab test
MACRO( ADD_MATLAB_TEST EXEFILE ${ARGN} )
    IF ( NOT MATLAB_EXE )
        MESSAGE( FATAL_ERROR "MATLAB_EXE not set, did you call CREATE_MATLAB_WRAPPER()" )
    ENDIF()
    CONFIGURE_FILE( ${EXEFILE}.m ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.m )
    CREATE_TEST_NAME( MATLAB ${EXEFILE} ${ARGN} )
    IF ( USING_MSVC )
        SET( MATLAB_OPTIONS "-logfile" "log_${EXEFILE}" )
    ENDIF()
    SET( MATLAB_COMMAND "try, ${EXEFILE}, catch ME, disp(getReport(ME)), clear all global, exit(1), end, disp('ALL TESTS PASSED'); exit(0)" )
    SET( MATLAB_DEBUGGER_OPTIONS )
    IF ( MATLAB_DEBUGGER )
        SET( MATLAB_DEBUGGER_OPTIONS -D${MATLAB_DEBUGGER} )
        STRING(REPLACE "\"" "" MATLAB_DEBUGGER_OPTIONS "${MATLAB_DEBUGGER_OPTIONS}")
    ENDIF()
    ADD_TEST( NAME ${TESTNAME} COMMAND "${MATLAB_EXE}" ${MATLAB_OPTIONS} -r "${MATLAB_COMMAND}" ${MATLAB_DEBUGGER_OPTIONS} )
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES PASS_REGULAR_EXPRESSION "ALL TESTS PASSED" FAIL_REGULAR_EXPRESSION "FAILED" )
    SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY ENVIRONMENT RATES_DIRECTORY=${RATES_DIRECTORY} )
ENDMACRO()


# Create a script to start matlab preloading libraries
FUNCTION( CREATE_MATLAB_WRAPPER )
    SET( tmp_libs ${MEX_LIBCXX} ${MEX_FILES} )
    STRING(REGEX REPLACE ";" ":" tmp_libs "${tmp_libs}")
    STRING(REGEX REPLACE ";" ":" tmp_path "${MATLABPATH}")
    IF ( USING_MSVC )
        # Create a matlab wrapper for windows
        SET( MATLAB_GUI "${CMAKE_CURRENT_BINARY_DIR}/tmp/matlab-gui.bat" )
        SET( MATLAB_CMD "${CMAKE_CURRENT_BINARY_DIR}/tmp/matlab-cmd.bat" )
        SET( MATLAB_INSTALL_CMD "matlab-cmd.bat" )
        FILE( WRITE "${MATLAB_GUI}" "@echo off\n")
        FILE( WRITE "${MATLAB_CMD}" "@echo off\n")
        FILE( APPEND "${MATLAB_GUI}" "matlab  -singleCompThread -nosplash %*\n")
        FILE( APPEND "${MATLAB_CMD}" "matlab  -singleCompThread -nosplash -nodisplay -nodesktop -wait %*\n")
        FILE( WRITE "${MATLAB_GUI}" "@echo on\n")
        FILE( WRITE "$${MATLAB_CMD}" "@echo on\n")
    ELSE()
        # Create a matlab wrapper for linux/mac
        SET( MATLAB_GUI "${CMAKE_CURRENT_BINARY_DIR}/tmp/matlab-gui" )
        SET( MATLAB_CMD "${CMAKE_CURRENT_BINARY_DIR}/tmp/matlab-cmd" )
        SET( MATLAB_INSTALL_CMD "matlab-cmd" )
        FILE( WRITE "${MATLAB_GUI}" "LD_PRELOAD=\"${tmp_libs}\" MKL_NUM_THREADS=1 MATLABPATH=\"${tmp_path}\" \"${Matlab_MAIN_PROGRAM}\" -singleCompThread -nosplash \"$@\"\n")
        FILE( WRITE "${MATLAB_CMD}" "LD_PRELOAD=\"${tmp_libs}\" MKL_NUM_THREADS=1 MATLABPATH=\"${tmp_path}\" \"${Matlab_MAIN_PROGRAM}\" -singleCompThread -nosplash -nodisplay -nojvm \"$@\"\n")
    ENDIF()
    FILE( COPY "${MATLAB_GUI}" DESTINATION "${${PROJ}_INSTALL_DIR}/mex"
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )
    FILE( COPY "${MATLAB_CMD}" DESTINATION "${${PROJ}_INSTALL_DIR}/mex"
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )
    SET( MATLAB_EXE "${${PROJ}_INSTALL_DIR}/mex/${MATLAB_INSTALL_CMD}" CACHE INTERNAL "" )
ENDFUNCTION()


# Macro to change the classification of a package
MACRO( SET_PACKAGE_CLASSIFICATION  PACKAGE_LIST  PACKAGE_NAME  CLASS )
    LIST(FIND ${PACKAGE_LIST} ${PACKAGE_NAME} PACKAGE_NAME_IDX)
    IF (PACKAGE_NAME_IDX EQUAL -1)
        MESSAGE(FATAL_ERROR "Package ${PACKAGE_NAME} not found in list of packages!")
    ELSE()
        MATH(EXPR PACKAGE_CLASSIFICATION_IDX "${PACKAGE_NAME_IDX}+2")
        LIST(INSERT ${PACKAGE_LIST} ${PACKAGE_CLASSIFICATION_IDX} ${CLASS})
        MATH(EXPR PACKAGE_CLASSIFICATION_IDX "${PACKAGE_CLASSIFICATION_IDX} + 1")
        LIST(REMOVE_AT ${PACKAGE_LIST} ${PACKAGE_CLASSIFICATION_IDX})
    ENDIF()
ENDMACRO()


# Macro to "disable" a package on the given platform (this mearly changes it to experimental)
MACRO( PACKAGE_DISABLE_ON_PLATFORMS  PACKAGE_LIST  PACKAGE_NAME )
    FOREACH(HOSTTYPE ${ARGN})
        IF (${PROJECT_NAME}_HOSTTYPE STREQUAL ${HOSTTYPE})
            SET_PACKAGE_CLASSIFICATION(${PACKAGE_LIST} ${PACKAGE_NAME} EX)
            IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
                MESSAGE(
                  "\n***"
                  "\n*** WARNING: User has set ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON but the"
                  "\n*** package ${PACKAGE_NAME} is not supported on this platform type '${HOSTTYPE}'!"
                  "\n***\n"
               )
            ENDIF()
        ENDIF()
    ENDFOREACH()
ENDMACRO()


# Append a list to a file
FUNCTION( APPEND_LIST FILENAME VARS PREFIX POSTFIX )
    FOREACH( tmp ${VARS} )
        FILE( APPEND "${FILENAME}" "${PREFIX}" )
        FILE( APPEND "${FILENAME}" "${tmp}" )
        FILE( APPEND "${FILENAME}" "${POSTFIX}" )
    ENDFOREACH ()
ENDFUNCTION()


# add custom target distclean
# cleans and removes cmake generated files etc.
MACRO( ADD_DISTCLEAN ${ARGN} )
    SET(DISTCLEANED
        cmake.depends
        cmake.check_depends
        CMakeCache.txt
        CMakeFiles
        CMakeTmp
        CMakeDoxy*
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
        docs
        examples
        latex_docs
        lib
        Makefile.config
        install_manifest.txt
        test
        matlab
        Matlab
        mex
        tmp
        #tmp#
        bin
        cmake
        cppclean
        compile_commands.json
        ${ARGN}
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
            Debug
        )
        SET( DISTCLEAN_FILE "${CMAKE_CURRENT_BINARY_DIR}/distclean.bat" )
        FILE( WRITE  "${DISTCLEAN_FILE}" "del /s /q /f " )
        APPEND_LIST( "${DISTCLEAN_FILE}" "${DISTCLEANED}" " " " " )
        FILE( APPEND "${DISTCLEAN_FILE}" "\n" )
        APPEND_LIST( "${DISTCLEAN_FILE}" "${DISTCLEANED}" "for /d %%x in ("   ") do rd /s /q \"%%x\"\n" )
        ADD_CUSTOM_COMMAND(
            DEPENDS clean
            COMMENT "distribution clean"
            COMMAND distclean.bat & del /s/q/f distclean.bat
            TARGET  distclean
        )
    ENDIF()
ENDMACRO()


# add custom target mex_clean
MACRO( ADD_MEXCLEAN )
    IF (UNIX)
        ADD_CUSTOM_TARGET( mexclean 
            COMMENT "mex clean"
            COMMAND rm
            ARGS    -Rf libmatlab.* *.mex* test/*.mex*
        )
    ENDIF(UNIX)
ENDMACRO()


# Print the current repo version and create target to write to a file
SET( WriteRepoVersionCmakeFile "${CMAKE_CURRENT_LIST_DIR}/WriteRepoVersion.cmake" )
FUNCTION( WRITE_REPO_VERSION FILENAME )
    SET( CMD ${CMAKE_COMMAND} -Dfilename="${FILENAME}" -Dsrc_dir="${${PROJ}_SOURCE_DIR}" 
             -Dtmp_file="${CMAKE_CURRENT_BINARY_DIR}/tmp/version.h" -DPROJ=${PROJ} 
             -P "${WriteRepoVersionCmakeFile}" )
    EXECUTE_PROCESS( COMMAND ${CMD} )
    ADD_CUSTOM_TARGET( write_repo_version  COMMENT "Write repo version"  COMMAND ${CMD} )
ENDFUNCTION()


# Check if we can include a python module
FUNCTION( FIND_PYTHON_MODULE MODULE)
	STRING(TOUPPER ${MODULE} MODULE2)
    SET( PY_${MODULE}_FOUND FALSE )
	IF ( NOT PY_${MODULE2} )
		IF ( ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED" )
			SET(PY_${MODULE2}_FIND_REQUIRED TRUE)
		ENDIF()
        SET( CMD "import ${MODULE}; print('Success')" )
		EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}" "-c" "${CMD}"
			RESULT_VARIABLE _${MODULE2}_status
			OUTPUT_VARIABLE _${MODULE2}_output
			ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE )
        IF ( "${_${MODULE2}_output}" STREQUAL "Success" )
            SET( PY_${MODULE2}_FOUND TRUE )
        ENDIF()
    ELSE()
        SET( PY_${MODULE2}_FOUND TRUE )
	ENDIF()
    IF ( PY_${MODULE2}_FOUND )
        MESSAGE( STATUS "Performing Test PYTHON_${MODULE2} - Success" )
    ELSE()
        MESSAGE( STATUS "Performing Test PYTHON_${MODULE2} - Failure" )
    ENDIF()
    IF ( NOT PY_${MODULE2}_FOUND AND PY_${MODULE2}_FIND_REQUIRED )
        MESSAGE( FATAL_ERROR "Unable to find required python module: ${MODULE2}" )
    ENDIF()
	SET( PY_${MODULE2}_FOUND ${PY_${MODULE2}_FOUND} PARENT_SCOPE )
ENDFUNCTION(FIND_PYTHON_MODULE)
