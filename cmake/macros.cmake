INCLUDE(CheckCCompilerFlag)
INCLUDE(CheckCSourceCompiles)
INCLUDE(CheckCXXCompilerFlag)
INCLUDE(CheckCXXSourceCompiles)
IF ( NOT TEST_FAIL_REGULAR_EXPRESSION )
    SET( TEST_FAIL_REGULAR_EXPRESSION "(FAILED)|(leaked context IDs detected)|(handles are still allocated)" )
ENDIF()


# Check that the PROJ and ${PROJ}_INSTALL_DIR variables are set 
# These variables are used to generate the ADD_PROJ_TEST macros
IF ( NOT PROJ )
    MESSAGE(FATAL_ERROR "PROJ must be set before including macros.cmake")
ENDIF()
IF ( NOT ${PROJ}_INSTALL_DIR )
    MESSAGE(FATAL_ERROR "${PROJ}_INSTALL_DIR must be set before including macros.cmake")
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
    ADD_SUBDIRECTORY( ${PACKAGE} )
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
    SET( SOURCES ${SOURCES} ${T_CXXSOURCES} ${T_CSOURCES} ${T_FSOURCES} )
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
    SET( COPY_TARGET copy-${PROJ}-include )
    IF( NOT TARGET ${COPY_TARGET} )
        ADD_CUSTOM_TARGET( ${COPY_TARGET} ALL )
        ADD_DEPENDENCIES( copy-${PROJ}-include ${COPY_TARGET} )
    ENDIF()
    # Copy the header files to the include path
    FILE( GLOB HFILES RELATIVE "${${PROJ}_SOURCE_DIR}" ${HEADERS} )
    FOREACH( HFILE ${HFILES} )
        SET( SRC_FILE "${${PROJ}_SOURCE_DIR}/${HFILE}" )
        SET( DST_FILE "${${PROJ}_INSTALL_DIR}/include/${HFILE}" )
        # Only copy the headers if the exisit in the project source directory
        IF ( EXISTS "${SRC_FILE}" )
            ADD_CUSTOM_COMMAND(TARGET ${COPY_TARGET} 
                PRE_BUILD 
                COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SRC_FILE}" "${DST_FILE}"
                DEPENDS "${SRC_FILE}"
            )
        ENDIF()
    ENDFOREACH()
    # Add the library and install the package
    IF ( NOT ONLY_BUILD_DOCS AND ( SOURCES OR CUDASOURCES ) )
        IF( USE_CUDA )
            CUDA_COMPILE( CUBINS ${CUDASOURCES} )
        ENDIF()
        ADD_LIBRARY( ${PACKAGE} ${LIB_TYPE} ${SOURCES} ${CUBINS} )
        TARGET_LINK_LIBRARIES( ${PACKAGE} ${SYSTEM_LDFLAGS} )
        IF ( TARGET write_repo_version )
            ADD_DEPENDENCIES( ${PACKAGE} write_repo_version )
        ENDIF()
        ADD_DEPENDENCIES ( ${PACKAGE} copy-${PROJ}-include )
        INSTALL( TARGETS ${PACKAGE} DESTINATION ${${PROJ}_INSTALL_DIR}/lib )
    ELSE()
        ADD_CUSTOM_TARGET( ${PACKAGE} ALL )
    ENDIF()
    # Clear the sources
    SET( HEADERS "" )
    SET( CSOURCES "" )
    SET( CXXSOURCES "" )
ENDMACRO()


# Macro to verify that a variable has been set
MACRO( VERIFY_VARIABLE VARIABLE_NAME )
    IF ( NOT ${VARIABLE_NAME} )
        MESSAGE( FATAL_ERROR "PLease set: " ${VARIABLE_NAME} )
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
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
    set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
    set(CMAKE_SHARED_LINKER_FLAGS)
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
    # Add the static flag if necessary
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static") # Add static flag
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -static ") 
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static ")
    SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")                # Add static flag
    SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-static")              # Add static flag
ENDMACRO()


# Macro to identify the compiler
MACRO( IDENTIFY_COMPILER )
    # SET the C/C++ compiler
    IF ( CMAKE_C_COMPILER_WORKS OR CMAKE_C_COMPILER_WORKS )
        IF( CMAKE_COMPILE_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
            SET( USING_GCC TRUE )
            ADD_DEFINITIONS( -DUSING_GCC )
            MESSAGE("Using gcc")
        ELSEIF( MSVC OR MSVC_IDE OR MSVC60 OR MSVC70 OR MSVC71 OR MSVC80 OR CMAKE_COMPILER_2005 OR MSVC90 OR MSVC10 )
            IF( NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows" )
                MESSAGE( FATAL_ERROR "Using microsoft compilers on non-windows system?" )
            ENDIF()
            SET( USING_MSVC TRUE )
            ADD_DEFINITIONS( -DUSING_MSVC )
            MESSAGE("Using Microsoft")
        ELSEIF( (${CMAKE_C_COMPILER_ID} MATCHES "Intel") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Intel") ) 
            SET(USING_ICC TRUE)
            ADD_DEFINITIONS( -DUSING_ICC )
            MESSAGE("Using icc")
        ELSEIF( (${CMAKE_C_COMPILER_ID} MATCHES "PGI") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "PGI") )
            SET(USING_PGCC TRUE)
            ADD_DEFINITIONS( -DUSING_PGCC )
            MESSAGE("Using pgCC")
        ELSEIF( (${CMAKE_C_COMPILER_ID} MATCHES "CRAY") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "CRAY") OR
                (${CMAKE_C_COMPILER_ID} MATCHES "Cray") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Cray") )
            SET(USING_CRAY TRUE)
            ADD_DEFINITIONS( -DUSING_CRAY )
            MESSAGE("Using Cray")
        ELSEIF( (${CMAKE_C_COMPILER_ID} MATCHES "CLANG") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "CLANG") OR
                (${CMAKE_C_COMPILER_ID} MATCHES "Clang") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang") )
            SET(USING_CLANG TRUE)
            ADD_DEFINITIONS( -DUSING_CLANG )
            MESSAGE("Using Clang")
        ELSE()
            SET(USING_DEFAULT TRUE)
            MESSAGE("${CMAKE_C_COMPILER_ID}")
            MESSAGE(WARNING "Unknown C/C++ compiler, default flags will be used")
        ENDIF()
    ENDIF()
    # SET the Fortran++ compiler
    IF ( CMAKE_Fortran_COMPILER_WORKS )
        IF( CMAKE_COMPILE_IS_GFORTRAN OR (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU") )
            SET( USING_GFORTRAN TRUE )
            MESSAGE("Using gfortran")
        ELSEIF ( (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel") ) 
            SET(USING_IFORT TRUE)
            MESSAGE("Using ifort")
        ELSEIF ( ${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
            SET(USING_PGF90 TRUE)
            MESSAGE("Using pgf90")
        ELSE()
            SET(USING_DEFAULT TRUE)
            MESSAGE("${CMAKE_Fortran_COMPILER_ID}")
            MESSAGE("Unknown Fortran compiler, default flags will be used")
        ENDIF()
    ENDIF()
ENDMACRO()


# Macro to set the proper warnings
MACRO( SET_WARNINGS )
  IF ( USING_GCC )
    # Add gcc specific compiler options
    # Note: adding -Wlogical-op causes a wierd linking error on Titan using the nvcc wrapper:
    #    /usr/bin/ld: cannot find gical-op: No such file or directory
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -Wall -Wextra -Wcast-align -Wlogical-op -Wno-char-subscripts -Wno-unused-parameter -Wno-unused-variable") 
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wcast-align -Wlogical-op -Wno-char-subscripts -Wno-unused-parameter -Wno-unused-variable")
  ELSEIF ( USING_MSVC )
    # Add Microsoft specifc compiler options
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} /D _SCL_SECURE_NO_WARNINGS /D _CRT_SECURE_NO_WARNINGS /D _ITERATOR_DEBUG_LEVEL=0" )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _SCL_SECURE_NO_WARNINGS /D _CRT_SECURE_NO_WARNINGS /D _ITERATOR_DEBUG_LEVEL=0" )
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
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -Wall -Wno-char-subscripts -Wno-unused-private-field -Wno-unused-variable")
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall -Wno-char-subscripts -Wno-unused-private-field -Wno-unused-variable")
  ELSEIF ( USING_DEFAULT )
    # Add default compiler options
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS}")
  ENDIF()
ENDMACRO()


# Macro to add user compile flags
MACRO( ADD_USER_FLAGS )
    STRING( STRIP "${CMAKE_C_FLAGS} ${CFLAGS} ${CFLAGS_EXTRA}" CMAKE_C_FLAGS )
    STRING( STRIP "${CMAKE_CXX_FLAGS} ${CXXFLAGS} ${CXXFLAGS_EXTRA}" CMAKE_CXX_FLAGS )
    STRING( STRIP "${CMAKE_Fortran_FLAGS} ${FFLAGS} ${FFLAGS_EXTRA}" CMAKE_Fortran_FLAGS )
ENDMACRO()


# Macro to add user c++ std 
MACRO( ADD_CXX_STD )
    IF ( NOT CXX_STD )
        MESSAGE( FATAL_ERROR "The desired c++ standard must be specified: CXX_STD=(98,11,14,NONE)")
    ENDIF()
    IF ( ${CXX_STD} STREQUAL "NONE" )
        # Do nothing
        return()
    ELSEIF ( (NOT ${CXX_STD} STREQUAL "98") AND (NOT ${CXX_STD} STREQUAL "11") AND (NOT ${CXX_STD} STREQUAL "14") )
        MESSAGE( FATAL_ERROR "Unknown c++ standard (98,11,14,NONE)" )
    ENDIF()
    # Add the flags
    SET( CMAKE_CXX_STANDARD ${CXX_STD} )
    MESSAGE( "C++ standard: ${CXX_STD}" )
    SET( CXX_STD_FLAG )
    IF ( USING_GCC )
        # GNU: -std=
        IF ( ${CXX_STD} STREQUAL "98" )
            SET( CXX_STD_FLAG -std=c++98 )
        ELSEIF ( ${CXX_STD} STREQUAL "11" )
            SET( CXX_STD_FLAG -std=c++11 )
        ELSEIF ( ${CXX_STD} STREQUAL "14" )
            SET( CXX_STD_FLAG -std=c++1y )
        ELSE()
            MESSAGE(FATAL_ERROR "Unknown standard")
        ENDIF()
    ELSEIF ( USING_MSVC )
        # Microsoft: Does not support this level of control
    ELSEIF ( USING_ICC )
        # ICC: -std=
        SET( CXX_STD_FLAG -std=c++${CXX_STD} )
    ELSEIF ( USING_CRAY )
        # Cray: Does not seem to support controlling the std?
    ELSEIF ( USING_PGCC )
        # PGI: -std=
        IF ( ${CXX_STD} STREQUAL "98" )
            SET( CXX_STD_FLAG --c++0x )
        ELSEIF ( ${CXX_STD} STREQUAL "11" )
            SET( CXX_STD_FLAG --c++11 )
        ELSEIF ( ${CXX_STD} STREQUAL "14" )
            MESSAGE( FATAL_ERROR "C++14 features are not availible yet for PGI" )
        ELSE()
            MESSAGE(FATAL_ERROR "Unknown standard")
        ENDIF()
    ELSEIF ( USING_CLANG )
        # Clang: -std=
        IF ( ( ${CXX_STD} STREQUAL "98") OR ( ${CXX_STD} STREQUAL "11" ) )
            SET( CXX_STD_FLAG -std=c++${CXX_STD} )
        ELSEIF ( ${CXX_STD} STREQUAL "14" )
            SET( CXX_STD_FLAG -std=c++1y )
        ELSE()
            MESSAGE(FATAL_ERROR "Unknown standard")
        ENDIF()
    ELSEIF ( USING_DEFAULT )
        # Default: do nothing
    ENDIF()
    ADD_DEFINITIONS( -DCXX_STD=${CXX_STD} )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_STD_FLAG}")
ENDMACRO()


# Macro to set the compile/link flags
MACRO( SET_COMPILER_FLAGS )
    # Initilaize the compiler
    IDENTIFY_COMPILER()
    # Set the default flags for each build type
    IF ( USING_MSVC )
        SET(CMAKE_C_FLAGS_DEBUG   "-D_DEBUG /DEBUG /Od /EHsc /MDd /Zi /Z7" )
        SET(CMAKE_CXX_FLAGS_DEBUG "-D_DEBUG /DEBUG /Od /EHsc /MDd /Zi /Z7" )
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
    IF ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
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
    # Add the c++ standard flags
    ADD_CXX_STD()
    # Set the warnings to use
    SET_WARNINGS()
    # Test the compile flags
    CHECK_C_COMPILER_FLAG( "${CMAKE_C_FLAGS}" CHECK_C_FLAGS )
    CHECK_CXX_COMPILER_FLAG( "${CMAKE_CXX_FLAGS}" CHECK_CXX_FLAGS )
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
    TARGET_LINK_LIBRARIES( ${EXE} ${${PROJ}_LIBS} ${${PROJ}_LIBS} )
    TARGET_LINK_LIBRARIES( ${EXE} ${${PROJECT_NAME}_LIBRARIES} )
    # Add external libraries
    SET_TARGET_PROPERTIES( ${EXE} PROPERTIES LINK_FLAGS "${LDFLAGS}" )
    SET_TARGET_PROPERTIES( ${EXE} PROPERTIES LINK_FLAGS "${LDFLAGS_EXTRA}" )
    TARGET_LINK_LIBRARIES( ${EXE} ${LINK_LIBRARIES} )
    TARGET_LINK_LIBRARIES( ${EXE} ${CUDA_LIBS} )
    TARGET_LINK_LIBRARIES( ${EXE} ${EXTERNAL_LIBS} )
    IF ( USE_MPI OR USE_EXT_MPI OR HAVE_MPI )
        TARGET_LINK_LIBRARIES( ${EXE} ${MPI_LINK_FLAGS} ${MPI_LIBRARIES} )
    ENDIF()
    TARGET_LINK_LIBRARIES( ${EXE} ${LAPACK_LIBS} ${BLAS_LIBS} ${BLAS_LAPACK_LIBS} )
    TARGET_LINK_LIBRARIES( ${EXE} ${COVERAGE_LIBS} ${LDLIBS} ${LDLIBS_EXTRA} )
    TARGET_LINK_LIBRARIES( ${EXE} ${SYSTEM_LIBS} ${SYSTEM_LDFLAGS} )
ENDMACRO()


# Check if we want to keep the test
FUNCTION( KEEP_TEST RESULT )
    SET( ${RESULT} 1 PARENT_SCOPE )
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
    # Check if test has already been added
    SET( tmp )
    IF ( TARGET ${EXEFILE} )
        GET_TARGET_PROPERTY(tmp ${EXEFILE} LOCATION)
        STRING(REGEX REPLACE "//" "/" tmp "${tmp}" )        
    ENDIF()
    IF ( NOT tmp )
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
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/$(Configuration)/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/$(CONFIGURATION)$(EFFECTIVE_PLATFORM_NAME)/${EXEFILE}" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/$(OutDir)/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSE()
        # We are trying to add 2 different tests with the same name
        MESSAGE( "Existing test: ${tmp}" )
        MESSAGE( "New test:      ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}" )
        MESSAGE( FATAL_ERROR "Trying to add 2 different tests with the same name" )
    ENDIF()
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
    IF ( NOT ARGN )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES RESOURCE_LOCK ${EXEFILE} )
    ELSE()
        FOREACH( tmp ${ARGN} )
            SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES RESOURCE_LOCK ${tmp} )
        ENDFOREACH()
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
    ADD_PROJ_PROVISIONAL_TEST ( ${EXEFILE} )
    CREATE_TEST_NAME( ${EXEFILE} ${ARGN} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    IF ( USE_MPI_FOR_SERIAL_TESTS )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${EXE} ${ARGN} )
    ELSE()
        ADD_TEST( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
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
    ADD_PROJ_PROVISIONAL_TEST ( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
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
            ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${EXE} ${ARGN} )
        ELSE()
            ADD_TEST( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
        ENDIF()
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS 1 )
        ADD_RESOURCE_LOCK( ${TESTNAME} ${EXEFILE} ${ARGN} )
    ELSEIF( (USE_MPI OR USE_EXT_MPI) AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_WEEKLY" ${ARGN} )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${PROCS} )
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
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs" ${ARGN} )
    IF ( NOT ( USE_MPI OR USE_EXT_MPI ) )
        MESSAGE("Disabling test ${TESTNAME} (configured without MPI)")
    ELSEIF ( ${PROCS} GREATER ${TEST_MAX_PROCS} )
        MESSAGE("Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS})")
    ELSE()
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${PROCS} )
        ADD_RESOURCE_LOCK( ${TESTNAME} ${EXEFILE} ${ARGN} )
    ENDIF()
ENDFUNCTION()


# Add a parallel test that may use both MPI and threads
# This allows us to correctly compute the number of processors used by the test
MACRO( ADD_${PROJ}_TEST_THREAD_MPI EXEFILE PROCS THREADS ${ARGN} )
    ADD_PROJ_PROVISIONAL_TEST( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_${THREADS}threads" ${ARGN} )
    MATH( EXPR TOT_PROCS "${PROCS} * ${THREADS}" )
    IF ( ${TOT_PROCS} GREATER ${TEST_MAX_PROCS} )
        MESSAGE("Disabling test ${TESTNAME} (exceeds maximum number of processors ${TEST_MAX_PROCS}")
    ELSEIF ( ( ${PROCS} STREQUAL "1" ) AND NOT USE_MPI_FOR_SERIAL_TESTS )
        ADD_TEST ( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
        SET_TESTS_PROPERTIES ( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}" PROCESSORS ${TOT_PROCS} )
        ADD_RESOURCE_LOCK( ${TESTNAME} ${EXEFILE} ${ARGN} )
    ELSEIF ( USE_MPI OR USE_EXT_MPI )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
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
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    ADD_DEPENDENCIES( build-examples ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE2 ${EXEFILE} LOCATION)
    ADD_CUSTOM_COMMAND( TARGET ${EXEFILE} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${EXE2} "${EXAMPLE_INSTALL_DIR}/${EXEFILE}"
    )
    IF( ${PROCS} STREQUAL "1" AND (NOT USE_MPI_FOR_SERIAL_TESTS) )
        CREATE_TEST_NAME( "example--${EXEFILE}" ${ARGN} )
        ADD_TEST( ${TESTNAME} ${EXE} ${ARGN} )
    ELSEIF ( USE_EXT_MPI AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "example--${EXEFILE}_${PROCS}procs" ${ARGN} )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
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
    ELSEIF( ( ${${FLAG}} STREQUAL "false" ) OR ( ${${FLAG}} STREQUAL "0" ) OR ( ${${FLAG}} STREQUAL "OFF" ) )
        SET( ${FLAG} 0 )
    ELSEIF( ( ${${FLAG}} STREQUAL "true" ) OR ( ${${FLAG}} STREQUAL "1" ) OR ( ${${FLAG}} STREQUAL "ON" ) )
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
FUNCTION( ADD_MATLAB_MEX MEXFILE )
    # Set the MEX compiler and default link flags
    IF ( USING_MSVC )
        SET( MEX "\"${MATLAB_DIRECTORY}/bin/mex.bat\"" )
        SET( MEX_LDFLAGS ${MEX_LDFLAGS} -L${CMAKE_CURRENT_BINARY_DIR}/.. 
            -L${CMAKE_CURRENT_BINARY_DIR}/../Debug -L${CMAKE_CURRENT_BINARY_DIR} )
    ELSE()
        SET( MEX "${MATLAB_DIRECTORY}/bin/mex" )
        SET( MEX_LDFLAGS ${MEX_LDFLAGS} )
    ENDIF()
    SET( MEX_LDFLAGS ${MEX_LDFLAGS} -L${CMAKE_CURRENT_BINARY_DIR}/.. ${SYSTEM_LDFLAGS} )
    SET( MEX_LIBS    ${MEX_LIBS}    ${SYSTEM_LIBS}    )
    IF ( NOT USING_MSVC )
        FOREACH( rpath ${CMAKE_INSTALL_RPATH} )
            SET( MEX_LDFLAGS ${MEX_LDFLAGS} "-Wl,-rpath,${rpath},--no-undefined" )
        ENDFOREACH()
    ENDIF()
    # Create the mex comamnd
    STRING(REGEX REPLACE "[.]cpp" "" TARGET ${MEXFILE})
    STRING(REGEX REPLACE "[.]c"   "" TARGET ${TARGET})
    FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/temp/${TARGET}" )
    IF ( USING_MSVC )
        SET( MEX_BAT_FILE "${CMAKE_CURRENT_BINARY_DIR}/temp/${TARGET}/compile_mex.bat" )
        FILE( WRITE  "${MEX_BAT_FILE}" "${MEX} \"${CMAKE_CURRENT_SOURCE_DIR}/${MEXFILE}\" " )
        APPEND_LIST( "${MEX_BAT_FILE}" "${MEX_FLAGS}"   "\"" "\" " )
        APPEND_LIST( "${MEX_BAT_FILE}" "${COMPFLAGS}"   "\"" "\" " )
        APPEND_LIST( "${MEX_BAT_FILE}" "${MEX_INCLUDE}" "\"" "\" " )
        APPEND_LIST( "${MEX_BAT_FILE}" "${MEX_LDFLAGS}" "\"" "\" " )
        APPEND_LIST( "${MEX_BAT_FILE}" "${MEX_LIBS}"    "\"" "\" " )
        APPEND_LIST( "${MEX_BAT_FILE}" "${COVERAGE_MATLAB_LIBS}" "\"" "\" " )
        SET( MEX_COMMAND "${MEX_BAT_FILE}" )
    ELSE()
        STRING(REPLACE " " ";" MEX_CFLAGS "$$CFLAGS ${CMAKE_C_FLAGS}")
        STRING(REPLACE " " ";" MEX_CXXFLAGS "$$CXXFLAGS ${CMAKE_CXX_FLAGS}")
        SET( MEX_COMMAND ${MEX} ${CMAKE_CURRENT_SOURCE_DIR}/${MEXFILE} 
            ${MEX_FLAGS} ${MEX_INCLUDE} 
            CFLAGS="${MEX_CFLAGS}" 
            CXXFLAGS="${MEX_CXXFLAGS}" 
            LDFLAGS="${MEX_LDFLAGS}" 
            ${MEX_LIBS} ${COVERAGE_MATLAB_LIBS} 
        )
    ENDIF()
    ADD_CUSTOM_COMMAND( 
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.${MEX_EXTENSION}
        COMMAND ${MEX_COMMAND}
        COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_BINARY_DIR}/temp/${TARGET}/${TARGET}.${MEX_EXTENSION}" "${CMAKE_CURRENT_BINARY_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_BINARY_DIR}/temp/${TARGET}/${TARGET}.${MEX_EXTENSION}" "${${PROJ}_INSTALL_DIR}/matlab"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/temp/${TARGET}
        DEPENDS ${${PROJ}_LIBS} ${MATLAB_TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/${MEXFILE}
    )
    ADD_CUSTOM_TARGET( ${TARGET}
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.${MEX_EXTENSION}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${MEXFILE}
    )
    ADD_DEPENDENCIES( ${TARGET} ${${PROJ}_LIBS} ${MATLAB_TARGET} )
    ADD_DEPENDENCIES( mex ${TARGET} )
    INSTALL( FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.${MEX_EXTENSION} DESTINATION ${${PROJ}_INSTALL_DIR}/mex )
ENDFUNCTION()


# Add a matlab test
MACRO( ADD_MATLAB_TEST EXEFILE ${ARGN} )
    FIND_PROGRAM( MATLAB_EXE NAME matlab ) 
    CONFIGURE_FILE( ${EXEFILE}.m ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.m )
    CREATE_TEST_NAME( MATLAB ${EXEFILE} ${ARGN} )
    IF ( USING_MSVC )
        SET( MATLAB_OPTIONS "-nodisplay" "-nodesktop" "-wait" "-logfile" "log_${EXEFILE}" )
    ELSE()
        SET( MATLAB_OPTIONS "-nojvm -nosplash -nodisplay" )
    ENDIF()
    SET( MATLAB_COMMAND "addpath('${${PROJ}_BINARY_DIR}/matlab'); try, ${EXEFILE}, catch ME, ME, exit(1), end, disp('ALL TESTS PASSED'); exit(0)" )
    SET( MATLAB_DEBUGGER_OPTIONS )
    IF ( MATLAB_DEBUGGER )
        SET( MATLAB_DEBUGGER_OPTIONS -D${MATLAB_DEBUGGER} )
        STRING(REPLACE "\"" "" MATLAB_DEBUGGER_OPTIONS "${MATLAB_DEBUGGER_OPTIONS}")
    ENDIF()
    IF ( USING_MSVC )
        FILE( WRITE  "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.bat" "@echo off\n")
        FILE( APPEND "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.bat" "matlab ")
        APPEND_LIST( "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.bat" "${MATLAB_OPTIONS}" "" " " )
        FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.bat " -r \"${MATLAB_COMMAND}\" ${MATLAB_DEBUGGER_OPTIONS}\n" )
        FILE(APPEND  ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.bat "@echo on\n")
        FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.bat "type log_${EXEFILE}\n" )
        ADD_TEST( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.bat )
    ELSE()
        ADD_TEST( ${TESTNAME} matlab ${MATLAB_OPTIONS} -r "${MATLAB_COMMAND}" ${MATLAB_DEBUGGER_OPTIONS} )
    ENDIF()
    SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES PASS_REGULAR_EXPRESSION "ALL TESTS PASSED" FAIL_REGULAR_EXPRESSION "FAILED" )
    SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY ENVIRONMENT RATES_DIRECTORY=${RATES_DIRECTORY} )
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
        latex_docs
        lib
        Makefile.config
        install_manifest.txt
        test
        matlab
        mex
        tmp
        #tmp#
        bin
        cmake
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


