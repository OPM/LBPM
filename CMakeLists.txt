# Set some CMake properties    
CMAKE_MINIMUM_REQUIRED( VERSION 3.9 )
if( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.20.0") 
    CMAKE_POLICY( SET CMP0115 OLD )
endif()


MESSAGE("====================")
MESSAGE("Configuring LBPM-WIA")
MESSAGE("====================")


# Set the project name
SET( PROJ LBPM )          # Set the project name for CMake
SET( LBPM_LIB lbpm-wia )  # Set the final library name
SET( LBPM_INC  )          # Set an optional subfolder for includes (e.g. include/name/...)


# Initialize the project
PROJECT( ${PROJ} LANGUAGES CXX)

# Prevent users from building in place
IF ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}" )
    MESSAGE( FATAL_ERROR "Building code in place is a bad idea" )
ENDIF()


# Set the default C++ standard
SET( CMAKE_CXX_EXTENSIONS OFF )
IF ( NOT CMAKE_CXX_STANDARD )
    IF ( CXX_STD )
        MESSAGE( FATAL_ERROR "CXX_STD is obsolete, please set CMAKE_CXX_STANDARD" )
    ENDIF()
    SET( CMAKE_CXX_STANDARD 14 )
ENDIF()
IF ( ( "${CMAKE_CXX_STANDARD}" GREATER "90" ) OR ( "${CMAKE_CXX_STANDARD}" LESS "14" ) )
    MESSAGE( FATAL_ERROR "C++14 or newer required" )
ENDIF()


# Set source/install paths
SET( ${PROJ}_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )
SET( ${PROJ}_BUILD_DIR  "${CMAKE_CURRENT_BINARY_DIR}" )
IF( ${PROJ}_INSTALL_DIR )
    SET( ${PROJ}_INSTALL_DIR "${${PROJ}_INSTALL_DIR}" )
ELSEIF( PREFIX )
    SET( ${PROJ}_INSTALL_DIR "${PREFIX}" )
ELSEIF( NOT ${PROJ}_INSTALL_DIR )
    SET( ${PROJ}_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}" )
ENDIF()
INCLUDE_DIRECTORIES( "${${PROJ}_INSTALL_DIR}/include" )
SET( CMAKE_MODULE_PATH ${${PROJ}_SOURCE_DIR} ${${PROJ}_SOURCE_DIR}/cmake )


# Include macros
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/macros.cmake" )
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/libraries.cmake" )
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/LBPM-macros.cmake" )


# Check if we are only compiling docs
CHECK_ENABLE_FLAG( ONLY_BUILD_DOCS 0 )


# Set testing paramaters
SET( DROP_METHOD "http" )
SET( DROP_SITE "" )
SET( DROP_LOCATION "/CDash/submit.php?project=LBPM-WIA" )
SET( TRIGGER_SITE "" )
SET( DROP_SITE_CDASH TRUE )
ENABLE_TESTING()
INCLUDE( CTest )


# Check the compile mode and compile flags
IF ( NOT ONLY_BUILD_DOCS )
    CONFIGURE_SYSTEM()
ENDIF()


# Add some directories to include
INCLUDE_DIRECTORIES( "${${PROJ}_INSTALL_DIR}/include" )


# Create the target for documentation
ADD_CUSTOM_TARGET( doc )
ADD_CUSTOM_TARGET( latex_docs )
CHECK_ENABLE_FLAG( USE_DOXYGEN 1 )
CHECK_ENABLE_FLAG( USE_LATEX 1 )
FILE( MAKE_DIRECTORY "${${PROJ}_INSTALL_DIR}/doc" )
IF ( USE_DOXYGEN )
    SET( DOXYFILE_LATEX NO )
    SET( DOXYFILE_IN "${${PROJ}_SOURCE_DIR}/doxygen/Doxyfile.in" )
    SET( DOXY_HEADER_FILE "${${PROJ}_SOURCE_DIR}/doxygen/html/header.html" )
    SET( DOXY_FOOTER_FILE "${${PROJ}_SOURCE_DIR}/doxygen/html/footer.html" )
    SET( DOXYFILE_OUTPUT_DIR "${${PROJ}_INSTALL_DIR}/doc" )
    SET( DOXYFILE_SRC_HTML_DIR "${${PROJ}_SOURCE_DIR}/doxygen/html" )
    SET( DOXYFILE_SOURCE_DIR "${${PROJ}_SOURCE_DIR}" )
    SET( REL_PACKAGE_HTML "" )
    SET( DOXYGEN_MACROS "" )
    MESSAGE("DOXYGEN_MACROS = ${DOXYGEN_MACROS}")
    INCLUDE( "${${PROJ}_SOURCE_DIR}/cmake/UseDoxygen.cmake" )
    IF ( DOXYGEN_FOUND )
        ADD_DEPENDENCIES( doxygen latex_docs )
        ADD_DEPENDENCIES( doc latex_docs doxygen )
    ELSE()
        SET( USE_DOXYGEN 0 )
    ENDIF()
ENDIF()


# Create custom targets for build-test, check, and distclean
ADD_CUSTOM_TARGET( build-test )
ADD_CUSTOM_TARGET( build-examples )
ADD_CUSTOM_TARGET( check COMMAND  make test  )
ADD_DISTCLEAN( analysis null_timer tests liblbpm-wia.* cpu gpu cuda hip example common IO threadpool StackTrace )


# Check for CUDA
CHECK_ENABLE_FLAG( USE_CUDA 0 )
CHECK_ENABLE_FLAG( USE_HIP 0 )
NULL_USE( CMAKE_CUDA_FLAGS )
IF ( USE_CUDA )
    ADD_DEFINITIONS( -DUSE_CUDA )
    ENABLE_LANGUAGE( CUDA )
ELSEIF ( USE_HIP )
    ENABLE_LANGUAGE( HIP )
    ADD_DEFINITIONS( -DUSE_HIP )
ENDIF()



# Configure external packages
IF ( NOT ONLY_BUILD_DOCS )
    CONFIGURE_MPI()     # MPI must be before other libraries
    CONFIGURE_MIC()
    CONFIGURE_HDF5()
    CONFIGURE_NETCDF()
    CONFIGURE_SILO()
    CONFIGURE_LBPM()
    CONFIGURE_TIMER( 0 "${${PROJ}_INSTALL_DIR}/null_timer" FALSE )
    CONFIGURE_LINE_COVERAGE()
    # Set the external library link list
    SET( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${TIMER_LIBS} )
ENDIF()



# Macro to create 1,2,4 processor tests
MACRO( ADD_LBPM_TEST_1_2_4 EXENAME ${ARGN} )
    ADD_LBPM_TEST( ${EXENAME} ${ARGN} )
    ADD_LBPM_TEST_PARALLEL( ${EXENAME} 2 ${ARGN} )
    ADD_LBPM_TEST_PARALLEL( ${EXENAME} 4 ${ARGN} )
ENDMACRO()


# Add the src directories
IF ( NOT ONLY_BUILD_DOCS )
    BEGIN_PACKAGE_CONFIG( lbpm-wia-library )
    ADD_PACKAGE_SUBDIRECTORY( common )
    ADD_PACKAGE_SUBDIRECTORY( analysis )
    ADD_PACKAGE_SUBDIRECTORY( IO )
    ADD_PACKAGE_SUBDIRECTORY( threadpool )
    ADD_PACKAGE_SUBDIRECTORY( StackTrace )
    ADD_PACKAGE_SUBDIRECTORY( models )
    IF ( USE_CUDA )
        ADD_PACKAGE_SUBDIRECTORY( cuda )
    ELSEIF ( USE_HIP )
        ADD_PACKAGE_SUBDIRECTORY( hip )
    ELSE()
        ADD_PACKAGE_SUBDIRECTORY( cpu )
    ENDIF()
    INSTALL_LBPM_TARGET( lbpm-wia-library  )
    ADD_SUBDIRECTORY( tests )
    ADD_SUBDIRECTORY( example )
    #ADD_SUBDIRECTORY( workflows )
    INSTALL_PROJ_LIB()
    CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/ValgrindSuppresionFile ${CMAKE_CURRENT_BINARY_DIR}/test/ValgrindSuppresionFile COPYONLY )
ENDIF()

