# Set some CMake properties    
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
CMAKE_POLICY(SET CMP0014 OLD)
IF ( ${CMAKE_MAJOR_VERSION} EQUAL 3 )
    CMAKE_POLICY(SET CMP0026 OLD)
ENDIF()


MESSAGE("====================")
MESSAGE("Configuring LBPM-WIA")
MESSAGE("====================")


# Prevent users from building in place
IF ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}" )
    MESSAGE( FATAL_ERROR "Building code in place is a bad idea" )
ENDIF()


# Set the project name
PROJECT( LBPM-WIA )
IF ( NOT CXX_STD )
    SET( CXX_STD 11 )
ENDIF()


# Set some common paths
SET( LBPM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
SET( LBPM_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} )
IF( PREFIX )
    SET( LBPM_INSTALL_DIR ${PREFIX} )
ELSEIF( NOT LBPM_INSTALL_DIR )
    SET( LBPM_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR} )
ENDIF()
SET( CMAKE_MODULE_PATH ${LBPM_SOURCE_DIR} ${LBPM_SOURCE_DIR}/cmake )


# Include macros
SET( PROJ LBPM )
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


# Create custom targets for build-test, check, and distclean
ADD_CUSTOM_TARGET( doc )
ADD_CUSTOM_TARGET( latex_docs )
ADD_CUSTOM_TARGET( build-test )
ADD_CUSTOM_TARGET( build-examples )
ADD_CUSTOM_TARGET( check COMMAND  make test  )
ADD_DISTCLEAN( analysis null_timer tests liblbpm-wia.* cpu gpu example common visit IO threadpool )


# Check the compile mode and compile flags
IF ( NOT ONLY_BUILD_DOCS )
    CONFIGURE_SYSTEM()
ENDIF()


# Add some directories to include
INCLUDE_DIRECTORIES( "${LBPM_INSTALL_DIR}/include" )


# Set doxygen info
CHECK_ENABLE_FLAG( USE_DOXYGEN 1 )
CHECK_ENABLE_FLAG( USE_LATEX 1 )
FILE( MAKE_DIRECTORY "${LBPM_INSTALL_DIR}/doc" )
IF ( USE_DOXYGEN )
    SET( DOXYFILE_IN ${LBPM_SOURCE_DIR}/doxygen/Doxyfile.in)
    SET( DOXY_HEADER_FILE ${LBPM_SOURCE_DIR}/doxygen/html/header.html )
    SET( DOXY_FOOTER_FILE ${LBPM_SOURCE_DIR}/doxygen/html/footer.html )
    SET( DOXYFILE_OUTPUT_DIR ${LBPM_BUILD_DIR}/doc/ )
    SET( DOXYFILE_SRC_HTML_DIR ${LBPM_SOURCE_DIR}/doxygen/html )
    SET( DOXYFILE_SOURCE_DIR ${LBPM_SOURCE_DIR} )
    SET( REL_PACKAGE_HTML "" )
    SET( DOXYGEN_MACROS "" )
    MESSAGE("DOXYGEN_MACROS = ${DOXYGEN_MACROS}")
    INCLUDE( ${LBPM_SOURCE_DIR}/cmake/UseDoxygen.cmake )
    SET(DOXYFILE_LATEX "YES")
    IF ( DOXYGEN_FOUND )
        ADD_DEPENDENCIES( doxygen latex_docs )
        ADD_DEPENDENCIES( doc latex_docs doxygen )
    ELSE()
        SET( USE_DOXYGEN 0 )
    ENDIF()
ENDIF()


# Configure external packages
IF ( NOT ONLY_BUILD_DOCS )
    CONFIGURE_MPI()     # MPI must be before other libraries
    CONFIGURE_CUDA()
    CONFIGURE_MIC()
    CONFIGURE_NETCDF()
    CONFIGURE_LBPM()
    CONFIGURE_TIMER( 0 "${${PROJ}_INSTALL_DIR}/null_timer" )
    CONFIGURE_LINE_COVERAGE()
    INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SharedPtr.cmake" )
    CONFIGURE_SHARED_PTR( "${LBPM_INSTALL_DIR}/include" "std" )
    # Set the external library link list
    SET( EXTERNAL_LIBS ${EXTERNAL_LIBS} ${TIMER_LIBS} )
ENDIF()


# Configure internal libraries (order matters)
SET( LBPM_LIBS  lbpm-wia )


# Macro to create 1,2,4 processor tests
MACRO( ADD_LBPM_TEST_1_2_4 EXENAME ${ARGN} )
    ADD_LBPM_TEST( ${EXENAME} ${ARGN} )
    ADD_LBPM_TEST_PARALLEL( ${EXENAME} 2 ${ARGN} )
    ADD_LBPM_TEST_PARALLEL( ${EXENAME} 4 ${ARGN} )
ENDMACRO()


# Add the src directories
IF ( NOT ONLY_BUILD_DOCS )
    BEGIN_PACKAGE_CONFIG( lbpm-wia )
    ADD_PACKAGE_SUBDIRECTORY( common )
    ADD_PACKAGE_SUBDIRECTORY( analysis )
    ADD_PACKAGE_SUBDIRECTORY( IO )
    ADD_PACKAGE_SUBDIRECTORY( threadpool )
    IF ( USE_CUDA )
        ADD_PACKAGE_SUBDIRECTORY( gpu )
    ELSE()
        ADD_PACKAGE_SUBDIRECTORY( cpu )
    ENDIF()
    INSTALL_LBPM_TARGET( lbpm-wia )
    ADD_SUBDIRECTORY( tests )
    ADD_SUBDIRECTORY( threadpool/test )
    ADD_SUBDIRECTORY( example )
    ADD_SUBDIRECTORY( workflows )
ENDIF()


# Create the visit plugin
IF ( USE_VISIT )
    INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindVisIt.cmake )
    VISIT_PLUGIN( visit LBM )
ENDIF()


