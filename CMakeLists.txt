# Set some CMake properties    
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
CMAKE_POLICY(SET CMP0011 OLD)
CMAKE_POLICY(SET CMP0014 OLD)


MESSAGE("====================")
MESSAGE("Configuring LBPM-WIA")
MESSAGE("====================")


# Set the project name
PROJECT( LBPM-WIA )


# Check if we are only compiling docs
INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/macros.cmake )
INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/libraries.cmake )
CHECK_ENABLE_FLAG( ONLY_BUILD_DOCS 0 )


# Set testing paramaters
SET( DROP_METHOD "http" )
SET( DROP_SITE "oblivion.engr.colostate.edu" )
SET( DROP_LOCATION "/CDash/submit.php?project=LBPM-WIA" )
SET( TRIGGER_SITE "" )
SET( DROP_SITE_CDASH TRUE )
ENABLE_TESTING()
INCLUDE( CTest )


# Set some common paths
SET( LBPM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
SET( LBPM_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} )
IF( PREFIX )
    SET( LBPM_INSTALL_DIR ${PREFIX} )
ELSEIF( NOT LBPM_INSTALL_DIR )
    SET( LBPM_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR} )
ENDIF()
SET( CMAKE_MODULE_PATH ${LBPM_SOURCE_DIR} ${LBPM_SOURCE_DIR}/cmake )


# Create custom targets for build-test, check, and distclean
ADD_CUSTOM_TARGET( doc )
ADD_CUSTOM_TARGET( latex_docs )
ADD_CUSTOM_TARGET( build-test )
ADD_CUSTOM_TARGET( check COMMAND  make test  )
ADD_DISTCLEAN()


# Check the compile mode and compile flags
IF ( NOT ONLY_BUILD_DOCS )
    CONFIGURE_SYSTEM()
ENDIF()


# Add some directories to include
INCLUDE_DIRECTORIES( ${LBPM_INSTALL_DIR}/include )


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
    CONFIGURE_LBPM()
    CONFIGURE_LINE_COVERAGE()
ENDIF()


# Configure internal libraries (order matters)
SET( LBPM_LIBS  lbpm-wia )


# Add the src directories
IF ( NOT ONLY_BUILD_DOCS )
    BEGIN_PACKAGE_CONFIG( lbpm-wia )
    ADD_PACKAGE_SUBDIRECTORY( include )
    IF ( USE_CUDA )
        ADD_PACKAGE_SUBDIRECTORY( gpu )
    ELSE()
        ADD_PACKAGE_SUBDIRECTORY( cpu )
    ENDIF()
    INSTALL_LBPM_TARGET( lbpm-wia )
    ADD_SUBDIRECTORY( tests )
    INSTALL( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/example DESTINATION ${LBPM_INSTALL_DIR} )
ENDIF()


