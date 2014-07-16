# Tool for building visit plugins
#
# The script will prompt the user to specify VISIT_ROOT_DIR if the prefix
# cannot be determined by the location of xml2cmake in the system path.
# Users can set the environment variable VISIT_BIN_PATH before running cmake
# (e.g. VISIT_BIN_PATH=/usr/local/bin instead of VISIT_ROOT_DIR)



# Find the xml2cmake executable and set VISIT_XML_CMAKE
SET( VISIT_XML_CMAKE )
FIND_PROGRAM( VISIT_XML_CMAKE
    NAMES xml2cmake
    PATHS 
        "${VISIT_ROOT_DIR}"
        "${VISIT_BIN_PATH}"
        "$ENV{VISIT_ROOT_DIR}"
        "$ENV{VISIT_BIN_PATH}"
        PATH_SUFFIXES bin bin64
        NO_DEFAULT_PATH
)
IF( NOT VISIT_XML_CMAKE )
    MESSAGE( FATAL_ERROR "xml2cmake not found in:\n"
        "${VISIT_ROOT_DIR}/bin\n"
        "${VISIT_BIN_PATH}\n"
        "$ENV{VISIT_ROOT_DIR}/bin\n"
        "$ENV{VISIT_BIN_PATH}\n"
    )
ELSE()
    MESSAGE( "VISIT_XML_CMAKE=${VISIT_XML_CMAKE}" )
ENDIF()


# Install plugin 
MACRO( VISIT_PLUGIN SRC_DIR TARGET )
    FILE( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${SRC_DIR}" )
    ADD_CUSTOM_TARGET( 
        ${SRC_DIR}
        ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/${SRC_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${SRC_DIR}"
        COMMAND ${VISIT_XML_CMAKE} ${TARGET}.xml
        COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} .
        COMMAND make
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${SRC_DIR}"
        SOURCES ${SRC_DIR}
    )
   
ENDMACRO()
