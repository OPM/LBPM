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
    CONFIGURE_FILE( "${CMAKE_CURRENT_SOURCE_DIR}/${SRC_DIR}/${TARGET}.xml" "${CMAKE_CURRENT_BINARY_DIR}/${SRC_DIR}/${TARGET}.xml" )
    FILE( GLOB ConfigFiles RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/${SRC_DIR}" 
        "${CMAKE_CURRENT_SOURCE_DIR}/${SRC_DIR}/*.C" "${CMAKE_CURRENT_SOURCE_DIR}/${SRC_DIR}/*.h" )
    ADD_CUSTOM_TARGET(copy-${SRC_DIR})
    FOREACH( ConfigFile ${ConfigFiles} )
        ADD_CUSTOM_COMMAND(TARGET copy-${SRC_DIR} PRE_BUILD 
            COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${SRC_DIR}/${ConfigFile}" 
                "${CMAKE_CURRENT_BINARY_DIR}/${SRC_DIR}/${ConfigFile}"
        )
    ENDFOREACH()
    ADD_CUSTOM_TARGET( 
        ${SRC_DIR}
        COMMAND ${VISIT_XML_CMAKE} ${TARGET}.xml
        COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}  .
        COMMAND make
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${SRC_DIR}"
        SOURCES ${SRC_DIR}
        DEPENDS lbpm-wia copy-${SRC_DIR}
    )
ENDMACRO()


