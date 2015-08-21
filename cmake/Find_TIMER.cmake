# Configure the Profiler App for the current project
# The 'CONFIGURE_TIMER' macro searches for the timer library if included,
# or creates a dummy null timer if it is not used:
#    CONFIGURE_TIMER( DEFAULT_USE_TIMER NULL_TIMER_DIR )
# This function assumes that USE_TIMER is set to indicate if the timer should be used
# If USE_TIMER is set, TIMER_DIRECTORY specifies the install path for the timer
# If USE_TIMER is not set we will create a summy timer that does nothing.
# The input argument DEFAULT_USE_TIMER specifies if the timer library is included by default.
# The input argument NULL_TIMER_DIR specifies the location to install the dummy timer.  
#    If it is an empty string, the default install path "${CMAKE_CURRENT_BINARY_DIR}/null_timer" is used.
# This function will set the following variables and add the appropriate paths to the include list
#    TIMER_INCLUDE  - Path to the timer headers
#    TIMER_CXXFLAGS - C++ flags for the timer library
#    TIMER_LDFLAGS  - Linker flags to link the timer library
#    TIMER_LDLIBS   - Linker libraries to link the timer library
FUNCTION( CONFIGURE_TIMER DEFAULT_USE_TIMER NULL_TIMER_DIR )
    # Determine if we want to use the timer utility
    CHECK_ENABLE_FLAG( USE_TIMER ${DEFAULT_USE_TIMER} )
    SET( TIMER_INCLUDE )
    SET( TIMER_CXXFLAGS )
    SET( TIMER_LDFLAGS )
    SET( TIMER_LDLIBS )
    IF ( USE_TIMER )
        # Check if we specified the timer directory
        EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E remove -f "${NULL_TIMER_DIR}/ProfilerApp.h" )
        IF ( NOT TIMER_DIRECTORY AND TIMER_INSTALL_DIR )
            SET( TIMER_DIRECTORY ${TIMER_INSTALL_DIR} )
        ENDIF()
        IF ( TIMER_DIRECTORY )
            VERIFY_PATH( ${TIMER_DIRECTORY} )
            VERIFY_PATH( ${TIMER_DIRECTORY}/include )
            VERIFY_PATH( ${TIMER_DIRECTORY}/lib )
            FIND_LIBRARY( TIMER_LIBS  NAMES timerutility  PATHS ${TIMER_DIRECTORY}/lib  NO_DEFAULT_PATH )
            SET( TIMER_INCLUDE ${TIMER_DIRECTORY}/include )
            SET( TIMER_CXXFLAGS "-DUSE_TIMER -I${TIMER_DIRECTORY}/include" )
            SET( TIMER_LDFLAGS -L${TIMER_DIRECTORY}/lib )
            SET( TIMER_LDLIBS -ltimerutility )
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for TIMER is not yet supported.  Use -D TIMER_DIRECTORY=" )
        ENDIF()
        SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} "${TIMER_DIRECTORY}/lib" PARENT_SCOPE )
        INCLUDE_DIRECTORIES( "${TIMER_INCLUDE}" )
        ADD_DEFINITIONS( -DUSE_TIMER )
        MESSAGE( "Using timer utility" )
        MESSAGE( "  TIMER_LIBRARIES = ${TIMER_LIBS}" )
    ELSE()
        IF ( "${NULL_TIMER_DIR}" STREQUAL "" )
            SET( NULL_TIMER_DIR "${CMAKE_CURRENT_BINARY_DIR}/null_timer" )
        ENDIF()
        FILE(WRITE  "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_START(...)       do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_STOP(...)        do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_START2(...)      do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_STOP2(...)       do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_SCOPED(...)      do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_SYNCHRONIZE()    do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_SAVE(...)        do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_STORE_TRACE(X)   do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_ENABLE(...)      do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_DISABLE()        do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_ENABLE_TRACE()   do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_DISABLE_TRACE()  do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_ENABLE_MEMORY()  do {} while(0)\n" )
        FILE(APPEND "${NULL_TIMER_DIR}/ProfilerApp.h" "#define PROFILE_DISABLE_MEMORY() do {} while(0)\n" )
        SET( TIMER_INCLUDE  "${NULL_TIMER_DIR}" )
        INCLUDE_DIRECTORIES( "${TIMER_INCLUDE}" )
        MESSAGE( "Disabling timer utility" )
    ENDIF()
    SET( TIMER_INCLUDE  "${TIMER_INCLUDE}"  PARENT_SCOPE )
    SET( TIMER_CXXFLAGS "${TIMER_CXXFLAGS}" PARENT_SCOPE )
    SET( TIMER_LDFLAGS  "${TIMER_LDFLAGS}"  PARENT_SCOPE )
    SET( TIMER_LDLIBS   "${TIMER_LDLIBS}"   PARENT_SCOPE )
    SET( USE_TIMER      "${USE_TIMER}"      PARENT_SCOPE )
ENDFUNCTION()

# Check that a path is valid
FUNCTION( VERIFY_PATH PATH_NAME )
    IF ("${PATH_NAME}" STREQUAL "")
        MESSAGE( FATAL_ERROR "Path is not set: ${PATH_NAME}" )
    ENDIF()
    IF ( NOT EXISTS ${PATH_NAME} )
        MESSAGE( FATAL_ERROR "Path does not exist: ${PATH_NAME}" )
    ENDIF()
ENDFUNCTION()

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
    ENDIF ()
ENDMACRO()

