# This script compares the output of TEST against GOLD, 
#    ensuring that all lines within GOLD are int TEST.
# Note that TEST may have additional lines that are not checked
CMAKE_POLICY(SET CMP0007 OLD)

FILE(READ "${TEST}" output )
FILE(READ "${GOLD}" sol )

macro(LIST_REPLACE LIST INDEX NEWVALUE)
    list(INSERT ${LIST} ${INDEX} ${NEWVALUE})
    MATH(EXPR __INDEX "${INDEX} + 1")
    list(REMOVE_AT ${LIST} ${__INDEX})
endmacro(LIST_REPLACE)

# Convert file contents into a CMake list (where each element in the list is one line of the file)
STRING(REGEX REPLACE ";" "\\\\;" data "${output}")
STRING(REGEX REPLACE ";" "\\\\;" sol  "${sol}")
STRING(REGEX REPLACE "\n" ";" data "${data}")
STRING(REGEX REPLACE "\n" ";" sol "${sol}")
LIST( LENGTH data N_data )
LIST( LENGTH sol N_sol )
MATH( EXPR N_data "${N_data}-1" )
MATH( EXPR N_sol "${N_sol}-1" )
FOREACH( index RANGE ${N_data} )
    LIST(GET data ${index} tmp )
    STRING(REGEX REPLACE "(\n|\r)" "" tmp "${tmp}")
    STRING(STRIP "${tmp}" tmp )
    LIST_REPLACE( data ${index} "${tmp}")
ENDFOREACH()
FOREACH( index RANGE ${N_sol} )
    LIST( GET sol ${index} tmp )
    STRING(REGEX REPLACE "(\n|\r)" "" tmp "${tmp}")
    STRING(STRIP "${tmp}" tmp )
    LIST_REPLACE( sol ${index} "${tmp}")
ENDFOREACH()

# Check that each line of sol is present in data (and delete it)
FOREACH( tmp ${sol} )
    LIST(FIND data "${tmp}" result )
    IF ( ${result} EQUAL -1 )
        MESSAGE("Test output:\n${output}\n\n")
        MESSAGE(FATAL_ERROR "Did not find '${tmp}' in test output\n" )
    ELSE()
        LIST(REMOVE_AT data ${result} )
    ENDIF()
ENDFOREACH()

# Finished
MESSAGE( "All lines in ${GOLD} were found in ${TEST}")

