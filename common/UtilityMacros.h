// This file contains useful macros including ERROR, WARNING, INSIST, ASSERT, etc.
#ifndef included_UtilityMacros
#define included_UtilityMacros

#include "common/Utilities.h"

#include <iostream>
#include <sstream>
#include <stdexcept>


/*! \defgroup Macros Set of utility macro functions 
 *  \details  These functions are a list of C++ macros that are used 
 *     for common operations, including checking for errors.
 *  \addtogroup Macros
 *  @{
 */


/*! \def NULL_STATEMENT
 *  \brief    A null statement
 *  \details  A statement that does nothing, for insure++ make it something 
 * more complex than a simple C null statement to avoid a warning.
 */
#ifdef __INSURE__
    #define NULL_STATEMENT do{if(0) int nullstatement=0 }}while(0)
#else
    #define NULL_STATEMENT
#endif


/*! \def NULL_USE(variable)
 *  \brief    A null use of a variable
 *  \details  A null use of a variable, use to avoid GNU compiler warnings about unused variables.
 *  \param variable  Variable to pretend to use
 */
#define NULL_USE(variable) do {                         \
    if(0) {char *temp = (char *)&variable; temp++;}     \
}while(0)


/*! \def ERROR(MSG)
 *  \brief      Throw error
 *  \details    Throw an error exception from within any C++ source code.  The 
 *     macro argument may be any standard ostream expression.  The file and
 *     line number of the abort are also printed.
 *  \param MSG  Error message to print
 */
#define ERROR(MSG) do {                                 \
    Utilities::abort(MSG,__FILE__,__LINE__);            \
}while(0)


/*! \def WARNING(MSG)
 *  \brief   Print a warning
 *  \details Print a warning without exit.  Print file and line number of the warning.
 *  \param MSG  Warning message to print
 */
#define WARNING(MSG) do {                                           \
    std::stringstream tboxos;                                       \
    tboxos << MSG << std::ends;                                     \
    printf("WARNING: %s\n   Warning called in %s on line %i\n",tboxos.str().c_str(),__FILE__,__LINE__); \
}while(0)


/*! \def ASSERT(EXP)
 *  \brief Assert error
 *  \details Throw an error exception from within any C++ source code if the
 *     given expression is not true.  This is a parallel-friendly version
 *     of assert.
 *     The file and line number of the abort are printed along with the stack trace (if availible).
 *  \param EXP  Expression to evaluate
 */
#define ASSERT(EXP) do {                                            \
    if ( !(EXP) ) {                                                 \
        std::stringstream tboxos;                                   \
        tboxos << "Failed assertion: " << #EXP << std::ends;        \
        Utilities::abort(tboxos.str(), __FILE__, __LINE__);         \
    }                                                               \
}while(0)


/*! \def INSIST(EXP,MSG)
 *  \brief Insist error
 *  \details Throw an error exception from within any C++ source code if the
 *     given expression is not true.  This will also print the given message.
 *     This is a parallel-friendly version of assert.
 *     The file and line number of the abort are printed along with the stack trace (if availible).
 *  \param EXP  Expression to evaluate
 *  \param MSG  Debug message to print
 */
#define INSIST(EXP,MSG) do {                                        \
    if ( !(EXP) ) {                                                 \
        std::stringstream tboxos;                                   \
        tboxos << "Failed insist: " << #EXP << std::endl;           \
        tboxos << "Message: " << MSG << std::ends;                  \
        Utilities::abort(tboxos.str(), __FILE__, __LINE__);         \
    }                                                               \
}while(0)



/**
 * Macro for use when assertions are to be included
 * only when debugging.
 */
/*! \def CHECK_ASSERT(EXP)
 *  \brief Assert error (debug only)
 *  \details Throw an error exception from within any C++ source code if the
 *     given expression is not true.  This only runs if DEBUG_CHECK_ASSERTIONS
 *     is enabled.  If enabled, this is the same as a call to ASSERT.
 *  \param EXP  Expression to evaluate
 */
#ifdef DEBUG_CHECK_ASSERTIONS
    #define CHECK_ASSERT(EXP) ASSERT(EXP)
#else
    #define CHECK_ASSERT(EXP) 
#endif


/*! @} */


#endif
