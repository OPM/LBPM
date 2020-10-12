/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
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
#ifndef NULL_STATEMENT
#ifdef __INSURE__
#define NULL_STATEMENT            \
    do {                          \
        if ( 0 )                  \
            int nullstatement = 0 \
    } while ( 0 )
#else
#define NULL_STATEMENT
#endif
#endif


/*! \def NULL_USE(variable)
 *  \brief    A null use of a variable
 *  \details  A null use of a variable, use to avoid GNU compiler warnings about unused variables.
 *  \param variable  Variable to pretend to use
 */
#ifndef NULL_USE
#define NULL_USE( variable )                \
    do {                                    \
        if ( 0 ) {                          \
            auto temp = (char *) &variable; \
            temp++;                         \
        }                                   \
    } while ( 0 )
#endif


/*! \def ERROR(MSG)
 *  \brief      Throw error
 *  \details    Throw an error exception from within any C++ source code.  The 
 *     macro argument may be any standard ostream expression.  The file and
 *     line number of the abort are also printed.
 *  \param MSG  Error message to print
 */
#define ERROR(MSG)                                                \
    do {                                                          \
        ::Utilities::abort( MSG, __FILE__, __LINE__ );            \
    } while ( 0 )


/*! \def WARNING(MSG)
 *  \brief   Print a warning
 *  \details Print a warning without exit.  Print file and line number of the warning.
 *  \param MSG  Warning message to print
 */
#define WARNING(MSG)                                                    \
    do {                                                                \
        std::stringstream tboxos;                                       \
        tboxos << MSG << std::ends;                                     \
        printf("WARNING: %s\n   Warning called in %s on line %i\n",     \
            tboxos.str().c_str(),__FILE__,__LINE__);                    \
    }while(0)


/*! \def ASSERT(EXP)
 *  \brief Assert error
 *  \details Throw an error exception from within any C++ source code if the
 *     given expression is not true.  This is a parallel-friendly version
 *     of assert.
 *     The file and line number of the abort are printed along with the stack trace (if availible).
 *  \param EXP  Expression to evaluate
 */
#define ASSERT(EXP)                                                     \
    do {                                                                \
        if ( !(EXP) ) {                                                 \
            std::stringstream tboxos;                                   \
            tboxos << "Failed assertion: " << #EXP << std::ends;        \
            ::Utilities::abort(tboxos.str(), __FILE__, __LINE__);       \
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
        ::Utilities::abort(tboxos.str(), __FILE__, __LINE__);         \
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


/*! \def DISABLE_WARNINGS
 *  \brief Reenable warnings
 *  \details This will re-enable warnings after a call to DIASABLE_WARNINGS
 */
/*! \def ENABLE_WARNINGS
 *  \brief Supress all warnings
 *  \details This will start to supress all compile warnings.
 *      Be sure to follow with ENABLE_WARNINGS
 */
// clang-format off
#ifdef DISABLE_WARNINGS
    // Macros previously defined
#elif defined( USING_MSVC )
    #define DISABLE_WARNINGS __pragma( warning( push, 0 ) )
    #define ENABLE_WARNINGS __pragma( warning( pop ) )
#elif defined( USING_CLANG )
    #define DISABLE_WARNINGS                                                \
        _Pragma( "clang diagnostic push" ) _Pragma( "clang diagnostic ignored \"-Wall\"" ) \
        _Pragma( "clang diagnostic ignored \"-Wextra\"" )                   \
        _Pragma( "clang diagnostic ignored \"-Wunused-private-field\"" )    \
        _Pragma( "clang diagnostic ignored \"-Wmismatched-new-delete\"" )
    #define ENABLE_WARNINGS _Pragma( "clang diagnostic pop" )
#elif defined( USING_GCC )
    // Note: We cannot disable the -Wliteral-suffix message with this macro because the
    // pragma command cannot suppress warnings from the C++ preprocessor.  See gcc bug #53431.
    #define DISABLE_WARNINGS                                                \
        _Pragma( "GCC diagnostic push" ) _Pragma( "GCC diagnostic ignored \"-Wall\"" ) \
        _Pragma( "GCC diagnostic ignored \"-Wextra\"" )                     \
        _Pragma( "GCC diagnostic ignored \"-Wpragmas\"" )                     \
        _Pragma( "GCC diagnostic ignored \"-Wunused-local-typedefs\"" )     \
        _Pragma( "GCC diagnostic ignored \"-Woverloaded-virtual\"" )        \
        _Pragma( "GCC diagnostic ignored \"-Wunused-parameter\"" )          \
        _Pragma( "GCC diagnostic ignored \"-Warray-bounds\"" )              \
        _Pragma( "GCC diagnostic ignored \"-Wterminate\"" )
    #define ENABLE_WARNINGS _Pragma( "GCC diagnostic pop" )
#else
    #define DISABLE_WARNINGS
    #define ENABLE_WARNINGS
#endif
// clang-format on



/*! @} */


#endif
