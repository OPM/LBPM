#ifndef included_StackTrace_Utilities
#define included_StackTrace_Utilities

#include <stdexcept>
#include <string>
#include <thread>
#include <typeinfo>

#include "StackTrace/StackTrace.h"
#include "StackTrace/string_view.h"


namespace StackTrace {
namespace Utilities {


/*!
 * Aborts the run after printing an error message with file and
 * line number information.
 */
void abort( const std::string &message, const std::string &filename, const int line );


/*!
 * Set the behavior of abort
 * @param throwException    Throw an exception instead of MPI_Abort (default is false)
 * @param stackType         Type of stack to get (1: thread local stack, 2: all threads, 3: global)
 */
void setAbortBehavior( bool throwException, int stackType = 2 );


//! Function to terminate the application
void terminate( const StackTrace::abort_error &err );


//! Function to set the error handlers
void setErrorHandlers();


//! Function to clear the error handlers
void clearErrorHandlers();


/*!
 * Function to get the memory availible.
 * This function will return the total memory availible
 * Note: depending on the implimentation, this number may be rounded to
 * to a multiple of the page size.
 * If this function fails, it will return 0.
 */
size_t getSystemMemory();


/*!
 * Function to get the memory usage.
 * This function will return the total memory used by the application.
 * Note: depending on the implimentation, this number may be rounded to
 * to a multiple of the page size.
 * If this function fails, it will return 0.
 */
size_t getMemoryUsage();


//! Function to get an arbitrary point in time
double time();


//! Function to get the resolution of time
double tick();


/*!
 * Sleep for X ms
 * @param N         Time to sleep (ms)
 */
inline void sleep_ms( int N ) { std::this_thread::sleep_for( std::chrono::milliseconds( N ) ); }


/*!
 * Sleep for X s
 * @param N         Time to sleep (s)
 */
inline void sleep_s( int N ) { std::this_thread::sleep_for( std::chrono::seconds( N ) ); }


//! Cause a segfault
void cause_segfault();


/*!
 * @brief  Call system command
 * @details  This function calls a system command, waits for the program
 *   to execute, captures and returns the output and exit code.
 * @param[in] cmd           Command to execute
 * @param[out] exit_code    Exit code returned from child process
 * @return                  Returns string containing the output
 */
std::string exec( const StackTrace::string_view &cmd, int &exit_code );


//! Return the hopefully demangled name of the given type
std::string getTypeName( const std::type_info &id );


//! Return the hopefully demangled name of the given type
template<class TYPE>
inline std::string getTypeName()
{
    return getTypeName( typeid( TYPE ) );
}


} // namespace Utilities
} // namespace StackTrace


#endif
