#ifndef included_Utilities
#define included_Utilities

#include <chrono>
#include <cstdarg>
#include <iostream>
#include <mutex>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <thread>
#include <vector>


namespace Utilities {


/*!
 * Aborts the run after printing an error message with file and
 * linenumber information.
 */
void abort( const std::string &message, const std::string &filename, const int line );


/*!
 * Set the behavior of abort
 * @param printMemory       Print the current memory usage (default is true)
 * @param printStack        Print the current call stack (default is true)
 * @param throwException    Throw an exception instead of MPI_Abort (default is false)
 */
void setAbortBehavior( bool printMemory, bool printStack, bool throwException );

//! Function to set the error handlers
void setErrorHandlers();


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


//! std::string version of sprintf
inline std::string stringf( const char *format, ... );


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


//! Factor a number into it's prime factors
std::vector<int> factor(size_t number);

//! Print AMP Banner
void nullUse( void* );

} // namespace Utilities


#include "common/UtilityMacros.h"


// stringf
inline std::string Utilities::stringf( const char *format, ... )
{
    va_list ap;
    va_start( ap, format );
    char tmp[4096];
    vsprintf( tmp, format, ap );
    va_end( ap );
    return std::string( tmp );
}


#endif
