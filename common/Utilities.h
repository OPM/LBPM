#ifndef included_Utilities
#define included_Utilities

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>



/*!
 * Utilities is a Singleton class containing basic routines for error 
 * reporting, file manipulations, etc.  Included are a set of \ref Macros "macros" that are commonly used.
 */
namespace Utilities
{

    /*!
     * Aborts the run after printing an error message with file and
     * linenumber information.
     */
    void abort(const std::string &message, const std::string &filename, const int line);


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

    //! Factor a number into it's prime factors
    std::vector<int> factor(size_t number);

    //! Print AMP Banner
    void nullUse( void* );

} // namespace Utilities


#include "common/UtilityMacros.h"

#endif


