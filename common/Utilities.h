#ifndef included_Utilities
#define included_Utilities

#include <cstdarg>
#include <vector>

#include "StackTrace/Utilities.h"


namespace Utilities {


// Functions inherited from StackTrace::Utilities
using StackTrace::Utilities::abort;
using StackTrace::Utilities::cause_segfault;
using StackTrace::Utilities::clearErrorHandlers;
using StackTrace::Utilities::exec;
using StackTrace::Utilities::getMemoryUsage;
using StackTrace::Utilities::getSystemMemory;
using StackTrace::Utilities::setAbortBehavior;
using StackTrace::Utilities::setErrorHandlers;
using StackTrace::Utilities::tick;
using StackTrace::Utilities::time;
using StackTrace::Utilities::sleep_ms;
using StackTrace::Utilities::sleep_s;


//! std::string version of sprintf
inline std::string stringf( const char *format, ... );


//! Factor a number into it's prime factors
std::vector<int> factor(size_t number);


//! Null use function
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
