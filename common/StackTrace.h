#ifndef included_AtomicStackTrace
#define included_AtomicStackTrace

#include <functional>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <thread>
#include <memory>
#include <set>


// Check for and include MPI
// clang-format off
#if defined(USE_MPI) || defined(USE_EXT_MPI)
    #include "mpi.h"
#elif defined(__has_include)
    #if __has_include("mpi.h")
        #include "mpi.h"
    #else
        typedef int MPI_Comm;
    #endif
#else
    typedef int MPI_Comm;
#endif
// clang-format on


namespace StackTrace {


struct stack_info {
    void *address;
    void *address2;
    std::string object;
    std::string function;
    std::string filename;
    int line;
    //! Default constructor
    stack_info() : address( nullptr ), address2( nullptr ), line( 0 ) {}
    //! Operator==
    bool operator==( const stack_info& rhs ) const;
    //! Operator!=
    bool operator!=( const stack_info& rhs ) const;
    //! Print the stack info
    std::string print() const;
    //! Compute the number of bytes needed to store the object
    size_t size() const;
    //! Pack the data to a byte array, returning a pointer to the end of the data
    char* pack( char* ptr ) const;
    //! Unpack the data from a byte array, returning a pointer to the end of the data
    const char* unpack( const char* ptr );
    //! Pack a vector of data to a memory block
    static std::vector<char> packArray( const std::vector<stack_info>& data );
    //! Unpack a vector of data from a memory block
    static std::vector<stack_info> unpackArray( const char* data );
};


struct multi_stack_info {
    int N;
    stack_info stack;
    std::vector<multi_stack_info> children;
    //! Default constructor
    multi_stack_info() : N( 0 ) {}
    //! Add the given stack to the multistack
    void add( size_t N, const stack_info *stack );
    //! Print the stack info
    std::vector<std::string> print( const std::string& prefix=std::string() ) const;
};


/*!
 * @brief  Get the current call stack
 * @details  This function returns the current call stack for the current thread
 * @return      Returns vector containing the stack
 */
std::vector<stack_info> getCallStack();


/*!
 * @brief  Get the current call stack for a thread
 * @details  This function returns the current call stack for the given thread
 * @param[in] id    The thread id of the stack we want to return
 * @return          Returns vector containing the stack
 */
std::vector<stack_info> getCallStack( std::thread::native_handle_type id );


/*!
 * @brief  Get the current call stack for all threads
 * @details  This function returns the current call stack for all threads
 *    in the current process.
 *    Note: This functionality may not be availible on all platforms
 * @return          Returns vector containing the stack
 */
multi_stack_info getAllCallStacks( );


/*!
 * @brief  Get the current call stack for all threads/processes
 * @details  This function returns the current call stack for all threads
 *    for all processes in the current process.  This function requires
 *    the user to call globalCallStackInitialize() before calling this
 *    routine, and globalCallStackFinalize() before exiting.
 *    Note: This functionality may not be availible on all platforms
 * @return          Returns vector containing the stack
 */
multi_stack_info getGlobalCallStacks( );


//! Function to return the current call stack for the current thread
std::vector<void *> backtrace();

//! Function to return the current call stack for the given thread
std::vector<void *> backtrace( std::thread::native_handle_type id );

//! Function to return the current call stack for all threads
std::vector<std::vector<void *>> backtraceAll();


//! Function to return the stack info for a given address
stack_info getStackInfo( void *address );


//! Function to return the stack info for a given address
std::vector<stack_info> getStackInfo( const std::vector<void *> &address );


//! Function to return the signal name
std::string signalName( int signal );


/*!
 * Return the symbols from the current executable (not availible for all platforms)
 * @return      Returns 0 if sucessful
 */
int getSymbols(
    std::vector<void *> &address, std::vector<char> &type, std::vector<std::string> &obj );


/*!
 * Return the name of the executable
 * @return      Returns the name of the executable (usually the full path)
 */
std::string getExecutable();


/*!
 * Return the search path for the symbols
 * @return      Returns the search path for the symbols
 */
std::string getSymPaths();


//!< Terminate type
enum class terminateType { signal, exception };

/*!
 * Set the error handlers
 * @param[in]   Function to terminate the program: abort(msg,type)
 */
void setErrorHandlers( std::function<void( std::string, terminateType )> abort );


/*!
 * Set the given signals to the handler
 * @param[in]   Function to terminate the program: abort(msg,type)
 */
void setSignals( const std::vector<int>& signals, void (*handler) (int) );


//! Clear a signal set by setSignals
void clearSignal( int signal );


//! Clear all signals set by setSignals
void clearSignals( );


//! Return a list of all signals that can be caught
std::vector<int> allSignalsToCatch( );

//! Return a default list of signals to catch
std::vector<int> defaultSignalsToCatch( );


//! Get a list of the active threads
std::set<std::thread::native_handle_type> activeThreads( );

//! Get a handle to this thread
std::thread::native_handle_type thisThread( );


//! Initialize globalCallStack functionallity
void globalCallStackInitialize( MPI_Comm comm );

//! Clean up globalCallStack functionallity
void globalCallStackFinalize( );


/*!
 * @brief  Call system command
 * @details  This function calls a system command, waits for the program
 *   to execute, captures and returns the output and exit code.
 * @param[in] cmd           Command to execute
 * @param[out] exit_code    Exit code returned from child process
 * @return                  Returns string containing the output
 */
std::string exec( const std::string& cmd, int& exit_code );


} // namespace StackTrace

#endif
