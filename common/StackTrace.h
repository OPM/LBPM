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
#ifndef included_StackTrace
#define included_StackTrace

#include <functional>
#include <iostream>
#include <set>
#include <thread>
#include <vector>


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
    //! Reset the stack
    void clear();
    //! Operator==
    bool operator==( const stack_info &rhs ) const;
    //! Operator!=
    bool operator!=( const stack_info &rhs ) const;
    //! Get the minimum width to print the addresses
    int getAddressWidth() const;
    //! Print the stack info
    std::string print( int widthAddress = 16, int widthObject = 20, int widthFunction = 32 ) const;
    //! Compute the number of bytes needed to store the object
    size_t size() const;
    //! Pack the data to a byte array, returning a pointer to the end of the data
    char *pack( char *ptr ) const;
    //! Unpack the data from a byte array, returning a pointer to the end of the data
    const char *unpack( const char *ptr );
    //! Pack a vector of data to a memory block
    static std::vector<char> packArray( const std::vector<stack_info> &data );
    //! Unpack a vector of data from a memory block
    static std::vector<stack_info> unpackArray( const char *data );
};


struct multi_stack_info {
    int N;                                  // Number of threads/processes
    stack_info stack;                       // Current stack item
    std::vector<multi_stack_info> children; // Children
    //! Default constructor
    multi_stack_info() : N( 0 ) {}
    //! Construct from a simple call stack
    explicit multi_stack_info( const std::vector<stack_info> & );
    //! Copy constructor from a simple call stack
    multi_stack_info &operator=( const std::vector<stack_info> & );
    //! Reset the stack
    void clear();
    //! Add the given stack to the multistack
    void add( size_t len, const stack_info *stack );
    //! Print the stack info
    std::vector<std::string> print( const std::string &prefix = std::string() ) const;

private:
    void print2( const std::string &prefix, int w[3], std::vector<std::string> &text ) const;
    int getAddressWidth() const;
    int getObjectWidth() const;
    int getFunctionWidth() const;
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
multi_stack_info getAllCallStacks();


/*!
 * @brief  Get the current call stack for all threads/processes
 * @details  This function returns the current call stack for all threads
 *    for all processes in the current process.  This function requires
 *    the user to call globalCallStackInitialize() before calling this
 *    routine, and globalCallStackFinalize() before exiting.
 *    Note: This functionality may not be availible on all platforms
 * @return          Returns vector containing the stack
 */
multi_stack_info getGlobalCallStacks();


/*!
 * @brief  Clean up the stack trace
 * @details  This function modifies the stack trace to remove entries
 *    related to acquiring the stack trace in an attempt to make it
 *    more useful for display/users.
 * @param[in,out] stack     The stack trace to modify
 */
void cleanupStackTrace( multi_stack_info &stack );


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
int getSymbols( std::vector<void *> &address,
                std::vector<char> &type,
                std::vector<std::string> &obj );


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
 * @param[in] abort     Function to terminate the program: abort(msg,type)
 */
void setErrorHandlers( std::function<void( std::string, terminateType )> abort );


/*!
 * Set the given signals to the handler
 * @param[in] signals   Signals to handle
 * @param[in] handler   Function to terminate the program: abort(msg,type)
 */
void setSignals( const std::vector<int> &signals, void ( *handler )( int ) );


//! Clear a signal set by setSignals
void clearSignal( int signal );


//! Clear all signals set by setSignals
void clearSignals();


//! Return a list of all signals that can be caught
std::vector<int> allSignalsToCatch();

//! Return a default list of signals to catch
std::vector<int> defaultSignalsToCatch();


//! Get a list of the active threads
std::set<std::thread::native_handle_type> activeThreads();

//! Get a handle to this thread
std::thread::native_handle_type thisThread();


//! Initialize globalCallStack functionallity
void globalCallStackInitialize( MPI_Comm comm );

//! Clean up globalCallStack functionallity
void globalCallStackFinalize();


/*!
 * @brief  Call system command
 * @details  This function calls a system command, waits for the program
 *   to execute, captures and returns the output and exit code.
 * @param[in] cmd           Command to execute
 * @param[out] exit_code    Exit code returned from child process
 * @return                  Returns string containing the output
 */
std::string exec( const std::string &cmd, int &exit_code );


} // namespace StackTrace


#endif
