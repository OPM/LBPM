#ifndef included_StackTrace
#define included_StackTrace

#include <array>
#include <functional>
#include <iostream>
#include <set>
#include <thread>
#include <vector>

#include "StackTrace/string_view.h"


namespace StackTrace {

//! Class to contain stack trace info for a single thread/process
struct stack_info {
    uint32_t line;
    void *address;
    void *address2;
    std::array<char, 56> object;
    std::array<char, 48> objectPath;
    std::array<char, 64> filename;
    std::array<char, 64> filenamePath;
    std::array<char, 256> function;
    //! Default constructor
    stack_info();
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
    //! Print the stack info
    static void print( std::ostream &out, const std::vector<stack_info> &stack,
        const StackTrace::string_view &prefix = "" );
    //! Print the stack info
    void print2(
        char *txt, int widthAddress = 16, int widthObject = 20, int widthFunction = 32 ) const;
    //! Compute the number of bytes needed to store the object
    size_t size() const;
    //! Pack the data to a byte array, returning a pointer to the end of the data
    char *pack( char *ptr ) const;
    //! Unpack the data from a byte array, returning a pointer to the end of the data
    const char *unpack( const char *ptr );
};


//! Class to contain stack trace info for multiple threads/processes
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
    //! Is the stack empty
    bool empty() const { return N == 0; }
    //! Add the given stack to the multistack
    void add( size_t len, const stack_info *stack );
    //! Add the given stack to the multistack
    void add( const multi_stack_info &stack );
    //! Compute the number of bytes needed to store the object
    size_t size() const;
    //! Pack the data to a byte array, returning a pointer to the end of the data
    char *pack( char *ptr ) const;
    //! Unpack the data from a byte array, returning a pointer to the end of the data
    const char *unpack( const char *ptr );
    //! Print the stack info
    std::vector<std::string> print( const StackTrace::string_view &prefix = "" ) const;
    //! Print the stack info
    void print( std::ostream &out, const StackTrace::string_view &prefix = "" ) const;
    //! Print the stack info
    std::string printString( const StackTrace::string_view &prefix = "" ) const;

private:
    template<class FUN>
    void print2( int Np, char *prefix, int w[3], bool c, FUN &fun ) const;
    int getAddressWidth() const;
    int getObjectWidth() const;
    int getFunctionWidth() const;
};


//!< Terminate type
enum class terminateType : uint8_t { signal, exception, abort, MPI, unknown };
enum class printStackType : uint8_t { local = 1, threaded = 2, global = 3 };

//!< Class to contain exception info from abort
class abort_error : public std::exception
{
public:
    std::string message;       //!< Abort message
    std::string filename;      //!< File where abort was called
    terminateType type;        //!< What caused the termination
    printStackType stackType;  //!< Print the local stack, all threads, or global call stack
    uint8_t signal;            //!< Signal number
    int line;                  //!< Line number where abort was called
    size_t bytes;              //!< Memory in use during abort
    std::vector<void *> stack; //!< Local call stack for abort
public:
    virtual const char *what() const noexcept override;
    abort_error();
    virtual ~abort_error() {}

private:
    mutable std::string d_msg;
};


//!< Class to contain symbol information
struct symbols_struct {
    char type;
    void *address;
    std::array<char, 56> obj;
    std::array<char, 56> objPath;
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
const char *signalName( int signal );


/*!
 * Return the symbols from the current executable (not availible for all platforms)
 * @return      Returns the symbols loaded
 */
std::vector<symbols_struct> getSymbols();


//! Clear internal symbol data
void clearSymbols();


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


/*!
 * Set the given signals to the handler
 * @param[in] signals   Signals to handle
 * @param[in] handler   Function to terminate the program: abort(signal)
 */
void setSignals( const std::vector<int> &signals, void ( *handler )( int ) );


//! Clear a signal set by setSignals
void clearSignal( int signal );


//! Clear a signal set by setSignals
void clearSignals( const std::vector<int> &signals );


//! Clear all signals set by setSignals
void clearSignals();


//! Raise a signal
void raiseSignal( int signal );


//! Return a list of all signals that can be caught
std::vector<int> allSignalsToCatch();

//! Return a default list of signals to catch
std::vector<int> defaultSignalsToCatch();


//! Get a list of the active threads
std::vector<std::thread::native_handle_type> activeThreads();

//! Get a handle to this thread
std::thread::native_handle_type thisThread();


/*!
 * @brief  Call system command
 * @details  This function calls a system command, waits for the program
 *   to execute, captures and returns the output and exit code.
 * @param[in] cmd           Command to execute
 * @param[out] exit_code    Exit code returned from child process
 * @return                  Returns string containing the output
 */
std::string exec( const string_view &cmd, int &exit_code );


/*!
 * @brief  Create stack from string
 * @details  This function creates the call stack from the string generated by print
 * @param[in] str           Vector of strings containing call stack
 * @return                  Returns the call stack
 */
multi_stack_info generateFromString( const std::vector<std::string> &str );


/*!
 * @brief  Create stack from string
 * @details  This function creates the call stack from the string
 * @param[in] str           String containing call stack
 * @return                  Returns the call stack
 */
multi_stack_info generateFromString( const std::string &str );


} // namespace StackTrace


#endif
