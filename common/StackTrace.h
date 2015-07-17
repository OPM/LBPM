#ifndef included_StackTrace
#define included_StackTrace

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>



namespace StackTrace {


struct stack_info {
    void *address;
    void *address2;
    std::string object;
    std::string function;
    std::string filename;
    int line;
    //! Default constructor
    stack_info(): address(NULL), address2(NULL), line(0) {}
    //! Print the stack info
    std::string print() const;
};


//! Function to return the current call stack
std::vector<stack_info> getCallStack();


//! Function to return the stack info for a given address
stack_info getStackInfo( void* address );


/*!
 * Return the symbols from the current executable (not availible for all platforms)
 * @return      Returns 0 if sucessful
 */
int getSymbols( std::vector<void*>& address, std::vector<char>& type, std::vector<std::string>& obj );


} // namespace StackTrace


#endif

