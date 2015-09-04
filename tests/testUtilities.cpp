#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <math.h>
#include <stdexcept>
#include <string.h>
#include <stdint.h>

#include "common/Utilities.h"
#include "common/StackTrace.h"
#include "common/UnitTest.h"
#include "common/MPI_Helpers.h"


// Detect the OS (defines which tests we allow to fail)
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64) || defined(_MSC_VER)
    #define USE_WINDOWS
#elif defined(__APPLE__)
    #define USE_MAC
#elif defined(__linux) || defined(__unix) || defined(__posix)
    #define USE_LINUX
#else
    #error Unknown OS
#endif


// Function to return the call stack
std::vector<std::string> get_call_stack() 
{
    std::vector<StackTrace::stack_info> stack = StackTrace::getCallStack();
    std::vector<std::string> stack2(stack.size());
    for (size_t i=0; i<stack.size(); i++)
        stack2[i] = stack[i].print();
    // Trick compiler to skip inline for this function with fake recursion
    if ( stack.size() > 10000 ) { stack2 = get_call_stack(); } 
    return stack2;
}


// The main function
int main(int argc, char *argv[]) 
{
    int rank = 0;
    MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    UnitTest ut;
    Utilities::setAbortBehavior( true, true, true );

    // Limit the scope of variables
    { 
        // Test the memory usage
        double t0 = Utilities::time();
        size_t n_bytes1 = Utilities::getMemoryUsage();
        double time1 = Utilities::time() - t0;
        uint64_t *tmp = new uint64_t[0x100000];
        memset(tmp,0xAA,0x100000*sizeof(uint64_t));
        NULL_USE(tmp);
        t0 = Utilities::time();
        size_t n_bytes2 = Utilities::getMemoryUsage();
        double time2 = Utilities::time() - t0;
        delete [] tmp;
        t0 = Utilities::time();
        size_t n_bytes3 = Utilities::getMemoryUsage();
        double time3 = Utilities::time() - t0;
        std::cout << "Number of bytes used for a basic test: " << n_bytes1 << ", " << n_bytes2 << ", " << n_bytes3 << std::endl;
        std::cout << "   Time to query: " << time1*1e6 << " us, " << time2*1e6 << " us, " << time3*1e6 << " us" << std::endl;
        if ( n_bytes1==0 ) {
            ut.failure("getMemoryUsage returns 0");
        } else {
            ut.passes("getMemoryUsage returns non-zero");
            if ( n_bytes2>n_bytes1 ) {
                ut.passes("getMemoryUsage increases size");
            } else {
                #if defined(USE_MAC)
                    ut.expected_failure("getMemoryUsage does not increase size");
                #else
                    ut.failure("getMemoryUsage increases size");
                #endif
            }
            if ( n_bytes1==n_bytes3 ) {
                ut.passes("getMemoryUsage decreases size properly");
            } else {
                #if defined(USE_MAC) || defined(USE_WINDOWS)
                    ut.expected_failure("getMemoryUsage does not decrease size properly");
                #else
                    ut.failure("getMemoryUsage does not decrease size properly");
                #endif
            }
        }

        // Test getting the current call stack
        std::vector<std::string> call_stack = get_call_stack();
        if ( rank==0 ) {
            std::cout << "Call stack:" << std::endl;
            for (size_t i=0; i<call_stack.size(); i++)
                std::cout << "   " << call_stack[i] << std::endl;
        }
        if ( !call_stack.empty() ) {
            ut.passes("non empty call stack");
            if ( call_stack[0].find("get_call_stack()") != std::string::npos )
                ut.passes("call stack decoded function symbols");
            else
                ut.expected_failure("call stack was unable to decode function symbols");
        } else {
            ut.failure("non empty call stack");
        }

        // Test catching an error
        try {
            ERROR("Test");
            ut.failure("Failed to catch RAY_ERROR");
        } catch (...) {
            ut.passes("Caught RAY_ERROR");
        }
        try {
            throw std::logic_error("test");
            ut.failure("Failed to catch exception");
        } catch (...) {
            ut.passes("Caught exception");
        }
        
        // Test time/tick
        double time = Utilities::time();
        double res = Utilities::tick();
        if ( time==0 || res==0 )
            ut.failure("time/tick");
        else
            ut.passes("time/tick");

    }

    // Finished
    ut.report();
    size_t N_errors = ut.NumFailGlobal();
    if ( N_errors==0 )
        printf("All tests passed\n");
    MPI_Finalize();
    return (int) N_errors;
}


