#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <sys/stat.h>
#include <math.h>
#include <stdexcept>
#include <vector>
#include <string>

#include "common/Utilities.h"
#include "common/UnitTest.h"

#ifdef USE_MPI
    #include "mpi.h"
#endif


// Function to return the call stack
std::vector<std::string> get_call_stack() 
{
    std::vector<std::string> stack = Utilities::getCallStack();
    // Trick compiler to skip inline for this function with fake recursion
    if ( stack.size() > 10000 ) { stack = get_call_stack(); } 
    return stack;
}


// The main function
int main(int argc, char *argv[]) 
{
    int rank = 0;
    #ifdef USE_MPI
        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    #endif
    UnitTest ut;
    Utilities::setAbortBehavior( true, true, true );

    // Test the memory usage
    double t0 = Utilities::time();
    size_t n_bytes1 = Utilities::getMemoryUsage();
    double time1 = Utilities::time() - t0;
    double *tmp = new double[0x100000];
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
        if ( n_bytes2>n_bytes1 )
            ut.passes("getMemoryUsage increases size");
        else
            ut.failure("getMemoryUsage increases size");
        if ( n_bytes1==n_bytes3 )
            ut.passes("getMemoryUsage decreases size properly");
        else
            ut.expected_failure("getMemoryUsage does not decrease size properly");
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

    // Finished
    ut.report();
    size_t N_errors = ut.NumFailGlobal();
    if ( N_errors==0 )
        printf("All tests passed\n");
    #ifdef USE_MPI
        MPI_Finalize();
    #endif
    return (int) N_errors;
}


