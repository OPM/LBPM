#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "threadpool/atomic_helpers.h"
#include "common/Utilities.h"
#include "common/UnitTest.h"

#define perr std::cerr
#define pout std::cout
#define printp printf


#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    // Using windows
    #define USE_WINDOWS
    #define NOMINMAX
    #include <stdlib.h>
    #include <windows.h>
    #include <process.h>
#elif defined(__APPLE__)
    // Using MAC
    #define USE_MAC
    #include <unistd.h>
    #include <mach/mach_init.h>
    #include <mach/thread_policy.h>
#elif defined(__linux) || defined(__unix) || defined(__posix)
    // Using Linux
    #define USE_LINUX
    #include <pthread.h>
    #include <unistd.h>
#else
    #error Unknown OS
#endif


#ifdef USE_WINDOWS
    #include <windows.h>
    #define TIME_TYPE LARGE_INTEGER
    #define get_time(x) QueryPerformanceCounter(x)
    #define get_diff(start,end,f) (((double)(end.QuadPart-start.QuadPart))/((double)f.QuadPart))
    #define get_frequency(f) QueryPerformanceFrequency(f)
    #define sleep(x) Sleep(x*1000)
#elif defined(USE_LINUX) || defined(USE_MAC)
    #include <sys/time.h>
    #define TIME_TYPE timeval
    #define get_time(x) gettimeofday(x,NULL);
    #define get_diff(start,end,f) (((double)end.tv_sec-start.tv_sec)+1e-6*((double)end.tv_usec-start.tv_usec))
    #define get_frequency(f) (*f=timeval())
#else
    #error Unknown OS
#endif


// Function to increment/decrement a counter N times
struct counter_data {
    AtomicOperations::counter_t *counter;
    int N;
};
void modify_counter( counter_data *data ) {
    int N = data->N;
    AtomicOperations::counter_t &counter = *(data->counter);
    if ( N > 0 ) {
        for (int i=0; i<N; i++)
            counter.increment();
    } else if ( N < 0 ) {
        for (int i=0; i<-N; i++)
            counter.decrement();
    }
}


// Define the thread handle type
#ifdef USE_WINDOWS
    typedef HANDLE thread_handle;
#elif defined(USE_LINUX) || defined(USE_MAC)
    typedef pthread_t* thread_handle;
#else
    #error Unknown OS
#endif

// Create a thread
#ifdef USE_WINDOWS
    static thread_handle create_thread( void (*routine)(void*), void* data ) {
        return (HANDLE)_beginthread( routine, 0, data);
    }
#elif defined(USE_LINUX) || defined(USE_MAC)
    static thread_handle create_thread( void (*routine)(void*), void* data ) {
        pthread_t *id = new pthread_t;
        pthread_create( id, NULL, (void*(*)(void*)) routine, data );
        return id;
    }
#else
    #error Unknown OS
#endif

// Destroy a thread
#ifdef USE_WINDOWS
    static void destroy_thread( thread_handle id ) {
        WaitForMultipleObjects( 1, &id, 1, 10000 );
    }
#elif defined(USE_LINUX) || defined(USE_MAC)
    static void destroy_thread( thread_handle id ) {
        pthread_join(*id,NULL);
        delete id;
    }
#else
    #error Unknown OS
#endif


/******************************************************************
* The main program                                                *
******************************************************************/
#ifdef USE_WINDOWS
    int __cdecl main(int, char **) {
#elif defined(USE_LINUX) || defined(USE_MAC)
    int main(int, char*[]) {
#else
    #error Unknown OS
#endif
    UnitTest ut;

    int N_threads = 64;     // Number of threads
    int N_count = 1000000;  // Number of work items

    TIME_TYPE start, end, f;
    get_frequency(&f);

    // Ensure we are using all processors
    #ifdef __USE_GNU
        int N_procs = sysconf( _SC_NPROCESSORS_ONLN );
        cpu_set_t mask;
        CPU_ZERO(&mask);
        for (int i=0; i<N_procs; i++)
            CPU_SET(i,&mask);
        sched_setaffinity(getpid(), sizeof(cpu_set_t), &mask );
    #endif

    // Create the counter we want to test
    AtomicOperations::counter_t count;
    counter_data data;
    data.counter = &count;
    data.N = 0;
    if ( count.increment() == 1 )
        ut.passes("increment count");
    else
        ut.failure("increment count");
    if ( count.decrement() == 0 )
        ut.passes("decrement count");
    else
        ut.failure("decrement count");
    count.setCount(3);
    if ( count.getCount() == 3 )
        ut.passes("set count");
    else
        ut.failure("set count");
    count.setCount(0);

    // Increment the counter in serial
    data.N = N_count;
    get_time(&start);
    modify_counter( &data );
    get_time(&end);
    double time_inc_serial = get_diff(start,end,f)/N_count;
    int val = count.getCount();
    if ( val != N_count ) {
        char tmp[100];
        sprintf(tmp,"Count of %i did not match expected count of %i",val,N_count);
        ut.failure(tmp);
    }
    printp("Time to increment (serial) = %0.1f ns\n",1e9*time_inc_serial);

    // Decrement the counter in serial
    data.N = -N_count;
    get_time(&start);
    modify_counter( &data );
    get_time(&end);
    double time_dec_serial = get_diff(start,end,f)/N_count;
    val = count.getCount();
    if ( val != 0 ) {
        char tmp[100];
        sprintf(tmp,"Count of %i did not match expected count of %i",val,0);
        ut.failure(tmp);
    }
    printp("Time to decrement (serial) = %0.1f ns\n",1e9*time_dec_serial);

    // Increment the counter in parallel
    data.N = N_count;
    std::vector<thread_handle> thread_ids(N_threads);
    get_time(&start);
    for (int i=0; i<N_threads; i++) {
        thread_ids[i] = create_thread( (void (*)(void*)) modify_counter, (void*) &data );
    }
    for (int i=0; i<N_threads; i++) {
        destroy_thread( thread_ids[i] );
    }
    get_time(&end);
    double time_inc_parallel = get_diff(start,end,f)/(N_count*N_threads);
    val = count.getCount();
    if ( val != N_count*N_threads ) {
        char tmp[100];
        sprintf(tmp,"Count of %i did not match expected count of %i",val,N_count*N_threads);
        ut.failure(tmp);
    }
    printp("Time to increment (parallel) = %0.1f ns\n",1e9*time_inc_parallel);

    // Decrement the counter in parallel
    data.N = -N_count;
    get_time(&start);
    for (int i=0; i<N_threads; i++) {
        thread_ids[i] = create_thread( (void (*)(void*)) modify_counter, (void*) &data );
    }
    for (int i=0; i<N_threads; i++) {
        destroy_thread( thread_ids[i] );
    }
    get_time(&end);
    double time_dec_parallel = get_diff(start,end,f)/(N_count*N_threads);
    val = count.getCount();
    if ( val != 0 ) {
        char tmp[100];
        sprintf(tmp,"Count of %i did not match expected count of %i",val,0);
        ut.failure(tmp);
    }
    printp("Time to decrement (parallel) = %0.1f ns\n",1e9*time_dec_parallel);

    // Check the time to increment/decrement
    if ( time_inc_serial>100e-9 || time_dec_serial>100e-9 || time_inc_parallel>100e-9 || time_dec_serial>100e-9 ) {
        #if USE_GCOV
            ut.expected_failure("Time to increment/decrement count is too expensive");
        #else
            ut.failure("Time to increment/decrement count is too expensive");
        #endif
    } else {
        ut.passes("Time to increment/decrement passed");
    }

    // Finished
    ut.report();
    int N_errors = ut.NumFailGlobal();
    return N_errors;
}



