#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "threadpool/thread_pool.h"
#include "common/UnitTest.h"
#include "common/Utilities.h"
#include "ProfilerApp.h"
#include "math.h"


#define perr std::cerr
#define pout std::cout
#define printp printf

#define MAX(x,y) ((x) > (y) ? (x) : (y))

#ifdef USE_WINDOWS
    #include <windows.h>
    #define TIME_TYPE LARGE_INTEGER
    #define get_time(x) QueryPerformanceCounter(x)
    #define get_diff(start,end,f) (((double)(end.QuadPart-start.QuadPart))/((double)f.QuadPart))
    #define get_frequency(f) QueryPerformanceFrequency(f)
    #define sleep(x) Sleep(x*1000)
    #define sleepMs(x) Sleep(x)
#elif defined(USE_LINUX) || defined(USE_MAC)
    #include <sys/time.h>
    #define TIME_TYPE timeval
    #define get_time(x) gettimeofday(x,NULL);
    #define get_diff(start,end,f) (((double)end.tv_sec-start.tv_sec)+1e-6*((double)end.tv_usec-start.tv_usec))
    #define get_frequency(f) (*f=timeval())
    #define sleepMs(x) usleep(1000*x)
#else
    #error Unknown OS
#endif

#ifdef USE_MPI
    #include "mpi.h"
#endif


// Wrapper function for mpi barrier
static inline void barrier() {
    #ifdef USE_MPI
        MPI_Barrier( MPI_COMM_WORLD );
    #endif
}

// Function to waste CPU cycles
void waste_cpu(int N) {
    if ( N > 10000 ) { PROFILE_START("waste_cpu",2); }
    double pi = 3.141592653589793;
    double x = 1.0;
    N = std::max(10,N);
    { for (int i=0; i<N; i++) x = sqrt(x*exp(pi/x)); } // style to limit gcov hits
    if ( fabs(x-2.926064057273157) > 1e-12 ) { abort(); }
    if ( N > 10000 ) { PROFILE_STOP("waste_cpu",2); }
}


// Function to sleep for N seconds then increment a global count
static volatile int global_sleep_count = 0;
void sleep_inc(int N) {
    PROFILE_START("sleep_inc");
    sleep(N);
    ++global_sleep_count;
    PROFILE_STOP("sleep_inc");
}
void sleep_inc2(double x) {
    sleepMs(static_cast<int>(round(x*1000)));
    ++global_sleep_count;
}
void sleep_msg( double x, std::string msg ) {
    PROFILE_START(msg);
    sleepMs(static_cast<int>(round(x*1000)));
    PROFILE_STOP(msg);
}
bool check_inc(int N) {
    return global_sleep_count==N;
}


// Function to return the processor for the given thread
void print_processor( ThreadPool* tpool ) 
{
    int rank = 0;
    #ifdef USE_MPI
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    #endif
    int thread = tpool->getThreadNumber();
    int processor = ThreadPool::getCurrentProcessor();
    char tmp[100];
    sprintf(tmp,"%i:  Thread,proc = %i,%i\n",rank,thread,processor);
    std::cout << tmp;
    sleepMs(100);
}


// Function to test how a member thread interacts with the thread pool
int test_member_thread( ThreadPool *tpool ) {
    int N_errors = 0;
    // Member threads are not allowed to wait for the pool to finish
    try {
        tpool->wait_pool_finished();
        N_errors++;
    } catch (...) {
    }
    // Member threads are not allowed to change the size of the pool
    try {
        tpool->wait_pool_finished();
        N_errors++;
    } catch (...) {
    }
    return N_errors;
}


// Function to test creating and locking a mutex
int test_mutex(bool recursive) {
    int N_errors = 0;
    Mutex lock(recursive);      // Create a lock
    Mutex lock2 = lock;         // Copy the lock
    // Test getting and releasing the lock
    lock.lock();
    lock.unlock();
    lock2.lock();
    lock2.unlock();
    bool own1 = lock.ownLock();
    bool own2 = lock2.ownLock();
    lock.lock();
    bool own3 = lock.ownLock();
    bool own4 = lock2.ownLock();
    lock.unlock();
    bool own5 = lock.ownLock();
    if ( own1 || own2 || !own3 || !own4 || own5 )
        return 1;
    if ( recursive ) {
        // Test the behavior of a recursive lock
        lock.lock();
        if ( !lock.tryLock() )
            return 1;
        lock.unlock();
        lock.lock();
        lock.unlock();
    } else {
        // Test the behavior of a non-recursive lock
        lock.lock();
        if ( lock.tryLock() )
            return 1;
        lock.unlock();
        try {
            lock.unlock();
            N_errors++;
        } catch (...) {
        }
        try {
            lock.lock();
            lock.lock();
            N_errors++;
        } catch (...) {
            lock.unlock();
        }
        try {
            lock.lock();
            lock2.lock();
            N_errors++;
            lock.unlock();
            lock2.unlock();
        } catch (...) {
            lock.unlock();
        }
    }
    return N_errors;
}


// Functions to test the templates
int myfun0() { return 0; }
int myfun1(int) { return 1; }
int myfun2(int,float) { return 2; }
int myfun3(int,float,double) { return 3; }
int myfun4(int,float,double,char) { return 4; }
int myfun5(int,float,double,char,std::string) { return 5; }
int myfun6(int,float,double,char,std::string,int) { return 6; }
int myfun7(int,float,double,char,std::string,int,int) { return 7; }


// Function to test instantiation of functions with different number of arguments
void vfunarg00() { }
void vfunarg01(int) { }
void vfunarg02(int,char) { }
void vfunarg03(int,char,double) { }
void vfunarg04(int,char,double,int) { }
void vfunarg05(int,char,double,int,char) { }
void vfunarg06(int,char,double,int,char,double) { }
void vfunarg07(int,char,double,int,char,double,int) { }
void vfunarg08(int,char,double,int,char,double,int,char) { }
void vfunarg09(int,char,double,int,char,double,int,char,double) { }
void vfunarg10(int,char,double,int,char,double,int,char,double,int) { }
void vfunarg11(int,char,double,int,char,double,int,char,double,int,char) { }
void vfunarg12(int,char,double,int,char,double,int,char,double,int,char,double) { }
void vfunarg13(int,char,double,int,char,double,int,char,double,int,char,double,int) { }
void vfunarg14(int,char,double,int,char,double,int,char,double,int,char,double,int,char) { }
void vfunarg15(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double) { }
void vfunarg16(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int) { }
void vfunarg17(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char) { }
void vfunarg18(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double) { }
void vfunarg19(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int) { }
void vfunarg20(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char) { }
void vfunarg21(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double) { }
void vfunarg22(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int) { }
void vfunarg23(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char) { }
void vfunarg24(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double) { }
int funarg00() { return 0; }
int funarg01(int) { return 1; }
int funarg02(int,char) { return 2; }
int funarg03(int,char,double) { return 3; }
int funarg04(int,char,double,int) { return 4; }
int funarg05(int,char,double,int,char) { return 5; }
int funarg06(int,char,double,int,char,double) { return 6; }
int funarg07(int,char,double,int,char,double,int) { return 7; }
int funarg08(int,char,double,int,char,double,int,char) { return 8; }
int funarg09(int,char,double,int,char,double,int,char,double) { return 9; }
int funarg10(int,char,double,int,char,double,int,char,double,int) { return 10; }
int funarg11(int,char,double,int,char,double,int,char,double,int,char) { return 11; }
int funarg12(int,char,double,int,char,double,int,char,double,int,char,double) { return 12; }
int funarg13(int,char,double,int,char,double,int,char,double,int,char,double,int) { return 13; }
int funarg14(int,char,double,int,char,double,int,char,double,int,char,double,int,char) { return 14; }
int funarg15(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double) { return 15; }
int funarg16(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int) { return 16; }
int funarg17(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char) { return 17; }
int funarg18(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double) { return 18; }
int funarg19(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int) { return 19; }
int funarg20(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char) { return 20; }
int funarg21(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double) { return 21; }
int funarg22(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int) { return 22; }
int funarg23(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char) { return 23; }
int funarg24(int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double,int,char,double) { return 24; }
int test_function_arguements( ThreadPool* tpool ) {
    int N_errors = 0;
    // Test some basic types of instantiations
    ThreadPool::thread_id_t id0 = TPOOL_ADD_WORK( tpool, myfun0, (NULL) );
    ThreadPool::thread_id_t id1 = TPOOL_ADD_WORK( tpool, myfun1, ( (int) 1 ) );
    ThreadPool::thread_id_t id2 = TPOOL_ADD_WORK( tpool, myfun2, ( (int) 1, (float) 2 ) );
    ThreadPool::thread_id_t id3 = TPOOL_ADD_WORK( tpool, myfun3, ( (int) 1, (float) 2, (double) 3 ) );
    ThreadPool::thread_id_t id4 = TPOOL_ADD_WORK( tpool, myfun4, ( (int) 1, (float) 2, (double) 3, (char) 4 ) );
    ThreadPool::thread_id_t id5 = TPOOL_ADD_WORK( tpool, myfun5, ( (int) 1, (float) 2, (double) 3, (char) 4, std::string("test") ) );
    ThreadPool::thread_id_t id52= TPOOL_ADD_WORK( tpool, myfun5, ( (int) 1, (float) 2, (double) 3, (char) 4, std::string("test") ), -1 );
    ThreadPool::thread_id_t id6 = TPOOL_ADD_WORK( tpool, myfun6, ( (int) 1, (float) 2, (double) 3, (char) 4, std::string("test"), (int) 1 ) );
    ThreadPool::thread_id_t id7 = TPOOL_ADD_WORK( tpool, myfun7, ( (int) 1, (float) 2, (double) 3, (char) 4, std::string("test"), (int) 1, (int) 1 ) );
    tpool->wait_pool_finished();
    if ( !tpool->isFinished(id0) ) { N_errors++; }
    if ( tpool->getFunctionRet<int>(id0) != 0 ) { N_errors++; }
    if ( tpool->getFunctionRet<int>(id1) != 1 ) { N_errors++; }
    if ( tpool->getFunctionRet<int>(id2) != 2 ) { N_errors++; }
    if ( tpool->getFunctionRet<int>(id3) != 3 ) { N_errors++; }
    if ( tpool->getFunctionRet<int>(id4) != 4 ) { N_errors++; }
    if ( tpool->getFunctionRet<int>(id5) != 5 ) { N_errors++; }
    if ( tpool->getFunctionRet<int>(id52)!= 5 ) { N_errors++; }
    if ( tpool->getFunctionRet<int>(id6) != 6 ) { N_errors++; }
    if ( tpool->getFunctionRet<int>(id7) != 7 ) { N_errors++; }
    // Test all the different numbers of arguments allowed
    TPOOL_ADD_WORK( tpool, vfunarg00, (NULL) );
    TPOOL_ADD_WORK( tpool, vfunarg01, ( 1) );
    TPOOL_ADD_WORK( tpool, vfunarg02, ( 1, 'a' ) );
    TPOOL_ADD_WORK( tpool, vfunarg03, ( 1, 'a', 3.0 ) );
    TPOOL_ADD_WORK( tpool, vfunarg04, ( 1, 'a', 3.0, 4 ) );
    TPOOL_ADD_WORK( tpool, vfunarg05, ( 1, 'a', 3.0, 4, 'e' ) );
    TPOOL_ADD_WORK( tpool, vfunarg06, ( 1, 'a', 3.0, 4, 'e', 6.0 ) );
    TPOOL_ADD_WORK( tpool, vfunarg07, ( 1, 'a', 3.0, 4, 'e', 6.0, 7 ) );
    TPOOL_ADD_WORK( tpool, vfunarg08, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h' ) );
    TPOOL_ADD_WORK( tpool, vfunarg09, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0 ) );
    TPOOL_ADD_WORK( tpool, vfunarg10, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10 ) );
    TPOOL_ADD_WORK( tpool, vfunarg11, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k' ) );
    TPOOL_ADD_WORK( tpool, vfunarg12, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0 ) );
    TPOOL_ADD_WORK( tpool, vfunarg13, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13 ) );
    TPOOL_ADD_WORK( tpool, vfunarg14, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n' ) );
    TPOOL_ADD_WORK( tpool, vfunarg15, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0 ) );
    TPOOL_ADD_WORK( tpool, vfunarg16, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16 ) );
    TPOOL_ADD_WORK( tpool, vfunarg17, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q' ) );
    TPOOL_ADD_WORK( tpool, vfunarg18, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0 ) );
    TPOOL_ADD_WORK( tpool, vfunarg19, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19 ) );
    TPOOL_ADD_WORK( tpool, vfunarg20, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19, 't' ) );
    TPOOL_ADD_WORK( tpool, vfunarg21, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19, 't', 21.0 ) );
    TPOOL_ADD_WORK( tpool, vfunarg22, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19, 't', 21.0, 22 ) );
    TPOOL_ADD_WORK( tpool, vfunarg23, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19, 't', 21.0, 22, 'w' ) );
    TPOOL_ADD_WORK( tpool, vfunarg24, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19, 't', 21.0, 22, 'w', 24.0 ) );
    std::vector<ThreadPool::thread_id_t> ids(25);
    ids[0]  = TPOOL_ADD_WORK( tpool, funarg00, (NULL) );
    ids[1]  = TPOOL_ADD_WORK( tpool, funarg01, ( 1 ) );
    ids[2]  = TPOOL_ADD_WORK( tpool, funarg02, ( 1, 'a' ) );
    ids[3]  = TPOOL_ADD_WORK( tpool, funarg03, ( 1, 'a', 3.0 ) );
    ids[4]  = TPOOL_ADD_WORK( tpool, funarg04, ( 1, 'a', 3.0, 4 ) );
    ids[5]  = TPOOL_ADD_WORK( tpool, funarg05, ( 1, 'a', 3.0, 4, 'e' ) );
    ids[6]  = TPOOL_ADD_WORK( tpool, funarg06, ( 1, 'a', 3.0, 4, 'e', 6.0 ) );
    ids[7]  = TPOOL_ADD_WORK( tpool, funarg07, ( 1, 'a', 3.0, 4, 'e', 6.0, 7) );
    ids[8]  = TPOOL_ADD_WORK( tpool, funarg08, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h' ) );
    ids[9]  = TPOOL_ADD_WORK( tpool, funarg09, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0 ) );
    ids[10] = TPOOL_ADD_WORK( tpool, funarg10, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10 ) );
    ids[11] = TPOOL_ADD_WORK( tpool, funarg11, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k' ) );
    ids[12] = TPOOL_ADD_WORK( tpool, funarg12, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0 ) );
    ids[13] = TPOOL_ADD_WORK( tpool, funarg13, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13 ) );
    ids[14] = TPOOL_ADD_WORK( tpool, funarg14, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'h' ) );
    ids[15] = TPOOL_ADD_WORK( tpool, funarg15, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'h', 15.0 ) );
    ids[16] = TPOOL_ADD_WORK( tpool, funarg16, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16 ) );
    ids[17] = TPOOL_ADD_WORK( tpool, funarg17, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q' ) );
    ids[18] = TPOOL_ADD_WORK( tpool, funarg18, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0 ) );
    ids[19] = TPOOL_ADD_WORK( tpool, funarg19, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19 ) );
    ids[20] = TPOOL_ADD_WORK( tpool, funarg20, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19, 't' ) );
    ids[21] = TPOOL_ADD_WORK( tpool, funarg21, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19, 't', 21.0 ) );
    ids[22] = TPOOL_ADD_WORK( tpool, funarg22, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19, 't', 21.0, 22 ) );
    ids[23] = TPOOL_ADD_WORK( tpool, funarg23, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19, 't', 21.0, 22, 'w' ) );
    ids[24] = TPOOL_ADD_WORK( tpool, funarg24, ( 1, 'a', 3.0, 4, 'e', 6.0, 7, 'h', 9.0, 10, 'k', 12.0, 13, 'n', 15.0, 16, 'q', 18.0, 19, 't', 21.0, 22, 'w', 24.0 ) );
    tpool->wait_all( ids );
    for (size_t i=0; i<ids.size(); i++) {
        if ( tpool->getFunctionRet<int>(ids[i]) != static_cast<int>(i) ) 
            N_errors++; 
    }
    return N_errors;
}


/******************************************************************
* Examples to derive a user work item                             *
******************************************************************/
class UserWorkItemVoid: public ThreadPool::WorkItem {
public:
    // User defined constructor (does not need to match any intrefaces)
    UserWorkItemVoid( int dummy )
    {
        // User initialized variables
        NULL_USE(dummy);
        // Set class variables
        ThreadPool::WorkItem::d_has_result = false;
        ThreadPool::WorkItem::d_state = 0;
    }
    // User defined run (can do anything)
    void run()
    {
        // Set the state (always do this first)
        ThreadPool::WorkItem::d_state = 1;
        // Perform the tasks
        printf("Hello work from UserWorkItem (void)");
        // Set the state (always do this last)
        ThreadPool::WorkItem::d_state = 2;
    }
    // User defined destructor
    virtual ~UserWorkItemVoid() 
    {
    }
};
class UserWorkItemInt: public ThreadPool::WorkItemRet<int> {
public:
    // User defined constructor (does not need to match any intrefaces)
    UserWorkItemInt( int dummy )
    {
        // User initialized variables
        NULL_USE(dummy);
        // Set class variables
        ThreadPool::WorkItem::d_has_result = true;
        ThreadPool::WorkItem::d_state = 0;
    }
    // User defined run (can do anything)
    void run()
    {
        // Set the state (always do this first)
        ThreadPool::WorkItem::d_state = 1;
        // Perform the tasks
        printf("Hello work from UserWorkItem (int)");
        // Store the results (it's type will match the template)
        ThreadPool::WorkItemRet<int>::d_result = 1;
        // Set the state (always do this last)
        ThreadPool::WorkItem::d_state = 2;
    }
    // User defined destructor
    virtual ~UserWorkItemInt() 
    {
    }
};


/******************************************************************
* test the time to run N tasks in parallel                        *
******************************************************************/
inline double run_parallel( ThreadPool *tpool, int N_tasks, int N_work )
{
    // Make sure the thread pool is empty
    tpool->wait_pool_finished();
    // Add the work
    TIME_TYPE start, end, f;
    get_frequency(&f);
    std::vector<ThreadPool::thread_id_t> ids;
    ids.reserve(N_tasks);
    get_time(&start);
    for (int i=0; i<N_tasks; i++)
        ids.push_back( TPOOL_ADD_WORK( tpool, waste_cpu, (N_work) ) );
    // Wait for the thread pool to finish
    tpool->wait_pool_finished();
    // Compute the time spent running the tasks
    get_time(&end);
    return get_diff(start,end,f);
}



/******************************************************************
* The main program                                                *
******************************************************************/
#ifdef USE_WINDOWS
    int __cdecl main(int argc, char **argv) {
#elif defined(USE_LINUX) || defined(USE_MAC)
    int main(int argc, char* argv[]) {
#else
    #error Unknown OS
#endif

    int N_threads = 4;      // Number of threads
    int N_work = 2000;      // Number of work items
    int N_it = 10;          // Number of cycles to run
    int N_problem = 4;      // Problem size
    PROFILE_ENABLE(3);
    PROFILE_ENABLE_TRACE();
    UnitTest ut;


    // Initialize MPI and set the error handlers
    int rank = 0;
    int size = 1;
    #ifdef USE_MPI
        int provided_thread_support=-1;
        MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided_thread_support);
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        Utilities::setErrorHandlers();
    #endif
    NULL_USE(size);
    

    // Disable OS specific warnings for all non-root ranks
    if ( rank > 0 )
        ThreadPool::set_OS_warnings(1);


    // Initialize the data
    std::vector<int> data1(N_work,0);
    std::vector<int> priority(N_work,0);
    for (int i=0; i<N_work; i++) {
        data1[i] = N_problem;
        priority[i] = i%128;
    }
    TIME_TYPE start, end, f, start2, end2;
    get_frequency(&f);


    // Print the size of the thread pool class
    printp("Size of ThreadPool = %i\n",(int)sizeof(ThreadPool));


    // Create and test a mutex
    barrier();
    printp("Testing mutex\n");
    int N_errors_mutex = test_mutex(false);
    N_errors_mutex    += test_mutex(true);
    if ( N_errors_mutex == 0 )
        ut.passes("test mutex");
    else
        ut.failure("Errors found testing mutex");


    // Get the number of processors availible
    barrier();
    int N_procs = 0;
    try {
        N_procs = ThreadPool::getNumberOfProcessors();
    } catch (...) {
    }
    if ( N_procs>0 )
        ut.passes("getNumberOfProcessors");
    else
        ut.failure("getNumberOfProcessors");
    printp("%i processors availible\n",N_procs);


    // Get the processor affinities for the process
    barrier();
    std::vector<int> cpus;
    try {
        cpus = ThreadPool::getProcessAffinity();
        printp("%i cpus for current process: ",(int)cpus.size());
        for (size_t i=0; i<cpus.size(); i++)
            printp("%i ",cpus[i]);
        printp("\n");
    } catch (...) {
    }
    if ( !cpus.empty() ) {
        ut.passes("getProcessAffinity");
    } else {
        #ifdef __APPLE__
            ut.expected_failure("getProcessAffinity");
        #else
            ut.failure("getProcessAffinity");
        #endif
    }


    // Test setting the process affinities
    barrier();
    bool pass = false;
    if ( !cpus.empty() && N_procs>0 ) {
        if ( cpus.size()==1 ) {
            cpus.resize(N_procs);
            for (int i=0; i<N_procs; i++)
                cpus.push_back( i );
            try {
                ThreadPool::setProcessAffinity( cpus );
            } catch (...) {
            }
            cpus = ThreadPool::getProcessAffinity();
            std::vector<int> cpus = ThreadPool::getProcessAffinity();
            printp("%i cpus for current process (updated): ",(int)cpus.size());
            for (size_t i=0; i<cpus.size(); i++)
                printp("%i ",cpus[i]);
            printp("\n");
            pass = cpus.size() > 1;
        } else {
            std::vector<int> cpus_orig = cpus;
            std::vector<int> cpus_tmp(1,cpus[0]);
            try {
                ThreadPool::setProcessAffinity( cpus_tmp );
            } catch (...) {
            }
            cpus = ThreadPool::getProcessAffinity();
            if ( cpus.size() == 1 )
                pass = true;
            try {
                ThreadPool::setProcessAffinity( cpus_orig );
            } catch (...) {
            }
            cpus = ThreadPool::getProcessAffinity();
            if ( cpus.size() != cpus_orig.size() )
                pass = false;
        }
    }
    if ( pass ) {
        ut.passes("setProcessAffinity");
    } else {
        #ifdef __APPLE__
            ut.expected_failure("setProcessAffinity");
        #else
            ut.failure("setProcessAffinity");
        #endif
    }
    int N_procs_used = std::min<int>(N_procs,N_threads);
    printp("%i processors used\n",N_procs_used);


    // Create the thread pool
    barrier();
    printp("Creating thread pool\n");
    ThreadPool tpool0;
    ThreadPool tpool;
    ThreadPool::thread_id_t id;
    id = TPOOL_ADD_WORK( &tpool, waste_cpu, ( data1[0] ) );
    if ( id==ThreadPool::thread_id_t() || !tpool.isValid(id) )
        ut.failure("Errors with id");
    tpool.setNumThreads(N_threads);
    if ( tpool.getNumThreads()==N_threads )
        ut.passes("Created thread pool");
    else
        ut.failure("Failed to create tpool with desired number of threads");

    // Test creating/destroying a thread pool using new
    barrier();
    pass = true;
    try {
        ThreadPool *tpool2 = new ThreadPool(MAX_NUM_THREADS-1);
        if ( tpool2->getNumThreads() != MAX_NUM_THREADS-1 )
            pass = false;
        if ( !ThreadPool::is_valid(tpool2) )
            pass = false;
        delete tpool2;
        // Check that tpool2 is invalid
        // Note: valgrind will report this as an invalid memory read, but we want to keep the test)
        if ( ThreadPool::is_valid(tpool2) )
            pass = false;
    } catch(...) {
        pass = false;
    }
    if ( tpool.getNumThreads()==N_threads )
        ut.passes("Created/destroyed thread pool with new");
    else
        ut.failure("Created/destroyed thread pool with new");


    // Test setting the thread affinities
    barrier();
    if ( cpus.size()>1 ) {
        sleepMs(50);
        // First make sure we can get the thread affinities
        std::vector<int> procs = ThreadPool::getThreadAffinity( );
        if ( procs == cpus ) {
            ut.passes("getThreadAffinity() matches procs");
        } else {
            char msg[100];
            sprintf(msg,"getThreadAffinity() does not match procs (%i,%i)",
                static_cast<int>(procs.size()), static_cast<int>(cpus.size()));
            ut.failure(msg);
        }
        pass = true;
        for (int i=0; i<N_threads; i++) {
            std::vector<int> procs_thread = tpool.getThreadAffinity( i );
            if ( procs_thread != procs ) {
                printp("%i: Initial thread affinity: ",rank);
                for (size_t i=0; i<procs_thread.size(); i++)
                    printp("%i ",procs_thread[i]);
                printp("\n");
                pass = false;
            }
        }
        if ( pass )
            ut.passes("getThreadAffinity(thread) matches procs");
        else
            ut.failure("getThreadAffinity(thread) does not match procs");
        // Try to set the thread affinities
        pass = true;
        if ( !procs.empty() ) {
            int N_procs_thread = std::max<int>((int)cpus.size()/N_threads,1);
            for (int i=0; i<N_threads; i++) {
                std::vector<int> procs_thread(N_procs_thread,-1);
                for (int j=0; j<N_procs_thread; j++)
                    procs_thread[j] = procs[(i*N_procs_thread+j)%procs.size()];
                tpool.setThreadAffinity( i, procs_thread );
                sleepMs(10);    // Give time for OS to update thread affinities
                std::vector<int> procs_thread2 = tpool.getThreadAffinity( i );
                if ( procs_thread2 != procs_thread ) {
                    printp("%i: Final thread affinity: ",rank);
                    for (size_t i=0; i<procs_thread.size(); i++)
                        printp("%i ",procs_thread[i]);
                    printp("\n");
                    pass = false;
                }
            }
        }
        if ( pass )
            ut.passes("setThreadAffinity passes");
        else
            ut.failure("setThreadAffinity failed to change affinity");
    }


    // Reset the thread affinities 
    barrier();
    tpool.setNumThreads(tpool.getNumThreads(),"none");
    //tpool.setNumThreads(tpool.getNumThreads(),"independent");
    for (int i=0; i<N_threads; i++) {
        std::vector<int> procs_thread = tpool.getThreadAffinity( i );
        printp("Thread affinity: ");
        for (size_t i=0; i<procs_thread.size(); i++)
            printp("%i ",procs_thread[i]);
        printp("\n");
    }

    // Print the current processors by thread id
    barrier();
    print_processor(&tpool);
    for (int i=0; i<N_threads; i++)
        TPOOL_ADD_WORK( &tpool, print_processor, ( &tpool ) );
    tpool.wait_pool_finished();

    // Run some basic tests
    barrier();
    get_time(&start);
    for (int n=0; n<N_it; n++) {
        for (int i=0; i<N_work; i++)
            waste_cpu(data1[i]);
    }
    get_time(&end);
    double time = get_diff(start,end,f);
    printp("Time for serial cycle = %0.0f us\n",1e6*time/N_it);
    printp("Time for serial item = %0.0f ns\n",1e9*time/(N_it*N_work));
    id = TPOOL_ADD_WORK( &tpool, waste_cpu, ( data1[0] ) );
    tpool.wait(id);
    std::vector<ThreadPool::thread_id_t> ids2;
    ids2.push_back( TPOOL_ADD_WORK( &tpool, waste_cpu, ( data1[0] ) ) );
    tpool.wait(ids2[0]);

    // Test calling functions with different number of arguments
    barrier();
    printp("Testing arguments:\n");
    int N_errors_args = test_function_arguements( &tpool );
    if ( N_errors_args == 0 )
        ut.passes("Calling function with default arguments");
    else
        ut.failure("Error calling function with default arguments");


    // Check that the threads can sleep in parallel (this does not depend on the number of processors)
    barrier();
    tpool.wait_pool_finished();
    get_time(&start);
    sleep_inc(1);
    get_time(&end);
    double sleep_serial = get_diff(start,end,f);
    ids2.clear();
    get_time(&start);
    for (int i=0; i<N_threads; i++)
        ids2.push_back( TPOOL_ADD_WORK( &tpool,sleep_inc, (1) ) );
    tpool.wait_all(N_procs_used,&ids2[0]);
    ids2.clear();
    get_time(&end);
    double sleep_parallel = get_diff(start,end,f);
    double sleep_speedup = N_procs_used*sleep_serial/sleep_parallel;
    printf("%i:  Speedup on %i sleeping threads: %0.3f\n",rank,N_procs_used,sleep_speedup);
    printf("%i:    ts = %0.3f, tp = %0.3f\n",rank,sleep_serial,sleep_parallel);
    if ( fabs(sleep_serial-1.0)<0.05 && fabs(sleep_parallel-1.0)<0.075 )
        ut.passes("Passed thread sleep");
    else 
        ut.failure("Failed thread sleep");


    // Check that the threads are actually working in parallel
    barrier();
    if ( N_procs_used>1 ) {
        #ifdef USE_MPI
            // Use a non-blocking serialization of the MPI processes 
            // if we do not have a sufficient number of processors
            bool serialize_mpi = N_procs < N_threads*size;
            int buf;
            MPI_Request request;
            MPI_Status status;
            if ( serialize_mpi && rank>0 ) {
                MPI_Irecv( &buf, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &request );
                int flag = false;
                while ( !flag ) {
                    MPI_Test( &request, &flag, &status );
                    sleep(1);
                }
            }
        #endif
        int N = 20000000;    // Enough work to keep the processor busy for ~ 1 s
        // Run in serial
        get_time(&start);
        waste_cpu(N);
        get_time(&end);
        double time_serial = get_diff(start,end,f);
        // Run in parallel
        double time_parallel2 = run_parallel( &tpool, N_procs_used, N/1000 );
        double time_parallel  = run_parallel( &tpool, N_procs_used, N );
        double speedup = N_procs_used*time_serial/time_parallel;
        printf("%i:  Speedup on %i procs: %0.3f\n",rank,N_procs_used,speedup);
        printf("%i:    ts = %0.3f, tp = %0.3f, tp2 = %0.3f\n",rank,time_serial,time_parallel,time_parallel2);
        if ( speedup > 1.4 ) {
            ut.passes("Passed speedup test");
        } else {
            #ifdef USE_GCOV
                ut.expected_failure("Times do not indicate tests are running in parallel (gcov)");
            #else
                ut.failure("Times do not indicate tests are running in parallel");
            #endif
        }
        #ifdef USE_MPI
            if ( serialize_mpi ) {
                if ( rank<size-1 )
                    MPI_Send( &N, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD );
                if ( rank==size-1 ) {
                    for (int i=0; i<size-1; i++)
                        MPI_Send( &N, 1, MPI_INT, i, 1, MPI_COMM_WORLD );
                } else {
                    MPI_Irecv( &buf, 1, MPI_INT, size-1, 1, MPI_COMM_WORLD, &request );
                    int flag = false;
                    MPI_Status status;
                    while ( !flag ) {
                        MPI_Test( &request, &flag, &status );
                        sleep(1);
                    }
                }
            }
        #endif
    } else {
        ut.expected_failure("Testing thread performance with less than 1 processor");
    }


    // Test adding a work item with a dependency
    barrier();
    {
        // Test that we sucessfully wait on the work items
        std::vector<ThreadPool::thread_id_t> ids;
        ids.reserve(5);
        global_sleep_count = 0;     // Reset the count before this test
        ThreadPool::thread_id_t id1 = TPOOL_ADD_WORK( &tpool, sleep_inc, ( 1 ) );
        ThreadPool::thread_id_t id2 = TPOOL_ADD_WORK( &tpool, sleep_inc, ( 2 ) );
        ThreadPool::WorkItem *wait1 = new WorkItemFull<bool,int>( check_inc, 1 );
        ThreadPool::WorkItem *wait2 = new WorkItemFull<bool,int>( check_inc, 2 );
        wait1->add_dependency(id1);
        wait2->add_dependency(id1);  wait2->add_dependency(id2);
        ids.clear();
        ids.push_back( tpool.add_work(wait1) );
        ids.push_back( tpool.add_work(wait2) );
        tpool.wait_all(ids.size(),&ids[0]);
        if ( !tpool.getFunctionRet<bool>(ids[0]) || !tpool.getFunctionRet<bool>(ids[1]) )
            ut.failure("Failed to wait on required dependency");
        else
            ut.passes("Dependencies");
        tpool.wait_pool_finished();
        // Check that we can handle more complex dependencies
        id1 = TPOOL_ADD_WORK( &tpool, sleep_inc2, ( 0.5 ) );
        for (int i=0; i<10; i++) {
            wait1 = new WorkItemFull<bool,int>( check_inc, 1 );
            wait1->add_dependency(id1);
            tpool.add_work(wait1);
        }
        tpool.wait_pool_finished();
        ids.clear();
        for (int i=0; i<5; i++)
            ids.push_back( TPOOL_ADD_WORK( &tpool, sleep_inc2, (0.5) ) );
        sleep_inc2(0.002);
        ThreadPool::WorkItem *work = new WorkItemFull<void,int>( waste_cpu, 100 );
        work->add_dependencies(ids);
        id = tpool.add_work(work,10);
        tpool.wait(id);
    }

    // Test the timing adding a single item
    barrier();
    for (int it=0; it<2; it++) {
        ThreadPool *tpool_ptr = NULL;
        if ( it==0 ) {
            printp("Testing timmings (adding a single item to empty tpool):\n");
            tpool_ptr = &tpool0;
        } else if ( it==1 ) {
            printp("Testing timmings (adding a single item):\n");
            tpool_ptr = &tpool;
        }
        std::vector<ThreadPool::thread_id_t> ids(N_work);
        double time_add = 0.0;
        double time_wait = 0.0;
        get_time(&start);
        for (int n=0; n<N_it; n++) {
            get_time(&start2);
            for (int i=0; i<N_work; i++)
                ids[i] = TPOOL_ADD_WORK( tpool_ptr, waste_cpu, ( data1[i] ), priority[i] );
            get_time(&end2);
            time_add += get_diff(start2,end2,f);
            get_time(&start2);
            tpool_ptr->wait_all(N_work,&ids[0]);
            //tpool_ptr->wait_pool_finished();
            get_time(&end2);
            time_wait += get_diff(start2,end2,f);
            if ( (n+1)%100 == 0 )
                printp("Cycle %i of %i finished\n",n+1,N_it);
        }
        get_time(&end);
        time = get_diff(start,end,f);
        printp("  time = %0.0f ms\n",1e3*time);
        printp("  time / cycle = %0.0f us\n",1e6*time/N_it);
        printp("  average time / item = %0.0f ns\n",1e9*time/(N_it*N_work));
        printp("     create and add = %0.0f ns\n",1e9*time_add/(N_it*N_work));
        printp("     wait = %0.0f us\n",1e9*time_wait/(N_it*N_work));
    }

    // Test the timing pre-creating the work items and adding multiple at a time
    barrier();
    for (int it=0; it<2; it++) {
        ThreadPool *tpool_ptr = NULL;
        if ( it==0 ) {
            printp("Testing timmings (adding a block of items to empty tpool):\n");
            tpool_ptr = &tpool0;
        } else if ( it==1 ) {
            printp("Testing timmings (adding a block of items):\n");
            tpool_ptr = &tpool;
        }
        double time_create_work = 0.0;
        double time_add_work = 0.0;
        double time_wait_work = 0.0;
        std::vector<ThreadPool::WorkItem*> work(N_work);
        get_time(&start);
        for (int n=0; n<N_it; n++) {
            get_time(&start2);
            for (int i=0; i<N_work; i++)
                work[i] = new WorkItemFull<void,int>( waste_cpu, data1[i] );
            get_time(&end2);
            time_create_work += get_diff(start2,end2,f);
            get_time(&start2);
            std::vector<ThreadPool::thread_id_t> ids = tpool_ptr->add_work( work, priority );
            get_time(&end2);
            time_add_work += get_diff(start2,end2,f);
            get_time(&start2);
            tpool_ptr->wait_all(ids);
            get_time(&end2);
            time_wait_work += get_diff(start2,end2,f);
            if ( (n+1)%100 == 0 )
                printp("Cycle %i of %i finished\n",n+1,N_it);
        }
        get_time(&end);
        time = get_diff(start,end,f);
        printp("  time = %0.0f ms\n",1e3*time);
        printp("  time / cycle = %0.0f us\n",1e6*time/N_it);
        printp("  average time / item = %0.0f ns\n",1e9*time/(N_it*N_work));
        printp("     create = %0.0f ns\n",1e9*time_create_work/(N_it*N_work));
        printp("     add = %0.0f ns\n",1e9*time_add_work/(N_it*N_work));
        printp("     wait = %0.0f ns\n",1e9*time_wait_work/(N_it*N_work));
    }

    // Run a dependency test that tests a simple case that should keep the thread pool busy
    // Note: Checking the results requires looking at the trace data
    tpool.wait_pool_finished();
    PROFILE_START("Dependency test");
    for (int i=0; i<10; i++) {
        char msg[3][100];
        sprintf(msg[0],"Item %i-%i",i,0);
        sprintf(msg[1],"Item %i-%i",i,1);
        sprintf(msg[2],"Item %i-%i",i,2);
        ThreadPool::WorkItem *work  = new WorkItemFull<void,double,std::string>(sleep_msg,0.5,msg[0]);
        ThreadPool::WorkItem *work1 = new WorkItemFull<void,double,std::string>(sleep_msg,0.1,msg[1]);
        ThreadPool::WorkItem *work2 = new WorkItemFull<void,double,std::string>(sleep_msg,0.1,msg[2]);
        ThreadPool::thread_id_t id = tpool.add_work(work);
        work1->add_dependency(id);
        work2->add_dependency(id);
        tpool.add_work(work1);
        tpool.add_work(work2);
    }
    tpool.wait_pool_finished();
    PROFILE_STOP("Dependency test");

    tpool.wait_pool_finished();
    barrier();
    ut.report();
    int N_errors = static_cast<int>(ut.NumFailGlobal());


    // Shudown MPI
    PROFILE_SAVE("test_thread_pool");
    pout << "Shutting down\n";
    barrier();
    #ifdef USE_MPI
        MPI_Finalize( );
        sleepMs(10);
    #endif
    return N_errors;
}



