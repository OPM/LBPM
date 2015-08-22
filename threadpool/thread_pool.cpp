#define _CRT_NONSTDC_NO_DEPRECATE 
#include "thread_pool.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <typeinfo>
#include <stdexcept>
#include <climits>

#include "ProfilerApp.h"
#include "common/Utilities.h"
#define perr std::cerr
#define pout std::cout
#define printp printf

#define MONITOR_THREADPOOL_PERFORMANCE 0

#if 0
    #define PROFILE_THREAD_START(X)   PROFILE_START(X,3)
    #define PROFILE_THREAD_START2(X)  PROFILE_START2(X,3)
    #define PROFILE_THREAD_STOP(X)    PROFILE_STOP(X,3)
    #define PROFILE_THREAD_STOP2(X)   PROFILE_STOP2(X,3)
#else
    #define PROFILE_THREAD_START(X)   do {} while(0)
    #define PROFILE_THREAD_START2(X)  do {} while(0)
    #define PROFILE_THREAD_STOP(X)    do {} while(0)
    #define PROFILE_THREAD_STOP2(X)   do {} while(0)
#endif


// Include system dependent headers and define some functions
#ifdef __WORDSIZE
    #define ARCH_SIZE __WORDSIZE
#elif defined(_WIN64)
    #define ARCH_SIZE 64
#elif defined(_WIN32) // Note: WIN64 also defines WIN32
    #define ARCH_SIZE 32
#endif
#ifdef USE_WINDOWS
    #include <windows.h>
    #define get_time(x) QueryPerformanceCounter(x)
    #define get_frequency(f) QueryPerformanceFrequency(f)
    #define get_diff(start,end,f) \
        static_cast<double>(end.QuadPart-start.QuadPart)/static_cast<double>(f.QuadPart)
    #define TIME_TYPE LARGE_INTEGER
#elif defined(USE_LINUX)
    #include <sys/time.h>
    #include <errno.h>
    #define Sleep(x) usleep(1000*x)
    #define get_time(x) gettimeofday(x,NULL);
    #define get_frequency(f) (*f=timeval())
    #define get_diff(start,end,f) 1e-6*static_cast<double>( \
        0xF4240*(static_cast<int64_t>(end.tv_sec)-static_cast<int64_t>(start.tv_sec)) + \
                (static_cast<int64_t>(end.tv_usec)-static_cast<int64_t>(start.tv_usec)) )
    #define TIME_TYPE timeval
#elif defined(USE_MAC)
    #include <sys/time.h>
    #include <mach/mach.h>
    #include <errno.h>
    #define Sleep(x) usleep(1000*x)
    #define get_time(x) gettimeofday(x,NULL);
    #define get_frequency(f) (*f=timeval())
    #define get_diff(start,end,f) 1e-6*static_cast<double>( \
        0xF4240*(static_cast<int64_t>(end.tv_sec)-static_cast<int64_t>(start.tv_sec)) + \
                (static_cast<int64_t>(end.tv_usec)-static_cast<int64_t>(start.tv_usec)) )
    #define TIME_TYPE timeval
    #ifndef ARCH_SIZE
        #ifdef __LP64__
            #define ARCH_SIZE 64
        #else
            #define ARCH_SIZE 32
        #endif
    #endif
#else
    #error Unknown OS
#endif


// Check the ARCH_SIZE and set macros
// Note: ARCH_SIZE must match the number of bits in size_t
#if ARCH_SIZE == 64
    // 32-bit macros
#elif ARCH_SIZE == 32
    // 64-bit macros
#else
    #error Cannot identify 32 vs 64-bit
#endif


#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))


#if MONITOR_THREADPOOL_PERFORMANCE==1
    static TIME_TYPE frequency;            // Clock frequency (only used for windows)
    static double total_add_work_time[3] = {0,0,0};
#endif



// Helper functions
template <class T>  void quicksort(std::vector<T> &x);
static inline bool find_id(const std::vector<ThreadPool::thread_id_t> &x_in, const ThreadPool::thread_id_t &id );


// Function to generate a random size_t number (excluding 0 and ~0)
static size_t rand_size_t() 
{
    size_t key = 0;
    double tmp = 1;
    if ( sizeof(size_t)==4 ) {
        while ( tmp < 4e9 ) {
            key ^= rand()*0x9E3779B9;   // 2^32*0.5*(sqrt(5)-1)
            tmp *= RAND_MAX;
        }
    } else if ( sizeof(size_t)==8 ) {
        while ( tmp < 1.8e19 ) {
            key ^= rand()*0x9E3779B97F4A7C15;  // 2^64*0.5*(sqrt(5)-1)
            tmp *= RAND_MAX;
        }
    } else {
        throw std::logic_error("Unhandled case");
    }
    if ( key==0 || (~key)==0 )
        key = rand_size_t();
    return key;
}


/******************************************************************
* Run some basic compile-time checks                              *
******************************************************************/
#if MAX_NUM_THREADS%64 != 0
    // We use a bit array for d_active and d_cancel
    #error MAX_NUM_THREADS must be a multiple of 64
#endif
#if MAX_NUM_THREADS >= 65535
    // We store N_threads as a short int 
    #error MAX_NUM_THREADS must < 65535
#endif
#if MAX_QUEUED >= 65535
    // We store the indicies to the queue list as short ints
    #error MAX_QUEUED must < 65535
#endif


/******************************************************************
* Convert a string to binary                                      *
******************************************************************/
template<class T>
static inline std::string convert_binary(T x) {
    char buffer[65];
    T mask = ((size_t)1)<<(8*sizeof(T)-1);
    for (size_t i=0; i<8*sizeof(T); i++) {
        if ( ( x & mask ) == 0 )
            buffer[i] = '0';
        else
            buffer[i] = '1';
        mask >>= 1;
    }
    buffer[8*sizeof(T)] = 0;
    return std::string(buffer);
}


/******************************************************************
* Get/Set a bit                                                   *
******************************************************************/
static inline void set_bit( volatile ThreadPool::uint64* x, size_t index, bool val )
{
    ThreadPool::uint64 mask = 0x01;
    mask <<= index%64;
    if ( val )
        x[index/64] |= mask;
    else
        x[index/64] &= ~mask;
}
static inline bool get_bit( const volatile ThreadPool::uint64* x, size_t index )
{
    ThreadPool::uint64 mask = 0x01;
    mask <<= index%64;
    return (x[index/64]&mask)!=0;
}



/******************************************************************
* Some mutex helper functions                                     *
******************************************************************/
#if defined(USE_LINUX) || defined(USE_MAC)
    // Store a set of global attributes for the thread pool
    static pthread_mutexattr_t threadpool_global_attr;
    static int initialize_threadpool_global_attr() {
        pthread_mutexattr_init(&threadpool_global_attr);
        #ifdef __USE_UNIX98
            pthread_mutexattr_settype( &threadpool_global_attr, PTHREAD_MUTEX_ERRORCHECK );
        #endif
        return 1;
    }
    static int threadpool_global_attr_dummy = 0;
    static inline void throw_pthread_error( std::string msg, int value ) {
        std::string code;
        if ( value==0 ) {
            code = "SUCCESS";
        } else if ( value==EINVAL ) {
            code = "EINVAL";
        } else if ( value==EBUSY ) {
            code = "EBUSY";
        } else if ( value==EAGAIN ) {
            code = "EAGAIN";
        } else if ( value==EDEADLK ) {
            code = "EDEADLK";
        } else if ( value==EPERM ) {
            code = "EPERM";
        } else {
            char tmp[100];
            sprintf(tmp,"Unknown (%i)",value);
            code = std::string(tmp);
        }
        throw std::logic_error(msg+code);
    }
#endif
#ifdef USE_WINDOWS
    static inline void lock_mutex( CRITICAL_SECTION *lock ) {
        EnterCriticalSection(lock); 
    }
    static inline void unlock_mutex( CRITICAL_SECTION *lock ) {
        LeaveCriticalSection(lock);
    }
    static CRITICAL_SECTION* create_mutex( ) {
        CRITICAL_SECTION *lock = new CRITICAL_SECTION;
        if (!InitializeCriticalSectionAndSpinCount(lock,0x00000400) ) 
            throw std::exception();
        return lock;
    }
    static void destroy_mutex( CRITICAL_SECTION *lock ) {
        DeleteCriticalSection(lock);
        delete lock;
    }
#elif defined(USE_LINUX) || defined(USE_MAC)
    static inline void lock_mutex( pthread_mutex_t *lock ) {
        int retval = pthread_mutex_lock(lock);
        if ( retval != 0 )
            throw_pthread_error("Error locking mutex: ",retval);
    }
    static inline void unlock_mutex( pthread_mutex_t *lock ) {
        int retval = pthread_mutex_unlock(lock);
        if ( retval != 0 )
            throw_pthread_error("Error unlocking mutex: ",retval);
    }
    static pthread_mutex_t* create_mutex( ) {
        pthread_mutex_t* lock = NULL;
        #if defined(USE_LINUX) || defined(USE_MAC)
            if (threadpool_global_attr_dummy!=1)
                threadpool_global_attr_dummy = initialize_threadpool_global_attr();
        #endif
        // We are creating a new mutex
        lock = new pthread_mutex_t;
        int error = pthread_mutex_init(lock,&threadpool_global_attr);
        if ( error != 0 )
            throw_pthread_error("Error initializing mutex: ",error);
        return lock;
    }
    static void destroy_mutex( pthread_mutex_t* lock ) {
        pthread_mutex_destroy(lock);
        delete lock;
    }
#else
    #error Unknown OS
#endif


/******************************************************************
* Mutex class                                                     *
******************************************************************/
Mutex::Mutex()
{
    d_lock = create_mutex();
    d_recursive = false;
    d_count = new int;
    d_lock_count = new int;
    d_thread = new size_t;
    *d_count = 1;
    *d_lock_count = 0;
    *d_thread = 0;
}
Mutex::Mutex(bool recursive)
{
    d_lock = create_mutex();
    d_recursive = recursive;
    d_count = new int;
    d_lock_count = new int;
    d_thread = new size_t;
    *d_count = 1;
    *d_lock_count = 0;
    *d_thread = 0;
}
Mutex::Mutex(const Mutex& rhs)
{
    rhs.lock();
    d_lock = rhs.d_lock;
    d_count = rhs.d_count;
    d_recursive = rhs.d_recursive;
    d_lock_count = rhs.d_lock_count;
    d_thread = rhs.d_thread;
    ++(*d_count);
    rhs.unlock();
}
Mutex& Mutex::operator=(const Mutex& rhs)
{
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    rhs.lock();
    this->d_lock = rhs.d_lock;
    this->d_count = rhs.d_count;
    this->d_recursive = rhs.d_recursive;
    this->d_lock_count = rhs.d_lock_count;
    this->d_thread = rhs.d_thread;
    ++(*this->d_count);
    rhs.unlock();
    return *this;
}
Mutex::~Mutex()
{
    lock();
    bool destroy = (*d_count)==1;
    (*d_count)--;
    unlock();
    if ( destroy ) {
        delete d_count;
        delete d_lock_count;
        delete d_thread;
        destroy_mutex(d_lock);
    }
}
void Mutex::lock() const
{
    // Check if we already own the lock
    size_t id = ThreadPool::getThreadId();
    if ( *d_lock_count>0 && *d_thread==id ) {
        if ( !d_recursive )
            throw std::logic_error("Lock is already locked and non-recursive");
        // Increment the lock count and return
        ++(*d_lock_count);
        return;
    }
    // Acquire the lock
    lock_mutex(d_lock);
    if ( *d_lock_count != 0 )   // If we are getting the lock, the count must be 0
        throw std::logic_error("Internal error");
    *d_lock_count = 1;  // Change lock count after acquiring mutex
    *d_thread = id;
}
bool Mutex::tryLock() const
{
    // Check if we already own the lock
    size_t id = ThreadPool::getThreadId();
    if ( *d_lock_count>0 && *d_thread==id ) {
        if ( !d_recursive )
            return false;
        // Increment the lock count and return
        ++(*d_lock_count);
        return true;
    }
    // Try and acquire the lock
    #ifdef USE_WINDOWS
        bool success = TryEnterCriticalSection(d_lock)!=0;
    #elif defined(USE_LINUX) || defined(USE_MAC)
        bool success = pthread_mutex_trylock(const_cast<pthread_mutex_t*>(d_lock))==0;
    #else
        #error Unknown OS
    #endif
    if ( success ) {
        if ( *d_lock_count != 0 )   // If we are getting the lock, the count must be 0
            throw std::logic_error("Internal error");
        *d_lock_count = 1;  // Chage lock count after acquiring mutex
        *d_thread = id;
    }
    return success;
}
void Mutex::unlock() const
{
    // Check if we already own the lock
    size_t id = ThreadPool::getThreadId();
    if ( *d_lock_count <= 0 )
        throw std::logic_error("Trying to release a lock that has not been locked");
    if ( *d_thread != id )
        throw std::logic_error("Thread that does not own lock is attempting to release");
    // Release the lock
    --(*d_lock_count);  // Change lock count before releasing mutex
    if ( *d_lock_count == 0 ) {
        *d_thread = 0;
        unlock_mutex(d_lock);
    }
}
bool Mutex::ownLock() const
{
    size_t id = ThreadPool::getThreadId();
    if ( *d_lock_count>0 && *d_thread==id )
        return true;
    return false;
}


/******************************************************************
* Functions to deal with the signaling                            *
******************************************************************/
#ifdef USE_WINDOWS
    static inline bool SIGNAL_EVENT(HANDLE event) {
        SetEvent(event);
         return false;
    }
#elif defined(USE_LINUX) || defined(USE_MAC)
    static inline bool SIGNAL_EVENT(pthread_cond_t *event) {
        int retval = pthread_cond_signal(event);
        if ( retval == -1 ) {
            perr << "Error signaling event\n";
            return true;
        }
        return false;
    }
#else
    #error Not programmed
#endif


/******************************************************************
* Simple function to check if the parity is odd (true) or even    *
******************************************************************/
static inline bool is_odd8(size_t x) {  // This only works for 64-bit integers
    x ^= (x >> 1);
    x ^= (x >> 2);
    x ^= (x >> 4);
    x ^= (x >> 8);
    x ^= (x >> 16);
    x ^= (x >> 32);
    return (x & 0x01) > 0;
}
template<class int_type>
static inline int count_bits(int_type x) {
    int count = 0;
    for (size_t i=0; i<8*sizeof(int_type); i++) {
        if ( (x>>i)&0x01 )
            ++count;
    }
    return count;
}


/******************************************************************
* Set the bahvior of OS warnings                                  *
******************************************************************/
static int global_OS_behavior = 0;
void ThreadPool::set_OS_warnings( int behavior )
{
    ASSERT(behavior>=0&&behavior<=2);
    global_OS_behavior = behavior;
}
static void OS_warning( const std::string& message )
{
    if ( global_OS_behavior==0 ) {
        pout << "Warning: " << message << std::endl;
    } else if ( global_OS_behavior==2 ) {
        perr << "Error: " << message << std::endl;
    }
}


/******************************************************************
* Function to return the number of prcessors availible            *
******************************************************************/
int ThreadPool::getNumberOfProcessors()
{
    #if defined(USE_LINUX) || defined(USE_MAC)
        return sysconf( _SC_NPROCESSORS_ONLN );
    #elif defined(USE_WINDOWS)
        SYSTEM_INFO sysinfo;
        GetSystemInfo( &sysinfo );
        return static_cast<int>(sysinfo.dwNumberOfProcessors);
    #else
        #error Unknown OS
    #endif
}


/******************************************************************
* Function to return the processor number of the current thread   *
******************************************************************/
int ThreadPool::getCurrentProcessor()
{
    #if defined(USE_LINUX) 
        return sched_getcpu()+1;
    #elif defined(USE_MAC)
        OS_warning("MAC does not support getCurrentProcessor");
        return 0;
    #elif defined(USE_WINDOWS)
        return GetCurrentProcessorNumber()+1;
    #else
        #error Unknown OS
    #endif
}


/******************************************************************
* Function to get/set the affinity of the current process         *
******************************************************************/
std::vector<int> ThreadPool::getProcessAffinity()
{
    std::vector<int> procs;
    #ifdef USE_LINUX
        #ifdef _GNU_SOURCE
            cpu_set_t mask;
            int error = sched_getaffinity(getpid(), sizeof(cpu_set_t), &mask );
            if ( error!=0 )
                throw std::logic_error("Error getting process affinity");
            for (int i=0; i<(int)sizeof(cpu_set_t)*CHAR_BIT; i++) {
                if ( CPU_ISSET(i,&mask) )
                    procs.push_back(i);
            }
        #else
            #warning sched_getaffinity is not supported for this compiler/OS
            OS_warning("sched_getaffinity is not supported for this compiler/OS");
            procs.clear();
        #endif
    #elif defined(USE_MAC)
        // MAC does not support getting or setting the affinity
        OS_warning("MAC does not support getting the process affinity");
        procs.clear();
    #elif defined(USE_WINDOWS)
        HANDLE hProc = GetCurrentProcess();
        size_t procMask;
        size_t sysMask;
        PDWORD_PTR procMaskPtr = reinterpret_cast<PDWORD_PTR>(&procMask);
        PDWORD_PTR sysMaskPtr  = reinterpret_cast<PDWORD_PTR>(&sysMask);
        GetProcessAffinityMask(hProc,procMaskPtr,sysMaskPtr);
        for (int i=0; i<(int)sizeof(size_t)*CHAR_BIT; i++) {
            if ( (procMask&0x1) != 0 )
                procs.push_back(i);
            procMask >>= 1;
        }
    #else
        #error Unknown OS
    #endif
    return procs;
}
void ThreadPool::setProcessAffinity( std::vector<int> procs )
{
    #ifdef USE_LINUX
        #ifdef _GNU_SOURCE
            cpu_set_t mask;
            CPU_ZERO(&mask);
            for (size_t i=0; i<procs.size(); i++)
                CPU_SET(procs[i],&mask);
            int error = sched_setaffinity(getpid(), sizeof(cpu_set_t), &mask );
            if ( error!=0 )
                throw std::logic_error("Error setting process affinity");
        #else
            #warning sched_setaffinity is not supported for this compiler/OS
            OS_warning("sched_setaffinity is not supported for this compiler/OS");
            procs.clear();
        #endif
    #elif defined(USE_MAC)
        // MAC does not support getting or setting the affinity
        OS_warning("Warning: MAC does not support setting the process affinity");
        procs.clear();
    #elif defined(USE_WINDOWS)
        DWORD mask = 0;
        for (size_t i=0; i<procs.size(); i++)
            mask |= ((DWORD)1) << procs[i];
        HANDLE hProc = GetCurrentProcess();
        SetProcessAffinityMask( hProc, mask );
    #else
        #error Unknown OS
    #endif
}


/******************************************************************
* Function to get the thread affinities                           *
******************************************************************/
#ifdef USE_WINDOWS
    DWORD GetThreadAffinityMask(HANDLE thread)
    {
        DWORD mask = 1;
        DWORD old = 0;
        // try every CPU one by one until one works or none are left
        while(mask)
        {
            old = SetThreadAffinityMask(thread, mask);
            if(old)
            {   // this one worked
                SetThreadAffinityMask(thread, old); // restore original
                return old;
            }
            else
            {
                if(GetLastError() != ERROR_INVALID_PARAMETER)
                    return 0; // fatal error, might as well throw an exception
            }
            mask <<= 1;
        }

        return 0;
    }
#endif
std::vector<int> ThreadPool::getThreadAffinity()
{
    std::vector<int> procs;
    #ifdef USE_LINUX
        #ifdef _GNU_SOURCE
            cpu_set_t mask;
            int error = pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &mask );
            if ( error!=0 )
                throw std::logic_error("Error getting thread affinity");
            for (int i=0; i<(int)sizeof(cpu_set_t)*CHAR_BIT; i++) {
                if ( CPU_ISSET(i,&mask) )
                    procs.push_back(i);
            }
        #else
            #warning pthread_getaffinity_np is not supported
            OS_warning("pthread does not support pthread_getaffinity_np");
            procs.clear();
        #endif
    #elif defined(USE_MAC)
        // MAC does not support getting or setting the affinity
        OS_warning("MAC does not support getting the thread affinity");
        procs.clear();
    #elif defined(USE_WINDOWS)
        size_t procMask = GetThreadAffinityMask(GetCurrentThread());
        for (int i=0; i<(int)sizeof(size_t)*CHAR_BIT; i++) {
            if ( (procMask&0x1) != 0 )
                procs.push_back(i);
            procMask >>= 1;
        }
    #else
        #error Unknown OS
    #endif
    return procs;
}
std::vector<int> ThreadPool::getThreadAffinity( int thread ) const
{
    if ( thread >= getNumThreads() )
        std::logic_error("Invalid thread number");
    std::vector<int> procs;
    #ifdef USE_LINUX
        #ifdef _GNU_SOURCE
            cpu_set_t mask;
            int error = pthread_getaffinity_np(d_hThread[thread], sizeof(cpu_set_t), &mask );
            if ( error!=0 )
                throw std::logic_error("Error getting thread affinity");
            for (int i=0; i<(int)sizeof(cpu_set_t)*CHAR_BIT; i++) {
                if ( CPU_ISSET(i,&mask) )
                    procs.push_back(i);
            }
        #else
            #warning pthread_getaffinity_np is not supported
            OS_warning("pthread does not support pthread_getaffinity_np");
            procs.clear();
        #endif
    #elif defined(USE_MAC)
        // MAC does not support getting or setting the affinity
        OS_warning("MAC does not support getting the thread affinity");
        procs.clear();
    #elif defined(USE_WINDOWS)
        size_t procMask = GetThreadAffinityMask(d_hThread[thread]);
        for (int i=0; i<(int)sizeof(size_t)*CHAR_BIT; i++) {
            if ( (procMask&0x1) != 0 )
                procs.push_back(i);
            procMask >>= 1;
        }
    #else
        #error Unknown OS
    #endif
    return procs;
}


/******************************************************************
* Function to set the thread affinity                             *
******************************************************************/
void ThreadPool::setThreadAffinity( std::vector<int> procs )
{
    #ifdef USE_LINUX
        #ifdef _GNU_SOURCE
            cpu_set_t mask;
            CPU_ZERO(&mask);
            for (size_t i=0; i<procs.size(); i++)
                CPU_SET(procs[i],&mask);
            int error = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &mask );
            if ( error!=0 )
                throw std::logic_error("Error setting thread affinity");
        #else
            #warning pthread_getaffinity_np is not supported
            OS_warning("pthread does not support pthread_setaffinity_np");
            procs.clear();
        #endif
    #elif defined(USE_MAC)
        // MAC does not support getting or setting the affinity
        NULL_USE(procs);
        OS_warning("MAC does not support setting the thread affinity");
    #elif defined(USE_WINDOWS)
        DWORD mask = 0;
        for (size_t i=0; i<procs.size(); i++)
            mask |= ((DWORD)1) << procs[i];
        SetThreadAffinityMask( GetCurrentThread(), mask );
    #else
        #error Unknown OS
    #endif
}
void ThreadPool::setThreadAffinity( int thread, std::vector<int> procs ) const
{
    if ( thread >= getNumThreads() )
        std::logic_error("Invalid thread number");
    #ifdef USE_LINUX
        #ifdef __USE_GNU
            cpu_set_t mask;
            CPU_ZERO(&mask);
            for (size_t i=0; i<procs.size(); i++)
                CPU_SET(procs[i],&mask);
            int error = pthread_setaffinity_np(d_hThread[thread], sizeof(cpu_set_t), &mask );
            if ( error!=0 )
                throw std::logic_error("Error setting thread affinity");
        #else
            #warning pthread_getaffinity_np is not supported
            OS_warning("pthread does not support pthread_setaffinity_np");
            procs.clear();
        #endif
    #elif defined(USE_MAC)
        // MAC does not support getting or setting the affinity
        NULL_USE(procs);
        OS_warning("MAC does not support getting the process affinity");
    #elif defined(USE_WINDOWS)
        DWORD mask = 0;
        for (size_t i=0; i<procs.size(); i++)
            mask |= ((DWORD)1) << procs[i];
        SetThreadAffinityMask( d_hThread[thread], mask );
    #else
        #error Unknown OS
    #endif
}


/******************************************************************
* Function to perform some basic checks before we start           *
******************************************************************/
void ThreadPool::check_startup(size_t size0) 
{
    // Check the size of the class to make sure that we don't have any
    // byte alignment problems between a library implimentation and a calling pacakge
    size_t size1 = sizeof(ThreadPool);
    size_t size2 = ((size_t)&d_NULL_HEAD)-((size_t)this)+sizeof(size_t);
    size_t size3 = ((size_t)&d_NULL_TAIL)-((size_t)this)+sizeof(size_t);
    if ( size0!=size1 || size1<size2 || size1<size3 )
        throw std::logic_error("Internal data format problem");
    // Check the size of variables 
    if ( sizeof(ThreadPool::uint64)!=8 )
        throw std::logic_error("uint64 is not 64 bits");
    if ( sizeof(AtomicOperations::int32_atomic)!=4 )
        throw std::logic_error("AtomicOperations::int32_atomic is not 32 bits");
    if ( sizeof(AtomicOperations::int64_atomic)!=8 )
        throw std::logic_error("AtomicOperations::int32_atomic is not 64 bits");
    // Check getting/setting a bit
    uint64 x[2] = {0x0,0x7};
    set_bit(x,2,true);
    set_bit(x,66,false);
    if ( x[0]!=4 || x[1]!=3 || !get_bit(x,2) || get_bit(x,66) )
        throw std::logic_error("Getting/setting a bit failed");
    // Check the thread id
    bool pass = true;
    ThreadPool::thread_id_t id;
    if ( id.getPriority()!=-128 )
        pass = false;
    id.reset(3,564,NULL);
    if ( id.getPriority()!=3 || id.getLocalID()!=564 )
        pass = false;
    if ( count_bits(0x0)!=0 || count_bits(0x03)!=2 )
        pass = false;
    if ( count_bits(~((size_t)0)) != 8*sizeof(size_t) )
        pass = false;
    if ( sizeof(size_t)==8 ) {
        if ( is_odd8(0x0) || !is_odd8(0x02) || is_odd8(0x03) )
            pass = false;
        if ( is_odd8(~((size_t)0)) || !is_odd8(MAXID64))
            pass = false;
        for (size_t i=0; i<1024; i++) {
            if ( (count_bits(MAXID64-i)%2==1) != is_odd8(MAXID64-i) ) {
                printp("%i %i %s\n",count_bits(MAXID64-i),is_odd8(MAXID64-i)?1:0,
                    convert_binary<unsigned long long int>(MAXID64-i).c_str());
                pass = false;
            }
        }
    }
    initialize_id();
    advance_id(); advance_id(); 
    ThreadPool::thread_id_t id2;
    id2.reset(3,d_id_assign,NULL);
    if ( isValid(id) || !isValid(id2) )
        pass = false;
    if ( !pass ) {
        throw std::logic_error("Thread pool failed to initialize");
    }
}


/******************************************************************
* Function to initialize the thread pool                          *
******************************************************************/
void ThreadPool::initialize( const int N, const char* affinity, int N_procs, const int* procs ) 
{
    // Get the clock frequency
    #if MONITOR_THREADPOOL_PERFORMANCE==1
        get_frequency( &frequency );
    #endif
    // Initialize the header/tail
    d_NULL_HEAD = rand_size_t();
    d_NULL_TAIL = d_NULL_HEAD;
    for (int i=0; i<MAX_NUM_THREADS; i++)
        d_hThread[i] = 0;
    // Initialize the variables to NULL values
    d_id_assign = 0;
    d_signal_empty = false;
    d_signal_count = 0;
    d_N_threads = 0;
    d_num_active = 0;
    d_queue_size = 0;
    d_N_wait = 0;
    for (int i=0; i<MAX_NUM_THREADS; i++)
        d_ThreadId[i] = ~((size_t)0);
    memset((void*)d_active,0,MAX_NUM_THREADS/8);
    memset((void*)d_cancel,0,MAX_NUM_THREADS/8);
    for (int i=0; i<MAX_QUEUED; i++) {
        d_queue_ids[i].reset();
        d_queue_list[i].reset();
        d_queue_list[i].position = i;
        d_queue_list[i].prev = i-1;
        d_queue_list[i].next = i+1;
    }
    d_queue_head = -1;
    d_queue_free = 0;
    for (int i=0; i<MAX_WAIT; i++)
        d_wait[i] = NULL;
    d_wait_finished = 0;
    d_lock_queue = 0;
    for (int i=0; i<MAX_NUM_THREADS; i++)
        d_hThread[i] = 0;
    #if defined(USE_LINUX) || defined(USE_MAC)
        d_queue_not_empty = 0;
    #endif
    // Initialize the id
    initialize_id();
    // Create the mutex lock and signal variables
    d_lock_queue = create_mutex();
    #ifdef USE_WINDOWS
        d_wait_finished = CreateEvent(NULL,FALSE,FALSE,NULL);
    #elif defined(USE_LINUX) || defined(USE_MAC)
        d_queue_not_empty = new pthread_cond_t;
        d_wait_finished = new pthread_cond_t;
        int error = pthread_cond_init(d_queue_not_empty,NULL);
        if ( error == -1 )
            perr << "Error creating d_queue_not_empty\n";
        error = pthread_cond_init(d_wait_finished,NULL);
        if ( error == -1 )
            perr << "Error creating d_wait_finished\n";
    #else
        #error Not programmed
    #endif
    // Create the threads
    setNumThreads(N,affinity,N_procs,procs);
}



/******************************************************************
* This is the de-constructor                                      *
******************************************************************/
ThreadPool::~ThreadPool() {
    if ( !is_valid(this) )
        throw std::logic_error("Thread pool is not valid");
    // Destroy the threads
    setNumThreads(0);
    // Delete all remaining data
    destroy_mutex(d_lock_queue);
    #ifdef USE_WINDOWS
        CloseHandle(d_wait_finished);
    #elif defined(USE_LINUX) || defined(USE_MAC)
        pthread_cond_destroy(d_wait_finished);
        pthread_cond_destroy(d_queue_not_empty);
        delete d_queue_not_empty;     d_queue_not_empty=NULL;
        delete d_wait_finished;       d_wait_finished=NULL;
    #else
        #error Not programmed
    #endif
    d_N_threads = -1;
    d_NULL_HEAD = 0;
    d_NULL_TAIL = 0;
    // Print the performance metrics
    #if MONITOR_THREADPOOL_PERFORMANCE==1
        printp("ThreadPool Performance:\n");
        printp("add_work: %e %e %e\n",total_add_work_time[0],total_add_work_time[1],total_add_work_time[2]);
    #endif
}


/******************************************************************
* Check if the pointer points to a valid thread pool object       *
******************************************************************/
bool ThreadPool::is_valid( const ThreadPool* tpool )
{
    if ( tpool == NULL )
        return false;
    if ( tpool->d_N_threads<0 || tpool->d_N_threads>MAX_NUM_THREADS )
        return false;
    if ( tpool->d_NULL_HEAD==0 || tpool->d_NULL_HEAD!=tpool->d_NULL_TAIL )
        return false;
    return true;
}


/******************************************************************
* This function creates the threads in the thread pool            *
******************************************************************/
void ThreadPool::setNumThreads( int num_worker_threads, 
    const char* affinity2, int N_procs, const int* procs ) 
{
    // Check if we are a member thread
    if ( isMemberThread() )
        throw std::logic_error("Member threads are not allowed to change the number of threads in the pool");
    // Determing the number of threads we need to create or destroy
    if ( num_worker_threads > MAX_NUM_THREADS ) {
        printp("Warning: Maximum Number of Threads is %i\n",MAX_NUM_THREADS);
        printp("         Only that number will be created\n");
        num_worker_threads = MAX_NUM_THREADS;
    } else if ( num_worker_threads < 0 ) {
        printp("Error: cannot have a negitive number of threads\n");
        printp("       Setting the number of threads to 0\n");
        num_worker_threads = 0;
    } 
    int d_N_threads_diff = num_worker_threads-d_N_threads;
    if ( d_N_threads_diff > 0 ) {
        // Create new threads
        lock_mutex(d_lock_queue);
        // Check that no threads are in the process of being deleted
        for (int i=0; i<MAX_NUM_THREADS/64; i++) {
            if ( d_cancel[i] != 0 )
                throw std::logic_error("Threads are being created and destroyed at the same time");
        }
        // Create the thread attributes (linux only)
        #if defined(USE_LINUX) || defined(USE_MAC)
            pthread_attr_t attr;
            pthread_attr_init(&attr);
            //int ptmp;
            //pthread_attr_setstacksize(&attr,2097152);     // Default stack size is 8MB
            //pthread_attr_setschedpolicy(&attr,1);
            //pthread_attr_getschedpolicy(&attr,&ptmp);
            //pout << "getschedpolicy = " << ptmp << std::endl;
        #endif
        // Create the threads
        void **tmp = new void*[2*d_N_threads_diff];
        int j = d_N_threads;
        for (int i=0; i<d_N_threads_diff; i++) {
            d_N_threads++;
            tmp[0+2*i] = this;
            tmp[1+2*i] = reinterpret_cast<void*>(static_cast<size_t>(j));
            bool error = false;
            set_bit(d_cancel,j,true);
            #ifdef USE_WINDOWS
                d_hThread[j] = (HANDLE)_beginthread( create_new_thread, 0, (void *) &tmp[2*i]);
                error = d_hThread==(HANDLE)(-1);
            #elif defined(USE_LINUX) || defined(USE_MAC)
                int rtn = pthread_create( &d_hThread[j], &attr, (void *(*)(void*)) create_new_thread, (void *) &tmp[2*i]);
                error = rtn!=0;
            #else
                #error Not programmed
            #endif
            if ( error ) {
                pout << "Warning: Only able to create " << i << " threads\n";
                break;
            }
            j++;
        }
        // Wait for all of the threads to finish initialization
        while ( 1 ) {
            unlock_mutex(d_lock_queue);
            Sleep(25);
            lock_mutex(d_lock_queue);
            bool wait = false;
            for (int i=0; i<MAX_NUM_THREADS/64; i++) {
                if ( d_cancel[i] != 0 )
                    wait = true;
            }
            if ( !wait ) 
                break;
        }
        // Delete the thread attributes (linux only)
        #if defined(USE_LINUX) || defined(USE_MAC)
            pthread_attr_destroy(&attr);
        #endif
        // Release the lock
        unlock_mutex(d_lock_queue);
        Sleep(25);
        delete [] tmp;
    } else if ( d_N_threads_diff < 0 ) {
        // Reduce the number of threads
        if ( num_worker_threads==0 ) {
            // Special case if we want to delete all of the threads
            wait_pool_finished();
        }
        // Lock the mutex for the deletion of existing threads
        lock_mutex(d_lock_queue);
        // Tell the threads to shutdown
        for (int i=0; i>d_N_threads_diff; i--)
            set_bit(d_cancel,d_N_threads-1+i,true);
        #ifdef USE_WINDOWS
            // Release the lock
            unlock_mutex(d_lock_queue);
            // Wake all threads to process the shutdown (Doesn't require blocking)
            for (int i=0; i<d_N_threads; i++) {
                ResumeThread(d_hThread[i]);
            }
        #elif defined(USE_LINUX) || defined(USE_MAC)
            // Wake all threads to process the shutdown
            int error = pthread_cond_broadcast(d_queue_not_empty);
            if ( error != 0 )
                perr << "Error in signaling thread";
            // Release the lock
            unlock_mutex(d_lock_queue);
        #else
            #error Not programmed
        #endif
        Sleep(25);
        // Wait for all of the threads to close
        #ifdef USE_WINDOWS
            int j = d_N_threads+d_N_threads_diff;
            WaitForMultipleObjects( -d_N_threads_diff, &d_hThread[j], 1, 10000 );
        #elif defined(USE_LINUX) || defined(USE_MAC)
            for (int i=0; i>d_N_threads_diff; i--) {
                int rtn = pthread_join(d_hThread[d_N_threads-1+i],NULL);
                if ( rtn != 0 ) {
                    perr << "error\n";
                    perr << "Error joining threads";
                }
            }
        #else
            #error Not programmed
        #endif
        for (int i=0; i>d_N_threads_diff; i--) {
            set_bit(d_cancel,d_N_threads-1+i,false);
            d_hThread[d_N_threads-1+i] = 0;
            d_ThreadId[d_N_threads-1+i] = ~((size_t)0);
        }
        d_N_threads += d_N_threads_diff;
    }
    if ( d_N_threads == 0 )
        return;
    // Get the default thread affinity to use
    std::vector<int> cpus;
    int tmp = global_OS_behavior;
    global_OS_behavior = 1;
    OS_warning("Dummy message (should not print)");
    try {
        cpus = ThreadPool::getProcessAffinity();
    } catch(...) {
        pout << "Warning: Unable to get default cpus for thread affinities\n";
    }
    if ( !cpus.empty() && N_procs>0 ) {
        cpus.resize(N_procs);
        for (int i=0; i<N_procs; i++)
            cpus[i] = procs[i];
    }
    // Set the affinity model and the associated thread affinities
    // Note: not all OS's support setting the thread affinities
    std::vector<std::vector<int> > t_procs(d_N_threads);
    std::string affinity(affinity2);
    if ( cpus.empty() ) {
        // We do not have a list of cpus to use, do nothing (OS not supported)
    } else if ( affinity=="none" ) {
        // We are using the default thread affinities (all threads get all procs of the program)
        for (int i=0; i<d_N_threads; i++)
            t_procs[i] = cpus;
    } else if ( affinity=="independent" ) {
        // We want to use an independent set of processors for each thread
        if ( (int) cpus.size() == d_N_threads ) {
            // The number of cpus matches the number of threads
            for (int i=0; i<d_N_threads; i++)
                t_procs[i] = std::vector<int>(1,cpus[i]);
        } else if ( (int) cpus.size() > d_N_threads ) {
            // There are more cpus than threads, threads will use more the one processor
            int N_procs_thread = (cpus.size()+d_N_threads-1)/d_N_threads;
            size_t k = 0;
            for (int i=0; i<d_N_threads; i++) {
                for (int j=0; j<N_procs_thread && k<cpus.size(); j++) {
                    t_procs[i].push_back( cpus[k] );
                    k++;
                }
            }
        } else {
            // There are fewer cpus than threads, threads will share a processor
            int N_threads_proc = (cpus.size()+d_N_threads-1)/cpus.size();
            for (int i=0; i<d_N_threads; i++)
                t_procs[i].push_back( cpus[i/N_threads_proc] );
        }
    } else {
        global_OS_behavior = tmp;
        throw std::logic_error("Unknown affinity model");
    }
    try {
        for (int i=0; i<d_N_threads; i++) {
            ThreadPool::setThreadAffinity( i, t_procs[i] );
            std::vector<int> cpus2 = getThreadAffinity( i );
            if ( cpus2 != t_procs[i] )
                pout << "Warning: error setting affinities (failed to set)\n";
        }
    } catch (...) {
        pout << "Warning: error setting affinities (exception)\n";
    }
    global_OS_behavior = tmp;
}


/******************************************************************
* Get an item in the work queue that is ready to be processed     *
******************************************************************/
int ThreadPool::getThreadNumber() const
{
    size_t id = getThreadId();
    int index = 0;
    for (int i=0; i<d_N_threads; i++) {
        if ( d_ThreadId[i]==id )
            index = i+1;
    }
    return index;
}


/******************************************************************
* Get an item in the work queue that is ready to be processed     *
******************************************************************/
short int ThreadPool::get_work_item( )
{
    const thread_id_t *ids = const_cast<const thread_id_t*>(d_queue_ids);
    const queue_list_struct *list = const_cast<const queue_list_struct*>(d_queue_list);
    short int index = d_queue_head;
    short int index2 = check_dependecies(list,ids,index);
    while ( index2==-1 && index!=-1 ) {
        index = d_queue_list[index].next;
        index2 = index==-1 ? -1:check_dependecies(list,ids,index);
    }
    return index2;
}
inline short int ThreadPool::check_dependecies( const ThreadPool::queue_list_struct *list, 
    const thread_id_t *queue, short int index )
{
    if ( index==-1 )
        return -1;
    WorkItem* work = reinterpret_cast<WorkItem*>(queue[index].d_work);
    // Loop through the dependencies, removing any that have finished,
    // and search for any that have not started (keeping the one with the fewest dependencies)
    size_t N_active = 0;
    thread_id_t* ids = work->d_ids;
    short int index2 = index;
    int N_dependencies = static_cast<int>(work->d_N_ids);
    for (int i=N_dependencies-1; i>=0; i--) {
        WorkItem* work2 = reinterpret_cast<WorkItem*>(ids[i].d_work);
        char state = work2->d_state;
        if ( state==0 ) {
            // We found a new potential item to process
            index2 = work2->d_tpool_index;
            index2 = check_dependecies(list,queue,index2);
            if ( index2 != -1 )
                break;
        } else if ( state==1 || state==-1 ) {
            // We found an item that is processing
            N_active++;
        } else if ( state==2 ) {
            // The item has finished
            ids[i].reset();
            std::swap(ids[i],ids[work->d_N_ids-1]);
            work->d_N_ids--;
            continue;
        }
    }
    if ( N_active>0 ) {
        // Some dependencies are working, choose a different work item
        index2 = -1;
    }
    return index2;
}


/******************************************************************
* This is the function that controls the individual thread and    *
* allows it to do work.                                           *
******************************************************************/
void ThreadPool::tpool_thread(int thread_id) 
{
    if ( getThreadId()==0 )
        throw std::logic_error("Invalid thread id");
    bool shutdown = false;
    bool printInfo = false;
    d_ThreadId[thread_id] = getThreadId();
    // Acquire mutex 
    lock_mutex(d_lock_queue);
    if ( get_bit(d_active,thread_id) )
        throw std::logic_error("Thread cannot already be active");
    d_num_active++;
    set_bit(d_active,thread_id,true);
    set_bit(d_cancel,thread_id,false);
    if ( printInfo ) {
        // Print the pid
        printp("pid = %i\n",(int)getpid());
        // Print the processor affinities for the process
        try {
            std::vector<int> cpus = ThreadPool::getProcessAffinity();
            printp("%i cpus for current thread: ",(int)cpus.size());
            for (size_t i=0; i<cpus.size(); i++)
                printp("%i ",cpus[i]);
            printp("\n");
        } catch(...) {
            printp("Unable to get process affinity\n");
        }
    }
    // Check for shutdown
    shutdown = false;
    //pout << "Thread initialized\n";
    PROFILE_THREAD_START("thread active");
    while ( !shutdown ) {
        // Check if there is work to do
        if ( d_queue_size>0 ) {
            // Get next work item to process
            short int work_index = ThreadPool::get_work_item();
            if ( work_index==-1 ) {
                unlock_mutex(d_lock_queue);
                Sleep(0);
                lock_mutex(d_lock_queue);
                continue;
            }
            // Remove the work item from the queue
            #ifdef D_DEBUG
                short int cur = d_queue_list[work_index].position;
            #endif
            short int next = d_queue_list[work_index].next;
            short int prev = d_queue_list[work_index].prev;
            if ( prev==-1 ) {
                d_queue_head = next;
            } else {
                d_queue_list[prev].next = next;
            }
            if ( next!=-1 ) {
                d_queue_list[next].prev = prev;
            }
            --d_queue_size;
            #ifdef D_DEBUG
                if ( cur!=work_index || ( d_queue_size>0 && d_queue_head==-1 ) )
                    throw std::logic_error("Internal error with threadpool");
            #endif
            thread_id_t work_id = const_cast<thread_id_t&>(d_queue_ids[work_index]);
            d_queue_ids[work_index].reset();
            d_queue_list[work_index].reset();
            d_queue_list[work_index].next = d_queue_free;
            d_queue_free = work_index;
            WorkItem* work = reinterpret_cast<WorkItem*>(work_id.d_work);
            work->d_state = -1;
            // Release mutex
            unlock_mutex(d_lock_queue);
            // Start work here 
            PROFILE_THREAD_START("thread working");
            work->run();
            if ( work->d_state!=2 ) { throw std::logic_error("Work item is not changing state"); }
            PROFILE_THREAD_STOP("thread working");
            // Work finished, acquire mutex and remove it from the active list
            lock_mutex(d_lock_queue);
            // Check if any threads are waiting on the current work item
            for (int i=0; i<d_N_wait; i++) {
                wait_event_struct* wait = const_cast<wait_event_struct*>(d_wait[i]);
                bool found = false;
                if ( wait->ids.empty() ) {
                    // Special case where we just want to wait for any work items to finish
                    found = true;
                } else {
                    found = find_id( wait->ids, work_id );
                }
                if ( found ) {
                    wait_type event = 0;
                    volatile int* count = &(wait->count);
                    if ( *count == 1 )
                        event = const_cast<wait_type>(wait->wait_event);
                    --(*count);
                    if ( event != 0 )
                        SIGNAL_EVENT(event);
                }
            }
            // Check the signal count and signal if desired
            if ( d_signal_count > 0 ) {
                --d_signal_count;
                if ( d_signal_count == 0 )
                    SIGNAL_EVENT(d_wait_finished);
            }
        } else {
            int N_active = --d_num_active;
            set_bit(d_active,thread_id,false);
            // Alert main thread that a thread finished processing 
            if ( N_active==0 ) { 
                if ( d_signal_empty ) {
                    SIGNAL_EVENT(d_wait_finished);
                    d_signal_empty = false;
                }
            }
            // Wait for work
            PROFILE_THREAD_STOP2("thread active");
            #ifdef USE_WINDOWS
                unlock_mutex(d_lock_queue);
                SuspendThread(d_hThread[thread_id]);
                lock_mutex(d_lock_queue);
            #elif defined(USE_LINUX) || defined(USE_MAC)
                pthread_cond_wait(d_queue_not_empty,d_lock_queue);
            #endif
            PROFILE_THREAD_START2("thread active");
            ++d_num_active;
            set_bit(d_active,thread_id,true);
        }
        // Check if there is a shutdown requested
        shutdown = get_bit(d_cancel,thread_id);
    }
    PROFILE_THREAD_STOP("thread active");
    d_num_active--;
    set_bit(d_active,thread_id,false);
    // Release mutex
    unlock_mutex(d_lock_queue);
    return;
}



/******************************************************************
* This is the function that adds work to the thread pool          *
* Note: this version uses a last in - first out work scheduling.  *
******************************************************************/
void ThreadPool::add_work( size_t N, ThreadPool::WorkItem* work[], 
    const int* priority, ThreadPool::thread_id_t* ids ) 
{
    #if MONITOR_THREADPOOL_PERFORMANCE
        TIME_TYPE start_time_local;
        get_time(&start_time_local);
    #endif
    // If we have a very long list, break it up into smaller pieces to keep the threads busy
    const size_t block_size = MAX_QUEUED/4;
    if ( N > block_size ) {
        size_t N_sets = (N+block_size-1)/block_size;
        for (size_t i=0; i<N_sets; i++) {
            size_t index = i*block_size;
            size_t N2 = std::min<size_t>(block_size,N-index);
            add_work( N2, &work[index], &priority[index], &ids[index] );
        }
        return;
    }
    // Create the thread ids (can be done without blocking)
    for (size_t i=0; i<N; i++) {
        ids[i].reset(priority[i],advance_id(),work[i]);
        work[i]->d_tpool_index = -2;
    }
    // If there are no threads, perform the work immediately
    if ( d_N_threads < 1 ) {
        for (size_t i=0; i<N; i++) {
            work[i]->run();
        }
        return;
    }
    // Wait for enough room in the queue (doesn't need blocking since it isn't that precise)
    if ( N > static_cast<size_t>(MAX_QUEUED-d_queue_size) ) {
        int N_wait = static_cast<int>( N - (MAX_QUEUED-d_queue_size) );
        while ( N_wait > 0 ) {
            d_signal_count = static_cast<unsigned char>(std::min(N_wait,255));
            #ifdef USE_WINDOWS
                DWORD ret = WaitForSingleObject( d_wait_finished, INFINITE );
            #elif defined(USE_LINUX) || defined(USE_MAC)
                lock_mutex(d_lock_queue);
                if ( d_signal_count > 0 )
                    pthread_cond_wait(d_wait_finished,d_lock_queue);
                unlock_mutex(d_lock_queue);
            #else
                #error Not programmed
            #endif
            N_wait = static_cast<int>( N - (MAX_QUEUED-d_queue_size) );
        }
    }
    // Get the lock and add the work items
    lock_mutex(d_lock_queue);
    #if MONITOR_THREADPOOL_PERFORMANCE
        TIME_TYPE stop_time_local;
        get_time(&stop_time_local);
        total_add_work_time[0] += get_diff(start_time_local,stop_time_local,frequency);
    #endif
    // Next create the work items and add them to the queue
    for (size_t i=0; i<N; i++) {
        queue_list_struct *work_item = const_cast<queue_list_struct*>(&d_queue_list[d_queue_free]);
        d_queue_free = work_item->next;
        work_item->next = -1;
        work_item->prev = -1;
        d_queue_ids[work_item->position] = ids[i];
        reinterpret_cast<WorkItem*>(ids[i].d_work)->d_tpool_index = work_item->position;
        if ( d_queue_head==-1 ) {
            d_queue_head = work_item->position;
        } else if ( ids[i] > d_queue_ids[d_queue_list[d_queue_head].position] ) {
            work_item->next = d_queue_head;
            d_queue_list[d_queue_head].prev = work_item->position;
            d_queue_head = work_item->position;
        } else {
            short int prev = d_queue_head;
            short int cur = d_queue_list[prev].next;
            while ( cur!=-1 ) {
                if ( d_queue_ids[cur] < ids[i] )
                    break;
                prev = cur;
                cur = d_queue_list[prev].next;
            }
            work_item->prev = prev;
            work_item->next = cur;
            if ( cur != -1 )
                d_queue_list[cur].prev = work_item->position;
            d_queue_list[prev].next = work_item->position;
        }
        ++d_queue_size;
    }
    int num_active2 = d_num_active;       // Copy the number of active threads to a local variable
    unlock_mutex(d_lock_queue);
    #if MONITOR_THREADPOOL_PERFORMANCE
        get_time(&stop_time_local);
        total_add_work_time[1] += get_diff(start_time_local,stop_time_local,frequency);
    #endif
    // Activate sleeping threads
    #ifdef USE_WINDOWS
        for (int i=0; i<d_N_threads; i++) {
            if ( num_active2 == d_N_threads ) {
                // All threads are active, no need to activate
                break;
            } else if ( d_queue_size == 0 ) {
                // Queue is empty, no need to activate
                break;
            } else if ( !get_bit(d_active,i) ) {
                // Thread is inactive, wake it
                ResumeThread(d_hThread[i]);
            }
        }
    #elif defined(USE_LINUX) || defined(USE_MAC)
        if ( num_active2 == d_N_threads ) {
            // All threads are active, no need to wake anybody
        } else if ( d_queue_size == 0 ) {
            // Queue is empty, no need to activate
        } else if ( N == 1 ) {
            // Added 1 item to the queue, wake 1 worker
            int error = pthread_cond_signal(d_queue_not_empty);
            if ( error != 0 )
                perr << "Error in signaling thread";
        } else {
            // Added multple items in the queue, wake all workers
            int error = pthread_cond_broadcast(d_queue_not_empty);
            if ( error != 0 )
                perr << "Error in signaling thread";
        }
    #endif
    #if MONITOR_THREADPOOL_PERFORMANCE
        get_time(&stop_time_local);
        total_add_work_time[2] += get_diff(start_time_local,stop_time_local,frequency);
    #endif
}




/******************************************************************
* This function checks if the work item has finished              *
******************************************************************/
bool ThreadPool::isFinished(ThreadPool::thread_id_t id) const
{
    if ( !isValid(id) ) {
        // The thread id is not valid
        return false;
    }
    return reinterpret_cast<WorkItem*>(id.d_work)->d_state==2;
}



/******************************************************************
* This function removes a finished work item                      *
******************************************************************/
ThreadPool::WorkItem* ThreadPool::getFinishedWorkItem(ThreadPool::thread_id_t id) const
{
    if ( !isValid(id) ) 
        return NULL;
    if ( reinterpret_cast<WorkItem*>(id.d_work)->d_state!=2 )
        return NULL;
    // Return the result
    WorkItem* work = reinterpret_cast<WorkItem*>(id.d_work);
    return work;
}



/******************************************************************
* This function waits for a some of the work items to finish      *
******************************************************************/
static inline void check_finished( size_t N_work, const ThreadPool::thread_id_t *ids, size_t& N_finished, bool* finished)
{
    for (size_t k=0; k<N_work; k++) {
        if ( !finished[k] && ids[k].finished() ) {
            N_finished++;
            finished[k] = true;
        }
    }
}
int ThreadPool::wait_some(size_t N_work, const ThreadPool::thread_id_t *ids, size_t N_wait, bool* finished) const
{
    // Check the inputs
    if ( N_wait<=0 || N_wait>N_work ) {
        printp("Invalid arguments in thread pool wait (%i,%i)\n",(int)N_work,(int)N_wait);
        return -1;
    }
    size_t N_finished = 0;
    memset(finished,0,N_work*sizeof(bool));
    // Check that all the ids are valid
    size_t next_id = d_id_assign-1;
    for (size_t k=0; k<N_work; k++) {
        if ( !ids[k].initialized() ) {
            finished[k] = true;
            N_finished++;
        }
        size_t local_id = ids[k].getLocalID();
        bool test = local_id==0 || local_id>MAXID64 || local_id<=next_id;
        test = test && !finished[k];
        if ( test ) {
            printp("Invalid ids for wait\n");
            return -1;
        }
    }
    // Check which ids have finished
    check_finished(N_work,ids,N_finished,finished);
    // If enough ids have finished return
    if ( N_finished >= N_wait ) {
        return 0;
    }
    // Acquire the lock and update the finished list
    // It is possible that in the time required to acquire the lock, the work items may finish
    lock_mutex(d_lock_queue);
    check_finished(N_work,ids,N_finished,finished);
    if ( N_finished >= N_wait ) {
        unlock_mutex(d_lock_queue);
        return 0;
    }
    // Create the wait event struct
    wait_event_struct* tmp = new wait_event_struct(&wait_pool);
    wait_type event = tmp->wait_event;
    tmp->count = static_cast<int>(N_wait-N_finished);
    tmp->ids.reserve(N_wait-N_finished);
    for (size_t k=0; k<N_work; k++) {
        if ( !finished[k] )
            tmp->ids.push_back(ids[k]);
    }
    quicksort(tmp->ids);
    d_wait[d_N_wait] = tmp;
    d_N_wait++;
    // Wait for a signal indicating that a thread has finished
    #ifdef USE_WINDOWS
        unlock_mutex(d_lock_queue);
        DWORD ret = WaitForSingleObject( event, INFINITE );
        lock_mutex(d_lock_queue);
    #elif defined(USE_LINUX) || defined(USE_MAC)
        pthread_cond_wait(event,d_lock_queue);
    #endif
    // Check for remaining references to the wait struct and delete the structure
    for (int k=0; k<d_N_wait; k++) {
        if ( d_wait[k] == tmp ) {
            for (int m=k+1; m<d_N_wait; m++)
                d_wait[m-1] = d_wait[m];
            d_wait[d_N_wait-1] = NULL;
            break;
        }
    }
    d_N_wait--;
    delete tmp;
    unlock_mutex(d_lock_queue);
    // Update the ids that have finished
    check_finished(N_work,ids,N_finished,finished);
    if ( N_finished<N_wait && N_work!=0 ) {
        throw std::logic_error("Internal error: failed to wait");
    }
    return 0;
}



/******************************************************************
* This function waits for all of the threads to finish their work *
******************************************************************/
void ThreadPool::wait_pool_finished() const 
{
    // First check that we are not one of the threads
    if ( isMemberThread() ) {
        throw std::logic_error("Member thread attempted to call wait_pool_finished");
    }
    lock_mutex(d_lock_queue);
    // Wait for all threads to finish their work
    while ( d_num_active>0 || d_queue_size>0 ) {
        d_signal_empty = true;
        #ifdef USE_WINDOWS
            unlock_mutex(d_lock_queue);
            DWORD ret = WaitForSingleObject( d_wait_finished, INFINITE );
            lock_mutex(d_lock_queue);
        #elif defined(USE_LINUX) || defined(USE_MAC)
            pthread_cond_wait(d_wait_finished,d_lock_queue);
        #else
            #error Not programmed
        #endif
    }
    d_signal_empty = false;
    unlock_mutex(d_lock_queue);
}



/******************************************************************
* These functions create the unique id to assign each work item   *
* If id is a 32-bit number we have 4e9 possible work items        *
* If id is a 64-bit number we have 9e19 possible work items and   *
*    we have some checking that will catch some invalid ids       *
******************************************************************/
inline void ThreadPool::initialize_id() 
{
    // Note that the best option is to use a 64-bit integer
    if ( sizeof(size_t)==8 ) {
        // Set the starting value to 2^56-3
        d_id_assign = MAXID64;
    } else if ( sizeof(size_t)==4 ) {
        // Set the starting value to 2^32-3
        d_id_assign = MAXID32;
    } else {
        throw std::logic_error("Internal error: failed to initialize ids");
    }
}
inline size_t ThreadPool::advance_id() 
{
    size_t id = AtomicOperations::atomic_decrement( &d_id_assign );
    if ( id==0 )
        throw std::logic_error("Ran out of valid ids");
    return id;
}


/******************************************************************
* Function to check if the current thread is a member thread      *
******************************************************************/
inline bool ThreadPool::isMemberThread() const
{
    size_t id = getThreadId();
    for (int i=0; i<d_N_threads; i++) {
        if ( id==d_ThreadId[i] )
            return true;
    }
    return false;
}


/******************************************************************
* Member functions of wait_event_struct                           *
******************************************************************/
ThreadPool::wait_event_struct::wait_event_struct( wait_pool_struct* wait_pool )
{
    count = 0;
    ThreadId = getThreadId();
    d_wait_pool = wait_pool;
    wait_event = d_wait_pool->pop();
}
ThreadPool::wait_event_struct::~wait_event_struct( )
{
    d_wait_pool->push(wait_event);
}


/******************************************************************
* Member functions of wait_pool_struct                            *
******************************************************************/
ThreadPool::wait_pool_struct::wait_pool_struct( )
{
    d_size = 16;
    d_count = 0;
    d_pool = new wait_type[d_size];
    memset(const_cast<wait_type*>(d_pool),0,d_size*sizeof(wait_type));
    d_lock = create_mutex( );
}
ThreadPool::wait_pool_struct::~wait_pool_struct( )
{
    for (size_t i=0; i<d_count; i++) {
        #ifdef USE_WINDOWS
            CloseHandle(d_pool[i]);
        #elif defined(USE_LINUX) || defined(USE_MAC)
            pthread_cond_destroy(d_pool[i]);
            delete d_pool[i];
        #else
            #error Not programmed
        #endif
    }
    delete [] d_pool;
    destroy_mutex( d_lock );
    d_size = 0;
    d_count = 0;
    d_pool = 0;
    d_lock = 0;
}
void ThreadPool::wait_pool_struct::push( ThreadPool::wait_type event )
{
    lock_mutex(d_lock);
    if ( d_count >= d_size ) {
        volatile wait_type* tmp = d_pool;
        d_pool = new wait_type[2*d_size];
        memset((void*)d_pool,0,2*d_size*sizeof(wait_type));
        memcpy((void*)d_pool,(void*)tmp,d_size*sizeof(wait_type));
        delete [] d_pool;
        d_size = 2*d_size;
    }
    d_pool[d_count] = event;
    ++d_count;
    unlock_mutex(d_lock);
}
ThreadPool::wait_type ThreadPool::wait_pool_struct::pop( )
{
    lock_mutex(d_lock);
    wait_type event = 0;
    if ( d_count == 0 ) {
        #ifdef USE_WINDOWS
            event = CreateEvent(NULL,FALSE,FALSE,NULL);
        #elif defined(USE_LINUX) || defined(USE_MAC)
            event = new pthread_cond_t;
            int error = pthread_cond_init(event,NULL);
            if ( error == -1 )
                std::logic_error("Error creating wait_event");
        #else
            #error Not programmed
        #endif
    } else {
        event = d_pool[d_count-1];
        --d_count;
    }
    unlock_mutex(d_lock);
    return event;
}


/******************************************************************
* templated quicksort routine                                     *
******************************************************************/
template <class T>
void quicksort(std::vector<T> &x)
{
    int n = (int) x.size();
    if ( n <= 1 )
        return;
    T *arr = &x[0];
    bool test;
    int i, ir, j, jstack, k, l, istack[100];
    T a, tmp_a;
    jstack = 0;
    l = 0;
    ir = n-1;
    while (1) {
        if ( ir-l < 7 ) {             // Insertion sort when subarray small enough.
            for ( j=l+1; j<=ir; j++ ) {
                a = arr[j];
                test = true;
                for (i=j-1; i>=0; i--) {
                    if ( arr[i] < a ) {
                        arr[i+1] = a;
                        test = false;
                        break;
                    }
                    arr[i+1] = arr[i];
                }
                if ( test ) {
                    i = l-1;
                    arr[i+1] = a;
                }
            }
            if ( jstack==0 )
                return;
            ir = istack[jstack];    // Pop stack and begin a new round of partitioning.
            l = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l+ir)/2;           // Choose median of left, center and right elements as partitioning
                                    // element a. Also rearrange so that a(l) < a(l+1) < a(ir).
            tmp_a = arr[k];
            arr[k] = arr[l+1];
            arr[l+1] = tmp_a;
            if ( arr[l]>arr[ir] ) {
                tmp_a = arr[l];
                arr[l] = arr[ir];
                arr[ir] = tmp_a;
            }
            if ( arr[l+1] > arr[ir] ) {
                tmp_a = arr[l+1];
                arr[l+1] = arr[ir];
                arr[ir] = tmp_a;
            }
            if ( arr[l] > arr[l+1] ) {
                tmp_a = arr[l];
                arr[l] = arr[l+1];
                arr[l+1] = tmp_a;
            }
            // Scan up to find element > a
            j = ir;
            a = arr[l+1];           // Partitioning element.
            for (i=l+2; i<=ir; i++) { 
                if ( arr[i]<a ) 
                    continue;
                while ( arr[j]>a )  // Scan down to find element < a.
                    j--;
                if ( j < i )
                    break;          // Pointers crossed. Exit with partitioning complete.
                tmp_a = arr[i];     // Exchange elements of both arrays.
                arr[i] = arr[j];
                arr[j] = tmp_a;
            }
            arr[l+1] = arr[j];      // Insert partitioning element in both arrays.
            arr[j] = a;
            jstack += 2;
            // Push pointers to larger subarray on stack, process smaller subarray immediately.
            if ( ir-i+1 >= j-l ) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack] = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
}


/************************************************************************
* Function to find the id in a sorted vector                            *
************************************************************************/
inline bool find_id(const std::vector<ThreadPool::thread_id_t> &x_in, const ThreadPool::thread_id_t &id ) 
{
    if ( x_in.empty() )
        return false;
    size_t n = x_in.size();
    const ThreadPool::thread_id_t *x = &x_in[0];   // Use the pointer for speed
    if ( n<4 ) {
        for (size_t i=0; i<n; i++) {
            if ( x[i] == id ) 
                return true;
        }
    }
    // Check if value is within the range of x
    if ( id == x[0] )
        return true;
    if ( id < x[0] )
        return false;
    if ( id == x[n-1] )
        return true;
    if ( id > x[n-1] )
        return false;
    // Perform the search
    size_t lower = 0;
    size_t upper = n-1;
    size_t index;
    while ( (upper-lower) != 1 ) {
        index = (upper+lower)/2;
        if ( x[index] == id )
            return true;
        if ( x[index] >= id )
            upper = index;
        else
            lower = index;
    }
    return false;
}


/************************************************************************
* Function to add dependencies to the work item                         *
************************************************************************/
void ThreadPool::WorkItem::add_dependencies( size_t N, const ThreadPool::thread_id_t* ids)
{
    if ( d_tpool_index != -1 ) {
        // The item has already been added to the threadpool, 
        // we are not allowed to add dependencies
        throw std::logic_error("Cannot add dependency to work item once it has been added the the threadpool");
    }
    if ( static_cast<size_t>(d_N_ids)+N > 0xFFFF ) {
        throw std::logic_error("Cannot add more than 65000 dependencies");
    }
    for (size_t i=0; i<N; i++) {
        if ( ids[i].initialized() ) {
            if ( d_N_ids >= d_size ) {
                thread_id_t* tmp = d_ids;
                unsigned int N2 = d_size;
                if ( N2 == 0 ) { N2 = 8; }
                while ( N2 <= d_N_ids )
                    N2 *= 2;
                d_ids = new thread_id_t[N2];
                for (size_t i=0; i<d_N_ids; i++)
                    std::swap(d_ids[i],tmp[i]);
                delete [] tmp;
                d_size = N2;
            }
            d_ids[d_N_ids+i] = ids[i];
            d_N_ids++;
        }
    }
}


