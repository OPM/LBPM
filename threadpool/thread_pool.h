// Copyright © 2004 Mark Berrill. All Rights Reserved. This work is distributed with permission,
// but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#ifndef included_ThreadPool
#define included_ThreadPool
#include <stdio.h>
#include <typeinfo>
#include <iostream>
#include <stdarg.h>
#include <string.h>
#include <vector>
#include <stdexcept>
#include <map>

#include "threadpool/atomic_helpers.h"


// Choose the OS 
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    // Using windows
    #define USE_WINDOWS
    #include <stdlib.h>
    #include <windows.h>
    #include <process.h>
    #define NOMINMAX
    // Disable warning: the inline specifier cannot be used when a friend 
    // declaration refers to a specialization of a function template
    #pragma warning(disable:4396)
#elif defined(__APPLE__)
    // Using MAC
    //   https://developer.apple.com/library/mac/#releasenotes/Performance/RN-AffinityAPI
    //   http://plugins.svn.wordpress.org/wp-xhprof-profiler/trunk/facebook-xhprof/extension/xhprof..c
    #define USE_MAC
    #include <unistd.h>
    #include <mach/mach_init.h>
    #include <mach/thread_policy.h>
    #define cpu_set_t thread_affinity_policy_data_t
    #define CPU_SET(cpu_id, new_mask) \
        *new_mask.affinity_tag = (cpu_id + 1)
    #define CPU_ZERO(new_mask)                 \
        (*(new_mask)).affinity_tag = THREAD_AFFINITY_TAG_NULL
    #define sched_setaffinity(pid, size, mask)       \
        thread_policy_set(mach_thread_self(), THREAD_AFFINITY_POLICY, mask, \
                          THREAD_AFFINITY_POLICY_COUNT)
    #define sched_getaffinity(pid, size, mask) \
        thread_policy_get(mach_thread_self(), THREAD_AFFINITY_POLICY, mask, \
                          THREAD_AFFINITY_POLICY_COUNT)
    /*
    #define CPU_ZERO(new_mask) \
        *new_mask.affinity_tag == THREAD_AFFINITY_TAG_NULL
    #define SET_AFFINITY(pid, size, mask) \
        thread_policy_set(mach_thread_self(), THREAD_AFFINITY_POLICY, mask, THREAD_AFFINITY_POLICY_COUNT)
    #define GET_AFFINITY(pid, size, mask) \
        thread_policy_get(mach_thread_self(), THREAD_AFFINITY_POLICY, mask, THREAD_AFFINITY_POLICY_COUNT)
    */
#elif defined(__linux) || defined(__unix) || defined(__posix)
    // Using Linux
    #define USE_LINUX
    #include <pthread.h>
    #include <unistd.h>
#else
    #error Unknown OS
#endif


// Set some definitions
#define MAX_NUM_THREADS 128         // The maximum number of threads (must be a multiple of 64)
#define MAX_QUEUED 1024             // The maximum number of items in the work queue at any moment
#define MAX_WAIT  128               // The maximum number of active waits at any given time




/** \class Mutex
 * \brief Functions for locking/unlocking a mutex
 * \details This class provides basic routines for creating, 
 *    locking, and unlocking a mutex <BR>
 *    The lock may be recursive, meaning that the same thread
 *    may lock and unlock the lock multiple times before releasing it.
 *    In this case unlock must be called the same number of times before
 *    another thread may lock the mutex.
 */
class Mutex {
public:
    //! Empty constructor (equivilent to Mutex(false) )
    Mutex();
    /** Default constructor
     * \param recursive     If set to true a thread may repeated lock a mutex.
     *                      If set to false an attept to repeatedly lock will throw an error.*/
    Mutex(bool recursive);
    //! Destructor
    ~Mutex();
    //! Copy constructor
    Mutex(const Mutex &);
    //! Assignment operator
    Mutex& operator=(const Mutex&);
    //! Lock the mutex
    void lock() const;
    //! Unlock the mutex
    void unlock() const;
    //! Try to lock the mutex and return true if successful
    bool tryLock() const;
    //! Return true if we already own the lock
    bool ownLock() const;
private:
    bool d_recursive;               // Is the lock recursive (this attribute cannot be changed)
    volatile int* d_count;          // Number of copies of the mutex
    volatile int* d_lock_count;     // Number of times a thread has locked the mutex
    volatile size_t* d_thread;      // Pointer to the thread id that owns the lock
    #ifdef USE_WINDOWS
        CRITICAL_SECTION *d_lock;
    #elif defined(USE_LINUX) || defined(USE_MAC)
        pthread_mutex_t *d_lock;
    #else
        #error Unknown OS
    #endif
friend class ThreadPool;
};


/** \class ThreadPool
 *
 * \brief This is a concrete class that provides for a basic thread pool. 
 * \details This class implements a basic thread pool that can be used for a wide variety of applications.
 * An example call usage is provided below.  The ability to return a value is provided.  Note that there
 * is a small overhead to using this functionality. <BR>
 * <pre>Example: <BR>
 *    Existing function call:
 *       double x = myfun_1(a,b);
 *       double y = myfun_2(c,d); <BR>
 *    Threaded call (processing in parallel): 
 *       thread_id_t ids[2]; 
 *       ids[0] = TPOOL_ADD_WORK( tpool, myfun_1, (a,b) ); 
 *       ids[1] = TPOOL_ADD_WORK( tpool, myfun_2, (c,d) ); 
 *       int error = wait_all(2,ids);
 *       double x = getFunctionRet(ids[0]);
 *       double y = getFunctionRet(ids[1]); <BR>
 *   </pre>
 */
class ThreadPool {

public:

    //! Convience typedef
    typedef unsigned long long int uint64;

    //! Function to get a unique id for the current thread
    static inline size_t getThreadId();


public:
    ///// Member classes

    /** \class thread_id_t
     *
     * \brief This a class to hold the work item id
     * \details This class hold the id of the work item that is being processed by the thread pool.
     *      It is created when a work item is added to the thread pool and is used by various routines within the thread pool.
     */
    class thread_id_t {
        public:
            //! Empty constructor
            inline thread_id_t( );
            //! Destructor
            inline ~thread_id_t( );
            //! Copy constructors
            inline thread_id_t( const thread_id_t& rhs );
            inline thread_id_t& operator=( const thread_id_t& rhs ) volatile;
            #ifndef USE_WINDOWS
                inline thread_id_t( const volatile thread_id_t& rhs );
                inline thread_id_t& operator=( const thread_id_t& rhs );
                inline thread_id_t& operator=( const volatile thread_id_t& rhs );
                inline thread_id_t& operator=( const volatile thread_id_t& rhs ) volatile;
            #endif
            // Overload key operators
            inline bool operator==(const thread_id_t& rhs ) const { return d_id==rhs.d_id; }
            inline bool operator!=(const thread_id_t& rhs ) const { return d_id!=rhs.d_id; }
            inline bool operator>=(const thread_id_t& rhs ) const { return d_id>=rhs.d_id; }
            inline bool operator<=(const thread_id_t& rhs ) const { return d_id<=rhs.d_id; }
            inline bool operator> (const thread_id_t& rhs ) const { return d_id>rhs.d_id;  }
            inline bool operator< (const thread_id_t& rhs ) const { return d_id<rhs.d_id;  }
            inline bool operator==(const volatile thread_id_t& rhs ) const volatile { return d_id==rhs.d_id; }
            inline bool operator!=(const volatile thread_id_t& rhs ) const volatile { return d_id!=rhs.d_id; }
            inline bool operator>=(const volatile thread_id_t& rhs ) const volatile { return d_id>=rhs.d_id; }
            inline bool operator<=(const volatile thread_id_t& rhs ) const volatile { return d_id<=rhs.d_id; }
            inline bool operator> (const volatile thread_id_t& rhs ) const volatile { return d_id>rhs.d_id;  }
            inline bool operator< (const volatile thread_id_t& rhs ) const volatile { return d_id<rhs.d_id;  }
            //! Reset the id back to a NULL id
            inline void reset() volatile;
            inline void reset();
            //! Check if the work has finished
            inline bool finished( ) const;
        private:
            // Default constructor
            inline void reset( int priority, size_t local_id, void* work );
            // Get the local id
            inline size_t getLocalID() const;
            // Get the priority
            inline int getPriority() const; 
            // Check if the id is initialized
            inline bool initialized() const volatile { return d_id!=0x0FFFFFFFFFFFFFFF; }
            // Friends
            friend class ThreadPool;
            template<typename T> friend void std::swap(T&, T&);
            // Data
            uint64 d_id;                                        // 64-bit data to store id
            AtomicOperations::int32_atomic* volatile d_count;   // Reference count
            void* d_work;                                       // Pointer to the work item
    };


    //! Base class for the work item (users should derive from WorkItemRet)
    class WorkItem {
        public:
            //! Function to run the routine
            virtual void run()=0;
            //! Will the routine return a result
            bool has_result() const { return d_has_result; }
            //! Empty deconstructor
            virtual ~WorkItem() { delete [] d_ids; d_ids=NULL; d_N_ids=0; d_size=0; }
            //! Get the number of work ids that this work item depends on
            inline size_t get_N_dependencies() const { return d_N_ids; }
            //! Return the list of work ids that we depend on
            std::vector<ThreadPool::thread_id_t> get_dependencies() const;
            /*!
             * \brief Add a work item to the list of dependencies
             * \param id    Id of the work item to add
             */
            void add_dependency( const ThreadPool::thread_id_t& id ) { add_dependencies(1,&id); }
            /*!
             * \brief Add a list of work item to the list of dependencies
             * \param ids   Ids of the work item to add
             */
            inline void add_dependencies( const std::vector<ThreadPool::thread_id_t>& ids ) { 
                if ( !ids.empty() ) { add_dependencies(ids.size(),&ids[0]); }
            }
            /*!
             * \brief Add a list of work item to the list of dependencies
             * \param N     Number of items to add
             * \param ids   Ids of the work item to add
             */
            void add_dependencies( size_t N, const ThreadPool::thread_id_t* ids);
        protected:
            friend class ThreadPool;
            inline WorkItem(): d_has_result(false), d_state(0), d_tpool_index(-1), d_N_ids(0), d_size(0), d_ids(NULL) {}
            bool d_has_result;          // Derived classes must set the result flag (true: has a result)
            volatile char d_state;      // Derived classes must set the state (0: not scheduled, -1: scheduled, 1: started, 2: finished)
            short int d_tpool_index;    // Index of the item in the thread pool (-1: not added)
        private:
            WorkItem(const WorkItem&);          // Private copy constructor
            WorkItem& operator=(const WorkItem&); // Private assignment operator
            short unsigned int d_N_ids;         // Number of dependencies
            short unsigned int d_size;          // Size of d_ids
            thread_id_t* d_ids;                 // Pointer to id list
    };


    /*!
     * \brief   Class to define a work item returning a variable
     * \details This is the class that defines a work item to be processed.  Users may derive their own 
     * class and add work using the add_work routine, or can use the TPOOL_ADD_WORK macro.
     * Note: this class is templated on the return argument type and may be a void type.
     */
    template <typename return_type> 
    class WorkItemRet: public ThreadPool::WorkItem {
        public:
            //! Run the work item
            virtual void run()=0;
            //! Return the results
            return_type get_results() const { return d_result; }
            //! Virtual destructor
            virtual ~WorkItemRet() {}
        protected:
            return_type d_result;
            inline WorkItemRet(): WorkItem() { d_has_result = true; }
        private:
            WorkItemRet(const WorkItemRet&);            // Private copy constructor
            WorkItemRet& operator=(const WorkItemRet&); // Private assignment operator
    };


public:
    ///// Member functions

    //! Empty constructor
    ThreadPool() 
    {
        // Note: we need the constructor in the header to ensure that check_startup
        //       is able to check for changes in the byte alignment
        check_startup(sizeof(ThreadPool));
        initialize(0,"none",0,NULL);
        if ( !is_valid(this) )
            throw std::logic_error("Thread pool is not valid");
    }


    /*!
     *  Constructor that initialize the thread pool with N threads
     * @param N    The desired number of worker threads
     * @param affinity          The affinity scheduler to use:
     *                          none - Let the OS handle the affinities (default)
     *                          independent - Give each thread an independent set of processors
     * @param procs             The processors to use (defaults to the process affinitiy list)
     */
    ThreadPool( const int N, const std::string& affinity="none", const std::vector<int>& procs=std::vector<int>() )
    {
        // Note: we need the constructor in the header to ensure that check_startup
        //       is able to check for changes in the byte alignment
        check_startup(sizeof(ThreadPool));
        const int* procs2 = procs.empty() ? NULL:(&procs[0]);
        initialize(N,affinity.c_str(),procs.size(),procs2);
        if ( !is_valid(this) )
            throw std::logic_error("Thread pool is not valid");
    }


    //! Destructor
    ~ThreadPool();


    //! Function to return the number of processors availible
    static int getNumberOfProcessors();


    //! Function to return the processor number that the current thread is running on
    static int getCurrentProcessor();


    //! Function to return the affinity of the current process
    static std::vector<int> getProcessAffinity();


    //! Function to set the affinity of the current process
    static void setProcessAffinity( std::vector<int> procs );


    //! Function to return the affinity of the current thread
    static std::vector<int> getThreadAffinity();


    /*!
     *  Function to return the affinity of the given thread
     *  @param thread   The index of the thread
     */
    std::vector<int> getThreadAffinity( int thread ) const;


    /*!
     *  Function to set the affinity of the current thread
     *  @param procs    The processors to use
     */
    static void setThreadAffinity( std::vector<int> procs );


    /*!
     *  Set the given thread to have the specified affinity
     *  @param thread   The index of the thread
     *  @param procs    The processors to use
     */
    void setThreadAffinity( int thread, std::vector<int> procs ) const;


    //! Function to return the number of threads in the thread pool
    int getNumThreads() const { return d_N_threads; }


    /*!
     * \brief   Function to set the number of threads in the thread pool
     * \details  This function will change the number of worker threads in the ThreadPool
     *   to the number specified.  This function will immediately change the number of threads
     *   in the ThreadPool without checking the existing work unless the desired number of 
     *   threads is 0.  In this case, the function will wait for all work items to finish
     *   before deleting the existing work threads.
     *   Member threads may not call this function.
     * @param N                 The desired number of worker threads
     * @param affinity          The affinity scheduler to use:
     *                          none - Let the OS handle the affinities (default)
     *                          independent - Give each thread an independent set of processors
     * @param procs             The processors to use (defaults to the process affinitiy list)
     */
    inline void setNumThreads( const int N, const std::string& affinity="none", 
        const std::vector<int>& procs=std::vector<int>() )
    {
        const int* procs2 = procs.empty() ? NULL:(&procs[0]);
        setNumThreads(N,affinity.c_str(),procs.size(),procs2);
    }


    /*!
     * \brief   Function to return the current thread number
     * \details  This function will return the thread number of current active thread.
     *    If the thread is not a member of the thread pool, this function will return 0.
     */
    int getThreadNumber() const;


    //! Function to check if the work item is valid
    /*!
     * This function checks if the work item has a valid id.  
     *   Note: this function does not require blocking and will return immediately.  
     * @param id                The id of the work item
     */
    inline bool isValid(const thread_id_t& id) const;


    /*!
     * \brief    Function to check if the work item has finished processing
     * \details  This function checks if the work item has finished processing. 
     * @param id                The id of the work item
     */
    bool isFinished(thread_id_t id) const;


    /*!
     * \brief   Function to get the returned function value
     * \details This is the function returns the value that was returned from the working function.
     *   If the work item has not finished or was not found it will return 0.  
     * @param id                The id of the work item
     */
    template <class return_type> 
    inline return_type getFunctionRet(const thread_id_t& id) const;

    
    /*!
     * \brief   Function to add a work item
     * \details This function adds a work item to the queue
     *   Note: any thread may call this routine.
     * @param work              Pointer to the work item to add
     *                          Note that the threadpool will automatically destroy the item when finished
     * @param priority          A value indicating the priority of the work item (0-default)
     */
    inline thread_id_t add_work( ThreadPool::WorkItem* work, int priority=0);


    /*!
     * \brief   Function to add multiple work items
     * \details This function adds multiple work item to the queue
     *   Note: any thread may call this routine.
     * @param work              Vector of pointers to the work items to add
     *                          Note that the threadpool will automatically destroy the item when finished
     * @param priority          Vector of values indicating the priority of the work items
     */
    inline std::vector<thread_id_t> add_work( const std::vector<ThreadPool::WorkItem*>& work, 
        const std::vector<int>& priority=std::vector<int>() );


    /*!
     * \brief   Function to wait until a specific work item has finished
     * \details This is the function waits for a specific work item to finished.  It returns 0 if successful.
     *   Note: any thread may call this routine, but they will block until finished.
     *   For worker threads this may eventually lead to a deadlock.
     * @param id                The work item to wait for
     */
    inline int wait(thread_id_t id) const;


    /*!
     * \brief   Function to wait until any of the given work items have finished their work
     * \details This is the function waits for any of the given work items to finish. 
     *   If successful it returns the index of a finished work item (the index in the array ids).
     *   If unseccessful it will return -1.
     *   Note: any thread may call this routine, but they will block until finished.
     *   For worker threads this may eventually lead to a deadlock.
     * @param N_work            The number of work items
     * @param ids               Array of work items to wait for
     */
    inline int wait_any(size_t N_work, const thread_id_t *ids);


    /*!
     * \brief   Function to wait until any of the given work items have finished their work
     * \details This is the function waits for any of the given work items to finish. 
     *   If successful it returns the index of a finished work item (the index in the array ids).
     *   If unseccessful it will return -1.
     *   Note: any thread may call this routine, but they will block until finished.
     *   For worker threads this may eventually lead to a deadlock.
     * @param ids               Vector of work items to wait for
     */
    inline int wait_any(const std::vector<thread_id_t>& ids) const;


    /*!
     * \brief   Function to wait until all of the given work items have finished their work
     * \details This is the function waits for all given of the work items to finish.  It returns 0 if successful.
     *   Note: any thread may call this routine, but they will block until finished.
     *   For worker threads this may eventually lead to a deadlock.
     * @param N_work            The number of work items
     * @param ids               Array of work items to wait for
     */
    inline int wait_all(size_t N_work, const thread_id_t *ids) const;


    /*!
     * \brief   Function to wait until all of the given work items have finished their work
     * \details This is the function waits for all given of the work items to finish.  It returns 0 if successful.
     *   Note: any thread may call this routine, but they will block until finished.
     *   For worker threads this may eventually lead to a deadlock.
     * @param ids               Vector of work items to wait for
     */
    inline int wait_all(const std::vector<thread_id_t>& ids) const;


    /*!
     * \brief   Function to wait until all work items in the thread pool have finished their work
     * \details This function will wait until all work has finished.  
     *   Note: member threads may not call this function.  
     *   Only one non-member thread should call this routine at a time.
     */
    void wait_pool_finished() const;


    /*!
     * \brief   Function to check if the thread pool is valid
     * \details Sometimes it is necessary to work with raw pointers for the thread pool.
     *    If the thread pool is invalid and used, the program will likely fail catastrophically.
     *    This function checks if the thread pool is valid is a relatively safe manner.
     *    If the thread pool is pointing to an invalid memory address, because it has been 
     *    freed, never allocated, or otherwise corrupted, this function will return false.
     * @param tpool         Pointer to the ThreadPool to check
     */
    static bool is_valid( const ThreadPool* tpool );


    /*!
     * \brief   Function to enable/disable OS warning messages
     * \details Some of the functions such as setting/getting the thread affinities
     *      are not supported on all platforms.  This function controls the behavior
     *      of these functions on systems where they are not supported.  The default
     *      behavior is to print a warning message.  Other options include ignoring 
     *      the messages (the functions will return empty sets), or throwing an exception.
     *      Note: this is a global property and will affect all thread pools in an application.
     * @param behavior      The behavior of OS specific messages/errors
     *                      0: Print a warning message
     *                      1: Ignore the messages
     *                      2: Throw an error
     */
    static void set_OS_warnings( int behavior=0 );


private:

    friend class ThreadPoolData;

    // Convience typedefs
    #ifdef USE_WINDOWS
        typedef HANDLE wait_type;
    #elif defined(USE_LINUX) || defined(USE_MAC)
        typedef pthread_cond_t* wait_type;
    #else
        #error Unknown OS
    #endif    


private:
    ///// Member data structures

    // Structure to store properties for each work item (linked list)
    struct queue_list_struct {
        short int position;             // Position of the work item in the list
        short int prev;                 // Next item in the list
        short int next;                 // Next item in the list
        queue_list_struct(): position(-1), prev(-1), next(-1) {}
        inline void reset() volatile { prev=-1; next=-1; }
        inline void reset() { prev=-1; next=-1; }
        private:
            queue_list_struct( const queue_list_struct& );
            queue_list_struct& operator=( const queue_list_struct& );
    };

    // Structure to store a pool of wait events (thread safe)
    struct wait_pool_struct {
        wait_pool_struct( );
        ~wait_pool_struct( );
        void push( wait_type event );
        wait_type pop();
        private:
            volatile unsigned int d_count;
            volatile unsigned int d_size;
            volatile wait_type *d_pool;
            #ifdef USE_WINDOWS
                CRITICAL_SECTION *d_lock;
            #elif defined(USE_LINUX) || defined(USE_MAC)
                pthread_mutex_t *d_lock;
            #else
                #error Unknown OS
            #endif
            wait_pool_struct& operator=( const wait_pool_struct& );
            wait_pool_struct( const wait_pool_struct& );
    };

    // Structure to store wait events (note: both the constructor and destructor are NOT thread safe and must be blocked)
    struct wait_event_struct {
        int count;                          // The number of work items that must finish before we alert the thread
        size_t ThreadId;                    // Id of the waiting thread
        std::vector<thread_id_t> ids;       // The ids we are waiting on
        wait_type wait_event;               // Handle to a wait event
        wait_event_struct( wait_pool_struct* wait_pool );
        ~wait_event_struct( );
        private:
            wait_pool_struct* d_wait_pool;
            wait_event_struct( );
            wait_event_struct( const wait_event_struct& );
    };


private:
    ///// Member functions

    // Copy constructors ( we do not want the user to be able to copy the thread pool)
    ThreadPool(const ThreadPool&);
    ThreadPool& operator=(const ThreadPool&);

    // Function to initialize the thread pool
    void setNumThreads( int N, const char* affinity, int N_procs, const int* procs );
    void initialize(int N, const char* affinity, int N_procs, const int* procs);
    void check_startup(size_t size0);

    // Function to add an array of work items
    void add_work(size_t N, ThreadPool::WorkItem* work[], const int* priority, ThreadPool::thread_id_t* id);
        
    // Function to get a work item that has finished
    WorkItem* getFinishedWorkItem(ThreadPool::thread_id_t id) const;
        
    // This function provides a wrapper (needed for the threads)
    static void create_new_thread(void *arglist) {
        void **tmp = (void **) arglist;
        ThreadPool *call = reinterpret_cast<ThreadPool*>(tmp[0]);
        int id = static_cast<int>(reinterpret_cast<size_t>(tmp[1]));
        call->tpool_thread(id);
    }

    /* This is the function that controls the individual thread and allows it to do work.
     * Note: this version uses a last in - first out work scheduling.
     * param thread_init - Structure address contining the startup information for the thread */
    void tpool_thread( int id );

    // Some functions/variables used to get/test the unique work ids 
    inline void initialize_id();                    // A simple function to initialize the id (should only be called once)
    inline size_t advance_id();                     // A simple function to advance the return the id and advance (thread-safe)

    // Function to check if the current thread is a member of the thread pool
    inline bool isMemberThread() const;

    // Function to wait for some work items to finish
    int wait_some(size_t N_work, const thread_id_t *ids, size_t N_wait, bool *finished) const;

    // Helper functions to get the next availible item in the work queue
    inline short int get_work_item( );
    static inline short int check_dependecies( const ThreadPool::queue_list_struct *list,
        const thread_id_t *ids, short int index );


private:
    ///// Member data
    // Note: We want to store the variables in a certain order to optimize storage 
    //   and ensure consistent packing / object size
    size_t d_NULL_HEAD;                                 // Null data buffer to check memory bounds
    volatile AtomicOperations::int64_atomic d_id_assign; // An internal variable used to store the current id to assign
    volatile mutable bool d_signal_empty;               // Do we want to send a signal when the queue is empty
    volatile mutable unsigned char d_signal_count;      // Do we want to send a signal when the count drops to zero
    short int d_N_threads;                              // Number of threads
    volatile short int d_num_active;                    // Number of threads that are currently active
    volatile short int d_queue_head;                    // Index to work queue head
    volatile short int d_queue_free;                    // Index to free queue item
    volatile int d_queue_size;                          // Number of items in the work queue
    volatile mutable int d_N_wait;                      // The number of threads waiting
    size_t d_ThreadId[MAX_NUM_THREADS];                 // Unique id for each thread
    volatile uint64 d_active[MAX_NUM_THREADS/64];       // Which threads are currently active
    volatile uint64 d_cancel[MAX_NUM_THREADS/64];       // Which threads should be deleted
    thread_id_t volatile d_queue_ids[MAX_QUEUED];       // List of ids in the work queue
    queue_list_struct volatile d_queue_list[MAX_QUEUED]; // Work queue list
    volatile mutable wait_event_struct* d_wait[MAX_WAIT]; // The wait events to check
    wait_type d_wait_finished;                          // Handle to a wait event that indicates all threads have finished work
    mutable wait_pool_struct wait_pool;                 // Pool of wait events that we can use
    #ifdef USE_WINDOWS
        CRITICAL_SECTION *d_lock_queue;                 // Mutex lock for changing the queue
        HANDLE d_hThread[MAX_NUM_THREADS];              // Handles to the threads
    #elif defined(USE_LINUX) || defined(USE_MAC)
        pthread_mutex_t *d_lock_queue;                  // Mutex lock for changing the queue
        pthread_t d_hThread[MAX_NUM_THREADS];           // Handles to the threads 
        wait_type d_queue_not_empty;                    // Event condition 
    #else
        #error Unknown OS
    #endif
    size_t d_NULL_TAIL;                                 // Null data buffer to check memory bounds
};



// Swap the contents of the two ids
namespace std {
    template<> inline void swap<ThreadPool::thread_id_t>( 
        ThreadPool::thread_id_t& a, ThreadPool::thread_id_t& b )
    { 
        std::swap(a.d_id,b.d_id);
        std::swap(a.d_count,b.d_count);
        std::swap(a.d_work,b.d_work);
    }
}

#include "thread_pool.hpp"


#endif


