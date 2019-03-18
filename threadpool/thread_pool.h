// Copyright © 2004 Mark Berrill. All Rights Reserved. This work is distributed with permission,
// but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.
#ifndef included_AtomicModelThreadPool
#define included_AtomicModelThreadPool

#include <condition_variable>
#include <functional>
#include <map>
#include <mutex>
#include <stdarg.h>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <thread>
#include <typeinfo>
#include <vector>


#include "threadpool/atomic_helpers.h"
#include "threadpool/atomic_list.h"


// clang-format off


/** \class ThreadPool
 *
 * \brief This is a concrete class that provides for a basic thread pool.
 * \details This class implements a basic thread pool that can be used for a wide variety of
 * applications.
 * An example call usage is provided below.  The ability to return a value is provided.  Note that
 * there
 * is a small overhead to using this functionality. <BR>
 * <pre>Example: <BR>
 *    Existing function call:
 *       double x = myfun_1(a,b);
 *       double y = myfun_2(c,d); <BR>
 *    Threaded call (processing in parallel):
 *       thread_id_t ids[2];
 *       ids[0] = TPOOL_ADD_WORK( tpool, myfun_1, (a,b) );
 *       ids[1] = TPOOL_ADD_WORK( tpool, myfun_2, (c,d) );
 *       wait_all(2,ids);
 *       double x = getFunctionRet(ids[0]);
 *       double y = getFunctionRet(ids[1]); <BR>
 *   </pre>
 */
class ThreadPool
{
public:
    ///// Set some global properties
    constexpr static uint16_t MAX_THREADS = 128; // The maximum number of threads (must be a multiple of 64)
    constexpr static uint16_t MAX_WAIT = 16;     // The maximum number of active waits at any given time

public:
    ///// Member classes
    class WorkItem;

    /** \class thread_id_t
     *
     * \brief This a class to hold the work item id
     * \details This class hold the id of the work item that is being processed by the thread pool.
     *      It is created when a work item is added to the thread pool and is used by various
     * routines within the thread pool.
     */
    class thread_id_t
    {
    public:
        // nullID definitins
        static constexpr uint64_t nullThreadID = 0x0FFFFFFFFFFFFFFF;
        static constexpr uint64_t maxThreadID  = 0x00FFFFFFFFFFFFFD;
        //! Empty constructor
        inline thread_id_t();
        //! Destructor
        inline ~thread_id_t();
        //! Copy constructors
        inline thread_id_t( const volatile thread_id_t &rhs );
        inline thread_id_t( volatile thread_id_t &&rhs );
        inline thread_id_t &operator=( const thread_id_t &rhs ) volatile;
        inline thread_id_t &operator=( volatile thread_id_t &&rhs ) volatile;
#if !defined( WIN32 ) && !defined( _WIN32 ) && !defined( WIN64 ) && !defined( _WIN64 )
        inline thread_id_t( const thread_id_t &rhs );
        inline thread_id_t &operator=( thread_id_t &&rhs );
        inline thread_id_t &operator=( const thread_id_t &rhs );
        inline thread_id_t &operator=( const volatile thread_id_t &rhs );
        inline thread_id_t &operator=( const volatile thread_id_t &rhs ) volatile;
#endif
        // Overload key operators
        inline bool operator==( const thread_id_t &rhs ) const { return !((d_id^rhs.d_id)&nullThreadID); }
        inline bool operator!=( const thread_id_t &rhs ) const { return (d_id^rhs.d_id)&nullThreadID; }
        inline bool operator>=( const thread_id_t &rhs ) const { return d_id >= rhs.d_id; }
        inline bool operator<=( const thread_id_t &rhs ) const { return d_id <= rhs.d_id; }
        inline bool operator>(  const thread_id_t &rhs ) const { return d_id  > rhs.d_id; }
        inline bool operator<(  const thread_id_t &rhs ) const { return d_id  < rhs.d_id; }
        inline bool operator==( const volatile thread_id_t &rhs ) const volatile { return !((d_id^rhs.d_id)&nullThreadID); }
        inline bool operator!=( const volatile thread_id_t &rhs ) const volatile { return (d_id^rhs.d_id)&nullThreadID; }
        inline bool operator>=( const volatile thread_id_t &rhs ) const volatile { return d_id >= rhs.d_id; }
        inline bool operator<=( const volatile thread_id_t &rhs ) const volatile { return d_id <= rhs.d_id; }
        inline bool operator>(  const volatile thread_id_t &rhs ) const volatile { return d_id  > rhs.d_id; }
        inline bool operator<(  const volatile thread_id_t &rhs ) const volatile { return d_id  < rhs.d_id; }
        //! Reset the id back to a NULL id
        inline void reset() volatile;
        inline void reset();
        //! Check if the work has started (will return true if it has started or finished)
        inline bool started() const;
        //! Check if the work has finished
        inline bool finished() const;
        //! swap with rhs
        inline void swap( thread_id_t &rhs )
        {
            std::swap( this->d_id, rhs.d_id );
            std::swap( this->d_count, rhs.d_count );
            std::swap( this->d_work, rhs.d_work );
        }
        //! Check if thread id is null
        inline bool isNull( ) const { return d_id==nullThreadID; }
        //! Check if thread id is null
        inline WorkItem* getWork( ) const { return reinterpret_cast<WorkItem *>( d_work ); }

    private:
        // Reset the internal data to the given values
        inline void reset( int priority, uint64_t local_id, void *work );
        static inline uint64_t createId( int priority, uint64_t local_id );
        // Get the local id
        inline uint64_t getLocalID() const;
        // Get the priority
        inline int getPriority() const;
        // Increase the priority
        inline void setPriority( int priority );
        // Check if the id is initialized
        inline bool initialized() const volatile { return d_id != 0x0FFFFFFFFFFFFFFF; }
        // Get a pointer to the work structure
        inline WorkItem* work() const { return reinterpret_cast<WorkItem *>( d_work ); }
        // Is the id ready to process
        inline bool ready() const;
        // Friends
        friend class ThreadPool;
        // Data
        uint64_t d_id;                                    // 64-bit data to store id
        volatile AtomicOperations::int32_atomic *d_count; // Reference count
        void *d_work;                                     // Pointer to the work item
    };


    //! Base class for the work item (users should derive from WorkItemRet)
    class WorkItem
    {
    public:
        //! Function to run the routine
        virtual void run() = 0;
        //! Will the routine return a result
        virtual bool has_result() const = 0;
        //! Empty deconstructor
        virtual ~WorkItem()
        {
            delete[] d_ids;
            d_ids   = nullptr;
            d_N_ids = 0;
            d_size  = 0;
        }
        //! Get the number of work ids that this work item depends on
        inline size_t get_N_dependencies() const { return d_N_ids; }
        //! Return the list of work ids that we depend on
        std::vector<ThreadPool::thread_id_t> get_dependencies() const;
        /*!
         * \brief Add a work item to the list of dependencies
         * \param id    Id of the work item to add
         */
        void add_dependency( const ThreadPool::thread_id_t &id ) { add_dependencies( 1, &id ); }
        /*!
         * \brief Add a list of work item to the list of dependencies
         * \param ids   Ids of the work item to add
         */
        inline void add_dependencies( const std::vector<ThreadPool::thread_id_t> &ids )
        {
            if ( !ids.empty() )
                add_dependencies( ids.size(), &ids[0] );
        }
        /*!
         * \brief Add a list of work item to the list of dependencies
         *    Note: this function is thread-safe for the threadpool and does not need blocking.
         * \param N     Number of items to add
         * \param ids   Ids of the work item to add
         */
        void add_dependencies( size_t N, const ThreadPool::thread_id_t *ids );

    protected:
        friend class ThreadPool;
        inline WorkItem():
              d_state( 0 ),
              d_N_ids( 0 ),
              d_size( 0 ),
              d_count( 0 ),
              d_ids( nullptr )
        {
        }

    private:
        WorkItem( const WorkItem & );            // Private copy constructor
        WorkItem &operator=( const WorkItem & ); // Private assignment operator
        volatile char d_state;                   // Current state (0: not added to threadpool, 1: queued, 2: started, 3: finished)
        uint16_t d_N_ids;                        // Number of dependencies
        uint16_t d_size;                         // Size of d_ids
        AtomicOperations::int32_atomic d_count;  // Count used by a thread_id
        thread_id_t *d_ids;                      // Pointer to id list
        // Friends
        friend class ThreadPool::thread_id_t;
    };


    /*!
     * \brief   Class to define a work item returning a variable
     * \details This is the class that defines a work item to be processed.  Users may derive their
     * own
     * class and add work using the add_work routine, or can use the TPOOL_ADD_WORK macro.
     * Note: this class is templated on the return argument type and may be a void type.
     */
    template <typename return_type>
    class WorkItemRet : public ThreadPool::WorkItem
    {
    public:
        //! Run the work item
        virtual void run() override = 0;
        //! Will the routine return a result
        virtual bool has_result() const override final { return !std::is_same<return_type,void>::value; }
        //! Return the results
        return_type get_results() const { return d_result; }
        //! Virtual destructor
        virtual ~WorkItemRet() {}
    protected:
        return_type d_result;
    protected:
        inline WorkItemRet() : d_result( return_type() ) { }
    private:
        WorkItemRet( const WorkItemRet & );            // Private copy constructor
        WorkItemRet &operator=( const WorkItemRet & ); // Private assignment operator
    };


public:
    ///// Member functions

    /*!
     *  Constructor that initialize the thread pool with N threads
     * @param N                 The desired number of worker threads
     * @param affinity          The affinity scheduler to use:
     *                          none - Let the OS handle the affinities (default)
     *                          independent - Give each thread an independent set of processors
     * @param procs             The processors to use (defaults to the process affinitiy list)
     * @param queueSize         The maximum number of items in the queue before forcing a wait
     */
    ThreadPool( const int N = 0, const std::string &affinity = "none",
        const std::vector<int> &procs = std::vector<int>(), int queueSize = 1024 );


    //! Destructor
    ~ThreadPool();


    //! Function to return the number of processors availible
    static int getNumberOfProcessors();


    //! Function to return the processor number that the current thread is running on
    static int getCurrentProcessor();


    //! Function to return the affinity of the current process
    static std::vector<int> getProcessAffinity();


    //! Function to set the affinity of the current process
    static void setProcessAffinity( const std::vector<int>& procs );


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
    static void setThreadAffinity( const std::vector<int>& procs );


    /*!
     *  Set the given thread to have the specified affinity
     *  @param thread   The index of the thread
     *  @param procs    The processors to use
     */
    void setThreadAffinity( int thread, const std::vector<int>& procs ) const;


    //! Function to return the number of threads in the thread pool
    inline int getNumThreads() const { return d_N_threads; }


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
    void setNumThreads( const int N, const std::string &affinity = "none",
        const std::vector<int> &procs = std::vector<int>() );


    /*!
     * \brief   Function to set the maximum wait time
     * \details  This function sets the maximum time the thread pool will
     *    wait before warning about a possible hung thread.
     *    Default is to wait 10 minutes.
     * @param time              The number of seconds to wait (seconds)
     */
    inline void setMaxWaitTimeDebug( const int time ) { d_max_wait_time = time; }


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
    inline bool isValid( const thread_id_t &id ) const;


    /*!
     * \brief    Function to check if the work item has finished processing
     * \details  This function checks if the work item has finished processing.
     * @param id                The id of the work item
     */
    inline bool isFinished( thread_id_t& id ) const { return id.finished(); }


    /*!
     * \brief   Function to get the returned function value
     * \details This is the function returns the value that was returned from the working function.
     *   If the work item has not finished or was not found it will return 0.
     * @param id                The id of the work item
     */
    template <class return_type>
    static inline return_type getFunctionRet( const thread_id_t &id );


    /*!
     * \brief   Function to create a work item
     * \details This function creates a work item that can be added to the queue
     * @param routine           Function to call from the thread pool
     * @param args              Function arguments to pass
     */
    template<class Ret, class... Args>
    static inline WorkItem* createWork( std::function<Ret(Args...)> routine, std::tuple<Args...> &&args );


    /*!
     * \brief   Function to create a work item
     * \details This function creates a work item that can be added to the queue
     * @param routine           Function to call from the thread pool
     * @param args              Function arguments to pass
     */
    template<class Ret, class... Args>
    static inline WorkItem* createWork( Ret( *routine )( Args... ), std::tuple<Args...> &&args );


    /*!
     * \brief   Function to create a work item
     * \details This function creates a work item that can be added to the queue
     * @param routine           Function to call from the thread pool
     * @param args              Function arguments to pass
     */
    template<class Ret, class... Args>
    static inline WorkItem* createWork( std::function<Ret(Args...)> routine, Args... args );


    /*!
     * \brief   Function to create a work item
     * \details This function creates a work item that can be added to the queue
     * @param routine           Function to call from the thread pool
     * @param args              Function arguments to pass
     */
    template <class Ret, class... Args>
    static inline WorkItem* createWork( Ret( *routine )( Args... ), Args... args );


    /*!
     * \brief   Function to add a work item
     * \details This function adds a work item to the queue
     *   Note: any thread may call this routine.
     * @param work              Pointer to the work item to add
     *                          Note that the threadpool will automatically destroy the item when
     * finished
     * @param priority          A value indicating the priority of the work item (0-default)
     */
    inline thread_id_t add_work( ThreadPool::WorkItem *work, int priority = 0 );


    /*!
     * \brief   Function to add multiple work items
     * \details This function adds multiple work item to the queue
     *   Note: any thread may call this routine.
     * @param work              Vector of pointers to the work items to add
     *                          Note that the threadpool will automatically destroy the item when
     * finished
     * @param priority          Vector of values indicating the priority of the work items
     */
    inline std::vector<thread_id_t> add_work( const std::vector<ThreadPool::WorkItem *> &work,
        const std::vector<int> &priority = std::vector<int>() );


    /*!
     * \brief   Function to wait until a specific work item has finished
     * \details This is the function waits for a specific work item to finished.
     *   Note: any thread may call this routine, but they will block until finished.
     *   For worker threads this may eventually lead to a deadlock.
     * @param id                The work item to wait for
     */
    inline void wait( thread_id_t id ) const;


    /*!
     * \brief   Function to wait until any of the given work items have finished their work
     * \details This is the function waits for any of the given work items to finish.
     *   If successful it returns the index of a finished work item (the index in the array ids).
     *   Note: any thread may call this routine, but they will block until finished.
     *   For worker threads this may eventually lead to a deadlock.
     * @param ids               Vector of work items to wait for
     */
    inline size_t wait_any( const std::vector<thread_id_t> &ids ) const;


    /*!
     * \brief   Function to wait until all of the given work items have finished their work
     * \details This is the function waits for all given of the work items to finish.
     *   Note: any thread may call this routine, but they will block until finished.
     *   For worker threads this may eventually lead to a deadlock.
     * @param ids               Vector of work items to wait for
     */
    inline void wait_all( const std::vector<thread_id_t> &ids ) const;


    /*!
     * \brief   Function to wait until some of the given work items have finished their work
     * \details This is the function waits for some of the given work items to finish.
     *   If successful it returns the indicies of the finished work items (the index in the array ids).
     *   Note: any thread may call this routine, but they will block until finished.
     *   For worker threads this may eventually lead to a deadlock.
     * @param N_wait            Number of work items to wait for
     * @param ids               Vector of work items to wait for
     * @param max_wait          Maximum time to wait (seconds)
     */
    inline std::vector<int> wait_some( int N_wait, const std::vector<thread_id_t> &ids, int max_wait = 10000000 ) const;


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
    static bool is_valid( const ThreadPool *tpool );


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
    static void set_OS_warnings( int behavior = 0 );


    //! Return the number of items queued
    int N_queued( ) const { return d_queue_list.size(); }


    //! Set the error handler for threads
    void setErrorHandler( std::function<void(const std::string&)> fun );


public: // Static interface

    /*!
     * \brief   Function to return the number of work threads
     * \details This function returns the number of threads in the thread pool,
     *    or 0 if the thread pool is empty or does not exist
     * @param tpool         Threadpool to add work to (may be null)
     */
    static inline int numThreads( const ThreadPool* tpool ) { return tpool ? tpool->getNumThreads() : 0; }

    /*!
     * \brief   Function to add a work item
     * \details This function adds a work item to the queue
     *   Note: any thread may call this routine.
     * @param tpool         Threadpool to add work to (may be null)
     * @param work          Pointer to the work item to add
     *                      Note that the threadpool will automatically destroy the item when finished
     * @param priority      A value indicating the priority of the work item (0-default)
     */
    static inline thread_id_t add_work( ThreadPool* tpool, ThreadPool::WorkItem *work, int priority = 0 );


    /*!
     * \brief   Function to add multiple work items
     * \details This function adds multiple work item to the queue
     *   Note: any thread may call this routine.
     * @param tpool         Threadpool to add work to (may be null)
     * @param work          Vector of pointers to the work items to add
     *                      Note that the threadpool will automatically destroy the item when finished
     * @param priority      Vector of values indicating the priority of the work items
     */
    static inline std::vector<thread_id_t> add_work( ThreadPool* tpool, const std::vector<ThreadPool::WorkItem *> &work,
        const std::vector<int> &priority = std::vector<int>() );


    /*!
     * \brief   Function to wait until all of the given work items have finished their work
     * \details This is the function waits for all given of the work items to finish.
     *   Note: any thread may call this routine, but they will block until finished.
     *   For worker threads this may eventually lead to a deadlock.
     * @param tpool         Threadpool containing work (must match call to add_work)
     * @param ids           Vector of work items to wait for
     */
    static inline void wait_all( const ThreadPool* tpool, const std::vector<thread_id_t> &ids );


    /*!
     * \brief   Function to wait until all work items in the thread pool have finished their work
     * \details This function will wait until all work has finished.
     *   Note: member threads may not call this function.
     *   Only one non-member thread should call this routine at a time.
     * @param tpool         Threadpool containing work (must match call to add_work)
     */
    static inline void wait_pool_finished( const ThreadPool* tpool ) { if ( tpool ) { tpool->wait_pool_finished(); } }


private:
    ///// Member data structures

   
    // Implimentation of condition_variable which does not require a lock
    class condition_variable
    {
      public:
        condition_variable() { }
        ~condition_variable() { }
        inline void wait() const { std::unique_lock<std::mutex> lock(d_mutex); d_cv.wait(lock); }
        inline void wait_for( double seconds ) const
        {
            std::unique_lock<std::mutex> lock(d_mutex);
            if ( seconds < 4e-6 )
                d_cv.wait_for(lock,std::chrono::nanoseconds(static_cast<int>(1e9*seconds)));
            else if ( seconds < 4e-3 )
                d_cv.wait_for(lock,std::chrono::microseconds(static_cast<int>(1e6*seconds)));
            else if ( seconds < 4 )
                d_cv.wait_for(lock,std::chrono::milliseconds(static_cast<int>(1e3*seconds)));
            else
                d_cv.wait_for(lock,std::chrono::seconds(static_cast<int>(seconds)));
        }
        inline void notify_one() const { d_cv.notify_one(); }
        inline void notify_all() const { d_cv.notify_all(); }
      private:
        mutable std::condition_variable d_cv;
        mutable std::mutex d_mutex;
    };


    // Structure to wait on multiple ids
    // Note: this is thread safe without blocking as long as it is added to the wait list
    //    before calling wait
    class wait_ids_struct {
      public:
        wait_ids_struct() = delete;
        wait_ids_struct( const wait_ids_struct& ) = delete;
        wait_ids_struct& operator=( const wait_ids_struct & ) = delete;
        wait_ids_struct( size_t N, const ThreadPool::thread_id_t *ids, size_t N_wait,
            AtomicOperations::pool<condition_variable,128>& cv_pool, int N_wait_list, volatile wait_ids_struct **list );
        ~wait_ids_struct( );
        void id_finished( const ThreadPool::thread_id_t& id ) const;
        bool wait_for( double total_time, double recheck_time );
      private:
        mutable int d_wait;                     // The number of work items that must finish before we alert the thread
        mutable int d_N;                        // The number of ids we are waiting on
        mutable thread_id_t *d_ids;             // The ids we are waiting on
        AtomicOperations::pool<condition_variable,128>& d_cv_pool;
        condition_variable *d_wait_event;       // Handle to a wait event
        volatile mutable bool *d_finished;      // Has each id finished
        volatile mutable wait_ids_struct **d_ptr;
      private:
        inline bool check();
    };


private:
    ///// Member functions

    // Copy constructors ( we do not want the user to be able to copy the thread pool)
    ThreadPool( const ThreadPool & );
    ThreadPool &operator=( const ThreadPool & );

    // Function to check the startup
    void check_startup( );

    // Function to add an array of work items
    void add_work(
        size_t N, ThreadPool::WorkItem *work[], const int *priority, ThreadPool::thread_id_t *id );
    inline void add_work( const ThreadPool::thread_id_t& id );

    // Function to get a work item that has finished
    static inline WorkItem *getFinishedWorkItem( const ThreadPool::thread_id_t& id )
    {
        return id.finished() ? id.work():nullptr;
    }

    // This function provides a wrapper (needed for the threads)
    static inline void create_new_thread( ThreadPool *tpool, int id )
    {
        tpool->tpool_thread( id );
    }

    /* This is the function that controls the individual thread and allows it to do work.
     * Note: this version uses a last in - first out work scheduling.
     * param thread_init - Structure address contining the startup information for the thread */
    void tpool_thread( int id );

    // Function to check if the current thread is a member of the thread pool
    inline bool isMemberThread() const { return getThreadNumber()>=0; }

    // Function to wait for some work items to finish
    int wait_some( size_t N_work, const thread_id_t *ids, size_t N_wait, bool *finished, int max_wait ) const;
    
    // Check if we are waiting too long and pring debug info
    void print_wait_warning( ) const;


private:
    ///// Member data

    // Typedefs
    typedef volatile AtomicOperations::int32_atomic vint32_t;
    typedef volatile AtomicOperations::int64_atomic vint64_t;
    typedef volatile wait_ids_struct vwait_t;
    typedef AtomicOperations::pool<condition_variable,128> cond_t;
    typedef AtomicList<thread_id_t,std::greater<thread_id_t>> queue_type;

    // Internal data
    uint32_t d_NULL_HEAD;                 // Null data buffer to check memory bounds
    volatile mutable bool d_signal_empty; // Do we want to send a signal when the queue is empty
    uint16_t d_N_threads;                 // Number of threads
    int d_max_wait_time;                  // The maximum time in a wait command before printing a warning message
    vint32_t d_signal_count;              // Signal count
    vint32_t d_num_active;                // Number of threads that are currently active
    vint64_t d_id_assign;                 // An internal variable used to store the current id to assign
    vint64_t d_active[MAX_THREADS/64];    // Which threads are currently active
    vint64_t d_cancel[MAX_THREADS/64];    // Which threads should be deleted
    vint64_t d_N_added;                   // Number of items added to the work queue
    vint64_t d_N_started;                 // Number of items started
    vint64_t d_N_finished;                // Number of items finished
    mutable vwait_t *d_wait[MAX_WAIT];    // The wait events to check
    mutable wait_ids_struct *d_wait_last; // A cached copy of the last completed wait event (in case a thread still has a reference)
    condition_variable d_wait_finished;   // Condition variable to signal when all work is finished
    condition_variable d_wait_work;       // Condition variable to signal when there is new work
    mutable cond_t d_cond_pool;           // Condition pool
    std::thread d_thread[MAX_THREADS];    // Handles to the threads
    std::thread::id d_threadId[MAX_THREADS]; // Unique id for each thread
    queue_type d_queue_list;              // The work queue
    std::function<void(const std::string&)> d_errorHandler; // Error handler
    uint32_t d_NULL_TAIL;                 // Null data buffer to check memory bounds
};


#include "threadpool/thread_pool.hpp"


// clang-format on
#endif
