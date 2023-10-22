// This file contains the template functions for the thread pool
#ifndef included_ThreadPoolTmpl
#define included_ThreadPoolTmpl
#include "threadpool/thread_pool.h"
#include <functional>
#include <stdexcept>
#include <tuple>



/*! \addtogroup Macros
 *  @{
 */


/*! \def id = TPOOL_ADD_WORK(tpool,function,args,priority)
 *  \brief Add an item to the thread pool
 *  \details This a macro to automatically create and add a work item to the thread pool.
 *  \param tpool        Pointer to the thread pool to use
 *  \param function     Pointer to the function to use
 *  \param args         The arguments to pass to the function in the form (arg1,arg2,...)
 *  \param priority     Optional argument specifying the priority of the work item
 */
#define TPOOL_ADD_WORK2( TPOOL, FUNCTION, ARGS, PRIORITY, ... ) \
    ThreadPool_add_work( TPOOL, PRIORITY, FUNCTION, std::make_tuple ARGS )
#define TPOOL_ADD_WORK( TPOOL, FUNCTION, ... ) TPOOL_ADD_WORK2( TPOOL, FUNCTION, __VA_ARGS__, 0, 0 )


/*! @} */

// \cond HIDDEN_SYMBOLS


// Unpack a tuple and call a function
template<int...>
struct index_tuple {
};
template<int I, typename IndexTuple, typename... Types>
struct make_indexes_impl;
template<int I, int... Indexes, typename T, typename... Types>
struct make_indexes_impl<I, index_tuple<Indexes...>, T, Types...> {
    typedef typename make_indexes_impl<I + 1, index_tuple<Indexes..., I>, Types...>::type type;
};
template<int I, int... Indexes>
struct make_indexes_impl<I, index_tuple<Indexes...>> {
    typedef index_tuple<Indexes...> type;
};
template<typename... Types>
struct make_indexes : make_indexes_impl<0, index_tuple<>, Types...> {
};
template<class Ret, class... Args, int... Indexes>
inline Ret apply_helper(
    std::function<Ret( Args... )> &pf, index_tuple<Indexes...>, std::tuple<Args...> &&tup )
{
    return pf( std::forward<Args>( std::get<Indexes>( tup ) )... );
}
template<class Ret, class... Args>
inline Ret apply( std::function<Ret( Args... )> &pf, const std::tuple<Args...> &tup )
{
    return apply_helper( pf, typename make_indexes<Args...>::type(), std::tuple<Args...>( tup ) );
}
template<class Ret, class... Args>
inline Ret apply( std::function<Ret( Args... )> &pf, std::tuple<Args...> &&tup )
{
    return apply_helper(
        pf, typename make_indexes<Args...>::type(), std::forward<std::tuple<Args...>>( tup ) );
}


// Specialization for no return argument
template<>
class ThreadPool::WorkItemRet<void> : public ThreadPool::WorkItem
{
public:
    virtual void run() override = 0;
    void get_results() {}
    virtual ~WorkItemRet() {}
    virtual bool has_result() const override final { return false; }
};


// Final class for the work item
template<class Ret, class... Args>
class WorkItemFull;
template<class... Args>
class WorkItemFull<void, Args...> final : public ThreadPool::WorkItemRet<void>
{
private:
    std::function<void( Args... )> routine;
    std::tuple<Args...> args;
    WorkItemFull();

public:
    WorkItemFull( std::function<void( Args... )> &&routine2, Args... ts )
        : ThreadPool::WorkItemRet<void>(), routine( std::move( routine2 ) ), args( ts... )
    {
    }
    WorkItemFull( std::function<void( Args... )> &&routine2, std::tuple<Args...> &&ts )
        : ThreadPool::WorkItemRet<void>(), routine( std::move( routine2 ) ), args( ts )
    {
    }
    virtual void run() override { apply( routine, args ); }
    virtual ~WorkItemFull() {}
};
template<class Ret, class... Args>
class WorkItemFull final : public ThreadPool::WorkItemRet<Ret>
{
private:
    std::function<Ret( Args... )> routine;
    std::tuple<Args...> args;
    WorkItemFull();

public:
    WorkItemFull( std::function<Ret( Args... )> &&routine2, Args... ts )
        : ThreadPool::WorkItemRet<Ret>(), routine( std::move( routine2 ) ), args( ts... )
    {
    }
    WorkItemFull( std::function<Ret( Args... )> &&routine2, std::tuple<Args...> &&ts )
        : ThreadPool::WorkItemRet<Ret>(), routine( std::move( routine2 ) ), args( ts )
    {
    }
    virtual void run() override { this->d_result = apply( routine, args ); }
    virtual ~WorkItemFull() {}
};


// Functions to add work to the thread pool
// clang-format off
template<class Ret, class... Args>
inline ThreadPool::thread_id_t ThreadPool_add_work(
    ThreadPool *tpool, int priority, std::function<Ret( Args... )> routine, std::tuple<Args...> &&args )
{
    auto work = new WorkItemFull<Ret, Args...>( routine, std::move( args ) );
    return ThreadPool::add_work( tpool, work, priority );
}
template<class Ret, class... Args>
inline ThreadPool::thread_id_t ThreadPool_add_work(
    ThreadPool *tpool, int priority, Ret ( *routine )( Args... ), std::tuple<Args...> &&args )
{
    auto work = new WorkItemFull<Ret, Args...>( routine, std::move( args ) );
    return ThreadPool::add_work( tpool, work, priority );
}
template<class Ret, class... Args>
inline ThreadPool::thread_id_t ThreadPool_add_work(
    ThreadPool *tpool, int priority, Ret ( *routine )(), std::tuple<std::nullptr_t>&& )
{
    auto work = new WorkItemFull<Ret>( routine );
    return ThreadPool::add_work( tpool, work, priority );
}
template<class Ret, class... Args>
inline ThreadPool::thread_id_t ThreadPool_add_work(
ThreadPool *tpool, int priority, std::function<Ret( Args... )> routine, Args... args )
{
    auto work = new WorkItemFull<Ret, Args...>( routine, std::forward_as_tuple( args... ) );
    return ThreadPool::add_work( tpool, work, priority );
}
template<class Ret, class... Args>
inline ThreadPool::thread_id_t ThreadPool_add_work(
    ThreadPool *tpool, int priority, Ret ( *routine )( Args... ), Args... args )
{
    auto work = new WorkItemFull<Ret, Args...>( routine, std::forward_as_tuple( args... ) );
    return ThreadPool::add_work( tpool, work, priority );
}
template<class Ret>
inline ThreadPool::thread_id_t ThreadPool_add_work(
    ThreadPool *tpool, int priority, Ret ( *routine )(), void * )
{
    auto work = new WorkItemFull<Ret>( routine );
    return ThreadPool::add_work( tpool, work, priority );
}
template<class Ret, class... Args>
inline ThreadPool::WorkItem *ThreadPool::createWork(
    std::function<Ret( Args... )> routine, Args... args )
{
    return new WorkItemFull<Ret, Args...>( routine, std::forward_as_tuple( args... ) );
}
template<class Ret, class... Args>
inline ThreadPool::WorkItem *ThreadPool::createWork( Ret ( *routine )( Args... ), Args... args )
{
    return new WorkItemFull<Ret, Args...>( routine, std::forward_as_tuple( args... ) );
}
template<class Ret, class... Args>
inline ThreadPool::WorkItem *ThreadPool::createWork(
    std::function<Ret( Args... )> routine, std::tuple<Args...> &&args )
{
    return new WorkItemFull<Ret, Args...>( routine, std::move( args ) );
}
template<class Ret, class... Args>
inline ThreadPool::WorkItem *ThreadPool::createWork(
    Ret ( *routine )( Args... ), std::tuple<Args...> &&args )
{
    return new WorkItemFull<Ret, Args...>( routine, std::move( args ) );
}
// clang-format on


/******************************************************************
 * Function to get the returned function value                     *
 ******************************************************************/
// clang-format off
template<class T> inline constexpr T zeroConstructor();
template<> inline constexpr bool zeroConstructor<bool>() { return false; }
template<> inline constexpr char zeroConstructor<char>() { return 0; }
template<> inline constexpr unsigned char zeroConstructor<unsigned char>() { return 0; }
template<> inline constexpr int zeroConstructor<int>() { return 0; }
template<> inline constexpr unsigned int zeroConstructor<unsigned int>() { return 0; }
template<> inline constexpr long zeroConstructor<long>() { return 0; }
template<> inline constexpr unsigned long zeroConstructor<unsigned long>() { return 0; }
template<> inline constexpr float zeroConstructor<float>() { return 0; }
template<> inline constexpr double zeroConstructor<double>() { return 0; }
template<class T> inline constexpr T zeroConstructor() { return T(); }
template<class Ret>
inline Ret ThreadPool::getFunctionRet( const ThreadPool::thread_id_t &id )
{
    auto work = dynamic_cast<WorkItemRet<Ret> *>( getFinishedWorkItem( id ) );
    return work == nullptr ? zeroConstructor<Ret>() : work->get_results();
}
// clang-format on


/******************************************************************
 * Inline functions to wait for the work items to finish           *
 ******************************************************************/
inline void ThreadPool::wait( ThreadPool::thread_id_t id ) const
{
    bool finished;
    int N = wait_some( 1, &id, 1, &finished, 10000000 );
    if ( N != 1 )
        throw std::logic_error( "Failed to wait for id" );
}
inline size_t ThreadPool::wait_any( const std::vector<thread_id_t> &ids ) const
{
    if ( ids.empty() )
        return 0;
    auto finished = new bool[ids.size()];
    int N         = wait_some( ids.size(), &ids[0], 1, finished, 10000000 );
    if ( N < 1 )
        throw std::logic_error( "Failed to wait for any id" );
    for ( size_t i = 0; i < ids.size(); i++ ) {
        if ( finished[i] ) {
            delete[] finished;
            return i;
        }
    }
    throw std::logic_error( "wait_any failed" );
}
inline void ThreadPool::wait_all( const std::vector<thread_id_t> &ids ) const
{
    if ( ids.empty() )
        return;
    auto finished = new bool[ids.size()];
    int N         = wait_some( ids.size(), ids.data(), ids.size(), finished, 10000000 );
    if ( N != (int) ids.size() )
        throw std::logic_error( "Failed to wait for all ids" );
    delete[] finished;
}
inline void ThreadPool::wait_all( const ThreadPool *tpool, const std::vector<thread_id_t> &ids )
{
    if ( tpool )
        return tpool->wait_all( ids );
}
inline std::vector<int> ThreadPool::wait_some(
    int N_wait, const std::vector<thread_id_t> &ids, int max_wait ) const
{
    auto finished  = new bool[ids.size()];
    int N_finished = wait_some( ids.size(), ids.data(), N_wait, finished, max_wait );
    std::vector<int> index( N_finished, -1 );
    for ( size_t i = 0, j = 0; i < ids.size(); i++ ) {
        if ( finished[i] ) {
            index[j] = i;
            j++;
        }
    }
    delete[] finished;
    return index;
}


/******************************************************************
 * Functions to add work items.                                    *
 ******************************************************************/
inline ThreadPool::thread_id_t ThreadPool::add_work( WorkItem *work, int priority )
{
    ThreadPool::thread_id_t id;
    add_work( 1, &work, &priority, &id );
    return id;
}
inline std::vector<ThreadPool::thread_id_t> ThreadPool::add_work(
    const std::vector<ThreadPool::WorkItem *> &work, const std::vector<int> &priority )
{
    size_t N = work.size();
    if ( N == 0 )
        return std::vector<ThreadPool::thread_id_t>();
    if ( priority.size() != N && !priority.empty() )
        throw std::logic_error( "size of work and priority do not match" );
    const int *priority2 = nullptr;
    if ( priority.empty() ) {
        priority2 = new int[N];
        memset( const_cast<int *>( priority2 ), 0, N * sizeof( int ) );
    } else {
        priority2 = &priority[0];
    }
    std::vector<ThreadPool::thread_id_t> ids( N );
    add_work( N, const_cast<ThreadPool::WorkItem **>( &work[0] ), priority2, &ids[0] );
    if ( priority.empty() )
        delete[] priority2;
    return ids;
}
inline ThreadPool::thread_id_t ThreadPool::add_work(
    ThreadPool *tpool, ThreadPool::WorkItem *work, int priority )
{
    ThreadPool::thread_id_t id;
    if ( tpool ) {
        id = tpool->add_work( work, priority );
    } else {
        id.reset( priority, std::rand(), work );
        work->d_state = 2;
        work->run();
        work->d_state = 3;
    }
    return id;
}
inline std::vector<ThreadPool::thread_id_t> ThreadPool::add_work( ThreadPool *tpool,
    const std::vector<ThreadPool::WorkItem *> &work, const std::vector<int> &priority )
{
    if ( tpool ) {
        return tpool->add_work( work, priority );
    } else {
        std::vector<ThreadPool::thread_id_t> ids( work.size() );
        for ( size_t i = 0; i < work.size(); i++ )
            ids[i] = add_work( tpool, work[i], priority[i] );
        return ids;
    }
}


/******************************************************************
 * Class functions to for the thread id                            *
 ******************************************************************/
inline ThreadPool::thread_id_t::thread_id_t()
    : d_id( nullThreadID ), d_count( nullptr ), d_work( nullptr )
{
}
inline ThreadPool::thread_id_t::~thread_id_t() { reset(); }
inline ThreadPool::thread_id_t::thread_id_t( volatile ThreadPool::thread_id_t &&rhs )
    : d_id( std::move( rhs.d_id ) ),
      d_count( std::move( rhs.d_count ) ),
      d_work( std::move( rhs.d_work ) )
{
    rhs.d_count = nullptr;
    rhs.d_work  = nullptr;
    rhs.d_id    = nullThreadID;
}
inline ThreadPool::thread_id_t &ThreadPool::thread_id_t::operator=(
    const ThreadPool::thread_id_t &rhs ) volatile
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return const_cast<ThreadPool::thread_id_t &>( *this );
    this->reset();
    d_id    = rhs.d_id;
    d_count = rhs.d_count;
    d_work  = rhs.d_work;
    if ( d_count != nullptr )
        AtomicOperations::atomic_increment( d_count );
    return const_cast<ThreadPool::thread_id_t &>( *this );
}
inline ThreadPool::thread_id_t &ThreadPool::thread_id_t::operator=(
    volatile ThreadPool::thread_id_t &&rhs ) volatile
{
    std::swap( d_id, rhs.d_id );
    std::swap( d_work, rhs.d_work );
    std::swap( d_count, rhs.d_count );
    return const_cast<ThreadPool::thread_id_t &>( *this );
}
inline ThreadPool::thread_id_t::thread_id_t( const volatile ThreadPool::thread_id_t &rhs )
    : d_id( rhs.d_id ), d_count( rhs.d_count ), d_work( rhs.d_work )
{
    if ( d_count != nullptr )
        AtomicOperations::atomic_increment( d_count );
}
#if !defined( WIN32 ) && !defined( _WIN32 ) && !defined( WIN64 ) && !defined( _WIN64 )
inline ThreadPool::thread_id_t::thread_id_t( const thread_id_t &rhs )
    : d_id( rhs.d_id ), d_count( rhs.d_count ), d_work( rhs.d_work )
{
    if ( d_count != nullptr )
        AtomicOperations::atomic_increment( d_count );
}
inline ThreadPool::thread_id_t &ThreadPool::thread_id_t::operator=( ThreadPool::thread_id_t &&rhs )
{
    std::swap( d_id, rhs.d_id );
    std::swap( d_work, rhs.d_work );
    std::swap( d_count, rhs.d_count );
    return const_cast<ThreadPool::thread_id_t &>( *this );
}
inline ThreadPool::thread_id_t &ThreadPool::thread_id_t::operator=(
    const ThreadPool::thread_id_t &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return const_cast<ThreadPool::thread_id_t &>( *this );
    this->reset();
    d_id    = rhs.d_id;
    d_count = rhs.d_count;
    d_work  = rhs.d_work;
    if ( d_count != nullptr )
        AtomicOperations::atomic_increment( d_count );
    return const_cast<ThreadPool::thread_id_t &>( *this );
}
inline ThreadPool::thread_id_t &ThreadPool::thread_id_t::operator=(
    const volatile ThreadPool::thread_id_t &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return const_cast<ThreadPool::thread_id_t &>( *this );
    this->reset();
    d_id    = rhs.d_id;
    d_count = rhs.d_count;
    d_work  = rhs.d_work;
    if ( d_count != nullptr )
        AtomicOperations::atomic_increment( d_count );
    return const_cast<ThreadPool::thread_id_t &>( *this );
}
inline ThreadPool::thread_id_t &ThreadPool::thread_id_t::operator=(
    const volatile ThreadPool::thread_id_t &rhs ) volatile
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return const_cast<ThreadPool::thread_id_t &>( *this );
    this->reset();
    d_id    = rhs.d_id;
    d_count = rhs.d_count;
    d_work  = rhs.d_work;
    if ( d_count != nullptr )
        AtomicOperations::atomic_increment( d_count );
    return const_cast<ThreadPool::thread_id_t &>( *this );
}
#endif
inline void ThreadPool::thread_id_t::reset() volatile
{
    if ( d_count != nullptr ) {
        int count = AtomicOperations::atomic_decrement( d_count );
        if ( count == 0 ) {
            WorkItem *tmp = reinterpret_cast<ThreadPool::WorkItem *>( d_work );
            delete tmp;
        }
    }
    d_id    = nullThreadID;
    d_count = nullptr;
    d_work  = nullptr;
}
inline void ThreadPool::thread_id_t::reset()
{
    if ( d_count != nullptr ) {
        int count = AtomicOperations::atomic_decrement( d_count );
        if ( count == 0 ) {
            WorkItem *tmp = reinterpret_cast<ThreadPool::WorkItem *>( d_work );
            delete tmp;
        }
    }
    d_id    = nullThreadID;
    d_count = nullptr;
    d_work  = nullptr;
}
inline uint64_t ThreadPool::thread_id_t::createId( int priority, uint64_t local_id )
{
    if ( priority < -127 || priority > 127 || local_id > maxThreadID )
        throw std::logic_error( "Invalid priority or local id" );
    uint64_t id = priority + 128;
    id          = ( id << 56 ) + local_id;
    return id;
}
inline void ThreadPool::thread_id_t::reset( int priority, uint64_t local_id, void *work )
{
    if ( d_count != nullptr ) {
        int count = AtomicOperations::atomic_decrement( d_count );
        if ( count == 0 ) {
            WorkItem *tmp = reinterpret_cast<ThreadPool::WorkItem *>( d_work );
            delete tmp;
        }
    }
    // Create the id
    d_id = createId( priority, local_id );
    // Create the work and counter
    d_count = nullptr;
    d_work  = work;
    if ( d_work != nullptr ) {
        d_count  = &( reinterpret_cast<WorkItem *>( work )->d_count );
        *d_count = 1;
    }
}
inline uint64_t ThreadPool::thread_id_t::getLocalID() const
{
    if ( d_id == nullThreadID )
        return ~( (uint64_t) 0 );
    uint64_t tmp = d_id & 0x00FFFFFFFFFFFFFF;
    return static_cast<size_t>( tmp );
}
inline int ThreadPool::thread_id_t::getPriority() const
{
    if ( d_id == nullThreadID )
        return -128;
    uint64_t tmp = d_id >> 56;
    return static_cast<int>( tmp ) - 128;
}
inline void ThreadPool::thread_id_t::setPriority( int priority )
{
    if ( d_id == nullThreadID )
        return;
    d_id = createId( priority, getLocalID() );
}
inline bool ThreadPool::thread_id_t::started() const
{
    return d_id == nullThreadID ? true : reinterpret_cast<WorkItem *>( d_work )->d_state >= 2;
}
inline bool ThreadPool::thread_id_t::finished() const
{
    return d_id == nullThreadID ? true : reinterpret_cast<WorkItem *>( d_work )->d_state == 3;
}
inline bool ThreadPool::thread_id_t::ready() const
{
    bool ready = true;
    if ( !isNull() ) {
        auto tmp = work();
        for ( size_t i = 0; i < tmp->d_N_ids; i++ )
            ready = ready && tmp->d_ids[i].finished();
    }
    return ready;
}


/******************************************************************
 * This function checks if the id is valid                         *
 ******************************************************************/
inline bool ThreadPool::isValid( const ThreadPool::thread_id_t &id ) const
{
    uint64_t local_id = id.getLocalID();
    uint64_t next_id  = d_id_assign - 1;
    return local_id != 0 && id.initialized() && local_id <= thread_id_t::maxThreadID &&
           local_id > next_id;
}


/******************************************************************
 * Function to get the thread number                               *
 * (-1 if it is not a member thread)                               *
 ******************************************************************/
inline int ThreadPool::getThreadNumber() const
{
    std::thread::id id = std::this_thread::get_id();
    for ( int i = 0; i < d_N_threads; i++ ) {
        if ( id == d_threadId[i] )
            return i;
    }
    return -1;
}


// \endcond


#endif
