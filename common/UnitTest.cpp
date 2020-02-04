#include "common/UnitTest.h"
#include "common/Utilities.h"
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


#define pout std::cout
#define printp printf


/********************************************************************
 *  Constructor/Destructor                                           *
 ********************************************************************/
UnitTest::UnitTest() : d_verbose( false ), d_comm( MPI_COMM_SELF )
{
    if ( Utilities::MPI::MPI_active() )
        d_comm = MPI_COMM_WORLD;
}
UnitTest::~UnitTest() { reset(); }
void UnitTest::reset()
{
    d_mutex.lock();
    // Clear the data forcing a reallocation
    std::vector<std::string>().swap( d_pass );
    std::vector<std::string>().swap( d_fail );
    std::vector<std::string>().swap( d_expected );
    d_mutex.unlock();
}


/********************************************************************
 *  Add a pass, fail, expected failure message in a thread-safe way  *
 ********************************************************************/
void UnitTest::passes( std::string in )
{
    d_mutex.lock();
    if ( d_verbose )
        printf( "UnitTest: %i passes: %s\n", d_comm.getRank(), in.data() );
    d_pass.emplace_back( std::move( in ) );
    d_mutex.unlock();
}
void UnitTest::failure( std::string in )
{
    d_mutex.lock();
    if ( d_verbose )
        printf( "UnitTest: %i failed: %s\n", d_comm.getRank(), in.data() );
    d_fail.emplace_back( std::move( in ) );
    d_mutex.unlock();
}
void UnitTest::expected_failure( std::string in )
{
    d_mutex.lock();
    if ( d_verbose )
        printf( "UnitTest: %i expected_failure: %s\n", d_comm.getRank(), in.data() );
    d_expected.emplace_back( std::move( in ) );
    d_mutex.unlock();
}


/********************************************************************
 *  Print a global report                                            *
 *  Note: only rank 0 will print, all messages will be aggregated    *
 ********************************************************************/
static inline void print_messages( const std::vector<std::vector<std::string>> &messages )
{
    if ( messages.size() > 1 ) {
        for ( size_t i = 0; i < messages.size(); i++ ) {
            if ( !messages[i].empty() ) {
                printp( "     Proccessor %i:\n", static_cast<int>( i ) );
                for ( const auto &j : messages[i] )
                    pout << "        " << j << std::endl;
            }
        }
    } else {
        for ( const auto &j : messages[0] )
            pout << "    " << j << std::endl;
    }
}
void UnitTest::report( const int level0 ) const
{
    d_mutex.lock();
    int size = d_comm.getSize();
    int rank = d_comm.getRank();
    // Give all processors a chance to print any remaining messages
    d_comm.barrier();
    Utilities::sleep_ms( 10 );
    // Broadcast the print level from rank 0
    int level = d_comm.bcast( level0, 0 );
    if ( level < 0 || level > 2 )
        ERROR( "Invalid print level" );
    // Perform a global all gather to get the number of failures per processor
    auto N_pass        = d_comm.allGather<int>( d_pass.size() );
    auto N_fail        = d_comm.allGather<int>( d_fail.size() );
    auto N_expected    = d_comm.allGather<int>( d_expected.size() );
    int N_pass_tot     = 0;
    int N_fail_tot     = 0;
    int N_expected_tot = 0;
    for ( int i = 0; i < size; i++ ) {
        N_pass_tot += N_pass[i];
        N_fail_tot += N_fail[i];
        N_expected_tot += N_expected[i];
    }
    // Send all messages to rank 0 (if needed)
    std::vector<std::vector<std::string>> pass_messages_rank( size );
    std::vector<std::vector<std::string>> fail_messages_rank( size );
    std::vector<std::vector<std::string>> expected_fail_rank( size );
    // Get the pass messages
    if ( ( level == 1 && N_pass_tot <= 20 ) || level == 2 )
        pass_messages_rank = UnitTest::gatherMessages( d_pass, 1 );
    // Get the fail messages
    if ( level == 1 || level == 2 )
        fail_messages_rank = UnitTest::gatherMessages( d_fail, 2 );
    // Get the expected_fail messages
    if ( ( level == 1 && N_expected_tot <= 50 ) || level == 2 )
        expected_fail_rank = UnitTest::gatherMessages( d_expected, 2 );
    // Print the results of all messages (only rank 0 will print)
    if ( rank == 0 ) {
        pout << std::endl;
        // Print the passed tests
        pout << "Tests passed" << std::endl;
        if ( level == 0 || ( level == 1 && N_pass_tot > 20 ) ) {
            // We want to print a summary
            if ( size > 8 ) {
                // Print 1 summary for all processors
                printp( "     %i tests passed (use report level 2 for more detail)\n", N_pass_tot );
            } else {
                // Print a summary for each processor
                for ( int i = 0; i < size; i++ )
                    printp( "     %i tests passed (proc %i) (use report level 2 for more detail)\n",
                        N_pass[i], i );
            }
        } else {
            // We want to print all messages
            for ( int i = 0; i < size; i++ )
                ASSERT( (int) pass_messages_rank[i].size() == N_pass[i] );
            print_messages( pass_messages_rank );
        }
        pout << std::endl;
        // Print the tests that failed
        pout << "Tests failed" << std::endl;
        if ( level == 0 ) {
            // We want to print a summary
            if ( size > 8 ) {
                // Print 1 summary for all processors
                printp( "     %i tests failed (use report level 2 for more detail)\n", N_fail_tot );
            } else {
                // Print a summary for each processor
                for ( int i = 0; i < size; i++ )
                    printp( "     %i tests failed (proc %i) (use report level 2 for more detail)\n",
                        N_fail[i], i );
            }
        } else {
            // We want to print all messages
            for ( int i = 0; i < size; i++ )
                ASSERT( (int) fail_messages_rank[i].size() == N_fail[i] );
            print_messages( fail_messages_rank );
        }
        pout << std::endl;
        // Print the tests that expected failed
        pout << "Tests expected failed" << std::endl;
        if ( level == 0 || ( level == 1 && N_expected_tot > 50 ) ) {
            // We want to print a summary
            if ( size > 8 ) {
                // Print 1 summary for all processors
                printp( "     %i tests expected failed (use report level 2 for more detail)\n",
                    N_expected_tot );
            } else {
                // Print a summary for each processor
                for ( int i = 0; i < size; i++ )
                    printp( "     %i tests expected failed (proc %i) (use report level 2 for more "
                            "detail)\n",
                        N_expected[i], i );
            }
        } else {
            // We want to print all messages
            for ( int i = 0; i < size; i++ )
                ASSERT( (int) expected_fail_rank[i].size() == N_expected[i] );
            print_messages( expected_fail_rank );
        }
        pout << std::endl;
    }
    // Add a barrier to synchronize all processors (rank 0 is much slower)
    d_comm.barrier();
    Utilities::sleep_ms( 10 ); // Need a brief pause to allow any printing to finish
    d_mutex.unlock();
}


/********************************************************************
 *  Gather the messages to rank 0                                    *
 ********************************************************************/
std::vector<std::vector<std::string>> UnitTest::gatherMessages(
    const std::vector<std::string> &local_messages, int tag ) const
{
    const int rank = d_comm.getRank();
    const int size = d_comm.getSize();
    std::vector<std::vector<std::string>> messages( size );
    if ( rank == 0 ) {
        // Rank 0 should receive all messages
        for ( int i = 0; i < size; i++ ) {
            if ( i == 0 )
                messages[i] = local_messages;
            else
                messages[i] = unpack_message_stream( i, tag );
        }
    } else {
        // All other ranks send their message (use non-blocking communication)
        pack_message_stream( local_messages, 0, tag );
    }
    return messages;
}


/********************************************************************
 *  Pack and send the given messages                                 *
 ********************************************************************/
void UnitTest::pack_message_stream(
    const std::vector<std::string> &messages, const int rank, const int tag ) const
{
    // Get the size of the messages
    auto N_messages  = (int) messages.size();
    auto *msg_size   = new int[N_messages];
    int msg_size_tot = 0;
    for ( int i = 0; i < N_messages; i++ ) {
        msg_size[i] = (int) messages[i].size();
        msg_size_tot += msg_size[i];
    }
    // Allocate space for the message stream
    size_t size_data = ( N_messages + 1 ) * sizeof( int ) + msg_size_tot;
    auto *data       = new char[size_data];
    // Pack the message stream
    memcpy( data, &N_messages, sizeof( int ) );
    memcpy( &data[sizeof( int )], msg_size, N_messages * sizeof( int ) );
    size_t k = ( N_messages + 1 ) * sizeof( int );
    for ( int i = 0; i < N_messages; i++ ) {
        messages[i].copy( &data[k], msg_size[i] );
        k += msg_size[i];
    }
    // Send the message stream (using a non-blocking send)
    auto request = d_comm.Isend( data, size_data, rank, tag );
    // Wait for the communication to send and free the temporary memory
    d_comm.wait( request );
    delete[] data;
    delete[] msg_size;
}


/********************************************************************
 *  Receive and unpack a message stream                              *
 ********************************************************************/
std::vector<std::string> UnitTest::unpack_message_stream( const int rank, const int tag ) const
{
    // Probe the message to get the message size
    int size_data = d_comm.probe( rank, tag );
    ASSERT( size_data >= 0 );
    // Allocate memory to receive the data
    auto *data = new char[size_data];
    // receive the data (using a non-blocking receive)
    auto request = d_comm.Irecv( data, size_data, rank, tag );
    // Wait for the communication to be received
    d_comm.wait( request );
    // Unpack the message stream
    int N_messages = 0;
    memcpy( &N_messages, data, sizeof( int ) );
    if ( N_messages == 0 ) {
        delete[] data;
        return std::vector<std::string>();
    }
    std::vector<int> msg_size( N_messages );
    std::vector<std::string> messages( N_messages );
    memcpy( msg_size.data(), &data[sizeof( int )], N_messages * sizeof( int ) );
    int k = ( N_messages + 1 ) * sizeof( int );
    for ( int i = 0; i < N_messages; i++ ) {
        messages[i] = std::string( &data[k], msg_size[i] );
        k += msg_size[i];
    }
    // Delete the temporary memory
    delete[] data;
    return messages;
}


/********************************************************************
 *  Other functions                                                  *
 ********************************************************************/
size_t UnitTest::NumPassGlobal() const { return d_comm.sumReduce( d_pass.size() ); }
size_t UnitTest::NumFailGlobal() const { return d_comm.sumReduce( d_fail.size() ); }
size_t UnitTest::NumExpectedFailGlobal() const { return d_comm.sumReduce( d_expected.size() ); }

