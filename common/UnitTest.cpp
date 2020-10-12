/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
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
UnitTest::UnitTest()
{
#ifdef USE_MPI
    comm = MPI_COMM_WORLD;
#endif
}
UnitTest::~UnitTest() { reset(); }
void UnitTest::reset()
{
    mutex.lock();
    // Clear the data forcing a reallocation
    std::vector<std::string>().swap( pass_messages );
    std::vector<std::string>().swap( fail_messages );
    std::vector<std::string>().swap( expected_fail_messages );
    mutex.unlock();
}


/********************************************************************
 *  Add a pass, fail, expected failure message in a thread-safe way  *
 ********************************************************************/
void UnitTest::passes( const std::string &in )
{
    mutex.lock();
    pass_messages.push_back( in );
    mutex.unlock();
}
void UnitTest::failure( const std::string &in )
{
    mutex.lock();
    fail_messages.push_back( in );
    mutex.unlock();
}
void UnitTest::expected_failure( const std::string &in )
{
    mutex.lock();
    expected_fail_messages.push_back( in );
    mutex.unlock();
}


/********************************************************************
 *  Print a global report                                            *
 *  Note: only rank 0 will print, all messages will be aggregated    *
 ********************************************************************/
inline std::vector<int> UnitTest::allGather( int value ) const
{
    int size = getSize();
    std::vector<int> data( size, value );
#ifdef USE_MPI
    if ( size > 1 )
        MPI_Allgather( &value, 1, MPI_INT, data.data(), 1, MPI_INT, comm );
#endif
    return data;
}
inline void UnitTest::barrier() const
{
#ifdef USE_MPI
    if ( getSize() > 1 )
        MPI_Barrier( comm );
#endif
}
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
    mutex.lock();
    int size = getSize();
    int rank = getRank();
    // Broadcast the print level from rank 0
    int level = level0;
#ifdef USE_MPI
    if ( getSize() > 1 )
        MPI_Bcast( &level, 1, MPI_INT, 0, comm );
#endif
    if ( level < 0 || level > 2 )
        ERROR( "Invalid print level" );
    // Perform a global all gather to get the number of failures per processor
    auto N_pass             = allGather( pass_messages.size() );
    auto N_fail             = allGather( fail_messages.size() );
    auto N_expected_fail    = allGather( expected_fail_messages.size() );
    int N_pass_tot          = 0;
    int N_fail_tot          = 0;
    int N_expected_fail_tot = 0;
    for ( int i = 0; i < size; i++ ) {
        N_pass_tot += N_pass[i];
        N_fail_tot += N_fail[i];
        N_expected_fail_tot += N_expected_fail[i];
    }
    // Send all messages to rank 0 (if needed)
    std::vector<std::vector<std::string>> pass_messages_rank( size );
    std::vector<std::vector<std::string>> fail_messages_rank( size );
    std::vector<std::vector<std::string>> expected_fail_rank( size );
    // Get the pass messages
    if ( ( level == 1 && N_pass_tot <= 20 ) || level == 2 )
        pass_messages_rank = UnitTest::gatherMessages( pass_messages, 1 );
    // Get the fail messages
    if ( level == 1 || level == 2 )
        fail_messages_rank = UnitTest::gatherMessages( fail_messages, 2 );
    // Get the expected_fail messages
    if ( ( level == 1 && N_expected_fail_tot <= 50 ) || level == 2 )
        expected_fail_rank = UnitTest::gatherMessages( expected_fail_messages, 2 );
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
        if ( level == 0 || ( level == 1 && N_expected_fail_tot > 50 ) ) {
            // We want to print a summary
            if ( size > 8 ) {
                // Print 1 summary for all processors
                printp( "     %i tests expected failed (use report level 2 for more detail)\n",
                    N_expected_fail_tot );
            } else {
                // Print a summary for each processor
                for ( int i = 0; i < size; i++ )
                    printp( "     %i tests expected failed (proc %i) (use report level 2 for more "
                            "detail)\n",
                        N_expected_fail[i], i );
            }
        } else {
            // We want to print all messages
            for ( int i = 0; i < size; i++ )
                ASSERT( (int) expected_fail_rank[i].size() == N_expected_fail[i] );
            print_messages( expected_fail_rank );
        }
        pout << std::endl;
    }
    // Add a barrier to synchronize all processors (rank 0 is much slower)
    barrier();
    Utilities::sleep_ms( 10 ); // Need a brief pause to allow any printing to finish
    mutex.unlock();
}


/********************************************************************
 *  Gather the messages to rank 0                                    *
 ********************************************************************/
std::vector<std::vector<std::string>> UnitTest::gatherMessages(
    const std::vector<std::string> &local_messages, int tag ) const
{
    const int rank = getRank();
    const int size = getSize();
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
#ifdef USE_MPI
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
    MPI_Request request;
    MPI_Isend( data, size_data, MPI_CHAR, rank, tag, comm, &request );
    // Wait for the communication to send and free the temporary memory
    MPI_Status status;
    MPI_Wait( &request, &status );
    delete[] data;
    delete[] msg_size;
#else
    NULL_USE( messages );
    NULL_USE( rank );
    NULL_USE( tag );
#endif
}


/********************************************************************
 *  Receive and unpack a message stream                              *
 ********************************************************************/
std::vector<std::string> UnitTest::unpack_message_stream( const int rank, const int tag ) const
{
#ifdef USE_MPI
    // Probe the message to get the message size
    MPI_Status status;
    MPI_Probe( rank, tag, comm, &status );
    int size_data = -1;
    MPI_Get_count( &status, MPI_BYTE, &size_data );
    ASSERT( size_data >= 0 );
    // Allocate memory to receive the data
    auto *data = new char[size_data];
    // receive the data (using a non-blocking receive)
    MPI_Request request;
    MPI_Irecv( data, size_data, MPI_CHAR, rank, tag, comm, &request );
    // Wait for the communication to be received
    MPI_Wait( &request, &status );
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
    delete[] data;
    return messages;
#else
    NULL_USE( rank );
    NULL_USE( tag );
    return std::vector<std::string>();
#endif
}


/********************************************************************
 *  Other functions                                                  *
 ********************************************************************/
int UnitTest::getRank() const
{
    int rank = 0;
#ifdef USE_MPI
    int flag = 0;
    MPI_Initialized( &flag );
    if ( flag )
        MPI_Comm_rank( comm, &rank );
#endif
    return rank;
}
int UnitTest::getSize() const
{
    int size = 1;
#ifdef USE_MPI
    int flag = 0;
    MPI_Initialized( &flag );
    if ( flag )
        MPI_Comm_size( comm, &size );
#endif
    return size;
}
size_t UnitTest::NumPassGlobal() const
{
    size_t num = pass_messages.size();
#ifdef USE_MPI
    if ( getSize() > 1 ) {
        auto send = static_cast<int>( num );
        int sum   = 0;
        MPI_Allreduce( &send, &sum, 1, MPI_INT, MPI_SUM, comm );
        num = static_cast<size_t>( sum );
    }
#endif
    return num;
}
size_t UnitTest::NumFailGlobal() const
{
    size_t num = fail_messages.size();
#ifdef USE_MPI
    if ( getSize() > 1 ) {
        auto send = static_cast<int>( num );
        int sum   = 0;
        MPI_Allreduce( &send, &sum, 1, MPI_INT, MPI_SUM, comm );
        num = static_cast<size_t>( sum );
    }
#endif
    return num;
}
size_t UnitTest::NumExpectedFailGlobal() const
{
    size_t num = expected_fail_messages.size();
#ifdef USE_MPI
    if ( getSize() > 1 ) {
        auto send = static_cast<int>( num );
        int sum   = 0;
        MPI_Allreduce( &send, &sum, 1, MPI_INT, MPI_SUM, comm );
        num = static_cast<size_t>( sum );
    }
#endif
    return num;
}
