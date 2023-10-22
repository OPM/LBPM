#include "IO/PIO.h"
#include "common/MPI.h"
#include "common/Utilities.h"

#include <cstring>
#include <fstream>
#include <string>


namespace IO {


static ParallelStreamBuffer pout_buffer;
static ParallelStreamBuffer perr_buffer;
static ParallelStreamBuffer plog_buffer;


std::ostream pout( &pout_buffer );
std::ostream perr( &perr_buffer );
std::ostream plog( &plog_buffer );


/****************************************************************************
 *  Functions to control logging                                             *
 ****************************************************************************/
std::ofstream *global_filestream = NULL;
static void shutdownFilestream()
{
    if ( global_filestream != NULL ) {
        global_filestream->flush();
        global_filestream->close();
        delete global_filestream;
        global_filestream = NULL;
    }
}
void Utilities::logOnlyNodeZero( const std::string &filename )
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif
    if ( rank == 0 )
        logAllNodes( filename, true );
}
void Utilities::logAllNodes( const std::string &filename, bool singleStream )
{
    if ( singleStream )
        ERROR( "Not implimented yet" );

    // If the filestream was open, then close it and reset streams
    shutdownFilestream();

    // Open the log stream and redirect output
    std::string full_filename = filename;
    if ( !singleStream ) {
        int rank = 0;
#ifdef USE_MPI
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif
        char tmp[100];
        sprintf( tmp, ".%04i", rank );
        full_filename += std::string( tmp );
    }
    global_filestream = new std::ofstream( full_filename.c_str() );

    if ( !( *global_filestream ) ) {
        delete global_filestream;
        global_filestream = NULL;
        perr << "PIO: Could not open log file ``" << full_filename << "''\n";
    } else {
        pout_buffer.setOutputStream( global_filestream );
        pout_buffer.setOutputStream( &std::cout );
        perr_buffer.setOutputStream( global_filestream );
        perr_buffer.setOutputStream( &std::cerr );
        plog_buffer.setOutputStream( global_filestream );
    }
}


/****************************************************************************
 *  ParallelStreamBuffer class                                               *
 ****************************************************************************/
void Utilities::stopLogging()
{
    pout_buffer.reset();
    perr_buffer.reset();
    plog_buffer.reset();
    shutdownFilestream();
    delete global_filestream;
    global_filestream = NULL;
}


/****************************************************************************
 *  ParallelStreamBuffer class                                               *
 ****************************************************************************/
ParallelStreamBuffer::ParallelStreamBuffer()
    : d_rank( 0 ), d_size( 0 ), d_buffer_size( 0 ), d_buffer( NULL )
{
}
ParallelStreamBuffer::~ParallelStreamBuffer() { delete[] d_buffer; }
void ParallelStreamBuffer::setOutputStream( std::ostream *stream ) { d_stream.push_back( stream ); }
int ParallelStreamBuffer::sync()
{
    for ( size_t i = 0; i < d_stream.size(); i++ ) {
        std::ostream &stream = *d_stream[i];
        stream << d_buffer;
    }
    d_size = 0;
    memset( d_buffer, 0, d_buffer_size );
    return 0;
}
void ParallelStreamBuffer::reserve( size_t size )
{
    if ( size > d_buffer_size ) {
        if ( d_buffer_size == 0 ) {
            d_buffer_size = 1024;
            d_buffer      = new char[d_buffer_size];
            memset( d_buffer, 0, d_buffer_size );
        }
        while ( size > d_buffer_size ) {
            char *tmp = d_buffer;
            d_buffer_size *= 2;
            d_buffer = new char[d_buffer_size];
            memset( d_buffer, 0, d_buffer_size );
            memcpy( d_buffer, tmp, d_size );
            delete[] tmp;
        }
    }
}
std::streamsize ParallelStreamBuffer::xsputn( const char *text, std::streamsize n )
{
    reserve( d_size + n );
    memcpy( &d_buffer[d_size], text, n );
    d_size += n;
    if ( text[n - 1] == 0 || text[n - 1] == 10 ) {
        sync();
    }
    return n;
}
int ParallelStreamBuffer::overflow( int ch )
{
    reserve( d_size + 1 );
    d_buffer[d_size] = ch;
    d_size++;
    if ( ch == 0 || ch == 10 ) {
        sync();
    }
    return std::char_traits<char>::to_int_type( ch );
}
int ParallelStreamBuffer::underflow() { return -1; }
void ParallelStreamBuffer::reset()
{
    sync();
    d_stream.clear();
    delete[] d_buffer;
    d_buffer      = NULL;
    d_buffer_size = 0;
}


} // namespace IO
