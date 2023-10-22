#ifndef included_PIO
#define included_PIO

#include <iostream>
#include <vector>


namespace IO {


/*!
 * Parallel output stream pout writes to the standard output from node zero
 * only.  Output from other nodes is ignored.  If logging is enabled, then
 * output is mirrored to the log stream, as well.
 */
extern std::ostream pout;

/*!
 * Parallel output stream perr writes to the standard error from all nodes.
 * Output is prepended with the processor number.
 */
extern std::ostream perr;

/*!
 * Parallel output stream plog writes output to the log file.  When logging
 * from multiple processors, the processor number is appended to the filename.
 */
extern std::ostream plog;

/*!
 * Parallel output printp pout writes to the standard output from node zero
 * only.  Output from other nodes is ignored.  If logging is enabled, then
 * output is mirrored to the log stream, as well.
 * The format matches the format for printf
 */
inline int printp( const char *format, ... );


/*!
 * Class ParallelBuffer is a simple I/O stream utility that
 * intercepts output from an ostream and redirects the output as necessary
 * for parallel I/O.  This class defines a stream buffer class for an
 * ostream class.
 */
class ParallelStreamBuffer : public std::streambuf
{
public:
    /*!
     * Create a parallel buffer class.  The object will require further
     * initialization to set up the I/O streams and prefix string.
     */
    ParallelStreamBuffer();

    /*!
     * Set the output file stream (multiple output streams are supported)
     * @param stream    Output stream
     */
    void setOutputStream( std::ostream *stream );

    /*!
     * The destructor simply deallocates any internal data
     * buffers.  It does not modify the output streams.
     */
    virtual ~ParallelStreamBuffer();

    /*!
     * Synchronize the parallel buffer (called from streambuf).
     */
    virtual int sync();

    /**
     * Write the specified number of characters into the output stream (called
     * from streambuf).
     */
    virtual std::streamsize xsputn( const char *text, std::streamsize n );

    /*!
     * Write an overflow character into the parallel buffer (called from
     * streambuf).
     */
    virtual int overflow( int ch );

    /*!
     * Read an overflow character from the parallel buffer (called from
     * streambuf).  This is not implemented.  It is needed by the
     * MSVC++ stream implementation.
     */
    virtual int underflow();

    /*!
     * Clear the internal buffer's memory
     */
    virtual void reset();

private:
    int d_rank;
    size_t d_size;
    size_t d_buffer_size;
    char *d_buffer;
    std::vector<std::ostream *> d_stream;
    inline void reserve( size_t size );
};


namespace Utilities {

/*!
 * Log messages for node zero only to the specified filename.  All output
 * to pout, perr, and plog on node zero will go to the log file.
 */
void logOnlyNodeZero( const std::string &filename );

/*!
 * Log messages from all nodes.  The diagnostic data for processor XXXXX
 * will be sent to a file with the name filename.XXXXX, where filename is
 * the function argument.
 */
void logAllNodes( const std::string &filename, bool singleStream = false );

/*!
 * Stop logging messages, flush buffers, and reset memory.
 */
void stopLogging();


} // namespace Utilities


} // namespace IO


#include "IO/PIO.hpp"

#endif
