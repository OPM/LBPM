#ifndef included_UnitTest
#define included_UnitTest

#include "common/MPI.h"

#include <mutex>
#include <string>
#include <vector>


/*!
 * @brief Class UnitTest is simple utility for running unit tests.
 * It provides basic routines for tracing success or failure of tests,
 * and reporting the results.
 * \par Code Sample:
 * \code
    try {
        std::cout << "Testing tstOne" << std::endl;
        tstOne(&ut);
        ut.passes("Test XXX passed");
    } catch( ... ) {
        ut.failure("An unknown exception was thrown");
    }
    ut.report();
    return ut.NumFail();
 * \endcode

 */
class UnitTest final
{
public:
    //! Constructor
    UnitTest();

    //! Destructor
    ~UnitTest();

    // Copy constructor
    UnitTest( const UnitTest & ) = delete;

    // Assignment operator
    UnitTest &operator=( const UnitTest & ) = delete;

    //! Indicate a passed test (thread-safe)
    void passes( std::string in );

    //! Indicate a failed test (thread-safe)
    void failure( std::string in );

    //! Indicate an expected failed test (thread-safe)
    void expected_failure( std::string in );

    //! Return the number of passed tests locally
    inline size_t NumPassLocal() const { return d_pass.size(); }

    //! Return the number of failed tests locally
    inline size_t NumFailLocal() const { return d_fail.size(); }

    //! Return the number of expected failed tests locally
    inline size_t NumExpectedFailLocal() const { return d_expected.size(); }

    //! Return the number of passed tests locally
    size_t NumPassGlobal() const;

    //! Return the number of failed tests locally
    size_t NumFailGlobal() const;

    //! Return the number of expected failed tests locally
    size_t NumExpectedFailGlobal() const;

    /*!
     * Print a report of the passed and failed tests.
     * Note: This is a blocking call that all processors must execute together.
     * Note: Only rank 0 will print the messages (this is necessary as other ranks may not be able
     * to print correctly).
     * @param level     Optional integer specifying the level of reporting (default: 1)
     *                  0: Report the number of tests passed, failed, and expected failures.
     *                  1: Report the passed tests (if <=20) or number passed,
     *                     Report all failures,
     *                     Report the expected failed tests (if <=50) or the number passed.
     *                  2: Report all passed, failed, and expected failed tests.
     */
    void report( const int level = 1 ) const;

    //! Clear the messages
    void reset();

    //! Make the unit test operator verbose?
    void verbose( bool verbose = true ) { d_verbose = verbose; }

private:
    std::vector<std::string> d_pass;
    std::vector<std::string> d_fail;
    std::vector<std::string> d_expected;
    bool d_verbose;
    mutable std::mutex d_mutex;
    Utilities::MPI d_comm;

private:
    // Function to pack the messages into a single data stream and send to the given processor
    // Note: This function does not return until the message stream has been sent
    void pack_message_stream(
        const std::vector<std::string> &messages, const int rank, const int tag ) const;

    // Function to unpack the messages from a single data stream
    // Note: This function does not return until the message stream has been received
    std::vector<std::string> unpack_message_stream( const int rank, const int tag ) const;

    // Gather the messages
    inline std::vector<std::vector<std::string>> gatherMessages(
        const std::vector<std::string> &local_messages, int tag ) const;
};


#endif
