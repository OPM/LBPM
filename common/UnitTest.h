#ifndef included_UnitTest
#define included_UnitTest

#include <sstream>
#include <vector>
#include <string>
#ifdef USE_MPI
    #include "mpi.h"
#endif


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
class UnitTest {
public:

    //! Constructor
    UnitTest();

    //! Indicate a passed test
    virtual void passes (const std::string &in) { pass_messages.push_back(in); }

    //! Indicate a failed test
    virtual void failure (const std::string &in) { fail_messages.push_back(in); }

    //! Indicate an expected failed test
    virtual void expected_failure (const std::string &in) { expected_fail_messages.push_back(in); }

    //! Return the number of passed tests locally
    virtual size_t NumPassLocal () { return pass_messages.size(); }

    //! Return the number of failed tests locally
    virtual size_t NumFailLocal () { return fail_messages.size(); }

    //! Return the number of expected failed tests locally
    virtual size_t NumExpectedFailLocal () { return expected_fail_messages.size(); }

    //! Return the number of passed tests locally
    virtual size_t NumPassGlobal ();

    //! Return the number of failed tests locally
    virtual size_t NumFailGlobal ();

    //! Return the number of expected failed tests locally
    virtual size_t NumExpectedFailGlobal ();

    //! Return the rank of the current processor
    int getRank ();

    //! Return the number of processors
    int getSize ();

    /*!
     * Print a report of the passed and failed tests.
     * Note: This is a blocking call that all processors must execute together.
     * Note: Only rank 0 will print the messages (this is necessary as other ranks may not be able to print correctly).
     * @param level     Optional integer specifying the level of reporting (default: 1)
     *                  0: Report the number of tests passed, failed, and expected failures.
     *                  1: Report the number of passed tests (if <=20) or the number passed otherwise,
     *                     report all failures,
     *                     report the number of expected failed tests (if <=50) or the number passed otherwise.
     *                  2: Report all passed, failed, and expected failed tests.
     */
     virtual void report(const int level=1);

protected:
    std::vector<std::string> pass_messages;
    std::vector<std::string> fail_messages;
    std::vector<std::string> expected_fail_messages;
    #ifdef USE_MPI
        MPI_Comm comm;
    #endif

private:
    // Make the copy constructor private
    UnitTest(const UnitTest& p) {}

    // Function to pack the messages into a single data stream and send to the given processor
    // Note: This function does not return until the message stream has been sent
    void pack_message_stream(const std::vector<std::string>& messages, const int rank, const int tag);

    // Function to unpack the messages from a single data stream
    // Note: This function does not return until the message stream has been received
    std::vector<std::string> unpack_message_stream(const int rank, const int tag);

};


#endif
