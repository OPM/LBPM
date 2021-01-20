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
#ifndef included_UnitTest
#define included_UnitTest

#include <mutex>
#include <sstream>
#include <string>
#include <vector>
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
class UnitTest
{
public:
    //! Constructor
    UnitTest();

    //! Destructor
    virtual ~UnitTest();

    //! Indicate a passed test (thread-safe)
    virtual void passes( const std::string &in );

    //! Indicate a failed test (thread-safe)
    virtual void failure( const std::string &in );

    //! Indicate an expected failed test (thread-safe)
    virtual void expected_failure( const std::string &in );

    //! Return the number of passed tests locally
    virtual size_t NumPassLocal() const { return pass_messages.size(); }

    //! Return the number of failed tests locally
    virtual size_t NumFailLocal() const { return fail_messages.size(); }

    //! Return the number of expected failed tests locally
    virtual size_t NumExpectedFailLocal() const { return expected_fail_messages.size(); }

    //! Return the number of passed tests locally
    virtual size_t NumPassGlobal() const;

    //! Return the number of failed tests locally
    virtual size_t NumFailGlobal() const;

    //! Return the number of expected failed tests locally
    virtual size_t NumExpectedFailGlobal() const;

    //! Return the rank of the current processor
    int getRank() const;

    //! Return the number of processors
    int getSize() const;

    /*!
     * Print a report of the passed and failed tests.
     * Note: This is a blocking call that all processors must execute together.
     * Note: Only rank 0 will print the messages (this is necessary as other ranks may not be able
     * to print correctly).
     * @param level     Optional integer specifying the level of reporting (default: 1)
     *                  0: Report the number of tests passed, failed, and expected failures.
     *                  1: Report the number of passed tests (if <=20) or the number passed
     *                     otherwise, report all failures, report the number of expected
     *                     failed tests (if <=50) or the number passed otherwise.
     *                  2: Report all passed, failed, and expected failed tests.
     */
    virtual void report( const int level = 1 ) const;

    //! Clear the messages
    void reset();

protected:
    std::vector<std::string> pass_messages;
    std::vector<std::string> fail_messages;
    std::vector<std::string> expected_fail_messages;
    mutable std::mutex mutex;
#ifdef USE_MPI
    MPI_Comm comm;
#endif

private:
    // Make the copy constructor private
    UnitTest( const UnitTest & ) {}

    // Function to pack the messages into a single data stream and send to the given processor
    // Note: This function does not return until the message stream has been sent
    void pack_message_stream(
        const std::vector<std::string> &messages, const int rank, const int tag ) const;

    // Function to unpack the messages from a single data stream
    // Note: This function does not return until the message stream has been received
    std::vector<std::string> unpack_message_stream( const int rank, const int tag ) const;

    // Helper functions
    inline void barrier() const;
    inline std::vector<int> allGather( int value ) const;
    inline std::vector<std::vector<std::string>> gatherMessages(
        const std::vector<std::string> &local_messages, int tag ) const;
};


#endif
