#ifndef included_StackTraceErrorHandlers
#define included_StackTraceErrorHandlers


#include "StackTrace/StackTrace.h"
#include "common/MPI.h"

#include <functional>


namespace StackTrace
{


    /*!
     * Set the error handler
     * @param[in] abort     Function to terminate the program: abort(msg,type)
     */
    void setErrorHandler( std::function<void( const StackTrace::abort_error& )> abort );

    //! Clear the error handler
    void clearErrorHandler();


    //! Set an error handler for MPI
    void setMPIErrorHandler( MPI_Comm comm );

    //! Clear an error handler for MPI
    void clearMPIErrorHandler( MPI_Comm comm );


    //! Initialize globalCallStack functionallity
    void globalCallStackInitialize( MPI_Comm comm );

    //! Clean up globalCallStack functionallity
    void globalCallStackFinalize();


} // namespace StackTrace

#endif
