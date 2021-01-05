#include <iostream>
#include "common/MPI.h"
#include "common/Utilities.h"


int main (int argc, char **argv)
{
    Utilities::startup( argc, argv );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();

    for (int i=0; i<nprocs; i++) {
        if ( rank==i )
            printf("%i of %i: Hello world\n",rank,nprocs);
        comm.barrier();
    }

    // Create a memory leak for valgrind to find
    if ( nprocs==1 ) {
        double *x = new double[1];
        ASSERT(x!=NULL);
    }

    // set the error code
    // Note: the error code should be consistent across all processors
    int error = 0;
    
    // Finished
    Utilities::shutdown();
    return error; 
}
