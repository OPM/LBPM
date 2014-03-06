#include <iostream>
#include "mpi.h"


int main (int argc, char **argv)
{
	int rank,nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    for (int i=0; i<nprocs; i++) {
        if ( rank==i )
            printf("%i of %i: Hello world\n",rank,nprocs);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // set the error code
    // Note: the error code should be consistent across all processors
    int error = 0;
    
    // Finished
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return error; 
}
