#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "common/Utilities.h"
#include "models/FreeLeeModel.h"

inline void Initialize_Mask(ScaLBL_FreeLeeModel &LeeModel){
	// initialize a bubble
	int i,j,k,n;
	int rank = LeeModel.Mask->rank();
	int Nx = LeeModel.Mask->Nx;
	int Ny = LeeModel.Mask->Ny;
	int Nz = LeeModel.Mask->Nz;
	if (rank == 0) printf(" initialize mask...\n");

	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;
                LeeModel.Mask->id[n]=1;
                LeeModel.id[n] = LeeModel.Mask->id[n];
			}
		}
	}
}

inline void Initialize_DummyPhaseField(ScaLBL_FreeLeeModel &LeeModel, double ax, double ay, double az){
	// initialize a bubble
	int i,j,k,n;
	int rank = LeeModel.Mask->rank();
	int Nx = LeeModel.Mask->Nx;
	int Ny = LeeModel.Mask->Ny;
	int Nz = LeeModel.Mask->Nz;
	if (rank == 0) printf("Setting up dummy phase field with gradient {x,y,z} = {%f , %f , %f}...\n",ax,ay,az);
	
	double * Dummy;
	int Nh = (Nx+2)*(Ny+2)*(Nz+2);
	Dummy = new double [(Nx+2)*(Ny+2)*(Nz+2)];
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;
                LeeModel.Mask->id[n]=1;
                LeeModel.id[n] = LeeModel.Mask->id[n];
                int nh = (k+1)*(Nx+2)*(Ny+2) + (j+1)*(Nx+2) + i+1;
                Dummy[nh] = ax*double(i) + ay*double(j) + az*double(k);
			}
		}
	}
	ScaLBL_CopyToDevice(LeeModel.Phi, Dummy, sizeof(double)*Nh);
	
	LeeModel.MGTest();
}

int main( int argc, char **argv )
{

    // Initialize
    Utilities::startup( argc, argv );

    // Load the input database
    auto db = std::make_shared<Database>( argv[1] );

    { // Limit scope so variables that contain communicators will free before MPI_Finialize

        Utilities::MPI comm( MPI_COMM_WORLD );
        int rank   = comm.getRank();
        int nprocs = comm.getSize();

        if ( rank == 0 ) {
            printf( "********************************************************\n" );
            printf( "Running Mixed Gradient Test	\n" );
            printf( "********************************************************\n" );
        }
        // Initialize compute device
        int device = ScaLBL_SetDevice( rank );
        NULL_USE( device );
        ScaLBL_DeviceBarrier();
        comm.barrier();

        PROFILE_ENABLE( 1 );
        // PROFILE_ENABLE_TRACE();
        // PROFILE_ENABLE_MEMORY();
        PROFILE_SYNCHRONIZE();
        PROFILE_START( "Main" );
        Utilities::setErrorHandlers();

        auto filename = argv[1];
        ScaLBL_FreeLeeModel LeeModel( rank, nprocs, comm );
        LeeModel.ReadParams( filename );
        LeeModel.SetDomain();
        Initialize_Mask(LeeModel);
        //LeeModel.Create_DummyPhase_MGTest(); 
        LeeModel.Create_TwoFluid(); 
        
        Initialize_DummyPhaseField(LeeModel,1.0, 2.0, 3.0);
        LeeModel.WriteDebug_TwoFluid();

        PROFILE_STOP( "Main" );
        PROFILE_SAVE( file, level );
        // ****************************************************


    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();
    return 0;

}
