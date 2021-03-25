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

inline int MultiHaloNeighborCheck(ScaLBL_FreeLeeModel &LeeModel){
	int i,j,k,iq,stride,nread;
	int Nxh = LeeModel.Nxh;
	int Nyh = LeeModel.Nyh;
	int Np = LeeModel.Np;
	
	int *TmpMap;
	TmpMap = new int[Np];
	ScaLBL_CopyToHost(TmpMap, LeeModel.dvcMap, Np*sizeof(int));
	
	int *neighborList;
	neighborList = new int[18*Np];
	ScaLBL_CopyToHost(neighborList, LeeModel.NeighborList, 18*Np*sizeof(int));
	printf("Check stride for interior neighbors \n");
	int count = 0;
	for (int n=LeeModel.ScaLBL_Comm->FirstInterior(); n<LeeModel.ScaLBL_Comm->LastInterior(); n++){
		// q=0
			int idx = TmpMap[n];
			k = idx/Nxh/Nyh;
			j = (idx-k*Nxh*Nyh)/Nxh;
			i = (idx-k*Nxh*Nyh -j*Nxh);

			// q=1
			nread = neighborList[n]; 
			iq = TmpMap[nread%Np];
			stride = idx - iq;
			if (stride != 1){
				printf(" %i, %i, %i   q = 1 stride=%i \n ",i,j,k,stride);
				count++;
			}
			
			// q=2
			nread = neighborList[n+Np]; 
			iq = TmpMap[nread%Np];
			stride = iq - idx;
			if (stride != 1){
				printf(" %i, %i, %i   q = 2 stride=%i \n ",i,j,k,stride);
				count++;
			}
			

			// q=3
			nread = neighborList[n+2*Np]; 
			iq = TmpMap[nread%Np];
			stride = idx - iq;
			if (stride != Nxh){
				printf(" %i, %i, %i   q = 3 stride=%i \n ",i,j,k,stride);
				count++;
			}

			// q = 4
			nread = neighborList[n+3*Np]; 
			iq = TmpMap[nread%Np];
			stride = iq-idx;
			if (stride != Nxh){
				printf(" %i, %i, %i   q = 4 stride=%i \n ",i,j,k,stride);
				count++;
			}
			

			// q=5
			nread = neighborList[n+4*Np];
			iq = TmpMap[nread%Np];
			stride = idx - iq;
			if (stride != Nxh*Nyh){
				count++;
				printf(" %i, %i, %i   q = 5 stride=%i \n ",i,j,k,stride);
			}
			
			// q = 6
			nread = neighborList[n+5*Np];
			iq = TmpMap[nread%Np];
			stride = iq - idx;
			if (stride != Nxh*Nyh){
				count++;
				printf(" %i, %i, %i   q = 6 stride=%i \n ",i,j,k,stride);
			}

	}
	return count;
}

int main( int argc, char **argv )
{

    // Initialize
    Utilities::startup( argc, argv );
    int errors = 0;
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
        
        errors=MultiHaloNeighborCheck(LeeModel);
        
        Initialize_DummyPhaseField(LeeModel,1.0, 2.0, 3.0);
        LeeModel.WriteDebug_TwoFluid();

        PROFILE_STOP( "Main" );
        PROFILE_SAVE( file, level );
        // ****************************************************


    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();
    return errors;

}
