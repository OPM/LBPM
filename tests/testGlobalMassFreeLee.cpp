#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "common/Utilities.h"
#include "models/FreeLeeModel.h"

//*******************************************************************
// Implementation of Free-Energy Two-Phase LBM (Lee model)
//*******************************************************************

int main( int argc, char **argv )
{
	  
    // Initialize
    Utilities::startup( argc, argv );

    // Load the input database
    auto db = std::make_shared<Database>( argv[1] );
  
    { // Limit scope so variables that contain communicators will free before MPI_Finialize

        Utilities::MPI comm( MPI_COMM_WORLD );
        int rank = comm.getRank();
        int nprocs = comm.getSize();

        if (rank == 0){
            printf("********************************************************\n");
            printf("Running Free Energy Lee LBM	\n");
            printf("********************************************************\n");
        }
        // Initialize compute device
        int device=ScaLBL_SetDevice(rank);
        NULL_USE( device );
        ScaLBL_DeviceBarrier();
        comm.barrier();

        PROFILE_ENABLE(1);
        //PROFILE_ENABLE_TRACE();
        //PROFILE_ENABLE_MEMORY();
        PROFILE_SYNCHRONIZE();
        PROFILE_START("Main");
        Utilities::setErrorHandlers();

        auto filename = argv[1];
        ScaLBL_FreeLeeModel LeeModel( rank,nprocs,comm );
        LeeModel.ReadParams( filename );
        LeeModel.SetDomain();    
        LeeModel.ReadInput();    
        LeeModel.Create_TwoFluid();       
        LeeModel.Initialize_TwoFluid(); 
        /* check neighbors */
        
        
        /* Copy the initial density to test that global mass is conserved */
        int Nx = LeeModel.Dm->Nx;
        int Ny = LeeModel.Dm->Ny;
        int Nz = LeeModel.Dm->Nz;
        DoubleArray DensityInit(Nx,Ny,Nz);
        LeeModel.ScaLBL_Comm->RegularLayout(LeeModel.Map,LeeModel.Den,DensityInit);
        
        double MLUPS = LeeModel.Run_TwoFluid(LeeModel.timestepMax);	
        
        DoubleArray DensityFinal(Nx,Ny,Nz);
        LeeModel.ScaLBL_Comm->RegularLayout(LeeModel.Map,LeeModel.Den,DensityFinal);
        
        DoubleArray DensityChange(Nx,Ny,Nz);
        double totalChange=0.0;
        for (int k=1; k<Nz-1; k++){
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++){
                	double change = DensityFinal(i,j,k)-DensityInit(i,j,k);
                	DensityChange(i,j,k) = change;
                	totalChange += change;
                }
            }
        }
        printf("Density change, %f\n", totalChange);
        
    	FILE *OUTFILE;
        char LocalRankFilename[40];
    	sprintf(LocalRankFilename,"DensityChange.%05i.raw",rank);
    	OUTFILE = fopen(LocalRankFilename,"wb");
    	fwrite(DensityChange.data(),8,Nx*Ny*Nz,OUTFILE);
    	fclose(OUTFILE);
   
        //LeeModel.WriteDebug_TwoFluid();

        PROFILE_STOP("Main");
        // ****************************************************


    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();
    return 0;
}
