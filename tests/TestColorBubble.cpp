
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/MPI.h"
#include "models/ColorModel.h"

using namespace std;
inline void InitializeBubble(ScaLBL_ColorModel &ColorModel, double BubbleRadius){
	// initialize a bubble
	int i,j,k,n;
	int rank = ColorModel.Dm->rank();
	int nprocx = ColorModel.Dm->nprocx();
	int nprocy = ColorModel.Dm->nprocy();
	int nprocz = ColorModel.Dm->nprocz();
	int Nx = ColorModel.Dm->Nx;
	int Ny = ColorModel.Dm->Ny;
	int Nz = ColorModel.Dm->Nz;
	if (rank == 0) cout << "Setting up bubble..." << endl;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;
				ColorModel.Averages->SDs(i,j,k) = 100.f;
			}
		}
	}
	// Initialize the bubble
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nx + i;
				double iglobal= double(i+(Nx-2)*ColorModel.Dm->iproc())-double((Nx-2)*nprocx)*0.5;
				double jglobal= double(j+(Ny-2)*ColorModel.Dm->jproc())-double((Ny-2)*nprocy)*0.5;
				double kglobal= double(k+(Nz-2)*ColorModel.Dm->kproc())-double((Nz-2)*nprocz)*0.5;
				// Initialize phase position field for parallel bubble test
				if ((iglobal*iglobal)+(jglobal*jglobal)+(kglobal*kglobal) < BubbleRadius*BubbleRadius){
					ColorModel.Mask->id[n] = 2;
					ColorModel.Mask->id[n] = 2;
				}
				else{
					ColorModel.Mask->id[n]=1;
					ColorModel.Mask->id[n]=1;
				}
				ColorModel.id[n] = ColorModel.Mask->id[n];
				ColorModel.Dm->id[n] = ColorModel.Mask->id[n];
			}
		}
	}
	
	FILE *OUTFILE;
	char LocalRankFilename[40];
	sprintf(LocalRankFilename,"Bubble.%05i.raw",rank);
	OUTFILE = fopen(LocalRankFilename,"wb");
	fwrite(ColorModel.id,1,Nx*Ny*Nz,OUTFILE);
	fclose(OUTFILE);
	// initialize the phase indicator field
}
//***************************************************************************************
int main(int argc, char **argv)
{
	// Initialize MPI
    Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();
	int check=0;
	{
		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Color Model: TestColor	\n");
			printf("********************************************************\n");
			if ( argc < 2 ) {
				std::cerr << "Invalid number of arguments, no input file specified\n";
				return -1;
			}
		}
		auto filename = argv[1];
		ScaLBL_ColorModel ColorModel(rank,nprocs,comm);
		ColorModel.ReadParams(filename);
		ColorModel.SetDomain();    
		//ColorModel.ReadInput(); 
		double radius=15.5;
		InitializeBubble(ColorModel,radius);
		ColorModel.Create();       // creating the model will create data structure to match the pore structure and allocate variables
		ColorModel.Initialize();   // initializing the model will set initial conditions for variables
		ColorModel.Run();	       
		ColorModel.WriteDebug();
	}
    Utilities::shutdown();

	return check;
}

