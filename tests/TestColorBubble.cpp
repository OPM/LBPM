
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/MPI_Helpers.h"
#include "models/ColorModel.h"

using namespace std;
inline void InitializeBubble(ScaLBL_ColorModel &ColorModel, double BubbleRadius){
		// initialize a bubble
  int i,j,k,n;
  int Np=0;
  double sum=0.f;
  int rank = ColorModel.Mask->rank();
		int nprocx = ColorModel.Mask->nprocx();
		int nprocy = ColorModel.Mask->nprocy();
		int nprocz = ColorModel.Mask->nprocz();
		int Nx = ColorModel.Mask->Nx;
		int Ny = ColorModel.Mask->Ny;
		int Nz = ColorModel.Mask->Nz;
		if (rank == 0) cout << "Setting up bubble..." << endl;
		sum=0; Np=0;
		for (k=0;k<Nz;k++){
		  for (j=0;j<Ny;j++){
		    for (i=0;i<Nx;i++){
		      n = k*Nx*Ny + j*Nz + i;
		      ColorModel.Averages->SDs(i,j,k) = 100.f;
		      // Initialize phase positions field
		      if (ColorModel.Averages->SDs(i,j,k) < 0.0){
			ColorModel.Mask->id[n] = 0;
		      }
		      else {
			sum++;
			Np++;
		      }
		    }
		  }
		}
		// Initialize the bubble
		for (k=0;k<Nz;k++){
		  for (j=0;j<Ny;j++){
		    for (i=0;i<Nx;i++){
		      n = k*Nx*Ny + j*Nz + i;
		      int iglobal= i+(Nx-2)*ColorModel.Dm->iproc();
		      int jglobal= j+(Ny-2)*ColorModel.Dm->jproc();
		      int kglobal= k+(Nz-2)*ColorModel.Dm->kproc();
		      // Initialize phase position field for parallel bubble test
		      if ((iglobal-0.5*(Nx-2)*nprocx)*(iglobal-0.5*(Nx-2)*nprocx)
			  +(jglobal-0.5*(Ny-2)*nprocy)*(jglobal-0.5*(Ny-2)*nprocy)
			  +(kglobal-0.5*(Nz-2)*nprocz)*(kglobal-0.5*(Nz-2)*nprocz) < BubbleRadius*BubbleRadius){
                        ColorModel.Mask->id[n] = 2;
		      }
		      else{
                        ColorModel.Mask->id[n]=1;
		      }
		    }
		  }
		}

}
//***************************************************************************************
int main(int argc, char **argv)
{
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	int check;
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
		ScaLBL_ColorModel ColorModel;
		ColorModel.ReadParams(filename);
		//ColorModel.ReadInput();
		double radius=15.5;
		InitializeBubble(ColorModel,radius);
		ColorModel.Create();
		ColorModel.Initialize();
		ColorModel.Run();
	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************

	return check;
}

