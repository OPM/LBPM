
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

std::shared_ptr<Database> loadInputs( int nprocs )
{
    auto db = std::make_shared<Database>( "Domain.in" );
    const int dim = 50;
    db->putScalar<int>( "BC", 0 );
    if ( nprocs == 1 ){
        db->putVector<int>( "nproc", { 1, 1, 1 } );
        db->putVector<int>( "n", { 3, 1, 1 } );
        db->putScalar<int>( "nspheres", 0 );
        db->putVector<double>( "L", { 1, 1, 1 } );
    } else if ( nprocs == 2 ) {
        db->putVector<int>( "nproc", { 2, 1, 1 } );
        db->putVector<int>( "n", { dim, dim, dim } );
        db->putScalar<int>( "nspheres", 0 );
        db->putVector<double>( "L", { 1, 1, 1 } );
    } else if ( nprocs == 4 ) {
        db->putVector<int>( "nproc", { 2, 2, 1 } );
        db->putVector<int>( "n", { dim, dim, dim } );
        db->putScalar<int>( "nspheres", 0 );
        db->putVector<double>( "L", { 1, 1, 1 } );
    } else if (nprocs==8){
        db->putVector<int>( "nproc", { 2, 2, 2 } );
        db->putVector<int>( "n", { dim, dim, dim } );
        db->putScalar<int>( "nspheres", 0 );
        db->putVector<double>( "L", { 1, 1, 1 } );
    }
    return db;
}

void InitializeSquareTube(ScaLBL_ColorModel &ColorModel){
  int i,j,k,n;
  int rank = ColorModel.Mask->rank();
  int Nx = ColorModel.Mask->Nx;
  int Ny = ColorModel.Mask->Ny;
  int Nz = ColorModel.Mask->Nz;
  int nprocx = ColorModel.Mask->rank_info.nx;
  int nprocy = ColorModel.Mask->rank_info.ny;
  int iproc = ColorModel.Mask->rank_info.ix;
  int jproc = ColorModel.Mask->rank_info.jy;
  int kproc = ColorModel.Mask->rank_info.kz;

	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
                n = k*Nx*Ny + j*Nx + i;
				ColorModel.Mask->id[n]=0;
			}
		}
	}
	printf("rank=%i, %i,%i,%i \n",rank,iproc,jproc,kproc);
	// Initialize a square tube
	for (k=1;k<Nz-1;k++){
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				n = k*Nx*Ny + j*Nx + i;
				int iglobal= i+(Nx-2)*iproc;
				int jglobal= j+(Ny-2)*jproc;
				int kglobal= k+(Nz-2)*kproc;

				// Initialize phase position field for parallel bubble test
				if (iglobal < 2)						ColorModel.Mask->id[n]=0;
				else if (iglobal > (Nx-2)*nprocx-2)	    ColorModel.Mask->id[n]=0;
				else if (jglobal < 2)					ColorModel.Mask->id[n]=0;
				else if (jglobal > (Ny-2)*nprocy-2)	    ColorModel.Mask->id[n]=0;
				else if (kglobal < 20)					ColorModel.Mask->id[n]=1;
				else									ColorModel.Mask->id[n]=2;
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
	int check=0;
	{
		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Color Model: TestColor	\n");
			printf("********************************************************\n");
		}
		auto filename = argv[1];
		ScaLBL_ColorModel ColorModel(rank,nprocs,comm);
		ColorModel.ReadParams(filename);
		ColorModel.SetDomain();    
		//ColorModel.ReadInput(); 
		InitializeSquareTube(ColorModel);
		ColorModel.Create();       // creating the model will create data structure to match the pore structure and allocate variables
		ColorModel.Initialize();   // initializing the model will set initial conditions for variables
		ColorModel.Run();	       
		ColorModel.WriteDebug(); 
 
	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************

	return check;
}

