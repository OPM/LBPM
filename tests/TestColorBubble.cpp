
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

inline void AssignComponentLabels(char *id, double *phase, int Nx, int Ny, int Nz, int rank, MPI_Comm comm)
{
	int NLABELS=0;
	char VALUE=0;
	double AFFINITY=0.f;
	
	vector <char> Label;
	vector <double> Affinity;
	// Read the labels
	if (rank==0){
		printf("Component labels:\n");
		ifstream iFILE("ComponentLabels.csv");\
		if (iFILE.good()){
			while (!iFILE.eof()){
				iFILE>>VALUE;
				iFILE>>AFFINITY;
				Label.push_back(VALUE);
				Affinity.push_back(AFFINITY);
				NLABELS++;
				printf("%i %f\n",VALUE,AFFINITY);
			}
		}
		else{
			printf("Using default labels: Solid (0 --> -1.0), NWP (1 --> 1.0), WP (2 --> -1.0)\n");
			// Set default values
			VALUE=0; AFFINITY=-1.0;
			Label.push_back(VALUE);
			Affinity.push_back(AFFINITY);
			NLABELS++;
			printf("%i %f\n",VALUE,AFFINITY);
			VALUE=1; AFFINITY=1.0;
			Label.push_back(VALUE);
			Affinity.push_back(AFFINITY);
			NLABELS++;
			printf("%i %f\n",VALUE,AFFINITY);
			VALUE=2; AFFINITY=-1.0;
			Label.push_back(VALUE);
			Affinity.push_back(AFFINITY);
			NLABELS++;
			printf("%i %f\n",VALUE,AFFINITY);
		}
	}
	// Broadcast the list
	MPI_Bcast(&NLABELS,1,MPI_INT,0,comm);
	
	// Copy into contiguous buffers
	char *LabelList;
	double * AffinityList;
	LabelList=new char[NLABELS];
	AffinityList=new double[NLABELS];
	MPI_Bcast(&LabelList,NLABELS,MPI_CHAR,0,comm);
	MPI_Bcast(&AffinityList,NLABELS,MPI_DOUBLE,0,comm);

	// Assign the labels
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				int n = k*Nx*Ny+j*Nx+i;
				VALUE=id[n];
				// Assign the affinity from the paired list
				for (int idx=0; idx < NLABELS; idx++){
					if (VALUE == LabelList[idx]){
						AFFINITY=AffinityList[idx];
						idx = NLABELS;
					}
				}
				phase[n] = AFFINITY;
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
	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************

	return check;
}

