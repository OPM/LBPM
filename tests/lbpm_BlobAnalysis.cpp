/*
This code computes TCAT averages on a blob-by-blob basis in parallel
It requires that the blobs be labeled using BlobIdentify.cpp 
James E. McClure 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "Domain.h"
#include "TwoPhase.h"
#include "common/MPI_Helpers.h"

int main(int argc, char **argv)
{
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	// parallel domain size (# of sub-domains)
	int nprocx,nprocy,nprocz;
	int iproc,jproc,kproc;
	int Nx,Ny,Nz,N,nspheres;
	double Lx,Ly,Lz;

	int BC,nblobs;

	if (rank==0){
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		ifstream domain("Domain.in");
		domain >> nprocx;
		domain >> nprocy;
		domain >> nprocz;
		domain >> Nx;
		domain >> Ny;
		domain >> Nz;
		domain >> nspheres;
		domain >> Lx;
		domain >> Ly;
		domain >> Lz;
		//.......................................................................
	}
	//.................................................
	MPI_Barrier(MPI_COMM_WORLD);
	// Computational domain
	MPI_Bcast(&Nx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Ny,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Nz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocy,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nspheres,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	//.................................................
	Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
	TwoPhase Averages(Dm);
	//.......................................................................
	// Filenames used
	char LocalRankString[8];
	char LocalRankFilename[40];
	char LocalRestartFile[40];
	char tmpstr[10];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
	//...........................................................................
	if (rank == 0) cout << "Reading in domain from signed distance function..." << endl;
	//.......................................................................
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
	ReadBinaryFile(LocalRankFilename, Averages.SDs.data, N);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) cout << "Domain set." << endl;
	//.......................................................................
	//copies of data needed to perform checkpointing from cpu
	double *Den, *DistEven, *DistOdd;
	Den = new double[2*N];
	DistEven = new double[10*N];
	DistOdd = new double[9*N];
	//.........................................................................
	if (rank==0) printf("Reading restart file! \n");
	// Read in the restart file to CPU buffers
	ReadCheckpoint(LocalRestartFile, Den, DistEven, DistOdd, N);
	MPI_Barrier(MPI_COMM_WORLD);
	//.........................................................................
	// Populate the arrays needed to perform averaging
	Averages.Phase();
	Averages.LocalBlobID();


	/*	Averages.Initialize();
	Averages.ComputeDelPhi();
	Averages.ColorToSignedDistance(beta,Averages.Phase.data,Averages.SDn.data);
	Averages.UpdateMeshValues();
	Averages.ComputeLocal();
	Averages.Reduce();
	Averages.PrintAll(timestep);
*/
	// ****************************************************
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	// ****************************************************
}
