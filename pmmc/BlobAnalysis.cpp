// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "pmmc.h"
#include "Domain.h"

using namespace std;

inline void ReadFromAllRanks(){
	
}

int main(int argc, char **argv)
{
	//.......................................................................
	int nprocx,nprocy,nprocz,nprocs;
	int Nx, Ny, Nz;
	int nx,ny,nz;
	int nspheres;
	double Lx,Ly,Lz;
	//.......................................................................
	int i,j,k;
	//.......................................................................
	// Reading the domain information file
	//.......................................................................
	ifstream domain("Domain.in");
	domain >> nprocx;
	domain >> nprocy;
	domain >> nprocz;
	domain >> nx;
	domain >> ny;
	domain >> nz;
	domain >> nspheres;
	domain >> Lx;
	domain >> Ly;
	domain >> Lz;
	//.......................................................................

	nprocs = nprocx*nproy*nprocz;
	printf("Number of MPI ranks: %i \n", nprocs);
	Nx = (nx-2)*nprocx;
	Ny = (ny-2)*nprocy;
	Nz = (nz-2)*nprocz;
	printf("Full domain size: %i x %i x %i  \n", Nx,Ny,Nz);
	
	// Filenames used
	char LocalRankString[8];
	char LocalRankFilename[40];
	char LocalRestartFile[40];
	char BaseFilename[20];
	char tmpstr[10];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(BaseFilename,"%s","Phase");
	sprintf(LocalRankFilename,"%s%s",BaseFilename,LocalRankString);

	                   
	IntArray LocalBlobID(Nx,Ny,Nz);
	DoubleArray Phase(Nx,Ny,Nz);
	DoubleArray SignDist(Nx,Ny,Nz);
	DoubleArray Press(Nx,Ny,Nz);
	DoubleArray Vel_x(Nx,Ny,Nz);			// Velocity
	DoubleArray Vel_y(Nx,Ny,Nz);
	DoubleArray Vel_z(Nx,Ny,Nz);
	DoubleArray MeanCurvature(Nx,Ny,Nz);
	DoubleArray GaussCurvature(Nx,Ny,Nz);
	DoubleArray SignDist_x(Nx,Ny,Nz);		// Gradient of the signed distance
	DoubleArray SignDist_y(Nx,Ny,Nz);
	DoubleArray SignDist_z(Nx,Ny,Nz);
	DoubleArray Phase_x(Nx,Ny,Nz);			// Gradient of the phase indicator field
	DoubleArray Phase_y(Nx,Ny,Nz);
	DoubleArray Phase_z(Nx,Ny,Nz);
	
	double *Temp;
	Temp = new double [3*nx*ny*nz];
	
	// read the files and populate main arrays
	for (int proc=0; proc<nprocs; proc++){
		
		sprintf(LocalRankString,"%05d",proc);
		sprintf(LocalRankFilename,"%s%s",BaseFilename,LocalRankString);
		printf("Reading file %s \n",LocalRankFilename);
		ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);	
		
		for (k=1; k<nz-1; k++){
			for (j=1; j<ny-1; j++){
				for (i=1; i<nz-1; i++){
					n = k*nx*ny+j*nx+i;
					Phase(i-1,j-1,k-1) = Temp[n];
				}
			}
		}
		
	}
	
	printf("Read %i ranks of %s \n",nprocs,BaseFilename);
	
}

