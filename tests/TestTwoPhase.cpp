// Unit test for TwoPhase averaging class

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "TwoPhase.h"
#include "Extras.h"
#include "D3Q19.h"
#include "D3Q7.h"
#include "Color.h"
#include "common/MPI.h"
#include "Communication.h"
#include "IO/Mesh.h"
#include "IO/Writer.h"
#include "ProfilerApp.h"

int main(int argc, char **argv)
{
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	printf("Running two-phase averaging test on %i processors \n",nprocs);

	int npx,npy,npz;
	int i,j,k,n;
	int Nx,Ny,Nz;
	double Lx,Ly,Lz;
	Nx=Ny=Nz=40;
	npx=npy=npz=2;
//	npz=nprocs;
	Lx=Ly=Lz=1.0;
	int BC=0;	// periodic boundary condition

	Domain Dm(Nx,Ny,Nz,rank,npx,npy,npz,Lx,Ly,Lz,BC);

	for (i=0; i<Dm.Nx*Dm.Ny*Dm.Nz; i++) Dm.id[i] = 1;

	Dm.CommInit(MPI_COMM_WORLD);

	TwoPhase Averages(Dm);
	int timestep=0;

	Nx = Dm.Nx;
	Ny = Dm.Ny;
	Nz = Dm.Nz;

	int size_x = npx*(Nx-2);
	int size_y = npy*(Ny-2);
	int size_z = npz*(Nz-2);

	int radius = 0.3*min(min(size_y,size_z),size_x);
	radius = 15;
//	printf("rank=%i,iproc= %i,jproc=%i,kproc=%i \n",rank,Dm.iproc,Dm.jproc,Dm.kproc);
	printf("Sphere radius = %i \n",radius);
//	printf("sendcount_x = %i \n",Dm.sendCount_X);

	// Initializing a simple case
	double x,y,z,distance;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){

				x = 1.0*(Dm.iproc*(Nx-2)+i-1)-0.5*size_x;
				y = 1.0*(Dm.jproc*(Ny-2)+j-1)-0.5*size_y;
				z = 1.0*(Dm.kproc*(Nz-2)+k-1)-0.5*size_z;
				distance=sqrt(x*x+y*y+z*z)-1.0*radius;
				Averages.Phase(i,j,k) = distance;
				Averages.SDs(i,j,k) = 100.0;
				Averages.SDn(i,j,k) = Averages.Phase(i,j,k);
				Averages.Phase_tplus(i,j,k) = Averages.Phase(i,j,k);
				Averages.Phase_tminus(i,j,k) = Averages.Phase(i,j,k);
				Averages.Press(i,j,k) = 0.0;
				Averages.Vel_x(i,j,k) = 0.0;
				Averages.Vel_y(i,j,k) = 0.0;
				Averages.Vel_z(i,j,k) = 0.0;
			}
		}
	}

	//....................................................................
	// The following only need to be done once
	Averages.SetupCubes(Dm);
	Averages.UpdateSolid(); 	// unless the solid is deformable!
	//....................................................................
	// The following need to be called each time new averages are computed
	Averages.Initialize();
	Averages.UpdateMeshValues();
	Averages.ComputeLocal();
	Averages.Reduce();
	Averages.PrintAll(timestep);
	//....................................................................


	if (rank==0){
		FILE *PHASE;
		PHASE = fopen("Phase.00000","wb");
		fwrite(Averages.MeanCurvature.data,8,Nx*Ny*Nz,PHASE);
		fclose(PHASE);
	}
	// ****************************************************
	MPI_Barrier(MPI_COMM_WORLD);
	return 0;
	MPI_Finalize();
	// ****************************************************
}
