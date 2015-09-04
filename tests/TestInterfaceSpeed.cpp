#include <iostream>
#include <math.h>

#include "common/TwoPhase.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "IO/Mesh.h"
#include "IO/Writer.h"
#include "ProfilerApp.h"

#define RADIUS 15
#define CAPRAD 20
#define HEIGHT 15.5
#define N 60
#define SPEED -1
#define PI 3.14159

int main (int argc, char *argv[])
{
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);

	int npx,npy,npz;
	int i,j,k,n;
	int Nx,Ny,Nz;
	double Lx,Ly,Lz;
	Nx=Ny=Nz=N;
	npx=npy=npz=1;
	Lx=Ly=Lz=1.0;
	int BC=0;	// periodic boundary condition

	Domain Dm(Nx,Ny,Nz,rank,npx,npy,npz,Lx,Ly,Lz,BC);

	for (i=0; i<Dm.Nx*Dm.Ny*Dm.Nz; i++) Dm.id[i] = 1;

	Dm.CommInit(comm);

	TwoPhase Averages(Dm);
	int timestep=0;

	Nx = Dm.Nx;
	Ny = Dm.Ny;
	Nz = Dm.Nz;

	double Cx,Cy,Cz;
	double dist1,dist2;

	Cx = Cy = Cz = N*0.5;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;
				dist2 = fabs(Cz-k)-HEIGHT;

				Averages.Phase_tminus(i,j,k) = dist2;
			}
		} 
	}
	Cz += SPEED;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				
				dist1 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)) - RADIUS;
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;
				dist2 = fabs(Cz-k)-HEIGHT;

				Averages.SDs(i,j,k) = -dist1;
				Averages.Phase(i,j,k) = dist2;
				Averages.SDn(i,j,k) = dist2;
			}
		}   
	}
	Cz += SPEED;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;
				dist2 = fabs(Cz-k)-HEIGHT;

				Averages.Phase_tplus(i,j,k) = dist2;
			}
		}   
	}
	
	//...........................................................................

	//....................................................................
	// The following only need to be done once
	//Averages.SetupCubes(Dm);
	Averages.UpdateSolid(); 	// unless the solid is deformable!
	//....................................................................
	// The following need to be called each time new averages are computed
	Averages.Initialize();
	Averages.UpdateMeshValues();
	Averages.ComputeLocal();
	Averages.Reduce();
	Averages.PrintAll(timestep);
	//....................................................................
	
	printf("-------------------------------- \n");
	printf("NWP volume = %f \n", Averages.nwp_volume);
	printf("Area wn = %f, Analytical = %f \n", Averages.awn,2*PI*RADIUS*RADIUS);
	printf("Area ns = %f, Analytical = %f \n", Averages.ans, 2*PI*RADIUS*(N-2)-4*PI*RADIUS*HEIGHT);
	printf("Area ws = %f, Analytical = %f \n", Averages.aws, 4*PI*RADIUS*HEIGHT);
	printf("Area s = %f, Analytical = %f \n", Averages.As, 2*PI*RADIUS*(N-2));
	printf("Length wns = %f, Analytical = %f \n", Averages.lwns, 4*PI*RADIUS);
	printf("Geodesic curvature (wns) = %f, Analytical = %f \n", Averages.KGwns_global, 0.0);
	printf("Normal curvature (wns) = %f, Analytical = %f \n", Averages.KNwns_global, 1.0/RADIUS);
//	printf("Cos(theta_wns) = %f, Analytical = %f \n",efawns/lwns,1.0*RADIUS/CAPRAD);
	printf("Interface Velocity = %f,%f,%f \n",Averages.vawn_global(0),Averages.vawn_global(1),Averages.vawn_global(2));
	printf("Common Curve Velocity = %f,%f,%f \n",Averages.vawns_global(0),Averages.vawns_global(1),Averages.vawns_global(2));
	printf("-------------------------------- \n");	
	//.........................................................................	
	
	int toReturn = 0;

	if ( fabs(Averages.vawn_global(2)+0.2) > 0.01){
		printf("TestInterfaceSpeed: Error too high for kinematic velocity of wn interface \n");
		toReturn=1;
	}
	if ( fabs(Averages.vawns_global(2)+0.2) > 0.01){
		printf("TestInterfaceSpeed: Error too high for kinematic velocity of common curve \n");
		toReturn=2;
	}

	return toReturn;

	// ****************************************************
	MPI_Barrier(comm);
	return 0;
	MPI_Finalize();
	// ****************************************************
}
