
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/MPI.h"
#include "models/MRTModel.h"

void ParallelPlates(ScaLBL_MRTModel &MRT){
	// initialize empty domain
  int i,j,k,n;
	int Nx = MRT.Nx;
	int Ny = MRT.Ny;
	int Nz = MRT.Nz;
	Array<char> id_solid(Nx,Ny,Nz);

	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				if (i<2) MRT.Mask->id[n] = 0;
				else if (i>Nx-3) MRT.Mask->id[n] = 0;
				else MRT.Mask->id[n]=1;
				if (MRT.Mask->id[n] == 0) 	id_solid(i,j,k) = 1;
				else 						id_solid(i,j,k) = 0;
			}
		}
	}
	// Initialize the signed distance function
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				// Initialize distance to +/- 1
				MRT.Distance(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
			}
		}
	}
	if (MRT.rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
	CalcDist(MRT.Distance,id_solid,*MRT.Dm);
    if (MRT.rank == 0) cout << "Domain set." << endl;
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
			printf("Running Unit Test: TestPoiseuille	\n");
			printf("********************************************************\n");
		}
		
		int i,j,k,n;
		ScaLBL_MRTModel MRT(rank,nprocs,comm);
		auto filename = argv[1];
		MRT.ReadParams(filename);
		MRT.SetDomain();    // this reads in the domain 
		ParallelPlates(MRT);
		MRT.Create();       // creating the model will create data structure to match the pore structure and allocate variables
		MRT.Initialize();   // initializing the model will set initial conditions for variables
		MRT.Run();	 
		double *Vz;  	Vz= new double [3*MRT.Np];

		int SIZE=MRT.Np*sizeof(double);
		ScaLBL_D3Q19_Momentum(MRT.fq,MRT.Velocity, MRT.Np);
		ScaLBL_DeviceBarrier(); comm.barrier();
		ScaLBL_CopyToHost(&Vz[0],&MRT.Velocity[0],3*SIZE);
		
		if (rank == 0) printf("Force: %f,%f,%f \n",MRT.Fx,MRT.Fy,MRT.Fz);
		double mu = MRT.mu;
		int Nx = MRT.Nx;
		int Ny = MRT.Ny;
		int Nz = MRT.Nz;
		double Fz = MRT.Fz;
		double vz;
		double W = 1.f*Nx-4.f;
		j=Ny/2; k=Nz/2;
		if (rank == 0) printf("Channel width=%f \n",W);
		if (rank == 0) printf("ID flag vz       analytical\n");
		comm.barrier();

		if (rank == 0) {
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				printf("%i ",MRT.Mask->id[n]);
				n = MRT.Map(i,j,k);
				//printf("%i,%i,%i; %i :",i,j,k,n);
				if (n<0) {vz =0.f; printf(" b    "); }
				else { vz=Vz[n+2*MRT.Np]; printf(" a    "); }
				printf("%f ",vz);
				//Analytical solution
				double x=1.f*i-1.5;
				if (n<0) vz=0.f;
				else vz=Fz*x*(W-x)/(2.f*mu);
				printf("%f\n",vz);
			}
			printf("\n");
		}
		if (rank == 1) {
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				printf("%i ",MRT.Mask->id[n]);
				n = MRT.Map(i,j,k);
				//printf("%i,%i,%i; %i :",i,j,k,n);
				if (n<0) {vz =0.f; printf(" b    "); }
				else { vz=Vz[n+2*MRT.Np]; printf(" a    "); }
				printf("%f ",vz);
				//Analytical solution
				double x=1.f*i-1.5;
				if (n<0) vz=0.f;
				else vz=Fz*x*(W-x)/(2.f*mu);
				printf("%f\n",vz);
			}
			printf("\n");
		}
	}

    Utilities::shutdown();

	return check;
}
