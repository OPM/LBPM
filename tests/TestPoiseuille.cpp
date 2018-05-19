
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/MPI_Helpers.h"
#include "models/MRTModel.h"

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
			printf("Running Unit Test: TestPoiseuille	\n");
			printf("********************************************************\n");
		}
		
		int i,j,k,n;
		ScaLBL_MRTModel MRT(rank,nprocs,comm);
		auto filename = argv[1];
		MRT.ReadParams(filename);
		MRT.SetDomain();    // this reads in the domain 
		MRT.Create();       // creating the model will create data structure to match the pore structure and allocate variables
		MRT.Initialize();   // initializing the model will set initial conditions for variables
		MRT.Run();	 
		double *Vz;  	Vz= new double [MRT.Np];
		MRT.VelocityField(Vz);

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
		MPI_Barrier(comm);

		if (rank == 0) {
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				printf("%i ",MRT.Mask->id[n]);
				n = MRT.Map(i,j,k);
				//printf("%i,%i,%i; %i :",i,j,k,n);
				if (n<0) {vz =0.f; printf(" b    "); }
				else { vz=Vz[n]; printf(" a    "); }
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
				else { vz=Vz[n]; printf(" a    "); }
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

	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************

	return check;
}
