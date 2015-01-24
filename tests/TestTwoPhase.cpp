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
	int rank,npx,npy,npz;
	int i,j,k;
	int Nx,Ny,Nz;
	double Lx,Ly,Lz;
	Nx=Ny=Nz=40;
	rank=0;
	npx=npy=npz=1;
	Lx=Ly=Lz=1.0;

	FILE *TIMELOG;
	if (rank==0)	TIMELOG = fopen("timelog.tcat","a+");

	Domain Dm(Nx,Ny,Nz,rank,npx,npy,npz,Lx,Ly,Lz);

	TwoPhase Averages(Dm);
	int timestep=0;

	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				Averages.Phase(i,j,k) = sqrt((1.0*i-0.5*Nx)*(1.0*i-0.5*Nx)+(1.0*j-0.5*Ny)*(1.0*j-0.5*Ny)+(1.0*k-0.5*Nz)*(1.0*k-0.5*Nz))-0.3*Nx;
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

	printf("Nx=%i, \n",Averages.Nx);

	Averages.SetupCubes(Dm);
	Averages.Initialize();
	Averages.UpdateMeshValues();
	Averages.ComputeLocal();
	Averages.Reduce();
	Averages.PrintAll(timestep,TIMELOG);

	return 0;
}
