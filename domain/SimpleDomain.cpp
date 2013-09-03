#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "Domain.h"

int main(int argc, char **argv)
{
	//.......................................................................
	// Variable declaration
	//.......................................................................
	int i,j,k,N;
	int Nx,Ny,Nz;
	int nspheres;
	double Lx,Ly,Lz;
	// parallel domain size (# of sub-domains)
	int nprocx,nprocy,nprocz;
//	int iproc,jproc,kproc;
	//.......................................................................

	//.......................................................................
	// Reading the input file
	//.......................................................................
	ifstream domain("Domain.in");
	domain >> nspheres;
	domain >> Lx;
	domain >> Ly;
	domain >> Lz;
	domain >> Nx;
	domain >> Ny;
	domain >> Nz;
	domain >> nprocx;
	domain >> nprocy;
	domain >> nprocz;
	//.......................................................................
	printf("********************************************************\n");
	printf("Sub-domain size = %i x %i x %i\n",Nz,Nz,Nz);
	printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
	printf("Number of spheres = %i \n", nspheres);
	printf("Domain Length (x) = %f \n", Lx);
	printf("Domain Length (y) = %f \n", Ly);
	printf("Domain Length (z) = %f \n", Lz);
	printf("********************************************************\n");
	//.......................................................................
	N = Nx*Ny*Nz;
	printf("Read input media... \n");

	char *id;
	id = new char[N];
	//.......................................................................
	// Read from a sphere packing
	//.......................................................................
	double *cx,*cy,*cz,*r;
	cx = new double[nspheres];
	cy = new double[nspheres];
	cz = new double[nspheres];
	r = new double[nspheres];
	//.......................................................................
	printf("Reading the sphere packing \n");
	ReadSpherePacking(nspheres,cx,cy,cz,r);
	//.......................................................................

	//.......................................................................
	// Write out the files
	//.......................................................................
	int rank;
	char LocalRankString[8];
	char LocalRankFilename[40];
	for (k=0; k<nprocz; k++){
		for (j=0; j<nprocy; j++){
			for (i=0; i<nprocx; i++){
				//.......................................................................
				rank = k*nprocx*nprocy + j*nprocx + i;
				//.......................................................................
				sprintf(LocalRankString,"%05d",rank);
				sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
				//.......................................................................
//				printf("Assigning the local phase ID \n");
				AssignLocalSolidID(id,nspheres,cx,cy,cz,r,Lx,Ly,Lz,Nx,Ny,Nz,
								   i,j,k,nprocx,nprocy,nprocz);
				//.......................................................................
//				printf("Generating residual NWP \n");
				GenerateResidual(id,Nx,Ny,Nz,0.3);
				//.......................................................................
				WriteLocalSolidID(LocalRankFilename, id, N);
				//.......................................................................
				printf("Finished rank = %i \n",rank);
			}
		}
	}
	//.......................................................................

}
