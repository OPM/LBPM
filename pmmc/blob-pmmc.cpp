#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

#ifndef ARRAY_H_INC
#include "Array.h"
#define ARRAY_H_INC
#endif

#include "PointList.h"

using namespace std;

int main(int argc, char *argv[])
{

	double d;
//	double *phase;
//	double *distance;

	int i,j,k,p,q,r, n, Nx,Ny,Nz;
	double vS,vF;
	double readval;

	vS = 0.f;
	vF = 0.f;
	printf("Solid surface: S(x) = %f \n",vS);
	printf("Fluid surface: F(x) = %f \n",vF);

	Nx = Ny = Nz = 64;
	printf("Domain size is: %i x %i x %i \n",Nx,Ny,Nz);


	DoubleArray F(Nx,Ny,Nz);
	DoubleArray S(Nx,Ny,Nz);

	printf("Reading distance from a file... \n");
	ifstream DISTANCE("Distance.in",ios::binary);
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				DISTANCE.read((char *) (&readval), sizeof(readval));
				n = k*Nx*Ny+j*Nx+i;
				S(i,j,k) = readval;
			}
		}
	}
	DISTANCE.close();

	printf("Reading phase from a file... \n");
	ifstream PHASE("Phase.in",ios::binary);
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				PHASE.read((char *) (&readval), sizeof(readval));
				n = k*Nx*Ny+j*Nx+i;
				F(i,j,k) = -readval;
				//phase[n] = d;

			}
		}
	}
	PHASE.close();

	IntArray indicator(Nx,Ny,Nz); // indicates if (i,j,k) has been assigned yet
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				// Initialize indicator fucntion to -1
				indicator(i,j,k) = -1;
			}
		}
	}

	printf("Execute blob identification algorithm... \n");

	/* ****************************************************************
				IDENTIFY ALL BLOBS: F > vF, S > vS
	****************************************************************** */
	// Find blob domains, number of blobs
	int nblobs = 0;					// number of blobs
	int ncubes = 0;					// total number of nodes in any blob
	int N = (Nx-1)*(Ny-1)*(Nz-1);		// total number of nodes
	IntArray blobs(3,N);	// store indices for blobs (cubes)
	IntArray temp(3,N);	// temporary storage array
	IntArray  b(50);		// number of nodes in each blob

	// Loop over z=0 first -> blobs attached to this end considered "connected" for LB simulation
	i=0;
	int number=0;
	for (k=0;k<1;k++){
		for (j=0;j<Ny;j++){
			if ( F(i,j,k) > vF ){
				if ( S(i,j,k) > vS ){
					// node i,j,k is in the porespace
					number = number+ComputeBlob(blobs,nblobs,ncubes,indicator,F,S,vF,vS,i,j,k,temp);
				}
			}
		}
	}
	// Specify the blob on the z axis
	if (ncubes > 0){
		b(nblobs) = number;
		printf("Number of blobs is: %i \n",nblobs);
		nblobs++;
	}
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=1;i<Nx;i++){
				if ( indicator(i,j,k) == -1 ){
					if ( F(i,j,k) > vF ){
						if ( S(i,j,k) > vS ){
							// node i,j,k is in the porespace
							b(nblobs) = ComputeBlob(blobs,nblobs,ncubes,indicator,F,S,vF,vS,i,j,k,temp);
							nblobs++;
						}
					}
				}
				// Otherwise, this point has already been assigned - ignore

				// Make sure list blob_nodes is large enough
				if ( nblobs > b.Length-1){
					printf("Increasing size of blob list \n");
					b = IncreaseSize(b,b.Length);
				}
			}
		}
	}
	// Go over all cubes again -> add any that do not contain nw phase
	bool add=1;			// Set to false if any corners contain nw-phase ( F > vF)
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	int count_in=0,count_out=0;
	int nodx,nody,nodz;
	for (k=0;k<Nz-1;k++){
		for (j=0;j<Ny-1;j++){
			for (i=0;i<Nx-1;i++){
				// Loop over cube corners
				add=1;				// initialize to true - add unless corner occupied by nw-phase
				for (p=0;p<8;p++){
					nodx=i+cube[p][0];
					nody=j+cube[p][1];
					nodz=k+cube[p][2];
					if ( indicator(nodx,nody,nodz) > -1 ){
						// corner occupied by nw-phase  -> do not add
						add = 0;
					}
				}
				if ( add == 1 ){
					blobs(0,ncubes) = i;
					blobs(1,ncubes) = j;
					blobs(2,ncubes) = k;
					ncubes++;
					count_in++;
				}
				else { count_out++; }
			}
		}
	}
	b(nblobs) = count_in;
	nblobs++;
	/* ****************************************************************
		     COMPUTE GRADIENT OF F,S    CURVATURE OF F
	****************************************************************** */

//	DTVectorField3D gradF = Gradient(Fluid);
//	DTVectorField3D gradS = Gradient(Solid);

	/* ****************************************************************
		     VARIABLES FOR THE PMMC ALGORITHM
	****************************************************************** */
	double area,awn,aws,ans,lwns;
	
	// Initialize arrays for return variables (areas, volumes, etc. for all blobs)
	DoubleArray nw_areas(nblobs);
	DoubleArray ns_areas(nblobs);
	DoubleArray ws_areas(nblobs);
	DoubleArray nws_length(nblobs);
	DoubleArray volume(nblobs);		// Volumes of every blob
	
	for (int a=0; a<nblobs; a++){
		nw_areas(a) = 0.f;
		ns_areas(a) = 0.f;
		ws_areas(a) = 0.f;
		nws_length(a) = 0.f;
		volume(a) = 0.f;
	}

	/* ****************************************************************
			RUN PMMC ON EACH BLOB
	****************************************************************** */
	printf("Running the PMMC Algorithm \n");
	printf("The number of blobs is %i \n",nblobs);
	int n_nw_tris_beg, n_ns_tris_beg, n_ws_tris_beg, n_nws_seg_beg;
	int start=0,finish;
	int a,c;
	int newton_steps = 0;
	double blob_volume;
	for (a=0;a<nblobs;a++){
		
		finish = start+b(a);

		ComputeAreasPMMC(blobs, start, finish, F, S, vF, vS, 
						 blob_volume, ans, aws, awn, lwns, Nx, Ny, Nz);
		
		start = finish;

		volume(a) = blob_volume;
		ws_areas(a) = aws;
		nw_areas(a) = awn;
		ns_areas(a) = ans;
		
		// Last "blob" is just the ws interface
		if (a+1 < nblobs){
			printf("-------------------------------- \n");
			printf("Blob id = %i \n", a);
			printf("Blob volume = %f \n", blob_volume);
			printf("Area wn = %f \n", nw_areas(a));
			printf("Area ns = %f \n", ns_areas(a));
			printf("-------------------------------- \n");
		}
		
	}  // End of the blob loop
	
	area = 0.0;
	for (a=0;a<nblobs;a++)	area+= ws_areas(a);
	printf("Area ws = %f \n", area);

	printf("Done. \n");
}