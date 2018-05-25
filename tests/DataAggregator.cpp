#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

using namespace std;

int main(int argc, char **argv){
	
  printf("Aggregating block data into single file \n");
	unsigned int Bx,By,Bz;
	uint64_t Nx,Ny,Nz;
	uint64_t N;
	uint64_t NX,NY,NZ;
	uint64_t i,j,k;
	uint64_t x,y,z;
	uint64_t x0,y0,z0;
	uint64_t N_full;

	//Read the block size
	if (argc>9){
	  printf("Input arguments accepted \n");
	Nx = atoi(argv[1]);
	Ny = atoi(argv[2]);
	Nz = atoi(argv[3]);
	x0 = atoi(argv[4]);
	y0 = atoi(argv[5]);
	z0 = atoi(argv[6]);
	NX = atol(argv[7]);
	NY = atol(argv[8]);
	NZ = atol(argv[9]);
	printf("Size %i X %i X %i \n",NX,NY,NZ);
	fflush(stdout);
	}
	else{
	  printf("setting defaults \n");
	  Nx=Ny=Nz=1024;
	  x0=y0=z0=0;
	  NX=NY=8640;
	  NZ=6480;
	}
	//Bx = By = Bz = 9;
	//Nx = Ny = Nz = 1024;

	// compute the number of blocks to read
	Bx=NX/Nx+1;
	By=NY/Ny+1;
	Bz=NZ/Nz+1;
	if (Bx>8) Bx=8;
	if (By>8) By=8;
	if (Bz>8) Bz=8;

	printf("System size (output) is: %i x %i x %i \n",NX,NY,NZ);
	printf("Block size (read) is: %i x %i x %i \n",Nx,Ny,Nz);
	printf("Starting location (read) is: %i, %i, %i \n", x0,y0,z0);
	printf("Block number (read): %i x %i x %i \n",Bx,By,Bz);
	fflush(stdout);
	
	// Filenames used
	//char LocalRankString[8];
	char LocalRankFilename[40];
	char sx[2];
	char sy[2];
	char sz[2];
	char tmpstr[10];

	//sprintf(LocalRankString,"%05d",rank);
	N = Nx*Ny*Nz;
	N_full=NX*NY*NZ;
	
	char *id;
	id = new char [N];
	char *ID;
	ID = new char [N_full];
	
	for (unsigned int  bz=0; bz<Bz; bz++){
		for (unsigned int by=0; by<By; by++){
			for (unsigned int bx=0; bx<Bx; bx++){
				sprintf(sx,"%d",bx);
				sprintf(sy,"%d",by);
				sprintf(LocalRankFilename,"%s%s%s%s%s%s%s","a2_x",sx,"_y",sy,"_z",sz,".gbd");
				printf("Reading file: %s \n", LocalRankFilename);
				fflush(stdout);

				// Read the file
				size_t readID;
				FILE *IDFILE = fopen(LocalRankFilename,"rb");
				readID=fread(id,1,N,IDFILE);
				fclose(IDFILE);
				printf("Loading data ... \n");
				// Unpack the data into the main array
				for ( k=0;k<Nz;k++){
					for ( j=0;j<Ny;j++){
						for ( i=0;i<Nx;i++){
							x = bx*Nx + i;
							y = by*Ny + j;
							z = bz*Nz + k;
							if ( x<NX && y<NY && z<NZ){
							  ID[z*NX*NY+y*NX+x] = id[k*Nx*Ny+j*Nx+i];
							}
						}
					}
				}			
			}
		}
	}

	// Compute porosity 
	uint64_t count=0;
	for (k=0; k<NZ; k++){
	  for (j=0; j<NY; j++){
	    for (i=0; i<NX; i++){
	      if (ID[k*NX*NY+j*NX+i] < 215) count++; 
	    }
	  }
	}
	printf("Porosity is %f \n",double(count)/double(NX*NY*NZ));

	printf("Done getting data -- writing main file \n");
	FILE *OUT = fopen("FullData.raw","wb");
	fwrite(ID,1,N_full,OUT);
	fclose(OUT);
	printf("Completed! \n");
}
