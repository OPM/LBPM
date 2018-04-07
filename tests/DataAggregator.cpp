#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

int main(int argc, char **argv){
	
	unsigned long int Bx,By,Bz;
	unsigned long int Nx,Ny,Nz;
	unsigned long int NX,NY,NZ;
	unsigned long int i,j,k;
	unsigned long int x,y,z;
	unsigned long int x0,y0,z0;
	unsigned long int N, N_full;
	unsigned long int Z0,Zf;

	//Read the block size
	Nx = atoi(argv[1]);
	Ny = atoi(argv[2]);
	Nz = atoi(argv[3]);
	x0 = atoi(argv[4]);
	y0 = atoi(argv[5]);
	z0 = atoi(argv[6]);
	NX = atoi(argv[7]);
	NY = atoi(argv[8]);
	NZ = atoi(argv[9]);

	//Bx = By = Bz = 9;
	//Nx = Ny = Nz = 1024;
	N = Nx*Ny*Nz;
	N_full=NX*NY*NZ;

	// compute the number of blocks to read
	Bx=NX/Nx+1;
	By=NY/Ny+1;
	Bz=NZ/Nz+1;

	printf("System size (read) is: %i x %i x %i \n",NX,NY,NZ);
	printf("Block size (read) is: %i x %i x %i \n",Nx,Ny,Nz);
	printf("Starting block (read) is: %i, %i, %i \n", x0,y0,z0);
	printf("First z slice (write) is: %i \n", Z0);
	printf("Last z slice (write) is: %i \n", Zf);
	
	// Filenames used
	//char LocalRankString[8];
	char LocalRankFilename[40];
	char sx[2];
	char sy[2];
	char sz[2];
	char tmpstr[10];
	//sprintf(LocalRankString,"%05d",rank);
	
	char *id;
	id = new char [N];
	char *ID;
	ID = new char [N_full];
	
	for (unsigned long int bz=z0/Nz; bz<Bz; bz++){
		for (unsigned long int by=y0/Ny; by<By; by++){
			for (unsigned long int bx=x0/Nx; bx<Bx; bx++){
				sprintf(sx,"%d",bx);
				sprintf(sy,"%d",by);
				sprintf(sz,"%d",bz);
				sprintf(LocalRankFilename,"%s%s%s%s%s%s%s","a2_x",sx,"_y",sy,"_z",sz,".gbd");
				printf("Reading file: %s \n", LocalRankFilename);
				// Read the file
				size_t readID;
				FILE *IDFILE = fopen(LocalRankFilename,"rb");
				readID=fread(id,1,N,IDFILE);
				fclose(IDFILE);
				
				// Unpack the data into the main array
				for ( k=0;k<Nz;k++){
					for ( j=0;j<Ny;j++){
						for ( i=0;i<Nx;i++){
							x = bx*Nx + i -x0;
							y = by*Ny + j -y0;
							z = bz*Nz + k -z0;
							if (!(x<0) && x<x0+NX && !(y<0) && y<y0+NY &&
							    !(z<0) && z<z0+NZ)
							  ID[z*NX*NY+y*NX+x] = ID[k*Nx*Ny+j*Nx+i];
						}
					}
				}			
			}
		}
	}

	// Compute porosity 
	unsigned long int count=0;
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
