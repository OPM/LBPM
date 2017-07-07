/*
 * Pre-processor to generate signed distance function from segmented data
 * segmented data should be stored in a raw binary file as 1-byte integer (type char)
 * will output distance functions for phases
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "common/Array.h"
#include "common/Domain.h"

int main(int argc, char **argv)
{
	
	int rank=0;

		bool MULTINPUT=false;

		int NWP,SOLID,rank_offset;
		SOLID=atoi(argv[1]);
		NWP=atoi(argv[2]);

		if (rank==0){
			printf("Solid Label: %i \n",SOLID);
			printf("NWP Label: %i \n",NWP);
		}
		if (argc > 3){
			rank_offset = atoi(argv[3]);
		}
		else{
			MULTINPUT=true;
			rank_offset=0;
		}

		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		int nprocs, nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
		double Lx, Ly, Lz;
		int Nx,Ny,Nz;
		int i,j,k,n;
		int BC=0;
		char Filename[40];
		int xStart,yStart,zStart;
		//  char fluidValue,solidValue;

		std::vector<char> solidValues;
		std::vector<char> nwpValues;
		std::string line;

		if (rank==0){
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

			ifstream image("Segmented.in");
			image >> Filename; 	// Name of data file containing segmented data
			image >> Nx;   		// size of the binary file
			image >> Ny;
			image >> Nz;
			image >> xStart;	// offset for the starting voxel
			image >> yStart;
			image >> zStart;

		}
		nprocs=nprocx*nprocy*nprocz;

		char *SegData = NULL;
		// Rank=0 reads the entire segmented data and distributes to worker processes
		if (rank==0){
			printf("Dimensions of segmented image: %i x %i x %i \n",Nx,Ny,Nz);
			SegData = new char[Nx*Ny*Nz];
			FILE *SEGDAT = fopen(Filename,"rb");
			if (SEGDAT==NULL) ERROR("Error reading segmented data");
			size_t ReadSeg;
			ReadSeg=fread(SegData,1,Nx*Ny*Nz,SEGDAT);
			if (ReadSeg != size_t(Nx*Ny*Nz)) printf("lbpm_segmented_decomp: Error reading segmented data (rank=%i)\n",rank);
			fclose(SEGDAT);
			printf("Read segmented data from %s \n",Filename);
		}

		// Get the rank info
		int N = (nx+2)*(ny+2)*(nz+2);

		// number of sites to use for periodic boundary condition transition zone
		int z_transition_size = (nprocz*nz - (Nz - zStart))/2;
		if (z_transition_size < 0) z_transition_size=0;

		char LocalRankFilename[40];
		char *loc_id;
		loc_id = new char [(nx+2)*(ny+2)*(nz+2)];

		// Set up the sub-domains
		if (rank==0){
			printf("Distributing subdomains across %i processors \n",nprocs);
			printf("Process grid: %i x %i x %i \n",nprocx,nprocy,nprocz);
			printf("Subdomain size: %i \n",N);
			printf("Size of transition region: %i \n", z_transition_size);

			for (int kp=0; kp<nprocz; kp++){
				for (int jp=0; jp<nprocy; jp++){
					for (int ip=0; ip<nprocx; ip++){
						// rank of the process that gets this subdomain
						int rnk = kp*nprocx*nprocy + jp*nprocx + ip;
						printf("extracting rank=%i \n",rnk);
						// Pack and send the subdomain for rnk
						for (k=0;k<nz+2;k++){
							for (j=0;j<ny+2;j++){
								for (i=0;i<nx+2;i++){
									int x = xStart + ip*nx + i-1;
									int y = yStart + jp*ny + j-1;
									//		int z = zStart + kp*nz + k-1;
									int z = zStart + kp*nz + k-1 - z_transition_size;
									if (z<zStart) 	z=zStart;
									if (!(z<Nz))	z=Nz-1;
									int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
									int nglobal = z*Nx*Ny+y*Nx+x;
									loc_id[nlocal] = SegData[nglobal];
								}
							}
						}
						// relabel the data
						for (k=0;k<nz+2;k++){
							for (j=0;j<ny+2;j++){
								for (i=0;i<nx+2;i++){
									n = k*(nx+2)*(ny+2) + j*(nx+2) + i;;
									if (loc_id[n]==char(SOLID))     loc_id[n] = 0;
									else if (loc_id[n]==char(NWP))  loc_id[n] = 1;
									else                     loc_id[n] = 2;

								}
							}
						}
						
						printf("writing file \n");
						// Write the data for this rank data 
						sprintf(LocalRankFilename,"ID.%05i",rnk+rank_offset);
						FILE *ID = fopen(LocalRankFilename,"wb");
						fwrite(loc_id,1,(nx+2)*(ny+2)*(nz+2),ID);
						fclose(ID);
					}
				}
			}
		}

		printf("Domain decomposition completed successfully \n");
		return 0;

}
