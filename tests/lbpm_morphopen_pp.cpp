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
#include "analysis/distance.h"
#include "analysis/morphology.h"

//*************************************************************************
// Morpohologica pre-processor
//   Initialize phase distribution using morphological approach
//   Signed distance function is used to determine fluid configuration
//*************************************************************************

int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	{
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
		double Lx, Ly, Lz;
		int i,j,k,n;
		int BC=0;
		//  char fluidValue,solidValue;
		int MAXTIME=1000;
		int READ_FROM_BLOCK=0;

		char LocalRankString[8];
		char LocalRankFilename[40];

		string filename;
		double Rcrit_new, SW;
		if (argc > 1){
			filename=argv[1];
			Rcrit_new=0.f; 
			//SW=strtod(argv[2],NULL);
		}
		else ERROR("No input database provided\n");
		// read the input database 
		auto db = std::make_shared<Database>( filename );
		auto domain_db = db->getDatabase( "Domain" );

		// Read domain parameters
		auto L = domain_db->getVector<double>( "L" );
		auto size = domain_db->getVector<int>( "n" );
		auto nproc = domain_db->getVector<int>( "nproc" );
		auto ReadValues = domain_db->getVector<int>( "ReadValues" );
		auto WriteValues = domain_db->getVector<int>( "WriteValues" );
		SW = domain_db->getScalar<double>("Sw");

		// Generate the NWP configuration
		//if (rank==0) printf("Initializing morphological distribution with critical radius %f \n", Rcrit);
		if (rank==0) printf("Performing morphological opening with target saturation %f \n", SW);
		//	GenerateResidual(id,nx,ny,nz,Saturation);
		
		nx = size[0];
		ny = size[1];
		nz = size[2];
		nprocx = nproc[0];
		nprocy = nproc[1];
		nprocz = nproc[2];

		int N = (nx+2)*(ny+2)*(nz+2);

		std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
		//		std::shared_ptr<Domain> Dm (new Domain(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC));
		for (n=0; n<N; n++) Dm->id[n]=1;
		Dm->CommInit();

		signed char *id;
		id = new signed char [N];
		sprintf(LocalRankFilename,"ID.%05i",rank);
		size_t readID;
		FILE *IDFILE = fopen(LocalRankFilename,"rb");
		if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
		readID=fread(id,1,N,IDFILE);
		if (readID != size_t(N)) printf("lbpm_morphopen_pp: Error reading ID (rank=%i) \n",rank);
		fclose(IDFILE);

		nx+=2; ny+=2; nz+=2;
		// Generate the signed distance map
		// Initialize the domain and communication
		Array<char> id_solid(nx,ny,nz);
		DoubleArray SignDist(nx,ny,nz);

		// Solve for the position of the solid phase
		for (int k=0;k<nz;k++){
			for (int j=0;j<ny;j++){
				for (int i=0;i<nx;i++){
					int n = k*nx*ny+j*nx+i;
					// Initialize the solid phase
					if (id[n] > 0)	id_solid(i,j,k) = 1;
					else	     	id_solid(i,j,k) = 0;
				}
			}
		}
		// Initialize the signed distance function
		for (int k=0;k<nz;k++){
			for (int j=0;j<ny;j++){
				for (int i=0;i<nx;i++){
					int n = k*nx*ny+j*nx+i;
					// Initialize distance to +/- 1
					SignDist(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
				}
			}
		}

		if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
		CalcDist(SignDist,id_solid,*Dm);

		MPI_Barrier(comm);

		// Run the morphological opening
		MorphOpen(SignDist, id, Dm, SW);
		
		// calculate distance to non-wetting fluid
		if (domain_db->keyExists( "HistoryLabels" )){
			if (rank==0) printf("Relabel solid components that touch fluid 1 \n");
			auto LabelList = domain_db->getVector<int>( "ComponentLabels" );
			auto HistoryLabels = domain_db->getVector<int>( "HistoryLabels" );
			size_t NLABELS=LabelList.size();
			if (rank==0){
				for (unsigned int idx=0; idx < NLABELS; idx++){
					signed char VALUE = LabelList[idx];
					signed char NEWVAL = HistoryLabels[idx];
					printf("    Relabel component %d as %d \n", VALUE, NEWVAL);
				}
			}
			for (int k=0;k<nz;k++){
				for (int j=0;j<ny;j++){
					for (int i=0;i<nx;i++){
						int n = k*nx*ny+j*nx+i;
						// Initialize the solid phase
						if (id[n] == 1)	id_solid(i,j,k) = 0;
						else	     	id_solid(i,j,k) = 1;
					}
				}
			}
			// Initialize the signed distance function
			for (int k=0;k<nz;k++){
				for (int j=0;j<ny;j++){
					for (int i=0;i<nx;i++){
						int n = k*nx*ny+j*nx+i;
						// Initialize distance to +/- 1
						SignDist(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
					}
				}
			}
			CalcDist(SignDist,id_solid,*Dm);
			// re-label IDs near the non-wetting fluid
			for (int k=0;k<nz;k++){
				for (int j=0;j<ny;j++){
					for (int i=0;i<nx;i++){
						int n = k*nx*ny+j*nx+i;
						signed char LOCVAL = id[n];
						for (unsigned int idx=0; idx < NLABELS; idx++){
							signed char VALUE=LabelList[idx];
							signed char NEWVALUE=HistoryLabels[idx];
							if (LOCVAL == VALUE){
								idx = NLABELS;
								if (SignDist(i,j,k) < 2.0){
									id[n] = NEWVALUE;
								}
							}
						}
					}
				}
			}
		}

		if (rank==0) printf("Writing ID file \n");
		sprintf(LocalRankFilename,"ID.%05i",rank);

		FILE *ID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,N,ID);
		fclose(ID);
	}

	MPI_Barrier(comm);
	MPI_Finalize();
}
