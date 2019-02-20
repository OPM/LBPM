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

//*************************************************************************
// Morpohological drainage pre-processor 
//   Generate states on primary drainage using morphological approach
//   Signed distance function is used to determine fluid configuration
//*************************************************************************
inline void PackID(int *list, int count, char *sendbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This packs up the values that need to be sent from one processor to another
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = ID[n];
	}
}
//***************************************************************************************

inline void UnpackID(int *list, int count, char *recvbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This unpacks the values once they have been recieved from neighbors
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		ID[n] = recvbuf[idx];
	}
}

//***************************************************************************************


int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);

	//.......................................................................
	// Reading the domain information file
	//.......................................................................
	int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
	double Lx, Ly, Lz;
	int i,j,k,n;
	int BC=0;

	string filename;
	if (argc > 1){
		filename=argv[1];
	}
	else ERROR("No input database provided\n");
	// read the input database 
	auto db = std::make_shared<Database>( filename );
	auto domain_db = db->getDatabase( "Domain" );

	// Read domain parameters
	auto L = domain_db->getVector<double>( "L" );
	auto size = domain_db->getVector<int>( "n" );
	auto nproc = domain_db->getVector<int>( "nproc" );
	double SW = domain_db->getScalar<double>( "Sw" );
	
	nx = size[0];
	ny = size[1];
	nz = size[2];
	nprocx = nproc[0];
	nprocy = nproc[1];
	nprocz = nproc[2];

	// Check that the number of processors >= the number of ranks
	if ( rank==0 ) {
		printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
		printf("Number of MPI ranks used: %i \n", nprocs);
		printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
	}
	if ( nprocs < nprocx*nprocy*nprocz ){
		ERROR("Insufficient number of processors");
	}

	char LocalRankFilename[40];

	int BoundaryCondition=1;
	Domain Dm(domain_db,comm);
	
	nx+=2; ny+=2; nz+=2;
	int N = nx*ny*nz;
	char *id;
	id = new char[N];


	// Define communication sub-domain -- everywhere
	for (int k=0; k<nz; k++){
		for (int j=0; j<ny; j++){
			for (int i=0; i<nx; i++){
				n = k*nx*ny+j*nx+i;
				Dm.id[n] = 1;
			}
		}
	}

	Dm.CommInit();

	sprintf(LocalRankFilename,"ID.%05i",rank);
	size_t readID;
	FILE *IDFILE = fopen(LocalRankFilename,"rb");
	if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
	readID=fread(id,1,N,IDFILE);
	if (readID != size_t(N)) printf("lbpm_morphdrain_pp: Error reading ID (rank=%i) \n",rank);
	fclose(IDFILE);

	int xdim,ydim,zdim;
	xdim=Dm.Nx-2;
	ydim=Dm.Ny-2;
	zdim=Dm.Nz-2;
	fillHalo<double> fillData(Dm.Comm, Dm.rank_info,{xdim,ydim,zdim},{1,1,1},0,1);
	
	// Generate the signed distance map
	// Initialize the domain and communication
	Array<char> id_solid(nx,ny,nz);
	DoubleArray SignDist(nx,ny,nz);
	DoubleArray phase(nx,ny,nz);
	IntArray phase_label(nx,ny,nz);

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
	CalcDist(SignDist,id_solid,Dm);

	// Run the morphological opening
	MorphOpen(SignDist, id, Dm, SW);

	for (int k=0;k<nz;k++){
		for (int j=0;j<ny;j++){
			for (int i=0;i<nx;i++){
				int n = k*nx*ny+j*nx+i;
				if (id[n] == 1){
					phase(i,j,k) = 1.0;
				}
				else
					phase(i,j,k) = -1.0;
			}
		}
	}
	
	// Extract only the connected part
	BlobIDstruct new_index;
	double vF=0.0; double vS=0.0;
	ComputeGlobalBlobIDs(nx-2,ny-2,nz-2,Dm.rank_info,phase,SignDist,vF,vS,phase_label,Dm.comm);
	MPI_Barrier(comm);
	
	for (int k=0;k<nz;k++){
		for (int j=0;j<ny;j++){
			for (int i=0;i<nx;i++){
				int n = k*nx*ny+j*nx+i;
				if (id[n] == 1 && phase_label(i,j,k) > 1){
					id[n] = 2;
				}
			}
		}
	}

	// calculate distance to non-wetting fluid
	if (domain_db->keyExists( "HistoryLabels" )){
		if (rank==0) printf("Relabel solid components that touch fluid 1 \n");
		auto LabelList = domain_db->getVector<char>( "ComponentLabels" );
		auto HistoryLabels = domain_db->getVector<char>( "HistoryLabels" );
		size_t NLABELS=LabelList.size();
		if (rank==0){
			for (unsigned int idx=0; idx < NLABELS; idx++){
				char VALUE = LabelList[idx];
				char NEWVAL = HistoryLabels[idx];
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
						char VALUE=LabelList[idx];
						char NEWVALUE=HistoryLabels[idx];
						if (LOCVAL == VALUE){
							idx = NLABELS;
							if (SignDist(i,j,k) < 1.0){
								id[n] = NEWVALUE;
							}
						}
					}
				}
			}
		}
	}

	sprintf(LocalRankFilename,"ID.%05i",rank);
	FILE *ID = fopen(LocalRankFilename,"wb");
	fwrite(id,1,N,ID);
	fclose(ID);

	MPI_Barrier(comm);
	MPI_Finalize();
}
