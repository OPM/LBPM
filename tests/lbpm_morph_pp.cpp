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
		if (rank==0) printf("Performing morphological imbibition with target saturation %f \n", SW);
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
		signed char *id_connected;
		id = new signed char [N];
		id_connected = new signed char [N];
		sprintf(LocalRankFilename,"ID.%05i",rank);
		size_t readID;
		FILE *IDFILE = fopen(LocalRankFilename,"rb");
		if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
		readID=fread(id,1,N,IDFILE);
		if (readID != size_t(N)) printf("lbpm_morph_pp: Error reading ID (rank=%i) \n",rank);
		fclose(IDFILE);

		nx+=2; ny+=2; nz+=2;
		// Generate the signed distance map
		// Initialize the domain and communication
		Array<char> id_solid(nx,ny,nz);
		Array<int> phase_label(nx,ny,nz);
		DoubleArray SignDist(nx,ny,nz);
		DoubleArray phase(nx,ny,nz);
		
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					n=k*nx*ny+j*nx+i;
					if (id[n] == 1){
						phase(i,j,k) = 1.0;
					}
					else
						phase(i,j,k) = -1.0;
				}
			}
		}
		
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
		
		// Extract only the connected part of NWP
		BlobIDstruct new_index;
		double vF=0.0; double vS=0.0;
		ComputeGlobalBlobIDs(nx-2,ny-2,nz-2,Dm->rank_info,phase,SignDist,vF,vS,phase_label,Dm->Comm);
		MPI_Barrier(Dm->Comm);
			
		int count_connected=0;
		int count_porespace=0;
		int count_water=0;
		for (int k=1; k<nz-1; k++){
			for (int j=1; j<ny-1; j++){
				for (int i=1; i<nx-1; i++){
					n=k*nx*ny+j*nx+i;
					// only apply opening to connected component 
					if ( phase_label(i,j,k) == 0){
						count_connected++;
					}
					if (id[n] > 0){
						count_porespace++;
					}
					if (id[n] == 2){
						count_water++;
					}
				}
			}
		}
		count_connected=sumReduce( Dm->Comm, count_connected);
		count_porespace=sumReduce( Dm->Comm, count_porespace);
		count_water=sumReduce( Dm->Comm, count_water);
		
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					n=k*nx*ny+j*nx+i;
					// only apply opening to connected component 
					if ( phase_label(i,j,k) == 0){
						id_solid(i,j,k) = 1;
						id_connected[n] = 2;
						id[n] = 2;
					}
					else{
						id_solid(i,j,k) = 0;
						id_connected[n] = 0;
					}
				}
			}
		}
		CalcDist(SignDist,id_solid,*Dm);

		// target water increase in voxels, normalized by connected volume
		double St = (SW*count_porespace - count_water)/count_porespace;  
		
		
		signed char water=2;
		signed char notwater=1;
		// Run the morphological opening
		if (St > 0.0)
			MorphOpen(SignDist, id_connected, Dm, St, water, notwater);
		else {
			if(rank==0) printf("Initial condition satisfies condition for saturation target \n");
		}
		
		// re-label 
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					n=k*nx*ny+j*nx+i;
					// only apply opening to connected component 
					if ( id_connected[n] == 1){
						id[n] = 1;
					}
				}
			}
		}
		
		count_water=0;
		for (int k=1; k<nz-1; k++){
			for (int j=1; j<ny-1; j++){
				for (int i=1; i<nx-1; i++){
					n=k*nx*ny+j*nx+i;
					if (id[n] == 2){
						count_water++;
					}
				}
			}
		}
		count_water=sumReduce( Dm->Comm, count_water);
		
		SW = double(count_water) / count_porespace;
		if(rank==0) printf("Final saturation: %f \n", SW);
		
		if (rank==0) printf("Writing ID file \n");
		sprintf(LocalRankFilename,"ID.%05i",rank);

		FILE *ID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,N,ID);
		fclose(ID);
	}

	MPI_Barrier(comm);
	MPI_Finalize();
}
