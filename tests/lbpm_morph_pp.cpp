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
    Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
	{
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
	        int n, nx, ny, nz;
		char LocalRankFilename[40];

		string filename;
		double SW;
		if (argc > 1){
			filename=argv[1];
			//Rcrit_new=0.f; 
			//SW=strtod(argv[2],NULL);
		}
		else ERROR("No input database provided\n");
		// read the input database 
		auto db = std::make_shared<Database>( filename );
		auto domain_db = db->getDatabase( "Domain" );

		// Read domain parameters
		auto size = domain_db->getVector<int>( "n" );
		auto nproc = domain_db->getVector<int>( "nproc" );
		auto ReadValues = domain_db->getVector<int>( "ReadValues" );
		auto WriteValues = domain_db->getVector<int>( "WriteValues" );
		SW = domain_db->getScalar<double>("Sw");
		auto READFILE = domain_db->getScalar<std::string>( "Filename" );

		// Generate the NWP configuration
		//if (rank==0) printf("Initializing morphological distribution with critical radius %f \n", Rcrit);
		if (rank==0) printf("Performing morphological imbibition with target saturation %f \n", SW);
		//	GenerateResidual(id,nx,ny,nz,Saturation);
		
		nx = size[0];
		ny = size[1];
		nz = size[2];

		int N = (nx+2)*(ny+2)*(nz+2);

		std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
		std::shared_ptr<Domain> Mask (new Domain(domain_db,comm));
		//		std::shared_ptr<Domain> Dm (new Domain(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC));
		for (n=0; n<N; n++) Dm->id[n]=1;
		Dm->CommInit();

		signed char *id;
		id = new signed char [N];
		signed char *id_connected;
		id_connected = new signed char [N];
		Mask->Decomp(READFILE);
		Mask->CommInit();

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
					id[n] = Mask->id[n];
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
					if (Mask->id[n] > 0){
						id_solid(i,j,k) = 1;
					}
					else	    
						id_solid(i,j,k) = 0;
				}
			}
		}
		// Initialize the signed distance function
		for (int k=0;k<nz;k++){
			for (int j=0;j<ny;j++){
				for (int i=0;i<nx;i++){
					// Initialize distance to +/- 1
					SignDist(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
				}
			}
		}

		if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
		CalcDist(SignDist,id_solid,*Dm);
		comm.barrier();
		
		// Extract only the connected part of NWP
		double vF=0.0; double vS=0.0;
		ComputeGlobalBlobIDs(nx-2,ny-2,nz-2,Dm->rank_info,phase,SignDist,vF,vS,phase_label,Dm->Comm);
		Dm->Comm.barrier();
			
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
		count_connected = Dm->Comm.sumReduce( count_connected );
		count_porespace = Dm->Comm.sumReduce( count_porespace );
		count_water = Dm->Comm.sumReduce( count_water );
		
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
		count_water = Dm->Comm.sumReduce( count_water );
		
		SW = double(count_water) / count_porespace;
		if(rank==0) printf("Final saturation: %f \n", SW);
		
		if (rank==0) printf("Writing ID file \n");
		sprintf(LocalRankFilename,"ID.%05i",rank);

		FILE *ID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,N,ID);
		fclose(ID);
		
		// write the geometry to a single file
		for (int k=0;k<nz;k++){
			for (int j=0;j<ny;j++){ 
				for (int i=0;i<nx;i++){
					int n = k*nx*ny+j*nx+i;
					Mask->id[n] = id[n];
				}
			}
		}
		comm.barrier();

        auto filename2 = READFILE + ".morph.raw";
		if (rank==0) printf("Writing file to: %s \n", filename2.c_str());
		Mask->AggregateLabels(filename2);
	}

    Utilities::shutdown();
}
