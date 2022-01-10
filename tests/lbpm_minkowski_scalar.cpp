// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>

#include "common/Array.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "common/MPI.h"
#include "IO/MeshDatabase.h"
#include "IO/Mesh.h"
#include "IO/Writer.h"
#include "IO/netcdf.h"
#include "analysis/analysis.h"
#include "analysis/filters.h"
#include "analysis/distance.h"
#include "analysis/Minkowski.h"

#include "ProfilerApp.h"

int main(int argc, char **argv)
{

	// Initialize MPI
        Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
        int rank = comm.getRank();
        //int nprocs = comm.getSize();
	{
		Utilities::setErrorHandlers();
		PROFILE_START("Main");

		//std::vector<std::string> filenames;
		if ( argc<2 ) {
			if ( rank == 0 ){
				printf("At least one filename must be specified\n");
			}
			return 1;
		}
		std::string filename = std::string(argv[1]);
		if ( rank == 0 ){
			printf("Input data file: %s\n",filename.c_str());
		}

		auto db = std::make_shared<Database>( filename );
		auto domain_db = db->getDatabase( "Domain" );

		// Read domain parameters
		auto Filename = domain_db->getScalar<std::string>( "Filename" );
		auto L = domain_db->getVector<double>( "L" );
		auto size = domain_db->getVector<int>( "n" );
		auto SIZE = domain_db->getVector<int>( "N" );
		auto nproc = domain_db->getVector<int>( "nproc" );
		auto ReadValues = domain_db->getVector<char>( "ReadValues" );
		auto WriteValues = domain_db->getVector<char>( "WriteValues" );
		auto nx = size[0];
		auto ny = size[1];
		auto nz = size[2];
		/*auto nprocx = nproc[0];
		auto nprocy = nproc[1];
		auto nprocz = nproc[2];
		auto Nx = SIZE[0];
		auto Ny = SIZE[1];
		auto Nz = SIZE[2];
		*/
		int i,j,k,n;

		std::shared_ptr<Domain> Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
		comm.barrier();
		Dm->CommInit();
		
		/* read the data */
		if (domain_db->keyExists( "Filename" )){
			auto Filename = domain_db->getScalar<std::string>( "Filename" );
			Dm->Decomp(Filename);
		}
		else{
			Dm->ReadIDs();
		}
		
		// Compute the Minkowski functionals
		comm.barrier();
		std::shared_ptr<Minkowski> Averages(new Minkowski(Dm));

		// Calculate the distance		
		// Initialize the domain and communication
		nx+=2; ny+=2; nz+=2;
		Array<char> id(nx,ny,nz);
		DoubleArray Distance(nx,ny,nz);
		
		//if (rank==0){
		//printf("ID: %i, %i, %i \n",Dm->Nx, Dm->Ny, Dm->Nz);
			//			printf("ID: %i, %i, %i \n",id.size(0),id.size(1),id.size(2));
		//	printf("SDn: %i, %i, %i \n",Averages->SDn.size(0),Averages->SDn.size(1),Averages->SDn.size(2));
		//}

		// Solve for the position of the solid phase
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n = k*nx*ny+j*nx+i;
					// Initialize the object
					if (Dm->id[n] == ReadValues[0])		id(i,j,k) = 1;
					else		      					id(i,j,k) = 0;
				}
			}
		}
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n=k*nx*ny+j*nx+i;
					// Initialize distance to +/- 1
					Distance(i,j,k) = 2.0*double(id(i,j,k))-1.0;
				}
			}
		}

		//std::array<bool> bc(3)={1,1,1};
		if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
		CalcDist(Distance,id,*Dm);

		if (rank==0) printf("Computing Minkowski functionals \n");
		Averages->ComputeScalar(Distance,0.f);
		Averages->PrintAll();
	}
	PROFILE_STOP("Main");
	PROFILE_SAVE("Minkowski",true);
        Utilities::shutdown();
	return 0;
}

