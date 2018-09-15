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
#include "common/MPI_Helpers.h"
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
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
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
		auto nprocx = nproc[0];
		auto nprocy = nproc[1];
		auto nprocz = nproc[2];
		auto Nx = SIZE[0];
		auto Ny = SIZE[1];
		auto Nz = SIZE[2];
		
		int i,j,k,n;

		char *SegData = NULL;
		// Rank=0 reads the entire segmented data and distributes to worker processes
		if (rank==0){
			printf("Dimensions of segmented image: %i x %i x %i \n",Nx,Ny,Nz);
			SegData = new char[Nx*Ny*Nz];
			FILE *SEGDAT = fopen(Filename.c_str(),"rb");
			if (SEGDAT==NULL) ERROR("Error reading segmented data");
			size_t ReadSeg;
			ReadSeg=fread(SegData,1,Nx*Ny*Nz,SEGDAT);
			if (ReadSeg != size_t(Nx*Ny*Nz)) printf("lbpm_segmented_decomp: Error reading segmented data (rank=%i)\n",rank);
			fclose(SEGDAT);
			printf("Read segmented data from %s \n",Filename.c_str());
		}
		MPI_Barrier(comm);

		// Get the rank info
		int N = (nx+2)*(ny+2)*(nz+2);
		
		std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
		for (k=0;k<nz+2;k++){
			for (j=0;j<ny+2;j++){
				for (i=0;i<nx+2;i++){
					n = k*(nx+2)*(ny+2)+j*(nx+2)+i;
					Dm->id[n] = 1;
				}
			}
		}
		Dm->CommInit();

		int z_transition_size = 0;
		int xStart = 0;
		int yStart = 0;
		int zStart = 0;
		// Set up the sub-domains
		if (rank==0){
			printf("Distributing subdomain across %i processors \n",nprocs);
			printf("Process grid: %i x %i x %i \n",Dm->nprocx(),Dm->nprocy(),Dm->nprocz());
			printf("Subdomain size: %i \n",N);
			//printf("Size of transition region: %i \n", z_transition_size);
			char *tmp;
			tmp = new char[N];
			for (int kp=0; kp<nprocz; kp++){
				for (int jp=0; jp<nprocy; jp++){
					for (int ip=0; ip<nprocx; ip++){
						// rank of the process that gets this subdomain
						int rnk = kp*Dm->nprocx()*Dm->nprocy() + jp*Dm->nprocx() + ip;
						// Pack and send the subdomain for rnk
						for (k=0;k<nz+2;k++){
							for (j=0;j<ny+2;j++){
								for (i=0;i<nx+2;i++){
									int x = xStart + ip*nx + i-1;
									int y = yStart + jp*ny + j-1;
									//		int z = zStart + kp*nz + k-1;
									int z = zStart + kp*nz + k-1 - z_transition_size;
									if (x<xStart) 	x=xStart;
									if (!(x<Nx))	x=Nx-1;
									if (y<yStart) 	y=yStart;
									if (!(y<Ny))	y=Ny-1;
									if (z<zStart) 	z=zStart;
									if (!(z<Nz))	z=Nz-1;
									int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
									int nglobal = z*Nx*Ny+y*Nx+x;
									tmp[nlocal] = SegData[nglobal];
								}
							}
						}
						if (rnk==0){
							for (k=0;k<nz+2;k++){
								for (j=0;j<ny+2;j++){
									for (i=0;i<nx+2;i++){
										int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
										Dm->id[nlocal] = tmp[nlocal];
									}
								}
							}
						}
						else{
							printf("Sending data to process %i \n", rnk);
							MPI_Send(tmp,N,MPI_CHAR,rnk,15,comm);
						}
					}
				}
			}
		}
		else{
			// Recieve the subdomain from rank = 0
			printf("Ready to recieve data %i at process %i \n", N,rank);
			MPI_Recv(Dm->id,N,MPI_CHAR,0,15,comm,MPI_STATUS_IGNORE);
		}
		MPI_Barrier(comm);
		
		// Compute the Minkowski functionals
		MPI_Barrier(comm);
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
	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;
}

