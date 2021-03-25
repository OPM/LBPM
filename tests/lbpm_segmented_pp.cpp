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
#include "analysis/TwoPhase.h"
#include "analysis/distance.h"

inline void MeanFilter(DoubleArray &Mesh){
	for (int k=1; k<(int)Mesh.size(2)-1; k++){
		for (int j=1; j<(int)Mesh.size(1)-1; j++){
			for (int i=1; i<(int)Mesh.size(0)-1; i++){
				double sum;
				sum=Mesh(i,j,k)+Mesh(i+1,j,k)+Mesh(i-1,j,k)+Mesh(i,j+1,k)+Mesh(i,j-1,k)+
						+Mesh(i,j,k+1)+Mesh(i,j,k-1);
				Mesh(i,j,k) = sum/7.0;
			}
		}
	}
}


double ReadFromBlock( char *ID, int iproc, int jproc, int kproc, int Nx, int Ny, int Nz){

	int x,y,z,i,j,k;
	// Get the list of blocks that contain the local rank
	int b0x =iproc*(Nx-2)/1024;
	int b0y =jproc*(Ny-2)/1024;
	int b0z =kproc*(Nz-2)/1024;
	int Bx =(iproc+1)*(Nx-2)/1024+1;
	int By =(jproc+1)*(Ny-2)/1024+1;
	int Bz =(kproc+1)*(Nz-2)/1024+1;

	// arrays to hold the strings 
    char LocalRankFilename[40];
    char sx[2];
    char sy[2];
    char sz[2];
    
	// array to store ids read from block
	char *id;
	id = new char [1024*1024*1024];

	for (int bz=b0z; bz<Bz; bz++){
		for (int by=b0y; by<By; by++){
			for (int bx=b0x; bx<Bx; bx++){
				sprintf(sx,"%d",bx);
				sprintf(sy,"%d",by);
				sprintf(sz,"%d",bz);
				sprintf(LocalRankFilename,"%s%s%s%s%s%s%s","a2_x",sx,"_y",sy,"_z",sz,".gbd");
				//sprintf(LocalRankFilename,"%s%s%s%s%s%s%s","dis_",sx,"x_",sy,"y_",sz,"z.gbd");
				//printf("Reading file: %s \n", LocalRankFilename);
				//fflush(stdout);
				
				// Read the file
				FILE *IDFILE = fopen(LocalRankFilename,"rb");
				fread(id,1,1024*1024*1024,IDFILE);
				fclose(IDFILE);
				//printf("Loading data ... \n");
				
				// loop over global index
				for ( k=0;k<Nz;k++){
					for ( j=0;j<Ny;j++){
						for ( i=0;i<Nx;i++){
							//  id for this process relative to current block
							x = iproc*(Nx-2) + i - 1 - 1024*bx;
							y = jproc*(Ny-2) + j - 1 - 1024*by;
							z = kproc*(Nz-2) + k - 1 - 1024*bz;
							if (!(x<0) && !(y<0) && !(z<0) && x < 1024  && y < 1024  && z < 1024){
							  //char value=id[z*1024*1024 + y*1024 + x];
							  char value=0;
								  // check porespace nodes
								if (id[z*1024*1024 + y*1024 + x] < 215)	value = 1;
								// set the ID
								ID[k*Nx*Ny + j*Nx + i] = value;
							}
						}
					}
				}
			}
		}
	}

	// Compute porosity 
	uint64_t count=0;
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				if (ID[k*Nx*Ny+j*Nx+i] == 1) count++; 
			}
		}
	}
	double porosity = double(count)/double((Nx-2)*(Ny-2)*(Nz-2));
	return porosity;
	/*printf("Porosity is %f \n",double(count)/double(NX*NY*NZ));

	printf("Done getting data -- writing main file \n");
	FILE *OUT = fopen("FullData.raw","wb");
	fwrite(ID,1,N,OUT);
	fclose(OUT);
	printf("Completed! \n");
	 */
}



int main(int argc, char **argv)
{
	// Initialize MPI
	Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
       int rank = comm.getRank();
       int nprocs = comm.getSize();
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
		if (argc > 1) filename=argv[1];
		else ERROR("No input database provided\n");
		// read the input database 
		auto db = std::make_shared<Database>( filename );
		auto domain_db = db->getDatabase( "Domain" );

		// Read domain parameters
		auto L = domain_db->getVector<double>( "L" );
		auto size = domain_db->getVector<int>( "n" );
		auto nproc = domain_db->getVector<int>( "nproc" );
		auto ReadValues = domain_db->getVector<char>( "ReadValues" );
		auto WriteValues = domain_db->getVector<char>( "WriteValues" );

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
		std::shared_ptr<TwoPhase> Averages( new TwoPhase(Dm) );
		nx+=2; ny+=2; nz+=2;

		// Read the phase ID
		double porosity;
		sprintf(LocalRankFilename,"ID.%05i",rank);
		if (READ_FROM_BLOCK == 0){
			size_t readID;
			FILE *ID = fopen(LocalRankFilename,"rb");
			readID=fread(Dm->id,1,N,ID);
			if (readID != size_t(N)) printf("lbpm_segmented_pp: Error reading ID \n");
			fclose(ID);
		}
		else{
			// Read from the large block and write the local ID file
			if (rank==0) printf("Reading ID file from blocks \n");
			fflush(stdout);
			porosity = ReadFromBlock(Dm->id,Dm->iproc(),Dm->jproc(),Dm->kproc(),nx,ny,nz);
			
			MPI_Barrier(MPI_COMM_WORLD);
			if (rank==0) printf("Writing local ID files (poros=%f) \n",porosity);
			fflush(stdout);
			FILE *ID = fopen(LocalRankFilename,"wb");
			fwrite(Dm->id,1,N,ID);
			fclose(ID);
			if (rank==0) printf("Succeeded! \n");
			fflush(stdout);

		}

		// Initialize the domain and communication
		Array<char> id(nx,ny,nz);
		//TwoPhase Averages(Dm);
		//	DoubleArray Distance(nx,ny,nz);
		//	DoubleArray Phase(nx,ny,nz);

		int count = 0;
		// Solve for the position of the solid phase
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n = k*nx*ny+j*nx+i;
					// Initialize the solid phase
					if (Dm->id[n] > 0)	id(i,j,k) = 1;
					else		      	id(i,j,k) = 0;
				}
			}
		}
		// Initialize the signed distance function
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n=k*nx*ny+j*nx+i;
					// Initialize distance to +/- 1
					Averages->SDs(i,j,k) = 2.0*double(id(i,j,k))-1.0;
				}
			}
		}
		MeanFilter(Averages->SDs);

		if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
		CalcDist(Averages->SDs,id,*Dm);

		sprintf(LocalRankFilename,"SignDist.%05i",rank);
		FILE *DIST = fopen(LocalRankFilename,"wb");
		fwrite(Averages->SDs.data(),8,N,DIST);
		fclose(DIST);

	}
	comm.barrier();
	Utilities::shutdown();
	return 0;

}
