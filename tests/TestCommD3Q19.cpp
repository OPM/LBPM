//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/MPI_Helpers.h"

using namespace std;


std::shared_ptr<Database> loadInputs( int nprocs )
{
    auto db = std::make_shared<Database>();
    const int dim = 8;
    db->putScalar<int>( "BC", 0 );
    if ( nprocs == 1 ){
        db->putVector<int>( "nproc", { 1, 1, 1 } );
        db->putVector<int>( "n", { 3, 1, 1 } );
        db->putScalar<int>( "nspheres", 0 );
        db->putVector<double>( "L", { 1, 1, 1 } );
    } else if ( nprocs == 2 ) {
        db->putVector<int>( "nproc", { 2, 1, 1 } );
        db->putVector<int>( "n", { dim, dim, dim } );
        db->putScalar<int>( "nspheres", 0 );
        db->putVector<double>( "L", { 1, 1, 1 } );
    } else if ( nprocs == 4 ) {
        db->putVector<int>( "nproc", { 2, 2, 1 } );
        db->putVector<int>( "n", { dim, dim, dim } );
        db->putScalar<int>( "nspheres", 0 );
        db->putVector<double>( "L", { 1, 1, 1 } );
    } else if (nprocs==8){
        db->putVector<int>( "nproc", { 2, 2, 2 } );
        db->putVector<int>( "n", { dim, dim, dim } );
        db->putScalar<int>( "nspheres", 0 );
        db->putVector<double>( "L", { 1, 1, 1 } );
    }
    return db;
}


extern void GlobalFlipScaLBL_D3Q19_Init(double *dist, IntArray Map, int Np, int Nx, int Ny, int Nz,
								int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz)
{
	// Set of Discrete velocities for the D3Q19 Model
	static int D3Q19[18][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
	{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},
	{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};
	
	int q,i,j,k,n,N;
	int Cqx,Cqy,Cqz; // Discrete velocity
	int x,y,z;		// Global indices
	int xn,yn,zn; 	// Global indices of neighbor 
	int X,Y,Z;		// Global size
	int idx;
	X = Nx*nprocx;
	Y = Ny*nprocy;
	Z = Nz*nprocz;
    NULL_USE(Z);
	N = (Nx+2)*(Ny+2)*(Nz+2);	// size of the array including halo


	for (k=0; k<Nz; k++){ 
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){

				//n = (k+1)*(Nx+2)*(Ny+2) + (j+1)*(Nx+2) + i+1;
				idx=Map(i,j,k);

				if (idx > 0){
					// Get the 'global' index
					x = iproc*Nx+i;
					y = jproc*Ny+j;
					z = kproc*Nz+k;
					for (q=0; q<18; q++){
						// Odd distribution
						Cqx = D3Q19[q][0];
						Cqy = D3Q19[q][1];
						Cqz = D3Q19[q][2];
						xn = x - Cqx;
						yn = y - Cqy;
						zn = z - Cqz;
						xn=x; yn=y;zn=z;
						if (xn < 0) xn += nprocx*Nx;
						if (yn < 0) yn += nprocy*Ny;
						if (zn < 0) zn += nprocz*Nz;
						if (!(xn < nprocx*Nx)) xn -= nprocx*Nx;
						if (!(yn < nprocy*Ny)) yn -= nprocy*Ny;
						if (!(zn < nprocz*Nz)) zn -= nprocz*Nz;	

						dist[(q+1)*Np+idx] = (zn*X*Y+yn*X+xn) + (q+1)*0.01;

					}
				}
			}
		}
	}
}

extern int GlobalCheckDebugDist(double *dist, IntArray Map, int Np, int Nx, int Ny, int Nz, 
		int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int start, int finish)
{

	int returnValue = 0;
	int q,i,j,k,n,N,idx;
	int Cqx,Cqy,Cqz; // Discrete velocity
	int x,y,z;		// Global indices
	int xn,yn,zn; 	// Global indices of neighbor 
	int X,Y,Z;		// Global size
	X = Nx*nprocx;
	Y = Ny*nprocy;
	Z = Nz*nprocz;
	NULL_USE(Z);
	N = (Nx+2)*(Ny+2)*(Nz+2);	// size of the array including halo
	for (k=0; k<Nz; k++){ 
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){

				idx=Map(i,j,k);

				if (idx > start && idx< finish){
					// Get the 'global' index
					x = iproc*Nx+i;
					y = jproc*Ny+j;
					z = kproc*Nz+k;
					for (q=0; q<18; q++){

						if (dist[(q+1)*Np+idx] != (z*X*Y+y*X+x) + (q+1)*0.01){
							printf("******************************************\n");
							printf("error in distribution q = %i \n", (q+1));
							printf("i,j,k= %i, %i, %i \n", x,y,z);
							printf("dist = %5.2f \n", dist[(q+1)*Np+idx]);
							printf("n= %i \n",z*X*Y+y*X+x);
							returnValue++;
						}

					}
				}
			}
		}
	}
	return returnValue;
}

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
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	int check;

	{
		MPI_Request req1[18],req2[18];
		MPI_Status stat1[18],stat2[18];

		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Unit Test for D3Q19 MPI Communication	\n");
			printf("********************************************************\n");
		}

		// BGK Model parameters
		string FILENAME;
		unsigned int nBlocks, nthreads;
		int timestepMax, interval;
		double tau,Fx,Fy,Fz,tol;
		// Domain variables
		int i,j,k,n;

        // Load inputs
        auto db = loadInputs( nprocs );
	/*		auto filename = argv[1];
		auto input_db = std::make_shared<Database>( filename );
		auto db = input_db->getDatabase( "Domain" );
	*/
        int Nx = db->getVector<int>( "n" )[0];
        int Ny = db->getVector<int>( "n" )[1];
        int Nz = db->getVector<int>( "n" )[2];
		
		auto Dm  = std::shared_ptr<Domain>(new Domain(db,comm));      // full domain for analysis

		Nx += 2;
		Ny += 2;
		Nz += 2;
		int N = Nx*Ny*Nz;
		int dist_mem_size = N*sizeof(double);

		//.......................................................................
		// Assign the phase ID field
		//.......................................................................
		char LocalRankString[8];
		sprintf(LocalRankString,"%05d",rank);
		char LocalRankFilename[40];
		sprintf(LocalRankFilename,"ID.%05i",rank);

		char *id;
		id = new char[Nx*Ny*Nz];

		/*		if (rank==0) printf("Assigning phase ID from file \n");
		if (rank==0) printf("Initialize from segmented data: solid=0, NWP=1, WP=2 \n");
		FILE *IDFILE = fopen(LocalRankFilename,"rb");
		if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
		fread(id,1,N,IDFILE);
		fclose(IDFILE);
		*/
		// Setup the domain
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					id[n] = 1;
					Dm->id[n] = id[n];
				}
			}
		}
		Dm->CommInit();
		int iproc,jproc,kproc;
		int nprocx,nprocy,nprocz;
		iproc = Dm->iproc();
		jproc = Dm->jproc();
		kproc = Dm->kproc();
		nprocx = Dm->nprocx();
		nprocy = Dm->nprocy();
		nprocz = Dm->nprocz();

		if (rank==0){
			printf("********************************************************\n");
			printf("Sub-domain size = %i x %i x %i\n",Nz,Nz,Nz);
			printf("Parallel domain size = %i x %i x %i\n",Dm->nprocx(),Dm->nprocy(),Dm->nprocz());
			printf("********************************************************\n");
		}

		//.......................................................................
		// Compute the media porosity
		//.......................................................................
		double sum;
		double sum_local=0.0, porosity;
		char component = 0; // solid phase
		int Np=0;
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (id[n] == component){
						sum_local+=1.0;
					}
					else Np++;
				}
			}
		}
		MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
		double iVol_global=1.f/double((Nx-2)*(Ny-2)*(Nz-2)*nprocx*nprocy*nprocz);
	    	porosity = 1.0-sum*iVol_global;
		if (rank==0) printf("Media porosity = %f \n",porosity);
		//.......................................................................

		//...........................................................................
		MPI_Barrier(comm);
		if (rank == 0) cout << "Domain set." << endl;
		//...........................................................................

		//...........................................................................
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");
		// Create a communicator for the device (will use optimized layout)
		ScaLBL_Communicator ScaLBL_Comm(Dm);
		
		int Npad=(Np/16 + 2)*16;
		if (rank==0)    printf ("Set up memory efficient layout, %i | %i | %i \n", Np, Npad, N);
		auto neighborList= new int[18*Npad];
		IntArray Map(Nx,Ny,Nz);
		Map.fill(-2);		
		Np = ScaLBL_Comm.MemoryOptimizedLayoutAA(Map,neighborList,Dm->id,Np);
		MPI_Barrier(comm);
		int neighborSize=18*Np*sizeof(int);
		//......................device distributions.................................
		dist_mem_size = Np*sizeof(double);
		if (rank==0)	printf ("Allocating distributions \n");
		
		int *NeighborList;
		int *dvcMap;
		double *fq;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);
		//...........................................................................
		
		double *fq_host;
		fq_host = new double [19*Np];
		
		// Update GPU data structures
		if (rank==0)	printf ("Setting up device map and neighbor list \n");
		int *TmpMap;
		TmpMap=new int[Np];
		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int idx=Map(i,j,k);
					if (!(idx < 0))
						TmpMap[idx] = k*Nx*Ny+j*Nx+i;
				}
			}
		}
		ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
		ScaLBL_DeviceBarrier();
		delete [] TmpMap;

		//...........................................................................

		/*	// Write the communcation structure into a file for debugging
	char LocalCommFile[40];
	sprintf(LocalCommFile,"%s%s","Comm.",LocalRankString);
	FILE *CommFile;
	CommFile = fopen(LocalCommFile,"w");
	fprintf(CommFile,"rank=%d, ",rank);
	fprintf(CommFile,"i=%d,j=%d,k=%d :",iproc,jproc,kproc);
	fprintf(CommFile,"x=%d, ",rank_x);
	fprintf(CommFile,"X=%d, ",rank_X);
	fprintf(CommFile,"y=%d, ",rank_y);
	fprintf(CommFile,"Y=%d, ",rank_Y);
	fprintf(CommFile,"z=%d, ",rank_z);
	fprintf(CommFile,"Z=%d, ",rank_Z);
	fprintf(CommFile,"xy=%d, ",rank_xy);
	fprintf(CommFile,"XY=%d, ",rank_XY);
	fprintf(CommFile,"xY=%d, ",rank_xY);
	fprintf(CommFile,"Xy=%d, ",rank_Xy);
	fprintf(CommFile,"xz=%d, ",rank_xz);
	fprintf(CommFile,"XZ=%d, ",rank_XZ);
	fprintf(CommFile,"xZ=%d, ",rank_xZ);
	fprintf(CommFile,"Xz=%d, ",rank_Xz);
	fprintf(CommFile,"yz=%d, ",rank_yz);
	fprintf(CommFile,"YZ=%d, ",rank_YZ);
	fprintf(CommFile,"yZ=%d, ",rank_yZ);
	fprintf(CommFile,"Yz=%d, ",rank_Yz);
	fprintf(CommFile,"\n");
	fclose(CommFile);
		 */
		if (rank==0)	printf("Setting the distributions, size = : %i\n", Np);
		//...........................................................................
		GlobalFlipScaLBL_D3Q19_Init(fq_host, Map, Np, Nx-2, Ny-2, Nz-2, iproc,jproc,kproc,nprocx,nprocy,nprocz);
		ScaLBL_CopyToDevice(fq, fq_host, 19*dist_mem_size);
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);
		//*************************************************************************
		// First timestep
		ScaLBL_Comm.SendD3Q19AA(fq); //READ FROM NORMAL
	
		ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		
		// Second timestep
		ScaLBL_Comm.SendD3Q19AA(fq); //READ FROM NORMAL
		
		ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		
		//...........................................................................
		ScaLBL_CopyToHost(fq_host,fq,19*Np*sizeof(double));
		check =	GlobalCheckDebugDist(fq_host, Map, Np, Nx-2, Ny-2, Nz-2,iproc,jproc,kproc,nprocx,nprocy,nprocz,0,ScaLBL_Comm.next);
		//...........................................................................

		int timestep = 0;
		if (rank==0) printf("********************************************************\n");
		if (rank==0)	printf("No. of timesteps for timing: %i \n", 100);

		//.......create and start timer............
		double starttime,stoptime,cputime;
		MPI_Barrier(comm);
		starttime = MPI_Wtime();
		//.........................................


		//************ MAIN ITERATION LOOP (timing communications)***************************************/
		while (timestep < 100){

			// First timestep
			ScaLBL_Comm.SendD3Q19AA(fq); //READ FROM NORMAL
			
			ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
			
			// Second timestep
			ScaLBL_Comm.SendD3Q19AA(fq); //READ FROM NORMAL
			
			ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
			//*********************************************

			ScaLBL_DeviceBarrier();
			MPI_Barrier(comm);
			// Iteration completed!
			timestep++;
			//...................................................................
		}
		//************************************************************************/
		stoptime = MPI_Wtime();
		//	cout << "CPU time: " << (stoptime - starttime) << " seconds" << endl;
		cputime = stoptime - starttime;
		//	cout << "Lattice update rate: "<< double(Nx*Ny*Nz*timestep)/cputime/1000000 <<  " MLUPS" << endl;
		double MLUPS = double(Np)*double(timestep)/cputime*1e-6;
		if (rank==0) printf("********************************************************\n");
		if (rank==0) printf("CPU time = %f \n", cputime);
		if (rank==0) printf("Lattice update rate (per process)= %f MLUPS \n", MLUPS);
		MLUPS *= nprocs;
		if (rank==0) printf("Lattice update rate (process)= %f MLUPS \n", MLUPS);
		if (rank==0) printf("********************************************************\n");

		// Number of memory references from the swap algorithm (per timestep)
		// 18 reads and 18 writes for each lattice site
		double MemoryRefs = double(Np)*36;
		// number of memory references for the swap algorithm - GigaBytes / second
		if (rank==0) printf("DRAM bandwidth (per process)= %f GB/sec \n",MemoryRefs*8*double(timestep)*1e-9);
		// Report bandwidth in Gigabits per second
		// communication bandwidth includes both send and recieve
		if (rank==0) printf("Communication bandwidth (per process)= %f Gbit/sec \n",ScaLBL_Comm.CommunicationCount*64*timestep/1e9);
		if (rank==0) printf("Aggregated communication bandwidth = %f Gbit/sec \n",nprocs*ScaLBL_Comm.CommunicationCount*64*timestep/1e9);
	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************

	return check;
}
