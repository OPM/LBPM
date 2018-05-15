
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/MPI_Helpers.h"



std::shared_ptr<Database> loadInputs( int nprocs )
{
    const int dim = 12;
	std::shared_ptr<Database> db;
    if ( exists( "Domain.in" ) ) {
    	db = std::make_shared<Database>( "Domain.in" );
	} else if (nprocs==1) {
        db->putVector<int>( "nproc", { 1, 1, 1 } );
        db->putVector<int>( "n", { dim, dim, dim } );
        db->putScalar<int>( "nspheres", 0 );
        db->putVector<double>( "L", { 1, 1, 1 } );
	} else if (nprocs==2) {
        db->putVector<int>( "nproc", { 1, 1, 1 } );
        db->putVector<int>( "n", { dim/2, dim/2, dim/2 } );
        db->putScalar<int>( "nspheres", 0 );
        db->putVector<double>( "L", { 1, 1, 1 } );
    }
    db->putScalar<int>( "BC", 0 );
    return db;
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
		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Unit Test: TestPoiseuille	\n");
			printf("********************************************************\n");
		}

		// BGK Model parameters
		string FILENAME;
		unsigned int nBlocks, nthreads;
		int timestepMax, interval;
		double tau,Fx,Fy,Fz,tol;
		// Domain variables
		int i,j,k,n;
		int timestep = 0;


		tau = 1.0;
		double mu=(tau-0.5)/3.0;
		double rlx_setA=1.0/tau;
		double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
		Fx = 0; Fy = 0;
		Fz = 1e-3; //1.f; // 1e-3;

        // Load inputs
        auto db = loadInputs( nprocs );
        int Nx = db->getVector<int>( "n" )[0];
        int Ny = db->getVector<int>( "n" )[1];
        int Nz = db->getVector<int>( "n" )[2];
        int nprocx = db->getVector<int>( "nproc" )[0];
        int nprocy = db->getVector<int>( "nproc" )[1];
        int nprocz = db->getVector<int>( "nproc" )[2];


		if (rank==0){
			printf("********************************************************\n");
			printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
			printf("********************************************************\n");
		}

		MPI_Barrier(comm);
		int kproc = rank/(nprocx*nprocy);
		int jproc = (rank-nprocx*nprocy*kproc)/nprocx;
		int iproc = rank-nprocx*nprocy*kproc-nprocz*jproc;

		if (rank == 0) {
			printf("i,j,k proc=%d %d %d \n",iproc,jproc,kproc);
		}
		MPI_Barrier(comm);
		if (rank == 1){
			printf("i,j,k proc=%d %d %d \n",iproc,jproc,kproc);
			printf("\n\n");
		}

		double iVol_global = 1.0/Nx/Ny/Nz/nprocx/nprocy/nprocz;

		Domain Dm(db);

		Nx += 2;
		Ny += 2;
		Nz += 2;
		int N = Nx*Ny*Nz;

		//.......................................................................
		// Assign the phase ID field
		//.......................................................................
		char LocalRankString[8];
		sprintf(LocalRankString,"%05d",rank);
		char LocalRankFilename[40];
		sprintf(LocalRankFilename,"ID.%05i",rank);
		/*
		FILE *IDFILE = fopen(LocalRankFilename,"rb");
		if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
		fread(Dm.id,1,N,IDFILE);
		fclose(IDFILE);
		 */

		// initialize empty domain
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (i<2) Dm.id[n] = 0;
					else if (i>Nx-3) Dm.id[n] = 0;
					else Dm.id[n]=1;
				}
			}
		}

		Dm.CommInit(comm);
		MPI_Barrier(comm);

		//.......................................................................
		// Compute the media porosity
		//.......................................................................
		double sum;
		double sum_local=0.0, porosity;
		int Np=0;  // number of local pore nodes
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (Dm.id[n] > 0){
						sum_local+=1.0;
						Np++;
					}
				}
			}
		}
		MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
		porosity = sum*iVol_global;
		if (rank==0) printf("Media porosity = %f \n",porosity);

		MPI_Barrier(comm);
		if (rank == 0) cout << "Domain set." << endl;
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");

		// Create a communicator for the device
		ScaLBL_Communicator ScaLBL_Comm(Dm);

		//...........device phase ID.................................................
		if (rank==0)	printf ("Copying phase ID to device \n");
		char *ID;
		ScaLBL_AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
		// Copy to the device
		ScaLBL_CopyToDevice(ID, Dm.id, N);
		//...........................................................................

		if (rank==0){
			printf("Total domain size = %i \n",N);
			printf("Reduced domain size = %i \n",Np);
		}


		// LBM variables
		if (rank==0)	printf ("Allocating distributions \n");
		if (rank==0)	printf ("Set up the neighborlist \n");
		int Npad=Np+32;
		int neighborSize=18*Npad*sizeof(int);
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Npad];
		Np = ScaLBL_Comm.MemoryOptimizedLayoutAA(Map,neighborList,Dm.id,Np);
		MPI_Barrier(comm);

		//......................device distributions.................................
		int dist_mem_size = Np*sizeof(double);

		int *NeighborList;
		//		double *f_even,*f_odd;
		double * dist;
		double * Velocity;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &dist, 19*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
		ScaLBL_CopyToDevice(NeighborList,     neighborList, neighborSize);
		//...........................................................................

		/*
		 *  AA Algorithm begins here
		 *
		 */
		ScaLBL_D3Q19_Init(dist, Np);

		//.......create and start timer............
		double starttime,stoptime,cputime;

		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		starttime = MPI_Wtime();

		/************ MAIN ITERATION LOOP (timing communications)***************************************/
		//		ScaLBL_Comm.SendD3Q19(dist, &dist[10*Np]);
		//		ScaLBL_Comm.RecvD3Q19(dist, &dist[10*Np]);
		//		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		//

		if (rank==0) printf("Beginning AA timesteps...\n");
		if (rank==0) printf("********************************************************\n");

		while (timestep < 2000) {

			ScaLBL_Comm.SendD3Q19AA(dist); //READ FROM NORMAL
			ScaLBL_D3Q19_AAodd_MRT(NeighborList, dist,  ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
			ScaLBL_Comm.RecvD3Q19AA(dist); //WRITE INTO OPPOSITE
			ScaLBL_D3Q19_AAodd_MRT(NeighborList, dist, 0, ScaLBL_Comm.next, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
			timestep++;

			ScaLBL_Comm.SendD3Q19AA(dist); //READ FORM NORMAL
			ScaLBL_D3Q19_AAeven_MRT(dist, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
			ScaLBL_Comm.RecvD3Q19AA(dist); //WRITE INTO OPPOSITE
			ScaLBL_D3Q19_AAeven_MRT(dist, 0, ScaLBL_Comm.next, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
			timestep++;
			//************************************************************************/

		}
		//************************************************************************/
		stoptime = MPI_Wtime();
		//	cout << "CPU time: " << (stoptime - starttime) << " seconds" << endl;
		cputime = stoptime - starttime;
		//	cout << "Lattice update rate: "<< double(Nx*Ny*Nz*timestep)/cputime/1000000 <<  " MLUPS" << endl;
		double MLUPS = double(Np*timestep)/cputime/1000000;
	//	if (rank==0) printf("********************************************************\n");
	//	if (rank==0) printf("CPU time = %f \n", cputime);
	//	if (rank==0) printf("Lattice update rate (per process)= %f MLUPS \n", MLUPS);
		MLUPS *= nprocs;
	//	if (rank==0) printf("Lattice update rate (process)= %f MLUPS \n", MLUPS);
	//	if (rank==0) printf("********************************************************\n");

		// Number of memory references from the swap algorithm (per timestep)
		// 18 reads and 18 writes for each lattice site
		double MemoryRefs = Np*38;
		// number of memory references for the swap algorithm - GigaBytes / second
	//	if (rank==0) printf("DRAM bandwidth (per process)= %f GB/sec \n",MemoryRefs*8*timestep/1e9/cputime);
		// Report bandwidth in Gigabits per second
		// communication bandwidth includes both send and recieve
		//if (rank==0) printf("Communication bandwidth (per process)= %f Gbit/sec \n",ScaLBL_Comm.CommunicationCount*64*timestep/1e9/cputime);
	//	if (rank==0) printf("Aggregated communication bandwidth = %f Gbit/sec \n",nprocs*ScaLBL_Comm.CommunicationCount*64*timestep/1e9/cputime);

		double *Vz;
		Vz= new double [Np];
		int SIZE=Np*sizeof(double);
    	ScaLBL_D3Q19_Momentum(dist,Velocity, Np);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		ScaLBL_CopyToHost(&Vz[0],&Velocity[2*Np],SIZE);

		if (rank == 0) printf("Force: %f,%f,%f \n",Fx,Fy,Fz);

		double vz;
		double W = 1.f*Nx-4.f;
		j=Ny/2; k=Nz/2;
		if (rank == 0) printf("Channel width=%f \n",W);
		if (rank == 0) printf("ID flag vz       analytical\n");

		MPI_Barrier(comm);
		if (rank == 0) {
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				printf("%i ",Dm.id[n]);

				n = Map(i,j,k);
				//printf("%i,%i,%i; %i :",i,j,k,n);
				if (n<0) {vz =0.f; printf(" b    "); }
				else { vz=Vz[n]; printf(" a    "); }
				printf("%f ",vz);
				//Analytical solution
				double x=1.f*i-1.5;
				if (n<0) vz=0.f;
				else vz=Fz*x*(W-x)/(2.f*mu);
				printf("%f\n",vz);
			}
			printf("\n");
		}


		if (rank == 1) {
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				printf("%i ",Dm.id[n]);

				n = Map(i,j,k);
				//printf("%i,%i,%i; %i :",i,j,k,n);
				if (n<0) {vz =0.f; printf(" b    "); }
				else { vz=Vz[n]; printf(" a    "); }
				printf("%f ",vz);
				//Analytical solution
				double x=1.f*i-1.5;
				if (n<0) vz=0.f;
				else vz=Fz*x*(W-x)/(2.f*mu);
				printf("%f\n",vz);
			}
			printf("\n");
		}



	}

	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************

	return check;
}
