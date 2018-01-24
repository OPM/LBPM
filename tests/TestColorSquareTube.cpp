
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
		// parallel domain size (# of sub-domains)
		int nprocx,nprocy,nprocz;
		int iproc,jproc,kproc;


		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Color Model: TestColor	\n");
			printf("********************************************************\n");
		}

		// BGK Model parameters
		string FILENAME;
		unsigned int nBlocks, nthreads;
		int timestepMax, interval;
		double Fx,Fy,Fz,tol;
		// Domain variables
		double Lx,Ly,Lz;
		int nspheres;
		int Nx,Ny,Nz;
		int i,j,k,n;
		int dim = 50;
		//if (rank == 0) printf("dim=%d\n",dim);
		int timestep = 1;
		int timesteps = 100;
		int centralNode = 2;

		double tauA = 1.0;
		double tauB = 1.0;
		double rhoA = 1.0;
		double rhoB = 1.0;
		double alpha = 0.001;
		double beta = 0.95;
		double tau = 1.0;
		double mu=(tau-0.5)/3.0;
		double rlx_setA=1.0/tau;
		double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);

		Fx = Fy = 0.f;
		Fz = 0.f;

		int typeBC;
		double din, dout, flux;
		double inletA,inletB,outletA,outletB;
		inletA=1.f;
		inletB=0.f;
		outletA=0.f;
		outletB=1.f;
		typeBC=4;
		flux = 10.f;
		dout=1.f;

		if (rank==0){
			//.......................................................................
			// Reading the domain information file
			//.......................................................................
			ifstream domain("Domain.in");
			if (domain.good()){
				domain >> nprocx;
				domain >> nprocy;
				domain >> nprocz;
				domain >> Nx;
				domain >> Ny;
				domain >> Nz;
				domain >> nspheres;
				domain >> Lx;
				domain >> Ly;
				domain >> Lz;
			}
			else if (nprocs==1){
				nprocx=nprocy=nprocz=1;
				Nx=3; Ny = 1;
				Nz = 1;
				nspheres=0;
				Lx=Ly=Lz=1;
			}
			else if (nprocs==2){
				nprocx=2; nprocy=1;
				nprocz=1;
				Nx=Ny=Nz=dim;
				Nx = dim; Ny = dim; Nz = dim;
				nspheres=0;
				Lx=Ly=Lz=1;
			}
			else if (nprocs==4){
				nprocx=nprocy=2;
				nprocz=1;
				Nx=Ny=Nz=dim;
				nspheres=0;
				Lx=Ly=Lz=1;
			}
			else if (nprocs==8){
				nprocx=nprocy=nprocz=2;
				Nx=Ny=Nz=dim;
				nspheres=0;
				Lx=Ly=Lz=1;
			}
			//.......................................................................
		}
		// **************************************************************
		// Broadcast simulation parameters from rank 0 to all other procs
		MPI_Barrier(comm);
		//.................................................
		MPI_Bcast(&Nx,1,MPI_INT,0,comm);
		MPI_Bcast(&Ny,1,MPI_INT,0,comm);
		MPI_Bcast(&Nz,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
		MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
		MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
		//.................................................
		MPI_Barrier(comm);
		// **************************************************************
		// **************************************************************

		if (nprocs != nprocx*nprocy*nprocz){
			printf("nprocx =  %i \n",nprocx);
			printf("nprocy =  %i \n",nprocy);
			printf("nprocz =  %i \n",nprocz);
			INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
		}

		if (rank==0){
			printf("********************************************************\n");
			printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
			printf("********************************************************\n");
		}
		MPI_Barrier(comm);

		double iVol_global = 1.0/Nx/Ny/Nz/nprocx/nprocy/nprocz;
		int BoundaryCondition=0;

		Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);
		
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


		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
                    n = k*Nx*Ny + j*Nx + i;
					Dm.id[n]=0;
				}
			}
		}

		kproc = rank/(nprocx*nprocy);
		jproc = (rank-nprocx*nprocy*kproc)/nprocx;
		iproc = rank-nprocx*nprocy*kproc-nprocx*jproc;
		printf("rank=%i, %i,%i,%i \n",rank,iproc,jproc,kproc);
		// Initialize a square tube
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny + j*Nx + i;
					int iglobal= i+(Nx-2)*iproc;
					int jglobal= j+(Ny-2)*jproc;
					int kglobal= k+(Nz-2)*kproc;

					// Initialize phase position field for parallel bubble test
					if (iglobal < 2)						Dm.id[n]=0;
					else if (iglobal > (Nx-2)*nprocx-2)	Dm.id[n]=0;
					else if (jglobal < 2)					Dm.id[n]=0;
					else if (jglobal > (Ny-2)*nprocy-2)	Dm.id[n]=0;
					else if (kglobal < 20)					Dm.id[n]=1;
					else									Dm.id[n]=2;
				}
			}
		}
		Dm.CommInit(comm);

		//.......................................................................
		// Compute the media porosity, assign phase labels and solid composition
		//.......................................................................
		double sum;
		double sum_local=0.0, porosity;
		int Np=0;  // number of local pore nodes
		double *PhaseLabel;
		PhaseLabel = new double[N];
		Dm.AssignComponentLabels(PhaseLabel);
		//.......................................................................
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
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");
		MPI_Barrier(comm);

		// Create a communicator for the device (will use optimized layout)
		ScaLBL_Communicator ScaLBL_Comm(Dm);
		//Create a second communicator based on the regular data layout
		ScaLBL_Communicator ScaLBL_Comm_Regular(Dm);

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
		if (rank==0)	printf ("Set up the neighborlist \n");

		int neighborSize=18*Np*sizeof(int);
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Np];

		ScaLBL_Comm.MemoryOptimizedLayoutAA(Map,neighborList,Dm.id,Np);
		MPI_Barrier(comm);

		//......................device distributions.................................
		int dist_mem_size = Np*sizeof(double);
		if (rank==0)	printf ("Allocating distributions \n");

		int *NeighborList;
		int *dvcMap;
		//		double *f_even,*f_odd;
		double *fq, *Aq, *Bq;
		double *Den, *Phi;
		double *ColorGrad;
		double *Vel;
		double *Pressure;
		
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);

		ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Aq, 7*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Bq, 7*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*Nx*Ny*Nz);		
		ScaLBL_AllocateDeviceMemory((void **) &Pressure, sizeof(double)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &Vel, 3*sizeof(double)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &ColorGrad, 3*sizeof(double)*Np);
		
		//...........................................................................
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
		//for (int idx=0; idx<Np; idx++) printf("Map=%i\n",TmpMap[idx]);

		ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
		ScaLBL_DeviceBarrier();
		delete [] TmpMap;
		
		// copy the neighbor list 
		ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
		// initialize phi based on PhaseLabel (include solid component labels)
		ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));
		//...........................................................................

		if (rank==0)	printf ("Initializing distributions \n");
		// Initialize the phase field and variables
		ScaLBL_D3Q19_Init(fq, Np);
		if (rank==0)	printf ("Initializing phase field \n");
		ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, Np);
		if (Dm.kproc==0){
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
		}
		if (Dm.kproc == nprocz-1){
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-1);
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-2);
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-3);
		}

		//************ MAIN ITERATION LOOP (timing communications)***************************************

		if (rank==0) printf("Beginning AA timesteps...\n");
		if (rank==0) printf("********************************************************\n");
		if (rank==0) printf("No. of timesteps for timing: %i \n", timesteps);

		//.......create and start timer............
		double starttime,stoptime,cputime;

		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		starttime = MPI_Wtime();
		//timesteps=20;
		//timestep=1;
		while (timestep < timesteps) {
			
			// *************ODD TIMESTEP*************
			// Compute the Phase indicator field
			// Read for Aq, Bq happens in this routine (requires communication)
			ScaLBL_Comm.BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
			ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm.next, Np, Np);
			ScaLBL_Comm.BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
			ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm.next, Np);

			// Halo exchange for phase field
			ScaLBL_Comm_Regular.SendHalo(Phi);
			ScaLBL_Comm_Regular.RecvHalo(Phi);
			
			// Perform the collision operation
			ScaLBL_Comm.SendD3Q19AA(fq); //READ FROM NORMAL
			ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm.next, Np, Np);
			ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
			// Set BCs
			if (typeBC > 0){
				ScaLBL_Comm.Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
				ScaLBL_Comm.Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
			}
			if (typeBC == 3){
				ScaLBL_Comm.D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			if (typeBC == 4){
				din = ScaLBL_Comm.D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm.next, Np);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
			timestep++;

			// *************EVEN TIMESTEP*************
			// Compute the Phase indicator field
			ScaLBL_Comm.BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
			ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm.next, Np, Np);
			ScaLBL_Comm.BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
			ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm.next, Np);

			// Halo exchange for phase field
			ScaLBL_Comm_Regular.SendHalo(Phi);
			ScaLBL_Comm_Regular.RecvHalo(Phi);

			// Perform the collision operation
			ScaLBL_Comm.SendD3Q19AA(fq); //READ FORM NORMAL
			ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz,  Nx, Nx*Ny, ScaLBL_Comm.next, Np, Np);
			ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
			// Set boundary conditions
			if (typeBC > 0){
				ScaLBL_Comm.Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
				ScaLBL_Comm.Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
			}
			if (typeBC == 3){
				ScaLBL_Comm.D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			else if (typeBC == 4){
				din = ScaLBL_Comm.D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm.next, Np);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
			timestep++;
			//************************************************************************

		}
		//************************************************************************
		stoptime = MPI_Wtime();
		//	cout << "CPU time: " << (stoptime - starttime) << " seconds" << endl;
		cputime = stoptime - starttime;
		//	cout << "Lattice update rate: "<< double(Nx*Ny*Nz*timestep)/cputime/1000000 <<  " MLUPS" << endl;
		double MLUPS = double(Np*timestep)/cputime/1000000;
		if (rank==0) printf("********************************************************\n");
		if (rank==0) printf("CPU time = %f \n", cputime);
		if (rank==0) printf("Lattice update rate (per process)= %f MLUPS \n", MLUPS);
		MLUPS *= nprocs;
		if (rank==0) printf("Lattice update rate (process)= %f MLUPS \n", MLUPS);
		if (rank==0) printf("********************************************************\n");

		// Number of memory references for color model
		double MemoryRefs = double(Np)*(77*8+(9+7+7)*4); // extra memory refs to read from neighborlist (every other timestep)
        // number of memory references for the swap algorithm - GigaBytes / second
        if (rank==0) printf("DRAM bandwidth (per process)= %f GB/sec \n",MemoryRefs*timestep/1e9/cputime);
        // Report bandwidth in Gigabits per second
        // communication bandwidth includes both send and recieve
        if (rank==0) printf("Communication bandwidth (per process)= %f Gbit/sec \n",ScaLBL_Comm.CommunicationCount*64*timestep/1e9/cputime);
        if (rank==0) printf("Aggregated communication bandwidth = %f Gbit/sec \n",nprocs*ScaLBL_Comm.CommunicationCount*64*timestep/1e9/cputime);
	
    	double *VEL;
    	VEL= new double [3*Np];
    	int SIZE=3*Np*sizeof(double);
    	ScaLBL_D3Q19_Momentum(fq,Vel,Np);
    	ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    	ScaLBL_CopyToHost(&VEL[0],&Vel[0],SIZE);

    	sum_local=0.f;
    	sum = 0.f;
    	for (k=1;k<Nz-1;k++){
    		for (j=1;j<Ny-1;j++){
    			for (i=1;i<Nx-1;i++){
    				n = k*Nx*Ny+j*Nx+i;
    				if (Dm.id[n] > 0){
    					int idx = Map(i,j,k);
    					sum_local+=VEL[2*Np+idx];
    				}
    			}
    		}
    	}
    	MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
    	double PoreVel = sum*iVol_global;
    	if (rank==0) printf("Average velocity = %f \n",PoreVel);
    	
    	if (rank==0){
    		printf("Printing inlet velocity for rank=0 \n");
    		k=1;
    		for (j=1;j<Ny-1;j++){
    			for (i=1;i<Nx-1;i++){ 
    				n = k*Nx*Ny+j*Nx+i;
    				if (Dm.id[n] > 0){
    					int idx = Map(i,j,k);
    					double vz = VEL[2*Np+idx];
    					printf("%f ",vz);
    				}
    			}
    			printf("\n");
    		}
    	}

    	double *PHASE;
    	PHASE= new double [Nx*Ny*Nz];
    	SIZE=Nx*Ny*Nz*sizeof(double);
    	ScaLBL_CopyToHost(&PHASE[0],&Phi[0],SIZE);
    	
    	FILE *OUTFILE;
		sprintf(LocalRankFilename,"Phase.%05i.raw",rank);
		OUTFILE = fopen(LocalRankFilename,"wb");
    	fwrite(PHASE,8,N,OUTFILE);
    	fclose(OUTFILE);
    	
    	double *DENA, *DENB, *TMPDAT;
    	SIZE=Np*sizeof(double);
    	TMPDAT = new double [Np];
    	DENA= new double [Nx*Ny*Nz];
    	DENB= new double [Nx*Ny*Nz];
    	ScaLBL_CopyToHost(&TMPDAT[0],&Den[0],SIZE);
    	ScaLBL_Comm.RegularLayout(Map,TMPDAT,DENA);
    	ScaLBL_CopyToHost(&TMPDAT[0],&Den[Np],SIZE);
    	ScaLBL_Comm.RegularLayout(Map,TMPDAT,DENB);

    	FILE *AFILE;
		sprintf(LocalRankFilename,"na.%05i.raw",rank);
    	AFILE = fopen(LocalRankFilename,"wb");
    	fwrite(DENA,8,N,AFILE);
    	fclose(AFILE);

    	FILE *BFILE;
		sprintf(LocalRankFilename,"nb.%05i.raw",rank);
    	BFILE = fopen(LocalRankFilename,"wb");
    	fwrite(DENB,8,N,BFILE);
    	fclose(BFILE);

    	double *CG;
    	CG= new double [3*Np];
    	ScaLBL_D3Q19_Gradient(dvcMap, Phi, ColorGrad, 0, Np, Np, Nx, Ny, Nz);
    	
    	ScaLBL_CopyToHost(&CG[0],&ColorGrad[0],3*SIZE);
    	for (int idx=0; idx<Np; idx++){
    		double C=CG[idx]*CG[idx]+CG[Np+idx]*CG[Np+idx]+CG[2*Np+idx]*CG[2*Np+idx];
    		TMPDAT[idx]=C;
    	}
    	ScaLBL_Comm.RegularLayout(Map,TMPDAT,DENB);
    	FILE *CGFILE;
		sprintf(LocalRankFilename,"cgrad.%05i.raw",rank);
		CGFILE = fopen(LocalRankFilename,"wb");
    	fwrite(DENB,8,N,CGFILE);
    	fclose(CGFILE);
    	
 
	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************

	return check;
}

