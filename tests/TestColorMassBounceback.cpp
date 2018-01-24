
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
		int dim = 3;
		//if (rank == 0) printf("dim=%d\n",dim);
		int timestep = 0;
		int timesteps = 100;
		int centralNode = 2;

		double tauA = 1.0;
		double tauB = 1.0;
		double rhoA = 1.0;
		double rhoB = 1.0;
		double alpha = 0.005;
		double beta = 0.95;
		
		double tau = 1.0;
		double mu=(tau-0.5)/3.0;
		double rlx_setA=1.0/tau;
		double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);

		Fx = Fy = 0.f;
		Fz = 0.f;

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
				Nx=Ny=Nz=3;
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

		int Np=0;  // number of local pore nodes
		double *PhaseLabel;
		PhaseLabel = new double[N];
		//.......................................................................
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					Dm.id[n]=0;
				}
			}
		}
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					Dm.id[n]=1;
					Np++;
					// constant color
					PhaseLabel[n]= -1.0;
				}
			}
		}
		Dm.CommInit(comm);
		MPI_Barrier(comm);
		if (rank == 0) cout << "Domain set." << endl;
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");

		//Create a second communicator based on the regular data layout
		ScaLBL_Communicator ScaLBL_Comm_Regular(Dm);
		ScaLBL_Communicator ScaLBL_Comm(Dm);

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
		TmpMap=new int[Np*sizeof(int)];
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

		// copy the neighbor list 
		ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
		//...........................................................................

		// Distributions / densities for checking
		double nA,nB;
		double *DIST;
		DIST= new double [7*Np];
		double *DENSITY;
		DENSITY= new double [2*Np];
		int SIZE;
		
		int errc_odd_a=0;
		int errc_even_a=0;
		int errc_odd_b=0;
		int errc_even_b=0;
		
		//*******************Component A*******************
		// initialize phi based on PhaseLabel (include solid component labels)
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					// constant color
					PhaseLabel[n]= 1.0;
				}
			}
		}
		ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));

		if (rank==0)	printf ("Initializing distributions \n");
		// Initialize the phase field and variables
		ScaLBL_D3Q19_Init(fq, Np);
		if (rank==0)	printf ("Initializing phase field \n");
		ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, Np);

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

		ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm.next, Np);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		timestep++;

		printf("Check after odd time \n");
		SIZE=2*Np*sizeof(double);
		ScaLBL_CopyToHost(&DENSITY[0],&Den[0],SIZE);

    	// Check the distributions
		SIZE=7*Np*sizeof(double);
		ScaLBL_CopyToHost(&DIST[0],&Aq[0],SIZE);

		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (Dm.id[n] > 0){
						int idx = Map(i,j,k);
    					nA=DENSITY[idx];
    					nB=DENSITY[Np+idx];
    					//printf("i,j,k=%i,%i,%i \n",i,j,k);
    					//printf("   nA=%f, nB=%f \n",nA,nB);
    					double val=DIST[idx];
    					double error = fabs(val - 0.3333333333333333*nA);
    					if (error > 1.0e-12) {
    						printf("   q=0, Aq=%f \n",val);
    						errc_odd_b++;
    					}
    					for (int q=1; q<7; q++){
    						val=DIST[q*Np+idx];
    						error = fabs(val - 0.1111111111111111*nA);
    						if (error > 1.0e-12) {
    							printf("   q=%i, Aq=%f \n",q,val);
    							errc_odd_b++;
    						}
    					}
					}
				}
			}
		}

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

		ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm.next, Np);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		timestep++;

		printf("Check after even time \n");

		SIZE=2*Np*sizeof(double);
		ScaLBL_CopyToHost(&DENSITY[0],&Den[0],SIZE);

    	// Check the distributions
		SIZE=7*Np*sizeof(double);
		ScaLBL_CopyToHost(&DIST[0],&Aq[0],SIZE);

		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (Dm.id[n] > 0){
						int idx = Map(i,j,k);
    					nA=DENSITY[idx];
    					nB=DENSITY[Np+idx];
    					//printf("i,j,k=%i,%i,%i \n",i,j,k);
    					//printf("   nA=%f, nB=%f \n",nA,nB);
    					double val=DIST[idx];
    					double error = fabs(val - 0.3333333333333333*nA);
    					if (error > 1.0e-12) {
    						printf("   q=0, Aq=%f \n",val);
    						errc_even_b++;
    					}
    					for (int q=1; q<7; q++){
    						val=DIST[q*Np+idx];
    						error = fabs(val - 0.1111111111111111*nA);
    						if (error > 1.0e-12) {
    							printf("   q=%i, Aq=%f \n",q,val);
    							errc_even_b++;
    						}
    					}
					}
				}
			}
		}

		//*******************Component B*******************
		// initialize phi based on PhaseLabel (include solid component labels)
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					// constant color
					PhaseLabel[n]= -1.0;
				}
			}
		}
		ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));

		if (rank==0)	printf ("Initializing distributions \n");
		// Initialize the phase field and variables
		ScaLBL_D3Q19_Init(fq, Np);
		if (rank==0)	printf ("Initializing phase field \n");
		ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, Np);

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

		ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm.next, Np);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		timestep++;


		printf("Check after odd time \n");
		SIZE=2*Np*sizeof(double);
		ScaLBL_CopyToHost(&DENSITY[0],&Den[0],SIZE);

    	// Check the distributions
		SIZE=7*Np*sizeof(double);
		ScaLBL_CopyToHost(&DIST[0],&Bq[0],SIZE);

		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (Dm.id[n] > 0){
						int idx = Map(i,j,k);
    					nA=DENSITY[idx];
    					nB=DENSITY[Np+idx];
    					//printf("i,j,k=%i,%i,%i \n",i,j,k);
    					//printf("   nA=%f, nB=%f \n",nA,nB);
    					double val=DIST[idx];
    					double error = fabs(val - 0.3333333333333333*nB);
    					if (error > 1.0e-12) {
    						printf("   q=0, Bq=%f \n",val);
    						errc_odd_b++;
    					}
    					for (int q=1; q<7; q++){
    						val=DIST[q*Np+idx];
    						error = fabs(val - 0.1111111111111111*nB);
    						if (error > 1.0e-12) {
    							printf("   q=%i, Bq=%f \n",q,val);
    							errc_odd_b++;
    						}
    					}
					}
				}
			}
		}


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

		ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm.next, Np);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		timestep++;
		printf("Check after even time \n");

		SIZE=2*Np*sizeof(double);
		ScaLBL_CopyToHost(&DENSITY[0],&Den[0],SIZE);

    	// Check the distributions
		SIZE=7*Np*sizeof(double);
		ScaLBL_CopyToHost(&DIST[0],&Bq[0],SIZE);

		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (Dm.id[n] > 0){
						int idx = Map(i,j,k);
    					nA=DENSITY[idx];
    					nB=DENSITY[Np+idx];
    					//printf("i,j,k=%i,%i,%i \n",i,j,k);
    					//printf("   nA=%f, nB=%f \n",nA,nB);
    					double val=DIST[idx];
    					double error = fabs(val - 0.3333333333333333*nB);
    					if (error > 1.0e-12) {
    						printf("   q=0, Bq=%f \n",val);
    						errc_even_b++;
    					}
    					for (int q=1; q<7; q++){
    						val=DIST[q*Np+idx];
    						error = fabs(val - 0.1111111111111111*nB);
    						if (error > 1.0e-12) {
    							printf("   q=%i, Bq=%f \n",q,val);
    							errc_even_b++;
    						}
    					}
					}
				}
			}
		}
		printf("Error counts: A even=%i, A odd=%i, B even=%i, B odd=%i \n",errc_even_a,errc_odd_a,errc_even_b,errc_odd_b);
		int errc_total=errc_even_a+errc_odd_a+errc_even_b+errc_odd_b;
		if (errc_total>0) check=1;
		else check=0;

	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************

	return check;
}

