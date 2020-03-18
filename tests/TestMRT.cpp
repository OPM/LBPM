
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


extern void AA1_ScaLBL_D3Q19_Init(double *dist_even, double *dist_odd, int Nx, int Ny, int Nz,
		int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz)
{
	// Set of Discrete velocities for the D3Q19 Model
	static int D3Q19[18][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
			{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},
			{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},
			{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};

	int q,i,j,k,n,N;
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

				n = (k+1)*(Nx+2)*(Ny+2) + (j+1)*(Nx+2) + i+1;

				// Get the 'global' index
				x = iproc*Nx+i;
				y = jproc*Ny+j;
				z = kproc*Nz+k;
				for (q=0; q<9; q++){
					// Even distribution
					Cqx = D3Q19[2*q][0];
					Cqy = D3Q19[2*q][1];
					Cqz = D3Q19[2*q][2];
					//					xn = x - Cqx;
					//					yn = y - Cqy;
					//					zn = z - Cqz;
					xn = x;
					yn = y;
					zn = z;
					if (xn < 0) xn += nprocx*Nx;
					if (yn < 0) yn += nprocy*Ny;
					if (zn < 0) zn += nprocz*Nz;
					if (!(xn < nprocx*Nx)) xn -= nprocx*Nx;
					if (!(yn < nprocy*Ny)) yn -= nprocy*Ny;
					if (!(zn < nprocz*Nz)) zn -= nprocz*Nz;

					dist_even[(q+1)*N+n] = (zn*X*Y+yn*X+xn) + (2*q+1)*0.01;

					// Odd distribution
					//					xn = x + Cqx;
					//					yn = y + Cqy;
					//					zn = z + Cqz;
					xn = x;
					yn = y;
					zn = z;
					if (xn < 0) xn += nprocx*Nx;
					if (yn < 0) yn += nprocy*Ny;
					if (zn < 0) zn += nprocz*Nz;
					if (!(xn < nprocx*Nx)) xn -= nprocx*Nx;
					if (!(yn < nprocy*Ny)) yn -= nprocy*Ny;
					if (!(zn < nprocz*Nz)) zn -= nprocz*Nz;

					dist_odd[q*N+n] =  (zn*X*Y+yn*X+xn) + 2*(q+1)*0.01;

				}
			}
		}
	}

}

extern void GlobalFlipScaLBL_D3Q19_Init(double *dist_even, double *dist_odd, int Nx, int Ny, int Nz, 
		int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz)
{
	// Set of Discrete velocities for the D3Q19 Model
	static int D3Q19[18][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
			{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},
			{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},
			{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};

	int q,i,j,k,n,N;
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

				n = (k+1)*(Nx+2)*(Ny+2) + (j+1)*(Nx+2) + i+1;

				// Get the 'global' index
				x = iproc*Nx+i;
				y = jproc*Ny+j;
				z = kproc*Nz+k;
				for (q=0; q<9; q++){
					// Even distribution
					Cqx = D3Q19[2*q][0];
					Cqy = D3Q19[2*q][1];
					Cqz = D3Q19[2*q][2];
					//					xn = x - Cqx;
					//					yn = y - Cqy;
					//					zn = z - Cqz;
					xn = x;
					yn = y;
					zn = z;
					if (xn < 0) xn += nprocx*Nx;
					if (yn < 0) yn += nprocy*Ny;
					if (zn < 0) zn += nprocz*Nz;
					if (!(xn < nprocx*Nx)) xn -= nprocx*Nx;
					if (!(yn < nprocy*Ny)) yn -= nprocy*Ny;
					if (!(zn < nprocz*Nz)) zn -= nprocz*Nz;	

					dist_even[(q+1)*N+n] = (zn*X*Y+yn*X+xn) + (2*q+1)*0.01;

					// Odd distribution
					//					xn = x + Cqx;
					//					yn = y + Cqy;
					//					zn = z + Cqz;
					xn = x;
					yn = y;
					zn = z;
					if (xn < 0) xn += nprocx*Nx;
					if (yn < 0) yn += nprocy*Ny;
					if (zn < 0) zn += nprocz*Nz;
					if (!(xn < nprocx*Nx)) xn -= nprocx*Nx;
					if (!(yn < nprocy*Ny)) yn -= nprocy*Ny;
					if (!(zn < nprocz*Nz)) zn -= nprocz*Nz;

					dist_odd[q*N+n] =  (zn*X*Y+yn*X+xn) + 2*(q+1)*0.01;

				}
			}
		}
	}

}


extern int GlobalCheckDebugDistInterior(double *dist_even, double *dist_odd, int Nx, int Ny, int Nz,
		int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz)
{

	int returnValue = 0;
	int q,i,j,k,n,N;
	int Cqx,Cqy,Cqz; // Discrete velocity
	int x,y,z;		// Global indices
	int xn,yn,zn; 	// Global indices of neighbor
	int X,Y,Z;		// Global size
	X = Nx*nprocx;
	Y = Ny*nprocy;
	Z = Nz*nprocz;
	NULL_USE(Z);
	N = (Nx+2)*(Ny+2)*(Nz+2);	// size of the array including halo
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){

				n = (k+1)*(Nx+2)*(Ny+2) + (j+1)*(Nx+2) + i+1;

				// Get the 'global' index
				x = iproc*Nx+i;
				y = jproc*Ny+j;
				z = kproc*Nz+k;
				for (q=0; q<9; q++){

					if (dist_even[(q+1)*N+n] != (z*X*Y+y*X+x) + 2*(q+1)*0.01){
						printf("******************************************\n");
						printf("error in even distribution q = %i \n", 2*(q+1));
						printf("i,j,k= %i, %i, %i \n", x,y,z);
						printf("dist = %5.2f, expect %5.2f \n", dist_even[(q+1)*N+n], (z*X*Y+y*X+x) + 2*(q+1)*0.01);
						printf("n= %i \n",z*X*Y+y*X+x);
						returnValue++;
					}


					if (dist_odd[q*N+n] !=  (z*X*Y+y*X+x) + (2*q+1)*0.01){
						printf("******************************************\n");
						printf("error in odd distribution q = %i \n", 2*q+1);
						printf("i,j,k= %i, %i, %i \n", x,y,z);
						printf("dist = %5.2f, expect %5.2f \n", dist_odd[q*N+n],(z*X*Y+y*X+x) + (2*q+1)*0.01);
						printf("n= %i \n",z*X*Y+y*X+x);
						returnValue++;
					}
				}
			}
		}
	}
	return returnValue;
}


extern void HostToGold(IntArray Map, double * f_even_host, double * f_odd_host, double * f_even_gold,
		double * f_odd_gold, int Nx, int Ny, int Nz,int Np,int N) {
	int n;
	for (int k=1;k<Nz-1;k++){
		for (int j=1;j<Ny-1;j++){
			for (int i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				int idx = Map(i,j,k);
				if (!(idx < 0)){
					for (int q=0; q<9; q++){
						f_even_gold[(q+1)*N + n] = f_even_host[(q+1)*Np + idx];
						f_odd_gold[q*N + n] = f_odd_host[q*Np + idx];
					}
				}
			}
		}
	}
}

extern void GoldToHost(IntArray Map, double * f_even_host, double * f_odd_host, double * f_even_gold,
		double * f_odd_gold, int Nx, int Ny, int Nz,int Np,int N) {
	int n;
	for (int k=1;k<Nz-1;k++){
		for (int j=1;j<Ny-1;j++){
			for (int i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				int idx = Map(i,j,k);
				if (!(idx < 0)){
					for (int q=0; q<9; q++){
						f_even_host[(q+1)*Np + idx] = f_even_gold[(q+1)*N + n];
						f_odd_host[q*Np + idx] = f_odd_gold[q*N + n];
					}
				}
			}
		}
	}
}



extern void PrintSpecificHost(int centralNode, IntArray Map, double * f_even_host, double * f_odd_host, double * f_even_gold,
		double * f_odd_gold, int Nx, int Ny, int Nz,int Np,int N) {
	int n;
	for (int k=1;k<Nz-1;k++){
		for (int j=1;j<Ny-1;j++){
			for (int i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				int idx = Map(i,j,k);
				if (!(idx < 0)){
					if ( i == centralNode && (j == centralNode && k == centralNode)) {
						printf("%i: ",idx);
					}

					for (int q=0; q<2; q++){
						if ( i == centralNode && (j == centralNode && k == centralNode)) {
							printf("%.02f,",f_odd_host[(q*Np+idx)]);
							printf("%.02f,",f_even_host[(q+1)*Np + idx]);
						}
					}
				}
			}
		}
	}
	printf("\n\n");
}



extern void PrintFullHost(IntArray Map, double * f_even_host, double * f_odd_host, double * f_even_gold,
		double * f_odd_gold, int Nx, int Ny, int Nz,int Np,int N) {
	int n;
	for (int k=1;k<Nz-1;k++){
		for (int j=1;j<Ny-1;j++){
			for (int i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				int idx = Map(i,j,k);
				if (!(idx < 0)){
					//printf("%i: ",idx);
					for (int q=0; q<1; q++){

						printf("  %.02f  <%d>  ",f_even_host[(q+1)*Np + idx],idx);
						printf("%.02f  |",f_odd_host[(q*Np+idx)]);
					}
				}

			}
		}
	}
	printf("\n");
}

extern void HostToUnobtainium(IntArray Map, double * f_even_host, double * f_odd_host, double * f_even_unobtainium,
		double * f_odd_unobtainium, int Nx, int Ny, int Nz,int Np,int N) {
	int n;
	for (int k=1;k<Nz-1;k++){
		for (int j=1;j<Ny-1;j++){
			for (int i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				int idx = Map(i,j,k);
				if (!(idx < 0)){
					for (int q=0; q<9; q++){
						f_odd_unobtainium[(q*Np+idx)] = f_odd_host[(q*Np+idx)];
						f_even_unobtainium[(q+1)*Np + idx] = f_even_host[(q+1)*Np + idx];
					}
				}
			}
		}
	}
}

extern void CheckDistrMatch(IntArray Map, double * f_even_host, double * f_odd_host, double * f_even_unobtainium,
		double * f_odd_unobtainium, int Nx, int Ny, int Nz,int Np,int N) {
	int n; int err = 0;
	for (int k=1;k<Nz-1;k++){
		for (int j=1;j<Ny-1;j++){
			for (int i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				int idx = Map(i,j,k);
				if (!(idx < 0)){
					for (int q=0; q<9; q++){
						if (f_odd_unobtainium[(q*Np+idx)] != f_odd_host[(q*Np+idx)]) {
							err++;
						}
						if (f_even_unobtainium[(q+1)*Np + idx] != f_even_host[(q+1)*Np + idx]) {
							err++;
						}
					}
				}
			}
		}
	}
	if (err == 0) {
		printf("CORRECT");
	} else {
		printf("DIFF=%d",err);
	}
	printf("\n\n");
}

extern double CheckDistrMatchDouble(IntArray Map, double * f_even_host, double * f_odd_host, double * f_even_unobtainium,
		double * f_odd_unobtainium, int Nx, int Ny, int Nz,int Np,int N) {
	int n; int err = 0;
	for (int k=1;k<Nz-1;k++){
		for (int j=1;j<Ny-1;j++){
			for (int i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				int idx = Map(i,j,k);
				if (!(idx < 0)){
					for (int q=0; q<9; q++){
						if (f_odd_unobtainium[(q*Np+idx)] != f_odd_host[(q*Np+idx)]) {
							err++;
						}
						if (f_even_unobtainium[(q+1)*Np + idx] != f_even_host[(q+1)*Np + idx]) {
							err++;
						}
					}
				}
			}
		}
	}
	return err;
}


extern void PrintNeighborList(int * neighborList, int Np, int rank) {
	if (rank == 0) {
		int n;
		int neighbor;

		for (int i = 0; i < Np; i++) {
			printf("idx=%d: ",i);
			for (int l = 0; l < 10; l++) {   // was 18
				neighbor = neighborList[l*Np + i];
				printf("%d ",neighbor);
			}
			printf("\n");
		}
		printf("\n\n");
	}
}

extern void PrintSpecificNeighborList(int specificNode, int * neighborList, int Np, int rank) {

	if (rank == 0) {
		int n;
		int neighbor;

		for (int i = 0; i < Np; i++) {
			if ( i == specificNode) {
				printf("idx=%d: ",i);
			}
			for (int l = 0; l < 10; l++) {  // was 18
				neighbor = neighborList[l*Np + i];
				if (i == specificNode) {
					printf("%d ",neighbor);
				}
			}
		}
		printf("\n\n");
	}
}


extern int GlobalCheckDebugDist(double *dist_even, double *dist_odd, int Nx, int Ny, int Nz, 
		int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz) {

	int returnValue = 0;
	int q,i,j,k,n,N;
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

				n = (k+1)*(Nx+2)*(Ny+2) + (j+1)*(Nx+2) + i+1;

				// Get the 'global' index
				x = iproc*Nx+i;
				y = jproc*Ny+j;
				z = kproc*Nz+k;
				for (q=0; q<9; q++){

					if (dist_even[(q+1)*N+n] != (z*X*Y+y*X+x) + 2*(q+1)*0.01){
						printf("******************************************\n");
						printf("error in even distribution q = %i \n", 2*(q+1));
						printf("i,j,k= %i, %i, %i \n", x,y,z);
						printf("dist = %5.2f, expect %5.2f \n", dist_even[(q+1)*N+n], (z*X*Y+y*X+x) + 2*(q+1)*0.01);
						printf("n= %i \n",z*X*Y+y*X+x);
						returnValue++;
					}


					if (dist_odd[q*N+n] !=  (z*X*Y+y*X+x) + (2*q+1)*0.01){
						printf("******************************************\n");
						printf("error in odd distribution q = %i \n", 2*q+1);
						printf("i,j,k= %i, %i, %i \n", x,y,z);
						printf("dist = %5.2f, expect %5.2f \n", dist_odd[q*N+n],(z*X*Y+y*X+x) + (2*q+1)*0.01);
						printf("n= %i \n",z*X*Y+y*X+x);
						returnValue++;
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
		// parallel domain size (# of sub-domains)
		int nprocx,nprocy,nprocz;
		int iproc,jproc,kproc;


		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Permeability Test: TestMRT	\n");
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
		int timesteps = 1000;
		int centralNode = 2;

		double tau = 1.0;
		double rlx_setA = 1.f/tau;
		double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);

		Fx = Fy = 0.f;
		Fz = 1.0e-6;

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
			printf("Process grid = %i x %i x %i\n",nprocx,nprocy,nprocz);
			printf("********************************************************\n");
		}

		MPI_Barrier(comm);
		kproc = rank/(nprocx*nprocy);
		jproc = (rank-nprocx*nprocy*kproc)/nprocx;
		iproc = rank-nprocx*nprocy*kproc-nprocz*jproc;

		if (rank == 0) {
			printf("i,j,k proc=%d %d %d \n",iproc,jproc,kproc);
		}
		MPI_Barrier(comm);
		if (rank == 1){
			printf("i,j,k proc=%d %d %d \n",iproc,jproc,kproc);
			printf("\n\n");
		}

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

		FILE *IDFILE = fopen(LocalRankFilename,"rb");
		if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
		fread(Dm.id,1,N,IDFILE);
		fclose(IDFILE);

		MPI_Barrier(comm);
		Dm.CommInit();

		//.......................................................................
		// Compute the media porosity
		//.......................................................................
		double sum;
		double sum_local=0.0, porosity;
		int Np=0;  // number of local pore nodes
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
		MPI_Barrier(comm);
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
		//		double *f_even,*f_odd;
		double * dist;
		double * Velocity;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &dist, 19*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
		ScaLBL_CopyToDevice(NeighborList,     neighborList, neighborSize);
		//...........................................................................

		ScaLBL_D3Q19_Init(dist, Np);

		/************ MAIN ITERATION LOOP (timing communications)***************************************/

		if (rank==0) printf("Beginning AA timesteps...\n");
		if (rank==0) printf("********************************************************\n");
		if (rank==0) printf("No. of timesteps for timing: %i \n", timesteps);

		//.......create and start timer............
		double starttime,stoptime,cputime;

		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		starttime = MPI_Wtime();

		while (timestep < timesteps) {
			
			ScaLBL_Comm.SendD3Q19AA(dist); //READ FROM NORMAL
			ScaLBL_D3Q19_AAodd_MRT(NeighborList, dist, ScaLBL_Comm.next, Np, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
			ScaLBL_Comm.RecvD3Q19AA(dist); //WRITE INTO OPPOSITE
			ScaLBL_D3Q19_AAodd_MRT(NeighborList, dist, 0, ScaLBL_Comm.next, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
			timestep++;

			ScaLBL_Comm.SendD3Q19AA(dist); //READ FORM NORMAL
			ScaLBL_D3Q19_AAeven_MRT(dist, ScaLBL_Comm.next, Np, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
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
		if (rank==0) printf("********************************************************\n");
		if (rank==0) printf("CPU time = %f \n", cputime);
		if (rank==0) printf("Lattice update rate (per process)= %f MLUPS \n", MLUPS);
		MLUPS *= nprocs;
		if (rank==0) printf("Lattice update rate (process)= %f MLUPS \n", MLUPS);
		if (rank==0) printf("********************************************************\n");

		// Number of memory references from the swap algorithm (per timestep)
		// 18 reads and 18 writes for each lattice site
		double MemoryRefs = Np*38;
        // number of memory references for the swap algorithm - GigaBytes / second
        if (rank==0) printf("DRAM bandwidth (per process)= %f GB/sec \n",MemoryRefs*8*timestep/1e9/cputime);
        // Report bandwidth in Gigabits per second
        // communication bandwidth includes both send and recieve
        if (rank==0) printf("Communication bandwidth (per process)= %f Gbit/sec \n",ScaLBL_Comm.CommunicationCount*64*timestep/1e9/cputime);
        if (rank==0) printf("Aggregated communication bandwidth = %f Gbit/sec \n",nprocs*ScaLBL_Comm.CommunicationCount*64*timestep/1e9/cputime);


    	double *VEL;
    	VEL= new double [3*Np];
    	int SIZE=3*Np*sizeof(double);
    	ScaLBL_D3Q19_Momentum(dist,Velocity, Np);
    	ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    	ScaLBL_CopyToHost(&VEL[0],&Velocity[0],SIZE);

    	sum_local=0.f;
    	sum = 0.f;
    	for (k=1;k<Nz-1;k++){
    		for (j=1;j<Ny-1;j++){
    			for (i=1;i<Nx-1;i++){
    				n = k*Nx*Ny+j*Nx+i;
    				if (Dm.id[n] > 0){\
    					int idx = Map(i,j,k);
    					sum_local+=VEL[2*Np+idx];
    				}
    			}
    		}
    	}
    	MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
    	double PoreVel = sum*iVol_global;
    	if (rank==0) printf("Velocity = %f \n",PoreVel);

	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************

	return check;
}

