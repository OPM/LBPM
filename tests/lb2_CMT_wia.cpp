#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "pmmc.h"
#include "Domain.h"
#include "Extras.h"
#include "D3Q19.h"
#include "D3Q7.h"
#include "Color.h"
//#include "common/MPI_Helpers.h"
//#include "Communication.h"

//#define CBUB

#define NC 2

using namespace std;

inline double NormProb(short int value, short int *mu, short int *sigma, int k)
{
	double sum,m,s;
	for (int i=0; i<NC; i++){
		m = double(mu[i]);
		s = double(sigma[i]);
		sum += exp(-(value-m)*(value-m)/(2.0*s*s));
	}
	m = double(mu[k]);
	s = double(sigma[k]);
	return exp(-(value-m)*(value-m)/(2.0*s*s)) / sum;
}

extern "C" void CMT_ScaLBL_D3Q7_ColorCollideMass(char *ID, double *A_even, double *A_odd, double *B_even, double *B_odd, 
		double *Den, double *Phi, double *ColorGrad, double beta, int N)
{
	char id;

	int idx,n,q,Cqx,Cqy,Cqz;
	//	int sendLoc;

	double f0,f1,f2,f3,f4,f5,f6;
	double na,nb;		// density values
	double ux,uy,uz;	// flow velocity
	double nx,ny,nz,C;	// color gradient components
	double a1,a2,b1,b2;
	double sp,delta;
	double feq[6];		// equilibrium distributions
	// Set of Discrete velocities for the D3Q19 Model
	int D3Q7[3][3]={{1,0,0},{0,1,0},{0,0,1}};

	for (n=0; n<N; n++){
		id = ID[n];
		if (id > 0 ){

			//.....Load the Color gradient.........
			nx = ColorGrad[n];
			ny = ColorGrad[N+n];
			nz = ColorGrad[2*N+n];
			C = sqrt(nx*nx+ny*ny+nz*nz);
			nx = nx/C;
			ny = ny/C;
			nz = nz/C;
			//........................................................................
			//					READ THE DISTRIBUTIONS
			//		(read from opposite array due to previous swap operation)
			//........................................................................
			f2 = A_odd[n];
			f4 = A_odd[N+n];
			f6 = A_odd[2*N+n];
			f0 = A_even[n];
			f1 = A_even[N+n];
			f3 = A_even[2*N+n];
			f5 = A_even[3*N+n];
			na = f0+f1+f2+f3+f4+f5+f6;
			//........................................................................
			f2 = B_odd[n];
			f4 = B_odd[N+n];
			f6 = B_odd[2*N+n];
			f0 = B_even[n];
			f1 = B_even[N+n];
			f3 = B_even[2*N+n];
			f5 = B_even[3*N+n];
			nb = f0+f1+f2+f3+f4+f5+f6;
			//........................................................................
			//....Instantiate the density distributions
			// Generate Equilibrium Distributions and stream
			// Stationary value - distribution 0
			A_even[n] = 0.3333333333333333*na;
			B_even[n] = 0.3333333333333333*nb;
			// Non-Stationary equilibrium distributions
			feq[0] = 0.1111111111111111*(1+3*ux);
			feq[1] = 0.1111111111111111*(1-3*ux);
			feq[2] = 0.1111111111111111*(1+3*uy);
			feq[3] = 0.1111111111111111*(1-3*uy);
			feq[4] = 0.1111111111111111*(1+3*uz);
			feq[5] = 0.1111111111111111*(1-3*uz);
			// Construction and streaming for the components
			for (idx=0; idx<3; idx++){
				//...............................................
				// Distribution index
				q = 2*idx;
				// Associated discrete velocity
				Cqx = D3Q7[idx][0];
				Cqy = D3Q7[idx][1];
				Cqz = D3Q7[idx][2];
				// Generate the Equilibrium Distribution
				a1 = na*feq[q];
				b1 = nb*feq[q];
				a2 = na*feq[q+1];
				b2 = nb*feq[q+1];
				// Recolor the distributions
				if (C > 0.0){
					sp = nx*double(Cqx)+ny*double(Cqy)+nz*double(Cqz);
					//if (idx > 2)	sp = 0.7071067811865475*sp;
					//delta = sp*min( min(a1,a2), min(b1,b2) );
					delta = na*nb/(na+nb)*0.1111111111111111*sp;
					//if (a1>0 && b1>0){
					a1 += beta*delta;
					a2 -= beta*delta;
					b1 -= beta*delta;
					b2 += beta*delta;
				}
				// Save the re-colored distributions
				A_odd[N*idx+n] 		= a1;
				A_even[N*(idx+1)+n] = a2;
				B_odd[N*idx+n] 		= b1;
				B_even[N*(idx+1)+n] = b2;
				//...............................................
			}
		}
	}
}


int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc,&argv);

	int n,N,Nx,Ny,Nz;

	Nx = Ny = Nz = 202;
	int rank = 12;
	char LocalRankString[8];
	char LocalRankFilename[40];
	char LocalRestartFile[40];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRestartFile,"%s%s","Solid.",LocalRankString);

	// Peaks of the standard normal distributions that approximate the data distribution
	double beta = 0.99;
	short int *mu;
	short int *sigma;
	mu = new short int [NC];
	sigma = new short int [NC];

 	mu[0] = 27200;	sigma[0] = 1500;
 	mu[1] = -29000;	sigma[1] = 1200;

	printf("Nx, Ny, Nz: %i, %i %i \n", Nx,Ny,Nz);
	printf("Filename = %s \n",LocalRestartFile);
	printf("Number of components = %i \n",NC);

	N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);

/*	//......................device distributions.................................
	double *f_even,*f_odd;
	double *A_even,*A_odd,*B_even,*B_odd;
	//...........................................................................
	ScaLBL_AllocateDeviceMemory((void **) &f_even, 10*dist_mem_size);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &f_odd, 9*dist_mem_size);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &A_even, 4*dist_mem_size);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &A_odd, 3*dist_mem_size);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &B_even, 4*dist_mem_size);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &B_odd, 3*dist_mem_size);	// Allocate device memory
*/	
	printf("Set up ID \n");
	char *ID;
	ScaLBL_AllocateDeviceMemory((void **) &ID, N);
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				n = k*Nx*Ny+j*Nx+i;
				ID[n] = 0;
			}
		}
	}
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				n = k*Nx*Ny+j*Nx+i;
				ID[n] = 1;	
			}
		}
	}

	printf("Allocate memory \n");

	double *Phi,*Den, *ColorGrad;
	ScaLBL_AllocateDeviceMemory((void **) &Phi, dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Den, NC*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &ColorGrad, dist_mem_size);
	
	double *packed_even,*packed_odd;
	ScaLBL_AllocateDeviceMemory((void **) &packed_even, 4*dist_mem_size*NC);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &packed_odd, 3*dist_mem_size*NC);	// Allocate device memory
	
	//..............................................
	// Read the input file
	//..............................................
	printf("Read files \n");
	short int value;
	short int *Data;
	Data = new short int [N];
	ifstream File(LocalRestartFile,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		File.read((char*) &value, sizeof(value));
		Data[n] = value;
	}
	File.close();
	//..............................................
	// Initialize the density from the input file
	//..............................................
	printf("Initialize density... \n");
	double m,s,val;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				n = k*Nx*Ny+j*Nx+i;
				short int img_val;
				img_val = Data[n];
				val = double(img_val);
				double sum = 0.0;
				for (int nc=0; nc<NC; nc++){
					m = double(mu[nc]);
					s = double(sigma[nc]);
					sum += exp(-(val-m)*(val-m)/(2.0*s*s));
				}
				for (int nc=0; nc<NC; nc++){
					m = double(mu[nc]);
					s = double(sigma[nc]);
					if (sum != 0.0)
						Den[N*nc+n] =  exp(-(val-m)*(val-m)/(2.0*s*s)) / sum;
					else Den[N*nc+n] = 0.0;

					//Den[N*nc+n] = 1.0*img_val;//NormProb(img_val, mu, sigma, nc);
				}
			}
		}
	}
	//..............................................
	printf("Initialize distribution \n");
	ScaLBL_D3Q7_Init(ID, &packed_even[0], &packed_odd[0], &Den[0], Nx, Ny, Nz);
	ScaLBL_D3Q7_Init(ID, &packed_even[4*N], &packed_odd[3*N], &Den[N], Nx, Ny, Nz);
	ScaLBL_ComputePhaseField(ID, Phi, Den, N);

	int timestep=0;
	int timestepMax=10;
	printf("# timesteps for the LBM = %i \n",timestepMax);

	while (timestep < timestepMax){
		
		ScaLBL_D3Q19_ColorGradient(ID,Phi,ColorGrad,Nx,Ny,Nz);

		CMT_ScaLBL_D3Q7_ColorCollideMass(ID, &packed_even[0], &packed_odd[0], &packed_even[4*N], &packed_odd[3*N], Den, Phi,
								ColorGrad, beta, N);
		
/*		//...................................................................................
		ScaLBL_D3Q19_Pack(1,dvcSendList_x,0,sendCount_x,sendbuf_x,A_even,N);
		ScaLBL_D3Q19_Pack(1,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,B_even,N);
		//...Packing for X face(1,7,9,11,13)................................
		ScaLBL_D3Q19_Pack(0,dvcSendList_X,0,sendCount_X,sendbuf_X,A_odd,N);
		ScaLBL_D3Q19_Pack(0,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,B_odd,N);
		//...Packing for y face(4,8,9,16,18).................................
		ScaLBL_D3Q19_Pack(2,dvcSendList_y,0,sendCount_y,sendbuf_y,A_even,N);
		ScaLBL_D3Q19_Pack(2,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,B_even,N);
		//...Packing for Y face(3,7,10,15,17).................................
		ScaLBL_D3Q19_Pack(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,A_odd,N);
		ScaLBL_D3Q19_Pack(1,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,B_odd,N);
		//...Packing for z face(6,12,13,16,17)................................
		ScaLBL_D3Q19_Pack(3,dvcSendList_z,0,sendCount_z,sendbuf_z,A_even,N);
		ScaLBL_D3Q19_Pack(3,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,B_even,N);
		//...Packing for Z face(5,11,14,15,18)................................
		ScaLBL_D3Q19_Pack(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,A_odd,N);
		ScaLBL_D3Q19_Pack(2,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,B_odd,N);
		//...................................................................................

		//...................................................................................
		// Send all the distributions
		MPI_Isend(sendbuf_x, 2*sendCount_x,MPI_DOUBLE,rank_x,sendtag,comm,&req1[0]);
		MPI_Irecv(recvbuf_X, 2*recvCount_X,MPI_DOUBLE,rank_X,recvtag,comm,&req2[0]);
		MPI_Isend(sendbuf_X, 2*sendCount_X,MPI_DOUBLE,rank_X,sendtag,comm,&req1[1]);
		MPI_Irecv(recvbuf_x, 2*recvCount_x,MPI_DOUBLE,rank_x,recvtag,comm,&req2[1]);
		MPI_Isend(sendbuf_y, 2*sendCount_y,MPI_DOUBLE,rank_y,sendtag,comm,&req1[2]);
		MPI_Irecv(recvbuf_Y, 2*recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,comm,&req2[2]);
		MPI_Isend(sendbuf_Y, 2*sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,comm,&req1[3]);
		MPI_Irecv(recvbuf_y, 2*recvCount_y,MPI_DOUBLE,rank_y,recvtag,comm,&req2[3]);
		MPI_Isend(sendbuf_z, 2*sendCount_z,MPI_DOUBLE,rank_z,sendtag,comm,&req1[4]);
		MPI_Irecv(recvbuf_Z, 2*recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,comm,&req2[4]);
		MPI_Isend(sendbuf_Z, 2*sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,comm,&req1[5]);
		MPI_Irecv(recvbuf_z, 2*recvCount_z,MPI_DOUBLE,rank_z,recvtag,comm,&req2[5]);
*/		//...................................................................................
		
		ScaLBL_D3Q7_Swap(ID, &packed_even[0], &packed_odd[0], Nx, Ny, Nz);
		ScaLBL_D3Q7_Swap(ID, &packed_even[4*N], &packed_odd[3*N], Nx, Ny, Nz);
		
/*		//...................................................................................
		// Wait for completion of D3Q19 communication
		MPI_Waitall(6,req1,stat1);
		MPI_Waitall(6,req2,stat2);
		//...................................................................................
		// Unpack the distributions on the device
		//...................................................................................
		//...Map recieve list for the X face: q=2,8,10,12,13 .................................
		ScaLBL_D3Q19_Unpack(0,-1,0,0,dvcRecvList_X,0,recvCount_X,recvbuf_X,A_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(0,-1,0,0,dvcRecvList_X,recvCount_X,recvCount_X,recvbuf_X,B_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the x face: q=1,7,9,11,13..................................
		ScaLBL_D3Q19_Unpack(1,1,0,0,dvcRecvList_x,0,recvCount_x,recvbuf_x,A_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(1,1,0,0,dvcRecvList_x,recvCount_x,recvCount_x,recvbuf_x,B_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
		ScaLBL_D3Q19_Unpack(1,0,-1,0,dvcRecvList_Y,0,recvCount_Y,recvbuf_Y,A_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(1,0,-1,0,dvcRecvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,B_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
		ScaLBL_D3Q19_Unpack(2,0,1,0,dvcRecvList_y,0,recvCount_y,recvbuf_y,A_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(2,0,1,0,dvcRecvList_y,recvCount_y,recvCount_y,recvbuf_y,B_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
		ScaLBL_D3Q19_Unpack(2,0,0,-1,dvcRecvList_Z,0,recvCount_Z,recvbuf_Z,A_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(2,0,0,-1,dvcRecvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,B_odd,Nx,Ny,Nz);
		//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
		ScaLBL_D3Q19_Unpack(3,0,0,1,dvcRecvList_z,0,recvCount_z,recvbuf_z,A_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(3,0,0,1,dvcRecvList_z,recvCount_z,recvCount_z,recvbuf_z,B_even,Nx,Ny,Nz);
		//..................................................................................
*/		
		//..................................................................................
		ScaLBL_D3Q7_Density(ID, &packed_even[0], &packed_odd[0], &Den[0], Nx, Ny, Nz);
		ScaLBL_D3Q7_Density(ID, &packed_even[4*N], &packed_odd[3*N], &Den[N], Nx, Ny, Nz);
		
		//*************************************************************************
		// 		Compute the phase indicator field 
		//*************************************************************************
		ScaLBL_ComputePhaseField(ID, Phi, Den, N);
		//*************************************************************************
		timestep++;
	}
	printf("Write density values \n");
	FILE *PHASE;
	PHASE = fopen("Density.out","wb");
	fwrite(Den,8,2*N,PHASE);
	fclose(PHASE);
    
    // Close MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
