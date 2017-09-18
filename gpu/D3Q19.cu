#include <stdio.h>

#define NBLOCKS 1024
#define NTHREADS 128


// functionality for parallel reduction in Flux BC routines -- probably should be re-factored to another location
// functions copied from https://devblogs.nvidia.com/parallelforall/faster-parallel-reductions-kepler/

__inline__ __device__
double warpReduceSum(double val) {
	for (int offset = warpSize/2; offset > 0; offset /= 2)
		val += __shfl_down(val, offset);
	return val;
}

__inline__ __device__
double blockReduceSum(double val) {

	static __shared__ double shared[32]; // Shared mem for 32 partial sums
	int lane = threadIdx.x % warpSize;
	int wid = threadIdx.x / warpSize;

	val = warpReduceSum(val);     // Each warp performs partial reduction

	if (lane==0) shared[wid]=val; // Write reduced value to shared memory

	__syncthreads();              // Wait for all partial reductions

	//read from shared memory only if that warp existed
	val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

	if (wid==0) val = warpReduceSum(val); //Final reduce within first warp

	return val;
}

__global__ void deviceReduceKernel(double *in, double* out, int N) {
	double sum = 0;
	//reduce multiple elements per thread
	for (int i = blockIdx.x * blockDim.x + threadIdx.x;
			i < N;
			i += blockDim.x * gridDim.x) {
		sum += in[i];
	}
	sum = blockReduceSum(sum);
	if (threadIdx.x==0)
		out[blockIdx.x]=sum;
}

__global__  void dvc_ScaLBL_D3Q19_Pack(int q, int *list, int start, int count, double *sendbuf, double *dist, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		n = list[idx];
		sendbuf[start+idx] = dist[q*N+n];
	}
}

__global__ void dvc_ScaLBL_D3Q19_Unpack(int q,  int *list,  int start, int count,
		double *recvbuf, double *dist, int N){
	//....................................................................................
	// Unpack distribution from the recv buffer
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int n,idx;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = list[start+idx];
		// unpack the distribution to the proper location
		if (!(n<0)) dist[q*N+n] = recvbuf[start+idx];
	}
}
/*
__global__ void dvc_ScaLBL_D3Q19_MapRecv(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
						   int *d3q19_recvlist, int Nx, int Ny, int Nz){
	//....................................................................................
	// Unack distribution from the recv buffer
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int i,j,k,n,nn,idx;
	int N = Nx*Ny*Nz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = list[idx];
		// Get the 3-D indices
		k = n/(Nx*Ny);
		j = (n-Nx*Ny*k)/Nx;
		i = n-Nx*Ny*k-Nx*j;
		// Streaming for the non-local distribution
		i += Cqx;
		j += Cqy;
		k += Cqz;
		nn = k*Nx*Ny+j*Nx+i;
		// unpack the distribution to the proper location
		d3q19_recvlist[start+idx]=nn;
	}
}
 */


__global__ void dvc_ScaLBL_D3Q19_Init(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz)
{
	int n,N;
	N = Nx*Ny*Nz;
	char id;
	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N ){
			id = ID[n];
			if (id > 0 ){
				f_even[n] = 0.3333333333333333;
				f_odd[n] = 0.055555555555555555;		//double(100*n)+1.f;
				f_even[N+n] = 0.055555555555555555;	//double(100*n)+2.f;
				f_odd[N+n] = 0.055555555555555555;	//double(100*n)+3.f;
				f_even[2*N+n] = 0.055555555555555555;	//double(100*n)+4.f;
				f_odd[2*N+n] = 0.055555555555555555;	//double(100*n)+5.f;
				f_even[3*N+n] = 0.055555555555555555;	//double(100*n)+6.f;
				f_odd[3*N+n] = 0.0277777777777778;   //double(100*n)+7.f;
				f_even[4*N+n] = 0.0277777777777778;   //double(100*n)+8.f;
				f_odd[4*N+n] = 0.0277777777777778;   //double(100*n)+9.f;
				f_even[5*N+n] = 0.0277777777777778;  //double(100*n)+10.f;
				f_odd[5*N+n] = 0.0277777777777778;  //double(100*n)+11.f;
				f_even[6*N+n] = 0.0277777777777778;  //double(100*n)+12.f;
				f_odd[6*N+n] = 0.0277777777777778;  //double(100*n)+13.f;
				f_even[7*N+n] = 0.0277777777777778;  //double(100*n)+14.f;
				f_odd[7*N+n] = 0.0277777777777778;  //double(100*n)+15.f;
				f_even[8*N+n] = 0.0277777777777778;  //double(100*n)+16.f;
				f_odd[8*N+n] = 0.0277777777777778;  //double(100*n)+17.f;
				f_even[9*N+n] = 0.0277777777777778;  //double(100*n)+18.f;
			}
			else{
				for(int q=0; q<9; q++){
					f_even[q*N+n] = -1.0;
					f_odd[q*N+n] = -1.0;
				}
				f_even[9*N+n] = -1.0;
			}
		}
	}
}

//*************************************************************************
__global__  void dvc_ScaLBL_D3Q19_Swap_Compact(int *neighborList, double *disteven, double *distodd, int Np, int q){
	int n,nn;
	double f1,f2;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<Np){
			nn = neighborList[q*Np+n];
			if (!(nn<0)){
				f1 = distodd[q*Np+n];
				f2 = disteven[(q+1)*Np+nn];
				disteven[(q+1)*Np+nn] = f1;
				distodd[q*Np+n] = f2;
			}
		}
	}
}

__global__ void dvc_ScaLBL_AAodd_Compact(char * ID, int *d_neighborList, double *dist, int Np) {

	int n, nn, nnn;
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	char id;
	int nread,nwrite,naccess;
	int S = Np/NBLOCKS/NTHREADS+1;

	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		id = ID[n];
		if (n<Np) {
			if (id > 0) {

				f0 = dist[n];

				nread = d_neighborList[n];        naccess = Np;
				if (nread<0)  { nread=n;          naccess = 10*Np;       }
				f2 = dist[nread + naccess];

				nread = d_neighborList[n+2*Np];   naccess = 2*Np;
				if (nread<0)  { nread=n;          naccess = 11*Np;     }
				f4 = dist[nread + naccess];

				nread = d_neighborList[n+4*Np];   naccess = 3*Np;
				if (nread<0)  { nread=n;          naccess = 12*Np;     }
				f6 = dist[nread + naccess];

				nread = d_neighborList[n+6*Np];   naccess = 4*Np;
				if (nread<0)  { nread=n;          naccess = 13*Np;     }
				f8 = dist[nread + naccess];

				nread = d_neighborList[n+8*Np];   naccess = 5*Np;
				if (nread<0)  { nread=n;          naccess = 14*Np;     }
				f10 = dist[nread + naccess];

				nread = d_neighborList[n+10*Np];  naccess = 6*Np;
				if (nread<0)  { nread=n;          naccess = 15*Np;     }
				f12 = dist[nread + naccess];

				nread = d_neighborList[n+12*Np];  naccess = 7*Np;
				if (nread<0)  { nread=n;          naccess = 16*Np;     }
				f14 = dist[nread + naccess];

				nread = d_neighborList[n+14*Np];  naccess = 8*Np;
				if (nread<0)  { nread=n;          naccess = 17*Np;     }
				f16 = dist[nread + naccess];

				nread = d_neighborList[n+16*Np];  naccess = 9*Np;
				if (nread<0)  { nread=n;          naccess = 18*Np;     }
				f18 = dist[nread + naccess];



				nread = d_neighborList[n+Np];     naccess = 10*Np;
				if (nread<0)  { nread=n;          naccess = Np;    }
				f1 = dist[nread + naccess];

				nread = d_neighborList[n+3*Np];   naccess = 11*Np;
				if (nread<0)  { nread=n;          naccess = 2*Np;    }
				f3 = dist[nread + naccess];

				nread = d_neighborList[n+5*Np];   naccess = 12*Np;
				if (nread<0)  { nread=n;          naccess = 3*Np;    }
				f5 = dist[nread + naccess];

				nread = d_neighborList[n+7*Np];   naccess = 13*Np;
				if (nread<0)  { nread=n;          naccess = 4*Np;    }
				f7 = dist[nread + naccess];

				nread = d_neighborList[n+9*Np];   naccess = 14*Np;
				if (nread<0)  { nread=n;          naccess = 5*Np;    }
				f9 = dist[nread + naccess];

				nread = d_neighborList[n+11*Np];   naccess = 15*Np;
				if (nread<0)  { nread=n;          naccess = 6*Np;    }
				f11 = dist[nread + naccess];

				nread = d_neighborList[n+13*Np];   naccess = 16*Np;
				if (nread<0)  { nread=n;          naccess = 7*Np;    }
				f13 = dist[nread + naccess];

				nread = d_neighborList[n+15*Np];   naccess = 17*Np;
				if (nread<0)  { nread=n;          naccess = 8*Np;    }
				f15 = dist[nread + naccess];

				nread = d_neighborList[n+17*Np];   naccess = 18*Np;
				if (nread<0)  { nread=n;          naccess = 9*Np;    }
				f17 = dist[nread + naccess];




				// ORIGINAL CORRECT WRITES
				//				nwrite = d_neighborList[n];      naccess = 10*Np;
				//				if (nwrite<0) { nwrite=n;        naccess = Np;  }
				//				dist[nwrite + naccess]   = f1;

				//				nwrite = d_neighborList[n+Np];   naccess = Np;
				//				if (nwrite<0) { nwrite=n;        naccess = 10*Np; }
				//				dist[nwrite + naccess]   = f2;

				nread = d_neighborList[n];        naccess = Np;
				if (nread<0)  { nread=n;          naccess = 10*Np;       }
				dist[nread + naccess] = f1;

				nread = d_neighborList[n+2*Np];   naccess = 2*Np;
				if (nread<0)  { nread=n;          naccess = 11*Np;     }
				dist[nread + naccess] = f3;

				nread = d_neighborList[n+4*Np];   naccess = 3*Np;
				if (nread<0)  { nread=n;          naccess = 12*Np;     }
				dist[nread + naccess] = f5;

				nread = d_neighborList[n+6*Np];   naccess = 4*Np;
				if (nread<0)  { nread=n;          naccess = 13*Np;     }
				dist[nread + naccess] = f7;

				nread = d_neighborList[n+8*Np];   naccess = 5*Np;
				if (nread<0)  { nread=n;          naccess = 14*Np;     }
				dist[nread + naccess] = f9;

				nread = d_neighborList[n+10*Np];  naccess = 6*Np;
				if (nread<0)  { nread=n;          naccess = 15*Np;     }
				dist[nread + naccess] = f11;

				nread = d_neighborList[n+12*Np];  naccess = 7*Np;
				if (nread<0)  { nread=n;          naccess = 16*Np;     }
				dist[nread + naccess] = f13;

				nread = d_neighborList[n+14*Np];  naccess = 8*Np;
				if (nread<0)  { nread=n;          naccess = 17*Np;     }
				dist[nread + naccess] = f15;

				nread = d_neighborList[n+16*Np];  naccess = 9*Np;
				if (nread<0)  { nread=n;          naccess = 18*Np;     }
				dist[nread + naccess] = f17;



				nread = d_neighborList[n+Np];     naccess = 10*Np;
				if (nread<0)  { nread=n;          naccess = Np;    }
				dist[nread + naccess] = f2;

				nread = d_neighborList[n+3*Np];   naccess = 11*Np;
				if (nread<0)  { nread=n;          naccess = 2*Np;    }
				dist[nread + naccess] = f4;

				nread = d_neighborList[n+5*Np];   naccess = 12*Np;
				if (nread<0)  { nread=n;          naccess = 3*Np;    }
				dist[nread + naccess] = f6;

				nread = d_neighborList[n+7*Np];   naccess = 13*Np;
				if (nread<0)  { nread=n;          naccess = 4*Np;    }
				dist[nread + naccess] = f8;

				nread = d_neighborList[n+9*Np];   naccess = 14*Np;
				if (nread<0)  { nread=n;          naccess = 5*Np;    }
				dist[nread + naccess] = f10;

				nread = d_neighborList[n+11*Np];   naccess = 15*Np;
				if (nread<0)  { nread=n;          naccess = 6*Np;    }
				dist[nread + naccess]= f12;

				nread = d_neighborList[n+13*Np];   naccess = 16*Np;
				if (nread<0)  { nread=n;          naccess = 7*Np;    }
				dist[nread + naccess] = f14;

				nread = d_neighborList[n+15*Np];   naccess = 17*Np;
				if (nread<0)  { nread=n;          naccess = 8*Np;    }
				dist[nread + naccess] = f16;

				nread = d_neighborList[n+17*Np];   naccess = 18*Np;
				if (nread<0)  { nread=n;          naccess = 9*Np;    }
				dist[nread + naccess] = f18;
			}
		}
	}
}


__global__ void dvc_ScaLBL_AAeven_Compact(char * ID, double *dist, int Np) {


	int q,n,nn;
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	char id;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		id = ID[n];

		if ( n<Np ){
			if (id > 0) {

				//........................................................................
				//					READ THE DISTRIBUTIONS
				//		(read from opposite array due to previous swap operation)
				//........................................................................
				// even
				f2 = dist[10*Np+n];
				f4 = dist[11*Np+n];
				f6 = dist[12*Np+n];
				f8 = dist[13*Np+n];
				f10 = dist[14*Np+n];
				f12 = dist[15*Np+n];
				f14 = dist[16*Np+n];
				f16 = dist[17*Np+n];
				f18 = dist[18*Np+n];
				f0 = dist[n];
				// odd
				f1 = dist[Np+n];
				f3 = dist[2*Np+n];
				f5 = dist[3*Np+n];
				f7 = dist[4*Np+n];
				f9 = dist[5*Np+n];
				f11 = dist[6*Np+n];
				f13 = dist[7*Np+n];
				f15 = dist[8*Np+n];
				f17 = dist[9*Np+n];

				//........................................................................
				//					WRITE THE DISTRIBUTIONS
				// even
				//disteven[n] = f0;
				dist[Np+n] = f2;
				dist[2*Np+n] = f4;
				dist[3*Np+n] = f6;
				dist[4*Np+n] = f8;
				dist[5*Np+n] = f10;
				dist[6*Np+n] = f12;
				dist[7*Np+n] = f14;
				dist[8*Np+n] = f16;
				dist[9*Np+n] = f18;
				// odd
				dist[10*Np+n] = f1;
				dist[11*Np+n] = f3;
				dist[12*Np+n] = f5;
				dist[13*Np+n] = f7;
				dist[14*Np+n] = f9;
				dist[15*Np+n] = f11;
				dist[16*Np+n] = f13;
				dist[17*Np+n] = f15;
				dist[18*Np+n] = f17;
				//........................................................................
			}
		}
	}
}

//__global__ void dvc_ScaLBL_AAodd_Compact_OLD(char * ID, int *d_neighborList, double *d_disteven, double *d_distodd, int Np) {
//
//	int n, nn, nnn;
//	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
//	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
//	char id;
//	int nread,nwrite,naccess;
//	int S = Np/NBLOCKS/NTHREADS+1;
//
//	for (int s=0; s<S; s++){
//		//........Get 1-D index for this thread....................
//		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
//		id = ID[n];
//		if (n<Np) {
//			if (id > 0) {
//				// Convention
//				// nread - read neighbor location for distribution q
//				// nwrite - write neighbor location for distribution q
//				// read f1
//
//				// Compact:
//				//				nn = d_neighborList[n+2*0*Np];    // nread
//				//				nnn=d_neighborList[n+(2*0+1)*Np]; // nwrite
//				//				//				if ((!(nn<0)) && (!(nnn<0))) {
//				//				//printf("n=%d: nn=%d, nnn=%d\n",n,nn,nnn);
//				//				f1 = d_distodd[nnn + 0*Np];
//				//				f2 = d_disteven[nn + (0+1)*Np];
//				//				d_distodd[nnn+0*Np] = f2;
//				//				d_disteven[nn+(0+1)*Np] = f1;
//
//				// still need f0
//				//f0 = d_disteven[n];  // ?
//
//
//				nread = d_neighborList[n];     naccess = 10*Np;    // 1   "neighbor 0"
//				if (nread<0) { nread=n;        naccess = Np;     }
//				f2 = dist[nread + naccess];
//
//				nread = d_neighborList[n+Np];  naccess = Np;      // 2
//				if (nread<0) { nread=n;        naccess = 10*Np;  }
//				f1 = dist[nread + nacess];
//
////				nread = d_neighborList[n+2*Np]; // 3
////				if (nread<0) nread=n;
////				f4 = d_disteven[2*Np+nread];
////
////				nread = d_neighborList[n+3*Np]; // 4
////				if (nread<0) nread=n;
////				f3 = d_distodd[Np+nread];
////
////				nread = d_neighborList[n+4*Np]; // 5
////				if (nread<0) nread=n;
////				f6 = d_disteven[3*Np+nread];
////
////				nread = d_neighborList[n+5*Np]; // 6
////				if (nread<0) nread=n;
////				f5 = d_distodd[2*Np+nread];
////
////				nread = d_neighborList[n+6*Np];  // 7
////				if (nread<0) nread=n;
////				f8 = d_disteven[4*Np+nread];
////
////				nread = d_neighborList[n+7*Np]; // 8
////				if (nread<0) nread=n;
////				f7 = d_distodd[3*Np+nread];
////
////				nread = d_neighborList[n+8*Np]; // 9
////				if (nread<0) nread=n;
////				f10 = d_disteven[5*Np+nread];
////
////				nread = d_neighborList[n+9*Np]; // 10
////				if (nread<0) nread=n;
////				f9 = d_distodd[4*Np+nread];
////
////				nread = d_neighborList[n+10*Np]; // 11
////				if (nread<0) nread=n;
////				f12	 = d_disteven[6*Np+nread];
////
////				nread = d_neighborList[n+11*Np]; // 12
////				if (nread<0) nread=n;
////				f11 = d_distodd[5*Np+nread];
////
////				nread = d_neighborList[n+12*Np]; // 13
////				if (nread<0) nread=n;
////				f14 = d_disteven[7*Np+nread];
////
////				nread = d_neighborList[n+13*Np]; // 14
////				if (nread<0) nread=n;
////				f13 = d_distodd[6*Np+nread];
////
////				nread = d_neighborList[n+14*Np]; //15
////				if (nread<0) nread=n;
////				f16 = d_disteven[8*Np+nread];
////
////				nread = d_neighborList[n+15*Np];  //16
////				if (nread<0) nread=n;
////				f15 = d_distodd[7*Np+nread];
////
////				nread = d_neighborList[n+16*Np]; //17
////				if (nread<0) nread=n;
////				f18 = d_disteven[9*Np+nread];
////
////				nread = d_neighborList[n+17*Np]; // 18
////				if (nread<0) nread=n;
////				f17 = d_distodd[8*Np+nread];
//
//				// At this point we need all f0, f1, f2, f3, f4 ..., f18
//				// COLLISION
//
//
//
//				// write
//				nwrite = d_neighborList[n];  // 1
//				if (nwrite<0) nwrite=n;
//				d_disteven[Np+nwrite]   = f1;
//
//				nwrite = d_neighborList[Np+n];  // 2
//				if (nwrite<0) nwrite=n;
//				d_distodd[nwrite]       = f2;
//
////				nwrite = d_neighborList[2*Np+n];  // 3
////				if (nwrite<0) nwrite=n;
////				d_disteven[2*Np+nwrite] = f3;
////
////				nwrite = d_neighborList[3*Np+n];  // 4
////				if (nwrite<0) nwrite=n;
////				d_distodd[Np+nwrite]    = f4;
////
////				nwrite = d_neighborList[4*Np+n];  // 5
////				if (nwrite<0) nwrite=n;
////				d_disteven[3*Np+nwrite] = f5;
////
////				nwrite = d_neighborList[5*Np+n];  // 6
////				if (nwrite<0) nwrite=n;
////				d_distodd[2*Np+nwrite]  = f6;
////
////				nwrite = d_neighborList[6*Np+n];  // 7
////				if (nwrite<0) nwrite=n;
////				d_disteven[4*Np+nwrite] = f7;
////
////				nwrite = d_neighborList[7*Np+n];  // 8
////				if (nwrite<0) nwrite=n;
////				d_distodd[3*Np+nwrite]  = f8;
////
////				nwrite = d_neighborList[8*Np+n];  // 9
////				if (nwrite<0) nwrite=n;
////				d_disteven[5*Np+nwrite] = f9;
////
////				nwrite = d_neighborList[9*Np+n];  // 10
////				if (nwrite<0) nwrite=n;
////				d_distodd[4*Np+nwrite]  = f10;
////
////				nwrite = d_neighborList[10*Np+n];  // 11
////				if (nwrite<0) nwrite=n;
////				d_disteven[6*Np+nwrite] = f11;
////
////				nwrite = d_neighborList[11*Np+n];  // 12
////				if (nwrite<0) nwrite=n;
////				d_distodd[5*Np+nwrite]  = f12;
////
////				nwrite = d_neighborList[12*Np+n];  // 13
////				if (nwrite<0) nwrite=n;
////				d_disteven[7*Np+nwrite] = f13;
////
////				nwrite = d_neighborList[13*Np+n];  // 14
////				if (nwrite<0) nwrite=n;
////				d_distodd[6*Np+nwrite]  = f14;
////
////				nwrite = d_neighborList[14*Np+n];  // 15
////				if (nwrite<0) nwrite=n;
////				d_disteven[8*Np+nwrite] = f15;
////
////				nwrite = d_neighborList[15*Np+n];  // 16
////				if (nwrite<0) nwrite=n;
////				d_distodd[7*Np+nwrite]  = f16;
////
////				nwrite = d_neighborList[16*Np+n];  // 17
////				if (nwrite<0) nwrite=n;
////				d_disteven[9*Np+nwrite] = f17;
////
////				nwrite = d_neighborList[17*Np+n];  // 18
////				if (nwrite<0) nwrite=n;
////				d_distodd[8*Np+nwrite]  = f18;
//
//			}
//		}
//	}
//}

__global__  void dvc_ScaLBL_D3Q19_Swap(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz)
{
	int i,j,k,n,nn,N;
	// distributions
	char id;
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	N = Nx*Ny*Nz;

	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N){ 
			id = ID[n];
			if (id > 0){
				//.......Back out the 3-D indices for node n..............
				k = n/(Nx*Ny);
				j = (n-Nx*Ny*k)/Nx;
				i = n-Nx*Ny*k-Nx*j;
				//........................................................................
				// Retrieve even distributions from the local node (swap convention)
				//		f0 = disteven[n];  // Does not particupate in streaming
				f1 = distodd[n];
				f3 = distodd[N+n];
				f5 = distodd[2*N+n];
				f7 = distodd[3*N+n];
				f9 = distodd[4*N+n];
				f11 = distodd[5*N+n];
				f13 = distodd[6*N+n];
				f15 = distodd[7*N+n];
				f17 = distodd[8*N+n];
				//........................................................................

				//........................................................................
				// Retrieve odd distributions from neighboring nodes (swap convention)
				//........................................................................
				nn = n+1;							// neighbor index (pull convention)
				if (!(i+1<Nx))	nn -= Nx;			// periodic BC along the x-boundary
				//if (i+1<Nx){
				f2 = disteven[N+nn];					// pull neighbor for distribution 2
				if (f2 > 0.0){
					distodd[n] = f2;
					disteven[N+nn] = f1;
				}
				//}
				//........................................................................
				nn = n+Nx;							// neighbor index (pull convention)
				if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
				//if (j+1<Ny){
				f4 = disteven[2*N+nn];				// pull neighbor for distribution 4
				if (f4 > 0.0){
					distodd[N+n] = f4;
					disteven[2*N+nn] = f3;
					//	}
				}
				//........................................................................
				nn = n+Nx*Ny;						// neighbor index (pull convention)
				if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
				//if (k+1<Nz){
				f6 = disteven[3*N+nn];				// pull neighbor for distribution 6
				if (f6 > 0.0){
					distodd[2*N+n] = f6;
					disteven[3*N+nn] = f5;
					//	}
				}
				//........................................................................
				nn = n+Nx+1;						// neighbor index (pull convention)
				if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
				if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
				//if ((i+1<Nx) && (j+1<Ny)){
				f8 = disteven[4*N+nn];				// pull neighbor for distribution 8
				if (f8 > 0.0){
					distodd[3*N+n] = f8;
					disteven[4*N+nn] = f7;
					//	}
				}
				//........................................................................
				nn = n-Nx+1;						// neighbor index (pull convention)
				if (!(i+1<Nx))	nn -= Nx;		// periodic BC along the x-boundary
				if (j-1<0)		nn += Nx*Ny;	// Perioidic BC along the y-boundary
				//if (!(i-1<0) && (j+1<Ny)){
				f10 = disteven[5*N+nn];					// pull neighbor for distribution 9
				if (f10 > 0.0){
					distodd[4*N+n] = f10;
					disteven[5*N+nn] = f9;
					//	}
				}
				//........................................................................
				nn = n+Nx*Ny+1;						// neighbor index (pull convention)
				if (!(i+1<Nx))	nn -= Nx;		// periodic BC along the x-boundary
				if (!(k+1<Nz))	nn -= Nx*Ny*Nz;	// Perioidic BC along the z-boundary
				//if ( !(i-1<0) && !(k-1<0)){
				f12 = disteven[6*N+nn];				// pull distribution 11
				if (f12 > 0.0){
					distodd[5*N+n] = f12;
					disteven[6*N+nn] = f11;
					//	}
				}
				//........................................................................
				nn = n-Nx*Ny+1;						// neighbor index (pull convention)
				if (!(i+1<Nx))	nn -= Nx;		// periodic BC along the x-boundary
				if (k-1<0)		nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
				//if (!(i-1<0) && (k+1<Nz)){
				f14 = disteven[7*N+nn];				// pull neighbor for distribution 13
				if (f14 > 0.0){
					distodd[6*N+n] = f14;
					disteven[7*N+nn] = f13;
					//	}
				}
				//........................................................................
				nn = n+Nx*Ny+Nx;					// neighbor index (pull convention)
				if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
				if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
				//if (!(j-1<0) && !(k-1<0)){
				f16 = disteven[8*N+nn];				// pull neighbor for distribution 15
				if (f16 > 0.0){
					distodd[7*N+n] = f16;
					disteven[8*N+nn] = f15;
					//	}
				}
				//........................................................................
				nn = n-Nx*Ny+Nx;					// neighbor index (pull convention)
				if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
				if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
				//if (!(j-1<0) && (k+1<Nz)){
				f18 = disteven[9*N+nn];				// pull neighbor for distribution 17
				if (f18 > 0.0){
					distodd[8*N+n] = f18;
					disteven[9*N+nn] = f17;
					//	}
				}
				//........................................................................

			}
		}
	}
}


__global__  void dvc_ScaLBL_D3Q19_Velocity(char *ID, double *disteven, double *distodd, double *vel, int Nx, int Ny, int Nz)
{
	int n,N;
	// distributions
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double vx,vy,vz;
	char id;
	N = Nx*Ny*Nz;

	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N){
			id = ID[n];
			if (id==0){
				vel[n] = 0.0; vel[N+n] = 0.0; vel[2*N+n]=0.0;
				for(int q=0; q<9; q++){
					disteven[q*N+n] = -1.0;
					distodd[q*N+n] = -1.0;
				}
				disteven[9*N+n] = -1.0;
			}
			else{
				//........................................................................
				// Registers to store the distributions
				//........................................................................
				f2 = disteven[N+n];
				f4 = disteven[2*N+n];
				f6 = disteven[3*N+n];
				f8 = disteven[4*N+n];
				f10 = disteven[5*N+n];
				f12 = disteven[6*N+n];
				f14 = disteven[7*N+n];
				f16 = disteven[8*N+n];
				f18 = disteven[9*N+n];
				//........................................................................
				f1 = distodd[n];
				f3 = distodd[1*N+n];
				f5 = distodd[2*N+n];
				f7 = distodd[3*N+n];
				f9 = distodd[4*N+n];
				f11 = distodd[5*N+n];
				f13 = distodd[6*N+n];
				f15 = distodd[7*N+n];
				f17 = distodd[8*N+n];
				//.................Compute the velocity...................................
				vx = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
				vy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
				vz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
				//..................Write the velocity.....................................
				vel[n] = vx;
				vel[N+n] = vy;
				vel[2*N+n] = vz;
				//........................................................................
			}
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_Pressure(const char *ID, const double *disteven, const double *distodd,
		double *Pressure, int Nx, int Ny, int Nz)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	char id;
	N = Nx*Ny*Nz;

	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N){
			id = ID[n];
			if (id == 0)  Pressure[n] = 0.0;
			else{
				//.......................................................................
				// Registers to store the distributions
				//........................................................................
				f0 = disteven[n];
				f2 = disteven[N+n];
				f4 = disteven[2*N+n];
				f6 = disteven[3*N+n];
				f8 = disteven[4*N+n];
				f10 = disteven[5*N+n];
				f12 = disteven[6*N+n];
				f14 = disteven[7*N+n];
				f16 = disteven[8*N+n];
				f18 = disteven[9*N+n];
				//........................................................................
				f1 = distodd[n];
				f3 = distodd[1*N+n];
				f5 = distodd[2*N+n];
				f7 = distodd[3*N+n];
				f9 = distodd[4*N+n];
				f11 = distodd[5*N+n];
				f13 = distodd[6*N+n];
				f15 = distodd[7*N+n];
				f17 = distodd[8*N+n];
				//.................Compute the velocity...................................
				Pressure[n] = 0.3333333333333333*(f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+
						f9+f12+f11+f14+f13+f16+f15+f18+f17);
			}
		}
	}
}

__global__  void dvc_D3Q19_Velocity_BC_z(double *disteven, double *distodd, double uz,
		int Nx, int Ny, int Nz)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double din;

	N = Nx*Ny*Nz;
	n = Nx*Ny +  blockIdx.x*blockDim.x + threadIdx.x;

	if (n < 2*Nx*Ny){
		//........................................................................
		// Read distributions from "opposite" memory convention
		//........................................................................
		//........................................................................
		f1 = distodd[n];
		f3 = distodd[N+n];
		f5 = distodd[2*N+n];
		f7 = distodd[3*N+n];
		f9 = distodd[4*N+n];
		f11 = distodd[5*N+n];
		f13 = distodd[6*N+n];
		f15 = distodd[7*N+n];
		f17 = distodd[8*N+n];
		//........................................................................
		f0 = disteven[n];
		f2 = disteven[N+n];
		f4 = disteven[2*N+n];
		f6 = disteven[3*N+n];
		f8 = disteven[4*N+n];
		f10 = disteven[5*N+n];
		f12 = disteven[6*N+n];
		f14 = disteven[7*N+n];
		f16 = disteven[8*N+n];
		f18 = disteven[9*N+n];
		//...................................................

		// Determine the outlet flow velocity
		//	uz = 1.0 - (f0+f4+f3+f2+f1+f8+f7+f9+f10 +
		//			2*(f5+f15+f18+f11+f14))/din;
		din = (f0+f4+f3+f2+f1+f8+f7+f9+f10+2*(f5+f15+f18+f11+f14))/(1.0-uz);
		// Set the unknown distributions:
		f6 = f5 + 0.3333333333333333*din*uz;
		f16 = f15 + 0.1666666666666667*din*uz;
		f17 = f16 + f4 - f3-f15+f18+f8-f7	+f9-f10;
		f12= (din*uz+f5+ f15+f18+f11+f14-f6-f16-f17-f2+f1-f14+f11-f8+f7+f9-f10)*0.5;
		f13= din*uz+f5+ f15+f18+f11+f14-f6-f16-f17-f12;

		//........Store in "opposite" memory location..........
		disteven[3*N+n] = f6;
		disteven[6*N+n] = f12;
		distodd[6*N+n] = f13;
		disteven[8*N+n] = f16;
		distodd[8*N+n] = f17;
		//...................................................
	}
}

__global__ void dvc_D3Q19_Velocity_BC_Z(double *disteven, double *distodd, double uz,
		int Nx, int Ny, int Nz, int outlet){
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double dout;

	N = Nx*Ny*Nz;
	n = outlet +  blockIdx.x*blockDim.x + threadIdx.x;

	// Loop over the boundary - threadblocks delineated by start...finish
	if ( n<N-Nx*Ny ){
		// Read distributions from "opposite" memory convention
		//........................................................................
		f1 = distodd[n];
		f3 = distodd[N+n];
		f5 = distodd[2*N+n];
		f7 = distodd[3*N+n];
		f9 = distodd[4*N+n];
		f11 = distodd[5*N+n];
		f13 = distodd[6*N+n];
		f15 = distodd[7*N+n];
		f17 = distodd[8*N+n];
		//........................................................................
		f0 = disteven[n];
		f2 = disteven[N+n];
		f4 = disteven[2*N+n];
		f6 = disteven[3*N+n];
		f8 = disteven[4*N+n];
		f10 = disteven[5*N+n];
		f12 = disteven[6*N+n];
		f14 = disteven[7*N+n];
		f16 = disteven[8*N+n];
		f18 = disteven[9*N+n];
		//uz = -1.0 + (f0+f4+f3+f2+f1+f8+f7+f9+f10 + 2*(f6+f16+f17+f12+f13))/dout;
		dout = (f0+f4+f3+f2+f1+f8+f7+f9+f10 + 2*(f6+f16+f17+f12+f13))/(1.0+uz);
		f5 = f6 - 0.33333333333333338*dout* uz;
		f15 = f16 - 0.16666666666666678*dout* uz;
		f18 = f15 - f4 + f3-f16+f17-f8+f7-f9+f10;
		f11 = (-dout*uz+f6+ f16+f17+f12+f13-f5-f15-f18+f2-f1-f13+f12+f8-f7-f9+f10)*0.5;
		f14 = -dout*uz+f6+ f16+f17+f12+f13-f5-f15-f18-f11;
		//........Store in "opposite" memory location..........
		distodd[2*N+n] = f5;
		distodd[5*N+n] = f11;
		disteven[7*N+n] = f14;
		distodd[7*N+n] = f15;
		disteven[9*N+n] = f18;
		//...................................................
	}
}

__global__ void dvc_D3Q19_Flux_BC_z(double *disteven, double *distodd, double flux, double *dvcsum,
		int Nx, int Ny, int Nz){
	// Note that this routine assumes the distributions are stored "opposite"
	// odd distributions in disteven and even distributions in distodd.
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f6,f7,f8,f9;
	double f10,f12,f13,f16,f17;

	//double A = 1.f*double(Nx*Ny);
	double factor = 1.f/(double(Nx*Ny)*(1.0-flux));

	double sum = 0.f;

	N = Nx*Ny*Nz;
	n = Nx*Ny +  blockIdx.x*blockDim.x + threadIdx.x;

	if (n < 2*Nx*Ny){

		//........................................................................
		f1 = distodd[n];
		f3 = distodd[N+n];
//		f5 = distodd[2*N+n];
		f7 = distodd[3*N+n];
		f9 = distodd[4*N+n];
//		f11 = distodd[5*N+n];
		f13 = distodd[6*N+n];
//		f15 = distodd[7*N+n];
		f17 = distodd[8*N+n];
		//........................................................................
		f0 = disteven[n];
		f2 = disteven[N+n];
		f4 = disteven[2*N+n];
		f6 = disteven[3*N+n];
		f8 = disteven[4*N+n];
		f10 = disteven[5*N+n];
		f12 = disteven[6*N+n];
//		f14 = disteven[7*N+n];
		f16 = disteven[8*N+n];
//		f18 = disteven[9*N+n];
		//...................................................
		// compute local sum to determine the density value to set pressure
		//sum = (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17))/(A*(1.0-flux));
		sum = factor*(f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));
		//localsum[n]=sum;
	}

	//sum = warpReduceSum(sum);
	//if (threadIdx.x & (warpSize-1) == 0 ){
	//   atomicAdd(dvcsum,sum);
	//}

	sum = blockReduceSum(sum);
	if (threadIdx.x==0)
	   atomicAdd(dvcsum, sum);
}

__global__ void dvc_D3Q19_Flux_BC_Z(double *disteven, double *distodd, double flux, double *dvcsum,
		int Nx, int Ny, int Nz, int outlet){
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f7,f8,f9;
	double f10,f11,f14,f15,f18;

	N = Nx*Ny*Nz;
	n = outlet +  blockIdx.x*blockDim.x + threadIdx.x;

	double factor = 1.f/(double(Nx*Ny)*(1.0+flux));
	double sum = 0.f;

	// Loop over the boundary - threadblocks delineated by start...finish
	if ( n<N-Nx*Ny ){
		//........................................................................
		// Read distributions from "opposite" memory convention
		//........................................................................
		f1 = distodd[n];
		f3 = distodd[N+n];
		f5 = distodd[2*N+n];
		f7 = distodd[3*N+n];
		f9 = distodd[4*N+n];
		f11 = distodd[5*N+n];
//		f13 = distodd[6*N+n];
		f15 = distodd[7*N+n];
//		f17 = distodd[8*N+n];
		//........................................................................
		f0 = disteven[n];
		f2 = disteven[N+n];
		f4 = disteven[2*N+n];
//		f6 = disteven[3*N+n];
		f8 = disteven[4*N+n];
		f10 = disteven[5*N+n];
//		f12 = disteven[6*N+n];
		f14 = disteven[7*N+n];
//		f16 = disteven[8*N+n];
		f18 = disteven[9*N+n];

		// Local sum (based on the consistency condition)
		//sum = (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18))/(A*(1.0+flux));
		sum = factor*(f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));
		//localsum[n]=sum;
	}

	sum = blockReduceSum(sum);
	if (threadIdx.x==0)
		atomicAdd(dvcsum, sum);

}


__global__  void dvc_ScaLBL_D3Q19_Pressure_BC_z(double *disteven, double *distodd, double din,
		int Nx, int Ny, int Nz)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz,Cyz,Cxz;

	N = Nx*Ny*Nz;
	n = Nx*Ny +  blockIdx.x*blockDim.x + threadIdx.x;

	if (n < 2*Nx*Ny){

		//........................................................................
		// Read distributions from "opposite" memory convention
		//........................................................................
		//........................................................................
		f1 = distodd[n];
		f3 = distodd[N+n];
		f5 = distodd[2*N+n];
		f7 = distodd[3*N+n];
		f9 = distodd[4*N+n];
		f11 = distodd[5*N+n];
		f13 = distodd[6*N+n];
		f15 = distodd[7*N+n];
		f17 = distodd[8*N+n];
		//........................................................................
		f0 = disteven[n];
		f2 = disteven[N+n];
		f4 = disteven[2*N+n];
		f6 = disteven[3*N+n];
		f8 = disteven[4*N+n];
		f10 = disteven[5*N+n];
		f12 = disteven[6*N+n];
		f14 = disteven[7*N+n];
		f16 = disteven[8*N+n];
		f18 = disteven[9*N+n];
		//...................................................
		//........Determine the inlet flow velocity.........
		//			uz = -1 + (f0+f3+f4+f1+f2+f7+f8+f10+f9
		//					   + 2*(f5+f15+f18+f11+f14))/din;
		//........Set the unknown distributions..............
		//			f6 = f5 - 0.3333333333333333*din*uz;
		//			f16 = f15 - 0.1666666666666667*din*uz;
		//			f17 = f16 - f3 + f4-f15+f18-f7+f8-f10+f9;
		//			f12= 0.5*(-din*uz+f5+f15+f18+f11+f14-f6-f16-
		//					  f17+f1-f2-f14+f11+f7-f8-f10+f9);
		//			f13= -din*uz+f5+f15+f18+f11+f14-f6-f16-f17-f12;
		// Determine the inlet flow velocity
		ux = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14);
		uy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18);
		uz = din - (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f5 = f6 + 0.33333333333333338*uz;
		f11 = f12 + 0.16666666666666678*(uz+ux)-Cxz;
		f14 = f13 + 0.16666666666666678*(uz-ux)+Cxz;
		f15 = f16 + 0.16666666666666678*(uy+uz)-Cyz;
		f18 = f17 + 0.16666666666666678*(uz-uy)+Cyz;
		//........Store in "opposite" memory location..........
		distodd[2*N+n] = f5;
		distodd[5*N+n] = f11;
		disteven[7*N+n] = f14;
		distodd[7*N+n] = f15;
		disteven[9*N+n] = f18;
	}
}

__global__ void dvc_ScaLBL_D3Q19_Init_Simple(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz)
{
	int n,N;
	N = Nx*Ny*Nz;
	char id;
	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N ){
			id = ID[n];
			if (id > 0 ){
				f_even[n] = 0 + 0.01*0;
				f_odd[n] = 0+ 0.01*1;		//double(100*n)+1.f;
				f_even[N+n] = 1+ 0.01*2;	//double(100*n)+2.f;
				f_odd[N+n] = 1+ 0.01*3;	//double(100*n)+3.f;
				f_even[2*N+n] = 2+ 0.01*4;	//double(100*n)+4.f;
				f_odd[2*N+n] = 2+ 0.01*5;	//double(100*n)+5.f;
				f_even[3*N+n] = 3+ 0.01*6;	//double(100*n)+6.f;
				f_odd[3*N+n] = 3+ 0.01*7;   //double(100*n)+7.f;
				f_even[4*N+n] = 4+ 0.01*8;   //double(100*n)+8.f;
				f_odd[4*N+n] = 4+ 0.01*9;   //double(100*n)+9.f;
				f_even[5*N+n] = 5+ 0.01*10;  //double(100*n)+10.f;
				f_odd[5*N+n] = 5+ 0.01*11;  //double(100*n)+11.f;
				f_even[6*N+n] = 6+ 0.01*12;  //double(100*n)+12.f;
				f_odd[6*N+n] = 6+ 0.01*13;  //double(100*n)+13.f;
				f_even[7*N+n] = 7+ 0.01*14;  //double(100*n)+14.f;
				f_odd[7*N+n] = 7+ 0.01*15;  //double(100*n)+15.f;
				f_even[8*N+n] = 8+ 0.01*16;  //double(100*n)+16.f;
				f_odd[8*N+n] = 8+ 0.01*17;  //double(100*n)+17.f;
				f_even[9*N+n] = 9+ 0.01*18;  //double(100*n)+18.f;
			}
			else{
				for(int q=0; q<9; q++){
					f_even[q*N+n] = -1.0;
					f_odd[q*N+n] = -1.0;
				}
				f_even[9*N+n] = -1.0;
			}
		}
	}
}

//__global__ void dvc_ScaLBL_AAeven_Compact(char * ID, double *disteven, double *distodd, int Np) {
//
//
//	int q,n,nn;
//	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
//	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
//
//	char id;
//	int S = Np/NBLOCKS/NTHREADS + 1;
//	for (int s=0; s<S; s++){
//		//........Get 1-D index for this thread....................
//		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
//		id = ID[n];
//
//		if ( n<Np ){
//			if (id > 0) {
//
//				//........................................................................
//				//					READ THE DISTRIBUTIONS
//				//		(read from opposite array due to previous swap operation)
//				//........................................................................
//				// even
//				f2 = distodd[n];
//				f4 = distodd[Np+n];
//				f6 = distodd[2*Np+n];
//				f8 = distodd[3*Np+n];
//				f10 = distodd[4*Np+n];
//				f12 = distodd[5*Np+n];
//				f14 = distodd[6*Np+n];
//				f16 = distodd[7*Np+n];
//				f18 = distodd[8*Np+n];
//				f0 = disteven[n];
//				// odd
//				f1 = disteven[Np+n];
//				f3 = disteven[2*Np+n];
//				f5 = disteven[3*Np+n];
//				f7 = disteven[4*Np+n];
//				f9 = disteven[5*Np+n];
//				f11 = disteven[6*Np+n];
//				f13 = disteven[7*Np+n];
//				f15 = disteven[8*Np+n];
//				f17 = disteven[9*Np+n];
//
//				//........................................................................
//				//					WRITE THE DISTRIBUTIONS
//				// even
//				//disteven[n] = f0;
//				disteven[Np+n] = f2;
//				disteven[2*Np+n] = f4;
//				disteven[3*Np+n] = f6;
//				disteven[4*Np+n] = f8;
//				disteven[5*Np+n] = f10;
//				disteven[6*Np+n] = f12;
//				disteven[7*Np+n] = f14;
//				disteven[8*Np+n] = f16;
//				disteven[9*Np+n] = f18;
//				// odd
//				distodd[n] = f1;
//				distodd[Np+n] = f3;
//				distodd[2*Np+n] = f5;
//				distodd[3*Np+n] = f7;
//				distodd[4*Np+n] = f9;
//				distodd[5*Np+n] = f11;
//				distodd[6*Np+n] = f13;
//				distodd[7*Np+n] = f15;
//				distodd[8*Np+n] = f17;
//				//........................................................................
//			}
//		}
//	}
//}



__global__  void dvc_ScaLBL_D3Q19_Pressure_BC_Z(double *disteven, double *distodd, double dout,
		int Nx, int Ny, int Nz, int outlet)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz,Cyz,Cxz;

	N = Nx*Ny*Nz;
	n = outlet +  blockIdx.x*blockDim.x + threadIdx.x;

	// Loop over the boundary - threadblocks delineated by start...finish
	if ( n<N-Nx*Ny ){

		//........................................................................
		// Read distributions from "opposite" memory convention
		//........................................................................
		f1 = distodd[n];
		f3 = distodd[N+n];
		f5 = distodd[2*N+n];
		f7 = distodd[3*N+n];
		f9 = distodd[4*N+n];
		f11 = distodd[5*N+n];
		f13 = distodd[6*N+n];
		f15 = distodd[7*N+n];
		f17 = distodd[8*N+n];
		//........................................................................
		f0 = disteven[n];
		f2 = disteven[N+n];
		f4 = disteven[2*N+n];
		f6 = disteven[3*N+n];
		f8 = disteven[4*N+n];
		f10 = disteven[5*N+n];
		f12 = disteven[6*N+n];
		f14 = disteven[7*N+n];
		f16 = disteven[8*N+n];
		f18 = disteven[9*N+n];
		//........Determine the outlet flow velocity.........
		//			uz = 1 - (f0+f3+f4+f1+f2+f7+f8+f10+f9+
		//					  2*(f6+f16+f17+f12+f13))/dout;
		//...................................................
		//........Set the Unknown Distributions..............
		//			f5 = f6 + 0.33333333333333338*dout*uz;
		//			f15 = f16 + 0.16666666666666678*dout*uz;
		//			f18 = f15+f3-f4-f16+f17+f7-f8+f10-f9;
		//			f11= 0.5*(dout*uz+f6+ f16+f17+f12+f13-f5
		//				  -f15-f18-f1+f2-f13+f12-f7+f8+f10-f9);
		//			f14= dout*uz+f6+ f16+f17+f12+f13-f5-f15-f18-f11;
		// Determine the outlet flow velocity
		//ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
		//uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		//uz = -1.0 + (f0+f4+f3+f2+f1+f8+f7+f9+f10 + 2*(f6+f16+f17+f12+f13))/dout;

		// Determine the inlet flow velocity
		ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
		uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		uz = -dout + (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f6 = f5 - 0.33333333333333338*uz;
		f12 = f11 - 0.16666666666666678*(uz+ux)+Cxz;
		f13 = f14 - 0.16666666666666678*(uz-ux)-Cxz;
		f16 = f15 - 0.16666666666666678*(uy+uz)+Cyz;
		f17 = f18 - 0.16666666666666678*(uz-uy)-Cyz;

		//........Store in "opposite" memory location..........
		disteven[3*N+n] = f6;
		disteven[6*N+n] = f12;
		distodd[6*N+n] = f13;
		disteven[8*N+n] = f16;
		distodd[8*N+n] = f17;
		//...................................................
	}
}
//*************************************************************************

//extern "C" void ScaLBL_D3Q19_MapRecv(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
//			int *d3q19_recvlist, int Nx, int Ny, int Nz){
//	int GRID = count / 512 + 1;
//	dvc_ScaLBL_D3Q19_Unpack <<<GRID,512 >>>(q, Cqx, Cqy, Cqz, list, start, count, d3q19_recvlist, Nx, Ny, Nz);
//}

extern "C" void ScaLBL_D3Q19_Pack(int q, int *list, int start, int count, double *sendbuf, double *dist, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_Pack <<<GRID,512 >>>(q, list, start, count, sendbuf, dist, N);
}
extern "C" void ScaLBL_D3Q19_Unpack(int q, int *list,  int start, int count, double *recvbuf, double *dist, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_Unpack <<<GRID,512 >>>(q, list, start, count, recvbuf, dist, N);
}
//*************************************************************************
extern "C" void ScaLBL_D3Q19_Init(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz){
	//	dvc_ScaLBL_D3Q19_Init<<<NBLOCKS,NTHREADS >>>(ID, f_even, f_odd, Nx, Ny, Nz);
	dvc_ScaLBL_D3Q19_Init_Simple<<<NBLOCKS,NTHREADS >>>(ID, f_even, f_odd, Nx, Ny, Nz);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Init: %s \n",cudaGetErrorString(err));
	}

}
extern "C" void ScaLBL_D3Q19_Swap(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz){
	dvc_ScaLBL_D3Q19_Swap<<<NBLOCKS,NTHREADS >>>(ID, disteven, distodd, Nx, Ny, Nz);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Swap: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_Swap_Compact(int *neighborList, double *disteven, double *distodd, int Np)
{

	const int Q = 9;
	//	cudaStream_t streams[Q];
	// Launch the swap operation as different kernels
	for (int q=0; q<Q; q++){
		dvc_ScaLBL_D3Q19_Swap_Compact<<<NBLOCKS,NTHREADS >>>(neighborList, disteven, distodd, Np, q);
	}
	// cpu should wait for all kernels to finish (to avoid launch of dependent kernels)
	//cudaDeviceSynchronize();
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Swap: %s \n",cudaGetErrorString(err));
	}
}

//extern "C" void ScaLBL_D3Q19_AAeven_Compact(char * ID, double *d_disteven, double *d_distodd, int Np) {
//	dvc_ScaLBL_AAeven_Compact<<<NBLOCKS,NTHREADS>>>(ID, d_disteven, d_distodd,Np);
//	cudaError_t err = cudaGetLastError();
//	if (cudaSuccess != err){
//		printf("CUDA error in ScaLBL_D3Q19_Init: %s \n",cudaGetErrorString(err));
//	}
//}

extern "C" void ScaLBL_D3Q19_AAeven_Compact(char * ID, double *d_dist,  int Np) {
	dvc_ScaLBL_AAeven_Compact<<<NBLOCKS,NTHREADS>>>(ID, d_dist, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Compact(char * ID, int *d_neighborList, double *d_dist, int Np) {
	dvc_ScaLBL_AAodd_Compact<<<NBLOCKS,NTHREADS>>>(ID,d_neighborList, d_dist,Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Init: %s \n",cudaGetErrorString(err));
	}
}


extern "C" void ScaLBL_D3Q19_Velocity(char *ID, double *disteven, double *distodd, double *vel, int Nx, int \
		Ny, int Nz){

	dvc_ScaLBL_D3Q19_Velocity<<<NBLOCKS,NTHREADS >>>(ID, disteven, distodd, vel, Nx, Ny, Nz);
}
extern "C" void ScaLBL_D3Q19_Pressure(char *ID, double *disteven, double *distodd, double *Pressure,
		int Nx, int Ny, int Nz){
	dvc_ScaLBL_D3Q19_Pressure<<< NBLOCKS,NTHREADS >>>(ID, disteven, distodd, Pressure, Nx, Ny, Nz);
}

extern "C" void ScaLBL_D3Q19_Velocity_BC_z(double *disteven, double *distodd, double uz,int Nx, int Ny, int Nz){
	int GRID = Nx*Ny / 512 + 1;
	dvc_D3Q19_Velocity_BC_z<<<GRID,512>>>(disteven,distodd, uz, Nx, Ny, Nz);
}

extern "C" void ScaLBL_D3Q19_Velocity_BC_Z(double *disteven, double *distodd, double uz, int Nx, int Ny, int Nz, int outlet){
	int GRID = Nx*Ny / 512 + 1;
	dvc_D3Q19_Velocity_BC_Z<<<GRID,512>>>(disteven, distodd, uz, Nx, Ny, Nz, outlet);
}




extern "C" double ScaLBL_D3Q19_Flux_BC_z(double *disteven, double *distodd, double flux,int Nx, int Ny, int Nz){

	int GRID = Nx*Ny / 512 + 1;

	// IMPORTANT -- this routine may fail if Nx*Ny > 512*512
	if (Nx*Ny > 512*512){
		printf("WARNING (ScaLBL_D3Q19_Flux_BC_Z): CUDA reduction operation may fail if Nx*Ny > 512*512");
	}

	// Allocate memory to store the sums
	double din;
	double sum[1];
 	double *dvcsum;
	cudaMalloc((void **)&dvcsum,sizeof(double)*Nx*Ny);
	cudaMemset(dvcsum,0,sizeof(double)*Nx*Ny);

	// compute the local flux and store the result
	dvc_D3Q19_Flux_BC_z<<<GRID,512>>>(disteven, distodd, flux, dvcsum, Nx, Ny, Nz);

	// Now read the total flux
	cudaMemcpy(&sum[0],dvcsum,sizeof(double),cudaMemcpyDeviceToHost);
	din=sum[0];

	// free the memory needed for reduction
	cudaFree(dvcsum);

	return din;
}

extern "C" double ScaLBL_D3Q19_Flux_BC_Z(double *disteven, double *distodd, double flux, int Nx, int Ny, int Nz, int outlet){

	int GRID = Nx*Ny / 512 + 1;

	// IMPORTANT -- this routine may fail if Nx*Ny > 512*512
	if (Nx*Ny > 512*512){
		printf("WARNING (ScaLBL_D3Q19_Flux_BC_Z): CUDA reduction operation may fail if Nx*Ny > 512*512");
	}

	// Allocate memory to store the sums
	double dout;
	double sum[1];
 	double *dvcsum;
	cudaMalloc((void **)&dvcsum,sizeof(double)*Nx*Ny);
	cudaMemset(dvcsum,0,sizeof(double)*Nx*Ny);

	// compute the local flux and store the result
	dvc_D3Q19_Flux_BC_Z<<<GRID,512>>>(disteven, distodd, flux, dvcsum, Nx, Ny, Nz, outlet);

	// Now read the total flux
	cudaMemcpy(&sum[0],dvcsum,sizeof(double),cudaMemcpyDeviceToHost);

	// free the memory needed for reduction

	dout = sum[0];

	cudaFree(dvcsum);

	return dout;

}


extern "C" double deviceReduce(double *in, double* out, int N) {
	int threads = 512;
	int blocks = min((N + threads - 1) / threads, 1024);

	double sum = 0.f;
	deviceReduceKernel<<<blocks, threads>>>(in, out, N);
	deviceReduceKernel<<<1, 1024>>>(out, out, blocks);
	return sum;
}

extern "C" void ScaLBL_D3Q19_Pressure_BC_z(double *disteven, double *distodd, double din, int Nx, int Ny, int Nz){
	int GRID = Nx*Ny / 512 + 1;
	dvc_ScaLBL_D3Q19_Pressure_BC_z<<<GRID,512>>>(disteven, distodd, din, Nx, Ny, Nz);
}
extern "C" void ScaLBL_D3Q19_Pressure_BC_Z(double *disteven, double *distodd, double dout,
		int Nx, int Ny, int Nz, int outlet){
	int GRID = Nx*Ny / 512 + 1;
	dvc_ScaLBL_D3Q19_Pressure_BC_Z<<<GRID,512>>>(disteven, distodd, dout, Nx, Ny, Nz, outlet);
}
