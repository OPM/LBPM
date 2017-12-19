#include <stdio.h>

#define NBLOCKS 1024
#define NTHREADS 128


// functionality for parallel reduction in Flux BC routines -- probably should be re-factored to another location
// functions copied from https://devblogs.nvidia.com/parallelforall/faster-parallel-reductions-kepler/

__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;

    unsigned long long int old = *address_as_ull, assumed;

    do{ assumed = old;
    	old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val +__longlong_as_double(assumed)));
    } while (assumed != old);

    return __longlong_as_double(old);
}

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
	// Unack distribution from the recv buffer
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

__global__ void dvc_D3Q19_Flux_BC_z(char *ID, double *disteven, double *distodd, double flux, double *dvcsum,
								  int Nx, int Ny, int Nz){
	// Note that this routine assumes the distributions are stored "opposite"
	// odd distributions in disteven and even distributions in distodd.
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f6,f7,f8,f9;
	double f10,f12,f13,f16,f17;

	//double A = 1.f*double(Nx*Ny);
	double factor = 1.f/((1.0-flux));

	double sum = 0.f;

	N = Nx*Ny*Nz;
	n = Nx*Ny +  blockIdx.x*blockDim.x + threadIdx.x;

	if (n < 2*Nx*Ny){
	     char id=ID[n];
	     if (id > 0){
		//........................................................................
		f2 = distodd[n];
		f4 = distodd[N+n];
		f6 = distodd[2*N+n];
		f8 = distodd[3*N+n];
		f10 = distodd[4*N+n];
		f12 = distodd[5*N+n];
		//f14 = distodd[6*N+n];
		f16 = distodd[7*N+n];
		//f18 = distodd[8*N+n];
		//........................................................................
		f0 = disteven[n];
		f1 = disteven[N+n];
		f3 = disteven[2*N+n];
		//f5 = disteven[3*N+n];
		f7 = disteven[4*N+n];
		f9 = disteven[5*N+n];
		//f11 = disteven[6*N+n];
		f13 = disteven[7*N+n];
		//f15 = disteven[8*N+n];
		f17 = disteven[9*N+n];
		//...................................................
		// compute local sum to determine the density value to set pressure
		//sum = (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17))/(A*(1.0-flux));
		sum = factor*(f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));
		//localsum[n]=sum;
		}
	}

	//atomicAdd(dvcsum, sum);


	//sum = warpReduceSum(sum);
	//if (threadIdx.x & (warpSize-1) == 0 ){
	//   atomicAdd(dvcsum,sum);
	//}

	sum = blockReduceSum(sum);
	if (threadIdx.x==0)
	   atomicAdd(dvcsum, sum);
}

__global__ void dvc_D3Q19_Flux_BC_Z(char *ID, double *disteven, double *distodd, double flux, double *dvcsum,
								   int Nx, int Ny, int Nz, int outlet){
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f7,f8,f9;
	double f10,f11,f14,f15,f18;

	N = Nx*Ny*Nz;
	n = outlet +  blockIdx.x*blockDim.x + threadIdx.x;

	double factor = 1.f/(1.0+flux);
	double sum = 0.f;

	// Loop over the boundary - threadblocks delineated by start...finish
	if ( n<N-Nx*Ny ){
	    char id=ID[n];
	    if (id>0){
   		//........................................................................
		// Read distributions from "opposite" memory convention
		//........................................................................
		f2 = distodd[n];
		f4 = distodd[N+n];
		//f6 = distodd[2*N+n];
		f8 = distodd[3*N+n];
		f10 = distodd[4*N+n];
		//f12 = distodd[5*N+n];
		f14 = distodd[6*N+n];
		//f16 = distodd[7*N+n];
		f18 = distodd[8*N+n];
		//........................................................................
		f0 = disteven[n];
		f1 = disteven[N+n];
		f3 = disteven[2*N+n];
		f5 = disteven[3*N+n];
		f7 = disteven[4*N+n];
		f9 = disteven[5*N+n];
		f11 = disteven[6*N+n];
		//f13 = disteven[7*N+n];
		f15 = disteven[8*N+n];
		//f17 = disteven[9*N+n];

		// Local sum (based on the consistency condition)
		//sum = (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18))/(A*(1.0+flux));
		sum = factor*(f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));
		//localsum[n]=sum;
		}
	}

	//sum = warpReduceSum(sum);
//	if (threadIdx.x & (warpSize-1) == 0 ){
	//   atomicAdd(dvcsum,sum);
//	}

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
            f2 = distodd[n];
            f4 = distodd[N+n];
            f6 = distodd[2*N+n];
            f8 = distodd[3*N+n];
            f10 = distodd[4*N+n];
            f12 = distodd[5*N+n];
            f14 = distodd[6*N+n];
            f16 = distodd[7*N+n];
            f18 = distodd[8*N+n];
            //........................................................................
            f0 = disteven[n];
            f1 = disteven[N+n];
            f3 = disteven[2*N+n];
            f5 = disteven[3*N+n];
            f7 = disteven[4*N+n];
            f9 = disteven[5*N+n];
            f11 = disteven[6*N+n];
            f13 = disteven[7*N+n];
            f15 = disteven[8*N+n];
            f17 = disteven[9*N+n];
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
            disteven[3*N+n] = f5;
            disteven[6*N+n] = f11;
            distodd[6*N+n] = f14;
            disteven[8*N+n] = f15;
            distodd[8*N+n] = f18;
		}
	}

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
		f2 = distodd[n];
		f4 = distodd[N+n];
		f6 = distodd[2*N+n];
		f8 = distodd[3*N+n];
		f10 = distodd[4*N+n];
		f12 = distodd[5*N+n];
		f14 = distodd[6*N+n];
		f16 = distodd[7*N+n];
		f18 = distodd[8*N+n];
		//........................................................................
		f0 = disteven[n];
		f1 = disteven[N+n];
		f3 = disteven[2*N+n];
		f5 = disteven[3*N+n];
		f7 = disteven[4*N+n];
		f9 = disteven[5*N+n];
		f11 = disteven[6*N+n];
		f13 = disteven[7*N+n];
		f15 = disteven[8*N+n];
		f17 = disteven[9*N+n];
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
		distodd[2*N+n] = f6;
		distodd[5*N+n] = f12;
		disteven[7*N+n] = f13;
		distodd[7*N+n] = f16;
		disteven[9*N+n] = f17;
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
	dvc_ScaLBL_D3Q19_Init<<<NBLOCKS,NTHREADS >>>(ID, f_even, f_odd, Nx, Ny, Nz);
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
