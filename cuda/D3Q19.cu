#include <stdio.h>
#include <cooperative_groups.h>

#define NBLOCKS 1024
#define NTHREADS 256

/*
1. constants that are known at compile time should be defined using preprocessor macros (e.g. #define) or via C/C++ const variables at global/file scope.
2. Usage of __constant__ memory may be beneficial for programs who use certain values that don't change for the duration of the kernel and for which certain access patterns are present (e.g. all threads access the same value at the same time). This is not better or faster than constants that satisfy the requirements of item 1 above.
3. If the number of choices to be made by a program are relatively small in number, and these choices affect kernel execution, one possible approach for additional compile-time optimization would be to use templated code/kernels
 */

__constant__ __device__ double mrt_V1=0.05263157894736842;
__constant__ __device__ double mrt_V2=0.012531328320802;
__constant__ __device__ double mrt_V3=0.04761904761904762;
__constant__ __device__ double mrt_V4=0.004594820384294068;
__constant__ __device__ double mrt_V5=0.01587301587301587;
__constant__ __device__ double mrt_V6=0.0555555555555555555555555;
__constant__ __device__ double mrt_V7=0.02777777777777778;
__constant__ __device__ double mrt_V8=0.08333333333333333;
__constant__ __device__ double mrt_V9=0.003341687552213868;
__constant__ __device__ double mrt_V10=0.003968253968253968;
__constant__ __device__ double mrt_V11=0.01388888888888889;
__constant__ __device__ double mrt_V12=0.04166666666666666;


// functionality for parallel reduction in Flux BC routines -- probably should be re-factored to another location
// functions copied from https://devblogs.nvidia.com/parallelforall/faster-parallel-reductions-kepler/

//__shared__ double Transform[722]=
//	   {};

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* address, double val) { 
   unsigned long long int* address_as_ull = (unsigned long long int*)address;
   unsigned long long int old = *address_as_ull, assumed;

   do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val+__longlong_as_double(assumed)));
   } while (assumed != old);
   return __longlong_as_double(old);
}
#endif

using namespace cooperative_groups;
__device__ double reduce_sum(thread_group g, double *temp, double val)
{
    int lane = g.thread_rank();

    // Each iteration halves the number of active threads
    // Each thread adds its partial sum[i] to sum[lane+i]
    for (int i = g.size() / 2; i > 0; i /= 2)
    {
        temp[lane] = val;
        g.sync(); // wait for all threads to store
        if(lane<i) val += temp[lane + i];
        g.sync(); // wait for all threads to load
    }
    return val; // note: only thread 0 will return full sum
}

__device__ double thread_sum(double *input, double n) 
{
    double sum = 0;

    for(int i = blockIdx.x * blockDim.x + threadIdx.x;
        i < n / 4; 
        i += blockDim.x * gridDim.x)
    {
        int4 in = ((int4*)input)[i];
        sum += in.x + in.y + in.z + in.w;
    }
    return sum;
}

__global__ void sum_kernel_block(double *sum, double *input, int n)
{
	double my_sum = thread_sum(input, n);

    extern __shared__ double temp[];
    thread_group g = this_thread_block();
    double block_sum = reduce_sum(g, temp, my_sum);

    if (g.thread_rank() == 0) atomicAdd(sum, block_sum);
}

__inline__ __device__
double warpReduceSum(double val) {
	for (int offset = warpSize/2; offset > 0; offset /= 2)
		val += __shfl_down_sync(0xFFFFFFFF, val, offset, 32);
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

__global__ void dvc_ScaLBL_D3Q19_Pack(int q, int *list, int start, int count, double *sendbuf, double *dist, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		n = list[idx];
		sendbuf[start+idx] = dist[q*N+n];
		//printf("%f \n",dist[q*N+n]);
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
		if (!(n<0)) { dist[q*N+n] = recvbuf[start+idx];
		//printf("%f \n",,dist[q*N+n]);
		}
	}
}

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

__global__ void dvc_ScaLBL_D3Q19_AA_Init(double *f_even, double *f_odd, int Np)
{
	int n;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<Np ){
			f_even[n] = 0.3333333333333333;
			f_odd[n] = 0.055555555555555555;                //double(100*n)+1.f;
			f_even[Np+n] = 0.055555555555555555;    //double(100*n)+2.f;
			f_odd[Np+n] = 0.055555555555555555;     //double(100*n)+3.f;
			f_even[2*Np+n] = 0.055555555555555555;  //double(100*n)+4.f;
			f_odd[2*Np+n] = 0.055555555555555555;   //double(100*n)+5.f;
			f_even[3*Np+n] = 0.055555555555555555;  //double(100*n)+6.f;
			f_odd[3*Np+n] = 0.0277777777777778;   //double(100*n)+7.f;
			f_even[4*Np+n] = 0.0277777777777778;   //double(100*n)+8.f;
			f_odd[4*Np+n] = 0.0277777777777778;   //double(100*n)+9.f;
			f_even[5*Np+n] = 0.0277777777777778;  //double(100*n)+10.f;
			f_odd[5*Np+n] = 0.0277777777777778;  //double(100*n)+11.f;
			f_even[6*Np+n] = 0.0277777777777778;  //double(100*n)+12.f;
			f_odd[6*Np+n] = 0.0277777777777778;  //double(100*n)+13.f;
			f_even[7*Np+n] = 0.0277777777777778;  //double(100*n)+14.f;
			f_odd[7*Np+n] = 0.0277777777777778;  //double(100*n)+15.f;
			f_even[8*Np+n] = 0.0277777777777778;  //double(100*n)+16.f;
			f_odd[8*Np+n] = 0.0277777777777778;  //double(100*n)+17.f;
			f_even[9*Np+n] = 0.0277777777777778;  //double(100*n)+18.f;
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_Init(double *dist, int Np)
{
	int n;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<Np ){
			dist[n] = 0.3333333333333333;
			dist[Np+n] = 0.055555555555555555;		//double(100*n)+1.f;
			dist[2*Np+n] = 0.055555555555555555;	//double(100*n)+2.f;
			dist[3*Np+n] = 0.055555555555555555;	//double(100*n)+3.f;
			dist[4*Np+n] = 0.055555555555555555;	//double(100*n)+4.f;
			dist[5*Np+n] = 0.055555555555555555;	//double(100*n)+5.f;
			dist[6*Np+n] = 0.055555555555555555;	//double(100*n)+6.f;
			dist[7*Np+n] = 0.0277777777777778;   //double(100*n)+7.f;
			dist[8*Np+n] = 0.0277777777777778;   //double(100*n)+8.f;
			dist[9*Np+n] = 0.0277777777777778;   //double(100*n)+9.f;
			dist[10*Np+n] = 0.0277777777777778;  //double(100*n)+10.f;
			dist[11*Np+n] = 0.0277777777777778;  //double(100*n)+11.f;
			dist[12*Np+n] = 0.0277777777777778;  //double(100*n)+12.f;
			dist[13*Np+n] = 0.0277777777777778;  //double(100*n)+13.f;
			dist[14*Np+n] = 0.0277777777777778;  //double(100*n)+14.f;
			dist[15*Np+n] = 0.0277777777777778;  //double(100*n)+15.f;
			dist[16*Np+n] = 0.0277777777777778;  //double(100*n)+16.f;
			dist[17*Np+n] = 0.0277777777777778;  //double(100*n)+17.f;
			dist[18*Np+n] = 0.0277777777777778;  //double(100*n)+18.f;
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

//__launch_bounds__(512,4)

__global__ void 
dvc_ScaLBL_AAodd_Compact(int *d_neighborList, double *dist, int Np) {

	int n;
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	int nread;
	int S = Np/NBLOCKS/NTHREADS+1;

	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<Np) {

			f0 = dist[n];

			nread = d_neighborList[n];
			f1 = dist[nread];

			nread = d_neighborList[n+2*Np];
			f3 = dist[nread];

			nread = d_neighborList[n+4*Np];
			f5 = dist[nread];

			nread = d_neighborList[n+6*Np];
			f7 = dist[nread];

			nread = d_neighborList[n+8*Np];
			f9 = dist[nread];

			nread = d_neighborList[n+10*Np];
			f11 = dist[nread];

			nread = d_neighborList[n+12*Np];
			f13 = dist[nread];

			nread = d_neighborList[n+14*Np];
			f15 = dist[nread];

			nread = d_neighborList[n+16*Np];
			f17 = dist[nread];


			nread = d_neighborList[n+Np];
			f2 = dist[nread];

			nread = d_neighborList[n+3*Np];
			f4 = dist[nread];

			nread = d_neighborList[n+5*Np];
			f6 = dist[nread];

			nread = d_neighborList[n+7*Np];
			f8 = dist[nread];

			nread = d_neighborList[n+9*Np];
			f10 = dist[nread];

			nread = d_neighborList[n+11*Np];
			f12 = dist[nread];

			nread = d_neighborList[n+13*Np];
			f14 = dist[nread];

			nread = d_neighborList[n+15*Np];
			f16 = dist[nread];

			nread = d_neighborList[n+17*Np];
			f18 = dist[nread];

			// ORIGINAL CORRECT WRITES
			//                              nwrite = d_neighborList[n];      naccess = 10*Np;
			//                              if (nwrite<0) { nwrite=n;        naccess = Np;  }
			//                              dist[nwrite + naccess]   = f1;

			//                              nwrite = d_neighborList[n+Np];   naccess = Np;
			//                              if (nwrite<0) { nwrite=n;        naccess = 10*Np; }
			//                              dist[nwrite + naccess]   = f2;

			nread = d_neighborList[n];
			dist[nread] = f2;

			nread = d_neighborList[n+2*Np];
			dist[nread] = f4;

			nread = d_neighborList[n+4*Np];
			dist[nread] = f6;

			nread = d_neighborList[n+6*Np];
			dist[nread] = f8;

			nread = d_neighborList[n+8*Np];
			dist[nread] = f10;

			nread = d_neighborList[n+10*Np];
			dist[nread] = f12;

			nread = d_neighborList[n+12*Np];
			dist[nread] = f14;

			nread = d_neighborList[n+14*Np];
			dist[nread] = f16;

			nread = d_neighborList[n+16*Np];
			dist[nread] = f18;


			nread = d_neighborList[n+Np];
			dist[nread] = f1;

			nread = d_neighborList[n+3*Np];
			dist[nread] = f3;

			nread = d_neighborList[n+5*Np];
			dist[nread] = f5;

			nread = d_neighborList[n+7*Np];
			dist[nread] = f7;

			nread = d_neighborList[n+9*Np];
			dist[nread] = f9;

			nread = d_neighborList[n+11*Np];
			dist[nread]= f11;

			nread = d_neighborList[n+13*Np];
			dist[nread] = f13;

			nread = d_neighborList[n+15*Np];
			dist[nread] = f15;

			nread = d_neighborList[n+17*Np];
			dist[nread] = f17;

		}
	}
}


__global__ void 
dvc_ScaLBL_AAodd_MRT(int *neighborList, double *dist, int start, int finish, int Np, double rlx_setA, double rlx_setB, double Fx, double Fy, double Fz) {

	int n;
	double fq;
	// conserved momemnts
	double rho,jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;

	int nread;
	int S = Np/NBLOCKS/NTHREADS+1;

	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			// q=0
			fq = dist[n];
			rho = fq;
			m1  = -30.0*fq;
			m2  = 12.0*fq;

			// q=1
			nread = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
			fq = dist[nread]; // reading the f1 data into register fq
			//fp = dist[10*Np+n];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jx = fq;
			m4 = -4.0*fq;
			m9 = 2.0*fq;
			m10 = -4.0*fq;

			// f2 = dist[10*Np+n];
			nread = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
			fq = dist[nread];  // reading the f2 data into register fq
			//fq = dist[Np+n];
			rho += fq;
			m1 -= 11.0*(fq);
			m2 -= 4.0*(fq);
			jx -= fq;
			m4 += 4.0*(fq);
			m9 += 2.0*(fq);
			m10 -= 4.0*(fq);

			// q=3
			nread = neighborList[n+2*Np]; // neighbor 4
			fq = dist[nread];
			//fq = dist[11*Np+n];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jy = fq;
			m6 = -4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 = fq;
			m12 = -2.0*fq;

			// q = 4
			nread = neighborList[n+3*Np]; // neighbor 3
			fq = dist[nread];
			//fq = dist[2*Np+n];
			rho+= fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jy -= fq;
			m6 += 4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 += fq;
			m12 -= 2.0*fq;

			// q=5
			nread = neighborList[n+4*Np];
			fq = dist[nread];
			//fq = dist[12*Np+n];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jz = fq;
			m8 = -4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 -= fq;
			m12 += 2.0*fq;


			// q = 6
			nread = neighborList[n+5*Np];
			fq = dist[nread];
			//fq = dist[3*Np+n];
			rho+= fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jz -= fq;
			m8 += 4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 -= fq;
			m12 += 2.0*fq;

			// q=7
			nread = neighborList[n+6*Np];
			fq = dist[nread];
			//fq = dist[13*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jy += fq;
			m6 += fq;
			m9  += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 = fq;
			m16 = fq;
			m17 = -fq;

			// q = 8
			nread = neighborList[n+7*Np];
			fq = dist[nread];
			//fq = dist[4*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jy -= fq;
			m6 -= fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 += fq;
			m16 -= fq;
			m17 += fq;

			// q=9
			nread = neighborList[n+8*Np];
			fq = dist[nread];
			//fq = dist[14*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jy -= fq;
			m6 -= fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 -= fq;
			m16 += fq;
			m17 += fq;

			// q = 10
			nread = neighborList[n+9*Np];
			fq = dist[nread];
			//fq = dist[5*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jy += fq;
			m6 += fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 -= fq;
			m16 -= fq;
			m17 -= fq;

			// q=11
			nread = neighborList[n+10*Np];
			fq = dist[nread];
			//fq = dist[15*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jz += fq;
			m8 += fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 = fq;
			m16 -= fq;
			m18 = fq;

			// q=12
			nread = neighborList[n+11*Np];
			fq = dist[nread];
			//fq = dist[6*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jz -= fq;
			m8 -= fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 += fq;
			m16 += fq;
			m18 -= fq;

			// q=13
			nread = neighborList[n+12*Np];
			fq = dist[nread];
			//fq = dist[16*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jz -= fq;
			m8 -= fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 -= fq;
			m16 -= fq;
			m18 -= fq;

			// q=14
			nread = neighborList[n+13*Np];
			fq = dist[nread];
			//fq = dist[7*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jz += fq;
			m8 += fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 -= fq;
			m16 += fq;
			m18 += fq;

			// q=15
			nread = neighborList[n+14*Np];
			fq = dist[nread];
			//fq = dist[17*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy += fq;
			m6 += fq;
			jz += fq;
			m8 += fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 = fq;
			m17 += fq;
			m18 -= fq;

			// q=16
			nread = neighborList[n+15*Np];
			fq = dist[nread];
			//fq = dist[8*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy -= fq;
			m6 -= fq;
			jz -= fq;
			m8 -= fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 += fq;
			m17 -= fq;
			m18 += fq;

			// q=17
			//fq = dist[18*Np+n];
			nread = neighborList[n+16*Np];
			fq = dist[nread];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy += fq;
			m6 += fq;
			jz -= fq;
			m8 -= fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 -= fq;
			m17 += fq;
			m18 += fq;

			// q=18
			nread = neighborList[n+17*Np];
			fq = dist[nread];
			//fq = dist[9*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy -= fq;
			m6 -= fq;
			jz += fq;
			m8 += fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 -= fq;
			m17 -= fq;
			m18 -= fq;

			//..............incorporate external force................................................
			//..............carry out relaxation process...............................................
			m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho - 11*rho) - m1);
			m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho) - m2);
			m4 = m4 + rlx_setB*((-0.6666666666666666*jx) - m4);
			m6 = m6 + rlx_setB*((-0.6666666666666666*jy) - m6);
			m8 = m8 + rlx_setB*((-0.6666666666666666*jz) - m8);
			m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho) - m9);
			m10 = m10 + rlx_setA*(-0.5*((2*jx*jx-jy*jy-jz*jz)/rho) - m10);
			m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho) - m11);
			m12 = m12 + rlx_setA*(-0.5*((jy*jy-jz*jz)/rho) - m12);
			m13 = m13 + rlx_setA*((jx*jy/rho) - m13);
			m14 = m14 + rlx_setA*((jy*jz/rho) - m14);
			m15 = m15 + rlx_setA*((jx*jz/rho) - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);
			//.......................................................................................................
			//.................inverse transformation......................................................

			// q=0
			fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
			dist[n] = fq;

			// q = 1
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10)+0.16666666*Fx;
			nread = neighborList[n+Np];
			dist[nread] = fq;

			// q=2
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10) -  0.16666666*Fx;
			nread = neighborList[n];
			dist[nread] = fq;

			// q = 3
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) + 0.16666666*Fy;
			nread = neighborList[n+3*Np];
			dist[nread] = fq;

			// q = 4
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) - 0.16666666*Fy;
			nread = neighborList[n+2*Np];
			dist[nread] = fq;

			// q = 5
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) + 0.16666666*Fz;
			nread = neighborList[n+5*Np];
			dist[nread] = fq;

			// q = 6
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) - 0.16666666*Fz;
			nread = neighborList[n+4*Np];
			dist[nread] = fq;

			// q = 7
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+mrt_V7*m9+mrt_V11*m10+
					mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17) + 0.08333333333*(Fx+Fy);
			
			nread = neighborList[n+7*Np];
			dist[nread] = fq;

			// q = 8
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
					+mrt_V12*m12+0.25*m13+0.125*(m17-m16) - 0.08333333333*(Fx+Fy);
			nread = neighborList[n+6*Np];
			dist[nread] = fq;

			// q = 9
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+mrt_V7*m9+mrt_V11*m10+
					mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17) + 0.08333333333*(Fx-Fy);
			nread = neighborList[n+9*Np];
			dist[nread] = fq;

			// q = 10
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+mrt_V7*m9+mrt_V11*m10+
					mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17)- 0.08333333333*(Fx-Fy);
			nread = neighborList[n+8*Np];
			dist[nread] = fq;

			// q = 11
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12+0.25*m15+0.125*(m18-m16) + 0.08333333333*(Fx+Fz);
			nread = neighborList[n+11*Np];
			dist[nread] = fq;

			// q = 12
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
					mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18) - 0.08333333333*(Fx+Fz);
			nread = neighborList[n+10*Np];
			dist[nread]= fq;

			// q = 13
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15-0.125*(m16+m18) + 0.08333333333*(Fx-Fz);
			nread = neighborList[n+13*Np];
			dist[nread] = fq;

			// q= 14
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15+0.125*(m16+m18) - 0.08333333333*(Fx-Fz);
			nread = neighborList[n+12*Np];
			dist[nread] = fq;


			// q = 15
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18) + 0.08333333333*(Fy+Fz);
			nread = neighborList[n+15*Np];
			dist[nread] = fq;

			// q = 16
			fq =  mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17)- 0.08333333333*(Fy+Fz);
			nread = neighborList[n+14*Np];
			dist[nread] = fq;


			// q = 17
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18) + 0.08333333333*(Fy-Fz);
			nread = neighborList[n+17*Np];
			dist[nread] = fq;

			// q = 18
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18) - 0.08333333333*(Fy-Fz);
			nread = neighborList[n+16*Np];
			dist[nread] = fq;

		}
	}
}


//__launch_bounds__(512,1)
__global__ void 
dvc_ScaLBL_AAeven_MRT(double *dist, int start, int finish, int Np, double rlx_setA, double rlx_setB, double Fx, double Fy, double Fz) {

	int n;
	double fq;
	// conserved momemnts
	double rho,jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){

			// q=0
			fq = dist[n];
			rho = fq;
			m1  = -30.0*fq;
			m2  = 12.0*fq;

			// q=1
			fq = dist[2*Np+n];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jx = fq;
			m4 = -4.0*fq;
			m9 = 2.0*fq;
			m10 = -4.0*fq;

			// q=2
			fq = dist[1*Np+n];
			rho += fq;
			m1 -= 11.0*(fq);
			m2 -= 4.0*(fq);
			jx -= fq;
			m4 += 4.0*(fq);
			m9 += 2.0*(fq);
			m10 -= 4.0*(fq);

			// q=3
			fq = dist[4*Np+n];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jy = fq;
			m6 = -4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 = fq;
			m12 = -2.0*fq;

			// q = 4
			fq = dist[3*Np+n];
			rho+= fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jy -= fq;
			m6 += 4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 += fq;
			m12 -= 2.0*fq;

			// q=5
			fq = dist[6*Np+n];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jz = fq;
			m8 = -4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 -= fq;
			m12 += 2.0*fq;

			// q = 6
			fq = dist[5*Np+n];
			rho+= fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jz -= fq;
			m8 += 4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 -= fq;
			m12 += 2.0*fq;

			// q=7
			fq = dist[8*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jy += fq;
			m6 += fq;
			m9  += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 = fq;
			m16 = fq;
			m17 = -fq;

			// q = 8
			fq = dist[7*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jy -= fq;
			m6 -= fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 += fq;
			m16 -= fq;
			m17 += fq;

			// q=9
			fq = dist[10*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jy -= fq;
			m6 -= fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 -= fq;
			m16 += fq;
			m17 += fq;

			// q = 10
			fq = dist[9*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jy += fq;
			m6 += fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 -= fq;
			m16 -= fq;
			m17 -= fq;

			// q=11
			fq = dist[12*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jz += fq;
			m8 += fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 = fq;
			m16 -= fq;
			m18 = fq;

			// q=12
			fq = dist[11*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jz -= fq;
			m8 -= fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 += fq;
			m16 += fq;
			m18 -= fq;

			// q=13
			fq = dist[14*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jz -= fq;
			m8 -= fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 -= fq;
			m16 -= fq;
			m18 -= fq;

			// q=14
			fq = dist[13*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jz += fq;
			m8 += fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 -= fq;
			m16 += fq;
			m18 += fq;

			// q=15
			fq = dist[16*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy += fq;
			m6 += fq;
			jz += fq;
			m8 += fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 = fq;
			m17 += fq;
			m18 -= fq;

			// q=16
			fq = dist[15*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy -= fq;
			m6 -= fq;
			jz -= fq;
			m8 -= fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 += fq;
			m17 -= fq;
			m18 += fq;

			// q=17
			fq = dist[18*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy += fq;
			m6 += fq;
			jz -= fq;
			m8 -= fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 -= fq;
			m17 += fq;
			m18 += fq;

			// q=18
			fq = dist[17*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy -= fq;
			m6 -= fq;
			jz += fq;
			m8 += fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 -= fq;
			m17 -= fq;
			m18 -= fq;

			//........................................................................
			//					READ THE DISTRIBUTIONS
			//		(read from opposite array due to previous swap operation)
			//........................................................................

			//..............incorporate external force................................................
			//..............carry out relaxation process...............................................
			m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho - 11*rho) - m1);
			m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho) - m2);
			m4 = m4 + rlx_setB*((-0.6666666666666666*jx) - m4);
			m6 = m6 + rlx_setB*((-0.6666666666666666*jy) - m6);
			m8 = m8 + rlx_setB*((-0.6666666666666666*jz) - m8);
			m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho) - m9);
			m10 = m10 + rlx_setA*(-0.5*((2*jx*jx-jy*jy-jz*jz)/rho) - m10);
			m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho) - m11);
			m12 = m12 + rlx_setA*(-0.5*((jy*jy-jz*jz)/rho) - m12);
			m13 = m13 + rlx_setA*((jx*jy/rho) - m13);
			m14 = m14 + rlx_setA*((jy*jz/rho) - m14);
			m15 = m15 + rlx_setA*((jx*jz/rho) - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);
			//.......................................................................................................
			//.................inverse transformation......................................................

			// q=0
			fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
			dist[n] = fq;

			// q = 1
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10) + 0.16666666*Fx;
			dist[1*Np+n] = fq;

			// q=2
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10) -  0.16666666*Fx;
			dist[2*Np+n] = fq;

			// q = 3
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) + 0.16666666*Fy;
			dist[3*Np+n] = fq;

			// q = 4
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) - 0.16666666*Fy;
			dist[4*Np+n] = fq;

			// q = 5
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) + 0.16666666*Fz;
			dist[5*Np+n] = fq;

			// q = 6
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) - 0.16666666*Fz;
			dist[6*Np+n] = fq;

			// q = 7
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17) + 
					0.08333333333*(Fx+Fy);
			dist[7*Np+n] = fq;


			// q = 8
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
					+mrt_V12*m12+0.25*m13+0.125*(m17-m16) - 0.08333333333*(Fx+Fy);
			dist[8*Np+n] = fq;

			// q = 9
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17)+
					0.08333333333*(Fx-Fy);
			dist[9*Np+n] = fq;

			// q = 10
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17)-
					0.08333333333*(Fx-Fy);
			dist[10*Np+n] = fq;


			// q = 11
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12+0.25*m15+0.125*(m18-m16) + 0.08333333333*(Fx+Fz);
			dist[11*Np+n] = fq;

			// q = 12
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
					mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18)-
					0.08333333333*(Fx+Fz);
			dist[12*Np+n] = fq;

			// q = 13
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15-0.125*(m16+m18) + 0.08333333333*(Fx-Fz);
			dist[13*Np+n] = fq;

			// q= 14
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15+0.125*(m16+m18) - 0.08333333333*(Fx-Fz);

			dist[14*Np+n] = fq;

			// q = 15
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18) + 0.08333333333*(Fy+Fz);
			dist[15*Np+n] = fq;

			// q = 16
			fq =  mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17)- 0.08333333333*(Fy+Fz);
			dist[16*Np+n] = fq;


			// q = 17
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18) + 0.08333333333*(Fy-Fz);
			dist[17*Np+n] = fq;

			// q = 18
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18) - 0.08333333333*(Fy-Fz);
			dist[18*Np+n] = fq;
			//........................................................................
		}
	}
}

//__launch_bounds__(512,4)

__global__ void dvc_ScaLBL_AAeven_Compact( double *dist, int Np) {

	int n;
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;

		if ( n<Np ){

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


__global__  void dvc_ScaLBL_D3Q19_Momentum(double *dist, double *vel, int N)
{
	int n;
	// distributions
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double vx,vy,vz;

	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N){

			f2 = dist[2*N+n];
			f4 = dist[4*N+n];
			f6 = dist[6*N+n];
			f8 = dist[8*N+n];
			f10 = dist[10*N+n];
			f12 = dist[12*N+n];
			f14 = dist[14*N+n];
			f16 = dist[16*N+n];
			f18 = dist[18*N+n];
			//........................................................................
			f1 = dist[N+n];
			f3 = dist[3*N+n];
			f5 = dist[5*N+n];
			f7 = dist[7*N+n];
			f9 = dist[9*N+n];
			f11 = dist[11*N+n];
			f13 = dist[13*N+n];
			f15 = dist[15*N+n];
			f17 = dist[17*N+n];			

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

__global__  void dvc_ScaLBL_D3Q19_Pressure(const double *dist, double *Pressure, int N)
{
	int n;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N){				//.......................................................................
			// Registers to store the distributions
			//........................................................................
			//........................................................................
			// Registers to store the distributions
			//........................................................................
			f0 = dist[n];
			f2 = dist[2*N+n];
			f4 = dist[4*N+n];
			f6 = dist[6*N+n];
			f8 = dist[8*N+n];
			f10 = dist[10*N+n];
			f12 = dist[12*N+n];
			f14 = dist[14*N+n];
			f16 = dist[16*N+n];
			f18 = dist[18*N+n];
			//........................................................................
			f1 = dist[N+n];
			f3 = dist[3*N+n];
			f5 = dist[5*N+n];
			f7 = dist[7*N+n];
			f9 = dist[9*N+n];
			f11 = dist[11*N+n];
			f13 = dist[13*N+n];
			f15 = dist[15*N+n];
			f17 = dist[17*N+n];
			//.................Compute the velocity...................................
			Pressure[n] = 0.3333333333333333*(f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+
					f9+f12+f11+f14+f13+f16+f15+f18+f17);
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAeven_Pressure_BC_z(int *list, double *dist, double din, int count, int Np)
{
	int idx, n;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz,Cyz,Cxz;
	ux = uy = 0.0;

	idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < count){

		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
		f7 = dist[8*Np+n];
		f8 = dist[7*Np+n];
		f9 = dist[10*Np+n];
		f10 = dist[9*Np+n];
		f12 = dist[11*Np+n];
		f13 = dist[14*Np+n];
		f16 = dist[15*Np+n];
		f17 = dist[18*Np+n];
		//...................................................
		// Determine the inlet flow velocity
		//ux = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14);
		//uy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18);
		uz = din - (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f5 = f6 + 0.33333333333333338*uz;
		f11 = f12 + 0.16666666666666678*(uz+ux)-Cxz;
		f14 = f13 + 0.16666666666666678*(uz-ux)+Cxz;
		f15 = f16 + 0.16666666666666678*(uy+uz)-Cyz;
		f18 = f17 + 0.16666666666666678*(uz-uy)+Cyz;
		//........Store in "opposite" memory location..........
		dist[6*Np+n] = f5;
		dist[12*Np+n] = f11;
		dist[13*Np+n] = f14;
		dist[16*Np+n] = f15;
		dist[17*Np+n] = f18;
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAeven_Pressure_BC_Z(int *list, double *dist, double dout, int count, int Np)
{
	int idx,n;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz,Cyz,Cxz;
	ux = uy = 0.0;

	idx = blockIdx.x*blockDim.x + threadIdx.x;

	// Loop over the boundary - threadblocks delineated by start...finish
	if ( idx < count ){

		n = list[idx];
		//........................................................................
		// Read distributions 
		//........................................................................
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
		f7 = dist[8*Np+n];
		f8 = dist[7*Np+n];
		f9 = dist[10*Np+n];
		f10 = dist[9*Np+n];
		f11 = dist[12*Np+n];
		f14 = dist[13*Np+n];
		f15 = dist[16*Np+n];
		f18 = dist[17*Np+n];
		
		// Determine the outlet flow velocity
		//ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
		//uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		uz = -dout + (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f6 = f5 - 0.33333333333333338*uz;
		f12 = f11 - 0.16666666666666678*(uz+ux)+Cxz;
		f13 = f14 - 0.16666666666666678*(uz-ux)-Cxz;
		f16 = f15 - 0.16666666666666678*(uy+uz)+Cyz;
		f17 = f18 - 0.16666666666666678*(uz-uy)-Cyz;

		dist[5*Np+n] = f6;
		dist[11*Np+n] = f12;
		dist[14*Np+n] = f13;
		dist[15*Np+n] = f16;
		dist[18*Np+n] = f17;
		//...................................................
	}
}
__global__  void dvc_ScaLBL_D3Q19_Reflection_BC_z(int *list, double *dist, int count, int Np){
	int idx, n;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		double f5 = 0.111111111111111111111111 - dist[6*Np+n];
		double f11 = 0.05555555555555555555556 - dist[12*Np+n];
		double f14 = 0.05555555555555555555556 - dist[13*Np+n];
		double f15 = 0.05555555555555555555556 - dist[16*Np+n];
		double f18 = 0.05555555555555555555556 - dist[17*Np+n];
		
		dist[6*Np+n] = f5;
		dist[12*Np+n] = f11;
		dist[13*Np+n] = f14;
		dist[16*Np+n] = f15;
		dist[17*Np+n] = f18;
	}
}

__global__  void dvc_ScaLBL_D3Q19_Reflection_BC_Z(int *list, double *dist, int count, int Np){
	int idx, n;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		double f6 = 0.111111111111111111111111 - dist[5*Np+n];
		double f12 = 0.05555555555555555555556 - dist[11*Np+n];
		double f13 = 0.05555555555555555555556 - dist[14*Np+n] ;
		double f16 = 0.05555555555555555555556 - dist[15*Np+n];
		double f17 = 0.05555555555555555555556 - dist[18*Np+n];
		
		dist[5*Np+n] = f6;
		dist[11*Np+n] = f12;
		dist[14*Np+n] = f13;
		dist[15*Np+n] = f16;
		dist[18*Np+n] = f17;
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAodd_Pressure_BC_z(int *d_neighborList, int *list, double *dist, double din, int count, int Np)
{
	int idx, n;
	int nread;
	int nr5,nr11,nr14,nr15,nr18;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz,Cyz,Cxz;
	ux = uy = 0.0;

	idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < count){
		
		n = list[idx];
		f0 = dist[n];
				
		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+6*Np];
		f7 = dist[nread];

		nread = d_neighborList[n+8*Np];
		f9 = dist[nread];

		nread = d_neighborList[n+12*Np];
		f13 = dist[nread];

		nread = d_neighborList[n+16*Np];
		f17 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

		nread = d_neighborList[n+7*Np];
		f8 = dist[nread];

		nread = d_neighborList[n+9*Np];
		f10 = dist[nread];

		nread = d_neighborList[n+11*Np];
		f12 = dist[nread];

		nread = d_neighborList[n+15*Np];
		f16 = dist[nread];

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		nr11 = d_neighborList[n+10*Np];
		nr15 = d_neighborList[n+14*Np];
		nr14 = d_neighborList[n+13*Np];
		nr18 = d_neighborList[n+17*Np];
		
		//...................................................
		//........Determine the inlet flow velocity.........
		//ux = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14);
		//uy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18);
		uz = din - (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f5 = f6 + 0.33333333333333338*uz;
		f11 = f12 + 0.16666666666666678*(uz+ux)-Cxz;
		f14 = f13 + 0.16666666666666678*(uz-ux)+Cxz;
		f15 = f16 + 0.16666666666666678*(uy+uz)-Cyz;
		f18 = f17 + 0.16666666666666678*(uz-uy)+Cyz;
		//........Store in "opposite" memory location..........
		dist[nr5] = f5;
		dist[nr11] = f11;
		dist[nr14] = f14;
		dist[nr15] = f15;
		dist[nr18] = f18;
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAodd_Pressure_BC_Z(int *d_neighborList, int *list, double *dist, double dout, int count, int Np)
{
	int idx,n,nread;
	int nr6,nr12,nr13,nr16,nr17;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz,Cyz,Cxz;
	ux = uy = 0.0;

	idx = blockIdx.x*blockDim.x + threadIdx.x;

	// Loop over the boundary - threadblocks delineated by start...finish
	if ( idx < count ){

		n = list[idx];
		//........................................................................
		// Read distributions 
		//........................................................................
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		f5 = dist[nread];

		nread = d_neighborList[n+6*Np];
		f7 = dist[nread];

		nread = d_neighborList[n+8*Np];
		f9 = dist[nread];

		nread = d_neighborList[n+10*Np];
		f11 = dist[nread];

		nread = d_neighborList[n+14*Np];
		f15 = dist[nread];


		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+7*Np];
		f8 = dist[nread];

		nread = d_neighborList[n+9*Np];
		f10 = dist[nread];

		nread = d_neighborList[n+13*Np];
		f14 = dist[nread];

		nread = d_neighborList[n+17*Np];
		f18 = dist[nread];
		
		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		nr12 = d_neighborList[n+11*Np];
		nr16 = d_neighborList[n+15*Np];
		nr17 = d_neighborList[n+16*Np];
		nr13 = d_neighborList[n+12*Np];

		
		//........Determine the outlet flow velocity.........
		//ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
		//uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		uz = -dout + (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f6 = f5 - 0.33333333333333338*uz;
		f12 = f11 - 0.16666666666666678*(uz+ux)+Cxz;
		f13 = f14 - 0.16666666666666678*(uz-ux)-Cxz;
		f16 = f15 - 0.16666666666666678*(uy+uz)+Cyz;
		f17 = f18 - 0.16666666666666678*(uz-uy)-Cyz;

		//........Store in "opposite" memory location..........
		dist[nr6] = f6;
		dist[nr12] = f12;
		dist[nr13] = f13;
		dist[nr16] = f16;
		dist[nr17] = f17;
		//...................................................
	}
}


__global__  void dvc_ScaLBL_D3Q19_AAeven_Flux_BC_z(int *list, double *dist, double flux, double Area, 
		double *dvcsum, int count, int Np)
{
	int idx, n;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double factor = 1.f/(Area);
	double sum = 0.f;

	idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < count){
		
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
		f7 = dist[8*Np+n];
		f8 = dist[7*Np+n];
		f9 = dist[10*Np+n];
		f10 = dist[9*Np+n];
		f12 = dist[11*Np+n];
		f13 = dist[14*Np+n];
		f16 = dist[15*Np+n];
		f17 = dist[18*Np+n];
		sum = factor*(f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));
	}

	//sum = blockReduceSum(sum);
	//if (threadIdx.x==0)
	//   atomicAdd(dvcsum, sum);
	
    extern __shared__ double temp[];
    thread_group g = this_thread_block();
    double block_sum = reduce_sum(g, temp, sum);

    if (g.thread_rank() == 0) atomicAdd(dvcsum, block_sum);
}


__global__  void dvc_ScaLBL_D3Q19_AAodd_Flux_BC_z(int *d_neighborList, int *list, double *dist, double flux, 
		double Area, double *dvcsum, int count, int Np)
{
	int idx, n;
	int nread;

	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double factor = 1.f/(Area);
	double sum = 0.f;

	idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < count){
		
		n = list[idx];
				
		f0 = dist[n];
		
		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+6*Np];
		f7 = dist[nread];

		nread = d_neighborList[n+8*Np];
		f9 = dist[nread];

		nread = d_neighborList[n+12*Np];
		f13 = dist[nread];

		nread = d_neighborList[n+16*Np];
		f17 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

		nread = d_neighborList[n+7*Np];
		f8 = dist[nread];

		nread = d_neighborList[n+9*Np];
		f10 = dist[nread];

		nread = d_neighborList[n+11*Np];
		f12 = dist[nread];

		nread = d_neighborList[n+15*Np];
		f16 = dist[nread];

		sum = factor*(f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));

	}

	//sum = blockReduceSum(sum);
	//if (threadIdx.x==0)
	//   atomicAdd(dvcsum, sum);
	
    extern __shared__ double temp[];
    thread_group g = this_thread_block();
    double block_sum = reduce_sum(g, temp, sum);

    if (g.thread_rank() == 0) atomicAdd(dvcsum, block_sum);
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

extern "C" void ScaLBL_D3Q19_AA_Init(double *f_even, double *f_odd, int Np){
	dvc_ScaLBL_D3Q19_AA_Init<<<NBLOCKS,NTHREADS >>>(f_even, f_odd, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AA_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_Init(double *dist, int Np){
	dvc_ScaLBL_D3Q19_Init<<<NBLOCKS,NTHREADS >>>(dist, Np);
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

extern "C" void ScaLBL_D3Q19_AAeven_Compact( double *d_dist,  int Np) {
        cudaFuncSetCacheConfig(dvc_ScaLBL_AAeven_Compact, cudaFuncCachePreferL1);
	dvc_ScaLBL_AAeven_Compact<<<NBLOCKS,NTHREADS>>>(d_dist, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Compact( int *d_neighborList, double *d_dist, int Np) {
        cudaFuncSetCacheConfig(dvc_ScaLBL_AAodd_Compact, cudaFuncCachePreferL1);
	dvc_ScaLBL_AAodd_Compact<<<NBLOCKS,NTHREADS>>>(d_neighborList, d_dist,Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_Momentum(double *dist, double *vel, int Np){

	dvc_ScaLBL_D3Q19_Momentum<<<NBLOCKS,NTHREADS >>>(dist, vel, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Velocity: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_Pressure(double *fq, double *Pressure, int Np){
	dvc_ScaLBL_D3Q19_Pressure<<< NBLOCKS,NTHREADS >>>(fq, Pressure, Np);
}

extern "C" void ScaLBL_D3Q19_Velocity_BC_z(double *disteven, double *distodd, double uz,int Nx, int Ny, int Nz){
	int GRID = Nx*Ny / 512 + 1;
	dvc_D3Q19_Velocity_BC_z<<<GRID,512>>>(disteven,distodd, uz, Nx, Ny, Nz);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Velocity_BC_z: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_Velocity_BC_Z(double *disteven, double *distodd, double uz, int Nx, int Ny, int Nz, int outlet){
	int GRID = Nx*Ny / 512 + 1;
	dvc_D3Q19_Velocity_BC_Z<<<GRID,512>>>(disteven, distodd, uz, Nx, Ny, Nz, outlet);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Velocity_BC_Z: %s \n",cudaGetErrorString(err));
	}
}

extern "C" double ScaLBL_D3Q19_Flux_BC_z(double *disteven, double *distodd, double flux,int Nx, int Ny, int Nz){

	int GRID = Nx*Ny / 512 + 1;

	// IMPORTANT -- this routine may fail if Nx*Ny > 512*512
	if (Nx*Ny > 512*512){
		printf("WARNING (ScaLBL_D3Q19_Flux_BC_z): CUDA reduction operation may fail if Nx*Ny > 512*512");
	}

	// Allocate memory to store the sums
	double din;
	double sum[1];
 	double *dvcsum;
	int sharedBytes = NTHREADS*sizeof(double);
	cudaMalloc((void **)&dvcsum,sizeof(double)*Nx*Ny);
	cudaMemset(dvcsum,0,sizeof(double)*Nx*Ny);
	
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Flux_BC_z (memory allocation): %s \n",cudaGetErrorString(err));
	}

	// compute the local flux and store the result
	dvc_D3Q19_Flux_BC_z<<<GRID,512,sharedBytes>>>(disteven, distodd, flux, dvcsum, Nx, Ny, Nz);
	
	err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Flux_BC_z (flux calculation, step 1): %s \n",cudaGetErrorString(err));
	}

	// Now read the total flux
	cudaMemcpy(&sum[0],dvcsum,sizeof(double),cudaMemcpyDeviceToHost);
	din=sum[0];
	
	err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Flux_BC_z (flux calculation, step 2): %s \n",cudaGetErrorString(err));
	}

	// free the memory needed for reduction
	cudaFree(dvcsum);

	return din;
}


extern "C" void ScaLBL_D3Q19_AAeven_Pressure_BC_z(int *list, double *dist, double din, int count, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_AAeven_Pressure_BC_z<<<GRID,512>>>(list, dist, din, count, N);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_Pressure_BC_z (kernel): %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_Pressure_BC_Z(int *list, double *dist, double dout, int count, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_AAeven_Pressure_BC_Z<<<GRID,512>>>(list, dist, dout, count, N);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_Pressure_BC_Z (kernel): %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Pressure_BC_z(int *neighborList, int *list, double *dist, double din, int count, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_AAodd_Pressure_BC_z<<<GRID,512>>>(neighborList, list, dist, din, count, N);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_Pressure_BC_z (kernel): %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Pressure_BC_Z(int *neighborList, int *list, double *dist, double dout, int count, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_AAodd_Pressure_BC_Z<<<GRID,512>>>(neighborList, list, dist, dout, count, N);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_Pressure_BC_Z (kernel): %s \n",cudaGetErrorString(err));
	}
}


extern "C" double ScaLBL_D3Q19_AAeven_Flux_BC_z(int *list, double *dist, double flux, double area, 
		 int count, int N){

	int GRID = count / 512 + 1;

	// IMPORTANT -- this routine may fail if Nx*Ny > 512*512
	if (count > 512*512){
		printf("WARNING (ScaLBL_D3Q19_Flux_BC_Z): CUDA reduction operation may fail if count > 512*512");
	}

	// Allocate memory to store the sums
	double din;
	double sum[1];
 	double *dvcsum;
	cudaMalloc((void **)&dvcsum,sizeof(double)*count);
	cudaMemset(dvcsum,0,sizeof(double)*count);
	int sharedBytes = 512*sizeof(double);
	
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_Flux_BC_z (memory allocation): %s \n",cudaGetErrorString(err));
	}

	// compute the local flux and store the result
	dvc_ScaLBL_D3Q19_AAeven_Flux_BC_z<<<GRID,512,sharedBytes>>>(list, dist, flux, area, dvcsum, count, N);
	err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_Flux_BC_z (kernel): %s \n",cudaGetErrorString(err));
	}

	// Now read the total flux
	cudaMemcpy(&sum[0],dvcsum,sizeof(double),cudaMemcpyDeviceToHost);
	din=sum[0];
	err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_Flux_BC_z (reduction): %s \n",cudaGetErrorString(err));
	}

	// free the memory needed for reduction
	cudaFree(dvcsum);

	return din;
}

extern "C" double ScaLBL_D3Q19_AAodd_Flux_BC_z(int *neighborList, int *list, double *dist, double flux, 
		double area, int count, int N){

	int GRID = count / 512 + 1;

	// IMPORTANT -- this routine may fail if Nx*Ny > 512*512
	if (count > 512*512){
		printf("WARNING (ScaLBL_D3Q19_AAodd_Flux_BC_z): CUDA reduction operation may fail if count > 512*512");
	}

	// Allocate memory to store the sums
	double din;
	double sum[1];
 	double *dvcsum;
	cudaMalloc((void **)&dvcsum,sizeof(double)*count);
	cudaMemset(dvcsum,0,sizeof(double)*count);
	int sharedBytes = 512*sizeof(double);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_Flux_BC_z (memory allocation): %s \n",cudaGetErrorString(err));
	}

	// compute the local flux and store the result
	dvc_ScaLBL_D3Q19_AAodd_Flux_BC_z<<<GRID,512,sharedBytes>>>(neighborList, list, dist, flux, area, dvcsum, count, N);
	err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_Flux_BC_z (kernel): %s \n",cudaGetErrorString(err));
	}
	// Now read the total flux
	cudaMemcpy(&sum[0],dvcsum,sizeof(double),cudaMemcpyDeviceToHost);
	din=sum[0];
	err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_Flux_BC_z (reduction): %s \n",cudaGetErrorString(err));
	}

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

extern "C" void ScaLBL_D3Q19_Reflection_BC_z(int *list, double *dist, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_Reflection_BC_z<<<GRID,512>>>(list, dist, count, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Reflection_BC_z (kernel): %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_Reflection_BC_Z(int *list, double *dist, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_Reflection_BC_Z<<<GRID,512>>>(list, dist, count, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Reflection_BC_Z (kernel): %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_MRT(double *dist, int start, int finish, int Np, double rlx_setA, double rlx_setB, double Fx,
       double Fy, double Fz){
       
       dvc_ScaLBL_AAeven_MRT<<<NBLOCKS,NTHREADS >>>(dist,start,finish,Np,rlx_setA,rlx_setB,Fx,Fy,Fz);

       cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_MRT: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_MRT(int *neighborlist, double *dist, int start, int finish, int Np, double rlx_setA, double rlx_setB, double Fx,
       double Fy, double Fz){
       
       dvc_ScaLBL_AAodd_MRT<<<NBLOCKS,NTHREADS >>>(neighborlist,dist,start,finish,Np,rlx_setA,rlx_setB,Fx,Fy,Fz);

       cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_MRT: %s \n",cudaGetErrorString(err));
	}
}

