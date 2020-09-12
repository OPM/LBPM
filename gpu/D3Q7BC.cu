#include <math.h>
#include <stdio.h>
#include <cuda_profiler_api.h>

#define NBLOCKS 1024
#define NTHREADS 256


__global__ void dvc_ScaLBL_Solid_Dirichlet_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count)
{

    int idx;
    int iq,ib;
    double value_b,value_q;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		iq = BounceBackDist_list[idx];
        ib = BounceBackSolid_list[idx];
		value_b = BoundaryValue[ib];//get boundary value from a solid site
        value_q = dist[iq];
		dist[iq] = -1.0*value_q + value_b*0.25;//NOTE 0.25 is the speed of sound for D3Q7 lattice
	}
}

__global__ void dvc_ScaLBL_Solid_Neumann_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count)
{

    int idx;
    int iq,ib;
    double value_b,value_q;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		iq = BounceBackDist_list[idx];
        ib = BounceBackSolid_list[idx];
		value_b = BoundaryValue[ib];//get boundary value from a solid site
        value_q = dist[iq];
		dist[iq] = value_q + value_b;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z(int *list, double *dist, double Vin, int count, int Np)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
		//...................................................
		f5 = Vin - (f0+f1+f2+f3+f4+f6);
		dist[6*Np+n] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z(int *list, double *dist, double Vout, int count, int Np)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
		//...................................................
		f6 = Vout - (f0+f1+f2+f3+f4+f5);
		dist[5*Np+n] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z(int *d_neighborList, int *list, double *dist, double Vin, int count, int Np)
{
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		f5 = Vin - (f0+f1+f2+f3+f4+f6);
		dist[nr5] = f5;
	}
}

__global__ void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np)
{
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		f5 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		f6 = Vout - (f0+f1+f2+f3+f4+f5);
		dist[nr6] = f6;
	}
}

__global__ void dvc_ScaLBL_Poisson_D3Q7_BC_z(int *list, int *Map, double *Psi, double Vin, int count)
{
	int idx,n,nm;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		nm = Map[n];
		Psi[nm] = Vin;
	}
}


__global__ void dvc_ScaLBL_Poisson_D3Q7_BC_Z(int *list, int *Map, double *Psi, double Vout, int count)
{
	int idx,n,nm;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		nm = Map[n];
		Psi[nm] = Vout;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z(int *list, double *dist, double Cin, int count, int Np)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
		//...................................................
		f5 = Cin - (f0+f1+f2+f3+f4+f6);
		dist[6*Np+n] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z(int *list, double *dist, double Cout, int count, int Np)
{
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
		//...................................................
		f6 = Cout - (f0+f1+f2+f3+f4+f5);
		dist[5*Np+n] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z(int *d_neighborList, int *list, double *dist, double Cin, int count, int Np)
{
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		f5 = Cin - (f0+f1+f2+f3+f4+f6);
		dist[nr5] = f5;
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, int count, int Np)
{
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		f5 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		f6 = Cout - (f0+f1+f2+f3+f4+f5);
		dist[nr6] = f6;
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

extern "C" void ScaLBL_D3Q19_GreyIMRT_Init(double *dist, int Np, double Den){
	dvc_ScaLBL_D3Q19_GreyIMRT_Init<<<NBLOCKS,NTHREADS >>>(dist, Np, Den);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_GreyIMRT_Init: %s \n",cudaGetErrorString(err));
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

extern "C" void ScaLBL_D3Q19_AAeven_Compact(char * ID, double *d_dist,  int Np) {
        cudaFuncSetCacheConfig(dvc_ScaLBL_AAeven_Compact, cudaFuncCachePreferL1);
	dvc_ScaLBL_AAeven_Compact<<<NBLOCKS,NTHREADS>>>(ID, d_dist, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Compact(char * ID, int *d_neighborList, double *d_dist, int Np) {
        cudaFuncSetCacheConfig(dvc_ScaLBL_AAodd_Compact, cudaFuncCachePreferL1);
	dvc_ScaLBL_AAodd_Compact<<<NBLOCKS,NTHREADS>>>(ID,d_neighborList, d_dist,Np);
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
//*******************************************************************************
//*******************************************************************************
//*******************************************************************************


//*******************************************************************************
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

