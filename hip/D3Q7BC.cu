#include <math.h>
#include <stdio.h>
#include "hip/hip_runtime.h"

#define NBLOCKS 560
#define NTHREADS 128

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

__global__ void dvc_ScaLBL_Solid_DirichletAndNeumann_D3Q7(double *dist, double *BoundaryValue,int *BoundaryLabel, int *BounceBackDist_list, int *BounceBackSolid_list, int count)
{

    int idx;
    int iq,ib;
    double value_b,value_b_label,value_q;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		iq = BounceBackDist_list[idx];
        ib = BounceBackSolid_list[idx];
		value_b = BoundaryValue[ib];//get boundary value from a solid site
		value_b_label = BoundaryLabel[ib];//get boundary label (i.e. type of BC) from a solid site
        value_q = dist[iq];
        if (value_b_label==1){//Dirichlet BC
		    dist[iq] = -1.0*value_q + value_b*0.25;//NOTE 0.25 is the speed of sound for D3Q7 lattice
        }
        if (value_b_label==2){//Neumann BC
		    dist[iq] = value_q + value_b;
        }
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

__global__ void dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np)
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
	if (idx < count){
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

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, int count, int Np)
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

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
        fsum_partial = f0+f1+f2+f3+f4+f6;
        uz = VelocityZ[n];
		//...................................................
        f5 =(FluxIn+(1.0-0.5/tau)*f6-0.5*uz*fsum_partial/tau)/(1.0-0.5/tau+0.5*uz/tau); 
		dist[6*Np+n] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
    int idx,n;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < count){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
        fsum_partial = f0+f1+f2+f3+f4+f5;
        uz = VelocityZ[n];
		//...................................................
        f6 =(FluxIn+(1.0-0.5/tau)*f5+0.5*uz*fsum_partial/tau)/(1.0-0.5/tau-0.5*uz/tau); 
		dist[5*Np+n] = f6;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr5;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
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

        fsum_partial = f0+f1+f2+f3+f4+f6;
        uz = VelocityZ[n];
		//...................................................
        f5 =(FluxIn+(1.0-0.5/tau)*f6-0.5*uz*fsum_partial/tau)/(1.0-0.5/tau+0.5*uz/tau); 

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		dist[nr5] = f5;
	}
}

__global__ void dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np)
{
    //NOTE: FluxIn is the inward flux
	int idx, n;
    int nread,nr6;
	double f0,f1,f2,f3,f4,f5,f6;
    double fsum_partial;
    double uz;
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

        fsum_partial = f0+f1+f2+f3+f4+f5;
        uz = VelocityZ[n];
		//...................................................
        f6 =(FluxIn+(1.0-0.5/tau)*f5+0.5*uz*fsum_partial/tau)/(1.0-0.5/tau-0.5*uz/tau); 

		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		dist[nr6] = f6;
	}
}
//*************************************************************************

extern "C" void ScaLBL_Solid_Dirichlet_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_Solid_Dirichlet_D3Q7<<<GRID,512>>>(dist, BoundaryValue, BounceBackDist_list, BounceBackSolid_list, count);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_Solid_Dirichlet_D3Q7 (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_Solid_Neumann_D3Q7(double *dist, double *BoundaryValue, int *BounceBackDist_list, int *BounceBackSolid_list, int count){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_Solid_Neumann_D3Q7<<<GRID,512>>>(dist, BoundaryValue, BounceBackDist_list, BounceBackSolid_list, count);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_Solid_Neumann_D3Q7 (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_Solid_DirichletAndNeumann_D3Q7(double *dist, double *BoundaryValue,int *BoundaryLabel, int *BounceBackDist_list, int *BounceBackSolid_list, int count){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_Solid_DirichletAndNeumann_D3Q7<<<GRID,512>>>(dist, BoundaryValue, BoundaryLabel, BounceBackDist_list, BounceBackSolid_list, count);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("hip error in ScaLBL_Solid_DirichletAndNeumann_D3Q7 (kernel): %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z(int *list, double *dist, double Vin, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z<<<GRID,512>>>(list, dist, Vin, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z(int *list, double *dist, double Vout, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z<<<GRID,512>>>(list, dist, Vout, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z(int *d_neighborList, int *list, double *dist, double Vin, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z<<<GRID,512>>>(d_neighborList, list, dist, Vin, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z<<<GRID,512>>>(d_neighborList, list, dist, Vout, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_Poisson_D3Q7_BC_z(int *list, int *Map, double *Psi, double Vin, int count){
	int GRID = count / 512 + 1;
    dvc_ScaLBL_Poisson_D3Q7_BC_z<<<GRID,512>>>(list, Map, Psi, Vin, count);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_Poisson_D3Q7_BC_z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_Poisson_D3Q7_BC_Z(int *list, int *Map, double *Psi, double Vout, int count){
	int GRID = count / 512 + 1;
    dvc_ScaLBL_Poisson_D3Q7_BC_Z<<<GRID,512>>>(list, Map, Psi, Vout, count);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_Poisson_D3Q7_BC_Z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z(int *list, double *dist, double Cin, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z<<<GRID,512>>>(list, dist, Cin, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z(int *list, double *dist, double Cout, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z<<<GRID,512>>>(list, dist, Cout, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z(int *d_neighborList, int *list, double *dist, double Cin, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z<<<GRID,512>>>(d_neighborList, list, dist, Cin, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z<<<GRID,512>>>(d_neighborList, list, dist, Cout, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_BC_z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_BC_z<<<GRID,512>>>(list, dist, FluxIn, tau, VelocityZ, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAeven_Ion_Flux_BC_z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_BC_Z(int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAeven_Ion_Flux_BC_Z<<<GRID,512>>>(list, dist, FluxIn, tau, VelocityZ, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAeven_Ion_Flux_BC_Z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_BC_z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_BC_z<<<GRID,512>>>(d_neighborList, list, dist, FluxIn, tau, VelocityZ, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAodd_Ion_Flux_BC_z (kernel): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_BC_Z(int *d_neighborList, int *list, double *dist, double FluxIn, double tau, double *VelocityZ, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_AAodd_Ion_Flux_BC_Z<<<GRID,512>>>(d_neighborList, list, dist, FluxIn, tau, VelocityZ, count, Np);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q7_AAodd_Ion_Flux_BC_Z (kernel): %s \n",hipGetErrorString(err));
	}
}
