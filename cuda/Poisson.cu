#include <stdio.h>
#include <math.h>
//#include <cuda_profiler_api.h>

#define NBLOCKS 1024
#define NTHREADS 256

__global__  void dvc_ScaLBL_D3Q7_AAodd_Poisson_ElectricPotential(int *neighborList,int *Map, double *dist, double *Psi, int start, int finish, int Np){
	int n;
	double psi;//electric potential
	double fq;
	int nread;
    int idx;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {

            // q=0
            fq = dist[n];
            psi = fq;

            // q=1
            nread = neighborList[n]; 
            fq = dist[nread];
            psi += fq;
            
            // q=2
            nread = neighborList[n+Np]; 
            fq = dist[nread];  
            psi += fq;

            // q=3
            nread = neighborList[n+2*Np]; 
            fq = dist[nread];
            psi += fq;

            // q = 4
            nread = neighborList[n+3*Np]; 
            fq = dist[nread];
            psi += fq;

            // q=5
            nread = neighborList[n+4*Np];
            fq = dist[nread];
            psi += fq;

            // q = 6
            nread = neighborList[n+5*Np];
            fq = dist[nread];
            psi += fq;
            
            idx=Map[n];
            Psi[idx] = psi;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAeven_Poisson_ElectricPotential(int *Map, double *dist, double *Psi, int start, int finish, int Np){
	int n;
	double psi;//electric potential
	double fq;
    int idx;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {

            // q=0
            fq = dist[n];
            psi = fq;
            
            // q=1
            fq = dist[2*Np+n];
            psi += fq;

            // q=2 
            fq = dist[1*Np+n];
            psi += fq;

            // q=3
            fq = dist[4*Np+n];
            psi += fq;

            // q=4
            fq = dist[3*Np+n];
            psi += fq;

            // q=5
            fq = dist[6*Np+n];
            psi += fq;

            // q=6
            fq = dist[5*Np+n];
            psi += fq;

            idx=Map[n];
            Psi[idx] = psi;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAodd_Poisson(int *neighborList, int *Map, double *dist, double *Den_charge, double *Psi, double *ElectricField, double tau, double epsilon_LB,bool UseSlippingVelBC,int start, int finish, int Np){

	int n;
	double psi;//electric potential
    double Ex,Ey,Ez;//electric field
    double rho_e;//local charge density
	double f0,f1,f2,f3,f4,f5,f6;
	int nr1,nr2,nr3,nr4,nr5,nr6;
    double rlx=1.0/tau;
    int idx;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {

            //Load data
            //When Helmholtz-Smoluchowski slipping velocity BC is used, the bulk fluid is considered as electroneutral
            //and thus the net space charge density is zero. 
            rho_e = (UseSlippingVelBC==1) ? 0.0 : Den_charge[n] / epsilon_LB;
            idx=Map[n];
            psi = Psi[idx];
            
            // q=0
            f0 = dist[n];
            // q=1
            nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
            f1 = dist[nr1]; // reading the f1 data into register fq

            nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
            f2 = dist[nr2];  // reading the f2 data into register fq

            // q=3
            nr3 = neighborList[n+2*Np]; // neighbor 4
            f3 = dist[nr3];

            // q = 4
            nr4 = neighborList[n+3*Np]; // neighbor 3
            f4 = dist[nr4];

            // q=5
            nr5 = neighborList[n+4*Np];
            f5 = dist[nr5];

            // q = 6
            nr6 = neighborList[n+5*Np];
            f6 = dist[nr6];
            
            Ex = (f1-f2)*rlx*4.0;//NOTE the unit of electric field here is V/lu
            Ey = (f3-f4)*rlx*4.0;//factor 4.0 is D3Q7 lattice speed of sound
            Ez = (f5-f6)*rlx*4.0;
            ElectricField[n+0*Np] = Ex;
            ElectricField[n+1*Np] = Ey;
            ElectricField[n+2*Np] = Ez;

            // q = 0
            dist[n] = f0*(1.0-rlx) + 0.25*(rlx*psi+rho_e);

            // q = 1
            dist[nr2] = f1*(1.0-rlx) + 0.125*(rlx*psi+rho_e);

            // q = 2
            dist[nr1] = f2*(1.0-rlx) + 0.125*(rlx*psi+rho_e);

            // q = 3
            dist[nr4] = f3*(1.0-rlx) + 0.125*(rlx*psi+rho_e);

            // q = 4
            dist[nr3] = f4*(1.0-rlx) + 0.125*(rlx*psi+rho_e);

            // q = 5
            dist[nr6] = f5*(1.0-rlx) + 0.125*(rlx*psi+rho_e);

            // q = 6
            dist[nr5] = f6*(1.0-rlx) + 0.125*(rlx*psi+rho_e);
            //........................................................................
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAeven_Poisson(int *Map, double *dist, double *Den_charge, double *Psi, double *ElectricField, double tau, double epsilon_LB,bool UseSlippingVelBC,int start, int finish, int Np){

	int n;
	double psi;//electric potential
    double Ex,Ey,Ez;//electric field
    double rho_e;//local charge density
	double f0,f1,f2,f3,f4,f5,f6;
    double rlx=1.0/tau;
    int idx;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {

            //Load data
            //When Helmholtz-Smoluchowski slipping velocity BC is used, the bulk fluid is considered as electroneutral
            //and thus the net space charge density is zero. 
            rho_e = (UseSlippingVelBC==1) ? 0.0 : Den_charge[n] / epsilon_LB;
            idx=Map[n];
            psi = Psi[idx];

            f0 = dist[n];
            f1 = dist[2*Np+n];
            f2 = dist[1*Np+n];
            f3 = dist[4*Np+n];
            f4 = dist[3*Np+n];
            f5 = dist[6*Np+n];
            f6 = dist[5*Np+n];


            Ex = (f1-f2)*rlx*4.0;//NOTE the unit of electric field here is V/lu
            Ey = (f3-f4)*rlx*4.0;//factor 4.0 is D3Q7 lattice speed of sound
            Ez = (f5-f6)*rlx*4.0;
            ElectricField[n+0*Np] = Ex;
            ElectricField[n+1*Np] = Ey;
            ElectricField[n+2*Np] = Ez;

            // q = 0
            dist[n] = f0*(1.0-rlx) + 0.25*(rlx*psi+rho_e);

            // q = 1
            dist[1*Np+n] = f1*(1.0-rlx) + 0.125*(rlx*psi+rho_e); 

            // q = 2
            dist[2*Np+n] = f2*(1.0-rlx) + 0.125*(rlx*psi+rho_e); 

            // q = 3
            dist[3*Np+n] = f3*(1.0-rlx) + 0.125*(rlx*psi+rho_e);

            // q = 4
            dist[4*Np+n] = f4*(1.0-rlx) + 0.125*(rlx*psi+rho_e); 

            // q = 5
            dist[5*Np+n] = f5*(1.0-rlx) + 0.125*(rlx*psi+rho_e); 

            // q = 6
            dist[6*Np+n] = f6*(1.0-rlx) + 0.125*(rlx*psi+rho_e); 
            //........................................................................
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_Poisson_Init(int *Map, double *dist, double *Psi, int start, int finish, int Np){

	int n;
    int ijk;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
            ijk = Map[n];
            dist[0*Np+n] = 0.25*Psi[ijk];
            dist[1*Np+n] = 0.125*Psi[ijk];		
            dist[2*Np+n] = 0.125*Psi[ijk];	
            dist[3*Np+n] = 0.125*Psi[ijk];	
            dist[4*Np+n] = 0.125*Psi[ijk];	
            dist[5*Np+n] = 0.125*Psi[ijk];	
            dist[6*Np+n] = 0.125*Psi[ijk];	
		}
	}
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_ElectricPotential(int *neighborList,int *Map, double *dist, double *Psi, int start, int finish, int Np){

	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_AAodd_Poisson_ElectricPotential<<<NBLOCKS,NTHREADS >>>(neighborList,Map,dist,Psi,start,finish,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAodd_Poisson_ElectricPotential: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_ElectricPotential(int *Map, double *dist, double *Psi, int start, int finish, int Np){

	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_AAeven_Poisson_ElectricPotential<<<NBLOCKS,NTHREADS >>>(Map,dist,Psi,start,finish,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAeven_Poisson_ElectricPotential: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAodd_Poisson(int *neighborList, int *Map, double *dist, double *Den_charge, double *Psi, double *ElectricField, double tau, double epsilon_LB,bool UseSlippingVelBC,int start, int finish, int Np){

	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_AAodd_Poisson<<<NBLOCKS,NTHREADS >>>(neighborList,Map,dist,Den_charge,Psi,ElectricField,tau,epsilon_LB,UseSlippingVelBC,start,finish,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAodd_Poisson: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAeven_Poisson(int *Map, double *dist, double *Den_charge, double *Psi, double *ElectricField, double tau, double epsilon_LB,bool UseSlippingVelBC,int start, int finish, int Np){

	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_AAeven_Poisson<<<NBLOCKS,NTHREADS >>>(Map,dist,Den_charge,Psi,ElectricField,tau,epsilon_LB,UseSlippingVelBC,start,finish,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAeven_Poisson: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_Poisson_Init(int *Map, double *dist, double *Psi, int start, int finish, int Np){

	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_Poisson_Init<<<NBLOCKS,NTHREADS >>>(Map,dist,Psi,start,finish,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_Poisson_Init: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}
