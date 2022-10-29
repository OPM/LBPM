#include <stdio.h>
#include <math.h>
//#include <cuda_profiler_api.h>

#define NBLOCKS 1024
#define NTHREADS 512

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

__global__  void dvc_ScaLBL_D3Q19_AAeven_Poisson_ElectricPotential(
    int *Map, double *dist, double *Den_charge, double *Psi, double epsilon_LB, bool UseSlippingVelBC, int start, int finish, int Np) {
    int n;
    double psi,sum;        //electric potential
    double rho_e;      //local charge density
    double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15,
        f16, f17, f18;
    double Gs;
    int idx;

    for (n = start; n < finish; n++) {
        rho_e = (UseSlippingVelBC==1) ? 0.0 : Den_charge[n] / epsilon_LB;

        //........................................................................
        // q=0
        f0 = dist[n];
        f1 = dist[2 * Np + n];
        f2 = dist[1 * Np + n];
        f3 = dist[4 * Np + n];
        f4 = dist[3 * Np + n];
        f5 = dist[6 * Np + n];
        f6 = dist[5 * Np + n];
        f7 = dist[8 * Np + n];
        f8 = dist[7 * Np + n];
        f9 = dist[10 * Np + n];
        f10 = dist[9 * Np + n];
        f11 = dist[12 * Np + n];
        f12 = dist[11 * Np + n];
        f13 = dist[14 * Np + n];
        f14 = dist[13 * Np + n];
        f15 = dist[16 * Np + n];
        f16 = dist[15 * Np + n];
        f17 = dist[18 * Np + n];
        f18 = dist[17 * Np + n];

        psi = f0 + f2 + f1 + f4 + f3 + f6 + f5 + f8 + f7 + f10 + f9 + f12 +
                f11 + f14 + f13 + f16 + f15 + f18 + f17;

        idx = Map[n];

        Psi[idx] = psi - 0.5*rho_e;
    }
}

__global__  void dvc_ScaLBL_D3Q19_AAodd_Poisson(int *neighborList, int *Map,
		double *dist, double *Den_charge,
		double *Psi, double *ElectricField,
		double tau, double epsilon_LB, bool UseSlippingVelBC,
		int start, int finish, int Np) {
	int n;
	double psi;        //electric potential
	double Ex, Ey, Ez; //electric field
	double rho_e;      //local charge density
	double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15,
	f16, f17, f18;
	int nr1, nr2, nr3, nr4, nr5, nr6, nr7, nr8, nr9, nr10, nr11, nr12, nr13,
	nr14, nr15, nr16, nr17, nr18;
	double sum_q;
	double rlx = 1.0 / tau;
	int idx;

	double W0 = 0.5;
	double W1 = 1.0/24.0;
	double W2 = 1.0/48.0;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			//Load data
			//When Helmholtz-Smoluchowski slipping velocity BC is used, the bulk fluid is considered as electroneutral
			//and thus the net space charge density is zero. 
			rho_e = (UseSlippingVelBC==1) ? 0.0 : Den_charge[n] / epsilon_LB;

			// q=0
			f0 = dist[n];
			// q=1
			nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
			f1 = dist[nr1];        // reading the f1 data into register fq

			nr2 = neighborList[n + Np]; // neighbor 1 ( < 10Np => even part of dist)
			f2 = dist[nr2];             // reading the f2 data into register fq

			// q=3
			nr3 = neighborList[n + 2 * Np]; // neighbor 4
			f3 = dist[nr3];

			// q = 4
			nr4 = neighborList[n + 3 * Np]; // neighbor 3
			f4 = dist[nr4];

			// q=5
			nr5 = neighborList[n + 4 * Np];
			f5 = dist[nr5];

			// q = 6
			nr6 = neighborList[n + 5 * Np];
			f6 = dist[nr6];

			// q=7
			nr7 = neighborList[n + 6 * Np];
			f7 = dist[nr7];

			// q = 8
			nr8 = neighborList[n + 7 * Np];
			f8 = dist[nr8];

			// q=9
			nr9 = neighborList[n + 8 * Np];
			f9 = dist[nr9];

			// q = 10
			nr10 = neighborList[n + 9 * Np];
			f10 = dist[nr10];

			// q=11
			nr11 = neighborList[n + 10 * Np];
			f11 = dist[nr11];

			// q=12
			nr12 = neighborList[n + 11 * Np];
			f12 = dist[nr12];

			// q=13
			nr13 = neighborList[n + 12 * Np];
			f13 = dist[nr13];

			// q=14
			nr14 = neighborList[n + 13 * Np];
			f14 = dist[nr14];

			// q=15
			nr15 = neighborList[n + 14 * Np];
			f15 = dist[nr15];

			// q=16
			nr16 = neighborList[n + 15 * Np];
			f16 = dist[nr16];

			// q=17
			//fq = dist[18*Np+n];
			nr17 = neighborList[n + 16 * Np];
			f17 = dist[nr17];

			// q=18
			nr18 = neighborList[n + 17 * Np];
			f18 = dist[nr18];

			sum_q = f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+f17+f18;
			//error = 8.0*(sum_q - f0) + rho_e; 

			psi = 2.0*(f0*(1.0 - rlx) + rlx*(sum_q + 0.125*rho_e));

			idx = Map[n];
			Psi[idx] = psi;

			Ex = (f1 - f2 + 0.5*(f7 - f8 + f9 - f10 + f11 - f12 + f13 - f14))*4.0; //NOTE the unit of electric field here is V/lu
			Ey = (f3 - f4 + 0.5*(f7 - f8 - f9 + f10 + f15 - f16 + f17 - f18))*4.0;
			Ez = (f5 - f6 + 0.5*(f11 - f12 - f13 + f14 + f15 - f16 - f17 + f18))*4.0;
			ElectricField[n + 0 * Np] = Ex;
			ElectricField[n + 1 * Np] = Ey;
			ElectricField[n + 2 * Np] = Ez;

			// q = 0
			dist[n] = W0*psi; //f0 * (1.0 - rlx) -  (1.0-0.5*rlx)*W0*rho_e;

			// q = 1
			dist[nr2] = W1*psi; //f1 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			// q = 2
			dist[nr1] = W1*psi; //f2 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			// q = 3
			dist[nr4] = W1*psi; //f3 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			// q = 4
			dist[nr3] = W1*psi; //f4 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			// q = 5
			dist[nr6] = W1*psi; //f5 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			// q = 6
			dist[nr5] = W1*psi; //f6 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;
			//........................................................................

			// q = 7
			dist[nr8] = W2*psi; //f7 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q = 8
			dist[nr7] = W2*psi; //f8 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q = 9
			dist[nr10] = W2*psi; //f9 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q = 10
			dist[nr9] = W2*psi; //f10 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q = 11
			dist[nr12] = W2*psi; //f11 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q = 12
			dist[nr11] = W2*psi; //f12 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q = 13
			dist[nr14] = W2*psi; //f13 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q= 14
			dist[nr13] = W2*psi; //f14 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q = 15
			dist[nr16] = W2*psi; //f15 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q = 16
			dist[nr15] = W2*psi; //f16 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q = 17
			dist[nr18] = W2*psi; //f17 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;

			// q = 18
			dist[nr17] = W2*psi; //f18 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAeven_Poisson(int *Map, double *dist,
		double *Den_charge, double *Psi,
		double *ElectricField, double *Error, double tau,
		double epsilon_LB, bool UseSlippingVelBC,
		int start, int finish, int Np) {
	int n;
	double psi;        //electric potential
	double Ex, Ey, Ez; //electric field
	double rho_e;      //local charge density
	double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15,
	f16, f17, f18;
	double error,sum_q;
	double rlx = 1.0 / tau;
	int idx;
	double W0 = 0.5;
	double W1 = 1.0/24.0;
	double W2 = 1.0/48.0;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			//Load data
			//When Helmholtz-Smoluchowski slipping velocity BC is used, the bulk fluid is considered as electroneutral
			//and thus the net space charge density is zero. 
			rho_e = (UseSlippingVelBC==1) ? 0.0 : Den_charge[n] / epsilon_LB;

			f0 = dist[n];
			f1 = dist[2 * Np + n];
			f2 = dist[1 * Np + n];
			f3 = dist[4 * Np + n];
			f4 = dist[3 * Np + n];
			f5 = dist[6 * Np + n];
			f6 = dist[5 * Np + n];

			f7 = dist[8 * Np + n];
			f8 = dist[7 * Np + n];
			f9 = dist[10 * Np + n];
			f10 = dist[9 * Np + n];
			f11 = dist[12 * Np + n];
			f12 = dist[11 * Np + n];
			f13 = dist[14 * Np + n];
			f14 = dist[13 * Np + n];
			f15 = dist[16 * Np + n];
			f16 = dist[15 * Np + n];
			f17 = dist[18 * Np + n];
			f18 = dist[17 * Np + n];

			/* Ex = (f1 - f2) * rlx *
	             4.0; //NOTE the unit of electric field here is V/lu
	        Ey = (f3 - f4) * rlx *
	             4.0; //factor 4.0 is D3Q7 lattice squared speed of sound
	        Ez = (f5 - f6) * rlx * 4.0;
			 */
			Ex = (f1 - f2 + 0.5*(f7 - f8 + f9 - f10 + f11 - f12 + f13 - f14))*4.0; //NOTE the unit of electric field here is V/lu
			Ey = (f3 - f4 + 0.5*(f7 - f8 - f9 + f10 + f15 - f16 + f17 - f18))*4.0;
			Ez = (f5 - f6 + 0.5*(f11 - f12 - f13 + f14 + f15 - f16 - f17 + f18))*4.0;
			ElectricField[n + 0 * Np] = Ex;
			ElectricField[n + 1 * Np] = Ey;
			ElectricField[n + 2 * Np] = Ez;

			sum_q = f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+f17+f18;
			error = 8.0*(sum_q - f0) + rho_e; 
			
			Error[n] = error;

			psi = 2.0*(f0*(1.0 - rlx) + rlx*(sum_q + 0.125*rho_e));
	        
	        idx = Map[n];
	        Psi[idx] = psi;
	        
			// q = 0
			dist[n] =  W0*psi;//

			// q = 1
			dist[1 * Np + n] =  W1*psi;//f1 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			// q = 2
			dist[2 * Np + n] =  W1*psi;//f2 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			// q = 3
			dist[3 * Np + n] =  W1*psi;//f3 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			// q = 4
			dist[4 * Np + n] =  W1*psi;//f4 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			// q = 5
			dist[5 * Np + n] =  W1*psi;//f5 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			// q = 6
			dist[6 * Np + n] =  W1*psi;//f6 * (1.0 - rlx) +W1* (rlx * psi) - (1.0-0.5*rlx)*0.05555555555555555*rho_e;

			dist[7 * Np + n] =  W2*psi;//f7 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[8 * Np + n] =  W2*psi;//f8* (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[9 * Np + n] =  W2*psi;//f9 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[10 * Np + n] = W2*psi;//f10 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[11 * Np + n] = W2*psi;//f11 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[12 * Np + n] = W2*psi;//f12 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[13 * Np + n] = W2*psi;//f13 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[14 * Np + n] = W2*psi;//f14 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[15 * Np + n] = W2*psi;//f15 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[16 * Np + n] = W2*psi;//f16 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[17 * Np + n] = W2*psi;//f17 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			dist[18 * Np + n] = W2*psi;//f18 * (1.0 - rlx) +W2* (rlx * psi) - (1.0-0.5*rlx)*0.02777777777777778*rho_e;
			//........................................................................
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_Poisson_Init(int *Map, double *dist, double *Psi,
		int start, int finish, int Np) {
	int n;
	int ijk;
	double W0 = 0.5;
	double W1 = 1.0/24.0;
	double W2 = 1.0/48.0;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			ijk = Map[n];
			dist[0 * Np + n] = W0 * Psi[ijk];//3333333333333333* Psi[ijk];
			dist[1 * Np + n] = W1 * Psi[ijk];
			dist[2 * Np + n] = W1 * Psi[ijk];
			dist[3 * Np + n] = W1 * Psi[ijk];
			dist[4 * Np + n] = W1 * Psi[ijk];
			dist[5 * Np + n] = W1 * Psi[ijk];
			dist[6 * Np + n] = W1 * Psi[ijk];
			dist[7 * Np + n] = W2* Psi[ijk];
			dist[8 * Np + n] = W2* Psi[ijk];
			dist[9 * Np + n] = W2* Psi[ijk];
			dist[10 * Np + n] = W2* Psi[ijk];
			dist[11 * Np + n] = W2* Psi[ijk];
			dist[12 * Np + n] = W2* Psi[ijk];
			dist[13 * Np + n] = W2* Psi[ijk];
			dist[14 * Np + n] = W2* Psi[ijk];
			dist[15 * Np + n] = W2* Psi[ijk];
			dist[16 * Np + n] = W2* Psi[ijk];
			dist[17 * Np + n] = W2* Psi[ijk];
			dist[18 * Np + n] = W2* Psi[ijk];
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAeven_Poisson_Potential_BC_z(int *list,  double *dist, double Vin, int count, int Np) {
	
	double W1 = 1.0/24.0;
	double W2 = 1.0/48.0;
	
	int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < count){
        int n = list[idx];
        
        dist[6 * Np + n]  = W1*Vin;
        dist[12 * Np + n] = W2*Vin;
        dist[13 * Np + n] = W2*Vin;
        dist[16 * Np + n] = W2*Vin;
        dist[17 * Np + n] = W2*Vin;
    }
}

__global__  void dvc_ScaLBL_D3Q19_AAeven_Poisson_Potential_BC_Z(int *list, double *dist, double Vout, int count, int Np) {
	
	double W1 = 1.0/24.0;
	double W2 = 1.0/48.0;

	int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < count){		
        int n = list[idx];
        dist[5 * Np + n]  = W1*Vout;
        dist[11 * Np + n] = W2*Vout;
        dist[14 * Np + n] = W2*Vout;
        dist[15 * Np + n] = W2*Vout;
        dist[18 * Np + n] = W2*Vout;
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAodd_Poisson_Potential_BC_z(int *d_neighborList, int *list, double *dist, double Vin, int count, int Np) {
	
	double W1 = 1.0/24.0;
	double W2 = 1.0/48.0;
    	int nr5, nr11, nr14, nr15, nr18;
    
	int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < count){
		int n = list[idx];

        
        // Unknown distributions
        nr5 = d_neighborList[n + 4 * Np];
        nr11 = d_neighborList[n + 10 * Np];
        nr15 = d_neighborList[n + 14 * Np];
        nr14 = d_neighborList[n + 13 * Np];
        nr18 = d_neighborList[n + 17 * Np];
        
        dist[nr5]  = W1*Vin;
        dist[nr11] = W2*Vin;
        dist[nr15] = W2*Vin;
        dist[nr14] = W2*Vin;
        dist[nr18] = W2*Vin;
    }
}

__global__  void dvc_ScaLBL_D3Q19_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np)  {
	
	double W1 = 1.0/24.0;
	double W2 = 1.0/48.0;
    	int nr6, nr12, nr13, nr16, nr17;

	int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < count){		
        int n = list[idx];
        // unknown distributions
        nr6 = d_neighborList[n + 5 * Np];
        nr12 = d_neighborList[n + 11 * Np];
        nr16 = d_neighborList[n + 15 * Np];
        nr17 = d_neighborList[n + 16 * Np];
        nr13 = d_neighborList[n + 12 * Np];
        
        dist[nr6]  = W1*Vout;
        dist[nr12] = W2*Vout;
        dist[nr16] = W2*Vout;
        dist[nr17] = W2*Vout;
        dist[nr13] = W2*Vout;
	}
}

/* wrapper functions to launch kernels */
extern "C" void ScaLBL_D3Q19_AAeven_Poisson_Potential_BC_z(int *list,  double *dist, double Vin, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_AAeven_Poisson_Potential_BC_z<<<GRID,512>>>(list, dist, Vin, count, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_Poisson_Potential_BC_z (kernel): %s \n",cudaGetErrorString(err));
	}
}
//
extern "C" void ScaLBL_D3Q19_AAeven_Poisson_Potential_BC_Z(int *list, double *dist,  double Vout, int count, int Np){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_AAeven_Poisson_Potential_BC_Z<<<GRID,512>>>(list, dist, Vout, count, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_Poisson_Potential_BC_Z (kernel): %s \n",cudaGetErrorString(err));
	}
}
extern "C" void ScaLBL_D3Q19_AAodd_Poisson_Potential_BC_z(int *d_neighborList, int *list, double *dist, double Vin, int count,int Np) {
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_AAodd_Poisson_Potential_BC_z<<<GRID,512>>>(d_neighborList, list, dist, Vin, count, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_Poisson_Potential_BC_z (kernel): %s \n",cudaGetErrorString(err));
	}
}
//
extern "C" void ScaLBL_D3Q19_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np)  {

	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q19_AAodd_Poisson_Potential_BC_Z<<<GRID,512>>>(d_neighborList, list, dist, Vout, count, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_Poisson_Potential_BC_Z (kernel): %s \n",cudaGetErrorString(err));
	}
}


extern "C" void ScaLBL_D3Q19_AAodd_Poisson(int *neighborList, int *Map,
		double *dist, double *Den_charge,
		double *Psi, double *ElectricField, 
		double tau, double epsilon_LB, bool UseSlippingVelBC,
		int start, int finish, int Np) {
	//cudaProfilerStart();
	dvc_ScaLBL_D3Q19_AAodd_Poisson<<<NBLOCKS,NTHREADS >>>(neighborList, Map,
			dist, Den_charge, Psi, ElectricField, tau, epsilon_LB, UseSlippingVelBC, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in dvc_ScaLBL_D3Q19_AAodd_Poisson: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_Poisson(int *Map, double *dist,
		double *Den_charge, double *Psi,
		double *ElectricField, double *Error, double tau,
		double epsilon_LB, bool UseSlippingVelBC,
		int start, int finish, int Np) {

	dvc_ScaLBL_D3Q19_AAeven_Poisson<<<NBLOCKS,NTHREADS >>>( Map, dist, Den_charge, Psi,
			ElectricField, Error, tau, epsilon_LB, UseSlippingVelBC, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in dvc_ScaLBL_D3Q19_AAeven_Poisson: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_Poisson_Init(int *Map, double *dist, double *Psi,
		int start, int finish, int Np){
	//cudaProfilerStart();

	dvc_ScaLBL_D3Q19_Poisson_Init<<<NBLOCKS,NTHREADS >>>(Map, dist, Psi, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Poisson_Init: %s \n",cudaGetErrorString(err));
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
