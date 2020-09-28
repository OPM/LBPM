#include <stdio.h>
#include <math.h>
//#include <cuda_profiler_api.h>

#define NBLOCKS 1024
#define NTHREADS 256

__global__  void dvc_ScaLBL_D3Q7_AAodd_IonConcentration(int *neighborList, double *dist, double *Den, int start, int finish, int Np){
    int n,nread;
    double fq,Ci;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {

            // q=0
            fq = dist[n];
            Ci = fq;

            // q=1
            nread = neighborList[n]; 
            fq = dist[nread]; 
            Ci += fq;
            
            // q=2
            nread = neighborList[n+Np]; 
            fq = dist[nread];  
            Ci += fq;
            
            // q=3
            nread = neighborList[n+2*Np]; 
            fq = dist[nread];
            Ci += fq;
            
            // q=4
            nread = neighborList[n+3*Np]; 
            fq = dist[nread];
            Ci += fq;
            
            // q=5
            nread = neighborList[n+4*Np];
            fq = dist[nread];
            Ci += fq;
            
            // q=6
            nread = neighborList[n+5*Np];
            fq = dist[nread];
            Ci += fq;

            Den[n]=Ci;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAeven_IonConcentration(double *dist, double *Den, int start, int finish, int Np){
    int n;
    double fq,Ci;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {

            // q=0
            fq = dist[n];
            Ci = fq;
            
            // q=1
            fq = dist[2*Np+n];
            Ci += fq;

            // q=2
            fq = dist[1*Np+n];
            Ci += fq;

            // q=3
            fq = dist[4*Np+n];
            Ci += fq;

            // q=4
            fq = dist[3*Np+n];
            Ci += fq;

            // q=5
            fq = dist[6*Np+n];
            Ci += fq;

            // q=6
            fq = dist[5*Np+n];
            Ci += fq;

            Den[n]=Ci;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAodd_Ion(int *neighborList, double *dist, double *Den, double *Velocity, double *ElectricField, 
                                           double Di, int zi, double rlx, double Vt, int start, int finish, int Np){
	int n;
	double Ci;
    double ux,uy,uz;
    double uEPx,uEPy,uEPz;//electrochemical induced velocity
    double Ex,Ey,Ez;//electrical field
	double f0,f1,f2,f3,f4,f5,f6;
	int nr1,nr2,nr3,nr4,nr5,nr6;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {

            //Load data
            Ci=Den[n];
            Ex=ElectricField[n+0*Np];
            Ey=ElectricField[n+1*Np];
            Ez=ElectricField[n+2*Np];
            ux=Velocity[n+0*Np];
            uy=Velocity[n+1*Np];
            uz=Velocity[n+2*Np];
            uEPx=zi*Di/Vt*Ex;
            uEPy=zi*Di/Vt*Ey;
            uEPz=zi*Di/Vt*Ez;

            // q=0
            f0 = dist[n];
            // q=1
            nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
            f1 = dist[nr1]; // reading the f1 data into register fq
            // q=2
            nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
            f2 = dist[nr2];  // reading the f2 data into register fq
            // q=3
            nr3 = neighborList[n+2*Np]; // neighbor 4
            f3 = dist[nr3];
            // q=4
            nr4 = neighborList[n+3*Np]; // neighbor 3
            f4 = dist[nr4];
            // q=5
            nr5 = neighborList[n+4*Np];
            f5 = dist[nr5];
            // q=6
            nr6 = neighborList[n+5*Np];
            f6 = dist[nr6];
            
            // q=0
            dist[n] = f0*(1.0-rlx)+rlx*0.25*Ci;
            //dist[n] = f0*(1.0-rlx)+rlx*0.25*Ci*(1.0 - 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q = 1
            dist[nr2] = f1*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(ux+uEPx));
            //dist[nr2] = f1*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(ux+uEPx)+8.0*(ux+uEPx)*(ux+uEPx)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q=2
            dist[nr1] = f2*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(ux+uEPx));
            //dist[nr1] = f2*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(ux+uEPx)+8.0*(ux+uEPx)*(ux+uEPx)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q = 3
            dist[nr4] = f3*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(uy+uEPy));
            //dist[nr4] = f3*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(uy+uEPy)+8.0*(uy+uEPy)*(uy+uEPy)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q = 4
            dist[nr3] = f4*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(uy+uEPy));
            //dist[nr3] = f4*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(uy+uEPy)+8.0*(uy+uEPy)*(uy+uEPy)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q = 5
            dist[nr6] = f5*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(uz+uEPz));
            //dist[nr6] = f5*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(uz+uEPz)+8.0*(uz+uEPz)*(uz+uEPz)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q = 6
            dist[nr5] = f6*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(uz+uEPz));
            //dist[nr5] = f6*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(uz+uEPz)+8.0*(uz+uEPz)*(uz+uEPz)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAeven_Ion(double *dist, double *Den, double *Velocity, double *ElectricField, 
                                            double Di, int zi, double rlx, double Vt, int start, int finish, int Np){
	int n;
	double Ci;
    double ux,uy,uz;
    double uEPx,uEPy,uEPz;//electrochemical induced velocity
    double Ex,Ey,Ez;//electrical field
	double f0,f1,f2,f3,f4,f5,f6;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {

            //Load data
            Ci=Den[n];
            Ex=ElectricField[n+0*Np];
            Ey=ElectricField[n+1*Np];
            Ez=ElectricField[n+2*Np];
            ux=Velocity[n+0*Np];
            uy=Velocity[n+1*Np];
            uz=Velocity[n+2*Np];
            uEPx=zi*Di/Vt*Ex;
            uEPy=zi*Di/Vt*Ey;
            uEPz=zi*Di/Vt*Ez;

            f0 = dist[n];
            f1 = dist[2*Np+n];
            f2 = dist[1*Np+n];
            f3 = dist[4*Np+n];
            f4 = dist[3*Np+n];
            f5 = dist[6*Np+n];
            f6 = dist[5*Np+n];
            
            // q=0
            dist[n] = f0*(1.0-rlx)+rlx*0.25*Ci;
            //dist[n] = f0*(1.0-rlx)+rlx*0.25*Ci*(1.0 - 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q = 1
            dist[1*Np+n] = f1*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(ux+uEPx));
            //dist[1*Np+n] = f1*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(ux+uEPx)+8.0*(ux+uEPx)*(ux+uEPx)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q=2
            dist[2*Np+n] = f2*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(ux+uEPx));
            //dist[2*Np+n] = f2*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(ux+uEPx)+8.0*(ux+uEPx)*(ux+uEPx)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q = 3
            dist[3*Np+n] = f3*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(uy+uEPy));
            //dist[3*Np+n] = f3*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(uy+uEPy)+8.0*(uy+uEPy)*(uy+uEPy)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q = 4
            dist[4*Np+n] = f4*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(uy+uEPy));
            //dist[4*Np+n] = f4*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(uy+uEPy)+8.0*(uy+uEPy)*(uy+uEPy)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q = 5
            dist[5*Np+n] = f5*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(uz+uEPz));
            //dist[5*Np+n] = f5*(1.0-rlx) + rlx*0.125*Ci*(1.0+4.0*(uz+uEPz)+8.0*(uz+uEPz)*(uz+uEPz)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));

            // q = 6
            dist[6*Np+n] = f6*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(uz+uEPz));
            //dist[6*Np+n] = f6*(1.0-rlx) + rlx*0.125*Ci*(1.0-4.0*(uz+uEPz)+8.0*(uz+uEPz)*(uz+uEPz)- 2.0*((ux+uEPx)*(ux+uEPx) + (uy+uEPy)*(uy+uEPy) + (uz+uEPz)*(uz+uEPz)));
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_Ion_Init(double *dist, double *Den, double DenInit, int Np){

	int n;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<Np) {
            dist[0*Np+n] = 0.25*DenInit;
            dist[1*Np+n] = 0.125*DenInit;		
            dist[2*Np+n] = 0.125*DenInit;	
            dist[3*Np+n] = 0.125*DenInit;	
            dist[4*Np+n] = 0.125*DenInit;	
            dist[5*Np+n] = 0.125*DenInit;	
            dist[6*Np+n] = 0.125*DenInit;	
            Den[n] = DenInit;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_Ion_Init_FromFile(double *dist, double *Den, int Np){

	int n;
    double DenInit;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<Np) {
            DenInit = Den[n];
            dist[0*Np+n] = 0.25*DenInit;
            dist[1*Np+n] = 0.125*DenInit;		
            dist[2*Np+n] = 0.125*DenInit;	
            dist[3*Np+n] = 0.125*DenInit;	
            dist[4*Np+n] = 0.125*DenInit;	
            dist[5*Np+n] = 0.125*DenInit;	
            dist[6*Np+n] = 0.125*DenInit;	
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_Ion_ChargeDensity(double *Den, double *ChargeDensity, int IonValence, int ion_component, int start, int finish, int Np){

    int n;
    double Ci;//ion concentration of species i
    double CD;//charge density
    double CD_tmp;
    double F = 96485.0;//Faraday's constant; unit[C/mol]; F=e*Na, where Na is the Avogadro constant

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
            Ci = Den[n+ion_component*Np];
            CD = ChargeDensity[n];
            CD_tmp = F*IonValence*Ci;
            ChargeDensity[n] = CD*(ion_component>0) + CD_tmp;
		}
	}
}


extern "C" void ScaLBL_D3Q7_AAodd_IonConcentration(int *neighborList, double *dist, double *Den, int start, int finish, int Np){

	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_AAodd_IonConcentration<<<NBLOCKS,NTHREADS >>>(neighborList,dist,Den,start,finish,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAodd_IonConcentration: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAeven_IonConcentration(double *dist, double *Den, int start, int finish, int Np){

	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_AAeven_IonConcentration<<<NBLOCKS,NTHREADS >>>(dist,Den,start,finish,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAeven_IonConcentration: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion(int *neighborList, double *dist, double *Den, double *Velocity, double *ElectricField,  
                                      double Di, int zi, double rlx, double Vt, int start, int finish, int Np){
	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_AAodd_Ion<<<NBLOCKS,NTHREADS >>>(neighborList,dist,Den,Velocity,ElectricField,Di,zi,rlx,Vt,start,finish,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAodd_Ion: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion(double *dist, double *Den, double *Velocity, double *ElectricField, 
                                       double Di, int zi, double rlx, double Vt, int start, int finish, int Np){
	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_AAeven_Ion<<<NBLOCKS,NTHREADS >>>(dist,Den,Velocity,ElectricField,Di,zi,rlx,Vt,start,finish,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAeven_Ion: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_Ion_Init(double *dist, double *Den, double DenInit, int Np){

	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_Ion_Init<<<NBLOCKS,NTHREADS >>>(dist,Den,DenInit,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_Ion_Init: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_Ion_Init_FromFile(double *dist, double *Den, int Np){

	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_Ion_Init_FromFile<<<NBLOCKS,NTHREADS >>>(dist,Den,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_Ion_Init_FromFile: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_Ion_ChargeDensity(double *Den, double *ChargeDensity, int IonValence, int ion_component, int start, int finish, int Np){

	//cudaProfilerStart();
	dvc_ScaLBL_D3Q7_Ion_ChargeDensity<<<NBLOCKS,NTHREADS >>>(Den,ChargeDensity,IonValence,ion_component,start,finish,Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_Ion_ChargeDensity: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}
