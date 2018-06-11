/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
// GPU Functions for D3Q7 Lattice Boltzmann Methods

#define NBLOCKS 560
#define NTHREADS 128

__global__  void dvc_ScaLBL_Scalar_Pack(int *list, int count, double *sendbuf, double *Data, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		n = list[idx];
		sendbuf[idx] = Data[n];
	}
}
__global__  void dvc_ScaLBL_Scalar_Unpack(int *list, int count, double *recvbuf, double *Data, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		n = list[idx];
		Data[n] = recvbuf[idx];
	}
}

__global__  void dvc_ScaLBL_PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N){
	//....................................................................................
	// Pack distribution into the send buffer for the listed lattice sites
	//....................................................................................
	int idx,n,component;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		for (component=0; component<number; component++){
			n = list[idx];
			sendbuf[idx*number+component] = Data[number*n+component];
			Data[number*n+component] = 0.0;	// Set the data value to zero once it's in the buffer!
		}
	}
}


__global__ void dvc_ScaLBL_UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N){
	//....................................................................................
	// Unack distribution from the recv buffer
	// Sum to the existing density value
	//....................................................................................
	int idx,n,component;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
			for (component=0; component<number; component++){
			n = list[idx];
			Data[number*n+component] += recvbuf[idx*number+component];
		}
	}
}

__global__ void dvc_ScaLBL_D3Q7_Unpack(int q,  int *list,  int start, int count,
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
		n = list[idx];
		// unpack the distribution to the proper location
		if (!(n<0)) { dist[q*N+n] = recvbuf[start+idx];
		//printf("%f \n",,dist[q*N+n]);
		}
	}
}

__global__ void dvc_ScaLBL_D3Q7_Init(char *ID, double *f_even, double *f_odd, double *Den, int Nx, int Ny, int Nz)
{
	int n,N;
	N = Nx*Ny*Nz;
	double value;
	char id;
	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N){
		   id = ID[n];
			if (id > 0){
				value = Den[n];
				f_even[n] = 0.3333333333333333*value;
				f_odd[n] = 0.1111111111111111*value;		//double(100*n)+1.f;
				f_even[N+n] = 0.1111111111111111*value;	//double(100*n)+2.f;
				f_odd[N+n] = 0.1111111111111111*value;	//double(100*n)+3.f;
				f_even[2*N+n] = 0.1111111111111111*value;	//double(100*n)+4.f;
				f_odd[2*N+n] = 0.1111111111111111*value;	//double(100*n)+5.f;
				f_even[3*N+n] = 0.1111111111111111*value;	//double(100*n)+6.f;
			}
			else{
				for(int q=0; q<3; q++){
					f_even[q*N+n] = -1.0;
					f_odd[q*N+n] = -1.0;
				}
				f_even[3*N+n] = -1.0;
			}
		}
	}
}

//*************************************************************************
__global__  void dvc_ScaLBL_D3Q7_Swap(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz)
{
	int i,j,k,n,nn,N;
	// distributions
	double f1,f2,f3,f4,f5,f6;
	char id;
	N = Nx*Ny*Nz;

	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;

		if (n<N ){
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
                //........................................................................

                //........................................................................
                // Retrieve odd distributions from neighboring nodes (swap convention)
                //........................................................................
                nn = n+1;							// neighbor index (pull convention)
                if (!(i+1<Nx))	nn -= Nx;			// periodic BC along the x-boundary
                //if (i+1<Nx){
                f2 = disteven[N+nn];					// pull neighbor for distribution 2
                if (!(f2 < 0.0)){
                    distodd[n] = f2;
                    disteven[N+nn] = f1;
                }
                //}
                //........................................................................
                nn = n+Nx;							// neighbor index (pull convention)
                if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
                //if (j+1<Ny){
                f4 = disteven[2*N+nn];				// pull neighbor for distribution 4
                if (!(f4 < 0.0)){
                    distodd[N+n] = f4;
                    disteven[2*N+nn] = f3;
                }
                //........................................................................
                nn = n+Nx*Ny;						// neighbor index (pull convention)
                if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
                //if (k+1<Nz){
                f6 = disteven[3*N+nn];				// pull neighbor for distribution 6
                if (!(f6 < 0.0)){
                    distodd[2*N+n] = f6;
                    disteven[3*N+nn] = f5;
                }
			}
		}
	}
}

//*************************************************************************
__global__  void dvc_ScaLBL_D3Q7_Density(char *ID, double *disteven, double *distodd, double *Den,
		int Nx, int Ny, int Nz)
{
	char id;
	int n;
	double f0,f1,f2,f3,f4,f5,f6;
	int N = Nx*Ny*Nz;

	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N){
            id = ID[n];
            if (id > 0 ){
                // Read the distributions
                f0 = disteven[n];
                f2 = disteven[N+n];
                f4 = disteven[2*N+n];
                f6 = disteven[3*N+n];
                f1 = distodd[n];
                f3 = distodd[N+n];
                f5 = distodd[2*N+n];
                // Compute the density
                Den[n] = f0+f1+f2+f3+f4+f5+f6;
            }
		}
	}
}

extern "C" void ScaLBL_D3Q7_Unpack(int q, int *list,  int start, int count, double *recvbuf, double *dist, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_D3Q7_Unpack <<<GRID,512 >>>(q, list, start, count, recvbuf, dist, N);
}

extern "C" void ScaLBL_Scalar_Pack(int *list, int count, double *sendbuf, double *Data, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_Scalar_Pack <<<GRID,512 >>>(list, count, sendbuf, Data, N);
}

extern "C" void ScaLBL_Scalar_Unpack(int *list, int count, double *recvbuf, double *Data, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_Scalar_Unpack <<<GRID,512 >>>(list, count, recvbuf, Data, N);
}
extern "C" void ScaLBL_PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_PackDenD3Q7 <<<GRID,512 >>>(list, count, sendbuf, number, Data, N);
}

extern "C" void ScaLBL_UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_UnpackDenD3Q7 <<<GRID,512 >>>(list, count, recvbuf, number, Data, N);
}

extern "C" void ScaLBL_D3Q7_Init(char *ID, double *f_even, double *f_odd, double *Den, int Nx, int Ny, int Nz){
	dvc_ScaLBL_D3Q7_Init <<<NBLOCKS,NTHREADS >>>(ID, f_even, f_odd, Den, Nx, Ny, Nz);
}

extern "C" void ScaLBL_D3Q7_Swap(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz){
	dvc_ScaLBL_D3Q7_Swap <<<NBLOCKS,NTHREADS >>>(ID, disteven, distodd, Nx, Ny, Nz);
}

extern "C" void ScaLBL_D3Q7_Density(char *ID, double *disteven, double *distodd, double *Den,
										int Nx, int Ny, int Nz){
	dvc_ScaLBL_D3Q7_Density <<<NBLOCKS,NTHREADS >>>(ID, disteven, distodd, Den, Nx, Ny,  Nz);
}

