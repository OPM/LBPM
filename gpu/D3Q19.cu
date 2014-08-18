#define NBLOCKS 32
#define NTHREADS 128

__global__  void dvc_PackDist(int q, int *list, int start, int count, double *sendbuf, double *dist, int N){
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

__global__ void dvc_UnpackDist(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
					   double *recvbuf, double *dist, int Nx, int Ny, int Nz){
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
		i = n-Nx*Ny*k-Nz*j;
		// Streaming for the non-local distribution
		i += Cqx;
		j += Cqy;
		k += Cqz;
/*		if (i < 0) i += Nx;
		if (j < 0) j += Ny;
		if (k < 0) k += Nz;
		if (!(i<Nx)) i -= Nx;
		if (!(j<Ny)) j -= Ny;
		if (!(k<Nz)) k -= Nz;
*/
		nn = k*Nx*Ny+j*Nx+i;
		// unpack the distribution to the proper location
	//	if (recvbuf[start+idx] != dist[q*N+nn]){
	//		printf("Stopping to check error \n");
	//		printf("recvbuf[start+idx] = %f \n",recvbuf[start+idx]);
	//		printf("dist[q*N+nn] = %f \n",dist[q*N+nn]);
	//		printf("A bug! Again? \n");
	//		idx = count;
	//	}
//		list[idx] = nn;
		dist[q*N+nn] = recvbuf[start+idx];
//		if (dist[q*N+nn] > 0.0) dist[q*N+nn] = recvbuf[start+idx];
	}
}

__global__ void dvc_InitD3Q19(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz)
{
	int n,N;
	N = Nx*Ny*Nz;
	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N && ID[n] > 0){
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

//*************************************************************************
__global__  void dvc_SwapD3Q19(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz)
{
	int i,j,k,n,nn,N;
	// distributions
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	
	N = Nx*Ny*Nz;
	
	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<N && ID[n] > 0){
			//.......Back out the 3-D indices for node n..............
			k = n/(Nx*Ny);
			j = (n-Nx*Ny*k)/Nx;
			i = n-Nx*Ny*k-Nz*j;
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
			if (f2 > 0){
				distodd[n] = f2;
				disteven[N+nn] = f1;
			}
			//}
			//........................................................................
			nn = n+Nx;							// neighbor index (pull convention)
			if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
			//if (j+1<Ny){
			f4 = disteven[2*N+nn];				// pull neighbor for distribution 4
			if (f4 > 0){
				distodd[N+n] = f4;
				disteven[2*N+nn] = f3;
				//	}
			}
			//........................................................................
			nn = n+Nx*Ny;						// neighbor index (pull convention)
			if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
			//if (k+1<Nz){
			f6 = disteven[3*N+nn];				// pull neighbor for distribution 6
			if (f6 > 0){
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
			if (f8 > 0){
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
			if (f10 > 0){
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
			if (f12 > 0){
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
			if (f14 > 0){
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
			if (f16 > 0){
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
			if (f18 > 0){
				distodd[8*N+n] = f18;
				disteven[9*N+nn] = f17;
				//	}
			}
			//........................................................................
			
		}
	}
}

extern "C" void PackDist(int q, int *list, int start, int count, double *sendbuf, double *dist, int N){
	int GRID = count / 512 + 1;
	dvc_PackDist <<<GRID,512 >>>(q, list, start, count, sendbuf, dist, N);
}
extern "C" void UnpackDist(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
			double *recvbuf, double *dist, int Nx, int Ny, int Nz){
	int GRID = count / 512 + 1;
	dvc_UnpackDist <<<GRID,512 >>>(q, Cqx, Cqy, Cqz, list, start, count, recvbuf, dist, Nx, Ny, Nz);
}
//*************************************************************************
extern "C" void InitD3Q19(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz){
	dvc_InitD3Q19<<<NBLOCKS,NTHREADS >>>(ID, f_even, f_odd, Nx, Ny, Nz);
}
extern "C" void SwapD3Q19(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz){
	dvc_SwapD3Q19<<<NBLOCKS,NTHREADS >>>(ID, disteven, distodd, Nx, Ny, Nz);

}

