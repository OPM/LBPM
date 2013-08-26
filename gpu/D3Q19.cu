#include <cuda.h>

__global__ void InitD3Q19(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz, int S)
{
	int n,N;
	N = Nx*Ny*Nz;

	for (int s=0; s<S; s++){

		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;

		if (n<N){
			if (ID[n] > 0){
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

__global__ void PackDist(int q, int *list, int start, int count, double *sendbuf, double *dist, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
//	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[start+idx] = dist[q*N+n];
	}
}


__global__ void MapRecvDist(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
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
//	for (idx=0; idx<count; idx++){
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
	}
}

//*************************************************************************
__global__ void SwapD3Q19(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz, int S)
{
	int n,nn,N;
	// distributions
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	N = Nx*Ny*Nz;

	// S - number of threadblocks per grid block
	for (int s=0; s<S; s++){

		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;

//	for (n=0; n<N; n++){
		//.......Back out the 3-D indices for node n..............
		int	k = n/(Nx*Ny);
		int j = (n-Nx*Ny*k)/Nx;
		int i = n-Nx*Ny*k-Nz*j;

		if (n<N){
			if (ID[n] > 0){
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
}
//*************************************************************************
extern "C" void dvc_PackD3Q19(int faceGrid, int edgeGrid, int threads,double *f_even, double *f_odd, int N,
		int *dvcSendList_x, int *dvcSendList_y, int *dvcSendList_z, int *dvcSendList_X, int *dvcSendList_Y, int *dvcSendList_Z,
		int *dvcSendList_xy, int *dvcSendList_XY, int *dvcSendList_xY, int *dvcSendList_Xy,
		int *dvcSendList_xz, int *dvcSendList_XZ, int *dvcSendList_xZ, int *dvcSendList_Xz,
		int *dvcSendList_yz, int *dvcSendList_YZ, int *dvcSendList_yZ, int *dvcSendList_Yz,
		int sendCount_x, int sendCount_y, int sendCount_z, int sendCount_X, int sendCount_Y, int sendCount_Z,
		int sendCount_xy, int sendCount_XY, int sendCount_xY, int sendCount_Xy,
		int sendCount_xz, int sendCount_XZ, int sendCount_xZ, int sendCount_Xz,
		int sendCount_yz, int sendCount_YZ, int sendCount_yZ, int sendCount_Yz,
		double *sendbuf_x, double *sendbuf_y, double *sendbuf_z, double *sendbuf_X, double *sendbuf_Y, double *sendbuf_Z,
		double *sendbuf_xy, double *sendbuf_XY, double *sendbuf_xY, double *sendbuf_Xy,
		double *sendbuf_xz, double *sendbuf_XZ, double *sendbuf_xZ, double *sendbuf_Xz,
		double *sendbuf_yz, double *sendbuf_YZ, double *sendbuf_yZ, double *sendbuf_Yz)
{
	//...................................................................................
	PackDist<<<faceGrid,threads>>>(1,dvcSendList_x,0,sendCount_x,sendbuf_x,f_even,N);
	PackDist<<<faceGrid,threads>>>(4,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	PackDist<<<faceGrid,threads>>>(5,dvcSendList_x,2*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	PackDist<<<faceGrid,threads>>>(6,dvcSendList_x,3*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	PackDist<<<faceGrid,threads>>>(7,dvcSendList_x,4*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	//...Packing for X face<<<faceGrid,threads>>>(1,7,9,11,13)................................
	PackDist<<<faceGrid,threads>>>(0,dvcSendList_X,0,sendCount_X,sendbuf_X,f_odd,N);
	PackDist<<<faceGrid,threads>>>(3,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	PackDist<<<faceGrid,threads>>>(4,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	PackDist<<<faceGrid,threads>>>(5,dvcSendList_X,3*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	PackDist<<<faceGrid,threads>>>(6,dvcSendList_X,4*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	//...Packing for y face<<<faceGrid,threads>>>(4,8,9,16,18).................................
	PackDist<<<faceGrid,threads>>>(2,dvcSendList_y,0,sendCount_y,sendbuf_y,f_even,N);
	PackDist<<<faceGrid,threads>>>(4,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	PackDist<<<faceGrid,threads>>>(4,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,f_odd,N);
	PackDist<<<faceGrid,threads>>>(8,dvcSendList_y,3*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	PackDist<<<faceGrid,threads>>>(9,dvcSendList_y,4*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	//...Packing for Y face<<<faceGrid,threads>>>(3,7,10,15,17).................................
	PackDist<<<faceGrid,threads>>>(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,f_odd,N);
	PackDist<<<faceGrid,threads>>>(3,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	PackDist<<<faceGrid,threads>>>(5,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,f_even,N);
	PackDist<<<faceGrid,threads>>>(7,dvcSendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	PackDist<<<faceGrid,threads>>>(8,dvcSendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	//...Packing for z face<<<faceGrid,threads>>>(6,12,13,16,17)................................
	PackDist<<<faceGrid,threads>>>(3,dvcSendList_z,0,sendCount_z,sendbuf_z,f_even,N);
	PackDist<<<faceGrid,threads>>>(6,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,f_even,N);
	PackDist<<<faceGrid,threads>>>(6,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
	PackDist<<<faceGrid,threads>>>(8,dvcSendList_z,3*sendCount_z,sendCount_z,sendbuf_z,f_even,N);
	PackDist<<<faceGrid,threads>>>(8,dvcSendList_z,4*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
	//...Packing for Z face<<<faceGrid,threads>>>(5,11,14,15,18)................................
	PackDist<<<faceGrid,threads>>>(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,f_odd,N);
	PackDist<<<faceGrid,threads>>>(5,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
	PackDist<<<faceGrid,threads>>>(7,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
	PackDist<<<faceGrid,threads>>>(7,dvcSendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
	PackDist<<<faceGrid,threads>>>(9,dvcSendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
	//...Pack the xy edge <<<edgeGrid,threads>>>(8)................................
	PackDist<<<edgeGrid,threads>>>(4,dvcSendList_xy,0,sendCount_xy,sendbuf_xy,f_even,N);
	//...Pack the Xy edge <<<edgeGrid,threads>>>(9)................................
	PackDist<<<edgeGrid,threads>>>(4,dvcSendList_Xy,0,sendCount_Xy,sendbuf_Xy,f_odd,N);
	//...Pack the xY edge <<<edgeGrid,threads>>>(10)................................
	PackDist<<<edgeGrid,threads>>>(5,dvcSendList_xY,0,sendCount_xY,sendbuf_xY,f_even,N);
	//...Pack the XY edge <<<edgeGrid,threads>>>(7)................................
	PackDist<<<edgeGrid,threads>>>(3,dvcSendList_XY,0,sendCount_XY,sendbuf_XY,f_odd,N);
	//...Pack the xz edge <<<edgeGrid,threads>>>(12)................................
	PackDist<<<edgeGrid,threads>>>(6,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
	//...Pack the xZ edge <<<edgeGrid,threads>>>(14)................................
	PackDist<<<edgeGrid,threads>>>(7,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
	//...Pack the Xz edge <<<edgeGrid,threads>>>(13)................................
	PackDist<<<edgeGrid,threads>>>(6,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
	//...Pack the XZ edge <<<edgeGrid,threads>>>(11)................................
	PackDist<<<edgeGrid,threads>>>(5,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
	//...Pack the xz edge <<<edgeGrid,threads>>>(12)................................
	PackDist<<<edgeGrid,threads>>>(6,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
	//...Pack the xZ edge <<<edgeGrid,threads>>>(14)................................
	PackDist<<<edgeGrid,threads>>>(7,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
	//...Pack the Xz edge <<<edgeGrid,threads>>>(13)................................
	PackDist<<<edgeGrid,threads>>>(6,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
	//...Pack the XZ edge <<<edgeGrid,threads>>>(11)................................
	PackDist<<<edgeGrid,threads>>>(5,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
	//...Pack the yz edge <<<edgeGrid,threads>>>(16)................................
	PackDist<<<edgeGrid,threads>>>(8,dvcSendList_yz,0,sendCount_yz,sendbuf_yz,f_even,N);
	//...Pack the yZ edge <<<edgeGrid,threads>>>(18)................................
	PackDist<<<edgeGrid,threads>>>(9,dvcSendList_yZ,0,sendCount_yZ,sendbuf_yZ,f_even,N);
	//...Pack the Yz edge <<<edgeGrid,threads>>>(17)................................
	PackDist<<<edgeGrid,threads>>>(8,dvcSendList_Yz,0,sendCount_Yz,sendbuf_Yz,f_odd,N);
	//...Pack the YZ edge <<<edgeGrid,threads>>>(15)................................
	PackDist<<<edgeGrid,threads>>>(7,dvcSendList_YZ,0,sendCount_YZ,sendbuf_YZ,f_odd,N);
}
//...................................................................................
//*************************************************************************
extern "C" void dvc_PackDist(int grid, int threads, int q, int *SendList, int start,
		int sendCount, double *sendbuf, double *Dist, int N)
{
	//...................................................................................
	PackDist<<<grid,threads>>>(q,SendList,start,sendCount,sendbuf,Dist,N);
}
//*************************************************************************
extern "C" void dvc_UnpackDist(int grid, int threads, int q, int Cqx, int Cqy, int Cqz, int *RecvList, int start,
		int recvCount, double *recvbuf, double *Dist, int Nx, int Ny, int Nz)
{
	//...................................................................................
	MapRecvDist<<<grid,threads>>>(q,Cqx,Cqy,Cqz,RecvList,start,recvCount,recvbuf,Dist,Nx,Ny,Nz);
}
//*************************************************************************
extern "C" void dvc_SwapD3Q19( int nblocks, int nthreads, int S,
		char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz)
{
	SwapD3Q19 <<< nblocks, nthreads >>> (ID, f_even, f_odd, Nx, Ny, Nz, S);
}
//*************************************************************************
extern "C" void dvc_InitD3Q19(int nblocks, int nthreads, int S, char *ID, double *f_even, double *f_odd, int Nx,
							  int Ny, int Nz)
{
	InitD3Q19 <<<nblocks, nthreads>>>  (ID, f_even, f_odd, Nx, Ny, Nz, S);
}
//*************************************************************************

