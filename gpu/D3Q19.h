
//*************************************************************************
//*************************************************************************
extern "C" void dvc_InitD3Q19(int nblocks, int nthreads, int S, char *ID, double *f_even, double *f_odd, int Nx,
							  int Ny, int Nz);
//*************************************************************************
extern "C" void dvc_SwapD3Q19(int nblocks, int nthreads, int S,
		char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz);
//*************************************************************************
extern "C" void dvc_PackDist(int grid, int threads, int q, int *SendList, int start,
		int sendCount, double *sendbuf, double *Dist, int N);
//*************************************************************************
extern "C" void dvc_UnpackDist(int grid, int threads, int q, int Cqx, int Cqy, int Cqz, int *RecvList, int start,
		int recvCount, double *recvbuf, double *Dist, int Nx, int Ny, int Nz);
//*************************************************************************
