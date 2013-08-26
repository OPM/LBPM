// 
//*************************************************************************
//*************************************************************************
extern "C" void dvc_InitDenColor( int nblocks, int nthreads, int S,
		char *ID, double *Den, double *Phi, double das, double dbs, int N);
//*************************************************************************
extern "C" void dvc_ComputeColorGradient(int nBlocks, int nthreads, int S,
		char *ID, double *Phi, double *ColorGrad, int Nx, int Ny, int Nz);
//*************************************************************************
extern "C" void dvc_ColorCollide(int nBlocks, int nthreads, int S,
		char *ID, double *f_even, double *f_odd, double *ColorGrad, double *Velocity,
		double rlxA, double rlxB,double alpha, double beta, double Fx, double Fy, double Fz,
		int Nx, int Ny, int Nz, bool pBC);
//*************************************************************************
extern "C" void dvc_DensityStreamD3Q7(int nBlocks, int nthreads, int S,
		char *ID, double *Den, double *Copy, double *Phi, double *ColorGrad, double *Velocity,
		double beta, int Nx, int Ny, int Nz, bool pBC);
//*************************************************************************
extern "C" void dvc_ComputePhi(int nBlocks, int nthreads, int S,
		char *ID, double *Phi, double *Copy, double *Den, int N);
//*************************************************************************
