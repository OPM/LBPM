extern "C" void PackDist(int q, int *list, int start, int count, double *sendbuf, double *dist, int N);

extern "C" void UnpackDist(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,double *recvbuf, double *dist, int Nx, int Ny, int Nz);
	       
extern "C" void InitD3Q19(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz);

extern "C" void SwapD3Q19(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz);

extern "C" void MRT(char *ID, double *f_even, double *f_odd, double rlxA, double rlxB, double Fx, double Fy, double Fz,int Nx, int Ny, int Nz);

extern "C" void InitD3Q19(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz);
