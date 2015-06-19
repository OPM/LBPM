extern "C" void PackDist(int q, int *list, int start, int count, double *sendbuf, double *dist, int N);

extern "C" void UnpackDist(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,double *recvbuf, double *dist, int Nx, int Ny, int Nz);
	       
extern "C" void InitD3Q19(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz);

extern "C" void SwapD3Q19(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz);

extern "C" void MRT(char *ID, double *f_even, double *f_odd, double rlxA, double rlxB, double Fx, double Fy, double Fz,int Nx, int Ny, int Nz);

extern "C" void ComputeVelocityD3Q19(char *ID, double *disteven, double *distodd, double *vel,
				int Nx, int Ny, int Nz);

extern "C" void ComputePressureD3Q19(char *ID, double *disteven, double *distodd, double *Pressure,
									int Nx, int Ny, int Nz);

extern "C" void PressureBC_inlet(double *disteven, double *distodd, double din,
								  int Nx, int Ny, int Nz);

extern "C" void PressureBC_outlet(double *disteven, double *distodd, double dout,
								   int Nx, int Ny, int Nz, int outlet);
