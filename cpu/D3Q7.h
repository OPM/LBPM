// CPU Functions for D3Q7 Lattice Boltzmann Methods

extern "C" void dvc_PackValues(int *list, int count, double *sendbuf, double *Data, int N);

extern "C" void dvc_UnpackValues(int *list, int count, double *recvbuf, double *Data, int N);

extern "C" void dvc_PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N);

extern "C" void dvc_UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N);

extern "C" void dvc_InitD3Q7(char *ID, double *f_even, double *f_odd, double *Den, int Nx, int Ny, int Nz, int S);

extern "C" void dvc_SwapD3Q7(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz, int S);

extern "C" void dvc_ComputeDensityD3Q7(char *ID, double *disteven, double *distodd, double *Den, 
										int Nx, int Ny, int Nz, int S);