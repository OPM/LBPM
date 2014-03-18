// CPU Functions for D3Q7 Lattice Boltzmann Methods

extern "C" void PackValues(int *list, int count, double *sendbuf, double *Data, int N);

extern "C" void UnpackValues(int *list, int count, double *recvbuf, double *Data, int N);

extern "C" void PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N);

extern "C" void UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N);

extern "C" void InitD3Q7(char *ID, double *f_even, double *f_odd, double *Den, int Nx, int Ny, int Nz);

extern "C" void SwapD3Q7(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz);

extern "C" void ComputeDensityD3Q7(char *ID, double *disteven, double *distodd, double *Den,
										int Nx, int Ny, int Nz);
