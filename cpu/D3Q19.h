
extern void PackDist(int q, int *list, int start, int count, double *sendbuf, double *dist, int N);
extern void MapRecvDist(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
			double *recvbuf, double *dist, int Nx, int Ny, int Nz);
//*************************************************************************
extern void SwapD3Q19(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz);
