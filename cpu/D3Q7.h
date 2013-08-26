// CPU Functions for D3Q7 Lattice Boltzmann Methods

extern void PackValues(int *list, int count, double *sendbuf, double *Data, int N);

extern void UnpackValues(int *list, int count, double *recvbuf, double *Data, int N);

extern void PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N);

extern void UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N);
