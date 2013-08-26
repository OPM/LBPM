// GPU Functions for D3Q7 Lattice Boltzmann Methods

__global__ void PackValues(int *list, int count, double *sendbuf, double *Data, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		n = list[idx];
		sendbuf[idx] = Data[n];
	}
}
__global__ void UnpackValues(int *list, int count, double *recvbuf, double *Data, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		n = list[idx];
		Data[n] = recvbuf[idx];
	}
}

__global__ void PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N){
	//....................................................................................
	// Pack distribution into the send buffer for the listed lattice sites
	//....................................................................................
	int idx,n,component;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
//	for (idx=0; idx<count; idx++){
		for (component=0; component<number; component++){
			n = list[idx];
			sendbuf[idx*number+component] = Data[number*n+component];
			Data[number*n+component] = 0.0;	// Set the data value to zero once it's in the buffer!
		}
	}
}


__global__ void UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N){
	//....................................................................................
	// Unack distribution from the recv buffer
	// Sum to the existing density value
	//....................................................................................
	int idx,n,component;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		//	for (idx=0; idx<count; idx++){
		for (component=0; component<number; component++){
			n = list[idx];
			Data[number*n+component] += recvbuf[idx*number+component];
		}
	}
}

//***************************************************************************************
extern "C" void dvc_PackDenD3Q7(int grid, int threads, int *list, int count, double *sendbuf,
		int number, double *Data, int N)
{
	//...................................................................................
	PackDenD3Q7<<<grid,threads>>>(list,count,sendbuf,number,Data,N);
}
//***************************************************************************************
extern "C" void dvc_UnpackDenD3Q7(int grid, int threads, int *list, int count, double *recvbuf,
		int number, double *Data, int N)
{
	//...................................................................................
	UnpackDenD3Q7<<<grid,threads>>>(list,count,recvbuf,number,Data,N);
}
//***************************************************************************************
extern "C" void dvc_PackValues(int grid, int threads, int *list, int count, double *sendbuf,
		double *Data, int N)
{
	//...................................................................................
	PackValues<<<grid,threads>>>(list,count,sendbuf,Data,N);
}
//***************************************************************************************
extern "C" void dvc_UnpackValues(int grid, int threads, int *list, int count, double *recvbuf,
		double *Data, int N)
{
	//...................................................................................
	UnpackValues<<<grid,threads>>>(list,count,recvbuf,Data,N);
}
//***************************************************************************************
