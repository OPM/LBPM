// CPU Functions for D3Q7 Lattice Boltzmann Methods

extern "C" void dvc_PackValues(int *list, int count, double *sendbuf, double *Data, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = Data[n];
	}
}
extern "C" void dvc_UnpackValues(int *list, int count, double *recvbuf, double *Data, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	for (idx=0; idx<count; idx++){
		n = list[idx];
		Data[n] = recvbuf[idx];
	}
}

extern "C" void dvc_PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N){
	//....................................................................................
	// Pack distribution into the send buffer for the listed lattice sites
	//....................................................................................
	int idx,n,component;
	for (idx=0; idx<count; idx++){
		for (component=0; component<number; component++){
			n = list[idx];
			sendbuf[idx*number+component] = Data[number*n+component];
			Data[number*n+component] = 0.0;	// Set the data value to zero once it's in the buffer!
		}
	}
}


extern  "C" void dvc_UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N){
	//....................................................................................
	// Unack distribution from the recv buffer
	// Sum to the existing density value
	//....................................................................................
	int idx,n,component;
	for (idx=0; idx<count; idx++){
		for (component=0; component<number; component++){
			n = list[idx];
			Data[number*n+component] += recvbuf[idx*number+component];
		}
	}
}

