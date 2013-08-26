// 
//***************************************************************************************
//***************************************************************************************
extern "C" void dvc_PackDenD3Q7(int grid, int threads, int *list, int count, double *sendbuf,
		int number, double *Data, int N);
//***************************************************************************************
extern "C" void dvc_UnpackDenD3Q7(int grid, int threads, int *list, int count, double *recvbuf,
		int number, double *Data, int N);
//***************************************************************************************
extern "C" void dvc_PackValues(int grid, int threads, int *list, int count, double *sendbuf,
		double *Data, int N);
//***************************************************************************************
extern "C" void dvc_UnpackValues(int grid, int threads, int *list, int count, double *recvbuf,
		double *Data, int N);
//***************************************************************************************
