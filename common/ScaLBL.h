/* ScaLBL.h 
 *  Header file for Scalable Lattice Boltzmann Library
 *  Separate implementations for GPU and CPU must both follow the conventions defined in this header
 *  This libarry contains the essential components of the LBM
 *     - streaming implementations
 *     - collision terms to model various physics
 *     - communication framework for the LBM
 *  Refer to Domain.h for setup of parallel domains
 */
#include "Domain.h"

extern "C" void ScaLBL_AllocateDeviceMemory(void** address, size_t size);

extern "C" void ScaLBL_FreeDeviceMemory(void* pointer);

extern "C" void ScaLBL_CopyToDevice(void* dest, const void* source, size_t size);

extern "C" void ScaLBL_CopyToHost(void* dest, const void* source, size_t size);

extern "C" void ScaLBL_DeviceBarrier();

extern "C" void ScaLBL_D3Q19_Pack(int q, int *list, int start, int count, double *sendbuf, double *dist, int N);

extern "C" void ScaLBL_D3Q19_Unpack(int q, int *list, int start, int count, double *recvbuf, double *dist, int N);

extern "C" void ScaLBL_Scalar_Pack(int *list, int count, double *sendbuf, double *Data, int N);

extern "C" void ScaLBL_Scalar_Unpack(int *list, int count, double *recvbuf, double *Data, int N);

extern "C" void ScaLBL_PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N);

extern "C" void ScaLBL_UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N);

extern "C" void ScaLBL_D3Q7_Init(char *ID, double *f_even, double *f_odd, double *Den, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q7_Swap(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q7_Density(char *ID, double *disteven, double *distodd, double *Den,
                                                                                int Nx, int Ny, int Nz);
extern "C" void ScaLBL_D3Q19_Init(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_Swap(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_Swap_Compact(int *neighborList, double *disteven, double *distodd, int Np);

extern "C" void ScaLBL_D3Q19_MRT(char *ID, double *f_even, double *f_odd, double rlxA, double rlxB,
		double Fx, double Fy, double Fz,int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_Velocity(char *ID, double *disteven, double *distodd, double *vel,
		int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_Pressure(const char *ID, const double *disteven, const double *distodd, 
        double *Pressure, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_Pressure_BC_z(double *disteven, double *distodd, double din,
								  int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_Pressure_BC_Z(double *disteven, double *distodd, double dout,
								   int Nx, int Ny, int Nz, int outlet);


extern "C" double ScaLBL_D3Q19_Flux_BC_z(double *disteven, double *distodd, double flux,
								  int Nx, int Ny, int Nz);

extern "C" double ScaLBL_D3Q19_Flux_BC_Z(double *disteven, double *distodd, double flux,
								   int Nx, int Ny, int Nz, int outlet);

extern "C" void ScaLBL_Color_Init(char *ID, double *Den, double *Phi, double das, double dbs, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_ColorDistance_Init(char *ID, double *Den, double *Phi, double *Distance,
								double das, double dbs, double beta, double xp, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_ColorGradient(char *ID, double *phi, double *ColorGrad, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_ColorCollide( char *ID, double *disteven, double *distodd, double *phi, double *ColorGrad,
								double *Velocity, int Nx, int Ny, int Nz,double rlx_setA, double rlx_setB,
								double alpha, double beta, double Fx, double Fy, double Fz);

extern "C" void ScaLBL_D3Q19_ColorCollide_gen( char *ID, double *disteven, double *distodd, double *phi, double *ColorGrad,
								double *Velocity, int Nx, int Ny, int Nz, double tau1, double tau2,
					       double alpha, double beta, double Fx, double Fy, double Fz);

extern "C" void ScaLBL_D3Q7_ColorCollideMass(char *ID, double *A_even, double *A_odd, double *B_even, double *B_odd, 
					     double *Den, double *Phi, double *ColorGrad, double *Velocity, double beta, int N, bool pBC);


extern "C" void ScaLBL_ComputePhaseField(char *ID, double *Phi, double *Den, int N);

extern "C" void ScaLBL_Color_BC_z(double *Phi, double *Den, double *A_even, double *A_odd,
								  double *B_even, double *B_odd, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_Color_BC_Z(double *Phi, double *Den, double *A_even, double *A_odd,
								  double *B_even, double *B_odd, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_Velocity_BC_z(double *disteven, double *distodd, double uz,
								   int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_Velocity_BC_Z(double *disteven, double *distodd, double uz,
								   int Nx, int Ny, int Nz, int outlet);

extern "C" void ScaLBL_SetSlice_z(double *Phi, double value, int Nx, int Ny, int Nz, int Slice);

class ScaLBL_Communicator{
public:
	//......................................................................................
	ScaLBL_Communicator(Domain &Dm);
	//ScaLBL_Communicator(Domain &Dm, IntArray &Map);
	~ScaLBL_Communicator();
	//......................................................................................
	unsigned long int CommunicationCount,SendCount,RecvCount;
	//......................................................................................
	//  Set up for D319 distributions
	// 		- determines how much memory is allocated
	//		- buffers are reused to send D3Q7 distributions and halo exchange as needed
	//......................................................................................
	// Buffers to store data sent and recieved by this MPI process
	double *sendbuf_x, *sendbuf_y, *sendbuf_z, *sendbuf_X, *sendbuf_Y, *sendbuf_Z;
	double *sendbuf_xy, *sendbuf_yz, *sendbuf_xz, *sendbuf_Xy, *sendbuf_Yz, *sendbuf_xZ;
	double *sendbuf_xY, *sendbuf_yZ, *sendbuf_Xz, *sendbuf_XY, *sendbuf_YZ, *sendbuf_XZ;
	double *recvbuf_x, *recvbuf_y, *recvbuf_z, *recvbuf_X, *recvbuf_Y, *recvbuf_Z;
	double *recvbuf_xy, *recvbuf_yz, *recvbuf_xz, *recvbuf_Xy, *recvbuf_Yz, *recvbuf_xZ;
	double *recvbuf_xY, *recvbuf_yZ, *recvbuf_Xz, *recvbuf_XY, *recvbuf_YZ, *recvbuf_XZ;
	//......................................................................................

	void MemoryOptimizedLayout(IntArray &Map, int *neighborList, char *id, int Np);
	void SendD3Q19(double *f_even, double *f_odd);
	void RecvD3Q19(double *f_even, double *f_odd);
	void BiSendD3Q7(double *A_even, double *A_odd, double *B_even, double *B_odd);
	void BiRecvD3Q7(double *A_even, double *A_odd, double *B_even, double *B_odd);
	void SendHalo(double *data);
	void RecvHalo(double *data);

	// Debugging and unit testing functions
	void PrintD3Q19();

private:
	void D3Q19_MapRecv(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count, int *d3q19_recvlist);

	bool Lock; 	// use Lock to make sure only one call at a time to protect data in transit
				// only one set of Send requests can be active at any time (per instance)
	int i,j,k,n;
	int sendtag,recvtag;
	// Give the object it's own MPI communicator
	RankInfoStruct rank_info;
	MPI_Group Group;	// Group of processors associated with this domain
	MPI_Comm MPI_COMM_SCALBL;		// MPI Communicator for this domain
	MPI_Request req1[18],req2[18];
	MPI_Status stat1[18],stat2[18];
	//......................................................................................
	// MPI ranks for all 18 neighbors
	//......................................................................................
	// These variables are all private to prevent external things from modifying them!!
	//......................................................................................
	int rank;
	int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
	int rank_xy,rank_XY,rank_xY,rank_Xy;
	int rank_xz,rank_XZ,rank_xZ,rank_Xz;
	int rank_yz,rank_YZ,rank_yZ,rank_Yz;
	//......................................................................................
	int Nx,Ny,Nz,N;
	//......................................................................................
	int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z;
	int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz, sendCount_xZ;
	int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ, sendCount_XZ;
	//......................................................................................
	int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z;
	int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz, recvCount_xZ;
	int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ, recvCount_XZ;
	//......................................................................................
	// Send buffers that reside on the compute device
	int *dvcSendList_x, *dvcSendList_y, *dvcSendList_z, *dvcSendList_X, *dvcSendList_Y, *dvcSendList_Z;
	int *dvcSendList_xy, *dvcSendList_yz, *dvcSendList_xz, *dvcSendList_Xy, *dvcSendList_Yz, *dvcSendList_xZ;
	int *dvcSendList_xY, *dvcSendList_yZ, *dvcSendList_Xz, *dvcSendList_XY, *dvcSendList_YZ, *dvcSendList_XZ;
	// Recieve buffers that reside on the compute device
	int *dvcRecvList_x, *dvcRecvList_y, *dvcRecvList_z, *dvcRecvList_X, *dvcRecvList_Y, *dvcRecvList_Z;
	int *dvcRecvList_xy, *dvcRecvList_yz, *dvcRecvList_xz, *dvcRecvList_Xy, *dvcRecvList_Yz, *dvcRecvList_xZ;
	int *dvcRecvList_xY, *dvcRecvList_yZ, *dvcRecvList_Xz, *dvcRecvList_XY, *dvcRecvList_YZ, *dvcRecvList_XZ;
	// Recieve buffers for the distributions
	int *dvcRecvDist_x, *dvcRecvDist_y, *dvcRecvDist_z, *dvcRecvDist_X, *dvcRecvDist_Y, *dvcRecvDist_Z;
	int *dvcRecvDist_xy, *dvcRecvDist_yz, *dvcRecvDist_xz, *dvcRecvDist_Xy, *dvcRecvDist_Yz, *dvcRecvDist_xZ;
	int *dvcRecvDist_xY, *dvcRecvDist_yZ, *dvcRecvDist_Xz, *dvcRecvDist_XY, *dvcRecvDist_YZ, *dvcRecvDist_XZ;
	//......................................................................................

};



ScaLBL_Communicator::ScaLBL_Communicator(Domain &Dm){
	//......................................................................................
	Lock=false; // unlock the communicator
	//......................................................................................
	// Create a separate copy of the communicator for the device
	MPI_Comm_group(Dm.Comm,&Group);
	MPI_Comm_create(Dm.Comm,Group,&MPI_COMM_SCALBL);
	//......................................................................................
	// Copy the domain size and communication information directly from Dm
	Nx = Dm.Nx;
	Ny = Dm.Ny;
	Nz = Dm.Nz;
	N = Nx*Ny*Nz;
	rank=Dm.rank;
	rank_x=Dm.rank_x;
	rank_y=Dm.rank_y;
	rank_z=Dm.rank_z;
	rank_X=Dm.rank_X;
	rank_Y=Dm.rank_Y;
	rank_Z=Dm.rank_Z;
	rank_xy=Dm.rank_xy;
	rank_XY=Dm.rank_XY;
	rank_xY=Dm.rank_xY;
	rank_Xy=Dm.rank_Xy;
	rank_xz=Dm.rank_xz;
	rank_XZ=Dm.rank_XZ;
	rank_xZ=Dm.rank_xZ;
	rank_Xz=Dm.rank_Xz;
	rank_yz=Dm.rank_yz;
	rank_YZ=Dm.rank_YZ;
	rank_yZ=Dm.rank_yZ;
	rank_Yz=Dm.rank_Yz;
	sendCount_x=Dm.sendCount_x;
	sendCount_y=Dm.sendCount_y;
	sendCount_z=Dm.sendCount_z;
	sendCount_X=Dm.sendCount_X;
	sendCount_Y=Dm.sendCount_Y;
	sendCount_Z=Dm.sendCount_Z;
	sendCount_xy=Dm.sendCount_xy;
	sendCount_yz=Dm.sendCount_yz;
	sendCount_xz=Dm.sendCount_xz;
	sendCount_Xy=Dm.sendCount_Xy;
	sendCount_Yz=Dm.sendCount_Yz;
	sendCount_xZ=Dm.sendCount_xZ;
	sendCount_xY=Dm.sendCount_xY;
	sendCount_yZ=Dm.sendCount_yZ;
	sendCount_Xz=Dm.sendCount_Xz;
	sendCount_XY=Dm.sendCount_XY;
	sendCount_YZ=Dm.sendCount_YZ;
	sendCount_XZ=Dm.sendCount_XZ;
	recvCount_x=Dm.recvCount_x;
	recvCount_y=Dm.recvCount_y;
	recvCount_z=Dm.recvCount_z;
	recvCount_X=Dm.recvCount_X;
	recvCount_Y=Dm.recvCount_Y;
	recvCount_Z=Dm.recvCount_Z;
	recvCount_xy=Dm.recvCount_xy;
	recvCount_yz=Dm.recvCount_yz;
	recvCount_xz=Dm.recvCount_xz;
	recvCount_Xy=Dm.recvCount_Xy;
	recvCount_Yz=Dm.recvCount_Yz;
	recvCount_xZ=Dm.recvCount_xZ;
	recvCount_xY=Dm.recvCount_xY;
	recvCount_yZ=Dm.recvCount_yZ;
	recvCount_Xz=Dm.recvCount_Xz;
	recvCount_XY=Dm.recvCount_XY;
	recvCount_YZ=Dm.recvCount_YZ;
	recvCount_XZ=Dm.recvCount_XZ;
	//......................................................................................
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_x, 5*sendCount_x*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_X, 5*sendCount_X*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_y, 5*sendCount_y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_Y, 5*sendCount_Y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_z, 5*sendCount_z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_Z, 5*sendCount_Z*sizeof(double));	// Allocatevoid * memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_xy, sendCount_xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_xY, sendCount_xY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_Xy, sendCount_Xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_XY, sendCount_XY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_xz, sendCount_xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_xZ, sendCount_xZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_Xz, sendCount_Xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_XZ, sendCount_XZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_yz, sendCount_yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_yZ, sendCount_yZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_Yz, sendCount_Yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_YZ, sendCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_x, 5*recvCount_x*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_X, 5*recvCount_X*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_y, 5*recvCount_y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_Y, 5*recvCount_Y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_z, 5*recvCount_z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_Z, 5*recvCount_Z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_xy, recvCount_xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_xY, recvCount_xY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_Xy, recvCount_Xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_XY, recvCount_XY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_xz, recvCount_xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_xZ, recvCount_xZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_Xz, recvCount_Xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_XZ, recvCount_XZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_yz, recvCount_yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_yZ, recvCount_yZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_Yz, recvCount_Yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_YZ, recvCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_x, sendCount_x*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_X, sendCount_X*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_y, sendCount_y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_Y, sendCount_Y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_z, sendCount_z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_Z, sendCount_Z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_xy, sendCount_xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_xY, sendCount_xY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_Xy, sendCount_Xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_XY, sendCount_XY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_xz, sendCount_xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_xZ, sendCount_xZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_Xz, sendCount_Xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_XZ, sendCount_XZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_yz, sendCount_yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_yZ, sendCount_yZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_Yz, sendCount_Yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_YZ, sendCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_x, recvCount_x*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_X, recvCount_X*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_y, recvCount_y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_Y, recvCount_Y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_z, recvCount_z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_Z, recvCount_Z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_xy, recvCount_xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_xY, recvCount_xY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_Xy, recvCount_Xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_XY, recvCount_XY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_xz, recvCount_xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_xZ, recvCount_xZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_Xz, recvCount_Xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_XZ, recvCount_XZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_yz, recvCount_yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_yZ, recvCount_yZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_Yz, recvCount_Yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_YZ, recvCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_x, 5*recvCount_x*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_X, 5*recvCount_X*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_y, 5*recvCount_y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_Y, 5*recvCount_Y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_z, 5*recvCount_z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_Z, 5*recvCount_Z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_xy, recvCount_xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_xY, recvCount_xY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_Xy, recvCount_Xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_XY, recvCount_XY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_xz, recvCount_xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_xZ, recvCount_xZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_Xz, recvCount_Xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_XZ, recvCount_XZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_yz, recvCount_yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_yZ, recvCount_yZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_Yz, recvCount_Yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvDist_YZ, recvCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	ScaLBL_CopyToDevice(dvcSendList_x,Dm.sendList_x,sendCount_x*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_X,Dm.sendList_X,sendCount_X*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_y,Dm.sendList_y,sendCount_y*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_Y,Dm.sendList_Y,sendCount_Y*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_z,Dm.sendList_z,sendCount_z*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_Z,Dm.sendList_Z,sendCount_Z*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_xy,Dm.sendList_xy,sendCount_xy*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_XY,Dm.sendList_XY,sendCount_XY*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_xY,Dm.sendList_xY,sendCount_xY*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_Xy,Dm.sendList_Xy,sendCount_Xy*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_xz,Dm.sendList_xz,sendCount_xz*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_XZ,Dm.sendList_XZ,sendCount_XZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_xZ,Dm.sendList_xZ,sendCount_xZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_Xz,Dm.sendList_Xz,sendCount_Xz*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_yz,Dm.sendList_yz,sendCount_yz*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_YZ,Dm.sendList_YZ,sendCount_YZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_yZ,Dm.sendList_yZ,sendCount_yZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_Yz,Dm.sendList_Yz,sendCount_Yz*sizeof(int));
	//......................................................................................
	ScaLBL_CopyToDevice(dvcRecvList_x,Dm.recvList_x,recvCount_x*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_X,Dm.recvList_X,recvCount_X*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_y,Dm.recvList_y,recvCount_y*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_Y,Dm.recvList_Y,recvCount_Y*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_z,Dm.recvList_z,recvCount_z*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_Z,Dm.recvList_Z,recvCount_Z*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_xy,Dm.recvList_xy,recvCount_xy*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_XY,Dm.recvList_XY,recvCount_XY*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_xY,Dm.recvList_xY,recvCount_xY*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_Xy,Dm.recvList_Xy,recvCount_Xy*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_xz,Dm.recvList_xz,recvCount_xz*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_XZ,Dm.recvList_XZ,recvCount_XZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_xZ,Dm.recvList_xZ,recvCount_xZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_Xz,Dm.recvList_Xz,recvCount_Xz*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_yz,Dm.recvList_yz,recvCount_yz*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_YZ,Dm.recvList_YZ,recvCount_YZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_yZ,Dm.recvList_yZ,recvCount_yZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_Yz,Dm.recvList_Yz,recvCount_Yz*sizeof(int));
	//......................................................................................

	MPI_Barrier(MPI_COMM_SCALBL);
	if (rank==0) printf("Mapping recieve distributions \n");

	//...................................................................................
	// Set up the recieve distribution lists
	//...................................................................................
	//...Map recieve list for the X face: q=2,8,10,12,13 .................................
	D3Q19_MapRecv(0,-1,0,0,Dm.recvList_X,0,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(3,-1,-1,0,Dm.recvList_X,recvCount_X,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(4,-1,1,0,Dm.recvList_X,2*recvCount_X,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(5,-1,0,-1,Dm.recvList_X,3*recvCount_X,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(6,-1,0,1,Dm.recvList_X,4*recvCount_X,recvCount_X,dvcRecvDist_X);
	//...................................................................................
	//...Map recieve list for the x face: q=1,7,9,11,13..................................
	D3Q19_MapRecv(1,1,0,0,Dm.recvList_x,0,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(4,1,1,0,Dm.recvList_x,recvCount_x,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(5,1,-1,0,Dm.recvList_x,2*recvCount_x,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(6,1,0,1,Dm.recvList_x,3*recvCount_x,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(7,1,0,-1,Dm.recvList_x,4*recvCount_x,recvCount_x,dvcRecvDist_x);
	//...................................................................................
	//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
	D3Q19_MapRecv(1,0,-1,0,Dm.recvList_Y,0,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(3,-1,-1,0,Dm.recvList_Y,recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(5,1,-1,0,Dm.recvList_Y,2*recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(7,0,-1,-1,Dm.recvList_Y,3*recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(8,0,-1,1,Dm.recvList_Y,4*recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	//...................................................................................
	//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
	D3Q19_MapRecv(2,0,1,0,Dm.recvList_y,0,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(4,1,1,0,Dm.recvList_y,recvCount_y,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(4,-1,1,0,Dm.recvList_y,2*recvCount_y,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(8,0,1,1,Dm.recvList_y,3*recvCount_y,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(9,0,1,-1,Dm.recvList_y,4*recvCount_y,recvCount_y,dvcRecvDist_y);
	//...................................................................................
	//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
	D3Q19_MapRecv(2,0,0,-1,Dm.recvList_Z,0,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(5,-1,0,-1,Dm.recvList_Z,recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(7,1,0,-1,Dm.recvList_Z,2*recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(7,0,-1,-1,Dm.recvList_Z,3*recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(9,0,1,-1,Dm.recvList_Z,4*recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
	D3Q19_MapRecv(3,0,0,1,Dm.recvList_z,0,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(6,1,0,1,Dm.recvList_z,recvCount_z,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(6,-1,0,1,Dm.recvList_z,2*recvCount_z,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(8,0,1,1,Dm.recvList_z,3*recvCount_z,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(8,0,-1,1,Dm.recvList_z,4*recvCount_z,recvCount_z,dvcRecvDist_z);
	//..................................................................................
	//...Map recieve list for the xy edge <<<8)................................
	D3Q19_MapRecv(3,-1,-1,0,Dm.recvList_XY,0,recvCount_XY,dvcRecvDist_XY);
	//...Map recieve list for the Xy edge <<<9)................................
	D3Q19_MapRecv(5,1,-1,0,Dm.recvList_xY,0,recvCount_xY,dvcRecvDist_xY);
	//...Map recieve list for the xY edge <<<10)................................
	D3Q19_MapRecv(4,-1,1,0,Dm.recvList_Xy,0,recvCount_Xy,dvcRecvDist_Xy);
	//...Map recieve list for the XY edge <<<7)................................
	D3Q19_MapRecv(4,1,1,0,Dm.recvList_xy,0,recvCount_xy,dvcRecvDist_xy);
	//...Map recieve list for the xz edge <<<12)................................
	D3Q19_MapRecv(5,-1,0,-1,Dm.recvList_XZ,0,recvCount_XZ,dvcRecvDist_XZ);
	//...Map recieve list for the xZ edge <<<14)................................
	D3Q19_MapRecv(6,-1,0,1,Dm.recvList_Xz,0,recvCount_Xz,dvcRecvDist_Xz);
	//...Map recieve list for the Xz edge <<<13)................................
	D3Q19_MapRecv(7,1,0,-1,Dm.recvList_xZ,0,recvCount_xZ,dvcRecvDist_xZ);
	//...Map recieve list for the XZ edge <<<11)................................
	D3Q19_MapRecv(6,1,0,1,Dm.recvList_xz,0,recvCount_xz,dvcRecvDist_xz);
	//...Map recieve list for the yz edge <<<16)................................
	D3Q19_MapRecv(7,0,-1,-1,Dm.recvList_YZ,0,recvCount_YZ,dvcRecvDist_YZ);
	//...Map recieve list for the yZ edge <<<18)................................
	D3Q19_MapRecv(8,0,-1,1,Dm.recvList_Yz,0,recvCount_Yz,dvcRecvDist_Yz);
	//...Map recieve list for the Yz edge <<<17)................................
	D3Q19_MapRecv(9,0,1,-1,Dm.recvList_yZ,0,recvCount_yZ,dvcRecvDist_yZ);
	//...Map recieve list for the YZ edge <<<15)................................
	D3Q19_MapRecv(8,0,1,1,Dm.recvList_yz,0,recvCount_yz,dvcRecvDist_yz);
	//...................................................................................

	//......................................................................................
	MPI_Barrier(MPI_COMM_SCALBL);
	ScaLBL_DeviceBarrier();
	//......................................................................................
	SendCount = sendCount_x+sendCount_X+sendCount_y+sendCount_Y+sendCount_z+sendCount_Z+
				sendCount_xy+sendCount_Xy+sendCount_xY+sendCount_XY+
				sendCount_xZ+sendCount_Xz+sendCount_xZ+sendCount_XZ+
				sendCount_yz+sendCount_Yz+sendCount_yZ+sendCount_YZ;

	RecvCount = recvCount_x+recvCount_X+recvCount_y+recvCount_Y+recvCount_z+recvCount_Z+
			recvCount_xy+recvCount_Xy+recvCount_xY+recvCount_XY+
			recvCount_xZ+recvCount_Xz+recvCount_xZ+recvCount_XZ+
			recvCount_yz+recvCount_Yz+recvCount_yZ+recvCount_YZ;

	CommunicationCount = SendCount+RecvCount;
	//......................................................................................

}

ScaLBL_Communicator::~ScaLBL_Communicator(){
	// destrutor does nothing (bad idea)
	// -- note that there needs to be a way to free memory allocated on the device!!!
}

void ScaLBL_Communicator::D3Q19_MapRecv(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
					   int *d3q19_recvlist){
	//....................................................................................
	// Map the recieve distributions to
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................

	int i,j,k,n,nn,idx;
	//int N = Nx*Ny*Nz;
	// Create work arrays to store the values from the device (so that host can do this work)
	int * ReturnDist;
	ReturnDist=new int [count];

	// Copy the list from the device
	//ScaLBL_CopyToHost(List, list, count*sizeof(int));

	for (idx=0; idx<count; idx++){
		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = list[idx];
		// Get the 3-D indices
		k = n/(Nx*Ny);
		j = (n-Nx*Ny*k)/Nx;
		i = n-Nx*Ny*k-Nx*j;
		// Streaming for the non-local distribution
		i += Cqx;
		j += Cqy;
		k += Cqz;
		// compute 1D index for the neighbor and save
		nn = k*Nx*Ny+j*Nx+i;
		//d3q19_recvlist[start+idx] = nn;
		ReturnDist[idx] = nn;
	}
	// Return updated version to the device
	ScaLBL_CopyToDevice(&d3q19_recvlist[start], ReturnDist, count*sizeof(int));

	// clean up the work arrays
	delete [] ReturnDist;

}

void ScaLBL_Communicator::MemoryOptimizedLayout(IntArray &Map, int *neighborList, char *id, int Np){
	/*
	 * Generate a memory optimized layout
	 *   id[n] == 0 implies that site n should be ignored (treat as a mask)
	 *   Map(i,j,k) = idx  <- this is the index for the memory optimized layout
	 *   neighborList(idx) <-stores the neighbors for the D3Q19 model
	 *   note that the number of communications remains the same
	 *   the index in the Send and Recv lists is also updated
	 *   this means that the commuincations are no longer valid for regular data structures
	 */
	int idx,i,j,k,n;

	// Check that Map has size matching sub-domain
	if (Map.size(0) != Nx)
		ERROR("ScaLBL_Communicator::MemoryOptimizedLayout: Map array dimensions do not match! \n");

	// Initialize Map
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				Map(i,j,k) = -2;
			}
		}
	}

	idx=0;
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				// domain interior
				Map(i,j,k) = -1;
				// Local index
				n = k*Nx*Ny+j*Nx+i;
				if (id[n] != 0){
					// Counts for the six faces
					if (i==1)       Map(n)=idx++;
					else if (j==1)  Map(n)=idx++;
					else if (k==1)  Map(n)=idx++;
					else if (i==Nx-2)  Map(n)=idx++;
					else if (j==Ny-2)  Map(n)=idx++;
					else if (k==Nz-2)  Map(n)=idx++;
				}
			}
		}
	}

	// Next loop over the domain interior in block-cyclic fashion
	int MemBlockSize=4;
	k=2;
	while (k<Nz-2){
		int kmax = min(k+MemBlockSize,Nz-2);
		j=2;
		while (j<Ny-2){
			int jmax = min(j+MemBlockSize,Ny-2);
			i=2;
			while (i<Nx-2){
				int imax = min(i+MemBlockSize,Nx-2);
				// Loop over the size of the local block
				//printf("Executing block [%i,%i]x[%i,%i] (idx=%i) \n",i,imax,k,kmax,idx);
				for (int kb=k; kb<kmax; kb++){
					for (int jb=j; jb<jmax; jb++){
						for (int ib=i; ib<imax; ib++){
							// Local index (regular layout)
							n = kb*Nx*Ny + jb*Nx + ib;
							if (id[n] != 0 ){
								Map(n) = idx;
								neighborList[idx++] = n; // index of self in regular layout
							}
						}
					}
				}
				i=imax;
			}
			j=jmax;
		}
		k=kmax;
	}


	if (idx > Np ){
		ERROR("ScaLBL_Communicator::MemoryOptimizedLayout: Failed to create memory efficient layout!\n");
	}
/*
	for (k=1;k<Nz-1;k++){
		printf("....k=%i .....\n",k);
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				n=k*Nx*Ny+j*Nx+i;
				idx=Map(i,j,k);
				printf("%i ",idx);
			}
			printf("\n");
		}
	}
*/

	// Now use Map to determine the neighbors for each lattice direction
	for (k=1;k<Nz-1;k++){
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				n=k*Nx*Ny+j*Nx+i;
				idx=Map(i,j,k);
				if (idx > Np) printf("ScaLBL_Communicator::MemoryOptimizedLayout: Map(%i,%i,%i) = %i > %i \n",i,j,k,Map(i,j,k),Np);
				else if (!(idx<0)){
					// store the idx associated with each neighbor
					// store idx for self if neighbor is in solid or out of domain
					//D3Q19 = {{1,0,0},{-1,0,0}
					//         {0,1,0},{0,-1,0}
					//         {0,0,1},{0,0,-1},
					//	       {1,1,0},{-1,-1,0},
					//         {1,-1,0},{-1,1,0},
					//         {1,0,1},{-1,0,-1},
					//         {1,0,-1},{-1,0,1},
					//	       {0,1,1},{0,-1,-1},
					//         {0,1,-1},{0,-1,1}};
					// note that only odd distributions need to be stored to execute the swap algorithm
					int neighbor;    // cycle through the neighbors of lattice site idx
					neighbor=Map(i+1,j,k);
					if (neighbor==-2)		neighborList[idx]=-1;
					else if (neighbor<0)	neighborList[idx]=idx;
					else					neighborList[idx]=neighbor;

					neighbor=Map(i,j+1,k);
					if (neighbor==-2)		neighborList[Np+idx]=-1;
					else if (neighbor<0)	neighborList[Np+idx]=idx;
					else					neighborList[Np+idx]=neighbor;

					neighbor=Map(i,j,k+1);
					if (neighbor==-2)		neighborList[2*Np+idx]=-1;
					else if (neighbor<0)	neighborList[2*Np+idx]=idx;
					else					neighborList[2*Np+idx]=neighbor;

					neighbor=Map(i+1,j+1,k);
					if (neighbor==-2)		neighborList[3*Np+idx]=-1;
					else if (neighbor<0)	neighborList[3*Np+idx]=idx;
					else					neighborList[3*Np+idx]=neighbor;

					neighbor=Map(i+1,j-1,k);
					if (neighbor==-2)		neighborList[4*Np+idx]=-1;
					else if (neighbor<0)	neighborList[4*Np+idx]=idx;
					else					neighborList[4*Np+idx]=neighbor;

					neighbor=Map(i+1,j,k+1);
					if (neighbor==-2)		neighborList[5*Np+idx]=-1;
					else if (neighbor<0)	neighborList[5*Np+idx]=idx;
					else					neighborList[5*Np+idx]=neighbor;

					neighbor=Map(i+1,j,k-1);
					if (neighbor==-2)		neighborList[6*Np+idx]=-1;
					else if (neighbor<0)	neighborList[6*Np+idx]=idx;
					else					neighborList[6*Np+idx]=neighbor;

					neighbor=Map(i,j+1,k+1);
					if (neighbor==-2)		neighborList[7*Np+idx]=-1;
					else if (neighbor<0)	neighborList[7*Np+idx]=idx;
					else					neighborList[7*Np+idx]=neighbor;

					neighbor=Map(i,j+1,k-1);
					if (neighbor==-2)		neighborList[8*Np+idx]=-1;
					else if (neighbor<0)	neighborList[8*Np+idx]=idx;
					else					neighborList[8*Np+idx]=neighbor;
				}
			}
		}
	}

	//for (idx=0; idx<Np; idx++)	printf("%i: %i %i\n", idx, neighborList[3*Np+idx],  neighborList[4*Np+idx]);
	//.......................................................................
	// Now map through  SendList and RecvList to update indices
	// First loop over the send lists

	int *TempBuffer;
	TempBuffer = new int [5*RecvCount];

	//.......................................................................
	// Re-index the send lists
	ScaLBL_CopyToHost(TempBuffer,dvcSendList_x,sendCount_x*sizeof(int));
	for (i=0; i<sendCount_x; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_x,TempBuffer,sendCount_x*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_y,sendCount_y*sizeof(int));
	for (i=0; i<sendCount_y; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_y,TempBuffer,sendCount_y*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_z,sendCount_z*sizeof(int));
	for (i=0; i<sendCount_z; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_z,TempBuffer,sendCount_z*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_X,sendCount_X*sizeof(int));
	for (i=0; i<sendCount_X; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_X,TempBuffer,sendCount_X*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_Y,sendCount_Y*sizeof(int));
	for (i=0; i<sendCount_Y; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_Y,TempBuffer,sendCount_Y*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_Z,sendCount_Z*sizeof(int));
	for (i=0; i<sendCount_Z; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_Z,TempBuffer,sendCount_Z*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_xy,sendCount_xy*sizeof(int));
	for (i=0; i<sendCount_xy; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_xy,TempBuffer,sendCount_xy*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_xY,sendCount_xY*sizeof(int));
	for (i=0; i<sendCount_xY; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_xY,TempBuffer,sendCount_xY*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_Xy,sendCount_Xy*sizeof(int));
	for (i=0; i<sendCount_Xy; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_Xy,TempBuffer,sendCount_Xy*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_XY,sendCount_XY*sizeof(int));
	for (i=0; i<sendCount_XY; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_XY,TempBuffer,sendCount_XY*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_xz,sendCount_xz*sizeof(int));
	for (i=0; i<sendCount_xz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_xz,TempBuffer,sendCount_xz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_xZ,sendCount_xZ*sizeof(int));
	for (i=0; i<sendCount_xZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_xZ,TempBuffer,sendCount_xZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_Xz,sendCount_Xz*sizeof(int));
	for (i=0; i<sendCount_Xz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_Xz,TempBuffer,sendCount_Xz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_XZ,sendCount_XZ*sizeof(int));
	for (i=0; i<sendCount_XZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_XZ,TempBuffer,sendCount_XZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_yz,sendCount_yz*sizeof(int));
	for (i=0; i<sendCount_yz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_yz,TempBuffer,sendCount_yz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_Yz,sendCount_Yz*sizeof(int));
	for (i=0; i<sendCount_Yz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_Yz,TempBuffer,sendCount_Yz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_yZ,sendCount_yZ*sizeof(int));
	for (i=0; i<sendCount_yZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_yZ,TempBuffer,sendCount_yZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_YZ,sendCount_YZ*sizeof(int));
	for (i=0; i<sendCount_YZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_YZ,TempBuffer,sendCount_YZ*sizeof(int));
	//.......................................................................
	// Re-index the recieve lists for the D3Q19 distributions
	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_x,5*recvCount_x*sizeof(int));
	for (i=0; i<5*recvCount_x; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
//		printf("dvcRecvList_x: Map(%i)=%i \n",n,idx);

	}
	ScaLBL_CopyToDevice(dvcRecvDist_x,TempBuffer,5*recvCount_x*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_y,5*recvCount_y*sizeof(int));
	for (i=0; i<5*recvCount_y; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_y,TempBuffer,5*recvCount_y*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_z,5*recvCount_z*sizeof(int));
	for (i=0; i<5*recvCount_z; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_z,TempBuffer,5*recvCount_z*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_X,5*recvCount_X*sizeof(int));
	for (i=0; i<5*recvCount_X; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_X,TempBuffer,5*recvCount_X*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_Y,5*recvCount_Y*sizeof(int));
	for (i=0; i<5*recvCount_Y; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_Y,TempBuffer,5*recvCount_Y*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_Z,5*recvCount_Z*sizeof(int));
	for (i=0; i<5*recvCount_Z; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_Z,TempBuffer,5*recvCount_Z*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_xy,recvCount_xy*sizeof(int));
	for (i=0; i<recvCount_xy; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_xy,TempBuffer,recvCount_xy*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_xY,recvCount_xY*sizeof(int));
	for (i=0; i<recvCount_xY; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_xY,TempBuffer,recvCount_xY*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_Xy,recvCount_Xy*sizeof(int));
	for (i=0; i<recvCount_Xy; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_Xy,TempBuffer,recvCount_Xy*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_XY,recvCount_XY*sizeof(int));
	for (i=0; i<recvCount_XY; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_XY,TempBuffer,recvCount_XY*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_xz,recvCount_xz*sizeof(int));
	for (i=0; i<recvCount_xz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_xz,TempBuffer,recvCount_xz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_xZ,recvCount_xZ*sizeof(int));
	for (i=0; i<recvCount_xZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_xZ,TempBuffer,recvCount_xZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_Xz,recvCount_Xz*sizeof(int));
	for (i=0; i<recvCount_Xz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_Xz,TempBuffer,recvCount_Xz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_XZ,recvCount_XZ*sizeof(int));
	for (i=0; i<recvCount_XZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_XZ,TempBuffer,recvCount_XZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_yz,recvCount_yz*sizeof(int));
	for (i=0; i<recvCount_yz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_yz,TempBuffer,recvCount_yz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_Yz,recvCount_Yz*sizeof(int));
	for (i=0; i<recvCount_Yz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_Yz,TempBuffer,recvCount_Yz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_yZ,recvCount_yZ*sizeof(int));
	for (i=0; i<recvCount_yZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_yZ,TempBuffer,recvCount_yZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_YZ,recvCount_YZ*sizeof(int));
	for (i=0; i<recvCount_YZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_YZ,TempBuffer,recvCount_YZ*sizeof(int));
	//.......................................................................


	// Reset the value of N to match the dense structure
	N = Np;

	// Clean up
	delete [] TempBuffer;
}




void ScaLBL_Communicator::SendD3Q19(double *f_even, double *f_odd){

	if (Lock==true){
		ERROR("ScaLBL Error (SendD3Q19): ScaLBL_Communicator is locked -- did you forget to match Send/Recv calls?");
	}
	else{
		Lock=true;
	}
	// assign tag of 19 to D3Q19 communication
	sendtag = recvtag = 19;
	ScaLBL_DeviceBarrier();
	// Pack the distributions
	ScaLBL_D3Q19_Pack(1,dvcSendList_x,0,sendCount_x,sendbuf_x,f_even,N);
	ScaLBL_D3Q19_Pack(4,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	ScaLBL_D3Q19_Pack(5,dvcSendList_x,2*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	ScaLBL_D3Q19_Pack(6,dvcSendList_x,3*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	ScaLBL_D3Q19_Pack(7,dvcSendList_x,4*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	//...Packing for X face(1,7,9,11,13)................................
	ScaLBL_D3Q19_Pack(0,dvcSendList_X,0,sendCount_X,sendbuf_X,f_odd,N);
	ScaLBL_D3Q19_Pack(3,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	ScaLBL_D3Q19_Pack(4,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	ScaLBL_D3Q19_Pack(5,dvcSendList_X,3*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	ScaLBL_D3Q19_Pack(6,dvcSendList_X,4*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	//...Packing for y face(4,8,9,16,18).................................
	ScaLBL_D3Q19_Pack(2,dvcSendList_y,0,sendCount_y,sendbuf_y,f_even,N);
	ScaLBL_D3Q19_Pack(4,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	ScaLBL_D3Q19_Pack(4,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,f_odd,N);
	ScaLBL_D3Q19_Pack(8,dvcSendList_y,3*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	ScaLBL_D3Q19_Pack(9,dvcSendList_y,4*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	//...Packing for Y face(3,7,10,15,17).................................
	ScaLBL_D3Q19_Pack(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,f_odd,N);
	ScaLBL_D3Q19_Pack(3,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	ScaLBL_D3Q19_Pack(5,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,f_even,N);
	ScaLBL_D3Q19_Pack(7,dvcSendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	ScaLBL_D3Q19_Pack(8,dvcSendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	//...Packing for z face(6,12,13,16,17)................................
	ScaLBL_D3Q19_Pack(3,dvcSendList_z,0,sendCount_z,sendbuf_z,f_even,N);
	ScaLBL_D3Q19_Pack(6,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,f_even,N);
	ScaLBL_D3Q19_Pack(6,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
	ScaLBL_D3Q19_Pack(8,dvcSendList_z,3*sendCount_z,sendCount_z,sendbuf_z,f_even,N);
	ScaLBL_D3Q19_Pack(8,dvcSendList_z,4*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
	//...Packing for Z face(5,11,14,15,18)................................
	ScaLBL_D3Q19_Pack(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,f_odd,N);
	ScaLBL_D3Q19_Pack(5,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
	ScaLBL_D3Q19_Pack(7,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
	ScaLBL_D3Q19_Pack(7,dvcSendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
	ScaLBL_D3Q19_Pack(9,dvcSendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
	//...Pack the xy edge (8)................................
	ScaLBL_D3Q19_Pack(4,dvcSendList_xy,0,sendCount_xy,sendbuf_xy,f_even,N);
	//...Pack the Xy edge (9)................................
	ScaLBL_D3Q19_Pack(4,dvcSendList_Xy,0,sendCount_Xy,sendbuf_Xy,f_odd,N);
	//...Pack the xY edge (10)................................
	ScaLBL_D3Q19_Pack(5,dvcSendList_xY,0,sendCount_xY,sendbuf_xY,f_even,N);
	//...Pack the XY edge (7)................................
	ScaLBL_D3Q19_Pack(3,dvcSendList_XY,0,sendCount_XY,sendbuf_XY,f_odd,N);
	//...Pack the xz edge (12)................................
	ScaLBL_D3Q19_Pack(6,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
	//...Pack the xZ edge (14)................................
	ScaLBL_D3Q19_Pack(7,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
	//...Pack the Xz edge (13)................................
	ScaLBL_D3Q19_Pack(6,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
	//...Pack the XZ edge (11)................................
	ScaLBL_D3Q19_Pack(5,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
	//...Pack the xz edge (12)................................
	//...Pack the yz edge (16)................................
	ScaLBL_D3Q19_Pack(8,dvcSendList_yz,0,sendCount_yz,sendbuf_yz,f_even,N);
	//...Pack the yZ edge (18)................................
	ScaLBL_D3Q19_Pack(9,dvcSendList_yZ,0,sendCount_yZ,sendbuf_yZ,f_even,N);
	//...Pack the Yz edge (17)................................
	ScaLBL_D3Q19_Pack(8,dvcSendList_Yz,0,sendCount_Yz,sendbuf_Yz,f_odd,N);
	//...Pack the YZ edge (15)................................
	ScaLBL_D3Q19_Pack(7,dvcSendList_YZ,0,sendCount_YZ,sendbuf_YZ,f_odd,N);
	//...................................................................................

	//...................................................................................
	// Send all the distributions
	MPI_Isend(sendbuf_x, 5*sendCount_x,MPI_DOUBLE,rank_x,sendtag,MPI_COMM_SCALBL,&req1[0]);
	MPI_Irecv(recvbuf_X, 5*recvCount_X,MPI_DOUBLE,rank_X,recvtag,MPI_COMM_SCALBL,&req2[0]);
	MPI_Isend(sendbuf_X, 5*sendCount_X,MPI_DOUBLE,rank_X,sendtag,MPI_COMM_SCALBL,&req1[1]);
	MPI_Irecv(recvbuf_x, 5*recvCount_x,MPI_DOUBLE,rank_x,recvtag,MPI_COMM_SCALBL,&req2[1]);
	MPI_Isend(sendbuf_y, 5*sendCount_y,MPI_DOUBLE,rank_y,sendtag,MPI_COMM_SCALBL,&req1[2]);
	MPI_Irecv(recvbuf_Y, 5*recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,MPI_COMM_SCALBL,&req2[2]);
	MPI_Isend(sendbuf_Y, 5*sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,MPI_COMM_SCALBL,&req1[3]);
	MPI_Irecv(recvbuf_y, 5*recvCount_y,MPI_DOUBLE,rank_y,recvtag,MPI_COMM_SCALBL,&req2[3]);
	MPI_Isend(sendbuf_z, 5*sendCount_z,MPI_DOUBLE,rank_z,sendtag,MPI_COMM_SCALBL,&req1[4]);
	MPI_Irecv(recvbuf_Z, 5*recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,MPI_COMM_SCALBL,&req2[4]);
	MPI_Isend(sendbuf_Z, 5*sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,MPI_COMM_SCALBL,&req1[5]);
	MPI_Irecv(recvbuf_z, 5*recvCount_z,MPI_DOUBLE,rank_z,recvtag,MPI_COMM_SCALBL,&req2[5]);
	MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_xy,sendtag,MPI_COMM_SCALBL,&req1[6]);
	MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_XY,recvtag,MPI_COMM_SCALBL,&req2[6]);
	MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_XY,sendtag,MPI_COMM_SCALBL,&req1[7]);
	MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_xy,recvtag,MPI_COMM_SCALBL,&req2[7]);
	MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_Xy,sendtag,MPI_COMM_SCALBL,&req1[8]);
	MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_xY,recvtag,MPI_COMM_SCALBL,&req2[8]);
	MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_xY,sendtag,MPI_COMM_SCALBL,&req1[9]);
	MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_Xy,recvtag,MPI_COMM_SCALBL,&req2[9]);
	MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_xz,sendtag,MPI_COMM_SCALBL,&req1[10]);
	MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_XZ,recvtag,MPI_COMM_SCALBL,&req2[10]);
	MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_XZ,sendtag,MPI_COMM_SCALBL,&req1[11]);
	MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_xz,recvtag,MPI_COMM_SCALBL,&req2[11]);
	MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_Xz,sendtag,MPI_COMM_SCALBL,&req1[12]);
	MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_xZ,recvtag,MPI_COMM_SCALBL,&req2[12]);
	MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_xZ,sendtag,MPI_COMM_SCALBL,&req1[13]);
	MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_Xz,recvtag,MPI_COMM_SCALBL,&req2[13]);
	MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_yz,sendtag,MPI_COMM_SCALBL,&req1[14]);
	MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_YZ,recvtag,MPI_COMM_SCALBL,&req2[14]);
	MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_YZ,sendtag,MPI_COMM_SCALBL,&req1[15]);
	MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_yz,recvtag,MPI_COMM_SCALBL,&req2[15]);
	MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_Yz,sendtag,MPI_COMM_SCALBL,&req1[16]);
	MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_yZ,recvtag,MPI_COMM_SCALBL,&req2[16]);
	MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_yZ,sendtag,MPI_COMM_SCALBL,&req1[17]);
	MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_Yz,recvtag,MPI_COMM_SCALBL,&req2[17]);
}

void ScaLBL_Communicator::RecvD3Q19(double *f_even, double *f_odd){
	//...................................................................................
	// Wait for completion of D3Q19 communication
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	ScaLBL_DeviceBarrier();

	//...................................................................................
	// Unpack the distributions on the device
	//...................................................................................
	//...Map recieve list for the X face: q=2,8,10,12,13 .................................
	ScaLBL_D3Q19_Unpack(0,dvcRecvDist_X,0,recvCount_X,recvbuf_X,f_odd,N);
	ScaLBL_D3Q19_Unpack(3,dvcRecvDist_X,recvCount_X,recvCount_X,recvbuf_X,f_odd,N);
	ScaLBL_D3Q19_Unpack(4,dvcRecvDist_X,2*recvCount_X,recvCount_X,recvbuf_X,f_odd,N);
	ScaLBL_D3Q19_Unpack(5,dvcRecvDist_X,3*recvCount_X,recvCount_X,recvbuf_X,f_odd,N);
	ScaLBL_D3Q19_Unpack(6,dvcRecvDist_X,4*recvCount_X,recvCount_X,recvbuf_X,f_odd,N);
	//...................................................................................
	//...Map recieve list for the x face: q=1,7,9,11,13..................................
	ScaLBL_D3Q19_Unpack(1,dvcRecvDist_x,0,recvCount_x,recvbuf_x,f_even,N);
	ScaLBL_D3Q19_Unpack(4,dvcRecvDist_x,recvCount_x,recvCount_x,recvbuf_x,f_even,N);
	ScaLBL_D3Q19_Unpack(5,dvcRecvDist_x,2*recvCount_x,recvCount_x,recvbuf_x,f_even,N);
	ScaLBL_D3Q19_Unpack(6,dvcRecvDist_x,3*recvCount_x,recvCount_x,recvbuf_x,f_even,N);
	ScaLBL_D3Q19_Unpack(7,dvcRecvDist_x,4*recvCount_x,recvCount_x,recvbuf_x,f_even,N);
	//...................................................................................
	//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
	ScaLBL_D3Q19_Unpack(1,dvcRecvDist_Y,0,recvCount_Y,recvbuf_Y,f_odd,N);
	ScaLBL_D3Q19_Unpack(3,dvcRecvDist_Y,recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,N);
	ScaLBL_D3Q19_Unpack(5,dvcRecvDist_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,f_even,N);
	ScaLBL_D3Q19_Unpack(7,dvcRecvDist_Y,3*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,N);
	ScaLBL_D3Q19_Unpack(8,dvcRecvDist_Y,4*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,N);
	//...................................................................................
	//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
	ScaLBL_D3Q19_Unpack(2,dvcRecvDist_y,0,recvCount_y,recvbuf_y,f_even,N);
	ScaLBL_D3Q19_Unpack(4,dvcRecvDist_y,recvCount_y,recvCount_y,recvbuf_y,f_even,N);
	ScaLBL_D3Q19_Unpack(4,dvcRecvDist_y,2*recvCount_y,recvCount_y,recvbuf_y,f_odd,N);
	ScaLBL_D3Q19_Unpack(8,dvcRecvDist_y,3*recvCount_y,recvCount_y,recvbuf_y,f_even,N);
	ScaLBL_D3Q19_Unpack(9,dvcRecvDist_y,4*recvCount_y,recvCount_y,recvbuf_y,f_even,N);
	//...................................................................................
	//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
	ScaLBL_D3Q19_Unpack(2,dvcRecvDist_Z,0,recvCount_Z,recvbuf_Z,f_odd,N);
	ScaLBL_D3Q19_Unpack(5,dvcRecvDist_Z,recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,N);
	ScaLBL_D3Q19_Unpack(7,dvcRecvDist_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,N);
	ScaLBL_D3Q19_Unpack(7,dvcRecvDist_Z,3*recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,N);
	ScaLBL_D3Q19_Unpack(9,dvcRecvDist_Z,4*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,N);
	//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
	ScaLBL_D3Q19_Unpack(3,dvcRecvDist_z,0,recvCount_z,recvbuf_z,f_even,N);
	ScaLBL_D3Q19_Unpack(6,dvcRecvDist_z,recvCount_z,recvCount_z,recvbuf_z,f_even,N);
	ScaLBL_D3Q19_Unpack(6,dvcRecvDist_z,2*recvCount_z,recvCount_z,recvbuf_z,f_odd,N);
	ScaLBL_D3Q19_Unpack(8,dvcRecvDist_z,3*recvCount_z,recvCount_z,recvbuf_z,f_even,N);
	ScaLBL_D3Q19_Unpack(8,dvcRecvDist_z,4*recvCount_z,recvCount_z,recvbuf_z,f_odd,N);
	//..................................................................................
	//...Map recieve list for the xy edge <<<8)................................
	ScaLBL_D3Q19_Unpack(3,dvcRecvDist_XY,0,recvCount_XY,recvbuf_XY,f_odd,N);
	//...Map recieve list for the Xy edge <<<9)................................
	ScaLBL_D3Q19_Unpack(5,dvcRecvDist_xY,0,recvCount_xY,recvbuf_xY,f_even,N);
	//...Map recieve list for the xY edge <<<10)................................
	ScaLBL_D3Q19_Unpack(4,dvcRecvDist_Xy,0,recvCount_Xy,recvbuf_Xy,f_odd,N);
	//...Map recieve list for the XY edge <<<7)................................
	ScaLBL_D3Q19_Unpack(4,dvcRecvDist_xy,0,recvCount_xy,recvbuf_xy,f_even,N);
	//...Map recieve list for the xz edge <<<12)................................
	ScaLBL_D3Q19_Unpack(5,dvcRecvDist_XZ,0,recvCount_XZ,recvbuf_XZ,f_odd,N);
	//...Map recieve list for the xZ edge <<<14)................................
	ScaLBL_D3Q19_Unpack(6,dvcRecvDist_Xz,0,recvCount_Xz,recvbuf_Xz,f_odd,N);
	//...Map recieve list for the Xz edge <<<13)................................
	ScaLBL_D3Q19_Unpack(7,dvcRecvDist_xZ,0,recvCount_xZ,recvbuf_xZ,f_even,N);
	//...Map recieve list for the XZ edge <<<11)................................
	ScaLBL_D3Q19_Unpack(6,dvcRecvDist_xz,0,recvCount_xz,recvbuf_xz,f_even,N);
	//...Map recieve list for the yz edge <<<16)................................
	ScaLBL_D3Q19_Unpack(7,dvcRecvDist_YZ,0,recvCount_YZ,recvbuf_YZ,f_odd,N);
	//...Map recieve list for the yZ edge <<<18)................................
	ScaLBL_D3Q19_Unpack(8,dvcRecvDist_Yz,0,recvCount_Yz,recvbuf_Yz,f_odd,N);
	//...Map recieve list for the Yz edge <<<17)................................
	ScaLBL_D3Q19_Unpack(9,dvcRecvDist_yZ,0,recvCount_yZ,recvbuf_yZ,f_even,N);
	//...Map recieve list for the YZ edge <<<15)................................
	ScaLBL_D3Q19_Unpack(8,dvcRecvDist_yz,0,recvCount_yz,recvbuf_yz,f_even,N);
	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................

}
void ScaLBL_Communicator::BiSendD3Q7(double *A_even, double *A_odd, double *B_even, double *B_odd){

	if (Lock==true){
		ERROR("ScaLBL Error (SendD3Q7): ScaLBL_Communicator is locked -- did you forget to match Send/Recv calls?");
	}
	else{
		Lock=true;
	}
	ScaLBL_DeviceBarrier();
	//...................................................................................
	sendtag = recvtag = 7;
	//...................................................................................
	ScaLBL_D3Q19_Pack(1,dvcSendList_x,0,sendCount_x,sendbuf_x,A_even,N);
	ScaLBL_D3Q19_Pack(1,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,B_even,N);
	//...Packing for X face(1,7,9,11,13)................................
	ScaLBL_D3Q19_Pack(0,dvcSendList_X,0,sendCount_X,sendbuf_X,A_odd,N);
	ScaLBL_D3Q19_Pack(0,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,B_odd,N);
	//...Packing for y face(4,8,9,16,18).................................
	ScaLBL_D3Q19_Pack(2,dvcSendList_y,0,sendCount_y,sendbuf_y,A_even,N);
	ScaLBL_D3Q19_Pack(2,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,B_even,N);
	//...Packing for Y face(3,7,10,15,17).................................
	ScaLBL_D3Q19_Pack(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,A_odd,N);
	ScaLBL_D3Q19_Pack(1,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,B_odd,N);
	//...Packing for z face(6,12,13,16,17)................................
	ScaLBL_D3Q19_Pack(3,dvcSendList_z,0,sendCount_z,sendbuf_z,A_even,N);
	ScaLBL_D3Q19_Pack(3,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,B_even,N);
	//...Packing for Z face(5,11,14,15,18)................................
	ScaLBL_D3Q19_Pack(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,A_odd,N);
	ScaLBL_D3Q19_Pack(2,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,B_odd,N);
	//...................................................................................

	//...................................................................................
	// Send all the distributions
	MPI_Isend(sendbuf_x, 2*sendCount_x,MPI_DOUBLE,rank_x,sendtag,MPI_COMM_SCALBL,&req1[0]);
	MPI_Irecv(recvbuf_X, 2*recvCount_X,MPI_DOUBLE,rank_X,recvtag,MPI_COMM_SCALBL,&req2[0]);
	MPI_Isend(sendbuf_X, 2*sendCount_X,MPI_DOUBLE,rank_X,sendtag,MPI_COMM_SCALBL,&req1[1]);
	MPI_Irecv(recvbuf_x, 2*recvCount_x,MPI_DOUBLE,rank_x,recvtag,MPI_COMM_SCALBL,&req2[1]);
	MPI_Isend(sendbuf_y, 2*sendCount_y,MPI_DOUBLE,rank_y,sendtag,MPI_COMM_SCALBL,&req1[2]);
	MPI_Irecv(recvbuf_Y, 2*recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,MPI_COMM_SCALBL,&req2[2]);
	MPI_Isend(sendbuf_Y, 2*sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,MPI_COMM_SCALBL,&req1[3]);
	MPI_Irecv(recvbuf_y, 2*recvCount_y,MPI_DOUBLE,rank_y,recvtag,MPI_COMM_SCALBL,&req2[3]);
	MPI_Isend(sendbuf_z, 2*sendCount_z,MPI_DOUBLE,rank_z,sendtag,MPI_COMM_SCALBL,&req1[4]);
	MPI_Irecv(recvbuf_Z, 2*recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,MPI_COMM_SCALBL,&req2[4]);
	MPI_Isend(sendbuf_Z, 2*sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,MPI_COMM_SCALBL,&req1[5]);
	MPI_Irecv(recvbuf_z, 2*recvCount_z,MPI_DOUBLE,rank_z,recvtag,MPI_COMM_SCALBL,&req2[5]);
	//...................................................................................
}
void ScaLBL_Communicator::BiRecvD3Q7(double *A_even, double *A_odd, double *B_even, double *B_odd){
	//...................................................................................
	// Wait for completion of D3Q19 communication
	MPI_Waitall(6,req1,stat1);
	MPI_Waitall(6,req2,stat2);
	ScaLBL_DeviceBarrier();
	//...................................................................................
	// Unpack the distributions on the device
/*	//...................................................................................
	//...Map recieve list for the X face: q=2,8,10,12,13 ................................
	ScaLBL_D3Q19_Unpack(0,-1,0,0,dvcRecvList_X,0,recvCount_X,recvbuf_X,A_odd,Nx,Ny,Nz);
	ScaLBL_D3Q19_Unpack(0,-1,0,0,dvcRecvList_X,recvCount_X,recvCount_X,recvbuf_X,B_odd,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the x face: q=1,7,9,11,13..................................
	ScaLBL_D3Q19_Unpack(1,1,0,0,dvcRecvList_x,0,recvCount_x,recvbuf_x,A_even,Nx,Ny,Nz);
	ScaLBL_D3Q19_Unpack(1,1,0,0,dvcRecvList_x,recvCount_x,recvCount_x,recvbuf_x,B_even,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the y face: q=4,8,9,16,18 .................................
	ScaLBL_D3Q19_Unpack(1,0,-1,0,dvcRecvList_Y,0,recvCount_Y,recvbuf_Y,A_odd,Nx,Ny,Nz);
	ScaLBL_D3Q19_Unpack(1,0,-1,0,dvcRecvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,B_odd,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the Y face: q=3,7,10,15,17 ................................
	ScaLBL_D3Q19_Unpack(2,0,1,0,dvcRecvList_y,0,recvCount_y,recvbuf_y,A_even,Nx,Ny,Nz);
	ScaLBL_D3Q19_Unpack(2,0,1,0,dvcRecvList_y,recvCount_y,recvCount_y,recvbuf_y,B_even,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the z face<<<6,12,13,16,17)................................
	ScaLBL_D3Q19_Unpack(2,0,0,-1,dvcRecvList_Z,0,recvCount_Z,recvbuf_Z,A_odd,Nx,Ny,Nz);
	ScaLBL_D3Q19_Unpack(2,0,0,-1,dvcRecvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,B_odd,Nx,Ny,Nz);
	//...Map recieve list for the Z face<<<5,11,14,15,18)................................
	ScaLBL_D3Q19_Unpack(3,0,0,1,dvcRecvList_z,0,recvCount_z,recvbuf_z,A_even,Nx,Ny,Nz);
	ScaLBL_D3Q19_Unpack(3,0,0,1,dvcRecvList_z,recvCount_z,recvCount_z,recvbuf_z,B_even,Nx,Ny,Nz);
	//..................................................................................
	*/
	//...Map recieve list for the X face: q=2,8,10,12,13 ................................
	ScaLBL_D3Q19_Unpack(0,dvcRecvDist_X,0,recvCount_X,recvbuf_X,A_odd,N);
	ScaLBL_D3Q19_Unpack(0,dvcRecvDist_X,recvCount_X,recvCount_X,recvbuf_X,B_odd,N);
	//...................................................................................
	//...Map recieve list for the x face: q=1,7,9,11,13..................................
	ScaLBL_D3Q19_Unpack(1,dvcRecvDist_x,0,recvCount_x,recvbuf_x,A_even,N);
	ScaLBL_D3Q19_Unpack(1,dvcRecvDist_x,recvCount_x,recvCount_x,recvbuf_x,B_even,N);
	//...................................................................................
	//...Map recieve list for the y face: q=4,8,9,16,18 .................................
	ScaLBL_D3Q19_Unpack(1,dvcRecvDist_Y,0,recvCount_Y,recvbuf_Y,A_odd,N);
	ScaLBL_D3Q19_Unpack(1,dvcRecvDist_Y,recvCount_Y,recvCount_Y,recvbuf_Y,B_odd,N);
	//...................................................................................
	//...Map recieve list for the Y face: q=3,7,10,15,17 ................................
	ScaLBL_D3Q19_Unpack(2,dvcRecvDist_y,0,recvCount_y,recvbuf_y,A_even,N);
	ScaLBL_D3Q19_Unpack(2,dvcRecvDist_y,recvCount_y,recvCount_y,recvbuf_y,B_even,N);
	//...................................................................................
	//...Map recieve list for the z face<<<6,12,13,16,17)................................
	ScaLBL_D3Q19_Unpack(2,dvcRecvDist_Z,0,recvCount_Z,recvbuf_Z,A_odd,N);
	ScaLBL_D3Q19_Unpack(2,dvcRecvDist_Z,recvCount_Z,recvCount_Z,recvbuf_Z,B_odd,N);
	//...Map recieve list for the Z face<<<5,11,14,15,18)................................
	ScaLBL_D3Q19_Unpack(3,dvcRecvDist_z,0,recvCount_z,recvbuf_z,A_even,N);
	ScaLBL_D3Q19_Unpack(3,dvcRecvDist_z,recvCount_z,recvCount_z,recvbuf_z,B_even,N);
	//..................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................
}
void ScaLBL_Communicator::SendHalo(double *data){
	//...................................................................................
	if (Lock==true){
		ERROR("ScaLBL Error (SendHalo): ScaLBL_Communicator is locked -- did you forget to match Send/Recv calls?");
	}
	else{
		Lock=true;
	}
	ScaLBL_DeviceBarrier();
	//...................................................................................
	sendtag = recvtag = 1;
	//...................................................................................
	ScaLBL_Scalar_Pack(dvcSendList_x, sendCount_x,sendbuf_x, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_y, sendCount_y,sendbuf_y, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_z, sendCount_z,sendbuf_z, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_X, sendCount_X,sendbuf_X, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_Y, sendCount_Y,sendbuf_Y, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_Z, sendCount_Z,sendbuf_Z, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_xy, sendCount_xy,sendbuf_xy, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_xY, sendCount_xY,sendbuf_xY, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_Xy, sendCount_Xy,sendbuf_Xy, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_XY, sendCount_XY,sendbuf_XY, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_xz, sendCount_xz,sendbuf_xz, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_xZ, sendCount_xZ,sendbuf_xZ, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_Xz, sendCount_Xz,sendbuf_Xz, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_XZ, sendCount_XZ,sendbuf_XZ, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_yz, sendCount_yz,sendbuf_yz, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_yZ, sendCount_yZ,sendbuf_yZ, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_Yz, sendCount_Yz,sendbuf_Yz, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_YZ, sendCount_YZ,sendbuf_YZ, data, N);
	//...................................................................................
	// Send / Recv all the phase indcator field values
	//...................................................................................
	MPI_Isend(sendbuf_x, sendCount_x,MPI_DOUBLE,rank_x,sendtag,MPI_COMM_SCALBL,&req1[0]);
	MPI_Irecv(recvbuf_X, recvCount_X,MPI_DOUBLE,rank_X,recvtag,MPI_COMM_SCALBL,&req2[0]);
	MPI_Isend(sendbuf_X, sendCount_X,MPI_DOUBLE,rank_X,sendtag,MPI_COMM_SCALBL,&req1[1]);
	MPI_Irecv(recvbuf_x, recvCount_x,MPI_DOUBLE,rank_x,recvtag,MPI_COMM_SCALBL,&req2[1]);
	MPI_Isend(sendbuf_y, sendCount_y,MPI_DOUBLE,rank_y,sendtag,MPI_COMM_SCALBL,&req1[2]);
	MPI_Irecv(recvbuf_Y, recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,MPI_COMM_SCALBL,&req2[2]);
	MPI_Isend(sendbuf_Y, sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,MPI_COMM_SCALBL,&req1[3]);
	MPI_Irecv(recvbuf_y, recvCount_y,MPI_DOUBLE,rank_y,recvtag,MPI_COMM_SCALBL,&req2[3]);
	MPI_Isend(sendbuf_z, sendCount_z,MPI_DOUBLE,rank_z,sendtag,MPI_COMM_SCALBL,&req1[4]);
	MPI_Irecv(recvbuf_Z, recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,MPI_COMM_SCALBL,&req2[4]);
	MPI_Isend(sendbuf_Z, sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,MPI_COMM_SCALBL,&req1[5]);
	MPI_Irecv(recvbuf_z, recvCount_z,MPI_DOUBLE,rank_z,recvtag,MPI_COMM_SCALBL,&req2[5]);
	MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_xy,sendtag,MPI_COMM_SCALBL,&req1[6]);
	MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_XY,recvtag,MPI_COMM_SCALBL,&req2[6]);
	MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_XY,sendtag,MPI_COMM_SCALBL,&req1[7]);
	MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_xy,recvtag,MPI_COMM_SCALBL,&req2[7]);
	MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_Xy,sendtag,MPI_COMM_SCALBL,&req1[8]);
	MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_xY,recvtag,MPI_COMM_SCALBL,&req2[8]);
	MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_xY,sendtag,MPI_COMM_SCALBL,&req1[9]);
	MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_Xy,recvtag,MPI_COMM_SCALBL,&req2[9]);
	MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_xz,sendtag,MPI_COMM_SCALBL,&req1[10]);
	MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_XZ,recvtag,MPI_COMM_SCALBL,&req2[10]);
	MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_XZ,sendtag,MPI_COMM_SCALBL,&req1[11]);
	MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_xz,recvtag,MPI_COMM_SCALBL,&req2[11]);
	MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_Xz,sendtag,MPI_COMM_SCALBL,&req1[12]);
	MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_xZ,recvtag,MPI_COMM_SCALBL,&req2[12]);
	MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_xZ,sendtag,MPI_COMM_SCALBL,&req1[13]);
	MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_Xz,recvtag,MPI_COMM_SCALBL,&req2[13]);
	MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_yz,sendtag,MPI_COMM_SCALBL,&req1[14]);
	MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_YZ,recvtag,MPI_COMM_SCALBL,&req2[14]);
	MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_YZ,sendtag,MPI_COMM_SCALBL,&req1[15]);
	MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_yz,recvtag,MPI_COMM_SCALBL,&req2[15]);
	MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_Yz,sendtag,MPI_COMM_SCALBL,&req1[16]);
	MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_yZ,recvtag,MPI_COMM_SCALBL,&req2[16]);
	MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_yZ,sendtag,MPI_COMM_SCALBL,&req1[17]);
	MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_Yz,recvtag,MPI_COMM_SCALBL,&req2[17]);
	//...................................................................................
}
void ScaLBL_Communicator::RecvHalo(double *data){

	//...................................................................................
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	ScaLBL_DeviceBarrier();
	//...................................................................................
	//...................................................................................
	ScaLBL_Scalar_Unpack(dvcRecvList_x, recvCount_x,recvbuf_x, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_y, recvCount_y,recvbuf_y, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_z, recvCount_z,recvbuf_z, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_X, recvCount_X,recvbuf_X, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Y, recvCount_Y,recvbuf_Y, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Z, recvCount_Z,recvbuf_Z, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xy, recvCount_xy,recvbuf_xy, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xY, recvCount_xY,recvbuf_xY, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Xy, recvCount_Xy,recvbuf_Xy, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_XY, recvCount_XY,recvbuf_XY, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xz, recvCount_xz,recvbuf_xz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xZ, recvCount_xZ,recvbuf_xZ, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Xz, recvCount_Xz,recvbuf_Xz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_XZ, recvCount_XZ,recvbuf_XZ, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_yz, recvCount_yz,recvbuf_yz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_yZ, recvCount_yZ,recvbuf_yZ, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Yz, recvCount_Yz,recvbuf_Yz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_YZ, recvCount_YZ,recvbuf_YZ, data, N);
	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................
}

void ScaLBL_Communicator::PrintD3Q19(){
	printf("Printing D3Q19 communication buffer contents \n");

	int i,n;
	double f;
	int *TempBuffer;
	TempBuffer = new int [5*recvCount_x];

	//.......................................................................
	// Re-index the send lists
	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_x,5*recvCount_x*sizeof(int));
	for (i=0; i<5*recvCount_x; i++){
		n = TempBuffer[i];
		f = recvbuf_x[i];
		printf("Receive %f to %i from buffer index %i \n", f, n, i);
	}

	delete [] TempBuffer;
}
