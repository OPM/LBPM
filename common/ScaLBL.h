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

extern "C" void AllocateDeviceMemory(void** address, size_t size);

//extern "C" void FreeDeviceMemory(void** address);

extern "C" void CopyToDevice(void* dest, void* source, size_t size);

extern "C" void CopyToHost(void* dest, void* source, size_t size);

extern "C" void DeviceBarrier();

extern "C" void PackDist(int q, int *list, int start, int count, double *sendbuf, double *dist, int N);

extern "C" void UnpackDist(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,double *recvbuf, double *dist, int Nx, int Ny, int Nz);

extern "C" void PackValues(int *list, int count, double *sendbuf, double *Data, int N);

extern "C" void UnpackValues(int *list, int count, double *recvbuf, double *Data, int N);

extern "C" void PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N);

extern "C" void UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N);

extern "C" void InitD3Q7(char *ID, double *f_even, double *f_odd, double *Den, int Nx, int Ny, int Nz);

extern "C" void SwapD3Q7(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz);

extern "C" void ComputeDensityD3Q7(char *ID, double *disteven, double *distodd, double *Den,
                                                                                int Nx, int Ny, int Nz);
extern "C" void InitD3Q19(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz);

extern "C" void SwapD3Q19(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz);

extern "C" void MRT(char *ID, double *f_even, double *f_odd, double rlxA, double rlxB,
		double Fx, double Fy, double Fz,int Nx, int Ny, int Nz);

extern "C" void ComputeVelocityD3Q19(char *ID, double *disteven, double *distodd, double *vel,
		int Nx, int Ny, int Nz);

extern "C" void ComputePressureD3Q19(char *ID, double *disteven, double *distodd, double *Pressure,
									int Nx, int Ny, int Nz);

extern "C" void PressureBC_inlet(double *disteven, double *distodd, double din,
								  int Nx, int Ny, int Nz);

extern "C" void PressureBC_outlet(double *disteven, double *distodd, double dout,
								   int Nx, int Ny, int Nz, int outlet);

extern "C" void InitDenColor(char *ID, double *Den, double *Phi, double das, double dbs, int Nx, int Ny, int Nz);

extern "C" void InitDenColorDistance(char *ID, double *Den, double *Phi, double *Distance,
								double das, double dbs, double beta, double xp, int Nx, int Ny, int Nz);

extern "C" void ComputeColorGradient(char *ID, double *phi, double *ColorGrad, int Nx, int Ny, int Nz);

extern "C" void ColorCollide( char *ID, double *disteven, double *distodd, double *ColorGrad,
								double *Velocity, int Nx, int Ny, int Nz,double rlx_setA, double rlx_setB,
								double alpha, double beta, double Fx, double Fy, double Fz, bool pBC);

extern "C" void ColorCollideOpt( char *ID, double *disteven, double *distodd, double *phi, double *ColorGrad,
								double *Velocity, int Nx, int Ny, int Nz,double rlx_setA, double rlx_setB,
								double alpha, double beta, double Fx, double Fy, double Fz);

extern "C" void DensityStreamD3Q7(char *ID, double *Den, double *Copy, double *Phi, double *ColorGrad, double *Velocity,
		double beta, int Nx, int Ny, int Nz, bool pBC);

extern "C" void ComputePhi(char *ID, double *Phi, double *Den, int N);

extern "C" void MassColorCollideD3Q7(char *ID, double *A_even, double *A_odd, double *B_even, double *B_odd,
		double *Den, double *Phi, double *ColorGrad, double *Velocity, double beta, int N, bool pBC);

extern "C" void ColorBC_inlet(double *Phi, double *Den, double *A_even, double *A_odd,
								  double *B_even, double *B_odd, int Nx, int Ny, int Nz);

extern "C" void ColorBC_outlet(double *Phi, double *Den, double *A_even, double *A_odd,
								  double *B_even, double *B_odd, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_Velocity_BC_z(double *disteven, double *distodd, double uz,
								   int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_Velocity_BC_Z(double *disteven, double *distodd, double uz,
								   int Nx, int Ny, int Nz, int outlet);

extern "C" void SetPhiSlice_z(double *Phi, double value, int Nx, int Ny, int Nz, int Slice);

class ScaLBL_Communicator{
public:
	//......................................................................................
	ScaLBL_Communicator(Domain &Dm);
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

	void SendD3Q19(double *f_even, double *f_odd);
	void RecvD3Q19(double *f_even, double *f_odd);
	void BiSendD3Q7(double *A_even, double *A_odd, double *B_even, double *B_odd);
	void BiRecvD3Q7(double *A_even, double *A_odd, double *B_even, double *B_odd);
	void SendHalo(double *data);
	void RecvHalo(double *data);

private:
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
	AllocateDeviceMemory((void **) &sendbuf_x, 5*sendCount_x*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_X, 5*sendCount_X*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_y, 5*sendCount_y*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_Y, 5*sendCount_Y*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_z, 5*sendCount_z*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_Z, 5*sendCount_Z*sizeof(double));	// Allocatevoid * memory
	AllocateDeviceMemory((void **) &sendbuf_xy, sendCount_xy*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_xY, sendCount_xY*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_Xy, sendCount_Xy*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_XY, sendCount_XY*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_xz, sendCount_xz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_xZ, sendCount_xZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_Xz, sendCount_Xz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_XZ, sendCount_XZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_yz, sendCount_yz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_yZ, sendCount_yZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_Yz, sendCount_Yz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_YZ, sendCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	AllocateDeviceMemory((void **) &recvbuf_x, 5*recvCount_x*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_X, 5*recvCount_X*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_y, 5*recvCount_y*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_Y, 5*recvCount_Y*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_z, 5*recvCount_z*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_Z, 5*recvCount_Z*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_xy, recvCount_xy*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_xY, recvCount_xY*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_Xy, recvCount_Xy*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_XY, recvCount_XY*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_xz, recvCount_xz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_xZ, recvCount_xZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_Xz, recvCount_Xz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_XZ, recvCount_XZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_yz, recvCount_yz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_yZ, recvCount_yZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_Yz, recvCount_Yz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_YZ, recvCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	AllocateDeviceMemory((void **) &dvcSendList_x, sendCount_x*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_X, sendCount_X*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_y, sendCount_y*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_Y, sendCount_Y*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_z, sendCount_z*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_Z, sendCount_Z*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_xy, sendCount_xy*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_xY, sendCount_xY*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_Xy, sendCount_Xy*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_XY, sendCount_XY*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_xz, sendCount_xz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_xZ, sendCount_xZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_Xz, sendCount_Xz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_XZ, sendCount_XZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_yz, sendCount_yz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_yZ, sendCount_yZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_Yz, sendCount_Yz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_YZ, sendCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	AllocateDeviceMemory((void **) &dvcRecvList_x, recvCount_x*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_X, recvCount_X*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_y, recvCount_y*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_Y, recvCount_Y*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_z, recvCount_z*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_Z, recvCount_Z*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_xy, recvCount_xy*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_xY, recvCount_xY*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_Xy, recvCount_Xy*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_XY, recvCount_XY*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_xz, recvCount_xz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_xZ, recvCount_xZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_Xz, recvCount_Xz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_XZ, recvCount_XZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_yz, recvCount_yz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_yZ, recvCount_yZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_Yz, recvCount_Yz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_YZ, recvCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	CopyToDevice(dvcSendList_x,Dm.sendList_x,sendCount_x*sizeof(int));
	CopyToDevice(dvcSendList_X,Dm.sendList_X,sendCount_X*sizeof(int));
	CopyToDevice(dvcSendList_y,Dm.sendList_y,sendCount_y*sizeof(int));
	CopyToDevice(dvcSendList_Y,Dm.sendList_Y,sendCount_Y*sizeof(int));
	CopyToDevice(dvcSendList_z,Dm.sendList_z,sendCount_z*sizeof(int));
	CopyToDevice(dvcSendList_Z,Dm.sendList_Z,sendCount_Z*sizeof(int));
	CopyToDevice(dvcSendList_xy,Dm.sendList_xy,sendCount_xy*sizeof(int));
	CopyToDevice(dvcSendList_XY,Dm.sendList_XY,sendCount_XY*sizeof(int));
	CopyToDevice(dvcSendList_xY,Dm.sendList_xY,sendCount_xY*sizeof(int));
	CopyToDevice(dvcSendList_Xy,Dm.sendList_Xy,sendCount_Xy*sizeof(int));
	CopyToDevice(dvcSendList_xz,Dm.sendList_xz,sendCount_xz*sizeof(int));
	CopyToDevice(dvcSendList_XZ,Dm.sendList_XZ,sendCount_XZ*sizeof(int));
	CopyToDevice(dvcSendList_xZ,Dm.sendList_xZ,sendCount_xZ*sizeof(int));
	CopyToDevice(dvcSendList_Xz,Dm.sendList_Xz,sendCount_Xz*sizeof(int));
	CopyToDevice(dvcSendList_yz,Dm.sendList_yz,sendCount_yz*sizeof(int));
	CopyToDevice(dvcSendList_YZ,Dm.sendList_YZ,sendCount_YZ*sizeof(int));
	CopyToDevice(dvcSendList_yZ,Dm.sendList_yZ,sendCount_yZ*sizeof(int));
	CopyToDevice(dvcSendList_Yz,Dm.sendList_Yz,sendCount_Yz*sizeof(int));
	//......................................................................................
	CopyToDevice(dvcRecvList_x,Dm.recvList_x,recvCount_x*sizeof(int));
	CopyToDevice(dvcRecvList_X,Dm.recvList_X,recvCount_X*sizeof(int));
	CopyToDevice(dvcRecvList_y,Dm.recvList_y,recvCount_y*sizeof(int));
	CopyToDevice(dvcRecvList_Y,Dm.recvList_Y,recvCount_Y*sizeof(int));
	CopyToDevice(dvcRecvList_z,Dm.recvList_z,recvCount_z*sizeof(int));
	CopyToDevice(dvcRecvList_Z,Dm.recvList_Z,recvCount_Z*sizeof(int));
	CopyToDevice(dvcRecvList_xy,Dm.recvList_xy,recvCount_xy*sizeof(int));
	CopyToDevice(dvcRecvList_XY,Dm.recvList_XY,recvCount_XY*sizeof(int));
	CopyToDevice(dvcRecvList_xY,Dm.recvList_xY,recvCount_xY*sizeof(int));
	CopyToDevice(dvcRecvList_Xy,Dm.recvList_Xy,recvCount_Xy*sizeof(int));
	CopyToDevice(dvcRecvList_xz,Dm.recvList_xz,recvCount_xz*sizeof(int));
	CopyToDevice(dvcRecvList_XZ,Dm.recvList_XZ,recvCount_XZ*sizeof(int));
	CopyToDevice(dvcRecvList_xZ,Dm.recvList_xZ,recvCount_xZ*sizeof(int));
	CopyToDevice(dvcRecvList_Xz,Dm.recvList_Xz,recvCount_Xz*sizeof(int));
	CopyToDevice(dvcRecvList_yz,Dm.recvList_yz,recvCount_yz*sizeof(int));
	CopyToDevice(dvcRecvList_YZ,Dm.recvList_YZ,recvCount_YZ*sizeof(int));
	CopyToDevice(dvcRecvList_yZ,Dm.recvList_yZ,recvCount_yZ*sizeof(int));
	CopyToDevice(dvcRecvList_Yz,Dm.recvList_Yz,recvCount_Yz*sizeof(int));
	//......................................................................................
	MPI_Barrier(MPI_COMM_SCALBL);
	DeviceBarrier();
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

void ScaLBL_Communicator::SendD3Q19(double *f_even, double *f_odd){

	if (Lock==true){
		ERROR("ScaLBL Error (SendD3Q19): ScaLBL_Communicator is locked -- did you forget to match Send/Recv calls?");
	}
	else{
		Lock=true;
	}
	// assign tag of 19 to D3Q19 communication
	sendtag = recvtag = 19;
	DeviceBarrier();
	// Pack the distributions
	PackDist(1,dvcSendList_x,0,sendCount_x,sendbuf_x,f_even,N);
	PackDist(4,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	PackDist(5,dvcSendList_x,2*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	PackDist(6,dvcSendList_x,3*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	PackDist(7,dvcSendList_x,4*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	//...Packing for X face(1,7,9,11,13)................................
	PackDist(0,dvcSendList_X,0,sendCount_X,sendbuf_X,f_odd,N);
	PackDist(3,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	PackDist(4,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	PackDist(5,dvcSendList_X,3*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	PackDist(6,dvcSendList_X,4*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	//...Packing for y face(4,8,9,16,18).................................
	PackDist(2,dvcSendList_y,0,sendCount_y,sendbuf_y,f_even,N);
	PackDist(4,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	PackDist(4,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,f_odd,N);
	PackDist(8,dvcSendList_y,3*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	PackDist(9,dvcSendList_y,4*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	//...Packing for Y face(3,7,10,15,17).................................
	PackDist(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,f_odd,N);
	PackDist(3,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	PackDist(5,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,f_even,N);
	PackDist(7,dvcSendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	PackDist(8,dvcSendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	//...Packing for z face(6,12,13,16,17)................................
	PackDist(3,dvcSendList_z,0,sendCount_z,sendbuf_z,f_even,N);
	PackDist(6,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,f_even,N);
	PackDist(6,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
	PackDist(8,dvcSendList_z,3*sendCount_z,sendCount_z,sendbuf_z,f_even,N);
	PackDist(8,dvcSendList_z,4*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
	//...Packing for Z face(5,11,14,15,18)................................
	PackDist(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,f_odd,N);
	PackDist(5,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
	PackDist(7,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
	PackDist(7,dvcSendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
	PackDist(9,dvcSendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
	//...Pack the xy edge (8)................................
	PackDist(4,dvcSendList_xy,0,sendCount_xy,sendbuf_xy,f_even,N);
	//...Pack the Xy edge (9)................................
	PackDist(4,dvcSendList_Xy,0,sendCount_Xy,sendbuf_Xy,f_odd,N);
	//...Pack the xY edge (10)................................
	PackDist(5,dvcSendList_xY,0,sendCount_xY,sendbuf_xY,f_even,N);
	//...Pack the XY edge (7)................................
	PackDist(3,dvcSendList_XY,0,sendCount_XY,sendbuf_XY,f_odd,N);
	//...Pack the xz edge (12)................................
	PackDist(6,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
	//...Pack the xZ edge (14)................................
	PackDist(7,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
	//...Pack the Xz edge (13)................................
	PackDist(6,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
	//...Pack the XZ edge (11)................................
	PackDist(5,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
	//...Pack the xz edge (12)................................
	//...Pack the yz edge (16)................................
	PackDist(8,dvcSendList_yz,0,sendCount_yz,sendbuf_yz,f_even,N);
	//...Pack the yZ edge (18)................................
	PackDist(9,dvcSendList_yZ,0,sendCount_yZ,sendbuf_yZ,f_even,N);
	//...Pack the Yz edge (17)................................
	PackDist(8,dvcSendList_Yz,0,sendCount_Yz,sendbuf_Yz,f_odd,N);
	//...Pack the YZ edge (15)................................
	PackDist(7,dvcSendList_YZ,0,sendCount_YZ,sendbuf_YZ,f_odd,N);
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
	DeviceBarrier();

	//...................................................................................
	// Unpack the distributions on the device
	//...................................................................................
	//...Map recieve list for the X face: q=2,8,10,12,13 .................................
	UnpackDist(0,-1,0,0,dvcRecvList_X,0,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
	UnpackDist(3,-1,-1,0,dvcRecvList_X,recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
	UnpackDist(4,-1,1,0,dvcRecvList_X,2*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
	UnpackDist(5,-1,0,-1,dvcRecvList_X,3*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
	UnpackDist(6,-1,0,1,dvcRecvList_X,4*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the x face: q=1,7,9,11,13..................................
	UnpackDist(1,1,0,0,dvcRecvList_x,0,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
	UnpackDist(4,1,1,0,dvcRecvList_x,recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
	UnpackDist(5,1,-1,0,dvcRecvList_x,2*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
	UnpackDist(6,1,0,1,dvcRecvList_x,3*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
	UnpackDist(7,1,0,-1,dvcRecvList_x,4*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
	UnpackDist(1,0,-1,0,dvcRecvList_Y,0,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
	UnpackDist(3,-1,-1,0,dvcRecvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
	UnpackDist(5,1,-1,0,dvcRecvList_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,f_even,Nx,Ny,Nz);
	UnpackDist(7,0,-1,-1,dvcRecvList_Y,3*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
	UnpackDist(8,0,-1,1,dvcRecvList_Y,4*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
	UnpackDist(2,0,1,0,dvcRecvList_y,0,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
	UnpackDist(4,1,1,0,dvcRecvList_y,recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
	UnpackDist(4,-1,1,0,dvcRecvList_y,2*recvCount_y,recvCount_y,recvbuf_y,f_odd,Nx,Ny,Nz);
	UnpackDist(8,0,1,1,dvcRecvList_y,3*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
	UnpackDist(9,0,1,-1,dvcRecvList_y,4*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
	UnpackDist(2,0,0,-1,dvcRecvList_Z,0,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
	UnpackDist(5,-1,0,-1,dvcRecvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
	UnpackDist(7,1,0,-1,dvcRecvList_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
	UnpackDist(7,0,-1,-1,dvcRecvList_Z,3*recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
	UnpackDist(9,0,1,-1,dvcRecvList_Z,4*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
	//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
	UnpackDist(3,0,0,1,dvcRecvList_z,0,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
	UnpackDist(6,1,0,1,dvcRecvList_z,recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
	UnpackDist(6,-1,0,1,dvcRecvList_z,2*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
	UnpackDist(8,0,1,1,dvcRecvList_z,3*recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
	UnpackDist(8,0,-1,1,dvcRecvList_z,4*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
	//..................................................................................
	//...Map recieve list for the xy edge <<<8)................................
	UnpackDist(3,-1,-1,0,dvcRecvList_XY,0,recvCount_XY,recvbuf_XY,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the Xy edge <<<9)................................
	UnpackDist(5,1,-1,0,dvcRecvList_xY,0,recvCount_xY,recvbuf_xY,f_even,Nx,Ny,Nz);
	//...Map recieve list for the xY edge <<<10)................................
	UnpackDist(4,-1,1,0,dvcRecvList_Xy,0,recvCount_Xy,recvbuf_Xy,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the XY edge <<<7)................................
	UnpackDist(4,1,1,0,dvcRecvList_xy,0,recvCount_xy,recvbuf_xy,f_even,Nx,Ny,Nz);
	//...Map recieve list for the xz edge <<<12)................................
	UnpackDist(5,-1,0,-1,dvcRecvList_XZ,0,recvCount_XZ,recvbuf_XZ,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the xZ edge <<<14)................................
	UnpackDist(6,-1,0,1,dvcRecvList_Xz,0,recvCount_Xz,recvbuf_Xz,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the Xz edge <<<13)................................
	UnpackDist(7,1,0,-1,dvcRecvList_xZ,0,recvCount_xZ,recvbuf_xZ,f_even,Nx,Ny,Nz);
	//...Map recieve list for the XZ edge <<<11)................................
	UnpackDist(6,1,0,1,dvcRecvList_xz,0,recvCount_xz,recvbuf_xz,f_even,Nx,Ny,Nz);
	//...Map recieve list for the yz edge <<<16)................................
	UnpackDist(7,0,-1,-1,dvcRecvList_YZ,0,recvCount_YZ,recvbuf_YZ,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the yZ edge <<<18)................................
	UnpackDist(8,0,-1,1,dvcRecvList_Yz,0,recvCount_Yz,recvbuf_Yz,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the Yz edge <<<17)................................
	UnpackDist(9,0,1,-1,dvcRecvList_yZ,0,recvCount_yZ,recvbuf_yZ,f_even,Nx,Ny,Nz);
	//...Map recieve list for the YZ edge <<<15)................................
	UnpackDist(8,0,1,1,dvcRecvList_yz,0,recvCount_yz,recvbuf_yz,f_even,Nx,Ny,Nz);
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
	DeviceBarrier();
	//...................................................................................
	sendtag = recvtag = 7;
	//...................................................................................
	PackDist(1,dvcSendList_x,0,sendCount_x,sendbuf_x,A_even,N);
	PackDist(1,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,B_even,N);
	//...Packing for X face(1,7,9,11,13)................................
	PackDist(0,dvcSendList_X,0,sendCount_X,sendbuf_X,A_odd,N);
	PackDist(0,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,B_odd,N);
	//...Packing for y face(4,8,9,16,18).................................
	PackDist(2,dvcSendList_y,0,sendCount_y,sendbuf_y,A_even,N);
	PackDist(2,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,B_even,N);
	//...Packing for Y face(3,7,10,15,17).................................
	PackDist(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,A_odd,N);
	PackDist(1,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,B_odd,N);
	//...Packing for z face(6,12,13,16,17)................................
	PackDist(3,dvcSendList_z,0,sendCount_z,sendbuf_z,A_even,N);
	PackDist(3,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,B_even,N);
	//...Packing for Z face(5,11,14,15,18)................................
	PackDist(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,A_odd,N);
	PackDist(2,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,B_odd,N);
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
	DeviceBarrier();
	//...................................................................................
	// Unpack the distributions on the device
	//...................................................................................
	//...Map recieve list for the X face: q=2,8,10,12,13 ................................
	UnpackDist(0,-1,0,0,dvcRecvList_X,0,recvCount_X,recvbuf_X,A_odd,Nx,Ny,Nz);
	UnpackDist(0,-1,0,0,dvcRecvList_X,recvCount_X,recvCount_X,recvbuf_X,B_odd,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the x face: q=1,7,9,11,13..................................
	UnpackDist(1,1,0,0,dvcRecvList_x,0,recvCount_x,recvbuf_x,A_even,Nx,Ny,Nz);
	UnpackDist(1,1,0,0,dvcRecvList_x,recvCount_x,recvCount_x,recvbuf_x,B_even,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the y face: q=4,8,9,16,18 .................................
	UnpackDist(1,0,-1,0,dvcRecvList_Y,0,recvCount_Y,recvbuf_Y,A_odd,Nx,Ny,Nz);
	UnpackDist(1,0,-1,0,dvcRecvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,B_odd,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the Y face: q=3,7,10,15,17 ................................
	UnpackDist(2,0,1,0,dvcRecvList_y,0,recvCount_y,recvbuf_y,A_even,Nx,Ny,Nz);
	UnpackDist(2,0,1,0,dvcRecvList_y,recvCount_y,recvCount_y,recvbuf_y,B_even,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the z face<<<6,12,13,16,17)................................
	UnpackDist(2,0,0,-1,dvcRecvList_Z,0,recvCount_Z,recvbuf_Z,A_odd,Nx,Ny,Nz);
	UnpackDist(2,0,0,-1,dvcRecvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,B_odd,Nx,Ny,Nz);
	//...Map recieve list for the Z face<<<5,11,14,15,18)................................
	UnpackDist(3,0,0,1,dvcRecvList_z,0,recvCount_z,recvbuf_z,A_even,Nx,Ny,Nz);
	UnpackDist(3,0,0,1,dvcRecvList_z,recvCount_z,recvCount_z,recvbuf_z,B_even,Nx,Ny,Nz);
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
	DeviceBarrier();
	//...................................................................................
	sendtag = recvtag = 1;
	//...................................................................................
	PackValues(dvcSendList_x, sendCount_x,sendbuf_x, data, N);
	PackValues(dvcSendList_y, sendCount_y,sendbuf_y, data, N);
	PackValues(dvcSendList_z, sendCount_z,sendbuf_z, data, N);
	PackValues(dvcSendList_X, sendCount_X,sendbuf_X, data, N);
	PackValues(dvcSendList_Y, sendCount_Y,sendbuf_Y, data, N);
	PackValues(dvcSendList_Z, sendCount_Z,sendbuf_Z, data, N);
	PackValues(dvcSendList_xy, sendCount_xy,sendbuf_xy, data, N);
	PackValues(dvcSendList_xY, sendCount_xY,sendbuf_xY, data, N);
	PackValues(dvcSendList_Xy, sendCount_Xy,sendbuf_Xy, data, N);
	PackValues(dvcSendList_XY, sendCount_XY,sendbuf_XY, data, N);
	PackValues(dvcSendList_xz, sendCount_xz,sendbuf_xz, data, N);
	PackValues(dvcSendList_xZ, sendCount_xZ,sendbuf_xZ, data, N);
	PackValues(dvcSendList_Xz, sendCount_Xz,sendbuf_Xz, data, N);
	PackValues(dvcSendList_XZ, sendCount_XZ,sendbuf_XZ, data, N);
	PackValues(dvcSendList_yz, sendCount_yz,sendbuf_yz, data, N);
	PackValues(dvcSendList_yZ, sendCount_yZ,sendbuf_yZ, data, N);
	PackValues(dvcSendList_Yz, sendCount_Yz,sendbuf_Yz, data, N);
	PackValues(dvcSendList_YZ, sendCount_YZ,sendbuf_YZ, data, N);
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
	DeviceBarrier();
	//...................................................................................
	//...................................................................................
	UnpackValues(dvcRecvList_x, recvCount_x,recvbuf_x, data, N);
	UnpackValues(dvcRecvList_y, recvCount_y,recvbuf_y, data, N);
	UnpackValues(dvcRecvList_z, recvCount_z,recvbuf_z, data, N);
	UnpackValues(dvcRecvList_X, recvCount_X,recvbuf_X, data, N);
	UnpackValues(dvcRecvList_Y, recvCount_Y,recvbuf_Y, data, N);
	UnpackValues(dvcRecvList_Z, recvCount_Z,recvbuf_Z, data, N);
	UnpackValues(dvcRecvList_xy, recvCount_xy,recvbuf_xy, data, N);
	UnpackValues(dvcRecvList_xY, recvCount_xY,recvbuf_xY, data, N);
	UnpackValues(dvcRecvList_Xy, recvCount_Xy,recvbuf_Xy, data, N);
	UnpackValues(dvcRecvList_XY, recvCount_XY,recvbuf_XY, data, N);
	UnpackValues(dvcRecvList_xz, recvCount_xz,recvbuf_xz, data, N);
	UnpackValues(dvcRecvList_xZ, recvCount_xZ,recvbuf_xZ, data, N);
	UnpackValues(dvcRecvList_Xz, recvCount_Xz,recvbuf_Xz, data, N);
	UnpackValues(dvcRecvList_XZ, recvCount_XZ,recvbuf_XZ, data, N);
	UnpackValues(dvcRecvList_yz, recvCount_yz,recvbuf_yz, data, N);
	UnpackValues(dvcRecvList_yZ, recvCount_yZ,recvbuf_yZ, data, N);
	UnpackValues(dvcRecvList_Yz, recvCount_Yz,recvbuf_Yz, data, N);
	UnpackValues(dvcRecvList_YZ, recvCount_YZ,recvbuf_YZ, data, N);
	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................
}

