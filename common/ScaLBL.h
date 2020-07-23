/* ScaLBL.h 
 *  Header file for Scalable Lattice Boltzmann Library
 *  Separate implementations for GPU and CPU must both follow the conventions defined in this header
 *  This libarry contains the essential components of the LBM
 *     - streaming implementations
 *     - collision terms to model various physics
 *     - communication framework for the LBM
 *  Refer to Domain.h for setup of parallel domains
 */
#ifndef ScalLBL_H
#define ScalLBL_H
#include "common/Domain.h"

extern "C" int ScaLBL_SetDevice(int rank);

extern "C" void ScaLBL_AllocateDeviceMemory(void** address, size_t size);

extern "C" void ScaLBL_FreeDeviceMemory(void* pointer);

extern "C" void ScaLBL_CopyToDevice(void* dest, const void* source, size_t size);

extern "C" void ScaLBL_CopyToHost(void* dest, const void* source, size_t size);

extern "C" void ScaLBL_AllocateZeroCopy(void** address, size_t size);

extern "C" void ScaLBL_CopyToZeroCopy(void* dest, const void* source, size_t size);

extern "C" void ScaLBL_DeviceBarrier();

extern "C" void ScaLBL_D3Q19_Pack(int q, int *list, int start, int count, double *sendbuf, double *dist, int N);

extern "C" void ScaLBL_D3Q19_Unpack(int q, int *list, int start, int count, double *recvbuf, double *dist, int N);

extern "C" void ScaLBL_D3Q7_Unpack(int q, int *list,  int start, int count, double *recvbuf, double *dist, int N);

extern "C" void ScaLBL_Scalar_Pack(int *list, int count, double *sendbuf, double *Data, int N);

extern "C" void ScaLBL_Scalar_Unpack(int *list, int count, double *recvbuf, double *Data, int N);

extern "C" void ScaLBL_Gradient_Unpack(double weight, double Cqx, double Cqy, double Cqz, 
		int *list, int start, int count, double *recvbuf, double *phi, double *grad, int N);

extern "C" void ScaLBL_PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N);

extern "C" void ScaLBL_UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N);

extern "C" void ScaLBL_D3Q19_Init(double *Dist, int Np);


extern "C" void ScaLBL_D3Q19_Momentum(double *dist, double *vel, int Np);

extern "C" void ScaLBL_D3Q19_Pressure(double *dist, double *press, int Np);

// BGK MODEL
extern "C" void ScaLBL_D3Q19_AAeven_BGK(double *dist, int start, int finish, int Np, double rlx, double Fx, double Fy, double Fz);

extern "C" void ScaLBL_D3Q19_AAodd_BGK(int *neighborList, double *dist, int start, int finish, int Np, double rlx, double Fx, double Fy, double Fz);

// GREYSCALE MODEL (Single-component)

extern "C" void ScaLBL_D3Q19_GreyIMRT_Init(double *Dist, int Np, double Den);

extern "C" void ScaLBL_D3Q19_AAeven_Greyscale(double *dist, int start, int finish, int Np, double rlx, double rlx_eff, double Fx, double Fy, double Fz,
                                              double *Poros,double *Perm, double *Velocity,double *Pressure);

extern "C" void ScaLBL_D3Q19_AAodd_Greyscale(int *neighborList, double *dist, int start, int finish, int Np, double rlx, double rlx_eff, double Fx, double Fy, double Fz, 
                                             double *Poros,double *Perm, double *Velocity,double *Pressure);

extern "C" void ScaLBL_D3Q19_AAeven_Greyscale_IMRT(double *dist, int start, int finish, int Np, double rlx, double rlx_eff, double Fx, double Fy, double Fz,
                                              double *Poros,double *Perm, double *Velocity,double Den,double *Pressure);

extern "C" void ScaLBL_D3Q19_AAodd_Greyscale_IMRT(int *neighborList, double *dist, int start, int finish, int Np, double rlx, double rlx_eff, double Fx, double Fy, double Fz, 
                                             double *Poros,double *Perm, double *Velocity,double Den,double *Pressure);

extern "C" void ScaLBL_D3Q19_AAeven_Greyscale_MRT(double *dist, int start, int finish, int Np, double rlx, double rlx_eff, double Fx, double Fy, double Fz,
                                              double *Poros,double *Perm, double *Velocity,double Den,double *Pressure);

extern "C" void ScaLBL_D3Q19_AAodd_Greyscale_MRT(int *neighborList, double *dist, int start, int finish, int Np, double rlx, double rlx_eff, double Fx, double Fy, double Fz, 
                                             double *Poros,double *Perm, double *Velocity,double Den,double *Pressure);
// GREYSCALE FREE-ENERGY MODEL (Two-component)

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleFE(double *dist, double *Aq, double *Bq, double *Den,
                double *DenGradA, double *DenGradB, double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double Gsc, double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure);

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleFE(int *neighborList, double *dist, double *Aq, double *Bq, double *Den,
                double *DenGradA, double *DenGradB, double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double Gsc, double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure);

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleFEChem(double *dist, double *Cq, double *Phi, double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double gamma,double kappaA,double kappaB,double lambdaA,double lambdaB,
                double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure,double *PressureGrad,double *PressTensorGrad,double *PhiLap);

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleFEChem(int *neighborList, double *dist, double *Cq, double *Phi, double *SolidForce, int start, int finish, int Np,
                double tauA,double tauB,double tauA_eff,double tauB_eff,double rhoA,double rhoB,double gamma,double kappaA,double kappaB,double lambdaA,double lambdaB,
                double Gx, double Gy, double Gz,
                double *Poros,double *Perm, double *Velocity,double *Pressure,double *PressureGrad,double *PressTensorGrad,double *PhiLap);

extern "C" void ScaLBL_D3Q7_GreyscaleFE_Init(double *Den, double *Cq, double *PhiLap, double gamma, double kappaA, double kappaB, double lambdaA, double lambdaB, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_GreyscaleFE_IMRT_Init(double *dist, double *Den, double rhoA, double rhoB, int Np);

extern "C" void ScaLBL_D3Q7_AAodd_GreyscaleFEDensity(int *NeighborList, double *Aq, double *Bq, double *Den, double *Phi, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_GreyscaleFEDensity(double *Aq, double *Bq, double *Den, double *Phi, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAodd_GreyscaleFEPhi(int *NeighborList, double *Cq, double *Phi, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_GreyscaleFEPhi(double *Cq, double *Phi, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_GreyscaleFE_Gradient(int *neighborList, double *Den, double *DenGrad, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_GreyscaleFE_Laplacian(int *neighborList, double *Den, double *DenLap, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_GreyscaleFE_Pressure(double *dist, double *Den, double *Porosity,double *Velocity,
                double *Pressure, double rhoA,double rhoB, int Np);

extern "C" void ScaLBL_D3Q19_GreyscaleFE_PressureTensor(int *neighborList, double *Phi,double *Pressure, double *PressTensor, double *PhiLap,
      		     double kappaA,double kappaB,double lambdaA,double lambdaB, int start, int finish, int Np);

// GREYSCALE SHAN-CHEN MODEL (Two-component)

extern "C" void ScaLBL_D3Q19_GreyscaleSC_Init(int *Map, double *distA, double *distB, double *DenA, double *DenB, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleSC_Density(int *NeighborList, int *Map, double *distA, double *distB, double *DenA, double *DenB, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleSC_Density(int *Map, double *distA, double *distB, double *DenA, double *DenB, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleSC_MRT(int *neighborList, int *Mpa, double *distA, double *distB, double *DenA,double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleSC_MRT(int *Map,double *distA, double *distB, double *DenA,double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleSC_BGK(int *neighborList, int *Map, double *distA, double *distB, double *DenA, double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleSC_BGK(int *Map, double *distA, double *distB, double *DenA, double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_GreyscaleSC_Gradient(int *neighborList, int *Map, double *Den, double *DenGrad, int strideY, int strideZ,int start, int finish, int Np);

extern "C" void ScaLBL_GreyscaleSC_BC_z(int *list, int *Map, double *DenA, double *DenB, double vA, double vB, int count);

extern "C" void ScaLBL_GreyscaleSC_BC_Z(int *list, int *Map, double *DenA, double *DenB, double vA, double vB, int count);

extern "C" void ScaLBL_GreyscaleSC_AAeven_Pressure_BC_z(int *list, double *distA, double *distB, double dinA, double dinB, int count, int N);

extern "C" void ScaLBL_GreyscaleSC_AAeven_Pressure_BC_Z(int *list, double *distA, double *distB, double doutA, double doutB, int count, int N);

extern "C" void ScaLBL_GreyscaleSC_AAodd_Pressure_BC_z(int *neighborList, int *list, double *distA, double *distB, double dinA, double dinB, int count, int N);

extern "C" void ScaLBL_GreyscaleSC_AAodd_Pressure_BC_Z(int *neighborList, int *list, double *distA, double *distB, double doutA, double doutB, int count, int N);

// GREYSCALE COLOR MODEL (Two-component)
//extern "C" void ScaLBL_D3Q19_GreyscaleColor_Init(double *dist, double *Porosity, int Np);

//extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColor(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
//        double *ColorGrad,double *Phi,double *GreySolidGrad, double *Poros,double *Perm,double *Vel, 
//        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
//		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np);
//
//extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColor(int *d_neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den, 
//		double *ColorGrad,double *Phi, double *GreySolidGrad, double *Poros,double *Perm,double *Vel, 
//        double rhoA, double rhoB, double tauA, double tauB, double tauA_eff,double tauB_eff, double alpha, double beta,
//		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColor(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
        double *Phi,double *GreySolidGrad, double *Poros,double *Perm,double *Vel, 
        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColor(int *d_neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *GreySolidGrad, double *Poros,double *Perm,double *Vel, 
        double rhoA, double rhoB, double tauA, double tauB, double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np);

// MRT MODEL
extern "C" void ScaLBL_D3Q19_AAeven_MRT(double *dist, int start, int finish, int Np, double rlx_setA, double rlx_setB, double Fx,
		double Fy, double Fz);

extern "C" void ScaLBL_D3Q19_AAodd_MRT(int *d_neighborList, double *dist, int start, int finish, int Np,
		double rlx_setA, double rlx_setB, double Fx, double Fy, double Fz);

// COLOR MODEL

extern "C" void ScaLBL_D3Q19_AAeven_Color(int *Map, double *dist, double *Aq, double *Bq, double *Den, double *Phi,
		double *Vel, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_Color(int *d_neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *Vel, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAodd_PhaseField(int *NeighborList, int *Map, double *Aq, double *Bq, 
			double *Den, double *Phi, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_PhaseField(int *Map, double *Aq, double *Bq, double *Den, double *Phi, 
			int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_Gradient(int *Map, double *Phi, double *ColorGrad, int start, int finish, int Np, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_PhaseField_Init(int *Map, double *Phi, double *Den, double *Aq, double *Bq, int start, int finish, int Np);

// Density functional hydrodynamics LBM
extern "C" void ScaLBL_DFH_Init(double *Phi, double *Den, double *Aq, double *Bq, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_DFH(int *neighborList, double *dist, double *Aq, double *Bq, double *Den, double *Phi,
		double *Gradient, double *SolidForce, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_DFH(int *neighborList, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *Gradient, double *SolidForce, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAodd_DFH(int *NeighborList, double *Aq, double *Bq, double *Den, double *Phi, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_DFH(double *Aq, double *Bq, double *Den, double *Phi, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_Gradient_DFH(int *NeighborList, double *Phi, double *ColorGrad, int start, int finish, int Np);

// BOUNDARY CONDITION ROUTINES

extern "C" void ScaLBL_D3Q19_AAodd_Pressure_BC_z(int *neighborList, int *list, double *dist, double din, int count, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_Pressure_BC_Z(int *neighborList, int *list, double *dist, double dout, int count, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_Pressure_BC_z(int *list, double *dist, double din, int count, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_Pressure_BC_Z(int *list, double *dist, double dout, int count, int Np);

extern "C" double ScaLBL_D3Q19_AAodd_Flux_BC_z(int *neighborList, int *list, double *dist, double flux, 
		double area, int count, int N);

extern "C" double ScaLBL_D3Q19_AAeven_Flux_BC_z(int *list, double *dist, double flux, double area, 
		 int count, int N);

extern "C" void ScaLBL_Color_BC_z(int *list, int *Map, double *Phi, double *Den, double vA, double vB, int count, int Np);

extern "C" void ScaLBL_Color_BC_Z(int *list, int *Map, double *Phi, double *Den, double vA, double vB, int count, int Np);

extern "C" void ScaLBL_SetSlice_z(double *Phi, double value, int Nx, int Ny, int Nz, int Slice);



class ScaLBL_Communicator{
public:
	//......................................................................................
	ScaLBL_Communicator(std::shared_ptr <Domain> Dm);

	//ScaLBL_Communicator(Domain &Dm, IntArray &Map);
	~ScaLBL_Communicator();
	//......................................................................................
	unsigned long int CommunicationCount,SendCount,RecvCount;
	int Nx,Ny,Nz,N;
	int BoundaryCondition;
	
	int next;
	int first_interior,last_interior;
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

	int LastExterior();
	int FirstInterior();
	int LastInterior();
	
	int MemoryOptimizedLayoutAA(IntArray &Map, int *neighborList, signed char *id, int Np);
//	void MemoryOptimizedLayout(IntArray &Map, int *neighborList, char *id, int Np);
//	void MemoryOptimizedLayoutFull(IntArray &Map, int *neighborList, char *id, int Np);
//	void MemoryDenseLayout(IntArray &Map, int *neighborList, char *id, int Np);
//	void MemoryDenseLayoutFull(IntArray &Map, int *neighborList, char *id, int Np);
//	void SendD3Q19(double *f_even, double *f_odd);
//	void RecvD3Q19(double *f_even, double *f_odd);
//	void SendD3Q19AA(double *f_even, double *f_odd);
//	void RecvD3Q19AA(double *f_even, double *f_odd);
	void SendD3Q19AA(double *dist);
	void RecvD3Q19AA(double *dist);
	void BiSendD3Q19AA(double *Aq, double *Bq);
	void BiRecvD3Q19AA(double *Aq, double *Bq);
//	void BiSendD3Q7(double *A_even, double *A_odd, double *B_even, double *B_odd);
//	void BiRecvD3Q7(double *A_even, double *A_odd, double *B_even, double *B_odd);
	void SendD3Q7AA(double *Aq);
	void RecvD3Q7AA(double *Aq);
	void BiSendD3Q7AA(double *Aq, double *Bq);
	void BiRecvD3Q7AA(double *Aq, double *Bq);
	void TriSendD3Q7AA(double *Aq, double *Bq, double *Cq);
	void TriRecvD3Q7AA(double *Aq, double *Bq, double *Cq);
	void SendHalo(double *data);
	void RecvHalo(double *data);
	void RecvGrad(double *Phi, double *Gradient);
	void RegularLayout(IntArray map, const double *data, DoubleArray &regdata);

	// Routines to set boundary conditions
	void Color_BC_z(int *Map, double *Phi, double *Den, double vA, double vB);
	void Color_BC_Z(int *Map, double *Phi, double *Den, double vA, double vB);
	void D3Q19_Pressure_BC_z(int *neighborList, double *fq, double din, int time);
	void D3Q19_Pressure_BC_Z(int *neighborList, double *fq, double dout, int time);
	double D3Q19_Flux_BC_z(int *neighborList, double *fq, double flux, int time);
	void GreyscaleSC_BC_z(int *Map, double *DenA, double *DenB, double vA, double vB);
	void GreyscaleSC_BC_Z(int *Map, double *DenA, double *DenB, double vA, double vB);
    void GreyscaleSC_Pressure_BC_z(int *neighborList, double *fqA, double *fqB, double dinA, double dinB, int time);
    void GreyscaleSC_Pressure_BC_Z(int *neighborList, double *fqA, double *fqB, double doutA, double doutB, int time);
//	void TestSendD3Q19(double *f_even, double *f_odd);
//	void TestRecvD3Q19(double *f_even, double *f_odd);

	// Debugging and unit testing functions
	void PrintD3Q19();

private:
	//void D3Q19_MapRecv_OLD(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count, int *d3q19_recvlist);
	void D3Q19_MapRecv(int Cqx, int Cqy, int Cqz, int *list,  int start, int count, int *d3q19_recvlist);

	bool Lock; 	// use Lock to make sure only one call at a time to protect data in transit
	// only one set of Send requests can be active at any time (per instance)
	int i,j,k,n;

	int iproc,jproc,kproc;
	int nprocx,nprocy,nprocz;
	int sendtag,recvtag;
	// Give the object it's own MPI communicator
	RankInfoStruct rank_info;
	Utilities::MPI MPI_COMM_SCALBL;		// MPI Communicator for this domain
	MPI_Request req1[18],req2[18];
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


#endif
