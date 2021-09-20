/** @file  ScaLBL.h */
/*  \details Header file for Scalable Lattice Boltzmann Library
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

/**
* \brief Set compute device
* @param rank      rank of MPI process 
*/
extern "C" int ScaLBL_SetDevice(int rank);

/**
* \brief Allocate memory 
* @param address      memory address
* @param size         size in bytes
*/
extern "C" void ScaLBL_AllocateDeviceMemory(void** address, size_t size);

/**
* \brief Free memory
* @param pointer         pointer to memory to free
*/
extern "C" void ScaLBL_FreeDeviceMemory(void* pointer);

/**
* \brief Copy memory from host to device 
* \details  Device memory should be close to simulation (based on NUMA cost)
*          Host memory may be a shared memory region (with possibly higher NUMA cost for simulation)
*          Analysis routine should minimize NUMA for host memory (based on process placement)
* @param dest    memory location to copy to
* @param source  memory region to copy from
* @param size    size of the region to copy in bytes        
*/
extern "C" void ScaLBL_CopyToDevice(void* dest, const void* source, size_t size);


/**
* \brief Copy memory from device to host
* \details  Device memory should be close to simulation (based on NUMA cost)
*          Host memory may be a shared memory region (with possibly higher NUMA cost for simulation)
*          Analysis routine should minimize NUMA for host memory (based on process placement)
* @param dest    memory location to copy to
* @param source  memory region to copy from
* @param size    size of the region to copy in bytes        
*/
extern "C" void ScaLBL_CopyToHost(void* dest, const void* source, size_t size);

/**
=* \brief Allocate zero copy memory buffer (i.e. shared memory)
* @param address      memory address
* @param size         size in bytes
*/
extern "C" void ScaLBL_AllocateZeroCopy(void** address, size_t size);

/**
* \brief Copy memory from host to zero copy buffer
* \details  Device memory should be close to simulation (based on NUMA cost)
*          Host memory may be a shared memory region (with possibly higher NUMA cost for simulation)
*          Analysis routine should minimize NUMA for host memory (based on process placement)
* @param dest    memory location to copy to
* @param source  memory region to copy from
* @param size    size of the region to copy in bytes        
*/
extern "C" void ScaLBL_CopyToZeroCopy(void* dest, const void* source, size_t size);

/**
* \brief  Device barrier routine     
*/
extern "C" void ScaLBL_DeviceBarrier();

/**
* \brief  Pack D3Q19 distributions for communication
* @param q  - index for distribution based on D3Q19 discrete velocity structure
* @param list - list of distributions to communicate
* @param start -  index to start parsing the list 
* @param count -  number of values to pack 
* @param sendbuf - memory buffer to hold values that will be sent
* @param dist - memory buffer to hold the distributions
* @param N - size of the distributions (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q19_Pack(int q, int *list, int start, int count, double *sendbuf, double *dist, int N);

/**
* \brief  Unpack D3Q19 distributions after communication
* @param q  - index for distribution based on D3Q19 discrete velocity structure
* @param list - list of distributions to communicate
* @param start -  index to start parsing the list 
* @param count -  number of values to unppack 
* @param recvbuf - memory buffer where recieved values have been stored
* @param dist - memory buffer to hold the distributions
* @param N - size of the distributions (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q19_Unpack(int q, int *list, int start, int count, double *recvbuf, double *dist, int N);

/**
* \brief Unpack D3Q7 distributions after communication
* @param q  - index for distribution based on D3Q19 discrete velocity structure
* @param list - list of distributions to communicate
* @param start -  index to start parsing the list 
* @param count -  number of values to unppack 
* @param recvbuf - memory buffer where recieved values have been stored
* @param dist - memory buffer to hold the distributions
* @param N - size of the distributions (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q7_Unpack(int q, int *list,  int start, int count, double *recvbuf, double *dist, int N);

/**
* \brief  Pack halo for scalar field to be prepare for communication
* @param list - list of distributions to communicate
* @param count -  number of values to ppack 
* @param sendbuf - memory buffer to pack values into
* @param Data - scalar field 
* @param N - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_Scalar_Pack(int *list, int count, double *sendbuf, double *Data, int N);

/**
* \brief  Pack halo for scalar field to be prepare for communication
* @param list - list of distributions to communicate
* @param count -  number of values to unppack 
* @param recvbuf - memory buffer where recieved values have been stored
* @param Data - scalar field 
* @param N - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_Scalar_Unpack(int *list, int count, double *recvbuf, double *Data, int N);

/**
* \brief  Unpack values and compute Shan-Chen type of gradient
* @param weight - weight value for gradient sum
* @param Cqx - contribution to x-part of gradient
* @param Cqy - contribution to y-part of gradient
* @param Cqz - contribution to z-part of gradient
* @param list - list of distributions to communicate
* @param start - location to start reading the list
* @param count -  number of values to unppack 
* @param recvbuf - memory buffer where recieved values have been stored
* @param phi - scalar field 
* @param grad - gradient 
* @param N - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_Gradient_Unpack(double weight, double Cqx, double Cqy, double Cqz, 
		int *list, int start, int count, double *recvbuf, double *phi, double *grad, int N);


/**
* \brief Initialize D3Q19 distributions
* @param Dist - D3Q19 distributions
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q19_Init(double *Dist, int Np);

/**
* \brief Compute momentum from D3Q19 distribution
* @param dist - D3Q19 distributions
* @param vel - memory buffer to store the momentum that is computed
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q19_Momentum(double *dist, double *vel, int Np);

/**
* \brief compute pressure from D3Q19 distribution
* @param dist - D3Q19 distributions
* @param press - memory buffer to store the pressure field that is computed
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q19_Pressure(double *dist, double *press, int Np);

// BGK MODEL
/**
* \brief BGK collision based on AA even access pattern for D3Q19 
* @param dist - D3Q19 distributions
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
* @param rlx - relaxation parameter
* @param Fx - force in x direction
* @param Fy - force in y direction
* @param Fz - force in z direction
*/
extern "C" void ScaLBL_D3Q19_AAeven_BGK(double *dist, int start, int finish, int Np, double rlx, double Fx, double Fy, double Fz);

/**
* \brief BGK collision based on AA odd access pattern for D3Q19 
* @param neighborList - neighbors based on D3Q19 lattice structure
* @param dist - D3Q19 distributions
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
* @param rlx - relaxation parameter
* @param Fx - force in x direction
* @param Fy - force in y direction
* @param Fz - force in z direction
*/
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

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColor(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
        double *Phi,double *GreySolidGrad, double *Poros,double *Perm,double *Vel,double *Pressure,
        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColor(int *d_neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *GreySolidGrad, double *Poros,double *Perm,double *Vel,double *Pressure, 
        double rhoA, double rhoB, double tauA, double tauB, double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColor_CP(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
        double *Phi, double *GreySolidW, double *GreySn, double *GreySw, double *GreyKn, double *GreyKw, double *Poros,double *Perm,double *Vel, double *Pressure,
        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, bool RecoloringOff, int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColor_CP(int *d_neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *GreySolidW, double *GreySn, double *GreySw, double *GreyKn, double *GreyKw, double *Poros, double *Perm,double *Vel,double *Pressure, 
        double rhoA, double rhoB, double tauA, double tauB, double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, bool RecoloringOff, int strideY, int strideZ, int start, int finish, int Np);

//extern "C" void ScaLBL_Update_GreyscalePotential(int *Map, double *Phi, double *Psi, double *Poro, double *Perm, double alpha, double W, 
//		int start, int finish, int Np);

// ION TRANSPORT MODEL

extern "C" void ScaLBL_D3Q7_AAodd_IonConcentration(int *neighborList, double *dist, double *Den, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_IonConcentration(double *dist, double *Den, int start, int finish, int Np);


extern "C" void ScaLBL_D3Q7_AAodd_Ion(int *neighborList, double *dist, double *Den, double *Velocity, double *ElectricField, 
                                      double Di, int zi, double rlx, double Vt, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_Ion(double *dist, double *Den, double *Velocity, double *ElectricField, 
                                       double Di, int zi, double rlx, double Vt, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_Ion_Init(double *dist, double *Den, double DenInit, int Np);
extern "C" void ScaLBL_D3Q7_Ion_Init_FromFile(double *dist, double *Den, int Np);

extern "C" void ScaLBL_D3Q7_Ion_ChargeDensity(double *Den, double *ChargeDensity, int IonValence, int ion_component, int start, int finish, int Np);

// LBM Poisson solver

/**
* \brief Poisson-Boltzmann collision based on AA odd access pattern for D3Q19 
* @param neighborList - neighbors based on D3Q19 lattice structure
* @param Map - mapping between sparse and dense representations
* @param dist - D3Q7 distributions
* @param Den_charge - charge density
* @param Psi - 
* @param ElectricField - electric field
* @param tau - relaxation time
* @param epsilon_LB - 
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q7_AAodd_Poisson(int *neighborList,int *Map, double *dist, double *Den_charge, double *Psi, double *ElectricField, double tau, double epsilon_LB,
        int start, int finish, int Np);

/**
* \brief Poisson-Boltzmann collision based on AA even access pattern for D3Q7 
* @param Map - mapping between sparse and dense representations
* @param dist - D3Q7 distributions
* @param Den_charge - charge density
* @param Psi - 
* @param ElectricField - electric field
* @param tau - relaxation time
* @param epsilon_LB - 
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q7_AAeven_Poisson(int *Map, double *dist, double *Den_charge, double *Psi, double *ElectricField, double tau,
					   double epsilon_LB, int start, int finish, int Np);
/**
* \brief Poisson-Boltzmann solver / solve electric potential based on AA odd access pattern for D3Q7 
* @param neighborList - neighbors based on D3Q19 lattice structure
* @param Map - mapping between sparse and dense representations
* @param dist - D3Q7 distributions
* @param Psi - 
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q7_AAodd_Poisson_ElectricPotential(int *neighborList,int *Map, double *dist, double *Psi, int start, int finish, int Np);

/**
* \brief Poisson-Boltzmann solver / solve electric potential based on AA odd access pattern for D3Q7
* @param Map - mapping between sparse and dense representations
* @param dist - D3Q7 distributions
* @param Psi - 
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q7_AAeven_Poisson_ElectricPotential(int *Map, double *dist, double *Psi, int start, int finish, int Np);

/**
* \brief Initialize Poisson-Boltzmann solver 
* @param Map - mapping between sparse and dense representations
* @param dist - D3Q7 distributions
* @param Psi - 
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q7_Poisson_Init(int *Map, double *dist, double *Psi, int start, int finish, int Np);

// LBM Stokes Model (adapted from MRT model)
extern "C" void ScaLBL_D3Q19_AAeven_StokesMRT(double *dist, double *Velocity, double *ChargeDensity, double *ElectricField, double rlx_setA, double rlx_setB, 
                double Gx, double Gy, double Gz,double rho0, double den_scale, double h, double time_conv, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_StokesMRT(int *neighborList, double *dist, double *Velocity, double *ChargeDensity, double *ElectricField, double rlx_setA, double rlx_setB, 
                double Gx, double Gy, double Gz, double rho0, double den_scale, double h, double time_conv,int start, int finish, int Np);

extern "C" void ScaLBL_PhaseField_InitFromRestart(double *Den, double *Aq, double *Bq, int start, int finish, int Np);

// MRT MODEL
/**
* \brief MRT collision based on AA even access pattern for D3Q19 
* @param dist - D3Q19 distributions
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
* @param rlx_setA - relaxation parameter for viscous modes
* @param rlx_setB - relaxation parameter for non-viscous modes
* @param Fx - force in x direction
* @param Fy - force in y direction
* @param Fz - force in z direction
*/
extern "C" void ScaLBL_D3Q19_AAeven_MRT(double *dist, int start, int finish, int Np, double rlx_setA, double rlx_setB, double Fx,
		double Fy, double Fz);

/**
* \brief MRT collision based on AA even access pattern for D3Q19 
* @param neighborList - index of neighbors based on D3Q19 lattice structure
* @param dist - D3Q19 distributions
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
* @param rlx_setA - relaxation parameter for viscous modes
* @param rlx_setB - relaxation parameter for non-viscous modes
* @param Fx - force in x direction
* @param Fy - force in y direction
* @param Fz - force in z direction
*/
extern "C" void ScaLBL_D3Q19_AAodd_MRT(int *neighborList, double *dist, int start, int finish, int Np,
		double rlx_setA, double rlx_setB, double Fx, double Fy, double Fz);

// COLOR MODEL
/**
* \brief Color model collision based on AA even access pattern for D3Q19 
* @param Map - mapping between sparse and dense data structures
* @param dist - D3Q19 distributions
* @param Aq - D3Q7 distribution for component A
* @param Bq - D3Q7 distribution for component B
* @param Den - density field
* @param Phi - phase indicator field
* @param Vel - velocity field
* @param rhoA - density of component A
* @param rhoB - density of component B
* @param tauA - relaxation time for component A
* @param tauB - relaxation time for component B
* @param alpha - parameter to control interfacial tension
* @param beta - parameter to control interface width
* @param Fx - force in x direction
* @param Fy - force in y direction
* @param Fz - force in z direction
* @param strideY - stride in y-direction for gradient computation
* @param strideZ - stride in z-direction for gradient computation
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q19_AAeven_Color(int *Map, double *dist, double *Aq, double *Bq, double *Den, double *Phi,
		double *Vel, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np);

/**
* \brief Color model collision based on AA even access pattern for D3Q19 
* @param NeighborList - neighbors based on D3Q19 lattice structure
* @param Map - mapping between sparse and dense data structures
* @param dist - D3Q19 distributions
* @param Aq - D3Q7 distribution for component A
* @param Bq - D3Q7 distribution for component B
* @param Den - density field
* @param Phi - phase indicator field
* @param Vel - velocity field
* @param rhoA - density of component A
* @param rhoB - density of component B
* @param tauA - relaxation time for component A
* @param tauB - relaxation time for component B
* @param alpha - parameter to control interfacial tension
* @param beta - parameter to control interface width
* @param Fx - force in x direction
* @param Fy - force in y direction
* @param Fz - force in z direction
* @param strideY - stride in y-direction for gradient computation
* @param strideZ - stride in z-direction for gradient computation
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q19_AAodd_Color(int *NeighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *Vel, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np);

/**
* \brief Compute phase field based on AA odd access pattern for D3Q19 
* @param NeighborList - neighbors based on D3Q19 lattice structure
* @param Map - mapping between sparse and dense data structures
* @param Aq - D3Q7 distribution for component A
* @param Bq - D3Q7 distribution for component B
* @param Den - density field
* @param Phi - phase indicator field
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q7_AAodd_PhaseField(int *NeighborList, int *Map, double *Aq, double *Bq, 
			double *Den, double *Phi, int start, int finish, int Np);

/**
* \brief Compute phase field based on AA even access pattern for D3Q19 
* @param Map - mapping between sparse and dense data structures
* @param Aq - D3Q7 distribution for component A
* @param Bq - D3Q7 distribution for component B
* @param Den - density field
* @param Phi - phase indicator field
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_D3Q7_AAeven_PhaseField(int *Map, double *Aq, double *Bq, double *Den, double *Phi, 
			int start, int finish, int Np);

/**
* \brief Initialize phase field for color model
* @param Map - mapping between sparse and dense data structures
* @param Phi - phase indicator field
* @param Den - density field
* @param Aq - D3Q7 distribution for component A
* @param Bq - D3Q7 distribution for component B
* @param start - lattice node to start loop
* @param finish - lattice node to finish loop
* @param Np - size of local sub-domain (derived from Domain structure)
*/
extern "C" void ScaLBL_PhaseField_Init(int *Map, double *Phi, double *Den, double *Aq, double *Bq, int start, int finish, int Np);


extern "C" void ScaLBL_D3Q7_AAodd_Color(int *neighborList, int *Map, double *Aq, double *Bq, double *Den, 
		double *Phi, double *ColorGrad, double *Vel, double rhoA, double rhoB, double beta, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_Color(int *Map, double *Aq, double *Bq, double *Den, 
		double *Phi, double *ColorGrad, double *Vel, double rhoA, double rhoB, double beta, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_Gradient(int *Map, double *Phi, double *ColorGrad, int start, int finish, int Np, int Nx, int Ny, int Nz);

extern "C" void ScaLBL_D3Q19_MixedGradient(int *Map, double *Phi, double *Gradient, int start, int finish, int Np, int Nx, int Ny, int Nz);


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

// FREE ENERGY LEE MODEL

extern "C" void ScaLBL_D3Q19_FreeLeeModel_TwoFluid_Init(double *gqbar, double *mu_phi, double *ColorGrad, double Fx, double Fy, double Fz, int Np);

extern "C" void ScaLBL_D3Q19_FreeLeeModel_SingleFluid_Init(double *gqbar, double Fx, double Fy, double Fz, int Np);

extern "C" void ScaLBL_FreeLeeModel_PhaseField_Init(int *Map, double *Phi, double *Den, double *hq, double *ColorGrad, 
                                                    double rhonA, double rhoB, double tauM, double W, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAodd_FreeLeeModel_PhaseField(int *neighborList, int *Map, double *hq, double *Den, double *Phi, 
                                                         double rhoA, double rhoB, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_FreeLeeModel_PhaseField(int *Map, double *hq, double *Den, double *Phi, 
			                                               double rhoA, double rhoB, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAodd_FreeLee_PhaseField(int *neighborList, int *Map, double *hq, double *Den, double *Phi, double *ColorGrad, double *Vel,
                                                          double rhoA, double rhoB, double tauM, double W, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_FreeLee_PhaseField( int *Map, double *hq, double *Den, double *Phi, double *ColorGrad, double *Vel,
		double rhoA, double rhoB, double tauM, double W, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q7_ComputePhaseField(int *Map,  double *hq, double *Den, double *Phi, double rhoA, double rhoB, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_FreeLeeModel(int *neighborList, int *Map, double *dist, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad, 
                                                double rhoA, double rhoB, double tauA, double tauB, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_FreeLeeModel(int *Map, double *dist, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
                                                double rhoA, double rhoB, double tauA, double tauB, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined(int *neighborList, int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
                                                double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined(int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
                                                double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_FreeLeeModel_Combined_HigherOrder(int *neighborList, int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
                                                double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_FreeLeeModel_Combined_HigherOrder(int *Map, double *dist, double *hq, double *Den,	double *Phi, double *mu_phi, double *Vel, double *Pressure, double *ColorGrad,
                                                double rhoA, double rhoB, double tauA, double tauB, double tauM, double kappa, double beta, double W, double Fx, double Fy, double Fz, 
                                                int strideY, int strideZ, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAodd_FreeLeeModel_SingleFluid_BGK(int *neighborList, double *dist, double *Vel, double *Pressure,  
                                                                double tau, double rho0, double Fx, double Fy, double Fz, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q19_AAeven_FreeLeeModel_SingleFluid_BGK(double *dist, double *Vel, double *Pressure, 
                                                                 double tau, double rho0, double Fx, double Fy, double Fz, int start, int finish, int Np);

extern "C" void ScaLBL_D3Q9_MGTest(int *Map, double *Phi,double *ColorGrad,int strideY, int strideZ, int start, int finish, int Np);

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

extern "C" void ScaLBL_D3Q19_Reflection_BC_z(int *list, double *dist, int count, int Np);

extern "C" void ScaLBL_D3Q19_Reflection_BC_Z(int *list, double *dist, int count, int Np);

extern "C" void ScaLBL_D3Q7_Reflection_BC_z(int *list, double *dist, int count, int Np);

extern "C" void ScaLBL_D3Q7_Reflection_BC_Z(int *list, double *dist, int count, int Np);

extern "C" void ScaLBL_SetSlice_z(double *Phi, double value, int Nx, int Ny, int Nz, int Slice);

extern "C" void ScaLBL_CopySlice_z(double *Phi, int Nx, int Ny, int Nz, int Source, int Destination);

extern "C" void ScaLBL_Solid_Dirichlet_D3Q7(double *dist,double *BoundaryValue,int *BounceBackDist_list,int *BounceBackSolid_list,int N);

extern "C" void ScaLBL_Solid_Neumann_D3Q7(double *dist,double *BoundaryValue,int *BounceBackDist_list,int *BounceBackSolid_list,int N);

extern "C" void ScaLBL_Solid_SlippingVelocityBC_D3Q19(double *dist, double *zeta_potential, double *ElectricField, double *SolidGrad,
                                               double epsilon_LB, double tau, double rho0,double den_scale, double h, double time_conv,
                                               int *BounceBackDist_list, int *BounceBackSolid_list, int *FluidBoundary_list,
                                               double *lattice_weight, float *lattice_cx, float *lattice_cy, float *lattice_cz,
                                               int count, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z(int *list, double *dist, double Vin, int count, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z(int *list, double *dist, double Vout, int count, int Np);

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z(int *d_neighborList, int *list, double *dist, double Vin, int count, int Np);

extern "C" void ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(int *d_neighborList, int *list, double *dist, double Vout, int count, int Np);

extern "C" void ScaLBL_Poisson_D3Q7_BC_z(int *list, int *Map, double *Psi, double Vin, int count);

extern "C" void ScaLBL_Poisson_D3Q7_BC_Z(int *list, int *Map, double *Psi, double Vout, int count);

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z(int *list, double *dist, double Cin, int count, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z(int *list, double *dist, double Cout, int count, int Np);

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z(int *d_neighborList, int *list, double *dist, double Cin, int count, int Np);

extern "C" void ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, int count, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_z(int *list, double *dist, double Cin, double tau, double *VelocityZ, int count, int Np);
extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_Z(int *list, double *dist, double Cout, double tau, double *VelocityZ, int count, int Np);
extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_z(int *d_neighborList, int *list, double *dist, double Cin, double tau, double *VelocityZ, int count, int Np);
extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, double tau, double *VelocityZ, int count, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_z(int *list, double *dist, double Cin, double tau, double *VelocityZ, int count, int Np);
extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_Z(int *list, double *dist, double Cout, double tau, double *VelocityZ, int count, int Np);
extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_z(int *d_neighborList, int *list, double *dist, double Cin, double tau, double *VelocityZ, int count, int Np);
extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, double tau, double *VelocityZ, int count, int Np);

extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_z(int *list, double *dist, double Cin, double tau, double *VelocityZ,double *ElectricField,double Di,double zi,double Vt,int count,int Np);
extern "C" void ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_Z(int *list, double *dist, double Cout, double tau, double *VelocityZ,double *ElectricField,double Di,double zi,double Vt,int count,int Np);
extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_z(int *d_neighborList, int *list, double *dist, double Cin, double tau, double *VelocityZ,double *ElectricField,double Di,double zi,double Vt, int count, int Np);
extern "C" void ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_Z(int *d_neighborList, int *list, double *dist, double Cout, double tau, double *VelocityZ,double *ElectricField,double Di,double zi,double Vt, int count, int Np);

/**
 * \class ScaLBL_Communicator
 *
 * @brief Generalized communication routines for lattice Boltzmann methods.
 *
*/
class ScaLBL_Communicator{
public:
	//......................................................................................
        /**
        *\brief Constructor
	* @param Dm - Domain information
        */
        ScaLBL_Communicator(std::shared_ptr <Domain> Dm);

        /**
        *\brief Destructor
        */
	~ScaLBL_Communicator();
	//......................................................................................
	unsigned long int CommunicationCount,SendCount,RecvCount;
	int Nx,Ny,Nz,N;
	int n_bb_d3q7, n_bb_d3q19; 
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
	
	double GetPerformance(int *NeighborList, double *fq, int Np);
	int MemoryOptimizedLayoutAA(IntArray &Map, int *neighborList, signed char *id, int Np, int width);
	void Barrier(){
		ScaLBL_DeviceBarrier();
		MPI_COMM_SCALBL.barrier();
	};
	void SendD3Q19AA(double *dist);
	void RecvD3Q19AA(double *dist);
	void SendD3Q7AA(double *fq, int Component);
	void RecvD3Q7AA(double *fq, int Component);
	void BiSendD3Q7AA(double *Aq, double *Bq);
	void BiRecvD3Q7AA(double *Aq, double *Bq);
	void TriSendD3Q7AA(double *Aq, double *Bq, double *Cq);
	void TriRecvD3Q7AA(double *Aq, double *Bq, double *Cq);
	void SendHalo(double *data);
	void RecvHalo(double *data);
	void RecvGrad(double *Phi, double *Gradient);
	void RegularLayout(IntArray map, const double *data, DoubleArray &regdata);
	void SetupBounceBackList(IntArray &Map, signed char *id, int Np, bool SlippingVelBC=false);
    void SolidDirichletD3Q7(double *fq, double *BoundaryValue);
    void SolidNeumannD3Q7(double *fq, double *BoundaryValue);
    void SolidSlippingVelocityBCD3Q19(double *fq, double *zeta_potential, double *ElectricField, double *SolidGrad,
                                      double epslion_LB, double tau, double rho0, double den_scale,double h, double time_conv);

    // Routines to set boundary conditions
    void Color_BC_z(int *Map, double *Phi, double *Den, double vA, double vB);
    void Color_BC_Z(int *Map, double *Phi, double *Den, double vA, double vB);
    void D3Q19_Pressure_BC_z(int *neighborList, double *fq, double din, int time);
    void D3Q19_Pressure_BC_Z(int *neighborList, double *fq, double dout, int time);
    void D3Q19_Reflection_BC_z(double *fq);
    void D3Q19_Reflection_BC_Z(double *fq);
    double D3Q19_Flux_BC_z(int *neighborList, double *fq, double flux, int time);
    void D3Q7_Poisson_Potential_BC_z(int *neighborList, double *fq, double Vin, int time);
    void D3Q7_Poisson_Potential_BC_Z(int *neighborList, double *fq, double Vout, int time);
    void Poisson_D3Q7_BC_z(int *Map, double *Psi, double Vin);
    void Poisson_D3Q7_BC_Z(int *Map, double *Psi, double Vout);
    void D3Q7_Ion_Concentration_BC_z(int *neighborList, double *fq, double Cin, int time);
    void D3Q7_Ion_Concentration_BC_Z(int *neighborList, double *fq, double Cout, int time);
    void D3Q7_Ion_Flux_Diff_BC_z(int *neighborList, double *fq, double Cin, double tau, double *VelocityZ, int time);
    void D3Q7_Ion_Flux_Diff_BC_Z(int *neighborList, double *fq, double Cout, double tau, double *VelocityZ, int time);
    void D3Q7_Ion_Flux_DiffAdvc_BC_z(int *neighborList, double *fq, double Cin, double tau, double *VelocityZ, int time);
    void D3Q7_Ion_Flux_DiffAdvc_BC_Z(int *neighborList, double *fq, double Cout, double tau, double *VelocityZ, int time);
    void D3Q7_Ion_Flux_DiffAdvcElec_BC_z(int *neighborList,double *fq,double Cin,double tau,double *VelocityZ,double *ElectricField_Z,double Di,double zi,double Vt, int time);
    void D3Q7_Ion_Flux_DiffAdvcElec_BC_Z(int *neighborList,double *fq,double Cout,double tau,double *VelocityZ,double *ElectricField_Z,double Di,double zi,double Vt, int time);
    void GreyscaleSC_BC_z(int *Map, double *DenA, double *DenB, double vA, double vB);
    void GreyscaleSC_BC_Z(int *Map, double *DenA, double *DenB, double vA, double vB);
    void GreyscaleSC_Pressure_BC_z(int *neighborList, double *fqA, double *fqB, double dinA, double dinB, int time);
    void GreyscaleSC_Pressure_BC_Z(int *neighborList, double *fqA, double *fqB, double doutA, double doutB, int time);
    // Debugging and unit testing functions
    void PrintD3Q19();

private:
	void D3Q19_MapRecv(int Cqx, int Cqy, int Cqz, const int *list,  int start, int count, int *d3q19_recvlist);

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
	int *bb_dist;
	int *bb_interactions;
    int *fluid_boundary;
    double *lattice_weight;
    float *lattice_cx, *lattice_cy, *lattice_cz;
	//......................................................................................

};


#endif
