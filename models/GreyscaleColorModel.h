/*
Implementation of multicomponent greyscale lattice boltzmann model
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/Communication.h"
#include "common/MPI.h"
#include "common/Database.h"
#include "common/ScaLBL.h"
#include "ProfilerApp.h"
#include "threadpool/thread_pool.h"

class ScaLBL_GreyscaleColorModel{
public:
	ScaLBL_GreyscaleColorModel(int RANK, int NP, MPI_Comm COMM);
	~ScaLBL_GreyscaleColorModel();	
	
	// functions in they should be run
	void ReadParams(string filename);
	void ReadParams(std::shared_ptr<Database> db0);
	void SetDomain();
	void ReadInput();
	void Create();
	void Initialize();
	void Run();
	void WriteDebug();
	void VelocityField();
	
	bool Restart,pBC;
	int timestep,timestepMax;
	int BoundaryCondition;
	double tauA,tauB;
    double tauA_eff,tauB_eff;
    double rhoA,rhoB;
	double tolerance;
	double Fx,Fy,Fz,flux;
	double din,dout;
    double dp;//solid particle diameter, unit in voxel
    double GreyPorosity;
    double Gsc;
	
	int Nx,Ny,Nz,N,Np;
	int rank,nprocx,nprocy,nprocz,nprocs;
	double Lx,Ly,Lz;

	std::shared_ptr<Domain> Dm;   // this domain is for analysis
	std::shared_ptr<Domain> Mask; // this domain is for lbm
	std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
    
    // input database
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> domain_db;
    std::shared_ptr<Database> greyscaleColor_db;
    std::shared_ptr<Database> analysis_db;
    std::shared_ptr<Database> vis_db;

    signed char *id;    
	int *NeighborList;
	double *fq,*Aq,*Bq;
    double *Den;
	double *Permeability;//grey voxel permeability
	double *Porosity;
	double *Velocity;
	double *Pressure_dvc;
    double *SolidForce;
    double *DenGradA;
    double *DenGradB;
    IntArray Map;
    DoubleArray SignDist;
    DoubleArray Velocity_x;
    DoubleArray Velocity_y;
    DoubleArray Velocity_z;
    DoubleArray PorosityMap;
    DoubleArray Pressure;
		
private:
	MPI_Comm comm;
    
	int dist_mem_size;
	int neighborSize;
	// filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
   
    void AssignComponentLabels(double *Porosity, double *Permeablity, double *SolidPotential);
    void AssignSolidForce(double *SolidPotential, double *SolidForce);
    void DensityField_Init();

};

