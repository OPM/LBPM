/*
 * Multi-relaxation time LBM Model
 */
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/MPI_Helpers.h"
#include "analysis/Minkowski.h"
#include "ProfilerApp.h"

class ScaLBL_Poisson{
public:
	ScaLBL_Poisson(int RANK, int NP, MPI_Comm COMM);
	~ScaLBL_Poisson();	
	
	// functions in they should be run
	void ReadParams(string filename);
	void ReadParams(std::shared_ptr<Database> db0);
	void SetDomain();
	void ReadInput();
	void Create();
	void Initialize();
	void Run(double *ChargeDensity);
    void getElectricPotential(int timestep);
    void getElectricPotential_debug(int timestep);
    void getElectricField(int timestep);
    void getElectricField_debug(int timestep);
    void DummyChargeDensity();//for debugging

	//bool Restart,pBC;
	int timestep,timestepMax;
    int analysis_interval;
	int BoundaryCondition;
    int BoundaryConditionSolid;
	double tau;
	double tolerance;
    double k2_inv;
    double epsilon0,epsilon0_LB,epsilonR,epsilon_LB;
    double Vin, Vout;
    double chargeDen_dummy;//for debugging
    short WriteLog;
	
	int Nx,Ny,Nz,N,Np;
	int rank,nprocx,nprocy,nprocz,nprocs;
	double Lx,Ly,Lz;
    double h;//image resolution

	std::shared_ptr<Domain> Dm;   // this domain is for analysis
	std::shared_ptr<Domain> Mask; // this domain is for lbm
	std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
	std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm_Regular;
    // input database
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> domain_db;
    std::shared_ptr<Database> electric_db;

    IntArray Map;
    DoubleArray Distance;
    DoubleArray Psi_host;
    int *NeighborList;
    int *dvcMap;
    //signed char *dvcID;
    double *fq;
    double *Psi; 
    double *ElectricField;
    double *ChargeDensityDummy;// for debugging

private:
	MPI_Comm comm;
	
	// filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
    char OutputFilename[200];
	FILE *TIMELOG;
   
    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);    	
    void AssignSolidBoundary(double *poisson_solid);
    void Potential_Init(double *psi_init);
    void ElectricField_LB_to_Phys(DoubleArray &Efield_reg);
    void SolveElectricPotentialAAodd();
    void SolveElectricPotentialAAeven();
    //void SolveElectricField();
    void SolvePoissonAAodd(double *ChargeDensity);
    void SolvePoissonAAeven(double *ChargeDensity);
    void getConvergenceLog(int timestep,double error);
    
};
