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
#include <cmath>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/MPI.h"
#include "analysis/Minkowski.h"
#include "ProfilerApp.h"

#define _USE_MATH_DEFINES
#ifndef ScaLBL_POISSON_INC
#define ScaLBL_POISSON_INC

class ScaLBL_Poisson{
public:
	ScaLBL_Poisson(int RANK, int NP, const Utilities::MPI& COMM);
	~ScaLBL_Poisson();	
	
	// functions in they should be run
	void ReadParams(string filename);
	void ReadParams(std::shared_ptr<Database> db0);
	void SetDomain();
	void ReadInput();
	void Create();
	void Initialize(double time_conv_from_Study);
	void Run(double *ChargeDensity,int timestep_from_Study);
    void getElectricPotential(DoubleArray &ReturnValues);
    void getElectricPotential_debug(int timestep);
    void getElectricField(DoubleArray &Values_x, DoubleArray &Values_y, DoubleArray &Values_z);
    void getElectricField_debug(int timestep);
    void DummyChargeDensity();//for debugging

	//bool Restart,pBC;
	int timestep,timestepMax;
    int analysis_interval;
	int BoundaryConditionInlet;
	int BoundaryConditionOutlet;
    int BoundaryConditionSolid;
	double tau;
	double tolerance;
    double k2_inv;
    double epsilon0,epsilon0_LB,epsilonR,epsilon_LB;
    double Vin, Vout;
    double chargeDen_dummy;//for debugging
    bool WriteLog;
    double Vin0,freqIn,t0_In,Vin_Type;
    double Vout0,freqOut,t0_Out,Vout_Type;
    bool   TestPeriodic;
    double TestPeriodicTime;//unit: [sec]
    double TestPeriodicTimeConv; //unit [sec/lt]
    double TestPeriodicSaveInterval; //unit [sec]
	
	int Nx,Ny,Nz,N,Np;
	int rank,nprocx,nprocy,nprocz,nprocs;
	double Lx,Ly,Lz;
    double h;//image resolution
    double time_conv;//phys to LB time converting factor; unit=[sec/lt]

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
	Utilities::MPI comm;
	
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
    void SolveElectricPotentialAAodd(int timestep_from_Study);
    void SolveElectricPotentialAAeven(int timestep_from_Study);
    //void SolveElectricField();
    void SolvePoissonAAodd(double *ChargeDensity);
    void SolvePoissonAAeven(double *ChargeDensity);
    void getConvergenceLog(int timestep,double error);
    double getBoundaryVoltagefromPeriodicBC(double V0,double freq,double t0,int V_type,int time_step);
    
};
#endif
