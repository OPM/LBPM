/*
 * Ion transporte LB Model
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <vector>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/MPI_Helpers.h"
#include "analysis/Minkowski.h"
#include "ProfilerApp.h"

class ScaLBL_IonModel{
public:
	ScaLBL_IonModel(int RANK, int NP, MPI_Comm COMM);
	~ScaLBL_IonModel();	
	
	// functions in they should be run
	void ReadParams(string filename,vector<int> &num_iter);
	void ReadParams(string filename);
	void ReadParams(std::shared_ptr<Database> db0);
	void SetDomain();
	void ReadInput();
	void Create();
	void Initialize();
	void Run(double *Velocity, double *ElectricField);
    void getIonConcentration(int timestep);
    void DummyFluidVelocity();
    void DummyElectricField();
    double CalIonDenConvergence(vector<double> &ci_avg_previous);

	//bool Restart,pBC;
	int timestep;
    vector<int> timestepMax;
	int BoundaryCondition;
	int BoundaryConditionSolid;
    double h;//domain resolution, unit [um/lu]
    double kb,electron_charge,T,Vt;
    double k2_inv;
    double tolerance;
    double fluidVelx_dummy,fluidVely_dummy,fluidVelz_dummy;
    double Ex_dummy,Ey_dummy,Ez_dummy;
	
	int number_ion_species;
    vector<double> IonDiffusivity;//User input unit [m^2/sec]
    vector<int> IonValence;
    vector<double> IonConcentration;//unit [mol/m^3]
    vector<double> Cin;//unit [mol/m^3]
    vector<double> Cout;//unit [mol/m^3]
	vector<double> tau;
	vector<double> time_conv;
	
	int Nx,Ny,Nz,N,Np;
	int rank,nprocx,nprocy,nprocz,nprocs;
	double Lx,Ly,Lz;

	std::shared_ptr<Domain> Dm;   // this domain is for analysis
	std::shared_ptr<Domain> Mask; // this domain is for lbm
	std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
    // input database
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> domain_db;
    std::shared_ptr<Database> ion_db;

    IntArray Map;
    DoubleArray Distance;
    int *NeighborList;
    double *fq;
    double *Ci; 
    double *ChargeDensity; 
    double *IonSolid;
    double *FluidVelocityDummy;
    double *ElectricFieldDummy;

private:
	MPI_Comm comm;
	
	// filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
    char OutputFilename[200];
   
    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);    	
    void AssignSolidBoundary(double *ion_solid);
    void AssignIonConcentration_FromFile(double *Ci,const vector<std::string> &File_ion);
    void IonConcentration_LB_to_Phys(DoubleArray &Den_reg);
};
