/*
 * Ion transporte LB Model
 */

#ifndef ScaLBL_IonModel_INC
#define ScaLBL_IonModel_INC

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
#include "common/MPI.h"
#include "analysis/Minkowski.h"
#include "ProfilerApp.h"

class ScaLBL_IonModel{
public:
	ScaLBL_IonModel(int RANK, int NP, const Utilities::MPI& COMM);
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
    void getIonConcentration(DoubleArray &IonConcentration, const size_t ic);
    void getIonConcentration_debug(int timestep);
    void DummyFluidVelocity();
    void DummyElectricField();
    double CalIonDenConvergence(vector<double> &ci_avg_previous);

	//bool Restart,pBC;
	int timestep;
    vector<int> timestepMax;
	int BoundaryConditionSolid;
    double h;//domain resolution, unit [um/lu]
    double kb,electron_charge,T,Vt;
    double k2_inv;
    double tolerance;
    double fluidVelx_dummy,fluidVely_dummy,fluidVelz_dummy;
    double Ex_dummy,Ey_dummy,Ez_dummy;
	
	size_t number_ion_species;
	vector<int> BoundaryConditionInlet;
	vector<int> BoundaryConditionOutlet;
    vector<double> IonDiffusivity;//User input unit [m^2/sec]
    vector<int> IonValence;
    vector<double> IonConcentration;//unit [mol/m^3]
    vector<double> Cin;//inlet boundary value, can be either concentration [mol/m^3] or flux [mol/m^2/sec]
    vector<double> Cout;//outlet boundary value, can be either concentration [mol/m^3] or flux [mol/m^2/sec]
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
	Utilities::MPI comm;
	
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
#endif
