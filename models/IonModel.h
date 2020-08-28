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
	void ReadParams(string filename,int num_iter,int num_iter_Stokes,double time_conv_Stokes);
	void ReadParams(std::shared_ptr<Database> db0);
	void SetDomain();
	void ReadInput();
	void Create();
	void Initialize();
	void Run(double *Velocity, double *ElectricField);
    void getIonConcentration(int timestep);
	
	//bool Restart,pBC;
	int timestep,timestepMax;
	int BoundaryCondition;
	int BoundaryConditionSolid;
    double h;//domain resolution, unit [um/lu]
    double time_conv;
    double kb,electron_charge,T,Vt;
    double k2_inv;
    double tolerance;
    double Ex,Ey,Ez;
	
	int number_ion_species;
    vector<double> IonDiffusivity;//User input unit [m^2/sec]
    vector<int> IonValence;
    vector<double> IonConcentration;//unit [mol/m^3]
    //vector<double> deltaT;
	vector<double> tau;
	
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

private:
	MPI_Comm comm;
	
	// filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
   
    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);    	
    void AssignSolidBoundary(double *ion_solid);
};
