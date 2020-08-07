/*
 * Multi-relaxation time LBM Model
 */
#include <stdio.h>
#include <stdlib.h>
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
	
	//bool Restart,pBC;
	int timestep,timestepMax;
    int analysis_interval;
	int BoundaryCondition;
	double tau;
	double tolerance;
    double k2_inv,deltaT;
    double epsilon0,epsilon0_LB,epsilonR,epsilon_LB;
	
	int Nx,Ny,Nz,N,Np;
	int rank,nprocx,nprocy,nprocz,nprocs;
	double Lx,Ly,Lz;
    double h;//image resolution

	std::shared_ptr<Domain> Dm;   // this domain is for analysis
	std::shared_ptr<Domain> Mask; // this domain is for lbm
	std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
    // input database
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> domain_db;
    std::shared_ptr<Database> electric_db;

    IntArray Map;
    DoubleArray Distance;
    DoubleArray Psi_host;
    int *NeighborList;
    double *fq;
    double *Psi; 
    double *ElectricField;

private:
	MPI_Comm comm;
	
	// filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
   
    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);    	
};
