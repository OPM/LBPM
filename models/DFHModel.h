/*
Implementation of color lattice boltzmann model
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/Communication.h"
#include "analysis/TwoPhase.h"
#include "analysis/runAnalysis.h"
#include "common/MPI_Helpers.h"
#include "ProfilerApp.h"
#include "threadpool/thread_pool.h"

class ScaLBL_DFHModel{
public:
	ScaLBL_DFHModel(int RANK, int NP, MPI_Comm COMM);
	~ScaLBL_DFHModel();	
	
	// functions in they should be run
	void ReadParams(string filename);
	void ReadParams(std::shared_ptr<Database> db0);
	void SetDomain();
	void ReadInput();
	void Create();
	void Initialize();
	void AssignSolidPotential();
	void Run();
	void WriteDebug();
	
	bool Restart,pBC;
	int timestep,timestepMax;
	int BoundaryCondition;
	double tauA,tauB,rhoA,rhoB,alpha,beta;
	double Fx,Fy,Fz,flux;
	double din,dout,inletA,inletB,outletA,outletB;
	
	int Nx,Ny,Nz,N,Np;
	int rank,nprocx,nprocy,nprocz,nprocs;
	double Lx,Ly,Lz;

	std::shared_ptr<Domain> Dm;   // this domain is for analysis
	std::shared_ptr<Domain> Mask; // this domain is for lbm
	std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
    std::shared_ptr<TwoPhase> Averages;
    
    // input database
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> domain_db;
    std::shared_ptr<Database> color_db;
    std::shared_ptr<Database> analysis_db;

    IntArray Map;
    char *id;
    int *NeighborList;
    int *dvcMap;
    double *fq, *Aq, *Bq;
    double *Den, *Phi;
    double *SolidPotential;
    double *Velocity;
    double *Gradient;
    double *Pressure;
		
private:
	MPI_Comm comm;
    
	int dist_mem_size;
	int neighborSize;
	// filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
   
    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);
    void AssignComponentLabels(double *phase);
    	
};

