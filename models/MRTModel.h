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

#include "common/Communication.h"
#include "common/MPI_Helpers.h"
#include "ProfilerApp.h"

class ScaLBL_MRTModel{
public:
	ScaLBL_MRTModel(int RANK, int NP, MPI_Comm COMM);
	~ScaLBL_MRTModel();	
	
	// functions in they should be run
	void ReadParams(string filename);
	void ReadParams(std::shared_ptr<Database> db0);
	void SetDomain();
	void ReadInput();
	void Create();
	void Initialize();
	void Run();
	void VelocityField(double *Vz);
	
	bool Restart,pBC;
	int timestep,timestepMax;
	int BoundaryCondition;
	double tau,mu;
	double Fx,Fy,Fz,flux;
	double din,dout;
	
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
    int *NeighborList;
    double *fq;
    double *Velocity;
    double *Pressure;
		
private:
	MPI_Comm comm;
    
	// filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
   
    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);    	
};