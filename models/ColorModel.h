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
#include "common/MPI.h"
#include "ProfilerApp.h"
#include "threadpool/thread_pool.h"


#ifndef ScaLBL_ColorModel_INC
#define ScaLBL_ColorModel_INC

class ScaLBL_ColorModel{
public:
	ScaLBL_ColorModel(int RANK, int NP, const Utilities::MPI& COMM);
	~ScaLBL_ColorModel();	
	
	// functions in they should be run
	void ReadParams(string filename);
	void ReadParams(std::shared_ptr<Database> db0);
	void SetDomain();
	void ReadInput();
	void Create();
	void Initialize();
	void Run();
	double Run(int returntime);
	void WriteDebug();
	void getPhaseField(DoubleArray &f);
	
	bool Restart,pBC;
	bool REVERSE_FLOW_DIRECTION;
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
	std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm_Regular;
    //std::shared_ptr<TwoPhase> Averages;
    std::shared_ptr<SubPhase> Averages;
    
    // input database
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> domain_db;
    std::shared_ptr<Database> color_db;
    std::shared_ptr<Database> analysis_db;
    std::shared_ptr<Database> vis_db;

    IntArray Map;
    signed char *id;    
	int *NeighborList;
	int *dvcMap;
	double *fq, *Aq, *Bq;
	double *Den, *Phi;
	double *ColorGrad;
	double *Velocity;
	double *Pressure;

	void AssignComponentLabels(double *phase);
		
private:
	Utilities::MPI comm;

	int dist_mem_size;
	int neighborSize;
	// filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
   
    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);
    double ImageInit(std::string filename);
    double MorphInit(const double beta, const double morph_delta);
    double SeedPhaseField(const double seed_water_in_oil);
    double MorphOpenConnected(double target_volume_change);
};

class FlowAdaptor{
public:
	FlowAdaptor(ScaLBL_ColorModel &M);
	~FlowAdaptor();
	double MoveInterface(ScaLBL_ColorModel &M);
	double ImageInit(ScaLBL_ColorModel &M, std::string Filename);
	double UpdateFractionalFlow(ScaLBL_ColorModel &M);
	void Flatten(ScaLBL_ColorModel &M);
	DoubleArray phi;
	DoubleArray phi_t;
private:
	int Nx, Ny, Nz;
	int timestep;
	int timestep_previous;
};
#endif

