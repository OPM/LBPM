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

class ScaLBL_ColorModel{
public:
	ScaLBL_ColorModel(string filename);
	~ScaLBL_ColorModel();	
	
	// functions in they should be run
	void ReadInput();
	void Create();
	void Initialize();
	void Run();
	
	bool Restart;
	int timestep,timestepMax;
	double tauA,tauB,rhoA,rhoB,alpha,beta;
	double Fx,Fy,Fz,flux;
	double din,dout,inletA,inletB,outletA,outletB;
	
	int Nx,Ny,Nz,N,Np;
	int nprocx,nprocy,nprocz,BC;
	double Lx,Ly,Lz;
		
private:
	MPI_Comm comm;
	Database db;
	Database domain_db;
	Database color_db;
	Database analysis_db;
      
	Domain Dm;   // this domain is for analysis
	Domain Mask; // this domain is for lbm

    ScaLBL_Communicator ScaLBL_Comm;

    std::shared_ptr<TwoPhase> Averages;

	// filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
    
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
    
    int rank;
    	
};

