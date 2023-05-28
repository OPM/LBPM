/*
Implementation of Lee et al JCP 2016 lattice boltzmann model
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/Communication.h"
#include "common/MPI.h"
#include "ProfilerApp.h"
#include "threadpool/thread_pool.h"
#include "common/ScaLBL.h"
#include "common/WideHalo.h"

#ifndef ScaLBL_FreeLeeModel_INC
#define ScaLBL_FreeLeeModel_INC

class ScaLBL_FreeLeeModel {
public:
    ScaLBL_FreeLeeModel(int RANK, int NP, const Utilities::MPI &COMM);
    ~ScaLBL_FreeLeeModel();

    // functions in they should be run
    void ReadParams(string filename);
    void ReadParams(std::shared_ptr<Database> db0);
    void SetDomain();
    void ReadInput();
    void Create_TwoFluid();
    void Initialize_TwoFluid();
    double Run_TwoFluid(int returntime);

    void WriteDebug_TwoFluid();
    void Create_SingleFluid();
    void Initialize_SingleFluid();
    void Run_SingleFluid();

    void WriteDebug_SingleFluid();
    // test utilities
    void Create_DummyPhase_MGTest();
    void MGTest();

    bool Restart, pBC;
    int timestep, timestepMax;
    int BoundaryCondition;
    double tauA, tauB, rhoA, rhoB;
    double tau, rho0; //only for single-fluid Lee model
    double tauM;      //relaxation time for phase field (or mass)
    double W, gamma, kappa, beta;
    double Fx, Fy, Fz, flux;
    double din, dout, inletA, inletB, outletA, outletB;

    int Nx, Ny, Nz, N, Np;
    int Nxh, Nyh, Nzh, Nh; // extra halo width
    int rank, nprocx, nprocy, nprocz, nprocs;
    double Lx, Ly, Lz;

    std::shared_ptr<Domain> Dm;   // this domain is for analysis
    std::shared_ptr<Domain> Mask; // this domain is for lbm
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm_Regular;
    std::shared_ptr<ScaLBLWideHalo_Communicator> ScaLBL_Comm_WideHalo;

    // input database
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> domain_db;
    std::shared_ptr<Database> freelee_db;
    std::shared_ptr<Database> analysis_db;
    std::shared_ptr<Database> vis_db;

    IntArray Map;
    signed char *id;
    int *NeighborList;
    int *dvcMap;
    double *gqbar, *hq;
    double *mu_phi, *Den, *Phi;
    double *ColorGrad;
    double *Velocity;
    double *Pressure;

    void getPhase(DoubleArray &PhaseValues);
    void getPotential(DoubleArray &PressureValues, DoubleArray &MuValues);
    void getVelocity(DoubleArray &Vx, DoubleArray &Vy, DoubleArray &Vz);
    void getData_RegularLayout(const double *data, DoubleArray &regdata);

    DoubleArray SignDist;

private:
    Utilities::MPI comm;

    size_t dist_mem_size;
    size_t neighborSize;
    // filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];

    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);
    void AssignComponentLabels_ChemPotential_ColorGrad();
};
#endif
