/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
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
#include "analysis/FlowAdaptor.h"
#include "analysis/TwoPhase.h"
#include "analysis/runAnalysis.h"
#include "common/MPI.h"
#include "ProfilerApp.h"
#include "threadpool/thread_pool.h"

#ifndef ScaLBL_ColorModel_INC
#define ScaLBL_ColorModel_INC

/**
 * \class ScaLBL_ColorModel
 *
 * @details
 * The ScaLBL_ColorModel class contains routines to initialize and run a two-component color lattice Boltzmann model
 * Momentum transport equations are described by a D3Q19 scheme
 * Mass transport equations are described by D3Q7 scheme
 */

class ScaLBL_ColorModel {
public:
    /**
    * \brief Constructor
    * @param RANK        processor rank 
    * @param NP        number of processors 
    * @param COMM        MPI communicator 
    */
    ScaLBL_ColorModel(int RANK, int NP, const Utilities::MPI &COMM);
    ~ScaLBL_ColorModel();

    /**
    * \brief Read simulation parameters
    * @param filename       input database file that includes "Color" section 
    */
    void ReadParams(string filename);

    /**
    * \brief Read simulation parameters
    * @param db0       input database that includes "Color" section 
    */
    void ReadParams(std::shared_ptr<Database> db0);

    /**
    * \brief Create domain data structures
    */
    void SetDomain();

    /**
    * \brief Read image data
    */
    void ReadInput();

    /**
    * \brief Create color model data structures
    */
    void Create();

    /**
    * \brief Initialize the simulation
    */
    void Initialize();

    /**
    * \brief Run the simulation
    */
    void Run();

    /**
    * \brief Run the simulation
    * @param returntime -  timestep at which the routine will return
    */
    double Run(int returntime);

    /**
    * \brief Debugging function to dump simulation state to disk
    */
    void WriteDebug();

    /**
    * \brief Copy the phase field for use by external methods 
    * @param f   - DoubleArray to hold the phase field
    */
    void getPhaseField(DoubleArray &f);

    bool Restart, pBC;
    bool REVERSE_FLOW_DIRECTION;
    int timestep, timestepMax;
    int BoundaryCondition;
    double tauA, tauB, rhoA, rhoB, alpha, beta;
    double Fx, Fy, Fz, flux;
    double din, dout, inletA, inletB, outletA, outletB;
    const double mDarcy_converter = 1013.0;

    int Nx, Ny, Nz, N, Np;
    int rank, nprocx, nprocy, nprocz, nprocs;
    double Lx, Ly, Lz;

    std::shared_ptr<Domain> Dm;   // this domain is for analysis
    std::shared_ptr<Domain> Mask; // this domain is for lbm
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm_Regular;
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

    /**
    * \brief Assign wetting affinity values 
    */
    void AssignComponentLabels(double *phase);

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
};

#endif
