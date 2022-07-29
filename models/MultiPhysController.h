/*
 * Multiphysics controller that coordinates the coupling between different models
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <vector>
#include <algorithm>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/MPI.h"
#include "analysis/Minkowski.h"
#include "ProfilerApp.h"

class ScaLBL_Multiphys_Controller {
public:
    ScaLBL_Multiphys_Controller(int RANK, int NP, const Utilities::MPI &COMM);
    ~ScaLBL_Multiphys_Controller();

    void ReadParams(string filename);
    void ReadParams(std::shared_ptr<Database> db0);
    int getStokesNumIter_PNP_coupling(double StokesTimeConv,
                                      const vector<double> &IonTimeConv);
    vector<int> getIonNumIter_PNP_coupling(double StokesTimeConv,
                                           const vector<double> &IonTimeConv);
    vector<int> getIonNumIter_NernstPlanck_coupling(const vector<double> &IonTimeConv);
    void getTimeConvMax_PNP_coupling(double StokesTimeConv,
                                     const vector<double> &IonTimeConv);

    bool Restart;
    int timestepMax;
    int num_iter_Stokes;
    vector<int> num_iter_Ion;
    int analysis_interval;
    int visualization_interval;
    int restart_interval;
    double tolerance;
    double time_conv_max;
    double time_conv_MainLoop;
    //double SchmidtNum;//Schmidt number = kinematic_viscosity/mass_diffusivity

    int rank, nprocs;

    // input database
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> study_db;

private:
    Utilities::MPI comm;

    // filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];

    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);
};
