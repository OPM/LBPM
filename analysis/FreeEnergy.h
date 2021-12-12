/*
 *  averaging tools for electrochemistry
 */

#ifndef FreeEnergyAnalyzer_INC
#define FreeEnergyAnalyzer_INC

#include <vector>
#include "common/Domain.h"
#include "common/Utilities.h"
#include "common/MPI.h"
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "analysis/distance.h"
#include "analysis/Minkowski.h"
#include "analysis/SubPhase.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"
#include "models/FreeLeeModel.h"

/**
 * \class FreeEnergyAnalyzer
 *
 * @brief 
 * The FreeEnergyAnalyzer class is constructed to analyze the LBPM free energy model for liquid-gas systems
 * 
 */

class FreeEnergyAnalyzer {
public:
    std::shared_ptr<Domain> Dm;
    double Volume;
    // input variables
    double rho_n, rho_w;
    double nu_n, nu_w;
    double gamma_wn, beta;
    double Fx, Fy, Fz;

    //...........................................................................
    int Nx, Ny, Nz;
    DoubleArray Rho;
    DoubleArray Phi;
    DoubleArray ChemicalPotential;
    DoubleArray Pressure;
    DoubleArray Vel_x;
    DoubleArray Vel_y;
    DoubleArray Vel_z;
    DoubleArray SDs;

    FreeEnergyAnalyzer(std::shared_ptr<Domain> Dm);
    ~FreeEnergyAnalyzer();

    void SetParams();
    void Basic(ScaLBL_FreeLeeModel &LeeModel, int timestep);
    void WriteVis(ScaLBL_FreeLeeModel &LeeModel,
                  std::shared_ptr<Database> input_db, int timestep);

private:
    FILE *TIMELOG;
};
#endif
