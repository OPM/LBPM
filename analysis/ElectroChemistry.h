/*
 *  averaging tools for electrochemistry
 */

#ifndef ElectroChem_INC
#define ElectroChem_INC

#include <vector>
#include "common/Domain.h"
#include "common/Utilities.h"
#include "common/MPI.h"
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "analysis/distance.h"
#include "analysis/Minkowski.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"
#include "models/IonModel.h"
#include "models/PoissonSolver.h"
#include "models/StokesModel.h"

class ElectroChemistryAnalyzer {
public:
    std::shared_ptr<Domain> Dm;
    double Volume;
    // input variables
    double rho_n, rho_w;
    double nu_n, nu_w;
    double gamma_wn, beta;
    double Fx, Fy, Fz;
    
    bool USE_MEMBRANE;

    //...........................................................................
    int Nx, Ny, Nz;
    DoubleArray Rho;                 // density field
    DoubleArray ChemicalPotential;   // density field
    DoubleArray ElectricalPotential; // density field
    DoubleArray ElectricalField_x;   // density field
    DoubleArray ElectricalField_y;   // density field
    DoubleArray ElectricalField_z;   // density field
    DoubleArray Pressure;            // pressure field
    DoubleArray Vel_x;               // velocity field
    DoubleArray Vel_y;
    DoubleArray Vel_z;
    DoubleArray SDs;
    DoubleArray IonFluxDiffusive_x; //ion diffusive flux components
    DoubleArray IonFluxDiffusive_y;
    DoubleArray IonFluxDiffusive_z;
    DoubleArray IonFluxAdvective_x; //ion advective flux components
    DoubleArray IonFluxAdvective_y;
    DoubleArray IonFluxAdvective_z;
    DoubleArray IonFluxElectrical_x; //ion electromigration flux components
    DoubleArray IonFluxElectrical_y;
    DoubleArray IonFluxElectrical_z;

    ElectroChemistryAnalyzer(std::shared_ptr<Domain> Dm);
    ElectroChemistryAnalyzer( ScaLBL_IonModel &IonModel);
    ~ElectroChemistryAnalyzer();

    void SetParams();
    void Basic(ScaLBL_IonModel &Ion, ScaLBL_Poisson &Poisson, ScaLBL_StokesModel &Stokes, int timestep);
    void Membrane(ScaLBL_IonModel &Ion, ScaLBL_Poisson &Poisson, int timestep);
    void WriteVis(ScaLBL_IonModel &Ion, ScaLBL_Poisson &Poisson, 
    		ScaLBL_StokesModel &Stokes,std::shared_ptr<Database> input_db, int timestep);
    void Basic(ScaLBL_IonModel &Ion, ScaLBL_Poisson &Poisson, int timestep);
    void WriteVis(ScaLBL_IonModel &Ion, ScaLBL_Poisson &Poisson,
                  std::shared_ptr<Database> input_db, int timestep);

private:
    FILE *TIMELOG;
};
#endif
