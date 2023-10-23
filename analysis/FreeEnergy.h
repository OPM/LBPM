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
