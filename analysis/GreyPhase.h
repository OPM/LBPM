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
 * Sub-phase averaging tools
 */

#ifndef GreyPhase_INC
#define GreyPhase_INC

#include <vector>
#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "common/Utilities.h"
#include "common/MPI.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"

/**
 * \class GreyPhase
 *
 * @brief 
 * The GreyPhase class tracks pressure, mass and momentum within a grey phase 
 * 
 */
class GreyPhase {
public:
    double p;
    double M, Px, Py, Pz;
    void reset() { p = M = Px = Py = Pz = 0.0; }

private:
};

/**
 * \class GreyPhaseAnalysis
 *
 * @brief 
 * The GreyPhaseAnalysis class is constructed to analyze the LBPM greyscale model
 * 
 */
class GreyPhaseAnalysis {
public:
    std::shared_ptr<Domain> Dm;
    double Volume;
    // input variables
    double rho_n, rho_w;
    double nu_n, nu_w;
    double gamma_wn, beta;
    double Fx, Fy, Fz;
    double grey_porosity;
    // outputs
    double saturation, water_flow_rate, oil_flow_rate;

    //simulation outputs (averaged values)
    GreyPhase Water, Oil;
    GreyPhase Water_local, Oil_local;
    //...........................................................................
    int Nx, Ny, Nz;
    //IntArray PhaseID;		// Phase ID array
    DoubleArray SDs;      // contains porosity map
    DoubleArray Porosity; // contains porosity map
    DoubleArray Rho_n;    // density field
    DoubleArray Rho_w;    // density field
    //DoubleArray Phi;		// phase indicator field
    //DoubleArray DelPhi;		// Magnitude of Gradient of the phase indicator field
    DoubleArray Pressure; // pressure field
    DoubleArray Vel_x;    // velocity field
    DoubleArray Vel_y;
    DoubleArray Vel_z;
    DoubleArray MobilityRatio;

    GreyPhaseAnalysis(std::shared_ptr<Domain> Dm);
    ~GreyPhaseAnalysis();

    void SetParams(double rhoA, double rhoB, double tauA, double tauB,
                   double force_x, double force_y, double force_z, double alpha,
                   double beta, double GreyPorosity);
    void Basic();
    void Write(int time);

private:
    FILE *TIMELOG;
};

#endif
