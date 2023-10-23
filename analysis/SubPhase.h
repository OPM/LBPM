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

#ifndef SubPhase_INC
#define SubPhase_INC

#include <vector>
#include "common/Domain.h"
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "analysis/distance.h"
#include "analysis/Minkowski.h"
#include "common/Utilities.h"
#include "common/MPI.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"

class phase {
public:
    int Nc;
    double p;
    double M, Px, Py, Pz, K, visc;
    double V, A, H, X;
    void reset() {
        p = M = Px = Py = Pz = K = 0.0;
        visc = 0.0;
        V = A = H = X = 0.0;
        Nc = 1;
    }

private:
};

class interface {
public:
    int Nc;
    double M, Px, Py, Pz, K;
    double Mw, Mn, Pnx, Pny, Pnz, Pwx, Pwy, Pwz, Kw, Kn;
    double V, A, H, X;
    void reset() {
        Nc = 0;
        M = Px = Py = Pz = K = 0.0;
        V = A = H = X = 0.0;
        Mw = Mn = Pnx = Pny = Pnz = Pwx = Pwy = Pwz = Kw = Kn = 0.0;
    }

private:
};

class SubPhase {
public:
    std::shared_ptr<Domain> Dm;
    double Volume;
    // input variables
    double rho_n, rho_w;
    double nu_n, nu_w;
    double gamma_wn, beta;
    double Fx, Fy, Fz;
    /*
	 * indices 
	 * 	    w - water phase
	 * 	    n - not water phase
	 * 		c - connected part
	 * 		d - disconnected part
	 * 		i - interface region
	 * 		b - bulk (total)
	 */
    // local entities
    phase wc, wd, wb, nc, nd, nb, solid;
    interface iwn, iwnc;
    interface ifs;

    // global entities
    phase gwc, gwd, gwb, gnc, gnd, gnb, gsolid;
    interface giwn, giwnc;
    interface gifs;
    /* fluid-solid wetting interaction */
    double total_wetting_interaction, count_wetting_interaction;
    double total_wetting_interaction_global, count_wetting_interaction_global;

    //...........................................................................
    int Nx, Ny, Nz;
    IntArray PhaseID;      // Phase ID array (solid=0, non-wetting=1, wetting=2)
    BlobIDArray Label_WP;  // Wetting phase label
    BlobIDArray Label_NWP; // Non-wetting phase label index (0:nblobs-1)
    std::vector<BlobIDType>
        Label_NWP_map;    // Non-wetting phase label for each index
    DoubleArray Rho_n;    // density field
    DoubleArray Rho_w;    // density field
    DoubleArray Phi;      // phase indicator field
    DoubleArray DelPhi;   // Magnitude of Gradient of the phase indicator field
    DoubleArray Pressure; // pressure field
    DoubleArray Vel_x;    // velocity field
    DoubleArray Vel_y;
    DoubleArray Vel_z;
    DoubleArray Dissipation;
    DoubleArray SDs;

    std::shared_ptr<Minkowski> morph_w;
    std::shared_ptr<Minkowski> morph_n;
    std::shared_ptr<Minkowski> morph_i;

    SubPhase(std::shared_ptr<Domain> Dm);
    ~SubPhase();

    void SetParams(double rhoA, double rhoB, double tauA, double tauB,
                   double force_x, double force_y, double force_z, double alpha,
                   double beta);
    void Basic();
    void Full();
    void Write(int time);
    void AggregateLabels(const std::string &filename);

private:
    FILE *TIMELOG;
    FILE *SUBPHASE;
};

#endif
