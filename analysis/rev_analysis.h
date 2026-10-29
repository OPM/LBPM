
#ifndef REV_ANALYSIS_H
#define REV_ANALYSIS_H
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
 * Representative elementary volume analysis
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <array>
#include <limits>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/MPI.h"
#include "analysis/Minkowski.h"
#include "ProfilerApp.h"
#include "common/Domain.h"



class ScaLBL_REV {
public:
    ScaLBL_REV(int RANK, int NP, const Utilities::MPI &COMM);
    ~ScaLBL_REV();

    // Standalone analysis: read settings, load the domain, then run analysis.
    // The standalone executable supplies false; explicit REV settings override it.
    void ReadParams(string filename, bool hydraulic_default = true);
    void SetDomain();
    void PoreSize();
    void DetRevAnalysis();
    void StatRevAnalysis();
    void sub_domain_select();

    // Captured at the first statistical size with CV < REV.gamma.
    // Keep these separate from the deterministic rev_x_* estimates.
    struct StatisticalREVResult {
        int size = 0;
        std::array<int, 3> dimensions{{0, 0, 0}};
        double mean = 0.0;
        double cv = -1.0;
        // Specific surface measured in the same samples at the porosity REV.
        double surface_mean = std::numeric_limits<double>::quiet_NaN();
        // Displayed hydraulic tortuosity (not tortuosity - 1), at the permeability REV.
        double tortuosity_mean = std::numeric_limits<double>::quiet_NaN();
    };
    StatisticalREVResult stat_poro_rev, stat_perm_rev;

    // Models supply dimensions, domain, distance, velocity, tau and force.
    template<class Model> void PoreSize(const Model& model) {
        PoreSizeImpl(GetModelData(model));
    }
    template<class Model> void DetRevAnalysis(const Model& model) {
        if (run_rev && (geometrical || hydraulic))
            DetRevAnalysisImpl(GetModelData(model));
    }
    template<class Model> void StatRevAnalysis(const Model& model) {
        if (run_rev && (geometrical || hydraulic))
            StatRevAnalysisImpl(GetModelData(model));
    }
    template<class Model> void sub_domain_select(const Model& model) {
        if (run_rev && select_subdomains && (geometrical || hydraulic))
            SubDomainSelectImpl(GetModelData(model));
    }

    double average_pore_size;
    double rev_x_poro, rev_x_perm, rev_x_tort, rev_x_surf;
    double size_step;
    int cc_interval, stat_samples;
    double gamma;
    bool geometrical = true;
    bool hydraulic = true;
    bool seed_anchor = true;

    DoubleArray Distance, Velocity_x, Velocity_y, Velocity_z;
    double tau, Fx, Fy, Fz;

    int Nx, Ny, Nz, N, Np;
    int rank, nprocx, nprocy, nprocz, nprocs;
    double Lx, Ly, Lz;

    std::shared_ptr<Domain> Dm;   // this domain is for analysis
    std::shared_ptr<Domain> Mask; // this domain is for lbm
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
    // input database
    std::shared_ptr<Database> db;

private:
    // Borrow model arrays so generic calls do not copy simulation fields.
    struct ModelData {
        int Nx, Ny, Nz;
        double voxel_length;
        const DoubleArray& Distance;
        const DoubleArray& Velocity_x;
        const DoubleArray& Velocity_y;
        const DoubleArray& Velocity_z;
        double tau, Fx, Fy, Fz;
    };

    template<class Model> static ModelData GetModelData(const Model& model) {
        INSIST(model.Dm, "REV analysis requires an initialized domain");
        return {model.Nx, model.Ny, model.Nz, model.Dm->voxel_length,
                model.Distance, model.Velocity_x, model.Velocity_y,
                model.Velocity_z, model.tau, model.Fx, model.Fy, model.Fz};
    }

    void PoreSizeImpl(const ModelData& model);
    void DetRevAnalysisImpl(const ModelData& model);
    void StatRevAnalysisImpl(const ModelData& model);
    void LoadVelocityFields();
    void SubDomainSelectImpl(const ModelData& model);
    std::shared_ptr<Database> selection_flow_params;
    Utilities::MPI comm;
    bool run_rev = false;
    bool standalone_ready = false;
    bool select_subdomains = false;
    bool save_subdomain_raw = false;
    // Zero means use each property's retained statistical REV dimensions.
    std::array<int, 3> subdomain_size{{0, 0, 0}};
};

#endif // REV_ANALYSIS_H
