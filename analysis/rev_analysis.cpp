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
#include "analysis/rev_analysis.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"
#include "IO/Reader.h"
#include <cmath>
#include <algorithm>
#include <random>

ScaLBL_REV::ScaLBL_REV(int RANK, int NP, const Utilities::MPI &COMM)
    : average_pore_size(0), rev_x_poro(-1), rev_x_perm(-1), rev_x_tort(-1),
      rev_x_surf(-1), size_step(0.06), cc_interval(10), stat_samples(50), gamma(0.1),
      tau(1.0), Fx(0), Fy(0), Fz(1.0e-5),
      Nx(0), Ny(0), Nz(0), N(0), Np(0), rank(RANK), nprocx(0), nprocy(0),
      nprocz(0), nprocs(NP), Lx(0), Ly(0), Lz(0), comm(COMM) {}

ScaLBL_REV::~ScaLBL_REV() = default;

void ScaLBL_REV::ReadParams(string filename, bool hydraulic_default) {

    db = std::make_shared<Database>(filename);
    run_rev = false;
    select_subdomains = false;
    subdomain_size = {{0, 0, 0}};
    stat_poro_rev = StatisticalREVResult{};
    stat_perm_rev = StatisticalREVResult{};
    if(db->keyExists("REV")){
        run_rev = true;
    }
    else{
        return;
    } 
    auto mrt_db = db->keyExists("MRT") ? db->getDatabase("MRT") : nullptr;

    size_step = 0.1;
    cc_interval = 10;
    gamma = 0.1;
    stat_samples = 50;
    geometrical = true;
    hydraulic = hydraulic_default;
    standalone_ready = false;
    tau = 1.0;
    Fx = Fy = 0.0;
    Fz = 1.0e-5;

    auto rev_db = db->getDatabase("REV");
    geometrical = rev_db->getWithDefault<bool>("Geometrical", true);
    hydraulic = rev_db->keyExists("Hydraulic")
        ? rev_db->getScalar<bool>("Hydraulic")
        : hydraulic_default;
    seed_anchor = rev_db->getWithDefault<bool>("Seed_Anchor", true);
    select_subdomains = rev_db->getWithDefault<bool>("SubdomainSelect", false);
    save_subdomain_raw = rev_db->getWithDefault<bool>("save_subdomain_raw", false);
    if (rev_db->keyExists("subdomain_size")) {
        const std::string message = "REV.subdomain_size requires one or three positive integer voxel dimensions";
        if (!rev_db->isType<double>("subdomain_size"))
            throw std::runtime_error(message);
        const auto size = rev_db->getVector<double>("subdomain_size");
        if (size.size() != 1 && size.size() != 3)
            throw std::runtime_error(message);
        for (int d = 0; d < 3; ++d) {
            const double value = size[size.size() == 1 ? 0 : d];
            if (!std::isfinite(value) || value < 1.0 ||
                value > std::numeric_limits<int>::max() || std::floor(value) != value)
                throw std::runtime_error(message);
            subdomain_size[d] = static_cast<int>(value);
        }
    }
    INSIST(!select_subdomains || comm.getSize() == 1,
           "REV.SubdomainSelect requires one MPI process");
    if (rev_db->keyExists("SizeStep"))
        size_step = rev_db->getScalar<double>("SizeStep");
    if (rev_db->keyExists("CCInterval"))
        cc_interval = rev_db->getScalar<int>("CCInterval");
    if (rev_db->keyExists("gamma"))
        gamma = rev_db->getScalar<double>("gamma");
    if (rev_db->keyExists("StatSamples"))
        stat_samples = rev_db->getScalar<int>("StatSamples");

    INSIST(size_step > 0.0, "REV.SizeStep must be positive");
    INSIST(cc_interval > 0, "REV.CCInterval must be positive");
    INSIST(stat_samples > 0, "REV.StatSamples must be positive");
    INSIST(std::isfinite(gamma) && gamma > 0.0, "REV.gamma must be positive and finite");
    auto simulation_params = mrt_db;
    if (rev_db->keyExists("SimulationParams")) {
        auto simulation_filename =
            rev_db->getScalar<std::string>("SimulationParams");
        INSIST(!simulation_filename.empty(), "REV.SimulationParams cannot be empty");
        const bool absolute = simulation_filename.front() == '/' ||
                              simulation_filename.front() == '\\' ||
                              (simulation_filename.size() > 1 &&
                               simulation_filename[1] == ':');
        if (!absolute) {
            auto separator = filename.find_last_of("/\\");
            if (separator != std::string::npos)
                simulation_filename = filename.substr(0, separator + 1) +
                                      simulation_filename;
        }

        auto simulation_db = std::make_shared<Database>(simulation_filename);
        INSIST(simulation_db->keyExists("MRT"),
               "REV.SimulationParams must reference a database with an MRT section");
        simulation_params = simulation_db->getDatabase("MRT");
    }
    if (simulation_params) {
        tau = simulation_params->getWithDefault<double>("tau", 1.0);
        if (simulation_params->keyExists("F")) {
            auto force = simulation_params->getVector<double>("F");
            INSIST(force.size() == 3, "MRT.F must have three components");
            Fx = force[0];
            Fy = force[1];
            Fz = force[2];
        }
    }
    selection_flow_params = simulation_params;
}

void ScaLBL_REV::SetDomain() {
    standalone_ready = false;
    stat_poro_rev = StatisticalREVResult{};
    stat_perm_rev = StatisticalREVResult{};
    if (!run_rev || (!geometrical && !hydraulic))
        return;
    INSIST(db, "Call REV.ReadParams before REV.SetDomain");
    INSIST(comm.getSize() == 1,
           "Standalone REV analysis currently requires one MPI process");

    auto domain_db = db->getDatabase("Domain");
    Dm = std::make_shared<Domain>(domain_db, comm);
    Mask = Dm;
    Nx = Dm->Nx;
    Ny = Dm->Ny;
    Nz = Dm->Nz;
    N = Nx * Ny * Nz;
    Lx = Dm->Lx;
    Ly = Dm->Ly;
    Lz = Dm->Lz;
    rank = Dm->rank();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();

    // Use the full domain when building the signed distance field.
    std::fill(Dm->id.begin(), Dm->id.end(), 1);
    Dm->CommInit();
    if (domain_db->keyExists("Filename")) {
        Dm->Decomp(domain_db->getScalar<std::string>("Filename"));
    } else if (domain_db->keyExists("GridFile")) {
        auto input_id = readMicroCT(*domain_db, comm);
        INSIST(input_id.size(0) == size_t(Nx - 2) &&
                   input_id.size(1) == size_t(Ny - 2) &&
                   input_id.size(2) == size_t(Nz - 2),
               "Input geometry does not match Domain.n");
        Array<signed char> id_view;
        id_view.viewRaw({size_t(Nx), size_t(Ny), size_t(Nz)}, Dm->id.data());
        fillHalo<signed char> fill(comm, Dm->rank_info,
                                  {Nx - 2, Ny - 2, Nz - 2}, {1, 1, 1}, 0, 1);
        fill.copy(input_id, id_view);
    } else {
        Dm->ReadIDs();
    }
    Np = Dm->PoreCount();

    Distance.resize(Nx, Ny, Nz);
    Array<char> pore(Nx, Ny, Nz);
    for (int i = 0; i < N; ++i)
        pore(i) = Dm->id[i] > 0 ? 1 : 0;
    CalcDist(Distance, pore, *Dm);

    Velocity_x.resize(Nx, Ny, Nz);
    Velocity_y.resize(Nx, Ny, Nz);
    Velocity_z.resize(Nx, Ny, Nz);
    Velocity_x.fill(0.0);
    Velocity_y.fill(0.0);
    Velocity_z.fill(0.0);
    if (hydraulic)
        LoadVelocityFields();

    average_pore_size = 0.0;
    rev_x_poro = rev_x_perm = rev_x_tort = rev_x_surf = -1.0;
    standalone_ready = true;
}

void ScaLBL_REV::LoadVelocityFields() {
    auto rev_db = db->getDatabase("REV");
    auto path = rev_db->getWithDefault<std::string>("VelocityPath", ".");
    std::string timestep;
    if (rev_db->keyExists("VelocityTimestep")) {
        timestep = rev_db->getScalar<std::string>("VelocityTimestep");
    } else {
        INSIST(exists(path + "/LBM.visit") || exists(path + "/summary.LBM"),
               "Standalone hydraulic REV requires saved velocity output: "
               "set REV.VelocityPath and REV.VelocityTimestep, "
               "or set REV.Hydraulic = false");
        auto timesteps = IO::readTimesteps(path);
        INSIST(!timesteps.empty(), "No saved velocity timesteps found");
        timestep = timesteps.back();
    }
    INSIST(exists(path + "/" + timestep + "/LBM.summary"),
           "Cannot find LBM.summary for REV velocity input");
    auto meshes = IO::getMeshList(path, timestep);
    auto entry = std::find_if(meshes.begin(), meshes.end(),
                              [](const IO::MeshDatabase& mesh) {
                                  return mesh.name == "domain";
                              });
    INSIST(entry != meshes.end(), "Saved velocity output has no domain mesh");
    INSIST(entry->domains.size() == 1,
           "Standalone REV requires saved output from one MPI domain");
    auto mesh = IO::getMesh(path, timestep, *entry, 0);
    auto domain_mesh = std::dynamic_pointer_cast<IO::DomainMesh>(mesh);
    INSIST(domain_mesh && domain_mesh->nx == Nx - 2 &&
               domain_mesh->ny == Ny - 2 && domain_mesh->nz == Nz - 2,
           "Saved velocity dimensions do not match the REV geometry");

    fillHalo<double> fill(comm, Dm->rank_info,
                          {Nx - 2, Ny - 2, Nz - 2}, {1, 1, 1}, 0, 1);
    const char* names[] = {"Velocity_x", "Velocity_y", "Velocity_z"};
    DoubleArray* fields[] = {&Velocity_x, &Velocity_y, &Velocity_z};
    for (int i = 0; i < 3; ++i) {
        auto variable = IO::getVariable(path, timestep, *entry, 0, names[i]);
        INSIST(variable && variable->dim == 1 &&
                   variable->type == IO::VariableType::VolumeVariable,
               std::string("Missing scalar velocity field: ") + names[i]);
        IO::reformatVariable(*mesh, *variable);
        INSIST(variable->data.size(0) == size_t(Nx - 2) &&
                   variable->data.size(1) == size_t(Ny - 2) &&
                   variable->data.size(2) == size_t(Nz - 2),
               "Saved velocity array dimensions do not match the REV geometry");
        fill.copy(variable->data, *fields[i]);
    }
}

void ScaLBL_REV::PoreSize() {
    INSIST(standalone_ready, "Call REV.SetDomain before REV.PoreSize");
    PoreSize(*this);
}

void ScaLBL_REV::DetRevAnalysis() {
    if (!run_rev || (!geometrical && !hydraulic))
        return;
    INSIST(standalone_ready, "Call REV.SetDomain before REV.DetRevAnalysis");
    DetRevAnalysis(*this);
}

void ScaLBL_REV::StatRevAnalysis() {
    if (!run_rev || (!geometrical && !hydraulic))
        return;
    INSIST(standalone_ready, "Call REV.SetDomain before REV.StatRevAnalysis");
    if (average_pore_size == 0.0)
        PoreSize();
    StatRevAnalysis(*this);
}

void ScaLBL_REV::sub_domain_select() {
    if (!run_rev || !select_subdomains || (!geometrical && !hydraulic))
        return;
    INSIST(standalone_ready, "Call REV.SetDomain before REV.sub_domain_select");
    sub_domain_select(*this);
}

struct det_convergence_metrics {
    double RE;
    double CC;
};

struct stat_convergence_metrics {
    double CV, mean, max, min;
};

struct VoxelDist {
    int x, y, z;
    double dist;
};

det_convergence_metrics det_compute_convergence(const std::vector<double>& data,
                                               const std::vector<double>& volumes,
                                               int iter_step, int cc_interval) {
    det_convergence_metrics m = {-1, -1};
    
    if (iter_step < 1) 
        return m;

    double mean = 0.0;
    double std_dev = 0.0;
    int valid_count = 0;
    int start_idx = (iter_step >= cc_interval - 1) ? (iter_step + 1 - cc_interval) : 0;
    for (int l = start_idx; l <= iter_step; l++) {
        if (data[l] != -1) {
            mean += data[l];
            valid_count++;
        }
    }

    if (valid_count == 0) return m;
    
    mean /= valid_count;

    for (int l = start_idx; l <= iter_step; l++) {
        if (data[l] != -1) {
            std_dev += (data[l] - mean) * (data[l] - mean);
        }
    }
    std_dev = sqrt(std_dev / valid_count);

    if (mean != 0) {
        m.CC = std_dev / mean;
    }

    if (data[iter_step] != -1 && data[iter_step - 1] != -1) {
        double sum_prev_curr = data[iter_step] + data[iter_step - 1];
        if (sum_prev_curr != 0 && volumes[iter_step] > 0) {
            m.RE = fabs((data[iter_step] - data[iter_step - 1]) / sum_prev_curr *
                        (volumes[iter_step - 1] / volumes[iter_step]));
        }
    }

    return m;
}

stat_convergence_metrics stat_compute_convergence(const std::vector<double>& data, int samples) {
    stat_convergence_metrics res = {-1.0, 0, -1, 1e100};
    int valid_count = 0;

    for (int l = 0; l < samples; l++) {
        if (data[l] != -1) {
            res.mean += data[l];
            if (data[l] > res.max) res.max = data[l];
            if (data[l] < res.min) res.min = data[l];
            valid_count++;
        }
    }

    if (valid_count == 0) return res;

    res.mean /= valid_count;

    double std_dev = 0.0;
    for (int l = 0; l < samples; l++) {
        if (data[l] != -1) {
            std_dev += (data[l] - res.mean) * (data[l] - res.mean);
        }
    }
    std_dev = std::sqrt(std_dev / valid_count);
    
    if (res.mean != 0) res.CV = std_dev / std::abs(res.mean);
    else if (std_dev == 0.0) res.CV = 0.0;

    return res;
}

void ScaLBL_REV::PoreSizeImpl(const ModelData &model) {
    int Nx = model.Nx, Ny = model.Ny, Nz = model.Nz;
    double h = model.voxel_length;

    DoubleArray distance_updated;
    distance_updated.resize(Nx, Ny, Nz);
    
    std::vector<VoxelDist> pore_voxels;

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                double val = model.Distance(i, j, k);
                distance_updated(i, j, k) = val;
                
                if (i > 0 && i < Nx - 1 && j > 0 && j < Ny - 1 && k > 0 && k < Nz - 1) {
                    if (val > 0) {
                        pore_voxels.push_back({i, j, k, val});
                    }
                }
            }
        }
    }

    std::sort(pore_voxels.begin(), pore_voxels.end(), [](const VoxelDist& a, const VoxelDist& b) {
        return a.dist > b.dist;
    });

    for (const auto& voxel : pore_voxels) {
        int i = voxel.x;
        int j = voxel.y;
        int k = voxel.z;
        double current_max = voxel.dist;

        int start_x = std::max(1, static_cast<int>(i - current_max - 1));
        int end_x   = std::min(Nx - 1, static_cast<int>(i + current_max + 2));
        int start_y = std::max(1, static_cast<int>(j - current_max - 1));
        int end_y   = std::min(Ny - 1, static_cast<int>(j + current_max + 2));
        int start_z = std::max(1, static_cast<int>(k - current_max - 1));
        int end_z   = std::min(Nz - 1, static_cast<int>(k + current_max + 2));

        double r_sq = current_max * current_max;

        for (int n = start_z; n < end_z; n++) {
            for (int m = start_y; m < end_y; m++) {
                for (int l = start_x; l < end_x; l++) {
                    double dist_sq = (i - l) * (i - l) + (j - m) * (j - m) + (k - n) * (k - n);
                    
                    if (dist_sq <= r_sq) {
                        double value2 = distance_updated(l, m, n);
                        if (value2 > 0 && value2 < current_max) {
                            distance_updated(l, m, n) = current_max;
                        }
                    }
                }
            }
        }
    }

    double sum = 0.0;
    int total_pore_voxels = 0;
    
    std::map<double, int> psd_histogram; 

    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                double val = distance_updated(i, j, k);
                if (val > 0) {
                    sum += val;
                    total_pore_voxels++;
                    psd_histogram[val]++;
                }
            }
        }
    }

    average_pore_size = (total_pore_voxels > 0) ? (2.0 * sum / total_pore_voxels) : 0.0;

    std::ofstream outfile("pore_size_distribution.csv");
    if (outfile.is_open()) {
        outfile << "radius_voxels diameter_um voxel_count volume_fraction\n";
        
        for (const auto& pair : psd_histogram) {
            double radius_vx = pair.first;
            int count = pair.second;
            
            double diameter_um = 2.0 * radius_vx * h;
            double volume_fraction = static_cast<double>(count) / total_pore_voxels;
            
            outfile << radius_vx << " "
                    << diameter_um << " " 
                    << count << " " 
                    << volume_fraction << "\n";
        }
        
        outfile << "\nAverage_Pore_Size_um " << average_pore_size * h << "\n";
        outfile.close();
    }

    std::cout << "\n\n\nAverage Pore Size: " << average_pore_size * h << " micrometers\n\n\n" << std::endl;
}

void ScaLBL_REV::DetRevAnalysisImpl(const ModelData &model) {
    if (!run_rev)
        return;

    PoreSizeImpl(model);

    double current_size_x_double = 2.0 * average_pore_size;
    int x_size = static_cast<int>(std::ceil(current_size_x_double)); //could be a problem if the geometry is in the yz plane.


    int Nx = model.Nx, Ny = model.Ny, Nz = model.Nz;
    int iter_step = 0; 
    double h = model.voxel_length;
    double mu = (model.tau - 0.5) / 3.f;
    
    rev_x_poro = -1;
    rev_x_perm = -1;
    rev_x_tort = -1;
    rev_x_surf = -1;

    int vector_size_x = 1000;
    std::vector<double> sub_volume(vector_size_x);
    std::vector<double> perm(vector_size_x);
    std::vector<double> poro(vector_size_x);
    std::vector<double> tort(vector_size_x);
    std::vector<double> surf(vector_size_x);

    std::ofstream geometric_log;
    if (geometrical) {
        geometric_log.open("det_rev_geometric.csv", std::ios::app);
        if (geometric_log.tellp() == 0) {
            geometric_log << "iter_step x_size_um x_size_vx "
                          << "porosity rev_poro_um rev_poro_vx RE_poro CC_poro "
                          << "surface rev_surf_um rev_surf_vx RE_surf CC_surf\n";
        }
    }

    std::ofstream hydraulic_log;
    if (hydraulic) {
        hydraulic_log.open("det_rev_hydraulic.csv", std::ios::app);
        if (hydraulic_log.tellp() == 0) {
            hydraulic_log << "iter_step x_size_um x_size_vx "
                          << "permeability rev_perm_um rev_perm_vx RE_perm CC_perm "
                          << "tortuosity rev_tort_um rev_tort_vx RE_tort CC_tort\n";
        }
    }


    const int max_x_size = Nx - 2;
    while (true) {
        x_size = std::min(x_size, max_x_size);

        //(still have to test for 2d geometries)
        int y_size = (Ny - 2) * x_size / (Nx - 2);
        int z_size = (Nz - 2) * x_size / (Nx - 2);

        int start_x = (Nx - 2 - x_size) / 2 + 1;
        int start_y = (Ny - 2 - y_size) / 2 + 1;
        int start_z = (Nz - 2 - z_size) / 2 + 1;

        int end_x = (Nx - 2 + x_size) / 2 + 1;
        int end_y = (Ny - 2 + y_size) / 2 + 1;
        int end_z = (Nz - 2 + z_size) / 2 + 1;

        start_x = std::max(1, start_x);
        end_x = std::min(Nx-2, end_x);
        start_y = std::max(1, start_y);
        end_y = std::min(Ny-2, end_y);
        start_z = std::max(1, start_z);
        end_z = std::min(Nz-2, end_z);

        // The common voxel volume cancels in V_(i-1) / V_i.
        sub_volume[iter_step] = static_cast<double>(end_x - start_x + 1) *
                               (end_y - start_y + 1) * (end_z - start_z + 1);

        double current_surf = 0;
        double count = 0.0;
        double vxb = 0.0, vyb = 0.0, vzb = 0.0;
        double vax = 0.0, vay = 0.0, vaz = 0.0;
        double v_tort=0.0;

        for (int k = start_z; k <= end_z; k++) { 
            for (int j = start_y; j <= end_y; j++) {
                for (int i = start_x; i <= end_x; i++) {
                    if (model.Distance(i, j, k) > 0) {
                        if (hydraulic) {
                            vxb = model.Velocity_x(i, j, k);
                            vyb = model.Velocity_y(i, j, k);
                            vzb = model.Velocity_z(i, j, k);
                            vax += vxb;
                            vay += vyb;
                            vaz += vzb;
                            v_tort += sqrt(vxb*vxb + vyb*vyb + vzb*vzb);
                        }
                        count += 1.0;
                    }
                    if (geometrical) {
                        if ((model.Distance(i, j, k) * model.Distance(i, j, k + 1) < 0) && (k != end_z))
                            current_surf++;

                        if ((model.Distance(i, j, k) * model.Distance(i, j + 1, k) < 0) && (j != end_y))
                            current_surf++;

                        if ((model.Distance(i, j, k) * model.Distance(i + 1, j, k) < 0) && (i != end_x))
                            current_surf++;
                    }
                }
            }
        }

        if (hydraulic) {
            if (count == 0)
            {
                vax = 0.0;
                vay = 0.0;
                vaz = 0.0;
            } else {
            vax /= count;
            vay /= count;
            vaz /= count;}

            double force_mag = sqrt(model.Fx * model.Fx + model.Fy * model.Fy + model.Fz * model.Fz);
            double dir_x = model.Fx / force_mag;
            double dir_y = model.Fy / force_mag;
            double dir_z = model.Fz / force_mag;
            if (force_mag == 0.0) {
                dir_x = 0.0;
                dir_y = 0.0;
                dir_z = 1.0;
                force_mag = 1.0;
            }

            double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);
            double current_poro = count / ((end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1));
            double current_perm = 1013 * h * h * mu * current_poro * flow_rate / force_mag;
            double current_tort = (vaz != 0) ? (v_tort / (vaz * count) - 1.0) : -1;

            perm[iter_step] = current_perm;
            tort[iter_step] = current_tort;

            //Permeability deterministic approach
            det_convergence_metrics perm_m = det_compute_convergence(perm, sub_volume, iter_step, cc_interval);
            if (perm_m.RE >= 0 && perm_m.RE < gamma && perm_m.CC >= 0 && perm_m.CC < gamma && rev_x_perm == -1) {
                rev_x_perm = x_size;
            }

            //Tortuosity deterministic approach
            det_convergence_metrics tort_m = det_compute_convergence(tort, sub_volume, iter_step, cc_interval);
            if (tort_m.RE >= 0 && tort_m.RE < gamma && tort_m.CC >= 0 && tort_m.CC < gamma && rev_x_tort == -1) {
                rev_x_tort = x_size;
            }

            hydraulic_log << iter_step << " " << x_size * h << " " << x_size << " "
                          << perm[iter_step] << " " << rev_x_perm*h << " " << rev_x_perm << " " << perm_m.RE << " " << perm_m.CC << " "
                          << tort[iter_step]+1 << " " << rev_x_tort*h << " " << rev_x_tort << " " << tort_m.RE << " " << tort_m.CC << "\n";
        }

        if (geometrical) {
            double vol = (end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1);
            current_surf = (current_surf / vol) / h;
            double current_poro = count / vol;

            poro[iter_step] = current_poro;
            surf[iter_step] = current_surf;

            //Porosity deterministic approach
            det_convergence_metrics poro_m = det_compute_convergence(poro, sub_volume, iter_step, cc_interval);
            if (poro_m.RE >= 0 && poro_m.RE < gamma && poro_m.CC >= 0 && poro_m.CC < gamma && rev_x_poro == -1) {
                rev_x_poro = x_size;
            }

            //Surface Deterministic approach
            det_convergence_metrics surf_m = det_compute_convergence(surf, sub_volume, iter_step, cc_interval);
            if (surf_m.RE >= 0 && surf_m.RE < gamma && surf_m.CC >= 0 && surf_m.CC < gamma && rev_x_surf == -1) {
                rev_x_surf = x_size;
            }

            geometric_log << iter_step << " " << x_size * h << " " << x_size << " "
                          << poro[iter_step] << " " << rev_x_poro*h << " " << rev_x_poro << " " << poro_m.RE << " " << poro_m.CC << " "
                          << surf[iter_step] << " " << rev_x_surf*h << " " << rev_x_surf << " " << surf_m.RE << " " << surf_m.CC << "\n";
        }

        // Stop after processing the full domain once.
        if (x_size == max_x_size) {
            break;
        }

        double min_step = std::max(1.0, average_pore_size);
        double proportional_step = current_size_x_double * size_step;
        current_size_x_double += std::max(min_step, proportional_step);
        x_size = static_cast<int>(std::ceil(current_size_x_double)); 

        iter_step++;
    }

    if (geometrical) {
        geometric_log.close();
    }

    if (hydraulic) {
        hydraulic_log.close();
    }

    printf("\n\nEnded Deterministic Analysis \n\n");
}

void ScaLBL_REV::StatRevAnalysisImpl(const ModelData &model) {
    if (!run_rev)
        return;

    stat_poro_rev = StatisticalREVResult{};
    stat_perm_rev = StatisticalREVResult{};

    std::mt19937 rng(0);
    if (!seed_anchor) {
        std::random_device rd;
        rng.seed(rd());
    }

    int Nx = model.Nx, Ny = model.Ny, Nz = model.Nz;
    
    double h = model.voxel_length;
    double mu = (model.tau - 0.5) / 3.0f;

    std::vector<double> stat_poro(stat_samples);
    std::vector<double> stat_perm(stat_samples);
    std::vector<double> stat_tort(stat_samples);
    std::vector<double> stat_surf(stat_samples);

    std::ofstream geometric_log;
    if (geometrical) {
        geometric_log.open("stat_rev_geometric.csv", std::ios::app);
        if (geometric_log.tellp() == 0) {
            geometric_log << "x_size_um x_size_vx "
                          << "stat_poro stat_poro_mean stat_poro_max stat_poro_min "
                          << "stat_surf stat_surf_mean stat_surf_max stat_surf_min\n";
        }
    }

    std::ofstream hydraulic_log;
    if (hydraulic) {
        hydraulic_log.open("stat_rev_hydraulic.csv", std::ios::app);
        if (hydraulic_log.tellp() == 0) {
            hydraulic_log << "x_size_um x_size_vx "
                          << "stat_perm stat_perm_mean stat_perm_max stat_perm_min "
                          << "stat_tort stat_tort_mean stat_tort_max stat_tort_min\n";
        }
    }

    double current_size_x_double = 1e100;
    if (rev_x_poro > 0) current_size_x_double = std::min(current_size_x_double, (double)rev_x_poro);
    if (rev_x_perm > 0) current_size_x_double = std::min(current_size_x_double, (double)rev_x_perm);
    if (rev_x_tort > 0) current_size_x_double = std::min(current_size_x_double, (double)rev_x_tort);
    if (rev_x_surf > 0) current_size_x_double = std::min(current_size_x_double, (double)rev_x_surf);
    if (current_size_x_double == 1e100) {
        current_size_x_double = 2.0 * average_pore_size;
    }
    if (current_size_x_double == -1) current_size_x_double = 2.0 * average_pore_size;
    int x_size = static_cast<int>(std::ceil(current_size_x_double));

    while (x_size > 0 && x_size < Nx - 2) {
        
        // Match the reported edge length: a cubic image uses L x L x L.
        int sub_y = std::max(1, int(int64_t(x_size) * (Ny - 2) / (Nx - 2)));
        int sub_z = std::max(1, int(int64_t(x_size) * (Nz - 2) / (Nx - 2)));

        int range_x = (Nx - 2) - x_size + 1;
        int range_y = (Ny - 2) - sub_y + 1;
        int range_z = (Nz - 2) - sub_z + 1;

        if (range_x < 1) range_x = 1;
        if (range_y < 1) range_y = 1;
        if (range_z < 1) range_z = 1;

        std::uniform_int_distribution<int> dist_x(1, range_x);
        std::uniform_int_distribution<int> dist_y(1, range_y);
        std::uniform_int_distribution<int> dist_z(1, range_z);

        for (int l = 0; l < stat_samples; l++) {
            int start_x = dist_x(rng);
            int start_y = dist_y(rng);
            int start_z = dist_z(rng);

            int end_x = std::min(Nx - 2, start_x + x_size - 1);
            int end_y = std::min(Ny - 2, start_y + sub_y - 1);
            int end_z = std::min(Nz - 2, start_z + sub_z - 1);

            double vax = 0.0, vay = 0.0, vaz = 0.0, count = 0.0, v_tort = 0.0, current_surf = 0.0;

            for (int k = start_z; k <= end_z; k++) {
                for (int j = start_y; j <= end_y; j++) {
                    for (int i = start_x; i <= end_x; i++) {
                        double dist_val = model.Distance(i, j, k);
                        if (dist_val > 0) {
                            if (hydraulic) {
                                vax += model.Velocity_x(i, j, k);
                                vay += model.Velocity_y(i, j, k);
                                vaz += model.Velocity_z(i, j, k);
                                v_tort += std::sqrt(std::pow(model.Velocity_x(i, j, k), 2) +
                                                    std::pow(model.Velocity_y(i, j, k), 2) +
                                                    std::pow(model.Velocity_z(i, j, k), 2));
                            }
                            count += 1.0;
                        }
                        
                        if (geometrical) {
                            if (k != end_z && (dist_val * model.Distance(i, j, k + 1) < 0)) current_surf++;
                            if (j != end_y && (dist_val * model.Distance(i, j + 1, k) < 0)) current_surf++;
                            if (i != end_x && (dist_val * model.Distance(i + 1, j, k) < 0)) current_surf++;
                        }
                    }
                }
            }

            double vol = (end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1);
            double current_poro = count / vol;

            if (hydraulic) {
                if (count > 0) {
                    vax /= count; vay /= count; vaz /= count;
                }

                double force_mag = std::sqrt(model.Fx * model.Fx + model.Fy * model.Fy + model.Fz * model.Fz);
                double dir_x = (force_mag > 0) ? model.Fx / force_mag : 0.0;
                double dir_y = (force_mag > 0) ? model.Fy / force_mag : 0.0;
                double dir_z = (force_mag > 0) ? model.Fz / force_mag : 1.0;
                if (force_mag == 0.0) force_mag = 1.0;

                double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);
                double current_perm = 1013.0 * h * h * mu * current_poro * flow_rate / force_mag;
                double current_tort = (vaz != 0) ? (v_tort / (vaz * count) - 1.0) : -1;

                stat_perm[l] = current_perm;
                stat_tort[l] = current_tort;
            }

            if (geometrical) {
                current_surf /= vol;
                current_surf /= h;
                stat_poro[l] = current_poro;
                stat_surf[l] = current_surf;
            }
        }

        if (hydraulic) {
            stat_convergence_metrics k_stats = stat_compute_convergence(stat_perm, stat_samples);
            stat_convergence_metrics t_stats = stat_compute_convergence(stat_tort, stat_samples);
            const bool has_tortuosity_samples =
                std::any_of(stat_tort.begin(), stat_tort.end(), [](double value) {
                    return value != -1.0;
                });
            if (stat_perm_rev.size == 0 && k_stats.CV >= 0.0 && k_stats.CV < gamma &&
                std::isfinite(k_stats.mean)) {
                stat_perm_rev = {x_size, {{x_size, sub_y, sub_z}}, k_stats.mean, k_stats.CV};
                // The existing statistical samples store T - 1, using -1 for
                // undefined flow. Do not turn an empty sample set into T = 1.
                if (std::isfinite(t_stats.mean) &&
                    std::any_of(stat_tort.begin(), stat_tort.end(), [](double value) {
                        return value != -1.0 && std::isfinite(value);
                    }))
                    stat_perm_rev.tortuosity_mean = t_stats.mean + 1.0;
            }
            hydraulic_log << x_size * h << " " << x_size << " "
                          << k_stats.CV << " " << k_stats.mean << " " << k_stats.max << " " << k_stats.min << " "
                          << t_stats.CV << " " << (has_tortuosity_samples ? t_stats.mean + 1 : 0) << " "
                          << t_stats.max + 1 << " " << (has_tortuosity_samples ? t_stats.min + 1 : 0) << "\n";
        }

        if (geometrical) {
            stat_convergence_metrics p_stats = stat_compute_convergence(stat_poro, stat_samples);
            stat_convergence_metrics s_stats = stat_compute_convergence(stat_surf, stat_samples);
            if (stat_poro_rev.size == 0 && p_stats.CV >= 0.0 && p_stats.CV < gamma &&
                std::isfinite(p_stats.mean)) {
                stat_poro_rev = {x_size, {{x_size, sub_y, sub_z}}, p_stats.mean, p_stats.CV,
                                 s_stats.mean};
            }
            geometric_log << x_size * h << " " << x_size << " "
                          << p_stats.CV << " " << p_stats.mean << " " << p_stats.max << " " << p_stats.min << " "
                          << s_stats.CV << " " << s_stats.mean  << " " << s_stats.max  << " " << s_stats.min  << "\n";
        }

        double min_step = std::max(1.0, average_pore_size);
        double proportional_step = current_size_x_double * size_step;
        current_size_x_double += std::max(min_step, proportional_step);
        x_size = static_cast<int>(std::ceil(current_size_x_double));
    }
    
    if (geometrical)
        geometric_log.close();
    if (hydraulic)
        hydraulic_log.close();
    printf("\n\nEnded Statistical Analysis\n\n");
}
