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
/* Flow adaptor class for multiphase flow methods */

#include "analysis/FlowAdaptor.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"

FlowAdaptor::FlowAdaptor(ScaLBL_ColorModel &M) {
    Nx = M.Dm->Nx;
    Ny = M.Dm->Ny;
    Nz = M.Dm->Nz;
    timestep = -1;
    timestep_previous = -1;

    phi.resize(Nx, Ny, Nz);
    phi.fill(0); // phase indicator field
    phi_t.resize(Nx, Ny, Nz);
    phi_t.fill(0); // time derivative for the phase indicator field
}

FlowAdaptor::~FlowAdaptor() {}

double FlowAdaptor::ImageInit(ScaLBL_ColorModel &M, std::string Filename) {
    int rank = M.rank;
    int Nx = M.Nx;
    int Ny = M.Ny;
    int Nz = M.Nz;
    if (rank == 0)
        printf("Re-initializing fluids from file: %s \n", Filename.c_str());
    M.Mask->Decomp(Filename);
    for (int i = 0; i < Nx * Ny * Nz; i++)
        M.id[i] = M.Mask->id[i]; // save what was read
    for (int i = 0; i < Nx * Ny * Nz; i++)
        M.Dm->id[i] = M.Mask->id[i]; // save what was read

    double *PhaseLabel;
    PhaseLabel = new double[Nx * Ny * Nz];
    M.AssignComponentLabels(PhaseLabel);

    double Count = 0.0;
    double PoreCount = 0.0;
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                if (M.id[Nx * Ny * k + Nx * j + i] == 2) {
                    PoreCount++;
                    Count++;
                } else if (M.id[Nx * Ny * k + Nx * j + i] == 1) {
                    PoreCount++;
                }
            }
        }
    }

    Count = M.Dm->Comm.sumReduce(Count);
    PoreCount = M.Dm->Comm.sumReduce(PoreCount);

    if (rank == 0)
        printf("   new saturation: %f (%f / %f) \n", Count / PoreCount, Count,
               PoreCount);
    ScaLBL_CopyToDevice(M.Phi, PhaseLabel, Nx * Ny * Nz * sizeof(double));
    M.Dm->Comm.barrier();

    ScaLBL_D3Q19_Init(M.fq, M.Np);
    ScaLBL_PhaseField_Init(M.dvcMap, M.Phi, M.Den, M.Aq, M.Bq, 0,
                           M.ScaLBL_Comm->LastExterior(), M.Np);
    ScaLBL_PhaseField_Init(M.dvcMap, M.Phi, M.Den, M.Aq, M.Bq,
                           M.ScaLBL_Comm->FirstInterior(),
                           M.ScaLBL_Comm->LastInterior(), M.Np);
    M.Dm->Comm.barrier();

    ScaLBL_CopyToHost(M.Averages->Phi.data(), M.Phi,
                      Nx * Ny * Nz * sizeof(double));

    delete PhaseLabel;
    double saturation = Count / PoreCount;
    return saturation;
}

double FlowAdaptor::UpdateFractionalFlow(ScaLBL_ColorModel &M) {

    double MASS_FRACTION_CHANGE = 0.006;
    double FRACTIONAL_FLOW_EPSILON = 5e-6;
    if (M.db->keyExists("FlowAdaptor")) {
        auto flow_db = M.db->getDatabase("FlowAdaptor");
        MASS_FRACTION_CHANGE =
            flow_db->getWithDefault<double>("mass_fraction_factor", 0.006);
        FRACTIONAL_FLOW_EPSILON =
            flow_db->getWithDefault<double>("fractional_flow_epsilon", 5e-6);
    }
    int Np = M.Np;
    double dA, dB, phi;
    double vx, vy, vz;
    double mass_a, mass_b, mass_a_global, mass_b_global;

    double *Aq_tmp, *Bq_tmp;
    double *Vel_x, *Vel_y, *Vel_z, *Phase;

    Aq_tmp = new double[7 * Np];
    Bq_tmp = new double[7 * Np];
    Phase = new double[Np];
    Vel_x = new double[Np];
    Vel_y = new double[Np];
    Vel_z = new double[Np];

    ScaLBL_CopyToHost(Aq_tmp, M.Aq, 7 * Np * sizeof(double));
    ScaLBL_CopyToHost(Bq_tmp, M.Bq, 7 * Np * sizeof(double));
    ScaLBL_CopyToHost(Vel_x, &M.Velocity[0], Np * sizeof(double));
    ScaLBL_CopyToHost(Vel_y, &M.Velocity[Np], Np * sizeof(double));
    ScaLBL_CopyToHost(Vel_z, &M.Velocity[2 * Np], Np * sizeof(double));

    int Nx = M.Nx;
    int Ny = M.Ny;
    int Nz = M.Nz;

    mass_a = mass_b = 0.0;
    double maxSpeed = 0.0;
    double localMaxSpeed = 0.0;
    /* compute mass change based on weights */
    double sum_weights_A = 0.0;
    double sum_weights_B = 0.0;
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int n = M.Map(i, j, k);
                //double distance = M.Averages->SDs(i,j,k);
                if (!(n < 0)) {
                    dA = Aq_tmp[n] + Aq_tmp[n + Np] + Aq_tmp[n + 2 * Np] +
                         Aq_tmp[n + 3 * Np] + Aq_tmp[n + 4 * Np] +
                         Aq_tmp[n + 5 * Np] + Aq_tmp[n + 6 * Np];
                    dB = Bq_tmp[n] + Bq_tmp[n + Np] + Bq_tmp[n + 2 * Np] +
                         Bq_tmp[n + 3 * Np] + Bq_tmp[n + 4 * Np] +
                         Bq_tmp[n + 5 * Np] + Bq_tmp[n + 6 * Np];
                    phi = (dA - dB) / (dA + dB);
                    Phase[n] = phi;
                    mass_a += dA;
                    mass_b += dB;
                    vx = Vel_x[n];
                    vy = Vel_y[n];
                    vz = Vel_z[n];
                    double local_momentum = sqrt(vx * vx + vy * vy + vz * vz);
                    double local_weight =
                        (FRACTIONAL_FLOW_EPSILON + local_momentum);
                    if (phi > 0.0) {
                        sum_weights_A += local_weight * dA;
                    } else {
                        sum_weights_B += local_weight * dB;
                    }
                    if (local_momentum > localMaxSpeed) {
                        localMaxSpeed = local_momentum;
                    }
                }
            }
        }
    }
    maxSpeed = M.Dm->Comm.maxReduce(localMaxSpeed);
    mass_a_global = M.Dm->Comm.sumReduce(mass_a);
    mass_b_global = M.Dm->Comm.sumReduce(mass_b);
    double sum_weights_A_global = M.Dm->Comm.sumReduce(sum_weights_A);
    double sum_weights_B_global = M.Dm->Comm.sumReduce(sum_weights_B);
    sum_weights_A_global /= (FRACTIONAL_FLOW_EPSILON + maxSpeed);
    sum_weights_B_global /= (FRACTIONAL_FLOW_EPSILON + maxSpeed);

    //double  total_momentum_A = sqrt(vax_global*vax_global+vay_global*vay_global+vaz_global*vaz_global);
    //double  total_momentum_B = sqrt(vbx_global*vbx_global+vby_global*vby_global+vbz_global*vbz_global);
    /* compute the total mass change */
    double TOTAL_MASS_CHANGE =
        MASS_FRACTION_CHANGE * (mass_a_global + mass_b_global);
    if (fabs(TOTAL_MASS_CHANGE) > 0.1 * mass_a_global)
        TOTAL_MASS_CHANGE = 0.1 * mass_a_global;
    if (fabs(TOTAL_MASS_CHANGE) > 0.1 * mass_b_global)
        TOTAL_MASS_CHANGE = 0.1 * mass_b_global;

    double MASS_FACTOR_A = TOTAL_MASS_CHANGE / sum_weights_A_global;
    double MASS_FACTOR_B = TOTAL_MASS_CHANGE / sum_weights_B_global;

    double LOCAL_MASS_CHANGE = 0.0;
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int n = M.Map(i, j, k);
                if (!(n < 0)) {
                    phi = Phase[n];
                    vx = Vel_x[n];
                    vy = Vel_y[n];
                    vz = Vel_z[n];
                    double local_momentum = sqrt(vx * vx + vy * vy + vz * vz);
                    double local_weight =
                        (FRACTIONAL_FLOW_EPSILON + local_momentum) /
                        (FRACTIONAL_FLOW_EPSILON + maxSpeed);
                    /* impose ceiling for spurious currents */
                    //if (local_momentum > maxSpeed) local_momentum =  maxSpeed;
                    if (phi > 0.0) {
                        LOCAL_MASS_CHANGE = MASS_FACTOR_A * local_weight;
                        Aq_tmp[n] -= 0.3333333333333333 * LOCAL_MASS_CHANGE;
                        Aq_tmp[n + Np] -=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        Aq_tmp[n + 2 * Np] -=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        Aq_tmp[n + 3 * Np] -=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        Aq_tmp[n + 4 * Np] -=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        Aq_tmp[n + 5 * Np] -=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        Aq_tmp[n + 6 * Np] -=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        //DebugMassA[n] = (-1.0)*LOCAL_MASS_CHANGE;
                    } else {
                        LOCAL_MASS_CHANGE = MASS_FACTOR_B * local_weight;
                        Bq_tmp[n] += 0.3333333333333333 * LOCAL_MASS_CHANGE;
                        Bq_tmp[n + Np] +=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        Bq_tmp[n + 2 * Np] +=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        Bq_tmp[n + 3 * Np] +=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        Bq_tmp[n + 4 * Np] +=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        Bq_tmp[n + 5 * Np] +=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        Bq_tmp[n + 6 * Np] +=
                            0.1111111111111111 * LOCAL_MASS_CHANGE;
                        //DebugMassB[n] = LOCAL_MASS_CHANGE;
                    }
                }
            }
        }
    }

    if (M.rank == 0)
        printf("Update Fractional Flow: change mass of fluid B by %f \n",
               TOTAL_MASS_CHANGE / mass_b_global);

    // Need to initialize Aq, Bq, Den, Phi directly
    //ScaLBL_CopyToDevice(Phi,phase.data(),7*Np*sizeof(double));
    ScaLBL_CopyToDevice(M.Aq, Aq_tmp, 7 * Np * sizeof(double));
    ScaLBL_CopyToDevice(M.Bq, Bq_tmp, 7 * Np * sizeof(double));

    delete Aq_tmp;
    delete Bq_tmp;
    delete Vel_x;
    delete Vel_y;
    delete Vel_z;
    delete Phase;

    return (TOTAL_MASS_CHANGE);
}

void FlowAdaptor::Flatten(ScaLBL_ColorModel &M) {

    ScaLBL_D3Q19_Init(M.fq, M.Np);
    ScaLBL_PhaseField_Init(M.dvcMap, M.Phi, M.Den, M.Aq, M.Bq, 0,
                           M.ScaLBL_Comm->LastExterior(), M.Np);
    ScaLBL_PhaseField_Init(M.dvcMap, M.Phi, M.Den, M.Aq, M.Bq,
                           M.ScaLBL_Comm->FirstInterior(),
                           M.ScaLBL_Comm->LastInterior(), M.Np);
}

double FlowAdaptor::MoveInterface(ScaLBL_ColorModel &M) {

    double INTERFACE_CUTOFF =
        M.color_db->getWithDefault<double>("move_interface_cutoff", 0.1);
    double MOVE_INTERFACE_FACTOR =
        M.color_db->getWithDefault<double>("move_interface_factor", 10.0);

    ScaLBL_CopyToHost(phi.data(), M.Phi, Nx * Ny * Nz * sizeof(double));
    /* compute the local derivative of phase indicator field */
    double beta = M.beta;
    double factor = 0.5 / beta;
    double total_interface_displacement = 0.0;
    double total_interface_sites = 0.0;
    for (int n = 0; n < Nx * Ny * Nz; n++) {
        /* compute the distance to the interface */
        double value1 = M.Averages->Phi(n);
        double dist1 = factor * log((1.0 + value1) / (1.0 - value1));
        double value2 = phi(n);
        double dist2 = factor * log((1.0 + value2) / (1.0 - value2));
        phi_t(n) = value2;
        if (value1 < INTERFACE_CUTOFF && value1 > -1 * INTERFACE_CUTOFF &&
            value2 < INTERFACE_CUTOFF && value2 > -1 * INTERFACE_CUTOFF) {
            /* time derivative of distance */
            double dxdt = 0.125 * (dist2 - dist1);
            /* extrapolate to move the distance further */
            double dist3 = dist2 + MOVE_INTERFACE_FACTOR * dxdt;
            /* compute the new phase interface */
            phi_t(n) = (2.f * (exp(-2.f * beta * (dist3))) /
                            (1.f + exp(-2.f * beta * (dist3))) -
                        1.f);
            total_interface_displacement += fabs(MOVE_INTERFACE_FACTOR * dxdt);
            total_interface_sites += 1.0;
        }
    }
    ScaLBL_CopyToDevice(M.Phi, phi_t.data(), Nx * Ny * Nz * sizeof(double));
    return total_interface_sites;
}

double FlowAdaptor::ShellAggregation(ScaLBL_ColorModel &M,
                                     const double target_delta_volume) {

    const RankInfoStruct rank_info(M.rank, M.nprocx, M.nprocy, M.nprocz);
    auto rank = M.rank;
    auto Nx = M.Nx;
    auto Ny = M.Ny;
    auto Nz = M.Nz;
    auto N = Nx * Ny * Nz;
    double vF = 0.f;
    double vS = 0.f;
    double delta_volume;
    double WallFactor = 1.0;
    bool USE_CONNECTED_NWP = false;

    DoubleArray phase(Nx, Ny, Nz);
    IntArray phase_label(Nx, Ny, Nz);
    ;
    DoubleArray phase_distance(Nx, Ny, Nz);
    Array<char> phase_id(Nx, Ny, Nz);
    fillHalo<double> fillDouble(M.Dm->Comm, M.Dm->rank_info,
                                {Nx - 2, Ny - 2, Nz - 2}, {1, 1, 1}, 0, 1);

    // Basic algorithm to
    // 1. Copy phase field to CPU
    ScaLBL_CopyToHost(phase.data(), M.Phi, N * sizeof(double));

    double count = 0.f;
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                if (phase(i, j, k) > 0.f && M.Averages->SDs(i, j, k) > 0.f)
                    count += 1.f;
            }
        }
    }
    double volume_initial = M.Dm->Comm.sumReduce(count);
    double PoreVolume = M.Dm->Volume * M.Dm->Porosity();
    /*ensure target isn't an absurdly small fraction of pore volume */
    if (volume_initial < target_delta_volume * PoreVolume) {
        volume_initial = target_delta_volume * PoreVolume;
    }

    // 2. Identify connected components of phase field -> phase_label

    double volume_connected = 0.0;
    double second_biggest = 0.0;
    if (USE_CONNECTED_NWP) {
        ComputeGlobalBlobIDs(Nx - 2, Ny - 2, Nz - 2, rank_info, phase,
                             M.Averages->SDs, vF, vS, phase_label, M.Dm->Comm);
        M.Dm->Comm.barrier();

        // only operate on component "0"ScaLBL_ColorModel &M,
        count = 0.0;

        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int label = phase_label(i, j, k);
                    if (label == 0) {
                        phase_id(i, j, k) = 0;
                        count += 1.0;
                    } else
                        phase_id(i, j, k) = 1;
                    if (label == 1) {
                        second_biggest += 1.0;
                    }
                }
            }
        }
        volume_connected = M.Dm->Comm.sumReduce(count);
        second_biggest = M.Dm->Comm.sumReduce(second_biggest);
    } else {
        // use the whole NWP
        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    if (M.Averages->SDs(i, j, k) > 0.f) {
                        if (phase(i, j, k) > 0.f) {
                            phase_id(i, j, k) = 0;
                        } else {
                            phase_id(i, j, k) = 1;
                        }
                    } else {
                        phase_id(i, j, k) = 1;
                    }
                }
            }
        }
    }

    // 3. Generate a distance map to the largest object -> phase_distance
    CalcDist(phase_distance, phase_id, *M.Dm);

    double temp, value;
    double factor = 0.5 / M.beta;
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                if (phase_distance(i, j, k) < 3.f) {
                    value = phase(i, j, k);
                    if (value > 1.f)
                        value = 1.f;
                    if (value < -1.f)
                        value = -1.f;
                    // temp -- distance based on analytical form McClure, Prins et al, Comp. Phys. Comm.
                    temp = -factor * log((1.0 + value) / (1.0 - value));
                    /// use this approximation close to the object
                    if (fabs(value) < 0.8 && M.Averages->SDs(i, j, k) > 1.f) {
                        phase_distance(i, j, k) = temp;
                    }
                    // erase the original object
                    phase(i, j, k) = -1.0;
                }
            }
        }
    }
    if (rank == 0)
        printf("Pathway volume / next largest ganglion %f \n",
               volume_connected / second_biggest);

    if (rank == 0)
        printf("MorphGrow with target volume fraction change %f \n",
               target_delta_volume / volume_initial);
    double target_delta_volume_incremental = target_delta_volume;
    if (fabs(target_delta_volume) > 0.01 * volume_initial)
        target_delta_volume_incremental = 0.01 * volume_initial *
                                          target_delta_volume /
                                          fabs(target_delta_volume);

    delta_volume =
        MorphGrow(M.Averages->SDs, phase_distance, phase_id, M.Averages->Dm,
                  target_delta_volume_incremental, WallFactor);

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                if (phase_distance(i, j, k) < 0.0)
                    phase_id(i, j, k) = 0;
                else
                    phase_id(i, j, k) = 1;
                //if (phase_distance(i,j,k) < 0.0 ) phase(i,j,k) = 1.0;
            }
        }
    }

    CalcDist(phase_distance, phase_id, *M.Dm); // re-calculate distance

    // 5. Update phase indicator field based on new distnace
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                double d = phase_distance(i, j, k);
                if (M.Averages->SDs(i, j, k) > 0.f) {
                    if (d < 3.f) {
                        //phase(i,j,k) = -1.0;
                        phase(i, j, k) = (2.f * (exp(-2.f * M.beta * d)) /
                                              (1.f + exp(-2.f * M.beta * d)) -
                                          1.f);
                    }
                }
            }
        }
    }
    fillDouble.fill(phase);

    count = 0.f;
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                if (phase(i, j, k) > 0.f && M.Averages->SDs(i, j, k) > 0.f) {
                    count += 1.f;
                }
            }
        }
    }
    double volume_final = M.Dm->Comm.sumReduce(count);

    delta_volume = (volume_final - volume_initial);
    if (rank == 0)
        printf("Shell Aggregation: change fluid volume fraction by %f \n",
               delta_volume / volume_initial);
    if (rank == 0)
        printf("   new saturation =  %f \n",
               volume_final /
                   (M.Mask->Porosity() *
                    double((Nx - 2) * (Ny - 2) * (Nz - 2) * M.nprocs)));

    // 6. copy back to the device
    //if (rank==0)  printf("MorphInit: copy data  back to device\n");
    ScaLBL_CopyToDevice(M.Phi, phase.data(), N * sizeof(double));

    // 7. Re-initialize phase field and density
    ScaLBL_PhaseField_Init(M.dvcMap, M.Phi, M.Den, M.Aq, M.Bq, 0,
                           M.ScaLBL_Comm->LastExterior(), M.Np);
    ScaLBL_PhaseField_Init(M.dvcMap, M.Phi, M.Den, M.Aq, M.Bq,
                           M.ScaLBL_Comm->FirstInterior(),
                           M.ScaLBL_Comm->LastInterior(), M.Np);
    auto BoundaryCondition = M.BoundaryCondition;
    if (BoundaryCondition == 1 || BoundaryCondition == 2 ||
        BoundaryCondition == 3 || BoundaryCondition == 4) {
        if (M.Dm->kproc() == 0) {
            ScaLBL_SetSlice_z(M.Phi, 1.0, Nx, Ny, Nz, 0);
            ScaLBL_SetSlice_z(M.Phi, 1.0, Nx, Ny, Nz, 1);
            ScaLBL_SetSlice_z(M.Phi, 1.0, Nx, Ny, Nz, 2);
        }
        if (M.Dm->kproc() == M.nprocz - 1) {
            ScaLBL_SetSlice_z(M.Phi, -1.0, Nx, Ny, Nz, Nz - 1);
            ScaLBL_SetSlice_z(M.Phi, -1.0, Nx, Ny, Nz, Nz - 2);
            ScaLBL_SetSlice_z(M.Phi, -1.0, Nx, Ny, Nz, Nz - 3);
        }
    }
    return delta_volume;
}

double FlowAdaptor::SeedPhaseField(ScaLBL_ColorModel &M,
                                   const double seed_water_in_oil) {
    srand(time(NULL));
    auto rank = M.rank;
    auto Np = M.Np;
    double mass_loss = 0.f;
    double count = 0.f;
    double *Aq_tmp, *Bq_tmp;

    Aq_tmp = new double[7 * Np];
    Bq_tmp = new double[7 * Np];

    ScaLBL_CopyToHost(Aq_tmp, M.Aq, 7 * Np * sizeof(double));
    ScaLBL_CopyToHost(Bq_tmp, M.Bq, 7 * Np * sizeof(double));

    for (int n = 0; n < M.ScaLBL_Comm->LastExterior(); n++) {
        double random_value = seed_water_in_oil * double(rand()) / RAND_MAX;
        double dA = Aq_tmp[n] + Aq_tmp[n + Np] + Aq_tmp[n + 2 * Np] +
                    Aq_tmp[n + 3 * Np] + Aq_tmp[n + 4 * Np] +
                    Aq_tmp[n + 5 * Np] + Aq_tmp[n + 6 * Np];
        double dB = Bq_tmp[n] + Bq_tmp[n + Np] + Bq_tmp[n + 2 * Np] +
                    Bq_tmp[n + 3 * Np] + Bq_tmp[n + 4 * Np] +
                    Bq_tmp[n + 5 * Np] + Bq_tmp[n + 6 * Np];
        double phase_id = (dA - dB) / (dA + dB);
        if (phase_id > 0.0) {
            Aq_tmp[n] -= 0.3333333333333333 * random_value;
            Aq_tmp[n + Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 2 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 3 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 4 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 5 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 6 * Np] -= 0.1111111111111111 * random_value;

            Bq_tmp[n] += 0.3333333333333333 * random_value;
            Bq_tmp[n + Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 2 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 3 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 4 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 5 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 6 * Np] += 0.1111111111111111 * random_value;
        }
        mass_loss += random_value * seed_water_in_oil;
    }

    for (int n = M.ScaLBL_Comm->FirstInterior();
         n < M.ScaLBL_Comm->LastInterior(); n++) {
        double random_value = seed_water_in_oil * double(rand()) / RAND_MAX;
        double dA = Aq_tmp[n] + Aq_tmp[n + Np] + Aq_tmp[n + 2 * Np] +
                    Aq_tmp[n + 3 * Np] + Aq_tmp[n + 4 * Np] +
                    Aq_tmp[n + 5 * Np] + Aq_tmp[n + 6 * Np];
        double dB = Bq_tmp[n] + Bq_tmp[n + Np] + Bq_tmp[n + 2 * Np] +
                    Bq_tmp[n + 3 * Np] + Bq_tmp[n + 4 * Np] +
                    Bq_tmp[n + 5 * Np] + Bq_tmp[n + 6 * Np];
        double phase_id = (dA - dB) / (dA + dB);
        if (phase_id > 0.0) {
            Aq_tmp[n] -= 0.3333333333333333 * random_value;
            Aq_tmp[n + Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 2 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 3 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 4 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 5 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 6 * Np] -= 0.1111111111111111 * random_value;

            Bq_tmp[n] += 0.3333333333333333 * random_value;
            Bq_tmp[n + Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 2 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 3 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 4 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 5 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 6 * Np] += 0.1111111111111111 * random_value;
        }
        mass_loss += random_value * seed_water_in_oil;
    }

    count = M.Dm->Comm.sumReduce(count);
    mass_loss = M.Dm->Comm.sumReduce(mass_loss);
    if (rank == 0)
        printf("Remove mass %f from %f voxels \n", mass_loss, count);

    // Need to initialize Aq, Bq, Den, Phi directly
    //ScaLBL_CopyToDevice(Phi,phase.data(),7*Np*sizeof(double));
    ScaLBL_CopyToDevice(M.Aq, Aq_tmp, 7 * Np * sizeof(double));
    ScaLBL_CopyToDevice(M.Bq, Bq_tmp, 7 * Np * sizeof(double));

    delete Aq_tmp;
    delete Bq_tmp;
    return (mass_loss);
}
