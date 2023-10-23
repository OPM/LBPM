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
#include "analysis/SubPhase.h"

// Constructor
SubPhase::SubPhase(std::shared_ptr<Domain> dm) : Dm(dm) {
    Nx = dm->Nx;
    Ny = dm->Ny;
    Nz = dm->Nz;
    Volume = (Nx - 2) * (Ny - 2) * (Nz - 2) * Dm->nprocx() * Dm->nprocy() *
             Dm->nprocz() * 1.0;

    morph_w = std::shared_ptr<Minkowski>(new Minkowski(Dm));
    morph_n = std::shared_ptr<Minkowski>(new Minkowski(Dm));
    morph_i = std::shared_ptr<Minkowski>(new Minkowski(Dm));

    // Global arrays
    PhaseID.resize(Nx, Ny, Nz);
    PhaseID.fill(0);
    Label_WP.resize(Nx, Ny, Nz);
    Label_WP.fill(0);
    Label_NWP.resize(Nx, Ny, Nz);
    Label_NWP.fill(0);
    Rho_n.resize(Nx, Ny, Nz);
    Rho_n.fill(0);
    Rho_w.resize(Nx, Ny, Nz);
    Rho_w.fill(0);
    Pressure.resize(Nx, Ny, Nz);
    Pressure.fill(0);
    Phi.resize(Nx, Ny, Nz);
    Phi.fill(0);
    DelPhi.resize(Nx, Ny, Nz);
    DelPhi.fill(0);
    Vel_x.resize(Nx, Ny, Nz);
    Vel_x.fill(0); // Gradient of the phase indicator field
    Vel_y.resize(Nx, Ny, Nz);
    Vel_y.fill(0);
    Vel_z.resize(Nx, Ny, Nz);
    Vel_z.fill(0);
    Dissipation.resize(Nx, Ny, Nz);
    Dissipation.fill(0);
    SDs.resize(Nx, Ny, Nz);
    SDs.fill(0);
    //.........................................

    //.........................................
    if (Dm->rank() == 0) {
        bool WriteHeader = false;
        SUBPHASE = fopen("subphase.csv", "r");
        if (SUBPHASE != NULL)
            fclose(SUBPHASE);
        else
            WriteHeader = true;

        SUBPHASE = fopen("subphase.csv", "a+");
        if (WriteHeader) {
            // If timelog is empty, write a short header to list the averages
            fprintf(SUBPHASE, "time rn rw nun nuw Fx Fy Fz iftwn wet ");
            fprintf(SUBPHASE, "pwc pwd pnc pnd ");                 // pressures
            fprintf(SUBPHASE, "Mwc Mwd Mwi Mnc Mnd Mni Msw Msn "); // mass
            fprintf(
                SUBPHASE,
                "Pwc_x Pwd_x Pwi_x Pnc_x Pnd_x Pni_x Psw_x Psn_x "); // momentum
            fprintf(SUBPHASE,
                    "Pwc_y Pwd_y Pwi_y Pnc_y Pnd_y Pni_y Psw_y Psn_y ");
            fprintf(SUBPHASE,
                    "Pwc_z Pwd_z Pwi_z Pnc_z Pnd_z Pni_z Psw_z Psn_z ");
            fprintf(SUBPHASE, "Kwc Kwd Kwi Knc Knd Kni "); // kinetic energy
            fprintf(SUBPHASE, "Dwc Dwd Dnc Dnd ");      // viscous dissipation
            fprintf(SUBPHASE, "Vwc Awc Hwc Xwc ");      // wc region
            fprintf(SUBPHASE, "Vwd Awd Hwd Xwd Nwd ");  // wd region
            fprintf(SUBPHASE, "Vnc Anc Hnc Xnc ");      // nc region
            fprintf(SUBPHASE, "Vnd And Hnd Xnd Nnd ");  // nd regionin
            fprintf(SUBPHASE, "Vi Ai Hi Xi ");          // interface region
            fprintf(SUBPHASE, "Vic Aic Hic Xic Nic\n"); // interface region

            // stress tensor?
        }

    } else {
        char LocalRankString[8];
        sprintf(LocalRankString, "%05d", Dm->rank());
        char LocalRankFilename[40];
        sprintf(LocalRankFilename, "%s%s", "subphase.csv.", LocalRankString);
        SUBPHASE = fopen(LocalRankFilename, "a+");

        fprintf(SUBPHASE, "time rn rw nun nuw Fx Fy Fz iftwn wet ");
        fprintf(SUBPHASE, "pwc pwd pnc pnd ");                 // pressures
        fprintf(SUBPHASE, "Mwc Mwd Mwi Mnc Mnd Mni Msw Msn "); // mass
        fprintf(SUBPHASE,
                "Pwc_x Pwd_x Pwi_x Pnc_x Pnd_x Pni_x Psw_x Psn_x "); // momentum
        fprintf(SUBPHASE, "Pwc_y Pwd_y Pwi_y Pnc_y Pnd_y Pni_y Psw_y Psn_y ");
        fprintf(SUBPHASE, "Pwc_z Pwd_z Pwi_z Pnc_z Pnd_z Pni_z Psw_z Psn_z ");
        fprintf(SUBPHASE, "Kwc Kwd Kwi Knc Knd Kni "); // kinetic energy
        fprintf(SUBPHASE, "Dwc Dwd Dnc Dnd ");         // viscous dissipation
        fprintf(SUBPHASE, "Vwc Awc Hwc Xwc ");         // wc region
        fprintf(SUBPHASE, "Vwd Awd Hwd Xwd Nwd ");     // wd region
        fprintf(SUBPHASE, "Vnc Anc Hnc Xnc ");         // nc region
        fprintf(SUBPHASE, "Vnd And Hnd Xnd Nnd ");     // nd region
        fprintf(SUBPHASE, "Vi Ai Hi Xi ");             // interface region
        fprintf(SUBPHASE, "Vic Aic Hic Xic Nic\n");    // interface region
    }

    if (Dm->rank() == 0) {
        bool WriteHeader = false;
        TIMELOG = fopen("timelog.csv", "r");
        if (TIMELOG != NULL)
            fclose(TIMELOG);
        else
            WriteHeader = true;

        TIMELOG = fopen("timelog.csv", "a+");
        if (WriteHeader) {
            // If timelog is empty, write a short header to list the averages
            fprintf(TIMELOG,
                    "sw krw krn krwf krnf vw vn force pw pn wet peff\n");
        }
    }
}

// Destructor
SubPhase::~SubPhase() {
    if (SUBPHASE != NULL) {
        fclose(SUBPHASE);
    }
}

void SubPhase::Write(int timestep) {
    if (Dm->rank() == 0) {
        fprintf(SUBPHASE, "%i %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ",
                timestep, rho_n, rho_w, nu_n, nu_w, Fx, Fy, Fz, gamma_wn,
                total_wetting_interaction_global);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g ", gwc.p, gwd.p, gnc.p, gnd.p);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ", gwc.M,
                gwd.M, giwn.Mw, gnc.M, gnd.M, giwn.Mn, gifs.Mw, gifs.Mn);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ", gwc.Px,
                gwd.Px, giwn.Pwx, gnc.Px, gnd.Px, giwn.Pnx, gifs.Pwx, gifs.Pnx);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ", gwc.Py,
                gwd.Py, giwn.Pwy, gnc.Py, gnd.Py, giwn.Pny, gifs.Pwy, gifs.Pny);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ", gwc.Pz,
                gwd.Pz, giwn.Pwz, gnc.Pz, gnd.Pz, giwn.Pnz, gifs.Pwz, gifs.Pnz);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %.8g %.8g ", gwc.K, gwd.K,
                giwn.Kw, gnc.K, gnd.K, giwn.Kn);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g ", gwc.visc, gwd.visc, gnc.visc,
                gnd.visc);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g ", gwc.V, gwc.A, gwc.H, gwc.X);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %i ", gwd.V, gwd.A, gwd.H, gwd.X,
                gwd.Nc);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g ", gnc.V, gnc.A, gnc.H, gnc.X);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %i ", gnd.V, gnd.A, gnd.H, gnd.X,
                gnd.Nc);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g ", giwn.V, giwn.A, giwn.H,
                giwn.X);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %i\n", giwnc.V, giwnc.A, giwnc.H,
                giwnc.X, giwnc.Nc);
        fflush(SUBPHASE);
    } else {
        fprintf(SUBPHASE, "%i %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ",
                timestep, rho_n, rho_w, nu_n, nu_w, Fx, Fy, Fz, gamma_wn,
                total_wetting_interaction);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g ", wc.p, wd.p, nc.p, nd.p);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ", wc.M,
                wd.M, iwn.Mw, nc.M, nd.M, iwn.Mn, ifs.Mw, ifs.Mn);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ", wc.Px,
                wd.Px, iwn.Pwx, nc.Px, nd.Px, iwn.Pnx, ifs.Pwx, ifs.Pnx);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ", wc.Py,
                wd.Py, iwn.Pwy, nc.Py, nd.Py, iwn.Pny, ifs.Pwy, ifs.Pny);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ", wc.Pz,
                wd.Pz, iwn.Pwz, nc.Pz, nd.Pz, iwn.Pnz, ifs.Pwz, ifs.Pnz);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %.8g %.8g ", wc.K, wd.K, iwn.Kw,
                nc.K, nd.K, iwn.Kn);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g ", wc.visc, wd.visc, nc.visc,
                nd.visc);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g ", wc.V, wc.A, wc.H, wc.X);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %i ", wd.V, wd.A, wd.H, wd.X,
                wd.Nc);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g ", nc.V, nc.A, nc.H, nc.X);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g %i ", nd.V, nd.A, nd.H, nd.X,
                nd.Nc);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g ", iwn.V, iwn.A, iwn.H, iwn.X);
        fprintf(SUBPHASE, "%.8g %.8g %.8g %.8g\n", iwnc.V, iwnc.A, iwnc.H,
                iwnc.X);
    }
}

void SubPhase::SetParams(double rhoA, double rhoB, double tauA, double tauB,
                         double force_x, double force_y, double force_z,
                         double alpha, double B) {
    Fx = force_x;
    Fy = force_y;
    Fz = force_z;
    rho_n = rhoA;
    rho_w = rhoB;
    nu_n = (tauA - 0.5) / 3.f;
    nu_w = (tauB - 0.5) / 3.f;
    gamma_wn = 5.796 * alpha;
    beta = B;
}

void SubPhase::Basic() {
    int i, j, k, n, imin, jmin, kmin, kmax;

    // If external boundary conditions are set, do not average over the inlet
    kmin = 1;
    kmax = Nz - 1;
    imin = jmin = 1;

    nb.reset();
    wb.reset();
    iwn.reset();

    double count_w = 0.0;
    double count_n = 0.0;

    /* compute the laplacian */
    Dm->CommunicateMeshHalo(Phi);
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                // Compute all of the derivatives using finite differences
                double fx = 0.5 * (Phi(i + 1, j, k) - Phi(i - 1, j, k));
                double fy = 0.5 * (Phi(i, j + 1, k) - Phi(i, j - 1, k));
                double fz = 0.5 * (Phi(i, j, k + 1) - Phi(i, j, k - 1));
                DelPhi(i, j, k) = sqrt(fx * fx + fy * fy + fz * fz);
            }
        }
    }
    Dm->CommunicateMeshHalo(DelPhi);

    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                // Compute volume averages
                if (Dm->id[n] > 0) {
                    // compute density
                    double nA = Rho_n(n);
                    double nB = Rho_w(n);
                    double phi = (nA - nB) / (nA + nB);
                    Phi(n) = phi;
                }
                if (Phi(n) != Phi(n)) {
                    // check for NaN
                    Phi(n) = 0.0;
                    //printf("Nan at %i %i %i \n",i,j,k);
                }
            }
        }
    }

    for (k = kmin; k < kmax; k++) {
        for (j = jmin; j < Ny - 1; j++) {
            for (i = imin; i < Nx - 1; i++) {
                n = k * Nx * Ny + j * Nx + i;
                // Compute volume averages
                if (Dm->id[n] > 0) {
                    // compute density
                    double nA = Rho_n(n);
                    double nB = Rho_w(n);
                    double phi = (nA - nB) / (nA + nB);
                    if (phi > 0.0) {
                        nA = 1.0;
                        nb.V += 1.0;
                        nb.M += nA * rho_n;
                        // velocity
                        nb.Px += rho_n * nA * Vel_x(n);
                        nb.Py += rho_n * nA * Vel_y(n);
                        nb.Pz += rho_n * nA * Vel_z(n);
                    } else {
                        nB = 1.0;
                        wb.M += nB * rho_w;
                        wb.V += 1.0;

                        // velocity
                        wb.Px += rho_w * nB * Vel_x(n);
                        wb.Py += rho_w * nB * Vel_y(n);
                        wb.Pz += rho_w * nB * Vel_z(n);
                    }
                    if (phi > 0.99) {
                        nb.p += Pressure(n);
                        count_n += 1.0;
                    } else if (phi < -0.99) {
                        wb.p += Pressure(n);
                        count_w += 1.0;
                    }
                    /* compute the film contribution */
                    else if (SDs(i, j, k) < 2.0) {
                        if (phi > 0.0) {
                            nA = 1.0;
                            iwn.V += 1.0;
                            iwn.Mn += nA * rho_n;
                            // velocity
                            iwn.Pnx += rho_n * nA * Vel_x(n);
                            iwn.Pny += rho_n * nA * Vel_y(n);
                            iwn.Pnz += rho_n * nA * Vel_z(n);
                        } else {
                            nB = 1.0;
                            iwn.Mw += nB * rho_w;
                            iwn.V += 1.0;

                            iwn.Pwx += rho_w * nB * Vel_x(n);
                            iwn.Pwy += rho_w * nB * Vel_y(n);
                            iwn.Pwz += rho_w * nB * Vel_z(n);
                        }
                    }
                }
            }
        }
    }

    total_wetting_interaction = count_wetting_interaction = 0.0;
    total_wetting_interaction_global = count_wetting_interaction_global = 0.0;
    for (k = kmin; k < kmax; k++) {
        for (j = jmin; j < Ny - 1; j++) {
            for (i = imin; i < Nx - 1; i++) {
                n = k * Nx * Ny + j * Nx + i;
                // compute contribution of wetting terms (within two voxels of solid)
                if (Dm->id[n] > 0 && SDs(i, j, k) < 2.0) {
                    count_wetting_interaction += 1.0;
                    total_wetting_interaction += DelPhi(i, j, k);
                }
            }
        }
    }

    total_wetting_interaction_global =
        Dm->Comm.sumReduce(total_wetting_interaction);
    count_wetting_interaction_global =
        Dm->Comm.sumReduce(count_wetting_interaction);

    gwb.V = Dm->Comm.sumReduce(wb.V);
    gnb.V = Dm->Comm.sumReduce(nb.V);
    gwb.M = Dm->Comm.sumReduce(wb.M);
    gnb.M = Dm->Comm.sumReduce(nb.M);
    gwb.Px = Dm->Comm.sumReduce(wb.Px);
    gwb.Py = Dm->Comm.sumReduce(wb.Py);
    gwb.Pz = Dm->Comm.sumReduce(wb.Pz);
    gnb.Px = Dm->Comm.sumReduce(nb.Px);
    gnb.Py = Dm->Comm.sumReduce(nb.Py);
    gnb.Pz = Dm->Comm.sumReduce(nb.Pz);

    giwn.Mw = Dm->Comm.sumReduce(iwn.Mw);
    giwn.Pwx = Dm->Comm.sumReduce(iwn.Pwx);
    giwn.Pwy = Dm->Comm.sumReduce(iwn.Pwy);
    giwn.Pwz = Dm->Comm.sumReduce(iwn.Pwz);

    giwn.Mn = Dm->Comm.sumReduce(iwn.Mn);
    giwn.Pnx = Dm->Comm.sumReduce(iwn.Pnx);
    giwn.Pny = Dm->Comm.sumReduce(iwn.Pny);
    giwn.Pnz = Dm->Comm.sumReduce(iwn.Pnz);

    count_w = Dm->Comm.sumReduce(count_w);
    count_n = Dm->Comm.sumReduce(count_n);
    if (count_w > 0.0)
        gwb.p = Dm->Comm.sumReduce(wb.p) / count_w;
    else
        gwb.p = 0.0;
    if (count_n > 0.0)
        gnb.p = Dm->Comm.sumReduce(nb.p) / count_n;
    else
        gnb.p = 0.0;

    // check for NaN
    bool err = false;
    if (gwb.V != gwb.V)
        err = true;
    if (gnb.V != gnb.V)
        err = true;
    if (gwb.p != gwb.p)
        err = true;
    if (gnb.p != gnb.p)
        err = true;
    if (gwb.Px != gwb.Px)
        err = true;
    if (gwb.Py != gwb.Py)
        err = true;
    if (gwb.Pz != gwb.Pz)
        err = true;
    if (gnb.Px != gnb.Px)
        err = true;
    if (gnb.Py != gnb.Py)
        err = true;
    if (gnb.Pz != gnb.Pz)
        err = true;

    if (Dm->rank() == 0) {
        /* align flow direction based on total mass flux */
        double dir_x = gwb.Px + gnb.Px;
        double dir_y = gwb.Py + gnb.Py;
        double dir_z = gwb.Pz + gnb.Pz;
        double flow_magnitude = dir_x * dir_x + dir_y * dir_y + dir_z * dir_z;
        double force_mag = sqrt(Fx * Fx + Fy * Fy + Fz * Fz);
        if (force_mag > 0.0) {
            dir_x = Fx / force_mag;
            dir_y = Fy / force_mag;
            dir_z = Fz / force_mag;
        } else {
            dir_x /= flow_magnitude;
            dir_y /= flow_magnitude;
            dir_z /= flow_magnitude;
        }
        if (Dm->BoundaryCondition == 1 || Dm->BoundaryCondition == 2 ||
            Dm->BoundaryCondition == 3 || Dm->BoundaryCondition == 4) {
            // compute the pressure drop
            double pressure_drop = (Pressure(Nx * Ny + Nx + 1) - 1.0) / 3.0;
            double length = ((Nz - 2) * Dm->nprocz());
            force_mag -= pressure_drop / length;
        }
        if (force_mag == 0.0 && flow_magnitude == 0.0) {
            // default to z direction
            dir_x = 0.0;
            dir_y = 0.0;
            dir_z = 1.0;
            force_mag = 1.0;
        }
        double Porosity = (gwb.V + gnb.V) / Dm->Volume;
        double saturation = gwb.V / (gwb.V + gnb.V);
        double water_flow_rate =
            gwb.V * (gwb.Px * dir_x + gwb.Py * dir_y + gwb.Pz * dir_z) / gwb.M /
            Dm->Volume;
        double not_water_flow_rate =
            gnb.V * (gnb.Px * dir_x + gnb.Py * dir_y + gnb.Pz * dir_z) / gnb.M /
            Dm->Volume;

        /* contribution from water films */
        double water_film_flow_rate =
            gwb.V * (giwn.Pwx * dir_x + giwn.Pwy * dir_y + giwn.Pwz * dir_z) /
            gwb.M / Dm->Volume;
        double not_water_film_flow_rate =
            gnb.V * (giwn.Pnx * dir_x + giwn.Pny * dir_y + giwn.Pnz * dir_z) /
            gnb.M / Dm->Volume;
        //double total_flow_rate = water_flow_rate + not_water_flow_rate;
        //double fractional_flow = water_flow_rate / total_flow_rate;
        double h = Dm->voxel_length;
        double krn = h * h * nu_n * Porosity * not_water_flow_rate / force_mag;
        double krw = h * h * nu_w * Porosity * water_flow_rate / force_mag;
        /* not counting films */
        double krnf = krn - h * h * nu_n * Porosity * not_water_film_flow_rate /
                                force_mag;
        double krwf =
            krw - h * h * nu_w * Porosity * water_film_flow_rate / force_mag;
        double eff_pressure = 1.0 / (krn + krw); // effective pressure drop

        fprintf(TIMELOG,
                "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
                saturation, krw, krn, krwf, krnf, h * water_flow_rate,
                h * not_water_flow_rate, force_mag, gwb.p, gnb.p,
                total_wetting_interaction_global, eff_pressure);
        fflush(TIMELOG);
    }
    if (err == true) {
        // exception if simulation produceds NaN
        printf("SubPhase.cpp: NaN encountered, may need to check simulation "
               "parameters \n");
    }
    ASSERT(err == false);
}

inline void InterfaceTransportMeasures(double beta, double rA, double rB,
                                       double nA, double nB, double nx,
                                       double ny, double nz, double ux,
                                       double uy, double uz, interface &I) {

    double A1, A2, A3, A4, A5, A6;
    double B1, B2, B3, B4, B5, B6;
    double nAB, delta;
    double phi = (nA - nB) / (nA + nB);
    // Instantiate mass transport distributions
    // Stationary value - distribution 0
    nAB = 1.0 / (nA + nB);
    //...............................................
    // q = 0,2,4
    // Cq = {1,0,0}, {0,1,0}, {0,0,1}
    delta = beta * nA * nB * nAB * 0.1111111111111111 * nx;
    if (!(nA * nB * nAB > 0))
        delta = 0;
    A1 = nA * (0.1111111111111111 * (1 + 4.5 * ux)) + delta;
    B1 = nB * (0.1111111111111111 * (1 + 4.5 * ux)) - delta;
    A2 = nA * (0.1111111111111111 * (1 - 4.5 * ux)) - delta;
    B2 = nB * (0.1111111111111111 * (1 - 4.5 * ux)) + delta;

    //...............................................
    // Cq = {0,1,0}
    delta = beta * nA * nB * nAB * 0.1111111111111111 * ny;
    if (!(nA * nB * nAB > 0))
        delta = 0;
    A3 = nA * (0.1111111111111111 * (1 + 4.5 * uy)) + delta;
    B3 = nB * (0.1111111111111111 * (1 + 4.5 * uy)) - delta;
    A4 = nA * (0.1111111111111111 * (1 - 4.5 * uy)) - delta;
    B4 = nB * (0.1111111111111111 * (1 - 4.5 * uy)) + delta;

    //...............................................
    // q = 4
    // Cq = {0,0,1}
    delta = beta * nA * nB * nAB * 0.1111111111111111 * nz;
    if (!(nA * nB * nAB > 0))
        delta = 0;
    A5 = nA * (0.1111111111111111 * (1 + 4.5 * uz)) + delta;
    B5 = nB * (0.1111111111111111 * (1 + 4.5 * uz)) - delta;
    A6 = nA * (0.1111111111111111 * (1 - 4.5 * uz)) - delta;
    B6 = nB * (0.1111111111111111 * (1 - 4.5 * uz)) + delta;

    double unx = (A1 - A2);
    double uny = (A3 - A4);
    double unz = (A5 - A6);
    double uwx = (B1 - B2);
    double uwy = (B3 - B4);
    double uwz = (B5 - B6);
    /*
	I.Mn += rA*nA;
	I.Mw += rB*nB;
	I.Pnx += rA*nA*unx;
	I.Pny += rA*nA*uny;
	I.Pnz += rA*nA*unz;
	I.Pwx += rB*nB*uwx;
	I.Pwy += rB*nB*uwy;
	I.Pwz += rB*nB*uwz;
	I.Kn += rA*nA*(unx*unx + uny*uny + unz*unz);
	I.Kw += rB*nB*(uwx*uwx + uwy*uwy + uwz*uwz);
	*/
    if (phi > 0.0) {
        I.Mn += rA;
        I.Pnx += rA * ux;
        I.Pny += rA * uy;
        I.Pnz += rA * uz;
    } else {
        I.Mw += rB;
        I.Pwx += rB * ux;
        I.Pwy += rB * uy;
        I.Pwz += rB * uz;
    }
    I.Kn += rA * nA * (unx * unx + uny * uny + unz * unz);
    I.Kw += rB * nB * (uwx * uwx + uwy * uwy + uwz * uwz);
}

void SubPhase::Full() {
    int i, j, k, n, imin, jmin, kmin, kmax;

    // If external boundary conditions are set, do not average over the inlet
    kmin = 1;
    kmax = Nz - 1;
    imin = jmin = 1;

    nd.reset();
    nc.reset();
    wd.reset();
    wc.reset();
    iwn.reset();
    iwnc.reset();
    ifs.reset();

    Dm->CommunicateMeshHalo(Phi);
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                // Compute all of the derivatives using finite differences
                double fx = 0.5 * (Phi(i + 1, j, k) - Phi(i - 1, j, k));
                double fy = 0.5 * (Phi(i, j + 1, k) - Phi(i, j - 1, k));
                double fz = 0.5 * (Phi(i, j, k + 1) - Phi(i, j, k - 1));
                DelPhi(i, j, k) = sqrt(fx * fx + fy * fy + fz * fz);
            }
        }
    }
    Dm->CommunicateMeshHalo(DelPhi);

    Dm->CommunicateMeshHalo(Vel_x);
    Dm->CommunicateMeshHalo(Vel_y);
    Dm->CommunicateMeshHalo(Vel_z);
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                // Compute velocity gradients using finite differences
                double phi = Phi(i, j, k);
                double nu = nu_n + 0.5 * (1.0 - phi) * (nu_w - nu_n);
                double rho = rho_n + 0.5 * (1.0 - phi) * (rho_w - rho_n);
                double ux = 0.5 * (Vel_x(i + 1, j, k) - Vel_x(i - 1, j, k));
                double uy = 0.5 * (Vel_x(i, j + 1, k) - Vel_x(i, j - 1, k));
                double uz = 0.5 * (Vel_x(i, j, k + 1) - Vel_x(i, j, k - 1));
                double vx = 0.5 * (Vel_y(i + 1, j, k) - Vel_y(i - 1, j, k));
                double vy = 0.5 * (Vel_y(i, j + 1, k) - Vel_y(i, j - 1, k));
                double vz = 0.5 * (Vel_y(i, j, k + 1) - Vel_y(i, j, k - 1));
                double wx = 0.5 * (Vel_z(i + 1, j, k) - Vel_z(i - 1, j, k));
                double wy = 0.5 * (Vel_z(i, j + 1, k) - Vel_z(i, j - 1, k));
                double wz = 0.5 * (Vel_z(i, j, k + 1) - Vel_z(i, j, k - 1));
                if (SDs(i, j, k) > 2.0) {
                    Dissipation(i, j, k) = 2 * rho * nu *
                                           (ux * ux + vy * vy + wz * wz +
                                            0.5 * (vx + uy) * (vx + uy) +
                                            0.5 * (vz + wy) * (vz + wy) +
                                            0.5 * (uz + wx) * (uz + wx));
                }
            }
        }
    }

    /*  Set up geometric analysis of each region */

    // non-wetting
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                if (SDs(n) <= 0.0) {
                    // Solid phase
                    morph_n->id(i, j, k) = 1;

                } else if (Phi(n) > 0.0) {
                    // non-wetting phase
                    morph_n->id(i, j, k) = 0;
                } else {
                    // wetting phase
                    morph_n->id(i, j, k) = 1;
                }
            }
        }
    }
    // measure the whole object
    morph_n->MeasureObject(); //0.5/beta,Phi);
    nd.V = morph_n->V();
    nd.A = morph_n->A();
    nd.H = morph_n->H();
    nd.X = morph_n->X();
    // measure only the connected part
    nd.Nc = morph_n->MeasureConnectedPathway(); //0.5/beta,Phi);
    nc.V = morph_n->V();
    nc.A = morph_n->A();
    nc.H = morph_n->H();
    nc.X = morph_n->X();
    // update disconnected part
    nd.V -= nc.V;
    nd.A -= nc.A;
    nd.H -= nc.H;
    nd.X -= nc.X;

    // compute global entities
    gnc.V = Dm->Comm.sumReduce(nc.V);
    gnc.A = Dm->Comm.sumReduce(nc.A);
    gnc.H = Dm->Comm.sumReduce(nc.H);
    gnc.X = Dm->Comm.sumReduce(nc.X);
    gnd.V = Dm->Comm.sumReduce(nd.V);
    gnd.A = Dm->Comm.sumReduce(nd.A);
    gnd.H = Dm->Comm.sumReduce(nd.H);
    gnd.X = Dm->Comm.sumReduce(nd.X);
    gnd.Nc = nd.Nc;
    // wetting
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                if (SDs(n) <= 0.0) {
                    // Solid phase
                    morph_w->id(i, j, k) = 1;

                } else if (Phi(n) < 0.0) {
                    // wetting phase
                    morph_w->id(i, j, k) = 0;
                } else {
                    // non-wetting phase
                    morph_w->id(i, j, k) = 1;
                }
            }
        }
    }
    morph_w->MeasureObject(); //-0.5/beta,Phi);
    wd.V = morph_w->V();
    wd.A = morph_w->A();
    wd.H = morph_w->H();
    wd.X = morph_w->X();
    // measure only the connected part
    wd.Nc = morph_w->MeasureConnectedPathway(); //-0.5/beta,Phi);
    wc.V = morph_w->V();
    wc.A = morph_w->A();
    wc.H = morph_w->H();
    wc.X = morph_w->X();
    // update disconnected part
    wd.V -= wc.V;
    wd.A -= wc.A;
    wd.H -= wc.H;
    wd.X -= wc.X;
    // compute global entities
    gwc.V = Dm->Comm.sumReduce(wc.V);
    gwc.A = Dm->Comm.sumReduce(wc.A);
    gwc.H = Dm->Comm.sumReduce(wc.H);
    gwc.X = Dm->Comm.sumReduce(wc.X);
    gwd.V = Dm->Comm.sumReduce(wd.V);
    gwd.A = Dm->Comm.sumReduce(wd.A);
    gwd.H = Dm->Comm.sumReduce(wd.H);
    gwd.X = Dm->Comm.sumReduce(wd.X);
    gwd.Nc = wd.Nc;

    /*  Set up geometric analysis of interface region */
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                if (SDs(n) <= 0.0) {
                    // Solid phase
                    morph_i->id(i, j, k) = 1;
                } else if (DelPhi(n) > 1e-4) {
                    // interface
                    morph_i->id(i, j, k) = 0;
                } else {
                    // not interface
                    morph_i->id(i, j, k) = 1;
                }
            }
        }
    }
    morph_i->MeasureObject();
    iwn.V = morph_i->V();
    iwn.A = morph_i->A();
    iwn.H = morph_i->H();
    iwn.X = morph_i->X();
    giwn.V = Dm->Comm.sumReduce(iwn.V);
    giwn.A = Dm->Comm.sumReduce(iwn.A);
    giwn.H = Dm->Comm.sumReduce(iwn.H);
    giwn.X = Dm->Comm.sumReduce(iwn.X);
    // measure only the connected part
    iwnc.Nc = morph_i->MeasureConnectedPathway();
    iwnc.V = morph_i->V();
    iwnc.A = morph_i->A();
    iwnc.H = morph_i->H();
    iwnc.X = morph_i->X();
    giwnc.V = Dm->Comm.sumReduce(iwnc.V);
    giwnc.A = Dm->Comm.sumReduce(iwnc.A);
    giwnc.H = Dm->Comm.sumReduce(iwnc.H);
    giwnc.X = Dm->Comm.sumReduce(iwnc.X);
    giwnc.Nc = iwnc.Nc;

    double vol_nc_bulk = 0.0;
    double vol_wc_bulk = 0.0;
    double vol_nd_bulk = 0.0;
    double vol_wd_bulk = 0.0;
    for (k = kmin; k < kmax; k++) {
        for (j = jmin; j < Ny - 1; j++) {
            for (i = imin; i < Nx - 1; i++) {
                n = k * Nx * Ny + j * Nx + i;
                // Compute volume averages
                if (SDs(n) > 0.0) {
                    // compute density
                    double nA = Rho_n(n);
                    double nB = Rho_w(n);
                    double phi = (nA - nB) / (nA + nB);
                    double ux = Vel_x(n);
                    double uy = Vel_y(n);
                    double uz = Vel_z(n);
                    double visc = Dissipation(n);

                    if (DelPhi(n) > 1e-3) {
                        // get the normal vector
                        double nx = 0.5 * (Phi(i + 1, j, k) - Phi(i - 1, j, k));
                        double ny = 0.5 * (Phi(i, j + 1, k) - Phi(i, j - 1, k));
                        double nz = 0.5 * (Phi(i, j, k + 1) - Phi(i, j, k - 1));
                        if (SDs(n) > 2.5) {
                            // not a film region
                            InterfaceTransportMeasures(beta, rho_w, rho_n, nA,
                                                       nB, nx, ny, nz, ux, uy,
                                                       uz, iwn);
                        } else {
                            // films that are close to the wetting fluid
                            if (morph_w->distance(i, j, k) < 2.5 && phi > 0.0) {
                                ifs.Mw += rho_w;
                                ifs.Pwx += rho_w * ux;
                                ifs.Pwy += rho_w * uy;
                                ifs.Pwz += rho_w * uz;
                            }
                            // films that are close to the NWP
                            if (morph_n->distance(i, j, k) < 2.5 && phi < 0.0) {
                                ifs.Mn += rho_n;
                                ifs.Pnx += rho_n * ux;
                                ifs.Pny += rho_n * uy;
                                ifs.Pnz += rho_n * uz;
                            }
                        }
                    } else if (phi > 0.0) {
                        if (morph_n->label(i, j, k) > 0) {
                            vol_nd_bulk += 1.0;
                            nd.p += Pressure(n);
                        } else {
                            vol_nc_bulk += 1.0;
                            nc.p += Pressure(n);
                        }
                    } else {
                        // water region
                        if (morph_w->label(i, j, k) > 0) {
                            vol_wd_bulk += 1.0;
                            wd.p += Pressure(n);
                        } else {
                            vol_wc_bulk += 1.0;
                            wc.p += Pressure(n);
                        }
                    }
                    if (phi > 0.0) {
                        if (morph_n->label(i, j, k) > 0) {
                            nA = 1.0;
                            nd.M += nA * rho_n;
                            nd.Px += nA * rho_n * ux;
                            nd.Py += nA * rho_n * uy;
                            nd.Pz += nA * rho_n * uz;
                            nd.K += nA * rho_n * (ux * ux + uy * uy + uz * uz);
                            nd.visc += visc;
                        } else {
                            nA = 1.0;
                            nc.M += nA * rho_n;
                            nc.Px += nA * rho_n * ux;
                            nc.Py += nA * rho_n * uy;
                            nc.Pz += nA * rho_n * uz;
                            nc.K += nA * rho_n * (ux * ux + uy * uy + uz * uz);
                            nc.visc += visc;
                        }
                    } else {
                        // water region
                        if (morph_w->label(i, j, k) > 0) {
                            nB = 1.0;
                            wd.M += nB * rho_w;
                            wd.Px += nB * rho_w * ux;
                            wd.Py += nB * rho_w * uy;
                            wd.Pz += nB * rho_w * uz;
                            wd.K += nB * rho_w * (ux * ux + uy * uy + uz * uz);
                            wd.visc += visc;
                        } else {
                            nB = 1.0;
                            wc.M += nB * rho_w;
                            wc.Px += nB * rho_w * ux;
                            wc.Py += nB * rho_w * uy;
                            wc.Pz += nB * rho_w * uz;
                            wc.K += nB * rho_w * (ux * ux + uy * uy + uz * uz);
                            wc.visc += visc;
                        }
                    }
                }
            }
        }
    }

    gnd.M = Dm->Comm.sumReduce(nd.M);
    gnd.Px = Dm->Comm.sumReduce(nd.Px);
    gnd.Py = Dm->Comm.sumReduce(nd.Py);
    gnd.Pz = Dm->Comm.sumReduce(nd.Pz);
    gnd.K = Dm->Comm.sumReduce(nd.K);
    gnd.visc = Dm->Comm.sumReduce(nd.visc);

    gwd.M = Dm->Comm.sumReduce(wd.M);
    gwd.Px = Dm->Comm.sumReduce(wd.Px);
    gwd.Py = Dm->Comm.sumReduce(wd.Py);
    gwd.Pz = Dm->Comm.sumReduce(wd.Pz);
    gwd.K = Dm->Comm.sumReduce(wd.K);
    gwd.visc = Dm->Comm.sumReduce(wd.visc);

    gnc.M = Dm->Comm.sumReduce(nc.M);
    gnc.Px = Dm->Comm.sumReduce(nc.Px);
    gnc.Py = Dm->Comm.sumReduce(nc.Py);
    gnc.Pz = Dm->Comm.sumReduce(nc.Pz);
    gnc.K = Dm->Comm.sumReduce(nc.K);
    gnc.visc = Dm->Comm.sumReduce(nc.visc);

    gwc.M = Dm->Comm.sumReduce(wc.M);
    gwc.Px = Dm->Comm.sumReduce(wc.Px);
    gwc.Py = Dm->Comm.sumReduce(wc.Py);
    gwc.Pz = Dm->Comm.sumReduce(wc.Pz);
    gwc.K = Dm->Comm.sumReduce(wc.K);
    gwc.visc = Dm->Comm.sumReduce(wc.visc);

    giwn.Mn = Dm->Comm.sumReduce(iwn.Mn);
    giwn.Pnx = Dm->Comm.sumReduce(iwn.Pnx);
    giwn.Pny = Dm->Comm.sumReduce(iwn.Pny);
    giwn.Pnz = Dm->Comm.sumReduce(iwn.Pnz);
    giwn.Kn = Dm->Comm.sumReduce(iwn.Kn);
    giwn.Mw = Dm->Comm.sumReduce(iwn.Mw);
    giwn.Pwx = Dm->Comm.sumReduce(iwn.Pwx);
    giwn.Pwy = Dm->Comm.sumReduce(iwn.Pwy);
    giwn.Pwz = Dm->Comm.sumReduce(iwn.Pwz);
    giwn.Kw = Dm->Comm.sumReduce(iwn.Kw);

    gifs.Mn = Dm->Comm.sumReduce(ifs.Mn);
    gifs.Pnx = Dm->Comm.sumReduce(ifs.Pnx);
    gifs.Pny = Dm->Comm.sumReduce(ifs.Pny);
    gifs.Pnz = Dm->Comm.sumReduce(ifs.Pnz);
    gifs.Mw = Dm->Comm.sumReduce(ifs.Mw);
    gifs.Pwx = Dm->Comm.sumReduce(ifs.Pwx);
    gifs.Pwy = Dm->Comm.sumReduce(ifs.Pwy);
    gifs.Pwz = Dm->Comm.sumReduce(ifs.Pwz);

    // pressure averaging
    gnc.p = Dm->Comm.sumReduce(nc.p);
    gnd.p = Dm->Comm.sumReduce(nd.p);
    gwc.p = Dm->Comm.sumReduce(wc.p);
    gwd.p = Dm->Comm.sumReduce(wd.p);

    if (vol_wc_bulk > 0.0)
        wc.p = wc.p / vol_wc_bulk;
    if (vol_nc_bulk > 0.0)
        nc.p = nc.p / vol_nc_bulk;
    if (vol_wd_bulk > 0.0)
        wd.p = wd.p / vol_wd_bulk;
    if (vol_nd_bulk > 0.0)
        nd.p = nd.p / vol_nd_bulk;

    vol_wc_bulk = Dm->Comm.sumReduce(vol_wc_bulk);
    vol_wd_bulk = Dm->Comm.sumReduce(vol_wd_bulk);
    vol_nc_bulk = Dm->Comm.sumReduce(vol_nc_bulk);
    vol_nd_bulk = Dm->Comm.sumReduce(vol_nd_bulk);

    if (vol_wc_bulk > 0.0)
        gwc.p = gwc.p / vol_wc_bulk;
    if (vol_nc_bulk > 0.0)
        gnc.p = gnc.p / vol_nc_bulk;
    if (vol_wd_bulk > 0.0)
        gwd.p = gwd.p / vol_wd_bulk;
    if (vol_nd_bulk > 0.0)
        gnd.p = gnd.p / vol_nd_bulk;
}

void SubPhase::AggregateLabels(const std::string &filename) {

    int nx = Dm->Nx;
    int ny = Dm->Ny;
    int nz = Dm->Nz;

    // assign the ID from the phase indicator field
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int n = k * nx * ny + j * nx + i;
                signed char local_id_val = Dm->id[n];
                if (local_id_val > 0) {
                    double value = Phi(i, j, k);
                    if (value > 0.0)
                        local_id_val = 1;
                    else
                        local_id_val = 2;
                }
                Dm->id[n] = local_id_val;
            }
        }
    }
    Dm->Comm.barrier();

    Dm->AggregateLabels(filename);
}
