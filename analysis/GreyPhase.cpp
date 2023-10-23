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
#include "analysis/GreyPhase.h"

// Constructor
GreyPhaseAnalysis::GreyPhaseAnalysis(std::shared_ptr<Domain> dm) : Dm(dm) {
    Nx = dm->Nx;
    Ny = dm->Ny;
    Nz = dm->Nz;
    Volume = (Nx - 2) * (Ny - 2) * (Nz - 2) * Dm->nprocx() * Dm->nprocy() *
             Dm->nprocz() * 1.0;

    // Global arrays
    SDs.resize(Nx, Ny, Nz);
    SDs.fill(0);
    Porosity.resize(Nx, Ny, Nz);
    Porosity.fill(0);
    //PhaseID.resize(Nx,Ny,Nz);       PhaseID.fill(0);
    Rho_n.resize(Nx, Ny, Nz);
    Rho_n.fill(0);
    Rho_w.resize(Nx, Ny, Nz);
    Rho_w.fill(0);
    Pressure.resize(Nx, Ny, Nz);
    Pressure.fill(0);
    //Phi.resize(Nx,Ny,Nz);         	Phi.fill(0);
    //DelPhi.resize(Nx,Ny,Nz);        DelPhi.fill(0);
    Vel_x.resize(Nx, Ny, Nz);
    Vel_x.fill(0); // Gradient of the phase indicator field
    Vel_y.resize(Nx, Ny, Nz);
    Vel_y.fill(0);
    Vel_z.resize(Nx, Ny, Nz);
    Vel_z.fill(0);
    MobilityRatio.resize(Nx, Ny, Nz);
    MobilityRatio.fill(0);
    //.........................................

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
            //fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
            fprintf(TIMELOG, "sw krw krn vw vn pw pn\n");
        }
    }
}

// Destructor
GreyPhaseAnalysis::~GreyPhaseAnalysis() {}

void GreyPhaseAnalysis::Write(int timestep) {}

void GreyPhaseAnalysis::SetParams(double rhoA, double rhoB, double tauA,
                                  double tauB, double force_x, double force_y,
                                  double force_z, double alpha, double B,
                                  double GreyPorosity) {
    Fx = force_x;
    Fy = force_y;
    Fz = force_z;
    rho_n = rhoA;
    rho_w = rhoB;
    nu_n = (tauA - 0.5) / 3.f;
    nu_w = (tauB - 0.5) / 3.f;
    gamma_wn = 6.0 * alpha;
    beta = B;
    grey_porosity = GreyPorosity;
}

void GreyPhaseAnalysis::Basic() {
    int i, j, k, n, imin, jmin, kmin, kmax;

    // If external boundary conditions are set, do not average over the inlet
    kmin = 1;
    kmax = Nz - 1;
    imin = jmin = 1;
    if (Dm->inlet_layers_z > 0 && Dm->kproc() == 0)
        kmin += Dm->inlet_layers_z;
    if (Dm->outlet_layers_z > 0 && Dm->kproc() == Dm->nprocz() - 1)
        kmax -= Dm->outlet_layers_z;

    Water_local.reset();
    Oil_local.reset();
    double count_w = 0.0;
    double count_n = 0.0;
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
                    double porosity = Porosity(n);
                    double mobility_ratio = MobilityRatio(n);

                    Water_local.M += nB * porosity;
                    Water_local.Px += porosity * (nA + nB) * Vel_x(n) * 0.5 *
                                      (1.0 - mobility_ratio);
                    Water_local.Py += porosity * (nA + nB) * Vel_y(n) * 0.5 *
                                      (1.0 - mobility_ratio);
                    Water_local.Pz += porosity * (nA + nB) * Vel_z(n) * 0.5 *
                                      (1.0 - mobility_ratio);

                    Oil_local.M += nA * porosity;
                    Oil_local.Px += porosity * (nA + nB) * Vel_x(n) * 0.5 *
                                    (1.0 + mobility_ratio);
                    Oil_local.Py += porosity * (nA + nB) * Vel_y(n) * 0.5 *
                                    (1.0 + mobility_ratio);
                    Oil_local.Pz += porosity * (nA + nB) * Vel_z(n) * 0.5 *
                                    (1.0 + mobility_ratio);

                    if (phi > 0.99) {
                        Oil_local.p += Pressure(n);
                        //Oil_local.p += pressure*(rho_n*nA)/(rho_n*nA+rho_w*nB);
                        count_n += 1.0;
                    } else if (phi < -0.99) {
                        Water_local.p += Pressure(n);
                        //Water_local.p += pressure*(rho_w*nB)/(rho_n*nA+rho_w*nB);
                        count_w += 1.0;
                    }
                }
            }
        }
    }
    Oil.M = Dm->Comm.sumReduce(Oil_local.M);
    Oil.Px = Dm->Comm.sumReduce(Oil_local.Px);
    Oil.Py = Dm->Comm.sumReduce(Oil_local.Py);
    Oil.Pz = Dm->Comm.sumReduce(Oil_local.Pz);

    Water.M = Dm->Comm.sumReduce(Water_local.M);
    Water.Px = Dm->Comm.sumReduce(Water_local.Px);
    Water.Py = Dm->Comm.sumReduce(Water_local.Py);
    Water.Pz = Dm->Comm.sumReduce(Water_local.Pz);

    //Oil.p /= Oil.M;
    //Water.p /= Water.M;
    count_w = Dm->Comm.sumReduce(count_w);
    count_n = Dm->Comm.sumReduce(count_n);
    if (count_w > 0.0)
        Water.p = Dm->Comm.sumReduce(Water_local.p) / count_w;
    else
        Water.p = 0.0;
    if (count_n > 0.0)
        Oil.p = Dm->Comm.sumReduce(Oil_local.p) / count_n;
    else
        Oil.p = 0.0;

    // check for NaN
    bool err = false;
    if (Water.M != Water.M)
        err = true;
    if (Water.p != Water.p)
        err = true;
    if (Water.Px != Water.Px)
        err = true;
    if (Water.Py != Water.Py)
        err = true;
    if (Water.Pz != Water.Pz)
        err = true;

    if (Oil.M != Oil.M)
        err = true;
    if (Oil.p != Oil.p)
        err = true;
    if (Oil.Px != Oil.Px)
        err = true;
    if (Oil.Py != Oil.Py)
        err = true;
    if (Oil.Pz != Oil.Pz)
        err = true;

    if (Dm->rank() == 0) {
        double force_mag = sqrt(Fx * Fx + Fy * Fy + Fz * Fz);
        double dir_x = 0.0;
        double dir_y = 0.0;
        double dir_z = 0.0;
        if (force_mag > 0.0) {
            dir_x = Fx / force_mag;
            dir_y = Fy / force_mag;
            dir_z = Fz / force_mag;
        } else {
            // default to z direction
            dir_x = 0.0;
            dir_y = 0.0;
            dir_z = 1.0;
        }
        if (Dm->BoundaryCondition == 1 || Dm->BoundaryCondition == 2 ||
            Dm->BoundaryCondition == 3 || Dm->BoundaryCondition == 4) {
            // compute the pressure drop
            double pressure_drop = (Pressure(Nx * Ny + Nx + 1) - 1.0) / 3.0;
            double length = ((Nz - 2) * Dm->nprocz());
            force_mag -= pressure_drop / length;
        }
        if (force_mag == 0.0) {
            // default to z direction
            dir_x = 0.0;
            dir_y = 0.0;
            dir_z = 1.0;
            force_mag = 1.0;
        }
        saturation = Water.M / (Water.M + Oil.M); // assume constant density
        water_flow_rate =
            grey_porosity * saturation *
            (Water.Px * dir_x + Water.Py * dir_y + Water.Pz * dir_z) / Water.M;
        oil_flow_rate = grey_porosity * (1.0 - saturation) *
                        (Oil.Px * dir_x + Oil.Py * dir_y + Oil.Pz * dir_z) /
                        Oil.M;

        double h = Dm->voxel_length;
        //TODO check if need greyporosity or domain porosity ? - compare to analytical solution
        double krn = h * h * nu_n * oil_flow_rate / force_mag;
        double krw = h * h * nu_w * water_flow_rate / force_mag;
        //printf("   water saturation = %f, fractional flow =%f \n",saturation,fractional_flow);
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g %.5g %.5g %.5g\n", saturation,
                krw, krn, h * water_flow_rate, h * oil_flow_rate, Water.p,
                Oil.p);
        fflush(TIMELOG);
    }

    if (err == true) {
        // exception if simulation produceds NaN
        printf("GreyPhaseAnalysis.cpp: NaN encountered, may need to check "
               "simulation parameters \n");
    }
    ASSERT(err == false);
}
/*
inline void InterfaceTransportMeasures( double beta, double rA, double rB, double nA, double nB, 
		double nx, double ny, double nz, double ux, double uy, double uz, interface &I){
	
	double A1,A2,A3,A4,A5,A6;
	double B1,B2,B3,B4,B5,B6;
	double nAB,delta;
	// Instantiate mass transport distributions
	// Stationary value - distribution 0
	nAB = 1.0/(nA+nB);
	//...............................................
	// q = 0,2,4
	// Cq = {1,0,0}, {0,1,0}, {0,0,1}
	delta = beta*nA*nB*nAB*0.1111111111111111*nx;
	if (!(nA*nB*nAB>0)) delta=0;
	A1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
	B1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
	A2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
	B2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;

	//...............................................
	// Cq = {0,1,0}
	delta = beta*nA*nB*nAB*0.1111111111111111*ny;
	if (!(nA*nB*nAB>0)) delta=0;
	A3 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
	B3 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
	A4 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
	B4 = nB*(0.1111111111111111*(1-4.5*uy))+delta;

	//...............................................
	// q = 4
	// Cq = {0,0,1}
	delta = beta*nA*nB*nAB*0.1111111111111111*nz;
	if (!(nA*nB*nAB>0)) delta=0;
	A5 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
	B5 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
	A6 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
	B6 = nB*(0.1111111111111111*(1-4.5*uz))+delta;
	
	double unx = (A1-A2);
	double uny = (A3-A4);
	double unz = (A5-A6);
	double uwx = (B1-B2);
	double uwy = (B3-B4);
	double uwz = (B5-B6);
	
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

}
*/
