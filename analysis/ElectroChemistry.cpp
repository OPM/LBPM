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
#include "analysis/ElectroChemistry.h"

ElectroChemistryAnalyzer::ElectroChemistryAnalyzer(std::shared_ptr<Domain> dm)
    : Dm(dm) {

    Nx = dm->Nx;
    Ny = dm->Ny;
    Nz = dm->Nz;
    Volume = (Nx - 2) * (Ny - 2) * (Nz - 2) * Dm->nprocx() * Dm->nprocy() *
             Dm->nprocz() * 1.0;

    ChemicalPotential.resize(Nx, Ny, Nz);
    ChemicalPotential.fill(0);
    ElectricalPotential.resize(Nx, Ny, Nz);
    ElectricalPotential.fill(0);
    ElectricalField_x.resize(Nx, Ny, Nz);
    ElectricalField_x.fill(0);
    ElectricalField_y.resize(Nx, Ny, Nz);
    ElectricalField_y.fill(0);
    ElectricalField_z.resize(Nx, Ny, Nz);
    ElectricalField_z.fill(0);
    Pressure.resize(Nx, Ny, Nz);
    Pressure.fill(0);
    Rho.resize(Nx, Ny, Nz);
    Rho.fill(0);
    Vel_x.resize(Nx, Ny, Nz);
    Vel_x.fill(0); // Gradient of the phase indicator field
    Vel_y.resize(Nx, Ny, Nz);
    Vel_y.fill(0);
    Vel_z.resize(Nx, Ny, Nz);
    Vel_z.fill(0);
    SDs.resize(Nx, Ny, Nz);
    SDs.fill(0);
    IonFluxDiffusive_x.resize(Nx, Ny, Nz);
    IonFluxDiffusive_x.fill(0);
    IonFluxDiffusive_y.resize(Nx, Ny, Nz);
    IonFluxDiffusive_y.fill(0);
    IonFluxDiffusive_z.resize(Nx, Ny, Nz);
    IonFluxDiffusive_z.fill(0);
    IonFluxAdvective_x.resize(Nx, Ny, Nz);
    IonFluxAdvective_x.fill(0);
    IonFluxAdvective_y.resize(Nx, Ny, Nz);
    IonFluxAdvective_y.fill(0);
    IonFluxAdvective_z.resize(Nx, Ny, Nz);
    IonFluxAdvective_z.fill(0);
    IonFluxElectrical_x.resize(Nx, Ny, Nz);
    IonFluxElectrical_x.fill(0);
    IonFluxElectrical_y.resize(Nx, Ny, Nz);
    IonFluxElectrical_y.fill(0);
    IonFluxElectrical_z.resize(Nx, Ny, Nz);
    IonFluxElectrical_z.fill(0);

    if (Dm->rank() == 0) {
        bool WriteHeader = false;
        TIMELOG = fopen("electrokinetic.csv", "r");
        if (TIMELOG != NULL)
            fclose(TIMELOG);
        else
            WriteHeader = true;

        TIMELOG = fopen("electrokinetic.csv", "a+");
        if (WriteHeader) {
            // If timelog is empty, write a short header to list the averages
            //fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
            fprintf(TIMELOG, "TBD TBD\n");
        }
    }
}

ElectroChemistryAnalyzer::ElectroChemistryAnalyzer(ScaLBL_IonModel &IonModel)
    : Dm(IonModel.Dm) {

    Nx = Dm->Nx;
    Ny = Dm->Ny;
    Nz = Dm->Nz;
    Volume = (Nx - 2) * (Ny - 2) * (Nz - 2) * Dm->nprocx() * Dm->nprocy() *
             Dm->nprocz() * 1.0;

    if (Dm->rank() == 0)
        printf("Analyze system with sub-domain size = %i x %i x %i \n", Nx, Ny,
               Nz);

    USE_MEMBRANE = IonModel.USE_MEMBRANE;

    ChemicalPotential.resize(Nx, Ny, Nz);
    ChemicalPotential.fill(0);
    ElectricalPotential.resize(Nx, Ny, Nz);
    ElectricalPotential.fill(0);
    ElectricalField_x.resize(Nx, Ny, Nz);
    ElectricalField_x.fill(0);
    ElectricalField_y.resize(Nx, Ny, Nz);
    ElectricalField_y.fill(0);
    ElectricalField_z.resize(Nx, Ny, Nz);
    ElectricalField_z.fill(0);
    Pressure.resize(Nx, Ny, Nz);
    Pressure.fill(0);
    Rho.resize(Nx, Ny, Nz);
    Rho.fill(0);
    Vel_x.resize(Nx, Ny, Nz);
    Vel_x.fill(0); // Gradient of the phase indicator field
    Vel_y.resize(Nx, Ny, Nz);
    Vel_y.fill(0);
    Vel_z.resize(Nx, Ny, Nz);
    Vel_z.fill(0);
    SDs.resize(Nx, Ny, Nz);
    SDs.fill(0);
    IonFluxDiffusive_x.resize(Nx, Ny, Nz);
    IonFluxDiffusive_x.fill(0);
    IonFluxDiffusive_y.resize(Nx, Ny, Nz);
    IonFluxDiffusive_y.fill(0);
    IonFluxDiffusive_z.resize(Nx, Ny, Nz);
    IonFluxDiffusive_z.fill(0);
    IonFluxAdvective_x.resize(Nx, Ny, Nz);
    IonFluxAdvective_x.fill(0);
    IonFluxAdvective_y.resize(Nx, Ny, Nz);
    IonFluxAdvective_y.fill(0);
    IonFluxAdvective_z.resize(Nx, Ny, Nz);
    IonFluxAdvective_z.fill(0);
    IonFluxElectrical_x.resize(Nx, Ny, Nz);
    IonFluxElectrical_x.fill(0);
    IonFluxElectrical_y.resize(Nx, Ny, Nz);
    IonFluxElectrical_y.fill(0);
    IonFluxElectrical_z.resize(Nx, Ny, Nz);
    IonFluxElectrical_z.fill(0);

    if (Dm->rank() == 0) {
        printf("Set up analysis routines for %lu ions \n",
               IonModel.number_ion_species);

        bool WriteHeader = false;
        TIMELOG = fopen("electrokinetic.csv", "r");
        if (TIMELOG != NULL)
            fclose(TIMELOG);
        else
            WriteHeader = true;

        TIMELOG = fopen("electrokinetic.csv", "a+");
        if (WriteHeader) {
            // If timelog is empty, write a short header to list the averages
            //fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
            fprintf(TIMELOG, "timestep voltage_out voltage_in ");
            fprintf(TIMELOG, "voltage_out_membrane voltage_in_membrane ");
            for (size_t i = 0; i < IonModel.number_ion_species; i++) {
                fprintf(TIMELOG, "rho_%lu_out rho_%lu_in ", i, i);
                fprintf(TIMELOG, "rho_%lu_out_membrane rho_%lu_in_membrane ", i,
                        i);
                fprintf(TIMELOG, "jx_%lu_out jx_%lu_in ", i, i);
                fprintf(TIMELOG, "jx_%lu_out_membrane jx_%lu_in_membrane ", i,
                        i);
                fprintf(TIMELOG, "jy_%lu_out jy_%lu_in ", i, i);
                fprintf(TIMELOG, "jy_%lu_out_membrane jy_%lu_in_membrane ", i,
                        i);
                fprintf(TIMELOG, "jz_%lu_out jz_%lu_in ", i, i);
                fprintf(TIMELOG, "jz_%lu_out_membrane jz_%lu_in_membrane ", i,
                        i);
            }
            fprintf(TIMELOG, "count_out count_in ");
            fprintf(TIMELOG, "count_out_membrane count_in_membrane\n");
        }
    }
}

ElectroChemistryAnalyzer::~ElectroChemistryAnalyzer() {
    if (Dm->rank() == 0) {
        fclose(TIMELOG);
    }
}

void ElectroChemistryAnalyzer::SetParams() {}

void ElectroChemistryAnalyzer::Membrane(ScaLBL_IonModel &Ion,
                                        ScaLBL_Poisson &Poisson, int timestep) {

    int i, j, k;

    Poisson.getElectricPotential(ElectricalPotential);

    if (Dm->rank() == 0)
        fprintf(TIMELOG, "%i ", timestep);

    /*   int iq, ip, nq, np, nqm, npm;
    Ion.MembraneDistance(i,j,k); // inside (-) or outside (+) the ion
    for (int link; link<Ion.IonMembrane->membraneLinkCount; link++){
        int iq = Ion.IonMembrane->membraneLinks[2*link];	
        int ip = Ion.IonMembrane->membraneLinks[2*link+1];	
        
		iq = membrane[2*link]; 	ip = membrane[2*link+1];
		nq = iq%Np;				np = ip%Np;
		nqm = Map[nq];			npm = Map[np];
    }
*/
    unsigned long int in_local_count, out_local_count;
    unsigned long int in_global_count, out_global_count;

    double value_in_local, value_out_local;
    double value_in_global, value_out_global;

    double value_membrane_in_local, value_membrane_out_local;
    double value_membrane_in_global, value_membrane_out_global;

    /* flux terms */
    double jx_membrane_in_local, jx_membrane_out_local;
    double jy_membrane_in_local, jy_membrane_out_local;
    double jz_membrane_in_local, jz_membrane_out_local;
    double jx_in_local, jx_out_local;
    double jy_in_local, jy_out_local;
    double jz_in_local, jz_out_local;

    double jx_membrane_in_global, jx_membrane_out_global;
    double jy_membrane_in_global, jy_membrane_out_global;
    double jz_membrane_in_global, jz_membrane_out_global;
    double jx_in_global, jx_out_global;
    double jy_in_global, jy_out_global;
    double jz_in_global, jz_out_global;

    unsigned long int membrane_in_local_count, membrane_out_local_count;
    unsigned long int membrane_in_global_count, membrane_out_global_count;

    double memdist, value, jx, jy, jz;
    in_local_count = 0;
    out_local_count = 0;
    membrane_in_local_count = 0;
    membrane_out_local_count = 0;

    value_membrane_in_local = 0.0;
    value_membrane_out_local = 0.0;
    value_in_local = 0.0;
    value_out_local = 0.0;
    for (k = Dm->inlet_layers_z; k < Nz; k++) {
        for (j = 1; j < Ny; j++) {
            for (i = 1; i < Nx; i++) {
                /* electric potential  */
                memdist = Ion.MembraneDistance(i, j, k);
                value = ElectricalPotential(i, j, k);
                if (memdist < 0.0) {
                    // inside the membrane
                    if (fabs(memdist) < 1.0) {
                        value_membrane_in_local += value;
                        membrane_in_local_count++;
                    }
                    value_in_local += value;
                    in_local_count++;

                } else {
                    // outside the membrane
                    if (fabs(memdist) < 1.0) {
                        value_membrane_out_local += value;
                        membrane_out_local_count++;
                    }
                    value_out_local += value;
                    out_local_count++;
                }
            }
        }
    }
    /* these only need to be computed the first time through */
    out_global_count = Dm->Comm.sumReduce(out_local_count);
    in_global_count = Dm->Comm.sumReduce(in_local_count);
    membrane_out_global_count = Dm->Comm.sumReduce(membrane_out_local_count);
    membrane_in_global_count = Dm->Comm.sumReduce(membrane_in_local_count);

    value_out_global = Dm->Comm.sumReduce(value_out_local);
    value_in_global = Dm->Comm.sumReduce(value_in_local);
    value_membrane_out_global = Dm->Comm.sumReduce(value_membrane_out_local);
    value_membrane_in_global = Dm->Comm.sumReduce(value_membrane_in_local);

    value_out_global /= out_global_count;
    value_in_global /= in_global_count;
    value_membrane_out_global /= membrane_out_global_count;
    value_membrane_in_global /= membrane_in_global_count;

    if (Dm->rank() == 0) {
        fprintf(TIMELOG, "%.8g ", value_out_global);
        fprintf(TIMELOG, "%.8g ", value_in_global);
        fprintf(TIMELOG, "%.8g ", value_membrane_out_global);
        fprintf(TIMELOG, "%.8g ", value_membrane_in_global);
    }

    value_membrane_in_local = 0.0;
    value_membrane_out_local = 0.0;
    value_in_local = 0.0;
    value_out_local = 0.0;
    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        Ion.getIonConcentration(Rho, ion);
        Ion.getIonFluxDiffusive(IonFluxDiffusive_x, IonFluxDiffusive_y,
                                IonFluxDiffusive_z, ion);
        Ion.getIonFluxAdvective(IonFluxAdvective_x, IonFluxAdvective_y,
                                IonFluxAdvective_z, ion);
        Ion.getIonFluxElectrical(IonFluxElectrical_x, IonFluxElectrical_y,
                                 IonFluxElectrical_z, ion);

        value_membrane_in_local = 0.0;
        value_membrane_out_local = 0.0;
        value_in_local = 0.0;
        value_out_local = 0.0;

        jx_membrane_in_local = jy_membrane_in_local = jz_membrane_in_local =
            0.0;
        jx_membrane_out_local = jy_membrane_out_local = jz_membrane_out_local =
            0.0;
        jx_in_local = jy_in_local = jz_in_local = 0.0;
        jx_out_local = jy_out_local = jz_out_local = 0.0;

        for (k = Dm->inlet_layers_z; k < Nz; k++) {
            for (j = 1; j < Ny; j++) {
                for (i = 1; i < Nx; i++) {
                    /* electric potential  */
                    memdist = Ion.MembraneDistance(i, j, k);
                    value = Rho(i, j, k);
                    jx = IonFluxDiffusive_x(i, j, k) +
                         IonFluxAdvective_x(i, j, k) +
                         IonFluxElectrical_x(i, j, k);
                    jy = IonFluxDiffusive_y(i, j, k) +
                         IonFluxAdvective_y(i, j, k) +
                         IonFluxElectrical_y(i, j, k);
                    jz = IonFluxDiffusive_z(i, j, k) +
                         IonFluxAdvective_z(i, j, k) +
                         IonFluxElectrical_z(i, j, k);

                    if (memdist < 0.0) {
                        // inside the membrane
                        if (fabs(memdist) < 1.0) {
                            value_membrane_in_local += value;
                            jx_membrane_in_local += jx;
                            jy_membrane_in_local += jy;
                            jz_membrane_in_local += jz;
                        }
                        value_in_local += value;
                        jx_in_local += jx;
                        jy_in_local += jy;
                        jz_in_local += jz;

                    } else {
                        // outside the membrane
                        if (fabs(memdist) < 1.0) {
                            value_membrane_out_local += value;
                            jx_membrane_out_local += jx;
                            jy_membrane_out_local += jy;
                            jz_membrane_out_local += jz;
                        }
                        value_out_local += value;
                        jx_out_local += jx;
                        jy_out_local += jy;
                        jz_out_local += jz;
                    }
                }
            }
        }
        value_out_global = Dm->Comm.sumReduce(value_out_local);
        value_in_global = Dm->Comm.sumReduce(value_in_local);
        value_membrane_out_global =
            Dm->Comm.sumReduce(value_membrane_out_local);
        value_membrane_in_global = Dm->Comm.sumReduce(value_membrane_in_local);

        value_out_global /= out_global_count;
        value_in_global /= in_global_count;
        value_membrane_out_global /= membrane_out_global_count;
        value_membrane_in_global /= membrane_in_global_count;

        jx_out_global = Dm->Comm.sumReduce(jx_out_local);
        jx_in_global = Dm->Comm.sumReduce(jx_in_local);
        jx_membrane_out_global = Dm->Comm.sumReduce(jx_membrane_out_local);
        jx_membrane_in_global = Dm->Comm.sumReduce(jx_membrane_in_local);

        jx_out_global /= out_global_count;
        jx_in_global /= in_global_count;
        jx_membrane_out_global /= membrane_out_global_count;
        jx_membrane_in_global /= membrane_in_global_count;

        jy_out_global = Dm->Comm.sumReduce(jy_out_local);
        jy_in_global = Dm->Comm.sumReduce(jy_in_local);
        jy_membrane_out_global = Dm->Comm.sumReduce(jy_membrane_out_local);
        jy_membrane_in_global = Dm->Comm.sumReduce(jy_membrane_in_local);

        jy_out_global /= out_global_count;
        jy_in_global /= in_global_count;
        jy_membrane_out_global /= membrane_out_global_count;
        jy_membrane_in_global /= membrane_in_global_count;

        jz_out_global = Dm->Comm.sumReduce(jz_out_local);
        jz_in_global = Dm->Comm.sumReduce(jz_in_local);
        jz_membrane_out_global = Dm->Comm.sumReduce(jz_membrane_out_local);
        jz_membrane_in_global = Dm->Comm.sumReduce(jz_membrane_in_local);

        jz_out_global /= out_global_count;
        jz_in_global /= in_global_count;
        jz_membrane_out_global /= membrane_out_global_count;
        jz_membrane_in_global /= membrane_in_global_count;

        if (Dm->rank() == 0) {
            fprintf(TIMELOG, "%.8g ", value_out_global);
            fprintf(TIMELOG, "%.8g ", value_in_global);
            fprintf(TIMELOG, "%.8g ", value_membrane_out_global);
            fprintf(TIMELOG, "%.8g ", value_membrane_in_global);

            fprintf(TIMELOG, "%.8g ", jx_out_global);
            fprintf(TIMELOG, "%.8g ", jx_in_global);
            fprintf(TIMELOG, "%.8g ", jx_membrane_out_global);
            fprintf(TIMELOG, "%.8g ", jx_membrane_in_global);

            fprintf(TIMELOG, "%.8g ", jy_out_global);
            fprintf(TIMELOG, "%.8g ", jy_in_global);
            fprintf(TIMELOG, "%.8g ", jy_membrane_out_global);
            fprintf(TIMELOG, "%.8g ", jy_membrane_in_global);

            fprintf(TIMELOG, "%.8g ", jz_out_global);
            fprintf(TIMELOG, "%.8g ", jz_in_global);
            fprintf(TIMELOG, "%.8g ", jz_membrane_out_global);
            fprintf(TIMELOG, "%.8g ", jz_membrane_in_global);
        }
    }

    if (Dm->rank() == 0) {
        fprintf(TIMELOG, "%lu ", out_global_count);
        fprintf(TIMELOG, "%lu ", in_global_count);
        fprintf(TIMELOG, "%lu ", membrane_out_global_count);
        fprintf(TIMELOG, "%lu\n", membrane_in_global_count);
        fflush(TIMELOG);
    }
}

void ElectroChemistryAnalyzer::Basic(ScaLBL_IonModel &Ion,
                                     ScaLBL_Poisson &Poisson,
                                     ScaLBL_StokesModel &Stokes, int timestep) {

    int i, j, k;
    double Vin = 0.0;
    double Vout = 0.0;
    Poisson.getElectricPotential(ElectricalPotential);

    /* local sub-domain averages */
    double *rho_avg_local;
    double *rho_mu_avg_local;
    double *rho_mu_fluctuation_local;
    double *rho_psi_avg_local;
    double *rho_psi_fluctuation_local;
    /* global averages */
    double *rho_avg_global;
    double *rho_mu_avg_global;
    double *rho_mu_fluctuation_global;
    double *rho_psi_avg_global;
    double *rho_psi_fluctuation_global;

    /* local sub-domain averages */
    rho_avg_local = new double[Ion.number_ion_species];
    rho_mu_avg_local = new double[Ion.number_ion_species];
    rho_mu_fluctuation_local = new double[Ion.number_ion_species];
    rho_psi_avg_local = new double[Ion.number_ion_species];
    rho_psi_fluctuation_local = new double[Ion.number_ion_species];
    /* global averages */
    rho_avg_global = new double[Ion.number_ion_species];
    rho_mu_avg_global = new double[Ion.number_ion_species];
    rho_mu_fluctuation_global = new double[Ion.number_ion_species];
    rho_psi_avg_global = new double[Ion.number_ion_species];
    rho_psi_fluctuation_global = new double[Ion.number_ion_species];

    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        rho_avg_local[ion] = 0.0;
        rho_mu_avg_local[ion] = 0.0;
        rho_psi_avg_local[ion] = 0.0;
        Ion.getIonConcentration(Rho, ion);
        /* Compute averages for each ion */
        for (k = 1; k < Nz; k++) {
            for (j = 1; j < Ny; j++) {
                for (i = 1; i < Nx; i++) {
                    rho_avg_local[ion] += Rho(i, j, k);
                    rho_mu_avg_local[ion] += Rho(i, j, k) * Rho(i, j, k);
                    rho_psi_avg_local[ion] +=
                        Rho(i, j, k) * ElectricalPotential(i, j, k);
                }
            }
        }
        rho_avg_global[ion] = Dm->Comm.sumReduce(rho_avg_local[ion]) / Volume;
        rho_mu_avg_global[ion] =
            Dm->Comm.sumReduce(rho_mu_avg_local[ion]) / Volume;
        rho_psi_avg_global[ion] =
            Dm->Comm.sumReduce(rho_psi_avg_local[ion]) / Volume;

        if (rho_avg_global[ion] > 0.0) {
            rho_mu_avg_global[ion] /= rho_avg_global[ion];
            rho_psi_avg_global[ion] /= rho_avg_global[ion];
        }
    }

    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        rho_mu_fluctuation_local[ion] = 0.0;
        rho_psi_fluctuation_local[ion] = 0.0;
        /* Compute averages for each ion */
        for (k = 1; k < Nz; k++) {
            for (j = 1; j < Ny; j++) {
                for (i = 1; i < Nx; i++) {
                    rho_mu_fluctuation_local[ion] +=
                        (Rho(i, j, k) * Rho(i, j, k) - rho_mu_avg_global[ion]);
                    rho_psi_fluctuation_local[ion] +=
                        (Rho(i, j, k) * ElectricalPotential(i, j, k) -
                         rho_psi_avg_global[ion]);
                }
            }
        }
        rho_mu_fluctuation_global[ion] =
            Dm->Comm.sumReduce(rho_mu_fluctuation_local[ion]);
        rho_psi_fluctuation_global[ion] =
            Dm->Comm.sumReduce(rho_psi_fluctuation_local[ion]);
    }

    if (Dm->rank() == 0) {
        fprintf(TIMELOG, "%i ", timestep);
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            fprintf(TIMELOG, "%.8g ", rho_avg_global[ion]);
            fprintf(TIMELOG, "%.8g ", rho_mu_avg_global[ion]);
            fprintf(TIMELOG, "%.8g ", rho_psi_avg_global[ion]);
            fprintf(TIMELOG, "%.8g ", rho_mu_fluctuation_global[ion]);
            fprintf(TIMELOG, "%.8g ", rho_psi_fluctuation_global[ion]);
        }
        fprintf(TIMELOG, "%.8g %.8g\n", Vin, Vout);
        fflush(TIMELOG);
    }
    /*	else{
		fprintf(TIMELOG,"%i ",timestep); 
		for (int ion=0; ion<Ion.number_ion_species; ion++){
			fprintf(TIMELOG,"%.8g ",rho_avg_local[ion]);
			fprintf(TIMELOG,"%.8g ",rho_mu_avg_local[ion]);
			fprintf(TIMELOG,"%.8g ",rho_psi_avg_local[ion]);
			fprintf(TIMELOG,"%.8g ",rho_mu_fluctuation_local[ion]);
			fprintf(TIMELOG,"%.8g ",rho_psi_fluctuation_local[ion]);
		}
		fflush(TIMELOG);
	} */
}

void ElectroChemistryAnalyzer::WriteVis(ScaLBL_IonModel &Ion,
                                        ScaLBL_Poisson &Poisson,
                                        ScaLBL_StokesModel &Stokes,
                                        std::shared_ptr<Database> input_db,
                                        int timestep) {

    auto vis_db = input_db->getDatabase("Visualization");
    char VisName[40];
    auto format = vis_db->getWithDefault<string>("format", "hdf5");

    if (Dm->rank() == 0) {
        printf("ElectroChemistryAnalyzer::WriteVis (format = %s)\n",
               format.c_str());
        if (vis_db->getWithDefault<bool>("save_electric_potential", true)) {
            printf("       save electric potential \n");
        }
        if (vis_db->getWithDefault<bool>("save_concentration", true)) {
            printf("       save concentration \n");
        }
        if (vis_db->getWithDefault<bool>("save_velocity", false)) {
            printf("       save velocity \n");
        }
        if (vis_db->getWithDefault<bool>("save_ion_flux_diffusive", false)) {
            printf("       save ion flux (diffusive) \n");
        }
        if (vis_db->getWithDefault<bool>("save_ion_flux_advective", false)) {
            printf("       save ion flux (advective) \n");
        }
        if (vis_db->getWithDefault<bool>("save_ion_flux_electrical", false)) {
            printf("       save ion flux (electrical) \n");
        }
        if (vis_db->getWithDefault<bool>("save_electric_field", false)) {
            printf("       save electric field \n");
        }
    }

    std::vector<IO::MeshDataStruct> visData;
    fillHalo<double> fillData(Dm->Comm, Dm->rank_info,
                              {Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2}, {1, 1, 1},
                              0, 1);

    IO::initialize("", format, "false");
    // Create the MeshDataStruct
    visData.resize(1);

    visData[0].meshName = "domain";
    visData[0].mesh =
        std::make_shared<IO::DomainMesh>(Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2,
                                         Dm->Nz - 2, Dm->Lx, Dm->Ly, Dm->Lz);
    //electric potential
    auto ElectricPotentialVar = std::make_shared<IO::Variable>();
    //electric field
    auto ElectricFieldVar_x = std::make_shared<IO::Variable>();
    auto ElectricFieldVar_y = std::make_shared<IO::Variable>();
    auto ElectricFieldVar_z = std::make_shared<IO::Variable>();

    //ion concentration
    std::vector<shared_ptr<IO::Variable>> IonConcentration;
    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        IonConcentration.push_back(std::make_shared<IO::Variable>());
    }
    //fluid velocity
    auto VxVar = std::make_shared<IO::Variable>();
    auto VyVar = std::make_shared<IO::Variable>();
    auto VzVar = std::make_shared<IO::Variable>();
    // diffusive ion flux
    std::vector<shared_ptr<IO::Variable>> IonFluxDiffusive;
    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        //push in x-,y-, and z-component for each ion species
        IonFluxDiffusive.push_back(std::make_shared<IO::Variable>());
        IonFluxDiffusive.push_back(std::make_shared<IO::Variable>());
        IonFluxDiffusive.push_back(std::make_shared<IO::Variable>());
    }
    // advective ion flux
    std::vector<shared_ptr<IO::Variable>> IonFluxAdvective;
    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        //push in x-,y-, and z-component for each ion species
        IonFluxAdvective.push_back(std::make_shared<IO::Variable>());
        IonFluxAdvective.push_back(std::make_shared<IO::Variable>());
        IonFluxAdvective.push_back(std::make_shared<IO::Variable>());
    }
    // electro-migrational ion flux
    std::vector<shared_ptr<IO::Variable>> IonFluxElectrical;
    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        //push in x-,y-, and z-component for each ion species
        IonFluxElectrical.push_back(std::make_shared<IO::Variable>());
        IonFluxElectrical.push_back(std::make_shared<IO::Variable>());
        IonFluxElectrical.push_back(std::make_shared<IO::Variable>());
    }
    //--------------------------------------------------------------------------------------------------------------------

    //-------------------------------------Create Names for Variables------------------------------------------------------
    if (vis_db->getWithDefault<bool>("save_electric_potential", true)) {
        ElectricPotentialVar->name = "ElectricPotential";
        ElectricPotentialVar->type = IO::VariableType::VolumeVariable;
        ElectricPotentialVar->dim = 1;
        ElectricPotentialVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(ElectricPotentialVar);
    }

    if (vis_db->getWithDefault<bool>("save_concentration", true)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            sprintf(VisName, "IonConcentration_%zu", ion + 1);
            IonConcentration[ion]->name = VisName;
            IonConcentration[ion]->type = IO::VariableType::VolumeVariable;
            IonConcentration[ion]->dim = 1;
            IonConcentration[ion]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                               Dm->Nz - 2);
            visData[0].vars.push_back(IonConcentration[ion]);
        }
    }

    if (vis_db->getWithDefault<bool>("save_velocity", false)) {
        VxVar->name = "Velocity_x";
        VxVar->type = IO::VariableType::VolumeVariable;
        VxVar->dim = 1;
        VxVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VxVar);
        VyVar->name = "Velocity_y";
        VyVar->type = IO::VariableType::VolumeVariable;
        VyVar->dim = 1;
        VyVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VyVar);
        VzVar->name = "Velocity_z";
        VzVar->type = IO::VariableType::VolumeVariable;
        VzVar->dim = 1;
        VzVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VzVar);
    }

    if (vis_db->getWithDefault<bool>("save_ion_flux_diffusive", false)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            // x-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_x", ion + 1);
            IonFluxDiffusive[3 * ion + 0]->name = VisName;
            IonFluxDiffusive[3 * ion + 0]->type =
                IO::VariableType::VolumeVariable;
            IonFluxDiffusive[3 * ion + 0]->dim = 1;
            IonFluxDiffusive[3 * ion + 0]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                       Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxDiffusive[3 * ion + 0]);
            // y-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_y", ion + 1);
            IonFluxDiffusive[3 * ion + 1]->name = VisName;
            IonFluxDiffusive[3 * ion + 1]->type =
                IO::VariableType::VolumeVariable;
            IonFluxDiffusive[3 * ion + 1]->dim = 1;
            IonFluxDiffusive[3 * ion + 1]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                       Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxDiffusive[3 * ion + 1]);
            // z-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_z", ion + 1);
            IonFluxDiffusive[3 * ion + 2]->name = VisName;
            IonFluxDiffusive[3 * ion + 2]->type =
                IO::VariableType::VolumeVariable;
            IonFluxDiffusive[3 * ion + 2]->dim = 1;
            IonFluxDiffusive[3 * ion + 2]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                       Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxDiffusive[3 * ion + 2]);
        }
    }

    if (vis_db->getWithDefault<bool>("save_ion_flux_advective", false)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            // x-component of advective flux
            sprintf(VisName, "Ion%zu_FluxAdvective_x", ion + 1);
            IonFluxAdvective[3 * ion + 0]->name = VisName;
            IonFluxAdvective[3 * ion + 0]->type =
                IO::VariableType::VolumeVariable;
            IonFluxAdvective[3 * ion + 0]->dim = 1;
            IonFluxAdvective[3 * ion + 0]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                       Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxAdvective[3 * ion + 0]);
            // y-component of advective flux
            sprintf(VisName, "Ion%zu_FluxAdvective_y", ion + 1);
            IonFluxAdvective[3 * ion + 1]->name = VisName;
            IonFluxAdvective[3 * ion + 1]->type =
                IO::VariableType::VolumeVariable;
            IonFluxAdvective[3 * ion + 1]->dim = 1;
            IonFluxAdvective[3 * ion + 1]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                       Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxAdvective[3 * ion + 1]);
            // z-component of advective flux
            sprintf(VisName, "Ion%zu_FluxAdvective_z", ion + 1);
            IonFluxAdvective[3 * ion + 2]->name = VisName;
            IonFluxAdvective[3 * ion + 2]->type =
                IO::VariableType::VolumeVariable;
            IonFluxAdvective[3 * ion + 2]->dim = 1;
            IonFluxAdvective[3 * ion + 2]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                       Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxAdvective[3 * ion + 2]);
        }
    }

    if (vis_db->getWithDefault<bool>("save_ion_flux_electrical", false)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            // x-component of electro-migrational flux
            sprintf(VisName, "Ion%zu_FluxElectrical_x", ion + 1);
            IonFluxElectrical[3 * ion + 0]->name = VisName;
            IonFluxElectrical[3 * ion + 0]->type =
                IO::VariableType::VolumeVariable;
            IonFluxElectrical[3 * ion + 0]->dim = 1;
            IonFluxElectrical[3 * ion + 0]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                        Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxElectrical[3 * ion + 0]);
            // y-component of electro-migrational flux
            sprintf(VisName, "Ion%zu_FluxElectrical_y", ion + 1);
            IonFluxElectrical[3 * ion + 1]->name = VisName;
            IonFluxElectrical[3 * ion + 1]->type =
                IO::VariableType::VolumeVariable;
            IonFluxElectrical[3 * ion + 1]->dim = 1;
            IonFluxElectrical[3 * ion + 1]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                        Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxElectrical[3 * ion + 1]);
            // z-component of electro-migrational flux
            sprintf(VisName, "Ion%zu_FluxElectrical_z", ion + 1);
            IonFluxElectrical[3 * ion + 2]->name = VisName;
            IonFluxElectrical[3 * ion + 2]->type =
                IO::VariableType::VolumeVariable;
            IonFluxElectrical[3 * ion + 2]->dim = 1;
            IonFluxElectrical[3 * ion + 2]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                        Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxElectrical[3 * ion + 2]);
        }
    }

    if (vis_db->getWithDefault<bool>("save_electric_field", false)) {
        ElectricFieldVar_x->name = "ElectricField_x";
        ElectricFieldVar_x->type = IO::VariableType::VolumeVariable;
        ElectricFieldVar_x->dim = 1;
        ElectricFieldVar_x->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(ElectricFieldVar_x);
        ElectricFieldVar_y->name = "ElectricField_y";
        ElectricFieldVar_y->type = IO::VariableType::VolumeVariable;
        ElectricFieldVar_y->dim = 1;
        ElectricFieldVar_y->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(ElectricFieldVar_y);
        ElectricFieldVar_z->name = "ElectricField_z";
        ElectricFieldVar_z->type = IO::VariableType::VolumeVariable;
        ElectricFieldVar_z->dim = 1;
        ElectricFieldVar_z->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(ElectricFieldVar_z);
    }
    //--------------------------------------------------------------------------------------------------------------------

    //------------------------------------Save All Variables--------------------------------------------------------------
    if (vis_db->getWithDefault<bool>("save_electric_potential", true)) {
        ASSERT(visData[0].vars[0]->name == "ElectricPotential");
        Poisson.getElectricPotential(ElectricalPotential);
        Array<double> &ElectricPotentialData = visData[0].vars[0]->data;
        fillData.copy(ElectricalPotential, ElectricPotentialData);
    }

    if (vis_db->getWithDefault<bool>("save_concentration", true)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            sprintf(VisName, "IonConcentration_%zu", ion + 1);
            //IonConcentration[ion]->name = VisName;
            ASSERT(visData[0].vars[1 + ion]->name == VisName);
            Array<double> &IonConcentrationData =
                visData[0].vars[1 + ion]->data;
            Ion.getIonConcentration(Rho, ion);
            fillData.copy(Rho, IonConcentrationData);
        }
    }

    if (vis_db->getWithDefault<bool>("save_velocity", false)) {
        ASSERT(visData[0].vars[1 + Ion.number_ion_species + 0]->name ==
               "Velocity_x");
        ASSERT(visData[0].vars[1 + Ion.number_ion_species + 1]->name ==
               "Velocity_y");
        ASSERT(visData[0].vars[1 + Ion.number_ion_species + 2]->name ==
               "Velocity_z");
        Stokes.getVelocity(Vel_x, Vel_y, Vel_z);
        Array<double> &VelxData =
            visData[0].vars[1 + Ion.number_ion_species + 0]->data;
        Array<double> &VelyData =
            visData[0].vars[1 + Ion.number_ion_species + 1]->data;
        Array<double> &VelzData =
            visData[0].vars[1 + Ion.number_ion_species + 2]->data;
        fillData.copy(Vel_x, VelxData);
        fillData.copy(Vel_y, VelyData);
        fillData.copy(Vel_z, VelzData);
    }

    if (vis_db->getWithDefault<bool>("save_ion_flux_diffusive", false)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {

            // x-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_x", ion + 1);
            //IonFluxDiffusive[3*ion+0]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species + 3 * ion + 0]
                       ->name == VisName);
            // y-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_y", ion + 1);
            //IonFluxDiffusive[3*ion+1]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species + 3 * ion + 1]
                       ->name == VisName);
            // z-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_z", ion + 1);
            //IonFluxDiffusive[3*ion+2]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species + 3 * ion + 2]
                       ->name == VisName);

            Array<double> &IonFluxData_x =
                visData[0].vars[4 + Ion.number_ion_species + 3 * ion + 0]->data;
            Array<double> &IonFluxData_y =
                visData[0].vars[4 + Ion.number_ion_species + 3 * ion + 1]->data;
            Array<double> &IonFluxData_z =
                visData[0].vars[4 + Ion.number_ion_species + 3 * ion + 2]->data;
            Ion.getIonFluxDiffusive(IonFluxDiffusive_x, IonFluxDiffusive_y,
                                    IonFluxDiffusive_z, ion);
            fillData.copy(IonFluxDiffusive_x, IonFluxData_x);
            fillData.copy(IonFluxDiffusive_y, IonFluxData_y);
            fillData.copy(IonFluxDiffusive_z, IonFluxData_z);
        }
    }

    if (vis_db->getWithDefault<bool>("save_ion_flux_advective", false)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {

            // x-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxAdvective_x", ion + 1);
            //IonFluxDiffusive[3*ion+0]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species * (1 + 3) + 3 * ion + 0]
                       ->name == VisName);
            // y-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxAdvective_y", ion + 1);
            //IonFluxDiffusive[3*ion+1]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species * (1 + 3) + 3 * ion + 1]
                       ->name == VisName);
            // z-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxAdvective_z", ion + 1);
            //IonFluxDiffusive[3*ion+2]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species * (1 + 3) + 3 * ion + 2]
                       ->name == VisName);

            Array<double> &IonFluxData_x =
                visData[0]
                    .vars[4 + Ion.number_ion_species * (1 + 3) + 3 * ion + 0]
                    ->data;
            Array<double> &IonFluxData_y =
                visData[0]
                    .vars[4 + Ion.number_ion_species * (1 + 3) + 3 * ion + 1]
                    ->data;
            Array<double> &IonFluxData_z =
                visData[0]
                    .vars[4 + Ion.number_ion_species * (1 + 3) + 3 * ion + 2]
                    ->data;
            Ion.getIonFluxAdvective(IonFluxAdvective_x, IonFluxAdvective_y,
                                    IonFluxAdvective_z, ion);
            fillData.copy(IonFluxAdvective_x, IonFluxData_x);
            fillData.copy(IonFluxAdvective_y, IonFluxData_y);
            fillData.copy(IonFluxAdvective_z, IonFluxData_z);
        }
    }

    if (vis_db->getWithDefault<bool>("save_ion_flux_electrical", false)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {

            // x-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxElectrical_x", ion + 1);
            //IonFluxDiffusive[3*ion+0]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 0]
                       ->name == VisName);
            // y-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxElectrical_y", ion + 1);
            //IonFluxDiffusive[3*ion+1]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 1]
                       ->name == VisName);
            // z-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxElectrical_z", ion + 1);
            //IonFluxDiffusive[3*ion+2]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 2]
                       ->name == VisName);

            Array<double> &IonFluxData_x =
                visData[0]
                    .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 0]
                    ->data;
            Array<double> &IonFluxData_y =
                visData[0]
                    .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 1]
                    ->data;
            Array<double> &IonFluxData_z =
                visData[0]
                    .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 2]
                    ->data;
            Ion.getIonFluxElectrical(IonFluxElectrical_x, IonFluxElectrical_y,
                                     IonFluxElectrical_z, ion);
            fillData.copy(IonFluxElectrical_x, IonFluxData_x);
            fillData.copy(IonFluxElectrical_y, IonFluxData_y);
            fillData.copy(IonFluxElectrical_z, IonFluxData_z);
        }
    }

    if (vis_db->getWithDefault<bool>("save_electric_field", false)) {
        ASSERT(
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 0]->name ==
            "ElectricField_x");
        ASSERT(
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 1]->name ==
            "ElectricField_y");
        ASSERT(
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 2]->name ==
            "ElectricField_z");
        Poisson.getElectricField(ElectricalField_x, ElectricalField_y,
                                 ElectricalField_z);
        Array<double> &ElectricalFieldxData =
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 0]->data;
        Array<double> &ElectricalFieldyData =
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 1]->data;
        Array<double> &ElectricalFieldzData =
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 2]->data;
        fillData.copy(ElectricalField_x, ElectricalFieldxData);
        fillData.copy(ElectricalField_y, ElectricalFieldyData);
        fillData.copy(ElectricalField_z, ElectricalFieldzData);
    }

    if (vis_db->getWithDefault<bool>("write_silo", true))
        IO::writeData(timestep, visData, Dm->Comm);
    //--------------------------------------------------------------------------------------------------------------------
    /*    if (vis_db->getWithDefault<bool>( "save_8bit_raw", true )){
    	char CurrentIDFilename[40];
    	sprintf(CurrentIDFilename,"id_t%d.raw",timestep);
    	Averages.AggregateLabels(CurrentIDFilename);
    }
*/
}

void ElectroChemistryAnalyzer::Basic(ScaLBL_IonModel &Ion,
                                     ScaLBL_Poisson &Poisson, int timestep) {

    int i, j, k;
    double Vin = 0.0;
    double Vout = 0.0;
    Poisson.getElectricPotential(ElectricalPotential);

    /* local sub-domain averages */
    double *rho_avg_local;
    double *rho_mu_avg_local;
    double *rho_mu_fluctuation_local;
    double *rho_psi_avg_local;
    double *rho_psi_fluctuation_local;
    /* global averages */
    double *rho_avg_global;
    double *rho_mu_avg_global;
    double *rho_mu_fluctuation_global;
    double *rho_psi_avg_global;
    double *rho_psi_fluctuation_global;

    /* Get the distance to the membrane */
    if (Ion.USE_MEMBRANE) {
        //Ion.MembraneDistance;
    }

    /* local sub-domain averages */
    rho_avg_local = new double[Ion.number_ion_species];
    rho_mu_avg_local = new double[Ion.number_ion_species];
    rho_mu_fluctuation_local = new double[Ion.number_ion_species];
    rho_psi_avg_local = new double[Ion.number_ion_species];
    rho_psi_fluctuation_local = new double[Ion.number_ion_species];
    /* global averages */
    rho_avg_global = new double[Ion.number_ion_species];
    rho_mu_avg_global = new double[Ion.number_ion_species];
    rho_mu_fluctuation_global = new double[Ion.number_ion_species];
    rho_psi_avg_global = new double[Ion.number_ion_species];
    rho_psi_fluctuation_global = new double[Ion.number_ion_species];

    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        rho_avg_local[ion] = 0.0;
        rho_mu_avg_local[ion] = 0.0;
        rho_psi_avg_local[ion] = 0.0;
        Ion.getIonConcentration(Rho, ion);
        /* Compute averages for each ion */
        for (k = 1; k < Nz; k++) {
            for (j = 1; j < Ny; j++) {
                for (i = 1; i < Nx; i++) {
                    rho_avg_local[ion] += Rho(i, j, k);
                    rho_mu_avg_local[ion] += Rho(i, j, k) * Rho(i, j, k);
                    rho_psi_avg_local[ion] +=
                        Rho(i, j, k) * ElectricalPotential(i, j, k);
                }
            }
        }
        rho_avg_global[ion] = Dm->Comm.sumReduce(rho_avg_local[ion]) / Volume;
        rho_mu_avg_global[ion] =
            Dm->Comm.sumReduce(rho_mu_avg_local[ion]) / Volume;
        rho_psi_avg_global[ion] =
            Dm->Comm.sumReduce(rho_psi_avg_local[ion]) / Volume;

        if (rho_avg_global[ion] > 0.0) {
            rho_mu_avg_global[ion] /= rho_avg_global[ion];
            rho_psi_avg_global[ion] /= rho_avg_global[ion];
        }
    }

    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        rho_mu_fluctuation_local[ion] = 0.0;
        rho_psi_fluctuation_local[ion] = 0.0;
        /* Compute averages for each ion */
        for (k = 1; k < Nz; k++) {
            for (j = 1; j < Ny; j++) {
                for (i = 1; i < Nx; i++) {
                    rho_mu_fluctuation_local[ion] +=
                        (Rho(i, j, k) * Rho(i, j, k) - rho_mu_avg_global[ion]);
                    rho_psi_fluctuation_local[ion] +=
                        (Rho(i, j, k) * ElectricalPotential(i, j, k) -
                         rho_psi_avg_global[ion]);
                }
            }
        }
        rho_mu_fluctuation_global[ion] =
            Dm->Comm.sumReduce(rho_mu_fluctuation_local[ion]);
        rho_psi_fluctuation_global[ion] =
            Dm->Comm.sumReduce(rho_psi_fluctuation_local[ion]);
    }

    if (Dm->rank() == 0) {
        fprintf(TIMELOG, "%i ", timestep);
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            fprintf(TIMELOG, "%.8g ", rho_avg_global[ion]);
            fprintf(TIMELOG, "%.8g ", rho_mu_avg_global[ion]);
            fprintf(TIMELOG, "%.8g ", rho_psi_avg_global[ion]);
            fprintf(TIMELOG, "%.8g ", rho_mu_fluctuation_global[ion]);
            fprintf(TIMELOG, "%.8g ", rho_psi_fluctuation_global[ion]);
        }
        fprintf(TIMELOG, "%.8g %.8g\n", Vin, Vout);
        fflush(TIMELOG);
    }
    /*	else{
		fprintf(TIMELOG,"%i ",timestep); 
		for (int ion=0; ion<Ion.number_ion_species; ion++){
			fprintf(TIMELOG,"%.8g ",rho_avg_local[ion]);
			fprintf(TIMELOG,"%.8g ",rho_mu_avg_local[ion]);
			fprintf(TIMELOG,"%.8g ",rho_psi_avg_local[ion]);
			fprintf(TIMELOG,"%.8g ",rho_mu_fluctuation_local[ion]);
			fprintf(TIMELOG,"%.8g ",rho_psi_fluctuation_local[ion]);
		}
		fflush(TIMELOG);
	} */
}

void ElectroChemistryAnalyzer::WriteVis(ScaLBL_IonModel &Ion,
                                        ScaLBL_Poisson &Poisson,
                                        std::shared_ptr<Database> input_db,
                                        int timestep) {

    auto vis_db = input_db->getDatabase("Visualization");
    char VisName[40];
    auto format = vis_db->getWithDefault<string>("format", "hdf5");

    std::vector<IO::MeshDataStruct> visData;
    fillHalo<double> fillData(Dm->Comm, Dm->rank_info,
                              {Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2}, {1, 1, 1},
                              0, 1);

    if (Dm->rank() == 0) {
        printf("ElectroChemistryAnalyzer::WriteVis (format = %s)\n",
               format.c_str());
        if (vis_db->getWithDefault<bool>("save_electric_potential", true)) {
            printf("       save electric potential \n");
        }
        if (vis_db->getWithDefault<bool>("save_concentration", true)) {
            printf("       save concentration \n");
        }
        if (vis_db->getWithDefault<bool>("save_velocity", false)) {
            printf("       save velocity \n");
        }
        if (vis_db->getWithDefault<bool>("save_ion_flux_diffusive", false)) {
            printf("       save ion flux (diffusive) \n");
        }
        if (vis_db->getWithDefault<bool>("save_ion_flux_advective", false)) {
            printf("       save ion flux (advective) \n");
        }
        if (vis_db->getWithDefault<bool>("save_ion_flux_electrical", false)) {
            printf("       save ion flux (electrical) \n");
        }
        if (vis_db->getWithDefault<bool>("save_electric_field", false)) {
            printf("       save electric field \n");
        }
    }

    IO::initialize("", format, "false");
    // Create the MeshDataStruct
    visData.resize(1);

    visData[0].meshName = "domain";
    visData[0].mesh =
        std::make_shared<IO::DomainMesh>(Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2,
                                         Dm->Nz - 2, Dm->Lx, Dm->Ly, Dm->Lz);
    //electric potential
    auto ElectricPotentialVar = std::make_shared<IO::Variable>();
    //electric field
    auto ElectricFieldVar_x = std::make_shared<IO::Variable>();
    auto ElectricFieldVar_y = std::make_shared<IO::Variable>();
    auto ElectricFieldVar_z = std::make_shared<IO::Variable>();

    //ion concentration
    std::vector<shared_ptr<IO::Variable>> IonConcentration;
    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        IonConcentration.push_back(std::make_shared<IO::Variable>());
    }

    // diffusive ion flux
    std::vector<shared_ptr<IO::Variable>> IonFluxDiffusive;
    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        //push in x-,y-, and z-component for each ion species
        IonFluxDiffusive.push_back(std::make_shared<IO::Variable>());
        IonFluxDiffusive.push_back(std::make_shared<IO::Variable>());
        IonFluxDiffusive.push_back(std::make_shared<IO::Variable>());
    }

    // electro-migrational ion flux
    std::vector<shared_ptr<IO::Variable>> IonFluxElectrical;
    for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
        //push in x-,y-, and z-component for each ion species
        IonFluxElectrical.push_back(std::make_shared<IO::Variable>());
        IonFluxElectrical.push_back(std::make_shared<IO::Variable>());
        IonFluxElectrical.push_back(std::make_shared<IO::Variable>());
    }
    //--------------------------------------------------------------------------------------------------------------------

    //-------------------------------------Create Names for Variables------------------------------------------------------
    if (vis_db->getWithDefault<bool>("save_electric_potential", true)) {
        ElectricPotentialVar->name = "ElectricPotential";
        ElectricPotentialVar->type = IO::VariableType::VolumeVariable;
        ElectricPotentialVar->dim = 1;
        ElectricPotentialVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(ElectricPotentialVar);
    }

    if (vis_db->getWithDefault<bool>("save_concentration", true)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            sprintf(VisName, "IonConcentration_%zu", ion + 1);
            IonConcentration[ion]->name = VisName;
            IonConcentration[ion]->type = IO::VariableType::VolumeVariable;
            IonConcentration[ion]->dim = 1;
            IonConcentration[ion]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                               Dm->Nz - 2);
            visData[0].vars.push_back(IonConcentration[ion]);
        }
    }

    if (vis_db->getWithDefault<bool>("save_ion_flux_diffusive", false)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            // x-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_x", ion + 1);
            IonFluxDiffusive[3 * ion + 0]->name = VisName;
            IonFluxDiffusive[3 * ion + 0]->type =
                IO::VariableType::VolumeVariable;
            IonFluxDiffusive[3 * ion + 0]->dim = 1;
            IonFluxDiffusive[3 * ion + 0]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                       Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxDiffusive[3 * ion + 0]);
            // y-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_y", ion + 1);
            IonFluxDiffusive[3 * ion + 1]->name = VisName;
            IonFluxDiffusive[3 * ion + 1]->type =
                IO::VariableType::VolumeVariable;
            IonFluxDiffusive[3 * ion + 1]->dim = 1;
            IonFluxDiffusive[3 * ion + 1]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                       Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxDiffusive[3 * ion + 1]);
            // z-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_z", ion + 1);
            IonFluxDiffusive[3 * ion + 2]->name = VisName;
            IonFluxDiffusive[3 * ion + 2]->type =
                IO::VariableType::VolumeVariable;
            IonFluxDiffusive[3 * ion + 2]->dim = 1;
            IonFluxDiffusive[3 * ion + 2]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                       Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxDiffusive[3 * ion + 2]);
        }
    }

    if (vis_db->getWithDefault<bool>("save_ion_flux_electrical", false)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            // x-component of electro-migrational flux
            sprintf(VisName, "Ion%zu_FluxElectrical_x", ion + 1);
            IonFluxElectrical[3 * ion + 0]->name = VisName;
            IonFluxElectrical[3 * ion + 0]->type =
                IO::VariableType::VolumeVariable;
            IonFluxElectrical[3 * ion + 0]->dim = 1;
            IonFluxElectrical[3 * ion + 0]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                        Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxElectrical[3 * ion + 0]);
            // y-component of electro-migrational flux
            sprintf(VisName, "Ion%zu_FluxElectrical_y", ion + 1);
            IonFluxElectrical[3 * ion + 1]->name = VisName;
            IonFluxElectrical[3 * ion + 1]->type =
                IO::VariableType::VolumeVariable;
            IonFluxElectrical[3 * ion + 1]->dim = 1;
            IonFluxElectrical[3 * ion + 1]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                        Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxElectrical[3 * ion + 1]);
            // z-component of electro-migrational flux
            sprintf(VisName, "Ion%zu_FluxElectrical_z", ion + 1);
            IonFluxElectrical[3 * ion + 2]->name = VisName;
            IonFluxElectrical[3 * ion + 2]->type =
                IO::VariableType::VolumeVariable;
            IonFluxElectrical[3 * ion + 2]->dim = 1;
            IonFluxElectrical[3 * ion + 2]->data.resize(Dm->Nx - 2, Dm->Ny - 2,
                                                        Dm->Nz - 2);
            visData[0].vars.push_back(IonFluxElectrical[3 * ion + 2]);
        }
    }

    if (vis_db->getWithDefault<bool>("save_electric_field", false)) {
        ElectricFieldVar_x->name = "ElectricField_x";
        ElectricFieldVar_x->type = IO::VariableType::VolumeVariable;
        ElectricFieldVar_x->dim = 1;
        ElectricFieldVar_x->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(ElectricFieldVar_x);
        ElectricFieldVar_y->name = "ElectricField_y";
        ElectricFieldVar_y->type = IO::VariableType::VolumeVariable;
        ElectricFieldVar_y->dim = 1;
        ElectricFieldVar_y->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(ElectricFieldVar_y);
        ElectricFieldVar_z->name = "ElectricField_z";
        ElectricFieldVar_z->type = IO::VariableType::VolumeVariable;
        ElectricFieldVar_z->dim = 1;
        ElectricFieldVar_z->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(ElectricFieldVar_z);
    }
    //--------------------------------------------------------------------------------------------------------------------

    //------------------------------------Save All Variables--------------------------------------------------------------
    if (vis_db->getWithDefault<bool>("save_electric_potential", true)) {
        ASSERT(visData[0].vars[0]->name == "ElectricPotential");
        Poisson.getElectricPotential(ElectricalPotential);
        Array<double> &ElectricPotentialData = visData[0].vars[0]->data;
        fillData.copy(ElectricalPotential, ElectricPotentialData);
    }

    if (vis_db->getWithDefault<bool>("save_concentration", true)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {
            sprintf(VisName, "IonConcentration_%zu", ion + 1);
            //IonConcentration[ion]->name = VisName;
            ASSERT(visData[0].vars[1 + ion]->name == VisName);
            Array<double> &IonConcentrationData =
                visData[0].vars[1 + ion]->data;
            Ion.getIonConcentration(Rho, ion);
            fillData.copy(Rho, IonConcentrationData);
        }
    }

    if (vis_db->getWithDefault<bool>("save_ion_flux_diffusive", false)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {

            // x-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_x", ion + 1);
            //IonFluxDiffusive[3*ion+0]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species + 3 * ion + 0]
                       ->name == VisName);
            // y-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_y", ion + 1);
            //IonFluxDiffusive[3*ion+1]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species + 3 * ion + 1]
                       ->name == VisName);
            // z-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxDiffusive_z", ion + 1);
            //IonFluxDiffusive[3*ion+2]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species + 3 * ion + 2]
                       ->name == VisName);

            Array<double> &IonFluxData_x =
                visData[0].vars[4 + Ion.number_ion_species + 3 * ion + 0]->data;
            Array<double> &IonFluxData_y =
                visData[0].vars[4 + Ion.number_ion_species + 3 * ion + 1]->data;
            Array<double> &IonFluxData_z =
                visData[0].vars[4 + Ion.number_ion_species + 3 * ion + 2]->data;
            Ion.getIonFluxDiffusive(IonFluxDiffusive_x, IonFluxDiffusive_y,
                                    IonFluxDiffusive_z, ion);
            fillData.copy(IonFluxDiffusive_x, IonFluxData_x);
            fillData.copy(IonFluxDiffusive_y, IonFluxData_y);
            fillData.copy(IonFluxDiffusive_z, IonFluxData_z);
        }
    }

    if (vis_db->getWithDefault<bool>("save_ion_flux_electrical", false)) {
        for (size_t ion = 0; ion < Ion.number_ion_species; ion++) {

            // x-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxElectrical_x", ion + 1);
            //IonFluxDiffusive[3*ion+0]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 0]
                       ->name == VisName);
            // y-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxElectrical_y", ion + 1);
            //IonFluxDiffusive[3*ion+1]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 1]
                       ->name == VisName);
            // z-component of diffusive flux
            sprintf(VisName, "Ion%zu_FluxElectrical_z", ion + 1);
            //IonFluxDiffusive[3*ion+2]->name = VisName;
            ASSERT(visData[0]
                       .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 2]
                       ->name == VisName);

            Array<double> &IonFluxData_x =
                visData[0]
                    .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 0]
                    ->data;
            Array<double> &IonFluxData_y =
                visData[0]
                    .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 1]
                    ->data;
            Array<double> &IonFluxData_z =
                visData[0]
                    .vars[4 + Ion.number_ion_species * (1 + 6) + 3 * ion + 2]
                    ->data;
            Ion.getIonFluxElectrical(IonFluxElectrical_x, IonFluxElectrical_y,
                                     IonFluxElectrical_z, ion);
            fillData.copy(IonFluxElectrical_x, IonFluxData_x);
            fillData.copy(IonFluxElectrical_y, IonFluxData_y);
            fillData.copy(IonFluxElectrical_z, IonFluxData_z);
        }
    }

    if (vis_db->getWithDefault<bool>("save_electric_field", false)) {
        ASSERT(
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 0]->name ==
            "ElectricField_x");
        ASSERT(
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 1]->name ==
            "ElectricField_y");
        ASSERT(
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 2]->name ==
            "ElectricField_z");
        Poisson.getElectricField(ElectricalField_x, ElectricalField_y,
                                 ElectricalField_z);
        Array<double> &ElectricalFieldxData =
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 0]->data;
        Array<double> &ElectricalFieldyData =
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 1]->data;
        Array<double> &ElectricalFieldzData =
            visData[0].vars[4 + Ion.number_ion_species * (1 + 9) + 2]->data;
        fillData.copy(ElectricalField_x, ElectricalFieldxData);
        fillData.copy(ElectricalField_y, ElectricalFieldyData);
        fillData.copy(ElectricalField_z, ElectricalFieldzData);
    }

    if (vis_db->getWithDefault<bool>("write_silo", true))
        IO::writeData(timestep, visData, Dm->Comm);
    //--------------------------------------------------------------------------------------------------------------------
    /*    if (vis_db->getWithDefault<bool>( "save_8bit_raw", true )){
    	char CurrentIDFilename[40];
    	sprintf(CurrentIDFilename,"id_t%d.raw",timestep);
    	Averages.AggregateLabels(CurrentIDFilename);
    }
*/
}
