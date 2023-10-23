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
#include "analysis/FreeEnergy.h"

FreeEnergyAnalyzer::FreeEnergyAnalyzer(std::shared_ptr<Domain> dm) : Dm(dm) {

    Nx = dm->Nx;
    Ny = dm->Ny;
    Nz = dm->Nz;
    Volume = (Nx - 2) * (Ny - 2) * (Nz - 2) * Dm->nprocx() * Dm->nprocy() *
             Dm->nprocz() * 1.0;

    ChemicalPotential.resize(Nx, Ny, Nz);
    ChemicalPotential.fill(0);
    Phi.resize(Nx, Ny, Nz);
    Phi.fill(0);
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

    if (Dm->rank() == 0) {
        bool WriteHeader = false;
        TIMELOG = fopen("free.csv", "r");
        if (TIMELOG != NULL)
            fclose(TIMELOG);
        else
            WriteHeader = true;

        TIMELOG = fopen("free.csv", "a+");
        if (WriteHeader) {
            // If timelog is empty, write a short header to list the averages
            //fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
            fprintf(TIMELOG, "timestep\n");
        }
    }
}

FreeEnergyAnalyzer::~FreeEnergyAnalyzer() {
    if (Dm->rank() == 0) {
        fclose(TIMELOG);
    }
}

void FreeEnergyAnalyzer::SetParams() {}

void FreeEnergyAnalyzer::Basic(ScaLBL_FreeLeeModel &LeeModel, int timestep) {

    if (Dm->rank() == 0) {
        fprintf(TIMELOG, "%i ", timestep);
        /*for (int ion=0; ion<Ion.number_ion_species; ion++){
			fprintf(TIMELOG,"%.8g ",rho_avg_global[ion]);
			fprintf(TIMELOG,"%.8g ",rho_mu_avg_global[ion]);
			fprintf(TIMELOG,"%.8g ",rho_psi_avg_global[ion]);
			fprintf(TIMELOG,"%.8g ",rho_mu_fluctuation_global[ion]);
			fprintf(TIMELOG,"%.8g ",rho_psi_fluctuation_global[ion]);
		}
		*/
        fprintf(TIMELOG, "\n");
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

void FreeEnergyAnalyzer::WriteVis(ScaLBL_FreeLeeModel &LeeModel,
                                  std::shared_ptr<Database> input_db,
                                  int timestep) {

    auto vis_db = input_db->getDatabase("Visualization");

    std::vector<IO::MeshDataStruct> visData;
    fillHalo<double> fillData(Dm->Comm, Dm->rank_info,
                              {Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2}, {1, 1, 1},
                              0, 1);

    IO::initialize("", "silo", "false");
    // Create the MeshDataStruct
    visData.resize(1);

    visData[0].meshName = "domain";
    visData[0].mesh =
        std::make_shared<IO::DomainMesh>(Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2,
                                         Dm->Nz - 2, Dm->Lx, Dm->Ly, Dm->Lz);
    auto VisPhase = std::make_shared<IO::Variable>();
    auto VisPressure = std::make_shared<IO::Variable>();
    auto VisChemicalPotential = std::make_shared<IO::Variable>();
    auto VxVar = std::make_shared<IO::Variable>();
    auto VyVar = std::make_shared<IO::Variable>();
    auto VzVar = std::make_shared<IO::Variable>();

    if (vis_db->getWithDefault<bool>("save_phase_field", true)) {
        VisPhase->name = "Phase";
        VisPhase->type = IO::VariableType::VolumeVariable;
        VisPhase->dim = 1;
        VisPhase->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VisPhase);
    }

    if (vis_db->getWithDefault<bool>("save_potential", true)) {

        VisPressure->name = "Pressure";
        VisPressure->type = IO::VariableType::VolumeVariable;
        VisPressure->dim = 1;
        VisPressure->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VisPressure);

        VisChemicalPotential->name = "ChemicalPotential";
        VisChemicalPotential->type = IO::VariableType::VolumeVariable;
        VisChemicalPotential->dim = 1;
        VisChemicalPotential->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VisChemicalPotential);
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

    if (vis_db->getWithDefault<bool>("save_phase", true)) {
        ASSERT(visData[0].vars[0]->name == "Phase");
        LeeModel.getPhase(Phi);
        Array<double> &PhaseData = visData[0].vars[0]->data;
        fillData.copy(Phi, PhaseData);
    }

    if (vis_db->getWithDefault<bool>("save_potential", true)) {
        ASSERT(visData[0].vars[1]->name == "Pressure");
        LeeModel.getPotential(Pressure, ChemicalPotential);
        Array<double> &PressureData = visData[0].vars[1]->data;
        fillData.copy(Pressure, PressureData);

        ASSERT(visData[0].vars[2]->name == "ChemicalPotential");
        Array<double> &ChemicalPotentialData = visData[0].vars[2]->data;
        fillData.copy(ChemicalPotential, ChemicalPotentialData);
    }

    if (vis_db->getWithDefault<bool>("save_velocity", false)) {
        ASSERT(visData[0].vars[3]->name == "Velocity_x");
        ASSERT(visData[0].vars[4]->name == "Velocity_y");
        ASSERT(visData[0].vars[5]->name == "Velocity_z");
        LeeModel.getVelocity(Vel_x, Vel_y, Vel_z);
        Array<double> &VelxData = visData[0].vars[3]->data;
        Array<double> &VelyData = visData[0].vars[4]->data;
        Array<double> &VelzData = visData[0].vars[5]->data;
        fillData.copy(Vel_x, VelxData);
        fillData.copy(Vel_y, VelyData);
        fillData.copy(Vel_z, VelzData);
    }

    if (vis_db->getWithDefault<bool>("write_silo", true))
        IO::writeData(timestep, visData, Dm->Comm);

    /*    if (vis_db->getWithDefault<bool>( "save_8bit_raw", true )){
    	char CurrentIDFilename[40];
    	sprintf(CurrentIDFilename,"id_t%d.raw",timestep);
    	Averages.AggregateLabels(CurrentIDFilename);
    }
*/
}
