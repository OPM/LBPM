#include "analysis/ElectroChemistry.h"

ElectroChemistryAnalyzer::ElectroChemistryAnalyzer(std::shared_ptr <Domain> dm):
	Dm(dm),
	fillData(dm->Comm,dm->rank_info,{dm->Nx-2,dm->Ny-2,dm->Nz-2},{1,1,1},0,1)
{
	
    MPI_Comm_dup(dm->Comm,&comm);
	Nx=dm->Nx; Ny=dm->Ny; Nz=dm->Nz;
	Volume=(Nx-2)*(Ny-2)*(Nz-2)*Dm->nprocx()*Dm->nprocy()*Dm->nprocz()*1.0;
	
	ChemicalPotential.resize(Nx,Ny,Nz);       ChemicalPotential.fill(0);
	ElectricalPotential.resize(Nx,Ny,Nz);      ElectricalPotential.fill(0);
	Pressure.resize(Nx,Ny,Nz);     	Pressure.fill(0);
	Rho.resize(Nx,Ny,Nz);       	Rho.fill(0);
	Vel_x.resize(Nx,Ny,Nz);         Vel_x.fill(0);	    // Gradient of the phase indicator field
	Vel_y.resize(Nx,Ny,Nz);         Vel_y.fill(0);
	Vel_z.resize(Nx,Ny,Nz);         Vel_z.fill(0);
	SDs.resize(Nx,Ny,Nz);         	SDs.fill(0);
	
	if (Dm->rank()==0){
		bool WriteHeader=false;
		TIMELOG = fopen("electrokinetic.csv","r");
		if (TIMELOG != NULL)
			fclose(TIMELOG);
		else
			WriteHeader=true;

		TIMELOG = fopen("electrokinetic.csv","a+");
		if (WriteHeader)
		{
			// If timelog is empty, write a short header to list the averages
			//fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
			fprintf(TIMELOG,"TBD TBD\n");				
		}
	}

}

ElectroChemistryAnalyzer::~ElectroChemistryAnalyzer(){
	if (Dm->rank()==0){
		fclose(TIMELOG);
	}
}

void ElectroChemistryAnalyzer::SetParams(){
	
}

void ElectroChemistryAnalyzer::Basic(ScaLBL_IonModel &Ion, ScaLBL_Poisson &Poisson, ScaLBL_StokesModel &Stokes, int timestep){

	int i,j,k;
	Poisson.getElectricPotential(ElectricalPotential);
	
	/* local sub-domain averages */
	double rho_avg_local[Ion.number_ion_species];
	double rho_mu_avg_local[Ion.number_ion_species];
	double rho_mu_fluctuation_local[Ion.number_ion_species];	
	double rho_psi_avg_local[Ion.number_ion_species];
	double rho_psi_fluctuation_local[Ion.number_ion_species];
	/* global averages */
	double rho_avg_global[Ion.number_ion_species];
	double rho_mu_avg_global[Ion.number_ion_species];
	double rho_mu_fluctuation_global[Ion.number_ion_species];
	double rho_psi_avg_global[Ion.number_ion_species];
	double rho_psi_fluctuation_global[Ion.number_ion_species];
	
	for (int ion=0; ion<Ion.number_ion_species; ion++){
		rho_avg_local[ion] = 0.0;
		rho_mu_avg_local[ion] = 0.0;
		rho_psi_avg_local[ion] = 0.0;
		Ion.getIonConcentration(Rho,ion);
		/* Compute averages for each ion */
		for (k=1; k<Nz; k++){
			for (j=1; j<Ny; j++){
				for (i=1; i<Nx; i++){
					rho_avg_local[ion] += Rho(i,j,k);
					rho_mu_avg_local[ion] += Rho(i,j,k)*Rho(i,j,k);
					rho_psi_avg_local[ion] += Rho(i,j,k)*ElectricalPotential(i,j,k);
				}
			}
		}
		rho_avg_global[ion]=sumReduce( Dm->Comm, rho_avg_local[ion]);
		rho_mu_avg_global[ion]=sumReduce( Dm->Comm, rho_mu_avg_local[ion]);
		rho_psi_avg_global[ion]=sumReduce( Dm->Comm, rho_psi_avg_local[ion]);
		
		rho_mu_avg_global[ion] /= rho_avg_global[ion];
		rho_psi_avg_global[ion] /= rho_avg_global[ion];
	}
	
	for (int ion=0; ion<Ion.number_ion_species; ion++){
		rho_mu_fluctuation_local[ion] = 0.0;
		rho_psi_fluctuation_local[ion] = 0.0;
		/* Compute averages for each ion */
		for (k=1; k<Nz; k++){
			for (j=1; j<Ny; j++){
				for (i=1; i<Nx; i++){
					rho_mu_fluctuation_local[ion] += (Rho(i,j,k)*Rho(i,j,k) - rho_mu_avg_global[ion]);
					rho_psi_fluctuation_local[ion] += (Rho(i,j,k)*ElectricalPotential(i,j,k) - rho_psi_avg_global[ion]);
				}
			}
		}
		rho_mu_fluctuation_global[ion]=sumReduce( Dm->Comm, rho_mu_fluctuation_local[ion]);
		rho_psi_fluctuation_global[ion]=sumReduce( Dm->Comm, rho_psi_fluctuation_local[ion]);
	}
	
	if (Dm->rank()==0){	
		fprintf(TIMELOG,"%i ",timestep); 
		for (int ion=0; ion<Ion.number_ion_species; ion++){
			fprintf(TIMELOG,"%.8g ",rho_avg_global[ion]);
			fprintf(TIMELOG,"%.8g ",rho_mu_avg_global[ion]);
			fprintf(TIMELOG,"%.8g ",rho_psi_avg_global[ion]);
			fprintf(TIMELOG,"%.8g ",rho_mu_fluctuation_global[ion]);
			fprintf(TIMELOG,"%.8g ",rho_psi_fluctuation_global[ion]);
		}
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

void ElectroChemistryAnalyzer::WriteVis( ScaLBL_IonModel &Ion, ScaLBL_Poisson &Poisson, ScaLBL_StokesModel &Stokes, std::shared_ptr<Database> input_db, int timestep){
	
	auto vis_db =  input_db->getDatabase( "Visualization" );
    char VisName[40];

    IO::initialize("","silo","false");
    // Create the MeshDataStruct    
    visData.resize(1);

    visData[0].meshName = "domain";
    visData[0].mesh = std::make_shared<IO::DomainMesh>( Dm->rank_info,Dm->Nx-2,Dm->Ny-2,Dm->Nz-2,Dm->Lx,Dm->Ly,Dm->Lz );
    auto ElectricPotential = std::make_shared<IO::Variable>();
    auto IonConcentration = std::make_shared<IO::Variable>();
    auto VxVar = std::make_shared<IO::Variable>();
    auto VyVar = std::make_shared<IO::Variable>();
    auto VzVar = std::make_shared<IO::Variable>();
    
    if (vis_db->getWithDefault<bool>( "save_electric_potential", true )){
    	ElectricPotential->name = "ElectricPotential";
    	ElectricPotential->type = IO::VariableType::VolumeVariable;
    	ElectricPotential->dim = 1;
    	ElectricPotential->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
        visData[0].vars.push_back(ElectricPotential);
    }

    if (vis_db->getWithDefault<bool>( "save_concentration", true )){
    	for (int ion=0; ion<Ion.number_ion_species; ion++){
    		sprintf(VisName,"IonConcentration_%i",ion);
    		IonConcentration->name = VisName;
    		IonConcentration->type = IO::VariableType::VolumeVariable;
    		IonConcentration->dim = 1;
    		IonConcentration->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    		visData[0].vars.push_back(IonConcentration);
    	}
    }
    if (vis_db->getWithDefault<bool>( "save_velocity", false )){
        VxVar->name = "Velocity_x";
        VxVar->type = IO::VariableType::VolumeVariable;
        VxVar->dim = 1;
        VxVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
        visData[0].vars.push_back(VxVar);
        VyVar->name = "Velocity_y";
        VyVar->type = IO::VariableType::VolumeVariable;
        VyVar->dim = 1;
        VyVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
        visData[0].vars.push_back(VyVar);
        VzVar->name = "Velocity_z";
        VzVar->type = IO::VariableType::VolumeVariable;
        VzVar->dim = 1;
        VzVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
        visData[0].vars.push_back(VzVar);
    }
    
    if (vis_db->getWithDefault<bool>( "save_electric_potential", true )){
    	ASSERT(visData[0].vars[0]->name=="ElectricPotential");
    	Poisson.getElectricPotential(ElectricalPotential);
    	Array<double>& ElectricPotentialData = visData[0].vars[0]->data;
    	fillData.copy(ElectricalPotential,ElectricPotentialData);
    }

    if (vis_db->getWithDefault<bool>( "save_concentration", true )){
    	for (int ion=0; ion<Ion.number_ion_species; ion++){
    		sprintf(VisName,"IonConcentration_%i",ion);
    		IonConcentration->name = VisName;
    		ASSERT(visData[0].vars[1]->name==VisName);
    		Array<double>& IonConcentrationData = visData[0].vars[1]->data;
    		Ion.getIonConcentration(Rho,ion);
    		fillData.copy(Rho,IonConcentrationData);
    	}
    }

    if (vis_db->getWithDefault<bool>( "save_velocity", false )){
    	ASSERT(visData[0].vars[2]->name=="Velocity_x");
    	ASSERT(visData[0].vars[3]->name=="Velocity_y");
    	ASSERT(visData[0].vars[4]->name=="Velocity_z");
    	Stokes.getVelocity(Vel_x,Vel_y,Vel_z);
    	Array<double>& VelxData = visData[0].vars[2]->data;
    	Array<double>& VelyData = visData[0].vars[3]->data;
    	Array<double>& VelzData = visData[0].vars[4]->data;
    	fillData.copy(Vel_x,VelxData);
    	fillData.copy(Vel_y,VelyData);
    	fillData.copy(Vel_z,VelzData);
    }
    
    if (vis_db->getWithDefault<bool>( "write_silo", true ))
    	IO::writeData( timestep, visData, comm );

/*    if (vis_db->getWithDefault<bool>( "save_8bit_raw", true )){
    	char CurrentIDFilename[40];
    	sprintf(CurrentIDFilename,"id_t%d.raw",timestep);
    	Averages.AggregateLabels(CurrentIDFilename);
    }
*/	
}
