#include "analysis/ElectroChemistry.h"

ElectroChemistryAnalyzer::ElectroChemistryAnalyzer(std::shared_ptr <Domain> dm):
	Dm(dm){
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
	
	DoubleArray Rho;	           // density field 
	DoubleArray ChemicalPotential;	   // density field 
	DoubleArray ElectricalPotential;	// density field 
	DoubleArray Pressure; 	// pressure field
	DoubleArray Vel_x;		// velocity field
	DoubleArray Vel_y;
	DoubleArray Vel_z;
	DoubleArray SDs;
	
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
			fprintf(TIMELOG,"sw krw krn vw vn pw pn\n");				
		}
	}

}

ElectroChemistryAnalyzer::~ElectroChemistryAnalyzer(){
	
}

void ElectroChemistryAnalyzer::SetParams(){
	
}

void ElectroChemistryAnalyzer::Basic(){
	
}

void ElectroChemistryAnalyzer::Write(int time){
	
}
