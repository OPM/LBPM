#include "analysis/ElectroChemistry.h"

void ElectroChemistryAnalyzer::ElectroChemistryAnalyzer(std::shared_ptr <Domain> Dm):
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
}

void ElectroChemistryAnalyzer::~ElectroChemistryAnalyzer(){
	
}
