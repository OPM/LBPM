/*
 * Sub-phase averaging tools
 */

#ifndef ElectroChem_INC
#define ElectroChem_INC

#include <vector>
#include "common/Domain.h"
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "analysis/distance.h"
#include "analysis/Minkowski.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"

class ElectroChemistryAnalyzer{
public:
	std::shared_ptr <Domain> Dm;
	double Volume;
	// input variables
	double rho_n, rho_w;
	double nu_n, nu_w;
	double gamma_wn, beta;
	double Fx, Fy, Fz;

    //...........................................................................
    int Nx,Ny,Nz;
	DoubleArray Rho;	           // density field 
	DoubleArray ChemicalPotential;	   // density field 
	DoubleArray ElectricalPotential;	// density field 
	DoubleArray Pressure; 	// pressure field
	DoubleArray Vel_x;		// velocity field
	DoubleArray Vel_y;
	DoubleArray Vel_z;
	DoubleArray SDs;

	ElectroChemistryAnalyzer(std::shared_ptr <Domain> Dm);
	~ElectroChemistryAnalyzer();
	
	void SetParams();
	void Basic();
	void Write(int time);

private:
	FILE *TIMELOG;
};
#endif

