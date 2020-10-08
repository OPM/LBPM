/*
 * Sub-phase averaging tools
 */

#ifndef GreyPhase_INC
#define GreyPhase_INC

#include <vector>
#include "common/Domain.h"
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"

class GreyPhase{
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
	IntArray PhaseID;		// Phase ID array 
	DoubleArray Rho_n;	// density field 
	DoubleArray Rho_w;	// density field 
	DoubleArray Phi;		// phase indicator field
	DoubleArray DelPhi;		// Magnitude of Gradient of the phase indicator field
	DoubleArray Pressure; 	// pressure field
	DoubleArray Vel_x;		// velocity field
	DoubleArray Vel_y;
	DoubleArray Vel_z;

	GreyPhase(std::shared_ptr <Domain> Dm);
	~GreyPhase();
	
	void SetParams(double rhoA, double rhoB, double tauA, double tauB, double force_x, double force_y, double force_z, double alpha, double beta);
	void Basic();
	void Write(int time);

private:
	FILE *TIMELOG;
};

#endif

