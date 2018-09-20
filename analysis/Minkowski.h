// Header file for two-phase averaging class
#ifndef Minkowski_INC
#define Minkowski_INC

#include <vector>

#include "analysis/dcel.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "analysis/analysis.h"

#include "shared_ptr.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"


class Minkowski{
	//...........................................................................
	int kstart,kfinish;

	double isovalue;
	double Volume;

	// CSV / text file where time history of averages is saved
	FILE *LOGFILE;

public:
	//...........................................................................
	std::shared_ptr <Domain> Dm;
	//...........................................................................
	// Averaging variables
	//...........................................................................
	// local averages (to each MPI process)
	double Ai,Ji,Xi,Vi;
	// Global averages (all processes)
	double Ai_global,Ji_global,Xi_global,Vi_global;

	//...........................................................................
	int Nx,Ny,Nz;
	double V();
	double A();
	double J();
	double X();
	
	//...........................................................................
	Minkowski();
	Minkowski(std::shared_ptr <Domain> Dm);
	~Minkowski();
	void ComputeScalar(const DoubleArray& Field, const double isovalue);
	void PrintAll();

};

#endif

