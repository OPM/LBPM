// Header file for two-phase averaging class
#ifndef Minkowski_INC
#define Minkowski_INC

#include <vector>

#include "analysis/decl.h"
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
	int n_obj_pts;
	int n_obj_tris;
	//...........................................................................
	int nc;
	int kstart,kfinish;

	double isovalue;
	double Volume;
	// initialize lists for vertices for surfaces, common line
	DTMutableList<Point> obj_pts;
	DTMutableList<Point> tmp;

	// initialize triangle lists for surfaces
	IntArray obj_tris;

	// Temporary storage arrays
	DoubleArray CubeValues;
	DoubleArray Values;
	DoubleArray NormalVector;

	DoubleArray RecvBuffer;

	char *TempID;

	// CSV / text file where time history of averages is saved
	FILE *LOGFILE;

public:
	//...........................................................................
	std::shared_ptr <Domain> Dm;
	//...........................................................................
	// Averaging variables
	//...........................................................................
	// local averages (to each MPI process)
	double vol_n;						// volumes the exclude the interfacial region
	// Global averages (all processes)
	double vol_n_global;			// volumes the exclude the interfacial region
	double euler,Kn,Jn,An;
	double euler_global,Kn_global,Jn_global,An_global;
	//...........................................................................
	int Nx,Ny,Nz;
	IntArray PhaseID;	// Phase ID array (solid=0, non-wetting=1, wetting=2)
	DoubleArray SDn;
	DoubleArray MeanCurvature;
	DoubleArray GaussCurvature;
	DoubleArray SDn_x;		// Gradient of the signed distance
	DoubleArray SDn_y;
	DoubleArray SDn_z;

	double V();
	double A();
	double J();
	double X();
	
	//...........................................................................
	Minkowski(std::shared_ptr <Domain> Dm);
	~Minkowski();
	void Initialize();
	void UpdateMeshValues();
	void ComputeLocal();
	void Reduce();
	void NonDimensionalize(double D);
	void PrintAll();
	int GetCubeLabel(int i, int j, int k, IntArray &BlobLabel);
	void SortBlobs();
	
};

#endif

