/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

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
// Header file for two-phase averaging class
#ifndef Minkowski_INC
#define Minkowski_INC

#include <vector>

#include "analysis/pmmc.h"
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

