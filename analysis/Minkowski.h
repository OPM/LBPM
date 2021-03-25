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
// Header file for two-phase averaging class
#ifndef Minkowski_INC
#define Minkowski_INC

#include <memory>
#include <vector>

#include "analysis/dcel.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "analysis/distance.h"
#include "analysis/filters.h"

#include "common/Utilities.h"
#include "common/MPI.h"
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
	Array <char> id;
	Array <int> label;
	Array <double> distance;
	//...........................................................................
	// Averaging variables
	//...........................................................................
	// local averages (to each MPI process)
	double Ai,Ji,Xi,Vi;
	// Global averages (all processes)
	double Ai_global,Ji_global,Xi_global,Vi_global;
	int n_connected_components;
	//...........................................................................
	int Nx,Ny,Nz;
	double V(){
		return Vi;
	}
	double A(){
		return Ai;
	}
	double H(){
		return Ji;
	}
	double X(){
		return Xi;
	}
		
	//..........................................................................
	Minkowski(){};//NULL CONSTRUCTOR
	Minkowski(std::shared_ptr <Domain> Dm);
	~Minkowski();
	void MeasureObject();
	void MeasureObject(double factor, const DoubleArray &Phi);
	int MeasureConnectedPathway();
	int MeasureConnectedPathway(double factor, const DoubleArray &Phi);
	void ComputeScalar(const DoubleArray& Field, const double isovalue);

	void PrintAll();

};

#endif

