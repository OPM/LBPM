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

/**
 * \class Minkowski
 *
 * @brief 
 * The Minkowski class is constructed to analyze the geometric properties of structures based on the Minkowski functionals
 * 
 */


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
    /**
    * \brief   Null constructor
    */
	Minkowski(){};//NULL CONSTRUCTOR
	
    /**
    * \brief   Constructor based on an existing Domain
    * @param Dm    - Domain structure 
    */
	Minkowski(std::shared_ptr <Domain> Dm);
	~Minkowski();
	
    /**
    * \brief   Compute scalar minkowski functionals	
	 *  step 1. compute the distance to an object 
	 *  step 2. construct dcel to represent the isosurface 
	 *  step 3. compute the scalar Minkowski functionals
	 * THIS ALGORITHM ASSUMES THAT id() is populated with phase id to distinguish objects
	 *    0 - labels the object
	 *    1 - labels everything else
	 */    
	void MeasureObject();
	
	void MeasureObject(double factor, const DoubleArray &Phi);
	
    /**
    * \details   Compute scalar minkowski functionals for connected part of a structure
    *   step 1. compute connected components and extract largest region by volume
    *  step 2. compute the distance to the connected part of the structure 
    *  step 3. construct dcel to represent the isosurface 
    *  step 4. compute the scalar Minkowski functionals
    * THIS ALGORITHM ASSUMES THAT id() is populated with phase id to distinguish objects
    *    0 - labels the object
    *    1 - labels everything else
    */    
	int MeasureConnectedPathway();

	int MeasureConnectedPathway(double factor, const DoubleArray &Phi);
	
    /**
    * \brief   Compute scalar minkowski functionals
    * \details Construct an isosurface and return the geometric invariants based on the triangulated list
    * @param isovalue        - threshold value to use to determine iso-surface
    * @param Field           - DoubleArray containing the field to threshold
    */
	void ComputeScalar(const DoubleArray& Field, const double isovalue);

    /**
    * \brief   print the scalar invariants
    */
	void PrintAll();

};

#endif

