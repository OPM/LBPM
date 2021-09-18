/* Flow adaptor class for multiphase flow methods */

#ifndef ScaLBL_FlowAdaptor_INC
#define ScaLBL_FlowAdaptor_INC
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "models/ColorModel.h"


/**
 * \class FlowAdaptor
 *
 * @brief 
 * The FlowAdaptor class operates on a lattice Boltzmann model to alter the flow conditions
 * 
 */

class FlowAdaptor{
public:
	

    /**
    * \brief Create a flow adaptor to operate on the LB model
    * @param M        ScaLBL_ColorModel 
    */
	FlowAdaptor(ScaLBL_ColorModel &M);
	
    /**
	* \brief Destructor
	*/
	~FlowAdaptor();
	
    /**
    * \brief   Fast-forward interface motion
    * \details Accelerate the movement of interfaces based on the time derivative
    *       Optional keys to control behavior can be specified in the input database:
    *          move_interface_cutoff  -- identifies the diffuse interface region
    *          move_interface_factor -- determines how much to ``fast forward" 
    * @param M        ScaLBL_ColorModel 
    */
	double MoveInterface(ScaLBL_ColorModel &M);
	
    /**
     * \brief image re-initialization
     * \details Re-initialize LB simulation from image data  
     * @param M        ScaLBL_ColorModel 
     * @param Filename    name of input file to be used to read image 
    */
	double ImageInit(ScaLBL_ColorModel &M, std::string Filename);
	
	double ShellAggregation(ScaLBL_ColorModel &M, const double delta_volume);
	double UpdateFractionalFlow(ScaLBL_ColorModel &M);
	double SeedPhaseField(ScaLBL_ColorModel &M, const double seed_water_in_oil);
	void Flatten(ScaLBL_ColorModel &M);
	DoubleArray phi;
	DoubleArray phi_t;
private:
	int Nx, Ny, Nz;
	int timestep;
	int timestep_previous;
};
#endif