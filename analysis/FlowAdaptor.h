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
 * @brief 
 * The FlowAdaptor class operates on a lattice Boltzmann model to alter the flow conditions
 * 
 */

class FlowAdaptor {
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
     * \brief Image re-initialization
     * \details Re-initialize LB simulation from image data
     * @param M        ScaLBL_ColorModel 
     * @param Filename    name of input file to be used to read image 
    */
    double ImageInit(ScaLBL_ColorModel &M, std::string Filename);

    /**
     * \details Update volume fraction based on morphological algorithm. Dilation / erosion algorithm will be applied to 
     * grow / shrink the phase regions
     * @param M        ScaLBL_ColorModel 
     * @param delta_volume    target change in volume fraction 
    */
    double ShellAggregation(ScaLBL_ColorModel &M, const double delta_volume);

    /**
     * \details  Update fractional flow condition. Mass will be preferentially added or removed from 
     * phase regions based on where flow is occurring
     * @param M        ScaLBL_ColorModel 
    */
    double UpdateFractionalFlow(ScaLBL_ColorModel &M);

    /**
     * \brief image re-initialization
     * \details Re-initialize LB simulation from image data  
     * @param M        ScaLBL_ColorModel 
     * @param seed_water_in_oil    controls amount of mass to randomly seed into fluids 
    */
    double SeedPhaseField(ScaLBL_ColorModel &M, const double seed_water_in_oil);

    /**
     * \brief Re-initialize LB simulation
     * @param M        ScaLBL_ColorModel 
    */
    void Flatten(ScaLBL_ColorModel &M);
    DoubleArray phi;
    DoubleArray phi_t;

private:
    int Nx, Ny, Nz;
    int timestep;
    int timestep_previous;
};
#endif