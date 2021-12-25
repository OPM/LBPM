/* Flow adaptor class for multiphase flow methods */

#ifndef ScaLBL_Membrane_INC
#define ScaLBL_Membrane_INC
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"

/**
 * \class Membrane
 * @brief 
 * The Membrane class operates on ScaLBL data structures to insert membrane
 * 
 */

class Membrane {
public:
    /**
    * \brief Create a flow adaptor to operate on the LB model
    * @param         ScaLBL - originating data structures 
    * @param 		 neighborList - list of neighbors for each site
    */
    Membrane(ScaLBL_Communicator &ScaLBL, int *neighborList);

    /**
	* \brief Destructor
	*/
    ~Membrane();

    /**
    * \brief   Fast-forward interface motion
    * \details Accelerate the movement of interfaces based on the time derivative
    *       Optional keys to control behavior can be specified in the input database:
    *          move_interface_cutoff  -- identifies the diffuse interface region
    *          move_interface_factor -- determines how much to ``fast forward" 
    * @param Distance       - signed distance to membrane 
    * @param Map       - mapping between regular layout and compact layout
    * @param neighborLIst - list of neighbors for each site
    * @param membrane - links that form the membrane
    * @param Np - number of sites in compact layout
    */
    int Create(DoubleArray &Distance, IntArray &Map, int *neighborList, int *membrane, int Np);


private:

};
#endif