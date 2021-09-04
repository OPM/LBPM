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

class FlowAdaptor{
public:
	FlowAdaptor(ScaLBL_ColorModel &M);
	~FlowAdaptor();
	double MoveInterface(ScaLBL_ColorModel &M);
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