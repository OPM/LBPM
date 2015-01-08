// Unit test for TwoPhase averaging class

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "TwoPhase.h"
#include "Extras.h"
#include "D3Q19.h"
#include "D3Q7.h"
#include "Color.h"
#include "common/MPI.h"
#include "Communication.h"
#include "IO/Mesh.h"
#include "IO/Writer.h"
#include "ProfilerApp.h"

int main(int argc, char **argv)
{
	int rank,npx,npy,npz;
	int Nx,Ny,Nz;
	double Lx,Ly,Lz;
	Nx=Ny=Nz=40;
	rank=0;
	npx=npy=npz=1;
	Lx=Ly=Lz=1.0;

	Domain Dm(Nx,Ny,Nz,rank,npx,npy,npz,Lx,Ly,Lz);

	TwoPhase Averages(Dm);
	Averages.SetupCubes(Dm);
	Averages.Initialize();
	Averages.Compute();

	return 0;
}
