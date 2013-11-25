#include <iostream>
#include <math.h>
#include "pmmc.h"
#include "Domain.h"

using namespace std;

/*
 *  Compare the measured and analytical curvature for a sphere
 *
 */

inline void MeshCurvature(DoubleArray &f, DoubleArray  &MeanCurvature, DoubleArray &GaussCurvature,
							int Nx, int Ny, int Nz);


int main(int argc, char **argv)
{
	int Nx,Ny,Nz;
	double Lx,Ly,Lz;

	Lx = Ly = Lz = 1.0;
	Nx = Ny = Nz = 102;

	int nspheres=1;
	// Set up the signed distance function for a sphere
	double *cx,*cy,*cz,*rad;
	cx = new double[nspheres];
	cy = new double[nspheres];
	cz = new double[nspheres];
	rad = new double[nspheres];

	// Sphere in the middle of the domain
	cx[0] = cy[0] = cz[0] = 0.5;

	DoubleArray SignDist(Nx,Ny,Nz);

	// Compute the signed distance function
	SignedDistance(SignDist.data,nspheres,cx,cy,cz,rad,Lx,Ly,Lz,Nx,Ny,Nz,0,0,0,1,1,1);



}
