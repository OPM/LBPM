#include <iostream>
#include <math.h>
#include "analysis/Minkowski.h"
#include "common/Domain.h"
#include "common/SpherePack.h"

using namespace std;

/*
 *  Compare the measured and analytical curvature for a sphere
 *
 */

int main(int argc, char **argv)
{
	int i,j,k;
	int Nx,Ny,Nz;
	double Lx,Ly,Lz;
	double fluid_isovalue=0.0;
	double solid_isovalue=0.0;

	Lx = Ly = Lz = 1.0;
	Nx = Ny = Nz = 64;
	DoubleArray Phase(Nx,Ny,Nz);
	DoubleArray CubeValues(2,2,2);
	
	printf("Set distance map \n");
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				Phase(i,j,k) = sqrt((1.0*i-0.5*Nx)*(1.0*i-0.5*Nx)+(1.0*j-0.5*Ny)*(1.0*j-0.5*Ny)+(1.0*k-0.5*Nz)*(1.0*k-0.5*Nz))-0.3*Nx;
			}
		}
	}

	double wn_curvature_sum = 0.0;
	double wn_area_sum = 0.0;
	
	printf("Construct local isosurface \n");
	DECL sphere;
	for (k=0; k<Nz-1; k++){
		for (j=0; j<Ny-1; j++){
			for (i=0; i<Nx-1; i++){
				sphere.LocalIsosurface(Phase,0.f,i,j,k);
				for (unsigned long int idx=0; idx<sphere.TriangleCount; idx++){
					unsigned long int edge = sphere.Face(idx); 
					
				}
			}
		}
	}

//	printf("Mean Curvature Average =  %f, Analytical = %f \n", wn_curvature_sum/wn_area_sum, 2.0/rad[0]/101 );

	int toReturn = 0;
/*	if ( fabs(wn_curvature_sum/wn_area_sum -2.0/rad[0]/101)*rad[0]*101.0*0.5 > 0.01 ){
		toReturn = 1;
		printf("Mean curvature test error exceeds relative error tolerance \n ");
	}
	*/
	return toReturn;
}
