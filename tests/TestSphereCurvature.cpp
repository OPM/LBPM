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
	
	printf("Construct local isosurface \n");
	DECL sphere;
	Point V1,V2,V3;
	unsigned long int e1,e2,e3;
	double s,s1,s2,s3;
	double area = 0.f;
	for (k=0; k<Nz-1; k++){
		for (j=0; j<Ny-1; j++){
			for (i=0; i<Nx-1; i++){
				sphere.LocalIsosurface(Phase,0.f,i,j,k);
				for (unsigned long int idx=0; idx<sphere.TriangleCount; idx++){
					e1 = sphere.Face(idx); 
					e2 = sphere.halfedge.next(e1);
					e3 = sphere.halfedge.next(e2);
					V1 = sphere.vertex.coords(sphere.halfedge.v1(e1));
					V2 = sphere.vertex.coords(sphere.halfedge.v1(e2));
					V3 = sphere.vertex.coords(sphere.halfedge.v1(e3));
					s1 = sqrt((V1.x-V2.x)*(V1.x-V2.x)+(V1.y-V2.y)*(V1.y-V2.y)+(V1.z-V2.z)*(V1.z-V2.z));
					s2 = sqrt((V1.x-V3.x)*(V1.x-V3.x)+(V1.y-V3.y)*(V1.y-V3.y)+(V1.z-V3.z)*(V1.z-V3.z));
					s3 = sqrt((V2.x-V3.x)*(V2.x-V3.x)+(V2.y-V3.y)*(V2.y-V3.y)+(V2.z-V3.z)*(V2.z-V3.z));
					s = 0.5*(s1+s2+s3);
					area += sqrt(s*(s-s1)*(s-s2)*(s-s3));
				}
			}
		}
	}

	printf("Surface area  = %f (analytical = %f) \n", area,4*3.14159*0.3*0.3*double(Nx*Nx));
	//	printf("Mean Curvature Average =  %f, Analytical = %f \n", wn_curvature_sum/wn_area_sum, 2.0/rad[0]/101 );

	int toReturn = 0;
/*	if ( fabs(wn_curvature_sum/wn_area_sum -2.0/rad[0]/101)*rad[0]*101.0*0.5 > 0.01 ){
		toReturn = 1;
		printf("Mean curvature test error exceeds relative error tolerance \n ");
	}
	*/
	return toReturn;
}
