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
	Point P1,P2,P3;
	unsigned long int e1,e2,e3;
	double s,s1,s2,s3;
	double area = 0.f;
	double Xi = 0.f;
	double Vx,Vy,Vz,Wx,Wy,Wz,Nx,Ny,Nz,norm;
	for (k=0; k<Nz-1; k++){
		for (j=0; j<Ny-1; j++){
			for (i=0; i<Nx-1; i++){
				sphere.LocalIsosurface(Phase,0.f,i,j,k);
				for (unsigned long int idx=0; idx<sphere.TriangleCount; idx++){
					e1 = sphere.Face(idx); 
					e2 = sphere.halfedge.next(e1);
					e3 = sphere.halfedge.next(e2);
					P1 = sphere.vertex.coords(sphere.halfedge.v1(e1));
					P2 = sphere.vertex.coords(sphere.halfedge.v1(e2));
					P3 = sphere.vertex.coords(sphere.halfedge.v1(e3));
					// compute the area
					s1 = sqrt((P1.x-P2.x)*(P1.x-P2.x)+(P1.y-P2.y)*(P1.y-P2.y)+(P1.z-P2.z)*(P1.z-P2.z));
					s2 = sqrt((P1.x-P3.x)*(P1.x-P3.x)+(P1.y-P3.y)*(P1.y-P3.y)+(P1.z-P3.z)*(P1.z-P3.z));
					s3 = sqrt((P2.x-P3.x)*(P2.x-P3.x)+(P2.y-P3.y)*(P2.y-P3.y)+(P2.z-P3.z)*(P2.z-P3.z));
					s = 0.5*(s1+s2+s3);
					area += sqrt(s*(s-s1)*(s-s2)*(s-s3));
					// compute the normal vector
					Vx=P2.x-P1.x;
					Vy=P2.y-P1.y;
					Vz=P2.z-P1.z;
					Wx=P3.x-P2.x;
					Wy=P3.y-P2.y;
					Wz=P3.z-P2.z;
					Nx = Vy*Wz-Vz*Wy;
					Ny = Vz*Wx-Vx*Wz;
					Nz = Vx*Wy-Vy*Wx;
					norm = 1.f/sqrt(Nx*Nx+Ny*Ny+Nz*Nz);
					Nx *= norm;
					Ny *= norm;
					Nz *= norm;
					// Euler characteristic (half edge rule: one face - 0.5*(three edges))
					Xi -= 0.5;
				}
				// Euler characteristic -- each vertex shared by four cubes
				Xi += 0.25*double(sphere.VertexCount);
			}
		}
	}

	printf("Surface area  = %f (analytical = %f) \n", area,4*3.14159*0.3*0.3*double(Nx*Nx));
	printf("Euler characteristic  = %f (analytical = 2.0) \n",Xi);

	//	printf("Mean Curvature Average =  %f, Analytical = %f \n", wn_curvature_sum/wn_area_sum, 2.0/rad[0]/101 );

	int toReturn = 0;
/*	if ( fabs(wn_curvature_sum/wn_area_sum -2.0/rad[0]/101)*rad[0]*101.0*0.5 > 0.01 ){
		toReturn = 1;
		printf("Mean curvature test error exceeds relative error tolerance \n ");
	}
	*/
	return toReturn;
}
