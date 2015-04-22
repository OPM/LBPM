#include <iostream>
#include <math.h>
#include "pmmc.h"
#include "Domain.h"

using namespace std;

/*
 *  Compare the measured and analytical curvature for a sphere
 *
 */

int main(int argc, char **argv)
{
	int i,j,k;
	int Nx,Ny,Nz,N;
	double Lx,Ly,Lz;
	double fluid_isovalue=0.0;
	double solid_isovalue=0.0;

	Lx = Ly = Lz = 1.0;
	Nx = Ny = Nz = 102;
	N = Nx*Ny*Nz;

	//...........................................................................
	// Set up the cube list
	//...........................................................................
	int ncubes = (Nx-2)*(Ny-2)*(Nz-2);	// Exclude the "upper" halo
	IntArray cubeList(3,ncubes);
	pmmc_CubeListFromMesh(cubeList, ncubes, Nx, Ny, Nz);
	//...........................................................................

	//****************************************************************************************
	// Create the structures needed to carry out the PMMC algorithm
	//****************************************************************************************
	DTMutableList<Point> nw_pts(20);
	DTMutableList<Point> ns_pts(20);
	DTMutableList<Point> ws_pts(20);
	DTMutableList<Point> nws_pts(20);
	// initialize triangle lists for surfaces
	IntArray nw_tris(3,20);
	IntArray ns_tris(3,20);
	IntArray ws_tris(3,20);
	// initialize list for line segments
	IntArray nws_seg(2,20);

	DTMutableList<Point> tmp(20);

	DoubleArray wn_curvature(20);

	int n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
	int n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;

	// Initialize arrays for local solid surface
	DTMutableList<Point> local_sol_pts(20);
	int n_local_sol_pts = 0;
	IntArray local_sol_tris(3,18);
	int n_local_sol_tris;
	DoubleArray values(20);
	DTMutableList<Point> local_nws_pts(20);
	int n_local_nws_pts;
	//****************************************************************************************

	int nspheres=1;
	// Set up the signed distance function for a sphere
	double *cx,*cy,*cz,*rad;
	cx = new double[nspheres];
	cy = new double[nspheres];
	cz = new double[nspheres];
	rad = new double[nspheres];

	// Sphere in the middle of the domain
	cx[0] = cy[0] = cz[0] = 0.5;
	rad[0] = 0.3;

	DoubleArray SignDist(Nx,Ny,Nz);
	DoubleArray Phase(Nx,Ny,Nz);
	DoubleArray GaussCurvature(Nx,Ny,Nz);
	DoubleArray MeanCurvature(Nx,Ny,Nz);
	DoubleArray CubeValues(2,2,2);

	// Compute the signed distance function
	SignedDistance(Phase.get(),nspheres,cx,cy,cz,rad,Lx,Ly,Lz,Nx,Ny,Nz,0,0,0,1,1,1);

	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				SignDist(i,j,k) = 100.0;
				Phase(i,j,k) = sqrt((1.0*i-0.5*Nx)*(1.0*i-0.5*Nx)+(1.0*j-0.5*Ny)*(1.0*j-0.5*Ny)+(1.0*k-0.5*Nz)*(1.0*k-0.5*Nz))-0.3*Nx;
			}
		}
	}
	SignedDistance(SignDist.get(),0,cx,cy,cz,rad,Lx,Ly,Lz,Nx,Ny,Nz,0,0,0,1,1,1);

	pmmc_MeshCurvature(Phase, MeanCurvature, GaussCurvature, Nx, Ny, Nz);

	double wn_curvature_sum = 0.0;
	double wn_area_sum = 0.0;


	for (int c=0;c<ncubes;c++){

		// Get cube from the list
		i = cubeList(0,c);
		j = cubeList(1,c);
		k = cubeList(2,c);

		// Run PMMC
		n_local_sol_tris = 0;
		n_local_sol_pts = 0;
		n_local_nws_pts = 0;

		n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
		n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;

		// Construct the interfaces and common curve
		pmmc_ConstructLocalCube(SignDist, Phase, solid_isovalue, fluid_isovalue,
				nw_pts, nw_tris, values, ns_pts, ns_tris, ws_pts, ws_tris,
				local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
				n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
				n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
				i, j, k, Nx, Ny, Nz);

		// Interpolate the curvature onto the surface
		wn_curvature_sum += pmmc_CubeSurfaceInterpValue(CubeValues, MeanCurvature, nw_pts, nw_tris,
									wn_curvature, i, j, k, n_nw_pts, n_nw_tris);

		wn_area_sum += pmmc_CubeSurfaceArea(nw_pts, nw_tris, n_nw_tris);

	}

	printf("Mean Curvature Average =  %f, Analytical = %f \n", wn_curvature_sum/wn_area_sum, 2.0/rad[0]/101 );
	

	int toReturn = 0;
	if ( fabs(wn_curvature_sum/wn_area_sum -2.0/rad[0]/101)*rad[0]*101.0*0.5 > 0.01 ){
		toReturn = 1;
		printf("Mean curvature test error exceeds relative error tolerance \n ");
	}
	return toReturn;
}
