#include <iostream>
#include <math.h>
#include "pmmc.h"
//#include "PointList.h"
//#include "Array.h"

#define RADIUS 15
#define CAPRAD 20
#define HEIGHT 15.5
#define N 60
#define PI 3.14159

int main (int argc, char *argv[])
{
	
	printf("Unit test for analysis framework \n");

	//  printf("Radius = %s \n,"RADIUS);  
	int SIZE = N*N*N;
	int Nx,Ny,Nz;
	Nx = Ny = Nz = N;
	int i,j,k,p,q,r;
	
//	double *Solid; // cylinder
//	double *Phase; // region of the cylinder	
//	Solid = new double [SIZE];
//	Phase = new double [SIZE];
	DoubleArray SignDist(Nx,Ny,Nz);
	DoubleArray Phase(Nx,Ny,Nz);
	DoubleArray Fx(Nx,Ny,Nz);
	DoubleArray Fy(Nx,Ny,Nz);
	DoubleArray Fz(Nx,Ny,Nz);
	DoubleArray Sx(Nx,Ny,Nz);
	DoubleArray Sy(Nx,Ny,Nz);
	DoubleArray Sz(Nx,Ny,Nz);
	DoubleArray GaussCurvature(Nx,Ny,Nz);
	DoubleArray MeanCurvature(Nx,Ny,Nz);
	double fluid_isovalue = 0.0;
	double solid_isovalue = 0.0;
	
	/* ****************************************************************
	 VARIABLES FOR THE PMMC ALGORITHM
	 ****************************************************************** */
	//...........................................................................
	// Averaging variables
	//...........................................................................
	double awn,ans,aws,lwns,nwp_volume;
	double As;
	double efawns,Jwn;
	double KNwns,KGwns;
	DoubleArray Gwns(6);

	//	bool add=1;			// Set to false if any corners contain nw-phase ( F > fluid_isovalue)
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	//	int count_in=0,count_out=0;
	//	int nodx,nody,nodz;
	// initialize lists for vertices for surfaces, common line
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
	
	DoubleArray CubeValues(2,2,2);
	DTMutableList<Point> tmp(20);
	DoubleArray ContactAngle(20);
	DoubleArray KGwns_values(20);
	DoubleArray KNwns_values(20);
	DoubleArray wn_curvature(20);
	//	IntArray store;
	
	int n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
	int n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
	
	double s,s1,s2,s3;		// Triangle sides (lengths)
	Point A,B,C,P;
	//	double area;
	
	// Initialize arrays for local solid surface
	DTMutableList<Point> local_sol_pts(20);
	int n_local_sol_pts = 0;
	IntArray local_sol_tris(3,18);
	int n_local_sol_tris;
	DoubleArray values(20);
	DTMutableList<Point> local_nws_pts(20);
	int n_local_nws_pts;
		
	int c;
	//...........................................................................
	int ncubes = (Nx-2)*(Ny-2)*(Nz-2);	// Exclude the "upper" halo
	IntArray cubeList(3,ncubes);
	pmmc_CubeListFromMesh(cubeList, ncubes, Nx, Ny, Nz);
	int nc=1;
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				cubeList(0,nc) = i;
				cubeList(1,nc) = j;
				cubeList(2,nc) = k;
				nc++;
			}
		}
	}
	ncubes = nc;

	//...........................................................................
	double Cx,Cy,Cz;
	double dist1,dist2;
	Cx = Cy = Cz = N*0.51;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				
				dist1 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)) - RADIUS;
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;

				SignDist(i,j,k) = -dist1;
				Phase(i,j,k) = dist2;

//				if (k<Cz && 1.0*(Cz-k) < Phase(i,j,k) )  Phase(i,j,k) = 1.0*(Cz-k);
			}
		}   
	}

	pmmc_MeshGradient(Phase,Fx,Fy,Fz,Nx,Ny,Nz);
	pmmc_MeshGradient(SignDist,Sx,Sy,Sz,Nx,Ny,Nz);
	pmmc_MeshCurvature(Phase, MeanCurvature, GaussCurvature, Nx, Ny, Nz);
	
	// End of the loop to set the values
	awn = aws = ans = lwns = 0.0;
	nwp_volume = 0.0;
	As = 0.0;
	for (i=0; i<6; i++) Gwns(i) = 0.0;
	
	for (c=0;c<ncubes;c++){
		// Get cube from the list
		i = cubeList(0,c);
		j = cubeList(1,c);
		k = cubeList(2,c);
		
		for (p=0;p<8;p++){
			if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 
				&&  SignDist(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){
				nwp_volume += 0.125;
			}
		}
		
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

		efawns += pmmc_CubeContactAngle(CubeValues,ContactAngle,Fx,Fy,Fz,Sx,Sy,Sz,local_nws_pts,i,j,k,n_local_nws_pts);

		Jwn += pmmc_CubeSurfaceInterpValue(CubeValues, MeanCurvature, nw_pts, nw_tris,
									wn_curvature, i, j, k, n_nw_pts, n_nw_tris);

		pmmc_CurveOrientation(Gwns, nws_pts, n_nws_pts, i,j,k);

		pmmc_CurveCurvature(Phase, SignDist, KNwns_values, KGwns_values, KNwns, KGwns, nws_pts, n_nws_pts, i, j, k);

	//	if (n_nw_pts>0) printf("speed %f \n",InterfaceSpeed(0));
		
		//*******************************************************************
		// Compute the Interfacial Areas, Common Line length
		awn += pmmc_CubeSurfaceArea(nw_pts,nw_tris,n_nw_tris);
		ans += pmmc_CubeSurfaceArea(ns_pts,ns_tris,n_ns_tris);
		aws += pmmc_CubeSurfaceArea(ws_pts,ws_tris,n_ws_tris);
		As += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);
		lwns +=  pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);

	}
	KGwns /= lwns;
	KNwns /= lwns;
	Jwn /= awn;
	efawns /= lwns;
	for (i=0; i<6; i++) Gwns(i) /= lwns;

	printf("Analysis complete. \n");

	double CAPHEIGHT = CAPRAD-sqrt(CAPRAD*CAPRAD-RADIUS*RADIUS); // height of the sphereical cap
	double analytical,RelError;
	printf("Height of sphereical cap = %f \n",CAPHEIGHT);
	printf("-------------------------------- \n");
	printf("NWP volume = %f \n", nwp_volume);
	analytical = 2*PI*(CAPHEIGHT*CAPHEIGHT+RADIUS*RADIUS);
	RelError = fabs(awn-analytical)/analytical;
	printf("Area wn = %f, Analytical = %f, Rel. Error = %f \n", awn,analytical,RelError);
	analytical = 2*PI*RADIUS*(N-2)-4*PI*RADIUS*(CAPRAD-CAPHEIGHT);
	RelError = fabs(ans-analytical)/analytical;
	printf("Area ns = %f, Analytical = %f, Rel. Error = %f \n", ans, analytical,RelError);
	analytical = 4*PI*RADIUS*(CAPRAD-CAPHEIGHT);
	RelError = fabs(aws-analytical)/analytical;
	printf("Area ws = %f, Analytical = %f, Rel. Error = %f \n", aws, analytical,RelError);
	analytical = 2*PI*RADIUS*(N-2);
	RelError = fabs(As-analytical)/analytical;
	printf("Area s = %f, Analytical = %f, Rel. Error = %f \n", As, analytical,RelError);
	analytical = 4*PI*RADIUS;
	RelError = fabs(lwns-analytical)/analytical;
	printf("Length wns = %f, Analytical = %f, Rel. Error = %f  \n", lwns, analytical,RelError);
	analytical = 1.0*RADIUS/CAPRAD;
	RelError = fabs(efawns-analytical)/analytical;
	printf("Cos(theta_wns) = %f, Analytical = %f, Rel. Error = %f  \n",efawns,analytical,RelError);
	analytical = 2.0/CAPRAD;
	RelError = fabs(Jwn-analytical)/analytical;
	printf("Mean curvature (wn) = %f, Analytical = %f, Rel. Error = %f  \n", Jwn, analytical,RelError);
	printf("Geodesic curvature (wns) = %f, Analytical, Rel. Error = %f  = %f \n", KGwns, 0.0, 0.0);
	analytical = 1.0/RADIUS;
	RelError = fabs(KNwns-analytical)/analytical;
	printf("Normal curvature (wns) = %f, Analytical = %f, Rel. Error = %f  \n", KNwns, analytical,RelError);
	printf("-------------------------------- \n");	
	//.........................................................................
	printf("Gwns=%f,%f,%f,%f,%f,%f \n",Gwns(0),Gwns(1),Gwns(2),Gwns(3),Gwns(4),Gwns(5));
	
	int toReturn = 0;
	if (fabs(efawns - 1.0*RADIUS/CAPRAD)/(1.0*RADIUS/CAPRAD) > 0.01){
		toReturn = 1;
		printf("tests/TestContactAngle.cpp: exceeded error tolerance for the contact angle \n");
	}
	else{
		printf("Passed test: contact angle");
	}

	if (fabs(ans - 2*PI*RADIUS*(N-2)-4*PI*RADIUS*(CAPRAD-CAPHEIGHT))/(4*PI*RADIUS*(CAPRAD-CAPHEIGHT)) > 0.01 ){
		printf("tests/TestContactAngle.cpp: exceeded error tolerance for ns area \n");
	}
	else{
		printf("Passed test: ans");
	}
	if (fabs(awn-2*PI*(CAPHEIGHT*CAPHEIGHT+RADIUS*RADIUS))/(2*PI*(CAPHEIGHT*CAPHEIGHT+RADIUS*RADIUS)) > 0.01){
		printf("tests/TestContactAngle.cpp: exceeded error tolerance for wn area \n");
	}
	else{
		printf("Passed test: awn");
	}


	return toReturn;

}
