#include <iostream>
#include <math.h>
#include "pmmc.h"
//#include "PointList.h"
//#include "Array.h"

#define RADIUS 15
#define CAPRAD 20
#define HEIGHT 15.5
#define N 60
#define SPEED -1
#define PI 3.14159

int main (int argc, char *argv[])
{
	
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
	DoubleArray Phase_x(Nx,Ny,Nz);
	DoubleArray Phase_y(Nx,Ny,Nz);
	DoubleArray Phase_z(Nx,Ny,Nz);
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
	double efawns,Jwn;
	double KNwns,KGwns;
	double As;
	double dEs,dAwn,dAns;			 // Global surface energy (calculated by rank=0)
	double awn_global,ans_global,aws_global,lwns_global,nwp_volume_global;	
	double As_global;
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
	
	DTMutableList<Point> tmp(20);
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
	
	DoubleArray CubeValues(2,2,2);
	DoubleArray ContactAngle(20);
	DoubleArray KGwns_values(20);
	DoubleArray KNwns_values(20);
	DoubleArray wn_curvature(20);
	DoubleArray InterfaceSpeed(20);
	DoubleArray NormalVector(60);
	DoubleArray vawn(6);
	DoubleArray vawns(3);
	
	int c;
	//...........................................................................
	int ncubes = (Nx-2)*(Ny-2)*(Nz-2);	// Exclude the "upper" halo
	IntArray cubeList(3,ncubes);
	pmmc_CubeListFromMesh(cubeList, ncubes, Nx, Ny, Nz);
	//...........................................................................
	double Cx,Cy,Cz;
	double dist1,dist2;
	// Extra copies of phase indicator needed to compute time derivatives on CPU
	DoubleArray Phase_tminus(Nx,Ny,Nz);
	DoubleArray Phase_tplus(Nx,Ny,Nz);
	DoubleArray dPdt(Nx,Ny,Nz);
	Cx = Cy = Cz = N*0.5;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;
				dist2 = fabs(Cz-k)-HEIGHT;

				Phase_tminus(i,j,k) = dist2;
			}
		} 
	}
	Cz += SPEED;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				
				dist1 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)) - RADIUS;
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;
				dist2 = fabs(Cz-k)-HEIGHT;

				SignDist(i,j,k) = -dist1;
				Phase(i,j,k) = dist2;
			}
		}   
	}
	Cz += SPEED;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;
				dist2 = fabs(Cz-k)-HEIGHT;

				Phase_tplus(i,j,k) = dist2;
			}
		}   
	}
	
	//...........................................................................
	// Calculate the time derivative of the phase indicator field
	for (int n=0; n<Nx*Ny*Nz; n++)	dPdt(n) = 0.5*(Phase_tplus(n) - Phase_tminus(n));
	
	pmmc_MeshGradient(Phase,Phase_x,Phase_y,Phase_z,Nx,Ny,Nz);
	pmmc_MeshGradient(SignDist,Sx,Sy,Sz,Nx,Ny,Nz);
	
	double norm;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				norm = Phase_x(i,j,k)*Phase_x(i,j,k)+Phase_y(i,j,k)*Phase_y(i,j,k)+Phase_z(i,j,k)*Phase_z(i,j,k);
			}
		}
	}
	
	// End of the loop to set the values
	awn = aws = ans = lwns = 0.0;
	nwp_volume = 0.0;
	As = 0.0;
	
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

		efawns += pmmc_CubeContactAngle(CubeValues,ContactAngle,Phase_x,Phase_y,Phase_z,Sx,Sy,Sz,local_nws_pts,i,j,k,n_local_nws_pts);
		
		Jwn += pmmc_CubeSurfaceInterpValue(CubeValues, MeanCurvature, nw_pts, nw_tris,
									wn_curvature, i, j, k, n_nw_pts, n_nw_tris);
		
		pmmc_InterfaceSpeed(dPdt, Phase_x, Phase_y, Phase_z, CubeValues, nw_pts, nw_tris,
							NormalVector, InterfaceSpeed, vawn, i, j, k, n_nw_pts, n_nw_tris);
		
		pmmc_CommonCurveSpeed(CubeValues, dPdt, vawns, Phase_x,Phase_y,Phase_z,Sx,Sy,Sz,
				local_nws_pts,i,j,k,n_local_nws_pts);

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
	for (i=0;i<6;i++)	vawn(i) /= awn;
	for (i=0;i<3;i++)	vawns(i) /= lwns;
	
	printf("-------------------------------- \n");
	printf("NWP volume = %f \n", nwp_volume);
	printf("Area wn = %f, Analytical = %f \n", awn,2*PI*RADIUS*RADIUS);
	printf("Area ns = %f, Analytical = %f \n", ans, 2*PI*RADIUS*(N-2)-4*PI*RADIUS*HEIGHT);
	printf("Area ws = %f, Analytical = %f \n", aws, 4*PI*RADIUS*HEIGHT);
	printf("Area s = %f, Analytical = %f \n", As, 2*PI*RADIUS*(N-2));
	printf("Geodesic curvature (wns) = %f, Analytical = %f \n", KGwns, 0.0);
	printf("Normal curvature (wns) = %f, Analytical = %f \n", KNwns, 1.0/RADIUS);
	printf("Length wns = %f, Analytical = %f \n", lwns, 4*PI*RADIUS);
//	printf("Cos(theta_wns) = %f, Analytical = %f \n",efawns/lwns,1.0*RADIUS/CAPRAD);
	printf("Interface Velocity = %f,%f,%f \n",vawn(0),vawn(1),vawn(2));
	printf("Common Curve Velocity = %f,%f,%f \n",vawns(0),vawns(1),vawns(2));
	printf("-------------------------------- \n");	
	//.........................................................................	
	
	int toReturn = 0;
	return toReturn;
}
