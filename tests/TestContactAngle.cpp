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
	double fluid_isovalue = 0.0;
	double solid_isovalue = 0.0;

	
	/* ****************************************************************
	 VARIABLES FOR THE PMMC ALGORITHM
	 ****************************************************************** */
	//...........................................................................
	// Averaging variables
	//...........................................................................
	double awn,ans,aws,lwns,nwp_volume;
	double efawns;
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
	
	int c;
	//...........................................................................
	int ncubes = (Nx-2)*(Ny-2)*(Nz-2);	// Exclude the "upper" halo
	IntArray cubeList(3,ncubes);
	pmmc_CubeListFromMesh(cubeList, ncubes, Nx, Ny, Nz);
	//...........................................................................
	double Cx,Cy,Cz;
	double dist1,dist2;
	Cx = Cy = Cz = N*0.51;
	for (k=0; k<N; k++){
		for (j=0; j<N; j++){
			for (i=0; i<N; i++){
				
				dist1 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)) - RADIUS;
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;

				SignDist(i,j,k) = -dist1;
				Phase(i,j,k) = dist2;
			}
		}   
	}

	pmmc_MeshGradient(Phase,Fx,Fy,Fz,Nx,Ny,Nz);
	pmmc_MeshGradient(SignDist,Sx,Sy,Sz,Nx,Ny,Nz);
	
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

		efawns += pmmc_CubeContactAngle(CubeValues,ContactAngle,Fx,Fy,Fz,Sx,Sy,Sz,local_nws_pts,i,j,k,n_local_nws_pts);
		
		//*******************************************************************
		// Compute the Interfacial Areas, Common Line length
		awn += pmmc_CubeSurfaceArea(nw_pts,nw_tris,n_nw_tris);
		ans += pmmc_CubeSurfaceArea(ns_pts,ns_tris,n_ns_tris);
		aws += pmmc_CubeSurfaceArea(ws_pts,ws_tris,n_ws_tris);
		As += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);
		lwns +=  pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
	}

	printf("-------------------------------- \n");
	printf("NWP volume = %f \n", nwp_volume);
	printf("Area wn = %f, Analytical = %f \n", awn,2*PI*RADIUS*RADIUS);
	printf("Area ns = %f, Analytical = %f \n", ans, 2*PI*RADIUS*(N-2)-4*PI*RADIUS*HEIGHT);
	printf("Area ws = %f, Analytical = %f \n", aws, 4*PI*RADIUS*HEIGHT);
	printf("Area s = %f, Analytical = %f \n", As, 2*PI*RADIUS*(N-2));
	printf("Length wns = %f, Analytical = %f \n", lwns, 4*PI*RADIUS);
	printf("Cos(theta_wns) = %f, Analytical = %f \n",efawns/lwns,1.0*RADIUS/CAPRAD);
	printf("-------------------------------- \n");	
	//.........................................................................
	
/*	FILE *PHASE;
	PHASE = fopen("Phase.in","wb");
	fwrite(Phase,8,SIZE,PHASE);
	fclose(PHASE);
	
	FILE *SOLID;
	SOLID = fopen("Distance.in","wb");
	fwrite(Solid,8,SIZE,SOLID);
	fclose(SOLID);
*/
}
