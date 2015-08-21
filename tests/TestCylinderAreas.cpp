#include <iostream>
#include <math.h>
#include "common/pmmc.h"
//#include "common/PointList.h"
//#include "common/Array.h"

#define RADIUS 15
#define HEIGHT 15.5
#define N 60
#define PI 3.14159

int main (int argc, char *argv[])
{
	
	//  printf("Radius = %s \n,"RADIUS);  
	int Nx,Ny,Nz;
	Nx = Ny = Nz = N;
	int i,j,k,p,r;
	
//	double *Solid; // cylinder
//	double *Phase; // region of the cylinder	
//	Solid = new double [N*N*N];
//	Phase = new double [N*N*N];
	DoubleArray SignDist(Nx,Ny,Nz);
	DoubleArray Phase(Nx,Ny,Nz);
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
	
	int n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0;
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
	//...........................................................................
	double Cx,Cy,Cz;
	double dist1,dist2;
	Cx = Cy = Cz = N*0.51;
	for (k=0; k<N; k++){
		for (j=0; j<N; j++){
			for (i=0; i<N; i++){
				dist1 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)) - RADIUS;
				dist2 = fabs(Cz-k)-HEIGHT;
				//	printf("distances = %f, %f \n",dist1,dist2);
				//Solid.data[k*Nx*Ny+j*Nx+i] = dist1;
				//Phase[k*Nx*Ny+j*Nx+i] = dist2;
				SignDist(i,j,k) = -dist1;
				Phase(i,j,k) = dist2;
			}
		}   
	}

	FILE *STRIS;
	STRIS = fopen("solid-triangles.out","w");
	
	FILE *WN_TRIS;
	WN_TRIS = fopen("wn-tris.out","w");
	
	FILE *NS_TRIS;
	NS_TRIS = fopen("ns-tris.out","w");
	
	FILE *WS_TRIS;
	WS_TRIS = fopen("ws-tris.out","w");
	
	FILE *WNS_PTS;
	WNS_PTS = fopen("wns-pts.out","w");
	
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
		
		n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0;
		n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;

		// Construct the interfaces and common curve
		pmmc_ConstructLocalCube(SignDist, Phase, solid_isovalue, fluid_isovalue,
				nw_pts, nw_tris, values, ns_pts, ns_tris, ws_pts, ws_tris,
				local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
				n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
				n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
				i, j, k, Nx, Ny, Nz);
		
		//*******************************************************************
		// Compute the Interfacial Areas, Common Line length for blob p
		// nw surface
		double temp;
		for (r=0;r<n_nw_tris;r++){
			A = nw_pts(nw_tris(0,r));
			B = nw_pts(nw_tris(1,r));
			C = nw_pts(nw_tris(2,r));
			// Compute length of sides (assume dx=dy=dz)
			s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
			s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
			s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
			s = 0.5*(s1+s2+s3);
			temp = s*(s-s1)*(s-s2)*(s-s3);
			if (temp > 0.0) awn += sqrt(temp);
			
		}
		for (r=0;r<n_ns_tris;r++){
			A = ns_pts(ns_tris(0,r));
			B = ns_pts(ns_tris(1,r));
			C = ns_pts(ns_tris(2,r));
			// Compute length of sides (assume dx=dy=dz)
			s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
			s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
			s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
			s = 0.5*(s1+s2+s3);
			//ans=ans+sqrt(s*(s-s1)*(s-s2)*(s-s3));
			temp = s*(s-s1)*(s-s2)*(s-s3);
			if (temp > 0.0) ans += sqrt(temp);
		}
		for (r=0;r<n_ws_tris;r++){
			A = ws_pts(ws_tris(0,r));
			B = ws_pts(ws_tris(1,r));
			C = ws_pts(ws_tris(2,r));
			// Compute length of sides (assume dx=dy=dz)
			s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
			s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
			s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
			s = 0.5*(s1+s2+s3);
			//aws=aws+sqrt(s*(s-s1)*(s-s2)*(s-s3));
			temp = s*(s-s1)*(s-s2)*(s-s3);
			if (temp > 0.0) aws += sqrt(temp);
		}
		for (r=0;r<n_local_sol_tris;r++){
			A = local_sol_pts(local_sol_tris(0,r));
			B = local_sol_pts(local_sol_tris(1,r));
			C = local_sol_pts(local_sol_tris(2,r));
			// Compute length of sides (assume dx=dy=dz)
			s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
			s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
			s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
			s = 0.5*(s1+s2+s3);
			//aws=aws+sqrt(s*(s-s1)*(s-s2)*(s-s3));
			temp = s*(s-s1)*(s-s2)*(s-s3);
			if (temp > 0.0) As += sqrt(temp);
		}
		for (p=0; p < n_local_nws_pts-1; p++){
			// Extract the line segment
			A = local_nws_pts(p);
			B = local_nws_pts(p+1);
			// Compute the length of the segment
			s = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
			// Add the length to the common line 
			lwns += s;
		}
		//.......................................................................................
		// Write the triangle lists to text file
		for (r=0;r<n_nw_tris;r++){
			A = nw_pts(nw_tris(0,r));
			B = nw_pts(nw_tris(1,r));
			C = nw_pts(nw_tris(2,r));
			fprintf(WN_TRIS,"%f %f %f %f %f %f %f %f %f \n",A.x,A.y,A.z,B.x,B.y,B.z,C.x,C.y,C.z);
		}		
		for (r=0;r<n_ws_tris;r++){
			A = ws_pts(ws_tris(0,r));
			B = ws_pts(ws_tris(1,r));
			C = ws_pts(ws_tris(2,r));
			fprintf(WS_TRIS,"%f %f %f %f %f %f %f %f %f \n",A.x,A.y,A.z,B.x,B.y,B.z,C.x,C.y,C.z);
		}
		for (r=0;r<n_ns_tris;r++){
			A = ns_pts(ns_tris(0,r));
			B = ns_pts(ns_tris(1,r));
			C = ns_pts(ns_tris(2,r));
			fprintf(NS_TRIS,"%f %f %f %f %f %f %f %f %f \n",A.x,A.y,A.z,B.x,B.y,B.z,C.x,C.y,C.z);
		}
		for (p=0; p < n_nws_pts; p++){
			P = nws_pts(p);
			fprintf(WNS_PTS,"%f %f %f \n",P.x, P.y, P.z);
		}
		
		//.......................................................................................
		
		for (r=0;r<n_local_sol_tris;r++){
			A = local_sol_pts(local_sol_tris(0,r));
			B = local_sol_pts(local_sol_tris(1,r));
			C = local_sol_pts(local_sol_tris(2,r));
			fprintf(STRIS,"%f %f %f %f %f %f %f %f %f \n",A.x,A.y,A.z,B.x,B.y,B.z,C.x,C.y,C.z);
		}
		//*******************************************************************
		// Reset the triangle counts to zero
		n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0;
		n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
		//	n_ns_tris_beg = 0;//n_ns_tris;
		//	n_ws_tris_beg = 0;//n_ws_tris;
		//	n_nws_seg_beg = n_nws_seg;
		//*******************************************************************
	}
	fclose(WN_TRIS);
	fclose(NS_TRIS);
	fclose(WS_TRIS);
	fclose(WNS_PTS);
	fclose(STRIS);

	printf("-------------------------------- \n");
	printf("NWP volume = %f \n", nwp_volume);
	printf("Area wn = %f, Analytical = %f \n", awn,2*PI*RADIUS*RADIUS);
	printf("Area ns = %f, Analytical = %f \n", ans, 2*PI*RADIUS*(N-2)-4*PI*RADIUS*HEIGHT);
	printf("Area ws = %f, Analytical = %f \n", aws, 4*PI*RADIUS*HEIGHT);
	printf("Area s = %f, Analytical = %f \n", As, 2*PI*RADIUS*(N-2));
	printf("Length wns = %f, Analytical = %f \n", lwns, 4*PI*RADIUS);
	printf("-------------------------------- \n");	
	//.........................................................................
	
	int toReturn = 0;
	if (fabs(awn - 2*PI*RADIUS*RADIUS)/(2*PI*RADIUS*RADIUS) > 0.02){	
		toReturn += 1;
		printf("TestCylinderArea.cpp: error tolerance exceeded for wn area \n");
	}
	if (fabs(ans - (2*PI*RADIUS*(N-2)-4*PI*RADIUS*HEIGHT))/(2*PI*RADIUS*(N-2)-4*PI*RADIUS*HEIGHT)> 0.02 ){
		toReturn += 2;
		printf("TestCylinderArea.cpp: error tolerance exceeded for ns area \n");
	}
	if (fabs(aws - 4*PI*RADIUS*HEIGHT)/(4*PI*RADIUS*HEIGHT) > 0.02 ){
		toReturn += 3;
		printf("TestCylinderArea.cpp: error tolerance exceeded for ws area \n");
	}
	if (fabs(As - 2*PI*RADIUS*(N-2))/(2*PI*RADIUS*(N-2)) > 0.02 ){
		toReturn += 4;
		printf("TestCylinderArea.cpp: error tolerance exceeded for solid area \n");
	}
	if (fabs(lwns - 4*PI*RADIUS)/(4*PI*RADIUS) > 0.02 ){
		toReturn += 5;
		printf("TestCylinderArea.cpp: error tolerance exceeded for common curve length \n");
	}

	return toReturn;
}
