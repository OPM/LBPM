#include <iostream>
#include <math.h>
#include "analysis/pmmc.h"
//#include "common/PointList.h"
//#include "common/Array.h"

#define RADIUS 15
#define HEIGHT 15.5
#define N 60
#define PI 3.14159

int main (int argc, char **argv)
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
	
	int c;
	int newton_steps = 0;
	//...........................................................................
	int ncubes = (Nx-2)*(Ny-2)*(Nz-2);	// Exclude the "upper" halo
	IntArray cubeList(3,ncubes);
	int nc=0;
	//...........................................................................
	// Set up the cube list (very regular in this case due to lack of blob-ID)
	for (k=0; k<Nz-2; k++){
		for (j=0; j<Ny-2; j++){
			for (i=0; i<Nx-2; i++){
				cubeList(0,nc) = i;
				cubeList(1,nc) = j;
				cubeList(2,nc) = k;
				nc++;
			}
		}
	}
	
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
		
		n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
		n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
		
		// if there is a solid phase interface in the grid cell
		if (Interface(SignDist,solid_isovalue,i,j,k) == 1){
			
			/////////////////////////////////////////
			/// CONSTRUCT THE LOCAL SOLID SURFACE ///
			/////////////////////////////////////////
			
			// find the local solid surface
			//	SOL_SURF(SignDist,0.0,Phase,fluid_isovalue,i,j,k, Nx,Ny,Nz,local_sol_pts,n_local_sol_pts,
			//			 local_sol_tris,n_local_sol_tris,values);
			
			// find the local solid surface using the regular Marching Cubes algorithm
			SolidMarchingCubes(SignDist,0.0,Phase,fluid_isovalue,i,j,k,Nx,Ny,Nz,local_sol_pts,n_local_sol_pts,
							   local_sol_tris,n_local_sol_tris,values);
			
			/////////////////////////////////////////
			//////// TRIM THE SOLID SURFACE /////////
			/////////////////////////////////////////
			/*					TRIM(local_sol_pts, n_local_sol_pts, fluid_isovalue,local_sol_tris, n_local_sol_tris,
			 ns_pts, n_ns_pts, ns_tris, n_ns_tris, ws_pts, n_ws_pts,
			 ws_tris, n_ws_tris, values, local_nws_pts, n_local_nws_pts,
			 Phase, SignDist, i, j, k, newton_steps);
			 */					
			TRIM(local_sol_pts, n_local_sol_pts, fluid_isovalue,local_sol_tris, n_local_sol_tris,
				 ns_pts, n_ns_pts, ns_tris, n_ns_tris, ws_pts, n_ws_pts,
				 ws_tris, n_ws_tris, values, local_nws_pts, n_local_nws_pts);
			
			/////////////////////////////////////////
			//////// WRITE COMMON LINE POINTS ///////
			////////      TO MAIN ARRAYS      ///////
			/////////////////////////////////////////
			// SORT THE LOCAL COMMON LINE POINTS
			/////////////////////////////////////////
			// Make sure the first common line point is on a face
			// Common curve points are located pairwise and must
			// be searched and rearranged accordingly
			for (p=0; p<n_local_nws_pts-1; p++){
				P = local_nws_pts(p);
				if ( P.x == 1.0*i || P.x ==1.0*(i+1)|| 
					P.y == 1.0*j || P.y == 1.0*(j+1) || 
					P.z == 1.0*k || P.z == 1.0*(k+1) ){
					if (p%2 == 0){
						// even points
						// Swap the pair of points
						local_nws_pts(p) = local_nws_pts(0);
						local_nws_pts(0) = P;
						P = local_nws_pts(p+1);
						local_nws_pts(p+1) = local_nws_pts(1);
						local_nws_pts(1) = P;
						p = n_local_nws_pts; 
						
					}
					else{
						// odd points - flip the order
						local_nws_pts(p) = local_nws_pts(p-1);
						local_nws_pts(p-1) = P;
						p-=2;
						
					}
					// guarantee exit from the loop
				}
 			}
			// Two common curve points per triangle
			// 0-(1=2)-(3=4)-...
			for (p=1; p<n_local_nws_pts-1; p+=2){
				A = local_nws_pts(p);
				for (q=p+1; q<n_local_nws_pts; q++){
					B = local_nws_pts(q);
					if ( A.x == B.x && A.y == B.y && A.z == B.z){
						if (q%2 == 0){
							// even points
							// Swap the pair of points
							local_nws_pts(q) = local_nws_pts(p+1);
							local_nws_pts(p+1) = B;
							B = local_nws_pts(q+1);
							local_nws_pts(q+1) = local_nws_pts(p+2);
							local_nws_pts(p+2) = B;
							q = n_local_nws_pts;
							
						}
						else{
							// odd points - flip the order
							local_nws_pts(q) = local_nws_pts(q-1);
							local_nws_pts(q-1) = B;
							q-=2;
						}
					}
				}
			}
			map = n_nws_pts = 0;			
			nws_pts(n_nws_pts++) = local_nws_pts(0);
			for (p=2; p < n_local_nws_pts; p+=2){
				nws_pts(n_nws_pts++) = local_nws_pts(p);
				
			}
			nws_pts(n_nws_pts++) = local_nws_pts(n_local_nws_pts-1);
			
			for (q=0; q < n_nws_pts-1; q++){
				nws_seg(0,n_nws_seg) = map+q;
				nws_seg(1,n_nws_seg) = map+q+1;
				n_nws_seg++;
			}
			// End of the common line sorting algorithm
			/////////////////////////////////////////
			
			
			/////////////////////////////////////////
			////// CONSTRUCT THE nw SURFACE /////////
			/////////////////////////////////////////
			if ( n_local_nws_pts > 0){
				n_nw_tris =0;
				EDGE(Phase, fluid_isovalue, SignDist, i,j,k, Nx, Ny, Nz, nw_pts, n_nw_pts, nw_tris, n_nw_tris,
					 nws_pts, n_nws_pts);
			}
			else {
				MC(Phase, fluid_isovalue, SignDist, i,j,k, nw_pts, n_nw_pts, nw_tris, n_nw_tris);
			}
		}
		
		/////////////////////////////////////////
		////// CONSTRUCT THE nw SURFACE /////////
		/////////////////////////////////////////
		
		else if (Fluid_Interface(Phase,SignDist,fluid_isovalue,i,j,k) == 1){
			MC(Phase, fluid_isovalue, SignDist, i,j,k, nw_pts, n_nw_pts, nw_tris, n_nw_tris);
		}
		//******END OF BLOB PMMC*********************************************
		
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
		n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
		n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
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
