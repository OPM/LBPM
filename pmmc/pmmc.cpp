#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "Array.h"

#include "PointList.h"
//#include "vecLib/clapack.h"

using namespace std;

//--------------------------------------------------------------------------------------------------------
int ComputeBlob(IntArray &blobs, int &nblobs, int &ncubes, IntArray &indicator,
					   DoubleArray &F, DoubleArray &S, double vf, double vs, int startx, int starty,
					   int startz, IntArray &temp)
{
	// Compute the blob (F>vf|S>vs) starting from (i,j,k) - oil blob
	// F>vf => oil phase S>vs => in porespace
	// update the list of blobs, indicator mesh
	int m = F.m;  // maxima for the meshes
	int n = F.n;
	int o = F.o;

	int cubes_in_blob=0;
	int nrecent = 1;						// number of nodes added at most recent sweep
	temp(0,0) = startx;				// Set the initial point as a "seed" for the sweeps
	temp(1,0) = starty;
	temp(2,0) = startz;
	int ntotal = 1;					// total number of nodes in blob
	indicator(startx,starty,startz) = nblobs;

	int p,s,x,y,z,start,finish,nodx,nody,nodz;
	int imin=startx,imax=startx,jmin=starty,jmax=starty;	// initialize maxima / minima
	int kmin=startz,kmax=startz;
	int d[26][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
	{1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},
	{1,0,-1},{-1,0,-1},{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
	{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},
	{-1,-1,1},{-1,-1,-1}};   // directions to neighbors
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	bool status = 1;						// status == true => continue to look for points
	while (status == 1){
		start = ntotal - nrecent;
		finish = ntotal;
		nrecent = 0;						// set recent points back to zero for next sweep through
		for (s=start;s<finish;s++){
			// Loop over recent points; look for new points
			x = temp(0,s);
			y = temp(1,s);
			z = temp(2,s);
			// Looop over the directions
			for (p=0;p<26;p++){
				nodx=x+d[p][0];
				if (nodx < 0 ){ nodx = m-1; }		// Periodic BC for x
				if (nodx > m-1 ){ nodx = 0; }
				nody=y+d[p][1];
				if (nody < 0 ){ nody = n-1; }	// Periodic BC for y
				if (nody > n-1 ){ nody = 0; }
				nodz=z+d[p][2];
				if (nodz < 0 ){ nodz = 0; }		// No periodic BC for z
				if (nodz > o-1 ){ nodz = o-1; }
				if ( F(nodx,nody,nodz) > vf && S(nodx,nody,nodz) > vs
					 && indicator(nodx,nody,nodz) == -1 ){
					// Node is a part of the blob - add it to the list
					temp(0,ntotal) = nodx;
					temp(1,ntotal) = nody;
					temp(2,ntotal) = nodz;
					ntotal++;
					nrecent++;
					// Update the indicator map
					indicator(nodx,nody,nodz) = nblobs;
					// Update the min / max for the cube loop
					if ( nodx < imin ){ imin = nodx; }
					if ( nodx > imax ){ imax = nodx; }
					if ( nody < jmin ){ jmin = nody; }
					if ( nody > jmax ){ jmax = nody; }
					if ( nodz < kmin ){ kmin = nodz; }
					if ( nodz > kmax ){ kmax = nodz; }
				}
				else if (F(nodx,nody,nodz) > vf && S(nodx,nody,nodz) > vs
						 && indicator(nodx,nody,nodz) > -1 &&  indicator(nodx,nody,nodz) != nblobs){
					// Some kind of error in algorithm
					printf("Error in blob search algorithm!");
				}
			}

		}
		if ( nrecent == 0){
			status = 0;
		}
	}
	// Use points in temporary storage array to add cubes to the list of blobs
	if ( imin > 0) { imin = imin-1; }
//	if ( imax < m-1) { imax = imax+1; }
	if ( jmin > 0) { jmin = jmin-1; }
//	if ( jmax < n-1) { jmax = jmax+1; }
	if ( kmin > 0) { kmin = kmin-1; }
//	if ( kmax < o-1) { kmax = kmax+1; }
	int i,j,k;
	bool add;
	for (k=kmin;k<kmax;k++){
		for (j=jmin;j<jmax;j++){
			for (i=imin;i<imax;i++){
				// If cube(i,j,k) has any nodes in blob, add it to the list
				// Loop over cube edges
				add = 0;
				for (p=0;p<8;p++){
					nodx = i+cube[p][0];
					nody = j+cube[p][1];
					nodz = k+cube[p][2];
					if ( indicator(nodx,nody,nodz) == nblobs ){
						// Cube corner is in this blob
						add = 1;
					}
				}
				if (add == 1){
					// Add cube to the list
					blobs(0,ncubes) = i;
					blobs(1,ncubes) = j;
					blobs(2,ncubes) = k;
					ncubes++;
					cubes_in_blob++;
					// Loop again to check for overlap
					for (p=0;p<8;p++){
						nodx = i+cube[p][0];
						nody = j+cube[p][1];
						nodz = k+cube[p][2];
						if (indicator(nodx,nody,nodz) > -1 && indicator(nodx,nody,nodz) != nblobs){
							printf("Overlapping cube!");
							cout << i << ", " << j << ", " << k << endl;
						}
					}
				}
			}
		}
	}

	return cubes_in_blob;
}
//--------------------------------------------------------------------------------------------------
/*
DoubleArray SOLVE( DoubleArray &A, DoubleArray &b)
{
	// solves the system A*x = b exactly

	DoubleArray AA = A.Copy();
	long int n = A.n;
	long int nrhs =1;
	IntArray piv(2*n);
	DoubleArray solution = b.Copy();

	long int info = 0;
	int ret = dgesv_(&n,&nrhs,AA.Pointer(),&n,(long int *)piv.Pointer(),solution.Pointer(),&n,&info);

	//	DTErrorMessage("Computation", string("Error=")+ DTInt2String(ret) + "info = " + DTInt2String(info));

    return solution;
}
//-----------------------------------------------------------------------------
Point NEWTON(Point x0,  DoubleArray &F, double &v,DoubleArray &S, int i, int j, int k, int &newton_steps)
{
	Point pt;
	// Performs a Newton iteration to compute the point x from initial guess x0

	double tol = 0.001;
	double dist = 1;

	// Map x0 into unit coordinates
	Point X0 = x0;
	x0.x = x0.x - i;
	x0.y = x0.y - j;
	x0.z = x0.z - k;

	// Compute coefficients for tri-linear approximations to the functions F and S
	// coefficients for F
	double Af, Bf, Cf, Df, Ef, Ff, Gf, Hf;
	Af = F(i,j,k);
	Bf = F(i+1,j,k)-Af;
	Cf = F(i,j+1,k)-Af;
	Df = F(i,j,k+1)-Af;
	Ef = F(i+1,j+1,k)-Af-Bf-Cf;
	Ff = F(i+1,j,k+1)-Af-Bf-Df;
	Gf = F(i,j+1,k+1)-Af-Cf-Df;
	Hf = F(i+1,j+1,k+1)-Af-Bf-Cf-Df-Ef-Ff-Gf;
	// coefficients for S
	double As, Bs, Cs, Ds, Es, Fs, Gs, Hs;
	As = S(i,j,k);
	Bs = S(i+1,j,k)-As;
	Cs = S(i,j+1,k)-As;
	Ds = S(i,j,k+1)-As;
	Es = S(i+1,j+1,k)-As-Bs-Cs;
	Fs = S(i+1,j,k+1)-As-Bs-Ds;
	Gs = S(i,j+1,k+1)-As-Cs-Ds;
	Hs = S(i+1,j+1,k+1)-As-Bs-Cs-Ds-Es-Fs-Gs;

	// Compute coefficients for plane in direction of grad F, grad S at point x0
	// Compute grad u = grad F , v = grad S
	int count; int number;
	count = 0;
	number = int(floor(newton_steps));
	//number = 2;
	pt = x0;

	bool pt_on_cube_leg = 0;

	// Compute the u = grad F, v = grad S and the normal vector N = u x v
	double u1, u2, u3, v1, v2, v3;
	u1 = Bf+Ef*x0.y+Ff*x0.z+Hf*x0.y*x0.z;
	u2 = Cf+Ef*x0.x+Gf*x0.z+Hf*x0.x*x0.z;
	u3 = Df+Ff*x0.x+Gf*x0.y+Hf*x0.x*x0.y;
	v1 = Bs+Es*x0.y+Fs*x0.z+Hs*x0.y*x0.z;
	v2 = Cs+Es*x0.x+Gs*x0.z+Hs*x0.x*x0.z;
	v3 = Ds+Fs*x0.x+Gs*x0.y+Hs*x0.x*x0.y;
	// Compute the normal vector N
	Point N;
	if (x0.x == 0 || x0.x == 1){
		if (x0.y == 0 || x0.y == 1){
			// on a cube leg
			pt_on_cube_leg = 1;
		}
		else if (x0.z == 0 || x0.z == 1){
			// on a cube leg
			pt_on_cube_leg = 1;
		}
		else {
			// use this cube face as the normal
			N.x = 1;
			N.y = 0;
			N.z = 0;
		}
	}
	else if (x0.y == 0 || x0.y == 1){
		if (x0.z == 0 || x0.z == 1){
			// on a cube leg
			pt_on_cube_leg = 1;
		}
		else {
			// use this cube face as the normal
			N.x = 0;
			N.y = 1;
			N.z = 0;
		}
	}
	else if (x0.z == 0 || x0.z == 1){
		// use this cube face as the normal
		N.x = 0;
		N.y = 0;
		N.z = 1;
	}
	else {
		// Compute N = u x v
		N.x = u2*v3 - u3*v2;
		N.y = u3*v1 - u1*v3;
		N.z = u1*v2 - u2*v1;
	}

	if ( pt_on_cube_leg == 1){
		// return x0
		pt = x0;
	}
	else {
		while (count < number ){

			// Construct the Jacobean Matrix J
			DoubleArray J(3,3);
			J(0,0) = u1;
			J(0,1) = u2;
			J(0,2) = u3;
			J(1,0) = v1;
			J(1,1) = v2;
			J(1,2) = v3;
			J(2,0) = N.x;
			J(2,1) = N.y;
			J(2,2) = N.z;
			// Evaluate the vector M(x0)
			DoubleArray M(3);
			M(0) = -(Af+Bf*x0.x+Cf*x0.y+Df*x0.z+Ef*x0.x*x0.y+Ff*x0.x*x0.z+Gf*x0.y*x0.z+Hf*x0.x*x0.y*x0.z - v);
			M(1) = -(As+Bs*x0.x+Cs*x0.y+Ds*x0.z+Es*x0.x*x0.y+Fs*x0.x*x0.z+Gs*x0.y*x0.z+Hs*x0.x*x0.y*x0.z);
			M(2) = 0;
			// Compute the update to x0 using Newton's method
			Point dX;
			DoubleArray d = SOLVE(J,M);
			dX.x = d(0);
			dX.y = d(1);
			dX.z = d(2);

			if ( dX.x > 0.1 || dX.y > 0.1 || dX.z > 0.1 ||
				 dX.x < -0.1 || dX.y < -0.1 || dX.z < -0.1){
				count = number;
			}
			else {
				pt = x0+dX;
				x0 = pt;
			}
			M(0) = -(Af+Bf*x0.x+Cf*x0.y+Df*x0.z+Ef*x0.x*x0.y+Ff*x0.x*x0.z+Gf*x0.y*x0.z+Hf*x0.x*x0.y*x0.z - v);
			M(1) = -(As+Bs*x0.x+Cs*x0.y+Ds*x0.z+Es*x0.x*x0.y+Fs*x0.x*x0.z+Gs*x0.y*x0.z+Hs*x0.x*x0.y*x0.z);
			M(2) = 0;

			dist = sqrt(M(0)*M(0)+M(1)*M(1)+M(2)*M(2));

			// Compute u, v at new x0

			u1 = Bf+Ef*x0.y+Ff*x0.z+Hf*x0.y*x0.z;
			u2 = Cf+Ef*x0.x+Gf*x0.z+Hf*x0.x*x0.z;
			u3 = Df+Ff*x0.x+Gf*x0.y+Hf*x0.x*x0.y;
			v1 = Bs+Es*x0.y+Fs*x0.z+Hs*x0.y*x0.z;
			v2 = Cs+Es*x0.x+Gs*x0.z+Hs*x0.x*x0.z;
			v3 = Ds+Fs*x0.x+Gs*x0.y+Hs*x0.x*x0.y;

			count++;
		}
	}
	// map pt back to original coordinates
	pt.x = pt.x+i;
	pt.y = pt.y+j;
	pt.z = pt.z+k;

	return pt;
}
 */
//--------------------------------------------------------------------------------------------------

bool vertexcheck(Point &P, int n, int pos, DTMutableList<Point> &cellvertices){

    // returns true if P is a new vertex (one previously unencountered
    bool V = 1;
    int i = pos-n;
    for (i = pos-n; i < pos; i++){
        if ( P == cellvertices(i)){
            V = 0;
        }
    }

    return V;
}

//--------------------------------------------------------------------------------------------------

bool ShareSide( Point &A,  Point &B)
{
    // returns true if points A and B share an x,y, or z coordinate
    bool l = 0;
    if ( A.x != B.x && A.y != B.y && A.z != B.z){
        l=0;
    }
    else{
        if (floor(A.x)==A.x && floor(B.x)==B.x && A.x==B.x) {
            l = 1;
        }
        if (floor(A.y)==A.y && floor(B.y)==B.y && A.y==B.y) {
            l = 1;
        }
        if (floor(A.z)==A.z && floor(B.z)==B.z && A.z==B.z) {
            l = 1;
        }
    }

    return l;
}
//--------------------------------------------------------------------------------------------------

bool Interface( DoubleArray &A, const double v, int i, int j, int k){
    // returns true if grid cell i, j, k contains a section of the interface
    bool Y = 0;

    if ((A(i,j,k)-v)*(A(i+1,j,k)-v) < 0 ){
		Y=1;

    }
    // 2
    if ((A(i+1,j,k)-v)*(A(i+1,j+1,k)-v) < 0){
		Y=1;
    }
    //3
    if ((A(i+1,j+1,k)-v)*(A(i,j+1,k)-v) < 0){
		Y=1;
    }
    //4
    if ((A(i,j+1,k)-v)*(A(i,j,k)-v) < 0 ){
		Y=1;
    }
    //5
    if ((A(i,j,k)-v)*(A(i,j,k+1)-v) < 0 ){
		Y=1;
    }
    //6
    if ((A(i+1,j,k)-v)*(A(i+1,j,k+1)-v) < 0 ){
		Y=1;
    }
    //7
    if ((A(i+1,j+1,k)-v)*(A(i+1,j+1,k+1)-v) < 0 ){
		Y=1;
    }
    //8
    if ((A(i,j+1,k)-v)*(A(i,j+1,k+1)-v) < 0 ){
		Y=1;
    }
    //9
    if ((A(i,j,k+1)-v)*(A(i+1,j,k+1)-v) < 0 ){
		Y=1;
    }
    //10
    if ((A(i+1,j,k+1)-v)*(A(i+1,j+1,k+1)-v) < 0 ){
		Y=1;
    }
    //11
    if ((A(i+1,j+1,k+1)-v)*(A(i,j+1,k+1)-v) < 0 ){
		Y=1;
    }
    //12
    if ((A(i,j+1,k+1)-v)*(A(i,j,k+1)-v) < 0 ){
		Y=1;
    }
    return Y;
}
//--------------------------------------------------------------------------------------------------

bool Fluid_Interface( DoubleArray &A,  DoubleArray &S, const double v, int i, int j, int k){
    // returns true if grid cell i, j, k contains a section of the interface
    bool Y = 0;

    if ((A(i,j,k)-v)*(A(i+1,j,k)-v) < 0 && S(i,j,k) > 0 && S(i+1,j,k) > 0){
		Y=1;
    }
    // 2
    if ((A(i+1,j,k)-v)*(A(i+1,j+1,k)-v) < 0 && S(i+1,j,k) > 0 && S(i+1,j+1,k) > 0){
		Y=1;
    }
    //3
    if ((A(i+1,j+1,k)-v)*(A(i,j+1,k)-v) < 0 && S(i+1,j+1,k) > 0 && S(i,j+1,k) > 0){
		Y=1;
    }
    //4
    if ((A(i,j+1,k)-v)*(A(i,j,k)-v) < 0 && S(i,j,k) > 0 && S(i,j+1,k) > 0){
		Y=1;
    }
    //5
    if ((A(i,j,k)-v)*(A(i,j,k+1)-v) < 0 && S(i,j,k) > 0 && S(i,j,k+1) > 0){
		Y=1;
    }
    //6
    if ((A(i+1,j,k)-v)*(A(i+1,j,k+1)-v) < 0 && S(i+1,j,k) > 0 && S(i+1,j,k+1) > 0){
		Y=1;
    }
    //7
    if ((A(i+1,j+1,k)-v)*(A(i+1,j+1,k+1)-v) < 0 && S(i+1,j+1,k) > 0 && S(i+1,j+1,k+1) > 0){
		Y=1;
    }
    //8
    if ((A(i,j+1,k)-v)*(A(i,j+1,k+1)-v) < 0 && S(i,j+1,k) > 0 && S(i,j+1,k+1) > 0){
		Y=1;
    }
    //9
    if ((A(i,j,k+1)-v)*(A(i+1,j,k+1)-v) < 0 && S(i,j,k+1) > 0 && S(i+1,j,k+1) > 0){
		Y=1;
    }
    //10
    if ((A(i+1,j,k+1)-v)*(A(i+1,j+1,k+1)-v) < 0 && S(i+1,j,k+1) > 0 && S(i+1,j+1,k+1) > 0){
		Y=1;
    }
    //11
    if ((A(i+1,j+1,k+1)-v)*(A(i,j+1,k+1)-v) < 0 && S(i+1,j+1,k+1) > 0 && S(i,j+1,k+1) > 0){
		Y=1;
    }
    //12
    if ((A(i,j+1,k+1)-v)*(A(i,j,k+1)-v) < 0 && S(i,j,k+1) > 0 && S(i,j+1,k+1) > 0){
		Y=1;
    }
    return Y;
}
//--------------------------------------------------------------------------------------------------
bool Solid( DoubleArray &A, int i, int j, int k){

    bool X = 0;

    // return 0 if there is no solid phase in grid cell i,j,k

    if (A(i,j,k) == 0){
        X = 1;
    }
    if (A(i+1,j,k) == 0){
        X = 1;
    }
    if (A(i,j+1,k) == 0){
        X = 1;
    }
    if (A(i,j,k+1) == 0){
        X = 1;
    }
    if (A(i+1,j+1,k) == 0){
        X = 1;
    }
    if (A(i+1,j,k+1) == 0){
        X = 1;
    }
    if (A(i,j+1,k+1) == 0){
        X = 1;
    }
    if (A(i+1,j+1,k+1) == 0){
        X = 1;
    }
    return X;
}
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
void SOL_SURF(DoubleArray &A, const double &v, DoubleArray &B, const double &isovalue,
					 int i,int j,int k, int m, int n, int o, DTMutableList<Point>
					 &cellvertices, int &lengthvertices, IntArray &Tlist, int &nTris,
					 DoubleArray &values){

	// THIS SUBROUTINE COMPUTES THE VERTICES FOR THE SOLID PHASE IN
	// A PARTICULAR GRID CELL, THEN ARRANGES THEM INTO TRIANGLES
	// ALSO ORGANIZES THE LIST OF VALUES TO CORRESPOND WITH VERTEX LIST

	int N = 0;
	Point P;
	Point PlaceHolder;
	int pos = lengthvertices;
	float temp;


	// int m; int n; int o;
	int p; int q;

	// for each leg of the triangle, see if a vertex exists
	// if so, find the vertex, and perform an extrapolation

	// Go over each corner -- check to see if the corners are themselves vertices
	//1
	if (A(i,j,k) == v){
		P.x = i;
		P.y = j;
		P.z = k;
		values(pos) = B(i,j,k);
		cellvertices(pos++) = P;
		N++;
	}
	//2
	if (A(i+1,j,k) == v){
		P.x = i+1;
		P.y = j;
		P.z = k;
		values(pos) = B(i+1,j,k);
		cellvertices(pos++) = P;
		N++;
	}
	//3
	if (A(i+1,j+1,k) == v){
		P.x = i+1;
		P.y = j+1;
		P.z = k;
		values(pos) = B(i+1,j+1,k);
		cellvertices(pos++) = P;
		N++;
	}
	//4
	if (A(i,j+1,k) == v){
		P.x = i;
		P.y = j+1;
		P.z = k;
		values(pos) = B(i,j+1,k);
		cellvertices(pos++) = P;
		N++;
	}
	//5
	if (A(i,j,k+1) == v){
		P.x = i;
		P.y = j;
		P.z = k+1;
		values(pos) = B(i,j,k+1);
		cellvertices(pos++) = P;
		N++;
	}
	//6
	if (A(i+1,j,k+1) == v){
		P.x = i+1;
		P.y = j;
		P.z = k+1;
		values(pos) = B(i+1,j,k+1);
		cellvertices(pos++) = P;
		N++;
	}
	//7
	if (A(i+1,j+1,k+1) == v){
		P.x = i+1;
		P.y = j+1;
		P.z = k+1;
		values(pos) = B(i+1,j+1,k+1);
		cellvertices(pos++) = P;
		N++;
	}
	//8
	if (A(i,j+1,k+1) == v){
		P.x = i;
		P.y = j+1;
		P.z = k+1;
		values(pos) = B(i,j+1,k+1);
		cellvertices(pos++) = P;
		N++;
	}

    // Go through each side, compute P for sides of box spiraling up


    if ((A(i,j,k)-v)*(A(i+1,j,k)-v) < 0)
	{
        P.x = i + (A(i,j,k)-v)/(A(i,j,k)-A(i+1,j,k));
        P.y = j;
        P.z = k;
		// compute extrapolated fluid value at P
		//		if ( A(i,j,k) > v){
		//			 values(pos) = EXTRAP(B, isovalue, i,j,k,4, P);
		//		}
		//		else{
		//			values(pos) = EXTRAP(B,isovalue, i+1,j,k, 1, P);
		//		}
		// Interpolate value at P using padded mesh B
		values(pos) = B(i,j,k)*(1-P.x+i)+B(i+1,j,k)*(P.x-i);

        cellvertices(pos++) =  P;
        N++;
	}
    // 2
    if ((A(i+1,j,k)-v)*(A(i+1,j+1,k)-v) < 0)
	{
        P.x = i+1;
        P.y = j + (A(i+1,j,k)-v)/(A(i+1,j,k)-A(i+1,j+1,k));
        P.z = k;
        if (vertexcheck(P, N, pos, cellvertices) == 1){ // P is a new vertex (not counted twice)
														// compute extrapolated fluid value at P
														//			if ( A(i+1,j,k) > v){
														//				values(pos) = EXTRAP(B,isovalue, i+1,j,k, 5, P);
														//			}
														//			else{
														//				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k, 2, P);
														//			}
														// Interpolate value at P using padded mesh B
			values(pos) = B(i+1,j,k)*(1-P.y+j)+B(i+1,j+1,k)*(P.y-j);
            cellvertices(pos++) =  P;
            N++;
        }
	}
    //3
    if ((A(i+1,j+1,k)-v)*(A(i,j+1,k)-v) < 0)
	{
        P.x = i + (A(i,j+1,k)-v) / (A(i,j+1,k)-A(i+1,j+1,k));
        P.y = j+1;
        P.z = k;
        if (vertexcheck(P, N, pos, cellvertices) == 1){ // P is a new vertex (not counted twice)
														// compute extrapolated fluid value at P
														//			if ( A(i+1,j+1,k) > v){
														//				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k, 1, P);
														//			}
														//			else{
														//				values(pos) = EXTRAP(B,isovalue, i,j+1,k, 4, P);
														//			}
														// Interpolate value at P using padded mesh B
			values(pos) = B(i,j+1,k)*(1-P.x+i)+B(i+1,j+1,k)*(P.x-i);
            cellvertices(pos++) =  P;
            N++;
        }
	}
    //4
    if ((A(i,j+1,k)-v)*(A(i,j,k)-v) < 0)
	{
        P.x = i;
        P.y = j + (A(i,j,k)-v) / (A(i,j,k)-A(i,j+1,k));
        P.z = k;
        if (vertexcheck(P, N, pos, cellvertices) == 1){ // P is a new vertex (not counted twice)
														// compute extrapolated fluid value at P
														//			if ( A(i,j,k) > v){
														//				values(pos) = EXTRAP(B,isovalue, i,j,k, 5, P);
														//			}
														//			else{
														//				values(pos) = EXTRAP(B,isovalue, i,j+1,k, 2, P);
														//			}
														// Interpolate value at P using padded mesh B
			values(pos) = B(i,j,k)*(1-P.y+j)+B(i,j+1,k)*(P.y-j);
            cellvertices(pos++) =  P;
            N++;
        }
	}
    //5
    if ((A(i,j,k)-v)*(A(i,j,k+1)-v) < 0)
	{
        P.x = i;
        P.y = j;
        P.z = k + (A(i,j,k)-v) / (A(i,j,k)-A(i,j,k+1));
		if (vertexcheck(P, N, pos, cellvertices) == 1){ // P is a new vertex (not counted twice)
														// compute extrapolated fluid value at P
														//			if ( A(i,j,k) > v){
														//				values(pos) = EXTRAP(B,isovalue, i,j,k, 6, P);
														//			}
														//			else{
														//				values(pos) = EXTRAP(B,isovalue, i,j,k+1, 3, P);
														//			}

			// Interpolate value at P using padded mesh B
			values(pos) = B(i,j,k)*(1-P.z+k)+B(i,j,k+1)*(P.z-k);
            cellvertices(pos++) =  P;
            N++;
        }
	}
    //6
    if ((A(i+1,j,k)-v)*(A(i+1,j,k+1)-v) < 0)
	{
        P.x = i+1;
        P.y = j;
        P.z = k + (A(i+1,j,k)-v) / (A(i+1,j,k)-A(i+1,j,k+1));
        if (vertexcheck(P, N, pos, cellvertices) == 1){  // P is a new vertex (not counted twice)
														 // compute extrapolated fluid value at P
														 //			if ( A(i+1,j,k) > v){
														 //				values(pos) = EXTRAP(B,isovalue, i+1,j,k, 6, P);
														 //			}
														 //			else{
														 //				values(pos) = EXTRAP(B,isovalue, i+1,j,k+1, 3, P);
														 //			}

			// Interpolate value at P using padded mesh B
			values(pos) = B(i+1,j,k)*(1-P.z+k)+B(i+1,j,k+1)*(P.z-k);
			cellvertices(pos++) =  P;
            N++;
        }
	}
    //7
    if ((A(i+1,j+1,k)-v)*(A(i+1,j+1,k+1)-v) < 0)
	{
        P.x = i+1;
        P.y = j+1;
        P.z = k + (A(i+1,j+1,k)-v) / (A(i+1,j+1,k)-A(i+1,j+1,k+1));
        if (vertexcheck(P, N, pos, cellvertices) == 1){  // P is a new vertex (not counted twice)
														 // compute extrapolated fluid value at P
														 //			if ( A(i+1,j+1,k) > v){
														 //				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k, 6, P);
														 //			}
														 //			else{
														 //				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k+1, 3, P);
														 //			}

			// Interpolate value at P using padded mesh B
			values(pos) = B(i+1,j+1,k)*(1-P.z+k)+B(i+1,j+1,k+1)*(P.z-k);
			cellvertices(pos++) =  P;
            N++;
        }
	}
    //8
    if ((A(i,j+1,k)-v)*(A(i,j+1,k+1)-v) < 0)
	{
        P.x = i;
        P.y = j+1;
        P.z = k + (A(i,j+1,k)-v) / (A(i,j+1,k)-A(i,j+1,k+1));
        if (vertexcheck(P, N, pos, cellvertices) == 1){  // P is a new vertex (not counted twice)
														 // compute extrapolated fluid value at P
														 //			if ( A(i,j+1,k) > v){
														 //				values(pos) = EXTRAP(B,isovalue, i,j+1,k, 6, P);
														 //			}
														 //			else{
														 //				values(pos) = EXTRAP(B,isovalue, i,j+1,k+1, 3, P);
														 //			}

			// Interpolate value at P using padded mesh B
			values(pos) = B(i,j+1,k)*(1-P.z+k)+B(i,j+1,k+1)*(P.z-k);
			cellvertices(pos++) =  P;
            N++;
        }
	}
    //9
    if ((A(i,j,k+1)-v)*(A(i+1,j,k+1)-v) < 0)
	{
        P.x = i + (A(i,j,k+1)-v) / (A(i,j,k+1)-A(i+1,j,k+1));
        P.y = j;
        P.z = k+1;
        if (vertexcheck(P, N, pos, cellvertices) == 1){  // P is a new vertex (not counted twice)
														 // compute extrapolated fluid value at P
														 //			if ( A(i,j,k+1) > v){
														 //				values(pos) = EXTRAP(B,isovalue, i,j,k+1, 4, P);
														 //			}
														 //			else{
														 //				values(pos) = EXTRAP(B,isovalue, i+1,j,k+1, 1, P);
														 //			}

			// Interpolate value at P using padded mesh B
			values(pos) = B(i,j,k+1)*(1-P.x+i)+B(i+1,j,k+1)*(P.x-i);
			cellvertices(pos++) =  P;
            N++;
        }
	}
    //10
    if ((A(i+1,j,k+1)-v)*(A(i+1,j+1,k+1)-v) < 0)
	{
        P.x = i+1;
        P.y = j + (A(i+1,j,k+1)-v) / (A(i+1,j,k+1)-A(i+1,j+1,k+1));
        P.z = k+1;
        if (vertexcheck(P, N, pos, cellvertices) == 1){  // P is a new vertex (not counted twice)
														 // compute extrapolated fluid value at P
														 //			if ( A(i+1,j,k+1) > v){
														 //				values(pos) = EXTRAP(B,isovalue, i+1,j,k+1, 5, P);
														 //			}
														 //			else{
														 //				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k+1, 2, P);
														 //			}

			// Interpolate value at P using padded mesh B
			values(pos) = B(i+1,j,k+1)*(1-P.y+j)+B(i+1,j+1,k+1)*(P.y-j);
			cellvertices(pos++) =  P;
			N++;
        }
	}
    //11
    if ((A(i+1,j+1,k+1)-v)*(A(i,j+1,k+1)-v) < 0)
	{
        P.x = i+(A(i,j+1,k+1)-v) / (A(i,j+1,k+1)-A(i+1,j+1,k+1));
        P.y = j+1;
        P.z = k+1;
        if (vertexcheck(P, N, pos, cellvertices) == 1){  // P is a new vertex (not counted twice)
														 // compute extrapolated fluid value at P
														 //			if ( A(i,j+1,k+1) > v){
														 //				values(pos) = EXTRAP(B,isovalue, i,j+1,k+1, 4, P);
														 //			}
														 //			else{
														 //				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k+1, 1, P);
														 //			}

			// Interpolate value at P using padded mesh B
			values(pos) = B(i,j+1,k+1)*(1-P.x+i)+B(i+1,j+1,k+1)*(P.x-i);
			cellvertices(pos++) =  P;
            N++;
        }
	}
    //12
    if ((A(i,j+1,k+1)-v)*(A(i,j,k+1)-v) < 0)
	{
        P.x = i;
        P.y = j + (A(i,j,k+1)-v) / (A(i,j,k+1)-A(i,j+1,k+1));
        P.z = k+1;
        if (vertexcheck(P, N, pos, cellvertices) == 1){  // P is a new vertex (not counted twice)
														 // compute extrapolated fluid value at P
														 //			if ( A(i,j,k+1) > v){
														 //				values(pos) = EXTRAP(B,isovalue, i,j,k+1, 5, P);
														 //			}
														 //			else{
														 //				values(pos) = EXTRAP(B,isovalue, i,j+1,k+1, 2, P);
														 //			}

			// Interpolate value at P using padded mesh B
			values(pos) = B(i,j,k+1)*(1-P.y+j)+B(i,j+1,k+1)*(P.y-j);
			cellvertices(pos++) =  P;
            N++;
        }
	}

	lengthvertices = pos;

	// *    *    *   ARRANGE VERTICES SO THAT NEIGHBORS SHARE A FACE    *    *    *
	// *    *    *    PERFORM SAME OPERATIONS TO THE LIST OF VALUES     *    *    *

    for (q = pos-N; q < pos-2; q++) {
        for (o = q+2; o < pos-1; o++) {
            if (ShareSide(cellvertices(q), cellvertices(o)) == 1) {
                PlaceHolder = cellvertices(q+1);
                cellvertices(q+1) = cellvertices(o);
                cellvertices(o) = PlaceHolder;

				temp = values(q+1);
				values(q+1) = values(o);
				values(o) = temp;
            }
        }

		// make sure other neighbor of vertex 1 is in last spot
        if (q == pos-N){
            for (p = q+2; p < pos-1; p++){
                if (ShareSide(cellvertices(q), cellvertices(p)) == 1){
                    PlaceHolder = cellvertices(pos-1);
                    cellvertices(pos-1) = cellvertices(p);
                    cellvertices(p) = PlaceHolder;

					temp = values(pos-1);
					values(pos-1) = values(p);
					values(p) = temp;
                }
            }
        }
        if ( ShareSide(cellvertices(pos-2), cellvertices(pos-3)) != 1 ){
            if (ShareSide( cellvertices(pos-3), cellvertices(pos-1)) == 1 && ShareSide( cellvertices(pos-N),
																						cellvertices(pos-2)) == 1 ){
                PlaceHolder = cellvertices(pos-2);
                cellvertices(pos-2) = cellvertices(pos-1);
                cellvertices(pos-1) = PlaceHolder;

				temp = values(pos-2);
				values(pos-2) = values(pos-1);
				values(pos-1) = temp;
            }
        }
        if ( ShareSide(cellvertices(pos-1), cellvertices(pos-2)) != 1 ){
            if (ShareSide( cellvertices(pos-3), cellvertices(pos-1)) == 1 && ShareSide(cellvertices(pos-4),
																					   cellvertices(pos-2)) == 1 ){
                PlaceHolder = cellvertices(pos-3);
                cellvertices(pos-3) = cellvertices(pos-2);
                cellvertices(pos-2) = PlaceHolder;

				temp = values(pos-3);
				values(pos-3) = values(pos-2);
				values(pos-2) = temp;
            }
            if (ShareSide( cellvertices(pos-N+1), cellvertices(pos-3)) == 1 && ShareSide(cellvertices(pos-1),
																						 cellvertices(pos-N+1)) == 1 ){
                PlaceHolder = cellvertices(pos-2);
                cellvertices(pos-2) = cellvertices(pos-N+1);
                cellvertices(pos-N+1) = PlaceHolder;

				temp = values(pos-2);
				values(pos-2) = values(pos-N+1);
				values(pos-N+1) = temp;
            }
        }
        if ( ShareSide(cellvertices(pos-N), cellvertices(pos-N+1)) != 1 ){
            if (ShareSide( cellvertices(pos-N), cellvertices(pos-2)) == 1 && ShareSide(cellvertices(pos-1),
																					   cellvertices(pos-N+1)) == 1){
                PlaceHolder = cellvertices(pos-1);
                cellvertices(pos-1) = cellvertices(pos-N);
                cellvertices(pos-N) = PlaceHolder;

				temp = values(pos-1);
				values(pos-1) = values(pos-N);
				values(pos-N) = temp;
            }
        }
    }


	// *    *    *   ESTABLISH TRIANGLE CONNECTIONS	   *    *    *

	for (p=pos-N+2; p<pos; p++){
		Tlist(0,nTris) = pos-N;
        Tlist(1,nTris) = p-1;
        Tlist(2,nTris) = p;
        nTris++;
	}
}
//-------------------------------------------------------------------------------
void TRIM(DTMutableList<Point> &local_sol_pts, int &n_local_sol_pts, double isovalue,
				 IntArray &local_sol_tris, int &n_local_sol_tris,
				 DTMutableList<Point> &ns_pts, int &n_ns_pts, IntArray &ns_tris,
				 int &n_ns_tris, DTMutableList<Point> &ws_pts, int &n_ws_pts,
				 IntArray &ws_tris, int &n_ws_tris, DoubleArray &values,
				 DTMutableList<Point> &local_nws_pts, int &n_local_nws_pts,
				 DoubleArray &fluid_pad, DoubleArray &S, int &i, int &j, int &k,
				 int &newton_steps){
	// Trim the local solid surface

	int map_ws;
	int map_ns;
	int pts_on_tri;
	int p; int q;

	int a; int b; int c;
	Point A; Point B; Point C;
	Point D; Point E; Point F;
	Point P;

	bool all_to_ns = 1;
	// Check to see if all triangles in the cell are in ns_surface
	for (q=0; q < n_local_sol_pts; q++){
		if ( values(q) <= isovalue && all_to_ns == 1){
			all_to_ns = 0;
		}
	}
	bool all_to_ws = 1;
	// Check to see if all triangles in the cell are in ws surface
	for (q=0; q < n_local_sol_pts; q++){
		if ( values(q) >= isovalue && all_to_ws == 1){
			all_to_ws = 0;
		}
	}

	if (all_to_ws == 1){
		map_ws = n_ws_pts;
		for ( p=0; p<n_local_sol_pts; p++){
			ws_pts(n_ws_pts++) = local_sol_pts(p);
		}
		for ( p=0; p<n_local_sol_tris; p++){
			ws_tris(0,n_ws_tris) = local_sol_tris(0,p) + map_ws;
			ws_tris(1,n_ws_tris) = local_sol_tris(1,p) + map_ws;
			ws_tris(2,n_ws_tris) = local_sol_tris(2,p) + map_ws;
			n_ws_tris++;
		}
	}
	else if (all_to_ns == 1){
		map_ns = n_ns_pts;
		for ( p=0; p<n_local_sol_pts; p++){
			ns_pts(n_ns_pts++) = local_sol_pts(p);
		}
		for ( p=0; p<n_local_sol_tris; p++){
			ns_tris(0,n_ns_tris) = local_sol_tris(0,p) + map_ns;
			ns_tris(1,n_ns_tris) = local_sol_tris(1,p) + map_ns;
			ns_tris(2,n_ns_tris) = local_sol_tris(2,p) + map_ns;
			n_ns_tris++;
		}
	}
	else {
		// this section of surface contains a common line
		map_ns = n_ns_tris;
		map_ws = n_ws_tris;
		// Go through all triangles
		for ( p=0; p<n_local_sol_tris; p++){
			a = local_sol_tris(0,p);
			b = local_sol_tris(1,p);
			c = local_sol_tris(2,p);
			A = local_sol_pts(a);
			B = local_sol_pts(b);
			C = local_sol_pts(c);
			if (values(a) >= isovalue && values(b) >= isovalue && values(c) >= isovalue ){
				// Triangle is in ns surface
				// Add points
				ns_pts(n_ns_pts++) = A;
				ns_pts(n_ns_pts++) = B;
				ns_pts(n_ns_pts++) = C;
				// Add triangles
				ns_tris(0,n_ns_tris) = n_ns_pts-3;
				ns_tris(1,n_ns_tris) = n_ns_pts-2;
				ns_tris(2,n_ns_tris) = n_ns_pts-1;
				n_ns_tris++;

			}
			else if (values(a) <= isovalue && values(b) <= isovalue && values(c) <= isovalue ){
				// Triangle is in ws surface
				// Add points
				ws_pts(n_ws_pts++) = A;
				ws_pts(n_ws_pts++) = B;
				ws_pts(n_ws_pts++) = C;
				// Add triangles
				ws_tris(0,n_ws_tris) = n_ws_pts-3;
				ws_tris(1,n_ws_tris) = n_ws_pts-2;
				ws_tris(2,n_ws_tris) = n_ws_pts-1;
				n_ws_tris++;
				// If two of values(a), values(b), values(c) --> this tri leg is on common line
				// Since this leg is part of two triangles, one in ws surface and one in ns surface
				// we only need to do this procedure one time
				if (values(a) == isovalue && values(b) == isovalue){
					if ( n_local_nws_pts == 0 ){
						local_nws_pts(n_local_nws_pts++) = A;
						local_nws_pts(n_local_nws_pts++) = B;
					}
					else if ( A != local_nws_pts(n_local_nws_pts-1) ){
						local_nws_pts(n_local_nws_pts++) = A;
						local_nws_pts(n_local_nws_pts++) = B;
					}
				}
				if (values(a) == isovalue && values(c) == isovalue){
					if ( n_local_nws_pts == 0 ){
						local_nws_pts(n_local_nws_pts++) = A;
						local_nws_pts(n_local_nws_pts++) = C;
					}
					else if ( A != local_nws_pts(n_local_nws_pts-1) ){
						local_nws_pts(n_local_nws_pts++) = A;
						local_nws_pts(n_local_nws_pts++) = C;
					}
				}
				if (values(b) == isovalue && values(c) == isovalue){
					if ( n_local_nws_pts == 0 ){
						local_nws_pts(n_local_nws_pts++) = B;
						local_nws_pts(n_local_nws_pts++) = C;
					}
					else if ( A != local_nws_pts(n_local_nws_pts-1) ){
						local_nws_pts(n_local_nws_pts++) = B;
						local_nws_pts(n_local_nws_pts++) = C;
					}
				}
			}
			else {
				// Triangle contains common line

				////////////////////////////////////////
				///////// FIND THE COMMON LINE /////////
				////////////////////////////////////////
				pts_on_tri = 0;
				if ( (values(a)-isovalue)*(values(b)-isovalue) < 0){
					// compute common line vertex
					P = A + (values(a) - isovalue)/(values(a)-values(b))*(B-A);
					if ( n_local_nws_pts == 0 ){
						local_nws_pts(n_local_nws_pts++) = P;
						pts_on_tri++;
					}
					else if ( P != local_nws_pts(n_local_nws_pts-1) ){
						local_nws_pts(n_local_nws_pts++) = P;
						pts_on_tri++;
					}
				}
				if ( (values(b)-isovalue)*(values(c)-isovalue) < 0){
					// compute common line vertex
					P = B + (values(b) - isovalue)/(values(b)-values(c))*(C-B);
					if ( n_local_nws_pts == 0 ){
						local_nws_pts(n_local_nws_pts++) = P;
						pts_on_tri++;
					}
					else if ( P != local_nws_pts(n_local_nws_pts-1) ){
						local_nws_pts(n_local_nws_pts++) = P;
						pts_on_tri++;
					}
				}
				if ( (values(a)-isovalue)*(values(c)-isovalue) < 0){
					// compute common line vertex
					P = A + (values(a) - isovalue)/(values(a)-values(c))*(C-A);
					local_nws_pts(n_local_nws_pts++) = P;
					pts_on_tri++;
				}
				// Also consider case in which common line points may be at vertices
				if ( values(a) == isovalue ){
					if ( n_local_nws_pts == 0 ){
						local_nws_pts(n_local_nws_pts++) = A;
						pts_on_tri++;
					}
					else if ( A != local_nws_pts(n_local_nws_pts-1) ){
						local_nws_pts(n_local_nws_pts++) = A;
						pts_on_tri++;
					}				}
				if ( values(b) == isovalue ){
					if ( n_local_nws_pts == 0 ){
						local_nws_pts(n_local_nws_pts++) = B;
						pts_on_tri++;
					}
					else if ( B != local_nws_pts(n_local_nws_pts-1) ){
						local_nws_pts(n_local_nws_pts++) = B;
						pts_on_tri++;
					}				}
				if ( values(c) == isovalue ){
					if ( n_local_nws_pts == 0 ){
						local_nws_pts(n_local_nws_pts++) = C;
						pts_on_tri++;
					}
					else if ( C != local_nws_pts(n_local_nws_pts-1) ){
						local_nws_pts(n_local_nws_pts++) = C;
						pts_on_tri++;
					}
				}
				// If there is a section of common line found
				// Use these common line points as an initial guess in Newton iteration
				D = local_nws_pts(n_local_nws_pts-2);
				E = local_nws_pts(n_local_nws_pts-1);
				F = 0.5*(local_nws_pts(n_local_nws_pts-2) + local_nws_pts(n_local_nws_pts-1));
//				D = NEWTON(local_nws_pts(n_local_nws_pts-2),fluid_pad,isovalue,S,i,j,k,newton_steps);
//				E = NEWTON(local_nws_pts(n_local_nws_pts-1),fluid_pad,isovalue,S,i,j,k,newton_steps);
//				F = NEWTON(F,fluid_pad,isovalue,S,i,j,k,newton_steps);
				// Write over the initial common line points
				// Note point F goes in between points D and E
				local_nws_pts(n_local_nws_pts-2) = D;
				local_nws_pts(n_local_nws_pts-1) = F;
				local_nws_pts(n_local_nws_pts++) = E;
				// Construct the new triangles
				if ( (values(a)-isovalue)*(values(b)-isovalue) < 0 &&
					 (values(b)-isovalue)*(values(c)-isovalue) < 0){
					if (values(b) > isovalue){
						// Points
						ns_pts(n_ns_pts++) = B;
						ns_pts(n_ns_pts++) = D;
						ns_pts(n_ns_pts++) = E;
						ns_pts(n_ns_pts++) = F;
						// Triangles
						ns_tris(0,n_ns_tris) = n_ns_pts-4; // B
						ns_tris(1,n_ns_tris) = n_ns_pts-1; // F
						ns_tris(2,n_ns_tris) = n_ns_pts-3; // D
						n_ns_tris++;
						ns_tris(0,n_ns_tris) = n_ns_pts-4; // B
						ns_tris(1,n_ns_tris) = n_ns_pts-2; // E
						ns_tris(2,n_ns_tris) = n_ns_pts-1; // F
						n_ns_tris++;
						// Points
						ws_pts(n_ws_pts++) = A; //-5
						ws_pts(n_ws_pts++) = C; //-4
						ws_pts(n_ws_pts++) = D; //-3
						ws_pts(n_ws_pts++) = E; //-2
						ws_pts(n_ws_pts++) = F; //-1
												// Triangles (A,D,F),(A,C,F),(C,E,F)
						ws_tris(0,n_ws_tris) = n_ws_pts-5; // A
						ws_tris(1,n_ws_tris) = n_ws_pts-1; // F
						ws_tris(2,n_ws_tris) = n_ws_pts-3; // D
						n_ws_tris++;
						ws_tris(0,n_ws_tris) = n_ws_pts-5; // A
						ws_tris(1,n_ws_tris) = n_ws_pts-4; // C
						ws_tris(2,n_ws_tris) = n_ws_pts-1; // F
						n_ws_tris++;
						ws_tris(0,n_ws_tris) = n_ws_pts-4; // C
						ws_tris(1,n_ws_tris) = n_ws_pts-2; // E
						ws_tris(2,n_ws_tris) = n_ws_pts-1; // F
						n_ws_tris++;
					}
					else {
						// Points
						ws_pts(n_ws_pts++) = B; //-4
						ws_pts(n_ws_pts++) = D; //-3
						ws_pts(n_ws_pts++) = E; //-2
						ws_pts(n_ws_pts++) = F; //-1
												// Triangles
						ws_tris(0,n_ws_tris) = n_ws_pts-4; // B
						ws_tris(1,n_ws_tris) = n_ws_pts-1; // F
						ws_tris(2,n_ws_tris) = n_ws_pts-3; // D
						n_ws_tris++;
						ws_tris(0,n_ws_tris) = n_ws_pts-4; // B
						ws_tris(1,n_ws_tris) = n_ws_pts-2; // E
						ws_tris(2,n_ws_tris) = n_ws_pts-1; // F
						n_ws_tris++;
						// Points
						ns_pts(n_ns_pts++) = A; //-5
						ns_pts(n_ns_pts++) = C; //-4
						ns_pts(n_ns_pts++) = D; //-3
						ns_pts(n_ns_pts++) = E; //-2
						ns_pts(n_ns_pts++) = F; //-1
												// Triangles (A,D,F),(A,C,F),(C,E,F)
						ns_tris(0,n_ns_tris) = n_ns_pts-5; // A
						ns_tris(1,n_ns_tris) = n_ns_pts-1; // F
						ns_tris(2,n_ns_tris) = n_ns_pts-3; // D
						n_ns_tris++;
						ns_tris(0,n_ns_tris) = n_ns_pts-5; // A
						ns_tris(1,n_ns_tris) = n_ns_pts-4; // C
						ns_tris(2,n_ns_tris) = n_ns_pts-1; // F
						n_ns_tris++;
						ns_tris(0,n_ns_tris) = n_ns_pts-4; // C
						ns_tris(1,n_ns_tris) = n_ns_pts-2; // E
						ns_tris(2,n_ns_tris) = n_ns_pts-1; // F
						n_ns_tris++;
					}
				}
				else if ( (values(a)-isovalue)*(values(b)-isovalue) < 0 &&
						  (values(a)-isovalue)*(values(c)-isovalue) < 0){
					if (values(a) > isovalue){
						// Points
						ns_pts(n_ns_pts++) = A; //-4
						ns_pts(n_ns_pts++) = D; //-3
						ns_pts(n_ns_pts++) = E; //-2
						ns_pts(n_ns_pts++) = F; //-1

						// Triangles
						ns_tris(0,n_ns_tris) = n_ns_pts-4; // A
						ns_tris(1,n_ns_tris) = n_ns_pts-1; // F
						ns_tris(2,n_ns_tris) = n_ns_pts-3; // D
						n_ns_tris++;
						ns_tris(0,n_ns_tris) = n_ns_pts-4; // A
						ns_tris(1,n_ns_tris) = n_ns_pts-2; // E
						ns_tris(2,n_ns_tris) = n_ns_pts-1; // F
						n_ns_tris++;
						// Points
						ws_pts(n_ws_pts++) = B; //-5
						ws_pts(n_ws_pts++) = C; //-4
						ws_pts(n_ws_pts++) = D; //-3
						ws_pts(n_ws_pts++) = E; //-2
						ws_pts(n_ws_pts++) = F; //-1
												// Triangles (B,D,F),(B,C,F),(C,E,F)
						ws_tris(0,n_ws_tris) = n_ws_pts-5; // B
						ws_tris(1,n_ws_tris) = n_ws_pts-1; // F
						ws_tris(2,n_ws_tris) = n_ws_pts-3; // D
						n_ws_tris++;
						ws_tris(0,n_ws_tris) = n_ws_pts-5; // B
						ws_tris(1,n_ws_tris) = n_ws_pts-4; // C
						ws_tris(2,n_ws_tris) = n_ws_pts-1; // F
						n_ws_tris++;
						ws_tris(0,n_ws_tris) = n_ws_pts-4; // C
						ws_tris(1,n_ws_tris) = n_ws_pts-2; // E
						ws_tris(2,n_ws_tris) = n_ws_pts-1; // F
						n_ws_tris++;
					}
					else {
						// Points
						ws_pts(n_ws_pts++) = A; //-4
						ws_pts(n_ws_pts++) = D; //-3
						ws_pts(n_ws_pts++) = E; //-2
						ws_pts(n_ws_pts++) = F; //-1
												// Triangles
						ws_tris(0,n_ws_tris) = n_ws_pts-4; // A
						ws_tris(1,n_ws_tris) = n_ws_pts-1; // F
						ws_tris(2,n_ws_tris) = n_ws_pts-3; // D
						n_ws_tris++;
						ws_tris(0,n_ws_tris) = n_ws_pts-4; // A
						ws_tris(1,n_ws_tris) = n_ws_pts-2; // E
						ws_tris(2,n_ws_tris) = n_ws_pts-1; // F
						n_ws_tris++;
						// Points
						ns_pts(n_ns_pts++) = B; //-5
						ns_pts(n_ns_pts++) = C; //-4
						ns_pts(n_ns_pts++) = D; //-3
						ns_pts(n_ns_pts++) = E; //-2
						ns_pts(n_ns_pts++) = F; //-1
												// Triangles (B,D,F),(B,C,F),(C,E,F)
						ns_tris(0,n_ns_tris) = n_ns_pts-5; // B
						ns_tris(1,n_ns_tris) = n_ns_pts-1; // F
						ns_tris(2,n_ns_tris) = n_ns_pts-3; // D
						n_ns_tris++;
						ns_tris(0,n_ns_tris) = n_ns_pts-5; // B
						ns_tris(1,n_ns_tris) = n_ns_pts-4; // C
						ns_tris(2,n_ns_tris) = n_ns_pts-1; // F
						n_ns_tris++;
						ns_tris(0,n_ns_tris) = n_ns_pts-4; // C
						ns_tris(1,n_ns_tris) = n_ns_pts-2; // E
						ns_tris(2,n_ns_tris) = n_ns_pts-1; // F
						n_ns_tris++;
					}
				}
				else {
					if (values(c) > isovalue){
						// Points
						ns_pts(n_ns_pts++) = C; //-4
						ns_pts(n_ns_pts++) = D; //-3
						ns_pts(n_ns_pts++) = E; //-2
						ns_pts(n_ns_pts++) = F; //-1
												// Triangles
						ns_tris(0,n_ns_tris) = n_ns_pts-4; // C
						ns_tris(1,n_ns_tris) = n_ns_pts-1; // F
						ns_tris(2,n_ns_tris) = n_ns_pts-3; // D
						n_ns_tris++;
						ns_tris(0,n_ns_tris) = n_ns_pts-4; // C
						ns_tris(1,n_ns_tris) = n_ns_pts-2; // E
						ns_tris(2,n_ns_tris) = n_ns_pts-1; // F
						n_ns_tris++;
						// Points
						ws_pts(n_ws_pts++) = A; //-5
						ws_pts(n_ws_pts++) = B; //-4
						ws_pts(n_ws_pts++) = D; //-3
						ws_pts(n_ws_pts++) = E; //-2
						ws_pts(n_ws_pts++) = F; //-1
												// Triangles (A,E,F),(A,B,F),(B,D,F)
						ws_tris(0,n_ws_tris) = n_ws_pts-5; // A
						ws_tris(1,n_ws_tris) = n_ws_pts-1; // F
						ws_tris(2,n_ws_tris) = n_ws_pts-2; // E
						n_ws_tris++;
						ws_tris(0,n_ws_tris) = n_ws_pts-5; // A
						ws_tris(1,n_ws_tris) = n_ws_pts-4; // B
						ws_tris(2,n_ws_tris) = n_ws_pts-1; // F
						n_ws_tris++;
						ws_tris(0,n_ws_tris) = n_ws_pts-4; // B
						ws_tris(1,n_ws_tris) = n_ws_pts-3; // D
						ws_tris(2,n_ws_tris) = n_ws_pts-1; // F
						n_ws_tris++;
					}
					else {
						// Points
						ws_pts(n_ws_pts++) = C; //-4
						ws_pts(n_ws_pts++) = D; //-3
						ws_pts(n_ws_pts++) = E; //-2
						ws_pts(n_ws_pts++) = F; //-1
												// Triangles
						ws_tris(0,n_ws_tris) = n_ws_pts-4; // C
						ws_tris(1,n_ws_tris) = n_ws_pts-1; // F
						ws_tris(2,n_ws_tris) = n_ws_pts-3; // D
						n_ws_tris++;
						ws_tris(0,n_ws_tris) = n_ws_pts-4; // C
						ws_tris(1,n_ws_tris) = n_ws_pts-2; // E
						ws_tris(2,n_ws_tris) = n_ws_pts-1; // F
						n_ws_tris++;
						// Points
						ns_pts(n_ns_pts++) = A; //-5
						ns_pts(n_ns_pts++) = B; //-4
						ns_pts(n_ns_pts++) = D; //-3
						ns_pts(n_ns_pts++) = E; //-2
						ns_pts(n_ns_pts++) = F; //-1
												// Triangles (A,E,F),(A,B,F),(B,D,F)
						ns_tris(0,n_ns_tris) = n_ns_pts-5; // A
						ns_tris(1,n_ns_tris) = n_ns_pts-1; // F
						ns_tris(2,n_ns_tris) = n_ns_pts-2; // E
						n_ns_tris++;
						ns_tris(0,n_ns_tris) = n_ns_pts-5; // A
						ns_tris(1,n_ns_tris) = n_ns_pts-4; // B
						ns_tris(2,n_ns_tris) = n_ns_pts-1; // F
						n_ns_tris++;
						ns_tris(0,n_ns_tris) = n_ns_pts-4; // B
						ns_tris(1,n_ns_tris) = n_ns_pts-3; // D
						ns_tris(2,n_ns_tris) = n_ns_pts-1; // F
						n_ns_tris++;
					}
				}
			}
		}
	}
}
//-------------------------------------------------------------------------------
void MC( DoubleArray &A, double &v, DoubleArray &solid, int &i, int &j, int &k,
				DTMutableList<Point> &nw_pts, int &n_nw_pts, IntArray &nw_tris,
				int &n_nw_tris)
{
	int N = 0;		// n will be the number of vertices in this grid cell only
    Point P;
	Point pt;
	DoubleArray TEST(3);

    Point PlaceHolder;
    int m;
    int o;
    int p;


	// Go over each corner -- check to see if the corners are themselves vertices
	//1
	if (A(i,j,k) == v){
		P.x = i;
		P.y = j;
		P.z = k;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//2
	if (A(i+1,j,k) == v){
		P.x = i+1;
		P.y = j;
		P.z = k;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//3
	if (A(i+1,j+1,k) == v){
		P.x = i+1;
		P.y = j+1;
		P.z = k;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//4
	if (A(i,j+1,k) == v){
		P.x = i;
		P.y = j+1;
		P.z = k;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//5
	if (A(i,j,k+1) == v){
		P.x = i;
		P.y = j;
		P.z = k+1;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//6
	if (A(i+1,j,k+1) == v){
		P.x = i+1;
		P.y = j;
		P.z = k+1;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//7
	if (A(i+1,j+1,k+1) == v){
		P.x = i+1;
		P.y = j+1;
		P.z = k+1;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//8
	if (A(i,j+1,k+1) == v){
		P.x = i;
		P.y = j+1;
		P.z = k+1;

		nw_pts(n_nw_pts++) = P;
		N++;
	}

    // Go through each side, compute P for sides of box spiraling up


//	float val;
	if ((A(i,j,k)-v)*(A(i+1,j,k)-v) < 0)
	{
		// If both points are in the fluid region
		if (A(i,j,k) != 0 && A(i+1,j,k) != 0){
			P.x = i + (A(i,j,k)-v)/(A(i,j,k)-A(i+1,j,k));
			P.y = j;
			P.z = k;
			// Evaluate the function S at the new point
			if (  solid(i,j,k)*(1-P.x+i) + solid(i+1,j,k)*(P.x-i) > 0 ){
				// This point is in the fluid region
				nw_pts(n_nw_pts++) =  P;
				N++;
			}
		}
	}
	// if point A(i,j,k) is in the solid phase
	//		else if ( A(i,j,k) == 0 ){
	//			pt.x = i;
	//			pt.y = j;
	//			pt.z = k;
	//			val =  EXTRAP(A, v, i+1,j,k, 1,pt);
	//			// If extrapolated value gives a vertex
	//			if ( (A(i+1,j,k)- v)*(val-v) < 0 ){
	//				P.x = i + (val-v)/(val-A(i+1,j,k));
	//				P.y = j;
	//				P.z = k;
	//				if (  INTERP(solid(i,j,k), solid(i+1,j,k), P.x-i) > 0 ){
	//					nw_pts(n_nw_pts++) =  P;
	//					N++;
	//				}
	//			}
	//		}
	//		// if point A(i+1,j,k) is in the solid phase
	//		else if ( A(i+1,j,k) == 0 ){
	//			pt.x = i+1;
	//			pt.y = j;
	//			pt.z = k;
	//			val =  EXTRAP(A, v, i,j,k, 4,pt);
	//			// If extrapolated value gives a vertex
	//			if ( (A(i,j,k)- v)*(val-v) < 0 ){
	//				P.x = i + (A(i,j,k)-v)/(A(i,j,k)-val);
	//				P.y = j;
	//				P.z = k;
	//				if (  INTERP(solid(i,j,k), solid(i+1,j,k), P.x-i) > 0 ){
	//					nw_pts(n_nw_pts++) =  P;
	//					N++;
	//				}
	//			}
	//		}
	//	}
    // 2
    if ((A(i+1,j,k)-v)*(A(i+1,j+1,k)-v) < 0)
	{
		if ( A(i+1,j,k) != 0 && A(i+1,j+1,k) != 0 ){
			P.x = i+1;
			P.y = j + (A(i+1,j,k)-v)/(A(i+1,j,k)-A(i+1,j+1,k));
			P.z = k;
			// Evaluate the function S at the new point
			if (  solid(i+1,j,k)*(1-P.y+j) + solid(i+1,j+1,k)*(P.y-j) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){ // P is a new vertex (not counted twice)
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
	//		else if ( A(i+1,j,k)  == 0 ){
	//			pt.x = i+1;
	//			pt.y = j;
	//			pt.z = k;
	//			val =  EXTRAP(A, v, i+1,j+1,k, 2,pt);
	//			// If extrapolated value gives a vertex
	//			if ( (A(i+1,j+1,k)- v)*(val-v) < 0 ){
	//				P.x = i+1;
	//				P.y = j + (val-v)/(val-A(i+1,j+1,k));
	//				P.z = k;
	//				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
	//					INTERP(solid(i+1,j,k), solid(i+1,j+1,k), P.y-j) > 0 ){ // P is a new vertex (not counted twice)
	//					nw_pts(n_nw_pts++) =  P;
	//					N++;
	//				}
	//			}
	//		}
	//		else if ( A(i+1,j+1,k) == 0){
	//			pt.x = i+1;
	//			pt.y = j+1;
	//			pt.z = k;
	//			val =  EXTRAP(A, v, i+1,j,k, 5,pt);
	//			// If extrapolated value gives a vertex
	//			if ( (A(i+1,j,k)- v)*(val-v) < 0 ){
	//				P.x = i+1;
	//				P.y = j + (A(i+1,j,k)-v)/(A(i+1,j,k)-val);
	//				P.z = k;
	//				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
	//					INTERP(solid(i+1,j,k), solid(i+1,j+1,k), P.y-j) > 0 ){
	//					nw_pts(n_nw_pts++) =  P;
	//					N++;
	//				}
	//			}
	//		}
	//	 }
    //3
    if ((A(i+1,j+1,k)-v)*(A(i,j+1,k)-v) < 0 )
	{
		if ( A(i+1,j+1,k) != 0 && A(i,j+1,k) != 0 ){
			P.x = i + (A(i,j+1,k)-v) / (A(i,j+1,k)-A(i+1,j+1,k));
			P.y = j+1;
			P.z = k;
			// Evaluate the function S at the new point
			if (  solid(i,j+1,k)*(1-P.x+i) + solid(i+1,j+1,k)*(P.x-i) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){ // P is a new vertex (not counted twice)
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
	/*		else if ( A(i,j+1,k)  == 0 ){
		pt.x = i;
	pt.y = j+1;
	pt.z = k;
	val =  EXTRAP(A, v, i+1,j+1,k, 1,pt);
	// If extrapolated value gives a vertex
	if ( (A(i+1,j+1,k)- v)*(val-v) < 0 ){
		P.x = i + (val-v) / (val-A(i+1,j+1,k));
		P.y = j+1;
		P.z = k;
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j+1,k), solid(i+1,j+1,k), P.x-i) > 0 ){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
	}
else if ( A(i+1,j+1,k) == 0){
	pt.x = i+1;
	pt.y = j+1;
	pt.z = k;
	val =  EXTRAP(A, v, i,j+1,k, 4, pt);
	// If extrapolated value gives a vertex
	if ( (A(i,j+1,k)- v)*(val-v) < 0 ){
		P.x = i + (A(i,j+1,k)-v) / (A(i,j+1,k)-val);
		P.y = j+1;
		P.z = k;
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j+1,k), solid(i+1,j+1,k), P.x-i) > 0 ){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
}
	 }
*/

    //4
    if ((A(i,j+1,k)-v)*(A(i,j,k)-v) < 0 )
	{
		if (A(i,j+1,k) != 0 && A(i,j,k) != 0 ){
			P.x = i;
			P.y = j + (A(i,j,k)-v) / (A(i,j,k)-A(i,j+1,k));
			P.z = k;
			// Evaluate the function S at the new point
			if (  solid(i,j,k)*(1-P.y+j) + solid(i,j+1,k)*(P.y-j) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){ // P is a new vertex (not counted twice)
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
	/*		else if ( A(i,j+1,k)  == 0 ){
		pt.x = i;
	pt.y = j+1;
	pt.z = k;
	val =  EXTRAP(A, v, i,j,k, 5,pt);
	// If extrapolated value gives a vertex
	if ( (A(i,j,k)- v)*(val-v) < 0 ){
		P.x = i;
		P.y = j + (A(i,j,k)-v) / (A(i,j,k)-val);
		P.z = k;
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j,k), solid(i,j+1,k), P.y-j) > 0 ){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
		}
	}
	}
else if ( A(i,j,k) == 0){
	pt.x = i;
	pt.y = j;
	pt.z = k;
	val =  EXTRAP(A, v, i,j+1,k, 2,pt);
	// If extrapolated value gives a vertex
	if ( (A(i,j+1,k)- v)*(val-v) < 0 ){
		P.x = i;
		P.y = j + (val-v) / (val-A(i,j+1,k));
		P.z = k;
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j,k), solid(i,j+1,k), P.y-j) > 0){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
}
	 }
*/
    //5
    if ((A(i,j,k)-v)*(A(i,j,k+1)-v) < 0 )
	{
		if ( A(i,j,k) != 0 && A(i,j,k+1) != 0 ){
			P.x = i;
			P.y = j;
			P.z = k + (A(i,j,k)-v) / (A(i,j,k)-A(i,j,k+1));
			// Evaluate the function S at the new point
			if (  solid(i,j,k)*(1-P.z+k) + solid(i,j,k+1)*(P.z-k) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){ // P is a new vertex (not counted twice)
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
	/*		else if ( A(i,j,k)  == 0 ){
		pt.x = i;
	pt.y = j;
	pt.z = k;
	val =  EXTRAP(A, v, i,j,k+1, 3,pt);
	// If extrapolated value gives a vertex
	if ( (A(i,j,k+1)- v)*(val-v) < 0 ){
		P.x = i;
		P.y = j;
		P.z = k + (val-v) / (val-A(i,j,k+1));
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j,k), solid(i,j,k+1), P.z-k) > 0){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
	}
else if ( A(i,j,k+1) == 0){
	pt.x = i;
	pt.y = j;
	pt.z = k+1;
	val =  EXTRAP(A, v, i,j,k, 6,pt);
	// If extrapolated value gives a vertex
	if ( (A(i,j,k)- v)*(val-v) < 0 ){
		P.x = i;
		P.y = j;
		P.z = k + (A(i,j,k)-v) / (A(i,j,k)-val);
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j,k), solid(i,j,k+1), P.z-k) > 0){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
}
	 }
*/
    //6
    if ((A(i+1,j,k)-v)*(A(i+1,j,k+1)-v) < 0 )
	{
		if ( A(i+1,j,k) != 0 && A(i+1,j,k+1) != 0 ){
			P.x = i+1;
			P.y = j;
			P.z = k + (A(i+1,j,k)-v) / (A(i+1,j,k)-A(i+1,j,k+1));
			// Evaluate the function S at the new point
			if (  solid(i+1,j,k)*(1-P.z+k) + solid(i+1,j,k+1)*(P.z-k) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
	/*		else if ( A(i+1,j,k)  == 0 ){
		pt.x = i+1;
	pt.y = j;
	pt.z = k;
	val =  EXTRAP(A, v, i+1,j,k+1, 3,pt);
	// If extrapolated value gives a vertex
	if ( (A(i+1,j,k+1)- v)*(val-v) < 0 ){
		P.x = i+1;
		P.y = j;
		P.z = k + (val-v) / (val-A(i+1,j,k+1));
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i+1,j,k), solid(i+1,j,k+1), P.z-k) > 0){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
	}
else if ( A(i+1,j,k+1) == 0){
	pt.x = i+1;
	pt.y = j;
	pt.z = k+1;
	val =  EXTRAP(A, v, i+1,j,k, 6,pt);
	// If extrapolated value gives a vertex
	if ( (A(i+1,j,k)- v)*(val-v) < 0 ){
		P.x = i+1;
		P.y = j;
		P.z = k + (A(i+1,j,k)-v) / (A(i+1,j,k)-val);
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i+1,j,k), solid(i+1,j,k+1), P.z-k) > 0){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
}

	 }
*/
    //7
    if ((A(i+1,j+1,k)-v)*(A(i+1,j+1,k+1)-v) < 0 )
	{
		if ( A(i+1,j+1,k) != 0 && A(i+1,j+1,k+1) != 0 ){
			P.x = i+1;
			P.y = j+1;
			P.z = k + (A(i+1,j+1,k)-v) / (A(i+1,j+1,k)-A(i+1,j+1,k+1));
			// Evaluate the function S at the new point
			if (  solid(i+1,j+1,k)*(1-P.z+k) + solid(i+1,j+1,k+1)*(P.z-k) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
	/*		else if ( A(i+1,j+1,k)  == 0 ){
		pt.x = i+1;
	pt.y = j+1;
	pt.z = k;
	val =  EXTRAP(A, v, i+1,j+1,k+1, 3,pt);
	// If extrapolated value gives a vertex
	if ( (A(i+1,j+1,k+1)- v)*(val-v) < 0 ){
		P.x = i+1;
		P.y = j+1;
		P.z = k + (val-v) / (val-A(i+1,j+1,k+1));
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i+1,j+1,k), solid(i+1,j+1,k+1), P.z-k) > 0 ){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
	}
else if ( A(i+1,j+1,k+1) == 0){
	pt.x = i+1;
	pt.y = j+1;
	pt.z = k+1;
	val =  EXTRAP(A, v, i+1,j+1,k, 6,pt);
	// If extrapolated value gives a vertex
	if ( (A(i+1,j+1,k)- v)*(val-v) < 0 ){
		P.x = i+1;
		P.y = j+1;
		P.z = k + (A(i+1,j+1,k)-v) / (A(i+1,j+1,k)-val);
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i+1,j+1,k), solid(i+1,j+1,k+1), P.z-k) > 0 ){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
}
	 }
*/
    //8
    if ((A(i,j+1,k)-v)*(A(i,j+1,k+1)-v) < 0 )
	{
		if ( A(i,j+1,k) != 0 && A(i,j+1,k+1) != 0 ){
			P.x = i;
			P.y = j+1;
			P.z = k + (A(i,j+1,k)-v) / (A(i,j+1,k)-A(i,j+1,k+1));
			// Evaluate the function S at the new point
			if (  solid(i,j+1,k)*(1-P.z+k) + solid(i,j+1,k+1)*(P.z-k) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
	/*		else if ( A(i,j+1,k)  == 0 ){
		pt.x = i;
	pt.y = j+1;
	pt.z = k;
	val =  EXTRAP(A, v, i,j+1,k+1, 3,pt);
	// If extrapolated value gives a vertex
	if ( (A(i,j+1,k+1)- v)*(val-v) < 0 ){
		P.x = i;
		P.y = j+1;
		P.z = k + (val-v) / (val-A(i,j+1,k+1));
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j+1,k), solid(i,j+1,k+1), P.z-k) > 0 ){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
	}
else if ( A(i,j+1,k+1) == 0){
	pt.x = i;
	pt.y = j+1;
	pt.z = k+1;
	val =  EXTRAP(A, v, i,j+1,k, 6,pt);
	// If extrapolated value gives a vertex
	if ( (A(i,j+1,k)- v)*(val-v) < 0 ){
		P.x = i;
		P.y = j+1;
		P.z = k + (A(i,j+1,k)-v) / (A(i,j+1,k)-val);
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j+1,k), solid(i,j+1,k+1), P.z-k) > 0 ){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
}
	 }
*/
    //9
    if ((A(i,j,k+1)-v)*(A(i+1,j,k+1)-v) < 0 )
	{
		if ( A(i,j,k+1) != 0 && A(i+1,j,k+1) != 0 ){
			P.x = i + (A(i,j,k+1)-v) / (A(i,j,k+1)-A(i+1,j,k+1));
			P.y = j;
			P.z = k+1;
			// Evaluate the function S at the new point
			if (  solid(i,j,k+1)*(1-P.x+i) + solid(i+1,j,k+1)*(P.x-i) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
	/*		else if ( A(i,j,k+1)  == 0 ){
		pt.x = i;
	pt.y = j;
	pt.z = k+1;
	val =  EXTRAP(A, v, i+1,j,k+1, 1,pt);
	// If extrapolated value gives a vertex
	if ( (A(i+1,j,k+1)- v)*(val-v) < 0 ){
		P.x = i + (val-v) / (val-A(i+1,j,k+1));
		P.y = j;
		P.z = k+1;
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j,k+1), solid(i+1,j,k+1), P.x-i) > 0 ){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
	}
else if ( A(i+1,j,k+1) == 0){
	pt.x = i+1;
	pt.y = j;
	pt.z = k+1;
	val =  EXTRAP(A, v, i,j,k+1, 4,pt);
	// If extrapolated value gives a vertex
	if ( (A(i,j,k+1)- v)*(val-v) < 0 ){
		P.x = i + (A(i,j,k+1)-v) / (A(i,j,k+1)-val);
		P.y = j;
		P.z = k+1;
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j,k+1), solid(i+1,j,k+1), P.x-i) > 0){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
}
	 }
*/
    //10
    if ((A(i+1,j,k+1)-v)*(A(i+1,j+1,k+1)-v) < 0 )
	{
		if ( A(i+1,j,k+1) != 0 && A(i+1,j+1,k+1) != 0 ){
			P.x = i+1;
			P.y = j + (A(i+1,j,k+1)-v) / (A(i+1,j,k+1)-A(i+1,j+1,k+1));
			P.z = k+1;
			// Evaluate the function S at the new point
			if (  solid(i+1,j,k+1)*(1-P.y+j) + solid(i+1,j+1,k+1)*(P.y-j) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
	/*		else if ( A(i+1,j,k+1)  == 0 ){
		pt.x = i+1;
	pt.y = j;
	pt.z = k+1;
	val =  EXTRAP(A, v, i+1,j+1,k+1, 2,pt);
	// If extrapolated value gives a vertex
	if ( (A(i+1,j+1,k+1)- v)*(val-v) < 0 ){
		P.x = i+1;
		P.y = j + (val-v) / (val-A(i+1,j+1,k+1));
		P.z = k+1;
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i+1,j,k+1), solid(i+1,j+1,k+1), P.y-j) > 0){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
	}
else if ( A(i+1,j+1,k+1) == 0){
	pt.x = i+1;
	pt.y = j+1;
	pt.z = k+1;
	val =  EXTRAP(A, v, i+1,j,k+1, 5,pt);
	// If extrapolated value gives a vertex
	if ( (A(i+1,j,k+1)- v)*(val-v) < 0 ){
		P.x = i+1;
		P.y = j + (A(i+1,j,k+1)-v) / (A(i+1,j,k+1)-val);
		P.z = k+1;
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i+1,j,k+1), solid(i+1,j+1,k+1), P.y-j) > 0){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
}
	 }
*/
    //11
    if ((A(i+1,j+1,k+1)-v)*(A(i,j+1,k+1)-v) < 0 )
	{
		if ( A(i+1,j+1,k+1) != 0 && A(i,j+1,k+1) != 0 ){
			P.x = i+(A(i,j+1,k+1)-v) / (A(i,j+1,k+1)-A(i+1,j+1,k+1));
			P.y = j+1;
			P.z = k+1;
			// Evaluate the function S at the new point
			if (  solid(i,j+1,k+1)*(1-P.x+i) + solid(i+1,j+1,k+1)*(P.x-i) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
	/*		else if ( A(i+1,j+1,k+1)  == 0 ){
		pt.x = i+1;
	pt.y = j+1;
	pt.z = k+1;
	val =  EXTRAP(A, v, i,j+1,k+1, 4,pt);
	// If extrapolated value gives a vertex
	if ( (A(i,j+1,k+1)- v)*(val-v) < 0 ){
		P.x = i+(A(i,j+1,k+1)-v) / (A(i,j+1,k+1)-val);
		P.y = j+1;
		P.z = k+1;
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j+1,k+1), solid(i+1,j+1,k+1), P.x-i) > 0){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
	}
else if ( A(i,j+1,k+1) == 0){
	pt.x = i;
	pt.y = j+1;
	pt.z = k+1;
	val =  EXTRAP(A, v, i+1,j+1,k+1, 1,pt);
	// If extrapolated value gives a vertex
	if ( (A(i+1,j+1,k+1)- v)*(val-v) < 0 ){
		P.x = i+(val-v) / (val-A(i+1,j+1,k+1));
		P.y = j+1;
		P.z = k+1;
		if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1 &&
			INTERP(solid(i,j+1,k+1), solid(i+1,j+1,k+1), P.x-i) > 0){ // P is a new vertex (not counted twice)
			nw_pts(n_nw_pts++) =  P;
			N++;
		}
	}
}
	 }
*/
    //12
    if ((A(i,j+1,k+1)-v)*(A(i,j,k+1)-v) < 0 )
	{
		if ( A(i,j+1,k+1) != 0 && A(i,j,k+1) != 0 ){
			P.x = i;
			P.y = j + (A(i,j,k+1)-v) / (A(i,j,k+1)-A(i,j+1,k+1));
			P.z = k+1;
			// Evaluate the function S at the new point
			if (  solid(i,j,k+1)*(1-P.y+j) + solid(i,j+1,k+1)*(P.y-j) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}

    // sort vertices so that they are connected to "neighbors"

    // DTMatlabDataFile("/tmp/Dump.mat",DTFile::NewReadWrite);

    // n_nw_pts = number of vertices in total (location n_nw_pts is first unfilled position)
    // n = number of vertices in this grid cell


    for (m = n_nw_pts-N; m < n_nw_pts-2; m++) {
        for (o = m+2; o < n_nw_pts-1; o++) {
            if (ShareSide(nw_pts(m), nw_pts(o)) == 1) {
                PlaceHolder = nw_pts(m+1);
                nw_pts(m+1) = nw_pts(o);
                nw_pts(o) = PlaceHolder;
            }
        }

		// make sure other neighbor of vertex 1 is in last spot
        if (m == n_nw_pts-N){
            for (p = m+2; p < n_nw_pts-1; p++){
                if (ShareSide(nw_pts(m), nw_pts(p)) == 1){
                    PlaceHolder = nw_pts(n_nw_pts-1);
                    nw_pts(n_nw_pts-1) = nw_pts(p);
                    nw_pts(p) = PlaceHolder;
                }
            }
        }
        if ( ShareSide(nw_pts(n_nw_pts-2), nw_pts(n_nw_pts-3)) != 1 ){
            if (ShareSide( nw_pts(n_nw_pts-3), nw_pts(n_nw_pts-1)) == 1 &&
				ShareSide( nw_pts(n_nw_pts-N),nw_pts(n_nw_pts-2)) == 1 ){
                PlaceHolder = nw_pts(n_nw_pts-2);
                nw_pts(n_nw_pts-2) = nw_pts(n_nw_pts-1);
                nw_pts(n_nw_pts-1) = PlaceHolder;
            }
        }
        if ( ShareSide(nw_pts(n_nw_pts-1), nw_pts(n_nw_pts-2)) != 1 ){
            if (ShareSide( nw_pts(n_nw_pts-3), nw_pts(n_nw_pts-1)) == 1 &&
				ShareSide(nw_pts(n_nw_pts-4),nw_pts(n_nw_pts-2)) == 1 ){
                PlaceHolder = nw_pts(n_nw_pts-3);
                nw_pts(n_nw_pts-3) = nw_pts(n_nw_pts-2);
                nw_pts(n_nw_pts-2) = PlaceHolder;
            }
            if (ShareSide( nw_pts(n_nw_pts-N+1), nw_pts(n_nw_pts-3)) == 1 &&
				ShareSide(nw_pts(n_nw_pts-1),nw_pts(n_nw_pts-N+1)) == 1 ){
                PlaceHolder = nw_pts(n_nw_pts-2);
                nw_pts(n_nw_pts-2) = nw_pts(n_nw_pts-N+1);
                nw_pts(n_nw_pts-N+1) = PlaceHolder;
            }
        }
        if ( ShareSide(nw_pts(n_nw_pts-N), nw_pts(n_nw_pts-N+1)) != 1 ){
            if (ShareSide( nw_pts(n_nw_pts-N), nw_pts(n_nw_pts-2)) == 1 &&
				ShareSide(nw_pts(n_nw_pts-1), nw_pts(n_nw_pts-N+1)) == 1){
                PlaceHolder = nw_pts(n_nw_pts-1);
                nw_pts(n_nw_pts-1) = nw_pts(n_nw_pts-N);
                nw_pts(n_nw_pts-N) = PlaceHolder;
            }
        }
    }

	// *    *    *   ESTABLISH TRIANGLE CONNECTIONS	   *    *    *

	for (p=n_nw_pts-N+2; p<n_nw_pts; p++){
		nw_tris(0,n_nw_tris) = n_nw_pts-N;
        nw_tris(1,n_nw_tris) = p-1;
        nw_tris(2,n_nw_tris) = p;
        n_nw_tris++;
	}

}
//-------------------------------------------------------------------------------
void EDGE(DoubleArray &A, double &v, DoubleArray &solid, int &i, int &j, int &k, int &m, int &n, int &o,
				 DTMutableList<Point> &nw_pts, int &n_nw_pts, IntArray &nw_tris, int &n_nw_tris,
				 DTMutableList<Point> &local_nws_pts, int &n_local_nws_pts)
{
	// FIND THE POINTS ON THE nw SURFACE THAT ARE ON THE EDGE (COMMON LINE WITH SOLID PHASE)
	// function A is the fluid data padded (so that it has values inside the solid phase)

	int N = 0;		// n will be the number of vertices in this grid cell only
    Point P;

    Point temp;

    int p; int q; int r;

	// Add common line points to nw_pts
	for (p=0;p<n_local_nws_pts;p++){
		nw_pts(n_nw_pts++) = local_nws_pts(p);
	}


	// Go over each corner -- check to see if the corners are themselves vertices
	//1
	if (A(i,j,k) == v){
		P.x = i;
		P.y = j;
		P.z = k;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//2
	if (A(i+1,j,k) == v){
		P.x = i+1;
		P.y = j;
		P.z = k;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//3
	if (A(i+1,j+1,k) == v){
		P.x = i+1;
		P.y = j+1;
		P.z = k;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//4
	if (A(i,j+1,k) == v){
		P.x = i;
		P.y = j+1;
		P.z = k;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//5
	if (A(i,j,k+1) == v){
		P.x = i;
		P.y = j;
		P.z = k+1;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//6
	if (A(i+1,j,k+1) == v){
		P.x = i+1;
		P.y = j;
		P.z = k+1;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//7
	if (A(i+1,j+1,k+1) == v){
		P.x = i+1;
		P.y = j+1;
		P.z = k+1;
		nw_pts(n_nw_pts++) = P;
		N++;
	}
	//8
	if (A(i,j+1,k+1) == v){
		P.x = i;
		P.y = j+1;
		P.z = k+1;

		nw_pts(n_nw_pts++) = P;
		N++;
	}

    // Go through each side, compute P for sides of box spiraling up
//	float val;
	Point pt;

	// 1
    if ((A(i,j,k)-v)*(A(i+1,j,k)-v) < 0)
	{
		// If both points are in the fluid region
		if (A(i,j,k) != 0 && A(i+1,j,k) != 0){
			P.x = i + (A(i,j,k)-v)/(A(i,j,k)-A(i+1,j,k));
			P.y = j;
			P.z = k;
			// Evaluate the function S at the new point
			if (  solid(i,j,k)*(1-P.x+i) + solid(i+1,j,k)*(P.x-i) > 0 ){
				// This point is in the fluid region
				nw_pts(n_nw_pts++) =  P;
				N++;
			}
		}
	}

    // 2
    if ((A(i+1,j,k)-v)*(A(i+1,j+1,k)-v) < 0)
	{
		if ( A(i+1,j,k) != 0 && A(i+1,j+1,k) != 0 ){
			P.x = i+1;
			P.y = j + (A(i+1,j,k)-v)/(A(i+1,j,k)-A(i+1,j+1,k));
			P.z = k;
			// Evaluate the function S at the new point
			if (  solid(i+1,j,k)*(1-P.y+j) + solid(i+1,j+1,k)*(P.y-j) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){ // P is a new vertex (not counted twice)
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}

    //3
    if ((A(i+1,j+1,k)-v)*(A(i,j+1,k)-v) < 0 )
	{
		if ( A(i+1,j+1,k) != 0 && A(i,j+1,k) != 0 ){
			P.x = i + (A(i,j+1,k)-v) / (A(i,j+1,k)-A(i+1,j+1,k));
			P.y = j+1;
			P.z = k;
			// Evaluate the function S at the new point
			if (  solid(i,j+1,k)*(1-P.x+i) + solid(i+1,j+1,k)*(P.x-i) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){ // P is a new vertex (not counted twice)
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}

    //4
    if ((A(i,j+1,k)-v)*(A(i,j,k)-v) < 0 )
	{
		if (A(i,j+1,k) != 0 && A(i,j,k) != 0 ){
			P.x = i;
			P.y = j + (A(i,j,k)-v) / (A(i,j,k)-A(i,j+1,k));
			P.z = k;
			// Evaluate the function S at the new point
			if (  solid(i,j,k)*(1-P.y+j) + solid(i,j+1,k)*(P.y-j) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){ // P is a new vertex (not counted twice)
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
    //5
    if ((A(i,j,k)-v)*(A(i,j,k+1)-v) < 0 )
	{
		if ( A(i,j,k) != 0 && A(i,j,k+1) != 0 ){
			P.x = i;
			P.y = j;
			P.z = k + (A(i,j,k)-v) / (A(i,j,k)-A(i,j,k+1));
			// Evaluate the function S at the new point
			if (  solid(i,j,k)*(1-P.z+k) + solid(i,j,k+1)*(P.z-k) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){ // P is a new vertex (not counted twice)
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}

    //6
    if ((A(i+1,j,k)-v)*(A(i+1,j,k+1)-v) < 0 )
	{
		if ( A(i+1,j,k) != 0 && A(i+1,j,k+1) != 0 ){
			P.x = i+1;
			P.y = j;
			P.z = k + (A(i+1,j,k)-v) / (A(i+1,j,k)-A(i+1,j,k+1));
			// Evaluate the function S at the new point
			if (  solid(i+1,j,k)*(1-P.z+k) + solid(i+1,j,k+1)*(P.z-k) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}
    //7
    if ((A(i+1,j+1,k)-v)*(A(i+1,j+1,k+1)-v) < 0 )
	{
		if ( A(i+1,j+1,k) != 0 && A(i+1,j+1,k+1) != 0 ){
			P.x = i+1;
			P.y = j+1;
			P.z = k + (A(i+1,j+1,k)-v) / (A(i+1,j+1,k)-A(i+1,j+1,k+1));
			// Evaluate the function S at the new point
			if (  solid(i+1,j+1,k)*(1-P.z+k) + solid(i+1,j+1,k+1)*(P.z-k) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}

    //8
    if ((A(i,j+1,k)-v)*(A(i,j+1,k+1)-v) < 0 )
	{
		if ( A(i,j+1,k) != 0 && A(i,j+1,k+1) != 0 ){
			P.x = i;
			P.y = j+1;
			P.z = k + (A(i,j+1,k)-v) / (A(i,j+1,k)-A(i,j+1,k+1));
			// Evaluate the function S at the new point
			if (  solid(i,j+1,k)*(1-P.z+k) + solid(i,j+1,k+1)*(P.z-k) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}

    //9
    if ((A(i,j,k+1)-v)*(A(i+1,j,k+1)-v) < 0 )
	{
		if ( A(i,j,k+1) != 0 && A(i+1,j,k+1) != 0 ){
			P.x = i + (A(i,j,k+1)-v) / (A(i,j,k+1)-A(i+1,j,k+1));
			P.y = j;
			P.z = k+1;
			// Evaluate the function S at the new point
			if (  solid(i,j,k+1)*(1-P.x+i) + solid(i+1,j,k+1)*(P.x-i) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}

    //10
    if ((A(i+1,j,k+1)-v)*(A(i+1,j+1,k+1)-v) < 0 )
	{
		if ( A(i+1,j,k+1) != 0 && A(i+1,j+1,k+1) != 0 ){
			P.x = i+1;
			P.y = j + (A(i+1,j,k+1)-v) / (A(i+1,j,k+1)-A(i+1,j+1,k+1));
			P.z = k+1;
			// Evaluate the function S at the new point
			if (  solid(i+1,j,k+1)*(1-P.y+j) + solid(i+1,j+1,k+1)*(P.y-j) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}

    //11
    if ((A(i+1,j+1,k+1)-v)*(A(i,j+1,k+1)-v) < 0 )
	{
		if ( A(i+1,j+1,k+1) != 0 && A(i,j+1,k+1) != 0 ){
			P.x = i+(A(i,j+1,k+1)-v) / (A(i,j+1,k+1)-A(i+1,j+1,k+1));
			P.y = j+1;
			P.z = k+1;
			// Evaluate the function S at the new point
			if (  solid(i,j+1,k+1)*(1-P.x+i) + solid(i+1,j+1,k+1)*(P.x-i) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}

    //12
    if ((A(i,j+1,k+1)-v)*(A(i,j,k+1)-v) < 0 )
	{
		if ( A(i,j+1,k+1) != 0 && A(i,j,k+1) != 0 ){
			P.x = i;
			P.y = j + (A(i,j,k+1)-v) / (A(i,j,k+1)-A(i,j+1,k+1));
			P.z = k+1;
			// Evaluate the function S at the new point
			if (  solid(i,j,k+1)*(1-P.y+j) + solid(i,j+1,k+1)*(P.y-j) > 0 ){
				// This point is in the fluid region
				if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1){
					nw_pts(n_nw_pts++) =  P;
					N++;
				}
			}
		}
	}

	////////////////////////////////////////////////
	////// SORT VERTICES SO THAT THEY CONNECT //////
	//////        TO ALL "NEIGHBORS"          //////
	////////////////////////////////////////////////

	//  First common line point should connect to last MC point
	for (q=n_nw_pts-N; q<n_nw_pts-1; q++){
		if ( ShareSide(nw_pts(n_nw_pts-N-n_local_nws_pts), nw_pts(n_nw_pts-1)) == 0 &&
			 ShareSide(nw_pts(n_nw_pts-N-n_local_nws_pts), nw_pts(q)) == 1 ){
			// Switch point q with point n_nw_pts-N-n_local_nws_pts
			temp = nw_pts(n_nw_pts-1);
			nw_pts(n_nw_pts-1) = nw_pts(q);
			nw_pts(q) = temp;
		}
	}

	//  Last common line point should connect to first MC point
	for (q=n_nw_pts-N+1; q<n_nw_pts-1; q++){
		if ( ShareSide(nw_pts(n_nw_pts-N-1), nw_pts(n_nw_pts-N)) == 0 &&
			 ShareSide(nw_pts(n_nw_pts-N-1), nw_pts(q)) == 1 ){
			// Switch these points
			temp = nw_pts(n_nw_pts-N);
			nw_pts(n_nw_pts-N) = nw_pts(q);
			nw_pts(q) = temp;
		}
	}
	// All MC points should connect to their neighbors
	for (q=n_nw_pts-N; q<n_nw_pts-2; q++){
		if ( ShareSide(nw_pts(q), nw_pts(q+1)) == 0){
			for (r=q+2; r < n_nw_pts-1; r++){
				if ( ShareSide(nw_pts(q), nw_pts(q+1)) == 0 &&
					 ShareSide(nw_pts(q), nw_pts(r)) == 1){
					// Switch r and q+1
					temp = nw_pts(q+1);
					nw_pts(q+1) = nw_pts(r);
					nw_pts(r) = temp;
				}
			}
		}
	}


	// *    *    *   ESTABLISH TRIANGLE CONNECTIONS	   *    *    *
	for (p=n_nw_pts-N-n_local_nws_pts; p<n_nw_pts-2; p++){
		nw_tris(0,n_nw_tris) = n_nw_pts-1;
        nw_tris(1,n_nw_tris) = p;
        nw_tris(2,n_nw_tris) = p+1;
        n_nw_tris++;
	}

	//	for (p=n_nw_pts-N-n_local_nws_pts+2; p<n_nw_pts; p++){
	//		nw_tris(0,n_nw_tris) = n_nw_pts-N-n_local_nws_pts;
	//		nw_tris(1,n_nw_tris) = p-1;
	//		nw_tris(2,n_nw_tris) = p;
	//        n_nw_tris++;
	//	}

}
//--------------------------------------------------------------------------------------------------------
void ComputeAreasPMMC(IntArray &cubeList, int start, int finish,
							 DoubleArray &F, DoubleArray &S, double vF, double vS, 
							 double &blob_volume, double &ans, double &aws, double &awn, double &lwns,
							 int Nx, int Ny, int Nz)
{
	/* ****************************************************************
	 VARIABLES FOR THE PMMC ALGORITHM
	 ****************************************************************** */
	awn = aws = ans = lwns = 0.0;
	
//	bool add=1;			// Set to false if any corners contain nw-phase ( F > vF)
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
	
	int i,j,k,p,q,r;
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
	
	int n_nw_tris_beg, n_ns_tris_beg, n_ws_tris_beg;
	int c;
	int newton_steps = 0;
//	double blob_volume;
	
	/* ****************************************************************
	 RUN PMMC ON EACH BLOB
	 ****************************************************************** */
//	printf("Running the PMMC Algorithm \n");
//	printf("The number of blobs is %i \n",nblobs);

	// Store beginning points for surfaces for blob p
	n_nw_tris_beg = n_nw_tris;
	n_ns_tris_beg = n_ns_tris;
	n_ws_tris_beg = n_ws_tris;
//	n_nws_seg_beg = n_nws_seg;
	// Loop over all cubes
	blob_volume = 0;	// Initialize the volume for blob a to zero
	for (c=start;c<finish;c++){
		// Get cube from the list
		i = cubeList(0,c);
		j = cubeList(1,c);
		k = cubeList(2,c);

	//	printf("somewhere inside %i \n",c);

		
/*		// Compute the volume
		for (p=0;p<8;p++){
			if ( indicator(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > -1){
				blob_volume += 0.125;
			}
		}
*/		
		for (p=0;p<8;p++){
			if ( F(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 
				&&  S(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){
				blob_volume += 0.125;
			}
		}
		
		// Run PMMC
		n_local_sol_tris = 0;
		n_local_sol_pts = 0;
		n_local_nws_pts = 0;
		
		// if there is a solid phase interface in the grid cell
		if (Interface(S,vS,i,j,k) == 1){
			
			/////////////////////////////////////////
			/// CONSTRUCT THE LOCAL SOLID SURFACE ///
			/////////////////////////////////////////
			
			// find the local solid surface
			SOL_SURF(S,0,F,vF,i,j,k, Nx,Ny,Nz,local_sol_pts,n_local_sol_pts,
					 local_sol_tris,n_local_sol_tris,values);
			
			/////////////////////////////////////////
			//////// TRIM THE SOLID SURFACE /////////
			/////////////////////////////////////////
			TRIM(local_sol_pts, n_local_sol_pts, vF,local_sol_tris, n_local_sol_tris,
				 ns_pts, n_ns_pts, ns_tris, n_ns_tris, ws_pts, n_ws_pts,
				 ws_tris, n_ws_tris, values, local_nws_pts, n_local_nws_pts,
				 F, S, i, j, k, newton_steps);
			
			/////////////////////////////////////////
			//////// WRITE COMMON LINE POINTS ///////
			////////      TO MAIN ARRAYS      ///////
			/////////////////////////////////////////
			map = n_nws_pts;
			for (p=0; p < n_local_nws_pts; p++){
				nws_pts(n_nws_pts++) = local_nws_pts(p);
			}
			for (q=0; q < n_local_nws_pts-1; q++){
				nws_seg(0,n_nws_seg) = map+q;
				nws_seg(1,n_nws_seg) = map+q+1;
				n_nws_seg++;
			}
			
			/////////////////////////////////////////
			////// CONSTRUCT THE nw SURFACE /////////
			/////////////////////////////////////////
			if ( n_local_nws_pts > 0){
				EDGE(F, vF, S, i,j,k, Nx, Ny, Nz, nw_pts, n_nw_pts, nw_tris, n_nw_tris,
					 local_nws_pts, n_local_nws_pts);
			}
			else {
				MC(F, vF, S, i,j,k, nw_pts, n_nw_pts, nw_tris, n_nw_tris);
			}
		}
		
		/////////////////////////////////////////
		////// CONSTRUCT THE nw SURFACE /////////
		/////////////////////////////////////////
		
		else if (Fluid_Interface(F,S,vF,i,j,k) == 1){
			MC(F, vF, S, i,j,k, nw_pts, n_nw_pts, nw_tris, n_nw_tris);
		}
		//******END OF BLOB PMMC*********************************************

		//*******************************************************************
		// Compute the Interfacial Areas, Common Line length for blob p
		// nw surface
		for (r=n_nw_tris_beg;r<n_nw_tris;r++){
			A = nw_pts(nw_tris(0,r));
			B = nw_pts(nw_tris(1,r));
			C = nw_pts(nw_tris(2,r));
			// Compute length of sides (assume dx=dy=dz)
			s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
			s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
			s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
			s = 0.5*(s1+s2+s3);
			awn=awn+sqrt(s*(s-s1)*(s-s2)*(s-s3));
		}
		for (r=n_ns_tris_beg;r<n_ns_tris;r++){
			A = ns_pts(ns_tris(0,r));
			B = ns_pts(ns_tris(1,r));
			C = ns_pts(ns_tris(2,r));
			// Compute length of sides (assume dx=dy=dz)
			s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
			s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
			s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
			s = 0.5*(s1+s2+s3);
			ans=ans+sqrt(s*(s-s1)*(s-s2)*(s-s3));
		}
		for (r=n_ws_tris_beg;r<n_ws_tris;r++){
			A = ws_pts(ws_tris(0,r));
			B = ws_pts(ws_tris(1,r));
			C = ws_pts(ws_tris(2,r));
			// Compute length of sides (assume dx=dy=dz)
			s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
			s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
			s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
			s = 0.5*(s1+s2+s3);
			aws=aws+sqrt(s*(s-s1)*(s-s2)*(s-s3));
		}
		//*******************************************************************
		// Reset the triangle counts to zero
		n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
		n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
		
		n_nw_tris_beg = n_nw_tris;
		n_ns_tris_beg = n_ns_tris;
		n_ws_tris_beg = n_ws_tris;
//		n_nws_seg_beg = n_nws_seg;
		//*******************************************************************
	}
}

//--------------------------------------------------------------------------------------------------------

