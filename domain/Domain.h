// Created by James McClure
// Copyright 2008-2013
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

using namespace std;

inline void ReadSpherePacking(int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad)
{
	// Read in the full sphere pack
	int count;  
	double x,y,z,r;
	//...... READ IN THE SPHERES...................................
	char * trsh;
	trsh = new char[100];
	cout << "Reading the packing file..." << endl;
	ifstream pack ("pack.out");
	//.........Trash the header lines..........
	pack.getline(trsh, 100);
	pack.getline(trsh, 100);
	pack.getline(trsh, 100);
	pack.getline(trsh, 100);
	pack.getline(trsh, 100);
	//........read the spheres..................
	count = 0;
	pack >> x;
	pack >> y;
	pack >> z;
	pack >> r;
	while (! pack.eof()){
		List_cx[count] = x;
		List_cy[count] = y;
		List_cz[count] = z;
		List_rad[count] = r;
		pack >> x;
		pack >> y;
		pack >> z;
		pack >> r;
		count++;
	}
	pack.close();
	cout << "Number of spheres extracted is: " << count << endl;
	if (count != nspheres){
		printf("Specified number of spheres is probably incorrect!");
	}
	// .............................................................
}

inline void AssignLocalSolidID(char *ID, int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad,
						  double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, 
						  int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz)
{
	// Use sphere lists to determine which nodes are in porespace
	// Write out binary file for nodes
	char value;
	int N = Nx*Ny*Nz; 	// Domain size, including the halo
	double hx,hy,hz;
	double x,y,z;
	double cx,cy,cz,r;
	int imin,imax,jmin,jmax,kmin,kmax;
	int p,i,j,k,n;
	//............................................
	double min_x,min_y,min_z;
//	double max_x,max_y,max_z;
	//............................................
	// Lattice spacing for the entire domain
	// It should generally be true that hx=hy=hz
	// Otherwise, you will end up with ellipsoids
	hx = Lx/(Nx*nprocx-1);
	hy = Ly/(Ny*nprocy-1);
	hz = Lz/(Nz*nprocz-1);	
	//............................................
	// Get maximum and minimum for this domain
	// Halo is included !
	min_x = double(iproc*Nx-1)*hx;
	min_y = double(jproc*Ny-1)*hy;
	min_z = double(kproc*Nz-1)*hz;
//	max_x = ((iproc+1)*Nx+1)*hx;
//	max_y = ((jproc+1)*Ny+1)*hy;
//	max_z = ((kproc+1)*Nz+1)*hz;
	//............................................

	//............................................
		// Pre-initialize local ID 
	for (n=0;n<N;n++){
		ID[n]=1;
	}
	//............................................

	//............................................
	// .........Loop over the spheres.............
	for (p=0;p<nspheres;p++){
		// Get the sphere from the list, map to local min
		cx = List_cx[p] - min_x;
		cy = List_cy[p] - min_y;
		cz = List_cz[p] - min_z;
		r = List_rad[p];
		// Check if
		// Range for this sphere in global indexing
		imin = int ((cx-r)/hx)-1;
		imax = int ((cx+r)/hx)+1;
		jmin = int ((cy-r)/hy)-1;
		jmax = int ((cy+r)/hy)+1;
		kmin = int ((cz-r)/hz)-1;
		kmax = int ((cz+r)/hz)+1;
		// Obviously we have to do something at the edges
		if (imin<0)		imin = 0;
		if (imin>Nx)	imin = Nx;
		if (imax<0)		imax = 0;
		if (imax>Nx)	imax = Nx;
		if (jmin<0)		jmin = 0;
		if (jmin>Ny)	jmin = Ny;
		if (jmax<0)		jmax = 0;
		if (jmax>Ny)	jmax = Ny;
		if (kmin<0)		kmin = 0;
		if (kmin>Nz)	kmin = Nz;
		if (kmax<0)		kmax = 0;
		if (kmax>Nz)	kmax = Nz;
		// Loop over the domain for this sphere (may be null)
		for (i=imin;i<imax;i++){
			for (j=jmin;j<jmax;j++){
				for (k=kmin;k<kmax;k++){
					// Initialize ID value to 'fluid (=1)'
					x = i*hx;
					y = j*hy;
					z = k*hz;
					value = 1;
					// if inside sphere, set to zero
					if ( (cx-x)*(cx-x)+(cy-y)*(cy-y)+(cz-z)*(cz-z) < r*r){
						value=0;
					}
					// get the position in the list
					n = k*Nx*Ny+j*Nx+i;
					if ( ID[n] != 0 ){
						ID[n] = value;
					}
				}
			}
		}
	}
}

inline void SignedDistance(double *Distance, int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad,
						  double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, 
						  int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz)
{
	// Use sphere lists to determine which nodes are in porespace
	// Write out binary file for nodes
	int N = Nx*Ny*Nz; 	// Domain size, including the halo
	double hx,hy,hz;
	double x,y,z;
	double cx,cy,cz,r;
	int imin,imax,jmin,jmax,kmin,kmax;
	int p,i,j,k,n;
	//............................................
	double min_x,min_y,min_z;
	double distance;
	//............................................
	// Lattice spacing for the entire domain
	// It should generally be true that hx=hy=hz
	// Otherwise, you will end up with ellipsoids
	hx = Lx/((Nx-2)*nprocx-1);
	hy = Ly/((Ny-2)*nprocy-1);
	hz = Lz/((Nz-2)*nprocz-1);	
	//............................................
	// Get maximum and minimum for this domain
	// Halo is included !
	min_x = double(iproc*(Nx-2)-1)*hx;
	min_y = double(jproc*(Ny-2)-1)*hy;
	min_z = double(kproc*(Nz-2)-1)*hz;
	//............................................

	//............................................
		// Pre-initialize Distance 
	for (n=0;n<N;n++){
		Distance[n]=100.0;
	}
	//............................................

	//............................................
	// .........Loop over the spheres.............
	for (p=0;p<nspheres;p++){
		// Get the sphere from the list, map to local min
		cx = List_cx[p] - min_x;
		cy = List_cy[p] - min_y;
		cz = List_cz[p] - min_z;
		r = List_rad[p];
		// Check if
		// Range for this sphere in global indexing
		imin = int ((cx-r)/hx);
		imax = int ((cx+r)/hx)+2;
		jmin = int ((cy-r)/hy);
		jmax = int ((cy+r)/hy)+2;
		kmin = int ((cz-r)/hz);
		kmax = int ((cz+r)/hz)+2;
		// Obviously we have to do something at the edges
		if (imin<0)		imin = 0;
		if (imin>Nx)	imin = Nx;
		if (imax<0)		imax = 0;
		if (imax>Nx)	imax = Nx;
		if (jmin<0)		jmin = 0;
		if (jmin>Ny)	jmin = Ny;
		if (jmax<0)		jmax = 0;
		if (jmax>Ny)	jmax = Ny;
		if (kmin<0)		kmin = 0;
		if (kmin>Nz)	kmin = Nz;
		if (kmax<0)		kmax = 0;
		if (kmax>Nz)	kmax = Nz;
		// Loop over the domain for this sphere (may be null)
		for (i=imin;i<imax;i++){
			for (j=jmin;j<jmax;j++){
				for (k=kmin;k<kmax;k++){
					// x,y,z is distance in physical units
					x = i*hx;
					y = j*hy;
					z = k*hz;
					// if inside sphere, set to zero
					// get the position in the list
					n = k*Nx*Ny+j*Nx+i;			
					// Compute the distance
					distance = sqrt((cx-x)*(cx-x)+(cy-y)*(cy-y)+(cz-z)*(cz-z)) - r;
					// Assign the minimum distance
					if (distance < Distance[n])		Distance[n] = distance;
				
				}
			}
		}
	}
	
	// Map the distance to lattice units
	for (n=0; n<N; n++)	Distance[n] = Distance[n]/hx;
}

inline void GenerateResidual(char *ID, int Nx, int Ny, int Nz, double Saturation)
{
	//.......................................................................
	int i,j,k,n,Number,N;
	int x,y,z,ii,jj,kk;
	int sizeX,sizeY,sizeZ;
	int *SizeX, *SizeY, *SizeZ;

	srand(10009);
//	float bin;
	//.......................................................................
	N = Nx*Ny*Nz;
	
	int binCount=1;
	SizeX = new int [binCount];
	SizeY = new int [binCount];
	SizeZ = new int [binCount];
	
	SizeX[0] = 32;
	SizeY[0] = 32;
	SizeZ[0] = 32;
	
/*	ifstream Dist("distribution.in");
	for (int bin=0; bin<binCount; bin++){
		Dist >> SizeX[bin];
		Dist >> SizeY[bin];
		Dist >> SizeZ[bin];
	}
	Dist.close();
*/	//.......................................................................
//	cout << "Generating blocks... " << endl;	
	// Count for the total number of oil nodes 
	int count = 0;
	// Count the total number of non-solid nodes
	int total = 0;
	for (i=0;i<N;i++){
		if (ID[i] != 0) total++;
	}
	
	sizeX = sizeY = sizeZ = 32;
	
	//	for (int feature = 0; feature < Number; feature++){
	float sat = 0.f;
	Number = 0;		// number of features
	while (sat < Saturation){
		Number++;
		// Randomly generate a point in the domain
		x = Nx*float(rand())/float(RAND_MAX);
		y = Ny*float(rand())/float(RAND_MAX);
		z = Nz*float(rand())/float(RAND_MAX);
		
//		bin = binCount*float(rand())/float(RAND_MAX);
//		sizeX = SizeX[int(floor(bin))];
//		sizeY = SizeY[int(floor(bin))];
//		sizeZ = SizeZ[int(floor(bin))];
		
//		cout << "Sampling from bin no. " << floor(bin) << endl; 
//		cout << "Feature size is: " << sizeX << "x" << sizeY << "x" << sizeZ << endl; 
		
		for (k=z;k<z+sizeZ;k++){
			for (j=y;j<y+sizeY;j++){
				for (i=x;i<x+sizeX;i++){
					// Identify nodes in the domain (periodic BC)
					ii = i;
					jj = j;
					kk = k;					
					if (ii < 1)			ii+=(Nx-2);
					if (jj < 1)			jj+=(Ny-2);
					if (kk < 1)			kk+=(Nz-2);
					if (!(ii < Nx-1))		ii-=(Nx-2);
					if (!(jj < Ny-1))		jj-=(Ny-2);
					if (!(kk < Nz-1))		kk-=(Nz-2);
					
					n = kk*Nx*Ny+jj*Nx+ii;
					
					if (ID[n] == 1){
						ID[n] = 2;
						count++;
					}
				}
			}
		}
		sat = float(count)/total;
	}
	//.......................................................................
}

inline void WriteLocalSolidID(char *FILENAME, char *ID, int N)
{
	char value;
	ofstream File(FILENAME,ios::binary);
	for (int n=0; n<N; n++){
		value = ID[n];
		File.write((char*) &value, sizeof(value));
	}
	File.close();
}

inline void WriteLocalSolidDistance(char *FILENAME, double *Distance, int N)
{
	double value;
	ofstream File(FILENAME,ios::binary);
	for (int n=0; n<N; n++){
		value = Distance[n];
		File.write((char*) &value, sizeof(value));
	}
	File.close();
}
