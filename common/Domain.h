// Created by James McClure
// Copyright 2008-2013
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <exception>      // std::exception
#include <stdexcept>
#include <stdio.h>
#include "common/Utilities.h"


using namespace std;

struct Domain{
	
	Domain(int nx, int ny, int nz, int rnk, int npx, int npy, int npz, 
			double lx, double ly, double lz){
		Nx = nx+2; Ny = ny+2; Nz = nz+2; 
		Lx = lx, Ly = ly, Lz = lz;
		rank = rnk;
		nprocx=npx; nprocy=npy; nprocz=npz;
		N = Nx*Ny*Nz;
		ID = new char [N];
		Blobs.New(Nx,Ny,Nz);
	}

	// Basic domain information
	int Nx,Ny,Nz,N;
	int iproc,jproc,kproc;
	int nprocx,nprocy,nprocz;
	double Lx,Ly,Lz;
	int rank;
	//**********************************
	// MPI ranks for all 18 neighbors
	//**********************************
	int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
	int rank_xy,rank_XY,rank_xY,rank_Xy;
	int rank_xz,rank_XZ,rank_xZ,rank_Xz;
	int rank_yz,rank_YZ,rank_yZ,rank_Yz;
	//**********************************
	
	// Solid indicator function
	char *ID;
	// Blob information
	IntArray Blobs;
	
	void InitializeRanks()
	{
		// map the rank to the block index
		iproc = rank%nprocx;
		jproc = (rank/nprocx)%nprocy;
		kproc = rank/(nprocx*nprocy);

		// set up the neighbor ranks
	    int i = iproc;
	    int j = jproc;
	    int k = kproc;
		rank_X = getRankForBlock(i+1,j,k);
		rank_x = getRankForBlock(i-1,j,k);
		rank_Y = getRankForBlock(i,j+1,k);
		rank_y = getRankForBlock(i,j-1,k);
		rank_Z = getRankForBlock(i,j,k+1);
		rank_z = getRankForBlock(i,j,k-1);
		rank_XY = getRankForBlock(i+1,j+1,k);
		rank_xy = getRankForBlock(i-1,j-1,k);
		rank_Xy = getRankForBlock(i+1,j-1,k);
		rank_xY = getRankForBlock(i-1,j+1,k);
		rank_XZ = getRankForBlock(i+1,j,k+1);
		rank_xz = getRankForBlock(i-1,j,k-1);
		rank_Xz = getRankForBlock(i+1,j,k-1);
		rank_xZ = getRankForBlock(i-1,j,k+1);
		rank_YZ = getRankForBlock(i,j+1,k+1);
		rank_yz = getRankForBlock(i,j-1,k-1);
		rank_Yz = getRankForBlock(i,j+1,k-1);
		rank_yZ = getRankForBlock(i,j-1,k+1);
	}
	
private:
	int getRankForBlock( int i, int j, int k )
	{
		int i2 = (i+nprocx)%nprocx;
		int j2 = (j+nprocy)%nprocy;
		int k2 = (k+nprocz)%nprocz;
		return i2 + j2*nprocx + k2*nprocx*nprocy;
	}
	
};

inline void ReadSpherePacking(int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad)
{
	// Read in the full sphere pack
	//...... READ IN THE SPHERES...................................
	cout << "Reading the packing file..." << endl;
	FILE *fid = fopen("pack.out","rb");
	INSIST(fid!=NULL,"Error opening pack.out");
	//.........Trash the header lines..........
	char * line = new char[100];
	fgets(line, 100, fid);
	fgets(line, 100, fid);
	fgets(line, 100, fid);
	fgets(line, 100, fid);
	fgets(line, 100, fid);
	//........read the spheres..................
    // We will read until a blank like or end-of-file is reached
	int count = 0;
	while ( !feof(fid) && fgets(line,100,fid)>0 ) {
		char* line2 = line;
		List_cx[count] = strtod(line2,&line2);
		List_cy[count] = strtod(line2,&line2);
		List_cz[count] = strtod(line2,&line2);
		List_rad[count] = strtod(line2,&line2);
		count++;
	}
	cout << "Number of spheres extracted is: " << count << endl;
    INSIST( count==nspheres, "Specified number of spheres is probably incorrect!" );
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

#ifdef NORANDOM
	srand(10009);
#else
	srand(time(NULL));
#endif
//	float bin;
	//.......................................................................
	N = Nx*Ny*Nz;
	
	int bin, binCount;	
	ifstream Dist("BlobSize.in");
	Dist >> binCount;
//	printf("Number of blob sizes: %i \n",binCount);
	SizeX = new int [binCount];
	SizeY = new int [binCount];
	SizeZ = new int [binCount];
	for (bin=0; bin<binCount; bin++){
		Dist >> SizeX[bin];
		Dist >> SizeY[bin];
		Dist >> SizeZ[bin];
	//	printf("Blob %i dimension: %i x %i x %i \n",bin, SizeX[bin], SizeY[bin], SizeZ[bin]);
	}
	Dist.close();
	//.......................................................................
//	cout << "Generating blocks... " << endl;	
	// Count for the total number of oil nodes 
	int count = 0;
	// Count the total number of non-solid nodes
	int total = 0;
	for (i=0;i<N;i++){
		if (ID[i] != 0) total++;
	}
	
	float sat = 0.f;
	Number = 0;		// number of features
	while (sat < Saturation){
		Number++;
		// Randomly generate a point in the domain
		x = Nx*float(rand())/float(RAND_MAX);
		y = Ny*float(rand())/float(RAND_MAX);
		z = Nz*float(rand())/float(RAND_MAX);
		
		bin = int(floor(binCount*float(rand())/float(RAND_MAX)));
		sizeX = SizeX[bin];
		sizeY = SizeY[bin];
		sizeZ = SizeZ[bin];
		
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
					
					if (ID[n] == 2){
						ID[n] = 1;
						count++;
					}
				}
			}
		}
		sat = float(count)/total;
	}
	//.......................................................................
}

inline void FlipID(char *ID, int N)
{
	for (int n=0; n<N; n++){
		if  	 (ID[n] == 1)	ID[n] = 2;
		else if  (ID[n] == 2)	ID[n] = 1;
	}
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


inline void WriteCheckpoint(char *FILENAME, double *cDen, double *cDistEven, double *cDistOdd, int N)
{
	int q,n;
	double value;
	ofstream File(FILENAME,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		value = cDen[2*n];
		File.write((char*) &value, sizeof(value));
		value = cDen[2*n+1];
		File.write((char*) &value, sizeof(value));
		// Write the even distributions
		for (q=0; q<10; q++){
			value = cDistEven[q*N+n];
			File.write((char*) &value, sizeof(value));
		}
		// Write the odd distributions
		for (q=0; q<9; q++){
			value = cDistOdd[q*N+n];
			File.write((char*) &value, sizeof(value));
		}
	}
	File.close();

}

inline void ReadCheckpoint(char *FILENAME, double *cDen, double *cDistEven, double *cDistOdd, int N)
{
	int q,n;
	double value;
	ifstream File(FILENAME,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		File.read((char*) &value, sizeof(value));
		cDen[2*n] = value;
	//	if (n== 66276)	printf("Density a  = %f \n",value);
		File.read((char*) &value, sizeof(value));
		cDen[2*n+1] = value;
	//	if (n== 66276)	printf("Density b  = %f \n",value);
		// Read the even distributions
		for (q=0; q<10; q++){
			File.read((char*) &value, sizeof(value));
			cDistEven[q*N+n] = value;
	//		if (n== 66276)	printf("dist even %i  = %f \n",q,value);
		}
		// Read the odd distributions
		for (q=0; q<9; q++){
			File.read((char*) &value, sizeof(value));
			cDistOdd[q*N+n] = value;
	//		if (n== 66276)	printf("dist even %i  = %f \n",q,value);
		}
	}
	File.close();
}

inline void ReadBinaryFile(char *FILENAME, double *Data, int N)
{
	int n;
	double value;
	ifstream File(FILENAME,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		File.read((char*) &value, sizeof(value));
		Data[n] = value;

	}
	File.close();
}

