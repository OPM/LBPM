#ifndef Domain_INC
#define Domain_INC
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

#include "common/Array.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"

using namespace std;


//! Read the domain information file
void read_domain( int rank, int nprocs, MPI_Comm comm, 
    int& nprocx, int& nprocy, int& nprocz, int& nx, int& ny, int& nz,
    int& nspheres, double& Lx, double& Ly, double& Lz );

//! Class to hold domain info
struct Domain{
    // Default constructor
	Domain(int nx, int ny, int nz, int rnk, int npx, int npy, int npz, 
			double lx, double ly, double lz, int BC);

    // Destructor
	~Domain();
    
	// Basic domain information
	int Nx,Ny,Nz,N;
	int iproc,jproc,kproc;
 	int nprocx,nprocy,nprocz;
    double Lx,Ly,Lz,Volume;
	int rank;
	int BoundaryCondition;
	RankInfoStruct rank_info;
	MPI_Group Group;	// Group of processors associated with this domain
	MPI_Comm Comm;		// MPI Communicator for this domain

	//**********************************
	// MPI ranks for all 18 neighbors
	//**********************************
	int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
	int rank_xy,rank_XY,rank_xY,rank_Xy;
	int rank_xz,rank_XZ,rank_xZ,rank_Xz;
	int rank_yz,rank_YZ,rank_yZ,rank_Yz;
	//**********************************
	//......................................................................................
	// Get the actual D3Q19 communication counts (based on location of solid phase)
	// Discrete velocity set symmetry implies the sendcount = recvcount
	//......................................................................................
	int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z;
	int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz, sendCount_xZ;
	int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ, sendCount_XZ;
	//......................................................................................
	int *sendList_x, *sendList_y, *sendList_z, *sendList_X, *sendList_Y, *sendList_Z;
	int *sendList_xy, *sendList_yz, *sendList_xz, *sendList_Xy, *sendList_Yz, *sendList_xZ;
	int *sendList_xY, *sendList_yZ, *sendList_Xz, *sendList_XY, *sendList_YZ, *sendList_XZ;
	//......................................................................................
	int *sendBuf_x, *sendBuf_y, *sendBuf_z, *sendBuf_X, *sendBuf_Y, *sendBuf_Z;
	int *sendBuf_xy, *sendBuf_yz, *sendBuf_xz, *sendBuf_Xy, *sendBuf_Yz, *sendBuf_xZ;
	int *sendBuf_xY, *sendBuf_yZ, *sendBuf_Xz, *sendBuf_XY, *sendBuf_YZ, *sendBuf_XZ;
	//......................................................................................
	int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z;
	int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz, recvCount_xZ;
	int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ, recvCount_XZ;
	//......................................................................................
	int *recvList_x, *recvList_y, *recvList_z, *recvList_X, *recvList_Y, *recvList_Z;
	int *recvList_xy, *recvList_yz, *recvList_xz, *recvList_Xy, *recvList_Yz, *recvList_xZ;
	int *recvList_xY, *recvList_yZ, *recvList_Xz, *recvList_XY, *recvList_YZ, *recvList_XZ;
	//......................................................................................
	int *recvBuf_x, *recvBuf_y, *recvBuf_z, *recvBuf_X, *recvBuf_Y, *recvBuf_Z;
	int *recvBuf_xy, *recvBuf_yz, *recvBuf_xz, *recvBuf_Xy, *recvBuf_Yz, *recvBuf_xZ;
	int *recvBuf_xY, *recvBuf_yZ, *recvBuf_Xz, *recvBuf_XY, *recvBuf_YZ, *recvBuf_XZ;
	//......................................................................................	
	double *sendData_x, *sendData_y, *sendData_z, *sendData_X, *sendData_Y, *sendData_Z;
	double *sendData_xy, *sendData_yz, *sendData_xz, *sendData_Xy, *sendData_Yz, *sendData_xZ;
	double *sendData_xY, *sendData_yZ, *sendData_Xz, *sendData_XY, *sendData_YZ, *sendData_XZ;
	double *recvData_x, *recvData_y, *recvData_z, *recvData_X, *recvData_Y, *recvData_Z;
	double *recvData_xy, *recvData_yz, *recvData_xz, *recvData_Xy, *recvData_Yz, *recvData_xZ;
	double *recvData_xY, *recvData_yZ, *recvData_Xz, *recvData_XY, *recvData_YZ, *recvData_XZ;

	// Solid indicator function
	char *id;
	// Blob information
	IntArray BlobLabel;
	IntArray BlobGraph;

	void InitializeRanks();
	void CommInit(MPI_Comm comm);
	void CommunicateMeshHalo(DoubleArray &Mesh);
	void BlobComm(MPI_Comm comm);

	void AssignBlobConnections();

private:

	inline int getRankForBlock( int i, int j, int k )
	{
		int i2 = (i+nprocx)%nprocx;
		int j2 = (j+nprocy)%nprocy;
		int k2 = (k+nprocz)%nprocz;
		return i2 + j2*nprocx + k2*nprocx*nprocy;
	}
	inline void PackBlobData(int *list, int count, int *sendbuf, int *data){
		// Fill in the phase ID values from neighboring processors
		// This packs up the values that need to be sent from one processor to another
		int idx,n;
		for (idx=0; idx<count; idx++){
			n = list[idx];
			sendbuf[idx] = data[n];
		}
	}
	inline void UnpackBlobData(int *list, int count, int *recvbuf, int *data){
		// Fill in the phase ID values from neighboring processors
		// This unpacks the values once they have been recieved from neighbors
		int idx,n;
		for (idx=0; idx<count; idx++){
			n = list[idx];
			data[n] = recvbuf[idx];
		}
	}
	
	int VoxelConnection(int n);
	
	void getBlobConnections(int * List, int count, int neighbor, int direction);
};



// Inline function to read line without a return argument
static inline void fgetl( char * str, int num, FILE * stream )
{
    char* ptr = fgets( str, num, stream );
    if ( 0 ) {char *temp = (char *)&ptr; temp++;}
}



inline double SSO(DoubleArray &Distance, char *ID, Domain &Dm, int timesteps){
    /*
     * This routine converts the data in the Distance array to a signed distance
     * by solving the equation df/dt = sign(1-|grad f|), where Distance provides
     * the values of f on the mesh associated with domain Dm
     * It has been tested with segmented data initialized to values [-1,1]
     * and will converge toward the signed distance to the surface bounding the associated phases
     */

    int Q=26;
    int q,i,j,k;
    double dt=0.1;
    int in,jn,kn;
    double Dqx,Dqy,Dqz,Dx,Dy,Dz,W;
    double nx,ny,nz,Cqx,Cqy,Cqz,sign,norm;
    double TotalVariation=0.0;

    const static int D3Q27[26][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
        {1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},
        {0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1},{1,1,1},{-1,-1,-1},{1,1,-1},{-1,-1,1},
        {-1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1}};

    double weights[26];
    // Compute the weights from the finite differences
    for (q=0; q<Q; q++){
        weights[q] = sqrt(1.0*(D3Q27[q][0]*D3Q27[q][0]) + 1.0*(D3Q27[q][1]*D3Q27[q][1]) + 1.0*(D3Q27[q][2]*D3Q27[q][2]));
    }

    int xdim,ydim,zdim;
    xdim=Dm.Nx-2;
    ydim=Dm.Ny-2;
    zdim=Dm.Nz-2;
    fillHalo<double> fillData(Dm.Comm, Dm.rank_info,xdim,ydim,zdim,1,1,1,0,1);

    int count = 0;
    while (count < timesteps){

        // Communicate the halo of values
        fillData.fill(Distance);

        TotalVariation=0.0;
        // Execute the next timestep
        for (k=1;k<Dm.Nz-1;k++){
            for (j=1;j<Dm.Ny-1;j++){
                for (i=1;i<Dm.Nx-1;i++){

                	int n = k*Dm.Nx*Dm.Ny + j*Dm.Nx + i;
			//sign = Distance(i,j,k) / fabs(Distance(i,j,k));
		        sign = -1;
		    	if (ID[n] == 1) sign = 1;
              /*
                    if (!(i+1<Nx))     nx=0.5*Distance(i,j,k);
                    else             nx=0.5*Distance(i+1,j,k);;
                    if (!(j+1<Ny))     ny=0.5*Distance(i,j,k);
                    else             ny=0.5*Distance(i,j+1,k);
                    if (!(k+1<Nz))     nz=0.5*Distance(i,j,k);
                    else             nz=0.5*Distance(i,j,k+1);
                    if (i<1)          nx-=0.5*Distance(i,j,k);
                    else             nx-=0.5*Distance(i-1,j,k);
                    if (j<1)         ny-=0.5*Distance(i,j,k);
                    else             ny-=0.5*Distance(i,j-1,k);
                    if (k<1)          nz-=0.5*Distance(i,j,k);
                    else             nz-=0.5*Distance(i,j,k-1);
                    */

                    //............Compute the Gradient...................................
                    nx = 0.5*(Distance(i+1,j,k) - Distance(i-1,j,k));
                    ny = 0.5*(Distance(i,j+1,k) - Distance(i,j-1,k));
                    nz = 0.5*(Distance(i,j,k+1) - Distance(i,j,k-1));

                    W = 0.0;    Dx = Dy = Dz = 0.0;
                    // also ignore places where the gradient is zero since this will not
                    // result in any local change to Distance
                    if (nx*nx+ny*ny+nz*nz > 0.0 ){
                        for (q=0; q<26; q++){
                            Cqx = 1.0*D3Q27[q][0];
                            Cqy = 1.0*D3Q27[q][1];
                            Cqz = 1.0*D3Q27[q][2];
                            // get the associated neighbor
                            in = i + D3Q27[q][0];
                            jn = j + D3Q27[q][1];
                            kn = k + D3Q27[q][2];

                            // make sure the neighbor is in the domain (periodic BC)
                            /*                if (in < 0 ) in +=Nx;
                             *     don't need this in parallel since MPI handles the halos
                             if (jn < 0 ) jn +=Ny;
                             if (kn < 0 ) kn +=Nz;
                             if (!(in < Nx) ) in -=Nx;
                             if (!(jn < Ny) ) jn -=Ny;
                             if (!(kn < Nz) ) kn -=Nz;
                                         // symmetric boundary
                            if (in < 0 ) in = i;
                            if (jn < 0 ) jn = j;
                            if (kn < 0 ) kn = k;
                            if (!(in < Nx) ) in = i;
                            if (!(jn < Ny) ) jn = k;
                            if (!(kn < Nz) ) kn = k;
                            */

                            // Compute the gradient using upwind finite differences
                            Dqx = weights[q]*(Distance(i,j,k) - Distance(in,jn,kn))*Cqx;
                            Dqy = weights[q]*(Distance(i,j,k) - Distance(in,jn,kn))*Cqy;
                            Dqz = weights[q]*(Distance(i,j,k) - Distance(in,jn,kn))*Cqz;

                            // Only include upwind derivatives
                            if (sign*(nx*Cqx + ny*Cqy + nz*Cqz) < 0.0 ){

                            	Dx += Dqx;
                            	Dy += Dqy;
                            	Dz += Dqz;
                            	W += weights[q];
                            }
                        }
                        // Normalize by the weight to get the approximation to the gradient
                        if (fabs(W) > 0.0){
                        	Dx /= W;
                        	Dy /= W;
                        	Dz /= W;
                        }
                        norm = sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
                    }
                    else{
                        norm = 0.0;
                    }
                    Distance(i,j,k) += dt*sign*(1.0 - norm);
                    TotalVariation +=  dt*sign*(1.0 - norm);
                    // Disallow any change in phase
                   // if (Distance(i,j,k)*2.0*(ID[n]-1.0) < 0) Distance(i,j,k) = -Distance(i,j,k);
                }
            }
        }
        TotalVariation /= (Dm.Nx-2)*(Dm.Ny-2)*(Dm.Nz-2);
        count++;
    }

    return TotalVariation;
}


inline void ReadSpherePacking(int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad)
{
    // Read in the full sphere pack
    //...... READ IN THE SPHERES...................................
    cout << "Reading the packing file..." << endl;
    FILE *fid = fopen("pack.out","rb");
    INSIST(fid!=NULL,"Error opening pack.out");
    //.........Trash the header lines..........
    char line[100];
    fgetl(line, 100, fid);
    fgetl(line, 100, fid);
    fgetl(line, 100, fid);
    fgetl(line, 100, fid);
    fgetl(line, 100, fid);
    //........read the spheres..................
    // We will read until a blank like or end-of-file is reached
    int count = 0;
    while ( !feof(fid) && fgets(line,100,fid)!=NULL ) {
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
    int N = Nx*Ny*Nz;     // Domain size, including the halo
    double hx,hy,hz;
    double x,y,z;
    double cx,cy,cz,r;
    int imin,imax,jmin,jmax,kmin,kmax;
    int p,i,j,k,n;
    //............................................
    double min_x,min_y,min_z;
//    double max_x,max_y,max_z;
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
//    max_x = ((iproc+1)*Nx+1)*hx;
//    max_y = ((jproc+1)*Ny+1)*hy;
//    max_z = ((kproc+1)*Nz+1)*hz;
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
        if (imin<0)        imin = 0;
        if (imin>Nx)    imin = Nx;
        if (imax<0)        imax = 0;
        if (imax>Nx)    imax = Nx;
        if (jmin<0)        jmin = 0;
        if (jmin>Ny)    jmin = Ny;
        if (jmax<0)        jmax = 0;
        if (jmax>Ny)    jmax = Ny;
        if (kmin<0)        kmin = 0;
        if (kmin>Nz)    kmin = Nz;
        if (kmax<0)        kmax = 0;
        if (kmax>Nz)    kmax = Nz;
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
    int N = Nx*Ny*Nz;     // Domain size, including the halo
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
        imin = int ((cx-2*r)/hx);
        imax = int ((cx+2*r)/hx)+2;
        jmin = int ((cy-2*r)/hy);
        jmax = int ((cy+2*r)/hy)+2;
        kmin = int ((cz-2*r)/hz);
        kmax = int ((cz+2*r)/hz)+2;
        // Obviously we have to do something at the edges
        if (imin<0)        imin = 0;
        if (imin>Nx)    imin = Nx;
        if (imax<0)        imax = 0;
        if (imax>Nx)    imax = Nx;
        if (jmin<0)        jmin = 0;
        if (jmin>Ny)    jmin = Ny;
        if (jmax<0)        jmax = 0;
        if (jmax>Ny)    jmax = Ny;
        if (kmin<0)        kmin = 0;
        if (kmin>Nz)    kmin = Nz;
        if (kmax<0)        kmax = 0;
        if (kmax>Nz)    kmax = Nz;
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
                    if (distance < Distance[n])        Distance[n] = distance;
                
                }
            }
        }
    }
    
    // Map the distance to lattice units
    for (n=0; n<N; n++)    Distance[n] = Distance[n]/hx;
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


inline void WriteCheckpoint(const char *FILENAME, const double *cDen, const double *cDistEven, const double *cDistOdd, int N)
{
    int q,n;
    double value;
    ofstream File(FILENAME,ios::binary);
    for (n=0; n<N; n++){
        // Write the two density values
        value = cDen[n];
        File.write((char*) &value, sizeof(value));
        value = cDen[N+n];
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
    int q=0, n=0;
    double value=0;
    ifstream File(FILENAME,ios::binary);
    for (n=0; n<N; n++){
        // Write the two density values
        File.read((char*) &value, sizeof(value));
        cDen[n] = value;
    //    if (n== 66276)    printf("Density a  = %f \n",value);
        File.read((char*) &value, sizeof(value));
        cDen[N+n] = value;
    //    if (n== 66276)    printf("Density b  = %f \n",value);
        // Read the even distributions
        for (q=0; q<10; q++){
            File.read((char*) &value, sizeof(value));
            cDistEven[q*N+n] = value;
    //        if (n== 66276)    printf("dist even %i  = %f \n",q,value);
        }
        // Read the odd distributions
        for (q=0; q<9; q++){
            File.read((char*) &value, sizeof(value));
            cDistOdd[q*N+n] = value;
    //        if (n== 66276)    printf("dist even %i  = %f \n",q,value);
        }
    }
    File.close();
}

inline void ReadBinaryFile(char *FILENAME, double *Data, int N)
{
    int n;
    double value;
    ifstream File(FILENAME,ios::binary);
    if (File.good()){
        for (n=0; n<N; n++){
            // Write the two density values
            File.read((char*) &value, sizeof(value));
            Data[n] = value;

        }
    }
    else {
        for (n=0; n<N; n++) Data[n] = 1.2e-34;
    }
    File.close();

}
inline double Eikonal(DoubleArray &Distance, char *ID, Domain &Dm, int timesteps);


#endif
