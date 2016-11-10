/*
 * Pre-processor to generate signed distance function from segmented data
 * segmented data should be stored in a raw binary file as 1-byte integer (type char)
 * will output distance functions for phases
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "common/Array.h"
#include "common/Domain.h"


//*************************************************************************
// Morpohologica pre-processor
//   Initialize phase distribution using morphological approach
//   Signed distance function is used to determine fluid configuration
//*************************************************************************
inline void PackID(int *list, int count, char *sendbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This packs up the values that need to be sent from one processor to another
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = ID[n];
	}
}
//***************************************************************************************

inline void UnpackID(int *list, int count, char *recvbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This unpacks the values once they have been recieved from neighbors
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		ID[n] = recvbuf[idx];
	}
}

//***************************************************************************************



int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);

 	//double Rcrit_new=1.f; // Hard-coded 'Rcrit' to avoid any calculations under resolutions.
 	double Rcrit_new=0.f; 
	double SW=strtod(argv[1],NULL);
	if (rank==0){
		//printf("Critical radius %f (voxels)\n",Rcrit);
		printf("Target saturation %f \n",SW);
	}
	//	}

	//.......................................................................
	// Reading the domain information file
	//.......................................................................
	int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
	double Lx, Ly, Lz;
	int i,j,k,n;
	int BC=0;

	if (rank==0){
		ifstream domain("Domain.in");
		domain >> nprocx;
		domain >> nprocy;
		domain >> nprocz;
		domain >> nx;
		domain >> ny;
		domain >> nz;
		domain >> nspheres;
		domain >> Lx;
		domain >> Ly;
		domain >> Lz;

	}
	MPI_Barrier(comm);
	// Computational domain
	MPI_Bcast(&nx,1,MPI_INT,0,comm);
	MPI_Bcast(&ny,1,MPI_INT,0,comm);
	MPI_Bcast(&nz,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
	MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
	//.................................................
	MPI_Barrier(comm);

	// Check that the number of processors >= the number of ranks
	if ( rank==0 ) {
		printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
		printf("Number of MPI ranks used: %i \n", nprocs);
		printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
	}
	if ( nprocs < nprocx*nprocy*nprocz ){
		ERROR("Insufficient number of processors");
	}

	char LocalRankFilename[40];

	int BoundaryCondition=1;
	Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);

	nx+=2; ny+=2; nz+=2;
	int N = nx*ny*nz;
	char *id;
	id = new char[N];

	// Define communication sub-domain -- everywhere
	for (int k=0; k<nz; k++){
		for (int j=0; j<ny; j++){
			for (int i=0; i<nx; i++){
				n = k*nx*ny+j*nx+i;
				Dm.id[n] = 1;
			}
		}
	}
	Dm.CommInit(comm);

	DoubleArray SignDist(nx,ny,nz);
	// Read the signed distance from file
	sprintf(LocalRankFilename,"SignDist.%05i",rank);
	FILE *DIST = fopen(LocalRankFilename,"rb");
	size_t ReadSignDist;
	ReadSignDist=fread(SignDist.data(),8,N,DIST);
	if (ReadSignDist != size_t(N)) printf("lbpm_morphdrain_pp: Error reading signed distance function (rank=%i)\n",rank);
	fclose(DIST);

	int count,countGlobal,totalGlobal;
	count = 0;
	double maxdist=0;
	double maxdistGlobal;
	for (int k=1; k<nz-1; k++){
		for (int j=1; j<ny-1; j++){
			for (int i=1; i<nx-1; i++){
				n = k*nx*ny+j*nx+i;
				if (SignDist(i,j,k) < 0.0)  id[n] = 0;
				else{
					// initially saturated with wetting phase
					id[n] = 2;
					count++;
					// extract maximum distance for critical radius
					if ( SignDist(i,j,k) > maxdist) maxdist=SignDist(i,j,k);
				}
			}
		}
	}
	// total Global is the number of nodes in the pore-space
	MPI_Allreduce(&count,&totalGlobal,1,MPI_INT,MPI_SUM,comm);
	MPI_Allreduce(&maxdist,&maxdistGlobal,1,MPI_DOUBLE,MPI_MAX,comm);
	float porosity=float(totalGlobal)/(nprocx*nprocy*nprocz*(nx-2)*(ny-2)*(nz-2));
	if (rank==0) printf("Media Porosity: %f \n",porosity);
	if (rank==0) printf("Maximum pore size: %f \n",maxdistGlobal);\

	// Generate a histogram of pore size distribution
	// Get all local pore sizes (local maxima)
	if (rank==0) printf("Generating a histogram of pore sizes \n");
	int NumBins=100;
	int *BinCounts;
	BinCounts = new int [NumBins];
	int *GlobalHistogram;
	GlobalHistogram = new int [NumBins];
	double GlobalValue;
	double BinWidth,MinPoreSize,MaxPoreSize;
	std::vector<double> PoreSize;
	for (int k=1; k<nz-1; k++){
		for (int j=1; j<ny-1; j++){
			for (int i=1; i<nx-1; i++){
				n = k*nx*ny+j*nx+i;
				if (SignDist(i,j,k) > 0.0){
					// Generate a list of all local maxima (each processor -- aggregate these later)
					if ( SignDist(i,j,k) > SignDist(i+1,j,k) && SignDist(i,j,k) > SignDist(i-1,j,k) &&
							SignDist(i,j,k) > SignDist(i,j+1,k) && SignDist(i,j,k) > SignDist(i,j-1,k) &&
							SignDist(i,j,k) > SignDist(i,j,k+1) && SignDist(i,j,k) > SignDist(i,j,k-1) &&
							SignDist(i,j,k) > SignDist(i+1,j+1,k) && SignDist(i,j,k) > SignDist(i-1,j+1,k) &&
							SignDist(i,j,k) > SignDist(i+1,j-1,k) && SignDist(i,j,k) > SignDist(i-1,j-1,k) &&
							SignDist(i,j,k) > SignDist(i+1,j,k+1) && SignDist(i,j,k) > SignDist(i-1,j,k+1) &&
							SignDist(i,j,k) > SignDist(i+1,j,k-1) && SignDist(i,j,k) > SignDist(i-1,j,k-1) &&
							SignDist(i,j,k) > SignDist(i,j+1,k+1) && SignDist(i,j,k) > SignDist(i,j-1,k+1) &&
							SignDist(i,j,k) > SignDist(i,j+1,k-1) && SignDist(i,j,k) > SignDist(i,j-1,k-1) &&
							SignDist(i,j,k) > SignDist(i+1,j+1,k+1) && SignDist(i,j,k) > SignDist(i-1,j-1,k-1) &&
							SignDist(i,j,k) > SignDist(i+1,j-1,k+1) && SignDist(i,j,k) > SignDist(i-1,j+1,k-1) &&
							SignDist(i,j,k) > SignDist(i-1,j+1,k+1) && SignDist(i,j,k) > SignDist(i+1,j-1,k-1) &&
							SignDist(i,j,k) > SignDist(i+1,j+1,k-1) && SignDist(i,j,k) > SignDist(i-1,j-1,k+1)){
						// save the size of each pore
					  PoreSize.push_back(SignDist(i,j,k));
					}
				}
			}
		}
	}
	// Compute min and max poresize
	if (rank==0) printf("    computing local minimum and maximum... \n");
	MinPoreSize=MaxPoreSize=PoreSize[0];
	for (int idx=0; idx<PoreSize.size(); idx++){
		if (PoreSize[idx] < MinPoreSize) MinPoreSize=PoreSize[idx];
		if (PoreSize[idx] > MaxPoreSize) MaxPoreSize=PoreSize[idx];
	}
	// reduce to get global minimum and maximum
	MPI_Allreduce(&MinPoreSize,&GlobalValue,1,MPI_DOUBLE,MPI_MIN,comm);
	MinPoreSize=GlobalValue;
	MPI_Allreduce(&MaxPoreSize,&GlobalValue,1,MPI_DOUBLE,MPI_MAX,comm);
	MaxPoreSize=GlobalValue;
	//if (rank==0) printf("    MaxPoreSize %f\n", MaxPoreSize);
	//if (rank==0) printf("    MinPoreSize %f\n", MinPoreSize);
	// Generate histogram counts
	if (rank==0) printf("    generating local bin counts... \n");
	BinWidth=(MaxPoreSize-MinPoreSize)/double(NumBins);
	for (int idx=0; idx<PoreSize.size(); idx++){
		double value = PoreSize[idx];
		int myBin = 0;
		while (MinPoreSize+myBin*BinWidth < value) myBin++;
		  //int((value-MinPoreSize)/double(BinWidth));
		BinCounts[myBin]+=1;
	}
	if (rank==0) printf("    summing global bin counts... \n");
	// Reduce the counts to generate the fhistogram at rank=0
	MPI_Reduce(BinCounts,GlobalHistogram,NumBins,MPI_INT,MPI_SUM,0,comm);

	FILE *PORESIZE;
	if (rank==0){
		PORESIZE=fopen("PoreSize.hist","w");
		printf("    writing PoreSize.hist \n");
		int PoreCount=0;
		int Count;
		for (int idx=0; idx<NumBins; idx++){
		  double BinCenter=MinPoreSize+idx*BinWidth;
		  Count=GlobalHistogram[idx];
		  PoreCount+=Count;
			fprintf(PORESIZE,"%i %f\n",Count,BinCenter);
		}
		fclose(PORESIZE);
		// Compute quartiles
		double Q1,Q2,Q3,Q4;
		int Qval=PoreCount/4;
		Q1=Q2=Q3=MinPoreSize;
		Q4=MaxPoreSize;
		Count=0;
		for (int idx=0; idx<NumBins; idx++){
		  double BinCenter=MinPoreSize+idx*BinWidth;
		  Count+=GlobalHistogram[idx];
		  if (Count<Qval) Q1+=BinWidth;
		  if (Count<2*Qval) Q2+=BinWidth;
		  if (Count<3*Qval) Q3+=BinWidth;
		}

		      printf("Quartiles for pore size distribution \n");
		      printf("Q1 %f\n",Q1);
		      printf("Q2 %f\n",Q2);
		      printf("Q3 %f\n",Q3);
		      printf("Q4 %f\n",Q4);
	}
	MPI_Barrier(comm);
	
	Dm.CommInit(comm);
	int iproc = Dm.iproc;
	int jproc = Dm.jproc;
	int kproc = Dm.kproc;

	// Generate the NWP configuration
	//if (rank==0) printf("Initializing morphological distribution with critical radius %f \n", Rcrit);
	if (rank==0) printf("Performing morphological opening with target saturation %f \n", SW);
	//	GenerateResidual(id,nx,ny,nz,Saturation);

	// Communication buffers
	char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
	char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
	char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
	char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
	char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
	char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
	// send buffers
	sendID_x = new char [Dm.sendCount_x];
	sendID_y = new char [Dm.sendCount_y];
	sendID_z = new char [Dm.sendCount_z];
	sendID_X = new char [Dm.sendCount_X];
	sendID_Y = new char [Dm.sendCount_Y];
	sendID_Z = new char [Dm.sendCount_Z];
	sendID_xy = new char [Dm.sendCount_xy];
	sendID_yz = new char [Dm.sendCount_yz];
	sendID_xz = new char [Dm.sendCount_xz];
	sendID_Xy = new char [Dm.sendCount_Xy];
	sendID_Yz = new char [Dm.sendCount_Yz];
	sendID_xZ = new char [Dm.sendCount_xZ];
	sendID_xY = new char [Dm.sendCount_xY];
	sendID_yZ = new char [Dm.sendCount_yZ];
	sendID_Xz = new char [Dm.sendCount_Xz];
	sendID_XY = new char [Dm.sendCount_XY];
	sendID_YZ = new char [Dm.sendCount_YZ];
	sendID_XZ = new char [Dm.sendCount_XZ];
	//......................................................................................
	// recv buffers
	recvID_x = new char [Dm.recvCount_x];
	recvID_y = new char [Dm.recvCount_y];
	recvID_z = new char [Dm.recvCount_z];
	recvID_X = new char [Dm.recvCount_X];
	recvID_Y = new char [Dm.recvCount_Y];
	recvID_Z = new char [Dm.recvCount_Z];
	recvID_xy = new char [Dm.recvCount_xy];
	recvID_yz = new char [Dm.recvCount_yz];
	recvID_xz = new char [Dm.recvCount_xz];
	recvID_Xy = new char [Dm.recvCount_Xy];
	recvID_xZ = new char [Dm.recvCount_xZ];
	recvID_xY = new char [Dm.recvCount_xY];
	recvID_yZ = new char [Dm.recvCount_yZ];
	recvID_Yz = new char [Dm.recvCount_Yz];
	recvID_Xz = new char [Dm.recvCount_Xz];
	recvID_XY = new char [Dm.recvCount_XY];
	recvID_YZ = new char [Dm.recvCount_YZ];
	recvID_XZ = new char [Dm.recvCount_XZ];
	//......................................................................................
	int sendtag,recvtag;
	sendtag = recvtag = 7;

	int x,y,z;
	int ii,jj,kk;
	int Nx = nx;
	int Ny = ny;
	int Nz = nz;
	
	double sw_old=0.0;
	double sw_new=0.0; 
    double sw_diff_old = 1.0;
    double sw_diff_new = 1.0;

	// Increase the critical radius until the target saturation is met
	double deltaR=0.01; // amount to change the radius in voxel units
    double Rcrit_old;

	int GlobalNumber = 1;
	int imin,jmin,kmin,imax,jmax,kmax;
    
    while (sw_new<SW)
    {

        Rcrit_old = Rcrit_new;
		Rcrit_new += deltaR;
	    int Window=round(Rcrit_new);
        if (Window == 0) Window = 1; // If Window = 0 at the begining, after the following process will have sw=1.0
                                     // and sw<Sw will be immediately broken
        int LocalNumber=0;
        // Initialization: saturate medium with wetting phase - need this for each iteraction before SW is met
        for (int k=1; k<nz-1; k++){
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++){
                    n = k*nx*ny+j*nx+i;
                    if (SignDist(i,j,k) < 0.0)  id[n] = 0;
                    else{
                        // initially saturated with wetting phase
                        id[n] = 2;
                    }
                }
            }
        }
        for(k=0; k<Nz; k++){
            for(j=0; j<Ny; j++){
                for(i=0; i<Nx; i++){
                    n = k*nx*ny + j*nx+i;
                    if (SignDist(i,j,k) > Rcrit_new){
                        // loop over the window and update
                        imin=max(1,i-Window);
                        jmin=max(1,j-Window);
                        kmin=max(1,k-Window);
                        imax=min(Nx-1,i+Window);
                        jmax=min(Ny-1,j+Window);
                        kmax=min(Nz-1,k+Window);
                        for (kk=kmin; kk<kmax; kk++){
                            for (jj=jmin; jj<jmax; jj++){
                                for (ii=imin; ii<imax; ii++){
                                    int nn = kk*nx*ny+jj*nx+ii;
                                    double dsq = double((ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k));
                                    if (id[nn] == 2 && dsq <= Rcrit_new*Rcrit_new){
                                        LocalNumber++;
                                        id[nn]=1;
                                    }
                                }
                            }
                        }
                    }
                    // move on
                }
            }
        }
     

        // Pack and send the updated ID values
        PackID(Dm.sendList_x, Dm.sendCount_x ,sendID_x, id);
        PackID(Dm.sendList_X, Dm.sendCount_X ,sendID_X, id);
        PackID(Dm.sendList_y, Dm.sendCount_y ,sendID_y, id);
        PackID(Dm.sendList_Y, Dm.sendCount_Y ,sendID_Y, id);
        PackID(Dm.sendList_z, Dm.sendCount_z ,sendID_z, id);
        PackID(Dm.sendList_Z, Dm.sendCount_Z ,sendID_Z, id);
        PackID(Dm.sendList_xy, Dm.sendCount_xy ,sendID_xy, id);
        PackID(Dm.sendList_Xy, Dm.sendCount_Xy ,sendID_Xy, id);
        PackID(Dm.sendList_xY, Dm.sendCount_xY ,sendID_xY, id);
        PackID(Dm.sendList_XY, Dm.sendCount_XY ,sendID_XY, id);
        PackID(Dm.sendList_xz, Dm.sendCount_xz ,sendID_xz, id);
        PackID(Dm.sendList_Xz, Dm.sendCount_Xz ,sendID_Xz, id);
        PackID(Dm.sendList_xZ, Dm.sendCount_xZ ,sendID_xZ, id);
        PackID(Dm.sendList_XZ, Dm.sendCount_XZ ,sendID_XZ, id);
        PackID(Dm.sendList_yz, Dm.sendCount_yz ,sendID_yz, id);
        PackID(Dm.sendList_Yz, Dm.sendCount_Yz ,sendID_Yz, id);
        PackID(Dm.sendList_yZ, Dm.sendCount_yZ ,sendID_yZ, id);
        PackID(Dm.sendList_YZ, Dm.sendCount_YZ ,sendID_YZ, id);
        //......................................................................................
        MPI_Sendrecv(sendID_x,Dm.sendCount_x,MPI_CHAR,Dm.rank_x,sendtag,
                recvID_X,Dm.recvCount_X,MPI_CHAR,Dm.rank_X,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_X,Dm.sendCount_X,MPI_CHAR,Dm.rank_X,sendtag,
                recvID_x,Dm.recvCount_x,MPI_CHAR,Dm.rank_x,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_y,Dm.sendCount_y,MPI_CHAR,Dm.rank_y,sendtag,
                recvID_Y,Dm.recvCount_Y,MPI_CHAR,Dm.rank_Y,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Y,Dm.sendCount_Y,MPI_CHAR,Dm.rank_Y,sendtag,
                recvID_y,Dm.recvCount_y,MPI_CHAR,Dm.rank_y,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_z,Dm.sendCount_z,MPI_CHAR,Dm.rank_z,sendtag,
                recvID_Z,Dm.recvCount_Z,MPI_CHAR,Dm.rank_Z,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Z,Dm.sendCount_Z,MPI_CHAR,Dm.rank_Z,sendtag,
                recvID_z,Dm.recvCount_z,MPI_CHAR,Dm.rank_z,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xy,Dm.sendCount_xy,MPI_CHAR,Dm.rank_xy,sendtag,
                recvID_XY,Dm.recvCount_XY,MPI_CHAR,Dm.rank_XY,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_XY,Dm.sendCount_XY,MPI_CHAR,Dm.rank_XY,sendtag,
                recvID_xy,Dm.recvCount_xy,MPI_CHAR,Dm.rank_xy,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Xy,Dm.sendCount_Xy,MPI_CHAR,Dm.rank_Xy,sendtag,
                recvID_xY,Dm.recvCount_xY,MPI_CHAR,Dm.rank_xY,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xY,Dm.sendCount_xY,MPI_CHAR,Dm.rank_xY,sendtag,
                recvID_Xy,Dm.recvCount_Xy,MPI_CHAR,Dm.rank_Xy,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xz,Dm.sendCount_xz,MPI_CHAR,Dm.rank_xz,sendtag,
                recvID_XZ,Dm.recvCount_XZ,MPI_CHAR,Dm.rank_XZ,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_XZ,Dm.sendCount_XZ,MPI_CHAR,Dm.rank_XZ,sendtag,
                recvID_xz,Dm.recvCount_xz,MPI_CHAR,Dm.rank_xz,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Xz,Dm.sendCount_Xz,MPI_CHAR,Dm.rank_Xz,sendtag,
                recvID_xZ,Dm.recvCount_xZ,MPI_CHAR,Dm.rank_xZ,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xZ,Dm.sendCount_xZ,MPI_CHAR,Dm.rank_xZ,sendtag,
                recvID_Xz,Dm.recvCount_Xz,MPI_CHAR,Dm.rank_Xz,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_yz,Dm.sendCount_yz,MPI_CHAR,Dm.rank_yz,sendtag,
                recvID_YZ,Dm.recvCount_YZ,MPI_CHAR,Dm.rank_YZ,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_YZ,Dm.sendCount_YZ,MPI_CHAR,Dm.rank_YZ,sendtag,
                recvID_yz,Dm.recvCount_yz,MPI_CHAR,Dm.rank_yz,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Yz,Dm.sendCount_Yz,MPI_CHAR,Dm.rank_Yz,sendtag,
                recvID_yZ,Dm.recvCount_yZ,MPI_CHAR,Dm.rank_yZ,recvtag,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_yZ,Dm.sendCount_yZ,MPI_CHAR,Dm.rank_yZ,sendtag,
                recvID_Yz,Dm.recvCount_Yz,MPI_CHAR,Dm.rank_Yz,recvtag,comm,MPI_STATUS_IGNORE);
        //......................................................................................
        UnpackID(Dm.recvList_x, Dm.recvCount_x ,recvID_x, id);
        UnpackID(Dm.recvList_X, Dm.recvCount_X ,recvID_X, id);
        UnpackID(Dm.recvList_y, Dm.recvCount_y ,recvID_y, id);
        UnpackID(Dm.recvList_Y, Dm.recvCount_Y ,recvID_Y, id);
        UnpackID(Dm.recvList_z, Dm.recvCount_z ,recvID_z, id);
        UnpackID(Dm.recvList_Z, Dm.recvCount_Z ,recvID_Z, id);
        UnpackID(Dm.recvList_xy, Dm.recvCount_xy ,recvID_xy, id);
        UnpackID(Dm.recvList_Xy, Dm.recvCount_Xy ,recvID_Xy, id);
        UnpackID(Dm.recvList_xY, Dm.recvCount_xY ,recvID_xY, id);
        UnpackID(Dm.recvList_XY, Dm.recvCount_XY ,recvID_XY, id);
        UnpackID(Dm.recvList_xz, Dm.recvCount_xz ,recvID_xz, id);
        UnpackID(Dm.recvList_Xz, Dm.recvCount_Xz ,recvID_Xz, id);
        UnpackID(Dm.recvList_xZ, Dm.recvCount_xZ ,recvID_xZ, id);
        UnpackID(Dm.recvList_XZ, Dm.recvCount_XZ ,recvID_XZ, id);
        UnpackID(Dm.recvList_yz, Dm.recvCount_yz ,recvID_yz, id);
        UnpackID(Dm.recvList_Yz, Dm.recvCount_Yz ,recvID_Yz, id);
        UnpackID(Dm.recvList_yZ, Dm.recvCount_yZ ,recvID_yZ, id);
        UnpackID(Dm.recvList_YZ, Dm.recvCount_YZ ,recvID_YZ, id);
        //......................................................................................

        MPI_Allreduce(&LocalNumber,&GlobalNumber,1,MPI_INT,MPI_SUM,comm);

        count = 0;
        for (int k=1; k<Nz-1; k++){
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++){
                    n=k*Nx*Ny+j*Nx+i;
                    if (id[n] == 2){
                        count++;
                    }
                }
            }
        }
        MPI_Allreduce(&count,&countGlobal,1,MPI_INT,MPI_SUM,comm);
        sw_new = float(countGlobal)/totalGlobal;
        sw_diff_new = abs(sw_new-SW);
        // for test only
        if (rank==0){
            printf("Final saturation=%f\n",sw_new);
            printf("Final critical radius=%f\n",Rcrit_new);

        }
    }

    if (sw_diff_new<sw_diff_old){
        if (rank==0){
            printf("Final saturation=%f\n",sw_new);
            printf("Final critical radius=%f\n",Rcrit_new);

        }
    }
    else{
        if (rank==0){
            printf("Final saturation=%f\n",sw_old);
            printf("Final critical radius=%f\n",Rcrit_old);

        }
    }




	sprintf(LocalRankFilename,"ID.%05i",rank);
	FILE *ID = fopen(LocalRankFilename,"wb");
	fwrite(id,1,N,ID);
	fclose(ID);

	MPI_Barrier(comm);
	MPI_Finalize();
}
