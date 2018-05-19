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

inline void FlipID(char *ID, int N)
{
	for (int n=0; n<N; n++){
		if  	 (ID[n] == 1)	ID[n] = 2;
		else if  (ID[n] == 2)	ID[n] = 1;
	}
}

//*************************************************************************
// Implementation of Two-Phase Immiscible LBM using CUDA
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

	int InitialWetting;
	double Saturation;
	//	if (argc == 3){
	//sscanf(argv[1],"%lf",&Saturation);
	//sscanf(argv[2],"%d",&InitialWetting);
	Saturation=strtod(argv[1],NULL);
	InitialWetting=atoi(argv[2]);
	if (rank==0){
		printf("Initializing wetting phase saturation of %f \n",Saturation);
		if (InitialWetting == 1)
			printf("Initial connected phase labeled (1) \n");
		else
			printf("Initial connected phase labeled (2) \n");
	}
	//	}

	if (InitialWetting == 1)	Saturation=1.0-Saturation;
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

	int BoundaryCondition=0;
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
	Dm.CommInit();

	DoubleArray SignDist(nx,ny,nz);
	// Read the signed distance from file
	sprintf(LocalRankFilename,"SignDist.%05i",rank);
	FILE *DIST = fopen(LocalRankFilename,"rb");
	size_t ReadSignDist;
	ReadSignDist=fread(SignDist.data(),8,N,DIST);
	if (ReadSignDist != size_t(N)) printf("lbpm_random_pp: Error reading signed distance function (rank=%i)\n",rank);
	fclose(DIST);

	int count,countGlobal,totalGlobal;
	count = 0;
	for (int k=1; k<nz-1; k++){
		for (int j=1; j<ny-1; j++){
			for (int i=1; i<nx-1; i++){
				n = k*nx*ny+j*nx+i;
				if (SignDist(i,j,k) < 0.0)  id[n] = 0;
				else{
					id[n] = 2;
					count++;
				}
			}
		}
	}
	// total Global is the number of nodes in the pore-space
	MPI_Allreduce(&count,&totalGlobal,1,MPI_INT,MPI_SUM,comm);
	float porosity=float(totalGlobal)/(nprocx*nprocy*nprocz*(nx-2)*(ny-2)*(nz-2));
	if (rank==0) printf("Media Porosity: %f \n",porosity);

	Dm.CommInit();
	int iproc = Dm.iproc();
	int jproc = Dm.jproc();
	int kproc = Dm.kproc();

	int bin, binCount;
	ifstream Dist("BlobSize.in");
	Dist >> binCount;
	int *SizeX, *SizeY, *SizeZ;
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

	// Generate the residual NWP
	if (rank==0) printf("Initializing with saturation (phase 1) = %f \n",Saturation);
	//	GenerateResidual(id,nx,ny,nz,Saturation);

	int x,y,z;
	int sizeX,sizeY,sizeZ;
	int ii,jj,kk;
	int Nx = nx;
	int Ny = ny;
	int Nz = nz;
	float sat = 0.f;
	int Number = 0;		// number of features
	while (sat < Saturation){
		if (rank==0){
			Number++;
			// Randomly generate a point in the domain
			x = (Nx-2)*nprocx*float(rand())/float(RAND_MAX);
			y = (Ny-2)*nprocy*float(rand())/float(RAND_MAX);
			z = (Nz-2)*nprocz*float(rand())/float(RAND_MAX);

			bin = int(floor(binCount*float(rand())/float(RAND_MAX)));
			sizeX = SizeX[bin];
			sizeY = SizeY[bin];
			sizeZ = SizeZ[bin];
		}
		MPI_Bcast(&x,1,MPI_INT,0,comm);
		MPI_Bcast(&y,1,MPI_INT,0,comm);
		MPI_Bcast(&z,1,MPI_INT,0,comm);
		MPI_Bcast(&sizeX,1,MPI_INT,0,comm);
		MPI_Bcast(&sizeY,1,MPI_INT,0,comm);
		MPI_Bcast(&sizeZ,1,MPI_INT,0,comm);

		//if (rank==0) printf("Broadcast block at %i,%i,%i \n",x,y,z);

		for (k=z;k<z+sizeZ;k++){
			for (j=y;j<y+sizeY;j++){
				for (i=x;i<x+sizeX;i++){

					// Identify nodes in the domain (periodic BC)
					ii = i;
					jj = j;
					kk = k;

					if (ii>nprocx*(Nx-2)) ii-=nprocx*(Nx-2);
					if (jj>nprocy*(Ny-2)) jj-=nprocy*(Ny-2);
					if (kk>nprocz*(Nz-2)) kk-=nprocz*(Nz-2);

					// Check if this is in the subdomain
					if (ii < (iproc+1)*(Nx-2)+1 && jj < (jproc+1)*(Ny-2)+1 && kk < (kproc+1)*(Nz-2)+1 &&
							ii  > iproc*(Nx-2) && jj > jproc*(Ny-2) && kk > kproc*(Nz-2) ){

						// Map from global to local coordinates
						ii -= iproc*(Nx-2);
						jj -= jproc*(Ny-2);
						kk -= kproc*(Nz-2);

						n = kk*Nx*Ny+jj*Nx+ii;

						if (id[n] == 2){
							id[n] = 1;
							//count++;
						}

					}
				}
			}
		}
		count = 0;
		for (int k=1; k<Nz-1; k++){
			for (int j=1; j<Ny-1; j++){
				for (int i=1; i<Nx-1; i++){
					n=k*Nx*Ny+j*Nx+i;
					if (id[n] == 1){
						count++;
					}
				}
			}
		}
		MPI_Allreduce(&count,&countGlobal,1,MPI_INT,MPI_SUM,comm);
		sat = float(countGlobal)/totalGlobal;
		//if (rank==0) printf("New count=%i\n",countGlobal);
		//if (rank==0) printf("New saturation=%f\n",sat);
	}

	if (InitialWetting == 1)	FlipID(id,nx*ny*nz);

	// Fill in the phase ID from neighboring processors
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
	MPI_Sendrecv(sendID_x,Dm.sendCount_x,MPI_CHAR,Dm.rank_x(),sendtag,
			recvID_X,Dm.recvCount_X,MPI_CHAR,Dm.rank_X(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_X,Dm.sendCount_X,MPI_CHAR,Dm.rank_X(),sendtag,
			recvID_x,Dm.recvCount_x,MPI_CHAR,Dm.rank_x(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_y,Dm.sendCount_y,MPI_CHAR,Dm.rank_y(),sendtag,
			recvID_Y,Dm.recvCount_Y,MPI_CHAR,Dm.rank_Y(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Y,Dm.sendCount_Y,MPI_CHAR,Dm.rank_Y(),sendtag,
			recvID_y,Dm.recvCount_y,MPI_CHAR,Dm.rank_y(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_z,Dm.sendCount_z,MPI_CHAR,Dm.rank_z(),sendtag,
			recvID_Z,Dm.recvCount_Z,MPI_CHAR,Dm.rank_Z(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Z,Dm.sendCount_Z,MPI_CHAR,Dm.rank_Z(),sendtag,
			recvID_z,Dm.recvCount_z,MPI_CHAR,Dm.rank_z(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xy,Dm.sendCount_xy,MPI_CHAR,Dm.rank_xy(),sendtag,
			recvID_XY,Dm.recvCount_XY,MPI_CHAR,Dm.rank_XY(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XY,Dm.sendCount_XY,MPI_CHAR,Dm.rank_XY(),sendtag,
			recvID_xy,Dm.recvCount_xy,MPI_CHAR,Dm.rank_xy(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xy,Dm.sendCount_Xy,MPI_CHAR,Dm.rank_Xy(),sendtag,
			recvID_xY,Dm.recvCount_xY,MPI_CHAR,Dm.rank_xY(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xY,Dm.sendCount_xY,MPI_CHAR,Dm.rank_xY(),sendtag,
			recvID_Xy,Dm.recvCount_Xy,MPI_CHAR,Dm.rank_Xy(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xz,Dm.sendCount_xz,MPI_CHAR,Dm.rank_xz(),sendtag,
			recvID_XZ,Dm.recvCount_XZ,MPI_CHAR,Dm.rank_XZ(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XZ,Dm.sendCount_XZ,MPI_CHAR,Dm.rank_XZ(),sendtag,
			recvID_xz,Dm.recvCount_xz,MPI_CHAR,Dm.rank_xz(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xz,Dm.sendCount_Xz,MPI_CHAR,Dm.rank_Xz(),sendtag,
			recvID_xZ,Dm.recvCount_xZ,MPI_CHAR,Dm.rank_xZ(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xZ,Dm.sendCount_xZ,MPI_CHAR,Dm.rank_xZ(),sendtag,
			recvID_Xz,Dm.recvCount_Xz,MPI_CHAR,Dm.rank_Xz(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yz,Dm.sendCount_yz,MPI_CHAR,Dm.rank_yz(),sendtag,
			recvID_YZ,Dm.recvCount_YZ,MPI_CHAR,Dm.rank_YZ(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_YZ,Dm.sendCount_YZ,MPI_CHAR,Dm.rank_YZ(),sendtag,
			recvID_yz,Dm.recvCount_yz,MPI_CHAR,Dm.rank_yz(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Yz,Dm.sendCount_Yz,MPI_CHAR,Dm.rank_Yz(),sendtag,
			recvID_yZ,Dm.recvCount_yZ,MPI_CHAR,Dm.rank_yZ(),recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yZ,Dm.sendCount_yZ,MPI_CHAR,Dm.rank_yZ(),sendtag,
			recvID_Yz,Dm.recvCount_Yz,MPI_CHAR,Dm.rank_Yz(),recvtag,comm,MPI_STATUS_IGNORE);
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
		count = 0;
		for (int k=1; k<Nz-1; k++){
			for (int j=1; j<Ny-1; j++){
				for (int i=1; i<Nx-1; i++){
					n=k*Nx*Ny+j*Nx+i;
					if (id[n] == 1){
						count++;
					}
				}
			}
		}
		MPI_Allreduce(&count,&countGlobal,1,MPI_INT,MPI_SUM,comm);
		sat = float(countGlobal)/totalGlobal;
	if (rank==0) printf("Final saturation=%f\n",sat);

	sprintf(LocalRankFilename,"ID.%05i",rank);
	FILE *ID = fopen(LocalRankFilename,"wb");
	fwrite(id,1,N,ID);
	fclose(ID);

	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;
}
