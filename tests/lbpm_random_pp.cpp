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
#include <Array.h>
#include <Domain.h>

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
		MPI_Bcast(&x,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&y,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&z,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&sizeX,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&sizeY,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&sizeZ,1,MPI_INT,0,MPI_COMM_WORLD);

		if (x+sizeX < (iproc+1)*(Nx-2) && y+sizeY < (jproc+1)*(Ny-2) && z+sizeZ < (kproc+1)*(Nz-2) &&
				x-sizeX  > iproc*(Nx-2) && y-sizeY > jproc*(Nz-2) && z-sizeZ > kproc*(Nz-2) ){

			x -= iproc*(Nx-2);
			y -= jproc*(Ny-2);
			z -= kproc*(Nz-2);

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


int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

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
				printf("Begin from connected wetting phase \n");
			else
				printf("Begin from connected non-wetting phase \n");
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
	MPI_Barrier(MPI_COMM_WORLD);
	// Computational domain
	MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocy,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nspheres,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//.................................................
	MPI_Barrier(MPI_COMM_WORLD);

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
	TwoPhase Averages(Dm);

	nx+=2; ny+=2; nz+=2;
	int N = nx*ny*nz;
    char *id;
    id = new char[N];

    DoubleArray SignDist(nx,ny,nz);
    // Read the signed distance from file
    sprintf(LocalRankFilename,"SignDist.%05i",rank);
    FILE *DIST = fopen(LocalRankFilename,"rb");
    fread(SignDist.get(),8,N,DIST);
    fclose(DIST);

	int count,countGlobal,totalGlobal;
    count = 0;
    for (int k=0; k<nz; k++){
        for (int j=0; j<ny; j++){
            for (int i=0; i<nx; i++){
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
	MPI_Allreduce(&count,&totalGlobal,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    Dm.CommInit(MPI_COMM_WORLD);
    int iproc = Dm.iproc;
    int jproc = Dm.jproc;
    int kproc = Dm.kproc;

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

	// Generate the residual NWP
	if (rank==0) printf("Initializing with NWP saturation = %f \n",Saturation);
	//	GenerateResidual(id,nx,ny,nz,Saturation);

	int x,y,z,bin;
	int sizeX,sizeY,sizeZ;
	float sat = 0.f;
	Number = 0;		// number of features
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
		MPI_Bcast(&x,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&y,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&z,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&sizeX,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&sizeY,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&sizeZ,1,MPI_INT,0,MPI_COMM_WORLD);

		if (x+sizeX < (iproc+1)*(Nx-2) && y+sizeY < (jproc+1)*(Ny-2) && z+sizeZ < (kproc+1)*(Nz-2) &&
				x-sizeX  > iproc*(Nx-2) && y-sizeY > jproc*(Nz-2) && z-sizeZ > kproc*(Nz-2) ){

			x -= iproc*(Nx-2);
			y -= jproc*(Ny-2);
			z -= kproc*(Nz-2);

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
							//count++;
						}

					}
				}
			}
		}
		count = 0;
	    for (int k=0; k<nz; k++){
	        for (int j=0; j<ny; j++){
	            for (int i=0; i<nx; i++){
					if (ID[n] == 1){
						count++;
					}
	            }
	        }
	    }
		MPI_Allreduce(&count,&countGlobal,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		sat = float(countGlobal)/totalGlobal;

	}

	if (InitialWetting == 1)	FlipID(id,nx*ny*nz);



    sprintf(LocalRankFilename,"ID.%05i",rank);
    FILE *ID = fopen(LocalRankFilename,"wb");
    fwrite(id,1,N,ID);
    fclose(ID);

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;
}
