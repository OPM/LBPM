// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "pmmc.h"
//#include "Domain.h"

using namespace std;

inline void ReadCheckpoint(char *FILENAME, double *cDen, double *cDistEven, double *cDistOdd, int N)
{
	int q,n;
	double value;
	ifstream File(FILENAME,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		File.read((char*) &value, sizeof(value));
		cDen[n] = value;
	//	if (n== 66276)	printf("Density a  = %f \n",value);
		File.read((char*) &value, sizeof(value));
		cDen[N+n] = value;
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

inline void SetPeriodicBC(DoubleArray &Scalar, int nx, int ny, int nz){
	
	int i,j,k,in,jn,kn;
	for (k=0; k<nz; k++){
		for (j=0; j<ny; j++){
			for (i=0; i<nx; i++){
				in = i; jn=j; kn=k;
				if (i==0) in = nx-2 ;
				else if (i==nx-1) in = 0;
				if (j==0) jn = ny-2;
				else if (j==ny-1) jn = 0;
				if (k==0) kn = nz-2;
				else if (k==nz-1) kn = 0;	
				Scalar(i,j,k) = Scalar(in,jn,kn);
			}
		}
	}
}
inline void ReadFromRank(char *FILENAME, DoubleArray &Phase, int nx, int ny, int nz, int iproc, int
							jproc, int kproc)
{
	int i,j,k,q,n,N;
	int iglobal,jglobal,kglobal;
	double value;
	double denA,denB;

	N = nx*ny*nz;
	
	double *Den;
	
	Den = new double[2*N];

	ifstream File(FILENAME,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		File.read((char*) &value, sizeof(value));
		Den[2*n] = value;
		//	if (n== 66276)	printf("Density a  = %f \n",value);
		File.read((char*) &value, sizeof(value));
		Den[2*n+1] = value;

		//	if (n== 66276)	printf("Density b  = %f \n",value);
		// Read the even distributions
		for (q=0; q<10; q++){
			File.read((char*) &value, sizeof(value));
		}
		// Read the odd distributions
		for (q=0; q<9; q++){
			File.read((char*) &value, sizeof(value));
		}
	}
	File.close();
	
	// Compute the phase field
	for (k=1; k<nz-1; k++){
		for (j=1; j<ny-1; j++){
			for (i=1; i<nz-1; i++){
				//........................................................................
				n = k*nx*ny+j*nx+i;
				//........................................................................
				denA = Den[n];
				denB = Den[N+n];
				//........................................................................
				// save values in global arrays
				//........................................................................
				iglobal = iproc*(nx-2)+i;
				jglobal = jproc*(ny-2)+j;
				kglobal = kproc*(nz-2)+k;
				//........................................................................				
				Phase(iglobal,jglobal,kglobal) = (denA-denB)/(denA+denB);
				//........................................................................
			}
		}
	}
	
	delete Den;
}

int main(int argc, char **argv)
{
	printf("-----------------------------------------------------------\n");
	printf("Labeling Blobs from Two-Phase Lattice Boltzmann Simulation \n");
	printf("-----------------------------------------------------------\n");

	//.......................................................................
	int nprocx,nprocy,nprocz,nprocs;
	int Nx, Ny, Nz;
	int nx,ny,nz;
	int nspheres;
	double Lx,Ly,Lz;
	//.......................................................................
	int i,j,k,n,p,idx;
	int iproc,jproc,kproc;
	//.......................................................................
	// Reading the domain information file
	//.......................................................................
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
	//.......................................................................

	nx+=2;
	ny+=2;
	nz+=2;
	
	nprocs = nprocx*nprocy*nprocz;
	printf("Number of MPI ranks: %i \n", nprocs);
	Nx = (nx-2)*nprocx+2;
	Ny = (ny-2)*nprocy+2;
	Nz = (nz-2)*nprocz+2;
	printf("Full domain size: %i x %i x %i  \n", Nx,Ny,Nz);
	
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners

	DoubleArray Phase(Nx,Ny,Nz);
	DoubleArray SignDist(Nx,Ny,Nz);
	
	// Filenames used
	char LocalRankString[8];
	char LocalRankFilename[40];
	char LocalRestartFile[40];
	char BaseFilename[20];
	char tmpstr[10];

	int proc,iglobal,kglobal,jglobal;

	double * Temp;
	Temp = new double[nx*ny*nz];
	
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				SignDist(i,j,k) = -100.0;
			}
		}
	}
	
	// read the files and populate main arrays
	for ( kproc=0; kproc<nprocz; kproc++){
		for ( jproc=0; jproc<nprocy; jproc++){
			for ( iproc=0; iproc<nprocx; iproc++){
				
				proc = kproc*nprocx*nprocy + jproc*nprocx + iproc;

				sprintf(LocalRankString,"%05d",proc);
				sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
				ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);	
				for (k=1; k<nz-1; k++){
					for (j=1; j<ny-1; j++){
						for (i=1; i<nz-1; i++){

							//........................................................................
							n = k*nx*ny+j*nx+i;
							//........................................................................
							iglobal = iproc*(nx-2)+i;
							jglobal = jproc*(ny-2)+j;
							kglobal = kproc*(nz-2)+k;
							//........................................................................
							SignDist(iglobal,jglobal,kglobal) = Temp[n];
							//........................................................................
						}
					}
				}
				
				sprintf(LocalRankFilename,"%s%s","Phase.",LocalRankString);
				ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);	
				for (k=1; k<nz-1; k++){
					for (j=1; j<ny-1; j++){
						for (i=1; i<nx-1; i++){

							//........................................................................
							n = k*nx*ny+j*nx+i;
							//........................................................................
							iglobal = iproc*(nx-2)+i;
							jglobal = jproc*(ny-2)+j;
							kglobal = kproc*(nz-2)+k;
							//........................................................................
							Phase(iglobal,jglobal,kglobal) = Temp[n];
							//........................................................................
						}
					}
				}
				
				
			}
		}
	}
	printf("Read %i ranks of Phase, SignDist \n",nprocs);
	
	delete [] Temp;
	
	IntArray GlobalBlobID(Nx,Ny,Nz);

	SetPeriodicBC(SignDist, Nx, Ny, Nz);
	SetPeriodicBC(Phase, Nx, Ny, Nz);
	
	//		FILE *PHASE;
	//PHASE = fopen("Phase.dat","wb");
	//fwrite(Phase.data,8,Nx*Ny*Nz,PHASE);
	//fclose(PHASE);
	
	// Initialize the local blob ID
	// Initializing the blob ID
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				if (SignDist(i,j,k) < 0.0){
					// Solid phase 
					GlobalBlobID(i,j,k) = -2;
				}
				else{
					GlobalBlobID(i,j,k) = -1;
				}
			}
		}
	}
	
	// Compute the porosity
	double porosity=0.0;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				if (SignDist(i,j,k) > 0.0){ 
					porosity += 1.0;
				}
			}
		}
	}
	int N=int(porosity*1.25);
	porosity /= (Nx*Ny*Nz*1.0);
	printf("Media porosity is %f \n",porosity);

	/* ****************************************************************
				IDENTIFY ALL BLOBS: F > vF, S > vS
	****************************************************************** */
	// Find blob domains, number of blobs
	int nblobs = 0;					// number of blobs
	int ncubes = 0;					// total number of nodes in any blob
	IntArray blobs(3,N);	// store indices for blobs (cubes)
	IntArray temp(3,N);	// temporary storage array
	IntArray b(N);          // number of nodes in each blob

	double vF=0.0;
	double vS=0.0;
	double trimdist=1.0;
	printf("Execute blob identification algorithm... \n");
	// Loop over z=0 first -> blobs attached to this end considered "connected" for LB simulation
	i=0;
	int number=0;
	/*	for (k=0;k<1;k++){
		for (j=0;j<Ny;j++){
			if ( Phase(i,j,k) > vF ){
				if ( SignDist(i,j,k) > vS ){
					// node i,j,k is in the porespace
					number = number+ComputeBlob(blobs,nblobs,ncubes,GlobalBlobID,Phase,SignDist,vF,vS,i,j,k,temp);
				}
			}
		}
	}
	// Specify the blob on the z axis
	if (ncubes > 0){
		b(nblobs) = number;
//		BlobList.push_back[number];
		printf("Number of non-wetting phase blobs is: %i \n",nblobs-1);
		nblobs++;
	}
	*/
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				if ( GlobalBlobID(i,j,k) == -1 ){
					if ( Phase(i,j,k) > vF ){
						if ( SignDist(i,j,k) > vS ){
							// node i,j,k is in the porespace
							b(nblobs) = ComputeBlob(blobs,nblobs,ncubes,GlobalBlobID,Phase,SignDist,vF,vS,i,j,k,temp);
							nblobs++;							  
						}
					}
				}
				// Otherwise, this point has already been assigned - ignore

				// Make sure list blob_nodes is large enough
				if ( nblobs > (int)b.length()-1){
					printf("Increasing size of blob list \n");
					b.resize(2*b.length());
				}
			}
		}
	}
	// Go over all cubes again -> add any that do not contain nw phase
	bool add=1;			// Set to false if any corners contain nw-phase ( F > vF)
//	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	int count_in=0,count_out=0;
	int nodx,nody,nodz;
	for (k=0;k<Nz-1;k++){
		for (j=0;j<Ny-1;j++){
			for (i=0;i<Nx-1;i++){
				// Loop over cube corners
				add=1;				// initialize to true - add unless corner occupied by nw-phase
				for (p=0;p<8;p++){
					nodx=i+cube[p][0];
					nody=j+cube[p][1];
					nodz=k+cube[p][2];
					if ( GlobalBlobID(nodx,nody,nodz) > -1 ){
						// corner occupied by nw-phase  -> do not add
						add = 0;
					}
				}
				if ( add == 1 ){
					blobs(0,ncubes) = i;
					blobs(1,ncubes) = j;
					blobs(2,ncubes) = k;
					ncubes++;
					count_in++;
				}
				else { count_out++; }
			}
		}
	}
	b(nblobs) = count_in;
	nblobs++;
	printf("Identified %i blobs. Writing per-process output files. \n",nblobs);

	int sizeLoc = nx*ny*nz;
	int *LocalBlobID;
	LocalBlobID = new int [sizeLoc];

	printf("File size (4 bytes per entry) %i, \n",sizeLoc);
	// read the files and populate main arrays
	for ( kproc=0; kproc<nprocz; kproc++){
		for ( jproc=0; jproc<nprocy; jproc++){
			for ( iproc=0; iproc<nprocx; iproc++){

				proc = kproc*nprocx*nprocy + jproc*nprocx + iproc;

				sprintf(LocalRankString,"%05d",proc);
				sprintf(LocalRankFilename,"%s%s","BlobLabel.",LocalRankString);

				for (k=0; k<nz; k++){
					for (j=0; j<ny; j++){
						for (i=0; i<nx; i++){
							//........................................................................
							n = k*nx*ny+j*nx+i;
							//........................................................................
							iglobal = iproc*(nx-2)+i;
							jglobal = jproc*(ny-2)+j;
							kglobal = kproc*(nz-2)+k;
							// periodic BC
							if (iglobal < 0 ) iglobal+=Nx;
							if (jglobal < 0 ) jglobal+=Ny;
							if (kglobal < 0 ) kglobal+=Nz;
							if (!(iglobal < Nx) ) iglobal-=Nx;
							if (!(jglobal < Ny) ) jglobal-=Ny;
							if (!(kglobal < Nz) ) kglobal-=Nz;
							//........................................................................
							LocalBlobID[n] = GlobalBlobID(iglobal,jglobal,kglobal);
							//........................................................................
						}
					}
				}

				FILE *BLOBLOCAL;
				BLOBLOCAL = fopen(LocalRankFilename,"wb");
				fwrite(&LocalBlobID[0],4,sizeLoc,BLOBLOCAL);
				fclose(BLOBLOCAL);
			}
		}
	}
	printf("Wrote %i ranks of BlobLabel.xxxxx \n",nprocs);

      	FILE *BLOBS;
	BLOBS = fopen("Blobs.dat","wb");
	fwrite(GlobalBlobID.get(),4,Nx*Ny*Nz,BLOBS);
	fclose(BLOBS);
	
}

