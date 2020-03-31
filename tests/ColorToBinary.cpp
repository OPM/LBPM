// Sequential component labeling for two phase systems
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2015

#include <iostream>
#include <math.h>
#include "analysis/analysis.h"
#include "analysis/TwoPhase.h"

using namespace std;

inline void ReadFromRank(char *FILENAME, DoubleArray &Phase, int nx, int ny, int nz, int iproc, int
							jproc, int kproc)
{
	int i,j,k,q,n,N;
	int iglobal,jglobal,kglobal;
	double value;
	double denA,denB;
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	
	N = nx*ny*nz;
	
	double *Den, *DistEven, *DistOdd;
	
	Den = new double[2*N];
	DistEven = new double[10*N];
	DistOdd = new double[9*N];

	ifstream File(FILENAME,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		File.read((char*) &value, sizeof(value));
		Den[n] = value;
		//	if (n== 66276)	printf("Density a  = %f \n",value);
		File.read((char*) &value, sizeof(value));
		Den[N+n] = value;

		//	if (n== 66276)	printf("Density b  = %f \n",value);
		// Read the even distributions
		for (q=0; q<10; q++){
			File.read((char*) &value, sizeof(value));
			DistEven[q*N+n] = value;
		}
		// Read the odd distributions
		for (q=0; q<9; q++){
			File.read((char*) &value, sizeof(value));
			DistOdd[q*N+n] = value;
		}
	}
	File.close();
	
	// Compute the phase field, pressure and velocity
	for (k=1; k<nz-1; k++){
		for (j=1; j<ny-1; j++){
			for (i=1; i<nz-1; i++){
				//........................................................................
				n = k*nx*ny+j*nx+i;
				//........................................................................
				denA = Den[n];
				denB = Den[N+n];
				//........................................................................
				f0 = DistEven[n];
				f2 = DistEven[N+n];
				f4 = DistEven[2*N+n];
				f6 = DistEven[3*N+n];
				f8 = DistEven[4*N+n];
				f10 = DistEven[5*N+n];
				f12 = DistEven[6*N+n];
				f14 = DistEven[7*N+n];
				f16 = DistEven[8*N+n];
				f18 = DistEven[9*N+n];
				//........................................................................
				f1 = DistOdd[n];
				f3 = DistOdd[1*N+n];
				f5 = DistOdd[2*N+n];
				f7 = DistOdd[3*N+n];
				f9 = DistOdd[4*N+n];
				f11 = DistOdd[5*N+n];
				f13 = DistOdd[6*N+n];
				f15 = DistOdd[7*N+n];
				f17 = DistOdd[8*N+n];
				//........................................................................
				//.................Compute the pressure....................................
				value = 0.3333333333333333*(f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17);
				//........................................................................
				//.................Compute the velocity...................................
				//vx = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
				//vy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
				//vz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
				//........................................................................
				// save values in global arrays
				//........................................................................
				iglobal = iproc*(nx-2)+i-1;
				jglobal = jproc*(ny-2)+j-1;
				kglobal = kproc*(nz-2)+k-1;
				//........................................................................				
				Phase(iglobal,jglobal,kglobal) = (denA-denB)/(denA+denB);
			//	Pressure(iglobal,jglobal,kglobal) = value;
			//	Vel_x(iglobal,jglobal,kglobal) = vx;
			//	Vel_y(iglobal,jglobal,kglobal) = vy;
			//	Vel_z(iglobal,jglobal,kglobal) = vz;
				//........................................................................
			}
		}
	}
	
	delete Den;
	delete DistEven;
	delete DistOdd;
}

int main(int argc, char **argv)
{
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);

	printf("----------------------------------------------------------\n");
	printf("Creating single Binary file from restart (8-bit integer)\n");
	printf("ID=0 (solid) ID=1 (non-wetting), ID=2 (wetting) \n");
	printf("----------------------------------------------------------\n");

	if (nprocs != 1) INSIST(nprocs == 1,"Error: ComponentLabel --serial case!");

	//.......................................................................
	int nprocx,nprocy,nprocz;
	int Nx, Ny, Nz;
	int nx,ny,nz;
	int nspheres;
	double Lx,Ly,Lz;
	//.......................................................................
	int i,j,k,n;
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
	int BoundaryCondition=0;
	Nx = (nx-2)*nprocx;
	Ny = (ny-2)*nprocy;
	Nz = (nz-2)*nprocz;
	Domain Dm(Nx,Ny,Nz,rank,1,1,1,Lx,Ly,Lz,BoundaryCondition);

	printf("Full domain size: %i x %i x %i  \n", Nx,Ny,Nz);
	
	DoubleArray Phase(Nx,Ny,Nz);
	DoubleArray SignDist(Nx,Ny,Nz);
	
	// Filenames used
	char LocalRankString[8];
	char LocalRankFilename[40];
	char BaseFilename[20];
	                 
	int proc,iglobal,kglobal,jglobal;

	double * Temp;
	Temp = new double[nx*ny*nz];

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
							iglobal = iproc*(nx-2)+i-1;
							jglobal = jproc*(ny-2)+j-1;
							kglobal = kproc*(nz-2)+k-1;
							//........................................................................
							SignDist(iglobal,jglobal,kglobal) = Temp[n];
							//........................................................................
						}
					}
				}
				
		/*		sprintf(LocalRankFilename,"%s%s","Phase.",LocalRankString);
				ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);
				for (k=1; k<nz-1; k++){
					for (j=1; j<ny-1; j++){
						for (i=1; i<nz-1; i++){

							//........................................................................
							n = k*nx*ny+j*nx+i;
							//........................................................................
							iglobal = iproc*(nx-2)+i-1;
							jglobal = jproc*(ny-2)+j-1;
							kglobal = kproc*(nz-2)+k-1;
							//........................................................................
							Phase(iglobal,jglobal,kglobal) = Temp[n];
							//........................................................................
						}
					}
				}
*/
				sprintf(LocalRankFilename,"%s%s","Restart.",LocalRankString);
				ReadFromRank(LocalRankFilename,Phase,nx,ny,nz,iproc,jproc,kproc);
			}
		}
	}
	printf("Read %i ranks of %s \n",nprocs,BaseFilename);
	delete Temp;
	
	// Initializing the blob ID
	IntArray PhaseID(Nx,Ny,Nz);
	char *ID;
	ID = new char [Nx*Ny*Nz];
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				n = k*Nx*Ny+j*Nx+i;
				if (SignDist(i,j,k) < 0.0){
					// Solid phase 
					ID[n] = 0;
					PhaseID(i,j,k) = 0;
				}
				else if (Phase(i,j,k) < 0.0){
					// wetting phase
					ID[n] = 2;
					PhaseID(i,j,k) = 2;
				}
				else {
					// non-wetting phase
					ID[n] = 1;
					PhaseID(i,j,k) = 1;
				}
			}
		}
	}

	FILE *OUTFILE;
	OUTFILE = fopen("ID.dat","wb");
	fwrite(ID,1,Nx*Ny*Nz,OUTFILE);
	fclose(OUTFILE);

	OUTFILE = fopen("Phase.dat","wb");
	fwrite(Phase.data(),8,Nx*Ny*Nz,OUTFILE);
	fclose(OUTFILE);

	/*	OUTFILE = fopen("SignDist.dat","wb");
	fwrite(SignDist.data(),8,Nx*Ny*Nz,OUTFILE);
	fclose(OUTFILE);
	*/
	
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}

