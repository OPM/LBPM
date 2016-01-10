// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "common/TwoPhase.h"

//#include "Domain.h"

using namespace std;



int main(int argc, char **argv)
{
  // Initialize MPI
  int rank, nprocs;
  MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&nprocs);
  { // Limit scope so variables that contain communicators will free before MPI_Finialize

    if ( rank==0 ) {
        printf("-----------------------------------------------------------\n");
        printf("Unit test for torus (Euler-Poincarie characteristic) \n");
        printf("-----------------------------------------------------------\n");
    }

    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
    double Lx, Ly, Lz;
    int Nx,Ny,Nz;
    int i,j,k,n;

    if (rank==0){
  /*  	ifstream domain("Domain.in");
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
    	*/
    	// Set the domain for single processor test
    	nprocx=nprocy=nprocz=1;
    	nx=ny=nz=100;
    	nspheres=1;
    	Lx=Ly=Lz=1;
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
    if ( nprocs < nprocx*nprocy*nprocz )
        ERROR("Insufficient number of processors");

    // Set up the domain
	int BC=0;
    // Get the rank info
	Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
 //   const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
	TwoPhase Averages(Dm);
	int N = (nx+2)*(ny+2)*(nz+2);
	Nx = nx+2;
	Ny = ny+2;
	Nz = nz+2;
	if (rank == 0) cout << "Domain set." << endl;
	//.......................................................................
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				Dm.id[n] = 1;
			}
		}
	}
	//.......................................................................
    Dm.CommInit(comm); // Initialize communications for domains
	//.......................................................................

	//.......................................................................
	// Assign the phase ID field based and the signed distance
	//.......................................................................
    double R1,R2;
    double CX,CY,CZ'//CY1,CY2;
    CX=Nx*nprocx*0.5;
    CY=Ny*nprocy*0.5;
    CZ=Nz*nprocz*0.5;
    R1 = Nx*nprocx*0.2; // middle radius
    R2 = Nx*nprocx*0.1; // donut thickness
    //
    //CY1=Nx*nprocx*0.5+R1;
    //CY2=Ny*nprocy*0.5-R1;

    double x,y,z;
	if (rank==0) printf("Initializing the system \n");
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;

				// global position relative to center
				x = Dm.iproc*Nx+i - CX;
				y = Dm.jproc*Ny+j - CY;
				z = Dm.kproc*Nz+k - CZ;

				// Shrink the sphere sizes by two voxels to make sure they don't touch
				Averages.SDs(i,j,k) = 100.0;
				//..............................................................................
				// Single torus
				Averages.Phase(i,j,k) = sqrt((sqrt(x*x+y*y) - R1)*(sqrt(x*x+y*y) - R1) + z*z) - R2;
				// Double torus
				/*		y = Dm.jproc*Ny+j - CY1;
				//z = Dm.kproc*Nz+k - CZ +R1;
				Averages.Phase(i,j,k) = sqrt((sqrt(x*x+y*y) - R1)*(sqrt(x*x+y*y) - R1) + z*z) - R2;
				
				y = Dm.jproc*Ny+j - CY2;
				//z = Dm.kproc*Nz+k - CZ-R1;
				Averages.Phase(i,j,k) = min(Averages.Phase(i,j,k),
						sqrt((sqrt(x*x+y*y) - R1)*(sqrt(x*x+y*y) - R1) + z*z) - R2);
				*///..............................................................................

				//Averages.Phase(i,j,k) = - Averages.Phase(i,j,k);
				if (Averages.Phase(i,j,k) > 0.0){
					Dm.id[n] = 2;
				}
				else{
					Dm.id[n] = 1;
				}
				Averages.SDn(i,j,k) = Averages.Phase(i,j,k);
				Averages.Phase(i,j,k) = Averages.SDn(i,j,k);
				Averages.Phase_tplus(i,j,k) = Averages.SDn(i,j,k);
				Averages.Phase_tminus(i,j,k) = Averages.SDn(i,j,k);
				Averages.DelPhi(i,j,k) = 0.0;
				Averages.Press(i,j,k) = 0.0;
				Averages.Vel_x(i,j,k) = 0.0;
				Averages.Vel_y(i,j,k) = 0.0;
				Averages.Vel_z(i,j,k) = 0.0;
			}
		}
	}

	double vS;
	vS = 0.0;

    double beta = 0.95;
	if (rank==0) printf("initializing the system \n");

	Averages.UpdateSolid();
    Dm.CommunicateMeshHalo(Averages.Phase);
    Dm.CommunicateMeshHalo(Averages.SDn);

    Averages.Initialize();
    Averages.UpdateMeshValues();

	if (rank==0) printf("computing local averages  \n");
	Averages.AssignComponentLabels();
    Averages.ComponentAverages();
    Averages.PrintComponents(int(5));
	if (rank==0) printf("reducing averages  \n");


    Averages.Initialize();
    Averages.ComputeLocal();
    Averages.Reduce();
    Averages.PrintAll(int(5));
   // Averages.Reduce();

  } // Limit scope so variables that contain communicators will free before MPI_Finialize
  MPI_Barrier(comm);
  MPI_Finalize();
  return 0;  
}

