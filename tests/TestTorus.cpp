// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "analysis/TwoPhase.h"


std::shared_ptr<Database> loadInputs( int nprocs )
{
  //auto db = std::make_shared<Database>( "Domain.in" );
    auto db = std::make_shared<Database>();
    db->putScalar<int>( "BC", 0 );
    db->putVector<int>( "nproc", { 1, 1, 1 } );
    db->putVector<int>( "n", { 100, 100, 100 } );
    db->putScalar<int>( "nspheres", 1 );
    db->putVector<double>( "L", { 1, 1, 1 } );
    return db;
}


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
    int i,j,k,n;

    // Load inputs
    auto db = loadInputs( nprocs );
    int Nx = db->getVector<int>( "n" )[0];
    int Ny = db->getVector<int>( "n" )[1];
    int Nz = db->getVector<int>( "n" )[2];
    int nprocx = db->getVector<int>( "nproc" )[0];
    int nprocy = db->getVector<int>( "nproc" )[1];
    int nprocz = db->getVector<int>( "nproc" )[2];

    if (rank==0){
    	printf("********************************************************\n");
    	printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
    	printf("********************************************************\n");
    }

    // Get the rank info
		std::shared_ptr<Domain> Dm(new Domain(db,comm));
 //   const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
		std::shared_ptr<TwoPhase> Averages(new TwoPhase(Dm));
	Nx += 2;
	Ny += 2;
	Nz += 2;
	int N = Nx*Ny*Nz;
	//.......................................................................
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				Dm->id[n] = 1;
			}
		}
	}
	//.......................................................................
    Dm->CommInit(); // Initialize communications for domains
	//.......................................................................

	//.......................................................................
	// Assign the phase ID field based and the signed distance
	//.......................................................................
    double R1,R2;
    double CX,CY,CZ; //CY1,CY2;
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
				x = Dm->iproc()*Nx+i - CX;
				y = Dm->jproc()*Ny+j - CY;
				z = Dm->kproc()*Nz+k - CZ;

				// Shrink the sphere sizes by two voxels to make sure they don't touch
				Averages->SDs(i,j,k) = 100.0;
				//..............................................................................
				// Single torus
				Averages->Phase(i,j,k) = sqrt((sqrt(x*x+y*y) - R1)*(sqrt(x*x+y*y) - R1) + z*z) - R2;
				// Double torus
				/*		y = Dm->jproc()*Ny+j - CY1;
				//z = Dm->kproc()*Nz+k - CZ +R1;
				Averages->Phase(i,j,k) = sqrt((sqrt(x*x+y*y) - R1)*(sqrt(x*x+y*y) - R1) + z*z) - R2;
				
				y = Dm->jproc()*Ny+j - CY2;
				//z = Dm->kproc()*Nz+k - CZ-R1;
				Averages->Phase(i,j,k) = min(Averages->Phase(i,j,k),
						sqrt((sqrt(x*x+y*y) - R1)*(sqrt(x*x+y*y) - R1) + z*z) - R2);
				*///..............................................................................

				//Averages->Phase(i,j,k) = - Averages->Phase(i,j,k);
				if (Averages->Phase(i,j,k) > 0.0){
					Dm->id[n] = 2;
				}
				else{
					Dm->id[n] = 1;
				}
				Averages->SDn(i,j,k) = Averages->Phase(i,j,k);
				Averages->Phase(i,j,k) = Averages->SDn(i,j,k);
				Averages->Phase_tplus(i,j,k) = Averages->SDn(i,j,k);
				Averages->Phase_tminus(i,j,k) = Averages->SDn(i,j,k);
				Averages->DelPhi(i,j,k) = 0.0;
				Averages->Press(i,j,k) = 0.0;
				Averages->Vel_x(i,j,k) = 0.0;
				Averages->Vel_y(i,j,k) = 0.0;
				Averages->Vel_z(i,j,k) = 0.0;
			}
		}
	}

    double beta = 0.95;
	if (rank==0) printf("initializing the system \n");

	Averages->UpdateSolid();
    Dm->CommunicateMeshHalo(Averages->Phase);
    Dm->CommunicateMeshHalo(Averages->SDn);

    Averages->Initialize();
    Averages->UpdateMeshValues();

	if (rank==0) printf("computing local averages  \n");
	Averages->AssignComponentLabels();
    Averages->ComponentAverages();
    Averages->PrintComponents(int(5));
	if (rank==0) printf("reducing averages  \n");


    Averages->Initialize();
    Averages->ComputeLocal();
    Averages->Reduce();
    Averages->PrintAll(int(5));
   // Averages->Reduce();

  } // Limit scope so variables that contain communicators will free before MPI_Finialize
  MPI_Barrier(comm);
  MPI_Finalize();
  return 0;  
}

