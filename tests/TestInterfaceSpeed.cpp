#include <iostream>
#include <math.h>

#include "analysis/TwoPhase.h"
#include "common/MPI.h"
#include "common/Communication.h"
#include "IO/Mesh.h"
#include "IO/Writer.h"
#include "ProfilerApp.h"

#define RADIUS 15
#define CAPRAD 20
#define HEIGHT 15.5
#define N 60
#define SPEED -1
#define PI 3.14159

int main (int argc, char *argv[])
{
	// Initialize MPI
    Utilities::startup( argc, argv );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int toReturn = 0;
    {
    // Load inputs
	string FILENAME = argv[1];
    // Load inputs
	if (rank==0)	printf("Loading input database \n");
    auto db = std::make_shared<Database>( FILENAME );
    auto domain_db = db->getDatabase( "Domain" );
    int Nx = domain_db->getVector<int>( "n" )[0];
    int Ny = domain_db->getVector<int>( "n" )[1];
    int Nz = domain_db->getVector<int>( "n" )[2];

    auto Dm = std::make_shared<Domain>(domain_db,comm);

    Nx+=2; Ny+=2; Nz+=2;

	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;

	Dm->CommInit();

	auto Averages = std::make_shared<TwoPhase>(Dm);
	int timestep=0;

	double Cx,Cy,Cz;
	double dist1,dist2;

	Cx = Cy = Cz = N*0.5;
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;
				dist2 = fabs(Cz-k)-HEIGHT;

				Averages->Phase_tminus(i,j,k) = dist2;
			}
		} 
	}
	Cz += SPEED;
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				
				dist1 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)) - RADIUS;
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;
				dist2 = fabs(Cz-k)-HEIGHT;

				Averages->SDs(i,j,k) = -dist1;
				Averages->Phase(i,j,k) = dist2;
				Averages->SDn(i,j,k) = dist2;
			}
		}   
	}
	Cz += SPEED;
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				dist2 = sqrt((i-Cx)*(i-Cx)+(j-Cy)*(j-Cy)+(k-Cz)*(k-Cz)) - CAPRAD;
				dist2 = fabs(Cz-k)-HEIGHT;

				Averages->Phase_tplus(i,j,k) = dist2;
			}
		}   
	}
	
	//...........................................................................

	//....................................................................
	// The following only need to be done once
	//Averages->SetupCubes(Dm);
	Averages->UpdateSolid(); 	// unless the solid is deformable!
	//....................................................................
	// The following need to be called each time new averages are computed
	Averages->Initialize();
	Averages->UpdateMeshValues();
	Averages->ComputeLocal();
	Averages->Reduce();
	Averages->PrintAll(timestep);
	//....................................................................
	
	printf("-------------------------------- \n");
	printf("NWP volume = %f \n", Averages->nwp_volume);
	printf("Area wn = %f, Analytical = %f \n", Averages->awn,2*PI*RADIUS*RADIUS);
	printf("Area ns = %f, Analytical = %f \n", Averages->ans, 2*PI*RADIUS*(Nz-2)-4*PI*RADIUS*HEIGHT);
	printf("Area ws = %f, Analytical = %f \n", Averages->aws, 4*PI*RADIUS*HEIGHT);
	printf("Area s = %f, Analytical = %f \n", Averages->As, 2*PI*RADIUS*(Nz-2));
	printf("Length wns = %f, Analytical = %f \n", Averages->lwns, 4*PI*RADIUS);
	printf("Geodesic curvature (wn) = %f, Analytical = %f \n", Averages->KGwns_global, 0.0);
	printf("Geodesic curvature (ws) = %f, Analytical = %f \n", Averages->KNwns_global, 4*PI);
//	printf("Cos(theta_wns) = %f, Analytical = %f \n",efawns/lwns,1.0*RADIUS/CAPRAD);
	printf("Interface Velocity = %f,%f,%f \n",Averages->vawn_global(0),Averages->vawn_global(1),Averages->vawn_global(2));
	printf("Common Curve Velocity = %f,%f,%f \n",Averages->vawns_global(0),Averages->vawns_global(1),Averages->vawns_global(2));
	printf("-------------------------------- \n");	
	//.........................................................................	
	
	if (fabs(Averages->awn - 2*PI*RADIUS*RADIUS)/(2*PI*RADIUS*RADIUS) > 0.02){
		toReturn = 1;
		printf("TestCylinderArea.cpp: error tolerance exceeded for wn area \n");
	}
	if (fabs(Averages->ans - (2*PI*RADIUS*(Nz-2)-4*PI*RADIUS*HEIGHT))/(2*PI*RADIUS*(Nz-2)-4*PI*RADIUS*HEIGHT)> 0.02 ){
		toReturn = 2;
		printf("TestCylinderArea.cpp: error tolerance exceeded for ns area \n");
	}
	if (fabs(Averages->aws - 4*PI*RADIUS*HEIGHT)/(4*PI*RADIUS*HEIGHT) > 0.02 ){
		toReturn = 3;
		printf("TestCylinderArea.cpp: error tolerance exceeded for ws area \n");
	}
	if (fabs(Averages->As - 2*PI*RADIUS*(Nz-2))/(2*PI*RADIUS*(Nz-2)) > 0.02 ){
		toReturn = 4;
		printf("TestCylinderArea.cpp: error tolerance exceeded for solid area \n");
	}
	if (fabs(Averages->lwns - 4*PI*RADIUS)/(4*PI*RADIUS) > 0.02 ){
		toReturn = 5;
		printf("TestCylinderArea.cpp: error tolerance exceeded for common curve length \n");
	}
	if ( fabs(Averages->vawn_global(2)+0.25) > 0.01){
		printf("TestInterfaceSpeed: Error too high for kinematic velocity of wn interface \n");
		toReturn = 6;
	}
	if ( fabs(Averages->vawns_global(2)+0.25) > 0.01){
		printf("TestInterfaceSpeed: Error too high for kinematic velocity of common curve \n");
		toReturn = 7;
	}

	comm.barrier();
    }
        Utilities::shutdown();
	return toReturn;
}
