// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "analysis/Minkowski.h"
#include "IO/MeshDatabase.h"

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
        Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
        int rank = comm.getRank();
        int nprocs = comm.getSize();
	{ // Limit scope so variables that contain communicators will free before MPI_Finialize

		if ( rank==0 ) {
			printf("-----------------------------------------------------------\n");
			printf("Unit test 3D topologies  \n");
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
		auto Dm = std::make_shared<Domain>(db,comm);

		Nx += 2;
		Ny += 2;
		Nz += 2;
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

		// Create visualization structure
		std::vector<IO::MeshDataStruct> visData;
		fillHalo<double> fillData(Dm->Comm,Dm->rank_info,{Dm->Nx-2,Dm->Ny-2,Dm->Nz-2},{1,1,1},0,1);;    

		IO::initialize("","silo","false");
		// Create the MeshDataStruct    
		visData.resize(1);
		visData[0].meshName = "domain";
		visData[0].mesh = std::make_shared<IO::DomainMesh>( Dm->rank_info,Dm->Nx-2,Dm->Ny-2,Dm->Nz-2,Dm->Lx,Dm->Ly,Dm->Lz );
		auto PhaseVar = std::make_shared<IO::Variable>();
		PhaseVar->name = "phase";
		PhaseVar->type = IO::VariableType::VolumeVariable;
		PhaseVar->dim = 1;
		PhaseVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
		visData[0].vars.push_back(PhaseVar);

		//.......................................................................
		// Assign the phase ID field based and the signed distance
		//.......................................................................
		double R1,R2,R;
		double CX,CY,CZ; //CY1,CY2;
		CX=Nx*nprocx*0.5;
		CY=Ny*nprocy*0.5;
		CZ=Nz*nprocz*0.5;
		R1 = (Nx-2)*nprocx*0.3; // middle radius
		R2 = (Nx-2)*nprocx*0.1; // donut thickness
		R = 0.4*nprocx*(Nx-2);

		Minkowski Object(Dm);

		int timestep = 0;
		double x,y,z;

		// partial torus
		timestep += 1;
		for ( k=1;k<Nz-1;k++){
			for ( j=1;j<Ny-1;j++){
				for ( i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;

					// global position relative to center
					x = Dm->iproc()*(Nx-2)+i - CX - 0.1;
					y = Dm->jproc()*(Ny-2)+j - CY - 0.1;
					z = Dm->kproc()*(Nz-2)+k - CZ -0.1;

					//..............................................................................
					if (x <= 0 || y<=0) {
						// Single torus
						Object.distance(i,j,k) = R2 - sqrt((sqrt(x*x+y*y) - R1)*(sqrt(x*x+y*y) - R1) + z*z);
					}
					else {
						double d1 = R2-sqrt(x*x +(y-R1)*(y-R1) + z*z);
						double d2 = R2-sqrt((x-R1)*(x-R1)+y*y + z*z);
						Object.distance(i,j,k) = max(d1,d2);
					}

					if (Object.distance(i,j,k) > 0.0){
						Dm->id[n] = 2;
						Object.id(i,j,k) = 2;
					}
					else{
						Dm->id[n] = 1;
						Object.id(i,j,k) = 1;
					}
				}
			}
		}

		ASSERT(visData[0].vars[0]->name=="phase");
		Array<double>& PhaseData = visData[0].vars[0]->data;
		fillData.copy(Object.distance,PhaseData);
		IO::writeData( timestep, visData, comm );

		//spherical shell
		timestep += 1;
		for ( k=1;k<Nz-1;k++){
			for ( j=1;j<Ny-1;j++){
				for ( i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;

					// global position relative to center
					x = Dm->iproc()*(Nx-2)+i - CX - 0.1;
					y = Dm->jproc()*(Ny-2)+j - CY - 0.1;
					z = Dm->kproc()*(Nz-2)+k - CZ - 0.1;

					//..............................................................................
					// Single torus
					double d1 = sqrt(x*x+y*y+z*z)-(R1-R2);
					double d2 = R-sqrt(x*x+y*y+z*z);
					Object.distance(i,j,k) = min(d1,d2);

					if (Object.distance(i,j,k) > 0.0){
						Dm->id[n] = 2;
						Object.id(i,j,k) = 2;
					}
					else{
						Dm->id[n] = 1;
						Object.id(i,j,k) = 1;
					}
				}
			}
		}
		ASSERT(visData[0].vars[0]->name=="phase");
		PhaseData = visData[0].vars[0]->data;
		fillData.copy(Object.distance,PhaseData);
		IO::writeData( timestep, visData, comm );

		// bowl
		timestep += 1;
		for ( k=1;k<Nz-1;k++){
			for ( j=1;j<Ny-1;j++){
				for ( i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;

					// global position relative to center
					x = Dm->iproc()*(Nx-2)+i - CX - 0.1;
					y = Dm->jproc()*(Ny-2)+j - CY - 0.1;
					z = Dm->kproc()*(Nz-2)+k - CZ - 0.1;

					//..............................................................................
					// Bowl
					if (z > 0 ){
						Object.distance(i,j,k) = R2-sqrt((sqrt(x*x+y*y) - R1)*(sqrt(x*x+y*y) - R1) + z*z);
					}
					else
					{
						double d1 = sqrt(x*x+y*y+z*z)-(R1-R2);
						double d2 = R-sqrt(x*x+y*y+z*z);
						Object.distance(i,j,k) = min(d1,d2);
					}
					if (Object.distance(i,j,k) > 0.0){
						Dm->id[n] = 2;
						Object.id(i,j,k) = 2;
					}
					else{
						Dm->id[n] = 1;
						Object.id(i,j,k) = 1;
					}
				}
			}
		}

		ASSERT(visData[0].vars[0]->name=="phase");
		PhaseData = visData[0].vars[0]->data;
		fillData.copy(Object.distance,PhaseData);
		IO::writeData( timestep, visData, comm );

	} // Limit scope so variables that contain communicators will free before MPI_Finialize
        Utilities::shutdown();
	return 0;  
}

