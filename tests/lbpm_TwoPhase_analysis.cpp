/*
 * Compute porous media marching cubes analysis (PMMC) based on input image data 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>

#include "common/Array.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "common/MPI.h"
#include "IO/MeshDatabase.h"
#include "IO/Mesh.h"
#include "IO/Writer.h"
#include "IO/netcdf.h"
#include "analysis/analysis.h"
#include "analysis/filters.h"
#include "analysis/distance.h"
#include "analysis/Minkowski.h"
#include "analysis/TwoPhase.h"

#include "ProfilerApp.h"

int main(int argc, char **argv)
{

	// Initialize MPI
        Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
        int rank = comm.getRank();
        //int nprocs = comm.getSize();
	{
		Utilities::setErrorHandlers();
		PROFILE_START("Main");

		//std::vector<std::string> filenames;
		if ( argc<2 ) {
			if ( rank == 0 ){
				printf("At least one filename must be specified\n");
			}
			return 1;
		}
		std::string filename = std::string(argv[1]);
		if ( rank == 0 ){
			printf("Input data file: %s\n",filename.c_str());
		}

		auto db = std::make_shared<Database>( filename );
		auto domain_db = db->getDatabase( "Domain" );
		auto vis_db = db->getDatabase( "Visualization" );

		// Read domain parameters
		auto Filename = domain_db->getScalar<std::string>( "Filename" );
		auto size = domain_db->getVector<int>( "n" );
		auto SIZE = domain_db->getVector<int>( "N" );
		auto nproc = domain_db->getVector<int>( "nproc" );
		auto ReadValues = domain_db->getVector<char>( "ReadValues" );
		auto WriteValues = domain_db->getVector<char>( "WriteValues" );

		auto nx = size[0];
		auto ny = size[1];
		auto nz = size[2];
		int i,j,k,n;

		std::shared_ptr<Domain> Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
		comm.barrier();
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n = k*nx*ny+j*nx+i;
					// Initialize the object
					Dm->id[n] = 1;
				}
			}
		}
		Dm->CommInit();
		
		/* read the data */
		if (domain_db->keyExists( "Filename" )){
			auto Filename = domain_db->getScalar<std::string>( "Filename" );
			Dm->Decomp(Filename);
		}
		else{
			Dm->ReadIDs();
		}
		
	    fillHalo<double> fillDouble(Dm->Comm, Dm->rank_info, {nx,ny,nz}, {1, 1, 1}, 0, 1);
	    		
		// Compute the Minkowski functionals
		comm.barrier();				
		std::shared_ptr<TwoPhase> Averages( new TwoPhase(Dm) );
		std::shared_ptr<Minkowski> Geometry(new Minkowski(Dm));

		// Calculate the distance		
		// Initialize the domain and communication
		nx+=2; ny+=2; nz+=2;
		Array<char> id(nx,ny,nz);
		DoubleArray Distance(nx,ny,nz);

		/* * * * * * * * * * * * * * * * * * * * * */
		// Solve for the position of the solid phase
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n = k*nx*ny+j*nx+i;
					// Initialize the object
					if (Dm->id[n] == ReadValues[0])		id(i,j,k) = 0;
					else		      					id(i,j,k) = 1;
				}
			}
		}
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n=k*nx*ny+j*nx+i;
					// Initialize distance to +/- 1
					Averages->SDs(i,j,k) = 2.0*double(id(i,j,k))-1.0;
				}
			}
		}
		fillDouble.fill(Averages->SDs);
		if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
		CalcDist(Averages->SDs,id,*Dm);
		
		/* move the solid surface by a voxel to improve contact angle measurement
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n=k*nx*ny+j*nx+i;
					Averages->SDs(i,j,k) -= 1.0;
				}
			}
		}*/
		/* * * * * * * * * * * * * * * * * * * * * */
		// Solve for the position of the fluid phase
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n = k*nx*ny+j*nx+i;
					// Initialize the object
					if (Dm->id[n] == ReadValues[1]){
						id(i,j,k) = 1;
						Averages->Phase(i,j,k) = 1.0;
					}
					else	{
						id(i,j,k) = 0;
						Averages->Phase(i,j,k) = -1.0;
					}
				}
			}
		}
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n=k*nx*ny+j*nx+i;
					// Initialize distance to +/- 1
					Distance(i,j,k) = 2.0*double(id(i,j,k))-1.0;
				}
			}
		}
		fillDouble.fill(Averages->Phase);
		fillDouble.fill(Distance);

		if (rank==0) printf("Initialized fluid phase -- Converting to Signed Distance function \n");
		CalcDist(Distance,id,*Dm);
		Mean3D(Distance,Averages->SDn);
		/* * * * * * * * * * * * * * * * * * * * * */

		if (rank==0) printf("Computing Minkowski functionals \n");
		Geometry->ComputeScalar(Distance,0.f);
		Geometry->PrintAll();

		Averages->Initialize();
		Averages->UpdateMeshValues();
		Averages->ComputeStatic();
		Averages->Reduce();
		Averages->PrintStatic();
		/*        
        Averages.Initialize();
        Averages.ComponentAverages();
        Averages.SortBlobs();
        Averages.PrintComponents(timestep);
		 */      

		auto format = vis_db->getWithDefault<string>("format", "hdf5");

		vis_db = db->getDatabase("Visualization");

		std::vector<IO::MeshDataStruct> visData;
		fillHalo<double> fillData(Dm->Comm, Dm->rank_info,
				{Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2},
				{1, 1, 1}, 0, 1);

		auto PhaseVar = std::make_shared<IO::Variable>();
		auto SignDistVar = std::make_shared<IO::Variable>();

		IO::initialize("", format, "false");
		// Create the MeshDataStruct
		visData.resize(1);
		visData[0].meshName = "domain";
		visData[0].mesh = std::make_shared<IO::DomainMesh>(
				Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2, Dm->Lx, Dm->Ly,
				Dm->Lz);
		SignDistVar->name = "SignDist";
		SignDistVar->type = IO::VariableType::VolumeVariable;
		SignDistVar->dim = 1;
		SignDistVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
		visData[0].vars.push_back(SignDistVar);

		PhaseVar->name = "Phase";
		PhaseVar->type = IO::VariableType::VolumeVariable;
		PhaseVar->dim = 1;
		PhaseVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
		visData[0].vars.push_back(PhaseVar);

		Array<double> &SignData = visData[0].vars[0]->data;
		Array<double> &PhaseData = visData[0].vars[1]->data;

		ASSERT(visData[0].vars[0]->name == "SignDist");
		ASSERT(visData[0].vars[1]->name == "Phase");

		fillData.copy(Averages->SDs, SignData);
		fillData.copy(Averages->SDn, PhaseData);

		int timestep = 0;
		IO::writeData(timestep, visData, Dm->Comm);
	}
	PROFILE_STOP("Main");
	PROFILE_SAVE("TwoPhase",true);
	Utilities::shutdown();
	return 0;
}


