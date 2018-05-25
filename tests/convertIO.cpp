#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "common/Utilities.h"
#include "IO/Mesh.h"
#include "IO/Reader.h"
#include "IO/Writer.h"
#include "ProfilerApp.h"


int main(int argc, char **argv)
{
  // Initialize MPI
  int rank,nprocs;
  MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&nprocs);
  Utilities::setErrorHandlers();
  PROFILE_ENABLE(2);
  PROFILE_ENABLE_TRACE();
  PROFILE_START("Main");
  { // Limit scope

    // Get inputs
    if ( argc != 3 ) {
        std::cerr << "Error calling convertIO:\n";
        std::cerr << "   convertIO input_file format\n";
        return -1;
    }
    std::string filename = argv[1];
    std::string format = argv[2];
    std::string path = IO::getPath( filename );

    // Read the timesteps
    auto timesteps = IO::readTimesteps( filename );

    // Loop through the timesteps, reading/writing the data
    IO::initialize( "", format, false );
    for ( auto timestep : timesteps ) {
        
        // Read the list of MeshDatabase
        auto databases = IO::getMeshList( path, timestep );

        // Build the MeshDataStruct
        std::vector<IO::MeshDataStruct> meshData(databases.size());

        // Loop through the database
        int i = 0;
        PROFILE_START("Read");
        for ( const auto& database : databases ) {
            
            // Read the appropriate mesh domain
            ASSERT( (int) database.domains.size() == nprocs );
            meshData[i].meshName = database.name;
            meshData[i].mesh = IO::getMesh( path, timestep, database, rank );

            // Read the variables
            for ( auto var : database.variables ) {
                auto varData = IO::getVariable( path, timestep, database, rank, var.name );
                IO::reformatVariable( *meshData[i].mesh, *varData );
                meshData[i].vars.push_back( varData );
            }

            i++;
        }
        MPI_Barrier(comm);
        PROFILE_STOP("Read");

        // Save the mesh data to a new file
        PROFILE_START("Write");
        IO::writeData( timestep, meshData, MPI_COMM_WORLD );
        MPI_Barrier(comm);
        PROFILE_STOP("Write");
    }

  } // Limit scope
  PROFILE_STOP("Main");
  PROFILE_SAVE("convertData",true);
  MPI_Barrier(comm);
  MPI_Finalize();
  return 0;
}

