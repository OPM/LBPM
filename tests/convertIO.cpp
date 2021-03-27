#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/MPI.h"
#include "common/Utilities.h"
#include "IO/Mesh.h"
#include "IO/Reader.h"
#include "IO/Writer.h"
#include "ProfilerApp.h"

int main(int argc, char **argv)
{
    // Initialize MPI
    Utilities::startup( argc, argv );
    Utilities::setErrorHandlers();
    PROFILE_ENABLE(2);
    PROFILE_ENABLE_TRACE();
    PROFILE_START("Main");

    { // Limit scope


        Utilities::MPI comm( MPI_COMM_WORLD );
        // Get inputs
        if ( argc != 5 ) {
            std::cerr << "Error calling convertIO:\n";
            std::cerr << "   convertIO <input_path> <input_format> <output_path> <output_format>\n";
            return -1;
        }
        std::string path_in = argv[1];
        std::string format_in = argv[2];
        std::string path_out = argv[3];
        std::string format_out = argv[4];

        // Check that we have enough ranks to load and write the data
        // This is really only a bottleneck for the writer
        int N_domains = IO::maxDomains( path_in, format_in, comm );
        ASSERT( comm.getSize() == N_domains );

        // Read the timesteps
        auto timesteps = IO::readTimesteps( path_in, format_in );

        // Loop through the timesteps, reading/writing the data
        IO::initialize( path_out, format_out, false );
        for ( auto timestep : timesteps ) {
            
            // Set the domain to read (needs to be the current rank for the writer to be valid)
            int domain = comm.getRank();

            // Get the maximum number of domains for the 
            auto data = IO::readData( path_in, timestep, domain );

            // Save the mesh data to a new file
            IO::writeData( timestep, data, comm );

        }

    } // Limit scope

    // shutdown
    PROFILE_STOP("Main");
    PROFILE_SAVE("convertData",true);
    Utilities::shutdown();
    return 0;
}

