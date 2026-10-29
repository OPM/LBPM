#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <string>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/Database.h"
#include "analysis/TwoPhase.h"
#include "common/MPI.h"
#include "analysis/rev_analysis.h"

//#define WRITE_SURFACES

/*
 * Standalone representative elementary volume analysis
 * James E. McClure 2013-2014
 */

using namespace std;

namespace {
std::string missingHydraulicInputs(const std::string &filename) {
    Database input_db(filename);
    if (!input_db.keyExists("REV"))
        return {};

    auto rev_db = input_db.getDatabase("REV");
    if (!rev_db->getWithDefault<bool>("Hydraulic", false))
        return {};

    std::string missing;
    for (const char *key : {"SimulationParams", "VelocityPath", "VelocityTimestep"}) {
        auto value = rev_db->getWithDefault<std::string>(key, "");
        if (value.find_first_not_of(" \t\r\n") == std::string::npos) {
            if (!missing.empty())
                missing += ", ";
            missing += std::string("REV.") + key;
        }
    }
    return missing;
}
}


int main(int argc, char **argv)
{
	if (argc < 2) {
		std::cerr << "Usage: lbpm_rev_analysis params.db\n";
		return 1;
	}
	// Initialize MPI
    Utilities::startup( argc, argv, false );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();
	int status = 0;
	{
		if (rank == 0){
			printf("********************************************************\n");
			printf("Running REV Analysis \n");
			printf("********************************************************\n");
		}
		comm.barrier();
		
		auto filename = argv[1];
		auto missing = missingHydraulicInputs(filename);
		if (!missing.empty()) {
			if (rank == 0) {
				std::cerr << "lbpm_rev_analysis: Hydraulic = True requires you to provide "
				          << "SimulationParams, VelocityPath, and VelocityTimestep in the REV block.\n"
				          << "Missing or empty settings: " << missing
				          << ". Analysis was not run.\n";
			}
			status = 1;
		} else try {
			ScaLBL_REV REV(rank,nprocs,comm);
			REV.ReadParams(filename, false);
			REV.SetDomain();    // this reads in the domain
			REV.DetRevAnalysis();
			REV.StatRevAnalysis();
			REV.sub_domain_select();
		} catch (const std::exception &error) {
			if (rank == 0)
				std::cerr << "lbpm_rev_analysis: " << error.what() << '\n';
			status = 1;
		}
		cout << flush;
	}
    Utilities::shutdown();
	return status;
}
