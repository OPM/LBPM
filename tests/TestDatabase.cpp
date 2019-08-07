#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <memory>

#include "common/UnitTest.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "common/Database.h"
#include "ProfilerApp.h"


// Main
int main(int argc, char **argv)
{
    int rank,nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);
    Utilities::setAbortBehavior(true,2);
    Utilities::setErrorHandlers();
    UnitTest ut;
    int err=0;
    
    int BC=2;
    int npx=1; int npy=2; int npz=4;
    int nx=32; int ny=34; int nz=35;
    
    // write a simple database and test that it can be read by LBPM
    auto db = std::make_shared<Database>( );
    db->putScalar<int>( "BC", BC );
    db->putVector<int>( "nproc", { npx, npy, npz } );
    db->putVector<int>( "n", { nx, ny, nz } );

    std::ofstream OutStream("test.db");
//    db->putDatabase();
    OutStream << "Domain { \n";
    db->print(OutStream, "   ");
    OutStream << "} \n";
	printf("TestDatbase: writing test file\n");
	OutStream.close();

	auto new_db = std::make_shared<Database>( "test.db" );
	auto domain_db = new_db->getDatabase( "Domain" );
	if (domain_db->keyExists( "BC" )){
		auto newBC = domain_db->getScalar<int>( "BC" );
		if (newBC != BC){
			err=1;
			printf("TestDatbase: getScalar failed! \n");
		}
	}
	if (err==0) printf("TestDatabase: succeeded!\n");
	else 		printf("TestDatabase: errors encountered \n");

    // Finished
    PROFILE_SAVE("TestDatabase",true);
    MPI_Barrier(comm);
    MPI_Finalize();
    return err;
}


