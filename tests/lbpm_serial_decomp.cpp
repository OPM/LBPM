/*
 * Pre-processor to generate signed distance function from segmented data
 * segmented data should be stored in a raw binary file as 1-byte integer (type char)
 * will output distance functions for phases
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "common/Array.h"
#include "common/Domain.h"

int main(int argc, char **argv)
{

	int rank=0;
	string filename;
	if (argc > 1)
		filename=argv[1];
	else{
		ERROR("lbpm_serial_decomp: no in put database provided \n");
	}
	int rank_offset=0;

	// read the input database 
	auto db = std::make_shared<Database>( filename );
	auto domain_db = db->getDatabase( "Domain" );
	
	std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));

	Dm->Decomp(domain_db);
	
}
