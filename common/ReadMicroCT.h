#ifndef READMICROCT_H
#define READMICROCT_H


#include "common/Array.h"
#include "common/Communication.h"
#include "common/Database.h"


Array<uint8_t> readMicroCT( const std::string& filename );

Array<uint8_t> readMicroCT( const Database& domain, MPI_Comm comm );


#endif
