#ifndef SpherePack_INC
#define SpherePack_INC

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <exception>
#include <stdexcept>

#include "common/Array.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "common/Database.h"

/*
Simple tools to work with sphere packs
 */

void WriteLocalSolidID(char *FILENAME, char *ID, int N);

void WriteLocalSolidDistance(char *FILENAME, double *Distance, int N);

void ReadSpherePacking(int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad);

void AssignLocalSolidID(char *ID, int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad,
			double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, 
			int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz);

void SignedDistance(double *Distance, int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad,
		    double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, 
		    int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz);

#endif
