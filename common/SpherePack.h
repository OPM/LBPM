/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
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
#include "common/MPI.h"
#include "common/Communication.h"
#include "common/Database.h"

/*
Simple tools to work with sphere packs
 */

void WriteLocalSolidID(char *FILENAME, char *ID, int N);

void WriteLocalSolidDistance(char *FILENAME, double *Distance, int N);

void ReadSpherePacking(int nspheres, double *List_cx, double *List_cy,
                       double *List_cz, double *List_rad);

void AssignLocalSolidID(char *ID, int nspheres, double *List_cx,
                        double *List_cy, double *List_cz, double *List_rad,
                        double Lx, double Ly, double Lz, int Nx, int Ny, int Nz,
                        int iproc, int jproc, int kproc, int nprocx, int nprocy,
                        int nprocz);

void SignedDistance(double *Distance, int nspheres, double *List_cx,
                    double *List_cy, double *List_cz, double *List_rad,
                    double Lx, double Ly, double Lz, int Nx, int Ny, int Nz,
                    int iproc, int jproc, int kproc, int nprocx, int nprocy,
                    int nprocz);

#endif
