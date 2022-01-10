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
#include "common/SpherePack.h"

// Inline function to read line without a return argument
static inline void fgetl(char *str, int num, FILE *stream) {
    char *ptr = fgets(str, num, stream);
    if (0) {
        char *temp = (char *)&ptr;
        temp++;
    }
}

void WriteLocalSolidID(char *FILENAME, char *ID, int N) {
    char value;
    ofstream File(FILENAME, ios::binary);
    for (int n = 0; n < N; n++) {
        value = ID[n];
        File.write((char *)&value, sizeof(value));
    }
    File.close();
}

void WriteLocalSolidDistance(char *FILENAME, double *Distance, int N) {
    double value;
    ofstream File(FILENAME, ios::binary);
    for (int n = 0; n < N; n++) {
        value = Distance[n];
        File.write((char *)&value, sizeof(value));
    }
    File.close();
}

void ReadSpherePacking(int nspheres, double *List_cx, double *List_cy,
                       double *List_cz, double *List_rad) {
    // Read in the full sphere pack
    //...... READ IN THE SPHERES...................................
    cout << "Reading the packing file..." << endl;
    FILE *fid = fopen("pack.out", "rb");
    INSIST(fid != NULL, "Error opening pack.out");
    //.........Trash the header lines..........
    char line[100];
    fgetl(line, 100, fid);
    fgetl(line, 100, fid);
    fgetl(line, 100, fid);
    fgetl(line, 100, fid);
    fgetl(line, 100, fid);
    //........read the spheres..................
    // We will read until a blank like or end-of-file is reached
    int count = 0;
    while (!feof(fid) && fgets(line, 100, fid) != NULL) {
        char *line2 = line;
        List_cx[count] = strtod(line2, &line2);
        List_cy[count] = strtod(line2, &line2);
        List_cz[count] = strtod(line2, &line2);
        List_rad[count] = strtod(line2, &line2);
        count++;
    }
    cout << "Number of spheres extracted is: " << count << endl;
    INSIST(count == nspheres,
           "Specified number of spheres is probably incorrect!");
    // .............................................................
}

void AssignLocalSolidID(char *ID, int nspheres, double *List_cx,
                        double *List_cy, double *List_cz, double *List_rad,
                        double Lx, double Ly, double Lz, int Nx, int Ny, int Nz,
                        int iproc, int jproc, int kproc, int nprocx, int nprocy,
                        int nprocz) {
    // Use sphere lists to determine which nodes are in porespace
    // Write out binary file for nodes
    char value;
    int N = Nx * Ny * Nz; // Domain size, including the halo
    double hx, hy, hz;
    double x, y, z;
    double cx, cy, cz, r;
    int imin, imax, jmin, jmax, kmin, kmax;
    int p, i, j, k, n;
    //............................................
    double min_x, min_y, min_z;
    //    double max_x,max_y,max_z;
    //............................................
    // Lattice spacing for the entire domain
    // It should generally be true that hx=hy=hz
    // Otherwise, you will end up with ellipsoids
    hx = Lx / (Nx * nprocx - 1);
    hy = Ly / (Ny * nprocy - 1);
    hz = Lz / (Nz * nprocz - 1);
    //............................................
    // Get maximum and minimum for this domain
    // Halo is included !
    min_x = double(iproc * Nx - 1) * hx;
    min_y = double(jproc * Ny - 1) * hy;
    min_z = double(kproc * Nz - 1) * hz;
    //    max_x = ((iproc+1)*Nx+1)*hx;
    //    max_y = ((jproc+1)*Ny+1)*hy;
    //    max_z = ((kproc+1)*Nz+1)*hz;
    //............................................

    //............................................
    // Pre-initialize local ID
    for (n = 0; n < N; n++) {
        ID[n] = 1;
    }
    //............................................

    //............................................
    // .........Loop over the spheres.............
    for (p = 0; p < nspheres; p++) {
        // Get the sphere from the list, map to local min
        cx = List_cx[p] - min_x;
        cy = List_cy[p] - min_y;
        cz = List_cz[p] - min_z;
        r = List_rad[p];
        // Check if
        // Range for this sphere in global indexing
        imin = int((cx - r) / hx) - 1;
        imax = int((cx + r) / hx) + 1;
        jmin = int((cy - r) / hy) - 1;
        jmax = int((cy + r) / hy) + 1;
        kmin = int((cz - r) / hz) - 1;
        kmax = int((cz + r) / hz) + 1;
        // Obviously we have to do something at the edges
        if (imin < 0)
            imin = 0;
        if (imin > Nx)
            imin = Nx;
        if (imax < 0)
            imax = 0;
        if (imax > Nx)
            imax = Nx;
        if (jmin < 0)
            jmin = 0;
        if (jmin > Ny)
            jmin = Ny;
        if (jmax < 0)
            jmax = 0;
        if (jmax > Ny)
            jmax = Ny;
        if (kmin < 0)
            kmin = 0;
        if (kmin > Nz)
            kmin = Nz;
        if (kmax < 0)
            kmax = 0;
        if (kmax > Nz)
            kmax = Nz;
        // Loop over the domain for this sphere (may be null)
        for (i = imin; i < imax; i++) {
            for (j = jmin; j < jmax; j++) {
                for (k = kmin; k < kmax; k++) {
                    // Initialize ID value to 'fluid (=1)'
                    x = i * hx;
                    y = j * hy;
                    z = k * hz;
                    value = 1;
                    // if inside sphere, set to zero
                    if ((cx - x) * (cx - x) + (cy - y) * (cy - y) +
                            (cz - z) * (cz - z) <
                        r * r) {
                        value = 0;
                    }
                    // get the position in the list
                    n = k * Nx * Ny + j * Nx + i;
                    if (ID[n] != 0) {
                        ID[n] = value;
                    }
                }
            }
        }
    }
}

void SignedDistance(double *Distance, int nspheres, double *List_cx,
                    double *List_cy, double *List_cz, double *List_rad,
                    double Lx, double Ly, double Lz, int Nx, int Ny, int Nz,
                    int iproc, int jproc, int kproc, int nprocx, int nprocy,
                    int nprocz) {
    // Use sphere lists to determine which nodes are in porespace
    // Write out binary file for nodes
    int N = Nx * Ny * Nz; // Domain size, including the halo
    double hx, hy, hz;
    double x, y, z;
    double cx, cy, cz, r;
    int imin, imax, jmin, jmax, kmin, kmax;
    int p, i, j, k, n;
    //............................................
    double min_x, min_y, min_z;
    double distance;
    //............................................
    // Lattice spacing for the entire domain
    // It should generally be true that hx=hy=hz
    // Otherwise, you will end up with ellipsoids
    hx = Lx / ((Nx - 2) * nprocx - 1);
    hy = Ly / ((Ny - 2) * nprocy - 1);
    hz = Lz / ((Nz - 2) * nprocz - 1);
    //............................................
    // Get maximum and minimum for this domain
    // Halo is included !
    min_x = double(iproc * (Nx - 2) - 1) * hx;
    min_y = double(jproc * (Ny - 2) - 1) * hy;
    min_z = double(kproc * (Nz - 2) - 1) * hz;
    //............................................

    //............................................
    // Pre-initialize Distance
    for (n = 0; n < N; n++) {
        Distance[n] = 100.0;
    }
    //............................................

    //............................................
    // .........Loop over the spheres.............
    for (p = 0; p < nspheres; p++) {
        // Get the sphere from the list, map to local min
        cx = List_cx[p] - min_x;
        cy = List_cy[p] - min_y;
        cz = List_cz[p] - min_z;
        r = List_rad[p];
        // Check if
        // Range for this sphere in global indexing
        imin = int((cx - 2 * r) / hx);
        imax = int((cx + 2 * r) / hx) + 2;
        jmin = int((cy - 2 * r) / hy);
        jmax = int((cy + 2 * r) / hy) + 2;
        kmin = int((cz - 2 * r) / hz);
        kmax = int((cz + 2 * r) / hz) + 2;
        // Obviously we have to do something at the edges
        if (imin < 0)
            imin = 0;
        if (imin > Nx)
            imin = Nx;
        if (imax < 0)
            imax = 0;
        if (imax > Nx)
            imax = Nx;
        if (jmin < 0)
            jmin = 0;
        if (jmin > Ny)
            jmin = Ny;
        if (jmax < 0)
            jmax = 0;
        if (jmax > Ny)
            jmax = Ny;
        if (kmin < 0)
            kmin = 0;
        if (kmin > Nz)
            kmin = Nz;
        if (kmax < 0)
            kmax = 0;
        if (kmax > Nz)
            kmax = Nz;
        // Loop over the domain for this sphere (may be null)
        for (i = imin; i < imax; i++) {
            for (j = jmin; j < jmax; j++) {
                for (k = kmin; k < kmax; k++) {
                    // x,y,z is distance in physical units
                    x = i * hx;
                    y = j * hy;
                    z = k * hz;
                    // if inside sphere, set to zero
                    // get the position in the list
                    n = k * Nx * Ny + j * Nx + i;
                    // Compute the distance
                    distance = sqrt((cx - x) * (cx - x) + (cy - y) * (cy - y) +
                                    (cz - z) * (cz - z)) -
                               r;
                    // Assign the minimum distance
                    if (distance < Distance[n])
                        Distance[n] = distance;
                }
            }
        }
    }

    // Map the distance to lattice units
    for (n = 0; n < N; n++)
        Distance[n] = Distance[n] / hx;
}
