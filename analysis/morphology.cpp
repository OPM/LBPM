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
#include <analysis/morphology.h>
// Implementation of morphological opening routine

inline void PackID(const int *list, int count, signed char *sendbuf,
                   signed char *ID) {
    // Fill in the phase ID values from neighboring processors
    // This packs up the values that need to be sent from one processor to another
    int idx, n;

    for (idx = 0; idx < count; idx++) {
        n = list[idx];
        sendbuf[idx] = ID[n];
    }
}
//***************************************************************************************

inline void UnpackID(const int *list, int count, signed char *recvbuf,
                     signed char *ID) {
    // Fill in the phase ID values from neighboring processors
    // This unpacks the values once they have been recieved from neighbors
    int idx, n;

    for (idx = 0; idx < count; idx++) {
        n = list[idx];
        ID[n] = recvbuf[idx];
    }
}

Morphology::Morphology() {

    /* MPI tags*/
    sendtag = recvtag = 1381;
}

Morphology::~Morphology() {}

void Morphology::Initialize(std::shared_ptr<Domain> Dm, DoubleArray &Distance) {
    /* Loop over all faces and determine overlaps */
    size_t Nx = Dm->Nx;
    size_t Ny = Dm->Ny;
    size_t Nz = Dm->Nz;
    size_t N = Nx * Ny * Nz;

    int *tmpShift_x, *tmpShift_y, *tmpShift_z;
    double *tmpDistance;
    tmpShift_x = new int[N];
    tmpShift_y = new int[N];
    tmpShift_z = new int[N];
    tmpDistance = new double[N];

    double distance, boundary_distance;

    /* Loop over the local sub-domain and create overlap lists for each neighboring sub-domain */
    int sendLoc = 0; // counter for the local sub-domain send values
    int recvLoc = 0; // counter for the local recv
    //...................................................
    /* x face */
    sendCount = recvCount = 0;
    for (size_t k = 1; k < Nz - 1; k++) {
        for (size_t j = 1; j < Ny - 1; j++) {
            for (size_t i = 1; i < Nx - 1; i++) {
                distance = Distance(i, j, k);
                // Distance to x boundary
                boundary_distance = double(i - 1);
                if (distance > boundary_distance) {
                    tmpShift_x[sendCount] = Nx + i - 2;
                    tmpShift_y[sendCount] = j;
                    tmpShift_z[sendCount] = k;
                    tmpDistance[sendCount++] = distance;
                    int n = k * Nx * Ny + j * Nx + i;
                    sendID.push_back(n);
                }
            }
        }
    }
    Dm->Comm.Irecv(&recvCount, 1, Dm->rank_X(), recvtag + 0);
    Dm->Comm.send(&sendCount, 1, Dm->rank_x(), sendtag + 0);
    Dm->Comm.barrier();
    sendOffset_x = sendLoc;
    recvOffset_X = recvLoc;
    sendLoc += sendCount;
    recvLoc += recvCount;
    sendCount_x = sendCount;
    recvCount_X = recvCount;
    /* grow the arrays */
    xShift.resize(recvLoc);
    yShift.resize(recvLoc);
    zShift.resize(recvLoc);
    morphRadius.resize(recvLoc);
    //..............................
    /* send the morphological radius */
    Dm->Comm.Irecv(&morphRadius[recvOffset_X], recvCount, Dm->rank_X(),
                   recvtag + 0);
    Dm->Comm.send(&tmpDistance[0], sendCount, Dm->rank_x(), sendtag + 0);
    /* send the shift values */
    Dm->Comm.Irecv(&xShift[recvOffset_X], recvCount, Dm->rank_X(), recvtag + 1);
    Dm->Comm.send(&tmpShift_x[0], sendCount, Dm->rank_x(), sendtag + 1);
    Dm->Comm.Irecv(&yShift[recvOffset_X], recvCount, Dm->rank_X(), recvtag + 2);
    Dm->Comm.send(&tmpShift_y[0], sendCount, Dm->rank_x(), sendtag + 2);
    Dm->Comm.Irecv(&zShift[recvOffset_X], recvCount, Dm->rank_X(), recvtag + 3);
    Dm->Comm.send(&tmpShift_z[0], sendCount, Dm->rank_x(), sendtag + 3);
    Dm->Comm.barrier();
    //...................................................
    //...................................................
    /* X face */
    sendCount = recvCount = 0;
    for (size_t k = 1; k < Nz - 1; k++) {
        for (size_t j = 1; j < Ny - 1; j++) {
            for (size_t i = 1; i < Nx - 1; i++) {
                distance = Distance(i, j, k);
                // Distance to x boundary
                boundary_distance = double(Nx - i - 1);
                if (distance > boundary_distance) {
                    tmpShift_x[sendCount] = (i) - (Nx - 2);
                    tmpShift_y[sendCount] = j;
                    tmpShift_z[sendCount] = k;
                    tmpDistance[sendCount++] = distance;
                    int n = k * Nx * Ny + j * Nx + i;
                    sendID.push_back(n);
                }
            }
        }
    }
    Dm->Comm.Irecv(&recvCount, 1, Dm->rank_x(), recvtag + 0);
    Dm->Comm.send(&sendCount, 1, Dm->rank_X(), sendtag + 0);
    Dm->Comm.barrier();
    sendOffset_X = sendLoc;
    recvOffset_x = recvLoc;
    sendLoc += sendCount;
    recvLoc += recvCount;
    sendCount_X = sendCount;
    recvCount_x = recvCount;
    /* grow the arrays */
    xShift.resize(recvLoc);
    yShift.resize(recvLoc);
    zShift.resize(recvLoc);
    morphRadius.resize(recvLoc);
    //..............................
    /* send the morphological radius */
    Dm->Comm.Irecv(&morphRadius[recvOffset_x], recvCount, Dm->rank_x(),
                   recvtag + 0);
    Dm->Comm.send(&tmpDistance[0], sendCount, Dm->rank_X(), sendtag + 0);
    /* send the shift values */
    Dm->Comm.Irecv(&xShift[recvOffset_x], recvCount, Dm->rank_x(), recvtag + 1);
    Dm->Comm.send(&tmpShift_x[0], sendCount, Dm->rank_X(), sendtag + 1);
    Dm->Comm.Irecv(&yShift[recvOffset_x], recvCount, Dm->rank_x(), recvtag + 2);
    Dm->Comm.send(&tmpShift_y[0], sendCount, Dm->rank_X(), sendtag + 2);
    Dm->Comm.Irecv(&zShift[recvOffset_x], recvCount, Dm->rank_x(), recvtag + 3);
    Dm->Comm.send(&tmpShift_z[0], sendCount, Dm->rank_X(), sendtag + 3);
    Dm->Comm.barrier();
    //...................................................

    //...................................................
    /* y face */
    sendCount = recvCount = 0;
    for (size_t k = 1; k < Nz - 1; k++) {
        for (size_t j = 1; j < Ny - 1; j++) {
            for (size_t i = 1; i < Nx - 1; i++) {
                distance = Distance(i, j, k);
                // Distance to y boundary
                boundary_distance = double(j - 1);
                if (distance > boundary_distance) {
                    tmpShift_x[sendCount] = i;
                    tmpShift_y[sendCount] = Ny + j - 2;
                    tmpShift_z[sendCount] = k;
                    tmpDistance[sendCount++] = distance;
                    int n = k * Nx * Ny + j * Nx + i;
                    sendID.push_back(n);
                }
            }
        }
    }
    Dm->Comm.Irecv(&recvCount, 1, Dm->rank_Y(), recvtag + 0);
    Dm->Comm.send(&sendCount, 1, Dm->rank_y(), sendtag + 0);
    Dm->Comm.barrier();
    sendOffset_y = sendLoc;
    recvOffset_Y = recvLoc;
    sendLoc += sendCount;
    recvLoc += recvCount;
    sendCount_y = sendCount;
    recvCount_Y = recvCount;
    /* grow the arrays */
    xShift.resize(recvLoc);
    yShift.resize(recvLoc);
    zShift.resize(recvLoc);
    morphRadius.resize(recvLoc);
    //..............................
    /* send the morphological radius */
    Dm->Comm.Irecv(&morphRadius[recvOffset_Y], recvCount, Dm->rank_Y(),
                   recvtag + 0);
    Dm->Comm.send(&tmpDistance[0], sendCount, Dm->rank_y(), sendtag + 0);
    /* send the shift values */
    Dm->Comm.Irecv(&xShift[recvOffset_Y], recvCount, Dm->rank_Y(), recvtag + 1);
    Dm->Comm.send(&tmpShift_x[0], sendCount, Dm->rank_y(), sendtag + 1);
    Dm->Comm.Irecv(&yShift[recvOffset_Y], recvCount, Dm->rank_Y(), recvtag + 2);
    Dm->Comm.send(&tmpShift_y[0], sendCount, Dm->rank_y(), sendtag + 2);
    Dm->Comm.Irecv(&zShift[recvOffset_Y], recvCount, Dm->rank_Y(), recvtag + 3);
    Dm->Comm.send(&tmpShift_z[0], sendCount, Dm->rank_y(), sendtag + 3);
    Dm->Comm.barrier();
    //...................................................
    //...................................................
    /* X face */
    sendCount = recvCount = 0;
    for (size_t k = 1; k < Nz - 1; k++) {
        for (size_t j = 1; j < Ny - 1; j++) {
            for (size_t i = 1; i < Nx - 1; i++) {
                distance = Distance(i, j, k);
                // Distance to x boundary
                boundary_distance = double(Ny - j - 1);
                if (distance > boundary_distance) {
                    tmpShift_x[sendCount] = i;
                    tmpShift_y[sendCount] = j - (Ny - 2);
                    tmpShift_z[sendCount] = k;
                    tmpDistance[sendCount++] = distance;
                    int n = k * Nx * Ny + j * Nx + i;
                    sendID.push_back(n);
                }
            }
        }
    }
    Dm->Comm.Irecv(&recvCount, 1, Dm->rank_y(), recvtag + 0);
    Dm->Comm.send(&sendCount, 1, Dm->rank_Y(), sendtag + 0);
    Dm->Comm.barrier();
    sendOffset_Y = sendLoc;
    recvOffset_y = recvLoc;
    sendLoc += sendCount;
    recvLoc += recvCount;
    sendCount_Y = sendCount;
    recvCount_y = recvCount;
    /* grow the arrays */
    xShift.resize(recvLoc);
    yShift.resize(recvLoc);
    zShift.resize(recvLoc);
    morphRadius.resize(recvLoc);
    //..............................
    /* send the morphological radius */
    Dm->Comm.Irecv(&morphRadius[recvOffset_y], recvCount, Dm->rank_y(),
                   recvtag + 0);
    Dm->Comm.send(&tmpDistance[0], sendCount, Dm->rank_Y(), sendtag + 0);
    /* send the shift values */
    Dm->Comm.Irecv(&xShift[recvOffset_y], recvCount, Dm->rank_y(), recvtag + 1);
    Dm->Comm.send(&tmpShift_x[0], sendCount, Dm->rank_Y(), sendtag + 1);
    Dm->Comm.Irecv(&yShift[recvOffset_y], recvCount, Dm->rank_y(), recvtag + 2);
    Dm->Comm.send(&tmpShift_y[0], sendCount, Dm->rank_Y(), sendtag + 2);
    Dm->Comm.Irecv(&zShift[recvOffset_y], recvCount, Dm->rank_y(), recvtag + 3);
    Dm->Comm.send(&tmpShift_z[0], sendCount, Dm->rank_Y(), sendtag + 3);
    Dm->Comm.barrier();
    //...................................................

    //...................................................
    /* z face */
    sendCount = recvCount = 0;
    for (size_t k = 1; k < Nz - 1; k++) {
        for (size_t j = 1; j < Ny - 1; j++) {
            for (size_t i = 1; i < Nx - 1; i++) {
                distance = Distance(i, j, k);
                // Distance to z boundary
                boundary_distance = double(k - 1);
                if (distance > boundary_distance) {
                    tmpShift_x[sendCount] = i;
                    tmpShift_y[sendCount] = j;
                    tmpShift_z[sendCount] = (Nz - 2) + k;
                    tmpDistance[sendCount++] = distance;
                    int n = k * Nx * Ny + j * Nx + i;
                    sendID.push_back(n);
                }
            }
        }
    }
    Dm->Comm.Irecv(&recvCount, 1, Dm->rank_Z(), recvtag + 0);
    Dm->Comm.send(&sendCount, 1, Dm->rank_z(), sendtag + 0);
    Dm->Comm.barrier();
    sendOffset_z = sendLoc;
    recvOffset_Z = recvLoc;
    sendLoc += sendCount;
    recvLoc += recvCount;
    sendCount_z = sendCount;
    recvCount_Z = recvCount;
    /* grow the arrays */
    xShift.resize(recvLoc);
    yShift.resize(recvLoc);
    zShift.resize(recvLoc);
    morphRadius.resize(recvLoc);
    //..............................
    /* send the morphological radius */
    Dm->Comm.Irecv(&morphRadius[recvOffset_Z], recvCount, Dm->rank_Z(),
                   recvtag + 0);
    Dm->Comm.send(&tmpDistance[0], sendCount, Dm->rank_z(), sendtag + 0);
    /* send the shift values */
    Dm->Comm.Irecv(&xShift[recvOffset_Z], recvCount, Dm->rank_Z(), recvtag + 1);
    Dm->Comm.send(&tmpShift_x[0], sendCount, Dm->rank_z(), sendtag + 1);
    Dm->Comm.Irecv(&yShift[recvOffset_Z], recvCount, Dm->rank_Z(), recvtag + 2);
    Dm->Comm.send(&tmpShift_y[0], sendCount, Dm->rank_z(), sendtag + 2);
    Dm->Comm.Irecv(&zShift[recvOffset_Z], recvCount, Dm->rank_Z(), recvtag + 3);
    Dm->Comm.send(&tmpShift_z[0], sendCount, Dm->rank_z(), sendtag + 3);
    Dm->Comm.barrier();
    //...................................................
    /* Z face */
    sendCount = recvCount = 0;
    for (size_t k = 1; k < Nz - 1; k++) {
        for (size_t j = 1; j < Ny - 1; j++) {
            for (size_t i = 1; i < Nx - 1; i++) {
                distance = Distance(i, j, k);
                // Distance to x boundary
                boundary_distance = double(Nz - k - 1);
                if (distance > boundary_distance) {
                    tmpShift_x[sendCount] = i;
                    tmpShift_y[sendCount] = j;
                    tmpShift_z[sendCount] = k - (Nz - 2);
                    tmpDistance[sendCount++] = distance;
                    int n = k * Nx * Ny + j * Nx + i;
                    sendID.push_back(n);
                }
            }
        }
    }
    Dm->Comm.Irecv(&recvCount, 1, Dm->rank_z(), recvtag + 0);
    Dm->Comm.send(&sendCount, 1, Dm->rank_Z(), sendtag + 0);
    Dm->Comm.barrier();
    sendOffset_Z = sendLoc;
    recvOffset_z = recvLoc;
    sendLoc += sendCount;
    recvLoc += recvCount;
    sendCount_Z = sendCount;
    recvCount_z = recvCount;
    /* grow the arrays */
    xShift.resize(recvLoc);
    yShift.resize(recvLoc);
    zShift.resize(recvLoc);
    morphRadius.resize(recvLoc);
    //..............................
    /* send the morphological radius */
    Dm->Comm.Irecv(&morphRadius[recvOffset_z], recvCount, Dm->rank_z(),
                   recvtag + 0);
    Dm->Comm.send(&tmpDistance[0], sendCount, Dm->rank_Z(), sendtag + 0);
    /* send the shift values */
    Dm->Comm.Irecv(&xShift[recvOffset_z], recvCount, Dm->rank_z(), recvtag + 1);
    Dm->Comm.send(&tmpShift_x[0], sendCount, Dm->rank_Z(), sendtag + 1);
    Dm->Comm.Irecv(&yShift[recvOffset_z], recvCount, Dm->rank_z(), recvtag + 2);
    Dm->Comm.send(&tmpShift_y[0], sendCount, Dm->rank_Z(), sendtag + 2);
    Dm->Comm.Irecv(&zShift[recvOffset_z], recvCount, Dm->rank_z(), recvtag + 3);
    Dm->Comm.send(&tmpShift_z[0], sendCount, Dm->rank_Z(), sendtag + 3);
    Dm->Comm.barrier();
    //...................................................

    /* resize the send / recv lists */
    sendCount = sendLoc;
    recvCount = recvLoc;
    sendList.resize(sendLoc);
    recvList.resize(recvLoc);
    localID.resize(sendCount);
    nonlocalID.resize(recvCount);

    /*printf("  offset %i for send (x) %i \n", sendOffset_x, sendCount_x);
	printf("  offset %i for send (X) %i \n", sendOffset_X, sendCount_X);
	printf("  offset %i for send (y) %i \n", sendOffset_y, sendCount_y);
	printf("  offset %i for send (Y) %i \n", sendOffset_Y, sendCount_Y);
	printf("  offset %i for send (z) %i \n", sendOffset_z, sendCount_z);
	printf("  offset %i for send (Z) %i \n", sendOffset_Z, sendCount_Z);
	*/
}

int Morphology::GetOverlaps(std::shared_ptr<Domain> Dm, signed char *id,
                            const signed char ErodeLabel,
                            const signed char NewLabel) {

    int Nx = Dm->Nx;
    int Ny = Dm->Ny;
    int Nz = Dm->Nz;
    int LocalNumber = 0;
    int i, j, k, ii, jj, kk;
    int imin, jmin, kmin, imax, jmax, kmax;

    for (int idx = 0; idx < sendCount; idx++) {
        int n = sendID[idx];
        localID[idx] = id[n];
    }
    //printf("send x -- offset: %i, count: %i \n",sendOffset_x,sendCount_x);
    Dm->Comm.Irecv(&nonlocalID[recvOffset_X], recvCount_X, Dm->rank_x(),
                   recvtag + 2);
    Dm->Comm.send(&localID[sendOffset_x], sendCount_x, Dm->rank_X(),
                  sendtag + 2);

    //printf("send X \n");
    Dm->Comm.Irecv(&nonlocalID[recvOffset_x], recvCount_x, Dm->rank_X(),
                   recvtag + 3);
    Dm->Comm.send(&localID[sendOffset_X], sendCount_X, Dm->rank_x(),
                  sendtag + 3);

    //printf("send y \n");
    Dm->Comm.Irecv(&nonlocalID[recvOffset_Y], recvCount_Y, Dm->rank_y(),
                   recvtag + 4);
    Dm->Comm.send(&localID[sendOffset_y], sendCount_y, Dm->rank_Y(),
                  sendtag + 4);

    //printf("send Y \n");
    Dm->Comm.Irecv(&nonlocalID[recvOffset_y], recvCount_y, Dm->rank_Y(),
                   recvtag + 5);
    Dm->Comm.send(&localID[sendOffset_Y], sendCount_Y, Dm->rank_y(),
                  sendtag + 5);

    //printf("send z \n");
    Dm->Comm.Irecv(&nonlocalID[recvOffset_Z], recvCount_Z, Dm->rank_z(),
                   recvtag + 6);
    Dm->Comm.send(&localID[sendOffset_z], sendCount_z, Dm->rank_Z(),
                  sendtag + 6);

    //printf("send Z \n");
    Dm->Comm.Irecv(&nonlocalID[recvOffset_z], recvCount_z, Dm->rank_Z(),
                   recvtag + 7);
    Dm->Comm.send(&localID[sendOffset_Z], sendCount_Z, Dm->rank_z(),
                  sendtag + 7);

    for (int idx = 0; idx < recvCount; idx++) {
        double radius = morphRadius[idx];
        signed char label = nonlocalID[idx];
        /* get the neighboring site index */
        i = xShift[idx];
        j = yShift[idx];
        k = zShift[idx];
        int Window = int(radius);
        // loop over the window and update
        if (label == NewLabel) {
            imin = max(1, i - Window);
            jmin = max(1, j - Window);
            kmin = max(1, k - Window);
            imax = min(Nx - 1, i + Window);
            jmax = min(Ny - 1, j + Window);
            kmax = min(Nz - 1, k + Window);
            for (kk = kmin; kk < kmax; kk++) {
                for (jj = jmin; jj < jmax; jj++) {
                    for (ii = imin; ii < imax; ii++) {
                        int nn = kk * Nx * Ny + jj * Nx + ii;
                        double dsq =
                            double((ii - i) * (ii - i) + (jj - j) * (jj - j) +
                                   (kk - k) * (kk - k));
                        if (id[nn] == ErodeLabel && dsq <= radius * radius) {
                            LocalNumber += 1.0;
                            id[nn] = NewLabel;
                        }
                    }
                }
            }
        }
    }
    Dm->Comm.barrier();

    return LocalNumber;
}

//***************************************************************************************
double MorphOpen(DoubleArray &SignDist, signed char *id,
                 std::shared_ptr<Domain> Dm, double VoidFraction,
                 signed char ErodeLabel, signed char NewLabel) {
    // SignDist is the distance to the object that you want to constaing the morphological opening
    // VoidFraction is the the empty space where the object inst
    // id is a labeled map
    // Dm contains information about the domain structure

    int nx = Dm->Nx;
    int ny = Dm->Ny;
    int nz = Dm->Nz;
    int nprocx = Dm->nprocx();
    int nprocy = Dm->nprocy();
    int nprocz = Dm->nprocz();
    int rank = Dm->rank();

    int n;
    double final_void_fraction;
    double count, countGlobal, totalGlobal;
    count = 0.f;
    double maxdist = -200.f;
    double maxdistGlobal;
    for (int k = 1; k < nz - 1; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                n = k * nx * ny + j * nx + i;
                // extract maximum distance for critical radius
                if (SignDist(i, j, k) > maxdist)
                    maxdist = SignDist(i, j, k);
                if (id[n] == ErodeLabel) {
                    count += 1.0;
                    //id[n]  = 2;
                }
            }
        }
    }
    Dm->Comm.barrier();

    Morphology Structure;
    Structure.Initialize(Dm, SignDist);

    // total Global is the number of nodes in the pore-space
    totalGlobal = Dm->Comm.sumReduce(count);
    maxdistGlobal = Dm->Comm.sumReduce(maxdist);
    double volume = double(nprocx * nprocy * nprocz) * double(nx - 2) *
                    double(ny - 2) * double(nz - 2);
    double volume_fraction = totalGlobal / volume;
    if (rank == 0)
        printf("Volume fraction for morphological opening: %f \n",
               volume_fraction);
    if (rank == 0)
        printf("Maximum pore size: %f \n", maxdistGlobal);
    final_void_fraction = volume_fraction; //initialize

    int ii, jj, kk;
    int imin, jmin, kmin, imax, jmax, kmax;
    int Nx = nx;
    int Ny = ny;
    int Nz = nz;

    double void_fraction_old = 1.0;
    double void_fraction_new = 1.0;
    double void_fraction_diff_old = 1.0;
    double void_fraction_diff_new = 1.0;
    if (ErodeLabel == 1) {
        VoidFraction = 1.0 - VoidFraction;
    }

    // Increase the critical radius until the target saturation is met
    double deltaR = 0.05; // amount to change the radius in voxel units
    double Rcrit_old = maxdistGlobal;
    double Rcrit_new = maxdistGlobal;

    int numTry = 0;
    int maxTry = 100;
    while (!(void_fraction_new < VoidFraction) && numTry < maxTry) {
        numTry++;
        void_fraction_diff_old = void_fraction_diff_new;
        void_fraction_old = void_fraction_new;
        Rcrit_old = Rcrit_new;
        Rcrit_new -= deltaR * Rcrit_old;
        if (rank == 0)
            printf("Try %i with radius %f \n", numTry, Rcrit_new);
        if (Rcrit_new < 0.5) {
            numTry = maxTry;
        }
        int Window = round(Rcrit_new);
        if (Window == 0)
            Window =
                1; // If Window = 0 at the begining, after the following process will have sw=1.0
        // and sw<Sw will be immediately broken
        double LocalNumber = 0.f;
        for (int k = 1; k < Nz - 1; k++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int i = 1; i < Nx - 1; i++) {
                    n = k * nx * ny + j * nx + i;
                    if (SignDist(i, j, k) > Rcrit_new) {
                        // loop over the window and update
                        //printf("Distance(%i %i %i) = %f \n",i,j,k, SignDist(i,j,k));
                        imin = max(1, i - Window);
                        jmin = max(1, j - Window);
                        kmin = max(1, k - Window);
                        imax = min(Nx - 1, i + Window);
                        jmax = min(Ny - 1, j + Window);
                        kmax = min(Nz - 1, k + Window);
                        for (kk = kmin; kk < kmax; kk++) {
                            for (jj = jmin; jj < jmax; jj++) {
                                for (ii = imin; ii < imax; ii++) {
                                    int nn = kk * nx * ny + jj * nx + ii;
                                    double dsq = double((ii - i) * (ii - i) +
                                                        (jj - j) * (jj - j) +
                                                        (kk - k) * (kk - k));
                                    if (id[nn] == ErodeLabel &&
                                        dsq <= Rcrit_new * Rcrit_new) {
                                        LocalNumber += 1.0;
                                        id[nn] = NewLabel;
                                    }
                                }
                            }
                        }
                    }
                    // move on
                }
            }
        }
        //LocalNumber += Structure.GetOverlaps(Dm, id, ErodeLabel, NewLabel);

        count = 0.f;
        for (int k = 1; k < Nz - 1; k++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int i = 1; i < Nx - 1; i++) {
                    n = k * Nx * Ny + j * Nx + i;
                    if (id[n] == ErodeLabel) {
                        count += 1.0;
                    }
                }
            }
        }
        countGlobal = Dm->Comm.sumReduce(count);
        void_fraction_new = countGlobal / totalGlobal;
        void_fraction_diff_new = abs(void_fraction_new - VoidFraction);
        if (rank == 0) {
            printf("     %f ", void_fraction_new);
            printf("     %f\n", Rcrit_new);
        }
    }

    if (void_fraction_diff_new < void_fraction_diff_old) {
        final_void_fraction = void_fraction_new;
        if (rank == 0) {
            printf("Final void fraction =%f\n", void_fraction_new);
            printf("Final critical radius=%f\n", Rcrit_new);
        }
    } else {
        final_void_fraction = void_fraction_old;
        if (rank == 0) {
            printf("Final void fraction=%f\n", void_fraction_old);
            printf("Final critical radius=%f\n", Rcrit_old);
        }
    }
    return final_void_fraction;
}

//***************************************************************************************
double MorphDrain(DoubleArray &SignDist, signed char *id,
                  std::shared_ptr<Domain> Dm, double VoidFraction,
                  double InitialRadius) {
    // SignDist is the distance to the object that you want to constaing the morphological opening
    // VoidFraction is the the empty space where the object inst
    // id is a labeled map
    // Dm contains information about the domain structure

    signed char ErodeLabel = 2;
    signed char NewLabel = 1;

    int nx = Dm->Nx;
    int ny = Dm->Ny;
    int nz = Dm->Nz;
    int nprocx = Dm->nprocx();
    int nprocy = Dm->nprocy();
    int nprocz = Dm->nprocz();
    int rank = Dm->rank();

    DoubleArray phase(nx, ny, nz);
    IntArray phase_label(nx, ny, nz);
    Array<char> ID(nx, ny, nz);
    fillHalo<char> fillChar(Dm->Comm, Dm->rank_info, {nx - 2, ny - 2, nz - 2},
                            {1, 1, 1}, 0, 1);

    Morphology Structure;
    Structure.Initialize(Dm, SignDist);

    int n;
    double final_void_fraction;
    double count, countGlobal, totalGlobal;
    count = 0.f;
    double maxdist = -200.f;
    double maxdistGlobal;
    for (int k = 1; k < nz - 1; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                n = k * nx * ny + j * nx + i;
                // extract maximum distance for critical radius
                if (SignDist(i, j, k) > maxdist)
                    maxdist = SignDist(i, j, k);
                if (SignDist(i, j, k) > 0.0) {
                    count += 1.0;
                    id[n] = ErodeLabel;
                }
                ID(i, j, k) = id[n];
            }
        }
    }
    fillChar.fill(ID);
    Dm->Comm.barrier();

    // total Global is the number of nodes in the pore-space
    totalGlobal = Dm->Comm.sumReduce(count);
    maxdistGlobal = Dm->Comm.sumReduce(maxdist);
    double volume = double(nprocx * nprocy * nprocz) * double(nx - 2) *
                    double(ny - 2) * double(nz - 2);
    double volume_fraction = totalGlobal / volume;
    if (rank == 0)
        printf("Volume fraction for morphological opening: %f \n",
               volume_fraction);
    if (rank == 0)
        printf("Maximum pore size: %f \n", maxdistGlobal);

    int ii, jj, kk;
    int imin, jmin, kmin, imax, jmax, kmax;
    int Nx = nx;
    int Ny = ny;
    int Nz = nz;

    double void_fraction_old = 1.0;
    double void_fraction_new = 1.0;
    double void_fraction_diff_old = 1.0;
    double void_fraction_diff_new = 1.0;

    // Increase the critical radius until the target saturation is met
    double deltaR = 0.05; // amount to change the radius in voxel units
    double Rcrit_old = maxdistGlobal;
    double Rcrit_new = maxdistGlobal;

    if (InitialRadius < maxdistGlobal) {
        Rcrit_old = InitialRadius;
        Rcrit_new = InitialRadius;
    }
    //if (argc>2){
    //	Rcrit_new = strtod(argv[2],NULL);
    //	if (rank==0) printf("Max. distance =%f, Initial critical radius = %f \n",maxdistGlobal,Rcrit_new);
    //}
    Dm->Comm.barrier();

    FILE *DRAIN = fopen("morphdrain.csv", "w");
    fprintf(DRAIN, "sw radius\n");

    while (void_fraction_new > VoidFraction && Rcrit_new > 0.5) {
        void_fraction_diff_old = void_fraction_diff_new;
        void_fraction_old = void_fraction_new;
        Rcrit_old = Rcrit_new;
        Rcrit_new -= deltaR * Rcrit_old;
        int Window = round(Rcrit_new);
        if (Window == 0)
            Window =
                1; // If Window = 0 at the begining, after the following process will have sw=1.0
        // and sw<Sw will be immediately broken
        double LocalNumber = 0.f;
        for (int k = 1; k < Nz - 1; k++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int i = 1; i < Nx - 1; i++) {
                    n = k * nx * ny + j * nx + i;
                    if (SignDist(i, j, k) > Rcrit_new) {
                        // loop over the window and update
                        imin = max(1, i - Window);
                        jmin = max(1, j - Window);
                        kmin = max(1, k - Window);
                        imax = min(Nx - 1, i + Window);
                        jmax = min(Ny - 1, j + Window);
                        kmax = min(Nz - 1, k + Window);
                        for (kk = kmin; kk < kmax; kk++) {
                            for (jj = jmin; jj < jmax; jj++) {
                                for (ii = imin; ii < imax; ii++) {
                                    double dsq = double((ii - i) * (ii - i) +
                                                        (jj - j) * (jj - j) +
                                                        (kk - k) * (kk - k));
                                    if (ID(ii, jj, kk) == ErodeLabel &&
                                        dsq <=
                                            (Rcrit_new + 1) * (Rcrit_new + 1)) {
                                        LocalNumber += 1.0;
                                        //id[nn]=1;
                                        ID(ii, jj, kk) = NewLabel;
                                        id[kk * Nx * Ny + jj * Nx + ii] =
                                            NewLabel;
                                    }
                                }
                            }
                        }
                    }
                    // move on
                }
            }
        }
        LocalNumber += Structure.GetOverlaps(Dm, id, ErodeLabel, NewLabel);

        for (int k = 1; k < Nz - 1; k++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int i = 1; i < Nx - 1; i++) {
                    ID(i, j, k) = id[k * Nx * Ny + j * Nx + i];
                }
            }
        }
        fillChar.fill(ID);
        Dm->Comm.barrier();

        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    if (ID(i, j, k) == NewLabel) {
                        phase(i, j, k) = 1.0;
                    } else
                        phase(i, j, k) = -1.0;
                }
            }
        }

        // Extract only the connected part of NWP
        double vF = 0.0;
        double vS = 0.0;
        ComputeGlobalBlobIDs(nx - 2, ny - 2, nz - 2, Dm->rank_info, phase,
                             SignDist, vF, vS, phase_label, Dm->Comm);
        Dm->Comm.barrier();

        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    n = k * nx * ny + j * nx + i;
                    if (ID(i, j, k) == 1 && phase_label(i, j, k) > 1) {
                        ID(i, j, k) = ErodeLabel;
                    }
                    id[n] = ID(i, j, k);
                }
            }
        }

        count = 0.f;
        for (int k = 1; k < nz - 1; k++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    n = k * nx * ny + j * nx + i;
                    if (id[n] > 1) {
                        count += 1.0;
                    }
                }
            }
        }
        countGlobal = Dm->Comm.sumReduce(count);
        void_fraction_new = countGlobal / totalGlobal;
        void_fraction_diff_new = abs(void_fraction_new - VoidFraction);
        if (rank == 0) {
            fprintf(DRAIN, "%f ", void_fraction_new);
            fprintf(DRAIN, "%f\n", Rcrit_new);
            printf("     %f ", void_fraction_new);
            printf("     %f\n", Rcrit_new);
        }
    }

    if (void_fraction_diff_new < void_fraction_diff_old) {
        final_void_fraction = void_fraction_new;
        if (rank == 0) {
            printf("Final void fraction =%f\n", void_fraction_new);
            printf("Final critical radius=%f\n", Rcrit_new);
        }
    } else {
        final_void_fraction = void_fraction_old;
        if (rank == 0) {
            printf("Final void fraction=%f\n", void_fraction_old);
            printf("Final critical radius=%f\n", Rcrit_old);
        }
    }
    // label all WP components as 2
    for (int k = 1; k < nz - 1; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                n = k * nx * ny + j * nx + i;
                if (id[n] > 1) {
                    id[n] = 2;
                }
            }
        }
    }

    return final_void_fraction;
}

//double MorphGrow(DoubleArray &BoundaryDist, DoubleArray &Dist, Array<char> &id, std::shared_ptr<Domain> Dm, double TargetGrowth)
double MorphGrow(DoubleArray &BoundaryDist, DoubleArray &Dist, Array<char> &id,
                 std::shared_ptr<Domain> Dm, double TargetGrowth,
                 double WallFactor) {
    int Nx = Dm->Nx;
    int Ny = Dm->Ny;
    int Nz = Dm->Nz;
    int rank = Dm->rank();

    double count = 0.0;
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                if (Dist(i, j, k) < 0.0) {
                    count += 1.0;
                }
            }
        }
    }
    double count_original = Dm->Comm.sumReduce(count);

    // Estimate morph_delta
    double morph_delta = 0.0;
    if (TargetGrowth > 0.0)
        morph_delta = 0.1;
    else
        morph_delta = -0.1;
    double morph_delta_previous = 0.0;
    double GrowthEstimate = 0.0;
    double GrowthPrevious = 0.0;
    int COUNT_FOR_LOOP = 0;
    double ERROR = 100.0;
    if (rank == 0)
        printf("Estimate delta for growth=%f \n", TargetGrowth);
    while (ERROR > 0.01 && COUNT_FOR_LOOP < 10) {
        COUNT_FOR_LOOP++;
        count = 0.0;
        double MAX_DISPLACEMENT = 0.0;
        for (int k = 1; k < Nz - 1; k++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int i = 1; i < Nx - 1; i++) {
                    double walldist = BoundaryDist(i, j, k);
                    //double wallweight = 1.0 / (1+exp(-5.f*(walldist-1.f)));
                    double wallweight =
                        WallFactor / (1 + exp(-5.f * (walldist - 1.f)));
                    //wallweight = 1.0;
                    if (fabs(wallweight * morph_delta) > MAX_DISPLACEMENT)
                        MAX_DISPLACEMENT = fabs(wallweight * morph_delta);

                    if (Dist(i, j, k) - wallweight * morph_delta < 0.0) {
                        count += 1.0;
                    }
                }
            }
        }
        count = Dm->Comm.sumReduce(count);
        MAX_DISPLACEMENT = Dm->Comm.maxReduce(MAX_DISPLACEMENT);
        GrowthEstimate = count - count_original;
        ERROR = fabs((GrowthEstimate - TargetGrowth) / TargetGrowth);

        if (rank == 0)
            printf("     delta=%f, growth=%f, max. displacement = %f \n",
                   morph_delta, GrowthEstimate, MAX_DISPLACEMENT);
        // Now adjust morph_delta
        if (fabs(GrowthEstimate - GrowthPrevious) > 0.0) {
            double step_size = (TargetGrowth - GrowthEstimate) *
                               (morph_delta - morph_delta_previous) /
                               (GrowthEstimate - GrowthPrevious);
            GrowthPrevious = GrowthEstimate;
            morph_delta_previous = morph_delta;
            morph_delta += step_size;
        }
        if (morph_delta / morph_delta_previous > 2.0)
            morph_delta = morph_delta_previous * 2.0;

        //MAX_DISPLACEMENT *= max(TargetGrowth/GrowthEstimate,1.25);

        if (morph_delta > 0.0) {
            // object is growing
            if (MAX_DISPLACEMENT > 3.0) {
                morph_delta = 3.0;
                COUNT_FOR_LOOP = 100; // exit loop if displacement is too large
            }
        } else {
            // object is shrinking
            if (MAX_DISPLACEMENT > 1.0) {
                morph_delta = -1.0;
                COUNT_FOR_LOOP = 100; // exit loop if displacement is too large
            }
        }
    }
    if (rank == 0)
        printf("Final delta=%f \n", morph_delta);

    count = 0.0;
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                double walldist = BoundaryDist(i, j, k);
                //double wallweight = 1.0 / (1+exp(-5.f*(walldist-1.f)));
                //wallweight = 1.0;
                double wallweight =
                    WallFactor / (1 + exp(-5.f * (walldist - 1.f)));
                Dist(i, j, k) -= wallweight * morph_delta;

                if (Dist(i, j, k) < 0.0)
                    count += 1.0;
            }
        }
    }
    count = Dm->Comm.sumReduce(count);

    return count;
}
