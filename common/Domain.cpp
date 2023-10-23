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
// Created by James McClure
// Copyright 2008-2020
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <exception> // std::exception
#include <stdexcept>

#include "common/Domain.h"
#include "common/Array.h"
#include "common/Utilities.h"
#include "common/MPI.h"
#include "common/Communication.h"

// Inline function to read line without a return argument
static inline void fgetl(char *str, int num, FILE *stream) {
    char *ptr = fgets(str, num, stream);
    if (0) {
        char *temp = (char *)&ptr;
        temp++;
    }
}

void Domain::read_swc(const std::string &Filename) {
    //...... READ IN SWC FILE...................................
    int count = 0;
    int number_of_lines = 0;
    if (rank() == 0) {
        cout << "Reading SWC file..." << endl;
        {
            std::string line;
            std::ifstream myfile(Filename);
            while (std::getline(myfile, line))
                ++number_of_lines;
            number_of_lines -= 1;
        }
        std::cout << "    Number of lines in SWC file: " << number_of_lines
                  << endl;
    }
    count = Comm.sumReduce(number_of_lines); // nonzero only for rank=0
    number_of_lines = count;

    // set up structures to read
    double *List_cx = new double[number_of_lines];
    double *List_cy = new double[number_of_lines];
    double *List_cz = new double[number_of_lines];
    double *List_rad = new double[number_of_lines];
    int *List_index = new int[number_of_lines];
    int *List_parent = new int[number_of_lines];
    int *List_type = new int[number_of_lines];

    if (rank() == 0) {
        FILE *fid = fopen(Filename.c_str(), "rb");
        INSIST(fid != NULL, "Error opening SWC file");
        //.........Trash the header lines (x 1)..........
        char line[100];
        fgetl(line, 100, fid);
        //........read the spheres..................
        // We will read until a blank like or end-of-file is reached
        count = 0;
        while (!feof(fid) && fgets(line, 100, fid) != NULL) {
            char *line2 = line;
            List_index[count] = int(strtod(line2, &line2));
            List_type[count] = int(strtod(line2, &line2));
            List_cx[count] = strtod(line2, &line2);
            List_cy[count] = strtod(line2, &line2);
            List_cz[count] = strtod(line2, &line2);
            List_rad[count] = strtod(line2, &line2);
            List_parent[count] = int(strtod(line2, &line2));
            count++;
        }
        fclose(fid);
        cout << "   Number of lines extracted is: " << count << endl;
        INSIST(count == number_of_lines, "Problem reading swc file!");

        double min_cx = List_cx[0] - List_rad[0];
        double min_cy = List_cy[0] - List_rad[0];
        double min_cz = List_cz[0] - List_rad[0];
        for (count = 1; count < number_of_lines; count++) {
            double value_x = List_cx[count] - List_rad[count];
            double value_y = List_cy[count] - List_rad[count];
            double value_z = List_cz[count] - List_rad[count];
            if (value_x < min_cx)
                min_cx = value_x;
            if (value_y < min_cy)
                min_cy = value_y;
            if (value_z < min_cz)
                min_cz = value_z;
        }
        /* shift the swc data */
        printf("   shift swc data by %f, %f, %f \n", min_cx, min_cy, min_cz);
        for (count = 0; count < number_of_lines; count++) {
            List_cx[count] -= offset_x * voxel_length;
            List_cy[count] -= offset_y * voxel_length;
            List_cz[count] -= offset_z * voxel_length;
        }
    }
    /* everybody gets the swc file */
    Comm.bcast(List_cx, number_of_lines, 0);
    Comm.bcast(List_cy, number_of_lines, 0);
    Comm.bcast(List_cz, number_of_lines, 0);
    Comm.bcast(List_rad, number_of_lines, 0);
    Comm.bcast(List_index, number_of_lines, 0);
    Comm.bcast(List_parent, number_of_lines, 0);
    Comm.bcast(List_type, number_of_lines, 0);

    /* units of swc file are in micron */
    double start_x, start_y, start_z;
    /* box owned by this rank */
    start_x = rank_info.ix * (Nx - 2) * voxel_length;
    start_y = rank_info.jy * (Ny - 2) * voxel_length;
    start_z = rank_info.kz * (Nz - 2) * voxel_length;
    //finish_x = (rank_info.ix+1)*(Nx-2)*voxel_length;
    //finish_y = (rank_info.jy+1)*(Ny-2)*voxel_length;
    //finish_z = (rank_info.kz+1)*(Nz-2)*voxel_length;

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                id[k * Nx * Ny + j * Nx + i] = 1;
            }
        }
    }

    /* Loop over SWC input and populate domain ID */
    for (int idx = 0; idx < number_of_lines; idx++) {
        /* get the object information */
        int parent = List_parent[idx] - 1;
        if (parent < 0)
            parent = idx;
        double xi = List_cx[idx];
        double yi = List_cy[idx];
        double zi = List_cz[idx];
        double xp = List_cx[parent];
        double yp = List_cy[parent];
        double zp = List_cz[parent];
        double ri = List_rad[idx];
        double rp = List_rad[parent];

        int radius_in_voxels = int(List_rad[idx] / voxel_length);
        signed char label = char(List_type[idx]);

        double xmin = min(((xi - start_x - List_rad[idx]) / voxel_length),
                          ((xp - start_x - List_rad[parent]) / voxel_length));
        double ymin = min(((yi - start_y - List_rad[idx]) / voxel_length),
                          ((yp - start_y - List_rad[parent]) / voxel_length));
        double zmin = min(((zi - start_z - List_rad[idx]) / voxel_length),
                          ((zp - start_z - List_rad[parent]) / voxel_length));
        double xmax = max(((xi - start_x + List_rad[idx]) / voxel_length),
                          ((xp - start_x + List_rad[parent]) / voxel_length));
        double ymax = max(((yi - start_y + List_rad[idx]) / voxel_length),
                          ((yp - start_y + List_rad[parent]) / voxel_length));
        double zmax = max(((zi - start_z + List_rad[idx]) / voxel_length),
                          ((zp - start_z + List_rad[parent]) / voxel_length));

        /*	if (rank()==1){
			printf("%i %f %f %f %f\n",label,xi,yi,zi,ri);
			printf("parent %i %f %f %f %f\n",parent,xp,yp,zp,rp);
		}
		*/
        double length = sqrt((xi - xp) * (xi - xp) + (yi - yp) * (yi - yp) +
                             (zi - zp) * (zi - zp));
        if (length == 0.0)
            length = 1.0;
        double alpha = (xi - xp) / length;
        double beta = (yi - yp) / length;
        double gamma = (zi - zp) / length;

        int start_idx = int(xmin);
        int start_idy = int(ymin);
        int start_idz = int(zmin);
        int finish_idx = int(xmax);
        int finish_idy = int(ymax);
        int finish_idz = int(zmax);
        /* get the little box to loop over
		int start_idx = int((List_cx[idx] - List_rad[idx] - start_x)/voxel_length) + 1;
		int start_idy = int((List_cy[idx] - List_rad[idx] - start_y)/voxel_length) + 1;
		int start_idz = int((List_cz[idx] - List_rad[idx] - start_z)/voxel_length) + 1;
		int finish_idx = int((List_cx[idx] + List_rad[idx] - start_x)/voxel_length) + 1;
		int finish_idy = int((List_cy[idx] + List_rad[idx] - start_y)/voxel_length) + 1;
		int finish_idz = int((List_cz[idx] + List_rad[idx] - start_z)/voxel_length) + 1;
		*/

        if (start_idx < 0)
            start_idx = 0;
        if (start_idy < 0)
            start_idy = 0;
        if (start_idz < 0)
            start_idz = 0;
        if (start_idx > Nx - 1)
            start_idx = Nx;
        if (start_idy > Ny - 1)
            start_idy = Ny;
        if (start_idz > Nz - 1)
            start_idz = Nz;
        if (finish_idx < 0)
            finish_idx = 0;
        if (finish_idy < 0)
            finish_idy = 0;
        if (finish_idz < 0)
            finish_idz = 0;
        if (finish_idx > Nx - 1)
            finish_idx = Nx;
        if (finish_idy > Ny - 1)
            finish_idy = Ny;
        if (finish_idz > Nz - 1)
            finish_idz = Nz;

        /*	if (rank()==1) printf(" alpha = %f, beta = %f, gamma= %f\n",alpha, beta,gamma);
		if (rank()==1) printf(" xi = %f, yi = %f, zi= %f, ri = %f \n",xi, yi, zi, ri);
		if (rank()==1) printf(" xp = %f, yp = %f, zp= %f, rp = %f \n",xp, yp, zp, rp);

		if (rank()==1) printf( "start: %i, %i, %i \n",start_idx,start_idy,start_idz);
		if (rank()==1) printf( "finish: %i, %i, %i \n",finish_idx,finish_idy,finish_idz);
		*/

        for (int k = start_idz; k < finish_idz; k++) {
            for (int j = start_idy; j < finish_idy; j++) {
                for (int i = start_idx; i < finish_idx; i++) {

                    double x = i * voxel_length + start_x;
                    double y = j * voxel_length + start_y;
                    double z = k * voxel_length + start_z;

                    double distance;
                    double s = ((x - xp) * alpha + (y - yp) * beta +
                                (z - zp) * gamma) /
                               (alpha * alpha + beta * beta + gamma * gamma);

                    double di =
                        ri - sqrt((x - xi) * (x - xi) + (y - yi) * (y - yi) +
                                  (z - zi) * (z - zi));
                    double dp =
                        rp - sqrt((x - xp) * (x - xp) + (y - yp) * (y - yp) +
                                  (z - zp) * (z - zp));

                    if (s > length) {
                        distance = di;
                    } else if (s < 0.0) {
                        distance = dp;
                    } else {
                        // linear variation for radius
                        double radius = rp + (ri - rp) * s / length;
                        distance =
                            radius -
                            sqrt((x - xp - alpha * s) * (x - xp - alpha * s) +
                                 (y - yp - beta * s) * (y - yp - beta * s) +
                                 (z - zp - gamma * s) * (z - zp - gamma * s));
                    }
                    if (distance < di)
                        distance = di;
                    if (distance < dp)
                        distance = dp;

                    if (distance > 0.0) {
                        /* label the voxel */
                        //id[k*Nx*Ny + j*Nx + i] = label;
                        id[k * Nx * Ny + j * Nx + i] = 2;
                    }
                }
            }
        }
        //if (rank()==0) printf( "next line..\n");
    }
    delete[] List_cx;
    delete[] List_cy;
    delete[] List_cz;
    delete[] List_rad;
    delete[] List_index;
    delete[] List_type;
    delete[] List_parent;
}

/********************************************************
 * Constructors                                          *
 ********************************************************/
Domain::Domain(int nx, int ny, int nz, int rnk, int npx, int npy, int npz,
               double lx, double ly, double lz, int BC)
    : database(nullptr), Nx(0), Ny(0), Nz(0), Lx(0), Ly(0), Lz(0), Volume(0),
      BoundaryCondition(0), voxel_length(1),
      Comm(Utilities::MPI(MPI_COMM_WORLD).dup()), inlet_layers_x(0),
      inlet_layers_y(0), inlet_layers_z(0), inlet_layers_phase(1),
      outlet_layers_phase(2) {
    NULL_USE(rnk);
    NULL_USE(npy);
    NULL_USE(npz);
    // set up the neighbor ranks
    int myrank = Comm.getRank();
    rank_info =
        RankInfoStruct(myrank, rank_info.nx, rank_info.ny, rank_info.nz);

    Comm.barrier();

    auto db = std::make_shared<Database>();
    db->putScalar<int>("BC", BC);
    db->putVector<int>("nproc", {npx, npx, npx});
    db->putVector<int>("n", {nx, ny, nz});
    db->putScalar<int>("nspheres", 0);
    db->putVector<double>("L", {lx, ly, lz});
    initialize(db);
}
Domain::Domain(std::shared_ptr<Database> db, const Utilities::MPI &Communicator)
    : database(db), Nx(0), Ny(0), Nz(0), Lx(0), Ly(0), Lz(0), Volume(0),
      BoundaryCondition(0), inlet_layers_x(0), inlet_layers_y(0),
      inlet_layers_z(0), outlet_layers_x(0), outlet_layers_y(0),
      outlet_layers_z(0), inlet_layers_phase(1), outlet_layers_phase(2) {
    Comm = Communicator.dup();

    // set up the neighbor ranks
    int myrank = Comm.getRank();
    initialize(db);
    rank_info =
        RankInfoStruct(myrank, rank_info.nx, rank_info.ny, rank_info.nz);
    Comm.barrier();
}

/********************************************************
 * Destructor                                            *
 ********************************************************/
Domain::~Domain() {}

/********************************************************
 * Initialization                                        *
 ********************************************************/
void Domain::initialize(std::shared_ptr<Database> db) {
    d_db = db;
    auto nproc = d_db->getVector<int>("nproc");
    auto n = d_db->getVector<int>("n");

    ASSERT(n.size() == 3u);
    ASSERT(nproc.size() == 3u);
    int nx = n[0];
    int ny = n[1];
    int nz = n[2];
    offset_x = offset_y = offset_z = 0;

    if (d_db->keyExists("InletLayers")) {
        auto InletCount = d_db->getVector<int>("InletLayers");
        inlet_layers_x = InletCount[0];
        inlet_layers_y = InletCount[1];
        inlet_layers_z = InletCount[2];
    }
    if (d_db->keyExists("OutletLayers")) {
        auto OutletCount = d_db->getVector<int>("OutletLayers");
        outlet_layers_x = OutletCount[0];
        outlet_layers_y = OutletCount[1];
        outlet_layers_z = OutletCount[2];
    }
    if (d_db->keyExists("InletLayersPhase")) {
        inlet_layers_phase = d_db->getScalar<int>("InletLayersPhase");
    }
    if (d_db->keyExists("OutletLayersPhase")) {
        outlet_layers_phase = d_db->getScalar<int>("OutletLayersPhase");
    }
    voxel_length = 1.0;
    if (d_db->keyExists("voxel_length")) {
        voxel_length = d_db->getScalar<double>("voxel_length");
    } else if (d_db->keyExists("L")) {
        auto Length = d_db->getVector<double>("L");
        Lx = Length[0];
        Ly = Length[1];
        Lz = Length[2];
        voxel_length = Lx / (nx * nproc[0]);
    }
    Lx = nx * nproc[0] * voxel_length;
    Ly = ny * nproc[1] * voxel_length;
    Lz = nz * nproc[2] * voxel_length;
    Nx = nx + 2;
    Ny = ny + 2;
    Nz = nz + 2;
    // Initialize ranks
    int myrank = Comm.getRank();
    rank_info = RankInfoStruct(myrank, nproc[0], nproc[1], nproc[2]);

    // Fill remaining variables
    N = Nx * Ny * Nz;
    Volume = nx * ny * nz * nproc[0] * nproc[1] * nproc[2] * 1.0;

    if (myrank == 0)
        printf("voxel length = %f micron \n", voxel_length);

    id = std::vector<signed char>(N, 0);
    BoundaryCondition = d_db->getScalar<int>("BC");
    int nprocs = Comm.getSize();
    INSIST(nprocs == nproc[0] * nproc[1] * nproc[2],
           "Fatal error in processor count!");
}

/********************************************************
 * Get send/recv lists                                   *
 ********************************************************/
const std::vector<int> &Domain::getRecvList(const char *dir) const {
    if (dir[0] == 'x') {
        if (dir[1] == 0)
            return recvList_x;
        else if (dir[1] == 'y')
            return recvList_xy;
        else if (dir[1] == 'Y')
            return recvList_xY;
        else if (dir[1] == 'z')
            return recvList_xz;
        else if (dir[1] == 'Z')
            return recvList_xZ;
    } else if (dir[0] == 'y') {
        if (dir[1] == 0)
            return recvList_y;
        else if (dir[1] == 'z')
            return recvList_yz;
        else if (dir[1] == 'Z')
            return recvList_yZ;
    } else if (dir[0] == 'z') {
        if (dir[1] == 0)
            return recvList_z;
    } else if (dir[0] == 'X') {
        if (dir[1] == 0)
            return recvList_X;
        else if (dir[1] == 'y')
            return recvList_Xy;
        else if (dir[1] == 'Y')
            return recvList_XY;
        else if (dir[1] == 'z')
            return recvList_Xz;
        else if (dir[1] == 'Z')
            return recvList_XZ;
    } else if (dir[0] == 'Y') {
        if (dir[1] == 0)
            return recvList_Y;
        else if (dir[1] == 'z')
            return recvList_Yz;
        else if (dir[1] == 'Z')
            return recvList_YZ;
    } else if (dir[0] == 'Z') {
        if (dir[1] == 0)
            return recvList_Z;
    }
    throw std::logic_error("Internal error");
}
const std::vector<int> &Domain::getSendList(const char *dir) const {
    if (dir[0] == 'x') {
        if (dir[1] == 0)
            return sendList_x;
        else if (dir[1] == 'y')
            return sendList_xy;
        else if (dir[1] == 'Y')
            return sendList_xY;
        else if (dir[1] == 'z')
            return sendList_xz;
        else if (dir[1] == 'Z')
            return sendList_xZ;
    } else if (dir[0] == 'y') {
        if (dir[1] == 0)
            return sendList_y;
        else if (dir[1] == 'z')
            return sendList_yz;
        else if (dir[1] == 'Z')
            return sendList_yZ;
    } else if (dir[0] == 'z') {
        if (dir[1] == 0)
            return sendList_z;
    } else if (dir[0] == 'X') {
        if (dir[1] == 0)
            return sendList_X;
        else if (dir[1] == 'y')
            return sendList_Xy;
        else if (dir[1] == 'Y')
            return sendList_XY;
        else if (dir[1] == 'z')
            return sendList_Xz;
        else if (dir[1] == 'Z')
            return sendList_XZ;
    } else if (dir[0] == 'Y') {
        if (dir[1] == 0)
            return sendList_Y;
        else if (dir[1] == 'z')
            return sendList_Yz;
        else if (dir[1] == 'Z')
            return sendList_YZ;
    } else if (dir[0] == 'Z') {
        if (dir[1] == 0)
            return sendList_Z;
    }
    throw std::logic_error("Internal error");
}

/********************************************************
 * Decomp                                                *
 ********************************************************/
void Domain::Decomp(const std::string &Filename) {
    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int rank_offset = 0;
    int RANK = rank();
    int nprocs, nprocx, nprocy, nprocz, nx, ny, nz;
    int64_t global_Nx, global_Ny, global_Nz;
    int64_t i, j, k, n;
    int64_t xStart, yStart, zStart;
    int checkerSize;
    bool USE_CHECKER = false;
    //int inlet_layers_x, inlet_layers_y, inlet_layers_z;
    //int outlet_layers_x, outlet_layers_y, outlet_layers_z;
    xStart = yStart = zStart = 0;
    inlet_layers_x = 0;
    inlet_layers_y = 0;
    inlet_layers_z = 0;
    outlet_layers_x = 0;
    outlet_layers_y = 0;
    outlet_layers_z = 0;
    inlet_layers_phase = 1;
    outlet_layers_phase = 2;
    checkerSize = 32;

    // Read domain parameters
    //auto Filename = database->getScalar<std::string>( "Filename" );
    //auto L = database->getVector<double>( "L" );
    auto size = database->getVector<int>("n");
    auto SIZE = database->getVector<int>("N");
    auto nproc = database->getVector<int>("nproc");
    if (database->keyExists("offset")) {
        auto offset = database->getVector<int>("offset");
        xStart = offset[0];
        yStart = offset[1];
        zStart = offset[2];
        offset_x = xStart;
        offset_y = yStart;
        offset_z = zStart;
    }
    if (database->keyExists("InletLayers")) {
        auto InletCount = database->getVector<int>("InletLayers");
        inlet_layers_x = InletCount[0];
        inlet_layers_y = InletCount[1];
        inlet_layers_z = InletCount[2];
    }
    if (database->keyExists("OutletLayers")) {
        auto OutletCount = database->getVector<int>("OutletLayers");
        outlet_layers_x = OutletCount[0];
        outlet_layers_y = OutletCount[1];
        outlet_layers_z = OutletCount[2];
    }
    if (database->keyExists("checkerSize")) {
        checkerSize = database->getScalar<int>("checkerSize");
        USE_CHECKER = true;
    } else {
        checkerSize = SIZE[0];
    }
    if (database->keyExists("InletLayersPhase")) {
        inlet_layers_phase = database->getScalar<int>("InletLayersPhase");
    }
    if (database->keyExists("OutletLayersPhase")) {
        outlet_layers_phase = database->getScalar<int>("OutletLayersPhase");
    }
    auto ReadValues = database->getVector<int>("ReadValues");
    auto WriteValues = database->getVector<int>("WriteValues");
    auto ReadType = database->getScalar<std::string>("ReadType");

    if (ReadType == "8bit") {
    } else if (ReadType == "16bit") {
    } else if (ReadType == "swc") {
    } else {
        //printf("INPUT ERROR: Valid ReadType are 8bit, 16bit \n");
        ReadType = "8bit";
    }

    /* swc format for neurons */
    if (ReadType == "swc") {
        read_swc(Filename);
    } else {
        nx = size[0];
        ny = size[1];
        nz = size[2];
        nprocx = nproc[0];
        nprocy = nproc[1];
        nprocz = nproc[2];
        global_Nx = SIZE[0];
        global_Ny = SIZE[1];
        global_Nz = SIZE[2];
        nprocs = nprocx * nprocy * nprocz;
        char *SegData = NULL;

        if (RANK == 0) {
            printf("Input media: %s\n", Filename.c_str());
            printf("Relabeling %lu values\n", ReadValues.size());
            for (size_t idx = 0; idx < ReadValues.size(); idx++) {
                int oldvalue = ReadValues[idx];
                int newvalue = WriteValues[idx];
                printf("oldvalue=%d, newvalue =%d \n", oldvalue, newvalue);
            }

            // Rank=0 reads the entire segmented data and distributes to worker processes
            printf("Dimensions of segmented image: %ld x %ld x %ld \n",
                   global_Nx, global_Ny, global_Nz);
            int64_t SIZE = global_Nx * global_Ny * global_Nz;
            SegData = new char[SIZE];
            if (ReadType == "8bit") {
                printf("Reading 8-bit input data \n");
                FILE *SEGDAT = fopen(Filename.c_str(), "rb");
                if (SEGDAT == NULL)
                    ERROR("Domain.cpp: Error reading segmented data");
                size_t ReadSeg;
                ReadSeg = fread(SegData, 1, SIZE, SEGDAT);
                if (ReadSeg != size_t(SIZE))
                    printf("Domain.cpp: Error reading segmented data \n");
                fclose(SEGDAT);
            } else if (ReadType == "16bit") {
                printf("Reading 16-bit input data \n");
                short int *InputData;
                InputData = new short int[SIZE];
                FILE *SEGDAT = fopen(Filename.c_str(), "rb");
                if (SEGDAT == NULL)
                    ERROR("Domain.cpp: Error reading segmented data");
                size_t ReadSeg;
                ReadSeg = fread(InputData, 2, SIZE, SEGDAT);
                if (ReadSeg != size_t(SIZE))
                    printf("Domain.cpp: Error reading segmented data \n");
                fclose(SEGDAT);
                for (int n = 0; n < SIZE; n++) {
                    SegData[n] = char(InputData[n]);
                }
            } else if (ReadType == "SWC") {
            }
            printf("Read segmented data from %s \n", Filename.c_str());

            // relabel the data
            std::vector<long int> LabelCount(ReadValues.size(), 0);
            for (int k = 0; k < global_Nz; k++) {
                for (int j = 0; j < global_Ny; j++) {
                    for (int i = 0; i < global_Nx; i++) {
                        n = k * global_Nx * global_Ny + j * global_Nx + i;
                        //char locval = loc_id[n];
                        signed char locval = SegData[n];
                        for (size_t idx = 0; idx < ReadValues.size(); idx++) {
                            signed char oldvalue = ReadValues[idx];
                            signed char newvalue = WriteValues[idx];
                            if (locval == oldvalue) {
                                SegData[n] = newvalue;
                                LabelCount[idx]++;
                                idx = ReadValues.size();
                            }
                        }
                    }
                }
            }
            for (size_t idx = 0; idx < ReadValues.size(); idx++) {
                long int label = ReadValues[idx];
                long int count = LabelCount[idx];
                printf("Label=%ld, Count=%ld \n", label, count);
            }
            if (USE_CHECKER) {
                if (inlet_layers_x > 0) {
                    // use checkerboard pattern
                    printf("Checkerboard pattern at x inlet for %i layers \n",
                           inlet_layers_x);
                    for (int k = 0; k < global_Nz; k++) {
                        for (int j = 0; j < global_Ny; j++) {
                            for (int i = xStart; i < xStart + inlet_layers_x;
                                 i++) {
                                if ((j / checkerSize + k / checkerSize) % 2 ==
                                    0) {
                                    // void checkers
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = 2;
                                } else {
                                    // solid checkers
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = 0;
                                }
                            }
                        }
                    }
                }

                if (inlet_layers_y > 0) {
                    printf("Checkerboard pattern at y inlet for %i layers \n",
                           inlet_layers_y);
                    // use checkerboard pattern
                    for (int k = 0; k < global_Nz; k++) {
                        for (int j = yStart; j < yStart + inlet_layers_y; j++) {
                            for (int i = 0; i < global_Nx; i++) {
                                if ((i / checkerSize + k / checkerSize) % 2 ==
                                    0) {
                                    // void checkers
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = 2;
                                } else {
                                    // solid checkers
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = 0;
                                }
                            }
                        }
                    }
                }

                if (inlet_layers_z > 0) {
                    printf("Checkerboard pattern at z inlet for %i layers, "
                           "saturated with phase label=%i \n",
                           inlet_layers_z, inlet_layers_phase);
                    // use checkerboard pattern
                    for (int k = zStart; k < zStart + inlet_layers_z; k++) {
                        for (int j = 0; j < global_Ny; j++) {
                            for (int i = 0; i < global_Nx; i++) {
                                if ((i / checkerSize + j / checkerSize) % 2 ==
                                    0) {
                                    // void checkers
                                    //SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 2;
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] =
                                        inlet_layers_phase;
                                } else {
                                    // solid checkers
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = 0;
                                }
                            }
                        }
                    }
                }

                if (outlet_layers_x > 0) {
                    // use checkerboard pattern
                    printf("Checkerboard pattern at x outlet for %i layers \n",
                           outlet_layers_x);
                    for (int k = 0; k < global_Nz; k++) {
                        for (int j = 0; j < global_Ny; j++) {
                            for (int i = xStart + nx * nprocx - outlet_layers_x;
                                 i < xStart + nx * nprocx; i++) {
                                if ((j / checkerSize + k / checkerSize) % 2 ==
                                    0) {
                                    // void checkers
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = 2;
                                } else {
                                    // solid checkers
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = 0;
                                }
                            }
                        }
                    }
                }

                if (outlet_layers_y > 0) {
                    printf("Checkerboard pattern at y outlet for %i layers \n",
                           outlet_layers_y);
                    // use checkerboard pattern
                    for (int k = 0; k < global_Nz; k++) {
                        for (int j = yStart + ny * nprocy - outlet_layers_y;
                             j < yStart + ny * nprocy; j++) {
                            for (int i = 0; i < global_Nx; i++) {
                                if ((i / checkerSize + k / checkerSize) % 2 ==
                                    0) {
                                    // void checkers
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = 2;
                                } else {
                                    // solid checkers
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = 0;
                                }
                            }
                        }
                    }
                }

                if (outlet_layers_z > 0) {
                    printf("Checkerboard pattern at z outlet for %i layers, "
                           "saturated with phase label=%i \n",
                           outlet_layers_z, outlet_layers_phase);
                    // use checkerboard pattern
                    for (int k = zStart + nz * nprocz - outlet_layers_z;
                         k < zStart + nz * nprocz; k++) {
                        for (int j = 0; j < global_Ny; j++) {
                            for (int i = 0; i < global_Nx; i++) {
                                if ((i / checkerSize + j / checkerSize) % 2 ==
                                    0) {
                                    // void checkers
                                    //SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 2;
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] =
                                        outlet_layers_phase;
                                } else {
                                    // solid checkers
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = 0;
                                }
                            }
                        }
                    }
                }
            } else {
                if (inlet_layers_z > 0) {
                    printf("Mixed reflection pattern at z inlet for %i layers, "
                           "saturated with phase label=%i \n",
                           inlet_layers_z, inlet_layers_phase);
                    for (int k = zStart; k < zStart + inlet_layers_z; k++) {
                        for (int j = 0; j < global_Ny; j++) {
                            for (int i = 0; i < global_Nx; i++) {
                                signed char local_id =
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i];
                                signed char reflection_id =
                                    SegData[(zStart + nz * nprocz - 1) *
                                                global_Nx * global_Ny +
                                            j * global_Nx + i];
                                if (local_id < 1 && reflection_id > 0) {
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = reflection_id;
                                }
                            }
                        }
                    }
                }
                if (outlet_layers_z > 0) {
                    printf(
                        "Mixed reflection pattern at z outlet for %i layers, "
                        "saturated with phase label=%i \n",
                        outlet_layers_z, outlet_layers_phase);
                    for (int k = zStart + nz * nprocz - outlet_layers_z;
                         k < zStart + nz * nprocz; k++) {
                        for (int j = 0; j < global_Ny; j++) {
                            for (int i = 0; i < global_Nx; i++) {
                                signed char local_id =
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i];
                                signed char reflection_id =
                                    SegData[zStart * global_Nx * global_Ny +
                                            j * global_Nx + i];
                                if (local_id < 1 && reflection_id > 0) {
                                    SegData[k * global_Nx * global_Ny +
                                            j * global_Nx + i] = reflection_id;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Get the rank info
        int64_t N = (nx + 2) * (ny + 2) * (nz + 2);

        // number of sites to use for periodic boundary condition transition zone
        int64_t z_transition_size = (nprocz * nz - (global_Nz - zStart)) / 2;
        if (z_transition_size < 0)
            z_transition_size = 0;

        // Set up the sub-domains
        if (RANK == 0) {
            printf("Distributing subdomains across %i processors \n", nprocs);
            printf("Process grid: %i x %i x %i \n", nprocx, nprocy, nprocz);
            printf("Subdomain size: %i x %i x %i \n", nx, ny, nz);
            printf("Size of transition region: %ld \n", z_transition_size);
            auto loc_id = new char[(nx + 2) * (ny + 2) * (nz + 2)];
            for (int kp = 0; kp < nprocz; kp++) {
                for (int jp = 0; jp < nprocy; jp++) {
                    for (int ip = 0; ip < nprocx; ip++) {
                        // rank of the process that gets this subdomain
                        int rnk = kp * nprocx * nprocy + jp * nprocx + ip;
                        // Pack and send the subdomain for rnk
                        for (k = 0; k < nz + 2; k++) {
                            for (j = 0; j < ny + 2; j++) {
                                for (i = 0; i < nx + 2; i++) {
                                    int64_t x = xStart + ip * nx + i - 1;
                                    int64_t y = yStart + jp * ny + j - 1;
                                    // int64_t z = zStart + kp*nz + k-1;
                                    int64_t z = zStart + kp * nz + k - 1 -
                                                z_transition_size;
                                    if (x < xStart)
                                        x = xStart;
                                    if (!(x < global_Nx))
                                        x = global_Nx - 1;
                                    if (y < yStart)
                                        y = yStart;
                                    if (!(y < global_Ny))
                                        y = global_Ny - 1;
                                    if (z < zStart)
                                        z = zStart;
                                    if (!(z < global_Nz))
                                        z = global_Nz - 1;
                                    int64_t nlocal = k * (nx + 2) * (ny + 2) +
                                                     j * (nx + 2) + i;
                                    int64_t nglobal =
                                        z * global_Nx * global_Ny +
                                        y * global_Nx + x;
                                    loc_id[nlocal] = SegData[nglobal];
                                }
                            }
                        }
                        if (rnk == 0) {
                            for (k = 0; k < nz + 2; k++) {
                                for (j = 0; j < ny + 2; j++) {
                                    for (i = 0; i < nx + 2; i++) {
                                        int nlocal = k * (nx + 2) * (ny + 2) +
                                                     j * (nx + 2) + i;
                                        id[nlocal] = loc_id[nlocal];
                                    }
                                }
                            }
                        } else {
                            //printf("Sending data to process %i \n", rnk);
                            Comm.send(loc_id, N, rnk, 15);
                        }
                        // Write the data for this rank data
                        char LocalRankFilename[40];
                        sprintf(LocalRankFilename, "ID.%05i",
                                rnk + rank_offset);
                        FILE *ID = fopen(LocalRankFilename, "wb");
                        fwrite(loc_id, 1, (nx + 2) * (ny + 2) * (nz + 2), ID);
                        fclose(ID);
                    }
                }
            }
            delete[] loc_id;
        } else {
            // Recieve the subdomain from rank = 0
            //printf("Ready to recieve data %i at process %i \n", N,rank);
            Comm.recv(id.data(), N, 0, 15);
        }
        delete[] SegData;
    }
    /************************/
    // inlet layers only apply to lower part of domain
    if (rank_info.ix > 0)
        inlet_layers_x = 0;
    if (rank_info.jy > 0)
        inlet_layers_y = 0;
    if (rank_info.kz > 0)
        inlet_layers_z = 0;
    // outlet layers only apply to top part of domain
    if (rank_info.ix < nproc[0] - 1)
        outlet_layers_x = 0;
    if (rank_info.jy < nproc[1] - 1)
        outlet_layers_y = 0;
    if (rank_info.kz < nproc[2] - 1)
        outlet_layers_z = 0;
    /************************/
    Comm.barrier();
    ComputePorosity();
}

void Domain::ComputePorosity() {
    // Compute the porosity
    double sum;
    double sum_local = 0.0;
    double iVol_global = 1.0 / (1.0 * (Nx - 2) * (Ny - 2) * (Nz - 2) *
                                nprocx() * nprocy() * nprocz());
    if (BoundaryCondition > 0 && BoundaryCondition != 5)
        iVol_global =
            1.0 / (1.0 * (Nx - 2) * nprocx() * (Ny - 2) * nprocy() *
                   ((Nz - 2) * nprocz() - inlet_layers_z - outlet_layers_z));
    //.........................................................
    for (int k = inlet_layers_z + 1; k < Nz - outlet_layers_z - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                if (id[n] > 0) {
                    sum_local += 1.0;
                }
            }
        }
    }
    sum = Comm.sumReduce(sum_local);
    porosity = sum * iVol_global;
    if (rank() == 0)
        printf("Media porosity = %f \n", porosity);
    //.........................................................
}

void Domain::AggregateLabels(const std::string &filename) {

    int nx = Nx;
    int ny = Ny;
    int nz = Nz;

    int npx = nprocx();
    int npy = nprocy();
    int npz = nprocz();

    int ipx = iproc();
    int ipy = jproc();
    int ipz = kproc();

    int nprocs = nprocx() * nprocy() * nprocz();

    int full_nx = npx * (nx - 2);
    int full_ny = npy * (ny - 2);
    int full_nz = npz * (nz - 2);
    int local_size = (nx - 2) * (ny - 2) * (nz - 2);
    long int full_size = long(full_nx) * long(full_ny) * long(full_nz);

    auto LocalID = new signed char[local_size];

    //printf("aggregate labels: local size=%i, global size = %i",local_size, full_size);
    // assign the ID for the local sub-region
    for (int k = 1; k < nz - 1; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                int n = k * nx * ny + j * nx + i;
                signed char local_id_val = id[n];
                LocalID[(k - 1) * (nx - 2) * (ny - 2) + (j - 1) * (nx - 2) + i -
                        1] = local_id_val;
            }
        }
    }
    Comm.barrier();

    // populate the FullID
    if (rank() == 0) {
        auto FullID = new signed char[full_size];
        // first handle local ID for rank 0
        for (int k = 1; k < nz - 1; k++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    int x = i - 1;
                    int y = j - 1;
                    int z = k - 1;
                    int n_local = (k - 1) * (nx - 2) * (ny - 2) +
                                  (j - 1) * (nx - 2) + i - 1;
                    int n_full = z * full_nx * full_ny + y * full_nx + x;
                    FullID[n_full] = LocalID[n_local];
                }
            }
        }
        // next get the local ID from the other ranks
        for (int rnk = 1; rnk < nprocs; rnk++) {
            ipz = rnk / (npx * npy);
            ipy = (rnk - ipz * npx * npy) / npx;
            ipx = (rnk - ipz * npx * npy - ipy * npx);
            //printf("ipx=%i ipy=%i ipz=%i\n", ipx, ipy, ipz);
            int tag = 15 + rnk;
            Comm.recv(LocalID, local_size, rnk, tag);
            for (int k = 1; k < nz - 1; k++) {
                for (int j = 1; j < ny - 1; j++) {
                    for (int i = 1; i < nx - 1; i++) {
                        int x = i - 1 + ipx * (nx - 2);
                        int y = j - 1 + ipy * (ny - 2);
                        int z = k - 1 + ipz * (nz - 2);
                        int n_local = (k - 1) * (nx - 2) * (ny - 2) +
                                      (j - 1) * (nx - 2) + i - 1;
                        int n_full = z * full_nx * full_ny + y * full_nx + x;
                        FullID[n_full] = LocalID[n_local];
                    }
                }
            }
        }
        // write the output
        FILE *OUTFILE = fopen(filename.c_str(), "wb");
        fwrite(FullID, 1, full_size, OUTFILE);
        fclose(OUTFILE);
        delete[] FullID;
    } else {
        // send LocalID to rank=0
        int tag = 15 + rank();
        int dstrank = 0;
        Comm.send(LocalID, local_size, dstrank, tag);
    }
    delete[] LocalID;
    Comm.barrier();
}

/********************************************************
 * Initialize communication                              *
 ********************************************************/
void Domain::CommInit() {
    int i, j, k, n;
    int sendtag = 21;
    int recvtag = 21;
    //......................................................................................
    int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y,
        sendCount_Z;
    int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz,
        sendCount_xZ;
    int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ,
        sendCount_XZ;
    sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y =
        sendCount_Z = 0;
    sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz =
        sendCount_xZ = 0;
    sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ =
        sendCount_XZ = 0;
    //......................................................................................
    for (k = 1; k < Nz - 1; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                // Check the phase ID
                if (id[k * Nx * Ny + j * Nx + i] > 0) {
                    // Counts for the six faces
                    if (i == 1)
                        sendCount_x++;
                    if (j == 1)
                        sendCount_y++;
                    if (k == 1)
                        sendCount_z++;
                    if (i == Nx - 2)
                        sendCount_X++;
                    if (j == Ny - 2)
                        sendCount_Y++;
                    if (k == Nz - 2)
                        sendCount_Z++;
                    // Counts for the twelve edges
                    if (i == 1 && j == 1)
                        sendCount_xy++;
                    if (i == 1 && j == Ny - 2)
                        sendCount_xY++;
                    if (i == Nx - 2 && j == 1)
                        sendCount_Xy++;
                    if (i == Nx - 2 && j == Ny - 2)
                        sendCount_XY++;

                    if (i == 1 && k == 1)
                        sendCount_xz++;
                    if (i == 1 && k == Nz - 2)
                        sendCount_xZ++;
                    if (i == Nx - 2 && k == 1)
                        sendCount_Xz++;
                    if (i == Nx - 2 && k == Nz - 2)
                        sendCount_XZ++;

                    if (j == 1 && k == 1)
                        sendCount_yz++;
                    if (j == 1 && k == Nz - 2)
                        sendCount_yZ++;
                    if (j == Ny - 2 && k == 1)
                        sendCount_Yz++;
                    if (j == Ny - 2 && k == Nz - 2)
                        sendCount_YZ++;
                }
            }
        }
    }

    // allocate send lists
    sendList_x.resize(sendCount_x, 0);
    sendList_y.resize(sendCount_y, 0);
    sendList_z.resize(sendCount_z, 0);
    sendList_X.resize(sendCount_X, 0);
    sendList_Y.resize(sendCount_Y, 0);
    sendList_Z.resize(sendCount_Z, 0);
    sendList_xy.resize(sendCount_xy, 0);
    sendList_yz.resize(sendCount_yz, 0);
    sendList_xz.resize(sendCount_xz, 0);
    sendList_Xy.resize(sendCount_Xy, 0);
    sendList_Yz.resize(sendCount_Yz, 0);
    sendList_xZ.resize(sendCount_xZ, 0);
    sendList_xY.resize(sendCount_xY, 0);
    sendList_yZ.resize(sendCount_yZ, 0);
    sendList_Xz.resize(sendCount_Xz, 0);
    sendList_XY.resize(sendCount_XY, 0);
    sendList_YZ.resize(sendCount_YZ, 0);
    sendList_XZ.resize(sendCount_XZ, 0);
    // Populate the send list
    sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y =
        sendCount_Z = 0;
    sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz =
        sendCount_xZ = 0;
    sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ =
        sendCount_XZ = 0;
    for (k = 1; k < Nz - 1; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                // Local value to send
                n = k * Nx * Ny + j * Nx + i;
                if (id[n] > 0) {
                    // Counts for the six faces
                    if (i == 1)
                        sendList_x[sendCount_x++] = n;
                    if (j == 1)
                        sendList_y[sendCount_y++] = n;
                    if (k == 1)
                        sendList_z[sendCount_z++] = n;
                    if (i == Nx - 2)
                        sendList_X[sendCount_X++] = n;
                    if (j == Ny - 2)
                        sendList_Y[sendCount_Y++] = n;
                    if (k == Nz - 2)
                        sendList_Z[sendCount_Z++] = n;
                    // Counts for the twelve edges
                    if (i == 1 && j == 1)
                        sendList_xy[sendCount_xy++] = n;
                    if (i == 1 && j == Ny - 2)
                        sendList_xY[sendCount_xY++] = n;
                    if (i == Nx - 2 && j == 1)
                        sendList_Xy[sendCount_Xy++] = n;
                    if (i == Nx - 2 && j == Ny - 2)
                        sendList_XY[sendCount_XY++] = n;

                    if (i == 1 && k == 1)
                        sendList_xz[sendCount_xz++] = n;
                    if (i == 1 && k == Nz - 2)
                        sendList_xZ[sendCount_xZ++] = n;
                    if (i == Nx - 2 && k == 1)
                        sendList_Xz[sendCount_Xz++] = n;
                    if (i == Nx - 2 && k == Nz - 2)
                        sendList_XZ[sendCount_XZ++] = n;

                    if (j == 1 && k == 1)
                        sendList_yz[sendCount_yz++] = n;
                    if (j == 1 && k == Nz - 2)
                        sendList_yZ[sendCount_yZ++] = n;
                    if (j == Ny - 2 && k == 1)
                        sendList_Yz[sendCount_Yz++] = n;
                    if (j == Ny - 2 && k == Nz - 2)
                        sendList_YZ[sendCount_YZ++] = n;
                }
            }
        }
    }

    //......................................................................................
    int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y,
        recvCount_Z;
    int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz,
        recvCount_xZ;
    int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ,
        recvCount_XZ;
    req1[0] = Comm.Isend(&sendCount_x, 1, rank_x(), sendtag + 0);
    req2[0] = Comm.Irecv(&recvCount_X, 1, rank_X(), recvtag + 0);
    req1[1] = Comm.Isend(&sendCount_X, 1, rank_X(), sendtag + 1);
    req2[1] = Comm.Irecv(&recvCount_x, 1, rank_x(), recvtag + 1);
    req1[2] = Comm.Isend(&sendCount_y, 1, rank_y(), sendtag + 2);
    req2[2] = Comm.Irecv(&recvCount_Y, 1, rank_Y(), recvtag + 2);
    req1[3] = Comm.Isend(&sendCount_Y, 1, rank_Y(), sendtag + 3);
    req2[3] = Comm.Irecv(&recvCount_y, 1, rank_y(), recvtag + 3);
    req1[4] = Comm.Isend(&sendCount_z, 1, rank_z(), sendtag + 4);
    req2[4] = Comm.Irecv(&recvCount_Z, 1, rank_Z(), recvtag + 4);
    req1[5] = Comm.Isend(&sendCount_Z, 1, rank_Z(), sendtag + 5);
    req2[5] = Comm.Irecv(&recvCount_z, 1, rank_z(), recvtag + 5);
    req1[6] = Comm.Isend(&sendCount_xy, 1, rank_xy(), sendtag + 6);
    req2[6] = Comm.Irecv(&recvCount_XY, 1, rank_XY(), recvtag + 6);
    req1[7] = Comm.Isend(&sendCount_XY, 1, rank_XY(), sendtag + 7);
    req2[7] = Comm.Irecv(&recvCount_xy, 1, rank_xy(), recvtag + 7);
    req1[8] = Comm.Isend(&sendCount_Xy, 1, rank_Xy(), sendtag + 8);
    req2[8] = Comm.Irecv(&recvCount_xY, 1, rank_xY(), recvtag + 8);
    req1[9] = Comm.Isend(&sendCount_xY, 1, rank_xY(), sendtag + 9);
    req2[9] = Comm.Irecv(&recvCount_Xy, 1, rank_Xy(), recvtag + 9);
    req1[10] = Comm.Isend(&sendCount_xz, 1, rank_xz(), sendtag + 10);
    req2[10] = Comm.Irecv(&recvCount_XZ, 1, rank_XZ(), recvtag + 10);
    req1[11] = Comm.Isend(&sendCount_XZ, 1, rank_XZ(), sendtag + 11);
    req2[11] = Comm.Irecv(&recvCount_xz, 1, rank_xz(), recvtag + 11);
    req1[12] = Comm.Isend(&sendCount_Xz, 1, rank_Xz(), sendtag + 12);
    req2[12] = Comm.Irecv(&recvCount_xZ, 1, rank_xZ(), recvtag + 12);
    req1[13] = Comm.Isend(&sendCount_xZ, 1, rank_xZ(), sendtag + 13);
    req2[13] = Comm.Irecv(&recvCount_Xz, 1, rank_Xz(), recvtag + 13);
    req1[14] = Comm.Isend(&sendCount_yz, 1, rank_yz(), sendtag + 14);
    req2[14] = Comm.Irecv(&recvCount_YZ, 1, rank_YZ(), recvtag + 14);
    req1[15] = Comm.Isend(&sendCount_YZ, 1, rank_YZ(), sendtag + 15);
    req2[15] = Comm.Irecv(&recvCount_yz, 1, rank_yz(), recvtag + 15);
    req1[16] = Comm.Isend(&sendCount_Yz, 1, rank_Yz(), sendtag + 16);
    req2[16] = Comm.Irecv(&recvCount_yZ, 1, rank_yZ(), recvtag + 16);
    req1[17] = Comm.Isend(&sendCount_yZ, 1, rank_yZ(), sendtag + 17);
    req2[17] = Comm.Irecv(&recvCount_Yz, 1, rank_Yz(), recvtag + 17);
    Comm.waitAll(18, req1);
    Comm.waitAll(18, req2);
    Comm.barrier();
    // allocate recv lists
    recvList_x.resize(recvCount_x, 0);
    recvList_y.resize(recvCount_y, 0);
    recvList_z.resize(recvCount_z, 0);
    recvList_X.resize(recvCount_X, 0);
    recvList_Y.resize(recvCount_Y, 0);
    recvList_Z.resize(recvCount_Z, 0);
    recvList_xy.resize(recvCount_xy, 0);
    recvList_yz.resize(recvCount_yz, 0);
    recvList_xz.resize(recvCount_xz, 0);
    recvList_Xy.resize(recvCount_Xy, 0);
    recvList_Yz.resize(recvCount_Yz, 0);
    recvList_xZ.resize(recvCount_xZ, 0);
    recvList_xY.resize(recvCount_xY, 0);
    recvList_yZ.resize(recvCount_yZ, 0);
    recvList_Xz.resize(recvCount_Xz, 0);
    recvList_XY.resize(recvCount_XY, 0);
    recvList_YZ.resize(recvCount_YZ, 0);
    recvList_XZ.resize(recvCount_XZ, 0);
    //......................................................................................
    req1[0] = Comm.Isend(sendList_x.data(), sendCount_x, rank_x(), sendtag);
    req2[0] = Comm.Irecv(recvList_X.data(), recvCount_X, rank_X(), recvtag);
    req1[1] = Comm.Isend(sendList_X.data(), sendCount_X, rank_X(), sendtag);
    req2[1] = Comm.Irecv(recvList_x.data(), recvCount_x, rank_x(), recvtag);
    req1[2] = Comm.Isend(sendList_y.data(), sendCount_y, rank_y(), sendtag);
    req2[2] = Comm.Irecv(recvList_Y.data(), recvCount_Y, rank_Y(), recvtag);
    req1[3] = Comm.Isend(sendList_Y.data(), sendCount_Y, rank_Y(), sendtag);
    req2[3] = Comm.Irecv(recvList_y.data(), recvCount_y, rank_y(), recvtag);
    req1[4] = Comm.Isend(sendList_z.data(), sendCount_z, rank_z(), sendtag);
    req2[4] = Comm.Irecv(recvList_Z.data(), recvCount_Z, rank_Z(), recvtag);
    req1[5] = Comm.Isend(sendList_Z.data(), sendCount_Z, rank_Z(), sendtag);
    req2[5] = Comm.Irecv(recvList_z.data(), recvCount_z, rank_z(), recvtag);
    req1[6] = Comm.Isend(sendList_xy.data(), sendCount_xy, rank_xy(), sendtag);
    req2[6] = Comm.Irecv(recvList_XY.data(), recvCount_XY, rank_XY(), recvtag);
    req1[7] = Comm.Isend(sendList_XY.data(), sendCount_XY, rank_XY(), sendtag);
    req2[7] = Comm.Irecv(recvList_xy.data(), recvCount_xy, rank_xy(), recvtag);
    req1[8] = Comm.Isend(sendList_Xy.data(), sendCount_Xy, rank_Xy(), sendtag);
    req2[8] = Comm.Irecv(recvList_xY.data(), recvCount_xY, rank_xY(), recvtag);
    req1[9] = Comm.Isend(sendList_xY.data(), sendCount_xY, rank_xY(), sendtag);
    req2[9] = Comm.Irecv(recvList_Xy.data(), recvCount_Xy, rank_Xy(), recvtag);
    req1[10] = Comm.Isend(sendList_xz.data(), sendCount_xz, rank_xz(), sendtag);
    req2[10] = Comm.Irecv(recvList_XZ.data(), recvCount_XZ, rank_XZ(), recvtag);
    req1[11] = Comm.Isend(sendList_XZ.data(), sendCount_XZ, rank_XZ(), sendtag);
    req2[11] = Comm.Irecv(recvList_xz.data(), recvCount_xz, rank_xz(), recvtag);
    req1[12] = Comm.Isend(sendList_Xz.data(), sendCount_Xz, rank_Xz(), sendtag);
    req2[12] = Comm.Irecv(recvList_xZ.data(), recvCount_xZ, rank_xZ(), recvtag);
    req1[13] = Comm.Isend(sendList_xZ.data(), sendCount_xZ, rank_xZ(), sendtag);
    req2[13] = Comm.Irecv(recvList_Xz.data(), recvCount_Xz, rank_Xz(), recvtag);
    req1[14] = Comm.Isend(sendList_yz.data(), sendCount_yz, rank_yz(), sendtag);
    req2[14] = Comm.Irecv(recvList_YZ.data(), recvCount_YZ, rank_YZ(), recvtag);
    req1[15] = Comm.Isend(sendList_YZ.data(), sendCount_YZ, rank_YZ(), sendtag);
    req2[15] = Comm.Irecv(recvList_yz.data(), recvCount_yz, rank_yz(), recvtag);
    req1[16] = Comm.Isend(sendList_Yz.data(), sendCount_Yz, rank_Yz(), sendtag);
    req2[16] = Comm.Irecv(recvList_yZ.data(), recvCount_yZ, rank_yZ(), recvtag);
    req1[17] = Comm.Isend(sendList_yZ.data(), sendCount_yZ, rank_yZ(), sendtag);
    req2[17] = Comm.Irecv(recvList_Yz.data(), recvCount_Yz, rank_Yz(), recvtag);
    Comm.waitAll(18, req1);
    Comm.waitAll(18, req2);
    //......................................................................................
    for (int idx = 0; idx < recvCount_x; idx++)
        recvList_x[idx] -= (Nx - 2);
    for (int idx = 0; idx < recvCount_X; idx++)
        recvList_X[idx] += (Nx - 2);
    for (int idx = 0; idx < recvCount_y; idx++)
        recvList_y[idx] -= (Ny - 2) * Nx;
    for (int idx = 0; idx < recvCount_Y; idx++)
        recvList_Y[idx] += (Ny - 2) * Nx;
    for (int idx = 0; idx < recvCount_z; idx++)
        recvList_z[idx] -= (Nz - 2) * Nx * Ny;
    for (int idx = 0; idx < recvCount_Z; idx++)
        recvList_Z[idx] += (Nz - 2) * Nx * Ny;
    for (int idx = 0; idx < recvCount_xy; idx++)
        recvList_xy[idx] -= (Nx - 2) + (Ny - 2) * Nx;
    for (int idx = 0; idx < recvCount_XY; idx++)
        recvList_XY[idx] += (Nx - 2) + (Ny - 2) * Nx;
    for (int idx = 0; idx < recvCount_xY; idx++)
        recvList_xY[idx] -= (Nx - 2) - (Ny - 2) * Nx;
    for (int idx = 0; idx < recvCount_Xy; idx++)
        recvList_Xy[idx] += (Nx - 2) - (Ny - 2) * Nx;
    for (int idx = 0; idx < recvCount_xz; idx++)
        recvList_xz[idx] -= (Nx - 2) + (Nz - 2) * Nx * Ny;
    for (int idx = 0; idx < recvCount_XZ; idx++)
        recvList_XZ[idx] += (Nx - 2) + (Nz - 2) * Nx * Ny;
    for (int idx = 0; idx < recvCount_xZ; idx++)
        recvList_xZ[idx] -= (Nx - 2) - (Nz - 2) * Nx * Ny;
    for (int idx = 0; idx < recvCount_Xz; idx++)
        recvList_Xz[idx] += (Nx - 2) - (Nz - 2) * Nx * Ny;
    for (int idx = 0; idx < recvCount_yz; idx++)
        recvList_yz[idx] -= (Ny - 2) * Nx + (Nz - 2) * Nx * Ny;
    for (int idx = 0; idx < recvCount_YZ; idx++)
        recvList_YZ[idx] += (Ny - 2) * Nx + (Nz - 2) * Nx * Ny;
    for (int idx = 0; idx < recvCount_yZ; idx++)
        recvList_yZ[idx] -= (Ny - 2) * Nx - (Nz - 2) * Nx * Ny;
    for (int idx = 0; idx < recvCount_Yz; idx++)
        recvList_Yz[idx] += (Ny - 2) * Nx - (Nz - 2) * Nx * Ny;
    //......................................................................................

    //......................................................................................
}

void Domain::ReadIDs() {
    // Read the IDs from input file
    int nprocs = nprocx() * nprocy() * nprocz();
    size_t readID;
    char LocalRankString[8];
    char LocalRankFilename[40];
    //.......................................................................
    if (rank() == 0)
        printf("Read input media... \n");
    //.......................................................................
    sprintf(LocalRankString, "%05d", rank());
    sprintf(LocalRankFilename, "%s%s", "ID.", LocalRankString);
    // .......... READ THE INPUT FILE .......................................
    if (rank() == 0)
        printf("Initialize from segmented data: solid=0, NWP=1, WP=2 \n");
    sprintf(LocalRankFilename, "ID.%05i", rank());
    FILE *IDFILE = fopen(LocalRankFilename, "rb");
    if (!IDFILE)
        ERROR("Domain::ReadIDs --  Error opening file: ID.xxxxx");
    readID = fread(id.data(), 1, N, IDFILE);
    if (readID != size_t(N))
        printf("Domain::ReadIDs -- Error reading ID (rank=%i) \n", rank());
    fclose(IDFILE);

    // Compute the porosity
    double sum;
    double sum_local = 0.0;
    double iVol_global = 1.0 / (1.0 * (Nx - 2) * (Ny - 2) * (Nz - 2) * nprocs);
    if (BoundaryCondition > 0)
        iVol_global = 1.0 / (1.0 * (Nx - 2) * nprocx() * (Ny - 2) * nprocy() *
                             ((Nz - 2) * nprocz() - 6));
    //.........................................................
    // If external boundary conditions are applied remove solid
    if (BoundaryCondition > 0 && kproc() == 0) {
        if (inlet_layers_z < 4)
            inlet_layers_z = 4;
        for (int k = 0; k < inlet_layers_z; k++) {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int n = k * Nx * Ny + j * Nx + i;
                    id[n] = 1;
                }
            }
        }
    }
    if (BoundaryCondition > 0 && kproc() == nprocz() - 1) {
        if (outlet_layers_z < 4)
            outlet_layers_z = 4;
        for (int k = Nz - outlet_layers_z; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int n = k * Nx * Ny + j * Nx + i;
                    id[n] = 2;
                }
            }
        }
    }
    for (int k = inlet_layers_z + 1; k < Nz - outlet_layers_z - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                if (id[n] > 0) {
                    sum_local += 1.0;
                }
            }
        }
    }
    sum = Comm.sumReduce(sum_local);
    porosity = sum * iVol_global;
    if (rank() == 0)
        printf("Media porosity = %f \n", porosity);
    //.........................................................
}
int Domain::PoreCount() {
    /*
	 * count the number of nodes occupied by mobile phases
	 */
    int Npore = 0; // number of local pore nodes
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                if (id[n] > 0) {
                    Npore++;
                }
            }
        }
    }
    return Npore;
}

void Domain::CommunicateMeshHalo(DoubleArray &Mesh) {
    int sendtag, recvtag;
    sendtag = recvtag = 7;
    double *MeshData = Mesh.data();
    // send buffers
    auto sendData_x = new double[sendCount("x")];
    auto sendData_y = new double[sendCount("y")];
    auto sendData_z = new double[sendCount("z")];
    auto sendData_X = new double[sendCount("X")];
    auto sendData_Y = new double[sendCount("Y")];
    auto sendData_Z = new double[sendCount("Z")];
    auto sendData_xy = new double[sendCount("xy")];
    auto sendData_yz = new double[sendCount("yz")];
    auto sendData_xz = new double[sendCount("xz")];
    auto sendData_Xy = new double[sendCount("Xy")];
    auto sendData_Yz = new double[sendCount("Yz")];
    auto sendData_xZ = new double[sendCount("xZ")];
    auto sendData_xY = new double[sendCount("xY")];
    auto sendData_yZ = new double[sendCount("yZ")];
    auto sendData_Xz = new double[sendCount("Xz")];
    auto sendData_XY = new double[sendCount("XY")];
    auto sendData_YZ = new double[sendCount("YZ")];
    auto sendData_XZ = new double[sendCount("XZ")];
    // recv buffers
    auto recvData_x = new double[recvCount("x")];
    auto recvData_y = new double[recvCount("y")];
    auto recvData_z = new double[recvCount("z")];
    auto recvData_X = new double[recvCount("X")];
    auto recvData_Y = new double[recvCount("Y")];
    auto recvData_Z = new double[recvCount("Z")];
    auto recvData_xy = new double[recvCount("xy")];
    auto recvData_yz = new double[recvCount("yz")];
    auto recvData_xz = new double[recvCount("xz")];
    auto recvData_Xy = new double[recvCount("Xy")];
    auto recvData_xZ = new double[recvCount("xZ")];
    auto recvData_xY = new double[recvCount("xY")];
    auto recvData_yZ = new double[recvCount("yZ")];
    auto recvData_Yz = new double[recvCount("Yz")];
    auto recvData_Xz = new double[recvCount("Xz")];
    auto recvData_XY = new double[recvCount("XY")];
    auto recvData_YZ = new double[recvCount("YZ")];
    auto recvData_XZ = new double[recvCount("XZ")];
    // Pack data
    PackMeshData(sendList("x"), sendCount("x"), sendData_x, MeshData);
    PackMeshData(sendList("X"), sendCount("X"), sendData_X, MeshData);
    PackMeshData(sendList("y"), sendCount("y"), sendData_y, MeshData);
    PackMeshData(sendList("Y"), sendCount("Y"), sendData_Y, MeshData);
    PackMeshData(sendList("z"), sendCount("z"), sendData_z, MeshData);
    PackMeshData(sendList("Z"), sendCount("Z"), sendData_Z, MeshData);
    PackMeshData(sendList("xy"), sendCount("xy"), sendData_xy, MeshData);
    PackMeshData(sendList("Xy"), sendCount("Xy"), sendData_Xy, MeshData);
    PackMeshData(sendList("xY"), sendCount("xY"), sendData_xY, MeshData);
    PackMeshData(sendList("XY"), sendCount("XY"), sendData_XY, MeshData);
    PackMeshData(sendList("xz"), sendCount("xz"), sendData_xz, MeshData);
    PackMeshData(sendList("Xz"), sendCount("Xz"), sendData_Xz, MeshData);
    PackMeshData(sendList("xZ"), sendCount("xZ"), sendData_xZ, MeshData);
    PackMeshData(sendList("XZ"), sendCount("XZ"), sendData_XZ, MeshData);
    PackMeshData(sendList("yz"), sendCount("yz"), sendData_yz, MeshData);
    PackMeshData(sendList("Yz"), sendCount("Yz"), sendData_Yz, MeshData);
    PackMeshData(sendList("yZ"), sendCount("yZ"), sendData_yZ, MeshData);
    PackMeshData(sendList("YZ"), sendCount("YZ"), sendData_YZ, MeshData);
    // send/recv
    Comm.sendrecv(sendData_x, sendCount("x"), rank_x(), sendtag, recvData_X,
                  recvCount("X"), rank_X(), recvtag);
    Comm.sendrecv(sendData_X, sendCount("X"), rank_X(), sendtag, recvData_x,
                  recvCount("x"), rank_x(), recvtag);
    Comm.sendrecv(sendData_y, sendCount("y"), rank_y(), sendtag, recvData_Y,
                  recvCount("Y"), rank_Y(), recvtag);
    Comm.sendrecv(sendData_Y, sendCount("Y"), rank_Y(), sendtag, recvData_y,
                  recvCount("y"), rank_y(), recvtag);
    Comm.sendrecv(sendData_z, sendCount("z"), rank_z(), sendtag, recvData_Z,
                  recvCount("Z"), rank_Z(), recvtag);
    Comm.sendrecv(sendData_Z, sendCount("Z"), rank_Z(), sendtag, recvData_z,
                  recvCount("z"), rank_z(), recvtag);
    Comm.sendrecv(sendData_xy, sendCount("xy"), rank_xy(), sendtag, recvData_XY,
                  recvCount("XY"), rank_XY(), recvtag);
    Comm.sendrecv(sendData_XY, sendCount("XY"), rank_XY(), sendtag, recvData_xy,
                  recvCount("xy"), rank_xy(), recvtag);
    Comm.sendrecv(sendData_Xy, sendCount("Xy"), rank_Xy(), sendtag, recvData_xY,
                  recvCount("xY"), rank_xY(), recvtag);
    Comm.sendrecv(sendData_xY, sendCount("xY"), rank_xY(), sendtag, recvData_Xy,
                  recvCount("Xy"), rank_Xy(), recvtag);
    Comm.sendrecv(sendData_xz, sendCount("xz"), rank_xz(), sendtag, recvData_XZ,
                  recvCount("XZ"), rank_XZ(), recvtag);
    Comm.sendrecv(sendData_XZ, sendCount("XZ"), rank_XZ(), sendtag, recvData_xz,
                  recvCount("xz"), rank_xz(), recvtag);
    Comm.sendrecv(sendData_Xz, sendCount("Xz"), rank_Xz(), sendtag, recvData_xZ,
                  recvCount("xZ"), rank_xZ(), recvtag);
    Comm.sendrecv(sendData_xZ, sendCount("xZ"), rank_xZ(), sendtag, recvData_Xz,
                  recvCount("Xz"), rank_Xz(), recvtag);
    Comm.sendrecv(sendData_yz, sendCount("yz"), rank_yz(), sendtag, recvData_YZ,
                  recvCount("YZ"), rank_YZ(), recvtag);
    Comm.sendrecv(sendData_YZ, sendCount("YZ"), rank_YZ(), sendtag, recvData_yz,
                  recvCount("yz"), rank_yz(), recvtag);
    Comm.sendrecv(sendData_Yz, sendCount("Yz"), rank_Yz(), sendtag, recvData_yZ,
                  recvCount("yZ"), rank_yZ(), recvtag);
    Comm.sendrecv(sendData_yZ, sendCount("yZ"), rank_yZ(), sendtag, recvData_Yz,
                  recvCount("Yz"), rank_Yz(), recvtag);
    // unpack data
    UnpackMeshData(recvList("x"), recvCount("x"), recvData_x, MeshData);
    UnpackMeshData(recvList("X"), recvCount("X"), recvData_X, MeshData);
    UnpackMeshData(recvList("y"), recvCount("y"), recvData_y, MeshData);
    UnpackMeshData(recvList("Y"), recvCount("Y"), recvData_Y, MeshData);
    UnpackMeshData(recvList("z"), recvCount("z"), recvData_z, MeshData);
    UnpackMeshData(recvList("Z"), recvCount("Z"), recvData_Z, MeshData);
    UnpackMeshData(recvList("xy"), recvCount("xy"), recvData_xy, MeshData);
    UnpackMeshData(recvList("Xy"), recvCount("Xy"), recvData_Xy, MeshData);
    UnpackMeshData(recvList("xY"), recvCount("xY"), recvData_xY, MeshData);
    UnpackMeshData(recvList("XY"), recvCount("XY"), recvData_XY, MeshData);
    UnpackMeshData(recvList("xz"), recvCount("xz"), recvData_xz, MeshData);
    UnpackMeshData(recvList("Xz"), recvCount("Xz"), recvData_Xz, MeshData);
    UnpackMeshData(recvList("xZ"), recvCount("xZ"), recvData_xZ, MeshData);
    UnpackMeshData(recvList("XZ"), recvCount("XZ"), recvData_XZ, MeshData);
    UnpackMeshData(recvList("yz"), recvCount("yz"), recvData_yz, MeshData);
    UnpackMeshData(recvList("Yz"), recvCount("Yz"), recvData_Yz, MeshData);
    UnpackMeshData(recvList("yZ"), recvCount("yZ"), recvData_yZ, MeshData);
    UnpackMeshData(recvList("YZ"), recvCount("YZ"), recvData_YZ, MeshData);
    // Free sendData
    delete[] sendData_x;
    delete[] sendData_y;
    delete[] sendData_z;
    delete[] sendData_X;
    delete[] sendData_Y;
    delete[] sendData_Z;
    delete[] sendData_xy;
    delete[] sendData_xY;
    delete[] sendData_Xy;
    delete[] sendData_XY;
    delete[] sendData_xz;
    delete[] sendData_xZ;
    delete[] sendData_Xz;
    delete[] sendData_XZ;
    delete[] sendData_yz;
    delete[] sendData_yZ;
    delete[] sendData_Yz;
    delete[] sendData_YZ;
    // Free recvData
    delete[] recvData_x;
    delete[] recvData_y;
    delete[] recvData_z;
    delete[] recvData_X;
    delete[] recvData_Y;
    delete[] recvData_Z;
    delete[] recvData_xy;
    delete[] recvData_xY;
    delete[] recvData_Xy;
    delete[] recvData_XY;
    delete[] recvData_xz;
    delete[] recvData_xZ;
    delete[] recvData_Xz;
    delete[] recvData_XZ;
    delete[] recvData_yz;
    delete[] recvData_yZ;
    delete[] recvData_Yz;
    delete[] recvData_YZ;
}

// Ideally stuff below here should be moved somewhere else -- doesn't really belong here
void WriteCheckpoint(const char *FILENAME, const double *cDen,
                     const double *cfq, size_t Np) {
    double value;
    ofstream File(FILENAME, ios::binary);
    for (size_t n = 0; n < Np; n++) {
        // Write the two density values
        value = cDen[n];
        File.write((char *)&value, sizeof(value));
        value = cDen[Np + n];
        File.write((char *)&value, sizeof(value));
        // Write the even distributions
        for (size_t q = 0; q < 19; q++) {
            value = cfq[q * Np + n];
            File.write((char *)&value, sizeof(value));
        }
    }
    File.close();
}

void ReadCheckpoint(char *FILENAME, double *cPhi, double *cfq, size_t Np) {
    double value = 0;
    ifstream File(FILENAME, ios::binary);
    for (size_t n = 0; n < Np; n++) {
        File.read((char *)&value, sizeof(value));
        cPhi[n] = value;
        // Read the distributions
        for (size_t q = 0; q < 19; q++) {
            File.read((char *)&value, sizeof(value));
            cfq[q * Np + n] = value;
        }
    }
    File.close();
}

void ReadBinaryFile(char *FILENAME, double *Data, size_t N) {
    double value;
    ifstream File(FILENAME, ios::binary);
    if (File.good()) {
        for (size_t n = 0; n < N; n++) {
            // Write the two density values
            File.read((char *)&value, sizeof(value));
            Data[n] = value;
        }
    } else {
        for (size_t n = 0; n < N; n++)
            Data[n] = 1.2e-34;
    }
    File.close();
}

void Domain::ReadFromFile(const std::string &Filename,
                          const std::string &Datatype, double *UserData) {
    //........................................................................................
    // Reading the user-defined input file
    // NOTE: so far it only supports BC=0 (periodic) and BC=5 (mixed reflection)
    //       because if checkerboard or inlet/outlet buffer layers are added, the
    //       value of the void space is undefined.
    // NOTE: if BC=5 is used, where the inlet and outlet layers of the domain are modified,
    //       user needs to modify the input file accordingly before LBPM simulator read
    //       the input file.
    //........................................................................................
    int RANK = rank();
    int nprocs, nprocx, nprocy, nprocz, nx, ny, nz;
    int64_t global_Nx, global_Ny, global_Nz;
    int64_t i, j, k;
    //TODO These offset we may still need them
    int64_t xStart, yStart, zStart;
    xStart = yStart = zStart = 0;

    // Read domain parameters
    // TODO currently the size of the data is still read from Domain{};
    //      but user may have a user-specified size
    auto size = database->getVector<int>("n");
    auto SIZE = database->getVector<int>("N");
    auto nproc = database->getVector<int>("nproc");
    //TODO currently the funcationality "offset" is disabled as the user-defined input data may have a different size from that of the input domain
    if (database->keyExists("offset")) {
        auto offset = database->getVector<int>("offset");
        xStart = offset[0];
        yStart = offset[1];
        zStart = offset[2];
    }

    nx = size[0];
    ny = size[1];
    nz = size[2];
    nprocx = nproc[0];
    nprocy = nproc[1];
    nprocz = nproc[2];
    global_Nx = SIZE[0];
    global_Ny = SIZE[1];
    global_Nz = SIZE[2];
    nprocs = nprocx * nprocy * nprocz;

    double *SegData = NULL;
    if (RANK == 0) {
        printf("User-defined input file: %s (data type: %s)\n",
               Filename.c_str(), Datatype.c_str());
        printf("NOTE: currently only BC=0 or 5 supports user-defined input "
               "file!\n");
        // Rank=0 reads the entire segmented data and distributes to worker processes
        printf("Dimensions of the user-defined input file: %ld x %ld x %ld \n",
               global_Nx, global_Ny, global_Nz);
        int64_t SIZE = global_Nx * global_Ny * global_Nz;

        if (Datatype == "double") {
            printf("Reading input data as double precision floating number\n");
            SegData = new double[SIZE];
            FILE *SEGDAT = fopen(Filename.c_str(), "rb");
            if (SEGDAT == NULL)
                ERROR("Domain.cpp: Error reading user-defined file!\n");
            size_t ReadSeg;
            ReadSeg = fread(SegData, 8, SIZE, SEGDAT);
            if (ReadSeg != size_t(SIZE))
                printf("Domain.cpp: Error reading file: %s\n",
                       Filename.c_str());
            fclose(SEGDAT);
        } else {
            ERROR("Error: User-defined input file only supports "
                  "double-precision floating number!\n");
        }
        printf("Read file successfully from %s \n", Filename.c_str());
    }

    // Get the rank info
    int64_t N = (nx + 2) * (ny + 2) * (nz + 2);

    // number of sites to use for periodic boundary condition transition zone
    //int64_t z_transition_size = (nprocz*nz - (global_Nz - zStart))/2;
    //if (z_transition_size < 0) z_transition_size=0;
    int64_t z_transition_size = 0;

    //char LocalRankFilename[1000];//just for debug
    double *loc_id;
    loc_id = new double[(nx + 2) * (ny + 2) * (nz + 2)];

    // Set up the sub-domains
    if (RANK == 0) {
        printf("Decomposing user-defined input file\n");
        printf("Distributing subdomains across %i processors \n", nprocs);
        printf("Process grid: %i x %i x %i \n", nprocx, nprocy, nprocz);
        printf("Subdomain size: %i x %i x %i \n", nx, ny, nz);
        printf("Size of transition region: %ld \n", z_transition_size);

        for (int kp = 0; kp < nprocz; kp++) {
            for (int jp = 0; jp < nprocy; jp++) {
                for (int ip = 0; ip < nprocx; ip++) {
                    // rank of the process that gets this subdomain
                    int rnk = kp * nprocx * nprocy + jp * nprocx + ip;
                    // Pack and send the subdomain for rnk
                    for (k = 0; k < nz + 2; k++) {
                        for (j = 0; j < ny + 2; j++) {
                            for (i = 0; i < nx + 2; i++) {
                                int64_t x = xStart + ip * nx + i - 1;
                                int64_t y = yStart + jp * ny + j - 1;
                                // int64_t z = zStart + kp*nz + k-1;
                                int64_t z = zStart + kp * nz + k - 1 -
                                            z_transition_size;
                                if (x < xStart)
                                    x = xStart;
                                if (!(x < global_Nx))
                                    x = global_Nx - 1;
                                if (y < yStart)
                                    y = yStart;
                                if (!(y < global_Ny))
                                    y = global_Ny - 1;
                                if (z < zStart)
                                    z = zStart;
                                if (!(z < global_Nz))
                                    z = global_Nz - 1;
                                int64_t nlocal =
                                    k * (nx + 2) * (ny + 2) + j * (nx + 2) + i;
                                int64_t nglobal = z * global_Nx * global_Ny +
                                                  y * global_Nx + x;
                                loc_id[nlocal] = SegData[nglobal];
                            }
                        }
                    }
                    if (rnk == 0) {
                        for (k = 0; k < nz + 2; k++) {
                            for (j = 0; j < ny + 2; j++) {
                                for (i = 0; i < nx + 2; i++) {
                                    int nlocal = k * (nx + 2) * (ny + 2) +
                                                 j * (nx + 2) + i;
                                    UserData[nlocal] = loc_id[nlocal];
                                }
                            }
                        }
                    } else {
                        //printf("Sending data to process %i \n", rnk);
                        Comm.send(loc_id, N, rnk, 15);
                    }
                }
            }
        }

    } else {
        // Recieve the subdomain from rank = 0
        //printf("Ready to recieve data %i at process %i \n", N,rank);
        Comm.recv(UserData, N, 0, 15);
    }
    Comm.barrier();
}

void Domain::AggregateLabels(const std::string &filename,
                             DoubleArray &UserData) {

    int nx = Nx;
    int ny = Ny;
    int nz = Nz;

    int npx = nprocx();
    int npy = nprocy();
    int npz = nprocz();

    int ipx = iproc();
    int ipy = jproc();
    int ipz = kproc();

    int nprocs = nprocx() * nprocy() * nprocz();

    int full_nx = npx * (nx - 2);
    int full_ny = npy * (ny - 2);
    int full_nz = npz * (nz - 2);
    int local_size = (nx - 2) * (ny - 2) * (nz - 2);
    unsigned long int full_size = long(full_nx) * long(full_ny) * long(full_nz);

    double *LocalID;
    LocalID = new double[local_size];

    //printf("aggregate labels: local size=%i, global size = %i",local_size, full_size);
    // assign the ID for the local sub-region
    for (int k = 1; k < nz - 1; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                double local_id_val = UserData(i, j, k);
                LocalID[(k - 1) * (nx - 2) * (ny - 2) + (j - 1) * (nx - 2) + i -
                        1] = local_id_val;
            }
        }
    }
    Comm.barrier();

    // populate the FullID
    if (rank() == 0) {
        double *FullID;
        FullID = new double[full_size];
        // first handle local ID for rank 0
        for (int k = 1; k < nz - 1; k++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    int x = i - 1;
                    int y = j - 1;
                    int z = k - 1;
                    int n_local = (k - 1) * (nx - 2) * (ny - 2) +
                                  (j - 1) * (nx - 2) + i - 1;
                    unsigned long int n_full =
                        z * long(full_nx) * long(full_ny) + y * long(full_nx) +
                        x;
                    FullID[n_full] = LocalID[n_local];
                }
            }
        }
        // next get the local ID from the other ranks
        for (int rnk = 1; rnk < nprocs; rnk++) {
            ipz = rnk / (npx * npy);
            ipy = (rnk - ipz * npx * npy) / npx;
            ipx = (rnk - ipz * npx * npy - ipy * npx);
            //printf("ipx=%i ipy=%i ipz=%i\n", ipx, ipy, ipz);
            int tag = 15 + rnk;
            //MPI_Recv(LocalID,local_size,MPI_DOUBLE,rnk,tag,Comm,MPI_STATUS_IGNORE);
            Comm.recv(LocalID, local_size, rnk, tag);

            for (int k = 1; k < nz - 1; k++) {
                for (int j = 1; j < ny - 1; j++) {
                    for (int i = 1; i < nx - 1; i++) {
                        int x = i - 1 + ipx * (nx - 2);
                        int y = j - 1 + ipy * (ny - 2);
                        int z = k - 1 + ipz * (nz - 2);
                        int n_local = (k - 1) * (nx - 2) * (ny - 2) +
                                      (j - 1) * (nx - 2) + i - 1;
                        unsigned long int n_full =
                            z * long(full_nx) * long(full_ny) +
                            y * long(full_nx) + x;
                        FullID[n_full] = LocalID[n_local];
                    }
                }
            }
        }
        // write the output
        FILE *OUTFILE = fopen(filename.c_str(), "wb");
        fwrite(FullID, 8, full_size, OUTFILE);
        fclose(OUTFILE);
    } else {
        // send LocalID to rank=0
        int tag = 15 + rank();
        int dstrank = 0;
        //MPI_Send(LocalID,local_size,MPI_DOUBLE,dstrank,tag,Comm);
        Comm.send(LocalID, local_size, dstrank, tag);
    }
    Comm.barrier();
}
