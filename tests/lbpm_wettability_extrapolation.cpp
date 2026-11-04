#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <queue>
#include <cmath>
#include <algorithm>
#include "common/Array.h"
#include "common/Domain.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"

//*************************************************************************
// Wettability Extrapolation
//   Given some points with known contact angle, this algorithm extrapolates 
//   values to solid-fluid interface voxels ("pore") OR into the grains ("grain").
//*************************************************************************

using namespace std;

static int wrapPhysical(int index, int extentWithHalo)
{
    const int lastPhysical = extentWithHalo - 2;
    if (index == 0)
        return lastPhysical;
    if (index == extentWithHalo - 1)
        return 1;
    return index;
}

int main(int argc, char **argv)
{
    Utilities::startup( argc, argv );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    {
        string filename;
        if (argc > 1){
            filename = argv[1];
        } else {
            ERROR("No input database provided\n");
        }

        auto db = std::make_shared<Database>( filename );
        auto domain_db = db->getDatabase( "Domain" );
        auto READFILE = domain_db->getScalar<std::string>( "Filename" );
        auto BC = domain_db->getScalar<int>( "BC" );
        
        std::string protocol = domain_db->getScalar<std::string>( "protocol" );
        if (protocol != "pore" && protocol != "grain") {
            ERROR("Unknown protocol '" + protocol + "'. Expected 'pore' or 'grain'.");
        }

        bool periodic = (BC == 0);

        if (rank == 0) {
            printf("BC = %d -> using %s boundaries\n", BC, periodic ? "periodic" : "open");
            printf("Protocol: %s\n", protocol.c_str());
        }

        std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
        std::shared_ptr<Domain> Mask (new Domain(domain_db,comm));

        Mask->Decomp(READFILE);
        Mask->CommInit();

        Dm->Decomp(READFILE);

        for (size_t n = 0; n < Dm->id.size(); n++) {
            Dm->id[n] = 1;
        }
        Dm->CommInit();

        int nx = Mask->Nx;
        int ny = Mask->Ny;
        int nz = Mask->Nz;

        Array<char> id_solid(nx, ny, nz);
        DoubleArray SignDist(nx, ny, nz);
        DoubleArray FinalMap(nx, ny, nz);

        std::queue<int64_t> Q;

        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    int n = k * nx * ny + j * nx + i;

                    FinalMap(i, j, k) = 0.0;

                    if (Mask->id[n] > 0) {
                        id_solid(i, j, k) = 1; 
                    } else {
                        id_solid(i, j, k) = 0;
                        if (Mask->id[n] < 0 &&
                            (!periodic || (i > 0 && i < nx - 1 &&
                                           j > 0 && j < ny - 1 &&
                                           k > 0 && k < nz - 1))) {
                            Q.push(n);
                            FinalMap(i, j, k) = (double)Mask->id[n];
                        }
                    }
                }
            }
        }

        if (rank == 0) printf("Initialized solid phase.\n");

        if(protocol == "pore") {
            if (periodic) {
                for (int k = 0; k < nz; k++) {
                    for (int j = 0; j < ny; j++) {
                        for (int i = 0; i < nx; i++) {
                            SignDist(i, j, k) = 2.0 * double(id_solid(i, j, k)) - 1.0;
                        }
                    }
                }
                if (rank == 0) printf("Converting to Signed Distance function\n");
                CalcDist(SignDist, id_solid, *Dm);
                
            } else {
                int pad = 1; 
                
                int global_nx = domain_db->getVector<int>("N")[0];
                int global_ny = domain_db->getVector<int>("N")[1];
                int global_nz = domain_db->getVector<int>("N")[2];
                
                int nx_p = global_nx + 2 * pad;
                int ny_p = global_ny + 2 * pad;
                int nz_p = global_nz + 2 * pad;

                std::shared_ptr<Domain> DmPad(new Domain(nx_p, ny_p, nz_p, rank, 1, 1, 1, (double)nx_p, (double)ny_p, (double)nz_p, 0));
                
                for (size_t n = 0; n < DmPad->id.size(); n++) {
                    DmPad->id[n] = 1;
                }
                DmPad->CommInit();

                int NXP = DmPad->Nx;
                int NYP = DmPad->Ny;
                int NZP = DmPad->Nz;

                Array<char> id_solid_pad(NXP, NYP, NZP);
                DoubleArray SignDistPad(NXP, NYP, NZP);

                for (int k = 0; k < NZP; k++) {
                    for (int j = 0; j < NYP; j++) {
                        for (int i = 0; i < NXP; i++) {
                            int orig_i = std::max(0, std::min(nx - 1, i - pad));
                            int orig_j = std::max(0, std::min(ny - 1, j - pad));
                            int orig_k = std::max(0, std::min(nz - 1, k - pad));
                            
                            id_solid_pad(i, j, k) = id_solid(orig_i, orig_j, orig_k);
                            
                            if(i == 0 || j == 0 || k == 0 || i == NXP-1 || j == NYP-1 || k == NZP-1) {
                                id_solid_pad(i, j, k) = 0;
                            }
                            
                            SignDistPad(i, j, k)  = 2.0 * double(id_solid_pad(i, j, k)) - 1.0;
                        }
                    }
                }

                if (rank == 0) printf("Converting to Signed Distance function\n");
                CalcDist(SignDistPad, id_solid_pad, *DmPad);

                for (int k = 0; k < nz; k++) {
                    for (int j = 0; j < ny; j++) {
                        for (int i = 0; i < nx; i++) {
                            SignDist(i, j, k) = SignDistPad(i + pad, j + pad, k + pad);
                        }
                    }
                }

                if (rank == 0) printf("Cropped SignDist back to original size.\n");
            }

            if (rank == 0) printf("Extrapolation algorithm running (pore).\n");

            while (!Q.empty()) {
                int64_t n = Q.front();
                Q.pop();

                int i = n % nx;
                int j = (n / nx) % ny;
                int k = n / (nx * ny);

                double current_val = FinalMap(i, j, k);

                for (int dk = -1; dk <= 1; dk++) {
                    for (int dj = -1; dj <= 1; dj++) {
                        for (int di = -1; di <= 1; di++) {
                            if (di == 0 && dj == 0 && dk == 0) continue;

                            int ni = i + di;
                            int nj = j + dj;
                            int nk = k + dk;

                            if (periodic) {
                                ni = wrapPhysical(ni, nx);
                                nj = wrapPhysical(nj, ny);
                                nk = wrapPhysical(nk, nz);
                            } else {
                                if (ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz)
                                    continue;
                            }

                            if (std::abs(SignDist(ni, nj, nk) - (-0.5)) < 1e-3 && FinalMap(ni, nj, nk) == 0.0) {
                                FinalMap(ni, nj, nk) = current_val;
                                int64_t next_n = (int64_t)nk * nx * ny + (int64_t)nj * nx + ni;
                                Q.push(next_n);
                            }
                        }
                    }
                }
            }
        }
        
        else if(protocol == "grain") {
            if (periodic) {
                for (int k = 0; k < nz; k++) {
                    for (int j = 0; j < ny; j++) {
                        for (int i = 0; i < nx; i++) {
                            SignDist(i, j, k) = 2.0 * double(id_solid(i, j, k)) - 1.0;
                        }
                    }
                }
                
                if (rank == 0) printf("Converting to Signed Distance function\n");
                CalcDist(SignDist, id_solid, *Dm);
                
                if (rank == 0) printf("Extrapolation algorithm running (grain - periodic).\n");

                while (!Q.empty()) {
                    int64_t n = Q.front();
                    Q.pop();

                    int i = n % nx;
                    int j = (n / nx) % ny;
                    int k = n / (nx * ny);

                    double current_val = FinalMap(i, j, k);

                    for (int dk = -1; dk <= 1; dk++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            for (int di = -1; di <= 1; di++) {
                                if (di == 0 && dj == 0 && dk == 0) continue;

                                int ni = wrapPhysical(i + di, nx);
                                int nj = wrapPhysical(j + dj, ny);
                                int nk = wrapPhysical(k + dk, nz);

                                if (SignDist(ni, nj, nk) < 0 && FinalMap(ni, nj, nk) == 0.0) {
                                    FinalMap(ni, nj, nk) = current_val;
                                    int64_t next_n = (int64_t)nk * nx * ny + (int64_t)nj * nx + ni;
                                    Q.push(next_n);
                                }
                            }
                        }
                    }
                }
                
            } else {
                int pad = 1; 
                
                int global_nx = domain_db->getVector<int>("N")[0];
                int global_ny = domain_db->getVector<int>("N")[1];
                int global_nz = domain_db->getVector<int>("N")[2];
                
                int nx_p = global_nx + 2 * pad;
                int ny_p = global_ny + 2 * pad;
                int nz_p = global_nz + 2 * pad;

                std::shared_ptr<Domain> DmPad(new Domain(nx_p, ny_p, nz_p, rank, 1, 1, 1, (double)nx_p, (double)ny_p, (double)nz_p, 0));
                
                for (size_t n = 0; n < DmPad->id.size(); n++) {
                    DmPad->id[n] = 1;
                }
                DmPad->CommInit();

                int NXP = DmPad->Nx;
                int NYP = DmPad->Ny;
                int NZP = DmPad->Nz;

                Array<char> id_solid_pad(NXP, NYP, NZP);
                DoubleArray SignDistPad(NXP, NYP, NZP);
                DoubleArray FinalMapPad(NXP, NYP, NZP);
                std::queue<int64_t> Q_pad;

                for (int k = 0; k < NZP; k++) {
                    for (int j = 0; j < NYP; j++) {
                        for (int i = 0; i < NXP; i++) {
                            FinalMapPad(i, j, k) = 0.0;
                            
                            if(i == 0 || j == 0 || k == 0 || i == NXP-1 || j == NYP-1 || k == NZP-1) {
                                id_solid_pad(i, j, k) = 1; 
                            } else {
                                id_solid_pad(i, j, k) = id_solid(i - pad, j - pad, k - pad);
                                FinalMapPad(i, j, k)  = FinalMap(i - pad, j - pad, k - pad);
                                
                                if (FinalMapPad(i, j, k) != 0.0) {
                                    int64_t n_pad = (int64_t)k * NXP * NYP + (int64_t)j * NXP + i;
                                    Q_pad.push(n_pad);
                                }
                            }
                            SignDistPad(i, j, k) = 2.0 * double(id_solid_pad(i, j, k)) - 1.0;
                        }
                    }
                }

                if (rank == 0) printf("Converting to Signed Distance function (padded grain)\n");
                CalcDist(SignDistPad, id_solid_pad, *DmPad);

                if (rank == 0) printf("Extrapolation algorithm running (grain - open).\n");

                while (!Q_pad.empty()) {
                    int64_t n = Q_pad.front();
                    Q_pad.pop();

                    int i = n % NXP;
                    int j = (n / NXP) % NYP;
                    int k = n / (NXP * NYP);

                    double current_val = FinalMapPad(i, j, k);

                    for (int dk = -1; dk <= 1; dk++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            for (int di = -1; di <= 1; di++) {
                                if (di == 0 && dj == 0 && dk == 0) continue;

                                int ni = i + di;
                                int nj = j + dj;
                                int nk = k + dk;

                                if (ni < 0 || ni >= NXP || nj < 0 || nj >= NYP || nk < 0 || nk >= NZP)
                                    continue;

                                if (SignDistPad(ni, nj, nk) < 0 && FinalMapPad(ni, nj, nk) == 0.0) {
                                    FinalMapPad(ni, nj, nk) = current_val;
                                    int64_t next_n = (int64_t)nk * NXP * NYP + (int64_t)nj * NXP + ni;
                                    Q_pad.push(next_n);
                                }
                            }
                        }
                    }
                }

                for (int k = 0; k < nz; k++) {
                    for (int j = 0; j < ny; j++) {
                        for (int i = 0; i < nx; i++) {
                            FinalMap(i, j, k) = FinalMapPad(i + pad, j + pad, k + pad);
                        }
                    }
                }
                if (rank == 0) printf("Cropped FinalMap back to original size.\n");
            }
        }

        for (int k = 0; k < nz; k++){
            for (int j = 0; j < ny; j++){
                for (int i = 0; i < nx; i++){
                    int n = k * nx * ny + j * nx + i;
                    if (FinalMap(i, j, k) != 0.0) {
                        Mask->id[n] = (signed char)FinalMap(i, j, k);
                    }
                }
            }
        }

        std::string filename2 = READFILE;
        std::string extension = ".raw";
        size_t pos = filename2.rfind(extension);
        if (pos != std::string::npos && pos == (filename2.length() - extension.length())) {
            filename2 = filename2.substr(0, pos);
        }
        filename2 += "_wettability.raw";

        if (rank == 0) printf("Writing file to: %s \n", filename2.c_str());

        Mask->AggregateLabels( filename2 );
    }

    Utilities::shutdown();
    return 0;
}
