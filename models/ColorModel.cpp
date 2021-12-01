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
color lattice boltzmann model
 */

#include "models/ColorModel.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"
#include "common/Communication.h"
#include "common/ReadMicroCT.h"
#include <stdlib.h>
#include <time.h>

ScaLBL_ColorModel::ScaLBL_ColorModel(int RANK, int NP,
                                     const Utilities::MPI &COMM)
    : rank(RANK), nprocs(NP), Restart(0), timestep(0), timestepMax(0), tauA(0),
      tauB(0), rhoA(0), rhoB(0), alpha(0), beta(0), Fx(0), Fy(0), Fz(0),
      flux(0), din(0), dout(0), inletA(0), inletB(0), outletA(0), outletB(0),
      Nx(0), Ny(0), Nz(0), N(0), Np(0), nprocx(0), nprocy(0), nprocz(0),
      BoundaryCondition(0), Lx(0), Ly(0), Lz(0), id(nullptr),
      NeighborList(nullptr), dvcMap(nullptr), fq(nullptr), Aq(nullptr),
      Bq(nullptr), Den(nullptr), Phi(nullptr), ColorGrad(nullptr),
      Velocity(nullptr), Pressure(nullptr), comm(COMM) {
    REVERSE_FLOW_DIRECTION = false;
}
ScaLBL_ColorModel::~ScaLBL_ColorModel() {
    delete[] id;
    ScaLBL_FreeDeviceMemory(NeighborList);
    ScaLBL_FreeDeviceMemory(dvcMap);
    ScaLBL_FreeDeviceMemory(fq);
    ScaLBL_FreeDeviceMemory(Aq);
    ScaLBL_FreeDeviceMemory(Bq);
    ScaLBL_FreeDeviceMemory(Den);
    ScaLBL_FreeDeviceMemory(Phi);
    ScaLBL_FreeDeviceMemory(Pressure);
    ScaLBL_FreeDeviceMemory(Velocity);
    ScaLBL_FreeDeviceMemory(ColorGrad);
}


void ScaLBL_ColorModel::ReadParams(string filename) {
    // read the input database
    db = std::make_shared<Database>(filename);
    domain_db = db->getDatabase("Domain");
    color_db = db->getDatabase("Color");
    analysis_db = db->getDatabase("Analysis");
    vis_db = db->getDatabase("Visualization");

    // set defaults
    timestepMax = 100000;
    tauA = tauB = 1.0;
    rhoA = rhoB = 1.0;
    Fx = Fy = Fz = 0.0;
    alpha = 1e-3;
    beta = 0.95;
    Restart = false;
    din = dout = 1.0;
    flux = 0.0;

    // Color Model parameters
    if (color_db->keyExists("timestepMax")) {
        timestepMax = color_db->getScalar<int>("timestepMax");
    }
    if (color_db->keyExists("tauA")) {
        tauA = color_db->getScalar<double>("tauA");
    }
    if (color_db->keyExists("tauB")) {
        tauB = color_db->getScalar<double>("tauB");
    }
    if (color_db->keyExists("rhoA")) {
        rhoA = color_db->getScalar<double>("rhoA");
    }
    if (color_db->keyExists("rhoB")) {
        rhoB = color_db->getScalar<double>("rhoB");
    }
    if (color_db->keyExists("F")) {
        Fx = color_db->getVector<double>("F")[0];
        Fy = color_db->getVector<double>("F")[1];
        Fz = color_db->getVector<double>("F")[2];
    }
    if (color_db->keyExists("alpha")) {
        alpha = color_db->getScalar<double>("alpha");
    }
    if (color_db->keyExists("beta")) {
        beta = color_db->getScalar<double>("beta");
    }
    if (color_db->keyExists("Restart")) {
        Restart = color_db->getScalar<bool>("Restart");
    }
    if (color_db->keyExists("din")) {
        din = color_db->getScalar<double>("din");
    }
    if (color_db->keyExists("dout")) {
        dout = color_db->getScalar<double>("dout");
    }
    if (color_db->keyExists("flux")) {
        flux = color_db->getScalar<double>("flux");
    }
    inletA = 1.f;
    inletB = 0.f;
    outletA = 0.f;
    outletB = 1.f;


    BoundaryCondition = 0;
    if (color_db->keyExists("BC")) {
        BoundaryCondition = color_db->getScalar<int>("BC");
    } else if (domain_db->keyExists("BC")) {
        BoundaryCondition = domain_db->getScalar<int>("BC");
    }
    if (domain_db->keyExists("InletLayersPhase")) {
        int inlet_layers_phase = domain_db->getScalar<int>("InletLayersPhase");
        if (inlet_layers_phase == 2) {
            inletA = 0.0;
            inletB = 1.0;
        }
    }
    if (domain_db->keyExists("OutletLayersPhase")) {
        int outlet_layers_phase =
            domain_db->getScalar<int>("OutletLayersPhase");
        if (outlet_layers_phase == 1) {
            inletA = 1.0;
            inletB = 0.0;
        }
    }

    // Override user-specified boundary condition for specific protocols
    auto protocol = color_db->getWithDefault<std::string>("protocol", "none");
    if (protocol == "seed water") {
        if (BoundaryCondition != 0 && BoundaryCondition != 5) {
            BoundaryCondition = 0;
            if (rank == 0)
                printf("WARNING: protocol (seed water) supports only full "
                       "periodic boundary condition \n");
        }
        domain_db->putScalar<int>("BC", BoundaryCondition);
    } else if (protocol == "open connected oil") {
        if (BoundaryCondition != 0 && BoundaryCondition != 5) {
            BoundaryCondition = 0;
            if (rank == 0)
                printf("WARNING: protocol (open connected oil) supports only "
                       "full periodic boundary condition \n");
        }
        domain_db->putScalar<int>("BC", BoundaryCondition);
    } else if (protocol == "shell aggregation") {
        if (BoundaryCondition != 0 && BoundaryCondition != 5) {
            BoundaryCondition = 0;
            if (rank == 0)
                printf("WARNING: protocol (shell aggregation) supports only "
                       "full periodic boundary condition \n");
        }
        domain_db->putScalar<int>("BC", BoundaryCondition);
    } else if (protocol == "fractional flow") {
        if (BoundaryCondition != 0 && BoundaryCondition != 5) {
            BoundaryCondition = 0;
            if (rank == 0)
                printf("WARNING: protocol (fractional flow) supports only full "
                       "periodic boundary condition \n");
        }
        domain_db->putScalar<int>("BC", BoundaryCondition);
    } else if (protocol == "centrifuge") {
        if (BoundaryCondition != 3) {
            BoundaryCondition = 3;
            if (rank == 0)
                printf("WARNING: protocol (centrifuge) supports only constant "
                       "pressure boundary condition \n");
        }
        domain_db->putScalar<int>("BC", BoundaryCondition);
    } else if (protocol == "core flooding") {
        if (rank == 0)
            printf("Using core flooding protocol \n");
        if (BoundaryCondition != 4) {
            BoundaryCondition = 4;
            if (rank == 0)
                printf("WARNING: protocol (core flooding) supports only "
                       "volumetric flux boundary condition \n");
        }
        domain_db->putScalar<int>("BC", BoundaryCondition);
    }
}


void ScaLBL_ColorModel::SetDomain() {
    Dm = std::shared_ptr<Domain>(
        new Domain(domain_db, comm)); // full domain for analysis
    Mask = std::shared_ptr<Domain>(
        new Domain(domain_db, comm)); // mask domain removes immobile phases
    // domain parameters
    Nx = Dm->Nx;
    Ny = Dm->Ny;
    Nz = Dm->Nz;
    Lx = Dm->Lx;
    Ly = Dm->Ly;
    Lz = Dm->Lz;
    N = Nx * Ny * Nz;
    id = new signed char[N];
    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = 1; // initialize this way
    //Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm) ); // TwoPhase analysis object
    Averages =
        std::shared_ptr<SubPhase>(new SubPhase(Dm)); // TwoPhase analysis object
    comm.barrier();
    Dm->CommInit();
    comm.barrier();
    // Read domain parameters
    rank = Dm->rank();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();
}

void ScaLBL_ColorModel::ReadInput() {

    sprintf(LocalRankString, "%05d", rank);
    sprintf(LocalRankFilename, "%s%s", "ID.", LocalRankString);
    sprintf(LocalRestartFile, "%s%s", "Restart.", LocalRankString);

    if (color_db->keyExists("image_sequence")) {
        auto ImageList = color_db->getVector<std::string>("image_sequence");
        int IMAGE_INDEX = color_db->getWithDefault<int>("image_index", 0);
        std::string first_image = ImageList[IMAGE_INDEX];
        Mask->Decomp(first_image);
        IMAGE_INDEX++;
    } else if (domain_db->keyExists("GridFile")) {
        // Read the local domain data
        auto input_id = readMicroCT(*domain_db, MPI_COMM_WORLD);
        // Fill the halo (assuming GCW of 1)
        array<int, 3> size0 = {(int)input_id.size(0), (int)input_id.size(1),
                               (int)input_id.size(2)};
        ArraySize size1 = {(size_t)Mask->Nx, (size_t)Mask->Ny,
                           (size_t)Mask->Nz};
        ASSERT((int)size1[0] == size0[0] + 2 && (int)size1[1] == size0[1] + 2 &&
               (int)size1[2] == size0[2] + 2);
        fillHalo<signed char> fill(MPI_COMM_WORLD, Mask->rank_info, size0,
                                   {1, 1, 1}, 0, 1);
        Array<signed char> id_view;
        id_view.viewRaw(size1, Mask->id.data());
        fill.copy(input_id, id_view);
        fill.fill(id_view);
    } else if (domain_db->keyExists("Filename")) {
        auto Filename = domain_db->getScalar<std::string>("Filename");
        Mask->Decomp(Filename);
    } else {
        Mask->ReadIDs();
    }
    for (int i = 0; i < Nx * Ny * Nz; i++)
        id[i] = Mask->id[i]; // save what was read

    // Generate the signed distance map
    // Initialize the domain and communication
    Array<char> id_solid(Nx, Ny, Nz);
    // Solve for the position of the solid phase
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                // Initialize the solid phase
                signed char label = Mask->id[n];
                if (label > 0)
                    id_solid(i, j, k) = 1;
                else
                    id_solid(i, j, k) = 0;
            }
        }
    }
    // Initialize the signed distance function
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                // Initialize distance to +/- 1
                Averages->SDs(i, j, k) = 2.0 * double(id_solid(i, j, k)) - 1.0;
            }
        }
    }
    //	MeanFilter(Averages->SDs);
    Minkowski Solid(Dm);
    if (rank == 0)
        printf("Initialized solid phase -- Converting to Signed Distance "
               "function \n");
    CalcDist(Averages->SDs, id_solid, *Mask);
    Solid.ComputeScalar(Averages->SDs, 0.0);
    /* save averages */
    Averages->solid.V = Solid.Vi;
    Averages->solid.A = Solid.Ai;
    Averages->solid.H = Solid.Ji;
    Averages->solid.X = Solid.Xi;
    Averages->gsolid.V = Solid.Vi_global;
    Averages->gsolid.A = Solid.Ai_global;
    Averages->gsolid.H = Solid.Ji_global;
    Averages->gsolid.X = Solid.Xi_global;
    /* write to file */
    if (rank == 0) {
        FILE *SOLID = fopen("solid.csv", "w");
        fprintf(SOLID, "Vs As Hs Xs\n");
        fprintf(SOLID, "%.8g %.8g %.8g %.8g\n", Solid.Vi_global,
                Solid.Ai_global, Solid.Ji_global, Solid.Xi_global);
        fclose(SOLID);
    }
    if (rank == 0)
        cout << "Domain set." << endl;

    Averages->SetParams(rhoA, rhoB, tauA, tauB, Fx, Fy, Fz, alpha, beta);
}

void ScaLBL_ColorModel::AssignComponentLabels(double *phase) {
    size_t NLABELS = 0;
    signed char VALUE = 0;
    double AFFINITY = 0.f;

    auto LabelList = color_db->getVector<int>("ComponentLabels");
    auto AffinityList = color_db->getVector<double>("ComponentAffinity");
    auto WettingConvention =
        color_db->getWithDefault<std::string>("WettingConvention", "none");

    NLABELS = LabelList.size();
    if (NLABELS != AffinityList.size()) {
        ERROR("Error: ComponentLabels and ComponentAffinity must be the same "
              "length! \n");
    }

    if (WettingConvention == "SCAL") {
        for (size_t idx = 0; idx < NLABELS; idx++)
            AffinityList[idx] *= -1.0;
    }

    double *label_count;
    double *label_count_global;
    label_count = new double[NLABELS];
    label_count_global = new double[NLABELS];
    // Assign the labels
    for (size_t idx = 0; idx < NLABELS; idx++)
        label_count[idx] = 0;

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                VALUE = id[n];
                // Assign the affinity from the paired list
                for (unsigned int idx = 0; idx < NLABELS; idx++) {
                    if (VALUE == LabelList[idx]) {
                        AFFINITY = AffinityList[idx];
                        label_count[idx] += 1.0;
                        idx = NLABELS;
                    }
                }
                // fluid labels are reserved
                if (VALUE == 1)
                    AFFINITY = 1.0;
                else if (VALUE == 2)
                    AFFINITY = -1.0;
                phase[n] = AFFINITY;
            }
        }
    }

    // Set Dm to match Mask
    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = Mask->id[i];

    for (size_t idx = 0; idx < NLABELS; idx++)
        label_count_global[idx] = Dm->Comm.sumReduce(label_count[idx]);

    if (rank == 0) {
        printf("Component labels: %lu \n", NLABELS);
        for (unsigned int idx = 0; idx < NLABELS; idx++) {
            VALUE = LabelList[idx];
            AFFINITY = AffinityList[idx];
            double volume_fraction =
                double(label_count_global[idx]) /
                double((Nx - 2) * (Ny - 2) * (Nz - 2) * nprocs);
            printf("   label=%d, affinity=%f, volume fraction==%f\n", VALUE,
                   AFFINITY, volume_fraction);
        }
    }
}

void ScaLBL_ColorModel::Create() {
    /*
     *  This function creates the variables needed to run a LBM 
    */
    //.........................................................
    // don't perform computations at the eight corners
    //id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
    //id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;

    //.........................................................
    // Initialize communication structures in averaging domain
    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = Mask->id[i];
    Mask->CommInit();
    Np = Mask->PoreCount();
    //...........................................................................
    if (rank == 0)
        printf("Create ScaLBL_Communicator \n");
    // Create a communicator for the device (will use optimized layout)
    // ScaLBL_Communicator ScaLBL_Comm(Mask); // original
    ScaLBL_Comm =
        std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));
    ScaLBL_Comm_Regular =
        std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));

    int Npad = (Np / 16 + 2) * 16;
    if (rank == 0)
        printf("Set up memory efficient layout, %i | %i | %i \n", Np, Npad, N);
    Map.resize(Nx, Ny, Nz);
    Map.fill(-2);
    auto neighborList = new int[18 * Npad];
    Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map, neighborList,
                                              Mask->id.data(), Np, 1);
    comm.barrier();

    //...........................................................................
    //                MAIN  VARIABLES ALLOCATED HERE
    //...........................................................................
    // LBM variables
    if (rank == 0)
        printf("Allocating distributions \n");
    //......................device distributions.................................
    dist_mem_size = Np * sizeof(double);
    neighborSize = 18 * (Np * sizeof(int));
    //...........................................................................
    ScaLBL_AllocateDeviceMemory((void **)&NeighborList, neighborSize);
    ScaLBL_AllocateDeviceMemory((void **)&dvcMap, sizeof(int) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&fq, 19 * dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **)&Aq, 7 * dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **)&Bq, 7 * dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **)&Den, 2 * dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **)&Phi, sizeof(double) * Nx * Ny * Nz);
    ScaLBL_AllocateDeviceMemory((void **)&Pressure, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&Velocity, 3 * sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&ColorGrad, 3 * sizeof(double) * Np);
    //...........................................................................
    // Update GPU data structures
    if (rank == 0)
        printf("Setting up device map and neighbor list \n");
    fflush(stdout);
    int *TmpMap;
    TmpMap = new int[Np];
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int idx = Map(i, j, k);
                if (!(idx < 0))
                    TmpMap[idx] = k * Nx * Ny + j * Nx + i;
            }
        }
    }
    // check that TmpMap is valid
    for (int idx = 0; idx < ScaLBL_Comm->LastExterior(); idx++) {
        auto n = TmpMap[idx];
        if (n > Nx * Ny * Nz) {
            printf("Bad value! idx=%i \n", n);
            TmpMap[idx] = Nx * Ny * Nz - 1;
        }
    }
    for (int idx = ScaLBL_Comm->FirstInterior();
         idx < ScaLBL_Comm->LastInterior(); idx++) {
        auto n = TmpMap[idx];
        if (n > Nx * Ny * Nz) {
            printf("Bad value! idx=%i \n", n);
            TmpMap[idx] = Nx * Ny * Nz - 1;
        }
    }
    ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int) * Np);
    ScaLBL_Comm->Barrier();
    delete[] TmpMap;

    // copy the neighbor list
    ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
    delete[] neighborList;
    // initialize phi based on PhaseLabel (include solid component labels)
    double *PhaseLabel;
    PhaseLabel = new double[N];
    AssignComponentLabels(PhaseLabel);
    ScaLBL_CopyToDevice(Phi, PhaseLabel, N * sizeof(double));
    delete[] PhaseLabel;
}

/********************************************************
 * AssignComponentLabels                                 *
 ********************************************************/

void ScaLBL_ColorModel::Initialize() {

    /* if both capillary number and flux BC are specified */
    if (color_db->keyExists("capillary_number") && BoundaryCondition == 4) {
        double capillary_number =
            color_db->getScalar<double>("capillary_number");
        if (rank == 0)
            printf("   set flux to achieve Ca=%f \n", capillary_number);
        double MuB = rhoB * (tauB - 0.5) / 3.0;
        double IFT = 6.0 * alpha;
        double CrossSectionalArea =
            (double)(nprocx * (Nx - 2) * nprocy * (Ny - 2));
        flux = Mask->Porosity() * CrossSectionalArea * (Ny - 2) * IFT *
               capillary_number / MuB;
        if (rank == 0)
            printf("   flux=%f \n", flux);
    }
    color_db->putScalar<double>("flux", flux);

    if (rank == 0)
        printf("Initializing distributions \n");
    ScaLBL_D3Q19_Init(fq, Np);
    /*
     * This function initializes model
     */
    if (Restart == true) {
        if (rank == 0) {
            printf("Reading restart file! \n");
        }

        // Read in the restart file to CPU buffers
        int *TmpMap;
        TmpMap = new int[Np];

        double *cPhi, *cDist, *cDen;
        cPhi = new double[N];
        cDen = new double[2 * Np];
        cDist = new double[19 * Np];
        ScaLBL_CopyToHost(TmpMap, dvcMap, Np * sizeof(int));
        ScaLBL_CopyToHost(cPhi, Phi, N * sizeof(double));

        ifstream File(LocalRestartFile, ios::binary);
        int idx;
        double value, va, vb;
        for (int n = 0; n < Np; n++) {
            File.read((char *)&va, sizeof(va));
            File.read((char *)&vb, sizeof(vb));
            cDen[n] = va;
            cDen[Np + n] = vb;
        }
        for (int n = 0; n < Np; n++) {
            // Read the distributions
            for (int q = 0; q < 19; q++) {
                File.read((char *)&value, sizeof(value));
                cDist[q * Np + n] = value;
            }
        }
        File.close();

        for (int n = 0; n < ScaLBL_Comm->LastExterior(); n++) {
            va = cDen[n];
            vb = cDen[Np + n];
            value = (va - vb) / (va + vb);
            idx = TmpMap[n];
            if (!(idx < 0) && idx < N)
                cPhi[idx] = value;
        }
        for (int n = ScaLBL_Comm->FirstInterior();
             n < ScaLBL_Comm->LastInterior(); n++) {
            va = cDen[n];
            vb = cDen[Np + n];
            value = (va - vb) / (va + vb);
            idx = TmpMap[n];
            if (!(idx < 0) && idx < N)
                cPhi[idx] = value;
        }

        // Copy the restart data to the GPU
        ScaLBL_CopyToDevice(Den, cDen, 2 * Np * sizeof(double));
        ScaLBL_CopyToDevice(fq, cDist, 19 * Np * sizeof(double));
        ScaLBL_CopyToDevice(Phi, cPhi, N * sizeof(double));
        ScaLBL_Comm->Barrier();

        comm.barrier();
    }

    if (rank == 0)
        printf("Initializing phase field \n");
    ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0,
                           ScaLBL_Comm->LastExterior(), Np);
    ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq,
                           ScaLBL_Comm->FirstInterior(),
                           ScaLBL_Comm->LastInterior(), Np);

    // establish reservoirs for external bC
    if (BoundaryCondition == 1 || BoundaryCondition == 2 ||
        BoundaryCondition == 3 || BoundaryCondition == 4) {
        if (Dm->kproc() == 0) {
            ScaLBL_SetSlice_z(Phi, 1.0, Nx, Ny, Nz, 0);
            ScaLBL_SetSlice_z(Phi, 1.0, Nx, Ny, Nz, 1);
            ScaLBL_SetSlice_z(Phi, 1.0, Nx, Ny, Nz, 2);
        }
        if (Dm->kproc() == nprocz - 1) {
            ScaLBL_SetSlice_z(Phi, -1.0, Nx, Ny, Nz, Nz - 1);
            ScaLBL_SetSlice_z(Phi, -1.0, Nx, Ny, Nz, Nz - 2);
            ScaLBL_SetSlice_z(Phi, -1.0, Nx, Ny, Nz, Nz - 3);
        }
    }
    ScaLBL_CopyToHost(Averages->Phi.data(), Phi, N * sizeof(double));
}

double ScaLBL_ColorModel::Run(int returntime) {
    int nprocs = nprocx * nprocy * nprocz;
    const RankInfoStruct rank_info(rank, nprocx, nprocy, nprocz);
    //************ MAIN ITERATION LOOP ***************************************/
    comm.barrier();
    PROFILE_START("Loop");
    bool Regular = false;
    bool RESCALE_FORCE = false;
    bool SET_CAPILLARY_NUMBER = false;
    bool TRIGGER_FORCE_RESCALE = false;
    double tolerance = 0.01;
    auto WettingConvention = color_db->getWithDefault<std::string>( "WettingConvention", "none" );
    auto current_db = db->cloneDatabase();
    auto flow_db = db->getDatabase("FlowAdaptor");
    int MIN_STEADY_TIMESTEPS =
        flow_db->getWithDefault<int>("min_steady_timesteps", 1000000);
    int MAX_STEADY_TIMESTEPS =
        flow_db->getWithDefault<int>("max_steady_timesteps", 1000000);
    int RESCALE_FORCE_AFTER_TIMESTEP = MAX_STEADY_TIMESTEPS * 2;
    int INITIAL_TIMESTEP = timestep;

    double capillary_number = 1.0e-5;
    double Ca_previous = 0.0;
    double minCa = 8.0e-6;
    double maxCa = 1.0;
    if (color_db->keyExists("capillary_number")) {
        capillary_number = color_db->getScalar<double>("capillary_number");
        SET_CAPILLARY_NUMBER = true;
        maxCa = 2.0 * capillary_number;
        minCa = 0.8 * capillary_number;
    }
    if (color_db->keyExists("rescale_force_after_timestep")) {
        RESCALE_FORCE_AFTER_TIMESTEP =
            color_db->getScalar<int>("rescale_force_after_timestep");
        RESCALE_FORCE = true;
    }
    if (analysis_db->keyExists("tolerance")) {
        tolerance = analysis_db->getScalar<double>("tolerance");
    }
    
    runAnalysis analysis(current_db, rank_info, ScaLBL_Comm, Dm, Np, Regular,
                         Map);
    auto t1 = std::chrono::system_clock::now();
    int CURRENT_TIMESTEP = 0;
    int EXIT_TIMESTEP = min(timestepMax, returntime);
    while (timestep < EXIT_TIMESTEP) {
      PROFILE_START("Update");
        // *************ODD TIMESTEP*************
        timestep++;
        // Compute the Phase indicator field
        // Read for Aq, Bq happens in this routine (requires communication)
        ScaLBL_Comm->BiSendD3Q7AA(Aq, Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi,
                                     ScaLBL_Comm->FirstInterior(),
                                     ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->BiRecvD3Q7AA(Aq, Bq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0,
                                     ScaLBL_Comm->LastExterior(), Np);

        // Perform the collision operation
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        if (BoundaryCondition > 0 && BoundaryCondition < 5) {
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
        }
        // Halo exchange for phase field
        ScaLBL_Comm_Regular->SendHalo(Phi);

        ScaLBL_D3Q19_AAodd_Color(
            NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB,
            tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx, Nx * Ny,
            ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm_Regular->RecvHalo(Phi);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        // Set BCs
        if (BoundaryCondition == 3) {
            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }
        if (BoundaryCondition == 4) {
            din =
                ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 5) {
            ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
            ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
        }
        ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi,
                                 Velocity, rhoA, rhoB, tauA, tauB, alpha, beta,
                                 Fx, Fy, Fz, Nx, Nx * Ny, 0,
                                 ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->Barrier();

        // *************EVEN TIMESTEP*************
        timestep++;
        // Compute the Phase indicator field
        ScaLBL_Comm->BiSendD3Q7AA(Aq, Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi,
                                      ScaLBL_Comm->FirstInterior(),
                                      ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->BiRecvD3Q7AA(Aq, Bq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0,
                                      ScaLBL_Comm->LastExterior(), Np);

        // Perform the collision operation
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        // Halo exchange for phase field
        if (BoundaryCondition > 0 && BoundaryCondition < 5) {
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
        }
        ScaLBL_Comm_Regular->SendHalo(Phi);
        ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA,
                                  rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx,
                                  Nx * Ny, ScaLBL_Comm->FirstInterior(),
                                  ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm_Regular->RecvHalo(Phi);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        // Set boundary conditions
        if (BoundaryCondition == 3) {
            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 4) {
            din =
                ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 5) {
            ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
            ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
        }
        ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA,
                                  rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx,
                                  Nx * Ny, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->Barrier();
        //************************************************************************
        analysis.basic(
            timestep, current_db, *Averages, Phi, Pressure, Velocity, fq,
            Den); // allow initial ramp-up to get closer to steady state

        CURRENT_TIMESTEP += 2;
        if (CURRENT_TIMESTEP > MIN_STEADY_TIMESTEPS && BoundaryCondition == 0) {
            analysis.finish();

            double volB = Averages->gwb.V;
            double volA = Averages->gnb.V;
            volA /= Dm->Volume;
            volB /= Dm->Volume;
            ;
            //initial_volume = volA*Dm->Volume;
            double vA_x = Averages->gnb.Px / Averages->gnb.M;
            double vA_y = Averages->gnb.Py / Averages->gnb.M;
            double vA_z = Averages->gnb.Pz / Averages->gnb.M;
            double vB_x = Averages->gwb.Px / Averages->gwb.M;
            double vB_y = Averages->gwb.Py / Averages->gwb.M;
            double vB_z = Averages->gwb.Pz / Averages->gwb.M;
            double muA = rhoA * (tauA - 0.5) / 3.f;
            double muB = rhoB * (tauB - 0.5) / 3.f;
            double force_mag = sqrt(Fx * Fx + Fy * Fy + Fz * Fz);
            double dir_x = Fx / force_mag;
            double dir_y = Fy / force_mag;
            double dir_z = Fz / force_mag;
            if (force_mag == 0.0) {
                // default to z direction
                dir_x = 0.0;
                dir_y = 0.0;
                dir_z = 1.0;
                force_mag = 1.0;
            }
            double current_saturation = volB / (volA + volB);
            double flow_rate_A =
                volA * (vA_x * dir_x + vA_y * dir_y + vA_z * dir_z);
            double flow_rate_B =
                volB * (vB_x * dir_x + vB_y * dir_y + vB_z * dir_z);
            double Ca =
                fabs(muA * flow_rate_A + muB * flow_rate_B) / (5.796 * alpha);

            bool isSteady = false;
            if ((fabs((Ca - Ca_previous) / Ca) < tolerance &&
                 CURRENT_TIMESTEP > MIN_STEADY_TIMESTEPS))
                isSteady = true;
            if (CURRENT_TIMESTEP >= MAX_STEADY_TIMESTEPS)
                isSteady = true;

            if (isSteady && (Ca > maxCa || Ca < minCa) &&
                SET_CAPILLARY_NUMBER) {
                /* re-run the point if the actual Ca is too far from the target Ca */
                isSteady = false;
                RESCALE_FORCE = true;
                t1 = std::chrono::system_clock::now();
                CURRENT_TIMESTEP = 0;
                timestep = INITIAL_TIMESTEP;
                TRIGGER_FORCE_RESCALE = true;
                if (rank == 0)
                    printf("    Capillary number missed target value = %f "
                           "(measured value was Ca = %f) \n ",
                           capillary_number, Ca);
            }

            if (RESCALE_FORCE == true && SET_CAPILLARY_NUMBER == true &&
                CURRENT_TIMESTEP > RESCALE_FORCE_AFTER_TIMESTEP) {
                TRIGGER_FORCE_RESCALE = true;
            }

            if (TRIGGER_FORCE_RESCALE) {
                RESCALE_FORCE = false;
                TRIGGER_FORCE_RESCALE = false;
                double RESCALE_FORCE_FACTOR = capillary_number / Ca;
                if (RESCALE_FORCE_FACTOR > 2.0)
                    RESCALE_FORCE_FACTOR = 2.0;
                if (RESCALE_FORCE_FACTOR < 0.5)
                    RESCALE_FORCE_FACTOR = 0.5;
                Fx *= RESCALE_FORCE_FACTOR;
                Fy *= RESCALE_FORCE_FACTOR;
                Fz *= RESCALE_FORCE_FACTOR;
                force_mag = sqrt(Fx * Fx + Fy * Fy + Fz * Fz);
                if (force_mag > 1e-3) {
                    Fx *= 1e-3 / force_mag; // impose ceiling for stability
                    Fy *= 1e-3 / force_mag;
                    Fz *= 1e-3 / force_mag;
                }
                if (rank == 0)
                    printf("    -- adjust force by factor %f \n ",
                           capillary_number / Ca);
                Averages->SetParams(rhoA, rhoB, tauA, tauB, Fx, Fy, Fz, alpha,
                                    beta);
                color_db->putVector<double>("F", {Fx, Fy, Fz});
            }
            if (isSteady) {
                Averages->Full();
                Averages->Write(timestep);
                analysis.WriteVisData(timestep, current_db, *Averages, Phi,
                                      Pressure, Velocity, fq, Den);
                analysis.finish();

                if (rank == 0) {
                    printf("** WRITE STEADY POINT *** ");
                    printf("Ca = %f, (previous = %f) \n", Ca, Ca_previous);
                    double h = Dm->voxel_length;
                    // pressures
                    double pA = Averages->gnb.p;
                    double pB = Averages->gwb.p;
                    double pAc = Averages->gnc.p;
                    double pBc = Averages->gwc.p;
                    double pAB = (pA - pB) / (h * 6.0 * alpha);
                    double pAB_connected = (pAc - pBc) / (h * 6.0 * alpha);
                    // connected contribution
                    double Vol_nc = Averages->gnc.V / Dm->Volume;
                    double Vol_wc = Averages->gwc.V / Dm->Volume;
                    double Vol_nd = Averages->gnd.V / Dm->Volume;
                    double Vol_wd = Averages->gwd.V / Dm->Volume;
                    double Mass_n = Averages->gnc.M + Averages->gnd.M;
                    double Mass_w = Averages->gwc.M + Averages->gwd.M;
                    double vAc_x = Averages->gnc.Px / Mass_n;
                    double vAc_y = Averages->gnc.Py / Mass_n;
                    double vAc_z = Averages->gnc.Pz / Mass_n;
                    double vBc_x = Averages->gwc.Px / Mass_w;
                    double vBc_y = Averages->gwc.Py / Mass_w;
                    double vBc_z = Averages->gwc.Pz / Mass_w;
                    // disconnected contribution
                    double vAd_x = Averages->gnd.Px / Mass_n;
                    double vAd_y = Averages->gnd.Py / Mass_n;
                    double vAd_z = Averages->gnd.Pz / Mass_n;
                    double vBd_x = Averages->gwd.Px / Mass_w;
                    double vBd_y = Averages->gwd.Py / Mass_w;
                    double vBd_z = Averages->gwd.Pz / Mass_w;

                    double flow_rate_A_connected =
                        Vol_nc *
                        (vAc_x * dir_x + vAc_y * dir_y + vAc_z * dir_z);
                    double flow_rate_B_connected =
                        Vol_wc *
                        (vBc_x * dir_x + vBc_y * dir_y + vBc_z * dir_z);
                    double flow_rate_A_disconnected =
                        (Vol_nd) *
                        (vAd_x * dir_x + vAd_y * dir_y + vAd_z * dir_z);
                    double flow_rate_B_disconnected =
                        (Vol_wd) *
                        (vBd_x * dir_x + vBd_y * dir_y + vBd_z * dir_z);

                    double kAeff_connected =
                        h * h * muA * flow_rate_A_connected / (force_mag);
                    double kBeff_connected =
                        h * h * muB * flow_rate_B_connected / (force_mag);

                    // Saturation normalized effective permeability to account for decoupled phases and
                    // effective porosity.
                    double kAeff_connected_low =
                        (1.0 - current_saturation) * h * h * muA *
                        flow_rate_A_connected / (force_mag);
                    double kBeff_connected_low = current_saturation * h * h *
                                                 muB * flow_rate_B_connected /
                                                 (force_mag);

                    double kAeff_disconnected =
                        h * h * muA * flow_rate_A_disconnected / (force_mag);
                    double kBeff_disconnected =
                        h * h * muB * flow_rate_B_disconnected / (force_mag);

                    double kAeff = h * h * muA * (flow_rate_A) / (force_mag);
                    double kBeff = h * h * muB * (flow_rate_B) / (force_mag);

                    // Saturation normalized effective permeability to account for decoupled phases and
                    // effective porosity.
                    double kAeff_low = (1.0 - current_saturation) * h * h *
                                       muA * (flow_rate_A) / (force_mag);
                    double kBeff_low = current_saturation * h * h * muB *
                                       (flow_rate_B) / (force_mag);

                    double viscous_pressure_drop =
                        (rhoA * volA + rhoB * volB) * force_mag;
                    double Mobility = muA / muB; // visc contrast
                    double eff_pres =
                        1.0 / (kAeff + kBeff); // effective pressure drop

                    bool WriteHeader = false;
                    FILE *kr_log_file = fopen("relperm.csv", "r");
                    if (kr_log_file != NULL)
                        fclose(kr_log_file);
                    else
                        WriteHeader = true;
                    kr_log_file = fopen("relperm.csv", "a");
                    if (WriteHeader) {
                        fprintf(kr_log_file, "timesteps sat.water ");
                        fprintf(kr_log_file, "eff.perm.oil.upper.bound "
                                             "eff.perm.water.upper.bound ");
                        fprintf(kr_log_file, "eff.perm.oil.lower.bound "
                                             "eff.perm.water.lower.bound ");
                        fprintf(kr_log_file,
                                "eff.perm.oil.connected.upper.bound "
                                "eff.perm.water.connected.upper.bound ");
                        fprintf(kr_log_file,
                                "eff.perm.oil.connected.lower.bound "
                                "eff.perm.water.connected.lower.bound ");
                        fprintf(kr_log_file, "eff.perm.oil.disconnected "
                                             "eff.perm.water.disconnected ");
                        fprintf(kr_log_file,
                                "cap.pressure cap.pressure.connected "
                                "pressure.drop Ca M eff.pressure\n");
                    }
                    fprintf(kr_log_file, "%i %.5g ", CURRENT_TIMESTEP,
                            current_saturation);
                    fprintf(kr_log_file, "%.5g %.5g ", kAeff, kBeff);
                    fprintf(kr_log_file, "%.5g %.5g ", kAeff_low, kBeff_low);
                    fprintf(kr_log_file, "%.5g %.5g ", kAeff_connected,
                            kBeff_connected);
                    fprintf(kr_log_file, "%.5g %.5g ", kAeff_connected_low,
                            kBeff_connected_low);
                    fprintf(kr_log_file, "%.5g %.5g ", kAeff_disconnected,
                            kBeff_disconnected);
                    fprintf(kr_log_file, "%.5g %.5g %.5g %.5g %.5g ", pAB,
                            pAB_connected, viscous_pressure_drop, Ca, Mobility);
                    fprintf(kr_log_file, "%.5g\n", eff_pres);
                    fclose(kr_log_file);

                    if (WettingConvention == "SCAL"){
                    	WriteHeader = false;
                    	FILE *scal_log_file = fopen("SCAL.csv", "r");
                    	if (scal_log_file != NULL)
                    		fclose(scal_log_file);
                    	else
                    		WriteHeader = true;
                    	scal_log_file = fopen("relperm.csv", "a");
                    	if (WriteHeader) {
                    		fprintf(scal_log_file, "timesteps sat.water ");
                    		fprintf(scal_log_file, "eff.perm.oil.upper.bound "
                    				"eff.perm.water.upper.bound ");
                    		fprintf(scal_log_file, "eff.perm.oil.lower.bound "
                    				"eff.perm.water.lower.bound ");
                    		fprintf(scal_log_file,
                    				"eff.perm.oil.connected.upper.bound "
                    				"eff.perm.water.connected.upper.bound ");
                    		fprintf(scal_log_file,
                    				"eff.perm.oil.connected.lower.bound "
                    				"eff.perm.water.connected.lower.bound ");
                    		fprintf(scal_log_file, "eff.perm.oil.disconnected "
                    				"eff.perm.water.disconnected ");
                    		fprintf(scal_log_file,
                    				"cap.pressure cap.pressure.connected "
                    				"pressure.drop Ca M eff.pressure\n");
                    	}
                    	fprintf(scal_log_file, "%i %.5g ", CURRENT_TIMESTEP,
                    			current_saturation);
                    	fprintf(scal_log_file, "%.5g %.5g ", kAeff, kBeff);
                    	fprintf(scal_log_file, "%.5g %.5g ", kAeff_low, kBeff_low);
                    	fprintf(scal_log_file, "%.5g %.5g ", kAeff_connected,
                    			kBeff_connected);
                    	fprintf(scal_log_file, "%.5g %.5g ", kAeff_connected_low,
                    			kBeff_connected_low);
                    	fprintf(scal_log_file, "%.5g %.5g ", kAeff_disconnected,
                    			kBeff_disconnected);
                    	fprintf(scal_log_file, "%.5g %.5g %.5g %.5g %.5g ", pAB,
                    			pAB_connected, viscous_pressure_drop, Ca, Mobility);
                    	fprintf(scal_log_file, "%.5g\n", eff_pres);
                    	fclose(scal_log_file);

                    }

                    printf("  Measured capillary number %f \n ", Ca);
                }
                if (SET_CAPILLARY_NUMBER) {
                    Fx *= capillary_number / Ca;
                    Fy *= capillary_number / Ca;
                    Fz *= capillary_number / Ca;
                    if (force_mag > 1e-3) {
                        Fx *= 1e-3 / force_mag; // impose ceiling for stability
                        Fy *= 1e-3 / force_mag;
                        Fz *= 1e-3 / force_mag;
                    }
                    if (rank == 0)
                        printf("    -- adjust force by factor %f \n ",
                               capillary_number / Ca);
                    Averages->SetParams(rhoA, rhoB, tauA, tauB, Fx, Fy, Fz,
                                        alpha, beta);
                    color_db->putVector<double>("F", {Fx, Fy, Fz});
                } else {
                    if (rank == 0) {
                        printf("** Continue to simulate steady *** \n ");
                        printf("Ca = %f, (previous = %f) \n", Ca, Ca_previous);
                    }
                }
            }
        }
    }
    analysis.finish();
    PROFILE_STOP("Update");

    PROFILE_STOP("Loop");
    PROFILE_SAVE("lbpm_color_simulator", 1);
    //************************************************************************
    // Compute the walltime per timestep
    auto t2 = std::chrono::system_clock::now();
    double cputime =
        std::chrono::duration<double>(t2 - t1).count() / CURRENT_TIMESTEP;
    // Performance obtained from each node
    double MLUPS = double(Np) / cputime / 1000000;

    if (rank == 0)
        printf("********************************************************\n");
    if (rank == 0)
        printf("CPU time = %f \n", cputime);
    if (rank == 0)
        printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
    return (MLUPS);
    MLUPS *= nprocs;
}

void ScaLBL_ColorModel::Run() {
    int nprocs = nprocx * nprocy * nprocz;
    const RankInfoStruct rank_info(rank, nprocx, nprocy, nprocz);
    int analysis_interval =
        1000; // number of timesteps in between in situ analysis
    if (analysis_db->keyExists("analysis_interval")) {
        analysis_interval = analysis_db->getScalar<int>("analysis_interval");
    }

    //************ MAIN ITERATION LOOP ***************************************/
    comm.barrier();
    PROFILE_START("Loop");
    //std::shared_ptr<Database> analysis_db;
    bool Regular = false;
    auto current_db = db->cloneDatabase();
    runAnalysis analysis(current_db, rank_info, ScaLBL_Comm, Dm, Np, Regular,
                         Map);
    //analysis.createThreads( analysis_method, 4 );
    auto t1 = std::chrono::system_clock::now();
    while (timestep < timestepMax) {
        PROFILE_START("Update");

        // *************ODD TIMESTEP*************
        timestep++;
        // Compute the Phase indicator field
        // Read for Aq, Bq happens in this routine (requires communication)
        ScaLBL_Comm->BiSendD3Q7AA(Aq, Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi,
                                     ScaLBL_Comm->FirstInterior(),
                                     ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->BiRecvD3Q7AA(Aq, Bq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0,
                                     ScaLBL_Comm->LastExterior(), Np);

        // Perform the collision operation
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        if (BoundaryCondition > 0 && BoundaryCondition < 5) {
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
        }
        // Halo exchange for phase field
        ScaLBL_Comm_Regular->SendHalo(Phi);

        ScaLBL_D3Q19_AAodd_Color(
            NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB,
            tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx, Nx * Ny,
            ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm_Regular->RecvHalo(Phi);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        // Set BCs
        if (BoundaryCondition == 3) {
            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }
        if (BoundaryCondition == 4) {
            din =
                ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 5) {
            ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
            ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
        }
        ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi,
                                 Velocity, rhoA, rhoB, tauA, tauB, alpha, beta,
                                 Fx, Fy, Fz, Nx, Nx * Ny, 0,
                                 ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->Barrier();

        // *************EVEN TIMESTEP*************
        timestep++;
        // Compute the Phase indicator field
        ScaLBL_Comm->BiSendD3Q7AA(Aq, Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi,
                                      ScaLBL_Comm->FirstInterior(),
                                      ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->BiRecvD3Q7AA(Aq, Bq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0,
                                      ScaLBL_Comm->LastExterior(), Np);

        // Perform the collision operation
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        // Halo exchange for phase field
        if (BoundaryCondition > 0 && BoundaryCondition < 5) {
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
        }
        ScaLBL_Comm_Regular->SendHalo(Phi);
        ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA,
                                  rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx,
                                  Nx * Ny, ScaLBL_Comm->FirstInterior(),
                                  ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm_Regular->RecvHalo(Phi);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        // Set boundary conditions
        if (BoundaryCondition == 3) {
            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 4) {
            din =
                ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 5) {
            ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
            ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
        }
        ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA,
                                  rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx,
                                  Nx * Ny, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->Barrier();
        //************************************************************************
        PROFILE_STOP("Update");

        if (rank == 0 && timestep % analysis_interval == 0 &&
            BoundaryCondition == 4) {
            printf("%i %f \n", timestep, din);
        }
        // Run the analysis
        analysis.basic(timestep, current_db, *Averages, Phi, Pressure, Velocity,
                       fq, Den);
    }
    analysis.finish();
    PROFILE_STOP("Loop");
    PROFILE_SAVE("lbpm_color_simulator", 1);
    //************************************************************************
    ScaLBL_Comm->Barrier();
    if (rank == 0)
        printf("---------------------------------------------------------------"
               "----\n");
    // Compute the walltime per timestep
    auto t2 = std::chrono::system_clock::now();
    double cputime = std::chrono::duration<double>(t2 - t1).count() / timestep;
    // Performance obtained from each node
    double MLUPS = double(Np) / cputime / 1000000;

    if (rank == 0)
        printf("********************************************************\n");
    if (rank == 0)
        printf("CPU time = %f \n", cputime);
    if (rank == 0)
        printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
    MLUPS *= nprocs;
    if (rank == 0)
        printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
    if (rank == 0)
        printf("********************************************************\n");

    // ************************************************************************
}

void ScaLBL_ColorModel::WriteDebug() {
    // Copy back final phase indicator field and convert to regular layout
    DoubleArray PhaseField(Nx, Ny, Nz);
    //ScaLBL_Comm->RegularLayout(Map,Phi,PhaseField);
    ScaLBL_CopyToHost(PhaseField.data(), Phi, sizeof(double) * N);

    FILE *OUTFILE;
    sprintf(LocalRankFilename, "Phase.%05i.raw", rank);
    OUTFILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, OUTFILE);
    fclose(OUTFILE);

    ScaLBL_Comm->RegularLayout(Map, &Den[0], PhaseField);
    FILE *AFILE;
    sprintf(LocalRankFilename, "A.%05i.raw", rank);
    AFILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, AFILE);
    fclose(AFILE);

    ScaLBL_Comm->RegularLayout(Map, &Den[Np], PhaseField);
    FILE *BFILE;
    sprintf(LocalRankFilename, "B.%05i.raw", rank);
    BFILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, BFILE);
    fclose(BFILE);

    ScaLBL_Comm->RegularLayout(Map, Pressure, PhaseField);
    FILE *PFILE;
    sprintf(LocalRankFilename, "Pressure.%05i.raw", rank);
    PFILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, PFILE);
    fclose(PFILE);

    ScaLBL_Comm->RegularLayout(Map, &Velocity[0], PhaseField);
    FILE *VELX_FILE;
    sprintf(LocalRankFilename, "Velocity_X.%05i.raw", rank);
    VELX_FILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, VELX_FILE);
    fclose(VELX_FILE);

    ScaLBL_Comm->RegularLayout(Map, &Velocity[Np], PhaseField);
    FILE *VELY_FILE;
    sprintf(LocalRankFilename, "Velocity_Y.%05i.raw", rank);
    VELY_FILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, VELY_FILE);
    fclose(VELY_FILE);

    ScaLBL_Comm->RegularLayout(Map, &Velocity[2 * Np], PhaseField);
    FILE *VELZ_FILE;
    sprintf(LocalRankFilename, "Velocity_Z.%05i.raw", rank);
    VELZ_FILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, VELZ_FILE);
    fclose(VELZ_FILE);
}
