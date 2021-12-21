/*
Two-fluid greyscale color lattice boltzmann model
 */
#include "models/GreyscaleColorModel.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"
#include "common/Communication.h"
#include "common/ReadMicroCT.h"
#include <stdlib.h>
#include <time.h>

template <class TYPE> void DeleteArray(const TYPE *p) { delete[] p; }

ScaLBL_GreyscaleColorModel::ScaLBL_GreyscaleColorModel(
    int RANK, int NP, const Utilities::MPI &COMM)
    : rank(RANK), nprocs(NP), Restart(0), timestep(0), timestepMax(0), tauA(0),
      tauB(0), tauA_eff(0), tauB_eff(0), rhoA(0), rhoB(0), alpha(0), beta(0),
      Fx(0), Fy(0), Fz(0), flux(0), din(0), dout(0), inletA(0), inletB(0),
      outletA(0), outletB(0), GreyPorosity(0), RecoloringOff(0), Nx(0), Ny(0),
      Nz(0), N(0), Np(0), nprocx(0), nprocy(0), nprocz(0), BoundaryCondition(0),
      Lx(0), Ly(0), Lz(0), comm(COMM) {
    REVERSE_FLOW_DIRECTION = false;
}
ScaLBL_GreyscaleColorModel::~ScaLBL_GreyscaleColorModel() {}
void ScaLBL_GreyscaleColorModel::ReadParams(string filename) {
    // read the input database
    db = std::make_shared<Database>(filename);
    domain_db = db->getDatabase("Domain");
    greyscaleColor_db = db->getDatabase("Color");
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
    RecoloringOff = false;
    //W=1.0;

    // Color Model parameters
    if (greyscaleColor_db->keyExists("timestepMax")) {
        timestepMax = greyscaleColor_db->getScalar<int>("timestepMax");
    }
    if (greyscaleColor_db->keyExists("tauA")) {
        tauA = greyscaleColor_db->getScalar<double>("tauA");
    }
    if (greyscaleColor_db->keyExists("tauB")) {
        tauB = greyscaleColor_db->getScalar<double>("tauB");
    }
    tauA_eff = greyscaleColor_db->getWithDefault<double>("tauA_eff", tauA);
    tauB_eff = greyscaleColor_db->getWithDefault<double>("tauB_eff", tauB);
    if (greyscaleColor_db->keyExists("rhoA")) {
        rhoA = greyscaleColor_db->getScalar<double>("rhoA");
    }
    if (greyscaleColor_db->keyExists("rhoB")) {
        rhoB = greyscaleColor_db->getScalar<double>("rhoB");
    }
    if (greyscaleColor_db->keyExists("F")) {
        Fx = greyscaleColor_db->getVector<double>("F")[0];
        Fy = greyscaleColor_db->getVector<double>("F")[1];
        Fz = greyscaleColor_db->getVector<double>("F")[2];
    }
    if (greyscaleColor_db->keyExists("alpha")) {
        alpha = greyscaleColor_db->getScalar<double>("alpha");
    }
    if (greyscaleColor_db->keyExists("beta")) {
        beta = greyscaleColor_db->getScalar<double>("beta");
    }
    if (greyscaleColor_db->keyExists("Restart")) {
        Restart = greyscaleColor_db->getScalar<bool>("Restart");
    }
    if (greyscaleColor_db->keyExists("din")) {
        din = greyscaleColor_db->getScalar<double>("din");
    }
    if (greyscaleColor_db->keyExists("dout")) {
        dout = greyscaleColor_db->getScalar<double>("dout");
    }
    if (greyscaleColor_db->keyExists("flux")) {
        flux = greyscaleColor_db->getScalar<double>("flux");
    }
    if (greyscaleColor_db->keyExists("RecoloringOff")) {
        RecoloringOff = greyscaleColor_db->getScalar<bool>("RecoloringOff");
    }
    inletA = 1.f;
    inletB = 0.f;
    outletA = 0.f;
    outletB = 1.f;
    //if (BoundaryCondition==4) flux *= rhoA; // mass flux must adjust for density (see formulation for details)

    BoundaryCondition = 0;
    if (domain_db->keyExists("BC")) {
        BoundaryCondition = domain_db->getScalar<int>("BC");
    }

    // Override user-specified boundary condition for specific protocols
    auto protocol =
        greyscaleColor_db->getWithDefault<std::string>("protocol", "none");
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
    }
}

void ScaLBL_GreyscaleColorModel::SetDomain() {
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
    Averages = std::shared_ptr<GreyPhaseAnalysis>(
        new GreyPhaseAnalysis(Dm)); // TwoPhase analysis object
    comm.barrier();
    Dm->CommInit();
    comm.barrier();
    // Read domain parameters
    rank = Dm->rank();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();
}

void ScaLBL_GreyscaleColorModel::ReadInput() {

    sprintf(LocalRankString, "%05d", rank);
    sprintf(LocalRankFilename, "%s%s", "ID.", LocalRankString);
    sprintf(LocalRestartFile, "%s%s", "Restart.", LocalRankString);

    if (greyscaleColor_db->keyExists("image_sequence")) {
        auto ImageList =
            greyscaleColor_db->getVector<std::string>("image_sequence");
        int IMAGE_INDEX =
            greyscaleColor_db->getWithDefault<int>("image_index", 0);
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
    if (rank == 0)
        printf("Initialized solid phase -- Converting to Signed Distance "
               "function \n");
    CalcDist(Averages->SDs, id_solid, *Mask);

    if (rank == 0)
        cout << "Domain set." << endl;
}

void ScaLBL_GreyscaleColorModel::AssignComponentLabels() {
    // Initialize impermeability solid nodes and grey nodes
    // Key input parameters:
    // 1. ComponentLabels
    //    labels for various impermeable minerals and grey nodes
    // 2. ComponentAffinity
    //    for impermeable minerals, this is same as the wettability phase field in the normal color model
    //    for grey nodes, this is effectively the initial phase field values
    // **Convention for ComponentLabels:
    //   (1) zero and negative integers are for impermeability minerals
    //   (2) positive integers > 2 are for grey nodes
    //   (3) label = 1 and 2 are always conserved for open node of non-wetting and wetting phase, respectively.
    double *phase;
    phase = new double[N];

    size_t NLABELS = 0;
    signed char VALUE = 0;
    double AFFINITY = 0.f;

    auto LabelList = greyscaleColor_db->getVector<int>("ComponentLabels");
    auto AffinityList =
        greyscaleColor_db->getVector<double>("ComponentAffinity");

    NLABELS = LabelList.size();
    if (NLABELS != AffinityList.size()) {
        ERROR("Error: ComponentLabels and ComponentAffinity must be the same "
              "length! \n");
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
                for (size_t idx = 0; idx < NLABELS; idx++) {
                    //printf("idx=%i, value=%i, %i, \n",idx, VALUE,LabelList[idx]);
                    if (VALUE == LabelList[idx]) {
                        AFFINITY = AffinityList[idx];

                        label_count[idx] += 1.0;
                        idx = NLABELS;
                        //Mask->id[n] = 0; // set mask to zero since this is an immobile component
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
        printf("Number of component labels: %lu \n", NLABELS);
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

    ScaLBL_CopyToDevice(Phi, phase, N * sizeof(double));
    ScaLBL_Comm->Barrier();
    delete[] phase;
}

void ScaLBL_GreyscaleColorModel::
    AssignGreySolidLabels() //apply capillary penalty wetting strength W
{
    // ONLY initialize grey nodes
    // Key input parameters:
    // 1. GreySolidLabels
    //    labels for grey nodes
    // 2. GreySolidAffinity
    //    ranges [-1,1]
    //    water-wet > 0
    //    oil-wet   < 0
    //    neutral   = 0 (i.e. no penalty)
    double *GreySolidW_host = new double[Np];
    double *GreySn_host = new double[Np];
    double *GreySw_host = new double[Np];
    double *GreyKn_host = new double[Np];
    double *GreyKw_host = new double[Np];

    size_t NLABELS = 0;
    signed char VALUE = 0;
    double AFFINITY = 0.f;
    double Sn, Sw; //end-point saturation of greynodes set by users
    double Kn, Kw; // endpoint effective permeability

    auto LabelList = greyscaleColor_db->getVector<int>("GreySolidLabels");
    auto AffinityList =
        greyscaleColor_db->getVector<double>("GreySolidAffinity");
    auto SnList = greyscaleColor_db->getVector<double>("grey_endpoint_A");
    auto SwList = greyscaleColor_db->getVector<double>("grey_endpoint_B");
    auto KnList =
        greyscaleColor_db->getVector<double>("grey_endpoint_permeability_A");
    auto KwList =
        greyscaleColor_db->getVector<double>("grey_endpoint_permeability_B");

    NLABELS = LabelList.size();
    if (NLABELS != AffinityList.size()) {
        ERROR("Error: GreySolidLabels and GreySolidAffinity must be the same "
              "length! \n");
    }
    if (NLABELS != SnList.size() || NLABELS != SwList.size()) {
        ERROR("Error: GreySolidLabels, grey_endpoint_A, and grey_endpoint_B "
              "must be the same length! \n");
    }
    if (NLABELS != KnList.size() || NLABELS != KwList.size()) {
        ERROR("Error: GreySolidLabels, grey_endpoint_permeability_A, and "
              "grey_endpoint_permeability_B must be the same length! \n");
    }

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                VALUE = id[n];
                AFFINITY =
                    0.f; //all nodes except the specified grey nodes have grey-solid affinity = 0.0
                Sn = 99.0;
                Sw = -99.0;
                Kn = 0.0;
                Kw = 0.0;
                // Assign the affinity from the paired list
                for (unsigned int idx = 0; idx < NLABELS; idx++) {
                    if (VALUE == LabelList[idx]) {
                        AFFINITY = AffinityList[idx];
                        Sn = SnList[idx];
                        Sw = SwList[idx];
                        Kn = KnList[idx];
                        Kw = KwList[idx];
                        idx = NLABELS;
                    }
                }
                int idx = Map(i, j, k);
                if (!(idx < 0)) {
                    GreySolidW_host[idx] = AFFINITY;
                    GreySn_host[idx] = Sn;
                    GreySw_host[idx] = Sw;
                    GreyKn_host[idx] = Kn;
                    GreyKw_host[idx] = Kw;
                }
            }
        }
    }

    if (rank == 0) {
        printf("Number of Grey-solid labels: %lu \n", NLABELS);
        for (unsigned int idx = 0; idx < NLABELS; idx++) {
            VALUE = LabelList[idx];
            AFFINITY = AffinityList[idx];
            Sn = SnList[idx];
            Sw = SwList[idx];
            //printf("   grey-solid label=%d, grey-solid affinity=%f\n",VALUE,AFFINITY);
            printf("   grey-solid label=%d, grey-solid affinity=%.3g, "
                   "grey-solid Sn=%.3g, grey-solid Sw=%.3g\n",
                   VALUE, AFFINITY, Sn, Sw);
        }
        printf("NOTE: grey-solid affinity>0: water-wet || grey-solid "
               "affinity<0: oil-wet \n");
    }

    ScaLBL_CopyToDevice(GreySolidW, GreySolidW_host, Np * sizeof(double));
    ScaLBL_CopyToDevice(GreySn, GreySn_host, Np * sizeof(double));
    ScaLBL_CopyToDevice(GreySw, GreySw_host, Np * sizeof(double));
    ScaLBL_CopyToDevice(GreyKn, GreySn_host, Np * sizeof(double));
    ScaLBL_CopyToDevice(GreyKw, GreySw_host, Np * sizeof(double));
    ScaLBL_Comm->Barrier();
    delete[] GreySolidW_host;
    delete[] GreySn_host;
    delete[] GreySw_host;
}
////----------------------------------------------------------------------------------------------------------//

void ScaLBL_GreyscaleColorModel::AssignGreyPoroPermLabels() {

    double *Porosity, *Permeability;
    Porosity = new double[Np];
    Permeability = new double[Np];

    size_t NLABELS = 0;
    signed char VALUE = 0;
    double POROSITY =
        1.f; //default: label 1 or 2, i.e. open nodes and porosity=1.0
    double PERMEABILITY = 1.f;

    auto LabelList = greyscaleColor_db->getVector<int>("GreySolidLabels");
    auto PorosityList = greyscaleColor_db->getVector<double>("PorosityList");
    auto PermeabilityList =
        greyscaleColor_db->getVector<double>("PermeabilityList");

    NLABELS = LabelList.size();
    if (LabelList.size() != PorosityList.size()) {
        ERROR("Error: GreySolidLabels and PorosityList must be the same "
              "length! \n");
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
                POROSITY =
                    1.f; //default: label 1 or 2, i.e. open nodes and porosity=1.0
                // Assign the affinity from the paired list
                for (size_t idx = 0; idx < NLABELS; idx++) {
                    //printf("idx=%i, value=%i, %i, \n",idx, VALUE,LabelList[idx]);
                    if (VALUE == LabelList[idx]) {
                        POROSITY = PorosityList[idx];
                        label_count[idx] += 1.0;
                        idx = NLABELS;
                        //Mask->id[n] = 0; // set mask to zero since this is an immobile component
                    }
                }
                int idx = Map(i, j, k);
                if (!(idx < 0)) {
                    if (POROSITY <= 0.0) {
                        ERROR("Error: Porosity for grey voxels must be 0.0 < "
                              "Porosity <= 1.0 !\n");
                    } else {
                        Porosity[idx] = POROSITY;
                    }
                }
            }
        }
    }

    if (NLABELS != PermeabilityList.size()) {
        ERROR("Error: GreySolidLabels and PermeabilityList must be the same "
              "length! \n");
    }
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                VALUE = id[n];
                PERMEABILITY = 1.f;
                // Assign the affinity from the paired list
                for (unsigned int idx = 0; idx < NLABELS; idx++) {
                    //printf("idx=%i, value=%i, %i, \n",idx, VALUE,LabelList[idx]);
                    if (VALUE == LabelList[idx]) {
                        PERMEABILITY = PermeabilityList[idx];
                        idx = NLABELS;
                        //Mask->id[n] = 0; // set mask to zero since this is an immobile component
                    }
                }
                int idx = Map(i, j, k);
                if (!(idx < 0)) {
                    if (PERMEABILITY <= 0.0) {
                        ERROR("Error: Permeability for grey voxel must be > "
                              "0.0 ! \n");
                    } else {
                        Permeability[idx] =
                            PERMEABILITY / Dm->voxel_length / Dm->voxel_length;
                    }
                }
            }
        }
    }

    // Set Dm to match Mask
    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = Mask->id[i];

    for (size_t idx = 0; idx < NLABELS; idx++)
        label_count_global[idx] = Dm->Comm.sumReduce(label_count[idx]);

    //Initialize a weighted porosity after considering grey voxels
    GreyPorosity = 0.0;
    for (unsigned int idx = 0; idx < NLABELS; idx++) {
        double volume_fraction =
            double(label_count_global[idx]) /
            double((Nx - 2) * (Ny - 2) * (Nz - 2) * nprocs);
        GreyPorosity += volume_fraction * PorosityList[idx];
    }

    if (rank == 0) {
        printf("Image resolution: %.5g [um/voxel]\n", Dm->voxel_length);
        printf("Number of Grey-fluid labels: %lu \n", NLABELS);
        for (unsigned int idx = 0; idx < NLABELS; idx++) {
            VALUE = LabelList[idx];
            POROSITY = PorosityList[idx];
            PERMEABILITY = PermeabilityList[idx];
            double volume_fraction =
                double(label_count_global[idx]) /
                double((Nx - 2) * (Ny - 2) * (Nz - 2) * nprocs);
            printf("   grey-fluid label=%d, porosity=%.3g, permeability=%.3g "
                   "[um^2] (=%.3g [voxel^2]), volume fraction=%.3g\n",
                   VALUE, POROSITY, PERMEABILITY,
                   PERMEABILITY / Dm->voxel_length / Dm->voxel_length,
                   volume_fraction);
            printf("                        effective porosity=%.3g\n",
                   volume_fraction * POROSITY);
        }
        printf("The weighted porosity, considering both open and grey voxels, "
               "is %.3g\n",
               GreyPorosity);
    }

    ScaLBL_CopyToDevice(Porosity_dvc, Porosity, Np * sizeof(double));
    ScaLBL_CopyToDevice(Permeability_dvc, Permeability, Np * sizeof(double));
    ScaLBL_Comm->Barrier();
    delete[] Porosity;
    delete[] Permeability;
}

void ScaLBL_GreyscaleColorModel::Create() {
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
    //ScaLBL_AllocateDeviceMemory((void **) &Psi, sizeof(double)*Nx*Ny*Nz);//greyscale potential
    ScaLBL_AllocateDeviceMemory((void **)&Pressure, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&Velocity, 3 * sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&MobilityRatio, sizeof(double) * Np);
    //ScaLBL_AllocateDeviceMemory((void **) &GreySolidPhi, sizeof(double)*Nx*Ny*Nz);
    //ScaLBL_AllocateDeviceMemory((void **) &GreySolidGrad, 3*sizeof(double)*Np);
    ScaLBL_AllocateDeviceMemory((void **)&GreySolidW, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&GreySn, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&GreySw, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&GreyKn, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&GreyKw, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&Porosity_dvc, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&Permeability_dvc,
                                sizeof(double) * Np);
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

    // initialize phi based on PhaseLabel (include solid component labels)
    AssignComponentLabels(); //do open/black/grey nodes initialization
    AssignGreySolidLabels();
    AssignGreyPoroPermLabels();
    Averages->SetParams(rhoA, rhoB, tauA, tauB, Fx, Fy, Fz, alpha, beta,
                        GreyPorosity);
    ScaLBL_Comm->RegularLayout(
        Map, Porosity_dvc,
        Averages->Porosity); //porosity doesn't change over time
}

void ScaLBL_GreyscaleColorModel::Initialize() {
    /*
	 * This function initializes model
	 */
    if (rank == 0)
        printf("Initializing distributions \n");
    ScaLBL_D3Q19_Init(fq, Np);
    //ScaLBL_D3Q19_GreyscaleColor_Init(fq, Porosity_dvc, Np);

    if (rank == 0)
        printf("Initializing phase field \n");
    ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0,
                           ScaLBL_Comm->LastExterior(), Np);
    ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq,
                           ScaLBL_Comm->FirstInterior(),
                           ScaLBL_Comm->LastInterior(), Np);

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

        if (rank == 0)
            printf("Initializing phase field from Restart\n");
        ScaLBL_PhaseField_InitFromRestart(Den, Aq, Bq, 0,
                                          ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_PhaseField_InitFromRestart(Den, Aq, Bq,
                                          ScaLBL_Comm->FirstInterior(),
                                          ScaLBL_Comm->LastInterior(), Np);
    }

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
    //ScaLBL_CopyToHost(Averages->Phi.data(),Phi,N*sizeof(double));
}

void ScaLBL_GreyscaleColorModel::Run() {
    int nprocs = nprocx * nprocy * nprocz;
    const RankInfoStruct rank_info(rank, nprocx, nprocy, nprocz);

    bool SET_CAPILLARY_NUMBER = false;
    bool RESCALE_FORCE = false;
    bool MORPH_ADAPT = false;
    bool USE_MORPH = false;
    bool USE_SEED = false;
    bool USE_DIRECT = false;
    int MAX_MORPH_TIMESTEPS =
        50000; // maximum number of LBM timesteps to spend in morphological adaptation routine
    int MIN_STEADY_TIMESTEPS = 100000;
    int MAX_STEADY_TIMESTEPS = 200000;
    int RESCALE_FORCE_AFTER_TIMESTEP = 0;
    int RAMP_TIMESTEPS =
        0; //50000;		 // number of timesteps to run initially (to get a reasonable velocity field before other pieces kick in)
    int CURRENT_MORPH_TIMESTEPS =
        0; // counter for number of timesteps spent in  morphological adaptation routine (reset each time)
    int CURRENT_STEADY_TIMESTEPS =
        0; // counter for number of timesteps spent in  morphological adaptation routine (reset each time)
    int morph_interval = 100000;
    int analysis_interval =
        1000; // number of timesteps in between in situ analysis
    int morph_timesteps = 0;
    double morph_delta = 0.0;
    double seed_water = 0.0;
    double capillary_number = 0.0;
    double tolerance = 0.01;
    double Ca_previous = 0.f;
    double initial_volume = 0.0;
    double delta_volume = 0.0;
    double delta_volume_target = 0.0;

    //TODO -------- For temporary use - should be included in the analysis framework later -------------
    int visualization_interval = 50000;
    int restart_interval = 100000;
    if (analysis_db->keyExists("visualization_interval")) {
        visualization_interval =
            analysis_db->getScalar<int>("visualization_interval");
    }
    if (analysis_db->keyExists("restart_interval")) {
        restart_interval = analysis_db->getScalar<int>("restart_interval");
    }
    //-------------------------------------------------------------------------------------------------

    /* history for morphological algoirthm */
    double KRA_MORPH_FACTOR = 0.5;
    double volA_prev = 0.0;
    double log_krA_prev = 1.0;
    double log_krA_target = 1.0;
    double log_krA = 1.0;
    double slope_krA_volume = 0.0;
    if (greyscaleColor_db->keyExists("vol_A_previous")) {
        volA_prev = greyscaleColor_db->getScalar<double>("vol_A_previous");
    }
    if (greyscaleColor_db->keyExists("log_krA_previous")) {
        log_krA_prev = greyscaleColor_db->getScalar<double>("log_krA_previous");
    }
    if (greyscaleColor_db->keyExists("krA_morph_factor")) {
        KRA_MORPH_FACTOR =
            greyscaleColor_db->getScalar<double>("krA_morph_factor");
    }

    /* defaults for simulation protocols */
    auto protocol =
        greyscaleColor_db->getWithDefault<std::string>("protocol", "none");
    if (protocol == "seed water") {
        morph_delta = -0.05;
        seed_water = 0.01;
        USE_SEED = true;
        USE_MORPH = true;
    }

    if (greyscaleColor_db->keyExists("capillary_number")) {
        capillary_number =
            greyscaleColor_db->getScalar<double>("capillary_number");
        SET_CAPILLARY_NUMBER = true;
    }
    if (greyscaleColor_db->keyExists("rescale_force_after_timestep")) {
        RESCALE_FORCE_AFTER_TIMESTEP =
            greyscaleColor_db->getScalar<int>("rescale_force_after_timestep");
        RESCALE_FORCE = true;
    }
    if (greyscaleColor_db->keyExists("timestep")) {
        timestep = greyscaleColor_db->getScalar<int>("timestep");
    }
    if (BoundaryCondition != 0 && BoundaryCondition != 5 &&
        SET_CAPILLARY_NUMBER == true) {
        if (rank == 0)
            printf("WARINING: capillary number target only supported for BC = "
                   "0 or 5 \n");
        SET_CAPILLARY_NUMBER = false;
    }
    if (analysis_db->keyExists("seed_water")) {
        seed_water = analysis_db->getScalar<double>("seed_water");
        if (rank == 0)
            printf("Seed water in oil %f (seed_water) \n", seed_water);
        USE_SEED = true;
    }
    if (analysis_db->keyExists("morph_delta")) {
        morph_delta = analysis_db->getScalar<double>("morph_delta");
        if (rank == 0)
            printf("Target volume change %f (morph_delta) \n", morph_delta);
    }
    if (analysis_db->keyExists("morph_interval")) {
        morph_interval = analysis_db->getScalar<int>("morph_interval");
        USE_MORPH = true;
    }
    if (analysis_db->keyExists("tolerance")) {
        tolerance = analysis_db->getScalar<double>("tolerance");
    }
    if (analysis_db->keyExists("analysis_interval")) {
        analysis_interval = analysis_db->getScalar<int>("analysis_interval");
    }
    if (analysis_db->keyExists("min_steady_timesteps")) {
        MIN_STEADY_TIMESTEPS =
            analysis_db->getScalar<int>("min_steady_timesteps");
    }
    if (analysis_db->keyExists("max_steady_timesteps")) {
        MAX_STEADY_TIMESTEPS =
            analysis_db->getScalar<int>("max_steady_timesteps");
    }
    if (analysis_db->keyExists("max_morph_timesteps")) {
        MAX_MORPH_TIMESTEPS =
            analysis_db->getScalar<int>("max_morph_timesteps");
    }

    if (rank == 0) {
        printf("********************************************************\n");
        if (protocol == "seed water") {
            printf("  using protocol =  seed water \n");
            printf("     min_steady_timesteps = %i \n", MIN_STEADY_TIMESTEPS);
            printf("     max_steady_timesteps = %i \n", MAX_STEADY_TIMESTEPS);
            printf("     tolerance = %f \n", tolerance);
            printf("     morph_delta = %f \n", morph_delta);
            printf("     seed_water = %f \n", seed_water);
        }
        printf("No. of timesteps: %i \n", timestepMax);
        fflush(stdout);
    }

    //.......create and start timer............
    ScaLBL_Comm->Barrier();
    comm.barrier();
    //.........................................

    //************ MAIN ITERATION LOOP ***************************************/
    PROFILE_START("Loop");
    //std::shared_ptr<Database> analysis_db;
    auto current_db = db->cloneDatabase();
    //runAnalysis analysis( current_db, rank_info, ScaLBL_Comm, Dm, Np, Regular, Map );
    //analysis.createThreads( analysis_method, 4 );
    auto t1 = std::chrono::system_clock::now();
    while (timestep < timestepMax) {
        //if ( rank==0 ) { printf("Running timestep %i (%i MB)\n",timestep+1,(int)(Utilities::getMemoryUsage()/1048576)); }
        PROFILE_START("Update");
        // *************ODD TIMESTEP*************
        timestep++;
        // Compute the Phase indicator field
        // Read for Aq, Bq happens in this routine (requires communication)
        ScaLBL_Comm->BiSendD3Q7AA(Aq, Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi,
                                     ScaLBL_Comm->FirstInterior(),
                                     ScaLBL_Comm->LastInterior(), Np);
        //ScaLBL_Update_GreyscalePotential(dvcMap,Phi,Psi,Porosity_dvc,Permeability_dvc,alpha,W,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->BiRecvD3Q7AA(Aq, Bq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0,
                                     ScaLBL_Comm->LastExterior(), Np);
        //ScaLBL_Update_GreyscalePotential(dvcMap,Phi,Psi,Porosity_dvc,Permeability_dvc,alpha,W,0,ScaLBL_Comm->LastExterior(), Np);

        // Perform the collision operation
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        if (BoundaryCondition > 0 && BoundaryCondition < 5) {
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
        }
        // Halo exchange for phase field
        ScaLBL_Comm_Regular->SendHalo(Phi);
        ScaLBL_D3Q19_AAodd_GreyscaleColor_CP(
            NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, GreySolidW, GreySn,
            GreySw, GreyKn, GreyKw, Porosity_dvc, Permeability_dvc, Velocity,
            MobilityRatio, Pressure, rhoA, rhoB, tauA, tauB, tauA_eff, tauB_eff,
            alpha, beta, Fx, Fy, Fz, RecoloringOff, Nx, Nx * Ny,
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

        ScaLBL_D3Q19_AAodd_GreyscaleColor_CP(
            NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, GreySolidW, GreySn,
            GreySw, GreyKn, GreyKw, Porosity_dvc, Permeability_dvc, Velocity,
            MobilityRatio, Pressure, rhoA, rhoB, tauA, tauB, tauA_eff, tauB_eff,
            alpha, beta, Fx, Fy, Fz, RecoloringOff, Nx, Nx * Ny, 0,
            ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->Barrier();

        // *************EVEN TIMESTEP*************
        timestep++;
        // Compute the Phase indicator field
        ScaLBL_Comm->BiSendD3Q7AA(Aq, Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi,
                                      ScaLBL_Comm->FirstInterior(),
                                      ScaLBL_Comm->LastInterior(), Np);
        //ScaLBL_Update_GreyscalePotential(dvcMap,Phi,Psi,Porosity_dvc,Permeability_dvc,alpha,W,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->BiRecvD3Q7AA(Aq, Bq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0,
                                      ScaLBL_Comm->LastExterior(), Np);
        //ScaLBL_Update_GreyscalePotential(dvcMap,Phi,Psi,Porosity_dvc,Permeability_dvc,alpha,W,0,ScaLBL_Comm->LastExterior(), Np);

        // Perform the collision operation
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        // Halo exchange for phase field
        if (BoundaryCondition > 0 && BoundaryCondition < 5) {
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
        }
        ScaLBL_Comm_Regular->SendHalo(Phi);
        ScaLBL_D3Q19_AAeven_GreyscaleColor_CP(
            dvcMap, fq, Aq, Bq, Den, Phi, GreySolidW, GreySn, GreySw, GreyKn,
            GreyKw, Porosity_dvc, Permeability_dvc, Velocity, MobilityRatio,
            Pressure, rhoA, rhoB, tauA, tauB, tauA_eff, tauB_eff, alpha, beta,
            Fx, Fy, Fz, RecoloringOff, Nx, Nx * Ny,
            ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
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

        ScaLBL_D3Q19_AAeven_GreyscaleColor_CP(
            dvcMap, fq, Aq, Bq, Den, Phi, GreySolidW, GreySn, GreySw, GreyKn,
            GreyKw, Porosity_dvc, Permeability_dvc, Velocity, MobilityRatio,
            Pressure, rhoA, rhoB, tauA, tauB, tauA_eff, tauB_eff, alpha, beta,
            Fx, Fy, Fz, RecoloringOff, Nx, Nx * Ny, 0,
            ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->Barrier();
        //************************************************************************
        PROFILE_STOP("Update");

        //TODO For temporary use - writing Restart and Vis files should be included in the analysis framework in the future
        if (timestep % restart_interval == 0) {
            //Use rank=0 write out Restart.db
            if (rank == 0) {
                greyscaleColor_db->putScalar<int>("timestep", timestep);
                greyscaleColor_db->putScalar<bool>("Restart", true);
                current_db->putDatabase("Color", greyscaleColor_db);
                std::ofstream OutStream("Restart.db");
                current_db->print(OutStream, "");
                OutStream.close();
            }
            //Write out Restart data.
            std::shared_ptr<double> cDen;
            std::shared_ptr<double> cfq;
            cDen = std::shared_ptr<double>(new double[2 * Np],
                                           DeleteArray<double>);
            cfq = std::shared_ptr<double>(new double[19 * Np],
                                          DeleteArray<double>);
            ScaLBL_CopyToHost(
                cDen.get(), Den,
                2 * Np * sizeof(double)); // Copy restart data to the CPU
            ScaLBL_CopyToHost(
                cfq.get(), fq,
                19 * Np * sizeof(double)); // Copy restart data to the CPU

            ofstream RESTARTFILE(LocalRestartFile, ios::binary);
            double value;
            for (int n = 0; n < Np; n++) {
                // Write the two density values
                value = cDen.get()[n];
                RESTARTFILE.write((char *)&value, sizeof(value));
                value = cDen.get()[Np + n];
                RESTARTFILE.write((char *)&value, sizeof(value));
            }
            for (int n = 0; n < Np; n++) {
                // Write the distributions
                for (int q = 0; q < 19; q++) {
                    value = cfq.get()[q * Np + n];
                    RESTARTFILE.write((char *)&value, sizeof(value));
                }
            }
            RESTARTFILE.close();
            comm.barrier();
        }
        if (timestep % visualization_interval == 0) {
            WriteVisFiles();
        }
        //-----------------------------------------------------------------------------------------------------------------

        if (rank == 0 && timestep % analysis_interval == 0 &&
            BoundaryCondition == 4) {
            printf("%i %.5g \n", timestep, din);
        }

        if (timestep % analysis_interval == 0) {
            ScaLBL_Comm->RegularLayout(Map, Pressure, Averages->Pressure);
            ScaLBL_Comm->RegularLayout(Map, MobilityRatio,
                                       Averages->MobilityRatio);
            ScaLBL_Comm->RegularLayout(Map, &Den[0], Averages->Rho_n);
            ScaLBL_Comm->RegularLayout(Map, &Den[Np], Averages->Rho_w);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[0], Averages->Vel_x);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[Np], Averages->Vel_y);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[2 * Np], Averages->Vel_z);

            Averages->Basic();
        }

        // allow initial ramp-up to get closer to steady state
        if (timestep > RAMP_TIMESTEPS && timestep % analysis_interval == 0 &&
            USE_MORPH) {
            //analysis.finish();
            CURRENT_STEADY_TIMESTEPS += analysis_interval;

            double muA = rhoA * (tauA - 0.5) / 3.f;
            double muB = rhoB * (tauB - 0.5) / 3.f;
            double force_mag = sqrt(Fx * Fx + Fy * Fy + Fz * Fz);
            if (force_mag == 0.0) {
                force_mag = 1.0;
            }
            double current_saturation = Averages->saturation;
            double volA = current_saturation * GreyPorosity;
            double volB = (1.0 - current_saturation) * GreyPorosity;
            double flow_rate_A = Averages->oil_flow_rate;
            double flow_rate_B = Averages->water_flow_rate;
            double Ca =
                fabs(muA * flow_rate_A + muB * flow_rate_B) / (6.0 * alpha);

            if (morph_timesteps > morph_interval) {

                bool isSteady = false;
                if ((fabs((Ca - Ca_previous) / Ca) < tolerance &&
                     CURRENT_STEADY_TIMESTEPS > MIN_STEADY_TIMESTEPS))
                    isSteady = true;
                if (CURRENT_STEADY_TIMESTEPS > MAX_STEADY_TIMESTEPS)
                    isSteady = true;
                if (RESCALE_FORCE == true && SET_CAPILLARY_NUMBER == true &&
                    CURRENT_STEADY_TIMESTEPS > RESCALE_FORCE_AFTER_TIMESTEP) {
                    RESCALE_FORCE = false;
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
                        printf("    -- adjust force by factor %.5g \n ",
                               capillary_number / Ca);
                    Averages->SetParams(rhoA, rhoB, tauA, tauB, Fx, Fy, Fz,
                                        alpha, beta, GreyPorosity);
                    greyscaleColor_db->putVector<double>("F", {Fx, Fy, Fz});
                }
                if (isSteady) {
                    MORPH_ADAPT = true;
                    CURRENT_MORPH_TIMESTEPS = 0;
                    delta_volume_target =
                        Dm->Volume * volA *
                        morph_delta; // set target volume change
                    //****** ENDPOINT ADAPTATION ********/
                    double krA_TMP = fabs(muA * flow_rate_A / force_mag);
                    double krB_TMP = fabs(muB * flow_rate_B / force_mag);
                    log_krA = log(krA_TMP);
                    if (krA_TMP < 0.0) {
                        // cannot do endpoint adaptation if kr is negative
                        log_krA = log_krA_prev;
                    } else if (krA_TMP < krB_TMP && morph_delta > 0.0) {
                        /** morphological target based on relative permeability for A **/
                        log_krA_target = log(KRA_MORPH_FACTOR * (krA_TMP));
                        slope_krA_volume = (log_krA - log_krA_prev) /
                                           (Dm->Volume * (volA - volA_prev));
                        delta_volume_target = min(
                            delta_volume_target,
                            Dm->Volume * (volA + (log_krA_target - log_krA) /
                                                     slope_krA_volume));
                        if (rank == 0) {
                            printf("    Enabling endpoint adaptation: krA = "
                                   "%.5g, krB = %.5g \n",
                                   krA_TMP, krB_TMP);
                            printf("    log(kr)=%.5g, volume=%.5g, TARGET "
                                   "log(kr)=%.5g, volume change=%.5g \n",
                                   log_krA, volA, log_krA_target,
                                   delta_volume_target / (volA * Dm->Volume));
                        }
                    }
                    log_krA_prev = log_krA;
                    volA_prev = volA;
                    //******************************** **/
                    /**  compute averages & write data **/
                    /*Averages->Full();
					Averages->Write(timestep);
					analysis.WriteVisData(timestep, current_db, *Averages, Phi, Pressure, Velocity, fq, Den );
					analysis.finish();
					*/
                    if (rank == 0) {
                        printf("** WRITE STEADY POINT *** ");
                        printf("Ca = %.5g, (previous = %.5g) \n", Ca,
                               Ca_previous);
                        double h = Dm->voxel_length;

                        // pressures
                        double pA = Averages->Oil.p;
                        double pB = Averages->Water.p;
                        double pAB = (pA - pB) / (h * 6.0 * alpha);

                        double kAeff =
                            h * h * muA * (flow_rate_A) / (force_mag);
                        double kBeff =
                            h * h * muB * (flow_rate_B) / (force_mag);

                        double viscous_pressure_drop =
                            (rhoA * volA + rhoB * volB) * force_mag;
                        double Mobility = muA / muB;

                        bool WriteHeader = false;
                        FILE *kr_log_file = fopen("relperm.csv", "r");
                        if (kr_log_file != NULL)
                            fclose(kr_log_file);
                        else
                            WriteHeader = true;
                        kr_log_file = fopen("relperm.csv", "a");
                        if (WriteHeader)
                            fprintf(kr_log_file,
                                    "timesteps sat.water eff.perm.oil "
                                    "eff.perm.water cap.pressure.norm "
                                    "pressure.drop Ca M\n");

                        fprintf(kr_log_file,
                                "%i %.5g %.5g %.5g %.5g %.5g %.5g %.5g\n",
                                CURRENT_STEADY_TIMESTEPS, current_saturation,
                                kAeff, kBeff, pAB, viscous_pressure_drop, Ca,
                                Mobility);
                        fclose(kr_log_file);

                        printf("  Measured capillary number %.5g \n ", Ca);
                    }
                    if (SET_CAPILLARY_NUMBER) {
                        Fx *= capillary_number / Ca;
                        Fy *= capillary_number / Ca;
                        Fz *= capillary_number / Ca;
                        if (force_mag > 1e-3) {
                            Fx *= 1e-3 /
                                  force_mag; // impose ceiling for stability
                            Fy *= 1e-3 / force_mag;
                            Fz *= 1e-3 / force_mag;
                        }
                        if (rank == 0)
                            printf("    -- adjust force by factor %.5g \n ",
                                   capillary_number / Ca);
                        Averages->SetParams(rhoA, rhoB, tauA, tauB, Fx, Fy, Fz,
                                            alpha, beta, GreyPorosity);
                        greyscaleColor_db->putVector<double>("F", {Fx, Fy, Fz});
                    }

                    CURRENT_STEADY_TIMESTEPS = 0;
                } else {
                    if (rank == 0) {
                        printf("** Continue to simulate steady *** \n ");
                        printf("Ca = %.5g, (previous = %.5g) \n", Ca,
                               Ca_previous);
                    }
                }
                morph_timesteps = 0;
                Ca_previous = Ca;
            }

            if (MORPH_ADAPT) {
                CURRENT_MORPH_TIMESTEPS += analysis_interval;
                if (USE_SEED) {
                    delta_volume = volA * Dm->Volume - initial_volume;
                    CURRENT_MORPH_TIMESTEPS += analysis_interval;
                    double massChange = SeedPhaseField(seed_water);
                    if (rank == 0)
                        printf("***Seed water in oil %.5g, volume change %.5g "
                               "/ %.5g ***\n",
                               massChange, delta_volume, delta_volume_target);
                }

                if ((delta_volume - delta_volume_target) / delta_volume_target >
                    0.0) {
                    MORPH_ADAPT = false;
                    CURRENT_STEADY_TIMESTEPS = 0;
                    initial_volume = volA * Dm->Volume;
                    delta_volume = 0.0;
                    if (RESCALE_FORCE_AFTER_TIMESTEP > 0)
                        RESCALE_FORCE = true;
                } else if (!(USE_DIRECT) &&
                           CURRENT_MORPH_TIMESTEPS > MAX_MORPH_TIMESTEPS) {
                    MORPH_ADAPT = false;
                    CURRENT_STEADY_TIMESTEPS = 0;
                    initial_volume = volA * Dm->Volume;
                    delta_volume = 0.0;
                    RESCALE_FORCE = true;
                    if (RESCALE_FORCE_AFTER_TIMESTEP > 0)
                        RESCALE_FORCE = true;
                }
            }
            morph_timesteps += analysis_interval;
        }
        ScaLBL_Comm->Barrier();
    }
    //analysis.finish();
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

double
ScaLBL_GreyscaleColorModel::SeedPhaseField(const double seed_water_in_oil) {
    srand(time(NULL));
    double mass_loss = 0.f;
    double count = 0.f;
    double *Aq_tmp, *Bq_tmp;

    Aq_tmp = new double[7 * Np];
    Bq_tmp = new double[7 * Np];

    ScaLBL_CopyToHost(Aq_tmp, Aq, 7 * Np * sizeof(double));
    ScaLBL_CopyToHost(Bq_tmp, Bq, 7 * Np * sizeof(double));

    for (int n = 0; n < ScaLBL_Comm->LastExterior(); n++) {
        double random_value = seed_water_in_oil * double(rand()) / RAND_MAX;
        double dA = Aq_tmp[n] + Aq_tmp[n + Np] + Aq_tmp[n + 2 * Np] +
                    Aq_tmp[n + 3 * Np] + Aq_tmp[n + 4 * Np] +
                    Aq_tmp[n + 5 * Np] + Aq_tmp[n + 6 * Np];
        double dB = Bq_tmp[n] + Bq_tmp[n + Np] + Bq_tmp[n + 2 * Np] +
                    Bq_tmp[n + 3 * Np] + Bq_tmp[n + 4 * Np] +
                    Bq_tmp[n + 5 * Np] + Bq_tmp[n + 6 * Np];
        double phase_id = (dA - dB) / (dA + dB);
        if (phase_id > 0.0) {
            Aq_tmp[n] -= 0.3333333333333333 * random_value;
            Aq_tmp[n + Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 2 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 3 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 4 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 5 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 6 * Np] -= 0.1111111111111111 * random_value;

            Bq_tmp[n] += 0.3333333333333333 * random_value;
            Bq_tmp[n + Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 2 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 3 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 4 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 5 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 6 * Np] += 0.1111111111111111 * random_value;
        }
        mass_loss += random_value * seed_water_in_oil;
    }

    for (int n = ScaLBL_Comm->FirstInterior(); n < ScaLBL_Comm->LastInterior();
         n++) {
        double random_value = seed_water_in_oil * double(rand()) / RAND_MAX;
        double dA = Aq_tmp[n] + Aq_tmp[n + Np] + Aq_tmp[n + 2 * Np] +
                    Aq_tmp[n + 3 * Np] + Aq_tmp[n + 4 * Np] +
                    Aq_tmp[n + 5 * Np] + Aq_tmp[n + 6 * Np];
        double dB = Bq_tmp[n] + Bq_tmp[n + Np] + Bq_tmp[n + 2 * Np] +
                    Bq_tmp[n + 3 * Np] + Bq_tmp[n + 4 * Np] +
                    Bq_tmp[n + 5 * Np] + Bq_tmp[n + 6 * Np];
        double phase_id = (dA - dB) / (dA + dB);
        if (phase_id > 0.0) {
            Aq_tmp[n] -= 0.3333333333333333 * random_value;
            Aq_tmp[n + Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 2 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 3 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 4 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 5 * Np] -= 0.1111111111111111 * random_value;
            Aq_tmp[n + 6 * Np] -= 0.1111111111111111 * random_value;

            Bq_tmp[n] += 0.3333333333333333 * random_value;
            Bq_tmp[n + Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 2 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 3 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 4 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 5 * Np] += 0.1111111111111111 * random_value;
            Bq_tmp[n + 6 * Np] += 0.1111111111111111 * random_value;
        }
        mass_loss += random_value * seed_water_in_oil;
    }

    count = Dm->Comm.sumReduce(count);
    mass_loss = Dm->Comm.sumReduce(mass_loss);
    if (rank == 0)
        printf("Remove mass %.5g from %.5g voxels \n", mass_loss, count);

    // Need to initialize Aq, Bq, Den, Phi directly
    //ScaLBL_CopyToDevice(Phi,phase.data(),7*Np*sizeof(double));
    ScaLBL_CopyToDevice(Aq, Aq_tmp, 7 * Np * sizeof(double));
    ScaLBL_CopyToDevice(Bq, Bq_tmp, 7 * Np * sizeof(double));

    return (mass_loss);
}

//TODO for temporary use - writing visualization files should be included in the analysis framework in the future
void ScaLBL_GreyscaleColorModel::WriteVisFiles() {
    //NOTE: write_silo is always true

    std::vector<IO::MeshDataStruct> visData;
    fillHalo<double> fillData(Dm->Comm, Dm->rank_info,
                              {Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2}, {1, 1, 1},
                              0, 1);

    auto VxVar = std::make_shared<IO::Variable>();
    auto VyVar = std::make_shared<IO::Variable>();
    auto VzVar = std::make_shared<IO::Variable>();
    auto SignDistVar = std::make_shared<IO::Variable>();
    auto PressureVar = std::make_shared<IO::Variable>();
    auto PhaseVar = std::make_shared<IO::Variable>();

    // Create the MeshDataStruct
    IO::initialize("", "silo", "false");
    visData.resize(1);
    visData[0].meshName = "domain";
    visData[0].mesh =
        std::make_shared<IO::DomainMesh>(Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2,
                                         Dm->Nz - 2, Dm->Lx, Dm->Ly, Dm->Lz);

    // create a temp data for copy from device
    DoubleArray DataTemp(Nx, Ny, Nz);

    if (vis_db->getWithDefault<bool>("save_phase_field", true)) {

        PhaseVar->name = "Phase";
        PhaseVar->type = IO::VariableType::VolumeVariable;
        PhaseVar->dim = 1;
        PhaseVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(PhaseVar);

        ASSERT(visData[0].vars[0]->name == "Phase");
        Array<double> &PhaseData = visData[0].vars[0]->data;
        ScaLBL_CopyToHost(DataTemp.data(), Phi, sizeof(double) * Nx * Ny * Nz);
        fillData.copy(DataTemp, PhaseData);
    }

    if (vis_db->getWithDefault<bool>("save_pressure", false)) {

        PressureVar->name = "Pressure";
        PressureVar->type = IO::VariableType::VolumeVariable;
        PressureVar->dim = 1;
        PressureVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(PressureVar);

        ASSERT(visData[0].vars[1]->name == "Pressure");
        Array<double> &PressData = visData[0].vars[1]->data;
        ScaLBL_Comm->RegularLayout(Map, Pressure, DataTemp);
        fillData.copy(DataTemp, PressData);
    }

    if (vis_db->getWithDefault<bool>("save_velocity", false)) {

        VxVar->name = "Velocity_x";
        VxVar->type = IO::VariableType::VolumeVariable;
        VxVar->dim = 1;
        VxVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VxVar);
        VyVar->name = "Velocity_y";
        VyVar->type = IO::VariableType::VolumeVariable;
        VyVar->dim = 1;
        VyVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VyVar);
        VzVar->name = "Velocity_z";
        VzVar->type = IO::VariableType::VolumeVariable;
        VzVar->dim = 1;
        VzVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(VzVar);

        ASSERT(visData[0].vars[2]->name == "Velocity_x");
        ASSERT(visData[0].vars[3]->name == "Velocity_y");
        ASSERT(visData[0].vars[4]->name == "Velocity_z");
        Array<double> &VelxData = visData[0].vars[2]->data;
        Array<double> &VelyData = visData[0].vars[3]->data;
        Array<double> &VelzData = visData[0].vars[4]->data;
        ScaLBL_Comm->RegularLayout(Map, &Velocity[0], DataTemp);
        fillData.copy(DataTemp, VelxData);
        ScaLBL_Comm->RegularLayout(Map, &Velocity[Np], DataTemp);
        fillData.copy(DataTemp, VelyData);
        ScaLBL_Comm->RegularLayout(Map, &Velocity[2 * Np], DataTemp);
        fillData.copy(DataTemp, VelzData);
    }

    if (vis_db->getWithDefault<bool>("save_distance", false)) {

        SignDistVar->name = "SignDist";
        SignDistVar->type = IO::VariableType::VolumeVariable;
        SignDistVar->dim = 1;
        SignDistVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(SignDistVar);

        ASSERT(visData[0].vars[5]->name == "SignDist");
        Array<double> &SignData = visData[0].vars[5]->data;
        fillData.copy(Averages->SDs, SignData);
    }

    if (vis_db->getWithDefault<bool>("write_silo", true)) {
        IO::writeData(timestep, visData, Dm->Comm);
    }

    if (vis_db->getWithDefault<bool>("save_8bit_raw", true)) {
        //TODO
        //char CurrentIDFilename[40];
        //sprintf(CurrentIDFilename,"id_t%d.raw",timestep);
        //Averages.AggregateLabels(CurrentIDFilename);
    }
}

void ScaLBL_GreyscaleColorModel::WriteDebug() {
    // Copy back final phase indicator field and convert to regular layout
    DoubleArray PhaseField(Nx, Ny, Nz);
    //ScaLBL_Comm->RegularLayout(Map,Phi,PhaseField);
    ScaLBL_CopyToHost(PhaseField.data(), Phi, sizeof(double) * N);

    FILE *OUTFILE;
    sprintf(LocalRankFilename, "Phase.%05i.raw", rank);
    OUTFILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, OUTFILE);
    fclose(OUTFILE);

    //ScaLBL_CopyToHost(PhaseField.data(), Psi, sizeof(double)*N);
    //FILE *PSIFILE;
    //sprintf(LocalRankFilename,"Psi.%05i.raw",rank);
    //PSIFILE = fopen(LocalRankFilename,"wb");
    //fwrite(PhaseField.data(),8,N,PSIFILE);
    //fclose(PSIFILE);

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

    ScaLBL_Comm->RegularLayout(Map, &Porosity_dvc[0], PhaseField);
    FILE *POROS_FILE;
    sprintf(LocalRankFilename, "Porosity.%05i.raw", rank);
    POROS_FILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, POROS_FILE);
    fclose(POROS_FILE);

    ScaLBL_Comm->RegularLayout(Map, &Permeability_dvc[0], PhaseField);
    FILE *PERM_FILE;
    sprintf(LocalRankFilename, "Permeability.%05i.raw", rank);
    PERM_FILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, PERM_FILE);
    fclose(PERM_FILE);

    //ScaLBL_Comm->RegularLayout(Map,&GreySolidGrad[0],PhaseField);
    //FILE *GreySG_X_FILE;
    //sprintf(LocalRankFilename,"GreySolidGrad_X.%05i.raw",rank);
    //GreySG_X_FILE = fopen(LocalRankFilename,"wb");
    //fwrite(PhaseField.data(),8,N,GreySG_X_FILE);
    //fclose(GreySG_X_FILE);

    //ScaLBL_Comm->RegularLayout(Map,&GreySolidGrad[Np],PhaseField);
    //FILE *GreySG_Y_FILE;
    //sprintf(LocalRankFilename,"GreySolidGrad_Y.%05i.raw",rank);
    //GreySG_Y_FILE = fopen(LocalRankFilename,"wb");
    //fwrite(PhaseField.data(),8,N,GreySG_Y_FILE);
    //fclose(GreySG_Y_FILE);

    //ScaLBL_Comm->RegularLayout(Map,&GreySolidGrad[2*Np],PhaseField);
    //FILE *GreySG_Z_FILE;
    //sprintf(LocalRankFilename,"GreySolidGrad_Z.%05i.raw",rank);
    //GreySG_Z_FILE = fopen(LocalRankFilename,"wb");
    //fwrite(PhaseField.data(),8,N,GreySG_Z_FILE);
    //fclose(GreySG_Z_FILE);

    /*	ScaLBL_Comm->RegularLayout(Map,&ColorGrad[0],PhaseField);
	FILE *CGX_FILE;
	sprintf(LocalRankFilename,"Gradient_X.%05i.raw",rank);
	CGX_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,CGX_FILE);
	fclose(CGX_FILE);

	ScaLBL_Comm->RegularLayout(Map,&ColorGrad[Np],PhaseField);
	FILE *CGY_FILE;
	sprintf(LocalRankFilename,"Gradient_Y.%05i.raw",rank);
	CGY_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,CGY_FILE);
	fclose(CGY_FILE);

	ScaLBL_Comm->RegularLayout(Map,&ColorGrad[2*Np],PhaseField);
	FILE *CGZ_FILE;
	sprintf(LocalRankFilename,"Gradient_Z.%05i.raw",rank);
	CGZ_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,CGZ_FILE);
	fclose(CGZ_FILE);
*/
}
