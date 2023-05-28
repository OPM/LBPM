/*
 * Multi-relaxation time LBM Model
 */
#include "models/StokesModel.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"

ScaLBL_StokesModel::ScaLBL_StokesModel(int RANK, int NP,
                                       const Utilities::MPI &COMM)
    : rank(RANK), nprocs(NP), Restart(0), timestep(0), timestepMax(0), tau(0),
      Fx(0), Fy(0), Fz(0), flux(0), din(0), dout(0), mu(0), h(0), nu_phys(0),
      rho_phys(0), rho0(0), den_scale(0), time_conv(0), tolerance(0),
      epsilon0(0), epsilon0_LB(0), epsilonR(0), epsilon_LB(0),
      UseSlippingVelBC(0), Nx(0), Ny(0), Nz(0), N(0), Np(0), nprocx(0),
      nprocy(0), nprocz(0), BoundaryCondition(0), Lx(0), Ly(0), Lz(0),
      comm(COMM) {}
ScaLBL_StokesModel::~ScaLBL_StokesModel() {}

void ScaLBL_StokesModel::ReadParams(string filename, int num_iter) {
    // read the input database
    db = std::make_shared<Database>(filename);
    domain_db = db->getDatabase("Domain");
    stokes_db = db->getDatabase("Stokes");

    //------ Load number of iteration from multiphysics controller ------//
    timestepMax = num_iter;
    //-------------------------------------------------------------------//

    //---------------------- Default model parameters --------------------------//
    rho_phys = 1000.0; //by default use water density; unit [kg/m^3]
    nu_phys =
        1.004e-6; //by default use water kinematic viscosity at 20C; unit [m^2/sec]
    h = 1.0;      //image resolution;[um]
    tau = 1.0;
    mu = (tau - 0.5) / 3.0; //LB kinematic viscosity;unit [lu^2/lt]
    time_conv =
        h * h * mu /
        nu_phys; //time conversion factor from physical to LB unit; [sec/lt]
    rho0 = 1.0;  //LB density
    den_scale =
        rho_phys / rho0 * (h * h * h * 1.0e-18); //scale factor for density
    tolerance = 1.0e-8;
    Fx = Fy = 0.0;
    Fz = 1.0e-5;
    //Stokes solver also needs the following parameters for slipping velocity BC
    epsilon0 = 8.85e-12; //electric permittivity of vaccum; unit:[C/(V*m)]
    epsilon0_LB = epsilon0 * (h * 1.0e-6); //unit:[C/(V*lu)]
    epsilonR = 78.4;                     //default dielectric constant of water
    epsilon_LB = epsilon0_LB * epsilonR; //electric permittivity
    UseSlippingVelBC = false;
    //--------------------------------------------------------------------------//

    // Read domain parameters
    if (domain_db->keyExists("voxel_length")) { //default unit: um/lu
        h = domain_db->getScalar<double>("voxel_length");
    }

    // Single-fluid Navier-Stokes Model parameters
    //if (stokes_db->keyExists( "timestepMax" )){
    //	timestepMax = stokes_db->getScalar<int>( "timestepMax" );
    //}
    BoundaryCondition = 0;
    if (stokes_db->keyExists("BC")) {
        BoundaryCondition = stokes_db->getScalar<int>("BC");
    }
    if (stokes_db->keyExists("tolerance")) {
        tolerance = stokes_db->getScalar<double>("tolerance");
    }
    if (stokes_db->keyExists("tau")) {
        tau = stokes_db->getScalar<double>("tau");
    }
    if (stokes_db->keyExists("rho0")) {
        rho0 = stokes_db->getScalar<double>("rho0");
    }
    if (stokes_db->keyExists("nu_phys")) {
        nu_phys = stokes_db->getScalar<double>("nu_phys");
    }
    if (stokes_db->keyExists("rho_phys")) {
        rho_phys = stokes_db->getScalar<double>("rho_phys");
    }
    if (stokes_db->keyExists("F")) {
        Fx = stokes_db->getVector<double>("F")[0];
        Fy = stokes_db->getVector<double>("F")[1];
        Fz = stokes_db->getVector<double>("F")[2];
    }
    if (stokes_db->keyExists("Restart")) {
        Restart = stokes_db->getScalar<bool>("Restart");
    }
    if (stokes_db->keyExists("din")) {
        din = stokes_db->getScalar<double>("din");
    }
    if (stokes_db->keyExists("dout")) {
        dout = stokes_db->getScalar<double>("dout");
    }
    if (stokes_db->keyExists("flux")) {
        flux = stokes_db->getScalar<double>("flux");
    }
    if (stokes_db->keyExists("UseElectroosmoticVelocityBC")) {
        UseSlippingVelBC =
            stokes_db->getScalar<bool>("UseElectroosmoticVelocityBC");
    }
    if (stokes_db->keyExists("epsilonR")) {
        epsilonR = stokes_db->getScalar<double>("epsilonR");
    }

    // Re-calculate model parameters due to parameter read
    mu = (tau - 0.5) / 3.0;
    time_conv =
        (h * h * 1.0e-12) * mu /
        nu_phys; //time conversion factor from physical to LB unit; [sec/lt]
    den_scale =
        rho_phys / rho0 * (h * h * h * 1.0e-18); //scale factor for density
    epsilon0_LB = epsilon0 * (h * 1.0e-6);       //unit:[C/(V*lu)]
    epsilon_LB = epsilon0_LB * epsilonR;         //electric permittivity
}

void ScaLBL_StokesModel::ReadParams(string filename) {
    //NOTE the max time step is left unspecified

    // read the input database
    db = std::make_shared<Database>(filename);
    domain_db = db->getDatabase("Domain");
    stokes_db = db->getDatabase("Stokes");

    //---------------------- Default model parameters --------------------------//
    rho_phys = 1000.0; //by default use water density; unit [kg/m^3]
    nu_phys =
        1.004e-6; //by default use water kinematic viscosity at 20C; unit [m^2/sec]
    h = 1.0;      //image resolution;[um]
    tau = 1.0;
    mu = (tau - 0.5) / 3.0; //LB kinematic viscosity;unit [lu^2/lt]
    time_conv =
        h * h * mu /
        nu_phys; //time conversion factor from physical to LB unit; [sec/lt]
    rho0 = 1.0;  //LB density
    den_scale =
        rho_phys / rho0 * (h * h * h * 1.0e-18); //scale factor for density
    tolerance = 1.0e-8;
    Fx = Fy = 0.0;
    Fz = 1.0e-5;
    //Stokes solver also needs the following parameters for slipping velocity BC
    epsilon0 = 8.85e-12; //electric permittivity of vaccum; unit:[C/(V*m)]
    epsilon0_LB = epsilon0 * (h * 1.0e-6); //unit:[C/(V*lu)]
    epsilonR = 78.4;                     //default dielectric constant of water
    epsilon_LB = epsilon0_LB * epsilonR; //electric permittivity
    UseSlippingVelBC = false;
    //--------------------------------------------------------------------------//

    // Read domain parameters
    if (domain_db->keyExists("voxel_length")) { //default unit: um/lu
        h = domain_db->getScalar<double>("voxel_length");
    }

    // Single-fluid Navier-Stokes Model parameters
    //if (stokes_db->keyExists( "timestepMax" )){
    //	timestepMax = stokes_db->getScalar<int>( "timestepMax" );
    //}
    BoundaryCondition = 0;
    if (stokes_db->keyExists("BC")) {
        BoundaryCondition = stokes_db->getScalar<int>("BC");
    }
    if (stokes_db->keyExists("tolerance")) {
        tolerance = stokes_db->getScalar<double>("tolerance");
    }
    if (stokes_db->keyExists("tau")) {
        tau = stokes_db->getScalar<double>("tau");
    }
    if (stokes_db->keyExists("rho0")) {
        rho0 = stokes_db->getScalar<double>("rho0");
    }
    if (stokes_db->keyExists("nu_phys")) {
        nu_phys = stokes_db->getScalar<double>("nu_phys");
    }
    if (stokes_db->keyExists("rho_phys")) {
        rho_phys = stokes_db->getScalar<double>("rho_phys");
    }
    if (stokes_db->keyExists("F")) {
        Fx = stokes_db->getVector<double>("F")[0];
        Fy = stokes_db->getVector<double>("F")[1];
        Fz = stokes_db->getVector<double>("F")[2];
    }
    if (stokes_db->keyExists("Restart")) {
        Restart = stokes_db->getScalar<bool>("Restart");
    }
    if (stokes_db->keyExists("din")) {
        din = stokes_db->getScalar<double>("din");
    }
    if (stokes_db->keyExists("dout")) {
        dout = stokes_db->getScalar<double>("dout");
    }
    if (stokes_db->keyExists("flux")) {
        flux = stokes_db->getScalar<double>("flux");
    }
    if (stokes_db->keyExists("UseElectroosmoticVelocityBC")) {
        UseSlippingVelBC =
            stokes_db->getScalar<bool>("UseElectroosmoticVelocityBC");
    }
    if (stokes_db->keyExists("epsilonR")) {
        epsilonR = stokes_db->getScalar<double>("epsilonR");
    }

    // Re-calculate model parameters due to parameter read
    mu = (tau - 0.5) / 3.0;
    time_conv =
        (h * h * 1.0e-12) * mu /
        nu_phys; //time conversion factor from physical to LB unit; [sec/lt]
    den_scale =
        rho_phys / rho0 * (h * h * h * 1.0e-18); //scale factor for density
    epsilon0_LB = epsilon0 * (h * 1.0e-6);       //unit:[C/(V*lu)]
    epsilon_LB = epsilon0_LB * epsilonR;         //electric permittivity
}

void ScaLBL_StokesModel::SetDomain() {
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
    Distance.resize(Nx, Ny, Nz);
    Velocity_x.resize(Nx, Ny, Nz);
    Velocity_y.resize(Nx, Ny, Nz);
    Velocity_z.resize(Nx, Ny, Nz);

    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = 1; // initialize this way
    //Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm) ); // TwoPhase analysis object
    comm.barrier();
    Dm->BoundaryCondition = BoundaryCondition;
    Mask->BoundaryCondition = BoundaryCondition;
    Dm->CommInit();
    comm.barrier();

    rank = Dm->rank();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();
}

void ScaLBL_StokesModel::ReadInput() {

    sprintf(LocalRankString, "%05d", Dm->rank());
    sprintf(LocalRankFilename, "%s%s", "ID.", LocalRankString);
    sprintf(LocalRestartFile, "%s%s", "Restart.", LocalRankString);

    if (domain_db->keyExists("Filename")) {
        auto Filename = domain_db->getScalar<std::string>("Filename");
        Mask->Decomp(Filename);
    } else if (domain_db->keyExists("GridFile")) {
        // Read the local domain data
        auto input_id = readMicroCT(*domain_db, comm);
        // Fill the halo (assuming GCW of 1)
        array<int, 3> size0 = {(int)input_id.size(0), (int)input_id.size(1),
                               (int)input_id.size(2)};
        ArraySize size1 = {(size_t)Mask->Nx, (size_t)Mask->Ny,
                           (size_t)Mask->Nz};
        ASSERT((int)size1[0] == size0[0] + 2 && (int)size1[1] == size0[1] + 2 &&
               (int)size1[2] == size0[2] + 2);
        fillHalo<signed char> fill(comm, Mask->rank_info, size0, {1, 1, 1}, 0,
                                   1);
        Array<signed char> id_view;
        id_view.viewRaw(size1, Mask->id.data());
        fill.copy(input_id, id_view);
        fill.fill(id_view);
    } else {
        Mask->ReadIDs();
    }

    // Generate the signed distance map
    // Initialize the domain and communication
    Array<char> id_solid(Nx, Ny, Nz);
    // Solve for the position of the solid phase
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                // Initialize the solid phase
                if (Mask->id[n] > 0)
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
                Distance(i, j, k) = 2.0 * double(id_solid(i, j, k)) - 1.0;
            }
        }
    }
    //	MeanFilter(Averages->SDs);
    if (rank == 0)
        printf("LB Single-Fluid Solver: initialized solid phase & converting "
               "to Signed Distance function \n");
    CalcDist(Distance, id_solid, *Dm);
    if (rank == 0)
        cout << "    Domain set." << endl;
}

void ScaLBL_StokesModel::AssignZetaPotentialSolid(
    double *zeta_potential_solid) {
    size_t NLABELS = 0;
    signed char VALUE = 0;
    double AFFINITY = 0.f;

    auto LabelList = stokes_db->getVector<int>("SolidLabels");
    auto AffinityList = stokes_db->getVector<double>("ZetaPotentialSolidList");

    NLABELS = LabelList.size();
    if (NLABELS != AffinityList.size()) {
        ERROR("Error: LB Single-Fluid Solver: SolidLabels and "
              "ZetaPotentialSolidList must be the same length! \n");
    }

    double *label_count;
    double *label_count_global;
    label_count = new double[NLABELS];
    label_count_global = new double[NLABELS];

    for (size_t idx = 0; idx < NLABELS; idx++)
        label_count[idx] = 0;

    // Assign the labels
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                VALUE = Mask->id[n];
                AFFINITY = 0.f;
                // Assign the affinity from the paired list
                for (unsigned int idx = 0; idx < NLABELS; idx++) {
                    if (VALUE == LabelList[idx]) {
                        AFFINITY = AffinityList
                            [idx]; //no need to convert unit for zeta potential (i.e. volt)
                        label_count[idx] += 1.0;
                        idx = NLABELS;
                    }
                }
                zeta_potential_solid[n] = AFFINITY;
            }
        }
    }

    for (size_t idx = 0; idx < NLABELS; idx++)
        label_count_global[idx] = Dm->Comm.sumReduce(label_count[idx]);

    if (rank == 0) {
        printf("LB Single-Fluid Solver: number of solid labels: %lu \n",
               NLABELS);
        for (unsigned int idx = 0; idx < NLABELS; idx++) {
            VALUE = LabelList[idx];
            AFFINITY = AffinityList[idx];
            double volume_fraction =
                double(label_count_global[idx]) /
                double((Nx - 2) * (Ny - 2) * (Nz - 2) * nprocs);
            printf(
                "   label=%d, zeta potential=%.3g [V], volume fraction=%.2g\n",
                VALUE, AFFINITY, volume_fraction);
        }
    }
}

void ScaLBL_StokesModel::AssignSolidGrad(double *solid_grad) {
    double *Dst;
    Dst = new double[3 * 3 * 3];
    for (int kk = 0; kk < 3; kk++) {
        for (int jj = 0; jj < 3; jj++) {
            for (int ii = 0; ii < 3; ii++) {
                int index = kk * 9 + jj * 3 + ii;
                Dst[index] = sqrt(double(ii - 1) * double(ii - 1) +
                                  double(jj - 1) * double(jj - 1) +
                                  double(kk - 1) * double(kk - 1));
            }
        }
    }
    //implement a D3Q19 lattice
    double w_face = 1.0 / 18.0;
    double w_edge = 0.5 * w_face;
    double w_corner = 0.0;
    //local
    Dst[13] = 0.f;
    //faces
    Dst[4] = w_face;
    Dst[10] = w_face;
    Dst[12] = w_face;
    Dst[14] = w_face;
    Dst[16] = w_face;
    Dst[22] = w_face;
    // corners
    Dst[0] = w_corner;
    Dst[2] = w_corner;
    Dst[6] = w_corner;
    Dst[8] = w_corner;
    Dst[18] = w_corner;
    Dst[20] = w_corner;
    Dst[24] = w_corner;
    Dst[26] = w_corner;
    // edges
    Dst[1] = w_edge;
    Dst[3] = w_edge;
    Dst[5] = w_edge;
    Dst[7] = w_edge;
    Dst[9] = w_edge;
    Dst[11] = w_edge;
    Dst[15] = w_edge;
    Dst[17] = w_edge;
    Dst[19] = w_edge;
    Dst[21] = w_edge;
    Dst[23] = w_edge;
    Dst[25] = w_edge;

    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int idx = Map(i, j, k);
                if (!(idx < 0)) {
                    double phi_x = 0.f;
                    double phi_y = 0.f;
                    double phi_z = 0.f;
                    for (int kk = 0; kk < 3; kk++) {
                        for (int jj = 0; jj < 3; jj++) {
                            for (int ii = 0; ii < 3; ii++) {

                                int index = kk * 9 + jj * 3 + ii;
                                double weight = Dst[index];

                                int idi = i + ii - 1;
                                int idj = j + jj - 1;
                                int idk = k + kk - 1;

                                if (idi < 0)
                                    idi = 0;
                                if (idj < 0)
                                    idj = 0;
                                if (idk < 0)
                                    idk = 0;
                                if (!(idi < Nx))
                                    idi = Nx - 1;
                                if (!(idj < Ny))
                                    idj = Ny - 1;
                                if (!(idk < Nz))
                                    idk = Nz - 1;

                                int nn = idk * Nx * Ny + idj * Nx + idi;
                                double vec_x = double(ii - 1);
                                double vec_y = double(jj - 1);
                                double vec_z = double(kk - 1);
                                double GWNS = double(Mask->id[nn]);
                                //Since the solid unit normal vector is wanted, treat
                                //wet node as 0.0 and solid node as 1.0
                                GWNS = (GWNS > 0.0) ? 0.0 : 1.0;
                                phi_x += GWNS * weight * vec_x;
                                phi_y += GWNS * weight * vec_y;
                                phi_z += GWNS * weight * vec_z;
                            }
                        }
                    }
                    //solid_grad normalization
                    double phi_mag =
                        sqrt(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z);
                    if (phi_mag == 0.0)
                        phi_mag = 1.0;
                    solid_grad[idx + 0 * Np] = phi_x / phi_mag;
                    solid_grad[idx + 1 * Np] = phi_y / phi_mag;
                    solid_grad[idx + 2 * Np] = phi_z / phi_mag;
                }
            }
        }
    }
}

void ScaLBL_StokesModel::Create() {
    /*
	 *  This function creates the variables needed to run a LBM 
	 */
    int rank = Mask->rank();
    //.........................................................
    // Initialize communication structures in averaging domain
    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = Mask->id[i];
    Mask->CommInit();
    Np = Mask->PoreCount();
    //...........................................................................
    if (rank == 0)
        printf("LB Single-Fluid Solver: Create ScaLBL_Communicator \n");
    // Create a communicator for the device (will use optimized layout)
    // ScaLBL_Communicator ScaLBL_Comm(Mask); // original
    ScaLBL_Comm =
        std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));

    int Npad = (Np / 16 + 2) * 16;
    if (rank == 0)
        printf("LB Single-Fluid Solver: Set up memory efficient layout \n");
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
        printf("LB Single-Fluid Solver: Allocating distributions \n");
    //......................device distributions.................................
    size_t dist_mem_size = Np * sizeof(double);
    size_t neighborSize = 18 * (Np * sizeof(int));
    //...........................................................................
    ScaLBL_AllocateDeviceMemory((void **)&NeighborList, neighborSize);
    ScaLBL_AllocateDeviceMemory((void **)&fq, 19 * dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **)&Pressure, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&Velocity, 3 * sizeof(double) * Np);
    //...........................................................................
    // Update GPU data structures
    if (rank == 0)
        printf("LB Single-Fluid Solver: Setting up device map and neighbor "
               "list \n");
    // copy the neighbor list
    ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
    comm.barrier();

    if (UseSlippingVelBC == true) {
        ScaLBL_Comm->SetupBounceBackList(Map, Mask->id.data(), Np, 1);
        comm.barrier();

        //For slipping velocity BC, need zeta potential and solid unit normal vector
        ScaLBL_AllocateDeviceMemory((void **)&ZetaPotentialSolid,
                                    sizeof(double) * Nx * Ny * Nz);
        ScaLBL_AllocateDeviceMemory((void **)&SolidGrad,
                                    sizeof(double) * 3 *
                                        Np); //unit normal vector of solid nodes

        double *ZetaPotentialSolid_host;
        ZetaPotentialSolid_host = new double[Nx * Ny * Nz];
        AssignZetaPotentialSolid(ZetaPotentialSolid_host);
        double *SolidGrad_host;
        SolidGrad_host = new double[3 * Np];
        AssignSolidGrad(SolidGrad_host);
        ScaLBL_CopyToDevice(ZetaPotentialSolid, ZetaPotentialSolid_host,
                            Nx * Ny * Nz * sizeof(double));
        ScaLBL_CopyToDevice(SolidGrad, SolidGrad_host, 3 * Np * sizeof(double));
        ScaLBL_Comm->Barrier();
        delete[] ZetaPotentialSolid_host;
        delete[] SolidGrad_host;
    }
}

void ScaLBL_StokesModel::Initialize() {
    /*
	 * This function initializes model
	 */
    if (rank == 0)
        printf("LB Single-Fluid Solver: Initializing distributions \n");
    if (rank == 0)
        printf("***************************************************************"
               "*\n");
    ScaLBL_D3Q19_Init(fq, Np);

    if (rank == 0)
        printf("*****************************************************\n");
    if (rank == 0)
        printf("LB Single-Fluid Navier-Stokes Solver: \n");
    if (rank == 0)
        printf("      Time conversion factor: %.5g [sec/lt]\n", time_conv);
    if (rank == 0)
        printf("      Internal iteration: %i [lt]\n", timestepMax);
    if (rank == 0)
        printf("*****************************************************\n");
}

void ScaLBL_StokesModel::Run_Lite(double *ChargeDensity,
                                  double *ElectricField) {
    double rlx_setA = 1.0 / tau;
    double rlx_setB = 8.f * (2.f - rlx_setA) / (8.f - rlx_setA);
    timestep = 0;
    while (timestep < timestepMax) {
        //************************************************************************/
        //**************ODD TIMESTEP*************//
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAodd_StokesMRT(
            NeighborList, fq, Velocity, ChargeDensity, ElectricField, rlx_setA,
            rlx_setB, Fx, Fy, Fz, rho0, den_scale, h, time_conv, UseSlippingVelBC,
            ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
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
        ScaLBL_D3Q19_AAodd_StokesMRT(NeighborList, fq, Velocity, ChargeDensity,
                                     ElectricField, rlx_setA, rlx_setB, Fx, Fy,
                                     Fz, rho0, den_scale, h, time_conv, UseSlippingVelBC, 
                                     0, ScaLBL_Comm->LastExterior(), Np);

        if (UseSlippingVelBC == true) {
            ScaLBL_Comm->SolidSlippingVelocityBCD3Q19(
                fq, ZetaPotentialSolid, ElectricField, SolidGrad, epsilon_LB,
                1.0 / rlx_setA, rho0, den_scale, h, time_conv);
        }
        ScaLBL_Comm->Barrier();
        comm.barrier();

        //**************EVEN TIMESTEP*************//
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        ScaLBL_D3Q19_AAeven_StokesMRT(
            fq, Velocity, ChargeDensity, ElectricField, rlx_setA, rlx_setB, Fx,
            Fy, Fz, rho0, den_scale, h, time_conv, UseSlippingVelBC, 
            ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
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
        ScaLBL_D3Q19_AAeven_StokesMRT(fq, Velocity, ChargeDensity,
                                      ElectricField, rlx_setA, rlx_setB, Fx, Fy,
                                      Fz, rho0, den_scale, h, time_conv, UseSlippingVelBC, 
                                      0, ScaLBL_Comm->LastExterior(), Np);
        if (UseSlippingVelBC == true) {
            ScaLBL_Comm->SolidSlippingVelocityBCD3Q19(
                fq, ZetaPotentialSolid, ElectricField, SolidGrad, epsilon_LB,
                1.0 / rlx_setA, rho0, den_scale, h, time_conv);
        }
        ScaLBL_Comm->Barrier();
        comm.barrier();
        //************************************************************************/
    }
}

void ScaLBL_StokesModel::getVelocity(DoubleArray &Vel_x, DoubleArray &Vel_y,
                                     DoubleArray &Vel_z) {
    //get velocity in physical unit [m/sec]
    ScaLBL_D3Q19_Momentum(fq, Velocity, Np);
    ScaLBL_Comm->Barrier();
    comm.barrier();

    ScaLBL_Comm->RegularLayout(Map, &Velocity[0], Vel_x);
    Velocity_LB_to_Phys(Vel_x);
    ScaLBL_Comm->Barrier();
    comm.barrier();

    ScaLBL_Comm->RegularLayout(Map, &Velocity[Np], Vel_y);
    Velocity_LB_to_Phys(Vel_y);
    ScaLBL_Comm->Barrier();
    comm.barrier();

    ScaLBL_Comm->RegularLayout(Map, &Velocity[2 * Np], Vel_z);
    Velocity_LB_to_Phys(Vel_z);
    ScaLBL_Comm->Barrier();
    comm.barrier();
}

void ScaLBL_StokesModel::getVelocity_debug(int timestep) {
    //get velocity in physical unit [m/sec]
    ScaLBL_D3Q19_Momentum(fq, Velocity, Np);
    ScaLBL_Comm->Barrier();
    comm.barrier();

    DoubleArray PhaseField(Nx, Ny, Nz);
    ScaLBL_Comm->RegularLayout(Map, &Velocity[0], PhaseField);
    Velocity_LB_to_Phys(PhaseField);
    FILE *VELX_FILE;
    sprintf(LocalRankFilename, "Velocity_X_Time_%i.%05i.raw", timestep, rank);
    VELX_FILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, VELX_FILE);
    fclose(VELX_FILE);

    ScaLBL_Comm->RegularLayout(Map, &Velocity[Np], PhaseField);
    Velocity_LB_to_Phys(PhaseField);
    FILE *VELY_FILE;
    sprintf(LocalRankFilename, "Velocity_Y_Time_%i.%05i.raw", timestep, rank);
    VELY_FILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, VELY_FILE);
    fclose(VELY_FILE);

    ScaLBL_Comm->RegularLayout(Map, &Velocity[2 * Np], PhaseField);
    Velocity_LB_to_Phys(PhaseField);
    FILE *VELZ_FILE;
    sprintf(LocalRankFilename, "Velocity_Z_Time_%i.%05i.raw", timestep, rank);
    VELZ_FILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, VELZ_FILE);
    fclose(VELZ_FILE);
}

void ScaLBL_StokesModel::Velocity_LB_to_Phys(DoubleArray &Vel_reg) {
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = Map(i, j, k);
                if (!(idx < 0)) {
                    Vel_reg(i, j, k) =
                        Vel_reg(i, j, k) * (h * 1.0e-6) / time_conv;
                }
            }
        }
    }
}

vector<double>
ScaLBL_StokesModel::computeElectricForceAvg(double *ChargeDensity,
                                            double *ElectricField) {

    double *Ex_host;
    double *Ey_host;
    double *Ez_host;
    Ex_host = new double[Np];
    Ey_host = new double[Np];
    Ez_host = new double[Np];

    double *rhoE_host;
    rhoE_host = new double[Np];

    ScaLBL_CopyToHost(Ex_host, &ElectricField[0 * Np], Np * sizeof(double));
    ScaLBL_CopyToHost(Ey_host, &ElectricField[1 * Np], Np * sizeof(double));
    ScaLBL_CopyToHost(Ez_host, &ElectricField[2 * Np], Np * sizeof(double));
    ScaLBL_CopyToHost(rhoE_host, ChargeDensity, Np * sizeof(double));

    double count_loc = 0;
    double count;
    double Fx_avg, Fy_avg, Fz_avg; //average electric field induced force
    double Fx_loc, Fy_loc, Fz_loc;
    Fx_loc = Fy_loc = Fz_loc = 0.0;

    for (int idx = 0; idx < ScaLBL_Comm->LastExterior(); idx++) {
        Fx_loc += rhoE_host[idx] * Ex_host[idx] * (time_conv * time_conv) /
                  (h * h * 1.0e-12) / den_scale;
        Fy_loc += rhoE_host[idx] * Ey_host[idx] * (time_conv * time_conv) /
                  (h * h * 1.0e-12) / den_scale;
        Fz_loc += rhoE_host[idx] * Ez_host[idx] * (time_conv * time_conv) /
                  (h * h * 1.0e-12) / den_scale;
        count_loc += 1.0;
    }
    for (int idx = ScaLBL_Comm->FirstInterior();
         idx < ScaLBL_Comm->LastInterior(); idx++) {
        Fx_loc += rhoE_host[idx] * Ex_host[idx] * (time_conv * time_conv) /
                  (h * h * 1.0e-12) / den_scale;
        Fy_loc += rhoE_host[idx] * Ey_host[idx] * (time_conv * time_conv) /
                  (h * h * 1.0e-12) / den_scale;
        Fz_loc += rhoE_host[idx] * Ez_host[idx] * (time_conv * time_conv) /
                  (h * h * 1.0e-12) / den_scale;
        count_loc += 1.0;
    }

    Fx_avg = Dm->Comm.sumReduce(Fx_loc);
    Fy_avg = Dm->Comm.sumReduce(Fy_loc);
    Fz_avg = Dm->Comm.sumReduce(Fz_loc);
    count = Dm->Comm.sumReduce(count_loc);

    Fx_avg /= count;
    Fy_avg /= count;
    Fz_avg /= count;

    vector<double> F_avg{Fx_avg, Fy_avg, Fz_avg};

    delete[] Ex_host;
    delete[] Ey_host;
    delete[] Ez_host;
    delete[] rhoE_host;

    return F_avg;
}

double ScaLBL_StokesModel::CalVelocityConvergence(double &flow_rate_previous,
                                                  double *ChargeDensity,
                                                  double *ElectricField) {

    //-----------------------------------------------------
    ScaLBL_D3Q19_Momentum(fq, Velocity, Np);
    ScaLBL_Comm->Barrier();
    comm.barrier();
    ScaLBL_Comm->RegularLayout(Map, &Velocity[0], Velocity_x);
    ScaLBL_Comm->RegularLayout(Map, &Velocity[Np], Velocity_y);
    ScaLBL_Comm->RegularLayout(Map, &Velocity[2 * Np], Velocity_z);

    double count_loc = 0;
    double count;
    double vax, vay, vaz;
    double vax_loc, vay_loc, vaz_loc;
    vax_loc = vay_loc = vaz_loc = 0.f;
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                if (Distance(i, j, k) > 0) {
                    vax_loc += Velocity_x(i, j, k);
                    vay_loc += Velocity_y(i, j, k);
                    vaz_loc += Velocity_z(i, j, k);
                    count_loc += 1.0;
                }
            }
        }
    }
    vax = Dm->Comm.sumReduce(vax_loc);
    vay = Dm->Comm.sumReduce(vay_loc);
    vaz = Dm->Comm.sumReduce(vaz_loc);
    count = Dm->Comm.sumReduce(count_loc);

    vax /= count;
    vay /= count;
    vaz /= count;

    vector<double> Eforce;
    Eforce = computeElectricForceAvg(ChargeDensity, ElectricField);
    double TFx = Fx + Eforce[0]; //TF: total body force
    double TFy = Fy + Eforce[1];
    double TFz = Fz + Eforce[2];
    double force_mag = sqrt(TFx * TFx + TFy * TFy + TFz * TFz);
    double dir_x = TFx / force_mag;
    double dir_y = TFy / force_mag;
    double dir_z = TFz / force_mag;
    if (force_mag == 0.0) {
        // default to z direction
        dir_x = 0.0;
        dir_y = 0.0;
        dir_z = 1.0;
        force_mag = 1.0;
    }
    double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);
    double error = fabs(flow_rate - flow_rate_previous) / fabs(flow_rate);
    flow_rate_previous = flow_rate;
    //----------------------------------------------------

    //for debugging
    if (rank == 0) {
        printf("StokesModel: error: %.5g\n", error);
    }
    return error;
}

void ScaLBL_StokesModel::Run() {
    double rlx_setA = 1.0 / tau;
    double rlx_setB = 8.f * (2.f - rlx_setA) / (8.f - rlx_setA);

    Minkowski Morphology(Mask);

    if (rank == 0) {
        bool WriteHeader = false;
        FILE *log_file = fopen("Permeability.csv", "r");
        if (log_file != NULL)
            fclose(log_file);
        else
            WriteHeader = true;

        if (WriteHeader) {
            log_file = fopen("Permeability.csv", "a+");
            fprintf(log_file, "time Fx Fy Fz mu Vs As Js Xs vx vy vz k\n");
            fclose(log_file);
        }
    }

    ScaLBL_Comm->Barrier();
    comm.barrier();
    if (rank == 0)
        printf("***************************************************************"
               "*\n");
    if (rank == 0)
        printf("LB Single-Fluid Navier-Stokes Solver: timestepMax = %i\n",
               timestepMax);
    if (rank == 0)
        printf("***************************************************************"
               "*\n");
    timestep = 0;
    double error = 1.0;
    double flow_rate_previous = 0.0;
    auto t1 = std::chrono::system_clock::now();
    while (timestep < timestepMax && error > tolerance) {
        //************************************************************************/
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq, ScaLBL_Comm->FirstInterior(),
                               ScaLBL_Comm->LastInterior(), Np, rlx_setA,
                               rlx_setB, Fx, Fy, Fz);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
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
        ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(),
                               Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
        ScaLBL_Comm->Barrier();
        comm.barrier();
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        ScaLBL_D3Q19_AAeven_MRT(fq, ScaLBL_Comm->FirstInterior(),
                                ScaLBL_Comm->LastInterior(), Np, rlx_setA,
                                rlx_setB, Fx, Fy, Fz);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
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
        ScaLBL_D3Q19_AAeven_MRT(fq, 0, ScaLBL_Comm->LastExterior(), Np,
                                rlx_setA, rlx_setB, Fx, Fy, Fz);
        ScaLBL_Comm->Barrier();
        comm.barrier();
        //************************************************************************/

        if (timestep % 1000 == 0) {
            ScaLBL_D3Q19_Momentum(fq, Velocity, Np);
            ScaLBL_Comm->Barrier();
            comm.barrier();
            ScaLBL_Comm->RegularLayout(Map, &Velocity[0], Velocity_x);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[Np], Velocity_y);
            ScaLBL_Comm->RegularLayout(Map, &Velocity[2 * Np], Velocity_z);

            double count_loc = 0;
            double count;
            double vax, vay, vaz;
            double vax_loc, vay_loc, vaz_loc;
            vax_loc = vay_loc = vaz_loc = 0.f;
            for (int k = 1; k < Nz - 1; k++) {
                for (int j = 1; j < Ny - 1; j++) {
                    for (int i = 1; i < Nx - 1; i++) {
                        if (Distance(i, j, k) > 0) {
                            vax_loc += Velocity_x(i, j, k);
                            vay_loc += Velocity_y(i, j, k);
                            vaz_loc += Velocity_z(i, j, k);
                            count_loc += 1.0;
                        }
                    }
                }
            }

            vax = Dm->Comm.sumReduce(vax_loc);
            vay = Dm->Comm.sumReduce(vay_loc);
            vaz = Dm->Comm.sumReduce(vaz_loc);
            count = Dm->Comm.sumReduce(count_loc);

            vax /= count;
            vay /= count;
            vaz /= count;

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
            double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);

            error = fabs(flow_rate - flow_rate_previous) / fabs(flow_rate);
            flow_rate_previous = flow_rate;

            //if (rank==0) printf("Computing Minkowski functionals \n");
            Morphology.ComputeScalar(Distance, 0.f);
            //Morphology.PrintAll();
            double mu = (tau - 0.5) / 3.f;
            double Vs = Morphology.V();
            double As = Morphology.A();
            double Hs = Morphology.H();
            double Xs = Morphology.X();
            Vs = Dm->Comm.sumReduce(Vs);
            As = Dm->Comm.sumReduce(As);
            Hs = Dm->Comm.sumReduce(Hs);
            Xs = Dm->Comm.sumReduce(Xs);
            double h = Dm->voxel_length;
            double absperm =
                h * h * mu * Mask->Porosity() * flow_rate / force_mag;
            if (rank == 0) {
                printf("     %f\n", absperm);
                FILE *log_file = fopen("Permeability.csv", "a");
                fprintf(log_file,
                        "%i %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g "
                        "%.8g %.8g\n",
                        timestep, Fx, Fy, Fz, mu, h * h * h * Vs, h * h * As,
                        h * Hs, Xs, vax, vay, vaz, absperm);
                fclose(log_file);
            }
        }
    }
    //************************************************************************/
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
}

void ScaLBL_StokesModel::VelocityField() {

    std::vector<IO::MeshDataStruct> visData;
    fillHalo<double> fillData(Dm->Comm, Dm->rank_info,
                              {Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2}, {1, 1, 1},
                              0, 1);

    auto VxVar = std::make_shared<IO::Variable>();
    auto VyVar = std::make_shared<IO::Variable>();
    auto VzVar = std::make_shared<IO::Variable>();
    auto SignDistVar = std::make_shared<IO::Variable>();

    IO::initialize("", "silo", "false");
    // Create the MeshDataStruct
    visData.resize(1);
    visData[0].meshName = "domain";
    visData[0].mesh =
        std::make_shared<IO::DomainMesh>(Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2,
                                         Dm->Nz - 2, Dm->Lx, Dm->Ly, Dm->Lz);
    SignDistVar->name = "SignDist";
    SignDistVar->type = IO::VariableType::VolumeVariable;
    SignDistVar->dim = 1;
    SignDistVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
    visData[0].vars.push_back(SignDistVar);

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

    Array<double> &SignData = visData[0].vars[0]->data;
    Array<double> &VelxData = visData[0].vars[1]->data;
    Array<double> &VelyData = visData[0].vars[2]->data;
    Array<double> &VelzData = visData[0].vars[3]->data;

    ASSERT(visData[0].vars[0]->name == "SignDist");
    ASSERT(visData[0].vars[1]->name == "Velocity_x");
    ASSERT(visData[0].vars[2]->name == "Velocity_y");
    ASSERT(visData[0].vars[3]->name == "Velocity_z");

    fillData.copy(Distance, SignData);
    fillData.copy(Velocity_x, VelxData);
    fillData.copy(Velocity_y, VelyData);
    fillData.copy(Velocity_z, VelzData);

    IO::writeData(timestep, visData, Dm->Comm);
}
