/*
color lattice boltzmann model
 */
#include "models/DFHModel.h"

ScaLBL_DFHModel::ScaLBL_DFHModel(int RANK, int NP, const Utilities::MPI &COMM)
    : rank(RANK), nprocs(NP), Restart(0), timestep(0), timestepMax(0), tauA(0),
      tauB(0), rhoA(0), rhoB(0), alpha(0), beta(0), Fx(0), Fy(0), Fz(0),
      flux(0), din(0), dout(0), inletA(0), inletB(0), outletA(0), outletB(0),
      Nx(0), Ny(0), Nz(0), N(0), Np(0), nprocx(0), nprocy(0), nprocz(0),
      BoundaryCondition(0), Lx(0), Ly(0), Lz(0), comm(COMM) {}
ScaLBL_DFHModel::~ScaLBL_DFHModel() {}

/*void ScaLBL_DFHModel::WriteCheckpoint(const char *FILENAME, const double *cPhi, const double *cfq, int Np)
{
    int q,n;
    double value;
    ofstream File(FILENAME,ios::binary);
    for (n=0; n<Np; n++){
        // Write the two density values
        value = cPhi[n];
        File.write((char*) &value, sizeof(value));
        // Write the even distributions
        for (q=0; q<19; q++){
            value = cfq[q*Np+n];
            File.write((char*) &value, sizeof(value));
        }
    }
    File.close();

}

void ScaLBL_DFHModel::ReadCheckpoint(char *FILENAME, double *cPhi, double *cfq, int Np)
{
    int q=0, n=0;
    double value=0;
    ifstream File(FILENAME,ios::binary);
    for (n=0; n<Np; n++){
        File.read((char*) &value, sizeof(value));
        cPhi[n] = value;
        // Read the distributions
        for (q=0; q<19; q++){
            File.read((char*) &value, sizeof(value));
            cfq[q*Np+n] = value;
        }
    }
    File.close();
}
 */

void ScaLBL_DFHModel::ReadParams(string filename) {
    // read the input database
    db = std::make_shared<Database>(filename);
    domain_db = db->getDatabase("Domain");
    color_db = db->getDatabase("Color");
    analysis_db = db->getDatabase("Analysis");

    // Color Model parameters
    timestepMax = color_db->getWithDefault<int>("timestepMax", 100);
    tauA = color_db->getWithDefault<double>("tauA", 1.0);
    tauB = color_db->getWithDefault<double>("tauB", 1.0);
    rhoA = color_db->getWithDefault<double>("rhoA", 1.0);
    rhoB = color_db->getWithDefault<double>("rhoB", 1.0);
    alpha = color_db->getWithDefault<double>("alpha", 0.001);
    beta = color_db->getWithDefault<double>("beta", 0.95);
    Restart = color_db->getWithDefault<bool>("Restart", true);
    din = color_db->getWithDefault<double>("din", 1.0);
    dout = color_db->getWithDefault<double>("dout", 1.0);
    flux = color_db->getWithDefault<double>("flux", 0.0);
    if (color_db->keyExists("F")) {
        Fx = color_db->getVector<double>("F")[0];
        Fy = color_db->getVector<double>("F")[1];
        Fz = color_db->getVector<double>("F")[2];
    }
    inletA = 1.f;
    inletB = 0.f;
    outletA = 0.f;
    outletB = 1.f;

    BoundaryCondition = domain_db->getScalar<int>("BC");
    if (color_db->keyExists("BC")) {
        BoundaryCondition = color_db->getScalar<int>("BC");
    } else if (domain_db->keyExists("BC")) {
        BoundaryCondition = domain_db->getScalar<int>("BC");
    }

    // Read domain parameters
    auto L = domain_db->getVector<double>("L");
    auto size = domain_db->getVector<int>("n");
    auto nproc = domain_db->getVector<int>("nproc");
    Nx = size[0];
    Ny = size[1];
    Nz = size[2];
    Lx = L[0];
    Ly = L[1];
    Lz = L[2];
    nprocx = nproc[0];
    nprocy = nproc[1];
    nprocz = nproc[2];

    if (BoundaryCondition == 4)
        flux =
            din *
            rhoA; // mass flux must adjust for density (see formulation for details)
}
void ScaLBL_DFHModel::SetDomain() {
    Dm = std::shared_ptr<Domain>(
        new Domain(domain_db, comm)); // full domain for analysis
    Mask = std::shared_ptr<Domain>(
        new Domain(domain_db, comm)); // mask domain removes immobile phases
    Nx += 2;
    Ny += 2;
    Nz += 2;
    N = Nx * Ny * Nz;
    id = new char[N];
    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = 1; // initialize this way
    Averages =
        std::shared_ptr<TwoPhase>(new TwoPhase(Dm)); // TwoPhase analysis object
    comm.barrier();
    Dm->CommInit();
    comm.barrier();
    rank = Dm->rank();
}

void ScaLBL_DFHModel::ReadInput() {
    //.......................................................................
    if (rank == 0)
        printf("Read input media... \n");
    //.......................................................................
    Mask->ReadIDs();
    for (int i = 0; i < Nx * Ny * Nz; i++)
        id[i] = Mask->id[i]; // save what was read

    sprintf(LocalRankString, "%05d", rank);
    sprintf(LocalRankFilename, "%s%s", "ID.", LocalRankString);
    sprintf(LocalRestartFile, "%s%s", "Restart.", LocalRankString);

    // .......... READ THE INPUT FILE .......................................
    //...........................................................................
    if (rank == 0)
        cout << "Reading in signed distance function..." << endl;
    //.......................................................................
    sprintf(LocalRankString, "%05d", rank);
    sprintf(LocalRankFilename, "%s%s", "SignDist.", LocalRankString);
    ReadBinaryFile(LocalRankFilename, Averages->SDs.data(), N);
    comm.barrier();
    if (rank == 0)
        cout << "Domain set." << endl;
}

void ScaLBL_DFHModel::AssignComponentLabels(double *phase) {
    size_t NLABELS = 0;
    char VALUE = 0;
    double AFFINITY = 0.f;

    auto LabelList = color_db->getVector<char>("ComponentLabels");
    auto AffinityList = color_db->getVector<double>("ComponentAffinity");

    NLABELS = LabelList.size();
    if (NLABELS != AffinityList.size()) {
        ERROR("Error: ComponentLabels and ComponentAffinity must be the same "
              "length! \n");
    }

    if (rank == 0) {
        printf("Components labels: %lu \n", NLABELS);
        for (unsigned int idx = 0; idx < NLABELS; idx++) {
            VALUE = LabelList[idx];
            AFFINITY = AffinityList[idx];
            printf("   label=%i, affinity=%f\n", int(VALUE), AFFINITY);
        }
    }
    // Assign the labels
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int n = k * Nx * Ny + j * Nx + i;
                VALUE = id[n];
                // Assign the affinity from the paired list
                for (unsigned int idx = 0; idx < NLABELS; idx++) {
                    //printf("rank=%i, idx=%i, value=%i, %i, \n",rank(),idx, VALUE,LabelList[idx]);
                    if (VALUE == LabelList[idx]) {
                        AFFINITY = AffinityList[idx];
                        idx = NLABELS;
                        Mask->id[n] =
                            0; // set mask to zero since this is an immobile component
                    }
                }
                phase[n] = AFFINITY;
            }
        }
    }
    // Set Dm to match Mask
    for (int i = 0; i < Nx * Ny * Nz; i++)
        Dm->id[i] = Mask->id[i];
}

void ScaLBL_DFHModel::Create() {
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

    int Npad = (Np / 16 + 2) * 16;
    if (rank == 0)
        printf("Set up memory efficient layout, %i | %i | %i \n", Np, Npad, N);
    Map.resize(Nx, Ny, Nz);
    Map.fill(-2);
    auto neighborList = new int[18 * Npad];
    Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map, neighborList,
                                              Mask->id.data(), Np, 1);
    ScaLBL_Comm->Barrier();
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
    ScaLBL_AllocateDeviceMemory((void **)&Phi, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&Pressure, sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&Velocity, 3 * sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&Gradient, 3 * sizeof(double) * Np);
    ScaLBL_AllocateDeviceMemory((void **)&SolidPotential,
                                3 * sizeof(double) * Np);

    //...........................................................................
    // Update GPU data structures
    if (rank == 0)
        printf("Setting up device map and neighbor list \n");
    // copy the neighbor list
    ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);

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
    ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int) * Np);
    ScaLBL_DeviceBarrier();
    delete[] TmpMap;
}

/********************************************************
 * AssignComponentLabels                                 *
 ********************************************************/
void ScaLBL_DFHModel::AssignSolidPotential() {
    if (rank == 0)
        printf("Computing solid interaction potential (Shan-Chen type) \n");
    double *PhaseLabel;
    PhaseLabel = new double[Nx * Ny * Nz];
    AssignComponentLabels(PhaseLabel);
    double *Tmp;
    Tmp = new double[3 * Np];
    //Averages->UpdateMeshValues(); // this computes the gradient of distance field (among other things)
    // Create the distance stencil
    // Compute solid forces based on mean field approximation
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
    double w_face = 1.0; //1.f/18.f;
    double w_edge = 0.5; //1.f/36.f;
    double w_corner = 0.f;
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
                                if (!(Mask->id[nn] > 0)) {
                                    double vec_x = double(ii - 1);
                                    double vec_y = double(jj - 1);
                                    double vec_z = double(kk - 1);
                                    double GWNS = PhaseLabel[nn];
                                    phi_x += GWNS * weight * vec_x;
                                    phi_y += GWNS * weight * vec_y;
                                    phi_z += GWNS * weight * vec_z;
                                    /*
									double GAMMA=-2.f;
									if (distval > 2.f) ALPHA=0.f; // symmetric cutoff distance                                    
									phi_x += ALPHA*exp(GAMMA*distval)*vec_x/distval;
									phi_y += ALPHA*exp(GAMMA*distval)*vec_y/distval;
									phi_z += ALPHA*exp(GAMMA*distval)*vec_z/distval;
									*/
                                }
                            }
                        }
                    }
                    Tmp[idx] = phi_x;
                    Tmp[idx + Np] = phi_y;
                    Tmp[idx + 2 * Np] = phi_z;

                    /*                        double d = Averages->SDs(n);
                                         double dx = Averages->SDs_x(n);
                                         double dy = Averages->SDs_y(n);
                                         double dz = Averages->SDs_z(n);
                                         double value=cns*exp(-bns*fabs(d))-cws*exp(-bns*fabs(d));

                 Tmp[idx] = value*dx;
                 Tmp[idx+Np] = value*dy;
                 Tmp[idx+2*Np] = value*dz;
					 */
                }
            }
        }
    }
    ScaLBL_CopyToDevice(SolidPotential, Tmp, 3 * sizeof(double) * Np);
    ScaLBL_DeviceBarrier();
    delete[] Tmp;
    delete[] Dst;

    /*
	DoubleArray Psx(Nx,Ny,Nz);
	DoubleArray Psy(Nx,Ny,Nz);
	DoubleArray Psz(Nx,Ny,Nz);
	DoubleArray Psnorm(Nx,Ny,Nz);
	ScaLBL_Comm->RegularLayout(Map,&SolidPotential[0],Psx);
	ScaLBL_Comm->RegularLayout(Map,&SolidPotential[Np],Psy);
	ScaLBL_Comm->RegularLayout(Map,&SolidPotential[2*Np],Psz);

	for (int n=0; n<N; n++) Psnorm(n) = Psx(n)*Psx(n)+Psy(n)*Psy(n)+Psz(n)*Psz(n);
	FILE *PFILE;
	sprintf(LocalRankFilename,"Potential.%05i.raw",rank);
	PFILE = fopen(LocalRankFilename,"wb");
	fwrite(Psnorm.data(),8,N,PFILE);
	fclose(PFILE);
	 */
}
void ScaLBL_DFHModel::Initialize() {
    /*
	 * This function initializes model
	 */

    AssignSolidPotential();
    int rank = Dm->rank();
    double count_wet = 0.f;
    double count_wet_global;
    double *PhaseLabel;
    PhaseLabel = new double[Nx * Ny * Nz];
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int idx = Map(i, j, k);
                int n = k * Nx * Ny + j * Nx + i;
                if (!(idx < 0)) {
                    if (Mask->id[n] == 1)
                        PhaseLabel[idx] = 1.0;
                    else {
                        PhaseLabel[idx] = -1.0;
                        count_wet += 1.f;
                    }
                }
            }
        }
    }
    count_wet_global = Dm->Comm.sumReduce(count_wet);

    if (rank == 0)
        printf("Wetting phase volume fraction =%f \n",
               count_wet_global / double(Nx * Ny * Nz * nprocs));
    // initialize phi based on PhaseLabel (include solid component labels)
    ScaLBL_CopyToDevice(Phi, PhaseLabel, Np * sizeof(double));
    //...........................................................................

    if (rank == 0)
        printf("Initializing distributions \n");
    ScaLBL_D3Q19_Init(fq, Np);

    if (Restart == true) {
        if (rank == 0) {
            printf("Reading restart file! \n");
            ifstream restart("Restart.txt");
            if (restart.is_open()) {
                restart >> timestep;
                printf("Restarting from timestep =%i \n", timestep);
            } else {
                printf("WARNING:No Restart.txt file, setting timestep=0 \n");
                timestep = 0;
            }
        }
        //MPI_Bcast(&timestep,1,MPI_INT,0,comm);
        // Read in the restart file to CPU buffers
        double *cPhi = new double[Np];
        double *cDist = new double[19 * Np];
        ifstream File(LocalRestartFile, ios::binary);
        double value;
        for (int n = 0; n < Np; n++) {
            File.read((char *)&value, sizeof(value));
            cPhi[n] = value;
            // Read the distributions
            for (int q = 0; q < 19; q++) {
                File.read((char *)&value, sizeof(value));
                cDist[q * Np + n] = value;
            }
        }
        File.close();
        // Copy the restart data to the GPU
        ScaLBL_CopyToDevice(fq, cDist, 19 * Np * sizeof(double));
        ScaLBL_CopyToDevice(Phi, cPhi, Np * sizeof(double));
        ScaLBL_DeviceBarrier();
        delete[] cPhi;
        delete[] cDist;
        comm.barrier();
    }

    if (rank == 0)
        printf("Initializing phase field \n");
    ScaLBL_DFH_Init(Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np);
    ScaLBL_DFH_Init(Phi, Den, Aq, Bq, ScaLBL_Comm->FirstInterior(),
                    ScaLBL_Comm->LastInterior(), Np);
}

void ScaLBL_DFHModel::Run() {
    int nprocs = nprocx * nprocy * nprocz;
    const RankInfoStruct rank_info(rank, nprocx, nprocy, nprocz);

    if (rank == 0)
        printf("********************************************************\n");
    if (rank == 0)
        printf("No. of timesteps: %i \n", timestepMax);
    ScaLBL_DeviceBarrier();
    comm.barrier();
    //************ MAIN ITERATION LOOP ***************************************/
    auto t1 = std::chrono::system_clock::now();
    bool Regular = true;
    PROFILE_START("Loop");
    runAnalysis analysis(analysis_db, rank_info, ScaLBL_Comm, Dm, Np, Regular,
                         Map);
    while (timestep < timestepMax) {
        //if ( rank==0 ) { printf("Running timestep %i (%i MB)\n",timestep+1,(int)(Utilities::getMemoryUsage()/1048576)); }
        PROFILE_START("Update");
        // *************ODD TIMESTEP*************
        timestep++;
        // Compute the Phase indicator field
        // Read for Aq, Bq happens in this routine (requires communication)
        ScaLBL_Comm->BiSendD3Q7AA(Aq, Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAodd_DFH(NeighborList, Aq, Bq, Den, Phi,
                              ScaLBL_Comm->FirstInterior(),
                              ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->BiRecvD3Q7AA(Aq, Bq); //WRITE INTO OPPOSITE
        ScaLBL_D3Q7_AAodd_DFH(NeighborList, Aq, Bq, Den, Phi, 0,
                              ScaLBL_Comm->LastExterior(), Np);

        // compute the gradient
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient,
                                  ScaLBL_Comm->FirstInterior(),
                                  ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->SendHalo(Phi);
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, 0,
                                  ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->RecvGrad(Phi, Gradient);

        // Perform the collision operation
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAodd_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient,
                               SolidPotential, rhoA, rhoB, tauA, tauB, alpha,
                               beta, Fx, Fy, Fz, ScaLBL_Comm->FirstInterior(),
                               ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        // Set BCs
        if (BoundaryCondition > 0) {
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
        }
        if (BoundaryCondition == 3) {
            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }
        if (BoundaryCondition == 4) {
            din =
                ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }
        ScaLBL_D3Q19_AAodd_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient,
                               SolidPotential, rhoA, rhoB, tauA, tauB, alpha,
                               beta, Fx, Fy, Fz, 0, ScaLBL_Comm->LastExterior(),
                               Np);
        ScaLBL_DeviceBarrier();
        comm.barrier();

        // *************EVEN TIMESTEP*************
        timestep++;
        // Compute the Phase indicator field
        ScaLBL_Comm->BiSendD3Q7AA(Aq, Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAeven_DFH(Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(),
                               ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->BiRecvD3Q7AA(Aq, Bq); //WRITE INTO OPPOSITE
        ScaLBL_D3Q7_AAeven_DFH(Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(),
                               Np);

        // compute the gradient
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient,
                                  ScaLBL_Comm->FirstInterior(),
                                  ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->SendHalo(Phi);
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, 0,
                                  ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->RecvGrad(Phi, Gradient);

        // Perform the collision operation
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        ScaLBL_D3Q19_AAeven_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient,
                                SolidPotential, rhoA, rhoB, tauA, tauB, alpha,
                                beta, Fx, Fy, Fz, ScaLBL_Comm->FirstInterior(),
                                ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        // Set boundary conditions
        if (BoundaryCondition > 0) {
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
        }
        if (BoundaryCondition == 3) {
            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        } else if (BoundaryCondition == 4) {
            din =
                ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }
        ScaLBL_D3Q19_AAeven_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient,
                                SolidPotential, rhoA, rhoB, tauA, tauB, alpha,
                                beta, Fx, Fy, Fz, 0,
                                ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_DeviceBarrier();
        comm.barrier();
        //************************************************************************
        comm.barrier();
        PROFILE_STOP("Update");

        // Run the analysis
        analysis.run(timestep, analysis_db, *Averages, Phi, Pressure, Velocity,
                     fq, Den);
    }
    analysis.finish();
    PROFILE_STOP("Loop");
    PROFILE_SAVE("lbpm_color_simulator", 1);
    //************************************************************************
    ScaLBL_DeviceBarrier();
    comm.barrier();
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

void ScaLBL_DFHModel::WriteDebug() {
    // Copy back final phase indicator field and convert to regular layout
    DoubleArray PhaseField(Nx, Ny, Nz);
    ScaLBL_Comm->RegularLayout(Map, Phi, PhaseField);
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
}
