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
 * Multi-relaxation time LBM Model
 */
#include "models/BGKModel.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"
ScaLBL_BGKModel::ScaLBL_BGKModel(int RANK, int NP, const Utilities::MPI &COMM)
    : rank(RANK), nprocs(NP), Restart(0), timestep(0), timestepMax(0), tau(0),
      Fx(0), Fy(0), Fz(0), flux(0), din(0), dout(0), mu(0), Nx(0), Ny(0), Nz(0),
      N(0), Np(0), nprocx(0), nprocy(0), nprocz(0), BoundaryCondition(0), Lx(0),
      Ly(0), Lz(0), comm(COMM) {}
ScaLBL_BGKModel::~ScaLBL_BGKModel() {}

void ScaLBL_BGKModel::ReadParams(string filename) {
    // read the input database
    db = std::make_shared<Database>(filename);
    domain_db = db->getDatabase("Domain");
    mrt_db = db->getDatabase("BGK");
    vis_db = db->getDatabase("Visualization");

    tau = 1.0;
    timestepMax = 100000;
    ANALYSIS_INTERVAL = 1000;
    tolerance = 1.0e-8;
    Fx = Fy = 0.0;
    Fz = 1.0e-5;
    dout = 1.0;
    din = 1.0;

    // Color Model parameters
    if (mrt_db->keyExists("timestepMax")) {
        timestepMax = mrt_db->getScalar<int>("timestepMax");
    }
    if (mrt_db->keyExists("analysis_interval")) {
        ANALYSIS_INTERVAL = mrt_db->getScalar<int>("analysis_interval");
    }
    if (mrt_db->keyExists("tolerance")) {
        tolerance = mrt_db->getScalar<double>("tolerance");
    }
    if (mrt_db->keyExists("tau")) {
        tau = mrt_db->getScalar<double>("tau");
    }
    if (mrt_db->keyExists("F")) {
        Fx = mrt_db->getVector<double>("F")[0];
        Fy = mrt_db->getVector<double>("F")[1];
        Fz = mrt_db->getVector<double>("F")[2];
    }
    if (mrt_db->keyExists("Restart")) {
        Restart = mrt_db->getScalar<bool>("Restart");
    }
    if (mrt_db->keyExists("din")) {
        din = mrt_db->getScalar<double>("din");
    }
    if (mrt_db->keyExists("dout")) {
        dout = mrt_db->getScalar<double>("dout");
    }
    if (mrt_db->keyExists("flux")) {
        flux = mrt_db->getScalar<double>("flux");
    }

    // Read domain parameters
    if (mrt_db->keyExists("BoundaryCondition")) {
        BoundaryCondition = mrt_db->getScalar<int>("BC");
    } else if (domain_db->keyExists("BC")) {
        BoundaryCondition = domain_db->getScalar<int>("BC");
    }

    mu = (tau - 0.5) / 3.0;
}
void ScaLBL_BGKModel::SetDomain() {
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
    Dm->CommInit();
    comm.barrier();

    rank = Dm->rank();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();
}

void ScaLBL_BGKModel::ReadInput() {

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
        printf("Initialized solid phase -- Converting to Signed Distance "
               "function \n");
    CalcDist(Distance, id_solid, *Dm);
    if (rank == 0)
        cout << "Domain set." << endl;
}

void ScaLBL_BGKModel::Create() {
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
        printf("Create ScaLBL_Communicator \n");
    // Create a communicator for the device (will use optimized layout)
    // ScaLBL_Communicator ScaLBL_Comm(Mask); // original
    ScaLBL_Comm =
        std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));

    int Npad = (Np / 16 + 2) * 16;
    if (rank == 0)
        printf("Set up memory efficient layout \n");
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
        printf("Setting up device map and neighbor list \n");
    // copy the neighbor list
    ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
    comm.barrier();
    double MLUPS = ScaLBL_Comm->GetPerformance(NeighborList, fq, Np);
    printf("  MLPUS=%f from rank %i\n", MLUPS, rank);
}

void ScaLBL_BGKModel::Initialize() {
    /*
	 * This function initializes model
	 */
    if (rank == 0)
        printf("Initializing distributions \n");
    ScaLBL_D3Q19_Init(fq, Np);
}

void ScaLBL_BGKModel::Run() {
    double rlx = 1.0 / tau;
 
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

    //.......create and start timer............
    ScaLBL_DeviceBarrier();
    comm.barrier();
    if (rank == 0)
        printf("Beginning AA timesteps, timestepMax = %i \n", timestepMax);
    if (rank == 0)
        printf("********************************************************\n");
    timestep = 0;
    double error = 1.0;
    double flow_rate_previous = 0.0;
    auto t1 = std::chrono::system_clock::now();
    while (timestep < timestepMax && error > tolerance) {
        //************************************************************************/
      /*      			timestep++;
			ScaLBL_Comm.SendD3Q19AA(dist); //READ FROM NORMAL
			ScaLBL_D3Q19_AAodd_BGK(NeighborList, dist, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np, rlx, Fx, Fy, Fz);
			ScaLBL_Comm.RecvD3Q19AA(dist); //WRITE INTO OPPOSITE
			ScaLBL_D3Q19_AAodd_BGK(NeighborList, dist, 0, ScaLBL_Comm.next, Np, rlx, Fx, Fy, Fz);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

			timestep++;
			ScaLBL_Comm.SendD3Q19AA(dist); //READ FORM NORMAL
			ScaLBL_D3Q19_AAeven_BGK(dist, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np, rlx, Fx, Fy, Fz);
			ScaLBL_Comm.RecvD3Q19AA(dist); //WRITE INTO OPPOSITE
			ScaLBL_D3Q19_AAeven_BGK(dist, 0, ScaLBL_Comm.next, Np, rlx, Fx, Fy, Fz);
			ScaLBL_DeviceBarrier(); MPI_Barrie
      */					  
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAodd_BGK(NeighborList, fq, ScaLBL_Comm->FirstInterior(),
                               ScaLBL_Comm->LastInterior(), Np, rlx, Fx, Fy, Fz);
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
        ScaLBL_D3Q19_AAodd_BGK(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(),
                               Np, rlx, Fx, Fy, Fz);
        ScaLBL_DeviceBarrier();
        comm.barrier();
        timestep++;
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        ScaLBL_D3Q19_AAeven_BGK(fq, ScaLBL_Comm->FirstInterior(),
                                ScaLBL_Comm->LastInterior(), Np, rlx, Fx, Fy, Fz);
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
        ScaLBL_D3Q19_AAeven_BGK(fq, 0, ScaLBL_Comm->LastExterior(), Np,
                                rlx, Fx, Fy, Fz);
        ScaLBL_DeviceBarrier();
        comm.barrier();
        //************************************************************************/

        if (timestep % ANALYSIS_INTERVAL == 0) {
            ScaLBL_D3Q19_Momentum(fq, Velocity, Np);
            ScaLBL_DeviceBarrier();
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

void ScaLBL_BGKModel::VelocityField() {

    auto format = vis_db->getWithDefault<string>("format", "silo");

    /*	memcpy(Morphology.SDn.data(), Distance.data(), Nx*Ny*Nz*sizeof(double));
	Morphology.Initialize();
	Morphology.UpdateMeshValues();
	Morphology.ComputeLocal();
	Morphology.Reduce();
	
	double count_loc=0;
	double count;
	double vax,vay,vaz;
	double vax_loc,vay_loc,vaz_loc;
	vax_loc = vay_loc = vaz_loc = 0.f;
	for (int n=0; n<ScaLBL_Comm->LastExterior(); n++){
		vax_loc += VELOCITY[n];
		vay_loc += VELOCITY[Np+n];
		vaz_loc += VELOCITY[2*Np+n];
		count_loc+=1.0;
	}
	
	for (int n=ScaLBL_Comm->FirstInterior(); n<ScaLBL_Comm->LastInterior(); n++){
		vax_loc += VELOCITY[n];
		vay_loc += VELOCITY[Np+n];
		vaz_loc += VELOCITY[2*Np+n];
		count_loc+=1.0;
	}
	MPI_Allreduce(&vax_loc,&vax,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	MPI_Allreduce(&vay_loc,&vay,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	MPI_Allreduce(&vaz_loc,&vaz,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	MPI_Allreduce(&count_loc,&count,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	
	vax /= count;
	vay /= count;
	vaz /= count;
	
	double mu = (tau-0.5)/3.f;
	if (rank==0) printf("Fx Fy Fz mu Vs As Js Xs vx vy vz\n");
	if (rank==0) printf("%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",Fx, Fy, Fz, mu, 
						Morphology.V(),Morphology.A(),Morphology.J(),Morphology.X(),vax,vay,vaz);
						*/
    vis_db = db->getDatabase("Visualization");
    if (vis_db->getWithDefault<bool>("write_silo", false)) {

        std::vector<IO::MeshDataStruct> visData;
        fillHalo<double> fillData(Dm->Comm, Dm->rank_info,
                                  {Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2},
                                  {1, 1, 1}, 0, 1);

        auto VxVar = std::make_shared<IO::Variable>();
        auto VyVar = std::make_shared<IO::Variable>();
        auto VzVar = std::make_shared<IO::Variable>();
        auto SignDistVar = std::make_shared<IO::Variable>();

        IO::initialize("", format, "false");
        // Create the MeshDataStruct
        visData.resize(1);
        visData[0].meshName = "domain";
        visData[0].mesh = std::make_shared<IO::DomainMesh>(
            Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2, Dm->Lx, Dm->Ly,
            Dm->Lz);
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
}
