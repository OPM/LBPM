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
// Run the analysis, blob identification, and write restart files
#include "analysis/runAnalysis.h"
#include "analysis/analysis.h"
#include "common/Array.h"
#include "common/Communication.h"
#include "common/MPI.h"
#include "common/ScaLBL.h"
#include "models/ColorModel.h"

#include "IO/MeshDatabase.h"
#include "threadpool/thread_pool.h"

#include "ProfilerApp.h"

AnalysisType &operator|=(AnalysisType &lhs, AnalysisType rhs) {
    lhs = static_cast<AnalysisType>(
        static_cast<std::underlying_type<AnalysisType>::type>(lhs) |
        static_cast<std::underlying_type<AnalysisType>::type>(rhs));
    return lhs;
}
bool matches(AnalysisType x, AnalysisType y) {
    return (static_cast<std::underlying_type<AnalysisType>::type>(x) &
            static_cast<std::underlying_type<AnalysisType>::type>(y)) != 0;
}

// Create a shared_ptr to an array of values
template <class TYPE>
static inline std::shared_ptr<TYPE> make_shared_array(size_t N) {
    return std::shared_ptr<TYPE>(new TYPE[N],
                                 [](const TYPE *p) { delete[] p; });
}

// Helper class to write the restart file from a seperate thread
class WriteRestartWorkItem : public ThreadPool::WorkItemRet<void> {
public:
    WriteRestartWorkItem(const std::string &filename_,
                         std::shared_ptr<double> cDen_,
                         std::shared_ptr<double> cfq_, int N_)
        : filename(filename_), cfq(cfq_), cDen(cDen_), N(N_) {}
    virtual void run() {
        PROFILE_START("Save Checkpoint", 1);
        double value;
        ofstream File(filename, ios::binary);
        for (int n = 0; n < N; n++) {
            // Write the two density values
            value = cDen.get()[n];
            File.write((char *)&value, sizeof(value));
            value = cDen.get()[N + n];
            File.write((char *)&value, sizeof(value));
        }
        for (int n = 0; n < N; n++) {
            // Write the distributions
            for (int q = 0; q < 19; q++) {
                value = cfq.get()[q * N + n];
                File.write((char *)&value, sizeof(value));
            }
        }
        File.close();
        PROFILE_STOP("Save Checkpoint", 1);
    };

private:
    WriteRestartWorkItem();
    const std::string filename;
    std::shared_ptr<double> cfq, cDen;
    const int N;
};

// Helper class to compute the blob ids
typedef std::shared_ptr<std::pair<int, IntArray>> BlobIDstruct;
typedef std::shared_ptr<std::vector<BlobIDType>> BlobIDList;
static const std::string id_map_filename = "lbpm_id_map.txt";
class BlobIdentificationWorkItem1 : public ThreadPool::WorkItemRet<void> {
public:
    BlobIdentificationWorkItem1(int timestep_, int Nx_, int Ny_, int Nz_,
                                const RankInfoStruct &rank_info_,
                                std::shared_ptr<const DoubleArray> phase_,
                                const DoubleArray &dist_, BlobIDstruct last_id_,
                                BlobIDstruct new_index_, BlobIDstruct new_id_,
                                BlobIDList new_list_,
                                runAnalysis::commWrapper &&comm_)
        : timestep(timestep_), Nx(Nx_), Ny(Ny_), Nz(Nz_), rank_info(rank_info_),
          phase(phase_), dist(dist_), last_id(last_id_), new_index(new_index_),
          new_id(new_id_), new_list(new_list_), comm(std::move(comm_)) {}
    ~BlobIdentificationWorkItem1() {}
    virtual void run() {
        // Compute the global blob id and compare to the previous version
        PROFILE_START("Identify blobs", 1);
        double vF = 0.0;
        double vS = -1.0; // one voxel buffer region around solid
        IntArray &ids = new_index->second;
        new_index->first =
            ComputeGlobalBlobIDs(Nx - 2, Ny - 2, Nz - 2, rank_info, *phase,
                                 dist, vF, vS, ids, comm.comm);
        PROFILE_STOP("Identify blobs", 1);
    }

private:
    BlobIdentificationWorkItem1();
    int timestep;
    int Nx, Ny, Nz;
    const RankInfoStruct rank_info;
    std::shared_ptr<const DoubleArray> phase;
    const DoubleArray &dist;
    BlobIDstruct last_id, new_index, new_id;
    BlobIDList new_list;
    runAnalysis::commWrapper comm;
};
class BlobIdentificationWorkItem2 : public ThreadPool::WorkItemRet<void> {
public:
    BlobIdentificationWorkItem2(int timestep_, int Nx_, int Ny_, int Nz_,
                                const RankInfoStruct &rank_info_,
                                std::shared_ptr<const DoubleArray> phase_,
                                const DoubleArray &dist_, BlobIDstruct last_id_,
                                BlobIDstruct new_index_, BlobIDstruct new_id_,
                                BlobIDList new_list_,
                                runAnalysis::commWrapper &&comm_)
        : timestep(timestep_), Nx(Nx_), Ny(Ny_), Nz(Nz_), rank_info(rank_info_),
          phase(phase_), dist(dist_), last_id(last_id_), new_index(new_index_),
          new_id(new_id_), new_list(new_list_), comm(std::move(comm_)) {}
    ~BlobIdentificationWorkItem2() {}
    virtual void run() {
        // Compute the global blob id and compare to the previous version
        PROFILE_START("Identify blobs maps", 1);
        const IntArray &ids = new_index->second;
        static int max_id = -1;
        new_id->first = new_index->first;
        new_id->second = new_index->second;
        if (last_id.get() != NULL) {
            // Compute the timestep-timestep map
            const IntArray &old_ids = last_id->second;
            ID_map_struct map =
                computeIDMap(Nx, Ny, Nz, old_ids, ids, comm.comm);
            // Renumber the current timestep's ids
            getNewIDs(map, max_id, *new_list);
            renumberIDs(*new_list, new_id->second);
            writeIDMap(map, timestep, id_map_filename);
        } else {
            max_id = -1;
            ID_map_struct map(new_id->first);
            getNewIDs(map, max_id, *new_list);
            writeIDMap(map, timestep, id_map_filename);
        }
        PROFILE_STOP("Identify blobs maps", 1);
    }

private:
    BlobIdentificationWorkItem2();
    int timestep;
    int Nx, Ny, Nz;
    const RankInfoStruct rank_info;
    std::shared_ptr<const DoubleArray> phase;
    const DoubleArray &dist;
    BlobIDstruct last_id, new_index, new_id;
    BlobIDList new_list;
    runAnalysis::commWrapper comm;
};

// Helper class to write the vis file from a thread
class WriteVisWorkItem : public ThreadPool::WorkItemRet<void> {
public:
    WriteVisWorkItem(int timestep_, std::vector<IO::MeshDataStruct> &visData_,
                     TwoPhase &Avgerages_, std::array<int, 3> n_,
                     RankInfoStruct rank_info_,
                     runAnalysis::commWrapper &&comm_)
        : timestep(timestep_), visData(visData_), Averages(Avgerages_),
          n(std::move(n_)), rank_info(std::move(rank_info_)),
          comm(std::move(comm_)) {}
    ~WriteVisWorkItem() {}
    virtual void run() {
        PROFILE_START("Save Vis", 1);

        fillHalo<double> fillData(comm.comm, rank_info, n, {1, 1, 1}, 0, 1);

        ASSERT(visData[0].vars[0]->name == "phase");
        Array<double> &PhaseData = visData[0].vars[0]->data;
        fillData.copy(Averages.SDn, PhaseData);

        ASSERT(visData[0].vars[5]->name == "SignDist");
        Array<double> &SignData = visData[0].vars[5]->data;
        fillData.copy(Averages.SDs, SignData);

        ASSERT(visData[0].vars[1]->name == "Pressure");
        Array<double> &PressData = visData[0].vars[1]->data;
        fillData.copy(Averages.Press, PressData);

        ASSERT(visData[0].vars[2]->name == "Velocity_x");
        ASSERT(visData[0].vars[3]->name == "Velocity_y");
        ASSERT(visData[0].vars[4]->name == "Velocity_z");
        Array<double> &VelxData = visData[0].vars[2]->data;
        Array<double> &VelyData = visData[0].vars[3]->data;
        Array<double> &VelzData = visData[0].vars[4]->data;
        fillData.copy(Averages.Vel_x, VelxData);
        fillData.copy(Averages.Vel_y, VelyData);
        fillData.copy(Averages.Vel_z, VelzData);

        ASSERT(visData[0].vars[6]->name == "BlobID");
        Array<double> &BlobData = visData[0].vars[6]->data;
        fillData.copy(Averages.Label_NWP, BlobData);
        IO::writeData(timestep, visData, comm.comm);

        PROFILE_STOP("Save Vis", 1);
    };

private:
    WriteVisWorkItem();
    int timestep;
    std::array<int, 3> n;
    RankInfoStruct rank_info;
    std::vector<IO::MeshDataStruct> &visData;
    TwoPhase &Averages;
    runAnalysis::commWrapper comm;
};

// Helper class to write the vis file from a thread
class IOWorkItem : public ThreadPool::WorkItemRet<void> {
public:
    IOWorkItem(int timestep_, std::shared_ptr<Database> input_db_,
               std::vector<IO::MeshDataStruct> &visData_, SubPhase &Averages_,
               std::array<int, 3> n_, RankInfoStruct rank_info_,
               runAnalysis::commWrapper &&comm_)
        : timestep(timestep_), input_db(input_db_), visData(visData_),
          Averages(Averages_), n(std::move(n_)),
          rank_info(std::move(rank_info_)), comm(std::move(comm_)) {}
    ~IOWorkItem() {}
    virtual void run() {
        PROFILE_START("Save Vis", 1);

        auto color_db = input_db->getDatabase("Color");
        auto vis_db = input_db->getDatabase("Visualization");
        // int timestep = color_db->getWithDefault<int>( "timestep", 0 );

        fillHalo<double> fillData(comm.comm, rank_info, n, {1, 1, 1}, 0, 1);

        if (vis_db->getWithDefault<bool>("save_phase_field", true)) {
            ASSERT(visData[0].vars[0]->name == "phase");
            Array<double> &PhaseData = visData[0].vars[0]->data;
            fillData.copy(Averages.Phi, PhaseData);
        }

        if (vis_db->getWithDefault<bool>("save_pressure", false)) {
            ASSERT(visData[0].vars[1]->name == "Pressure");
            Array<double> &PressData = visData[0].vars[1]->data;
            fillData.copy(Averages.Pressure, PressData);
        }

        if (vis_db->getWithDefault<bool>("save_velocity", false)) {
            ASSERT(visData[0].vars[2]->name == "Velocity_x");
            ASSERT(visData[0].vars[3]->name == "Velocity_y");
            ASSERT(visData[0].vars[4]->name == "Velocity_z");
            Array<double> &VelxData = visData[0].vars[2]->data;
            Array<double> &VelyData = visData[0].vars[3]->data;
            Array<double> &VelzData = visData[0].vars[4]->data;
            fillData.copy(Averages.Vel_x, VelxData);
            fillData.copy(Averages.Vel_y, VelyData);
            fillData.copy(Averages.Vel_z, VelzData);
        }

        if (vis_db->getWithDefault<bool>("save_dissipation", false)) {
            ASSERT(visData[0].vars[5]->name == "ViscousDissipation");
            Array<double> &ViscousDissipation = visData[0].vars[5]->data;
            fillData.copy(Averages.Dissipation, ViscousDissipation);
        }

        if (vis_db->getWithDefault<bool>("save_distance", false)) {
            ASSERT(visData[0].vars[6]->name == "SignDist");
            Array<double> &SignData = visData[0].vars[6]->data;
            fillData.copy(Averages.SDs, SignData);
        }

        if (vis_db->getWithDefault<bool>("save_connected_components", false)) {
            ASSERT(visData[0].vars[7]->name == "BlobID");
            Array<double> &BlobData = visData[0].vars[7]->data;
            fillData.copy(Averages.morph_n->label, BlobData);
        }

        if (vis_db->getWithDefault<bool>("write_silo", true))
            IO::writeData(timestep, visData, comm.comm);

        if (vis_db->getWithDefault<bool>("save_8bit_raw", true)) {
            char CurrentIDFilename[40];
            sprintf(CurrentIDFilename, "id_t%d.raw", timestep);
            Averages.AggregateLabels(CurrentIDFilename);
        }

        PROFILE_STOP("Save Vis", 1);
    };

private:
    IOWorkItem();
    int timestep;
    std::array<int, 3> n;
    RankInfoStruct rank_info;
    std::shared_ptr<Database> input_db;
    std::vector<IO::MeshDataStruct> &visData;
    SubPhase &Averages;
    runAnalysis::commWrapper comm;
};

// Helper class to run the analysis from within a thread
// Note: Averages will be modified after the constructor is called
class AnalysisWorkItem : public ThreadPool::WorkItemRet<void> {
public:
    AnalysisWorkItem(AnalysisType type_, int timestep_, TwoPhase &Averages_,
                     BlobIDstruct ids, BlobIDList id_list_, double beta_)
        : type(type_), timestep(timestep_), Averages(Averages_), blob_ids(ids),
          id_list(id_list_), beta(beta_) {}
    ~AnalysisWorkItem() {}
    virtual void run() {
        Averages.NumberComponents_NWP = blob_ids->first;
        Averages.Label_NWP = blob_ids->second;
        Averages.Label_NWP_map = *id_list;
        Averages.NumberComponents_WP = 1;
        Averages.Label_WP.fill(0.0);
        if (matches(type, AnalysisType::CopyPhaseIndicator)) {
            // Averages.ColorToSignedDistance(beta,Averages.Phase,Averages.Phase_tplus);
        }
        if (matches(type, AnalysisType::ComputeAverages)) {
            PROFILE_START("Compute dist", 1);
            Averages.Initialize();
            Averages.ComputeDelPhi();
            Averages.ColorToSignedDistance(beta, Averages.Phase, Averages.SDn);
            Averages.ColorToSignedDistance(beta, Averages.Phase_tminus,
                                           Averages.Phase_tminus);
            Averages.ColorToSignedDistance(beta, Averages.Phase_tplus,
                                           Averages.Phase_tplus);
            Averages.UpdateMeshValues();
            Averages.ComputeLocal();
            Averages.Reduce();
            Averages.PrintAll(timestep);
            Averages.Initialize();
            Averages.ComponentAverages();
            Averages.SortBlobs();
            Averages.PrintComponents(timestep);
            PROFILE_STOP("Compute dist", 1);
        }
    }

private:
    AnalysisWorkItem();
    AnalysisType type;
    int timestep;
    TwoPhase &Averages;
    BlobIDstruct blob_ids;
    BlobIDList id_list;
    double beta;
};

class TCATWorkItem : public ThreadPool::WorkItemRet<void> {
public:
    TCATWorkItem(AnalysisType type_, int timestep_, TwoPhase &Averages_,
                 BlobIDstruct ids, BlobIDList id_list_, double beta_)
        : type(type_), timestep(timestep_), Averages(Averages_), blob_ids(ids),
          id_list(id_list_), beta(beta_) {}
    ~TCATWorkItem() {}
    virtual void run() {
        Averages.NumberComponents_NWP = blob_ids->first;
        Averages.Label_NWP = blob_ids->second;
        Averages.Label_NWP_map = *id_list;
        Averages.NumberComponents_WP = 1;
        Averages.Label_WP.fill(0.0);
        if (matches(type, AnalysisType::CopyPhaseIndicator)) {
            // Averages.ColorToSignedDistance(beta,Averages.Phase,Averages.Phase_tplus);
        }
        if (matches(type, AnalysisType::ComputeAverages)) {
            PROFILE_START("Compute TCAT", 1);
            Averages.Initialize();
            Averages.ComputeDelPhi();
            Averages.ColorToSignedDistance(beta, Averages.Phase, Averages.SDn);
            Averages.ColorToSignedDistance(beta, Averages.Phase_tminus,
                                           Averages.Phase_tminus);
            Averages.ColorToSignedDistance(beta, Averages.Phase_tplus,
                                           Averages.Phase_tplus);
            Averages.UpdateMeshValues();
            Averages.ComputeLocal();
            Averages.Reduce();
            Averages.PrintAll(timestep);
            PROFILE_STOP("Compute TCAT", 1);
        }
    }

private:
    TCATWorkItem();
    AnalysisType type;
    int timestep;
    TwoPhase &Averages;
    BlobIDstruct blob_ids;
    BlobIDList id_list;
    double beta;
};

class GanglionTrackingWorkItem : public ThreadPool::WorkItemRet<void> {
public:
    GanglionTrackingWorkItem(AnalysisType type_, int timestep_,
                             TwoPhase &Averages_, BlobIDstruct ids,
                             BlobIDList id_list_, double beta_)
        : type(type_), timestep(timestep_), Averages(Averages_), blob_ids(ids),
          id_list(id_list_), beta(beta_) {}
    ~GanglionTrackingWorkItem() {}
    virtual void run() {
        Averages.NumberComponents_NWP = blob_ids->first;
        Averages.Label_NWP = blob_ids->second;
        Averages.Label_NWP_map = *id_list;
        Averages.NumberComponents_WP = 1;
        Averages.Label_WP.fill(0.0);
        if (matches(type, AnalysisType::CopyPhaseIndicator)) {
            // Averages.ColorToSignedDistance(beta,Averages.Phase,Averages.Phase_tplus);
        }
        if (matches(type, AnalysisType::ComputeAverages)) {
            PROFILE_START("Compute ganglion", 1);
            Averages.Initialize();
            Averages.ComputeDelPhi();
            Averages.ColorToSignedDistance(beta, Averages.Phase, Averages.SDn);
            Averages.ColorToSignedDistance(beta, Averages.Phase_tminus,
                                           Averages.Phase_tminus);
            Averages.ColorToSignedDistance(beta, Averages.Phase_tplus,
                                           Averages.Phase_tplus);
            Averages.UpdateMeshValues();
            Averages.ComponentAverages();
            Averages.SortBlobs();
            Averages.PrintComponents(timestep);
            PROFILE_STOP("Compute ganglion", 1);
        }
    }

private:
    GanglionTrackingWorkItem();
    AnalysisType type;
    int timestep;
    TwoPhase &Averages;
    BlobIDstruct blob_ids;
    BlobIDList id_list;
    double beta;
};

class BasicWorkItem : public ThreadPool::WorkItemRet<void> {
public:
    BasicWorkItem(AnalysisType type_, int timestep_, SubPhase &Averages_)
        : type(type_), timestep(timestep_), Averages(Averages_) {}
    ~BasicWorkItem() {}
    virtual void run() {

        if (matches(type, AnalysisType::CopyPhaseIndicator)) {
            // Averages.ColorToSignedDistance(beta,Averages.Phase,Averages.Phase_tplus);
        }
        if (matches(type, AnalysisType::ComputeAverages)) {
            PROFILE_START("Compute basic averages", 1);
            Averages.Basic();
            PROFILE_STOP("Compute basic averages", 1);
        }
    }

private:
    BasicWorkItem();
    AnalysisType type;
    int timestep;
    SubPhase &Averages;
    double beta;
};

class SubphaseWorkItem : public ThreadPool::WorkItemRet<void> {
public:
    SubphaseWorkItem(AnalysisType type_, int timestep_, SubPhase &Averages_)
        : type(type_), timestep(timestep_), Averages(Averages_) {}
    ~SubphaseWorkItem() {}
    virtual void run() {

        PROFILE_START("Compute subphase", 1);
        Averages.Full();
        Averages.Write(timestep);
        PROFILE_STOP("Compute subphase", 1);
    }

private:
    SubphaseWorkItem();
    AnalysisType type;
    int timestep;
    SubPhase &Averages;
    double beta;
};

/******************************************************************
 *  MPI comm wrapper for use with analysis                         *
 ******************************************************************/
runAnalysis::commWrapper::commWrapper(int tag_, const Utilities::MPI &comm_,
                                      runAnalysis *analysis_)
    : comm(comm_), tag(tag_), analysis(analysis_) {}
runAnalysis::commWrapper::commWrapper(commWrapper &&rhs)
    : comm(rhs.comm), tag(rhs.tag), analysis(rhs.analysis) {
    rhs.tag = -1;
}
runAnalysis::commWrapper::~commWrapper() {
    if (tag == -1)
        return;
    comm.barrier();
    analysis->d_comm_used[tag] = false;
}
runAnalysis::commWrapper runAnalysis::getComm() {
    // Get a tag from root
    int tag = -1;
    if (d_rank == 0) {
        for (int i = 0; i < 1024; i++) {
            if (!d_comm_used[i]) {
                tag = i;
                break;
            }
        }
        if (tag == -1)
            ERROR("Unable to get comm");
    }
    tag = d_comm.bcast(tag, 0);
    d_comm_used[tag] = true;
    if (d_comms[tag].isNull())
        d_comms[tag] = d_comm.dup();
    return commWrapper(tag, d_comms[tag], this);
}

/******************************************************************
 *  Constructor/Destructors                                        *
 ******************************************************************/
runAnalysis::runAnalysis(std::shared_ptr<Database> input_db,
                         const RankInfoStruct &rank_info,
                         std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm,
                         std::shared_ptr<Domain> Dm, int Np, bool Regular,
                         IntArray Map)
    : d_Np(Np), d_regular(Regular), d_rank_info(rank_info), d_Map(Map),
      d_comm(Dm->Comm.dup()), d_ScaLBL_Comm(ScaLBL_Comm) {

    auto db = input_db->getDatabase("Analysis");
    auto vis_db = input_db->getDatabase("Visualization");

    /* set the I/O format */
    format = vis_db->getWithDefault<string>("format", "silo");

    // Ids of work items to use for dependencies
    ThreadPool::thread_id_t d_wait_blobID;
    ThreadPool::thread_id_t d_wait_analysis;
    ThreadPool::thread_id_t d_wait_vis;
    ThreadPool::thread_id_t d_wait_restart;
    ThreadPool::thread_id_t d_wait_subphase;

    char rankString[20];
    sprintf(rankString, "%05d", Dm->rank());
    d_n[0] = Dm->Nx - 2;
    d_n[1] = Dm->Ny - 2;
    d_n[2] = Dm->Nz - 2;
    d_N[0] = Dm->Nx;
    d_N[1] = Dm->Ny;
    d_N[2] = Dm->Nz;

    d_restart_interval = db->getWithDefault<int>("restart_interval", 100000);
    d_analysis_interval = db->getWithDefault<int>("analysis_interval", 1000);
    d_subphase_analysis_interval = INT_MAX;
    d_visualization_interval = INT_MAX;
    d_blobid_interval = INT_MAX;
    if (db->keyExists("blobid_interval")) {
        d_blobid_interval = db->getScalar<int>("blobid_interval");
    }
    if (db->keyExists("visualization_interval")) {
        d_visualization_interval = db->getScalar<int>("visualization_interval");
    }
    if (db->keyExists("subphase_analysis_interval")) {
        d_subphase_analysis_interval =
            db->getScalar<int>("subphase_analysis_interval");
    }

    auto restart_file =
        db->getWithDefault<std::string>("restart_file", "Restart");
    d_restartFile = restart_file + "." + rankString;

    d_rank = d_comm.getRank();
    writeIDMap(ID_map_struct(), 0, id_map_filename);

    // Initialize IO for silo
    //std::string  format = "silo";
    format = vis_db->getWithDefault<string>("format", "silo");

    IO::initialize("", format, "false");
    // Create the MeshDataStruct
    d_meshData.resize(1);

    d_meshData[0].meshName = "domain";
    d_meshData[0].mesh = std::make_shared<IO::DomainMesh>(
        d_rank_info, d_n[0], d_n[1], d_n[2], Dm->Lx, Dm->Ly, Dm->Lz);
    auto PhaseVar = std::make_shared<IO::Variable>();
    auto PressVar = std::make_shared<IO::Variable>();
    auto VxVar = std::make_shared<IO::Variable>();
    auto VyVar = std::make_shared<IO::Variable>();
    auto VzVar = std::make_shared<IO::Variable>();
    auto ViscousDissipationVar = std::make_shared<IO::Variable>();
    auto SignDistVar = std::make_shared<IO::Variable>();
    auto BlobIDVar = std::make_shared<IO::Variable>();

    if (vis_db->getWithDefault<bool>("save_phase_field", true)) {
        PhaseVar->name = "phase";
        PhaseVar->type = IO::VariableType::VolumeVariable;
        PhaseVar->dim = 1;
        PhaseVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(PhaseVar);
    }

    if (vis_db->getWithDefault<bool>("save_pressure", false)) {
        PressVar->name = "Pressure";
        PressVar->type = IO::VariableType::VolumeVariable;
        PressVar->dim = 1;
        PressVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(PressVar);
    }

    if (vis_db->getWithDefault<bool>("save_velocity", false)) {
        VxVar->name = "Velocity_x";
        VxVar->type = IO::VariableType::VolumeVariable;
        VxVar->dim = 1;
        VxVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(VxVar);
        VyVar->name = "Velocity_y";
        VyVar->type = IO::VariableType::VolumeVariable;
        VyVar->dim = 1;
        VyVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(VyVar);
        VzVar->name = "Velocity_z";
        VzVar->type = IO::VariableType::VolumeVariable;
        VzVar->dim = 1;
        VzVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(VzVar);
    }

    if (vis_db->getWithDefault<bool>("save_dissipation", false)) {
        ViscousDissipationVar->name = "ViscousDissipation";
        ViscousDissipationVar->type = IO::VariableType::VolumeVariable;
        ViscousDissipationVar->dim = 1;
        ViscousDissipationVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(ViscousDissipationVar);
    }

    if (vis_db->getWithDefault<bool>("save_distance", false)) {
        SignDistVar->name = "SignDist";
        SignDistVar->type = IO::VariableType::VolumeVariable;
        SignDistVar->dim = 1;
        SignDistVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(SignDistVar);
    }

    if (vis_db->getWithDefault<bool>("save_connected_components", false)) {
        BlobIDVar->name = "BlobID";
        BlobIDVar->type = IO::VariableType::VolumeVariable;
        BlobIDVar->dim = 1;
        BlobIDVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(BlobIDVar);
    }

    // Initialize the comms
    for (int i = 0; i < 1024; i++)
        d_comm_used[i] = false;
    // Initialize the threads
    int N_threads = db->getWithDefault<int>("N_threads", 4);
    auto method = db->getWithDefault<std::string>("load_balance", "default");
    createThreads(method, N_threads);
}

runAnalysis::runAnalysis(ScaLBL_ColorModel &ColorModel)
/*	std::shared_ptr<Database> input_db, const RankInfoStruct &rank_info,
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm, std::shared_ptr<Domain> Dm, int Np,
    bool Regular, IntArray Map )
    : d_Np( Np ),
      d_regular( Regular ),
      d_rank_info( rank_info ),
      d_Map( Map ),
      d_comm( Dm->Comm.dup() ),
      d_ScaLBL_Comm( ScaLBL_Comm )*/
{

    d_comm = ColorModel.Dm->Comm.dup();
    d_Np = ColorModel.Np;

    auto input_db = ColorModel.db;
    auto db = input_db->getDatabase("Analysis");
    auto vis_db = input_db->getDatabase("Visualization");

    // Ids of work items to use for dependencies
    ThreadPool::thread_id_t d_wait_blobID;
    ThreadPool::thread_id_t d_wait_analysis;
    ThreadPool::thread_id_t d_wait_vis;
    ThreadPool::thread_id_t d_wait_restart;
    ThreadPool::thread_id_t d_wait_subphase;

    char rankString[20];
    sprintf(rankString, "%05d", ColorModel.Dm->rank());
    d_n[0] = ColorModel.Dm->Nx - 2;
    d_n[1] = ColorModel.Dm->Ny - 2;
    d_n[2] = ColorModel.Dm->Nz - 2;
    d_N[0] = ColorModel.Dm->Nx;
    d_N[1] = ColorModel.Dm->Ny;
    d_N[2] = ColorModel.Dm->Nz;

    d_restart_interval = db->getWithDefault<int>("restart_interval", 100000);
    d_analysis_interval = db->getWithDefault<int>("analysis_interval", 1000);
    d_subphase_analysis_interval = INT_MAX;
    d_visualization_interval = INT_MAX;
    d_blobid_interval = INT_MAX;
    if (db->keyExists("blobid_interval")) {
        d_blobid_interval = db->getScalar<int>("blobid_interval");
    }
    if (db->keyExists("visualization_interval")) {
        d_visualization_interval = db->getScalar<int>("visualization_interval");
    }
    if (db->keyExists("subphase_analysis_interval")) {
        d_subphase_analysis_interval =
            db->getScalar<int>("subphase_analysis_interval");
    }

    auto restart_file =
        db->getWithDefault<std::string>("restart_file", "Restart");
    d_restartFile = restart_file + "." + rankString;

    d_rank = d_comm.getRank();
    writeIDMap(ID_map_struct(), 0, id_map_filename);
    // Initialize IO for silo
    //std::string  format = "silo";

    format = vis_db->getWithDefault<string>("format", "silo");

    IO::initialize("", format, "false");
    // Create the MeshDataStruct
    d_meshData.resize(1);

    d_meshData[0].meshName = "domain";
    d_meshData[0].mesh = std::make_shared<IO::DomainMesh>(
        d_rank_info, d_n[0], d_n[1], d_n[2], ColorModel.Dm->Lx,
        ColorModel.Dm->Ly, ColorModel.Dm->Lz);
    auto PhaseVar = std::make_shared<IO::Variable>();
    auto PressVar = std::make_shared<IO::Variable>();
    auto VxVar = std::make_shared<IO::Variable>();
    auto VyVar = std::make_shared<IO::Variable>();
    auto VzVar = std::make_shared<IO::Variable>();
    auto SignDistVar = std::make_shared<IO::Variable>();
    auto BlobIDVar = std::make_shared<IO::Variable>();

    if (vis_db->getWithDefault<bool>("save_phase_field", true)) {
        PhaseVar->name = "phase";
        PhaseVar->type = IO::VariableType::VolumeVariable;
        PhaseVar->dim = 1;
        PhaseVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(PhaseVar);
    }

    if (vis_db->getWithDefault<bool>("save_pressure", false)) {
        PressVar->name = "Pressure";
        PressVar->type = IO::VariableType::VolumeVariable;
        PressVar->dim = 1;
        PressVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(PressVar);
    }

    if (vis_db->getWithDefault<bool>("save_velocity", false)) {
        VxVar->name = "Velocity_x";
        VxVar->type = IO::VariableType::VolumeVariable;
        VxVar->dim = 1;
        VxVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(VxVar);
        VyVar->name = "Velocity_y";
        VyVar->type = IO::VariableType::VolumeVariable;
        VyVar->dim = 1;
        VyVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(VyVar);
        VzVar->name = "Velocity_z";
        VzVar->type = IO::VariableType::VolumeVariable;
        VzVar->dim = 1;
        VzVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(VzVar);
    }

    if (vis_db->getWithDefault<bool>("save_distance", false)) {
        SignDistVar->name = "SignDist";
        SignDistVar->type = IO::VariableType::VolumeVariable;
        SignDistVar->dim = 1;
        SignDistVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(SignDistVar);
    }

    if (vis_db->getWithDefault<bool>("save_connected_components", false)) {
        BlobIDVar->name = "BlobID";
        BlobIDVar->type = IO::VariableType::VolumeVariable;
        BlobIDVar->dim = 1;
        BlobIDVar->data.resize(d_n[0], d_n[1], d_n[2]);
        d_meshData[0].vars.push_back(BlobIDVar);
    }

    // Initialize the comms
    for (int i = 0; i < 1024; i++)
        d_comm_used[i] = false;
    // Initialize the threads
    int N_threads = db->getWithDefault<int>("N_threads", 4);
    auto method = db->getWithDefault<std::string>("load_balance", "default");
    createThreads(method, N_threads);
}
runAnalysis::~runAnalysis() {
    // Finish processing analysis
    finish();
}
void runAnalysis::finish() {
    PROFILE_START("finish");
    // Wait for the work items to finish
    d_tpool.wait_pool_finished();
    // Clear the wait ids
    d_wait_blobID.reset();
    d_wait_analysis.reset();
    d_wait_vis.reset();
    d_wait_subphase.reset();
    d_wait_restart.reset();
    // Syncronize
    d_comm.barrier();
    PROFILE_STOP("finish");
}

/******************************************************************
 *  Set the thread affinities                                      *
 ******************************************************************/
void print(const std::vector<int> &ids) {
    if (ids.empty())
        return;
    printf("%i", ids[0]);
    for (size_t i = 1; i < ids.size(); i++)
        printf(", %i", ids[i]);
    printf("\n");
}
void runAnalysis::createThreads(const std::string &method, int N_threads) {
    // Check if we are not using analysis threads
    if (method == "none")
        return;
    // Check if we have thread support
    auto thread_support = Utilities::MPI::queryThreadSupport();
    if (thread_support != Utilities::MPI::ThreadSupport::MULTIPLE &&
        N_threads > 0)
        std::cerr << "Warning: Failed to start MPI with necessary thread "
                     "support, errors may occur\n";
    // Create the threads
    const auto cores = d_tpool.getProcessAffinity();
    if (N_threads == 0) {
        // Special case to serials the analysis for debugging
        d_tpool.setNumThreads(0);
    } else if (cores.empty()) {
        // We were not able to get the cores for the process
        d_tpool.setNumThreads(N_threads);
    } else if (method == "default") {
        // Create the given number of threads, but let the OS manage affinities
        d_tpool.setNumThreads(N_threads);
    } else if (method == "independent") {
        int N = cores.size() - 1;
        d_tpool.setNumThreads(N);
        d_tpool.setThreadAffinity({cores[0]});
        for (int i = 0; i < N; i++)
            d_tpool.setThreadAffinity(i, {cores[i + 1]});
    }
    // Print the current affinities
    if (d_rank == 0) {
        printf("Affinities - rank 0:\n");
        printf("Main: ");
        print(d_tpool.getProcessAffinity());
        for (int i = 0; i < d_tpool.getNumThreads(); i++) {
            printf("Thread %i: ", i + 1);
            print(d_tpool.getThreadAffinity(i));
        }
    }
}

/******************************************************************
 *  Check which analysis we want to perform                        *
 ******************************************************************/
AnalysisType runAnalysis::computeAnalysisType(int timestep) {
    AnalysisType type = AnalysisType::AnalyzeNone;
    if (timestep % d_analysis_interval + 8 == d_analysis_interval) {
        // Copy the phase indicator field for the earlier timestep
        // printf("Copy phase indicator,timestep=%i\n",timestep);
        type |= AnalysisType::CopyPhaseIndicator;
    }
    if (timestep % d_blobid_interval == 0) {
        // Identify blobs and update global ids in time
        type |= AnalysisType::IdentifyBlobs;
    }
    /*#ifdef USE_CUDA
        if ( tpool.getQueueSize()<=3 && tpool.getNumThreads()>0 && timestep%50==0 ) {
            // Keep a few blob identifications queued up to keep the processors busy,
            // allowing us to track the blobs as fast as possible
            // Add more detailed estimates of the update frequency required to track blobs
            type |= AnalysisType::IdentifyBlobs;
        }
    #endif */
    if (timestep % d_analysis_interval + 4 == d_analysis_interval) {
        // Copy the averages to the CPU (and identify blobs)
        // printf("Copy sim state, timestep=%i \n",timestep);
        type |= AnalysisType::CopySimState;
        type |= AnalysisType::IdentifyBlobs;
    }
    if (timestep % d_analysis_interval == 0) {
        // Run the analysis
        // printf("Compute averages, timestep=%i \n",timestep);
        type |= AnalysisType::ComputeAverages;
    }
    if (timestep % d_restart_interval == 0) {
        // Write the restart file
        type |= AnalysisType::CreateRestart;
    }
    if (timestep % d_visualization_interval == 0) {
        // Write the visualization data
        type |= AnalysisType::WriteVis;
        type |= AnalysisType::CopySimState;
        type |= AnalysisType::IdentifyBlobs;
    }
    return type;
}

/******************************************************************
 *  Run the analysis                                               *
 ******************************************************************/
void runAnalysis::run(int timestep, std::shared_ptr<Database> input_db,
                      TwoPhase &Averages, const double *Phi, double *Pressure,
                      double *Velocity, double *fq, double *Den) {
    int N = d_N[0] * d_N[1] * d_N[2];
    NULL_USE(N);
    NULL_USE(Phi);

    auto db = input_db->getDatabase("Analysis");
    // int timestep = db->getWithDefault<int>( "timestep", 0 );

    // Check which analysis steps we need to perform
    auto type = computeAnalysisType(timestep);
    if (type == AnalysisType::AnalyzeNone)
        return;

    // Check how may queued items we have
    if (d_tpool.N_queued() > 20) {
        std::cerr << "Analysis queue is getting behind, waiting ...\n";
        finish();
    }

    PROFILE_START("run");

    // Copy the appropriate variables to the host (so we can spawn new threads)
    ScaLBL_DeviceBarrier();
    PROFILE_START("Copy data to host", 1);
    std::shared_ptr<DoubleArray> phase;
    /*   if ( matches(type,AnalysisType::CopyPhaseIndicator) ||
            matches(type,AnalysisType::ComputeAverages) ||
            matches(type,AnalysisType::CopySimState) ||
            matches(type,AnalysisType::IdentifyBlobs) )
    {
        phase = std::make_shared<DoubleArray>(d_N[0],d_N[1],d_N[2]);
        //ScaLBL_CopyToHost(phase->data(),Phi,N*sizeof(double));
        //   try 2 d_ScaLBL_Comm.RegulLayout(d_Map,Phi,Averages.Phase);
        //   memcpy(Averages.Phase.data(),phase->data(),N*sizeof(double));
        int Nx = d_N[0];
        int Ny = d_N[1];
        int Nz = d_N[2];
        double *TmpDat;
        TmpDat = new double [d_Np];
        ScaLBL_CopyToHost(&TmpDat[0],&Phi[0], d_Np*sizeof(double));
        for (int k=0; k<Nz; k++){
            for (int j=0; j<Ny; j++){
                for (int i=0; i<Nx; i++){
                    int n=k*Nx*Ny+j*Nx+i;
                    int idx=d_Map(i,j,k);
                    if (!(idx<0)){
                        double value=TmpDat[idx];
                        //regdata(i,j,k)=value;
                        phase->data()[n]=value;
                    }
                }
            }
        }
        delete [] TmpDat;
    }
    */
    // if ( matches(type,AnalysisType::CopyPhaseIndicator) ) {
    if (timestep % d_analysis_interval + 8 == d_analysis_interval) {
        if (d_regular)
            d_ScaLBL_Comm->RegularLayout(d_Map, Phi, Averages.Phase_tplus);
        else
            ScaLBL_CopyToHost(Averages.Phase_tplus.data(), Phi,
                              N * sizeof(double));
        // memcpy(Averages.Phase_tplus.data(),phase->data(),N*sizeof(double));
    }
    if (timestep % d_analysis_interval == 0) {
        if (d_regular)
            d_ScaLBL_Comm->RegularLayout(d_Map, Phi, Averages.Phase_tminus);
        else
            ScaLBL_CopyToHost(Averages.Phase_tminus.data(), Phi,
                              N * sizeof(double));
        // memcpy(Averages.Phase_tminus.data(),phase->data(),N*sizeof(double));
    }
    // if ( matches(type,AnalysisType::CopySimState) ) {
    if (timestep % d_analysis_interval + 4 == d_analysis_interval) {
        // Copy the members of Averages to the cpu (phase was copied above)
        PROFILE_START("Copy-Pressure", 1);
        ScaLBL_D3Q19_Pressure(fq, Pressure, d_Np);
        // ScaLBL_D3Q19_Momentum(fq,Velocity,d_Np);
        ScaLBL_DeviceBarrier();
        PROFILE_STOP("Copy-Pressure", 1);
        PROFILE_START("Copy-Wait", 1);
        PROFILE_STOP("Copy-Wait", 1);
        PROFILE_START("Copy-State", 1);
        // memcpy(Averages.Phase.data(),phase->data(),N*sizeof(double));
        if (d_regular)
            d_ScaLBL_Comm->RegularLayout(d_Map, Phi, Averages.Phase);
        else
            ScaLBL_CopyToHost(Averages.Phase.data(), Phi, N * sizeof(double));
        // copy other variables
        d_ScaLBL_Comm->RegularLayout(d_Map, Pressure, Averages.Press);
        d_ScaLBL_Comm->RegularLayout(d_Map, &Velocity[0], Averages.Vel_x);
        d_ScaLBL_Comm->RegularLayout(d_Map, &Velocity[d_Np], Averages.Vel_y);
        d_ScaLBL_Comm->RegularLayout(d_Map, &Velocity[2 * d_Np],
                                     Averages.Vel_z);
        PROFILE_STOP("Copy-State", 1);
    }
    std::shared_ptr<double> cfq, cDen;
    // if ( matches(type,AnalysisType::CreateRestart) ) {
    if (timestep % d_restart_interval == 0) {
        // Copy restart data to the CPU
        cDen = make_shared_array<double>(2 * d_Np);
        cfq = make_shared_array<double>(19 * d_Np);
        ScaLBL_CopyToHost(cfq.get(), fq, 19 * d_Np * sizeof(double));
        ScaLBL_CopyToHost(cDen.get(), Den, 2 * d_Np * sizeof(double));
    }
    PROFILE_STOP("Copy data to host", 1);

    // Spawn threads to do blob identification work
    if (matches(type, AnalysisType::IdentifyBlobs)) {
        phase = std::make_shared<DoubleArray>(d_N[0], d_N[1], d_N[2]);
        if (d_regular)
            d_ScaLBL_Comm->RegularLayout(d_Map, Phi, *phase);
        else
            ScaLBL_CopyToHost(phase->data(), Phi, N * sizeof(double));

        auto new_index =
            std::make_shared<std::pair<int, IntArray>>(0, IntArray());
        auto new_ids =
            std::make_shared<std::pair<int, IntArray>>(0, IntArray());
        auto new_list = std::make_shared<std::vector<BlobIDType>>();
        auto work1 = new BlobIdentificationWorkItem1(
            timestep, d_N[0], d_N[1], d_N[2], d_rank_info, phase, Averages.SDs,
            d_last_ids, new_index, new_ids, new_list, getComm());
        auto work2 = new BlobIdentificationWorkItem2(
            timestep, d_N[0], d_N[1], d_N[2], d_rank_info, phase, Averages.SDs,
            d_last_ids, new_index, new_ids, new_list, getComm());
        work1->add_dependency(d_wait_blobID);
        work2->add_dependency(d_tpool.add_work(work1));
        d_wait_blobID = d_tpool.add_work(work2);
        d_last_index = new_index;
        d_last_ids = new_ids;
        d_last_id_map = new_list;
    }

    // Spawn threads to do the analysis work
    // if (timestep%d_restart_interval==0){
    // if ( matches(type,AnalysisType::ComputeAverages) ) {
    if (timestep % d_analysis_interval == 0) {
        auto work = new AnalysisWorkItem(type, timestep, Averages, d_last_index,
                                         d_last_id_map, d_beta);
        work->add_dependency(d_wait_blobID);
        work->add_dependency(d_wait_analysis);
        work->add_dependency(
            d_wait_vis); // Make sure we are done using analysis before modifying
        d_wait_analysis = d_tpool.add_work(work);
    }

    // Spawn a thread to write the restart file
    //    if ( matches(type,AnalysisType::CreateRestart) ) {
    if (timestep % d_restart_interval == 0) {

        if (d_rank == 0) {
            input_db->putScalar<bool>("Restart", true);
            std::ofstream OutStream("Restart.db");
            input_db->print(OutStream, "");
            OutStream.close();
        }
        // Write the restart file (using a seperate thread)
        auto work =
            new WriteRestartWorkItem(d_restartFile.c_str(), cDen, cfq, d_Np);
        work->add_dependency(d_wait_restart);
        d_wait_restart = d_tpool.add_work(work);
    }

    // Save the results for visualization
    //    if ( matches(type,AnalysisType::CreateRestart) ) {
    if (timestep % d_restart_interval == 0) {
        // Write the vis files
        auto work = new WriteVisWorkItem(timestep, d_meshData, Averages, d_n,
                                         d_rank_info, getComm());
        work->add_dependency(d_wait_blobID);
        work->add_dependency(d_wait_analysis);
        work->add_dependency(d_wait_vis);
        d_wait_vis = d_tpool.add_work(work);
    }
    PROFILE_STOP("run");
}

/******************************************************************
 *  Run the analysis                                               *
 ******************************************************************/
void runAnalysis::basic(int timestep, std::shared_ptr<Database> input_db,
                        SubPhase &Averages, const double *Phi, double *Pressure,
                        double *Velocity, double *fq, double *Den) {
    int Nx = d_N[0];
    int Ny = d_N[1];
    int Nz = d_N[2];
    int N = Nx * Ny * Nz;
    NULL_USE(N);
    // Check which analysis steps we need to perform
    auto color_db = input_db->getDatabase("Color");
    auto vis_db = input_db->getDatabase("Visualization");

    // int timestep = color_db->getWithDefault<int>( "timestep", 0 );
    auto type = computeAnalysisType(timestep);
    if (type == AnalysisType::AnalyzeNone)
        return;

    // Check how may queued items we have
    if (d_tpool.N_queued() > 20) {
        std::cerr << "Analysis queue is getting behind, waiting ...\n";
        finish();
    }

    PROFILE_START("basic");

    // Copy the appropriate variables to the host (so we can spawn new threads)
    ScaLBL_DeviceBarrier();
    PROFILE_START("Copy data to host", 1);

    // if ( matches(type,AnalysisType::CopySimState) ) {
    if (timestep % d_analysis_interval == 0) {
        finish(); // can't copy if threads are still working on data
        // Copy the members of Averages to the cpu (phase was copied above)
        PROFILE_START("Copy-Pressure", 1);
        ScaLBL_D3Q19_Pressure(fq, Pressure, d_Np);
        // ScaLBL_D3Q19_Momentum(fq,Velocity,d_Np);
        ScaLBL_DeviceBarrier();
        PROFILE_STOP("Copy-Pressure", 1);
        PROFILE_START("Copy-Wait", 1);
        PROFILE_STOP("Copy-Wait", 1);
        PROFILE_START("Copy-State", 1);
        /*if (d_regular)
            d_ScaLBL_Comm->RegularLayout(d_Map,Phi,Averages.Phi);
        else */
        ScaLBL_CopyToHost(Averages.Phi.data(), Phi, N * sizeof(double));
        // copy other variables
        d_ScaLBL_Comm->RegularLayout(d_Map, Pressure, Averages.Pressure);
        d_ScaLBL_Comm->RegularLayout(d_Map, &Den[0], Averages.Rho_n);
        d_ScaLBL_Comm->RegularLayout(d_Map, &Den[d_Np], Averages.Rho_w);
        d_ScaLBL_Comm->RegularLayout(d_Map, &Velocity[0], Averages.Vel_x);
        d_ScaLBL_Comm->RegularLayout(d_Map, &Velocity[d_Np], Averages.Vel_y);
        d_ScaLBL_Comm->RegularLayout(d_Map, &Velocity[2 * d_Np],
                                     Averages.Vel_z);
        PROFILE_STOP("Copy-State", 1);
    }
    PROFILE_STOP("Copy data to host");

    // Spawn threads to do the analysis work
    // if (timestep%d_restart_interval==0){
    // if ( matches(type,AnalysisType::ComputeAverages) ) {
    if (timestep % d_analysis_interval == 0) {
        auto work = new BasicWorkItem(type, timestep, Averages);
        work->add_dependency(
            d_wait_subphase); // Make sure we are done using analysis before modifying
        work->add_dependency(d_wait_analysis);
        work->add_dependency(d_wait_vis);
        d_wait_analysis = d_tpool.add_work(work);
    }

    if (timestep % d_subphase_analysis_interval == 0) {
        auto work = new SubphaseWorkItem(type, timestep, Averages);
        work->add_dependency(
            d_wait_subphase); // Make sure we are done using analysis before modifying
        work->add_dependency(d_wait_analysis);
        work->add_dependency(d_wait_vis);
        d_wait_subphase = d_tpool.add_work(work);
    }

    if (timestep % d_restart_interval == 0) {
        std::shared_ptr<double> cfq, cDen;
        // Copy restart data to the CPU
        cDen = make_shared_array<double>(2 * d_Np);
        cfq = make_shared_array<double>(19 * d_Np);
        ScaLBL_CopyToHost(cfq.get(), fq, 19 * d_Np * sizeof(double));
        ScaLBL_CopyToHost(cDen.get(), Den, 2 * d_Np * sizeof(double));

        if (d_rank == 0) {
            color_db->putScalar<int>("timestep", timestep);
            color_db->putScalar<bool>("Restart", true);
            input_db->putDatabase("Color", color_db);
            std::ofstream OutStream("Restart.db");
            input_db->print(OutStream, "");
            OutStream.close();
        }
        // Write the restart file (using a seperate thread)
        auto work1 =
            new WriteRestartWorkItem(d_restartFile.c_str(), cDen, cfq, d_Np);
        work1->add_dependency(d_wait_restart);
        d_wait_restart = d_tpool.add_work(work1);
    }

    if (timestep % d_visualization_interval == 0) {
        // Write the vis files
        auto work = new IOWorkItem(timestep, input_db, d_meshData, Averages,
                                   d_n, d_rank_info, getComm());
        work->add_dependency(d_wait_analysis);
        work->add_dependency(d_wait_subphase);
        work->add_dependency(d_wait_vis);
        d_wait_vis = d_tpool.add_work(work);
    }

    PROFILE_STOP("basic");
}

void runAnalysis::WriteVisData(int timestep, std::shared_ptr<Database> input_db,
                               SubPhase &Averages, const double *Phi,
                               double *Pressure, double *Velocity, double *fq,
                               double *Den) {
    auto color_db = input_db->getDatabase("Color");
    auto vis_db = input_db->getDatabase("Visualization");
    // int timestep = color_db->getWithDefault<int>( "timestep", 0 );

    // Check which analysis steps we need to perform
    auto type = computeAnalysisType(timestep);
    if (type == AnalysisType::AnalyzeNone)
        return;

    // Check how may queued items we have
    if (d_tpool.N_queued() > 20) {
        std::cerr << "Analysis queue is getting behind, waiting ...\n";
        finish();
    }

    // Copy the appropriate variables to the host (so we can spawn new threads)
    ScaLBL_DeviceBarrier();

    PROFILE_START("write vis", 1);

    // if (Averages.WriteVis == true){
    auto work2 = new IOWorkItem(timestep, input_db, d_meshData, Averages, d_n,
                                d_rank_info, getComm());
    work2->add_dependency(d_wait_vis);
    d_wait_vis = d_tpool.add_work(work2);

    // Averages.WriteVis = false;

    PROFILE_STOP("write vis");
}
