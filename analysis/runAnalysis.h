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
#ifndef RunAnalysis_H_INC
#define RunAnalysis_H_INC

#include "analysis/SubPhase.h"
#include "analysis/TwoPhase.h"
#include "analysis/analysis.h"
#include "common/Communication.h"
#include "common/ScaLBL.h"
#include "threadpool/thread_pool.h"
#include "models/ColorModel.h"
#include <limits.h>

// Types of analysis
enum class AnalysisType : uint64_t {
    AnalyzeNone = 0,
    IdentifyBlobs = 0x01,
    CopyPhaseIndicator = 0x02,
    CopySimState = 0x04,
    ComputeAverages = 0x08,
    CreateRestart = 0x10,
    WriteVis = 0x20,
    ComputeSubphase = 0x40
};

//! Class to run the analysis in multiple threads
class runAnalysis {
public:
    //! Constructor
    runAnalysis(std::shared_ptr<Database> db, const RankInfoStruct &rank_info,
                std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm,
                std::shared_ptr<Domain> dm, int Np, bool Regular, IntArray Map);

    runAnalysis(ScaLBL_ColorModel &ColorModel);

    //! Destructor
    ~runAnalysis();

    //! Run the next analysis
    void run(int timestep, std::shared_ptr<Database> db, TwoPhase &Averages,
             const double *Phi, double *Pressure, double *Velocity, double *fq,
             double *Den);

    void basic(int timestep, std::shared_ptr<Database> db, SubPhase &Averages,
               const double *Phi, double *Pressure, double *Velocity,
               double *fq, double *Den);
    void WriteVisData(int timestep, std::shared_ptr<Database> vis_db,
                      SubPhase &Averages, const double *Phi, double *Pressure,
                      double *Velocity, double *fq, double *Den);

    //! Finish all active analysis
    void finish();

    /*!
     *  \brief    Set the affinities
     *  \details  This function will create the analysis threads and set the affinity
     *      of this thread and all analysis threads.  If MPI_THREAD_MULTIPLE is not
     *      enabled, the analysis threads will be disabled and the analysis will run in the current
     * thread.
     * @param[in] method    Method used to control the affinities:
     *                      none - Don't use threads (runs all analysis in the current thread)
     *                      default - Create the specified number of threads, but don't load balance
     *                      independent - Create the necessary number of threads to fill all cpus,
     *                                and set the affinities based on the current process such
     *                                that all threads run on independent cores
     * @param[in] N_threads Number of threads, only used by some of the methods
     */
    void createThreads(const std::string &method = "default",
                       int N_threads = 4);

private:
    runAnalysis();

    // Determine the analysis to perform
    AnalysisType computeAnalysisType(int timestep);

public:
    class commWrapper {
    public:
        Utilities::MPI comm;
        int tag;
        runAnalysis *analysis;
        commWrapper(int tag, const Utilities::MPI &comm, runAnalysis *analysis);
        commWrapper() = delete;
        commWrapper(const commWrapper &rhs) = delete;
        commWrapper &operator=(const commWrapper &rhs) = delete;
        commWrapper(commWrapper &&rhs);
        ~commWrapper();
    };

    // Get a comm (not thread safe)
    commWrapper getComm();

private:
    std::array<int, 3> d_n; // Number of local cells
    std::array<int, 3> d_N; // Number of local cells with ghosts
    int d_Np;
    int d_rank;
    int d_restart_interval, d_analysis_interval, d_blobid_interval,
        d_visualization_interval;
    int d_subphase_analysis_interval;
    double d_beta;
    bool d_regular;
    std::string format; // IO format string "silo" or "hdf5"

    ThreadPool d_tpool;
    RankInfoStruct d_rank_info;
    IntArray d_Map;
    std::shared_ptr<std::pair<int, IntArray>> d_last_ids;
    std::shared_ptr<std::pair<int, IntArray>> d_last_index;
    std::shared_ptr<std::vector<BlobIDType>> d_last_id_map;
    std::vector<IO::MeshDataStruct> d_meshData;
    std::string d_restartFile;
    Utilities::MPI d_comm;
    Utilities::MPI d_comms[1024];
    volatile bool d_comm_used[1024];
    std::shared_ptr<ScaLBL_Communicator> d_ScaLBL_Comm;

    // Ids of work items to use for dependencies
    ThreadPool::thread_id_t d_wait_blobID;
    ThreadPool::thread_id_t d_wait_analysis;
    ThreadPool::thread_id_t d_wait_subphase;
    ThreadPool::thread_id_t d_wait_vis;
    ThreadPool::thread_id_t d_wait_restart;

    // Friends
    friend commWrapper::~commWrapper();
};

#endif
