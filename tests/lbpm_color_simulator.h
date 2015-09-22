// Run the analysis, blob identification, and write restart files
#include "common/Array.h"

#define ANALYSIS_INTERVAL=1000
#define BLOBID_INTERVAL=1000

enum AnalysisType{ AnalyzeNone=0, IdentifyBlobs=0x01, CopyPhaseIndicator=0x02, 
    CopyAverages=0x04, CalcDist=0x08, CreateRestart=0x10 };


// Structure used to store ids
struct AnalysisWaitIdStruct {
    ThreadPool::thread_id_t blobID;
    ThreadPool::thread_id_t analysis;
    ThreadPool::thread_id_t restart;
};


// Helper class to write the restart file from a seperate thread
class WriteRestartWorkItem: public ThreadPool::WorkItem
{
public:
    WriteRestartWorkItem( const char* filename_, std::shared_ptr<double> cDen_,
        std::shared_ptr<double> cDistEven_, std::shared_ptr<double>cDistOdd_, int N_ ):
        filename(filename_), cDen(cDen_), cDistEven(cDistEven_), cDistOdd(cDistOdd_), N(N_) {}
    virtual void run() {
        ThreadPool::WorkItem::d_state = 1;  // Change state to in progress
        PROFILE_START("Save Checkpoint",1);
        WriteCheckpoint(filename,cDen.get(),cDistEven.get(),cDistOdd.get(),N);
        PROFILE_STOP("Save Checkpoint",1);
        PROFILE_SAVE("lbpm_color_simulator",1);
        ThreadPool::WorkItem::d_state = 2;  // Change state to finished
    };
private:
    WriteRestartWorkItem();
    const char* filename;
    std::shared_ptr<double> cDen, cDistEven, cDistOdd;
    const int N;
};


// Helper class to compute the blob ids
static const std::string id_map_filename = "lbpm_id_map.txt";
typedef std::shared_ptr<std::pair<int,IntArray> > BlobIDstruct;
typedef std::shared_ptr<std::vector<BlobIDType> > BlobIDList;
class BlobIdentificationWorkItem1: public ThreadPool::WorkItem
{
public:
    BlobIdentificationWorkItem1( int timestep_, int Nx_, int Ny_, int Nz_, const RankInfoStruct& rank_info_, 
        std::shared_ptr<const DoubleArray> phase_, const DoubleArray& dist_,
        BlobIDstruct last_id_, BlobIDstruct new_index_, BlobIDstruct new_id_, BlobIDList new_list_ ):
        timestep(timestep_), Nx(Nx_), Ny(Ny_), Nz(Nz_), rank_info(rank_info_),
        phase(phase_), dist(dist_), last_id(last_id_), new_index(new_index_), new_id(new_id_), new_list(new_list_) { }
    virtual void run() {
        ThreadPool::WorkItem::d_state = 1;  // Change state to in progress
        // Compute the global blob id and compare to the previous version
        PROFILE_START("Identify blobs",1);
        MPI_Comm newcomm;
        MPI_Comm_dup(MPI_COMM_WORLD,&newcomm);
        double vF = 0.0;
        double vS = -1.0; // one voxel buffer region around solid
        IntArray& ids = new_index->second;
        new_index->first = ComputeGlobalBlobIDs(Nx-2,Ny-2,Nz-2,rank_info,*phase,dist,vF,vS,ids,newcomm);
        MPI_Comm_free(&newcomm);
        PROFILE_STOP("Identify blobs",1);
        ThreadPool::WorkItem::d_state = 2;  // Change state to finished
    }
private:
    BlobIdentificationWorkItem1();
    int timestep;
    int Nx, Ny, Nz;
    const RankInfoStruct& rank_info;
    std::shared_ptr<const DoubleArray> phase;
    const DoubleArray& dist;
    BlobIDstruct last_id, new_index, new_id;
    BlobIDList new_list;
};
class BlobIdentificationWorkItem2: public ThreadPool::WorkItem
{
public:
    BlobIdentificationWorkItem2( int timestep_, int Nx_, int Ny_, int Nz_, const RankInfoStruct& rank_info_, 
        std::shared_ptr<const DoubleArray> phase_, const DoubleArray& dist_,
        BlobIDstruct last_id_, BlobIDstruct new_index_, BlobIDstruct new_id_, BlobIDList new_list_ ):
        timestep(timestep_), Nx(Nx_), Ny(Ny_), Nz(Nz_), rank_info(rank_info_),
        phase(phase_), dist(dist_), last_id(last_id_), new_index(new_index_), new_id(new_id_), new_list(new_list_) { }
    virtual void run() {
        ThreadPool::WorkItem::d_state = 1;  // Change state to in progress
        // Compute the global blob id and compare to the previous version
        PROFILE_START("Identify blobs maps",1);
        MPI_Comm newcomm;
        MPI_Comm_dup(MPI_COMM_WORLD,&newcomm);
        const IntArray& ids = new_index->second;
        static int max_id = -1;
        new_id->first = new_index->first;
        new_id->second = new_index->second;
        if ( last_id!=NULL ) {
            // Compute the timestep-timestep map
            const IntArray& old_ids = last_id->second;
            ID_map_struct map = computeIDMap(old_ids,ids,newcomm);
            // Renumber the current timestep's ids
            getNewIDs(map,max_id,*new_list);
            renumberIDs(*new_list,new_id->second);
            writeIDMap(map,timestep,id_map_filename);
        } else {
            max_id = -1;
            ID_map_struct map(new_id->first);
            getNewIDs(map,max_id,*new_list);
            writeIDMap(map,timestep,id_map_filename);
        }
        MPI_Comm_free(&newcomm);
        PROFILE_STOP("Identify blobs maps",1);
        ThreadPool::WorkItem::d_state = 2;  // Change state to finished
    }
private:
    BlobIdentificationWorkItem2();
    int timestep;
    int Nx, Ny, Nz;
    const RankInfoStruct& rank_info;
    std::shared_ptr<const DoubleArray> phase;
    const DoubleArray& dist;
    BlobIDstruct last_id, new_index, new_id;
    BlobIDList new_list;
};


// Helper class to run the analysis from within a thread
// Note: Averages will be modified after the constructor is called
class AnalysisWorkItem: public ThreadPool::WorkItem
{
public:
    AnalysisWorkItem( AnalysisType type_, int timestep_, TwoPhase& Averages_, 
        BlobIDstruct ids, BlobIDList id_list_, double beta_ ):
        type(type_), timestep(timestep_), Averages(Averages_), 
        blob_ids(ids), id_list(id_list_), beta(beta_) { }
    virtual void run() {
        ThreadPool::WorkItem::d_state = 1;  // Change state to in progress
        Averages.NumberComponents_NWP = blob_ids->first;
        Averages.Label_NWP = blob_ids->second;
        Averages.Label_NWP_map = *id_list;
        Averages.NumberComponents_WP = 1;
        Averages.Label_WP.fill(0.0);
        if ( (type&CopyPhaseIndicator) != 0 ) {
            // Averages.ColorToSignedDistance(beta,Averages.Phase,Averages.Phase_tplus);
        }
        if ( (type&CalcDist) != 0 ) {
            PROFILE_START("Compute dist",1);
            // Averages.ColorToSignedDistance(beta,Averages.Phase_tminus,Averages.Phase_tminus);
            Averages.Initialize();
            Averages.ComputeDelPhi();
            Averages.ColorToSignedDistance(beta,Averages.Phase,Averages.SDn);
            Averages.UpdateMeshValues();
            Averages.ComputeLocal();
            Averages.Reduce();
            Averages.PrintAll(timestep);
            Averages.Initialize();
            Averages.ComponentAverages();
            Averages.SortBlobs();
            Averages.PrintComponents(timestep);
            PROFILE_STOP("Compute dist",1);
        }
        ThreadPool::WorkItem::d_state = 2;  // Change state to finished
    }
private:
    AnalysisWorkItem();
    AnalysisType type;
    int timestep;
    TwoPhase& Averages;
    BlobIDstruct blob_ids;
    BlobIDList id_list;
    double beta;
};

// Function to start the analysis
void run_analysis( int timestep, int restart_interval, 
    const RankInfoStruct& rank_info, TwoPhase& Averages,
    BlobIDstruct& last_ids, BlobIDstruct& last_index, BlobIDList& last_id_map,
    int Nx, int Ny, int Nz, bool pBC, double beta, double err,
    const double *Phi, double *Pressure, const double *Velocity, 
    const char *ID, const double *f_even, const double *f_odd, const double *Den, 
    const char *LocalRestartFile, ThreadPool& tpool, AnalysisWaitIdStruct& wait )
{
    int N = Nx*Ny*Nz;

    // Determin the analysis we want to perform
    AnalysisType type = AnalyzeNone;
    if ( timestep%ANALYSIS_INTERVAL + 5 == ANALYSIS_INTERVAL ) {
        // Copy the phase indicator field for the earlier timestep
        type = static_cast<AnalysisType>( type | CopyPhaseIndicator );
    }
    if ( timestep%BLOBID_INTERVAL == 0 ) {
        // Identify blobs and update global ids in time
        type = static_cast<AnalysisType>( type | IdentifyBlobs );
    }
/*    #ifdef USE_CUDA
        if ( tpool.getQueueSize()<=3 && tpool.getNumThreads()>0 && timestep%50==0 ) {
            // Keep a few blob identifications queued up to keep the processors busy,
            // allowing us to track the blobs as fast as possible
            // Add more detailed estimates of the update frequency required to track blobs
            type = static_cast<AnalysisType>( type | IdentifyBlobs );
        }
    #endif
 */
    if ( timestep%ANALYSIS_INTERVAL == 0 ) {
        // Copy the averages to the CPU (and identify blobs)
        type = static_cast<AnalysisType>( type | CopyAverages );
        type = static_cast<AnalysisType>( type | IdentifyBlobs );
    }
    if ( timestep%ANALYSIS_INTERVAL == 5 ) {
        // Run the analysis
        type = static_cast<AnalysisType>( type | CalcDist );
    }
    if (timestep%restart_interval == 0) {
        // Write the restart file
        type = static_cast<AnalysisType>( type | CreateRestart );
    }
    
    // Return if we are not doing anything
    if ( type == AnalyzeNone )
        return;

    PROFILE_START("start_analysis");

    // Copy the appropriate variables to the host (so we can spawn new threads)
    DeviceBarrier();
    PROFILE_START("Copy data to host",1);
    std::shared_ptr<DoubleArray> phase;
    if ( (type&CopyPhaseIndicator)!=0 || (type&CalcDist)!=0 ||
         (type&CopyAverages)!=0 || (type&IdentifyBlobs)!=0 )
    {
        phase = std::shared_ptr<DoubleArray>(new DoubleArray(Nx,Ny,Nz));
        CopyToHost(phase->get(),Phi,N*sizeof(double));
    }
    if ( (type&CopyPhaseIndicator)!=0 ) {
        memcpy(Averages.Phase_tplus.get(),phase->get(),N*sizeof(double));
    }
    if ( (type&CalcDist)!=0 ) {
        memcpy(Averages.Phase_tminus.get(),phase->get(),N*sizeof(double));
    }
    if ( (type&CopyAverages) != 0 ) {
        // Copy the members of Averages to the cpu (phase was copied above)
        ComputePressureD3Q19(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
        memcpy(Averages.Phase.get(),phase->get(),N*sizeof(double));
        DeviceBarrier();
        CopyToHost(Averages.Press.get(),Pressure,N*sizeof(double));
        CopyToHost(Averages.Vel_x.get(),&Velocity[0],N*sizeof(double));
        CopyToHost(Averages.Vel_y.get(),&Velocity[N],N*sizeof(double));
        CopyToHost(Averages.Vel_z.get(),&Velocity[2*N],N*sizeof(double));

    }
    std::shared_ptr<double> cDen, cDistEven, cDistOdd;
    if ( (type&CreateRestart) != 0 ) {
        // Copy restart data to the CPU
        cDen = std::shared_ptr<double>(new double[2*N],DeleteArray<double>);
        cDistEven = std::shared_ptr<double>(new double[10*N],DeleteArray<double>);
        cDistOdd = std::shared_ptr<double>(new double[9*N],DeleteArray<double>);
        CopyToHost(cDistEven.get(),f_even,10*N*sizeof(double));
        CopyToHost(cDistOdd.get(),f_odd,9*N*sizeof(double));
        CopyToHost(cDen.get(),Den,2*N*sizeof(double));
    }
    PROFILE_STOP("Copy data to host",1);

    // Spawn threads to do blob identification work
    if ( (type&IdentifyBlobs)!=0 ) {
        BlobIDstruct new_index(new std::pair<int,IntArray>(0,IntArray()));
        BlobIDstruct new_ids(new std::pair<int,IntArray>(0,IntArray()));
        BlobIDList new_list(new std::vector<BlobIDType>());
        ThreadPool::WorkItem *work1 = new BlobIdentificationWorkItem1(timestep,
            Nx,Ny,Nz,rank_info,phase,Averages.SDs,last_ids,new_index,new_ids,new_list);
        ThreadPool::WorkItem *work2 = new BlobIdentificationWorkItem2(timestep,
            Nx,Ny,Nz,rank_info,phase,Averages.SDs,last_ids,new_index,new_ids,new_list);
        work1->add_dependency(wait.blobID);
        work2->add_dependency(tpool.add_work(work1));
        wait.blobID = tpool.add_work(work2);
        last_index = new_index;
        last_ids = new_ids;
        last_id_map = new_list;
    }

    // Spawn threads to do the analysis work
    if ( (type&CalcDist) != 0 ) {
        ThreadPool::WorkItem *work = new AnalysisWorkItem(
            type,timestep,Averages,last_index,last_id_map,beta);
        work->add_dependency(wait.blobID);
        work->add_dependency(wait.analysis);
        wait.analysis = tpool.add_work(work);
    }

    // Spawn a thread to write the restart file
    if ( (type&CreateRestart) != 0 ) {
        int rank = MPI_WORLD_RANK();
        if (pBC) {
            //err = fabs(sat_w - sat_w_previous);
            //sat_w_previous = sat_w;
            if (rank==0) printf("Timestep %i: change in saturation since last checkpoint is %f \n",timestep,err);
        } else {
            // Not clear yet
        }
        // Write the restart file (using a seperate thread)
        WriteRestartWorkItem *work = new WriteRestartWorkItem(LocalRestartFile,cDen,cDistEven,cDistOdd,N);
        work->add_dependency(wait.restart);
        wait.restart = tpool.add_work(work);
    }
    PROFILE_STOP("start_analysis");
}


