// Run the analysis, blob identification, and write restart files


enum AnalysisType{ AnalyzeNone=0, IdentifyBlobs=0x01, CopyPhaseIndicator=0x02, 
    CopyAverages=0x02, CalcDist=0x02, CreateRestart=0x10 };


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
        PROFILE_START("Save Checkpoint",1);
        WriteCheckpoint(filename,cDen.get(),cDistEven.get(),cDistOdd.get(),N);
        PROFILE_STOP("Save Checkpoint",1);
        PROFILE_SAVE("lbpm_color_simulator",1);
    };
private:
    WriteRestartWorkItem();
    const char* filename;
    std::shared_ptr<double> cDen, cDistEven, cDistOdd;
    const int N;
};


// Helper class to compute the blob ids
typedef std::shared_ptr<std::pair<int,IntArray> > BlobIDstruct;
class BlobIdentificationWorkItem: public ThreadPool::WorkItem
{
public:
    BlobIdentificationWorkItem( int Nx_, int Ny_, int Nz_, const RankInfoStruct& rank_info_, 
        std::shared_ptr<const double> phase_, const DoubleArray& dist_,
        BlobIDstruct last_id_, BlobIDstruct new_id_ ):
        Nx(Nx_), Ny(Ny_), Nz(Nz_), rank_info(rank_info_), phase(phase_), 
        dist(dist_), last_id(last_id_), new_id(new_id_) { }
    virtual void run() {
        // Compute the global blob id and compare to the previous version
        PROFILE_START("Identify blobs and maps",1);
        double vF = 0.0;
        double vS = 0.0;
        new_id->first = ComputeGlobalBlobIDs(Nx-2,Ny-2,Nz-2,rank_info,
            *phase,dist,vF,vS,new_id->second);
        if ( last_id==NULL ) {
            // Compute the timestep-timestep map
            ID_map_struct map = computeIDMap(last_id->second,new_id->second);
            // Renumber the current timestep's ids
            
        }
        PROFILE_STOP("Identify blobs and maps",1);
    }
private:
    BlobIdentificationWorkItem();
    int Nx, Ny, Nz;
    const RankInfoStruct& rank_info;
    std::shared_ptr<const double> phase;
    const DoubleArray& dist;
    BlobIDstruct last_id, new_id;
};


// Helper class to run the analysis from within a thread
// Note: Averages will be modified after the constructor is called
class AnalysisWorkItem: public ThreadPool::WorkItem
{
public:
    AnalysisWorkItem( AnalysisType type_, int timestep_, TwoPhase& Averages_, BlobIDstruct ids, double beta_ ):
        type(type_), timestep(timestep_), Averages(Averages_), blob_ids(ids), beta(beta_) { }
    virtual void run() {
        Averages.NumberComponents_NWP = blob_ids->first;
        Averages.Label_NWP = blob_ids->second;
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
    }
private:
    AnalysisWorkItem();
    AnalysisType type;
    int timestep;
    TwoPhase& Averages;
    BlobIDstruct blob_ids;
    double beta;
};


// Function to start the analysis
void run_analysis( int timestep, int restart_interval, 
    const RankInfoStruct& rank_info, TwoPhase& Averages, BlobIDstruct& last_ids,
    int Nx, int Ny, int Nz, bool pBC, double beta, double err,
    const double *Phi, double *Pressure, const double *Velocity, 
    const char *ID, const double *f_even, const double *f_odd, const double *Den, 
    const char *LocalRestartFile, ThreadPool& tpool, AnalysisWaitIdStruct& wait )
{
    int N = Nx*Ny*Nz;

    // Determin the analysis we want to perform
    AnalysisType type = AnalyzeNone;
    if ( timestep%1000 == 995 ) {
        // Copy the phase indicator field for the earlier timestep
        type = static_cast<AnalysisType>( type | CopyPhaseIndicator );
    }
    if ( timestep%1000 == 0 ) {
        type = static_cast<AnalysisType>( type | CopyAverages );
        type = static_cast<AnalysisType>( type | IdentifyBlobs );
    }
    if ( timestep%1000 == 5 ) {
        type = static_cast<AnalysisType>( type | CalcDist );
    }
    if (timestep%restart_interval == 0) {
        type = static_cast<AnalysisType>( type | CreateRestart );
    }
    
    // Return if we are not doing anything
    if ( type == AnalyzeNone )
        return;

    PROFILE_START("start_analysis");

    // Copy the appropriate variables to the host (so we can spawn new threads)
    DeviceBarrier();
    PROFILE_START("Copy data to host",1);
    std::shared_ptr<double> phase;
    if ( (type&CopyPhaseIndicator)!=0 || (type&CalcDist)!=0 ||
         (type&CopyAverages)!=0 || (type&IdentifyBlobs)!=0 )
    {
        std::shared_ptr<double>(new double[N],DeleteArray<double>);
        CopyToHost(phase.get(),Phi,N*sizeof(double));
    }
    if ( (type&CopyPhaseIndicator)!=0 ) {
        memcpy(Averages.Phase_tplus.get(),phase.get(),N*sizeof(double));
    }
    if ( (type&CalcDist)!=0 ) {
        memcpy(Averages.Phase_tminus.get(),phase.get(),N*sizeof(double));
    }
    if ( (type&CopyAverages) != 0 ) {
        // Copy the members of Averages to the cpu (phase was copied above)
        ComputePressureD3Q19(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
        memcpy(Averages.Phase.get(),phase.get(),N*sizeof(double));
        DeviceBarrier();
        CopyToHost(Averages.Press.get(),Pressure,N*sizeof(double));
        CopyToHost(Averages.Vel_x.get(),&Velocity[0],N*sizeof(double));
        CopyToHost(Averages.Vel_y.get(),&Velocity[N],N*sizeof(double));
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
        BlobIDstruct new_ids;
        ThreadPool::WorkItem *work = new BlobIdentificationWorkItem(
            Nx,Ny,Nz,rank_info,phase,Averages.SDs,last_ids,new_ids);
        work->add_dependency(wait.blobID);
        last_ids = new_ids;
        wait.blobID = tpool.add_work(work);
    }

    // Spawn threads to do the analysis work
    if ( (type&CalcDist) != 0 ) {
        ThreadPool::WorkItem *work = new AnalysisWorkItem(type,timestep,Averages,last_ids,beta);
        work->add_dependency(wait.blobID);
        work->add_dependency(wait.analysis);
        wait.analysis = tpool.add_work(work);
    }

    // Spawn a thread to write the restart file
    if ( (type&CreateRestart) != 0 ) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
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


