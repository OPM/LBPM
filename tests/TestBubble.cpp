#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "analysis/pmmc.h"
#include "common/ScaLBL.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "IO/Mesh.h"
#include "IO/Writer.h"
#include "ProfilerApp.h"

using namespace std;


std::shared_ptr<Database> readInput( const std::string& file )
{
    if ( exists( file )
        return std::make_shared<Database>( file );
    auto db = std::make_shared<Database>();
    db.putData( "Color", std::make_shared<Database>() );
    db.putData( "Domain", std::make_shared<Database>() );
    return db;
}


int main(int argc, char **argv)
{
  // Initialize MPI
  int provided_thread_support = -1;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided_thread_support);
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD,&comm);
  int rank = comm_rank(comm);
  int nprocs = comm_size(comm);
  if ( rank==0 && provided_thread_support<MPI_THREAD_MULTIPLE )
    std::cerr << "Warning: Failed to start MPI with necessary thread support, thread support will be disabled" << std::endl;
  { // Limit scope so variables that contain communicators will free before MPI_Finialize

    // parallel domain size (# of sub-domains)
    int nprocx,nprocy,nprocz;

    MPI_Request req1[18],req2[18];
    MPI_Status stat1[18],stat2[18];

    if (rank == 0){
        printf("********************************************************\n");
        printf("Running Hybrid Implementation of Color LBM    \n");
        printf("********************************************************\n");
    }

    PROFILE_ENABLE(1);
    //PROFILE_ENABLE_TRACE();
    //PROFILE_ENABLE_MEMORY();
    PROFILE_SYNCHRONIZE();
    PROFILE_START("Main");
    Utilities::setErrorHandlers();
    // Variables that specify the computational domain  
    string FILENAME;
    unsigned int nBlocks, nthreads;
    int Nx,Ny,Nz;
    int nspheres;
    double Lx,Ly,Lz;
    // Color Model parameters
    int i,j,k,n;

    // pmmc threshold values
    double fluid_isovalue,solid_isovalue;
    fluid_isovalue = 0.0;
    solid_isovalue = 0.0;
    nBlocks = 32;
    nthreads = 128;
    
    int RESTART_INTERVAL=1000;

    auto db = readInput( "input.in" );

    // Read variables from Color
    auto color_db = db->getDatabase( "Color" );
    auto tau   = color_db->getScalarWithDefault<double>( "tau", 1.0 );  
    auto alpha = color_db->getScalarWithDefault<double>( "alpha", 1e-3 );
    auto beta  = color_db->getScalarWithDefault<double>( "beta", 0.95 );
    auto phi_s = color_db->getScalarWithDefault<double>( "phi_s", 0.0 );
    auto Fx    = color_db->getVectorWithDefault<double>( "F", { 0, 0, 0 } )[0];
    auto Fy    = color_db->getScalarWithDefault<double>( "F", { 0, 0, 0 } )[1];
    auto Fz    = color_db->getScalarWithDefault<double>( "F", { 0, 0, 0 } )[2];
    auto Restart = color_db->getScalarWithDefault<bool>( "Restart", 0 );
    auto pBC   = color_db->getScalarWithDefault<double>( "pBC", 0 );
    auto din   = color_db->getScalarWithDefault<double>( "din", 1.0 );
    auto dout  = color_db->getScalarWithDefault<double>( "dout", 1.0 );
    auto timestepMax = color_db->getScalarWithDefault<int>( "timestepMax", 3 );
    auto interval = color_db->getScalarWithDefault<int>( "interval", 100 );
    auto tol   = color_db->getScalarWithDefault<double>( "tol", 1e-6 );
    auto das   = color_db->getScalarWithDefault<double>( "das", 0.1 );
    auto dbs   = color_db->getScalarWithDefault<double>( "dab", 0.9 );

    // Read variables from Domain
    auto domain_db = db->getDatabase( "Domain" );
    auto n = domain_db->getVectorWithDefault<int>( "n", { 50, 50, 50 } );
    auto L = domain_db->getVectorWithDefault<double>( "L", { 1.0, 1.0, 1.0 } );
    auto nproc = domain_db->getVectorWithDefault<int>( "nproc", { 1, 1, 1 } );
    auto nspheres = domain_db->getScalarWithDefault<int>( "nspheres", 1 );
    int Nx = n[0];
    int Ny = n[1];
    int Nz = n[2];
    int nprocx = nproc[0];
    int nprocy = nproc[1];
    int nprocz = nproc[2];
    double Lx = L[0];
    double Ly = L[1];
    double Lz = L[2];

    // Get the rank info
    const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
    int iproc = rank_info.ix;
    int jproc = rank_info.jy;
    int kproc = rank_info.kz;

    MPI_Barrier(comm);
    // **************************************************************
    // **************************************************************
    double Ps = -(das-dbs)/(das+dbs);
    double rlxA = 1.f/tau;
    double rlxB = 8.f*(2.f-rlxA)/(8.f-rlxA);
    double xIntPos;
    xIntPos = log((1.0+phi_s)/(1.0-phi_s))/(2.0*beta);     
    
    if (nprocs != nprocx*nprocy*nprocz){
        printf("Fatal error in processor number! \n");
        printf("   nprocx =  %i \n",nprocx);
        printf("   nprocy =  %i \n",nprocy);
        printf("   nprocz =  %i \n",nprocz);
        return 1;
    }

    if (rank==0){
        printf("********************************************************\n");
        printf("tau = %f \n", tau);
        printf("alpha = %f \n", alpha);        
        printf("beta = %f \n", beta);
        printf("das = %f \n", das);
        printf("dbs = %f \n", dbs);
        printf("Value of phi at solid surface = %f \n", phi_s);
        printf("Distance to phi = 0.0: %f \n", xIntPos);
        printf("gamma_{wn} = %f \n", 5.796*alpha);
        // printf("cos theta_c = %f \n", 1.05332*Ps);
        printf("Force(x) = %f \n", Fx);
        printf("Force(y) = %f \n", Fy);
        printf("Force(z) = %f \n", Fz);
        printf("Sub-domain size = %i x %i x %i\n",Nz,Nz,Nz);
        printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
        printf("********************************************************\n");
    }

    // Full domain used for averaging (do not use mask for analysis)
    Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,pBC);
    Dm.CommInit();

    // Mask that excludes the solid phase
    Domain Mask(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,pBC);
     
     MPI_Barrier(comm);

    Nx+=2; Ny+=2; Nz += 2;

    int N = Nx*Ny*Nz;
    int dist_mem_size = N*sizeof(double);

    int S = N/nthreads/nBlocks+1;

    // unsigned int nBlocks = N/nthreads + (N%nthreads == 0?0:1);
    // dim3 grid(nBlocks,1,1);

    if (rank==0) printf("Number of blocks = %i \n", nBlocks);
    if (rank==0) printf("Threads per block = %i \n", nthreads);
    if (rank==0) printf("Sweeps per thread = %i \n", S);
    if (rank==0) printf("Number of nodes per side = %i \n", Nx);
    if (rank==0) printf("Total Number of nodes = %i \n", N);
    if (rank==0) printf("********************************************************\n");

    //.......................................................................
    if (rank == 0)    printf("Read input media... \n");
    //.......................................................................
    
    //.......................................................................
    // Filenames used
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
    sprintf(LocalRankString,"%05d",rank);
    sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
    sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
    
//    printf("Local File Name =  %s \n",LocalRankFilename);
    // .......... READ THE INPUT FILE .......................................
//    char value;
    char *id;
    id = new char[N];
    int sum = 0;
    double iVol_global = 1.0/(1.0*Nx*Ny*Nz*nprocs);
    double porosity = 0;

    DoubleArray SDs(Nx,Ny,Nz);
    DoubleArray SDn(Nx,Ny,Nz);
    //.......................................................................
    
    double BubbleRadius = 15.5; // Radius of the capillary tube
    sum=0;
    for (k=0;k<Nz;k++){
        for (j=0;j<Ny;j++){
            for (i=0;i<Nx;i++){
                n = k*Nx*Ny + j*Nz + i;
                // Cylindrical capillary tube aligned with the z direction
                SDs(i,j,k) = 100;
                // Initialize phase positions field
                if (SDs(i,j,k) < 0.0){
                    id[n] = 0;
                }
                else {
                    id[n] = 1;
                    sum++;
                }
            }
        }
    }
    
    // Set up MPI communication structurese
    if (rank==0)    printf ("Setting up communication control structures \n");

	// Initialize communication structures in averaging domain
	for (i=0; i<Mask.Nx*Mask.Ny*Mask.Nz; i++) Mask.id[i] = id[i];
	Mask.CommInit(comm);

    //......................................................................................
    // Create a communicator for the device
    ScaLBL_Communicator ScaLBL_Comm(Mask);

    //...........device phase ID.................................................
    if (rank==0)    printf ("Copying phase ID to device \n");
    char *ID;
    ScaLBL_AllocateDeviceMemory((void **) &ID, N);                        // Allocate device memory
    // Copy to the device
    ScaLBL_CopyToDevice(ID, id, N);
    ScaLBL_DeviceBarrier();
    //...........................................................................

    //...........................................................................
    //                MAIN  VARIABLES ALLOCATED HERE
    //...........................................................................
    // LBM variables
    if (rank==0)    printf ("Allocating distributions \n");
    //......................device distributions.................................
    double *f_even,*f_odd;
    double *A_even,*A_odd,*B_even,*B_odd;
    //...........................................................................
    ScaLBL_AllocateDeviceMemory((void **) &f_even, 10*dist_mem_size);    // Allocate device memory
    ScaLBL_AllocateDeviceMemory((void **) &f_odd, 9*dist_mem_size);    // Allocate device memory
    ScaLBL_AllocateDeviceMemory((void **) &A_even, 4*dist_mem_size);    // Allocate device memory
    ScaLBL_AllocateDeviceMemory((void **) &A_odd, 3*dist_mem_size);    // Allocate device memory
    ScaLBL_AllocateDeviceMemory((void **) &B_even, 4*dist_mem_size);    // Allocate device memory
    ScaLBL_AllocateDeviceMemory((void **) &B_odd, 3*dist_mem_size);    // Allocate device memory
    //...........................................................................
    double *Phi,*Den;
    double *ColorGrad, *Velocity, *Pressure;
    //...........................................................................
    ScaLBL_AllocateDeviceMemory((void **) &Phi, dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &Pressure, dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &ColorGrad, 3*dist_mem_size);
    //...........................................................................
    //...........................................................................
    // Phase indicator (in array form as needed by PMMC algorithm)
    DoubleArray Phase(Nx,Ny,Nz);

    // Extra copies of phi needed to compute time derivatives on CPU
    DoubleArray Phase_tminus(Nx,Ny,Nz);
    DoubleArray Phase_tplus(Nx,Ny,Nz);
    DoubleArray dPdt(Nx,Ny,Nz);

    //copies of data needed to perform checkpointing from cpu
    double *cDen, *cDistEven, *cDistOdd;
    cDen = new double[2*N];
    cDistEven = new double[10*N];
    cDistOdd = new double[9*N];

    // data needed to perform CPU-based averaging
    //double *Vel = new double[3*N];        // fluid velocity

    DoubleArray MeanCurvature(Nx,Ny,Nz);
    DoubleArray GaussCurvature(Nx,Ny,Nz);
    DoubleArray SDs_x(Nx,Ny,Nz);        // Gradient of the signed distance
    DoubleArray SDs_y(Nx,Ny,Nz);
    DoubleArray SDs_z(Nx,Ny,Nz);
    DoubleArray Phase_x(Nx,Ny,Nz);            // Gradient of the phase indicator field
    DoubleArray Phase_y(Nx,Ny,Nz);
    DoubleArray Phase_z(Nx,Ny,Nz);
    DoubleArray Press(Nx,Ny,Nz);
    MeanCurvature.fill(0);
    GaussCurvature.fill(0);
    SDs_x.fill(0);
    SDs_y.fill(0);
    SDs_z.fill(0);
    Phase_x.fill(0);
    Phase_y.fill(0);
    Phase_z.fill(0);
    Press.fill(0);
    
    /*****************************************************************
     VARIABLES FOR THE PMMC ALGORITHM
     ****************************************************************** */
    //...........................................................................
    // Averaging variables
    //...........................................................................
    // local averages (to each MPI process)
    double awn,ans,aws,lwns,nwp_volume;
    double As;
    double vol_w, vol_n;                        // volumes the exclude the interfacial region
    double sat_w;
    double pan,paw;                                // local phase averaged pressure
    // double vx_w,vy_w,vz_w,vx_n,vy_n,vz_n;          // local phase averaged velocity
    // Global averages (all processes)
    double vol_w_global, vol_n_global;            // volumes the exclude the interfacial region
    double awn_global,ans_global,aws_global;
    double lwns_global;
    double efawns,efawns_global;                // averaged contact angle
    double Jwn,Jwn_global;                        // average mean curavture - wn interface
    DoubleArray van(3);
    DoubleArray vaw(3);
    DoubleArray vawn(3);
    DoubleArray Gwn(6);
    DoubleArray Gns(6);
    DoubleArray Gws(6);
    
    double nwp_volume_global;                    // volume for the wetting phase (for saturation)
    // double p_n_global,p_w_global;                // global phase averaged pressure
    // double vx_w_global,vy_w_global,vz_w_global;    // global phase averaged velocity
    // double vx_n_global,vy_n_global,vz_n_global;    // global phase averaged velocity
    double As_global;
    double dEs,dAwn,dAns;                        // Global surface energy (calculated by rank=0)
    double pan_global,paw_global;                                // local phase averaged pressure
    DoubleArray van_global(3);
    DoubleArray vaw_global(3);
    DoubleArray vawn_global(3);
    DoubleArray Gwn_global(6);
    DoubleArray Gns_global(6);
    DoubleArray Gws_global(6);
    //...........................................................................

    // bool add=1;            // Set to false if any corners contain nw-phase ( F > fluid_isovalue)
    int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
    DoubleArray CubeValues(2,2,2);
    // int count_in=0,count_out=0;
    // int nodx,nody,nodz;
    // initialize lists for vertices for surfaces, common line
    DTMutableList<Point> nw_pts(20);
    DTMutableList<Point> ns_pts(20);
    DTMutableList<Point> ws_pts(20);
    DTMutableList<Point> nws_pts(20);
    // initialize triangle lists for surfaces
    IntArray nw_tris(3,20);
    IntArray ns_tris(3,20);
    IntArray ws_tris(3,20);
    // initialize list for line segments
    IntArray nws_seg(2,20);
    DTMutableList<Point> tmp(20);
    DoubleArray Values(20);
    DoubleArray ContactAngle(20);
    DoubleArray Curvature(20);
    DoubleArray InterfaceSpeed(20);
    DoubleArray NormalVector(60);
    
    // IntArray store;
    
    int n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0;
    int n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
    
    // double s,s1,s2,s3;        // Triangle sides (lengths)
    Point A,B,C,P;
    // double area;
    
    // Initialize arrays for local solid surface
    DTMutableList<Point> local_sol_pts(20);
    int n_local_sol_pts = 0;
    IntArray local_sol_tris(3,18);
    int n_local_sol_tris;
    DoubleArray values(20);
    DTMutableList<Point> local_nws_pts(20);
    int n_local_nws_pts;
    
    //int n_nw_tris_beg, n_ns_tris_beg, n_ws_tris_beg;
    int c;
    //int newton_steps = 0;
    //...........................................................................
    int ncubes = (Nx-2)*(Ny-2)*(Nz-2);    // Exclude the "upper" halo
    IntArray cubeList(3,ncubes);
    int nc=0;
    //...........................................................................
    // Set up the cube list (very regular in this case due to lack of blob-ID)
    for (k=0; k<Nz-2; k++){
        for (j=0; j<Ny-2; j++){
            for (i=0; i<Nx-2; i++){
                cubeList(0,nc) = i;
                cubeList(1,nc) = j;
                cubeList(2,nc) = k;
                nc++;
            }
        }
    }
    //...........................................................................
    // Grids used to pack faces on the GPU for MPI
    int faceGrid,edgeGrid,packThreads;
    packThreads=512;
    edgeGrid=1;
    faceGrid=Nx*Ny/packThreads;
    NULL_USE(faceGrid);
    NULL_USE(edgeGrid);
    //...........................................................................    

    
    //...........................................................................
    int timestep = 0;
    if (rank==0) printf("********************************************************\n");
    if (rank==0)    printf("No. of timesteps: %i \n", timestepMax);

    //.......create a stream for the LB calculation.......
    //    cudaStream_t stream;
    //    cudaStreamCreate(&stream);

    //.......create and start timer............
    double starttime,stoptime,cputime;
    MPI_Barrier(comm);
    starttime = MPI_Wtime();
    //.........................................
    //...........................................................................
    //                MAIN  VARIABLES INITIALIZED HERE
    //...........................................................................
    double BubRad[5];
    BubRad[0] = 8.0;
    BubRad[1] = 10.0;
    BubRad[2] = 12.0;
    BubRad[3] = 15.0;
    BubRad[4] = 20.0;
    //...........................................................................
    if (rank==0){
        printf("--------------------------------------------------------------------------------------\n");
        printf("radius ");                                // Timestep, Change in Surface Energy
        printf("sw pw pn awn Jwn ");                    // Scalar averages
        printf("Gwn [xx, yy, zz, xy, xz, yz] ");        // Orientation tensors
        printf("--------------------------------------------------------------------------------------\n");
    }

    for (int bubbleCount=0; bubbleCount<4; bubbleCount++){
        char bubbleCountName[20];
        sprintf(bubbleCountName,"bubbleCount-%i",bubbleCount);
        PROFILE_START(bubbleCountName);
        BubbleRadius = BubRad[bubbleCount]; // Radius of the current bubble

        // Initialize the bubble
        for (k=0;k<Nz;k++){
            for (j=0;j<Ny;j++){
                for (i=0;i<Nx;i++){
                    n = k*Nx*Ny + j*Nz + i;
                    int iglobal= i+(Nx-2)*iproc;
                    int jglobal= j+(Ny-2)*jproc;
                    int kglobal= k+(Nz-2)*kproc;
                    // Cylindrical capillary tube aligned with the z direction
                    SDs(i,j,k) = 100;
                    // Initialize phase positions field
                    // if ((i-0.5*Nx)*(i-0.5*Nx)+(j-0.5*Ny)*(j-0.5*Ny)+(k-0.5*Nz)*(k-0.5*Nz) < BubbleRadius*BubbleRadius){
                    //     id[n] = 2;
                    // }
                    // Initialize phase position field for parallel bubble test
                    if ((iglobal-0.5*(Nx-2)*nprocx)*(iglobal-0.5*(Nx-2)*nprocx)
                            +(jglobal-0.5*(Ny-2)*nprocy)*(jglobal-0.5*(Ny-2)*nprocy)
                            +(kglobal-0.5*(Nz-2)*nprocz)*(kglobal-0.5*(Nz-2)*nprocz) < BubbleRadius*BubbleRadius){
                        id[n] = 2;
                    }
                    else{
                        id[n]=1;
                    }
                }
            }
        }
        // Copy the bubble to the device and initialize
        ScaLBL_CopyToDevice(ID, id, N);

        //...........................................................................
        //...........................................................................
        ScaLBL_D3Q19_Init(ID, f_even, f_odd, Nx, Ny, Nz);
        //......................................................................
        //ScaLBL_Color_InitDistance(ID, Den, Phi, SDs.data(), das, dbs, beta, xIntPos, Nx, Ny, Nz);
        ScaLBL_Color_Init(ID, Den, Phi, das, dbs, Nx, Ny, Nz);
        ScaLBL_D3Q7_Init(ID, A_even, A_odd, &Den[0], Nx, Ny, Nz);
        ScaLBL_D3Q7_Init(ID, B_even, B_odd, &Den[N], Nx, Ny, Nz);
        //......................................................................
        // Once phase has been initialized, map solid to account for 'smeared' interface
        //......................................................................
        for (i=0; i<N; i++)    SDs(i) -= (1.0); //
        //......................................................................
        //.......................................................................
        sprintf(LocalRankString,"%05d",rank);
        sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
        WriteLocalSolidID(LocalRankFilename, id, N);
        sprintf(LocalRankFilename,"%s%s","SDs.",LocalRankString);
        WriteLocalSolidDistance(LocalRankFilename, SDs.data(), N);
        //.......................................................................
        if (Restart == true){
            if (rank==0) printf("Reading restart file! \n");
            // Read in the restart file to CPU buffers
            ReadCheckpoint(LocalRestartFile, cDen, cDistEven, cDistOdd, N);
            // Copy the restart data to the GPU
            ScaLBL_CopyToDevice(f_even,cDistEven,10*N*sizeof(double));
            ScaLBL_CopyToDevice(f_odd,cDistOdd,9*N*sizeof(double));
            ScaLBL_CopyToDevice(Den,cDen,2*N*sizeof(double));
            ScaLBL_DeviceBarrier();
            MPI_Barrier(comm);
        }

        //*************************************************************************
        //         Compute the phase indicator field and reset Copy, Den
        //*************************************************************************
        ScaLBL_ComputePhaseField(ID, Phi, Den, N);
        //*************************************************************************    
        ScaLBL_DeviceBarrier();
        ScaLBL_Comm.SendHalo(Phi);
        ScaLBL_Comm.RecvHalo(Phi);
        ScaLBL_DeviceBarrier();
        MPI_Barrier(comm);
        //*************************************************************************

        if (rank==0 && pBC){
            printf("Setting inlet pressure = %f \n", din);
            printf("Setting outlet pressure = %f \n", dout);
        }
        if (pBC && kproc == 0)    {
            ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);            
            ScaLBL_Color_BC_z(Phi,Den,Velocity,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
        }
            
        if (pBC && kproc == nprocz-1){
            ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
            ScaLBL_Color_BC_Z(Phi,Den,Velocity,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
        }

        //...........................................................................
        // Copy the phase indicator field for the earlier timestep
        ScaLBL_DeviceBarrier();
        ScaLBL_CopyToHost(Phase_tplus.data(),Phi,N*sizeof(double));
        //...........................................................................
        //...........................................................................
        // Copy the data for for the analysis timestep
        //...........................................................................
        // Copy the phase from the GPU -> CPU
        //...........................................................................
        ScaLBL_DeviceBarrier();
        ScaLBL_D3Q19_Pressure(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
        ScaLBL_CopyToHost(Phase.data(),Phi,N*sizeof(double));
        ScaLBL_CopyToHost(Press.data(),Pressure,N*sizeof(double));
        MPI_Barrier(comm);
        //...........................................................................
        
        int timestep = 0;

        //************ MAIN ITERATION LOOP ***************************************/
        PROFILE_START("Time-loop");
        while (timestep < timestepMax){


            //*************************************************************************
            // Fused Color Gradient and Collision 
            //*************************************************************************
            ScaLBL_D3Q19_ColorCollide( ID,f_even,f_odd,Phi,ColorGrad,
                                 Velocity,Nx,Ny,Nz,rlxA,rlxB,alpha,beta,Fx,Fy,Fz);
            //*************************************************************************

            ScaLBL_DeviceBarrier();
            //*************************************************************************
            // Pack and send the D3Q19 distributions
            ScaLBL_Comm.SendD3Q19(f_even, f_odd);
            //*************************************************************************
            
            //*************************************************************************
            //         Carry out the density streaming step for mass transport
            //*************************************************************************
            ScaLBL_D3Q7_ColorCollideMass(ID, A_even, A_odd, B_even, B_odd, Den, Phi, 
                                    ColorGrad, Velocity, beta, N, pBC);
            //*************************************************************************

            ScaLBL_DeviceBarrier();
            MPI_Barrier(comm);
            //*************************************************************************
            //         Swap the distributions for momentum transport
            //*************************************************************************
            ScaLBL_D3Q19_Swap(ID, f_even, f_odd, Nx, Ny, Nz);
            //*************************************************************************

            ScaLBL_DeviceBarrier();
            MPI_Barrier(comm);
            //*************************************************************************
            // Wait for communications to complete and unpack the distributions
            ScaLBL_Comm.RecvD3Q19(f_even, f_odd);
            //*************************************************************************

            ScaLBL_DeviceBarrier();
            //*************************************************************************
            // Pack and send the D3Q7 distributions
            ScaLBL_Comm.BiSendD3Q7(A_even, A_odd, B_even, B_odd);
            //*************************************************************************

            ScaLBL_DeviceBarrier();
            ScaLBL_D3Q7_Swap(ID, A_even, A_odd, Nx, Ny, Nz);
            ScaLBL_D3Q7_Swap(ID, B_even, B_odd, Nx, Ny, Nz);

            ScaLBL_DeviceBarrier();
            MPI_Barrier(comm);

            //*************************************************************************
            // Wait for communication and unpack the D3Q7 distributions
            ScaLBL_Comm.BiRecvD3Q7(A_even, A_odd, B_even, B_odd);
            //*************************************************************************

            ScaLBL_DeviceBarrier();
            //..................................................................................
            ScaLBL_D3Q7_Density(ID, A_even, A_odd, &Den[0], Nx, Ny, Nz);
            ScaLBL_D3Q7_Density(ID, B_even, B_odd, &Den[N], Nx, Ny, Nz);
            
            //*************************************************************************
            //         Compute the phase indicator field 
            //*************************************************************************
            // ScaLBL_ComputePhaseField(ID, Phi, Copy, Den, N);
            ScaLBL_DeviceBarrier();
            MPI_Barrier(comm);
        
            ScaLBL_ComputePhaseField(ID, Phi, Den, N);
            //*************************************************************************
            ScaLBL_Comm.SendHalo(Phi);
            ScaLBL_DeviceBarrier();
            ScaLBL_Comm.RecvHalo(Phi);
            //*************************************************************************

            ScaLBL_DeviceBarrier();

            // Pressure boundary conditions
            if (pBC && kproc == 0)    {
                ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);            
	            ScaLBL_D3Q19_Velocity(ID,f_even,f_odd,Velocity,Nx,Ny,Nz);
                ScaLBL_Color_BC_z(Phi,Den,Velocity,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
            }
                
            if (pBC && kproc == nprocz-1){
                ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
	            ScaLBL_D3Q19_Velocity(ID,f_even,f_odd,Velocity,Nx,Ny,Nz);
                ScaLBL_Color_BC_Z(Phi,Den,Velocity,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
            }
            
            //...................................................................................

            MPI_Barrier(comm);

            // Timestep completed!
            timestep++;

            if (timestep%RESTART_INTERVAL == 0){
                // Copy the data to the CPU
                ScaLBL_CopyToHost(cDistEven,f_even,10*N*sizeof(double));
                ScaLBL_CopyToHost(cDistOdd,f_odd,9*N*sizeof(double));
                ScaLBL_CopyToHost(cDen,Den,2*N*sizeof(double));
                // Read in the restart file to CPU buffers
                WriteCheckpoint(LocalRestartFile, cDen, cDistEven, cDistOdd, N);
            }
        }
        PROFILE_STOP("Time-loop");
        // End the bubble loop
        //...........................................................................
        // Copy the phase indicator field for the later timestep
        ScaLBL_DeviceBarrier();
        ScaLBL_D3Q19_Pressure(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
        ScaLBL_CopyToHost(Phase_tminus.data(),Phi,N*sizeof(double));
        ScaLBL_CopyToHost(Phase_tplus.data(),Phi,N*sizeof(double));
        ScaLBL_CopyToHost(Phase.data(),Phi,N*sizeof(double));
        ScaLBL_CopyToHost(Press.data(),Pressure,N*sizeof(double));

        double temp=0.5/beta;
        for (n=0; n<N; n++){
          double value = Phase.data()[n];
            SDn(n) = temp*log((1.0+value)/(1.0-value));
        }

        //...........................................................................
        // Calculate the time derivative of the phase indicator field
        for (n=0; n<N; n++)    dPdt(n) = 0.1*(Phase_tplus(n) - Phase_tminus(n));
        //...........................................................................

        // Compute the gradients of the phase indicator and signed distance fields
        pmmc_MeshGradient(Phase,Phase_x,Phase_y,Phase_z,Nx,Ny,Nz);
        pmmc_MeshGradient(SDs,SDs_x,SDs_y,SDs_z,Nx,Ny,Nz);
        //...........................................................................
        // Compute the mesh curvature of the phase indicator field
        pmmc_MeshCurvature(SDn, MeanCurvature, GaussCurvature, Nx, Ny, Nz);
        //...........................................................................
        // Fill in the halo region for the mesh gradients and curvature
        //...........................................................................
        // Pressure
        Dm.CommunicateMeshHalo(Press);
        // Mean Curvature
        Dm.CommunicateMeshHalo(MeanCurvature);
        // Gaussian Curvature
        Dm.CommunicateMeshHalo(GaussCurvature);
        // Gradient of the phase indicator field
        Dm.CommunicateMeshHalo(Phase_x);
        Dm.CommunicateMeshHalo(Phase_y);
        Dm.CommunicateMeshHalo(Phase_z);

        //...........................................................................
        // Compute areas using porous medium marching cubes algorithm
        // McClure, Adalsteinsson, et al. (2007)
        //...........................................................................
        awn = aws = ans = lwns = 0.0;
        nwp_volume = 0.0;
        As = 0.0;

        // Compute phase averages
        pan = paw = 0.0;
        vaw(0) = vaw(1) = vaw(2) = 0.0;
        van(0) = van(1) = van(2) = 0.0;
        vawn(0) = vawn(1) = vawn(2) = 0.0;
        Gwn(0) = Gwn(1) = Gwn(2) = 0.0;
        Gwn(3) = Gwn(4) = Gwn(5) = 0.0;
        Gws(0) = Gws(1) = Gws(2) = 0.0;
        Gws(3) = Gws(4) = Gws(5) = 0.0;
        Gns(0) = Gns(1) = Gns(2) = 0.0;
        Gns(3) = Gns(4) = Gns(5) = 0.0;
        vol_w = vol_n =0.0;
        Jwn = efawns = 0.0;
        
        for (c=0;c<ncubes;c++){
            // Get cube from the list
            i = cubeList(0,c);
            j = cubeList(1,c);
            k = cubeList(2,c);

            //...........................................................................
            // Compute volume averages
            for (int p=0;p<8;p++){

                double delphi;
                if ( SDs(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){
                    // 1-D index for this cube corner
                    n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
                    // compute the norm of the gradient of the phase indicator field
                    delphi = sqrt(Phase_x(n)*Phase_x(n)+Phase_y(n)*Phase_y(n)+Phase_z(n)*Phase_z(n));
                    // Compute the non-wetting phase volume contribution
                    if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){
                        nwp_volume += 0.125;
                        // volume the excludes the interfacial region
                        if (delphi < 1e-4){
                            vol_n += 0.125;
                            // pressure
                            pan += 0.125*Press(n);
                        }
                    }
                    else if (delphi < 1e-4){
                        // volume the excludes the interfacial region
                        vol_w += 0.125;
                        // pressure
                        paw += 0.125*Press(n);
                    }
                }
            }

            //...........................................................................
            // Construct the interfaces and common curve
            pmmc_ConstructLocalCube(SDs, SDn, solid_isovalue, fluid_isovalue,
                    nw_pts, nw_tris, values, ns_pts, ns_tris, ws_pts, ws_tris,
                    local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                    n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                    n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
                    i, j, k, Nx, Ny, Nz);

            // Integrate the contact angle
            efawns += pmmc_CubeContactAngle(CubeValues,Values,Phase_x,Phase_y,Phase_z,SDs_x,SDs_y,SDs_z,
                    local_nws_pts,i,j,k,n_local_nws_pts);

            // Integrate the mean curvature
            Jwn    += pmmc_CubeSurfaceInterpValue(CubeValues,MeanCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);

            pmmc_InterfaceSpeed(dPdt, Phase_x, Phase_y, Phase_z, CubeValues, nw_pts, nw_tris,
                    NormalVector, InterfaceSpeed, vawn, i, j, k, n_nw_pts, n_nw_tris);

            //...........................................................................
            // Compute the Interfacial Areas, Common Line length
            /*        awn += pmmc_CubeSurfaceArea(nw_pts,nw_tris,n_nw_tris);
            ans += pmmc_CubeSurfaceArea(ns_pts,ns_tris,n_ns_tris);
            aws += pmmc_CubeSurfaceArea(ws_pts,ws_tris,n_ws_tris);
             */
            As  += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);

            // Compute the surface orientation and the interfacial area
            awn += pmmc_CubeSurfaceOrientation(Gwn,nw_pts,nw_tris,n_nw_tris);
            ans += pmmc_CubeSurfaceOrientation(Gns,ns_pts,ns_tris,n_ns_tris);
            aws += pmmc_CubeSurfaceOrientation(Gws,ws_pts,ws_tris,n_ws_tris);
            lwns +=  pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
            //...........................................................................
        }
        //...........................................................................
        MPI_Barrier(comm);
        MPI_Allreduce(&nwp_volume,&nwp_volume_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&awn,&awn_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&ans,&ans_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&aws,&aws_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&lwns,&lwns_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&As,&As_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&Jwn,&Jwn_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&efawns,&efawns_global,1,MPI_DOUBLE,MPI_SUM,comm);
        // Phase averages
        MPI_Allreduce(&vol_w,&vol_w_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&vol_n,&vol_n_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&paw,&paw_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&pan,&pan_global,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&vaw(0),&vaw_global(0),3,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&van(0),&van_global(0),3,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&vawn(0),&vawn_global(0),3,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&Gwn(0),&Gwn_global(0),6,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&Gns(0),&Gns_global(0),6,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&Gws(0),&Gws_global(0),6,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Barrier(comm);
        //.........................................................................
        // Compute the change in the total surface energy based on the defined interval
        // See McClure, Prins and Miller (2013) 
        //.........................................................................
        dAwn += awn_global;
        dAns += ans_global;
        dEs = 6.01603*alpha*(dAwn + 1.05332*Ps*dAns);
        dAwn = -awn_global;        // Get ready for the next analysis interval
        dAns = -ans_global;

        // Normalize the phase averages 
        // (density of both components = 1.0)
        paw_global = paw_global / vol_w_global;
        vaw_global(0) = vaw_global(0) / vol_w_global;
        vaw_global(1) = vaw_global(1) / vol_w_global;
        vaw_global(2) = vaw_global(2) / vol_w_global;
        pan_global = pan_global / vol_n_global;
        van_global(0) = van_global(0) / vol_n_global;
        van_global(1) = van_global(1) / vol_n_global;
        van_global(2) = van_global(2) / vol_n_global;

        // Normalize surface averages by the interfacial area
        Jwn_global /= awn_global;
        efawns_global /= lwns_global;

        if (awn_global > 0.0)    for (i=0; i<3; i++)        vawn_global(i) /= awn_global;
        if (awn_global > 0.0)    for (i=0; i<6; i++)        Gwn_global(i) /= awn_global;
        if (ans_global > 0.0)    for (i=0; i<6; i++)        Gns_global(i) /= ans_global;
        if (aws_global > 0.0)    for (i=0; i<6; i++)        Gws_global(i) /= aws_global;

        sat_w = 1.0 - nwp_volume_global*iVol_global/porosity;
        // Compute the specific interfacial areas and common line length (per unit volume)
        awn_global = awn_global*iVol_global;
        ans_global = ans_global*iVol_global;
        aws_global = aws_global*iVol_global;
        lwns_global = lwns_global*iVol_global;
        dEs = dEs*iVol_global;

        //.........................................................................
        if (rank==0){
            printf("%.5g ",BubbleRadius);                                    // bubble radius
            printf("%.5g %.5g ",paw_global,pan_global);            // saturation and pressure
            printf("%.5g ",awn_global);                                        // interfacial area
            printf("%.5g ",Jwn_global);                                        // curvature of wn interface
            printf("%.5g %.5g %.5g %.5g %.5g %.5g \n",
                    Gwn_global(0),Gwn_global(1),Gwn_global(2),Gwn_global(3),Gwn_global(4),Gwn_global(5));        // orientation of wn interface
        }


        shared_ptr<IO::TriList> mesh( new IO::TriList() );
        mesh->A.reserve(8*ncubes);
        mesh->B.reserve(8*ncubes);
        mesh->C.reserve(8*ncubes);
        for (c=0;c<ncubes;c++){
            // Get cube from the list
            i = cubeList(0,c);
            j = cubeList(1,c);
            k = cubeList(2,c);
            //...........................................................................
            // Construct the interfaces and common curve
            pmmc_ConstructLocalCube(SDs, Phase, solid_isovalue, fluid_isovalue,
                    nw_pts, nw_tris, values, ns_pts, ns_tris, ws_pts, ws_tris,
                    local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                    n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                    n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
                    i, j, k, Nx, Ny, Nz);

            //.......................................................................................
            // Write the triangle lists to text file
            for (int r=0;r<n_nw_tris;r++){
                A = nw_pts(nw_tris(0,r));
                B = nw_pts(nw_tris(1,r));
                C = nw_pts(nw_tris(2,r));
                // compare the trianlge orientation against the color gradient
                // Orientation of the triangle
                double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                double normal_x = Phase_x(i,j,k);
                double normal_y = Phase_y(i,j,k);
                double normal_z = Phase_z(i,j,k);

                // If the normals don't point in the same direction, flip the orientation of the triangle
                // Right hand rule for triangle orientation is used to determine rendering for most software
                if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                    P = A;
                    A = C;
                    C = P;
                }
                mesh->A.push_back(A);
                mesh->B.push_back(B);
                mesh->C.push_back(C);
            }        
        }
        std::vector<IO::MeshDataStruct> meshData(1);
        meshData[0].meshName = "wn-tris";
        meshData[0].mesh = mesh;
        for (size_t k=0; k<meshData.size(); k++) {
            shared_ptr<IO::Variable> dist( new IO::Variable() );
            dist->name = "distance";
            dist->dim = 1;
            dist->type = IO::VariableType::NodeVariable;
            dist->data.resize(3*mesh->A.size());
            for (size_t i=0; i<mesh->A.size(); i++) {
                const Point& a = mesh->A[i];
                const Point& b = mesh->B[i];
                const Point& c = mesh->C[i];
                dist->data(3*i+0) = sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
                dist->data(3*i+1) = sqrt(b.x*b.x+b.y*b.y+b.z*b.z);
                dist->data(3*i+2) = sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
            }
            meshData[k].vars.push_back(dist);
        }
        IO::writeData( bubbleCount, meshData, comm );

        PROFILE_STOP(bubbleCountName);
    }
    NULL_USE(sat_w);

    //************************************************************************/
    ScaLBL_DeviceBarrier();
    MPI_Barrier(comm);
    stoptime = MPI_Wtime();
    if (rank==0) printf("-------------------------------------------------------------------\n");
    // Compute the walltime per timestep
    cputime = (stoptime - starttime)/timestep;
    // Performance obtained from each node
    double MLUPS = double(Nx*Ny*Nz)/cputime/1000000;
    
    if (rank==0) printf("********************************************************\n");
    if (rank==0) printf("CPU time = %f \n", cputime);
    if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
    MLUPS *= nprocs;
    if (rank==0) printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
    if (rank==0) printf("********************************************************\n");
    
    //************************************************************************/
    // Write out the phase indicator field 
    //************************************************************************/
    sprintf(LocalRankFilename,"%s%s","Phase.",LocalRankString);
    //    printf("Local File Name =  %s \n",LocalRankFilename);
    // ScaLBL_CopyToHost(Phase.data(),Phi,N*sizeof(double));

    FILE *PHASE;
    PHASE = fopen(LocalRankFilename,"wb");
    fwrite(Press.data(),8,N,PHASE);
    // fwrite(MeanCurvature.data(),8,N,PHASE);
    fclose(PHASE);
    
/*  double *DensityValues;
    DensityValues = new double [2*N];
    ScaLBL_CopyToHost(DensityValues,Copy,2*N*sizeof(double));
    FILE *PHASE;
    PHASE = fopen(LocalRankFilename,"wb");
    fwrite(DensityValues,8,2*N,PHASE);
    fclose(PHASE);
*/    //************************************************************************/
    PROFILE_SAVE("TestBubble");
    
    // ****************************************************
    MPI_Barrier(comm);
  } // Limit scope so variables that contain communicators will free before MPI_Finialize
  Utilities::shutdown();
  return 0;
}
