/*
 * Multi-relaxation time LBM Model
 */
#include "models/MRT.h"
#include "models/ElectroModel.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"

ScaLBL_Poisson::ScaLBL_Poisson(int RANK, int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tau(0),
Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),mu(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{

}
ScaLBL_Poisson::~ScaLBL_Poisson(){

}

void ScaLBL_Poisson::ReadParams(string filename){
	// read the input database 
	db = std::make_shared<Database>( filename );
	domain_db = db->getDatabase( "Domain" );
	mrt_db = db->getDatabase( "MRT" );
	electric_db = db->getDatabase( "Electrochemistry" );
	
	tau = 1.0;
	timestepMax = 100000;
	tolerance = 1.0e-8;
	Fx = Fy = 0.0;
	Fz = 1.0e-5;

	// Color Model parameters
	if (mrt_db->keyExists( "timestepMax" )){
		timestepMax = mrt_db->getScalar<int>( "timestepMax" );
	}
	if (mrt_db->keyExists( "tolerance" )){
		tolerance = mrt_db->getScalar<double>( "tolerance" );
	}
	if (mrt_db->keyExists( "tau" )){
		tau = mrt_db->getScalar<double>( "tau" );
	}
	if (mrt_db->keyExists( "F" )){
		Fx = mrt_db->getVector<double>( "F" )[0];
		Fy = mrt_db->getVector<double>( "F" )[1];
		Fz = mrt_db->getVector<double>( "F" )[2];
	}
	if (mrt_db->keyExists( "Restart" )){
		Restart = mrt_db->getScalar<bool>( "Restart" );
	}
	if (mrt_db->keyExists( "din" )){
		din = mrt_db->getScalar<double>( "din" );
	}
	if (mrt_db->keyExists( "dout" )){
		dout = mrt_db->getScalar<double>( "dout" );
	}
	if (mrt_db->keyExists( "flux" )){
		flux = mrt_db->getScalar<double>( "flux" );
	}	
	
	// Read domain parameters
	if (domain_db->keyExists( "BC" )){
		BoundaryCondition = domain_db->getScalar<int>( "BC" );
	}
	

	mu=(tau-0.5)/3.0;
}
void ScaLBL_Poisson::SetDomain(){
	Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
	Mask  = std::shared_ptr<Domain>(new Domain(domain_db,comm));    // mask domain removes immobile phases

	// domain parameters
	Nx = Dm->Nx;
	Ny = Dm->Ny;
	Nz = Dm->Nz;
	Lx = Dm->Lx;
	Ly = Dm->Ly;
	Lz = Dm->Lz;
	
	N = Nx*Ny*Nz;
	Distance.resize(Nx,Ny,Nz);
	Velocity_x.resize(Nx,Ny,Nz);
	Velocity_y.resize(Nx,Ny,Nz);
	Velocity_z.resize(Nx,Ny,Nz);
	
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;               // initialize this way
	//Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm) ); // TwoPhase analysis object
	MPI_Barrier(comm);
	Dm->CommInit();
	MPI_Barrier(comm);
	
	rank = Dm->rank();	
	nprocx = Dm->nprocx();
	nprocy = Dm->nprocy();
	nprocz = Dm->nprocz();
}

void ScaLBL_Poisson::ReadInput(){
    
    sprintf(LocalRankString,"%05d",Dm->rank());
    sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
    sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);

    
    if (domain_db->keyExists( "Filename" )){
    	auto Filename = domain_db->getScalar<std::string>( "Filename" );
    	Mask->Decomp(Filename);
    }
    else if (domain_db->keyExists( "GridFile" )){
    	// Read the local domain data
    	auto input_id = readMicroCT( *domain_db, comm );
    	// Fill the halo (assuming GCW of 1)
    	array<int,3> size0 = { (int) input_id.size(0), (int) input_id.size(1), (int) input_id.size(2) };
    	ArraySize size1 = { (size_t) Mask->Nx, (size_t) Mask->Ny, (size_t) Mask->Nz };
    	ASSERT( (int) size1[0] == size0[0]+2 && (int) size1[1] == size0[1]+2 && (int) size1[2] == size0[2]+2 );
    	fillHalo<signed char> fill( comm, Mask->rank_info, size0, { 1, 1, 1 }, 0, 1 );
    	Array<signed char> id_view;
    	id_view.viewRaw( size1, Mask->id );
    	fill.copy( input_id, id_view );
    	fill.fill( id_view );
    }
    else{
    	Mask->ReadIDs();
    }

    // Generate the signed distance map
	// Initialize the domain and communication
	Array<char> id_solid(Nx,Ny,Nz);
	// Solve for the position of the solid phase
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				int n = k*Nx*Ny+j*Nx+i;
				// Initialize the solid phase
				if (Mask->id[n] > 0)	id_solid(i,j,k) = 1;
				else	     	    id_solid(i,j,k) = 0;
			}
		}
	}
	// Initialize the signed distance function
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				// Initialize distance to +/- 1
				Distance(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
			}
		}
	}
//	MeanFilter(Averages->SDs);
	if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
	CalcDist(Distance,id_solid,*Dm);
    if (rank == 0) cout << "Domain set." << endl;
}

void ScaLBL_Poisson::Create(){
	/*
	 *  This function creates the variables needed to run a LBM 
	 */
	int rank=Mask->rank();
	//.........................................................
	// Initialize communication structures in averaging domain
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i];
	Mask->CommInit();
	Np=Mask->PoreCount();
	//...........................................................................
	if (rank==0)    printf ("Create ScaLBL_Communicator \n");
	// Create a communicator for the device (will use optimized layout)
	// ScaLBL_Communicator ScaLBL_Comm(Mask); // original
	ScaLBL_Comm  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));

	int Npad=(Np/16 + 2)*16;
	if (rank==0)    printf ("Set up memory efficient layout \n");
	Map.resize(Nx,Ny,Nz);       Map.fill(-2);
	auto neighborList= new int[18*Npad];
	Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Mask->id,Np);
	MPI_Barrier(comm);
	//...........................................................................
	//                MAIN  VARIABLES ALLOCATED HERE
	//...........................................................................
	// LBM variables
	if (rank==0)    printf ("Allocating distributions \n");
	//......................device distributions.................................
	int dist_mem_size = Np*sizeof(double);
	int neighborSize=18*(Np*sizeof(int));
	//...........................................................................
	ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
	ScaLBL_AllocateDeviceMemory((void **) &fq, 7*dist_mem_size);  
	ScaLBL_AllocateDeviceMemory((void **) &Psi, sizeof(double)*Np);
	//...........................................................................
	// Update GPU data structures
	if (rank==0)    printf ("Setting up device map and neighbor list \n");
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	MPI_Barrier(comm);
	
}        

void ScaLBL_Poisson::Initialize(){
	/*
	 * This function initializes model
	 */
    if (rank==0)    printf ("Initializing distributions \n");
    ScaLBL_D3Q19_Init(fq, Np);
}

void ScaLBL_Poisson::Run(double *ChargeDensity){

	//.......create and start timer............
	double starttime,stoptime,cputime;
	ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
	starttime = MPI_Wtime();
	if (rank==0) printf("Beginning AA timesteps, timestepMax = %i \n", timestepMax);
	if (rank==0) printf("********************************************************\n");
	timestep=0;
	double error = 1.0;
	double flow_rate_previous = 0.0;
	while (timestep < timestepMax && error > tolerance) {
		//************************************************************************/
		timestep++;
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
		ScaLBL_D3Q7_AAodd_Poisson(NeighborList, fq,  ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		// Set boundary conditions
		/* ... */
		ScaLBL_D3Q7_AAodd_Poisson(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		timestep++;
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
		ScaLBL_D3Q7_AAeven_Poisson(fq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		// Set boundary conditions
		/* ... */
		ScaLBL_D3Q7_AAeven_Poisson(fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		//************************************************************************/
	}
	//************************************************************************/
	stoptime = MPI_Wtime();
	if (rank==0) printf("-------------------------------------------------------------------\n");
	// Compute the walltime per timestep
	cputime = (stoptime - starttime)/timestep;
	// Performance obtained from each node
	double MLUPS = double(Np)/cputime/1000000;

	if (rank==0) printf("********************************************************\n");
	if (rank==0) printf("CPU time = %f \n", cputime);
	if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
	MLUPS *= nprocs;
	if (rank==0) printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
	if (rank==0) printf("********************************************************\n");

}
