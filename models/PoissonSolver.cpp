/*
 * Multi-relaxation time LBM Model
 */
#include "models/PoissonSolver.h"
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
	electric_db = db->getDatabase( "Poisson" );
	
    k2_inv = 4.5;//the inverse of 2nd-rank moment of D3Q7 lattice
    gamma = 0.3;//time step of LB-Poisson equation
	tau = 0.5+k2_inv*gamma;
	timestepMax = 100000;
	tolerance = 1.0e-6;//stopping criterion for obtaining steady-state electricla potential
    h = 1.0;//resolution; unit: um/lu
    epsilon0 = 8.85e-12;//electrical permittivity of vaccum; unit:[C/(V*m)]
    epsilon0_LB = epsilon0*(h*1.0e-6);//unit:[C/(V*lu)]
    epsilonR = 78.4;//default dielectric constant for water
    epsilon_LB = epsilon0_LB*epsilonR;//electrical permittivity 
    analysis_interval = 1000; 

	// LB-Poisson Model parameters
	if (electric_db->keyExists( "timestepMax" )){
		timestepMax = electric_db->getScalar<int>( "timestepMax" );
	}
	if (electric_db->keyExists( "analysis_interval" )){
		analysis_interval = electric_db->getScalar<int>( "analysis_interval" );
	}
	if (electric_db->keyExists( "tolerance" )){
		tolerance = electric_db->getScalar<double>( "tolerance" );
	}
	if (electric_db->keyExists( "gamma" )){
		gamma = electric_db->getScalar<double>( "gamma" );
	}
	if (electric_db->keyExists( "epsilonR" )){
		epsilonR = electric_db->getScalar<double>( "epsilonR" );
	}
	// Read domain parameters
	if (domain_db->keyExists( "voxel_length" )){//default unit: um/lu
		h = domain_db->getScalar<double>( "voxel_length" );
	}
	if (domain_db->keyExists( "BC" )){
		BoundaryCondition = domain_db->getScalar<int>( "BC" );
	}
    //Re-calcualte model parameters if user updates input
    epsilon0_LB = epsilon0*(h*1.0e-6);//unit:[C/(V*lu)]
    epsilon_LB = epsilon0_LB*epsilonR;//electrical permittivity 
	tau = 0.5+k2_inv*gamma;

	if (rank==0) printf("***********************************************************************************\n");
	if (rank==0) printf("LB-Poisson Solver: steady-state MaxTimeStep = %i; steady-state tolerance = %.3g \n", timestepMax,tolerance);
	if (rank==0) printf("                   LB relaxation tau = %.5g \n", tau);
	if (rank==0) printf("***********************************************************************************\n");
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
	Psi_host.resize(Nx,Ny,Nz);

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
	if (rank==0) printf("LB-Poisson Solver: Initialized solid phase & converting to Signed Distance function \n");
	CalcDist(Distance,id_solid,*Dm);
    if (rank == 0) cout << "    Domain set." << endl;
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
	if (rank==0)    printf ("LB-Poisson Solver: Create ScaLBL_Communicator \n");
	// Create a communicator for the device (will use optimized layout)
	// ScaLBL_Communicator ScaLBL_Comm(Mask); // original
	ScaLBL_Comm  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));

	int Npad=(Np/16 + 2)*16;
	if (rank==0)    printf ("LB-Poisson Solver: Set up memory efficient layout \n");
	Map.resize(Nx,Ny,Nz);       Map.fill(-2);
	auto neighborList= new int[18*Npad];
	Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Mask->id,Np);
	MPI_Barrier(comm);
	//...........................................................................
	//                MAIN  VARIABLES ALLOCATED HERE
	//...........................................................................
	// LBM variables
	if (rank==0)    printf ("LB-Poisson Solver: Allocating distributions \n");
	//......................device distributions.................................
	int dist_mem_size = Np*sizeof(double);
	int neighborSize=18*(Np*sizeof(int));
	//...........................................................................
	ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
	ScaLBL_AllocateDeviceMemory((void **) &fq, 7*dist_mem_size);  
	ScaLBL_AllocateDeviceMemory((void **) &Psi, sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &ElectricField, 3*sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &zeta, sizeof(double)*ScaLBL_Comm->n_bb_d3q7);  
	//...........................................................................
	// initialize the zeta function (example is zeta is constant on solid surface)
	double *tmpZeta = new double[ScaLBL_Comm->n_bb_d3q7];
	for int (i=0; i<ScaLBL_Comm->n_bb_d3q7; i++){
		tmpZeta[i] = 1.0/k2_inv; // this has to be read from input file
	}
	ScaLBL_CopyToDevice(zeta, tmpZeta, sizeof(double)*ScaLBL_Comm->n_bb_d3q7);
	delete [] tmpZeta;
	
	// Update GPU data structures
	if (rank==0)    printf ("LB-Poisson Solver: Setting up device map and neighbor list \n");
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	MPI_Barrier(comm);
	
}        

void ScaLBL_Poisson::Initialize(){
	/*
	 * This function initializes model
	 */
    if (rank==0)    printf ("LB-Poisson Solver: initializing D3Q7 distributions\n");
    ScaLBL_D3Q7_Poisson_Init(fq, Np);
}

void ScaLBL_Poisson::Run(double *ChargeDensity){
    
	//.......create and start timer............
	//double starttime,stoptime,cputime;
	//ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
	//starttime = MPI_Wtime();

	timestep=0;
	double error = 1.0;
	double psi_avg_previous = 0.0;
	while (timestep < timestepMax && error > tolerance) {
		//************************************************************************/
		// *************ODD TIMESTEP*************//
        timestep++;
		ScaLBL_Comm->SendD3Q7AA(fq, 0); //READ FROM NORMAL
		ScaLBL_D3Q7_AAodd_Poisson(NeighborList, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, gamma, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->RecvD3Q7AA(fq, 0); //WRITE INTO OPPOSITE
		// Set boundary conditions
		/* ... */
		ScaLBL_D3Q7_AAodd_Poisson(NeighborList, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, gamma, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->SolidDirichletD3Q7(fq, zeta);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

		// *************EVEN TIMESTEP*************//
		timestep++;
		ScaLBL_Comm->SendD3Q7AA(fq, 0); //READ FORM NORMAL
		ScaLBL_D3Q7_AAeven_Poisson(fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, gamma, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->RecvD3Q7AA(fq, 0); //WRITE INTO OPPOSITE
		// Set boundary conditions
		/* ... */
		ScaLBL_D3Q7_AAeven_Poisson(fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, gamma, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->SolidDirichletD3Q7(fq, zeta);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		//************************************************************************/

        // Check convergence of steady-state solution
        if (timestep%analysis_interval==0){
        
			ScaLBL_Comm->RegularLayout(Map,&Psi,Psi_host);
			double count_loc=0;
			double count;
            double psi_avg;
            double psi_loc=0.f;

			for (int k=1; k<Nz-1; k++){
				for (int j=1; j<Ny-1; j++){
					for (int i=1; i<Nx-1; i++){
						if (Distance(i,j,k) > 0){
							psi_loc += Psi_host(i,j,k);
							count_loc+=1.0;
						}
					}
				}
			}
			MPI_Allreduce(&psi_loc,&psi_avg,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
			MPI_Allreduce(&count_loc,&count,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
			
			psi_avg /= count;
            double psi_avg_mag=psi_avg;
		    if (psi_avg==0.0) psi_avg_mag=1.0;
            error = fabs(psi_avg-psi_avg_previous)/fabs(psi_avg_mag);
			psi_avg_previous = psi_avg;
        }
	}

	//************************************************************************/
	//stoptime = MPI_Wtime();
	////if (rank==0) printf("LB-Poission Solver: a steady-state solution is obtained\n");
	////if (rank==0) printf("---------------------------------------------------------------------------\n");
	//// Compute the walltime per timestep
	//cputime = (stoptime - starttime)/timestep;
	//// Performance obtained from each node
	//double MLUPS = double(Np)/cputime/1000000;

	//if (rank==0) printf("******************* LB-Poisson Solver Statistics ********************\n");
	//if (rank==0) printf("CPU time = %f \n", cputime);
	//if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
	//MLUPS *= nprocs;
	//if (rank==0) printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
	//if (rank==0) printf("*********************************************************************\n");

}

