/*
 * Multi-relaxation time LBM Model
 */
#include "models/IonModel.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"

ScaLBL_IonModel::ScaLBL_IonModel(int RANK, int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tau(0),
Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),mu(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{

}
ScaLBL_IonModel::~ScaLBL_IonModel(){

}

void ScaLBL_IonModel::ReadParams(string filename){
	
    //fundamental constant
    kb = 1.38e-23;//Boltzmann constant;unit [J/K]
    electron_charge = 1.6e-19;//electron charge;unit [C]
    
    // read the input database 
	db = std::make_shared<Database>( filename );
	domain_db = db->getDatabase( "Domain" );
	ion_db = db->getDatabase( "Ions" );
	
	// Default model parameters		
    T = 300.0;//temperature; unit [K]
    Vt = kb*T/electron_charge;//thermal voltage; unit [V]
    k2_inv = 4.5;//the inverse of 2nd-rank moment of D3Q7 lattice
    h = 1.0;//resolution; unit: um/lu
	timestepMax = 100000;
	tolerance = 1.0e-8;
	number_ion_species = 1;
    IonDiffusivity.push_back(1.0e-9);//User input unit [m^2/sec]
    //TODO needs to scale the unit of diffusivity!
    IonValence.push_back(1);
    IonConcentration.push_back(1.0e-3);//unit [mol/m^3]
    // TODO rescale ion concentration unit
    deltaT.push_back(1.0);
    tau.push_back(0.5+k2_inv*deltaT[0]*IonDiffusivisty[0]);

	// LB-Ion Model parameters		
	if (ion_db->keyExists( "timestepMax" )){
		timestepMax = ion_db->getScalar<int>( "timestepMax" );
	}
	if (ion_db->keyExists( "analysis_interval" )){
		analysis_interval = ion_db->getScalar<int>( "analysis_interval" );
	}
	if (ion_db->keyExists( "tolerance" )){
		tolerance = ion_db->getScalar<double>( "tolerance" );
	}
	if (ion_db->keyExists( "temperature" )){
		T = ion_db->getScalar<int>( "temperature" );
	}
	if (ion_db->keyExists( "epsilonR" )){
		epsilonR = ion_db->getScalar<double>( "epsilonR" );
	}
	if (ion_db->keyExists( "number_ion_species" )){
		number_ion_species = ion_db->getScalar<int>( "number_ion_species" );
	}

    //read ion related list
	if (ion_db->keyExists( "deltaT" )){
        deltaT.clear();
	    deltaT = ion_db->getVector<double>( "deltaT" );
        if (deltaT.size()!=number_ion_species){
            ERROR("Error: number_ion_species and deltaT must be the same length! \n");
        }
    }
    //NOTE: Ion diffusivity has unit: [m^2/sec]
	if (ion_db->keyExists("IonDiffusivityList")){
        IonDiffusivity.clear();
	    IonDiffusivity = ion_db->getVector<double>( "IonDiffusivityList" );
        if (IonDiffusivity.size()!=number_ion_species){
		    ERROR("Error: number_ion_species and IonDiffusivityList must be the same length! \n");
        }
    }
    //read ion algebric valence list
	if (ion_db->keyExists("IonValenceList")){
        IonValence.clear();
	    IonValence = ion_db->getVector<int>( "IonValenceList" );
        if (IonValence.size()!=number_ion_species){
		    ERROR("Error: number_ion_species and IonValenceList must be the same length! \n");
        }
    }
    //read initial ion concentration list; unit [mol/m^3]
	if (ion_db->keyExists("IonConcentrationList")){
        IonConcentration.clear();
	    IonConcentration = ion_db->getVector<double>( "IonConcentrationList" );
        if (IonConcentration.size()!=number_ion_species){
		    ERROR("Error: number_ion_species and IonConcentrationList must be the same length! \n");
        }
    }

	// Read domain parameters
	if (domain_db->keyExists( "voxel_length" )){//default unit: um/lu
		h = domain_db->getScalar<double>( "voxel_length" );
	}
	if (domain_db->keyExists( "BC" )){
		BoundaryCondition = domain_db->getScalar<int>( "BC" );
	}
    //Re-calcualte model parameters if user updates input
    //TODO ion diffusivity needs rescale unit to LB unit 
    //TODO rescale ion initial concentration unit to LB unit
    if (deltaT.size()>1){
        tau.clear();
        for (int i=0;i<IonDiffusivity.size();i++){
            tau.push_back(0.5+k2_inv*deltaT*IonDiffusivity[i]);
        }
    }
}

void ScaLBL_IonModel::SetDomain(){
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

void ScaLBL_IonModel::ReadInput(){
    
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

void ScaLBL_IonModel::Create(){
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
	ScaLBL_AllocateDeviceMemory((void **) &fq, number_ion_species*7*dist_mem_size);  
	ScaLBL_AllocateDeviceMemory((void **) &Ci, number_ion_species*sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &ChargeDensity, sizeof(double)*Np);
	//...........................................................................
	// Update GPU data structures
	if (rank==0)    printf ("Setting up device map and neighbor list \n");
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	MPI_Barrier(comm);
	
}        

void ScaLBL_IonModel::Initialize(){
	/*
	 * This function initializes model
	 */
    if (rank==0)    printf ("Initializing D3Q7 distributions for ion transport\n");
	for (int ic=0; ic<number_ion_species; ic++){
        ScaLBL_D3Q7_Ion_Init(&fq[ic*Np*7],&Ci[ic*Np],IonConcentration[ic],Np); 
    }
}

void ScaLBL_IonModel::Run(double *Velocity, double *ElectricField){

    //LB-related parameter
    vector<double> rlx(tau.begin(),tau.end());
    for (double item : rlx){
        item = 1.0/item; 
    }
	//.......create and start timer............
	double starttime,stoptime,cputime;
	ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
	starttime = MPI_Wtime();

	if (rank==0) printf("***************************************************\n");
	if (rank==0) printf("LB-Ion Transport: timestepMax = %i\n", timestepMax);
	if (rank==0) printf("***************************************************\n");
	timestep=0;
	while (timestep < timestepMax) {
		//************************************************************************/
		// *************ODD TIMESTEP*************//
		timestep++;
        //Update ion concentration and charge density
		for (int ic=0; ic<number_ion_species; ic++){
			ScaLBL_Comm->SendD3Q7AA(fq, ic); //READ FROM NORMAL
        }
		for (int ic=0; ic<number_ion_species; ic++){
            //TODO should this loop be merged with Send & Recv ?
            //Sum up distribution to get ion concentration
            ScaLBL_D3Q7_AAodd_IonConcentration(NeighborList, &fq[ic*Np*7],&Ci[ic*Np],ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        }
        ScaLBL_D3Q7_IonChargeDensity(Ci, ChargeDensity, IonValence, number_ion_species, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);

		for (int ic=0; ic<number_ion_species; ic++){
			ScaLBL_Comm->RecvD3Q7AA(fq, ic); //WRITE INTO OPPOSITE
        }
		for (int ic=0; ic<number_ion_species; ic++){
            //TODO should this loop be merged with Send & Recv ?
            //Sum up distribution to get ion concentration
            ScaLBL_D3Q7_AAodd_IonConcentration(NeighborList, &fq[ic*Np*7],&Ci[ic*Np], 0, ScaLBL_Comm->LastExterior(), Np);
        }
        ScaLBL_D3Q7_IonChargeDensity(Ci, ChargeDensity, number_ion_species, 0, ScaLBL_Comm->LastExterior(), Np);

        //LB-Ion collison
		for (int ic=0; ic<number_ion_species; ic++){
            ScaLBL_D3Q7_AAodd_Ion(NeighborList, &fq[ic*Np*7],&Ci[ic*Np],Velocity,ElectricField,IonDiffusivity[ic],IonValence[ic],
                                  rlx[ic],deltaT[ic],Vt,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        }

		// Set boundary conditions
		/* ... */

		for (int ic=0; ic<number_ion_species; ic++){
            ScaLBL_D3Q7_AAodd_Ion(NeighborList, &fq[ic*Np*7],&Ci[ic*Np],Velocity,ElectricField,IonDiffusivity[ic],IonValence[ic],
                                  rlx[ic],deltaT[ic],Vt,0, ScaLBL_Comm->LastExterior(), Np);
        }
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

		// *************EVEN TIMESTEP*************//
		timestep++;
        //Update ion concentration and charge density
		for (int ic=0; ic<number_ion_species; ic++){
			ScaLBL_Comm->SendD3Q7AA(fq, ic); //READ FORM NORMAL
        }
		for (int ic=0; ic<number_ion_species; ic++){
            //TODO should this loop be merged with Send & Recv ?
            //Sum up distribution to get ion concentration
            ScaLBL_D3Q7_AAeven_IonConcentration(&fq[ic*Np*7],&Ci[ic*Np],ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        }
        ScaLBL_D3Q7_IonChargeDensity(Ci, ChargeDensity, IonValence, number_ion_species, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);

		for (int ic=0; ic<number_ion_species; ic++){
			ScaLBL_Comm->RecvD3Q7AA(fq, ic); //WRITE INTO OPPOSITE
        }
		for (int ic=0; ic<number_ion_species; ic++){
            //TODO should this loop be merged with Send & Recv ?
            //Sum up distribution to get ion concentration
            ScaLBL_D3Q7_AAeven_IonConcentration(&fq[ic*Np*7],&Ci[ic*Np], 0, ScaLBL_Comm->LastExterior(), Np);
        }
        ScaLBL_D3Q7_IonChargeDensity(Ci, ChargeDensity, number_ion_species, 0, ScaLBL_Comm->LastExterior(), Np);

        //LB-Ion collison
		for (int ic=0; ic<number_ion_species; ic++){
            ScaLBL_D3Q7_AAeven_Ion(&fq[ic*Np*7],&Ci[ic*Np],Velocity,ElectricField,IonDiffusivity[ic],IonValence[ic],
                                  rlx[ic],deltaT[ic],Vt,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        }

		// Set boundary conditions
		/* ... */
		
		for (int ic=0; ic<number_ion_species; ic++){
            ScaLBL_D3Q7_AAeven_Ion(&fq[ic*Np*7],&Ci[ic*Np],Velocity,ElectricField,IonDiffusivity[ic],IonValence[ic],
                                  rlx[ic],deltaT[ic],Vt,0, ScaLBL_Comm->LastExterior(), Np);
        }
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

