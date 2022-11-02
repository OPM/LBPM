/*
 * Multi-relaxation time LBM Model
 */
#include "models/PoissonSolver.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"


static inline bool fileExists( const std::string &filename )
{
    std::ifstream ifile( filename.c_str() );
    return ifile.good();
}


ScaLBL_Poisson::ScaLBL_Poisson(int RANK, int NP, const Utilities::MPI& COMM):
    rank(RANK), TIMELOG(nullptr), nprocs(NP),timestep(0),timestepMax(0),tau(0),k2_inv(0),tolerance(0),h(0),
    epsilon0(0),epsilon0_LB(0),epsilonR(0),epsilon_LB(0),Vin(0),Vout(0),Nx(0),Ny(0),Nz(0),N(0),Np(0),analysis_interval(0),
    chargeDen_dummy(0),WriteLog(0),nprocx(0),nprocy(0),nprocz(0),
    BoundaryConditionInlet(0),BoundaryConditionOutlet(0),BoundaryConditionSolidList(0),Lx(0),Ly(0),Lz(0),
    Vin0(0),freqIn(0),PhaseShift_In(0),Vout0(0),freqOut(0),PhaseShift_Out(0),
    TestPeriodic(0),TestPeriodicTime(0),TestPeriodicTimeConv(0),TestPeriodicSaveInterval(0),
    comm(COMM)
{
    if ( rank == 0 ) {
	    bool WriteHeader = !fileExists( "PoissonSolver_Convergence.csv" );

	    TIMELOG = fopen("PoissonSolver_Convergence.csv","a+");
	    if (WriteHeader)
		    fprintf(TIMELOG,"Timestep Error\n");
    }
}
ScaLBL_Poisson::~ScaLBL_Poisson()
{
	ScaLBL_FreeDeviceMemory(NeighborList);
	ScaLBL_FreeDeviceMemory(dvcMap);
	ScaLBL_FreeDeviceMemory(Psi);
	ScaLBL_FreeDeviceMemory(Psi_BCLabel);
	ScaLBL_FreeDeviceMemory(ElectricField);
	ScaLBL_FreeDeviceMemory(ResidualError);
	ScaLBL_FreeDeviceMemory(fq);

    if ( TIMELOG )
        fclose( TIMELOG );
}

void ScaLBL_Poisson::ReadParams(string filename){
	// read the input database 
	db = std::make_shared<Database>( filename );
	domain_db = db->getDatabase( "Domain" );
	electric_db = db->getDatabase( "Poisson" );
	
    k2_inv = 3.0;//speed of sound for D3Q19 lattice 
	tau = 0.5+k2_inv;
	timestepMax = 100000;
	tolerance = 1.0e-6;//stopping criterion for obtaining steady-state electricla potential
    h = 1.0;//resolution; unit: um/lu
    epsilon0 = 8.85e-12;//electric permittivity of vaccum; unit:[C/(V*m)]
    epsilon0_LB = epsilon0*(h*1.0e-6);//unit:[C/(V*lu)]
    epsilonR = 78.4;//default dielectric constant of water
    epsilon_LB = epsilon0_LB*epsilonR;//electric permittivity 
    analysis_interval = 1000; 
    chargeDen_dummy = 1.0e-3;//For debugging;unit=[C/m^3]
    WriteLog = false;
    TestPeriodic = false;
    TestPeriodicTime = 1.0;//unit: [sec]
    TestPeriodicTimeConv = 0.01; //unit [sec/lt]
    TestPeriodicSaveInterval = 0.1; //unit [sec]
    Restart = "false";

	// LB-Poisson Model parameters
	if (electric_db->keyExists( "Restart" )){
        Restart = electric_db->getScalar<bool>("Restart");
	}
	if (electric_db->keyExists( "timestepMax" )){
		timestepMax = electric_db->getScalar<int>( "timestepMax" );
	}
	if (electric_db->keyExists( "tau" )){
		tau = electric_db->getScalar<double>( "tau" );
	}
	if (electric_db->keyExists( "analysis_interval" )){
		analysis_interval = electric_db->getScalar<int>( "analysis_interval" );
	}
	if (electric_db->keyExists( "tolerance" )){
		tolerance = electric_db->getScalar<double>( "tolerance" );
	}
    //'tolerance_method' can be {"MSE","MSE_max"}
	tolerance_method = electric_db->getWithDefault<std::string>( "tolerance_method", "MSE" );
	lattice_scheme = electric_db->getWithDefault<std::string>( "lattice_scheme", "D3Q19" );
	if (electric_db->keyExists( "epsilonR" )){
		epsilonR = electric_db->getScalar<double>( "epsilonR" );
	}
	if (electric_db->keyExists( "DummyChargeDen" )){
		chargeDen_dummy = electric_db->getScalar<double>( "DummyChargeDen" );
	}
	if (electric_db->keyExists( "WriteLog" )){
		WriteLog = electric_db->getScalar<bool>( "WriteLog" );
	}
	if (electric_db->keyExists( "TestPeriodic" )){
		TestPeriodic = electric_db->getScalar<bool>( "TestPeriodic" );
	}
	if (electric_db->keyExists( "TestPeriodicTime" )){
		TestPeriodicTime = electric_db->getScalar<double>( "TestPeriodicTime" );
	}
	if (electric_db->keyExists( "TestPeriodicTimeConv" )){
		TestPeriodicTimeConv = electric_db->getScalar<double>( "TestPeriodicTimeConv" );
	}
	if (electric_db->keyExists( "TestPeriodicSaveInterval" )){
		TestPeriodicSaveInterval = electric_db->getScalar<double>( "TestPeriodicSaveInterval" );
	}

    // Read solid boundary condition specific to Poisson equation
    // BC_solid=1: Dirichlet-type surfacen potential
    // BC_solid=2: Neumann-type surfacen charge density
    BoundaryConditionSolidList.push_back(1);
	if (electric_db->keyExists( "BC_SolidList" )){
        BoundaryConditionSolidList.clear();
		BoundaryConditionSolidList = electric_db->getVector<int>( "BC_SolidList" );
	}
    // Read boundary condition for electric potential
    // BC = 0: normal periodic BC
    // BC = 1: fixed electric potential
    // BC = 2: sine/cosine periodic electric potential (need extra input parameters)
    BoundaryConditionInlet = 0;
    BoundaryConditionOutlet = 0;
	if (electric_db->keyExists( "BC_Inlet" )){
		BoundaryConditionInlet = electric_db->getScalar<int>( "BC_Inlet" );
	}
	if (electric_db->keyExists( "BC_Outlet" )){
		BoundaryConditionOutlet = electric_db->getScalar<int>( "BC_Outlet" );
	}

	// Read domain parameters
	if (domain_db->keyExists( "voxel_length" )){//default unit: um/lu
		h = domain_db->getScalar<double>( "voxel_length" );
	}

    //Re-calcualte model parameters if user updates input
    epsilon0_LB = epsilon0*(h*1.0e-6);//unit:[C/(V*lu)]
    epsilon_LB = epsilon0_LB*epsilonR;//electric permittivity 
    
    /* restart string */
    sprintf(LocalRankString, "%05d", rank);
    sprintf(LocalRestartFile, "%s%s", "Psi.", LocalRankString);

	if (rank==0) printf("***********************************************************************************\n");
	if (rank==0) printf("LB-Poisson Solver: steady-state MaxTimeStep = %i; steady-state tolerance = %.3g \n", timestepMax,tolerance);
	if (rank==0) printf("                   LB relaxation tau = %.5g \n", tau);
	if (rank==0) printf("***********************************************************************************\n");
    if (tolerance_method.compare("MSE")==0){
        if (rank==0) printf("LB-Poisson Solver: Use averaged MSE to check solution convergence.\n");
    }
    else if (tolerance_method.compare("MSE_max")==0){
        if (rank==0) printf("LB-Poisson Solver: Use maximum MSE to check solution convergence.\n");
    }
    else{
        if (rank==0) printf("LB-Poisson Solver: tolerance_method=%s cannot be identified!\n",tolerance_method.c_str());
    }
    if (lattice_scheme.compare("D3Q7")==0){
        if (rank==0) printf("LB-Poisson Solver: Use D3Q7 lattice structure.\n");
    }
    else if (lattice_scheme.compare("D3Q19")==0){
        if (rank==0) printf("LB-Poisson Solver: Use D3Q19 lattice structure.\n");
    }
    else{
        if (rank==0) printf("LB-Poisson Solver: lattice_scheme=%s cannot be identified!\n",lattice_scheme.c_str());
    }
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
	Psi_previous.resize(Nx,Ny,Nz);

	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;               // initialize this way
	//Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm) ); // TwoPhase analysis object
	comm.barrier();
    if (BoundaryConditionInlet==0 && BoundaryConditionOutlet==0){
        Dm->BoundaryCondition   = 0;
        Mask->BoundaryCondition = 0;
    }
    else if (BoundaryConditionInlet>0 && BoundaryConditionOutlet>0){
        Dm->BoundaryCondition   = 1;
        Mask->BoundaryCondition = 1;
    }
    else {//i.e. non-periodic and periodic BCs are mixed
        ERROR("Error: check the type of inlet and outlet boundary condition! Mixed periodic and non-periodic BCs are found!\n");
    }
	Dm->CommInit();
	comm.barrier();
	
	rank = Dm->rank();	
	nprocx = Dm->nprocx();
	nprocy = Dm->nprocy();
	nprocz = Dm->nprocz();
}

void ScaLBL_Poisson::ReadInput(){
    
    sprintf(LocalRankString,"%05d",Dm->rank());
    sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
    sprintf(LocalRestartFile,"%s%s","Psi.",LocalRankString);

    
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
    	id_view.viewRaw( size1, Mask->id.data() );
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

void ScaLBL_Poisson::AssignSolidBoundary(double *poisson_solid, int *poisson_solid_BClabel)
{
	signed char VALUE=0;
	double AFFINITY=0.f;
    int BoundaryConditionSolid=0;

	auto LabelList = electric_db->getVector<int>( "SolidLabels" );
	auto AffinityList = electric_db->getVector<double>( "SolidValues" );

	size_t NLABELS = LabelList.size();
	if (NLABELS != AffinityList.size() || NLABELS != BoundaryConditionSolidList.size()){
		ERROR("Error: LB-Poisson Solver: BC_SolidList, SolidLabels and SolidValues all must be of the same length! \n");
	}

	std::vector<double> label_count( NLABELS, 0.0 );
	std::vector<double> label_count_global( NLABELS, 0.0 );
	// Assign the labels

	for (size_t idx=0; idx<NLABELS; idx++) label_count[idx]=0;

	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				int n = k*Nx*Ny+j*Nx+i;
				VALUE=Mask->id[n];
                AFFINITY=0.f;
                BoundaryConditionSolid=0;
				// Assign the affinity from the paired list
				for (unsigned int idx=0; idx < NLABELS; idx++){
					if (VALUE == LabelList[idx]){
						AFFINITY=AffinityList[idx];
                        BoundaryConditionSolid=BoundaryConditionSolidList[idx];
                        if (BoundaryConditionSolid!=1 && BoundaryConditionSolid!=2){
		                    ERROR("Error: LB-Poisson Solver: Note only BC_SolidList of 1 or 2 is supported!\n");
                        }
                        //NOTE need to convert the user input phys unit to LB unit
                        if (BoundaryConditionSolid==2){
                            //for BCS=1, i.e. Dirichlet-type, no need for unit conversion
                            AFFINITY = AFFINITY*(h*h*1.0e-12)/epsilon_LB; 
                        }
						label_count[idx] += 1.0;
						idx = NLABELS;
						//Mask->id[n] = 0; // set mask to zero since this is an immobile component
					}
				}
				poisson_solid[n] = AFFINITY;
                poisson_solid_BClabel[n] = BoundaryConditionSolid;
			}
		}
	}

	for (size_t idx=0; idx<NLABELS; idx++)
		label_count_global[idx]=Dm->Comm.sumReduce(  label_count[idx]);

	if (rank==0){
		printf("LB-Poisson Solver: number of Poisson solid labels: %lu \n",NLABELS);
		for (unsigned int idx=0; idx<NLABELS; idx++){
			VALUE=LabelList[idx];
			AFFINITY=AffinityList[idx];
            BoundaryConditionSolid=BoundaryConditionSolidList[idx];
			double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocs);
            if (BoundaryConditionSolid==1){
			    printf("   label=%d, surface potential=%.3g [V], volume fraction=%.2g\n",VALUE,AFFINITY,volume_fraction); 
            }
            else if (BoundaryConditionSolid==2){
			    printf("   label=%d, surface charge density=%.3g [C/m^2], volume fraction=%.2g\n",VALUE,AFFINITY,volume_fraction); 
            }
            else{
		        ERROR("Error: LB-Poisson Solver: Note only BC_SolidList of 1 or 2 is supported!\n");
            }
		}
	}
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
	ScaLBL_Comm_Regular  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));

	int Npad=(Np/16 + 2)*16;
	if (rank==0)    printf ("LB-Poisson Solver: Set up memory efficient layout \n");
	Map.resize(Nx,Ny,Nz);       Map.fill(-2);
	auto neighborList= new int[18*Npad];
	Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Mask->id.data(),Npad,1);
	comm.barrier();

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
	ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);
	//ScaLBL_AllocateDeviceMemory((void **) &dvcID, sizeof(signed char)*Nx*Ny*Nz);
	ScaLBL_AllocateDeviceMemory((void **) &Psi, sizeof(double)*Nx*Ny*Nz);
	ScaLBL_AllocateDeviceMemory((void **) &Psi_BCLabel, sizeof(int)*Nx*Ny*Nz);
	ScaLBL_AllocateDeviceMemory((void **) &ElectricField, 3*sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &ResidualError, sizeof(double)*Np);
    if (lattice_scheme.compare("D3Q7")==0){
	    ScaLBL_AllocateDeviceMemory((void **) &fq, 7*dist_mem_size);  
    }
    else if (lattice_scheme.compare("D3Q19")==0){
	    ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);  
    }
	//...........................................................................
	
	// Update GPU data structures
	if (rank==0)    printf ("LB-Poisson Solver: Setting up device map and neighbor list \n");
	fflush(stdout);
	int *TmpMap;
	TmpMap=new int[Np];
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				int idx=Map(i,j,k);
				if (!(idx < 0))
					TmpMap[idx] = k*Nx*Ny+j*Nx+i;
			}
		}
	}
	comm.barrier();
	if (rank==0)    printf (" .... LB-Poisson Solver: check  neighbor list \n");
	// check that TmpMap is valid
	for (int idx=0; idx<ScaLBL_Comm->LastExterior(); idx++){
		auto n = TmpMap[idx];
		if (n > Nx*Ny*Nz){
			printf("Bad value! idx=%i \n", n);
			TmpMap[idx] = Nx*Ny*Nz-1;
		}
	}
	for (int idx=ScaLBL_Comm->FirstInterior(); idx<ScaLBL_Comm->LastInterior(); idx++){
		auto n = TmpMap[idx];
		if ( n > Nx*Ny*Nz ){
			printf("Bad value! idx=%i \n",n);
			TmpMap[idx] = Nx*Ny*Nz-1;
		}
	}
	comm.barrier();
	if (rank==0)    printf (" .... LB-Poisson Solver: copy  neighbor list to GPU \n");
	ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
	ScaLBL_Comm->Barrier();
	delete [] TmpMap;
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	ScaLBL_Comm->Barrier();
	comm.barrier();
	delete [] neighborList;
    // copy node ID
	//ScaLBL_CopyToDevice(dvcID, Mask->id, sizeof(signed char)*Nx*Ny*Nz);
	//ScaLBL_Comm->Barrier();
	
    //Initialize solid boundary for electric potential
	// DON'T USE WITH CELLULAR SYSTEM (NO SOLID -- NEED Membrane SOLUTION)
    ScaLBL_Comm->SetupBounceBackList(Map, Mask->id.data(), Np);
	comm.barrier();
}        

void ScaLBL_Poisson::Potential_Init(double *psi_init){

    //set up default boundary input parameters
    Vin0 = Vout0 = 1.0; //unit: [V]
    freqIn = freqOut = 50.0; //unit: [Hz]
    PhaseShift_In = PhaseShift_Out = 0.0; //unit: [radian]
    Vin  = 1.0; //Boundary-z (inlet)  electric potential
    Vout = 1.0; //Boundary-Z (outlet) electric potential

    if (BoundaryConditionInlet==0 && BoundaryConditionOutlet==0){
    
        signed char VALUE=0;
        double AFFINITY=0.f;

        auto LabelList = electric_db->getVector<int>( "InitialValueLabels" );
        auto AffinityList = electric_db->getVector<double>( "InitialValues" );

        size_t NLABELS = LabelList.size();
        if (NLABELS != AffinityList.size()){
            ERROR("Error: LB-Poisson Solver: InitialValueLabels and InitialValues must be of the same length! \n");
        }

        std::vector<double> label_count( NLABELS, 0.0 );
        std::vector<double> label_count_global( NLABELS, 0.0 );

        for (int k=0;k<Nz;k++){
            for (int j=0;j<Ny;j++){
                for (int i=0;i<Nx;i++){
                    int n = k*Nx*Ny+j*Nx+i;
                    VALUE=Mask->id[n];
                    AFFINITY=0.f;
                    // Assign the affinity from the paired list
                    for (unsigned int idx=0; idx < NLABELS; idx++){
                        if (VALUE == LabelList[idx]){
                            AFFINITY=AffinityList[idx];
                            label_count[idx] += 1.0;
                            idx = NLABELS;
                        }
                    }
                    int idx=Map(i,j,k);
                    if (!(idx<0))  psi_init[n] = AFFINITY;
                }
            }
        }

        for (size_t idx=0; idx<NLABELS; idx++)
            label_count_global[idx]=Dm->Comm.sumReduce(  label_count[idx]);

        if (rank==0){
            printf("LB-Poisson Solver: number of Poisson initial-value labels: %lu \n",NLABELS);
            for (unsigned int idx=0; idx<NLABELS; idx++){
                VALUE=LabelList[idx];
                AFFINITY=AffinityList[idx];
                double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocs);
                printf("   label=%d, initial potential=%.3g [V], volume fraction=%.2g\n",VALUE,AFFINITY,volume_fraction); 
            }
        }
    }
    else if (BoundaryConditionInlet>0 && BoundaryConditionOutlet>0){
        //read input parameters for inlet
        switch (BoundaryConditionInlet){
            case 1:
                if (electric_db->keyExists( "Vin" )){
                    Vin = electric_db->getScalar<double>( "Vin" );
                }
                if (rank==0) printf("LB-Poisson Solver: inlet boundary; fixed electric potential Vin = %.3g [V]\n",Vin);
                break;
            case 2:
                if (electric_db->keyExists( "Vin0" )){//voltage amplitude; unit: Volt
                    Vin0 = electric_db->getScalar<double>( "Vin0" );
                }
                if (electric_db->keyExists( "freqIn" )){//unit: Hz
                    freqIn = electric_db->getScalar<double>( "freqIn" );
                }
                if (electric_db->keyExists( "PhaseShift_In" )){//phase shift, unit: radian
                    PhaseShift_In = electric_db->getScalar<double>( "PhaseShift_In" );
                }
                if (rank==0){
                    printf("LB-Poisson Solver: inlet boundary; periodic electric potential Vin = %.3g*Cos[2*pi*%.3g*t+%.3g] [V] \n",Vin0,freqIn,PhaseShift_In);
                    printf("                                   V0 = %.3g [V], frequency = %.3g [Hz], phase shift = %.3g [radian] \n",Vin0,freqIn,PhaseShift_In);
                } 
                break;
        }
        //read input parameters for outlet
        switch (BoundaryConditionOutlet){
            case 1:
                if (electric_db->keyExists( "Vout" )){
                    Vout = electric_db->getScalar<double>( "Vout" );
                }
                if (rank==0) printf("LB-Poisson Solver: outlet boundary; fixed electric potential Vout = %.3g [V] \n",Vout);
                break;
            case 2:
                if (electric_db->keyExists( "Vout0" )){//voltage amplitude; unit: Volt
                    Vout0 = electric_db->getScalar<double>( "Vout0" );
                }
                if (electric_db->keyExists( "freqOut" )){//unit: Hz
                    freqOut = electric_db->getScalar<double>( "freqOut" );
                }
                if (electric_db->keyExists( "PhaseShift_Out" )){//timestep shift, unit: lt
                    PhaseShift_Out = electric_db->getScalar<double>( "PhaseShift_Out" );
                }
                if (rank==0){
                    printf("LB-Poisson Solver: outlet boundary; periodic electric potential Vout = %.3g*Cos[2*pi*%.3g*t+%.3g] [V]\n",Vout0,freqOut,PhaseShift_Out);
                    printf("                                    V0 = %.3g [V], frequency = %.3g [Hz], timestep shift = %.3g [radian] \n",Vout0,freqOut,PhaseShift_Out);
                } 
                break;
        }
        //calcualte inlet/outlet voltage for the case of BCInlet/Outlet=2
        if (BoundaryConditionInlet==2)  Vin  = getBoundaryVoltagefromPeriodicBC(Vin0,freqIn,PhaseShift_In,0);
        if (BoundaryConditionOutlet==2) Vout = getBoundaryVoltagefromPeriodicBC(Vout0,freqOut,PhaseShift_Out,0);

        //initialize a linear electrical potential between inlet and outlet
        double slope = (Vout-Vin)/(Nz-2);
        double psi_linearized;
        for (int k=0;k<Nz;k++){
            if (k==0 || k==1){
                psi_linearized = Vin;
            }
            else if (k==Nz-1 || k==Nz-2){
                psi_linearized = Vout;
            }
            else{
                psi_linearized = slope*(k-1)+Vin;
            }
            for (int j=0;j<Ny;j++){
                for (int i=0;i<Nx;i++){
                    int n = k*Nx*Ny+j*Nx+i;
                    if (Mask->id[n]>0){
                        psi_init[n] = psi_linearized;
                    }
                }
            }
        }
    }
    else{//mixed periodic and non-periodic BCs are not supported!
        ERROR("Error: check the type of inlet and outlet boundary condition! Mixed periodic and non-periodic BCs are found!\n");
    }
    /** RESTART **/
     if (Restart == true) {
         if (rank == 0) {
             printf("   POISSON MODEL: Reading restart file! \n");
         }
         ifstream File(LocalRestartFile, ios::binary);
         double value;
         // Read the distributions
         for (int n = 0; n < Nx*Ny*Nz; n++) {
        	 File.read((char *)&value, sizeof(value));
        	 psi_init[n] = value;
         }
         File.close();
     }
     /** END RESTART **/
}

double ScaLBL_Poisson::getBoundaryVoltagefromPeriodicBC(double V0, double freq, double phase_shift, int time_step){
    return V0*cos(2.0*M_PI*freq*time_conv*time_step+phase_shift);
}

void ScaLBL_Poisson::Initialize(double time_conv_from_Study){
	/*
	 * This function initializes model
     * "time_conv_from_Study" is the phys to LB time conversion factor, unit=[sec/lt]
     * which is used for periodic voltage input for inlet and outlet boundaries
	 */
    if (lattice_scheme.compare("D3Q7")==0){
        if (rank==0)    printf ("LB-Poisson Solver: initializing D3Q7 distributions\n");
    }
    else if (lattice_scheme.compare("D3Q19")==0){
        if (rank==0)    printf ("LB-Poisson Solver: initializing D3Q19 distributions\n");
    }
    //NOTE the initialization involves two steps:
    //1. assign solid boundary value (surface potential or surface change density)
    //2. Initialize electric potential for pore nodes
    double *psi_host;
    int *psi_BCLabel_host;
    psi_host = new double [Nx*Ny*Nz];
    psi_BCLabel_host = new int [Nx*Ny*Nz];
    time_conv = time_conv_from_Study;
    AssignSolidBoundary(psi_host,psi_BCLabel_host);//step1
    Potential_Init(psi_host);//step2
	ScaLBL_CopyToDevice(Psi, psi_host, Nx*Ny*Nz*sizeof(double));
	ScaLBL_CopyToDevice(Psi_BCLabel, psi_BCLabel_host, Nx*Ny*Nz*sizeof(int));
	ScaLBL_Comm->Barrier();
    if (lattice_scheme.compare("D3Q7")==0){
        ScaLBL_D3Q7_Poisson_Init(dvcMap, fq, Psi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_D3Q7_Poisson_Init(dvcMap, fq, Psi, 0, ScaLBL_Comm->LastExterior(), Np);
    }
    else if (lattice_scheme.compare("D3Q19")==0){
        /* switch to d3Q19 model */
        ScaLBL_D3Q19_Poisson_Init(dvcMap, fq, Psi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_D3Q19_Poisson_Init(dvcMap, fq, Psi, 0, ScaLBL_Comm->LastExterior(), Np);
    }
    delete [] psi_host;
    delete [] psi_BCLabel_host;
    
    //extra treatment for halo layer
    //if (BoundaryCondition==1){
	//	if (Dm->kproc()==0){
	//		ScaLBL_SetSlice_z(Psi,Vin,Nx,Ny,Nz,0);
	//	}
	//	if (Dm->kproc() == nprocz-1){
	//		ScaLBL_SetSlice_z(Psi,Vout,Nx,Ny,Nz,Nz-1);
	//	}
    //}
}

//void ScaLBL_Poisson::Run(double *ChargeDensity, bool UseSlippingVelBC, int timestep_from_Study){
//    
//	//.......create and start timer............
//	//double starttime,stoptime,cputime;
//	//comm.barrier();
//    //auto t1 = std::chrono::system_clock::now();
//	double *host_Error;
//	host_Error = new double [Np];
//
//	timestep=0;
//	double error = 1.0;
//	while (timestep < timestepMax && error > tolerance) {
//		//************************************************************************/
//		// *************ODD TIMESTEP*************//
//        timestep++;
//        //SolveElectricPotentialAAodd(timestep_from_Study,ChargeDensity, UseSlippingVelBC);//update electric potential
//        SolvePoissonAAodd(ChargeDensity, UseSlippingVelBC);//perform collision
//		ScaLBL_Comm->Barrier(); comm.barrier();
//
//		// *************EVEN TIMESTEP*************//
//		timestep++;
//		//SolveElectricPotentialAAeven(timestep_from_Study,ChargeDensity, UseSlippingVelBC);//update electric potential
//        SolvePoissonAAeven(ChargeDensity, UseSlippingVelBC);//perform collision
//		ScaLBL_Comm->Barrier(); comm.barrier();
//		//************************************************************************/
//
//		
//        // Check convergence of steady-state solution
//        if (timestep==2){
//            //save electric potential for convergence check
//        }
//        if (timestep%analysis_interval==0){
//	  /* get the elecric potential */
//            ScaLBL_CopyToHost(Psi_host.data(),Psi,sizeof(double)*Nx*Ny*Nz);
//        	if (rank==0) printf("   ... getting Poisson solver error \n");
//        	double err = 0.0;
//        	double max_error = 0.0;
//        	ScaLBL_CopyToHost(host_Error,ResidualError,sizeof(double)*Np);
//        	for (int idx=0; idx<Np; idx++){
//        		err = host_Error[idx]*host_Error[idx];
//        		if (err > max_error ){
//        			max_error = err;
//        		}
//        	}
//        	error=Dm->Comm.maxReduce(max_error);
//
//        	/* compute the eletric field */
//        	//ScaLBL_D3Q19_Poisson_getElectricField(fq, ElectricField, tau, Np);
//
//        }
//	}
//    if(WriteLog==true){
//        getConvergenceLog(timestep,error);
//    }
//
//	//************************************************************************/
//	////if (rank==0) printf("LB-Poission Solver: a steady-state solution is obtained\n");
//	////if (rank==0) printf("---------------------------------------------------------------------------\n");
//	//// Compute the walltime per timestep
//	//auto t2 = std::chrono::system_clock::now();
//	//double cputime = std::chrono::duration<double>( t2 - t1 ).count() / timestep;
//	//// Performance obtained from each node
//	//double MLUPS = double(Np)/cputime/1000000;
//
//	//if (rank==0) printf("******************* LB-Poisson Solver Statistics ********************\n");
//	//if (rank==0) printf("CPU time = %f \n", cputime);
//	//if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
//	//MLUPS *= nprocs;
//	//if (rank==0) printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
//	//if (rank==0) printf("*********************************************************************\n");
//
//}

void ScaLBL_Poisson::Run(double *ChargeDensity, bool UseSlippingVelBC, int timestep_from_Study){
    
    double error = 1.0;
    if (lattice_scheme.compare("D3Q7")==0){
        timestep=0;
        while (timestep < timestepMax && error > tolerance) {
            //************************************************************************/
            // *************ODD TIMESTEP*************//
            timestep++;
            SolveElectricPotentialAAodd(timestep_from_Study);//update electric potential
            SolvePoissonAAodd(ChargeDensity, UseSlippingVelBC,timestep);//perform collision
            ScaLBL_Comm->Barrier(); comm.barrier();

            // *************EVEN TIMESTEP*************//
            timestep++;
            SolveElectricPotentialAAeven(timestep_from_Study);//update electric potential
            SolvePoissonAAeven(ChargeDensity, UseSlippingVelBC,timestep);//perform collision
            ScaLBL_Comm->Barrier(); comm.barrier();
            //************************************************************************/

            // Check convergence of steady-state solution
            if (timestep==2){
                //save electric potential for convergence check
                ScaLBL_CopyToHost(Psi_previous.data(),Psi,sizeof(double)*Nx*Ny*Nz);
            }
            if (timestep%analysis_interval==0){
                if (tolerance_method.compare("MSE")==0){
                    double count_loc=0;
                    double count;
                    double MSE_loc=0.0;
                    ScaLBL_CopyToHost(Psi_host.data(),Psi,sizeof(double)*Nx*Ny*Nz);
                    for (int k=1; k<Nz-1; k++){
                        for (int j=1; j<Ny-1; j++){
                            for (int i=1; i<Nx-1; i++){
                                if (Distance(i,j,k) > 0){
                                    MSE_loc += (Psi_host(i,j,k) - Psi_previous(i,j,k))*(Psi_host(i,j,k) - Psi_previous(i,j,k));
                                    count_loc+=1.0;
                                }
                            }
                        }
                    }
                    error=Dm->Comm.sumReduce(MSE_loc);
                    count=Dm->Comm.sumReduce(count_loc);
                    error /= count;
                }
                else if (tolerance_method.compare("MSE_max")==0){
                    vector<double>MSE_loc;
                    double MSE_loc_max;
                    ScaLBL_CopyToHost(Psi_host.data(),Psi,sizeof(double)*Nx*Ny*Nz);
                    for (int k=1; k<Nz-1; k++){
                        for (int j=1; j<Ny-1; j++){
                            for (int i=1; i<Nx-1; i++){
                                if (Distance(i,j,k) > 0){
                                    MSE_loc.push_back((Psi_host(i,j,k) - Psi_previous(i,j,k))*(Psi_host(i,j,k) - Psi_previous(i,j,k)));
                                }
                            }
                        }
                    }
                    vector<double>::iterator it_max = max_element(MSE_loc.begin(),MSE_loc.end());
                    unsigned int idx_max=distance(MSE_loc.begin(),it_max);
                    MSE_loc_max=MSE_loc[idx_max];
                    error=Dm->Comm.maxReduce(MSE_loc_max);
                }
                else{
                    ERROR("Error: user-specified tolerance_method cannot be identified; check you input database! \n");
                }
                ScaLBL_CopyToHost(Psi_previous.data(),Psi,sizeof(double)*Nx*Ny*Nz);
            }
        }
    }
    else if (lattice_scheme.compare("D3Q19")==0){

        double *host_Error;
        host_Error = new double [Np];

        timestep=0;
        auto t1 = std::chrono::system_clock::now();
        while (timestep < timestepMax && error > tolerance) {
            //************************************************************************/
            // *************ODD TIMESTEP*************//
            timestep++;
            //SolveElectricPotentialAAodd(timestep_from_Study);//,ChargeDensity, UseSlippingVelBC);//update electric potential
            SolvePoissonAAodd(ChargeDensity, UseSlippingVelBC,timestep);//perform collision
            ScaLBL_Comm->Barrier(); comm.barrier();

            // *************EVEN TIMESTEP*************//
            timestep++;
            //SolveElectricPotentialAAeven(timestep_from_Study);//,ChargeDensity, UseSlippingVelBC);//update electric potential
            SolvePoissonAAeven(ChargeDensity, UseSlippingVelBC,timestep);//perform collision
            ScaLBL_Comm->Barrier(); comm.barrier();
            //************************************************************************/

            // Check convergence of steady-state solution
            if (timestep==2){
                //save electric potential for convergence check
            }
            if (timestep%analysis_interval==0){
          /* get the elecric potential */
                ScaLBL_CopyToHost(Psi_host.data(),Psi,sizeof(double)*Nx*Ny*Nz);
                if (rank==0) printf("   ... getting Poisson solver error \n");
                double err = 0.0;
                double max_error = 0.0;
                ScaLBL_CopyToHost(host_Error,ResidualError,sizeof(double)*Np);
                for (int idx=0; idx<Np; idx++){
                    err = host_Error[idx]*host_Error[idx];
                    if (err > max_error ){
                        max_error = err;
                    }
                }
                error=Dm->Comm.maxReduce(max_error);

                /* compute the eletric field */
                //ScaLBL_D3Q19_Poisson_getElectricField(fq, ElectricField, tau, Np);

            }
        }
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
        delete [] host_Error;

    }
    //************************************************************************/

    if(WriteLog==true){
        getConvergenceLog(timestep,error);
    }
}

void ScaLBL_Poisson::Run(double *ChargeDensity, DoubleArray MembraneDistance, bool UseSlippingVelBC, int timestep_from_Study){

	double error = 1.0;
	double threshold = 10000000.0;
	bool SET_THRESHOLD = false;
	if (electric_db->keyExists( "rescale_at_distance" )){
		SET_THRESHOLD = true;
		threshold = electric_db->getScalar<double>( "rescale_at_distance" );
	}
	if (BoundaryConditionInlet > 0)    SET_THRESHOLD = false;
	if (BoundaryConditionOutlet > 0)   SET_THRESHOLD = false;
	
	double *host_Error;
	host_Error = new double [Np];

	timestep=0;
	auto t1 = std::chrono::system_clock::now();
	while (timestep < timestepMax && error > tolerance) {
		//************************************************************************/
		// *************ODD TIMESTEP*************//
		timestep++;
		//SolveElectricPotentialAAodd(timestep_from_Study,ChargeDensity, UseSlippingVelBC);//update electric potential
		SolvePoissonAAodd(ChargeDensity, UseSlippingVelBC,timestep);//perform collision
		ScaLBL_Comm->Barrier(); comm.barrier();

		// *************EVEN TIMESTEP*************//
		timestep++;
		//SolveElectricPotentialAAeven(timestep_from_Study,ChargeDensity, UseSlippingVelBC);//update electric potential
		SolvePoissonAAeven(ChargeDensity, UseSlippingVelBC,timestep);//perform collision
		ScaLBL_Comm->Barrier(); comm.barrier();
		//************************************************************************/

		// Check convergence of steady-state solution
		if (timestep==2){
			//save electric potential for convergence check
		}
		if (timestep%analysis_interval==0){
			/* get the elecric potential */
			ScaLBL_CopyToHost(Psi_host.data(),Psi,sizeof(double)*Nx*Ny*Nz);
			if (rank==0) printf("   ... getting Poisson solver error \n");
			double err = 0.0;
			double max_error = 0.0;
			ScaLBL_CopyToHost(host_Error,ResidualError,sizeof(double)*Np);
			for (int idx=0; idx<Np; idx++){
				err = host_Error[idx]*host_Error[idx];
				if (err > max_error ){
					max_error = err;
				}
			}
			error=Dm->Comm.maxReduce(max_error);

			if (error > tolerance && SET_THRESHOLD){
				/* don't use this with an external BC */
				// cpompute the far-field electric potential
				double inside_local = 0.0;
				double outside_local = 0.0;
				double inside_count_local = 0.0;
				double outside_count_local = 0.0;
				/* global values */
				double inside_global = 0.0;
				double outside_global = 0.0;
				double inside_count_global = 0.0;
				double outside_count_global = 0.0;
				for (int k=1; k<Nz; k++){
					for (int j=1; j<Ny; j++){
						for (int i=1; i<Nx; i++){
							int n = k*Nx*Ny + j*Nx + i;
							double distance = MembraneDistance(i,j,k);
							if  (distance > threshold  && distance < (threshold + 1.0)){
								outside_count_local += 1.0;
								outside_local += Psi_host(n);
							}
							else if  (distance < (-1.0)*threshold && distance > (-1.0)*(threshold + 1.0)){
								inside_count_local += 1.0;
								inside_local += Psi_host(n);
							}
						}
					}
				}
				inside_count_global = Dm->Comm.sumReduce(inside_count_local);
				outside_count_global = Dm->Comm.sumReduce(outside_count_local);
				outside_global = Dm->Comm.sumReduce(outside_local);
				inside_global = Dm->Comm.sumReduce(inside_local);
				outside_global /= outside_count_global;
				inside_global  /= inside_count_global;
				
				if (rank==0) printf("   Rescaling far-field electric potential to value (outside): %f \n",outside_global);
				if (rank==0) printf("   Rescaling far-field electric potential to value (inside): %f \n",inside_global);

				// rescale the far-field electric potential
				for (int k=1; k<Nz; k++){
					for (int j=1; j<Ny; j++){
						for (int i=1; i<Nx; i++){
							int n = k*Nx*Ny + j*Nx + i;
							double distance = MembraneDistance(i,j,k);
							if  ( distance > (threshold + 1.0)){
								Psi_host(n) = outside_global;
							}
							else if  ( distance < (-1.0)*(threshold + 1.0)){
								Psi_host(n) = inside_global;
							}
						}
					}
				}				
				ScaLBL_CopyToDevice(Psi,Psi_host.data(),sizeof(double)*Nx*Ny*Nz);
			}
			/* compute the eletric field */
			//ScaLBL_D3Q19_Poisson_getElectricField(fq, ElectricField, tau, Np);
		}
	}
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
	delete [] host_Error;

	//************************************************************************/
	if(WriteLog==true){
		getConvergenceLog(timestep,error);
	}
}

void ScaLBL_Poisson::getConvergenceLog(int timestep,double error){
    if ( rank == 0 ) {
        fprintf(TIMELOG,"%i %.5g\n",timestep,error); 
        fflush(TIMELOG);
    }
}

void ScaLBL_Poisson::SolveElectricPotentialAAodd(int timestep_from_Study){

    if (lattice_scheme.compare("D3Q7")==0){
        ScaLBL_Comm->SendD3Q7AA(fq, 0); //READ FROM NORMAL
        ScaLBL_D3Q7_AAodd_Poisson_ElectricPotential(NeighborList, dvcMap, fq, Psi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q7AA(fq, 0); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        // Set boundary conditions
        if (BoundaryConditionInlet > 0){
            switch (BoundaryConditionInlet){
                case 1:
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
                case 2:
                    Vin  = getBoundaryVoltagefromPeriodicBC(Vin0,freqIn,PhaseShift_In,timestep_from_Study);
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
            }
        }
        if (BoundaryConditionOutlet > 0){
            switch (BoundaryConditionOutlet){
                case 1:
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_Z(NeighborList, fq, Vout, timestep);
                    break;
                case 2:
                    Vout = getBoundaryVoltagefromPeriodicBC(Vout0,freqOut,PhaseShift_Out,timestep_from_Study);
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_Z(NeighborList, fq, Vout, timestep);
                    break;
            }
        }
        //-------------------------//
        ScaLBL_D3Q7_AAodd_Poisson_ElectricPotential(NeighborList, dvcMap, fq, Psi, 0, ScaLBL_Comm->LastExterior(), Np);
    }
    else if (lattice_scheme.compare("D3Q19")==0){
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        //ScaLBL_D3Q19_AAodd_Poisson_ElectricPotential(NeighborList, dvcMap, fq, ChargeDensity, Psi, epsilon_LB, UseSlippingVelBC, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        
        // Set boundary conditions
        if (BoundaryConditionInlet > 0){
            switch (BoundaryConditionInlet){
                case 1:
                    ScaLBL_Comm->D3Q19_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
                case 2:
                    Vin  = getBoundaryVoltagefromPeriodicBC(Vin0,freqIn,PhaseShift_In,timestep_from_Study);
                    ScaLBL_Comm->D3Q19_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
            }
        }
        if (BoundaryConditionOutlet > 0){
            switch (BoundaryConditionOutlet){
                case 1:
                    ScaLBL_Comm->D3Q19_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
                case 2:
                    Vout = getBoundaryVoltagefromPeriodicBC(Vout0,freqOut,PhaseShift_Out,timestep_from_Study);
                    ScaLBL_Comm->D3Q19_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
            }
        }
        //-------------------------//
        
        //ScaLBL_D3Q19_AAodd_Poisson_ElectricPotential(NeighborList, dvcMap, fq, ChargeDensity, Psi, epsilon_LB, UseSlippingVelBC, 0, ScaLBL_Comm->LastExterior(), Np);
    }
}

void ScaLBL_Poisson::SolveElectricPotentialAAeven(int timestep_from_Study){

    if (lattice_scheme.compare("D3Q7")==0){
        ScaLBL_Comm->SendD3Q7AA(fq, 0); //READ FORM NORMAL
        ScaLBL_D3Q7_AAeven_Poisson_ElectricPotential(dvcMap, fq, Psi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q7AA(fq, 0); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        // Set boundary conditions
        if (BoundaryConditionInlet > 0){
            switch (BoundaryConditionInlet){
                case 1:
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
                case 2:
                    Vin  = getBoundaryVoltagefromPeriodicBC(Vin0,freqIn,PhaseShift_In,timestep_from_Study);
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
            }
        }
        if (BoundaryConditionOutlet > 0){
            switch (BoundaryConditionOutlet){
                case 1:
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_Z(NeighborList, fq, Vout, timestep);
                    break;
                case 2:
                    Vout = getBoundaryVoltagefromPeriodicBC(Vout0,freqOut,PhaseShift_Out,timestep_from_Study);
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_Z(NeighborList, fq, Vout, timestep);
                    break;
            }
        }
        //-------------------------//
        ScaLBL_D3Q7_AAeven_Poisson_ElectricPotential(dvcMap, fq, Psi, 0, ScaLBL_Comm->LastExterior(), Np);
    }
    else if (lattice_scheme.compare("D3Q19")==0){
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        //ScaLBL_D3Q19_AAeven_Poisson_ElectricPotential(dvcMap, fq, ChargeDensity, Psi, epsilon_LB, UseSlippingVelBC,
        //		ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->Barrier();
        

        // Set boundary conditions
        if (BoundaryConditionInlet > 0){
            switch (BoundaryConditionInlet){
                case 1:
                    ScaLBL_Comm->D3Q19_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
                case 2:
                    Vin  = getBoundaryVoltagefromPeriodicBC(Vin0,freqIn,PhaseShift_In,timestep_from_Study);
                    ScaLBL_Comm->D3Q19_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
            }
        }
        if (BoundaryConditionOutlet > 0){
            switch (BoundaryConditionOutlet){
                case 1:
                    ScaLBL_Comm->D3Q19_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
                case 2:
                    Vout = getBoundaryVoltagefromPeriodicBC(Vout0,freqOut,PhaseShift_Out,timestep_from_Study);
                    ScaLBL_Comm->D3Q19_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
            }
        }
       
        //-------------------------//
        //ScaLBL_D3Q19_AAeven_Poisson_ElectricPotential(dvcMap, fq, ChargeDensity, Psi, epsilon_LB, UseSlippingVelBC, 0, ScaLBL_Comm->LastExterior(), Np);
    }
}

void ScaLBL_Poisson::SolvePoissonAAodd(double *ChargeDensity, bool UseSlippingVelBC, int timestep){
	
    if (lattice_scheme.compare("D3Q7")==0){
        ScaLBL_D3Q7_AAodd_Poisson(NeighborList, dvcMap, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, UseSlippingVelBC, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_D3Q7_AAodd_Poisson(NeighborList, dvcMap, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, UseSlippingVelBC, 0, ScaLBL_Comm->LastExterior(), Np);
        //TODO: perhaps add another ScaLBL_Comm routine to update Psi values on solid boundary nodes.
        //something like:
        //ScaLBL_Comm->SolidDirichletBoundaryUpdates(Psi, Psi_BCLabel, timestep);
        ScaLBL_Comm->SolidDirichletAndNeumannD3Q7(fq, Psi, Psi_BCLabel);
        //if (BoundaryConditionSolid==1){
        //    ScaLBL_Comm->SolidDirichletD3Q7(fq, Psi);
        //}
        //else if (BoundaryConditionSolid==2){
        //    ScaLBL_Comm->SolidNeumannD3Q7(fq, Psi);
        //}
    }
    else if (lattice_scheme.compare("D3Q19")==0){
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAodd_Poisson(NeighborList, dvcMap, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, UseSlippingVelBC, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        //ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        
        // Set boundary conditions
        if (BoundaryConditionInlet > 0){
            switch (BoundaryConditionInlet){
                case 1:
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
                case 2:
                    Vin  = getBoundaryVoltagefromPeriodicBC(Vin0,freqIn,PhaseShift_In,timestep);
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
            }
        }
        if (BoundaryConditionOutlet > 0){
            switch (BoundaryConditionOutlet){
                case 1:
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_Z(NeighborList, fq, Vout, timestep);
                    break;
                case 2:
                    Vout = getBoundaryVoltagefromPeriodicBC(Vout0,freqOut,PhaseShift_Out,timestep);
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_Z(NeighborList, fq, Vout, timestep);
                    break;
            }
        }

        ScaLBL_D3Q19_AAodd_Poisson(NeighborList, dvcMap, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, UseSlippingVelBC, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->Barrier();
        //TODO: perhaps add another ScaLBL_Comm routine to update Psi values on solid boundary nodes.
        //something like:
        //ScaLBL_Comm->SolidDirichletAndNeumannD3Q7(fq, Psi, Psi_BCLabel);
    }
}

void ScaLBL_Poisson::SolvePoissonAAeven(double *ChargeDensity, bool UseSlippingVelBC, int timestep){

    if (lattice_scheme.compare("D3Q7")==0){
        ScaLBL_D3Q7_AAeven_Poisson(dvcMap, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, UseSlippingVelBC, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_D3Q7_AAeven_Poisson(dvcMap, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, UseSlippingVelBC, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->SolidDirichletAndNeumannD3Q7(fq, Psi, Psi_BCLabel);
        //if (BoundaryConditionSolid==1){
        //    ScaLBL_Comm->SolidDirichletD3Q7(fq, Psi);
        //}
        //else if (BoundaryConditionSolid==2){
        //    ScaLBL_Comm->SolidNeumannD3Q7(fq, Psi);
        //}
    }
    else if (lattice_scheme.compare("D3Q19")==0){
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAeven_Poisson(dvcMap, fq, ChargeDensity, Psi, ElectricField, ResidualError, tau, epsilon_LB, UseSlippingVelBC, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        //	ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        
        // Set boundary conditions
        if (BoundaryConditionInlet > 0){
            switch (BoundaryConditionInlet){
                case 1:
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
                case 2:
                    Vin  = getBoundaryVoltagefromPeriodicBC(Vin0,freqIn,PhaseShift_In,timestep);
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
                    break;
            }
        }
        if (BoundaryConditionOutlet > 0){
            switch (BoundaryConditionOutlet){
                case 1:
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_Z(NeighborList, fq, Vout, timestep);
                    break;
                case 2:
                    Vout = getBoundaryVoltagefromPeriodicBC(Vout0,freqOut,PhaseShift_Out,timestep);
                    ScaLBL_Comm->D3Q7_Poisson_Potential_BC_Z(NeighborList, fq, Vout, timestep);
                    break;
            }
        }
        
        ScaLBL_D3Q19_AAeven_Poisson(dvcMap, fq, ChargeDensity, Psi, ElectricField, ResidualError, tau, epsilon_LB, UseSlippingVelBC, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->Barrier();
        
        //ScaLBL_Comm->SolidDirichletAndNeumannD3Q7(fq, Psi, Psi_BCLabel);
    }
}

void ScaLBL_Poisson::Checkpoint(){

	if (rank == 0) {
		printf("   POISSON MODEL: Writing restart file! \n");
	}
	double value;
	double *cPsi;
	cPsi = new double[Nx*Ny*Nz];
	ScaLBL_CopyToHost(cPsi, Psi, Nx*Ny*Nz *sizeof(double));

	ofstream File(LocalRestartFile, ios::binary);
	for (int n = 0; n < Nx*Ny*Nz; n++) {
		value = cPsi[n];
		File.write((char *)&value, sizeof(value));
	}

	File.close();

	delete[] cPsi;
}

void ScaLBL_Poisson::DummyChargeDensity(){
    double *ChargeDensity_host;
    ChargeDensity_host = new double[Np];

	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				int idx=Map(i,j,k);
				if (!(idx < 0))
                    ChargeDensity_host[idx] = chargeDen_dummy*(h*h*h*1.0e-18);
            }
        }
    }
	ScaLBL_AllocateDeviceMemory((void **) &ChargeDensityDummy, sizeof(double)*Np);
	ScaLBL_CopyToDevice(ChargeDensityDummy, ChargeDensity_host, sizeof(double)*Np);
	ScaLBL_Comm->Barrier();
	delete [] ChargeDensity_host;
}

void ScaLBL_Poisson::getElectricPotential_debug(int timestep){
    //This function write out decomposed data
    DoubleArray PhaseField(Nx,Ny,Nz);
	//ScaLBL_Comm->RegularLayout(Map,Psi,PhaseField);
    ScaLBL_CopyToHost(PhaseField.data(),Psi,sizeof(double)*Nx*Ny*Nz);
    //ScaLBL_Comm->Barrier(); comm.barrier();
    FILE *OUTFILE;
    sprintf(LocalRankFilename,"Electric_Potential_Time_%i.%05i.raw",timestep,rank);
    OUTFILE = fopen(LocalRankFilename,"wb");
    fwrite(PhaseField.data(),8,N,OUTFILE);
    fclose(OUTFILE);
}

void ScaLBL_Poisson::getElectricPotential(DoubleArray &ReturnValues){
    //This function wirte out the data in a normal layout (by aggregating all decomposed domains)
	//ScaLBL_Comm->RegularLayout(Map,Psi,PhaseField);
    ScaLBL_CopyToHost(ReturnValues.data(),Psi,sizeof(double)*Nx*Ny*Nz);
}

void ScaLBL_Poisson::getElectricField(DoubleArray &Values_x, DoubleArray &Values_y, DoubleArray &Values_z){

	    ScaLBL_Comm->RegularLayout(Map,&ElectricField[0*Np],Values_x);
        ElectricField_LB_to_Phys(Values_x);
        ScaLBL_Comm->Barrier(); comm.barrier();

	    ScaLBL_Comm->RegularLayout(Map,&ElectricField[1*Np],Values_y);
        ElectricField_LB_to_Phys(Values_y);
        ScaLBL_Comm->Barrier(); comm.barrier();

	    ScaLBL_Comm->RegularLayout(Map,&ElectricField[2*Np],Values_z);
        ElectricField_LB_to_Phys(Values_z);
        ScaLBL_Comm->Barrier(); comm.barrier();

}

void ScaLBL_Poisson::getElectricField_debug(int timestep){

        //ScaLBL_D3Q7_Poisson_getElectricField(fq,ElectricField,tau,Np);
        //ScaLBL_Comm->Barrier(); comm.barrier();

        DoubleArray PhaseField(Nx,Ny,Nz);
	    ScaLBL_Comm->RegularLayout(Map,&ElectricField[0*Np],PhaseField);
        ElectricField_LB_to_Phys(PhaseField);
        FILE *EX;
        sprintf(LocalRankFilename,"ElectricField_X_Time_%i.%05i.raw",timestep,rank);
        EX = fopen(LocalRankFilename,"wb");
        fwrite(PhaseField.data(),8,N,EX);
        fclose(EX);

	    ScaLBL_Comm->RegularLayout(Map,&ElectricField[1*Np],PhaseField);
        ElectricField_LB_to_Phys(PhaseField);
        FILE *EY;
        sprintf(LocalRankFilename,"ElectricField_Y_Time_%i.%05i.raw",timestep,rank);
        EY = fopen(LocalRankFilename,"wb");
        fwrite(PhaseField.data(),8,N,EY);
        fclose(EY);

	    ScaLBL_Comm->RegularLayout(Map,&ElectricField[2*Np],PhaseField);
        ElectricField_LB_to_Phys(PhaseField);
        FILE *EZ;
        sprintf(LocalRankFilename,"ElectricField_Z_Time_%i.%05i.raw",timestep,rank);
        EZ = fopen(LocalRankFilename,"wb");
        fwrite(PhaseField.data(),8,N,EZ);
        fclose(EZ);
}

void ScaLBL_Poisson::ElectricField_LB_to_Phys(DoubleArray &Efield_reg){
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
                int idx=Map(i,j,k);
				if (!(idx < 0)){
                    Efield_reg(i,j,k) = Efield_reg(i,j,k)/(h*1.0e-6); 
                }
            }
        }
    }
}

void ScaLBL_Poisson::WriteVis( int timestep) {

    auto vis_db = db->getDatabase("Visualization");
    auto format = vis_db->getWithDefault<string>( "format", "hdf5" );

    DoubleArray ElectricalPotential(Nx, Ny, Nz);
    std::vector<IO::MeshDataStruct> visData;
    fillHalo<double> fillData(Dm->Comm, Dm->rank_info,
                              {Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2}, {1, 1, 1},
                              0, 1);

    IO::initialize("",format,"false");
    // Create the MeshDataStruct    
    visData.resize(1);

    visData[0].meshName = "domain";
    visData[0].mesh =
        std::make_shared<IO::DomainMesh>(Dm->rank_info, Dm->Nx - 2, Dm->Ny - 2,
                                         Dm->Nz - 2, Dm->Lx, Dm->Ly, Dm->Lz);
    //electric potential
    auto ElectricPotentialVar = std::make_shared<IO::Variable>();

    //--------------------------------------------------------------------------------------------------------------------

    //-------------------------------------Create Names for Variables------------------------------------------------------
    if (vis_db->getWithDefault<bool>("save_electric_potential", true)) {
        ElectricPotentialVar->name = "ElectricPotential";
        ElectricPotentialVar->type = IO::VariableType::VolumeVariable;
        ElectricPotentialVar->dim = 1;
        ElectricPotentialVar->data.resize(Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2);
        visData[0].vars.push_back(ElectricPotentialVar);
    }
    //--------------------------------------------------------------------------------------------------------------------

    //------------------------------------Save All Variables--------------------------------------------------------------
    if (vis_db->getWithDefault<bool>("save_electric_potential", true)) {
        ASSERT(visData[0].vars[0]->name == "ElectricPotential");
        getElectricPotential(ElectricalPotential);
        Array<double> &ElectricPotentialData = visData[0].vars[0]->data;
        fillData.copy(ElectricalPotential, ElectricPotentialData);
    }

    if (vis_db->getWithDefault<bool>("write_silo", true))
        IO::writeData(timestep, visData, Dm->Comm);
    //--------------------------------------------------------------------------------------------------------------------
}

//void ScaLBL_Poisson::SolveElectricField(){
//    ScaLBL_Comm_Regular->SendHalo(Psi);
//    ScaLBL_D3Q7_Poisson_ElectricField(NeighborList, dvcMap, dvcID, Psi, ElectricField, BoundaryConditionSolid,
//                                      Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
//    ScaLBL_Comm_Regular->RecvHalo(Psi);
//    ScaLBL_Comm->Barrier();
//    if (BoundaryCondition == 1){
//        ScaLBL_Comm->Poisson_D3Q7_BC_z(dvcMap,Psi,Vin);
//        ScaLBL_Comm->Poisson_D3Q7_BC_Z(dvcMap,Psi,Vout);
//    }
//    ScaLBL_D3Q7_Poisson_ElectricField(NeighborList, dvcMap, dvcID, Psi, ElectricField, BoundaryConditionSolid, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np);
//
//}

//void ScaLBL_Poisson::getElectricPotential(){
//
//        DoubleArray PhaseField(Nx,Ny,Nz);
//	    ScaLBL_Comm->RegularLayout(Map,Psi,PhaseField);
//        //ScaLBL_Comm->Barrier(); comm.barrier();
//        FILE *OUTFILE;
//        sprintf(LocalRankFilename,"Electric_Potential.%05i.raw",rank);
//        OUTFILE = fopen(LocalRankFilename,"wb");
//        fwrite(PhaseField.data(),8,N,OUTFILE);
//        fclose(OUTFILE);
//}

//old version where Psi is of size Np
//void ScaLBL_Poisson::AssignSolidBoundary(double *poisson_solid)
//{
//	size_t NLABELS=0;
//	signed char VALUE=0;
//	double AFFINITY=0.f;
//
//	auto LabelList = electric_db->getVector<int>( "SolidLabels" );
//	auto AffinityList = electric_db->getVector<double>( "SolidValues" );
//
//	NLABELS=LabelList.size();
//	if (NLABELS != AffinityList.size()){
//		ERROR("Error: LB-Poisson Solver: SolidLabels and SolidValues must be the same length! \n");
//	}
//
//	double label_count[NLABELS];
//	double label_count_global[NLABELS];
//	// Assign the labels
//
//	for (size_t idx=0; idx<NLABELS; idx++) label_count[idx]=0;
//
//	for (int k=0;k<Nz;k++){
//		for (int j=0;j<Ny;j++){
//			for (int i=0;i<Nx;i++){
//				int n = k*Nx*Ny+j*Nx+i;
//				VALUE=Mask->id[n];
//                AFFINITY=0.f;
//				// Assign the affinity from the paired list
//				for (unsigned int idx=0; idx < NLABELS; idx++){
//				      //printf("idx=%i, value=%i, %i, \n",idx, VALUE,LabelList[idx]);
//					if (VALUE == LabelList[idx]){
//						AFFINITY=AffinityList[idx];
//                        //NOTE need to convert the user input phys unit to LB unit
//                        if (BoundaryConditionSolid==2){
//                            //for BCS=1, i.e. Dirichlet-type, no need for unit conversion
//                            //TODO maybe there is a factor of gamm missing here ?
//                            AFFINITY = AFFINITY*(h*h*1.0e-12)/epsilon_LB; 
//                        }
//						label_count[idx] += 1.0;
//						idx = NLABELS;
//						//Mask->id[n] = 0; // set mask to zero since this is an immobile component
//					}
//				}
//				poisson_solid[n] = AFFINITY;
//			}
//		}
//	}
//
//	for (size_t idx=0; idx<NLABELS; idx++)
//		label_count_global[idx]=Dm->Comm.sumReduce(  label_count[idx]);
//
//	if (rank==0){
//		printf("LB-Poisson Solver: number of Poisson solid labels: %lu \n",NLABELS);
//		for (unsigned int idx=0; idx<NLABELS; idx++){
//			VALUE=LabelList[idx];
//			AFFINITY=AffinityList[idx];
//			double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocs);
//            switch (BoundaryConditionSolid){
//                case 1:
//			        printf("   label=%d, surface potential=%.3g [V], volume fraction=%.2g\n",VALUE,AFFINITY,volume_fraction); 
//                    break;
//                case 2: 
//			        printf("   label=%d, surface charge density=%.3g [C/m^2], volume fraction=%.2g\n",VALUE,AFFINITY,volume_fraction); 
//                    break;
//                default:
//			        printf("   label=%d, surface potential=%.3g [V], volume fraction=%.2g\n",VALUE,AFFINITY,volume_fraction); 
//                    break;
//            }
//		}
//	}
//}

// old version where Psi is of size Np
//void ScaLBL_Poisson::Potential_Init(double *psi_init){
//
//    if (BoundaryCondition==1){
//	    if (electric_db->keyExists( "Vin" )){
//	    	Vin = electric_db->getScalar<double>( "Vin" );
//	    }
//	    if (electric_db->keyExists( "Vout" )){
//	    	Vout = electric_db->getScalar<double>( "Vout" );
//	    }
//    }
//    //By default only periodic BC is applied and Vin=Vout=1.0, i.e. there is no potential gradient along Z-axis
//    double slope = (Vout-Vin)/(Nz-2);
//    double psi_linearized;
//	for (int k=0;k<Nz;k++){
//        if (k==0 || k==1){
//            psi_linearized = Vin;
//        }
//        else if (k==Nz-1 || k==Nz-2){
//            psi_linearized = Vout;
//        }
//        else{
//            psi_linearized = slope*(k-1)+Vin;
//        }
//		for (int j=0;j<Ny;j++){
//			for (int i=0;i<Nx;i++){
//				int idx = Map(i,j,k);
//				if (!(idx < 0)){
//                    psi_init[idx] = psi_linearized;
//                }
//            }
//        }
//    }
//}
