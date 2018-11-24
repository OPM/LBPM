/*
color lattice boltzmann model
 */
#include "models/ColorModel.h"
#include "analysis/distance.h"

ScaLBL_ColorModel::ScaLBL_ColorModel(int RANK, int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tauA(0),tauB(0),rhoA(0),rhoB(0),alpha(0),beta(0),
Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),inletA(0),inletB(0),outletA(0),outletB(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{

}
ScaLBL_ColorModel::~ScaLBL_ColorModel(){

}

/*void ScaLBL_ColorModel::WriteCheckpoint(const char *FILENAME, const double *cPhi, const double *cfq, int Np)
{
    int q,n;
    double value;
    ofstream File(FILENAME,ios::binary);
    for (n=0; n<Np; n++){
        // Write the two density values
        value = cPhi[n];
        File.write((char*) &value, sizeof(value));
        // Write the even distributions
        for (q=0; q<19; q++){
            value = cfq[q*Np+n];
            File.write((char*) &value, sizeof(value));
        }
    }
    File.close();

}

void ScaLBL_ColorModel::ReadCheckpoint(char *FILENAME, double *cPhi, double *cfq, int Np)
{
    int q=0, n=0;
    double value=0;
    ifstream File(FILENAME,ios::binary);
    for (n=0; n<Np; n++){
        File.read((char*) &value, sizeof(value));
        cPhi[n] = value;
        // Read the distributions
        for (q=0; q<19; q++){
            File.read((char*) &value, sizeof(value));
            cfq[q*Np+n] = value;
        }
    }
    File.close();
}
 */


void ScaLBL_ColorModel::ReadParams(string filename){
	// read the input database 
	db = std::make_shared<Database>( filename );
	domain_db = db->getDatabase( "Domain" );
	color_db = db->getDatabase( "Color" );
	analysis_db = db->getDatabase( "Analysis" );

	// Color Model parameters
	timestepMax = color_db->getScalar<int>( "timestepMax" );
	tauA = color_db->getScalar<double>( "tauA" );
	tauB = color_db->getScalar<double>( "tauB" );
	rhoA = color_db->getScalar<double>( "rhoA" );
	rhoB = color_db->getScalar<double>( "rhoB" );
	Fx = color_db->getVector<double>( "F" )[0];
	Fy = color_db->getVector<double>( "F" )[1];
	Fz = color_db->getVector<double>( "F" )[2];
	alpha = color_db->getScalar<double>( "alpha" );
	beta = color_db->getScalar<double>( "beta" );
	Restart = color_db->getScalar<bool>( "Restart" );
	din = color_db->getScalar<double>( "din" );
	dout = color_db->getScalar<double>( "dout" );
	flux = color_db->getScalar<double>( "flux" );
	inletA=1.f;
	inletB=0.f;
	outletA=0.f;
	outletB=1.f;

	// Read domain parameters
	auto L = domain_db->getVector<double>( "L" );
	auto size = domain_db->getVector<int>( "n" );
	auto nproc = domain_db->getVector<int>( "nproc" );
	BoundaryCondition = domain_db->getScalar<int>( "BC" );
	Nx = size[0];
	Ny = size[1];
	Nz = size[2];
	Lx = L[0];
	Ly = L[1];
	Lz = L[2];
	nprocx = nproc[0];
	nprocy = nproc[1];
	nprocz = nproc[2];
	
	if (BoundaryCondition==4) flux *= rhoA; // mass flux must adjust for density (see formulation for details)

}
void ScaLBL_ColorModel::SetDomain(){
	Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
	Mask  = std::shared_ptr<Domain>(new Domain(domain_db,comm));    // mask domain removes immobile phases
	Nx+=2; Ny+=2; Nz += 2;
	N = Nx*Ny*Nz;
	id = new char [N];
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;               // initialize this way
	Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm) ); // TwoPhase analysis object
	MPI_Barrier(comm);
	Dm->CommInit();
	MPI_Barrier(comm);
	rank = Dm->rank();
}

void ScaLBL_ColorModel::ReadInput(){
	size_t readID;
	Mask->ReadIDs();
	for (int i=0; i<Nx*Ny*Nz; i++) id[i] = Mask->id[i];  // save what was read

	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
	
	// Generate the signed distance map
	// Initialize the domain and communication
	Array<char> id_solid(Nx,Ny,Nz);
	int count = 0;
	// Solve for the position of the solid phase
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				int n = k*Nx*Ny+j*Nx+i;
				// Initialize the solid phase
				if (Mask->id[n] > 0)	id_solid(i,j,k) = 1;
				else	     	      	id_solid(i,j,k) = 0;
			}
		}
	}
	// Initialize the signed distance function
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				int n=k*Nx*Ny+j*Nx+i;
				// Initialize distance to +/- 1
				Averages->SDs(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
			}
		}
	}
//	MeanFilter(Averages->SDs);
	if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
	CalcDist(Averages->SDs,id_solid,*Mask);
	
	if (rank == 0) cout << "Domain set." << endl;
	
	Averages->SetParams(rhoA,rhoB,tauA,tauB,Fx,Fy,Fz,alpha);
}

void ScaLBL_ColorModel::AssignComponentLabels(double *phase)
{
	size_t NLABELS=0;
	char VALUE=0;
	double AFFINITY=0.f;

	auto LabelList = color_db->getVector<char>( "ComponentLabels" );
	auto AffinityList = color_db->getVector<double>( "ComponentAffinity" );

	NLABELS=LabelList.size();
	if (NLABELS != AffinityList.size()){
		ERROR("Error: ComponentLabels and ComponentAffinity must be the same length! \n");
	}

	if (rank==0){
		printf("Components labels: %lu \n",NLABELS);
		for (unsigned int idx=0; idx<NLABELS; idx++){
			VALUE=LabelList[idx];
			AFFINITY=AffinityList[idx];
			printf("   label=%i, affinity=%f\n",int(VALUE),AFFINITY); 
		}
	}
	// Assign the labels
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				int n = k*Nx*Ny+j*Nx+i;
				VALUE=id[n];
				// Assign the affinity from the paired list
				for (unsigned int idx=0; idx < NLABELS; idx++){
					//printf("rank=%i, idx=%i, value=%i, %i, \n",rank(),idx, VALUE,LabelList[idx]);
					if (VALUE == LabelList[idx]){
						AFFINITY=AffinityList[idx];
						idx = NLABELS;
						Mask->id[n] = 0; // set mask to zero since this is an immobile component
					}
				}
				// fluid labels are reserved
				if (VALUE == 1) AFFINITY=1.0;
				else if (VALUE == 2) AFFINITY=-1.0;
				phase[n] = AFFINITY;
			}
		}
	}
	// Set Dm to match Mask
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i]; 
}


void ScaLBL_ColorModel::Create(){
	/*
	 *  This function creates the variables needed to run a LBM 
	 */
	//.........................................................
	// don't perform computations at the eight corners
	//id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
	//id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;

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
	ScaLBL_Comm_Regular  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));

	int Npad=(Np/16 + 2)*16;
	if (rank==0)    printf ("Set up memory efficient layout, %i | %i | %i \n", Np, Npad, N);
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
	dist_mem_size = Np*sizeof(double);
	neighborSize=18*(Np*sizeof(int));
	//...........................................................................
	ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
	ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Aq, 7*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Bq, 7*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*Nx*Ny*Nz);		
	ScaLBL_AllocateDeviceMemory((void **) &Pressure, sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &ColorGrad, 3*sizeof(double)*Np);
	//...........................................................................
	// Update GPU data structures
	if (rank==0)	printf ("Setting up device map and neighbor list \n");
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
	// check that TmpMap is valid
	for (int idx=0; idx<ScaLBL_Comm->LastExterior(); idx++){
		int n = TmpMap[idx];
		if (n > Nx*Ny*Nz){
			printf("Bad value! idx=%i \n");
			TmpMap[idx] = Nx*Ny*Nz-1;
		}
	}
	for (int idx=ScaLBL_Comm->FirstInterior(); idx<ScaLBL_Comm->LastInterior(); idx++){
		int n = TmpMap[idx];
		if (n > Nx*Ny*Nz){
			printf("Bad value! idx=%i \n");
			TmpMap[idx] = Nx*Ny*Nz-1;
		}
	}
	ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
	ScaLBL_DeviceBarrier();
	delete [] TmpMap;
	
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	// initialize phi based on PhaseLabel (include solid component labels)
	double *PhaseLabel;
	PhaseLabel = new double[N];
	AssignComponentLabels(PhaseLabel);
	ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));
}        

/********************************************************
 * AssignComponentLabels                                 *
 ********************************************************/

void ScaLBL_ColorModel::Initialize(){
	
	if (rank==0)	printf ("Initializing distributions \n");
	ScaLBL_D3Q19_Init(fq, Np);
	/*
	 * This function initializes model
	 */
	if (Restart == true){
		if (rank==0){
			printf("Reading restart file! \n");
			ifstream restart("Restart.txt");
			if (restart.is_open()){
				restart  >> timestep;
				printf("Restarting from timestep =%i \n",timestep);
			}
			else{
				printf("WARNING:No Restart.txt file, setting timestep=0 \n");
				timestep=0;
			}
		}
		MPI_Bcast(&timestep,1,MPI_INT,0,comm);
		// Read in the restart file to CPU buffers
		int *TmpMap;
		TmpMap = new int[Np];
		
		double *cPhi, *cDist, *cDen;
		cPhi = new double[N];
		cDen = new double[2*Np];
		cDist = new double[19*Np];
		ScaLBL_CopyToHost(TmpMap, dvcMap, Np*sizeof(int));
        ScaLBL_CopyToHost(cPhi, Phi, N*sizeof(double));
    	
		ifstream File(LocalRestartFile,ios::binary);
		int idx;
		double value,va,vb;
		for (int n=0; n<Np; n++){
			File.read((char*) &va, sizeof(va));
			File.read((char*) &vb, sizeof(vb));
			cDen[n]    = va;
			cDen[Np+n] = vb;
		}
		for (int n=0; n<Np; n++){
			// Read the distributions
			for (int q=0; q<19; q++){
				File.read((char*) &value, sizeof(value));
				cDist[q*Np+n] = value;
			}
		}
		File.close();
		
		for (int n=0; n<ScaLBL_Comm->LastExterior(); n++){
			va = cDen[n];
			vb = cDen[Np + n];
			value = (va-vb)/(va+vb);
			idx = TmpMap[n];
			if (!(idx < 0) && idx<N)
				cPhi[idx] = value;
		}
		for (int n=ScaLBL_Comm->FirstInterior(); n<ScaLBL_Comm->LastInterior(); n++){
		  va = cDen[n];
		  vb = cDen[Np + n];
		  	value = (va-vb)/(va+vb);
		  	idx = TmpMap[n];
		  	if (!(idx < 0) && idx<N)
		  		cPhi[idx] = value;
		}
		
		// Copy the restart data to the GPU
		ScaLBL_CopyToDevice(Den,cDen,2*Np*sizeof(double));
		ScaLBL_CopyToDevice(fq,cDist,19*Np*sizeof(double));
		ScaLBL_CopyToDevice(Phi,cPhi,N*sizeof(double));
		ScaLBL_DeviceBarrier();

		MPI_Barrier(comm);
	}

	if (rank==0)	printf ("Initializing phase field \n");
	ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np);
	ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);

	if (BoundaryCondition >0 ){
		if (Dm->kproc()==0){
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
		}
		if (Dm->kproc() == nprocz-1){
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-1);
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-2);
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-3);
		}
	}

}

void ScaLBL_ColorModel::Run(){
	int nprocs=nprocx*nprocy*nprocz;
	const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
	
	bool SET_CAPILLARY_NUMBER = false;
	bool MORPH_ADAPT = false;
	bool USE_MORPH = false;
	int morph_interval;
	double morph_delta;
	int morph_timesteps = 0;
	int ramp_timesteps = 50000;
	double capillary_number;
	double tolerance = 1.f;
	double Ca_previous = 0.f;

	int target_saturation_index=0;
	std::vector<double> target_saturation;
	double TARGET_SATURATION = 0.f;
	if (color_db->keyExists( "target_saturation" )){
		target_saturation = color_db->getVector<double>( "target_saturation" );
		TARGET_SATURATION = target_saturation[0];
	}
	if (color_db->keyExists( "capillary_number" )){
		capillary_number = color_db->getScalar<double>( "capillary_number" );
		SET_CAPILLARY_NUMBER=true;
	}
	else{
		capillary_number=0;
	}
	if (BoundaryCondition != 0 && SET_CAPILLARY_NUMBER==true){
		if (rank == 0) printf("WARINING: capillary number target only supported for BC = 0 \n");
		SET_CAPILLARY_NUMBER=false;
	}
	if (analysis_db->keyExists( "morph_delta" )){
		morph_delta = analysis_db->getScalar<double>( "morph_delta" );
	}
	else{
		morph_delta=0.5;
	}
	if (analysis_db->keyExists( "morph_interval" )){
		morph_interval = analysis_db->getScalar<int>( "morph_interval" );
		USE_MORPH = true;
	}
	else{
		morph_interval=1000000;
		USE_MORPH = false;
	}
	if (analysis_db->keyExists( "tolerance" )){
		tolerance = analysis_db->getScalar<double>( "tolerance" );
	}
	else{
		tolerance = 0.02;
	}
	int analysis_interval = analysis_db->getScalar<int>( "analysis_interval" );
	
	if (rank==0){
		printf("********************************************************\n");
		printf("No. of timesteps: %i \n", timestepMax);
		fflush(stdout);
	}
	//.......create and start timer............
	double starttime,stoptime,cputime;
	ScaLBL_DeviceBarrier();
	MPI_Barrier(comm);
	starttime = MPI_Wtime();
	//.........................................
	
	//************ MAIN ITERATION LOOP ***************************************/
	PROFILE_START("Loop");
    //std::shared_ptr<Database> analysis_db;
	bool Regular = false;
	runAnalysis analysis( analysis_db, rank_info, ScaLBL_Comm, Dm, Np, Regular, beta, Map );
	//analysis.createThreads( analysis_method, 4 );
	while (timestep < timestepMax ) {
		//if ( rank==0 ) { printf("Running timestep %i (%i MB)\n",timestep+1,(int)(Utilities::getMemoryUsage()/1048576)); }
		PROFILE_START("Update");
		// *************ODD TIMESTEP*************
		timestep++;
		// Compute the Phase indicator field
		// Read for Aq, Bq happens in this routine (requires communication)
		ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
		ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();
		ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);
		
		// Perform the collision operation
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
		if (BoundaryCondition > 0){
			ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
			ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
		}
		// Halo exchange for phase field
		ScaLBL_Comm_Regular->SendHalo(Phi);

		ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm_Regular->RecvHalo(Phi);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();
		// Set BCs
		if (BoundaryCondition == 3){
			ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		if (BoundaryCondition == 4){
			din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

		// *************EVEN TIMESTEP*************
		timestep++;
		// Compute the Phase indicator field
		ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
		ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();
		ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);

		// Perform the collision operation
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
		// Halo exchange for phase field
		if (BoundaryCondition > 0){
			ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
			ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
		}
		ScaLBL_Comm_Regular->SendHalo(Phi);
		ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz,  Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm_Regular->RecvHalo(Phi);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();
		// Set boundary conditions
		if (BoundaryCondition == 3){
			ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		else if (BoundaryCondition == 4){
			din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		//************************************************************************
		
		MPI_Barrier(comm);
		PROFILE_STOP("Update");

		// Run the analysis
		analysis.run( timestep, *Averages, Phi, Pressure, Velocity, fq, Den );
		
		// allow initial ramp-up to get closer to steady state
		if (timestep > ramp_timesteps && timestep%analysis_interval == analysis_interval-20 && USE_MORPH){
			if ( morph_timesteps > morph_interval ){

				double volB = Averages->Volume_w(); 
				double volA = Averages->Volume_n(); 
				double vA_x = Averages->van_global(0); 
				double vA_y = Averages->van_global(1); 
				double vA_z = Averages->van_global(2); 
				double vB_x = Averages->vaw_global(0); 
				double vB_y = Averages->vaw_global(1); 
				double vB_z = Averages->vaw_global(2);
				double muA = rhoA*(tauA-0.5)/3.f; 
				double muB = rhoB*(tauB-0.5)/3.f;				
				
				double flow_rate_A = sqrt(vA_x*vA_x + vA_y*vA_y + vA_z*vA_z);
				double flow_rate_B = sqrt(vB_x*vB_x + vB_y*vB_y + vB_z*vB_z);
				double current_saturation = volB/(volA+volB);
				double Ca = fabs(volA*muA*flow_rate_A + volB*muB*flow_rate_B)/(5.796*alpha*double(Nx*Ny*Nz*nprocs));

				double force_magnitude = sqrt(Fx*Fx + Fy*Fy + Fz*Fz);
				//double krA = muA*volA*flow_rate_A/force_magnitude/double(Nx*Ny*Nz*nprocs);
				//double krB = muB*volB*flow_rate_B/force_magnitude/double(Nx*Ny*Nz*nprocs);

				if (fabs((Ca - Ca_previous)/Ca) < tolerance ){
					MORPH_ADAPT = true;
					if (rank==0){
						printf("** WRITE STEADY POINT *** ");
						printf("Ca = %f, (previous = %f) \n",Ca,Ca_previous);
						volA /= double(Nx*Ny*Nz*nprocs);
						volB /= double(Nx*Ny*Nz*nprocs);
						FILE * kr_log_file = fopen("relperm.csv","a");
						fprintf(kr_log_file,"%i %.5g %.5g %.5g %.5g %.5g %.5g ",timestep-analysis_interval+20,muA,muB,5.796*alpha,Fx,Fy,Fz);
						fprintf(kr_log_file,"%.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g\n",volA,volB,vA_x,vA_y,vA_z,vB_x,vB_y,vB_z);
						fclose(kr_log_file);

						printf("  Measured capillary number %f \n ",Ca);
					}

					if (SET_CAPILLARY_NUMBER ){
						Fx *= capillary_number / Ca;
						Fy *= capillary_number / Ca;
						Fz *= capillary_number / Ca;

						if (force_magnitude > 1e-3){
							Fx *= 1e-3/force_magnitude;   // impose ceiling for stability
							Fy *= 1e-3/force_magnitude;   
							Fz *= 1e-3/force_magnitude;   
						}
						if (force_magnitude < 1e-6){
							Fx *= 1e-6/force_magnitude;   // impose floor
							Fy *= 1e-6/force_magnitude;   
							Fz *= 1e-6/force_magnitude;   
						}
						if (rank == 0) printf("    -- adjust force by factor %f \n ",capillary_number / Ca);
						Averages->SetParams(rhoA,rhoB,tauA,tauB,Fx,Fy,Fz,alpha);
					}

					if (morph_delta > 0.f){
						// wetting phase saturation will decrease
						while (current_saturation < TARGET_SATURATION && target_saturation_index < target_saturation.size() ){
							TARGET_SATURATION = target_saturation[target_saturation_index++];
						}
					}
					else{
						// wetting phase saturation will increase
						while (current_saturation > TARGET_SATURATION && target_saturation_index < target_saturation.size() ){
							TARGET_SATURATION = target_saturation[target_saturation_index++];
							if (rank==0) printf("   Set target saturation as %f (currently %f)\n",TARGET_SATURATION,current_saturation);
						}
					}
				}
				else{
					if (rank==0){
						printf("** Continue to simulate steady *** \n ");
						printf("Ca = %f, (previous = %f) \n",Ca,Ca_previous);
					}
					morph_timesteps=0;
				}
				Ca_previous = Ca;
			}
			if (MORPH_ADAPT ){
				if (rank==0) printf("***Morphological step with target saturation %f ***\n",TARGET_SATURATION);
				double volB = Averages->Volume_w(); 
				double volA = Averages->Volume_n(); 
				double delta_volume = MorphInit(beta,morph_delta);
				double delta_volume_target = volB - (volA + volB)*TARGET_SATURATION; // change in volume to A
				// update the volume
				volA += delta_volume;
				volB -= delta_volume;
				if ((delta_volume_target - delta_volume) / delta_volume > 0.f){
					morph_delta *= 1.01*min((delta_volume_target - delta_volume) / delta_volume, 2.0);
					if (morph_delta > 1.f) morph_delta = 1.f;
					if (morph_delta < -1.f) morph_delta = -1.f;
					if (fabs(morph_delta) < 0.05 ) morph_delta = 0.05*(morph_delta)/fabs(morph_delta); // set minimum
					if (rank==0) printf("  Adjust morph delta: %f \n", morph_delta);
				}
				//MORPH_ADAPT = false;
				if (volB/(volA + volB) > TARGET_SATURATION){
					MORPH_ADAPT = false;
					TARGET_SATURATION = target_saturation[target_saturation_index++];
				}
				MPI_Barrier(comm);
				morph_timesteps = 0;
			}
			morph_timesteps += analysis_interval;
		}
	}
	analysis.finish();
	PROFILE_STOP("Loop");
	PROFILE_SAVE("lbpm_color_simulator",1);
	//************************************************************************
	ScaLBL_DeviceBarrier();
	MPI_Barrier(comm);
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

	// ************************************************************************
}

double ScaLBL_ColorModel::MorphInit(const double beta, const double morph_delta){
	const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);

	double vF = 0.f;
	double vS = 0.f;

	DoubleArray phase(Nx,Ny,Nz);
	IntArray phase_label(Nx,Ny,Nz);;
	DoubleArray phase_distance(Nx,Ny,Nz);
	Array<char> phase_id(Nx,Ny,Nz);

	// Basic algorithm to 
	// 1. Copy phase field to CPU
	ScaLBL_CopyToHost(phase.data(), Phi, N*sizeof(double));

	double count,count_global,volume_initial,volume_final;
	count = 0.f;
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				if (phase(i,j,k) > 0.f && Averages->SDs(i,j,k) > 0.f) count+=1.f;
			}
		}
	}
	MPI_Allreduce(&count,&count_global,1,MPI_DOUBLE,MPI_SUM,comm);
	volume_initial = count_global;

	// 2. Identify connected components of phase field -> phase_label
	BlobIDstruct new_index;
	ComputeGlobalBlobIDs(Nx-2,Ny-2,Nz-2,rank_info,phase,Averages->SDs,vF,vS,phase_label,comm);
	MPI_Barrier(comm);
	// only operate on component "0"
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				int label = phase_label(i,j,k);
				if (label == 0 )     phase_id(i,j,k) = 0;
				else 		     phase_id(i,j,k) = 1;
			}
		}
	}	
	// 3. Generate a distance map to the largest object -> phase_distance
	CalcDist(phase_distance,phase_id,*Dm);

	double temp,value;
	double factor=0.5/beta;
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				if (phase_distance(i,j,k) < 3.f ){
					value = phase(i,j,k);
					if (value > 1.f)   value=1.f;
					if (value < -1.f)  value=-1.f;
					// temp -- distance based on analytical form McClure, Prins et al, Comp. Phys. Comm.
					temp = -factor*log((1.0+value)/(1.0-value));
					/// use this approximation close to the object
					if (fabs(value) < 0.8 && Averages->SDs(i,j,k) > 1.f ){
						phase_distance(i,j,k) = temp;
					}
				}
			}
		}
	}

	// 4. Apply erosion / dilation operation to phase_distance
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				double walldist=Averages->SDs(i,j,k);
				double wallweight = 1.f / (1+exp(-5.f*(walldist-1.f))); 
				phase_distance(i,j,k) -= wallweight*morph_delta;
			}
		}
	}

	// 5. Update phase indicator field based on new distnace
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				int n = k*Nx*Ny + j*Nx + i;
				double d = phase_distance(i,j,k);
				if (Averages->SDs(i,j,k) > 0.f){
					// only update phase field in immediate proximity of largest component
					if (d < 3.f){
						phase(i,j,k) = (2.f*(exp(-2.f*beta*d))/(1.f+exp(-2.f*beta*d))-1.f);
					}
				}
			} 
		}
	}

	count = 0.f;
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				if (phase(i,j,k) > 0.f && Averages->SDs(i,j,k) > 0.f) count+=1.f;
			}
		}
	}
	MPI_Allreduce(&count,&count_global,1,MPI_DOUBLE,MPI_SUM,comm);
	volume_final=count_global;

	double delta_volume = (volume_final-volume_initial);
	if (rank == 0)  printf("MorphInit: change fluid volume fraction by %f \n", delta_volume/volume_initial);

	// 6. copy back to the device
	//if (rank==0)  printf("MorphInit: copy data  back to device\n");
	ScaLBL_CopyToDevice(Phi,phase.data(),N*sizeof(double));

	// 7. Re-initialize phase field and density
	ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np);
	ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
	if (BoundaryCondition >0 ){
		if (Dm->kproc()==0){
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
		}
		if (Dm->kproc() == nprocz-1){
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-1);
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-2);
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-3);
		}
	}
	return delta_volume;
}

void ScaLBL_ColorModel::WriteDebug(){
	// Copy back final phase indicator field and convert to regular layout
	DoubleArray PhaseField(Nx,Ny,Nz);
	//ScaLBL_Comm->RegularLayout(Map,Phi,PhaseField);
	ScaLBL_CopyToHost(PhaseField.data(), Phi, sizeof(double)*N);

	FILE *OUTFILE;
	sprintf(LocalRankFilename,"Phase.%05i.raw",rank);
	OUTFILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,OUTFILE);
	fclose(OUTFILE);
}
