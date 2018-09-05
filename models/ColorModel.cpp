/*
color lattice boltzmann model
 */
#include "models/ColorModel.h"

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

	if (BoundaryCondition==4) flux = din*rhoA; // mass flux must adjust for density (see formulation for details)

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
	//.......................................................................
	if (rank == 0)    printf("Read input media... \n");
	//.......................................................................
	Mask->ReadIDs();
	for (int i=0; i<Nx*Ny*Nz; i++) id[i] = Mask->id[i];  // save what was read

	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);

	// .......... READ THE INPUT FILE .......................................
	//...........................................................................
	if (rank == 0) cout << "Reading in signed distance function..." << endl;
	//.......................................................................
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
	ReadBinaryFile(LocalRankFilename, Averages->SDs.data(), N);
	MPI_Barrier(comm);
	if (rank == 0) cout << "Domain set." << endl;
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
	for (int idx=0; idx<ScaLBL_Comm->LastInterior(); idx++){
		if (idx == ScaLBL_Comm->LastExterior()) idx = ScaLBL_Comm->FirstInterior();
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
		double *cPhi = new double[Np];
		double *cDist = new double[19*Np];
		ifstream File(LocalRestartFile,ios::binary);
		double value,va,vb;
		for (int n=0; n<Np; n++){
			File.read((char*) &value, sizeof(va));
			File.read((char*) &value, sizeof(vb));
			value = (va-vb)/(va+vb);
			cPhi[n] = value;
		}
		for (int n=0; n<Np; n++){
			// Read the distributions
			for (int q=0; q<19; q++){
				File.read((char*) &value, sizeof(value));
				cDist[q*Np+n] = value;
			}
		}
		File.close();
		// Copy the restart data to the GPU
		ScaLBL_CopyToDevice(fq,cDist,19*Np*sizeof(double));
		ScaLBL_CopyToDevice(Phi,cPhi,Np*sizeof(double));
		ScaLBL_DeviceBarrier();
		delete [] cPhi;
		delete [] cDist;
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
		ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);
		
		// Perform the collision operation
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
		// Halo exchange for phase field
		ScaLBL_Comm_Regular->SendHalo(Phi);

		ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm_Regular->RecvHalo(Phi);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		// Set BCs
		if (BoundaryCondition > 0){
			ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
			ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
		}
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
		ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);

		// Perform the collision operation
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
		// Halo exchange for phase field
		ScaLBL_Comm_Regular->SendHalo(Phi);
		ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz,  Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm_Regular->RecvHalo(Phi);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		// Set boundary conditions
		if (BoundaryCondition > 0){
			ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
			ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
		}
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
