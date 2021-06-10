 /*
color lattice boltzmann model
 */
#include "models/ColorModel.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"
#include "common/Communication.h"
#include "common/ReadMicroCT.h"
#include <stdlib.h>
#include <time.h>


ScaLBL_ColorModel::ScaLBL_ColorModel(int RANK, int NP, const Utilities::MPI& COMM):
    rank(RANK), nprocs(NP), Restart(0), timestep(0), timestepMax(0),
    tauA(0), tauB(0), rhoA(0), rhoB(0), alpha(0), beta(0),
    Fx(0), Fy(0), Fz(0), flux(0), din(0), dout(0),
    inletA(0), inletB(0), outletA(0), outletB(0),
    Nx(0), Ny(0), Nz(0), N(0), Np(0), nprocx(0), nprocy(0), nprocz(0),
    BoundaryCondition(0), Lx(0), Ly(0), Lz(0), id(nullptr),
    NeighborList(nullptr), dvcMap(nullptr), fq(nullptr), Aq(nullptr), Bq(nullptr),
    Den(nullptr), Phi(nullptr), ColorGrad(nullptr), Velocity(nullptr), Pressure(nullptr),
    comm(COMM)
{
	REVERSE_FLOW_DIRECTION = false;
}
ScaLBL_ColorModel::~ScaLBL_ColorModel()
{
    delete [] id;
	ScaLBL_FreeDeviceMemory( NeighborList );
	ScaLBL_FreeDeviceMemory( dvcMap );
	ScaLBL_FreeDeviceMemory( fq );
	ScaLBL_FreeDeviceMemory( Aq );
	ScaLBL_FreeDeviceMemory( Bq );
	ScaLBL_FreeDeviceMemory( Den );
	ScaLBL_FreeDeviceMemory( Phi );		
	ScaLBL_FreeDeviceMemory( Pressure );
	ScaLBL_FreeDeviceMemory( Velocity );
	ScaLBL_FreeDeviceMemory( ColorGrad );
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
	color_db =  db->getDatabase( "Color" );
	analysis_db = db->getDatabase( "Analysis" );
	vis_db = db->getDatabase( "Visualization" );

	// set defaults
	timestepMax = 100000;
	tauA = tauB = 1.0;
	rhoA = rhoB = 1.0;
	Fx = Fy = Fz = 0.0;
	alpha=1e-3;
	beta=0.95;
	Restart=false;
	din=dout=1.0;
	flux=0.0;
	
	// Color Model parameters
	if (color_db->keyExists( "timestepMax" )){
		timestepMax = color_db->getScalar<int>( "timestepMax" );
	}
	if (color_db->keyExists( "tauA" )){
		tauA = color_db->getScalar<double>( "tauA" );
	}
	if (color_db->keyExists( "tauB" )){
		tauB = color_db->getScalar<double>( "tauB" );
	}
	if (color_db->keyExists( "rhoA" )){
		rhoA = color_db->getScalar<double>( "rhoA" );
	}
	if (color_db->keyExists( "rhoB" )){
		rhoB = color_db->getScalar<double>( "rhoB" );
	}
	if (color_db->keyExists( "F" )){
		Fx = color_db->getVector<double>( "F" )[0];
		Fy = color_db->getVector<double>( "F" )[1];
		Fz = color_db->getVector<double>( "F" )[2];
	}
	if (color_db->keyExists( "alpha" )){
		alpha = color_db->getScalar<double>( "alpha" );
	}
	if (color_db->keyExists( "beta" )){
		beta = color_db->getScalar<double>( "beta" );
	}
	if (color_db->keyExists( "Restart" )){
		Restart = color_db->getScalar<bool>( "Restart" );
	}
	if (color_db->keyExists( "din" )){
		din = color_db->getScalar<double>( "din" );
	}
	if (color_db->keyExists( "dout" )){
		dout = color_db->getScalar<double>( "dout" );
	}
	if (color_db->keyExists( "flux" )){
		flux = color_db->getScalar<double>( "flux" );
	}
	inletA=1.f;
	inletB=0.f;
	outletA=0.f;
	outletB=1.f;
	//if (BoundaryCondition==4) flux *= rhoA; // mass flux must adjust for density (see formulation for details)

	BoundaryCondition = 0;
	if (color_db->keyExists( "BC" )){
		BoundaryCondition = color_db->getScalar<int>( "BC" );
	}
	else if (domain_db->keyExists( "BC" )){
		BoundaryCondition = domain_db->getScalar<int>( "BC" );
	}
	if (domain_db->keyExists( "InletLayersPhase" )){
		int inlet_layers_phase = domain_db->getScalar<int>( "InletLayersPhase" );
		if (inlet_layers_phase == 2 ) {
		  inletA = 0.0;
		  inletB = 1.0;
		}
	}
	if (domain_db->keyExists( "OutletLayersPhase" )){
		int outlet_layers_phase = domain_db->getScalar<int>( "OutletLayersPhase" );
		if (outlet_layers_phase == 1 ) {
		  inletA = 1.0;
		  inletB = 0.0;
		}
	}
	
	// Override user-specified boundary condition for specific protocols
	auto protocol = color_db->getWithDefault<std::string>( "protocol", "none" );
	if (protocol == "seed water"){
		if (BoundaryCondition != 0 && BoundaryCondition != 5){
			BoundaryCondition = 0;
			if (rank==0) printf("WARNING: protocol (seed water) supports only full periodic boundary condition \n");
		}
		domain_db->putScalar<int>( "BC", BoundaryCondition );
	}
	else if (protocol == "open connected oil"){
		if (BoundaryCondition != 0 && BoundaryCondition != 5){
			BoundaryCondition = 0;
			if (rank==0) printf("WARNING: protocol (open connected oil) supports only full periodic boundary condition \n");
		}
		domain_db->putScalar<int>( "BC", BoundaryCondition );
	}
	else if (protocol == "shell aggregation"){
		if (BoundaryCondition != 0 && BoundaryCondition != 5){
			BoundaryCondition = 0;
			if (rank==0) printf("WARNING: protocol (shell aggregation) supports only full periodic boundary condition \n");
		}
		domain_db->putScalar<int>( "BC", BoundaryCondition );
	}  
	else if (protocol == "fractional flow"){
		if (BoundaryCondition != 0 && BoundaryCondition != 5){
			BoundaryCondition = 0;
			if (rank==0) printf("WARNING: protocol (fractional flow) supports only full periodic boundary condition \n");
		}
		domain_db->putScalar<int>( "BC", BoundaryCondition );
	}  	
}

void ScaLBL_ColorModel::SetDomain(){
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
	id = new signed char [N];
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;               // initialize this way
	//Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm) ); // TwoPhase analysis object
	Averages = std::shared_ptr<SubPhase> ( new SubPhase(Dm) ); // TwoPhase analysis object
	comm.barrier();
	Dm->CommInit();
	comm.barrier();
	// Read domain parameters
	rank = Dm->rank();	
	nprocx = Dm->nprocx();
	nprocy = Dm->nprocy();
	nprocz = Dm->nprocz();
}

void ScaLBL_ColorModel::ReadInput(){
	
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
	
	if (color_db->keyExists( "image_sequence" )){
		auto ImageList = color_db->getVector<std::string>( "image_sequence");
		int IMAGE_INDEX = color_db->getWithDefault<int>( "image_index", 0 );
		std::string first_image = ImageList[IMAGE_INDEX];
		Mask->Decomp(first_image);
		IMAGE_INDEX++;
	}
	else if (domain_db->keyExists( "GridFile" )){
        // Read the local domain data
	    auto input_id = readMicroCT( *domain_db, MPI_COMM_WORLD );
        // Fill the halo (assuming GCW of 1)
        array<int,3> size0 = { (int) input_id.size(0), (int) input_id.size(1), (int) input_id.size(2) };
        ArraySize size1 = { (size_t) Mask->Nx, (size_t) Mask->Ny, (size_t) Mask->Nz };
        ASSERT( (int) size1[0] == size0[0]+2 && (int) size1[1] == size0[1]+2 && (int) size1[2] == size0[2]+2 );
        fillHalo<signed char> fill( MPI_COMM_WORLD, Mask->rank_info, size0, { 1, 1, 1 }, 0, 1 );
        Array<signed char> id_view;
        id_view.viewRaw( size1, Mask->id.data() );
        fill.copy( input_id, id_view );
        fill.fill( id_view );
	}
	else if (domain_db->keyExists( "Filename" )){
		auto Filename = domain_db->getScalar<std::string>( "Filename" );
		Mask->Decomp(Filename);
	}
	else{
		Mask->ReadIDs();
	}
	for (int i=0; i<Nx*Ny*Nz; i++) id[i] = Mask->id[i];  // save what was read
	
	// Generate the signed distance map
	// Initialize the domain and communication
	Array<char> id_solid(Nx,Ny,Nz);
	// Solve for the position of the solid phase
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				int n = k*Nx*Ny+j*Nx+i;
				// Initialize the solid phase
				signed char label = Mask->id[n];
				if (label > 0)		id_solid(i,j,k) = 1;
				else	     		id_solid(i,j,k) = 0;
			}
		}
	}
	// Initialize the signed distance function
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				// Initialize distance to +/- 1
				Averages->SDs(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
			}
		}
	}
//	MeanFilter(Averages->SDs);
	Minkowski Solid(Dm);
	if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
	CalcDist(Averages->SDs,id_solid,*Mask);
	Solid.ComputeScalar(Averages->SDs,0.0);
	/* save averages */
	Averages->solid.V = Solid.Vi;
	Averages->solid.A = Solid.Ai;
	Averages->solid.H = Solid.Ji;
	Averages->solid.X = Solid.Xi;
	Averages->gsolid.V = Solid.Vi_global;
	Averages->gsolid.A = Solid.Ai_global;
	Averages->gsolid.H = Solid.Ji_global;
	Averages->gsolid.X = Solid.Xi_global;
	/* write to file */
	if (rank == 0) {
		FILE *SOLID = fopen("solid.csv","w");
		fprintf(SOLID,"Vs As Hs Xs\n");
		fprintf(SOLID,"%.8g %.8g %.8g %.8g\n",Solid.Vi_global,Solid.Ai_global,Solid.Ji_global,Solid.Xi_global);
		fclose(SOLID);
	}
	if (rank == 0) cout << "Domain set." << endl;
	
	Averages->SetParams(rhoA,rhoB,tauA,tauB,Fx,Fy,Fz,alpha,beta);
}

void ScaLBL_ColorModel::AssignComponentLabels(double *phase)
{
	size_t NLABELS=0;
	signed char VALUE=0;
	double AFFINITY=0.f;

	auto LabelList = color_db->getVector<int>( "ComponentLabels" );
	auto AffinityList = color_db->getVector<double>( "ComponentAffinity" );
	auto WettingConvention = color_db->getWithDefault<std::string>( "WettingConvention", "none" );

	NLABELS=LabelList.size();
	if (NLABELS != AffinityList.size()){
		ERROR("Error: ComponentLabels and ComponentAffinity must be the same length! \n");
	}

	if (WettingConvention == "SCAL"){
	  for (size_t idx=0; idx<NLABELS; idx++) AffinityList[idx] *= -1.0;
	}
	
	double * label_count;
	double * label_count_global;
	label_count = new double [NLABELS];
	label_count_global = new double [NLABELS];
	// Assign the labels
	for (size_t idx=0; idx<NLABELS; idx++) label_count[idx]=0;

	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				int n = k*Nx*Ny+j*Nx+i;
				VALUE=id[n];
				// Assign the affinity from the paired list
				for (unsigned int idx=0; idx < NLABELS; idx++){
				      //printf("idx=%i, value=%i, %i, \n",idx, VALUE,LabelList[idx]);
					if (VALUE == LabelList[idx]){
						AFFINITY=AffinityList[idx];
						label_count[idx] += 1.0;
						idx = NLABELS;
						//Mask->id[n] = 0; // set mask to zero since this is an immobile component
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
	
	for (size_t idx=0; idx<NLABELS; idx++)
		label_count_global[idx]=Dm->Comm.sumReduce(  label_count[idx]);

	if (rank==0){
		printf("Component labels: %lu \n",NLABELS);
		for (unsigned int idx=0; idx<NLABELS; idx++){
			VALUE=LabelList[idx];
			AFFINITY=AffinityList[idx];
			double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocs);
			printf("   label=%d, affinity=%f, volume fraction==%f\n",VALUE,AFFINITY,volume_fraction); 
		}
	}

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
	Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Mask->id.data(),Np,1);
	comm.barrier();

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
	ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
	ScaLBL_Comm->Barrier();
	delete [] TmpMap;
	
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
    delete [] neighborList;
	// initialize phi based on PhaseLabel (include solid component labels)
	double *PhaseLabel;
	PhaseLabel = new double[N];
	AssignComponentLabels(PhaseLabel);
	ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));
    delete [] PhaseLabel;
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
		}

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
		ScaLBL_Comm->Barrier();

		comm.barrier();
	}

	if (rank==0)	printf ("Initializing phase field \n");
	ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np);
	ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);

	// establish reservoirs for external bC
	if (BoundaryCondition == 1 || BoundaryCondition == 2 ||  BoundaryCondition == 3 || BoundaryCondition == 4 ){
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
	ScaLBL_CopyToHost(Averages->Phi.data(),Phi,N*sizeof(double));
}

double ScaLBL_ColorModel::Run(int returntime){
	int nprocs=nprocx*nprocy*nprocz;
	const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
	//************ MAIN ITERATION LOOP ***************************************/
	comm.barrier();
	PROFILE_START("Loop");
	//std::shared_ptr<Database> analysis_db;
	bool Regular = false;
	bool RESCALE_FORCE = false;
	bool SET_CAPILLARY_NUMBER = false;
	double tolerance = 0.01;
	auto current_db = db->cloneDatabase();
    auto flow_db = db->getDatabase( "FlowAdaptor" );
	int MIN_STEADY_TIMESTEPS = flow_db->getWithDefault<int>( "min_steady_timesteps", 1000000 );
	int MAX_STEADY_TIMESTEPS = flow_db->getWithDefault<int>( "max_steady_timesteps", 1000000 );

	int RESCALE_FORCE_AFTER_TIMESTEP = MAX_STEADY_TIMESTEPS*2;

	double capillary_number = 1.0e-5;
	double Ca_previous = 0.0;
	if (color_db->keyExists( "capillary_number" )){
		capillary_number = color_db->getScalar<double>( "capillary_number" );
		SET_CAPILLARY_NUMBER=true;
	}
	if (color_db->keyExists( "rescale_force_after_timestep" )){
		RESCALE_FORCE_AFTER_TIMESTEP = color_db->getScalar<int>( "rescale_force_after_timestep" );
		RESCALE_FORCE = true;
	}
	if (analysis_db->keyExists( "tolerance" )){
		tolerance = analysis_db->getScalar<double>( "tolerance" );
	}

	runAnalysis analysis( current_db, rank_info, ScaLBL_Comm, Dm, Np, Regular, Map );
	auto t1 = std::chrono::system_clock::now();
	int CURRENT_TIMESTEP = 0;
	int START_TIMESTEP = timestep;
	int EXIT_TIMESTEP = min(timestepMax,returntime);
	while (timestep < EXIT_TIMESTEP ) {
		//if ( rank==0 ) { printf("Running timestep %i (%i MB)\n",timestep+1,(int)(Utilities::getMemoryUsage()/1048576)); }
		PROFILE_START("Update");
		// *************ODD TIMESTEP*************
		timestep++;
		// Compute the Phase indicator field
		// Read for Aq, Bq happens in this routine (requires communication)
		ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
		ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
		ScaLBL_Comm->Barrier();
		ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);

		// Perform the collision operation
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
		if (BoundaryCondition > 0 && BoundaryCondition < 5){
			ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
			ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
		}
		// Halo exchange for phase field
		ScaLBL_Comm_Regular->SendHalo(Phi);

		ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm_Regular->RecvHalo(Phi);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_Comm->Barrier();
		// Set BCs
		if (BoundaryCondition == 3){
			ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		if (BoundaryCondition == 4){
			din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		else if (BoundaryCondition == 5){
			ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
			ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
		}
		ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->Barrier(); 

		// *************EVEN TIMESTEP*************
		timestep++;
		// Compute the Phase indicator field
		ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
		ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
		ScaLBL_Comm->Barrier();
		ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);

		// Perform the collision operation
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
		// Halo exchange for phase field
		if (BoundaryCondition > 0 && BoundaryCondition < 5){
			ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
			ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
		}
		ScaLBL_Comm_Regular->SendHalo(Phi);
		ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz,  Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm_Regular->RecvHalo(Phi);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_Comm->Barrier();
		// Set boundary conditions
		if (BoundaryCondition == 3){
			ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		else if (BoundaryCondition == 4){
			din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		else if (BoundaryCondition == 5){
			ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
			ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
		}
		ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->Barrier(); 
		//************************************************************************
		analysis.basic(timestep, current_db, *Averages, Phi, Pressure, Velocity, fq, Den );		// allow initial ramp-up to get closer to steady state

		CURRENT_TIMESTEP += 2;
		if (CURRENT_TIMESTEP > MIN_STEADY_TIMESTEPS){
			analysis.finish();

			double volB = Averages->gwb.V; 
			double volA = Averages->gnb.V; 
			volA /= Dm->Volume;
			volB /= Dm->Volume;;
			//initial_volume = volA*Dm->Volume;
			double vA_x = Averages->gnb.Px/Averages->gnb.M; 
			double vA_y = Averages->gnb.Py/Averages->gnb.M; 
			double vA_z = Averages->gnb.Pz/Averages->gnb.M; 
			double vB_x = Averages->gwb.Px/Averages->gwb.M; 
			double vB_y = Averages->gwb.Py/Averages->gwb.M; 
			double vB_z = Averages->gwb.Pz/Averages->gwb.M;
			double muA = rhoA*(tauA-0.5)/3.f; 
			double muB = rhoB*(tauB-0.5)/3.f;				
			double force_mag = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
			double dir_x = Fx/force_mag;
			double dir_y = Fy/force_mag;
			double dir_z = Fz/force_mag;
			if (force_mag == 0.0){
				// default to z direction
				dir_x = 0.0;
				dir_y = 0.0;
				dir_z = 1.0;
				force_mag = 1.0;
			}
			double current_saturation = volB/(volA+volB);
			double flow_rate_A = volA*(vA_x*dir_x + vA_y*dir_y + vA_z*dir_z);
			double flow_rate_B = volB*(vB_x*dir_x + vB_y*dir_y + vB_z*dir_z);
			double Ca = fabs(muA*flow_rate_A + muB*flow_rate_B)/(5.796*alpha);


			bool isSteady = false;
			if ( (fabs((Ca - Ca_previous)/Ca) < tolerance &&  CURRENT_TIMESTEP > MIN_STEADY_TIMESTEPS))
				isSteady = true;
			if (CURRENT_TIMESTEP >= MAX_STEADY_TIMESTEPS)
				isSteady = true;
			if (RESCALE_FORCE == true && SET_CAPILLARY_NUMBER == true && CURRENT_TIMESTEP > RESCALE_FORCE_AFTER_TIMESTEP){
				RESCALE_FORCE = false;
				double RESCALE_FORCE_FACTOR = capillary_number / Ca;
				if (RESCALE_FORCE_FACTOR > 2.0) RESCALE_FORCE_FACTOR = 2.0;
				if (RESCALE_FORCE_FACTOR < 0.5) RESCALE_FORCE_FACTOR = 0.5;
				Fx *= RESCALE_FORCE_FACTOR;
				Fy *= RESCALE_FORCE_FACTOR;
				Fz *= RESCALE_FORCE_FACTOR;
				force_mag = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
				if (force_mag > 1e-3){
					Fx *= 1e-3/force_mag;   // impose ceiling for stability
					Fy *= 1e-3/force_mag;   
					Fz *= 1e-3/force_mag;   
				}
				if (rank == 0) printf("    -- adjust force by factor %f \n ",capillary_number / Ca);
				Averages->SetParams(rhoA,rhoB,tauA,tauB,Fx,Fy,Fz,alpha,beta);
				color_db->putVector<double>("F",{Fx,Fy,Fz});
			}
			if ( isSteady ){
				Averages->Full();
				Averages->Write(timestep);
				analysis.WriteVisData(timestep, current_db, *Averages, Phi, Pressure, Velocity, fq, Den );
				analysis.finish();

				if (rank==0){
					printf("** WRITE STEADY POINT *** ");
					printf("Ca = %f, (previous = %f) \n",Ca,Ca_previous);
					double h = Dm->voxel_length;		
					// pressures
					double pA = Averages->gnb.p;
					double pB = Averages->gwb.p;
					double pAc = Averages->gnc.p;
					double pBc = Averages->gwc.p;
					double pAB = (pA-pB)/(h*6.0*alpha);
					double pAB_connected = (pAc-pBc)/(h*6.0*alpha);
					// connected contribution
					double Vol_nc = Averages->gnc.V/Dm->Volume;
					double Vol_wc = Averages->gwc.V/Dm->Volume;
					double Vol_nd = Averages->gnd.V/Dm->Volume;
					double Vol_wd = Averages->gwd.V/Dm->Volume;
					double Mass_n = Averages->gnc.M + Averages->gnd.M;
					double Mass_w = Averages->gwc.M + Averages->gwd.M;
					double vAc_x = Averages->gnc.Px/Mass_n; 
					double vAc_y = Averages->gnc.Py/Mass_n; 
					double vAc_z = Averages->gnc.Pz/Mass_n; 
					double vBc_x = Averages->gwc.Px/Mass_w; 
					double vBc_y = Averages->gwc.Py/Mass_w; 
					double vBc_z = Averages->gwc.Pz/Mass_w;
					// disconnected contribution
					double vAd_x = Averages->gnd.Px/Mass_n; 
					double vAd_y = Averages->gnd.Py/Mass_n; 
					double vAd_z = Averages->gnd.Pz/Mass_n; 
					double vBd_x = Averages->gwd.Px/Mass_w; 
					double vBd_y = Averages->gwd.Py/Mass_w; 
					double vBd_z = Averages->gwd.Pz/Mass_w;

					double flow_rate_A_connected = Vol_nc*(vAc_x*dir_x + vAc_y*dir_y + vAc_z*dir_z);
					double flow_rate_B_connected = Vol_wc*(vBc_x*dir_x + vBc_y*dir_y + vBc_z*dir_z);
					double flow_rate_A_disconnected = (Vol_nd)*(vAd_x*dir_x + vAd_y*dir_y + vAd_z*dir_z);
					double flow_rate_B_disconnected = (Vol_wd)*(vBd_x*dir_x + vBd_y*dir_y + vBd_z*dir_z);

					double kAeff_connected = h*h*muA*flow_rate_A_connected/(force_mag);
					double kBeff_connected = h*h*muB*flow_rate_B_connected/(force_mag);

					double kAeff_disconnected = h*h*muA*flow_rate_A_disconnected/(force_mag);
					double kBeff_disconnected = h*h*muB*flow_rate_B_disconnected/(force_mag);

					double kAeff = h*h*muA*(flow_rate_A)/(force_mag);
					double kBeff = h*h*muB*(flow_rate_B)/(force_mag);

					double viscous_pressure_drop = (rhoA*volA + rhoB*volB)*force_mag;
					double Mobility = muA/muB;

					bool WriteHeader=false;
					FILE * kr_log_file = fopen("relperm.csv","r");
					if (kr_log_file != NULL)
						fclose(kr_log_file);
					else
						WriteHeader=true;
					kr_log_file = fopen("relperm.csv","a");
					if (WriteHeader)
						fprintf(kr_log_file,"timesteps sat.water eff.perm.oil eff.perm.water eff.perm.oil.connected eff.perm.water.connected eff.perm.oil.disconnected eff.perm.water.disconnected cap.pressure cap.pressure.connected pressure.drop Ca M\n");

					fprintf(kr_log_file,"%i %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g\n",CURRENT_TIMESTEP,current_saturation,kAeff,kBeff,kAeff_connected,kBeff_connected,kAeff_disconnected,kBeff_disconnected,pAB,pAB_connected,viscous_pressure_drop,Ca,Mobility);
					fclose(kr_log_file);

					printf("  Measured capillary number %f \n ",Ca);
				}
				if (SET_CAPILLARY_NUMBER ){
					Fx *= capillary_number / Ca;
					Fy *= capillary_number / Ca;
					Fz *= capillary_number / Ca;
					if (force_mag > 1e-3){
						Fx *= 1e-3/force_mag;   // impose ceiling for stability
						Fy *= 1e-3/force_mag;   
						Fz *= 1e-3/force_mag;   
					}
					if (rank == 0) printf("    -- adjust force by factor %f \n ",capillary_number / Ca);
					Averages->SetParams(rhoA,rhoB,tauA,tauB,Fx,Fy,Fz,alpha,beta);
					color_db->putVector<double>("F",{Fx,Fy,Fz});
				}
				else{
					if (rank==0){
						printf("** Continue to simulate steady *** \n ");
						printf("Ca = %f, (previous = %f) \n",Ca,Ca_previous);
					}
				}
			}
		}
	}
	analysis.finish();
	PROFILE_STOP("Update");

	PROFILE_STOP("Loop");
	PROFILE_SAVE("lbpm_color_simulator",1);
	//************************************************************************
	// Compute the walltime per timestep
	auto t2 = std::chrono::system_clock::now();
	double cputime = std::chrono::duration<double>( t2 - t1 ).count() / (timestep - START_TIMESTEP);
	// Performance obtained from each node
	double MLUPS = double(Np)/cputime/1000000;

	if (rank==0) printf("********************************************************\n");
	if (rank==0) printf("CPU time = %f \n", cputime);
	if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
	return(MLUPS);
	MLUPS *= nprocs;

}

void ScaLBL_ColorModel::Run(){
	int nprocs=nprocx*nprocy*nprocz;
	const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
	
	int IMAGE_INDEX = 0;
	int IMAGE_COUNT = 0;
	std::vector<std::string> ImageList;
	bool SET_CAPILLARY_NUMBER = false;
	bool RESCALE_FORCE = false;
	bool MORPH_ADAPT = false;
	bool USE_MORPH = false;
	bool USE_SEED = false;
	bool USE_DIRECT = false;
	bool USE_MORPHOPEN_OIL = false;
	int MAX_MORPH_TIMESTEPS = 50000; // maximum number of LBM timesteps to spend in morphological adaptation routine
	int MIN_STEADY_TIMESTEPS = 100000;
	int MAX_STEADY_TIMESTEPS = 200000;
	int RESCALE_FORCE_AFTER_TIMESTEP = 0;
	int RAMP_TIMESTEPS = 0;//50000;		 // number of timesteps to run initially (to get a reasonable velocity field before other pieces kick in)
	int CURRENT_MORPH_TIMESTEPS=0;   // counter for number of timesteps spent in  morphological adaptation routine (reset each time)
	int CURRENT_STEADY_TIMESTEPS=0;   // counter for number of timesteps spent in  morphological adaptation routine (reset each time)
	int morph_interval = 100000;
	int analysis_interval = 1000; 	// number of timesteps in between in situ analysis 
	int morph_timesteps = 0;
	double morph_delta = 0.0;
	double seed_water = 0.0;
	double capillary_number = 0.0;
	double tolerance = 0.01;
	double Ca_previous = 0.f;
	double initial_volume = 0.0;
	double delta_volume = 0.0;
	double delta_volume_target = 0.0;
	
	/* history for morphological algoirthm */
	double KRA_MORPH_FACTOR=0.5;
	double volA_prev = 0.0; 
	double log_krA_prev = 1.0;
	double log_krA_target = 1.0;
	double log_krA = 1.0;
	double slope_krA_volume = 0.0;
	if (color_db->keyExists( "vol_A_previous" )){
		volA_prev  = color_db->getScalar<double>( "vol_A_previous" );
	}
	if (color_db->keyExists( "log_krA_previous" )){
		log_krA_prev  = color_db->getScalar<double>( "log_krA_previous" );
	}
	if (color_db->keyExists( "krA_morph_factor" )){
		KRA_MORPH_FACTOR  = color_db->getScalar<double>( "krA_morph_factor" );
	}
	/* defaults for simulation protocols */
	auto protocol = color_db->getWithDefault<std::string>( "protocol", "none" );
	if (protocol == "image sequence"){
		// Get the list of images
		USE_DIRECT = true;
		ImageList = color_db->getVector<std::string>( "image_sequence");
		IMAGE_INDEX = color_db->getWithDefault<int>( "image_index", 0 );
		IMAGE_COUNT = ImageList.size();
		morph_interval = 10000;
		USE_MORPH = true;
	}
	else if (protocol == "seed water"){
		morph_delta = -0.05;
		seed_water = 0.01;
		USE_SEED = true;
		USE_MORPH = true;
	}
	else if (protocol == "open connected oil"){
		morph_delta = -0.05;
		USE_MORPH = true;
		USE_MORPHOPEN_OIL = true;
	}
	else if (protocol == "shell aggregation"){
		morph_delta = -0.05;
		USE_MORPH = true;
	}  
	if (color_db->keyExists( "capillary_number" )){
		capillary_number = color_db->getScalar<double>( "capillary_number" );
		SET_CAPILLARY_NUMBER=true;
	}
	if (color_db->keyExists( "rescale_force_after_timestep" )){
		RESCALE_FORCE_AFTER_TIMESTEP = color_db->getScalar<int>( "rescale_force_after_timestep" );
		RESCALE_FORCE = true;
	}
	if (color_db->keyExists( "timestep" )){
		timestep = color_db->getScalar<int>( "timestep" );
	}
	if (BoundaryCondition != 0 && BoundaryCondition != 5 && SET_CAPILLARY_NUMBER==true){
		if (rank == 0) printf("WARINING: capillary number target only supported for BC = 0 or 5 \n");
		SET_CAPILLARY_NUMBER=false;
	}
	if (analysis_db->keyExists( "seed_water" )){
		seed_water = analysis_db->getScalar<double>( "seed_water" );
		if (rank == 0) printf("Seed water in oil %f (seed_water) \n",seed_water);
		ASSERT(protocol == "seed water");
	}
	if (analysis_db->keyExists( "morph_delta" )){
		morph_delta = analysis_db->getScalar<double>( "morph_delta" );
		if (rank == 0) printf("Target volume change %f (morph_delta) \n",morph_delta);
	}
	if (analysis_db->keyExists( "morph_interval" )){
		morph_interval = analysis_db->getScalar<int>( "morph_interval" );
		USE_MORPH = true;
	}
	if (analysis_db->keyExists( "use_morphopen_oil" )){
		USE_MORPHOPEN_OIL = analysis_db->getScalar<bool>( "use_morphopen_oil" );
		if (rank == 0 && USE_MORPHOPEN_OIL) printf("Volume change by morphological opening \n");
		USE_MORPH = true;
	}
	if (analysis_db->keyExists( "tolerance" )){
		tolerance = analysis_db->getScalar<double>( "tolerance" );
	}
	if (analysis_db->keyExists( "analysis_interval" )){
		analysis_interval = analysis_db->getScalar<int>( "analysis_interval" );
	}
	if (analysis_db->keyExists( "min_steady_timesteps" )){
		MIN_STEADY_TIMESTEPS = analysis_db->getScalar<int>( "min_steady_timesteps" );
	}
	if (analysis_db->keyExists( "max_steady_timesteps" )){
		MAX_STEADY_TIMESTEPS = analysis_db->getScalar<int>( "max_steady_timesteps" );
	}
	if (analysis_db->keyExists( "max_morph_timesteps" )){
		MAX_MORPH_TIMESTEPS = analysis_db->getScalar<int>( "max_morph_timesteps" );
	}

	if (rank==0){
		printf("********************************************************\n");
		if (protocol == "image sequence"){
			printf("  using protocol = image sequence \n");
			printf("     min_steady_timesteps = %i \n",MIN_STEADY_TIMESTEPS);
			printf("     max_steady_timesteps = %i \n",MAX_STEADY_TIMESTEPS);
			printf("     tolerance = %f \n",tolerance);
			std::string first_image = ImageList[IMAGE_INDEX];
			printf("     first image in sequence: %s ***\n", first_image.c_str());
		}
		else if (protocol == "seed water"){
			printf("  using protocol =  seed water \n");
			printf("     min_steady_timesteps = %i \n",MIN_STEADY_TIMESTEPS);
			printf("     max_steady_timesteps = %i \n",MAX_STEADY_TIMESTEPS);
			printf("     tolerance = %f \n",tolerance);
			printf("     morph_delta = %f \n",morph_delta);
			printf("     seed_water = %f \n",seed_water);
		}
		else if (protocol == "open connected oil"){
			printf("  using protocol = open connected oil \n");
			printf("     min_steady_timesteps = %i \n",MIN_STEADY_TIMESTEPS);
			printf("     max_steady_timesteps = %i \n",MAX_STEADY_TIMESTEPS);
			printf("     tolerance = %f \n",tolerance);
			printf("     morph_delta = %f \n",morph_delta);
		}
		else if (protocol == "shell aggregation"){
			printf("  using protocol = shell aggregation \n"); 
			printf("     min_steady_timesteps = %i \n",MIN_STEADY_TIMESTEPS);
			printf("     max_steady_timesteps = %i \n",MAX_STEADY_TIMESTEPS);
			printf("     tolerance = %f \n",tolerance);
			printf("     morph_delta = %f \n",morph_delta);
		}  
		else if (protocol == "fractional flow"){
			printf("  using protocol = fractional flow \n"); 
			printf("     min_steady_timesteps = %i \n",MIN_STEADY_TIMESTEPS);
			printf("     max_steady_timesteps = %i \n",MAX_STEADY_TIMESTEPS);
			printf("     tolerance = %f \n",tolerance);
			printf("     morph_delta = %f \n",morph_delta);
		} 
		printf("No. of timesteps: %i \n", timestepMax);
		fflush(stdout);
	}

	//************ MAIN ITERATION LOOP ***************************************/
	comm.barrier();
	PROFILE_START("Loop");
    //std::shared_ptr<Database> analysis_db;
	bool Regular = false;
	auto current_db = db->cloneDatabase();
	runAnalysis analysis( current_db, rank_info, ScaLBL_Comm, Dm, Np, Regular, Map );
	//analysis.createThreads( analysis_method, 4 );
    auto t1 = std::chrono::system_clock::now();
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
		ScaLBL_Comm->Barrier();
		ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);

		// Perform the collision operation
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
		if (BoundaryCondition > 0 && BoundaryCondition < 5){
			ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
			ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
		}
		// Halo exchange for phase field
		ScaLBL_Comm_Regular->SendHalo(Phi);

		ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm_Regular->RecvHalo(Phi);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_Comm->Barrier();
		// Set BCs
		if (BoundaryCondition == 3){
			ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		if (BoundaryCondition == 4){
			din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		else if (BoundaryCondition == 5){
			ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
			ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
		}
		ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->Barrier(); 

		// *************EVEN TIMESTEP*************
		timestep++;
		// Compute the Phase indicator field
		ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
		ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
		ScaLBL_Comm->Barrier();
		ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);

		// Perform the collision operation
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
		// Halo exchange for phase field
		if (BoundaryCondition > 0 && BoundaryCondition < 5){
			ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
			ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
		}
		ScaLBL_Comm_Regular->SendHalo(Phi);
		ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz,  Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm_Regular->RecvHalo(Phi);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_Comm->Barrier();
		// Set boundary conditions
		if (BoundaryCondition == 3){
			ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		else if (BoundaryCondition == 4){
			din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
		else if (BoundaryCondition == 5){
			ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
			ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
		}
		ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
				alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->Barrier(); 
		//************************************************************************
		PROFILE_STOP("Update");

		if (rank==0 && timestep%analysis_interval == 0 && BoundaryCondition == 4){
			printf("%i %f \n",timestep,din);
		}
		// Run the analysis
		analysis.basic(timestep, current_db, *Averages, Phi, Pressure, Velocity, fq, Den );

		// allow initial ramp-up to get closer to steady state
		if (timestep > RAMP_TIMESTEPS && timestep%analysis_interval == 0 && USE_MORPH){
			analysis.finish();
			CURRENT_STEADY_TIMESTEPS += analysis_interval;

			double volB = Averages->gwb.V; 
			double volA = Averages->gnb.V; 
			volA /= Dm->Volume;
			volB /= Dm->Volume;;
			//initial_volume = volA*Dm->Volume;
			double vA_x = Averages->gnb.Px/Averages->gnb.M; 
			double vA_y = Averages->gnb.Py/Averages->gnb.M; 
			double vA_z = Averages->gnb.Pz/Averages->gnb.M; 
			double vB_x = Averages->gwb.Px/Averages->gwb.M; 
			double vB_y = Averages->gwb.Py/Averages->gwb.M; 
			double vB_z = Averages->gwb.Pz/Averages->gwb.M;
			double muA = rhoA*(tauA-0.5)/3.f; 
			double muB = rhoB*(tauB-0.5)/3.f;				
			double force_mag = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
			double dir_x = Fx/force_mag;
			double dir_y = Fy/force_mag;
			double dir_z = Fz/force_mag;
			if (force_mag == 0.0){
				// default to z direction
				dir_x = 0.0;
				dir_y = 0.0;
				dir_z = 1.0;
				force_mag = 1.0;
			}
			double current_saturation = volB/(volA+volB);
			double flow_rate_A = volA*(vA_x*dir_x + vA_y*dir_y + vA_z*dir_z);
			double flow_rate_B = volB*(vB_x*dir_x + vB_y*dir_y + vB_z*dir_z);
			double Ca = fabs(muA*flow_rate_A + muB*flow_rate_B)/(5.796*alpha);
			
			if ( morph_timesteps > morph_interval ){
				
				bool isSteady = false;
				if ( (fabs((Ca - Ca_previous)/Ca) < tolerance &&  CURRENT_STEADY_TIMESTEPS > MIN_STEADY_TIMESTEPS))
					isSteady = true;
				if (CURRENT_STEADY_TIMESTEPS > MAX_STEADY_TIMESTEPS)
					isSteady = true;
				if (RESCALE_FORCE == true && SET_CAPILLARY_NUMBER == true && CURRENT_STEADY_TIMESTEPS > RESCALE_FORCE_AFTER_TIMESTEP){
				  RESCALE_FORCE = false;
				  double RESCALE_FORCE_FACTOR = capillary_number / Ca;
				  if (RESCALE_FORCE_FACTOR > 2.0) RESCALE_FORCE_FACTOR = 2.0;
				  if (RESCALE_FORCE_FACTOR < 0.5) RESCALE_FORCE_FACTOR = 0.5;
				  Fx *= RESCALE_FORCE_FACTOR;
				  Fy *= RESCALE_FORCE_FACTOR;
				  Fz *= RESCALE_FORCE_FACTOR;
				  force_mag = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
				  if (force_mag > 1e-3){
				    Fx *= 1e-3/force_mag;   // impose ceiling for stability
				    Fy *= 1e-3/force_mag;   
				    Fz *= 1e-3/force_mag;   
				  }
				  if (rank == 0) printf("    -- adjust force by factor %f \n ",capillary_number / Ca);
				  Averages->SetParams(rhoA,rhoB,tauA,tauB,Fx,Fy,Fz,alpha,beta);
				  color_db->putVector<double>("F",{Fx,Fy,Fz});
				}
				if ( isSteady ){
					MORPH_ADAPT = true;
					CURRENT_MORPH_TIMESTEPS=0;
					delta_volume_target = Dm->Volume*volA *morph_delta; // set target volume change
					//****** ENDPOINT ADAPTATION ********/
					double krA_TMP= fabs(muA*flow_rate_A / force_mag);
					double krB_TMP= fabs(muB*flow_rate_B / force_mag);
					log_krA = log(krA_TMP);
					if (krA_TMP < 0.0){
						// cannot do endpoint adaptation if kr is negative
						log_krA = log_krA_prev;
					}
					else if (krA_TMP < krB_TMP && morph_delta > 0.0){
						/** morphological target based on relative permeability for A **/
						log_krA_target = log(KRA_MORPH_FACTOR*(krA_TMP));
						slope_krA_volume = (log_krA - log_krA_prev)/(Dm->Volume*(volA - volA_prev));
						delta_volume_target=min(delta_volume_target,Dm->Volume*(volA+(log_krA_target - log_krA)/slope_krA_volume));
						if (rank==0){
							printf("    Enabling endpoint adaptation: krA = %f, krB = %f \n",krA_TMP,krB_TMP);	
							printf("    log(kr)=%f, volume=%f, TARGET log(kr)=%f, volume change=%f \n",log_krA, volA, log_krA_target, delta_volume_target/(volA*Dm->Volume));							
						}
					}
					log_krA_prev = log_krA;
					volA_prev = volA;
					//******************************** **/
					/**  compute averages & write data **/
					Averages->Full();
					Averages->Write(timestep);
					analysis.WriteVisData(timestep, current_db, *Averages, Phi, Pressure, Velocity, fq, Den );
					analysis.finish();
					
					if (rank==0){
						printf("** WRITE STEADY POINT *** ");
						printf("Ca = %f, (previous = %f) \n",Ca,Ca_previous);
						double h = Dm->voxel_length;		
						// pressures
						double pA = Averages->gnb.p;
						double pB = Averages->gwb.p;
						double pAc = Averages->gnc.p;
						double pBc = Averages->gwc.p;
						double pAB = (pA-pB)/(h*6.0*alpha);
						double pAB_connected = (pAc-pBc)/(h*6.0*alpha);
						// connected contribution
						double Vol_nc = Averages->gnc.V/Dm->Volume;
						double Vol_wc = Averages->gwc.V/Dm->Volume;
						double Vol_nd = Averages->gnd.V/Dm->Volume;
						double Vol_wd = Averages->gwd.V/Dm->Volume;
						double Mass_n = Averages->gnc.M + Averages->gnd.M;
						double Mass_w = Averages->gwc.M + Averages->gwd.M;
						double vAc_x = Averages->gnc.Px/Mass_n; 
						double vAc_y = Averages->gnc.Py/Mass_n; 
						double vAc_z = Averages->gnc.Pz/Mass_n; 
						double vBc_x = Averages->gwc.Px/Mass_w; 
						double vBc_y = Averages->gwc.Py/Mass_w; 
						double vBc_z = Averages->gwc.Pz/Mass_w;
						// disconnected contribution
						double vAd_x = Averages->gnd.Px/Mass_n; 
						double vAd_y = Averages->gnd.Py/Mass_n; 
						double vAd_z = Averages->gnd.Pz/Mass_n; 
						double vBd_x = Averages->gwd.Px/Mass_w; 
						double vBd_y = Averages->gwd.Py/Mass_w; 
						double vBd_z = Averages->gwd.Pz/Mass_w;
						
						double flow_rate_A_connected = Vol_nc*(vAc_x*dir_x + vAc_y*dir_y + vAc_z*dir_z);
						double flow_rate_B_connected = Vol_wc*(vBc_x*dir_x + vBc_y*dir_y + vBc_z*dir_z);
						double flow_rate_A_disconnected = (Vol_nd)*(vAd_x*dir_x + vAd_y*dir_y + vAd_z*dir_z);
						double flow_rate_B_disconnected = (Vol_wd)*(vBd_x*dir_x + vBd_y*dir_y + vBd_z*dir_z);
						
						double kAeff_connected = h*h*muA*flow_rate_A_connected/(force_mag);
						double kBeff_connected = h*h*muB*flow_rate_B_connected/(force_mag);
						
						double kAeff_disconnected = h*h*muA*flow_rate_A_disconnected/(force_mag);
						double kBeff_disconnected = h*h*muB*flow_rate_B_disconnected/(force_mag);

						double kAeff = h*h*muA*(flow_rate_A)/(force_mag);
						double kBeff = h*h*muB*(flow_rate_B)/(force_mag);

						double viscous_pressure_drop = (rhoA*volA + rhoB*volB)*force_mag;
						double Mobility = muA/muB;
						
						bool WriteHeader=false;
						FILE * kr_log_file = fopen("relperm.csv","r");
						if (kr_log_file != NULL)
							fclose(kr_log_file);
						else
							WriteHeader=true;
						kr_log_file = fopen("relperm.csv","a");
						if (WriteHeader)
							fprintf(kr_log_file,"timesteps sat.water eff.perm.oil eff.perm.water eff.perm.oil.connected eff.perm.water.connected eff.perm.oil.disconnected eff.perm.water.disconnected cap.pressure cap.pressure.connected pressure.drop Ca M\n");

						fprintf(kr_log_file,"%i %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g\n",CURRENT_STEADY_TIMESTEPS,current_saturation,kAeff,kBeff,kAeff_connected,kBeff_connected,kAeff_disconnected,kBeff_disconnected,pAB,pAB_connected,viscous_pressure_drop,Ca,Mobility);
						fclose(kr_log_file);

						printf("  Measured capillary number %f \n ",Ca);
					}
					if (SET_CAPILLARY_NUMBER ){
						Fx *= capillary_number / Ca;
						Fy *= capillary_number / Ca;
						Fz *= capillary_number / Ca;
						if (force_mag > 1e-3){
							Fx *= 1e-3/force_mag;   // impose ceiling for stability
							Fy *= 1e-3/force_mag;   
							Fz *= 1e-3/force_mag;   
						}
						if (rank == 0) printf("    -- adjust force by factor %f \n ",capillary_number / Ca);
						Averages->SetParams(rhoA,rhoB,tauA,tauB,Fx,Fy,Fz,alpha,beta);
						color_db->putVector<double>("F",{Fx,Fy,Fz});
					}
					
					CURRENT_STEADY_TIMESTEPS = 0;
				}
				else{
					if (rank==0){
						printf("** Continue to simulate steady *** \n ");
						printf("Ca = %f, (previous = %f) \n",Ca,Ca_previous);
					}
				}
				morph_timesteps=0;
				Ca_previous = Ca;
			}

			if (MORPH_ADAPT ){
				CURRENT_MORPH_TIMESTEPS += analysis_interval;
				if (USE_DIRECT){
					// Use image sequence
					IMAGE_INDEX++;
					MORPH_ADAPT = false;
					if (IMAGE_INDEX < IMAGE_COUNT){
						std::string next_image = ImageList[IMAGE_INDEX];
						if (rank==0) printf("***Loading next image in sequence (%i) ***\n",IMAGE_INDEX);
						color_db->putScalar<int>("image_index",IMAGE_INDEX);
						ImageInit(next_image);
					}
					else{
						if (rank==0) printf("Finished simulating image sequence \n");
						timestep = timestepMax;
					}
				}
				else if (USE_SEED){
					delta_volume = volA*Dm->Volume - initial_volume;
					CURRENT_MORPH_TIMESTEPS += analysis_interval;
					double massChange = SeedPhaseField(seed_water);
					if (rank==0) printf("***Seed water in oil %f, volume change %f / %f ***\n", massChange, delta_volume, delta_volume_target);
				}
				else if (USE_MORPHOPEN_OIL){
					delta_volume = volA*Dm->Volume - initial_volume;
					if (rank==0) printf("***Morphological opening of connected oil, target volume change %f ***\n", delta_volume_target);
					MorphOpenConnected(delta_volume_target);
				}
				else {
					if (rank==0) printf("***Shell aggregation, target volume change %f ***\n", delta_volume_target);
					//double delta_volume_target = volB - (volA + volB)*TARGET_SATURATION; // change in volume to A
					delta_volume += MorphInit(beta,delta_volume_target-delta_volume);
				}

				if ( (delta_volume - delta_volume_target)/delta_volume_target > 0.0 ){
					MORPH_ADAPT = false;
					CURRENT_STEADY_TIMESTEPS=0;
					initial_volume = volA*Dm->Volume;
					delta_volume = 0.0;
					if (RESCALE_FORCE_AFTER_TIMESTEP > 0)
						RESCALE_FORCE = true;
				}
				else if (!(USE_DIRECT) && CURRENT_MORPH_TIMESTEPS > MAX_MORPH_TIMESTEPS) {
					MORPH_ADAPT = false;
					CURRENT_STEADY_TIMESTEPS=0;
					initial_volume = volA*Dm->Volume;
					delta_volume = 0.0;
					RESCALE_FORCE = true;
					if (RESCALE_FORCE_AFTER_TIMESTEP > 0)
						RESCALE_FORCE = true;
				}
			}
			morph_timesteps += analysis_interval;
		}
		comm.barrier();
	}
	analysis.finish();
	PROFILE_STOP("Loop");
	PROFILE_SAVE("lbpm_color_simulator",1);
	//************************************************************************
	ScaLBL_Comm->Barrier();
	if (rank==0) printf("-------------------------------------------------------------------\n");
	// Compute the walltime per timestep
    auto t2 = std::chrono::system_clock::now();
	double cputime = std::chrono::duration<double>( t2 - t1 ).count() / timestep;
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

double ScaLBL_ColorModel::ImageInit(std::string Filename){
	
	if (rank==0) printf("Re-initializing fluids from file: %s \n", Filename.c_str());
	Mask->Decomp(Filename);
	for (int i=0; i<Nx*Ny*Nz; i++) id[i] = Mask->id[i];  // save what was read
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i];  // save what was read

	double *PhaseLabel;
	PhaseLabel = new double[Nx*Ny*Nz];
	AssignComponentLabels(PhaseLabel);

	double Count = 0.0;
	double PoreCount = 0.0;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				if (id[Nx*Ny*k+Nx*j+i] == 2){
					PoreCount++;
					Count++;
				}
				else if (id[Nx*Ny*k+Nx*j+i] == 1){
					PoreCount++;						
				}
			}
		}
	}

	Count=Dm->Comm.sumReduce(  Count);
	PoreCount=Dm->Comm.sumReduce(  PoreCount);
	
	if (rank==0) printf("   new saturation: %f (%f / %f) \n", Count / PoreCount, Count, PoreCount);
	ScaLBL_CopyToDevice(Phi, PhaseLabel, Nx*Ny*Nz*sizeof(double));
	comm.barrier();
	
	ScaLBL_D3Q19_Init(fq, Np);
	ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np);
	ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
	comm.barrier();
	
	ScaLBL_CopyToHost(Averages->Phi.data(),Phi,Nx*Ny*Nz*sizeof(double));

	double saturation = Count/PoreCount;
	return saturation;

}

double ScaLBL_ColorModel::MorphOpenConnected(double target_volume_change){
	
	int nx = Nx;
	int ny = Ny;
	int nz = Nz;
	int n;
	int N = nx*ny*nz;
	double volume_change=0.0;
	
	if (target_volume_change < 0.0){
		Array<char> id_solid(nx,ny,nz);
		Array<int> phase_label(nx,ny,nz);
		DoubleArray distance(Nx,Ny,Nz);
		DoubleArray phase(nx,ny,nz);
		signed char *id_connected;
		id_connected = new signed char [nx*ny*nz];

		ScaLBL_CopyToHost(phase.data(), Phi, N*sizeof(double));

		// Extract only the connected part of NWP
		double vF=0.0; double vS=0.0;
		ComputeGlobalBlobIDs(nx-2,ny-2,nz-2,Dm->rank_info,phase,Averages->SDs,vF,vS,phase_label,Dm->Comm);
		comm.barrier();

		long long count_connected=0;
		long long count_porespace=0;
		long long count_water=0;
		for (int k=1; k<nz-1; k++){
			for (int j=1; j<ny-1; j++){
				for (int i=1; i<nx-1; i++){
					n=k*nx*ny+j*nx+i;
					// only apply opening to connected component 
					if ( phase_label(i,j,k) == 0){
						count_connected++;
					}
					if (id[n] > 0){
						count_porespace++;
					}
					if (id[n] == 2){
						count_water++;
					}
				}
			}
		}
		count_connected=Dm->Comm.sumReduce(  count_connected);
		count_porespace=Dm->Comm.sumReduce(  count_porespace);
		count_water=Dm->Comm.sumReduce(  count_water);

		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					n=k*nx*ny+j*nx+i;
					// only apply opening to connected component 
					if ( phase_label(i,j,k) == 0){
						id_solid(i,j,k) = 1;
						id_connected[n] = 2;
						id[n] = 2;
						/* delete the connected component */
						phase(i,j,k) = -1.0;
					}
					else{
						id_solid(i,j,k) = 0;
						id_connected[n] = 0;
					}
				}
			}
		}
		CalcDist(distance,id_solid,*Dm);

		signed char water=2;
		signed char notwater=1;
		double SW=-(target_volume_change)/count_connected;
		MorphOpen(distance, id_connected, Dm, SW, water, notwater);
		
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					n=k*nx*ny+j*nx+i;
					// only apply opening to connected component 
					if ( id_connected[n] == 1){
						id_solid(i,j,k) = 0;
					}
					else{
						id_solid(i,j,k) = 1;
					}
				}
			}
		}
		CalcDist(distance,id_solid,*Dm);
		
		// re-initialize
		double beta = 0.95;
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					n=k*nx*ny+j*nx+i;
					double d = distance(i,j,k);
					if (Averages->SDs(i,j,k) > 0.f){
						if (d < 3.f){
							phase(i,j,k) = (2.f*(exp(-2.f*beta*d))/(1.f+exp(-2.f*beta*d))-1.f);	
						}
					}
				} 
			}
		}
		
		int count_morphopen=0.0;
		for (int k=1; k<nz-1; k++){
			for (int j=1; j<ny-1; j++){
				for (int i=1; i<nx-1; i++){
					n=k*nx*ny+j*nx+i;
					// only apply opening to connected component 
					if ( id_connected[n] == 1){
						count_morphopen++;
					}
				}
			}
		}
		count_morphopen=Dm->Comm.sumReduce(  count_morphopen);
		volume_change = double(count_morphopen - count_connected);
		
		if (rank==0)  printf("   opening of connected oil %f \n",volume_change/count_connected);

		ScaLBL_CopyToDevice(Phi,phase.data(),N*sizeof(double));
		ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		if (BoundaryCondition == 1 || BoundaryCondition == 2 || BoundaryCondition == 3 || BoundaryCondition == 4){
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
	return(volume_change);
}

double ScaLBL_ColorModel::SeedPhaseField(const double seed_water_in_oil){
  srand(time(NULL));
  double mass_loss =0.f;
  double count =0.f;
  double *Aq_tmp, *Bq_tmp;
  
  Aq_tmp = new double [7*Np];
  Bq_tmp = new double [7*Np];

  ScaLBL_CopyToHost(Aq_tmp, Aq, 7*Np*sizeof(double));
  ScaLBL_CopyToHost(Bq_tmp, Bq, 7*Np*sizeof(double));
  
 
   for (int n=0; n < ScaLBL_Comm->LastExterior(); n++){
    double random_value = seed_water_in_oil*double(rand())/ RAND_MAX;
    double dA = Aq_tmp[n] + Aq_tmp[n+Np]  + Aq_tmp[n+2*Np] + Aq_tmp[n+3*Np] + Aq_tmp[n+4*Np] + Aq_tmp[n+5*Np] + Aq_tmp[n+6*Np];
    double dB = Bq_tmp[n] + Bq_tmp[n+Np]  + Bq_tmp[n+2*Np] + Bq_tmp[n+3*Np] + Bq_tmp[n+4*Np] + Bq_tmp[n+5*Np] + Bq_tmp[n+6*Np];
    double phase_id = (dA - dB) / (dA + dB);
    if (phase_id > 0.0){
      Aq_tmp[n] -= 0.3333333333333333*random_value;
      Aq_tmp[n+Np] -= 0.1111111111111111*random_value;
      Aq_tmp[n+2*Np] -= 0.1111111111111111*random_value;
      Aq_tmp[n+3*Np] -= 0.1111111111111111*random_value;
      Aq_tmp[n+4*Np] -= 0.1111111111111111*random_value;
      Aq_tmp[n+5*Np] -= 0.1111111111111111*random_value;
      Aq_tmp[n+6*Np] -= 0.1111111111111111*random_value;
      
      Bq_tmp[n] += 0.3333333333333333*random_value;
      Bq_tmp[n+Np] += 0.1111111111111111*random_value;
      Bq_tmp[n+2*Np] += 0.1111111111111111*random_value;
      Bq_tmp[n+3*Np] += 0.1111111111111111*random_value;
      Bq_tmp[n+4*Np] += 0.1111111111111111*random_value;
      Bq_tmp[n+5*Np] += 0.1111111111111111*random_value;
      Bq_tmp[n+6*Np] += 0.1111111111111111*random_value;
    }
    mass_loss += random_value*seed_water_in_oil;
  }

  for (int n=ScaLBL_Comm->FirstInterior(); n < ScaLBL_Comm->LastInterior(); n++){
    double random_value = seed_water_in_oil*double(rand())/ RAND_MAX;
    double dA = Aq_tmp[n] + Aq_tmp[n+Np]  + Aq_tmp[n+2*Np] + Aq_tmp[n+3*Np] + Aq_tmp[n+4*Np] + Aq_tmp[n+5*Np] + Aq_tmp[n+6*Np];
    double dB = Bq_tmp[n] + Bq_tmp[n+Np]  + Bq_tmp[n+2*Np] + Bq_tmp[n+3*Np] + Bq_tmp[n+4*Np] + Bq_tmp[n+5*Np] + Bq_tmp[n+6*Np];
    double phase_id = (dA - dB) / (dA + dB);
    if (phase_id > 0.0){
      Aq_tmp[n] -= 0.3333333333333333*random_value;
      Aq_tmp[n+Np] -= 0.1111111111111111*random_value;
      Aq_tmp[n+2*Np] -= 0.1111111111111111*random_value;
      Aq_tmp[n+3*Np] -= 0.1111111111111111*random_value;
      Aq_tmp[n+4*Np] -= 0.1111111111111111*random_value;
      Aq_tmp[n+5*Np] -= 0.1111111111111111*random_value;
      Aq_tmp[n+6*Np] -= 0.1111111111111111*random_value;
      
      Bq_tmp[n] += 0.3333333333333333*random_value;
      Bq_tmp[n+Np] += 0.1111111111111111*random_value;
      Bq_tmp[n+2*Np] += 0.1111111111111111*random_value;
      Bq_tmp[n+3*Np] += 0.1111111111111111*random_value;
      Bq_tmp[n+4*Np] += 0.1111111111111111*random_value;
      Bq_tmp[n+5*Np] += 0.1111111111111111*random_value;
      Bq_tmp[n+6*Np] += 0.1111111111111111*random_value;
    }
    mass_loss += random_value*seed_water_in_oil;
  }

  count= Dm->Comm.sumReduce(  count);
  mass_loss= Dm->Comm.sumReduce(  mass_loss);
  if (rank == 0) printf("Remove mass %f from %f voxels \n",mass_loss,count);

  // Need to initialize Aq, Bq, Den, Phi directly
  //ScaLBL_CopyToDevice(Phi,phase.data(),7*Np*sizeof(double));
  ScaLBL_CopyToDevice(Aq, Aq_tmp, 7*Np*sizeof(double));
  ScaLBL_CopyToDevice(Bq, Bq_tmp, 7*Np*sizeof(double));

  return(mass_loss);
}

double ScaLBL_ColorModel::MorphInit(const double beta, const double target_delta_volume){
	const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);

	double vF = 0.f;
	double vS = 0.f;
	double delta_volume;
	double WallFactor = 1.0;
	bool USE_CONNECTED_NWP = false;

	DoubleArray phase(Nx,Ny,Nz);
	IntArray phase_label(Nx,Ny,Nz);;
	DoubleArray phase_distance(Nx,Ny,Nz);
	Array<char> phase_id(Nx,Ny,Nz);
	fillHalo<double> fillDouble(Dm->Comm,Dm->rank_info,{Nx-2,Ny-2,Nz-2},{1,1,1},0,1);
	

	// Basic algorithm to 
	// 1. Copy phase field to CPU
	ScaLBL_CopyToHost(phase.data(), Phi, N*sizeof(double));

	double count = 0.f;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				if (phase(i,j,k) > 0.f && Averages->SDs(i,j,k) > 0.f) count+=1.f;
			}
		}
	}
	double volume_initial = Dm->Comm.sumReduce(  count);
	double PoreVolume = Dm->Volume*Dm->Porosity();
	/*ensure target isn't an absurdly small fraction of pore volume */
	if (volume_initial < target_delta_volume*PoreVolume){
		volume_initial = target_delta_volume*PoreVolume;
	}
	/*
	sprintf(LocalRankFilename,"phi_initial.%05i.raw",rank);
	FILE *INPUT = fopen(LocalRankFilename,"wb");
	fwrite(phase.data(),8,N,INPUT);
	fclose(INPUT);
	*/
	// 2. Identify connected components of phase field -> phase_label
	
	double volume_connected = 0.0;
       	double second_biggest = 0.0;
	if (USE_CONNECTED_NWP){
		ComputeGlobalBlobIDs(Nx-2,Ny-2,Nz-2,rank_info,phase,Averages->SDs,vF,vS,phase_label,comm);
		comm.barrier();

		// only operate on component "0"
		count = 0.0;

		for (int k=0; k<Nz; k++){
			for (int j=0; j<Ny; j++){
				for (int i=0; i<Nx; i++){
					int label = phase_label(i,j,k);
					if (label == 0 ){
						phase_id(i,j,k) = 0;
						count += 1.0;
					}
					else 		
						phase_id(i,j,k) = 1;
					if (label == 1 ){
						second_biggest += 1.0;
					}
				}
			}
		}	
		volume_connected = Dm->Comm.sumReduce(  count);
		second_biggest = Dm->Comm.sumReduce(  second_biggest);
	}
	else {
		// use the whole NWP 
		for (int k=0; k<Nz; k++){
			for (int j=0; j<Ny; j++){
				for (int i=0; i<Nx; i++){
					if (Averages->SDs(i,j,k) > 0.f){
						if (phase(i,j,k) > 0.f ){
							phase_id(i,j,k) = 0;
						}
						else {
							phase_id(i,j,k) = 1;
						}
					}
					else {
						phase_id(i,j,k) = 1;
					}
				}
			}
		}
	}

	/*int reach_x, reach_y, reach_z;
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
			}
		}
	}*/

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
					// erase the original object
					phase(i,j,k) = -1.0;
				}
			}
		}
	}

	if (USE_CONNECTED_NWP){
	if (volume_connected - second_biggest < 2.0*fabs(target_delta_volume) && target_delta_volume < 0.0){
		// if connected volume is less than 2% just delete the whole thing
		if (rank==0) printf("Connected region has shrunk! \n");
		REVERSE_FLOW_DIRECTION = true;
	}
	
/*	else{*/
		if (rank==0) printf("Pathway volume / next largest ganglion %f \n",volume_connected/second_biggest );
	}
		if (rank==0) printf("MorphGrow with target volume fraction change %f \n", target_delta_volume/volume_initial);
		double target_delta_volume_incremental = target_delta_volume;
		if (fabs(target_delta_volume) > 0.01*volume_initial)  
			target_delta_volume_incremental = 0.01*volume_initial*target_delta_volume/fabs(target_delta_volume);
		delta_volume = MorphGrow(Averages->SDs,phase_distance,phase_id,Averages->Dm, target_delta_volume_incremental, WallFactor);

		for (int k=0; k<Nz; k++){
			for (int j=0; j<Ny; j++){
				for (int i=0; i<Nx; i++){
					if (phase_distance(i,j,k) < 0.0 ) phase_id(i,j,k) = 0;
					else 		     				  phase_id(i,j,k) = 1;
					//if (phase_distance(i,j,k) < 0.0 ) phase(i,j,k) = 1.0;
				}
			}
		}	

		CalcDist(phase_distance,phase_id,*Dm); // re-calculate distance

		// 5. Update phase indicator field based on new distnace
		for (int k=0; k<Nz; k++){
			for (int j=0; j<Ny; j++){
				for (int i=0; i<Nx; i++){
					double d = phase_distance(i,j,k);
					if (Averages->SDs(i,j,k) > 0.f){
						if (d < 3.f){
							//phase(i,j,k) = -1.0;
							phase(i,j,k) = (2.f*(exp(-2.f*beta*d))/(1.f+exp(-2.f*beta*d))-1.f);	
						}
					}
				} 
			}
		}
		fillDouble.fill(phase);
	//}

	count = 0.f;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				if (phase(i,j,k) > 0.f && Averages->SDs(i,j,k) > 0.f){
					count+=1.f;
				}
			}
		}
	}
	double volume_final= Dm->Comm.sumReduce(  count);

	delta_volume = (volume_final-volume_initial);
	if (rank == 0)  printf("MorphInit: change fluid volume fraction by %f \n", delta_volume/volume_initial);
	if (rank == 0)  printf("   new saturation =  %f \n", volume_final/(Mask->Porosity()*double((Nx-2)*(Ny-2)*(Nz-2)*nprocs)));

	// 6. copy back to the device
	//if (rank==0)  printf("MorphInit: copy data  back to device\n");
	ScaLBL_CopyToDevice(Phi,phase.data(),N*sizeof(double));
	/*
	sprintf(LocalRankFilename,"dist_final.%05i.raw",rank);
	FILE *DIST = fopen(LocalRankFilename,"wb");
	fwrite(phase_distance.data(),8,N,DIST);
	fclose(DIST);
	
	sprintf(LocalRankFilename,"phi_final.%05i.raw",rank);
	FILE *PHI = fopen(LocalRankFilename,"wb");
	fwrite(phase.data(),8,N,PHI);
	fclose(PHI);
	*/
	// 7. Re-initialize phase field and density
	ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np);
	ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
	if (BoundaryCondition == 1 || BoundaryCondition == 2 || BoundaryCondition == 3 || BoundaryCondition == 4){
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

    ScaLBL_Comm->RegularLayout(Map,&Den[0],PhaseField);
	FILE *AFILE;
	sprintf(LocalRankFilename,"A.%05i.raw",rank);
	AFILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,AFILE);
	fclose(AFILE);

	ScaLBL_Comm->RegularLayout(Map,&Den[Np],PhaseField);
	FILE *BFILE;
	sprintf(LocalRankFilename,"B.%05i.raw",rank);
	BFILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,BFILE);
	fclose(BFILE);

	ScaLBL_Comm->RegularLayout(Map,Pressure,PhaseField);
	FILE *PFILE;
	sprintf(LocalRankFilename,"Pressure.%05i.raw",rank);
	PFILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,PFILE);
	fclose(PFILE);

	ScaLBL_Comm->RegularLayout(Map,&Velocity[0],PhaseField);
	FILE *VELX_FILE;
	sprintf(LocalRankFilename,"Velocity_X.%05i.raw",rank);
	VELX_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,VELX_FILE);
	fclose(VELX_FILE);

	ScaLBL_Comm->RegularLayout(Map,&Velocity[Np],PhaseField);
	FILE *VELY_FILE;
	sprintf(LocalRankFilename,"Velocity_Y.%05i.raw",rank);
	VELY_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,VELY_FILE);
	fclose(VELY_FILE);

	ScaLBL_Comm->RegularLayout(Map,&Velocity[2*Np],PhaseField);
	FILE *VELZ_FILE;
	sprintf(LocalRankFilename,"Velocity_Z.%05i.raw",rank);
	VELZ_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,VELZ_FILE);
	fclose(VELZ_FILE);

/*	ScaLBL_Comm->RegularLayout(Map,&ColorGrad[0],PhaseField);
	FILE *CGX_FILE;
	sprintf(LocalRankFilename,"Gradient_X.%05i.raw",rank);
	CGX_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,CGX_FILE);
	fclose(CGX_FILE);

	ScaLBL_Comm->RegularLayout(Map,&ColorGrad[Np],PhaseField);
	FILE *CGY_FILE;
	sprintf(LocalRankFilename,"Gradient_Y.%05i.raw",rank);
	CGY_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,CGY_FILE);
	fclose(CGY_FILE);

	ScaLBL_Comm->RegularLayout(Map,&ColorGrad[2*Np],PhaseField);
	FILE *CGZ_FILE;
	sprintf(LocalRankFilename,"Gradient_Z.%05i.raw",rank);
	CGZ_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,CGZ_FILE);
	fclose(CGZ_FILE);
*/
}

FlowAdaptor::FlowAdaptor(ScaLBL_ColorModel &M){
	Nx = M.Dm->Nx;
	Ny = M.Dm->Ny;
	Nz = M.Dm->Nz;
	timestep=-1;
	timestep_previous=-1;	
	
	phi.resize(Nx,Ny,Nz);         phi.fill(0);	    // phase indicator field
	phi_t.resize(Nx,Ny,Nz);       phi_t.fill(0);	// time derivative for the phase indicator field
}

FlowAdaptor::~FlowAdaptor(){
	
}

double FlowAdaptor::UpdateFractionalFlow(ScaLBL_ColorModel &M){	
	
  	  double MASS_FRACTION_CHANGE = 0.05;
	  if (M.db->keyExists( "FlowAdaptor" )){
		  auto flow_db = M.db->getDatabase( "FlowAdaptor" );
		  MASS_FRACTION_CHANGE = flow_db->getWithDefault<double>( "fractional_flow_increment", 0.05);
	  }

	  int Np = M.Np;
	  double dA, dB, phi;
	  double vx,vy,vz;
	  double vax,vay,vaz;
	  double vbx,vby,vbz;
	  double vax_global,vay_global,vaz_global;
	  double vbx_global,vby_global,vbz_global;

	  double mass_a, mass_b, mass_a_global, mass_b_global;
	  
	  double *Aq_tmp, *Bq_tmp;
	  double *Vel_x, *Vel_y, *Vel_z, *Phase;
	  
	  Aq_tmp = new double [7*Np];
	  Bq_tmp = new double [7*Np];
	  Phase = new double [Np];
	  Vel_x = new double [Np];
	  Vel_y = new double [Np];
	  Vel_z = new double [Np];

	  ScaLBL_CopyToHost(Aq_tmp, M.Aq, 7*Np*sizeof(double));
	  ScaLBL_CopyToHost(Bq_tmp, M.Bq, 7*Np*sizeof(double));
	  ScaLBL_CopyToHost(Vel_x, &M.Velocity[0], Np*sizeof(double));
	  ScaLBL_CopyToHost(Vel_y, &M.Velocity[Np], Np*sizeof(double));
	  ScaLBL_CopyToHost(Vel_z, &M.Velocity[2*Np], Np*sizeof(double));
	  
	  /* DEBUG STRUCTURES */
	  int Nx = M.Nx;  int Ny = M.Ny;  int Nz = M.Nz;
	  int N = Nx*Ny*Nz;
	  double * DebugMassA, *DebugMassB;
	  DebugMassA = new double[Np];
	  DebugMassB = new double[Np];
	  
	  /* compute the total momentum */
	  vax = vay = vaz = 0.0;
	  vbx = vby = vbz = 0.0;
	  mass_a = mass_b = 0.0;
	  for (int n=0; n < M.ScaLBL_Comm->LastExterior(); n++){
		  dA = Aq_tmp[n] + Aq_tmp[n+Np]  + Aq_tmp[n+2*Np] + Aq_tmp[n+3*Np] + Aq_tmp[n+4*Np] + Aq_tmp[n+5*Np] + Aq_tmp[n+6*Np];
		  dB = Bq_tmp[n] + Bq_tmp[n+Np]  + Bq_tmp[n+2*Np] + Bq_tmp[n+3*Np] + Bq_tmp[n+4*Np] + Bq_tmp[n+5*Np] + Bq_tmp[n+6*Np];
		  phi = (dA - dB) / (dA + dB);
		  Phase[n] = phi;
		  mass_a += dA;
		  mass_b += dB;
		  if (phi > 0.0){
			  vax += Vel_x[n];
			  vay += Vel_y[n];
			  vaz += Vel_z[n];
		  }
		  else {
			  vbx += Vel_x[n];
			  vby += Vel_y[n];
			  vbz += Vel_z[n];
		  }
	  }
	  for (int n=M.ScaLBL_Comm->FirstInterior(); n < M.ScaLBL_Comm->LastInterior(); n++){
		    dA = Aq_tmp[n] + Aq_tmp[n+Np]  + Aq_tmp[n+2*Np] + Aq_tmp[n+3*Np] + Aq_tmp[n+4*Np] + Aq_tmp[n+5*Np] + Aq_tmp[n+6*Np];
		    dB = Bq_tmp[n] + Bq_tmp[n+Np]  + Bq_tmp[n+2*Np] + Bq_tmp[n+3*Np] + Bq_tmp[n+4*Np] + Bq_tmp[n+5*Np] + Bq_tmp[n+6*Np];
		    phi = (dA - dB) / (dA + dB);
		    Phase[n] = phi;
			  mass_a += dA;
			  mass_b += dB;
			  if (phi > 0.0){
				  vax += Vel_x[n];
				  vay += Vel_y[n];
				  vaz += Vel_z[n];
			  }
			  else {
				  vbx += Vel_x[n];
				  vby += Vel_y[n];
				  vbz += Vel_z[n];
			  }
	   }
	   mass_a_global = M.Dm->Comm.sumReduce(mass_a);
	   mass_b_global = M.Dm->Comm.sumReduce(mass_b);
	   vax_global = M.Dm->Comm.sumReduce(vax);
	   vay_global = M.Dm->Comm.sumReduce(vay);
	   vaz_global = M.Dm->Comm.sumReduce(vaz);
	   vbx_global = M.Dm->Comm.sumReduce(vbx);
	   vby_global = M.Dm->Comm.sumReduce(vby);
	   vbz_global = M.Dm->Comm.sumReduce(vbz);	 
	   
	   double  total_momentum_A = sqrt(vax_global*vax_global+vay_global*vay_global+vaz_global*vaz_global);
	   double  total_momentum_B = sqrt(vbx_global*vbx_global+vby_global*vby_global+vbz_global*vbz_global);
	   double total_momentum = total_momentum_A + total_momentum_B;
	   /* compute the total mass change */
	   double TOTAL_MASS_CHANGE = MASS_FRACTION_CHANGE*(mass_a_global + mass_b_global);
	   if (fabs(TOTAL_MASS_CHANGE) > 0.1*mass_a_global )
		   TOTAL_MASS_CHANGE = 0.1*mass_a_global;
	   if (fabs(TOTAL_MASS_CHANGE) > 0.1*mass_b_global )
		   TOTAL_MASS_CHANGE = 0.1*mass_b_global;
	   
	   double LOCAL_MASS_CHANGE = 0.0;
	   for (int n=0; n < M.ScaLBL_Comm->LastExterior(); n++){	    
		   phi = Phase[n];
		   vx = Vel_x[n];
		   vy = Vel_y[n];
		   vz = Vel_z[n];
		   double local_momentum = sqrt(vx*vx+vy*vy+vz*vz);
		   if (phi > 0.0){
			   LOCAL_MASS_CHANGE = TOTAL_MASS_CHANGE*local_momentum/total_momentum_A;
			   Aq_tmp[n] -= 0.3333333333333333*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+2*Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+3*Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+4*Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+5*Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+6*Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   DebugMassA[n] = (-1.0)*LOCAL_MASS_CHANGE;
		   }
		   else{
			   LOCAL_MASS_CHANGE = TOTAL_MASS_CHANGE*local_momentum/total_momentum_B;
			   Bq_tmp[n] += 0.3333333333333333*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+2*Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+3*Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+4*Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+5*Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+6*Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   DebugMassB[n] = LOCAL_MASS_CHANGE;
		   }
	   }

	  for (int n=M.ScaLBL_Comm->FirstInterior(); n < M.ScaLBL_Comm->LastInterior(); n++){
		   phi = Phase[n];
		   vx = Vel_x[n];
		   vy = Vel_y[n];
		   vz = Vel_z[n];
		   double local_momentum = sqrt(vx*vx+vy*vy+vz*vz);
		   if (phi > 0.0){
			   LOCAL_MASS_CHANGE = TOTAL_MASS_CHANGE*local_momentum/total_momentum_A;
			   Aq_tmp[n] -= 0.3333333333333333*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+2*Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+3*Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+4*Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+5*Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Aq_tmp[n+6*Np] -= 0.1111111111111111*LOCAL_MASS_CHANGE;
			   DebugMassA[n] = (-1.0)*LOCAL_MASS_CHANGE;
		   }
		   else{
			   LOCAL_MASS_CHANGE = TOTAL_MASS_CHANGE*local_momentum/total_momentum_B;
			   Bq_tmp[n] += 0.3333333333333333*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+2*Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+3*Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+4*Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+5*Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   Bq_tmp[n+6*Np] += 0.1111111111111111*LOCAL_MASS_CHANGE;
			   DebugMassB[n] = LOCAL_MASS_CHANGE;
		   }
	  }

	  if (M.rank == 0) printf("Update Fractional Flow: change mass of fluid B by %f \n",TOTAL_MASS_CHANGE/mass_b_global);

	  /* Print out debugging info with mass update  */
	  // initialize the array
	  double value;
	  char LocalRankFilename[40];
	  DoubleArray regdata(Nx,Ny,Nz);
	  regdata.fill(0.f);
	  for (int k=0; k<Nz; k++){
		  for (int j=0; j<Ny; j++){
			  for (int i=0; i<Nx; i++){
				  int idx=M.Map(i,j,k);
				  if (!(idx<0)){
					  value=DebugMassA[idx];
					  regdata(i,j,k)=value;
				  }
			  }
		  }
	  }
	  
	  FILE *AFILE;
	  sprintf(LocalRankFilename,"dA.%05i.raw",M.rank);
	  AFILE = fopen(LocalRankFilename,"wb");
	  fwrite(regdata.data(),8,N,AFILE);
	  fclose(AFILE);

	  regdata.fill(0.f);
	  for (int k=0; k<Nz; k++){
		  for (int j=0; j<Ny; j++){
			  for (int i=0; i<Nx; i++){
				  int idx=M.Map(i,j,k);
				  if (!(idx<0)){
					  value=DebugMassB[idx];
					  regdata(i,j,k)=value;
				  }
			  }
		  }
	  }
	  FILE *BFILE;
	  sprintf(LocalRankFilename,"dB.%05i.raw",M.rank);
	  BFILE = fopen(LocalRankFilename,"wb");
	  fwrite(regdata.data(),8,N,BFILE);
	  fclose(BFILE);

	  // Need to initialize Aq, Bq, Den, Phi directly
	  //ScaLBL_CopyToDevice(Phi,phase.data(),7*Np*sizeof(double));
	  ScaLBL_CopyToDevice(M.Aq, Aq_tmp, 7*Np*sizeof(double));
	  ScaLBL_CopyToDevice(M.Bq, Bq_tmp, 7*Np*sizeof(double));

	  return(TOTAL_MASS_CHANGE);
}

void FlowAdaptor::Flatten(ScaLBL_ColorModel &M){	
	
	  int Np = M.Np;
	  double dA, dB;
	  
	  double *Aq_tmp, *Bq_tmp;
	  
	  Aq_tmp = new double [7*Np];
	  Bq_tmp = new double [7*Np];

	  ScaLBL_CopyToHost(Aq_tmp, M.Aq, 7*Np*sizeof(double));
	  ScaLBL_CopyToHost(Bq_tmp, M.Bq, 7*Np*sizeof(double));

	  for (int n=0; n < M.ScaLBL_Comm->LastExterior(); n++){
		  dA = Aq_tmp[n] + Aq_tmp[n+Np]  + Aq_tmp[n+2*Np] + Aq_tmp[n+3*Np] + Aq_tmp[n+4*Np] + Aq_tmp[n+5*Np] + Aq_tmp[n+6*Np];
		  dB = Bq_tmp[n] + Bq_tmp[n+Np]  + Bq_tmp[n+2*Np] + Bq_tmp[n+3*Np] + Bq_tmp[n+4*Np] + Bq_tmp[n+5*Np] + Bq_tmp[n+6*Np];
		  if (dA > 1.0){
			  double mass_change = dA - 1.0;
			  Aq_tmp[n] -= 0.333333333333333*mass_change;
			  Aq_tmp[n+Np] -= 0.111111111111111*mass_change;
			  Aq_tmp[n+2*Np] -= 0.111111111111111*mass_change;
			  Aq_tmp[n+3*Np] -= 0.111111111111111*mass_change;
			  Aq_tmp[n+4*Np] -= 0.111111111111111*mass_change;
			  Aq_tmp[n+5*Np] -= 0.111111111111111*mass_change;
			  Aq_tmp[n+6*Np] -= 0.111111111111111*mass_change;
		  }
		  if (dB > 1.0){
			  double mass_change = dB - 1.0;
			  Bq_tmp[n] -= 0.333333333333333*mass_change;
			  Bq_tmp[n+Np] -= 0.111111111111111*mass_change;
			  Bq_tmp[n+2*Np] -= 0.111111111111111*mass_change;
			  Bq_tmp[n+3*Np] -= 0.111111111111111*mass_change;
			  Bq_tmp[n+4*Np] -= 0.111111111111111*mass_change;
			  Bq_tmp[n+5*Np] -= 0.111111111111111*mass_change;
			  Bq_tmp[n+6*Np] -= 0.111111111111111*mass_change;
		  }
	  }
	  for (int n=M.ScaLBL_Comm->FirstInterior(); n < M.ScaLBL_Comm->LastInterior(); n++){
		  dA = Aq_tmp[n] + Aq_tmp[n+Np]  + Aq_tmp[n+2*Np] + Aq_tmp[n+3*Np] + Aq_tmp[n+4*Np] + Aq_tmp[n+5*Np] + Aq_tmp[n+6*Np];
		  dB = Bq_tmp[n] + Bq_tmp[n+Np]  + Bq_tmp[n+2*Np] + Bq_tmp[n+3*Np] + Bq_tmp[n+4*Np] + Bq_tmp[n+5*Np] + Bq_tmp[n+6*Np];
		  if (dA > 1.0){
			  double mass_change = dA - 1.0;
			  Aq_tmp[n] -= 0.333333333333333*mass_change;
			  Aq_tmp[n+Np] -= 0.111111111111111*mass_change;
			  Aq_tmp[n+2*Np] -= 0.111111111111111*mass_change;
			  Aq_tmp[n+3*Np] -= 0.111111111111111*mass_change;
			  Aq_tmp[n+4*Np] -= 0.111111111111111*mass_change;
			  Aq_tmp[n+5*Np] -= 0.111111111111111*mass_change;
			  Aq_tmp[n+6*Np] -= 0.111111111111111*mass_change;
		  }
		  if (dB > 1.0){
			  double mass_change = dB - 1.0;
			  Bq_tmp[n] -= 0.333333333333333*mass_change;
			  Bq_tmp[n+Np] -= 0.111111111111111*mass_change;
			  Bq_tmp[n+2*Np] -= 0.111111111111111*mass_change;
			  Bq_tmp[n+3*Np] -= 0.111111111111111*mass_change;
			  Bq_tmp[n+4*Np] -= 0.111111111111111*mass_change;
			  Bq_tmp[n+5*Np] -= 0.111111111111111*mass_change;
			  Bq_tmp[n+6*Np] -= 0.111111111111111*mass_change;
		  }
	   }

	  ScaLBL_CopyToDevice(M.Aq, Aq_tmp, 7*Np*sizeof(double));
	  ScaLBL_CopyToDevice(M.Bq, Bq_tmp, 7*Np*sizeof(double));
}

double FlowAdaptor::MoveInterface(ScaLBL_ColorModel &M){	
	
    double INTERFACE_CUTOFF = M.color_db->getWithDefault<double>( "move_interface_cutoff", 0.1 );
    double MOVE_INTERFACE_FACTOR = M.color_db->getWithDefault<double>( "move_interface_factor", 10.0 );

    ScaLBL_CopyToHost( phi.data(), M.Phi, Nx*Ny*Nz* sizeof( double ) );
    /* compute the local derivative of phase indicator field */
    double beta = M.beta;
    double factor = 0.5/beta;
    double total_interface_displacement = 0.0;
    double total_interface_sites = 0.0;
    for (int n=0; n<Nx*Ny*Nz; n++){
    	/* compute the distance to the interface */
    	double value1 = M.Averages->Phi(n);
    	double dist1 = factor*log((1.0+value1)/(1.0-value1));
    	double value2 = phi(n);
    	double dist2 = factor*log((1.0+value2)/(1.0-value2));
    	phi_t(n) = value2;
    	if (value1 < INTERFACE_CUTOFF && value1 > -1*INTERFACE_CUTOFF && value2 < INTERFACE_CUTOFF && value2 > -1*INTERFACE_CUTOFF ){
    		/* time derivative of distance */
    		double dxdt = 0.125*(dist2-dist1);
    		/* extrapolate to move the distance further */
    		double dist3 = dist2 + MOVE_INTERFACE_FACTOR*dxdt;
    		/* compute the new phase interface */
    		phi_t(n) = (2.f*(exp(-2.f*beta*(dist3)))/(1.f+exp(-2.f*beta*(dist3))) - 1.f);
    		total_interface_displacement += fabs(MOVE_INTERFACE_FACTOR*dxdt);
    		total_interface_sites += 1.0;
    	}
	}
    ScaLBL_CopyToDevice( M.Phi, phi_t.data(), Nx*Ny*Nz* sizeof( double ) );	
    return total_interface_sites;
}

