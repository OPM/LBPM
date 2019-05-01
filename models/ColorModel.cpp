/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
color lattice boltzmann model
 */
#include "models/ColorModel.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"
#include <stdlib.h>
#include <time.h>

ScaLBL_ColorModel::ScaLBL_ColorModel(int RANK, int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tauA(0),tauB(0),rhoA(0),rhoB(0),alpha(0),beta(0),
Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),inletA(0),inletB(0),outletA(0),outletB(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{
	REVERSE_FLOW_DIRECTION = false;
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

	BoundaryCondition = 0;
	if (domain_db->keyExists( "BC" )){
		BoundaryCondition = domain_db->getScalar<int>( "BC" );
	}
	if (BoundaryCondition==4) flux *= rhoA; // mass flux must adjust for density (see formulation for details)

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
	MPI_Barrier(comm);
	Dm->CommInit();
	MPI_Barrier(comm);
	// Read domain parameters
	rank = Dm->rank();	
	nprocx = Dm->nprocx();
	nprocy = Dm->nprocy();
	nprocz = Dm->nprocz();
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
	
	Averages->SetParams(rhoA,rhoB,tauA,tauB,Fx,Fy,Fz,alpha,beta);
}

void ScaLBL_ColorModel::AssignComponentLabels(double *phase)
{
	size_t NLABELS=0;
	signed char VALUE=0;
	double AFFINITY=0.f;

	auto LabelList = color_db->getVector<int>( "ComponentLabels" );
	auto AffinityList = color_db->getVector<double>( "ComponentAffinity" );

	NLABELS=LabelList.size();
	if (NLABELS != AffinityList.size()){
		ERROR("Error: ComponentLabels and ComponentAffinity must be the same length! \n");
	}

	double label_count[NLABELS];
	double label_count_global[NLABELS];
	// Assign the labels

	for (int idx=0; idx<NLABELS; idx++) label_count[idx]=0;

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

	for (int idx=0; idx<NLABELS; idx++)		label_count_global[idx]=sumReduce( Dm->Comm, label_count[idx]);

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
	ScaLBL_CopyToHost(Averages->Phi.data(),Phi,N*sizeof(double));
}

void ScaLBL_ColorModel::Run(){
	int nprocs=nprocx*nprocy*nprocz;
	const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
	
	bool SET_CAPILLARY_NUMBER = false;
	bool MORPH_ADAPT = false;
	bool USE_MORPH = false;
	bool USE_SEED = false;
	int analysis_interval = 1000; 	// number of timesteps in between in situ analysis 
	int MAX_MORPH_TIMESTEPS = 50000; // maximum number of LBM timesteps to spend in morphological adaptation routine
	int MIN_STEADY_TIMESTEPS = 100000;
	int MAX_STEADY_TIMESTEPS = 200000;
	int RAMP_TIMESTEPS = 0;//50000;		 // number of timesteps to run initially (to get a reasonable velocity field before other pieces kick in)
	int morph_interval = 1000000;
	int CURRENT_MORPH_TIMESTEPS=0;   // counter for number of timesteps spent in  morphological adaptation routine (reset each time)
	int CURRENT_STEADY_TIMESTEPS=0;   // counter for number of timesteps spent in  morphological adaptation routine (reset each time)
	int morph_timesteps = 0;
	double morph_delta = 0.0;
	double seed_water = 0.0;
	double capillary_number = 0.0;
	double tolerance = 0.01;
	double Ca_previous = 0.f;
	double initial_volume = 0.0;
	double delta_volume = 0.0;
	double delta_volume_target = 0.0;
	double RESIDUAL_ENDPOINT_THRESHOLD = 0.04;
	
	if (color_db->keyExists( "residual_endpoint_threshold" )){
		RESIDUAL_ENDPOINT_THRESHOLD  = color_db->getScalar<double>( "residual_endpoint_threshold" );
	}
	if (color_db->keyExists( "capillary_number" )){
		capillary_number = color_db->getScalar<double>( "capillary_number" );
		SET_CAPILLARY_NUMBER=true;
	}
	if (BoundaryCondition != 0 && SET_CAPILLARY_NUMBER==true){
		if (rank == 0) printf("WARINING: capillary number target only supported for BC = 0 \n");
		SET_CAPILLARY_NUMBER=false;
	}
	if (analysis_db->keyExists( "seed_water" )){
		seed_water = analysis_db->getScalar<double>( "seed_water" );
		USE_SEED = true;
	}
	if (analysis_db->keyExists( "morph_delta" )){
		morph_delta = analysis_db->getScalar<double>( "morph_delta" );
	}
	if (analysis_db->keyExists( "morph_interval" )){
		morph_interval = analysis_db->getScalar<int>( "morph_interval" );
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
		//analysis.run( timestep, *Averages, Phi, Pressure, Velocity, fq, Den );
		analysis.basic( timestep, *Averages, Phi, Pressure, Velocity, fq, Den );

		if (rank==0 && timestep%analysis_interval == 0 && BoundaryCondition > 0){
			printf("....inlet pressure=%f \n",din);
		}
		
		// allow initial ramp-up to get closer to steady state
		if (timestep > RAMP_TIMESTEPS && timestep%analysis_interval == 0 && USE_MORPH){
			analysis.finish();
			CURRENT_STEADY_TIMESTEPS += analysis_interval;

			double volB = Averages->gwb.V; 
			double volA = Averages->gnb.V; 
			volA /= Dm->Volume;
			volB /= Dm->Volume;;
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

				if ( isSteady ){
					MORPH_ADAPT = true;
					CURRENT_MORPH_TIMESTEPS=0;
					delta_volume_target = (volA )*morph_delta; // set target volume change
					Averages->Full();
					Averages->Write(timestep);
					analysis.WriteVisData( timestep, *Averages, Phi, Pressure, Velocity, fq, Den );
					analysis.finish();
					
					if (rank==0){
						printf("** WRITE STEADY POINT *** ");
						printf("Ca = %f, (previous = %f) \n",Ca,Ca_previous);
						double pA = Averages->gnb.p;
						double pB = Averages->gwb.p;

						double h = Dm->voxel_length;		
						double kAeff = h*h*muA*flow_rate_A/(rhoA*force_mag);
						double kBeff = h*h*muB*flow_rate_B/(rhoB*force_mag);
						double pAB = (pA-pB)/(h*5.796*alpha);
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
							fprintf(kr_log_file,"timesteps sat.water eff.perm.oil eff.perm.water cap.pressure pressure.drop Ca M\n",CURRENT_STEADY_TIMESTEPS,current_saturation,kAeff,kBeff,pAB,viscous_pressure_drop,Ca,Mobility);

						fprintf(kr_log_file,"%i %.5g %.5g %.5g %.5g %.5g %.5g %.5g\n",CURRENT_STEADY_TIMESTEPS,current_saturation,kAeff,kBeff,pAB,viscous_pressure_drop,Ca,Mobility);
						fclose(kr_log_file);

						printf("  Measured capillary number %f \n ",Ca);
						CURRENT_STEADY_TIMESTEPS = 0;
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
						if (force_mag < 1e-7){
							Fx *= 1e-7/force_mag;   // impose floor
							Fy *= 1e-7/force_mag;   
							Fz *= 1e-7/force_mag;   
						}
						if (rank == 0) printf("    -- adjust force by factor %f \n ",capillary_number / Ca);
						Averages->SetParams(rhoA,rhoB,tauA,tauB,Fx,Fy,Fz,alpha,beta);
					}
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
				if (USE_SEED){
					delta_volume = volA - initial_volume;
					CURRENT_MORPH_TIMESTEPS += analysis_interval;
					double massChange = SeedPhaseField(seed_water);
					if (rank==0) printf("***Seed water in oil %f, volume change %f / %f ***\n", seed_water, delta_volume, delta_volume_target);
				}
				else{
					if (rank==0) printf("***Morphological step with target volume change %f ***\n", delta_volume_target);
					//double delta_volume_target = volB - (volA + volB)*TARGET_SATURATION; // change in volume to A
					delta_volume += MorphInit(beta,delta_volume_target-delta_volume);
				}
				if ( (delta_volume - delta_volume_target)/delta_volume_target > 0.0 ){
					MORPH_ADAPT = false;
					initial_volume = volA;
					delta_volume = 0.0;
					CURRENT_STEADY_TIMESTEPS=0;
				}
				else if (CURRENT_MORPH_TIMESTEPS > MAX_MORPH_TIMESTEPS) {
					delta_volume = 0.0;
					initial_volume = volA;
					MORPH_ADAPT = false;
					CURRENT_STEADY_TIMESTEPS=0;
				}
				if ( REVERSE_FLOW_DIRECTION ){
					//if (rank==0) printf("*****REVERSE FLOW DIRECTION***** \n");
					delta_volume = 0.0;
					// flow direction will reverse after next steady point
					MORPH_ADAPT = false;
					CURRENT_STEADY_TIMESTEPS=0;
					//morph_delta *= (-1.0);
					REVERSE_FLOW_DIRECTION = false;
				}
				MPI_Barrier(comm);
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

double ScaLBL_ColorModel::MorphOpenConnected(double target_volume_change){
	
	int nx = Nx;
	int ny = Ny;
	int nz = Nz;
	int n;
	int N = nx*ny*nz;
	double volume_change;
	Array<char> id_solid(nx,ny,nz);
	Array<int> phase_label(nx,ny,nz);
	DoubleArray distance(Nx,Ny,Nz);
	DoubleArray phase(nx,ny,nz);
	signed char *id_connected;
	id_connected = new signed char [nx*ny*nz];
	
	ScaLBL_CopyToHost(phase.data(), Phi, N*sizeof(double));

	// Extract only the connected part of NWP
	BlobIDstruct new_index;
	double vF=0.0; double vS=0.0;
	ComputeGlobalBlobIDs(nx-2,ny-2,nz-2,Dm->rank_info,phase,Averages->SDs,vF,vS,phase_label,Dm->Comm);
	MPI_Barrier(Dm->Comm);
		
	int count_oil=0;
	int count_connected=0;
	int count_porespace=0;
	int count_water=0;
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
	count_connected=sumReduce( Dm->Comm, count_connected);
	count_porespace=sumReduce( Dm->Comm, count_porespace);
	count_water=sumReduce( Dm->Comm, count_water);
	
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

	double SW=0.5;
	// target water increase in voxels, normalized by connected volume
	double St = (SW*count_porespace - count_water)/count_porespace;  
	
	signed char water=2;
	signed char notwater=1;
	MorphOpen(distance, id_connected, Dm, St, water, notwater);
	
	for (int k=0; k<nz; k++){
		for (int j=0; j<ny; j++){
			for (int i=0; i<nx; i++){
				n=k*nx*ny+j*nx+i;
				// only apply opening to connected component 
				if ( id_connected[n] == 1){
					phase(i,j,k) = 1.0;
				}
			}
		}
	}
	
	
	
	ScaLBL_CopyToDevice(Phi,phase.data(),N*sizeof(double));
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
	return(volume_change);
}

double ScaLBL_ColorModel::SeedPhaseField(const double seed_water_in_oil){
	srand(time(NULL));
	double mass_loss =0.f;
	double count =0.f;
	DoubleArray phase(Nx,Ny,Nz);

	ScaLBL_CopyToHost(phase.data(), Phi, N*sizeof(double));
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				double random_value = double(rand())/ RAND_MAX;

				if (Averages->SDs(i,j,k) < 0.f){
					// skip
				}
				else if (phase(i,j,k) > 0.f ){
					phase(i,j,k) -= random_value*seed_water_in_oil;
					mass_loss += random_value*seed_water_in_oil;
					count++;
				}
				else {

				}
			}
		}
	}
	count= sumReduce( Dm->Comm, count);
	mass_loss= sumReduce( Dm->Comm, mass_loss);
	if (rank == 0) printf("Remove mass %f from %f voxels \n",mass_loss,count);
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
	return(mass_loss);
}

double ScaLBL_ColorModel::MorphInit(const double beta, const double target_delta_volume){
	const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);

	double vF = 0.f;
	double vS = 0.f;
	double delta_volume;

	DoubleArray phase(Nx,Ny,Nz);
	IntArray phase_label(Nx,Ny,Nz);;
	DoubleArray phase_distance(Nx,Ny,Nz);
	Array<char> phase_id(Nx,Ny,Nz);
	fillHalo<double> fillDouble(Dm->Comm,Dm->rank_info,{Nx-2,Ny-2,Nz-2},{1,1,1},0,1);
	

	// Basic algorithm to 
	// 1. Copy phase field to CPU
	ScaLBL_CopyToHost(phase.data(), Phi, N*sizeof(double));

	double count,count_global,volume_initial,volume_final,volume_connected;
	count = 0.f;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				if (phase(i,j,k) > 0.f && Averages->SDs(i,j,k) > 0.f) count+=1.f;
			}
		}
	}
	volume_initial = sumReduce( Dm->Comm, count);
	/*
	sprintf(LocalRankFilename,"phi_initial.%05i.raw",rank);
	FILE *INPUT = fopen(LocalRankFilename,"wb");
	fwrite(phase.data(),8,N,INPUT);
	fclose(INPUT);
	*/
	// 2. Identify connected components of phase field -> phase_label
	BlobIDstruct new_index;
	ComputeGlobalBlobIDs(Nx-2,Ny-2,Nz-2,rank_info,phase,Averages->SDs,vF,vS,phase_label,comm);
	MPI_Barrier(comm);
	
	// only operate on component "0"
	count = 0.0;
	double second_biggest = 0.0;

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
	volume_connected = sumReduce( Dm->Comm, count);
	second_biggest = sumReduce( Dm->Comm, second_biggest);

	int reach_x, reach_y, reach_z;
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
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
					// erase the original object
					phase(i,j,k) = -1.0;
				}
			}
		}
	}

	if (volume_connected - second_biggest < 2.0*fabs(target_delta_volume) && target_delta_volume < 0.0){
		// if connected volume is less than 2% just delete the whole thing
		if (rank==0) printf("Connected region has shrunk! \n");
		REVERSE_FLOW_DIRECTION = true;
	}
/*	else{*/
		if (rank==0) printf("Pathway volume / next largest ganglion %f \n",volume_connected/second_biggest );
		if (rank==0) printf("MorphGrow with target volume fraction change %f \n", target_delta_volume/volume_initial);
		double target_delta_volume_incremental = target_delta_volume;
		if (fabs(target_delta_volume) > 0.01*volume_initial)  
			target_delta_volume_incremental = 0.01*volume_initial*target_delta_volume/fabs(target_delta_volume);
		delta_volume = MorphGrow(Averages->SDs,phase_distance,phase_id,Averages->Dm, target_delta_volume_incremental);

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
					int n = k*Nx*Ny + j*Nx + i;
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
	volume_final= sumReduce( Dm->Comm, count);

	delta_volume = (volume_final-volume_initial);
	if (rank == 0)  printf("MorphInit: change fluid volume fraction by %f \n", delta_volume/volume_initial);
	if (rank == 0)  printf("   new saturation =  %f \n", volume_final/(0.238323*double((Nx-2)*(Ny-2)*(Nz-2)*nprocs)));

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
