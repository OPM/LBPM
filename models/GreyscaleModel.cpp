/*
  Copyright 2020 Equinor ASA

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

  Greyscale lattice boltzmann model
 */
#include "models/GreyscaleModel.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"
#include <stdlib.h>
#include <time.h>

template<class TYPE>
void DeleteArray( const TYPE *p )
{
    delete [] p;
}

ScaLBL_GreyscaleModel::ScaLBL_GreyscaleModel(int RANK, int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tau(0),tau_eff(0),Den(0),Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),GreyPorosity(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{
	SignDist.resize(Nx,Ny,Nz);           
    SignDist.fill(0);

}
ScaLBL_GreyscaleModel::~ScaLBL_GreyscaleModel(){

}

void ScaLBL_GreyscaleModel::ReadParams(string filename){
	// read the input database 
	db = std::make_shared<Database>( filename );
	domain_db = db->getDatabase( "Domain" );
	greyscale_db =  db->getDatabase( "Greyscale" );
	analysis_db = db->getDatabase( "Analysis" );
	vis_db = db->getDatabase( "Visualization" );

	// set defaults
	timestepMax = 100000;
	tau = 1.0;
    tau_eff = tau;
    Den = 1.0;//constant density
	tolerance = 0.01;
	Fx = Fy = Fz = 0.0;
	Restart=false;
	din=dout=1.0;
	flux=0.0;
    dp = 10.0; //unit of 'dp': voxel
    CollisionType = 1; //1: IMRT; 2: BGK; 3: MRT
	
	// ---------------------- Greyscale Model parameters -----------------------//
	if (greyscale_db->keyExists( "timestepMax" )){
		timestepMax = greyscale_db->getScalar<int>( "timestepMax" );
	}
	if (greyscale_db->keyExists( "tau" )){
		tau = greyscale_db->getScalar<double>( "tau" );
	}
	tau_eff = greyscale_db->getWithDefault<double>( "tau_eff", tau );
	if (greyscale_db->keyExists( "Den" )){
		Den = greyscale_db->getScalar<double>( "Den" );
	}
	if (greyscale_db->keyExists( "dp" )){
		dp = greyscale_db->getScalar<double>( "dp" );
	}
	if (greyscale_db->keyExists( "F" )){
		Fx = greyscale_db->getVector<double>( "F" )[0];
		Fy = greyscale_db->getVector<double>( "F" )[1];
		Fz = greyscale_db->getVector<double>( "F" )[2];
	}
	if (greyscale_db->keyExists( "Restart" )){
		Restart = greyscale_db->getScalar<bool>( "Restart" );
	}
	if (greyscale_db->keyExists( "din" )){
		din = greyscale_db->getScalar<double>( "din" );
	}
	if (greyscale_db->keyExists( "dout" )){
		dout = greyscale_db->getScalar<double>( "dout" );
	}
	if (greyscale_db->keyExists( "flux" )){
		flux = greyscale_db->getScalar<double>( "flux" );
	}
	if (greyscale_db->keyExists( "tolerance" )){
		tolerance = greyscale_db->getScalar<double>( "tolerance" );
	}
	auto collision = greyscale_db->getWithDefault<std::string>( "collision", "IMRT" );
	if (collision == "BGK"){
        CollisionType=2;
	}
    else if (collision == "MRT"){
        CollisionType=3;
    }
	// ------------------------------------------------------------------------//
    
    //------------------------ Other Domain parameters ------------------------//
	BoundaryCondition = 0;
	if (domain_db->keyExists( "BC" )){
		BoundaryCondition = domain_db->getScalar<int>( "BC" );
	}
	// ------------------------------------------------------------------------//
}

void ScaLBL_GreyscaleModel::SetDomain(){
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

	SignDist.resize(Nx,Ny,Nz);
	Velocity_x.resize(Nx,Ny,Nz);
	Velocity_y.resize(Nx,Ny,Nz);
	Velocity_z.resize(Nx,Ny,Nz);
	PorosityMap.resize(Nx,Ny,Nz);
	Pressure.resize(Nx,Ny,Nz);

	id = new signed char [N];
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;               // initialize this way
	MPI_Barrier(comm);
	Dm->CommInit();
	MPI_Barrier(comm);
	// Read domain parameters
	rank = Dm->rank();	
	nprocx = Dm->nprocx();
	nprocy = Dm->nprocy();
	nprocz = Dm->nprocz();
}

void ScaLBL_GreyscaleModel::ReadInput(){
	
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
	
	if (domain_db->keyExists( "Filename" )){
		auto Filename = domain_db->getScalar<std::string>( "Filename" );
		Mask->Decomp(Filename);
	}
	else{
        if (rank==0) printf("Filename of input image is not found, reading ID.0* instead.");
		Mask->ReadIDs();
	}
	for (int i=0; i<Nx*Ny*Nz; i++) id[i] = Mask->id[i];  // save what was read
	
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
				SignDist(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
			}
		}
	}
//	MeanFilter(SignDist);
	if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
	CalcDist(SignDist,id_solid,*Mask);
	
	if (rank == 0) cout << "Domain set." << endl;
}

/********************************************************
 * AssignComponentLabels                                 *
 ********************************************************/
void ScaLBL_GreyscaleModel::AssignComponentLabels(double *Porosity, double *Permeability)
{
	size_t NLABELS=0;
	signed char VALUE=0;
	double POROSITY=0.f;
	double PERMEABILITY=0.f;

	auto LabelList = greyscale_db->getVector<int>( "ComponentLabels" );
	auto PorosityList = greyscale_db->getVector<double>( "PorosityList" );
	auto PermeabilityList = greyscale_db->getVector<double>( "PermeabilityList" );

	NLABELS=LabelList.size();
	if (NLABELS != PorosityList.size()){
		ERROR("Error: ComponentLabels and PorosityList must be the same length! \n");
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
						POROSITY=PorosityList[idx];
						label_count[idx] += 1.0;
						idx = NLABELS;
						//Mask->id[n] = 0; // set mask to zero since this is an immobile component
					}
				}
				int idx = Map(i,j,k);
				if (!(idx < 0)){
                    if (POROSITY<=0.0){
                        ERROR("Error: Porosity for grey voxels must be 0.0 < Porosity <= 1.0 !\n");
                    }
                    else{
					    Porosity[idx] = POROSITY;
                    }
                }
			}
		}
	}

	if (NLABELS != PermeabilityList.size()){
		ERROR("Error: ComponentLabels and PermeabilityList must be the same length! \n");
	}
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				int n = k*Nx*Ny+j*Nx+i;
				VALUE=id[n];
				// Assign the affinity from the paired list
				for (unsigned int idx=0; idx < NLABELS; idx++){
					//printf("idx=%i, value=%i, %i, \n",idx, VALUE,LabelList[idx]);
					if (VALUE == LabelList[idx]){
						PERMEABILITY=PermeabilityList[idx];
						idx = NLABELS;
						//Mask->id[n] = 0; // set mask to zero since this is an immobile component
					}
				}
				int idx = Map(i,j,k);
				if (!(idx < 0)){
                    if (PERMEABILITY<=0.0){
                        ERROR("Error: Permeability for grey voxel must be > 0.0 ! \n");
                    }
                    else{
					    Permeability[idx] = PERMEABILITY/Dm->voxel_length/Dm->voxel_length;
                    }
                }
			}
		}
	}


	// Set Dm to match Mask
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i]; 
	
	for (int idx=0; idx<NLABELS; idx++)		label_count_global[idx]=sumReduce( Dm->Comm, label_count[idx]);
    //Initialize a weighted porosity after considering grey voxels
    GreyPorosity=0.0;
	for (unsigned int idx=0; idx<NLABELS; idx++){
		double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocs);
        GreyPorosity+=volume_fraction*PorosityList[idx];
    }

	if (rank==0){
        printf("Image resolution: %.5g [um/voxel]\n",Dm->voxel_length);
		printf("Number of component labels: %lu \n",NLABELS);
		for (unsigned int idx=0; idx<NLABELS; idx++){
			VALUE=LabelList[idx];
			POROSITY=PorosityList[idx];
			PERMEABILITY=PermeabilityList[idx];
			double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocs);
			printf("   label=%d: porosity=%.3g, permeability=%.3g [um^2] (=%.3g [voxel^2]), volume fraction=%.3g\n",
                    VALUE,POROSITY,PERMEABILITY,PERMEABILITY/Dm->voxel_length/Dm->voxel_length,volume_fraction); 
            printf("             effective porosity=%.3g\n",volume_fraction*POROSITY);
		}
        printf("The weighted porosity, considering both open and grey voxels, is %.3g\n",GreyPorosity);
	}
}


void ScaLBL_GreyscaleModel::Create(){
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
	ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Permeability, sizeof(double)*Np);		
	ScaLBL_AllocateDeviceMemory((void **) &Porosity, sizeof(double)*Np);		
	ScaLBL_AllocateDeviceMemory((void **) &Pressure_dvc, sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
	//...........................................................................
	// Update GPU data structures
	if (rank==0)	printf ("Setting up device neighbor list \n");
	fflush(stdout);
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	// initialize phi based on PhaseLabel (include solid component labels)
	double *Poros, *Perm;
	Poros = new double[Np];
	Perm  = new double[Np];
	AssignComponentLabels(Poros,Perm);
	ScaLBL_CopyToDevice(Porosity, Poros, Np*sizeof(double));
	ScaLBL_CopyToDevice(Permeability, Perm, Np*sizeof(double));
    delete [] Poros;
    delete [] Perm;
}        


void ScaLBL_GreyscaleModel::Initialize(){
	if (rank==0)	printf ("Initializing distributions \n");
    //TODO: for BGK, you need to consider voxel porosity
    //      for IMRT, the whole set of feq is different
    //      if in the future you have different collison mode, need to write two set of initialization functions
    if (CollisionType==1){
	    ScaLBL_D3Q19_GreyIMRT_Init(fq, Np, Den);
        if (rank==0) printf("Collision model: Incompressible MRT.\n");
    }
    else if (CollisionType==2){
	    ScaLBL_D3Q19_Init(fq, Np);
        if (rank==0) printf("Collision model: BGK.\n");
    }
    else if (CollisionType==3){
	    ScaLBL_D3Q19_Init(fq, Np);
        if (rank==0) printf("Collision model: MRT.\n");
    }
    else{
        if (rank==0) printf("Unknown collison type! IMRT collision is used.\n"); 
	    ScaLBL_D3Q19_GreyIMRT_Init(fq, Np, Den);
        CollisionType=1;
		greyscale_db->putScalar<std::string>( "collision", "IMRT" );
    }

	if (Restart == true){
		if (rank==0){
			printf("Initializing distributions from Restart! \n");
		}
		// Read in the restart file to CPU buffers
        std::shared_ptr<double> cfq;
        cfq = std::shared_ptr<double>(new double[19*Np],DeleteArray<double>);
        FILE *File;
        File=fopen(LocalRestartFile,"rb");
        fread(cfq.get(),sizeof(double),19*Np,File);
        fclose(File);

		// Copy the restart data to the GPU
		ScaLBL_CopyToDevice(fq,cfq.get(),19*Np*sizeof(double));
		ScaLBL_DeviceBarrier();

		MPI_Barrier(comm);
	}
}

void ScaLBL_GreyscaleModel::Run(){
	int nprocs=nprocx*nprocy*nprocz;
	const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
	
	int analysis_interval = 1000; 	// number of timesteps in between in situ analysis 
	int visualization_interval = 1000; 	 
	int restart_interval = 10000; 	// number of timesteps in between in saving distributions for restart 
	if (analysis_db->keyExists( "analysis_interval" )){
		analysis_interval = analysis_db->getScalar<int>( "analysis_interval" );
	}
	if (analysis_db->keyExists( "visualization_interval" )){
		visualization_interval = analysis_db->getScalar<int>( "visualization_interval" );
	}
	if (analysis_db->keyExists( "restart_interval" )){
		restart_interval = analysis_db->getScalar<int>( "restart_interval" );
	}
	if (greyscale_db->keyExists( "timestep" )){
		timestep = greyscale_db->getScalar<int>( "timestep" );
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
	
	Minkowski Morphology(Mask);

	//************ MAIN ITERATION LOOP ***************************************/
	PROFILE_START("Loop");
	auto current_db = db->cloneDatabase();
	double rlx = 1.0/tau;
    double rlx_eff = 1.0/tau_eff;
	double error = 1.0;
	double flow_rate_previous = 0.0;
	while (timestep < timestepMax && error > tolerance) {
		//************************************************************************/
		// *************ODD TIMESTEP*************//
		timestep++;
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        switch (CollisionType){
            case 1: 
                    ScaLBL_D3Q19_AAodd_Greyscale_IMRT(NeighborList, fq,  ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
            case 2: 
                    ScaLBL_D3Q19_AAodd_Greyscale(NeighborList, fq,  ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Pressure_dvc);
                    break;
            case 3: 
                    ScaLBL_D3Q19_AAodd_Greyscale_MRT(NeighborList, fq,  ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
            default: 
                    ScaLBL_D3Q19_AAodd_Greyscale_IMRT(NeighborList, fq,  ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
        }
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();
		// Set BCs
		if (BoundaryCondition == 3){
			ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
        switch (CollisionType){
            case 1: 
		            ScaLBL_D3Q19_AAodd_Greyscale_IMRT(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
            case 2: 
		            ScaLBL_D3Q19_AAodd_Greyscale(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Pressure_dvc);
                    break;
            case 3: 
		            ScaLBL_D3Q19_AAodd_Greyscale_MRT(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
            default: 
		            ScaLBL_D3Q19_AAodd_Greyscale_IMRT(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
        }
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

		// *************EVEN TIMESTEP*************//
		timestep++;
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
        switch (CollisionType){
            case 1: 
		            ScaLBL_D3Q19_AAeven_Greyscale_IMRT(fq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
            case 2: 
		            ScaLBL_D3Q19_AAeven_Greyscale(fq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Pressure_dvc);
                    break;
            case 3: 
		            ScaLBL_D3Q19_AAeven_Greyscale_MRT(fq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
            default: 
		            ScaLBL_D3Q19_AAeven_Greyscale_IMRT(fq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
        }
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();
		// Set BCs
		if (BoundaryCondition == 3){
			ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
		}
        switch (CollisionType){
            case 1: 
		            ScaLBL_D3Q19_AAeven_Greyscale_IMRT(fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
            case 2: 
		            ScaLBL_D3Q19_AAeven_Greyscale(fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Pressure_dvc);
                    break;
            case 3: 
		            ScaLBL_D3Q19_AAeven_Greyscale_MRT(fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
            default: 
		            ScaLBL_D3Q19_AAeven_Greyscale_IMRT(fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx, rlx_eff, Fx, Fy, Fz,Porosity,Permeability,Velocity,Den,Pressure_dvc);
                    break;
        }
        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		//************************************************************************/
		
		if (timestep%analysis_interval==0){
			ScaLBL_Comm->RegularLayout(Map,&Velocity[0],Velocity_x);
			ScaLBL_Comm->RegularLayout(Map,&Velocity[Np],Velocity_y);
			ScaLBL_Comm->RegularLayout(Map,&Velocity[2*Np],Velocity_z);
			//ScaLBL_Comm->RegularLayout(Map,Porosity,PorosityMap);
			//ScaLBL_Comm->RegularLayout(Map,Pressure_dvc,Pressure);
			
			double count_loc=0;
			double count;
			double vax,vay,vaz;
			double vax_loc,vay_loc,vaz_loc;
            //double px_loc,py_loc,pz_loc;
            //double px,py,pz;
            //double mass_loc,mass_glb;
            
            //parameters for domain average
	        int64_t i,j,k,n,imin,jmin,kmin,kmax;
            // If external boundary conditions are set, do not average over the inlet and outlet
            kmin=1; kmax=Nz-1;
            //In case user forgets to specify the inlet/outlet buffer layers for BC>0
            if (BoundaryCondition > 0 && Dm->kproc() == 0) kmin=4;
            if (BoundaryCondition > 0 && Dm->kproc() == Dm->nprocz()-1) kmax=Nz-4;

            imin=jmin=1;
            // If inlet/outlet layers exist use these as default
            //if (Dm->inlet_layers_x > 0) imin = Dm->inlet_layers_x;
            //if (Dm->inlet_layers_y > 0) jmin = Dm->inlet_layers_y;
            if (BoundaryCondition > 0 && Dm->inlet_layers_z > 0 && Dm->kproc() == 0) kmin = 1 + Dm->inlet_layers_z;//"1" indicates the halo layer
            if (BoundaryCondition > 0 && Dm->outlet_layers_z > 0 && Dm->kproc() == Dm->nprocz()-1) kmax = Nz-1 - Dm->outlet_layers_z; 

//			px_loc = py_loc = pz_loc = 0.f;
//            mass_loc = 0.f;
//			for (int k=kmin; k<kmax; k++){
//				for (int j=jmin; j<Ny-1; j++){
//					for (int i=imin; i<Nx-1; i++){
//						if (SignDist(i,j,k) > 0){
//							px_loc   += Velocity_x(i,j,k)*Den*PorosityMap(i,j,k);
//							py_loc   += Velocity_y(i,j,k)*Den*PorosityMap(i,j,k);
//							pz_loc   += Velocity_z(i,j,k)*Den*PorosityMap(i,j,k);
//							mass_loc += Den*PorosityMap(i,j,k);
//						}
//					}
//				}
//			}
//			MPI_Allreduce(&px_loc,  &px,      1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
//			MPI_Allreduce(&py_loc,  &py,      1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
//			MPI_Allreduce(&pz_loc,  &pz,      1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
//			MPI_Allreduce(&mass_loc,&mass_glb,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
//			
//			vax = px/mass_glb;
//			vay = py/mass_glb;
//			vaz = pz/mass_glb;
            
			vax_loc = vay_loc = vaz_loc = 0.f;
			for (int k=kmin; k<kmax; k++){
				for (int j=jmin; j<Ny-1; j++){
					for (int i=imin; i<Nx-1; i++){
						if (SignDist(i,j,k) > 0){
							vax_loc += Velocity_x(i,j,k);
							vay_loc += Velocity_y(i,j,k);
							vaz_loc += Velocity_z(i,j,k);
							count_loc+=1.0;
						}
					}
				}
			}
            vax = sumReduce( Mask->Comm, vax_loc);
            vay = sumReduce( Mask->Comm, vay_loc);
            vaz = sumReduce( Mask->Comm, vaz_loc);
            count = sumReduce( Mask->Comm, count_loc);

			vax /= count;
			vay /= count;
			vaz /= count;

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
			//double flow_rate = (px*dir_x + py*dir_y + pz*dir_z)/mass_glb;
			double flow_rate = (vax*dir_x + vay*dir_y + vaz*dir_z);
			
			error = fabs(flow_rate - flow_rate_previous) / fabs(flow_rate);
			flow_rate_previous = flow_rate;
			
			//if (rank==0) printf("Computing Minkowski functionals \n");
			Morphology.ComputeScalar(SignDist,0.f);
			//Morphology.PrintAll();
			double mu = (tau-0.5)/3.f;
			double Vs = Morphology.V();
			double As = Morphology.A();
			double Hs = Morphology.H();
			double Xs = Morphology.X();
			Vs = sumReduce( Dm->Comm, Vs);
			As = sumReduce( Dm->Comm, As);
			Hs = sumReduce( Dm->Comm, Hs);
			Xs = sumReduce( Dm->Comm, Xs);

			double h = Dm->voxel_length;
			//double absperm = h*h*mu*Mask->Porosity()*flow_rate / force_mag;
			double absperm = h*h*mu*GreyPorosity*flow_rate / force_mag;

            if (rank==0){
				printf("     AbsPerm = %.5g [micron^2]\n",absperm);
                bool WriteHeader=false;
                FILE * log_file = fopen("Permeability.csv","r");
                if (log_file != NULL)
                    fclose(log_file);
                else
                    WriteHeader=true;
                log_file = fopen("Permeability.csv","a");
                if (WriteHeader)
                    fprintf(log_file,"timestep Fx Fy Fz mu Vs As Hs Xs vax vay vaz AbsPerm \n",
                            timestep,Fx,Fy,Fz,mu,h*h*h*Vs,h*h*As,h*Hs,Xs,vax,vay,vaz,absperm);

                fprintf(log_file,"%i %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",timestep, Fx, Fy, Fz, mu, 
                        h*h*h*Vs,h*h*As,h*Hs,Xs,vax,vay,vaz, absperm);
                fclose(log_file);
            }
		}

		if (timestep%visualization_interval==0){
            VelocityField();
        }

		if (timestep%restart_interval==0){
            //Use rank=0 write out Restart.db
            if (rank==0) {
                greyscale_db->putScalar<int>("timestep",timestep);    		
                greyscale_db->putScalar<bool>( "Restart", true );
                current_db->putDatabase("Greyscale", greyscale_db);
                std::ofstream OutStream("Restart.db");
                current_db->print(OutStream, "");
                OutStream.close();
      
            }
            //Write out Restart data.
            std::shared_ptr<double> cfq;
            cfq = std::shared_ptr<double>(new double[19*Np],DeleteArray<double>);
            ScaLBL_CopyToHost(cfq.get(),fq,19*Np*sizeof(double));// Copy restart data to the CPU

            FILE *RESTARTFILE;
            RESTARTFILE=fopen(LocalRestartFile,"wb");
            fwrite(cfq.get(),sizeof(double),19*Np,RESTARTFILE);
            fclose(RESTARTFILE);
		    MPI_Barrier(comm);
        }
	}

	PROFILE_STOP("Loop");
	PROFILE_SAVE("lbpm_greyscale_simulator",1);
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

void ScaLBL_GreyscaleModel::VelocityField(){

/*	Minkowski Morphology(Mask);
	int SIZE=Np*sizeof(double);
	ScaLBL_D3Q19_Momentum(fq,Velocity, Np);
	ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
	ScaLBL_CopyToHost(&VELOCITY[0],&Velocity[0],3*SIZE);

	memcpy(Morphology.SDn.data(), Distance.data(), Nx*Ny*Nz*sizeof(double));
	Morphology.Initialize();
	Morphology.UpdateMeshValues();
	Morphology.ComputeLocal();
	Morphology.Reduce();
	
	double count_loc=0;
	double count;
	double vax,vay,vaz;
	double vax_loc,vay_loc,vaz_loc;
	vax_loc = vay_loc = vaz_loc = 0.f;
	for (int n=0; n<ScaLBL_Comm->LastExterior(); n++){
		vax_loc += VELOCITY[n];
		vay_loc += VELOCITY[Np+n];
		vaz_loc += VELOCITY[2*Np+n];
		count_loc+=1.0;
	}
	
	for (int n=ScaLBL_Comm->FirstInterior(); n<ScaLBL_Comm->LastInterior(); n++){
		vax_loc += VELOCITY[n];
		vay_loc += VELOCITY[Np+n];
		vaz_loc += VELOCITY[2*Np+n];
		count_loc+=1.0;
	}
	MPI_Allreduce(&vax_loc,&vax,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	MPI_Allreduce(&vay_loc,&vay,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	MPI_Allreduce(&vaz_loc,&vaz,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	MPI_Allreduce(&count_loc,&count,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
	
	vax /= count;
	vay /= count;
	vaz /= count;
	
	double mu = (tau-0.5)/3.f;
	if (rank==0) printf("Fx Fy Fz mu Vs As Js Xs vx vy vz\n");
	if (rank==0) printf("%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",Fx, Fy, Fz, mu, 
						Morphology.V(),Morphology.A(),Morphology.J(),Morphology.X(),vax,vay,vaz);
						*/
	
	std::vector<IO::MeshDataStruct> visData;
	fillHalo<double> fillData(Dm->Comm,Dm->rank_info,{Dm->Nx-2,Dm->Ny-2,Dm->Nz-2},{1,1,1},0,1);

	auto VxVar = std::make_shared<IO::Variable>();
	auto VyVar = std::make_shared<IO::Variable>();
	auto VzVar = std::make_shared<IO::Variable>();
	auto SignDistVar = std::make_shared<IO::Variable>();
	auto PressureVar = std::make_shared<IO::Variable>();

	IO::initialize("","silo","false");
	// Create the MeshDataStruct	
	visData.resize(1);
	visData[0].meshName = "domain";
	visData[0].mesh = std::make_shared<IO::DomainMesh>( Dm->rank_info,Dm->Nx-2,Dm->Ny-2,Dm->Nz-2,Dm->Lx,Dm->Ly,Dm->Lz );
	SignDistVar->name = "SignDist";
	SignDistVar->type = IO::VariableType::VolumeVariable;
	SignDistVar->dim = 1;
	SignDistVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
	visData[0].vars.push_back(SignDistVar);
	
	VxVar->name = "Velocity_x";
	VxVar->type = IO::VariableType::VolumeVariable;
	VxVar->dim = 1;
	VxVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
	visData[0].vars.push_back(VxVar);
	VyVar->name = "Velocity_y";
	VyVar->type = IO::VariableType::VolumeVariable;
	VyVar->dim = 1;
	VyVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
	visData[0].vars.push_back(VyVar);
	VzVar->name = "Velocity_z";
	VzVar->type = IO::VariableType::VolumeVariable;
	VzVar->dim = 1;
	VzVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
	visData[0].vars.push_back(VzVar);
	
	PressureVar->name = "Pressure";
	PressureVar->type = IO::VariableType::VolumeVariable;
	PressureVar->dim = 1;
	PressureVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
	visData[0].vars.push_back(PressureVar);

	Array<double>& SignData  = visData[0].vars[0]->data;
	Array<double>& VelxData = visData[0].vars[1]->data;
	Array<double>& VelyData = visData[0].vars[2]->data;
	Array<double>& VelzData = visData[0].vars[3]->data;
	Array<double>& PressureData = visData[0].vars[4]->data;
	
    ASSERT(visData[0].vars[0]->name=="SignDist");
    ASSERT(visData[0].vars[1]->name=="Velocity_x");
    ASSERT(visData[0].vars[2]->name=="Velocity_y");
    ASSERT(visData[0].vars[3]->name=="Velocity_z");
    ASSERT(visData[0].vars[4]->name=="Pressure");
	
	ScaLBL_Comm->RegularLayout(Map,&Velocity[0],Velocity_x);
	ScaLBL_Comm->RegularLayout(Map,&Velocity[Np],Velocity_y);
	ScaLBL_Comm->RegularLayout(Map,&Velocity[2*Np],Velocity_z);
	ScaLBL_Comm->RegularLayout(Map,Pressure_dvc,Pressure);

    fillData.copy(SignDist,SignData);
    fillData.copy(Velocity_x,VelxData);
    fillData.copy(Velocity_y,VelyData);
    fillData.copy(Velocity_z,VelzData);
    fillData.copy(Pressure,PressureData);
	
    IO::writeData( timestep, visData, Dm->Comm );

}

void ScaLBL_GreyscaleModel::WriteDebug(){
	// Copy back final phase indicator field and convert to regular layout
	DoubleArray PhaseField(Nx,Ny,Nz);

	//ScaLBL_CopyToHost(Porosity.data(), Poros, sizeof(double)*N);

//	FILE *OUTFILE;
//	sprintf(LocalRankFilename,"Phase.%05i.raw",rank);
//	OUTFILE = fopen(LocalRankFilename,"wb");
//	fwrite(PhaseField.data(),8,N,OUTFILE);
//	fclose(OUTFILE);
//
//    ScaLBL_Comm->RegularLayout(Map,&Den[0],PhaseField);
//	FILE *AFILE;
//	sprintf(LocalRankFilename,"A.%05i.raw",rank);
//	AFILE = fopen(LocalRankFilename,"wb");
//	fwrite(PhaseField.data(),8,N,AFILE);
//	fclose(AFILE);
//
//	ScaLBL_Comm->RegularLayout(Map,&Den[Np],PhaseField);
//	FILE *BFILE;
//	sprintf(LocalRankFilename,"B.%05i.raw",rank);
//	BFILE = fopen(LocalRankFilename,"wb");
//	fwrite(PhaseField.data(),8,N,BFILE);
//	fclose(BFILE);
//
//	ScaLBL_Comm->RegularLayout(Map,Pressure,PhaseField);
//	FILE *PFILE;
//	sprintf(LocalRankFilename,"Pressure.%05i.raw",rank);
//	PFILE = fopen(LocalRankFilename,"wb");
//	fwrite(PhaseField.data(),8,N,PFILE);
//	fclose(PFILE);

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

	ScaLBL_Comm->RegularLayout(Map,&Porosity[0],PhaseField);
	FILE *POROS_FILE;
	sprintf(LocalRankFilename,"Porosity.%05i.raw",rank);
	POROS_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,POROS_FILE);
	fclose(POROS_FILE);

	ScaLBL_Comm->RegularLayout(Map,&Permeability[0],PhaseField);
	FILE *PERM_FILE;
	sprintf(LocalRankFilename,"Permeability.%05i.raw",rank);
	PERM_FILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,PERM_FILE);
	fclose(PERM_FILE);
}
