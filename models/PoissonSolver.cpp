/*
 * Multi-relaxation time LBM Model
 */
#include "models/PoissonSolver.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"

ScaLBL_Poisson::ScaLBL_Poisson(int RANK, int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP),timestep(0),timestepMax(0),tau(0),k2_inv(0),tolerance(0),h(0),
epsilon0(0),epsilon0_LB(0),epsilonR(0),epsilon_LB(0),Vin(0),Vout(0),Nx(0),Ny(0),Nz(0),N(0),Np(0),analysis_interval(0),
chargeDen_dummy(0),WriteLog(0),
nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),BoundaryConditionSolid(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{

}
ScaLBL_Poisson::~ScaLBL_Poisson(){

}

void ScaLBL_Poisson::ReadParams(string filename){
	// read the input database 
	db = std::make_shared<Database>( filename );
	domain_db = db->getDatabase( "Domain" );
	electric_db = db->getDatabase( "Poisson" );
	
    k2_inv = 4.0;//speed of sound for D3Q7 lattice 
	tau = 0.5+k2_inv;
	timestepMax = 100000;
	tolerance = 1.0e-6;//stopping criterion for obtaining steady-state electricla potential
    h = 1.0;//resolution; unit: um/lu
    epsilon0 = 8.85e-12;//electric permittivity of vaccum; unit:[C/(V*m)]
    epsilon0_LB = epsilon0*(h*1.0e-6);//unit:[C/(V*lu)]
    epsilonR = 78.4;//default dielectric constant of water
    epsilon_LB = epsilon0_LB*epsilonR;//electric permittivity 
    analysis_interval = 1000; 
    Vin  = 1.0; //Boundary-z (inlet)  electric potential
    Vout = 1.0; //Boundary-Z (outlet) electric potential
    chargeDen_dummy = 1.0e-3;//For debugging;unit=[C/m^3]
    WriteLog = false;

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
	if (electric_db->keyExists( "epsilonR" )){
		epsilonR = electric_db->getScalar<double>( "epsilonR" );
	}
	if (electric_db->keyExists( "DummyChargeDen" )){
		chargeDen_dummy = electric_db->getScalar<double>( "DummyChargeDen" );
	}
	if (electric_db->keyExists( "WriteLog" )){
		WriteLog = electric_db->getScalar<bool>( "WriteLog" );
	}

    // Read solid boundary condition specific to Poisson equation
    BoundaryConditionSolid = 1;
	if (electric_db->keyExists( "BC_Solid" )){
		BoundaryConditionSolid = electric_db->getScalar<int>( "BC_Solid" );
	}
    // Read boundary condition for electric potential
    // BC = 0: normal periodic BC
    // BC = 1: fixed inlet and outlet potential
    BoundaryCondition = 0;
	if (electric_db->keyExists( "BC" )){
		BoundaryCondition = electric_db->getScalar<int>( "BC" );
	}

	// Read domain parameters
	if (domain_db->keyExists( "voxel_length" )){//default unit: um/lu
		h = domain_db->getScalar<double>( "voxel_length" );
	}


    //Re-calcualte model parameters if user updates input
    epsilon0_LB = epsilon0*(h*1.0e-6);//unit:[C/(V*lu)]
    epsilon_LB = epsilon0_LB*epsilonR;//electric permittivity 

	if (rank==0) printf("***********************************************************************************\n");
	if (rank==0) printf("LB-Poisson Solver: steady-state MaxTimeStep = %i; steady-state tolerance = %.3g \n", timestepMax,tolerance);
	if (rank==0) printf("                   LB relaxation tau = %.5g \n", tau);
	if (rank==0) printf("***********************************************************************************\n");

    switch (BoundaryConditionSolid){
        case 1:
          if (rank==0) printf("LB-Poisson Solver: solid boundary: Dirichlet-type surfacen potential is assigned\n");  
          break;
        case 2:
          if (rank==0) printf("LB-Poisson Solver: solid boundary: Neumann-type surfacen charge density is assigned\n");  
          break;
        default:
          if (rank==0) printf("LB-Poisson Solver: solid boundary: Dirichlet-type surfacen potential is assigned\n");  
          break;
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

	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;               // initialize this way
	//Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm) ); // TwoPhase analysis object
	MPI_Barrier(comm);
	Dm->BoundaryCondition = BoundaryCondition;
	Mask->BoundaryCondition = BoundaryCondition;
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

void ScaLBL_Poisson::AssignSolidBoundary(double *poisson_solid)
{
	size_t NLABELS=0;
	signed char VALUE=0;
	double AFFINITY=0.f;

	auto LabelList = electric_db->getVector<int>( "SolidLabels" );
	auto AffinityList = electric_db->getVector<double>( "SolidValues" );

	NLABELS=LabelList.size();
	if (NLABELS != AffinityList.size()){
		ERROR("Error: LB-Poisson Solver: SolidLabels and SolidValues must be the same length! \n");
	}

	double label_count[NLABELS];
	double label_count_global[NLABELS];
	// Assign the labels

	for (size_t idx=0; idx<NLABELS; idx++) label_count[idx]=0;

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
			}
		}
	}

	for (size_t idx=0; idx<NLABELS; idx++)
		label_count_global[idx]=sumReduce( Dm->Comm, label_count[idx]);

	if (rank==0){
		printf("LB-Poisson Solver: number of Poisson solid labels: %lu \n",NLABELS);
		for (unsigned int idx=0; idx<NLABELS; idx++){
			VALUE=LabelList[idx];
			AFFINITY=AffinityList[idx];
			double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocs);
            switch (BoundaryConditionSolid){
                case 1:
			        printf("   label=%d, surface potential=%.3g [V], volume fraction=%.2g\n",VALUE,AFFINITY,volume_fraction); 
                    break;
                case 2: 
			        printf("   label=%d, surface charge density=%.3g [C/m^2], volume fraction=%.2g\n",VALUE,AFFINITY,volume_fraction); 
                    break;
                default:
			        printf("   label=%d, surface potential=%.3g [V], volume fraction=%.2g\n",VALUE,AFFINITY,volume_fraction); 
                    break;
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
	ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);
	//ScaLBL_AllocateDeviceMemory((void **) &dvcID, sizeof(signed char)*Nx*Ny*Nz);
	ScaLBL_AllocateDeviceMemory((void **) &fq, 7*dist_mem_size);  
	ScaLBL_AllocateDeviceMemory((void **) &Psi, sizeof(double)*Nx*Ny*Nz);
	ScaLBL_AllocateDeviceMemory((void **) &ElectricField, 3*sizeof(double)*Np);
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
	ScaLBL_DeviceBarrier();
	delete [] TmpMap;
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	ScaLBL_DeviceBarrier();
	MPI_Barrier(comm);
	delete [] neighborList;
    // copy node ID
	//ScaLBL_CopyToDevice(dvcID, Mask->id, sizeof(signed char)*Nx*Ny*Nz);
	//ScaLBL_DeviceBarrier();
	
    //Initialize solid boundary for electric potential
    ScaLBL_Comm->SetupBounceBackList(Map, Mask->id, Np);
	MPI_Barrier(comm);
}        

void ScaLBL_Poisson::Potential_Init(double *psi_init){

    if (BoundaryCondition==1){
	    if (electric_db->keyExists( "Vin" )){
	    	Vin = electric_db->getScalar<double>( "Vin" );
	    }
	    if (electric_db->keyExists( "Vout" )){
	    	Vout = electric_db->getScalar<double>( "Vout" );
	    }
    }
    //By default only periodic BC is applied and Vin=Vout=1.0, i.e. there is no potential gradient along Z-axis
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


void ScaLBL_Poisson::Initialize(){
	/*
	 * This function initializes model
	 */
    if (rank==0)    printf ("LB-Poisson Solver: initializing D3Q7 distributions\n");
    //NOTE the initialization involves two steps:
    //1. assign solid boundary value (surface potential or surface change density)
    //2. Initialize electric potential for pore nodes
    double *psi_host;
    psi_host = new double [Nx*Ny*Nz];
    AssignSolidBoundary(psi_host);//step1
    Potential_Init(psi_host);//step2
	ScaLBL_CopyToDevice(Psi, psi_host, Nx*Ny*Nz*sizeof(double));
	ScaLBL_DeviceBarrier();
    ScaLBL_D3Q7_Poisson_Init(dvcMap, fq, Psi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
    ScaLBL_D3Q7_Poisson_Init(dvcMap, fq, Psi, 0, ScaLBL_Comm->LastExterior(), Np);
    delete [] psi_host;
    
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
        
        SolveElectricPotentialAAodd();//update electric potential
        SolvePoissonAAodd(ChargeDensity);//perform collision
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

		// *************EVEN TIMESTEP*************//
		timestep++;
		SolveElectricPotentialAAeven();//update electric potential
        SolvePoissonAAeven(ChargeDensity);//perform collision
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		//************************************************************************/

        // Check convergence of steady-state solution
        if (timestep%analysis_interval==0){
        
			//ScaLBL_Comm->RegularLayout(Map,Psi,Psi_host);
            ScaLBL_CopyToHost(Psi_host.data(),Psi,sizeof(double)*Nx*Ny*Nz);
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
    if(WriteLog==true){
        getConvergenceLog(timestep,error);
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

void ScaLBL_Poisson::getConvergenceLog(int timestep,double error){
    if (rank==0){
		bool WriteHeader=false;
		TIMELOG = fopen("PoissonSolver_Convergence.csv","r");
		if (TIMELOG != NULL)
			fclose(TIMELOG);
		else
			WriteHeader=true;

		TIMELOG = fopen("PoissonSolver_Convergence.csv","a+");
		if (WriteHeader)
		{
			fprintf(TIMELOG,"Timestep Error\n");				
            fprintf(TIMELOG,"%i %.5g\n",timestep,error); 
            fflush(TIMELOG);
		}
        else {
            fprintf(TIMELOG,"%i %.5g\n",timestep,error); 
            fflush(TIMELOG);
        }
    }
}

void ScaLBL_Poisson::SolveElectricPotentialAAodd(){
	ScaLBL_Comm->SendD3Q7AA(fq, 0); //READ FROM NORMAL
	ScaLBL_D3Q7_AAodd_Poisson_ElectricPotential(NeighborList, dvcMap, fq, Psi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
	ScaLBL_Comm->RecvD3Q7AA(fq, 0); //WRITE INTO OPPOSITE
    ScaLBL_DeviceBarrier();
	// Set boundary conditions
	if (BoundaryCondition == 1){
		ScaLBL_Comm->D3Q7_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
		ScaLBL_Comm->D3Q7_Poisson_Potential_BC_Z(NeighborList, fq, Vout, timestep);
	}
    //-------------------------//
	ScaLBL_D3Q7_AAodd_Poisson_ElectricPotential(NeighborList, dvcMap, fq, Psi, 0, ScaLBL_Comm->LastExterior(), Np);
}

void ScaLBL_Poisson::SolveElectricPotentialAAeven(){
	ScaLBL_Comm->SendD3Q7AA(fq, 0); //READ FORM NORMAL
	ScaLBL_D3Q7_AAeven_Poisson_ElectricPotential(dvcMap, fq, Psi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
	ScaLBL_Comm->RecvD3Q7AA(fq, 0); //WRITE INTO OPPOSITE
    ScaLBL_DeviceBarrier();
	// Set boundary conditions
	if (BoundaryCondition == 1){
		ScaLBL_Comm->D3Q7_Poisson_Potential_BC_z(NeighborList, fq,  Vin, timestep);
		ScaLBL_Comm->D3Q7_Poisson_Potential_BC_Z(NeighborList, fq, Vout, timestep);
	}
    //-------------------------//
	ScaLBL_D3Q7_AAeven_Poisson_ElectricPotential(dvcMap, fq, Psi, 0, ScaLBL_Comm->LastExterior(), Np);
}

void ScaLBL_Poisson::SolvePoissonAAodd(double *ChargeDensity){
	ScaLBL_D3Q7_AAodd_Poisson(NeighborList, dvcMap, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
	ScaLBL_D3Q7_AAodd_Poisson(NeighborList, dvcMap, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, 0, ScaLBL_Comm->LastExterior(), Np);
    if (BoundaryConditionSolid==1){
	    ScaLBL_Comm->SolidDirichletD3Q7(fq, Psi);
    }
    else if (BoundaryConditionSolid==2){
	    ScaLBL_Comm->SolidNeumannD3Q7(fq, Psi);
    }
}

void ScaLBL_Poisson::SolvePoissonAAeven(double *ChargeDensity){
	ScaLBL_D3Q7_AAeven_Poisson(dvcMap, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
	ScaLBL_D3Q7_AAeven_Poisson(dvcMap, fq, ChargeDensity, Psi, ElectricField, tau, epsilon_LB, 0, ScaLBL_Comm->LastExterior(), Np);
    if (BoundaryConditionSolid==1){
	    ScaLBL_Comm->SolidDirichletD3Q7(fq, Psi);
    }
    else if (BoundaryConditionSolid==2){
	    ScaLBL_Comm->SolidNeumannD3Q7(fq, Psi);
    }
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
	ScaLBL_DeviceBarrier();
	delete [] ChargeDensity_host;
}

void ScaLBL_Poisson::getElectricPotential_debug(int timestep){
    //This function write out decomposed data
    DoubleArray PhaseField(Nx,Ny,Nz);
	//ScaLBL_Comm->RegularLayout(Map,Psi,PhaseField);
    ScaLBL_CopyToHost(PhaseField.data(),Psi,sizeof(double)*Nx*Ny*Nz);
    //ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
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
        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

	    ScaLBL_Comm->RegularLayout(Map,&ElectricField[1*Np],Values_y);
        ElectricField_LB_to_Phys(Values_y);
        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

	    ScaLBL_Comm->RegularLayout(Map,&ElectricField[2*Np],Values_z);
        ElectricField_LB_to_Phys(Values_z);
        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

}

void ScaLBL_Poisson::getElectricField_debug(int timestep){

        //ScaLBL_D3Q7_Poisson_getElectricField(fq,ElectricField,tau,Np);
        //ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

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

//void ScaLBL_Poisson::SolveElectricField(){
//    ScaLBL_Comm_Regular->SendHalo(Psi);
//    ScaLBL_D3Q7_Poisson_ElectricField(NeighborList, dvcMap, dvcID, Psi, ElectricField, BoundaryConditionSolid,
//                                      Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
//    ScaLBL_Comm_Regular->RecvHalo(Psi);
//    ScaLBL_DeviceBarrier();
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
//        //ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
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
//		label_count_global[idx]=sumReduce( Dm->Comm, label_count[idx]);
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
