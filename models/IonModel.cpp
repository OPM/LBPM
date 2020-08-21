/*
 * Dilute Ion Transport LBM Model
 */
#include "models/IonModel.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"

ScaLBL_IonModel::ScaLBL_IonModel(int RANK, int NP, MPI_Comm COMM):
rank(RANK),nprocs(NP),timestep(0),timestepMax(0),time_conv(0),kb(0),electron_charge(0),T(0),Vt(0),k2_inv(0),h(0),
tolerance(0),number_ion_species(0),Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),
BoundaryCondition(0),BoundaryConditionSolid(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{

}
ScaLBL_IonModel::~ScaLBL_IonModel(){

}

void ScaLBL_IonModel::ReadParams(string filename,int num_iter,int num_iter_Stokes,double time_conv_Stokes){
    
    // read the input database 
	db = std::make_shared<Database>( filename );
	domain_db = db->getDatabase( "Domain" );
	ion_db = db->getDatabase( "Ions" );
	
    //------ Load number of iteration from multiphysics controller ------//
	timestepMax = num_iter;
    //compute time conversion factor for ion model
    time_conv = num_iter_Stokes*time_conv_Stokes/num_iter;
    //-------------------------------------------------------------------//

    // Universal constant
    kb = 1.38e-23;//Boltzmann constant;unit [J/K]
    electron_charge = 1.6e-19;//electron charge;unit [C]

	//---------------------- Default model parameters --------------------------//		
    T = 300.0;//temperature; unit [K]
    Vt = kb*T/electron_charge;//thermal voltage; unit [Vy]
    k2_inv = 4.5;//speed of sound for D3Q7 lattice
    h = 1.0;//resolution; unit: um/lu
	tolerance = 1.0e-8;
	number_ion_species = 1;
    IonDiffusivity.push_back(1.0e-9);//user-input diffusivity has physical unit [m^2/sec]
    IonValence.push_back(1);//algebraic valence charge
    IonConcentration.push_back(1.0e-3);//user-input ion concentration has physical unit [mol/m^3]
    //deltaT.push_back(1.0);
    //tau.push_back(0.5+k2_inv*deltaT[0]*IonDiffusivity[0]);
    tau.push_back(0.5+k2_inv*time_conv/(h*1.0e-6)/(h*1.0e-6)*IonDiffusivity[0]);
    //--------------------------------------------------------------------------//

	// Read domain parameters
	if (domain_db->keyExists( "voxel_length" )){//default unit: um/lu
		h = domain_db->getScalar<double>( "voxel_length" );
	}
    BoundaryCondition = 0;
	if (domain_db->keyExists( "BC" )){
		BoundaryCondition = domain_db->getScalar<int>( "BC" );
	}

	// LB-Ion Model parameters		
	//if (ion_db->keyExists( "timestepMax" )){
	//	timestepMax = ion_db->getScalar<int>( "timestepMax" );
	//}
	if (ion_db->keyExists( "tolerance" )){
		tolerance = ion_db->getScalar<double>( "tolerance" );
	}
	if (ion_db->keyExists( "temperature" )){
		T = ion_db->getScalar<int>( "temperature" );
        //re-calculate thermal voltage 
        Vt = kb*T/electron_charge;//thermal voltage; unit [Vy]
	}
	if (ion_db->keyExists( "number_ion_species" )){
		number_ion_species = ion_db->getScalar<int>( "number_ion_species" );
	}
    //read ion related list
    //NOTE: ion diffusivity has INPUT unit: [m^2/sec]
    //      it must be converted to LB unit: [lu^2/lt]
	if (ion_db->keyExists("IonDiffusivityList")){
        IonDiffusivity.clear();
	    IonDiffusivity = ion_db->getVector<double>( "IonDiffusivityList" );
        // time relaxation parameters tau also needs update
        tau.clear();
        if (IonDiffusivity.size()!=number_ion_species){
		    ERROR("Error: number_ion_species and IonDiffusivityList must be the same length! \n");
        }
        else{
            for (int i=0; i<IonDiffusivity.size();i++){
                //First, convert ion diffusivity in physical unit to LB unit
                IonDiffusivity[i] = IonDiffusivity[i]*time_conv/(h*h*1.0e-12);//LB diffusivity has unit [lu^2/lt]
                //Second, re-calculate tau
                tau.push_back(0.5+k2_inv*IonDiffusivity[i]);
            }
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
    //read initial ion concentration list; INPUT unit [mol/m^3]
    //it must be converted to LB unit [mol/lu^3]
	if (ion_db->keyExists("IonConcentrationList")){
        IonConcentration.clear();
	    IonConcentration = ion_db->getVector<double>( "IonConcentrationList" );
        if (IonConcentration.size()!=number_ion_species){
		    ERROR("Error: number_ion_species and IonConcentrationList must be the same length! \n");
        }
        else{
            for (int i=0; i<IonConcentration.size();i++){
                IonConcentration[i] = IonConcentration[i]*(h*h*h*1.0e-18);//LB ion concentration has unit [mol/lu^3]
            }
        }
    }

    //Read solid boundary condition specific to Ion model
    BoundaryConditionSolid = 0;
	if (ion_db->keyExists( "BC_Solid" )){
		BoundaryConditionSolid = ion_db->getScalar<int>( "BC_Solid" );
	}


	if (rank==0) printf("*****************************************************\n");
	if (rank==0) printf("LB Ion Transport Solver: \n");
	if (rank==0) printf("      Time conversion factor: %.5g [sec/lt]\n", time_conv);
	if (rank==0) printf("      Internal iteration: %i [lt]\n", timestepMax);
    for (int i=0; i<number_ion_species;i++){
	    if (rank==0) printf("      Ion %i: LB relaxation tau = %.5g\n", i+1,tau[i]);
    }
	if (rank==0) printf("*****************************************************\n");

    switch (BoundaryConditionSolid){
        case 0:
          if (rank==0) printf("LB Ion Solver: solid boundary: non-flux boundary is assigned\n");  
          break;
        case 1:
          if (rank==0) printf("LB Ion Solver: solid boundary: Neumann-type surfacen ion concentration is assigned\n");  
          break;
        default:
          if (rank==0) printf("LB Ion Solver: solid boundary: non-flux boundary is assigned\n");  
          break;
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
	if (rank==0) printf("LB Ion Solver: Initialized solid phase & converting to Signed Distance function \n");
	CalcDist(Distance,id_solid,*Dm);
    if (rank == 0) cout << "    Domain set." << endl;
}


void ScaLBL_IonModel::AssignSolidBoundary(double *ion_solid)
{
	size_t NLABELS=0;
	signed char VALUE=0;
	double AFFINITY=0.f;

	auto LabelList = ion_db->getVector<int>( "SolidLabels" );
	auto AffinityList = ion_db->getVector<double>( "SolidValues" );

	NLABELS=LabelList.size();
	if (NLABELS != AffinityList.size()){
		ERROR("Error: LB Ion Solver: SolidLabels and SolidValues must be the same length! \n");
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
				      //printf("idx=%i, value=%i, %i, \n",idx, VALUE,LabelList[idx]);
					if (VALUE == LabelList[idx]){
						AFFINITY=AffinityList[idx];
                        //NOTE need to convert the user input phys unit to LB unit
                        AFFINITY = AFFINITY*(h*h*1.0e-12); 
						label_count[idx] += 1.0;
						idx = NLABELS;
						//Mask->id[n] = 0; // set mask to zero since this is an immobile component
					}
				}
				ion_solid[n] = AFFINITY;
			}
		}
	}

	for (size_t idx=0; idx<NLABELS; idx++)
		label_count_global[idx]=sumReduce( Dm->Comm, label_count[idx]);

	if (rank==0){
		printf("LB Ion Solver: number of ion solid labels: %lu \n",NLABELS);
		for (unsigned int idx=0; idx<NLABELS; idx++){
			VALUE=LabelList[idx];
			AFFINITY=AffinityList[idx];
			double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocs);
			printf("   label=%d, surface ion concentration=%.3g [mol/m^2], volume fraction=%.2g\n",VALUE,AFFINITY,volume_fraction); 
		}
	}
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
	if (rank==0)    printf ("LB Ion Solver: Create ScaLBL_Communicator \n");
	// Create a communicator for the device (will use optimized layout)
	// ScaLBL_Communicator ScaLBL_Comm(Mask); // original
	ScaLBL_Comm  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));

	int Npad=(Np/16 + 2)*16;
	if (rank==0)    printf ("LB Ion Solver: Set up memory efficient layout \n");
	Map.resize(Nx,Ny,Nz);       Map.fill(-2);
	auto neighborList= new int[18*Npad];
	Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Mask->id,Np);
	MPI_Barrier(comm);
	//...........................................................................
	//                MAIN  VARIABLES ALLOCATED HERE
	//...........................................................................
	// LBM variables
	if (rank==0)    printf ("LB Ion Solver: Allocating distributions \n");
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
	if (rank==0)    printf ("LB Ion Solver: Setting up device map and neighbor list \n");
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	MPI_Barrier(comm);
	
    //Initialize solid boundary for electrical potential
    //if ion concentration at solid surface is specified
    if (BoundaryConditionSolid==1){

	    ScaLBL_AllocateDeviceMemory((void **) &IonSolid, sizeof(double)*Nx*Ny*Nz);
        ScaLBL_Comm->SetupBounceBackList(Map, Mask->id, Np);
        MPI_Barrier(comm);

        double *IonSolid_host;
        IonSolid_host = new double[Nx*Ny*Nz];
        AssignSolidBoundary(IonSolid_host);
        ScaLBL_CopyToDevice(IonSolid, IonSolid_host, Nx*Ny*Nz*sizeof(double));
        ScaLBL_DeviceBarrier();
        delete [] IonSolid_host;
    }


}        

void ScaLBL_IonModel::Initialize(){
	/*
	 * This function initializes model
	 */
    if (rank==0)    printf ("LB Ion Solver: initializing D3Q7 distributions\n");
	for (int ic=0; ic<number_ion_species; ic++){
        ScaLBL_D3Q7_Ion_Init(&fq[ic*Np*7],&Ci[ic*Np],IonConcentration[ic],Np); 
    }
    if (rank==0)    printf ("LB Ion Solver: initializing charge density\n");
	for (int ic=0; ic<number_ion_species; ic++){
        ScaLBL_D3Q7_Ion_ChargeDensity(Ci, ChargeDensity, IonValence[ic], ic, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_D3Q7_Ion_ChargeDensity(Ci, ChargeDensity, IonValence[ic], ic, 0, ScaLBL_Comm->LastExterior(), Np);
    }
}

void ScaLBL_IonModel::Run(double *Velocity, double *ElectricField){

    //LB-related parameter
    vector<double> rlx(tau.begin(),tau.end());
    for (double item : rlx){
        item = 1.0/item; 
    }
	//.......create and start timer............
	//double starttime,stoptime,cputime;
	//ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
	//starttime = MPI_Wtime();

	timestep=0;
	while (timestep < timestepMax) {
		//************************************************************************/
		// *************ODD TIMESTEP*************//
		timestep++;
        //Update ion concentration and charge density
		for (int ic=0; ic<number_ion_species; ic++){
			ScaLBL_Comm->SendD3Q7AA(fq, ic); //READ FROM NORMAL
            ScaLBL_D3Q7_AAodd_IonConcentration(NeighborList, &fq[ic*Np*7],&Ci[ic*Np],ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
			ScaLBL_Comm->RecvD3Q7AA(fq, ic); //WRITE INTO OPPOSITE
            ScaLBL_DeviceBarrier();
            ScaLBL_D3Q7_AAodd_IonConcentration(NeighborList, &fq[ic*Np*7],&Ci[ic*Np], 0, ScaLBL_Comm->LastExterior(), Np);
            ScaLBL_D3Q7_Ion_ChargeDensity(Ci, ChargeDensity, IonValence[ic], ic, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
            ScaLBL_D3Q7_Ion_ChargeDensity(Ci, ChargeDensity, IonValence[ic], ic, 0, ScaLBL_Comm->LastExterior(), Np);
        }

        //LB-Ion collison
		for (int ic=0; ic<number_ion_species; ic++){
            ScaLBL_D3Q7_AAodd_Ion(NeighborList, &fq[ic*Np*7],&Ci[ic*Np],Velocity,ElectricField,IonDiffusivity[ic],IonValence[ic],
                                  rlx[ic],Vt,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        }

		// Set boundary conditions
		/* ... */

		for (int ic=0; ic<number_ion_species; ic++){
            ScaLBL_D3Q7_AAodd_Ion(NeighborList, &fq[ic*Np*7],&Ci[ic*Np],Velocity,ElectricField,IonDiffusivity[ic],IonValence[ic],
                                  rlx[ic],Vt,0, ScaLBL_Comm->LastExterior(), Np);
        }
        if (BoundaryConditionSolid==1){
		    for (int ic=0; ic<number_ion_species; ic++){
                //TODO IonSolid may also be species-dependent
		        ScaLBL_Comm->SolidNeumannD3Q7(&fq[ic*Np*7], IonSolid);
		        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
            }
        }

		// *************EVEN TIMESTEP*************//
		timestep++;
        //Update ion concentration and charge density
		for (int ic=0; ic<number_ion_species; ic++){
			ScaLBL_Comm->SendD3Q7AA(fq, ic); //READ FORM NORMAL
            ScaLBL_D3Q7_AAeven_IonConcentration(&fq[ic*Np*7],&Ci[ic*Np],ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
			ScaLBL_Comm->RecvD3Q7AA(fq, ic); //WRITE INTO OPPOSITE
            ScaLBL_D3Q7_AAeven_IonConcentration(&fq[ic*Np*7],&Ci[ic*Np], 0, ScaLBL_Comm->LastExterior(), Np);
            ScaLBL_D3Q7_Ion_ChargeDensity(Ci, ChargeDensity, IonValence[ic], ic, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
            ScaLBL_D3Q7_Ion_ChargeDensity(Ci, ChargeDensity, IonValence[ic], ic, 0, ScaLBL_Comm->LastExterior(), Np);
        }

        //LB-Ion collison
		for (int ic=0; ic<number_ion_species; ic++){
            ScaLBL_D3Q7_AAeven_Ion(&fq[ic*Np*7],&Ci[ic*Np],Velocity,ElectricField,IonDiffusivity[ic],IonValence[ic],
                                  rlx[ic],Vt,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        }

		// Set boundary conditions
		/* ... */
		
		for (int ic=0; ic<number_ion_species; ic++){
            ScaLBL_D3Q7_AAeven_Ion(&fq[ic*Np*7],&Ci[ic*Np],Velocity,ElectricField,IonDiffusivity[ic],IonValence[ic],
                                  rlx[ic],Vt,0, ScaLBL_Comm->LastExterior(), Np);
        }
        if (BoundaryConditionSolid==1){
		    for (int ic=0; ic<number_ion_species; ic++){
                //TODO IonSolid may also be species-dependent
		        ScaLBL_Comm->SolidNeumannD3Q7(&fq[ic*Np*7], IonSolid);
		        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
            }
        }
		//************************************************************************/
	}
	//************************************************************************/
	//stoptime = MPI_Wtime();
	//if (rank==0) printf("-------------------------------------------------------------------\n");
	//// Compute the walltime per timestep
	//cputime = (stoptime - starttime)/timestep;
	//// Performance obtained from each node
	//double MLUPS = double(Np)/cputime/1000000;

	//if (rank==0) printf("********************************************************\n");
	//if (rank==0) printf("CPU time = %f \n", cputime);
	//if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
	//MLUPS *= nprocs;
	//if (rank==0) printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
	//if (rank==0) printf("********************************************************\n");

}

void ScaLBL_IonModel::getIonConcentration(){
	for (int ic=0; ic<number_ion_species; ic++){
        ScaLBL_IonConcentration_Phys(Ci, h, ic, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_IonConcentration_Phys(Ci, h, ic, 0, ScaLBL_Comm->LastExterior(), Np);
    }

    DoubleArray PhaseField(Nx,Ny,Nz);
	for (int ic=0; ic<number_ion_species; ic++){
	    ScaLBL_Comm->RegularLayout(Map,&Ci[ic*Np],PhaseField);
        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

        FILE *OUTFILE;
        sprintf(LocalRankFilename,"Ion%02i.%05i.raw",ic+1,rank);
        OUTFILE = fopen(LocalRankFilename,"wb");
        fwrite(PhaseField.data(),8,N,OUTFILE);
        fclose(OUTFILE);
    }

}

