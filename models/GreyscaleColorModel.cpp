/*
Greyscale lattice boltzmann model
 */
#include "models/GreyscaleColorModel.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"
#include <stdlib.h>
#include <time.h>

template<class TYPE>
void DeleteArray( const TYPE *p )
{
    delete [] p;
}

ScaLBL_GreyscaleColorModel::ScaLBL_GreyscaleColorModel(int RANK, int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tauA(0),tauB(0),tauA_eff(0),tauB_eff(0),
rhoA(0),rhoB(0),gamma(0),kappaA(0),kappaB(0),lambdaA(0),lambdaB(0),
Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),GreyPorosity(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{
	SignDist.resize(Nx,Ny,Nz);           
    SignDist.fill(0);

}
ScaLBL_GreyscaleColorModel::~ScaLBL_GreyscaleColorModel(){

}

void ScaLBL_GreyscaleColorModel::ReadParams(string filename){
	// read the input database 
	db = std::make_shared<Database>( filename );
	domain_db = db->getDatabase( "Domain" );
	greyscaleColor_db =  db->getDatabase( "GreyscaleColor" );
	analysis_db = db->getDatabase( "Analysis" );
	vis_db = db->getDatabase( "Visualization" );

	// set defaults
	timestepMax = 100000;
	tauA = 1.0;
	tauB = 1.0;
    tauA_eff = tauA;//the effective viscosity of the Darcy term
    tauB_eff = tauB;
    rhoA = 1.0;//constant molecular mass (after LB scaling)
    rhoB = 1.0;
	tolerance = 0.01;
	Fx = Fy = Fz = 0.0;
	Restart=false;
	din=dout=1.0;
	flux=0.0;
    dp = 10.0; //unit of 'dp': voxel
    //Gsc = 1.0;
    gamma = 1.0;//may also have gammaA and gammaB;
    kappaA = 1.0e-3;
    kappaB = 1.0e-3;
    lambdaA = 1.0e-3;
    lambdaB = 1.0e-3;
	
	// ---------------------- Greyscale Model parameters -----------------------//
	if (greyscaleColor_db->keyExists( "timestepMax" )){
		timestepMax = greyscaleColor_db->getScalar<int>( "timestepMax" );
	}
	if (greyscaleColor_db->keyExists( "tauA" )){
		tauA = greyscaleColor_db->getScalar<double>( "tauA" );
	}
	if (greyscaleColor_db->keyExists( "tauB" )){
		tauB = greyscaleColor_db->getScalar<double>( "tauB" );
	}
	tauA_eff = greyscaleColor_db->getWithDefault<double>( "tauA_eff", tauA);
	tauB_eff = greyscaleColor_db->getWithDefault<double>( "tauB_eff", tauB);
	if (greyscaleColor_db->keyExists( "rhoA" )){
		rhoA = greyscaleColor_db->getScalar<double>( "rhoA" );
	}
	if (greyscaleColor_db->keyExists( "rhoB" )){
		rhoB = greyscaleColor_db->getScalar<double>( "rhoB" );
	}
//	if (greyscaleColor_db->keyExists( "Gsc" )){
//		Gsc = greyscaleColor_db->getScalar<double>( "Gsc" );
//	}
	if (greyscaleColor_db->keyExists( "gamma" )){
		gamma = greyscaleColor_db->getScalar<double>( "gamma" );
	}
	if (greyscaleColor_db->keyExists( "kappaA" )){
		kappaA = greyscaleColor_db->getScalar<double>( "kappaA" );
	}
	if (greyscaleColor_db->keyExists( "kappaB" )){
		kappaB = greyscaleColor_db->getScalar<double>( "kappaB" );
	}
	if (greyscaleColor_db->keyExists( "lambdaA" )){
		lambdaA = greyscaleColor_db->getScalar<double>( "lambdaA" );
	}
	if (greyscaleColor_db->keyExists( "lambdaB" )){
		lambdaB = greyscaleColor_db->getScalar<double>( "lambdaB" );
	}
	if (greyscaleColor_db->keyExists( "dp" )){
		dp = greyscaleColor_db->getScalar<double>( "dp" );
	}
	if (greyscaleColor_db->keyExists( "F" )){
		Fx = greyscaleColor_db->getVector<double>( "F" )[0];
		Fy = greyscaleColor_db->getVector<double>( "F" )[1];
		Fz = greyscaleColor_db->getVector<double>( "F" )[2];
	}
	if (greyscaleColor_db->keyExists( "Restart" )){
		Restart = greyscaleColor_db->getScalar<bool>( "Restart" );
	}
	if (greyscaleColor_db->keyExists( "din" )){
		din = greyscaleColor_db->getScalar<double>( "din" );
	}
	if (greyscaleColor_db->keyExists( "dout" )){
		dout = greyscaleColor_db->getScalar<double>( "dout" );
	}
	if (greyscaleColor_db->keyExists( "flux" )){
		flux = greyscaleColor_db->getScalar<double>( "flux" );
	}
	if (greyscaleColor_db->keyExists( "tolerance" )){
		tolerance = greyscaleColor_db->getScalar<double>( "tolerance" );
	}
	//auto collision = greyscaleColor_db->getWithDefault<std::string>( "collision", "IMRT" );
	//if (collision == "BGK"){
    //    CollisionType=2;
	//}
	// ------------------------------------------------------------------------//
    
    //------------------------ Other Domain parameters ------------------------//
	BoundaryCondition = 0;
	if (domain_db->keyExists( "BC" )){
		BoundaryCondition = domain_db->getScalar<int>( "BC" );
	}
	// ------------------------------------------------------------------------//
}

void ScaLBL_GreyscaleColorModel::SetDomain(){
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

void ScaLBL_GreyscaleColorModel::ReadInput(){
	
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

void ScaLBL_GreyscaleColorModel::AssignSolidForce(double *SolidPotential, double *SolidForce){

	double *Dst;
	Dst = new double [3*3*3];
	for (int kk=0; kk<3; kk++){
		for (int jj=0; jj<3; jj++){
			for (int ii=0; ii<3; ii++){
				int index = kk*9+jj*3+ii;
				Dst[index] = sqrt(double(ii-1)*double(ii-1) + double(jj-1)*double(jj-1)+ double(kk-1)*double(kk-1));
			}
		}
	}
	double w_face = 1.f/18.f;
	double w_edge = 1.f/36.f;
	double w_corner = 0.f;
	//local 
	Dst[13] = 0.f;
	//faces
	Dst[4] = w_face;
	Dst[10] = w_face;
	Dst[12] = w_face;
	Dst[14] = w_face;
	Dst[16] = w_face;
	Dst[22] = w_face;
	// corners
	Dst[0] = w_corner;
	Dst[2] = w_corner;
	Dst[6] = w_corner;
	Dst[8] = w_corner;
	Dst[18] = w_corner;
	Dst[20] = w_corner;
	Dst[24] = w_corner;
	Dst[26] = w_corner;
	// edges
	Dst[1] = w_edge;
	Dst[3] = w_edge;
	Dst[5] = w_edge;
	Dst[7] = w_edge;
	Dst[9] = w_edge;
	Dst[11] = w_edge;
	Dst[15] = w_edge;
	Dst[17] = w_edge;
	Dst[19] = w_edge;
	Dst[21] = w_edge;
	Dst[23] = w_edge;
	Dst[25] = w_edge;

	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				int idx=Map(i,j,k);
				if (!(idx < 0)){
					double phi_x = 0.f;
					double phi_y = 0.f;
					double phi_z = 0.f;
					for (int kk=0; kk<3; kk++){
						for (int jj=0; jj<3; jj++){
							for (int ii=0; ii<3; ii++){

								int index = kk*9+jj*3+ii;
								double weight= Dst[index];

								int idi=i+ii-1;
								int idj=j+jj-1;
								int idk=k+kk-1;

								if (idi < 0) idi=0;
								if (idj < 0) idj=0;
								if (idk < 0) idk=0;
								if (!(idi < Nx)) idi=Nx-1;
								if (!(idj < Ny)) idj=Ny-1;
								if (!(idk < Nz)) idk=Nz-1;

								int nn = idk*Nx*Ny + idj*Nx + idi;
								if ((Mask->id[nn] <= 0)||(Mask->id[nn]>=3)){
									double vec_x = double(ii-1);
									double vec_y = double(jj-1);
									double vec_z = double(kk-1);
									double GWNS=SolidPotential[nn];
									phi_x += GWNS*weight*vec_x;
									phi_y += GWNS*weight*vec_y;
									phi_z += GWNS*weight*vec_z;
								}
							}
						}
					}
					SolidForce[idx] = phi_x;
					SolidForce[idx+Np] = phi_y;
					SolidForce[idx+2*Np] = phi_z;
				}
			}
		}
	}
	delete [] Dst;
}


void ScaLBL_GreyscaleColorModel::AssignComponentLabels(double *Porosity, double *Permeability, double *SolidPotential)
{
	size_t NLABELS=0;
	signed char VALUE=0;
	double POROSITY=0.f;
	double PERMEABILITY=0.f;
	double AFFINITY=0.f;

	auto LabelList = greyscaleColor_db->getVector<int>( "ComponentLabels" );
	auto AffinityList = greyscaleColor_db->getVector<double>( "ComponentAffinity" );
	auto PorosityList = greyscaleColor_db->getVector<double>( "PorosityList" );
	auto PermeabilityList = greyscaleColor_db->getVector<double>( "PermeabilityList" );

    //1. Requirement for "ComponentLabels":
    //   *labels can be a nagative integer, 0, 1, 2, or a positive integer >= 3
    //   *label = 1 and 2 are reserved for NW and W phase respectively.
    //2. Requirement for "ComponentAffinity":
    //   *should be in the same length as "ComponentLabels"
    //   *could leave Affinity=0.0 for label=1 and 2
    //3. Requirement for "PorosityList":
    //   *for ComponentLables <=0, put porosity value = 0.0;
    //   *for ComponentLabels >=3, put the corresponding sub-resolution porosity
    //   *for ComponentLabels =1, 2, put porosity=1 (or if users accidentally put other values it should still be fine)
    //4. Requirement for "PermeabilityList":
    //   *for ComponentLabels <=2, does not matter, can leave it as 1.0 

	NLABELS=LabelList.size();
	if (NLABELS != PorosityList.size() || NLABELS != AffinityList.size() || NLABELS != PermeabilityList.size()){
		ERROR("Error: ComponentLabels, ComponentAffinity, PorosityList and PermeabilityList must all be the same length! \n");
	}

	double label_count[NLABELS];
	double label_count_global[NLABELS];

	for (int idx=0; idx<NLABELS; idx++) label_count[idx]=0;

    //Populate the poroisty map, NOTE only for node_ID > 0, i.e. open or grey nodes
    //For node_ID <= 0: these are solid nodes of various wettability
	for (int k=1;k<Nz-1;k++){
		for (int j=1;j<Ny-1;j++){
			for (int i=1;i<Nx-1;i++){
				int n = k*Nx*Ny+j*Nx+i;
				VALUE=id[n];
				for (unsigned int idx=0; idx < NLABELS; idx++){
					if ((VALUE>0) && (VALUE == LabelList[idx])){
						POROSITY=PorosityList[idx];
						label_count[idx] += 1.0;
						idx = NLABELS;
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

    //Populate the permeability map, NOTE only for node_ID > 0, i.e. open or grey nodes
    //For node_ID <= 0: these are solid nodes of various wettability
	for (int k=1;k<Nz-1;k++){
		for (int j=1;j<Ny-1;j++){
			for (int i=1;i<Nx-1;i++){
				int n = k*Nx*Ny+j*Nx+i;
				VALUE=id[n];
				// Assign the affinity from the paired list
				for (unsigned int idx=0; idx < NLABELS; idx++){
					//printf("idx=%i, value=%i, %i, \n",idx, VALUE,LabelList[idx]);
					if ( (VALUE>0) && (VALUE == LabelList[idx])){
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

    //Populate the solid potential map, for ALL range of node_ID except node = 1,2, i.e. NW and W phase
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				int n = k*Nx*Ny+j*Nx+i;
				VALUE=id[n];
				// Assign the affinity from the paired list
				for (unsigned int idx=0; idx < NLABELS; idx++){
					if (VALUE == LabelList[idx]){
                        if (VALUE<=0){
                            AFFINITY=AffinityList[idx];
                        }
                        else if (VALUE>=3){
						    AFFINITY=AffinityList[idx]*(1.0-PorosityList[idx]);//BE CAREFUL! Requires for node_ID<=0, user puts porosity=0.0
                        }
                        else{//i.e. label = 1 or 2
                            AFFINITY=0.0;
                        } 
						idx = NLABELS;
					}
				}
				//NOTE: node_ID = 1 and 2 are reserved
				if ((VALUE == 1)||(VALUE == 2)) AFFINITY=0.0;//NOTE: still need this as users may forget to put label=1,2 in ComponentLabelLists
				SolidPotential[n] = AFFINITY;
			}
		}
	}


	// Set Dm to match Mask
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i]; 
	
	for (int idx=0; idx<NLABELS; idx++)		label_count_global[idx]=Dm->Comm.sumReduce(label_count[idx]);

    //Initialize a weighted porosity after considering grey voxels
    GreyPorosity=0.0;
	for (unsigned int idx=0; idx<NLABELS; idx++){
		double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocs);
        GreyPorosity+=volume_fraction*PorosityList[idx];
    }

	if (rank==0){
        printf("Image resolution: %.5g [um/voxel]\n",Dm->voxel_length);
		printf("Component labels: %lu \n",NLABELS);
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

void ScaLBL_GreyscaleColorModel::Density_and_Phase_Init(){

	size_t NLABELS=0;
	signed char VALUE=0;

    vector<int> LabelList{1,2};
    vector<double> SwList{0.0,1.0}; 

	if (greyscaleColor_db->keyExists( "GreyNodeLabels" )){
        LabelList.clear();
	    LabelList = greyscaleColor_db->getVector<int>( "GreyNodeLabels" );
	}
	if (greyscaleColor_db->keyExists( "GreyNodeSw" )){
        SwList.clear();
	    SwList = greyscaleColor_db->getVector<double>( "GreyNodeSw" );
	}

	NLABELS=LabelList.size();
	if (NLABELS != SwList.size()){
		ERROR("Error: GreyNodeLabels and GreyNodeSw must be the same length! \n");
	}
	
    double *Den_temp;
	Den_temp=new double [2*Np];
    double nA=0.5;//to prevent use may forget to specify all greynodes, then must initialize something to start with, givning just zeros is too risky.
    double nB=0.5;

    double *Phi_temp;
    Phi_temp=new double [Np];
    double phi = 0.0;

	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				int n = k*Nx*Ny+j*Nx+i;
				VALUE=Mask->id[n];
                if (VALUE>0){
                    for (unsigned int idx=0; idx < NLABELS; idx++){
                        if (VALUE == LabelList[idx]){
                            double Sw = SwList[idx];
                            if ((Sw<0.0) || (Sw>1.0)) ERROR("Error: Initial saturation for grey nodes must be between [0.0, 1.0]! \n");
                            nB=Sw;
                            nA=1.0-Sw;
                            phi = nA-nB;
                            idx = NLABELS;
                        }
                    }
                    if (VALUE==1){//label=1 reserved for NW phase
                        nA=1.0;
                        nB=0.0; 
                        phi = nA-nB;
                    }
                    else if(VALUE==2){//label=2 reserved for W phase
                        nA=0.0;
                        nB=1.0; 
                        phi = nA-nB;
                    }
                    int idx = Map(i,j,k);
                    Den_temp[idx+0*Np] = nA;
                    Den_temp[idx+1*Np] = nB;
                    Phi_temp[idx] = phi;
                }
			}
		}
	}
    //copy to device
	ScaLBL_CopyToDevice(Den, Den_temp, 2*Np*sizeof(double));
	ScaLBL_CopyToDevice(Phi, Phi_temp, 1*Np*sizeof(double));
	ScaLBL_DeviceBarrier();
	delete [] Den_temp;
}

void ScaLBL_GreyscaleColorModel::Create(){
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
	ScaLBL_AllocateDeviceMemory((void **) &PressureGrad, 3*sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
	//ScaLBL_AllocateDeviceMemory((void **) &Aq, 7*sizeof(double)*Np);
	//ScaLBL_AllocateDeviceMemory((void **) &Bq, 7*sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &Cq, 7*sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &Den, 2*sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &PhiLap, sizeof(double)*Np);//laplacian of phase field		
	ScaLBL_AllocateDeviceMemory((void **) &SolidForce, 3*sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &PressTensor, 6*sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &PressTensorGrad, 18*sizeof(double)*Np);
	//ScaLBL_AllocateDeviceMemory((void **) &DenGradA, 3*sizeof(double)*Np);
	//ScaLBL_AllocateDeviceMemory((void **) &DenGradB, 3*sizeof(double)*Np);
	//ScaLBL_AllocateDeviceMemory((void **) &DenLapA, sizeof(double)*Np);
	//ScaLBL_AllocateDeviceMemory((void **) &DenLapB, sizeof(double)*Np);
	//...........................................................................
	// Update GPU data structures
	if (rank==0)	printf ("Setting up device neighbor list \n");
	fflush(stdout);
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	// initialize phi based on PhaseLabel (include solid component labels)
	double *Poros, *Perm;
	Poros = new double[Np];
	Perm = new double[Np];
	double *SolidForce_host = new double[3*Np];
    double *SolidPotential_host = new double [Nx*Ny*Nz];
	AssignComponentLabels(Poros,Perm,SolidPotential_host);
    AssignSolidForce(SolidPotential_host,SolidForce_host);
	ScaLBL_CopyToDevice(Porosity, Poros, Np*sizeof(double));
	ScaLBL_CopyToDevice(Permeability, Perm, Np*sizeof(double));
	ScaLBL_CopyToDevice(SolidForce, SolidForce_host, 3*Np*sizeof(double));

	ScaLBL_DeviceBarrier();
    //TODO make the following smart pointers
	delete [] SolidForce_host;
	delete [] SolidPotential_host;
    delete [] Poros;
    delete [] Perm;
}        


void ScaLBL_GreyscaleColorModel::Initialize(){
	if (Restart == true){
		if (rank==0){
			printf("Initializing density field and distributions from Restart! \n");
		}
		// Read in the restart file to CPU buffers
        std::shared_ptr<double> cfq;
        cfq = std::shared_ptr<double>(new double[19*Np],DeleteArray<double>);
        std::shared_ptr<double> cDen;
        cDen = std::shared_ptr<double>(new double[2*Np],DeleteArray<double>);
        FILE *File;
        File=fopen(LocalRestartFile,"rb");
        fread(cfq.get(),sizeof(double),19*Np,File);
        fread(cDen.get(),sizeof(double),2*Np,File);
        fclose(File);

		// Copy the restart data to the GPU
		ScaLBL_CopyToDevice(fq,cfq.get(),19*Np*sizeof(double));
		ScaLBL_CopyToDevice(Den,cDen.get(),2*Np*sizeof(double));
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);

        //ScaLBL_D3Q7_GreyColorIMRT_Init(Den, Aq, Bq, Phi, 0, ScaLBL_Comm->LastExterior(), Np);//initialize D3Q7 density components
        //ScaLBL_D3Q7_GreyColorIMRT_Init(Den, Aq, Bq, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_D3Q19_GreyscaleColor_Laplacian(NeighborList, Phi, PhiLap, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_D3Q19_GreyscaleColor_Laplacian(NeighborList, Phi, PhiLap, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_D3Q7_GreyColorIMRT_Init(Den, Cq, PhiLap, gamma,kappaA,kappaB,lambdaA,lambdaB, 0, ScaLBL_Comm->LastExterior(), Np);//initialize D3Q7 density components
        ScaLBL_D3Q7_GreyColorIMRT_Init(Den, Cq, PhiLap, gamma,kappaA,kappaB,lambdaA,lambdaB, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);

        //TODO need to initialize velocity field !
        //this is required for calculating the pressure_dvc
        //can make a funciton to update velocity, such as ScaLBL_D3Q19_GreyColorIMRT_Velocity
	}
    else{
        if (rank==0)	printf ("Initializing density field \n");
        Density_and_Phase_Init();//initialize density field
        //ScaLBL_D3Q7_GreyColorIMRT_Init(Den, Aq, Bq, Phi, 0, ScaLBL_Comm->LastExterior(), Np);//initialize D3Q7 density components
        //ScaLBL_D3Q7_GreyColorIMRT_Init(Den, Aq, Bq, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_D3Q19_GreyscaleColor_Laplacian(NeighborList, Phi, PhiLap, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_D3Q19_GreyscaleColor_Laplacian(NeighborList, Phi, PhiLap, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_D3Q7_GreyColorIMRT_Init(Den, Cq, PhiLap, gamma,kappaA,kappaB,lambdaA,lambdaB, 0, ScaLBL_Comm->LastExterior(), Np);//initialize D3Q7 density components
        ScaLBL_D3Q7_GreyColorIMRT_Init(Den, Cq, PhiLap, gamma,kappaA,kappaB,lambdaA,lambdaB, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        
        if (rank==0)	printf ("Initializing distributions \n");
        ScaLBL_D3Q19_GreyColorIMRT_Init(fq, Den, rhoA, rhoB, Np);

        //Velocity also needs initialization
        if (rank==0)	printf ("Initializing velocity field \n");
        double *vel_init;
        vel_init = new double [3*Np];
        for (int i=0;i<3*Np;i++) vel_init[i]=0.0;
		ScaLBL_CopyToDevice(Velocity,vel_init,3*Np*sizeof(double));
		ScaLBL_DeviceBarrier();
        delete [] vel_init;
    }
}

void ScaLBL_GreyscaleColorModel::Run(){
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
	if (greyscaleColor_db->keyExists( "timestep" )){
		timestep = greyscaleColor_db->getScalar<int>( "timestep" );
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
	double error = 1.0;
	double flow_rate_previous = 0.0;
	while (timestep < timestepMax && error > tolerance) {
		//************************************************************************/
		// *************ODD TIMESTEP*************//
		timestep++;
		// Compute the density field
		// Read for Aq, Bq happens in this routine (requires communication)
		//ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
		ScaLBL_Comm->SendD3Q7AA(Cq); //READ FROM NORMAL
		//ScaLBL_D3Q7_AAodd_GreyscaleColorDensity(NeighborList, Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_D3Q7_AAodd_GreyscaleColorPhi(NeighborList, Cq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		//ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
		ScaLBL_Comm->RecvD3Q7AA(Cq); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();
		//ScaLBL_D3Q7_AAodd_GreyscaleColorDensity(NeighborList, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_D3Q7_AAodd_GreyscaleColorPhi(NeighborList, Cq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);

        // Update local pressure
        ScaLBL_D3Q19_GreyscaleColor_Pressure(fq, Den, Porosity, Velocity, Pressure_dvc, rhoA, rhoB, Np);
        // Compute pressure gradient
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, Pressure_dvc, PressureGrad, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(Pressure_dvc);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, Pressure_dvc, PressureGrad, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(Pressure_dvc,PressureGrad);
		ScaLBL_DeviceBarrier();
        // Compute Pressure Tensor
		ScaLBL_Comm->SendHalo(Phi);
        ScaLBL_D3Q19_GreyscaleColor_PressureTensor(NeighborList,Phi,PressTensor,PhiLap,kappaA,kappaB,lambdaA,lambdaB,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(),Np);
		ScaLBL_Comm->RecvHalo(Phi);
		ScaLBL_DeviceBarrier();
        ScaLBL_D3Q19_GreyscaleColor_PressureTensor(NeighborList,Phi,PressTensor,PhiLap,kappaA,kappaB,lambdaA,lambdaB,0,ScaLBL_Comm->LastExterior(),Np);
        /* Compute gradient of the pressure tensor */
		// call the recv Grad function once per tensor element
        // 1st tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[0*Np], &PressTensorGrad[0*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[0*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[0*Np], &PressTensorGrad[0*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[0*Np],&PressTensorGrad[0*Np]);
		// 2nd tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[1*Np], &PressTensorGrad[3*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[1*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[1*Np], &PressTensorGrad[3*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[1*Np],&PressTensorGrad[3*Np]);
		// 3rd tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[2*Np], &PressTensorGrad[6*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[2*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[2*Np], &PressTensorGrad[6*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[2*Np],&PressTensorGrad[6*Np]);
		// 4th tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[3*Np], &PressTensorGrad[9*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[3*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[3*Np], &PressTensorGrad[9*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[3*Np],&PressTensorGrad[9*Np]);
		// 5th tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[4*Np], &PressTensorGrad[12*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[4*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[4*Np], &PressTensorGrad[12*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[4*Np],&PressTensorGrad[12*Np]);
		// 6th tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[5*Np], &PressTensorGrad[15*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[5*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[5*Np], &PressTensorGrad[15*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[5*Np],&PressTensorGrad[15*Np]);

//        //compute Den gradients - component A
//		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &Den[0*Np], DenGradA, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
//		ScaLBL_Comm->SendHalo(&Den[0*Np]);
//		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &Den[0*Np], DenGradA, 0, ScaLBL_Comm->LastExterior(), Np);
//		ScaLBL_Comm->RecvGrad(&Den[0*Np],DenGradA);
//        //compute Den gradients - component B
//		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &Den[1*Np], DenGradB, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
//		ScaLBL_Comm->SendHalo(&Den[1*Np]);
//		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &Den[1*Np], DenGradB, 0, ScaLBL_Comm->LastExterior(), Np);
//		ScaLBL_Comm->RecvGrad(&Den[1*Np],DenGradB);


		ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
//        ScaLBL_D3Q19_AAodd_GreyscaleColor(NeighborList, fq, Aq, Bq, Den, DenGradA, DenGradB, SolidForce, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, 
//                tauA,tauB,tauA_eff,tauB_eff,rhoA,rhoB,Gsc,Fx,Fy,Fz,Porosity,Permeability,Velocity,Pressure_dvc);
        ScaLBL_D3Q19_AAodd_GreyscaleColorChem(NeighborList, fq, Cq, Phi, Den, SolidForce, 
                                              ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np,
                                              tauA, tauB, tauA_eff, tauB_eff, rhoA, rhoB, gamma,kappaA,kappaB,lambdaA,lambdaB, Fx, Fy, Fz, 
                                              Porosity, Permeability, Velocity, Pressure_dvc,PressureGrad,PressTensorGrad,PhiLap);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();

//		// Set BCs
//		if (BoundaryCondition == 3){
//			ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
//			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
//		}
//        ScaLBL_D3Q19_AAodd_GreyscaleColor(NeighborList, fq, Aq, Bq, Den, DenGradA, DenGradB, SolidForce, 0, ScaLBL_Comm->LastExterior(), Np, 
//                tauA,tauB,tauA_eff,tauB_eff,rhoA,rhoB,Gsc,Fx,Fy,Fz,Porosity,Permeability,Velocity,Pressure_dvc);
        ScaLBL_D3Q19_AAodd_GreyscaleColorChem(NeighborList, fq, Cq, Phi, Den, SolidForce, 
                                              0, ScaLBL_Comm->LastExterior(), Np,
                                              tauA, tauB, tauA_eff, tauB_eff, rhoA, rhoB, gamma,kappaA,kappaB,lambdaA,lambdaB, Fx, Fy, Fz, 
                                              Porosity, Permeability, Velocity, Pressure_dvc,PressureGrad,PressTensorGrad,PhiLap);
		ScaLBL_DeviceBarrier(); 
        MPI_Barrier(comm);

		// *************EVEN TIMESTEP*************//
		timestep++;
		// Compute the density field
		// Read for Aq, Bq happens in this routine (requires communication)
		//ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
		ScaLBL_Comm->SendD3Q7AA(Cq); //READ FROM NORMAL
		//ScaLBL_D3Q7_AAeven_GreyscaleColorDensity(Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_D3Q7_AAeven_GreyscaleColorPhi(Cq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		//ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
		ScaLBL_Comm->RecvD3Q7AA(Cq); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();
		//ScaLBL_D3Q7_AAeven_GreyscaleColorDensity(Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_D3Q7_AAeven_GreyscaleColorPhi(Cq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);

        // Update local pressure
        ScaLBL_D3Q19_GreyscaleColor_Pressure(fq, Den, Porosity, Velocity, Pressure_dvc, rhoA, rhoB, Np);
        // Compute pressure gradient
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, Pressure_dvc, PressureGrad, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(Pressure_dvc);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, Pressure_dvc, PressureGrad, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(Pressure_dvc,PressureGrad);
		ScaLBL_DeviceBarrier();
        // Compute Pressure Tensor
		ScaLBL_Comm->SendHalo(Phi);
        ScaLBL_D3Q19_GreyscaleColor_PressureTensor(NeighborList,Phi,PressTensor,PhiLap,kappaA,kappaB,lambdaA,lambdaB,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(),Np);
		ScaLBL_Comm->RecvHalo(Phi);
		ScaLBL_DeviceBarrier();
        ScaLBL_D3Q19_GreyscaleColor_PressureTensor(NeighborList,Phi,PressTensor,PhiLap,kappaA,kappaB,lambdaA,lambdaB,0,ScaLBL_Comm->LastExterior(),Np);
        /* Compute gradient of the pressure tensor */
		// call the recv Grad function once per tensor element
        // 1st tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[0*Np], &PressTensorGrad[0*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[0*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[0*Np], &PressTensorGrad[0*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[0*Np],&PressTensorGrad[0*Np]);
		// 2nd tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[1*Np], &PressTensorGrad[3*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[1*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[1*Np], &PressTensorGrad[3*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[1*Np],&PressTensorGrad[3*Np]);
		// 3rd tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[2*Np], &PressTensorGrad[6*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[2*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[2*Np], &PressTensorGrad[6*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[2*Np],&PressTensorGrad[6*Np]);
		// 4th tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[3*Np], &PressTensorGrad[9*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[3*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[3*Np], &PressTensorGrad[9*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[3*Np],&PressTensorGrad[9*Np]);
		// 5th tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[4*Np], &PressTensorGrad[12*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[4*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[4*Np], &PressTensorGrad[12*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[4*Np],&PressTensorGrad[12*Np]);
		// 6th tensor element
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[5*Np], &PressTensorGrad[15*Np], ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->SendHalo(&PressTensor[5*Np]);
		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &PressTensor[5*Np], &PressTensorGrad[15*Np], 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->RecvGrad(&PressTensor[5*Np],&PressTensorGrad[15*Np]);

//        //compute Den gradients - component A
//		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &Den[0*Np], DenGradA, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
//		ScaLBL_Comm->SendHalo(&Den[0*Np]);
//		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &Den[0*Np], DenGradA, 0, ScaLBL_Comm->LastExterior(), Np);
//		ScaLBL_Comm->RecvGrad(&Den[0*Np],DenGradA);
//        //compute Den gradients - component B
//		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &Den[1*Np], DenGradB, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
//		ScaLBL_Comm->SendHalo(&Den[1*Np]);
//		ScaLBL_D3Q19_GreyscaleColor_Gradient(NeighborList, &Den[1*Np], DenGradB, 0, ScaLBL_Comm->LastExterior(), Np);
//		ScaLBL_Comm->RecvGrad(&Den[1*Np],DenGradB);

		ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
//        ScaLBL_D3Q19_AAeven_GreyscaleColor(fq, Aq, Bq, Den, DenGradA, DenGradB, SolidForce, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, 
//                tauA,tauB,tauA_eff,tauB_eff,rhoA,rhoB,Gsc,Fx,Fy,Fz,Porosity,Permeability,Velocity,Pressure_dvc);
        ScaLBL_D3Q19_AAeven_GreyscaleColorChem(fq, Cq, Phi, Den, SolidForce, 
                                               ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np,
                                               tauA, tauB, tauA_eff, tauB_eff, rhoA, rhoB, gamma,kappaA,kappaB,lambdaA,lambdaB, Fx, Fy, Fz, 
                                               Porosity, Permeability, Velocity, Pressure_dvc,PressureGrad,PressTensorGrad,PhiLap);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();

//		// Set BCs
//		if (BoundaryCondition == 3){
//			ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
//			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
//		}
//        ScaLBL_D3Q19_AAeven_GreyscaleColor(fq, Aq, Bq, Den, DenGradA, DenGradB, SolidForce, 0, ScaLBL_Comm->LastExterior(), Np, 
//                tauA,tauB,tauA_eff,tauB_eff,rhoA,rhoB,Gsc,Fx,Fy,Fz,Porosity,Permeability,Velocity,Pressure_dvc);
        ScaLBL_D3Q19_AAeven_GreyscaleColorChem(fq, Cq, Phi, Den, SolidForce, 
                                               0, ScaLBL_Comm->LastExterior(), Np,
                                               tauA, tauB, tauA_eff, tauB_eff, rhoA, rhoB, gamma,kappaA,kappaB,lambdaA,lambdaB, Fx, Fy, Fz, 
                                               Porosity, Permeability, Velocity, Pressure_dvc,PressureGrad,PressTensorGrad,PhiLap);
        ScaLBL_DeviceBarrier(); 
        MPI_Barrier(comm);
		//************************************************************************/
		
//		if (timestep%analysis_interval==0){
//			ScaLBL_Comm->RegularLayout(Map,&Velocity[0],Velocity_x);
//			ScaLBL_Comm->RegularLayout(Map,&Velocity[Np],Velocity_y);
//			ScaLBL_Comm->RegularLayout(Map,&Velocity[2*Np],Velocity_z);
//			//ScaLBL_Comm->RegularLayout(Map,Porosity,PorosityMap);
//			//ScaLBL_Comm->RegularLayout(Map,Pressure_dvc,Pressure);
//			
//			double count_loc=0;
//			double count;
//			double vax,vay,vaz;
//			double vax_loc,vay_loc,vaz_loc;
//            //double px_loc,py_loc,pz_loc;
//            //double px,py,pz;
//            //double mass_loc,mass_glb;
//            
//            //parameters for domain average
//	        int64_t i,j,k,n,imin,jmin,kmin,kmax;
//            // If external boundary conditions are set, do not average over the inlet and outlet
//            kmin=1; kmax=Nz-1;
//            //In case user forgets to specify the inlet/outlet buffer layers for BC>0
//            if (BoundaryCondition > 0 && Dm->kproc() == 0) kmin=4;
//            if (BoundaryCondition > 0 && Dm->kproc() == Dm->nprocz()-1) kmax=Nz-4;
//
//            imin=jmin=1;
//            // If inlet/outlet layers exist use these as default
//            //if (Dm->inlet_layers_x > 0) imin = Dm->inlet_layers_x;
//            //if (Dm->inlet_layers_y > 0) jmin = Dm->inlet_layers_y;
//            if (BoundaryCondition > 0 && Dm->inlet_layers_z > 0 && Dm->kproc() == 0) kmin = 1 + Dm->inlet_layers_z;//"1" indicates the halo layer
//            if (BoundaryCondition > 0 && Dm->outlet_layers_z > 0 && Dm->kproc() == Dm->nprocz()-1) kmax = Nz-1 - Dm->outlet_layers_z; 
//
////			px_loc = py_loc = pz_loc = 0.f;
////            mass_loc = 0.f;
////			for (int k=kmin; k<kmax; k++){
////				for (int j=jmin; j<Ny-1; j++){
////					for (int i=imin; i<Nx-1; i++){
////						if (SignDist(i,j,k) > 0){
////							px_loc   += Velocity_x(i,j,k)*Den*PorosityMap(i,j,k);
////							py_loc   += Velocity_y(i,j,k)*Den*PorosityMap(i,j,k);
////							pz_loc   += Velocity_z(i,j,k)*Den*PorosityMap(i,j,k);
////							mass_loc += Den*PorosityMap(i,j,k);
////						}
////					}
////				}
////			}
////			MPI_Allreduce(&px_loc,  &px,      1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
////			MPI_Allreduce(&py_loc,  &py,      1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
////			MPI_Allreduce(&pz_loc,  &pz,      1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
////			MPI_Allreduce(&mass_loc,&mass_glb,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
////			
////			vax = px/mass_glb;
////			vay = py/mass_glb;
////			vaz = pz/mass_glb;
//            
//			vax_loc = vay_loc = vaz_loc = 0.f;
//			for (int k=kmin; k<kmax; k++){
//				for (int j=jmin; j<Ny-1; j++){
//					for (int i=imin; i<Nx-1; i++){
//						if (SignDist(i,j,k) > 0){
//							vax_loc += Velocity_x(i,j,k);
//							vay_loc += Velocity_y(i,j,k);
//							vaz_loc += Velocity_z(i,j,k);
//							count_loc+=1.0;
//						}
//					}
//				}
//			}
//			//MPI_Allreduce(&vax_loc,&vax,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
//			//MPI_Allreduce(&vay_loc,&vay,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
//			//MPI_Allreduce(&vaz_loc,&vaz,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
//			//MPI_Allreduce(&count_loc,&count,1,MPI_DOUBLE,MPI_SUM,Mask->Comm);
//			
//            vax = Mask->Comm.sumReduce( vax_loc );
//            vay = Mask->Comm.sumReduce( vay_loc );
//            vaz = Mask->Comm.sumReduce( vaz_loc );
//            count = Mask->Comm.sumReduce( count_loc );
//
//			vax /= count;
//			vay /= count;
//			vaz /= count;
//
//			double force_mag = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
//			double dir_x = Fx/force_mag;
//			double dir_y = Fy/force_mag;
//			double dir_z = Fz/force_mag;
//			if (force_mag == 0.0){
//				// default to z direction
//				dir_x = 0.0;
//				dir_y = 0.0;
//				dir_z = 1.0;
//				force_mag = 1.0;
//			}
//			//double flow_rate = (px*dir_x + py*dir_y + pz*dir_z)/mass_glb;
//			double flow_rate = (vax*dir_x + vay*dir_y + vaz*dir_z);
//			
//			error = fabs(flow_rate - flow_rate_previous) / fabs(flow_rate);
//			flow_rate_previous = flow_rate;
//			
//			//if (rank==0) printf("Computing Minkowski functionals \n");
//			Morphology.ComputeScalar(SignDist,0.f);
//			//Morphology.PrintAll();
//			double mu = (tau-0.5)/3.f;
//			double Vs = Morphology.V();
//			double As = Morphology.A();
//			double Hs = Morphology.H();
//			double Xs = Morphology.X();
//			Vs = Dm->Comm.sumReduce( Vs);
//			As = Dm->Comm.sumReduce( As);
//			Hs = Dm->Comm.sumReduce( Hs);
//			Xs = Dm->Comm.sumReduce( Xs);
//
//			double h = Dm->voxel_length;
//			//double absperm = h*h*mu*Mask->Porosity()*flow_rate / force_mag;
//			double absperm = h*h*mu*GreyPorosity*flow_rate / force_mag;
//
//            if (rank==0){
//				printf("     AbsPerm = %.5g [micron^2]\n",absperm);
//                bool WriteHeader=false;
//                FILE * log_file = fopen("Permeability.csv","r");
//                if (log_file != NULL)
//                    fclose(log_file);
//                else
//                    WriteHeader=true;
//                log_file = fopen("Permeability.csv","a");
//                if (WriteHeader)
//                    fprintf(log_file,"timestep Fx Fy Fz mu Vs As Hs Xs vax vay vaz AbsPerm \n",
//                            timestep,Fx,Fy,Fz,mu,h*h*h*Vs,h*h*As,h*Hs,Xs,vax,vay,vaz,absperm);
//
//                fprintf(log_file,"%i %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",timestep, Fx, Fy, Fz, mu, 
//                        h*h*h*Vs,h*h*As,h*Hs,Xs,vax,vay,vaz, absperm);
//                fclose(log_file);
//            }
//		}

		if (timestep%visualization_interval==0){
            VelocityField();
        }

		if (timestep%restart_interval==0){
            //Use rank=0 write out Restart.db
            if (rank==0) {
                greyscaleColor_db->putScalar<int>("timestep",timestep);    		
                greyscaleColor_db->putScalar<bool>( "Restart", true );
                current_db->putDatabase("GreyscaleColor", greyscaleColor_db);
                std::ofstream OutStream("Restart.db");
                current_db->print(OutStream, "");
                OutStream.close();
      
            }
            //Write out Restart data.
            std::shared_ptr<double> cfq;
            cfq = std::shared_ptr<double>(new double[19*Np],DeleteArray<double>);
            ScaLBL_CopyToHost(cfq.get(),fq,19*Np*sizeof(double));// Copy restart data to the CPU
            std::shared_ptr<double> cDen;
            cDen = std::shared_ptr<double>(new double[2*Np],DeleteArray<double>);
            ScaLBL_CopyToHost(cDen.get(),Den,2*Np*sizeof(double));// Copy restart data to the CPU

            FILE *RESTARTFILE;
            RESTARTFILE=fopen(LocalRestartFile,"wb");
            fwrite(cfq.get(),sizeof(double),19*Np,RESTARTFILE);
            fwrite(cDen.get(),sizeof(double),2*Np,RESTARTFILE);
            fclose(RESTARTFILE);
		    MPI_Barrier(comm);
        }
	}

	PROFILE_STOP("Loop");
	PROFILE_SAVE("lbpm_greyscaleColor_simulator",1);
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

void ScaLBL_GreyscaleColorModel::VelocityField(){

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

void ScaLBL_GreyscaleColorModel::WriteDebug(){
	// Copy back final phase indicator field and convert to regular layout
	DoubleArray PhaseField(Nx,Ny,Nz);

	//ScaLBL_CopyToHost(Porosity.data(), Poros, sizeof(double)*N);

//	FILE *OUTFILE;
//	sprintf(LocalRankFilename,"Phase.%05i.raw",rank);
//	OUTFILE = fopen(LocalRankFilename,"wb");
//	fwrite(PhaseField.data(),8,N,OUTFILE);
//	fclose(OUTFILE);
//
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

	ScaLBL_Comm->RegularLayout(Map,Pressure_dvc,PhaseField);
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

//	ScaLBL_Comm->RegularLayout(Map,&Porosity[0],PhaseField);
//	FILE *POROS_FILE;
//	sprintf(LocalRankFilename,"Porosity.%05i.raw",rank);
//	POROS_FILE = fopen(LocalRankFilename,"wb");
//	fwrite(PhaseField.data(),8,N,POROS_FILE);
//	fclose(POROS_FILE);
//
//	ScaLBL_Comm->RegularLayout(Map,&Permeability[0],PhaseField);
//	FILE *PERM_FILE;
//	sprintf(LocalRankFilename,"Permeability.%05i.raw",rank);
//	PERM_FILE = fopen(LocalRankFilename,"wb");
//	fwrite(PhaseField.data(),8,N,PERM_FILE);
//	fclose(PERM_FILE);
}
