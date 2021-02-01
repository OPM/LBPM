/*
color lattice boltzmann model
 */
#include "models/FreeLeeModel.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"
#include "common/Communication.h"
#include "common/ReadMicroCT.h"
#include <stdlib.h>
#include <time.h>

ScaLBL_FreeLeeModel::ScaLBL_FreeLeeModel(int RANK, int NP, const Utilities::MPI& COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tauA(0),tauB(0),tauM(0),rhoA(0),rhoB(0),W(0),gamma(0),kappa(0),beta(0),
Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),inletA(0),inletB(0),outletA(0),outletB(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{
	
}
ScaLBL_FreeLeeModel::~ScaLBL_FreeLeeModel(){

}
void ScaLBL_FreeLeeModel::ReadParams(string filename){
	// read the input database 
	db = std::make_shared<Database>( filename );
	domain_db = db->getDatabase( "Domain" );
	freelee_db =  db->getDatabase( "FreeLee" );
	analysis_db = db->getDatabase( "Analysis" );
	vis_db = db->getDatabase( "Visualization" );

	// set defaults
	timestepMax = 100000;
	tauA = tauB = 1.0;
    tauM = 1.0;//relaxation time for phase field
	rhoA = rhoB = 1.0;
	Fx = Fy = Fz = 0.0;
	gamma=1e-3;//surface tension
	W=5.0;//interfacial thickness
    beta = 12.0*gamma/W;
    kappa = 3.0*gamma*W/2.0;//beta and kappa are related to surface tension \gamma
	Restart=false;
	din=dout=1.0;
	flux=0.0;
	
	// Color Model parameters
	if (freelee_db->keyExists( "timestepMax" )){
		timestepMax = freelee_db->getScalar<int>( "timestepMax" );
	}
	if (freelee_db->keyExists( "tauA" )){
		tauA = freelee_db->getScalar<double>( "tauA" );
	}
	if (freelee_db->keyExists( "tauB" )){
		tauB = freelee_db->getScalar<double>( "tauB" );
	}
	if (freelee_db->keyExists( "rhoA" )){
		rhoA = freelee_db->getScalar<double>( "rhoA" );
	}
	if (freelee_db->keyExists( "rhoB" )){
		rhoB = freelee_db->getScalar<double>( "rhoB" );
	}
	if (freelee_db->keyExists( "F" )){
		Fx = freelee_db->getVector<double>( "F" )[0];
		Fy = freelee_db->getVector<double>( "F" )[1];
		Fz = freelee_db->getVector<double>( "F" )[2];
	}
	if (freelee_db->keyExists( "gamma" )){
		gamma = freelee_db->getScalar<double>( "gamma" );
	}
	if (freelee_db->keyExists( "W" )){
		W = freelee_db->getScalar<double>( "W" );
	}
	if (freelee_db->keyExists( "Restart" )){
		Restart = freelee_db->getScalar<bool>( "Restart" );
	}
	if (freelee_db->keyExists( "din" )){
		din = freelee_db->getScalar<double>( "din" );
	}
	if (freelee_db->keyExists( "dout" )){
		dout = freelee_db->getScalar<double>( "dout" );
	}
	if (freelee_db->keyExists( "flux" )){
		flux = freelee_db->getScalar<double>( "flux" );
	}
	inletA=1.f;
	inletB=0.f;
	outletA=0.f;
	outletB=1.f;
    //update secondary parameters
    beta = 12.0*gamma/W;
    kappa = 3.0*gamma*W/2.0;//beta and kappa are related to surface tension \gamma
	//if (BoundaryCondition==4) flux *= rhoA; // mass flux must adjust for density (see formulation for details)

	BoundaryCondition = 0;
	if (domain_db->keyExists( "BC" )){
		BoundaryCondition = domain_db->getScalar<int>( "BC" );
	}
}

void ScaLBL_FreeLeeModel::SetDomain(){
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
	Nxh = Nx+2;
	Nyh = Ny+2;
	Nzh = Nz+2;
	Nh = Nxh*Nyh*Nzh;
	id = new signed char [N];
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;               // initialize this way

	comm.barrier();
	Dm->CommInit();
	comm.barrier();
	// Read domain parameters
	rank = Dm->rank();	
	nprocx = Dm->nprocx();
	nprocy = Dm->nprocy();
	nprocz = Dm->nprocz();
}

void ScaLBL_FreeLeeModel::ReadInput(){
	
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
	
	if (freelee_db->keyExists( "image_sequence" )){
		auto ImageList = freelee_db->getVector<std::string>( "image_sequence");
		int IMAGE_INDEX = freelee_db->getWithDefault<int>( "image_index", 0 );
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
	SignDist.resize(Nx,Ny,Nz);
	// Initialize the signed distance function
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				// Initialize distance to +/- 1
				SignDist(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
			}
		}
	}
	if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
	CalcDist(SignDist,id_solid,*Mask);
	
	if (rank == 0) cout << "Domain set." << endl;

}

void ScaLBL_FreeLeeModel::Create(){
	/*
	 *  This function creates the variables needed to run a LBM 
	 */
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
	ScaLBL_Comm_WideHalo  = std::shared_ptr<ScaLBLWideHalo_Communicator>(new ScaLBLWideHalo_Communicator(Mask,2));

	// create the layout for the LBM
	int Npad=(Np/16 + 2)*16;
	if (rank==0)    printf ("Set up memory efficient layout, %i | %i | %i \n", Np, Npad, N);
	Map.resize(Nx,Ny,Nz);       Map.fill(-2);
	auto neighborList= new int[18*Npad];
	Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Mask->id.data(),Np,2);
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
	ScaLBL_AllocateDeviceMemory((void **) &gqbar, 19*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &hq, 7*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &mu_phi, dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Den, dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*Nh);		
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
					TmpMap[idx] = ScaLBL_Comm_WideHalo->Map(i,j,k);
			}
		}
	}
	// check that TmpMap is valid
	for (int idx=0; idx<ScaLBL_Comm->LastExterior(); idx++){
		auto n = TmpMap[idx];
		if (n > Nxh*Nyh*Nzh){
			printf("Bad value! idx=%i \n", n);
			TmpMap[idx] = Nxh*Nyh*Nzh-1;
		}
	}
	for (int idx=ScaLBL_Comm->FirstInterior(); idx<ScaLBL_Comm->LastInterior(); idx++){
		auto n = TmpMap[idx];
		if ( n > Nxh*Nyh*Nzh ){
			printf("Bad value! idx=%i \n",n);
			TmpMap[idx] = Nxh*Nyh*Nzh-1;
		}
	}
    // copy the device map
	ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	comm.barrier();
	delete [] TmpMap;
	delete [] neighborList;
}        

void ScaLBL_FreeLeeModel::AssignComponentLabels_ChemPotential_ColorGrad()
{
	double *phase;
	phase = new double[Nh];

	size_t NLABELS=0;
	signed char VALUE=0;
	double AFFINITY=0.f;

	auto LabelList = greyscaleColor_db->getVector<int>( "ComponentLabels" );
	auto AffinityList = greyscaleColor_db->getVector<double>( "ComponentAffinity" );

	NLABELS=LabelList.size();
	if (NLABELS != AffinityList.size()){
		ERROR("Error: ComponentLabels and ComponentAffinity must be the same length! \n");
	}

	double label_count[NLABELS];
	double label_count_global[NLABELS];

	// Assign the labels
	for (size_t idx=0; idx<NLABELS; idx++) label_count[idx]=0;
	for (int k=0;k<Nzh;k++){
		for (int j=0;j<Nyh;j++){
			for (int i=0;i<Nxh;i++){

                //idx for double-halo array 'phase'
				int nh = k*Nxh*Nyh+j*Nxh+i;

                //idx for single-halo array Mask->id[n]
                int x=i-1;
                int y=j-1;
                int z=k-1;
				if (x<0)   x=0;
				if (y<0)   y=0;
				if (z<0)   z=0;
				if (x>=Nx) x=Nx-1;
				if (y>=Ny) y=Ny-1;
				if (z>=Nz) z=Nz-1;
				int n = z*Nx*Ny+y*Nx+x;
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
				phase[nh] = AFFINITY;
			}
		}
	}

	// Set Dm to match Mask
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i]; 
	
	for (size_t idx=0; idx<NLABELS; idx++)
		label_count_global[idx] = sumReduce( Dm->Comm, label_count[idx]);

	if (rank==0){
		printf("Number of component labels: %lu \n",NLABELS);
		for (unsigned int idx=0; idx<NLABELS; idx++){
			VALUE=LabelList[idx];
			AFFINITY=AffinityList[idx];
			double volume_fraction  = double(label_count_global[idx])/double((Nx-2)*(Ny-2)*(Nz-2)*nprocs);
			printf("   label=%d, affinity=%f, volume fraction==%f\n",VALUE,AFFINITY,volume_fraction); 
		}
	}

    //compute color gradient and laplacian of phase field
	double *ColorGrad_host, mu_phi_host;
	ColorGrad_host = new double[3*Np];
	mu_phi_host    = new double[Np];

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
	double w_face = 1.0/18.0;
	double w_edge = 1.0/36.0;
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

    double cs2_inv = 3.0;//inverse of c_s^2 for D3Q19 lattice
    int width = 2;//For better readability: make halo width explicity wherever possible
	for (int k=width; k<Nzh-width; k++){
		for (int j=width; j<Nyh-width; j++){
			for (int i=width; i<Nxh-width; i++){

                //idx for double-halo array 'phase'
				int nh = k*Nxh*Nyh+j*Nxh+i;
				
                int idx=Map(i-width+1,j-width+1,k-width+1);
				if (!(idx < 0)){
					double phi_x  = 0.f;
					double phi_y  = 0.f;
					double phi_z  = 0.f;
					double phi_Lap = 0.f;//Laplacian of the phase field
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
								if (!(idi < Nxh)) idi=Nxh-1;
								if (!(idj < Nyh)) idj=Nyh-1;
								if (!(idk < Nzh)) idk=Nzh-1;

								int nn = idk*Nxh*Nyh + idj*Nxh + idi;
								double vec_x = double(ii-1);
								double vec_y = double(jj-1);
								double vec_z = double(kk-1);
								double GWNS=phase[nn];
								double GWNS_local=phase[nh];
								phi_x  += GWNS*weight*vec_x;
								phi_y  += GWNS*weight*vec_y;
								phi_z  += GWNS*weight*vec_z;
                                phi_Lap += weight*(GWNS-GWNS_local);//Laplacian of the phase field
							}
						}
					}
                    //store color gradient
                    ColorGrad_host[idx+0*Np] = cs2_inv*phi_x;
                    ColorGrad_host[idx+1*Np] = cs2_inv*phi_y;
                    ColorGrad_host[idx+2*Np] = cs2_inv*phi_z;
                    //compute chemical potential
                    phi_Lap = 2.0*cs2_inv*phi_Lap;
                    mu_phi_host[idx] = 4.0*beta*phase[nh]*(phase[nh]+1.0)*(phase[nh]-1.0) - kappa*phi_Lap;
				}
			}
		}
	}

    //copy all data to device
	ScaLBL_CopyToDevice(Phi, phase, Nh*sizeof(double));
	ScaLBL_CopyToDevice(ColorGrad, ColorGrad_host, 3*Np*sizeof(double));
	ScaLBL_CopyToDevice(mu_phi, mu_phi_host, Np*sizeof(double));
	comm.barrier();
    delete [] phase;
    delete [] ColorGrad_host;
    delete [] mu_phi_host;
    delete [] Dst;
}

void ScaLBL_FreeLeeModel::Initialize(){
	/*
	 * This function initializes model
	 */
	if (rank==0)	printf ("Initializing phase field, chemical potential and color gradient\n");
	AssignComponentLabels_ChemPotential_ColorGrad();//initialize phase field Phi

	if (rank==0)	printf ("Initializing distributions for momentum transport\n");
	ScaLBL_D3Q19_FreeLeeModel_Init(gqbar, mu_phi, ColorGrad, Fx, Fy, Fz, Np);

	if (rank==0)	printf ("Initializing density field and distributions for phase-field transport\n");
	ScaLBL_FreeLeeModel_PhaseField_Init(dvcMap, Phi, Den, hq, ColorGrad, rhoA, rhoB, tauM, W, 0, ScaLBL_Comm->LastExterior(), Np);
	ScaLBL_FreeLeeModel_PhaseField_Init(dvcMap, Phi, Den, hq, ColorGrad, rhoA, rhoB, tauM, W, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);

	if (Restart == true){
        //TODO need to revise this function
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
        //ScaLBL_CopyToHost(cPhi, Phi, N*sizeof(double));
    	
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

        if (rank==0)	printf ("Initializing phase and density fields on device from Restart\n");
        ScaLBL_FreeLeeModel_PhaseField_InitFromRestart(Den, hq, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_FreeLeeModel_PhaseField_InitFromRestart(Den, hq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
	}

	// establish reservoirs for external bC
    // TODO to be revised
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
	//ScaLBL_CopyToHost(Averages->Phi.data(),Phi,N*sizeof(double));
}

void ScaLBL_FreeLeeModel::Run(){
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
	comm.barrier();
	starttime = MPI_Wtime();
	//.........................................

	//************ MAIN ITERATION LOOP ***************************************/
	PROFILE_START("Loop");
	while (timestep < timestepMax ) {
		//if ( rank==0 ) { printf("Running timestep %i (%i MB)\n",timestep+1,(int)(Utilities::getMemoryUsage()/1048576)); }
		PROFILE_START("Update");
		// *************ODD TIMESTEP*************
		timestep++;
        //-------------------------------------------------------------------------------------------------------------------
		// Compute the Phase indicator field
		// Read for hq happens in this routine (requires communication)
		ScaLBL_Comm->SendD3Q7AA(hq,0); //READ FROM NORMAL
		ScaLBL_D3Q7_AAodd_FreeLeeModel_PhaseField(NeighborList, dvcMap, hq, Den, Phi, rhoA, rhoB, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->RecvD3Q7AA(hq,0); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();
		ScaLBL_D3Q7_AAodd_FreeLeeModel_PhaseField(NeighborList, dvcMap, hq, Den, Phi, rhoA, rhoB,  0, ScaLBL_Comm->LastExterior(), Np);

		// Perform the collision operation
		ScaLBL_Comm->SendD3Q19AA(gqbar); //READ FROM NORMAL
		if (BoundaryCondition > 0 && BoundaryCondition < 5){
            //TODO to be revised
			ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
			ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
		}
		// Halo exchange for phase field
		ScaLBL_Comm_WideHalo->Send(Phi);

		ScaLBL_D3Q19_AAodd_FreeLeeModel(NeighborList, dvcMap, gqbar, hq, Den, Phi, mu_phi, Velocity, Pressure, rhoA, rhoB, tauA, tauB, tauM,
				                        kappa, beta, W, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm_WideHalo->Recv(Phi);
		ScaLBL_Comm->RecvD3Q19AA(gqbar); //WRITE INTO OPPOSITE
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
		else if (BoundaryCondition == 5){
			ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
			ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
		}
		ScaLBL_D3Q19_AAodd_FreeLeeModel(NeighborList, dvcMap, gqbar, hq, Den, Phi, mu_phi, Velocity, Pressure, rhoA, rhoB, tauA, tauB, tauM,
				                        kappa, beta, W, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_DeviceBarrier(); 
		MPI_Barrier(ScaLBL_Comm->MPI_COMM_SCALBL);

		// *************EVEN TIMESTEP*************
		timestep++;
		// Compute the Phase indicator field
		ScaLBL_Comm->SendD3Q7AA(hq,0); //READ FROM NORMAL
		ScaLBL_D3Q7_AAeven_FreeLeeModel_PhaseField(dvcMap, hq, Den, Phi, rhoA, rhoB,  ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm->RecvD3Q7AA(hq,0); //WRITE INTO OPPOSITE
		ScaLBL_DeviceBarrier();
		ScaLBL_D3Q7_AAeven_FreeLeeModel_PhaseField(dvcMap, hq, Den, Phi, rhoA, rhoB,  0, ScaLBL_Comm->LastExterior(), Np);

		// Perform the collision operation
		ScaLBL_Comm->SendD3Q19AA(gqbar); //READ FORM NORMAL
		// Halo exchange for phase field
		if (BoundaryCondition > 0 && BoundaryCondition < 5){
			ScaLBL_Comm->Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
			ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
		}
		ScaLBL_Comm_WideHalo->Send(Phi);
		ScaLBL_D3Q19_AAeven_FreeLeeModel(dvcMap, gqbar, hq, Den, Phi, mu_phi, Velocity, Pressure, rhoA, rhoB, tauA, tauB, tauM,
				                        kappa, beta, W, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
		ScaLBL_Comm_WideHalo->Recv(Phi);
		ScaLBL_Comm->RecvD3Q19AA(gqbar); //WRITE INTO OPPOSITE
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
		else if (BoundaryCondition == 5){
			ScaLBL_Comm->D3Q19_Reflection_BC_z(fq);
			ScaLBL_Comm->D3Q19_Reflection_BC_Z(fq);
		}
		ScaLBL_D3Q19_AAeven_FreeLeeModel(dvcMap, gqbar, hq, Den, Phi, mu_phi, Velocity, Pressure, rhoA, rhoB, tauA, tauB, tauM,
				                        kappa, beta, W, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np);
		ScaLBL_Comm->Barrier();
		//************************************************************************
		PROFILE_STOP("Update");
	}
	PROFILE_STOP("Loop");
	PROFILE_SAVE("lbpm_color_simulator",1);
	//************************************************************************
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


void ScaLBL_FreeLeeModel::WriteDebug(){
	// Copy back final phase indicator field and convert to regular layout
	DoubleArray PhaseData(Nxh,Nyh,Nzh);
	//ScaLBL_Comm->RegularLayout(Map,Phi,PhaseField);
	ScaLBL_CopyToHost(PhaseData.data(), Phi, sizeof(double)*Nh);

	FILE *OUTFILE;
	sprintf(LocalRankFilename,"Phase.%05i.raw",rank);
	OUTFILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseData.data(),8,Nh,OUTFILE);
	fclose(OUTFILE);

	DoubleArray PhaseField(Nx,Ny,Nz);
    ScaLBL_Comm->RegularLayout(Map,Den,PhaseField);
	FILE *AFILE;
	sprintf(LocalRankFilename,"Density.%05i.raw",rank);
	AFILE = fopen(LocalRankFilename,"wb");
	fwrite(PhaseField.data(),8,N,AFILE);
	fclose(AFILE);

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
