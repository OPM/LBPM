/*
color lattice boltzmann model
 */
#include "models/ColorModel.h"

ScaLBL_ColorModel::ScaLBL_ColorModel(string filename){
	// read the input database 
    db = std::make_shared<Database>( filename );
    domain_db = db->getDatabase( "Domain" );
    color_db = db->getDatabase( "Color" );
    analysis_db = db->getDatabase( "Analysis" );
    
    // Color Model parameters
    timestepMax = domain_db->getScalar<int>( "timestepMax" );
    tauA = domain_db->getScalar<double>( "tauA" );
    tauB = domain_db->getScalar<double>( "tauB" );
    rhoA = domain_db->getScalar<double>( "rhoA" );
    rhoB = domain_db->getScalar<double>( "rhoB" );
    Fx = domain_db->getVector<double>( "F" )[0];
    Fy = domain_db->getVector<double>( "F" )[1];
    Fz = domain_db->getVector<double>( "F" )[2];
    alpha = domain_db->getScalar<double>( "alpha" );
    beta = domain_db->getScalar<double>( "beta" );
    Restart = domain_db->getScalar<int>( "Restart" );
    din = domain_db->getScalar<double>( "din" );
    dout = domain_db->getScalar<double>( "dout" );
    flux = domain_db->getScalar<double>( "flux" );;
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
    
    if (BoundaryCondition==4) flux = din*rhoA; // mass flux must adjust for density (see formulation for details)

    // Full domain used for analysis
    Domain Dm(domain_db);
    for (int i=0; i<Dm.Nx*Dm.Ny*Dm.Nz; i++) Dm.id[i] = 1;
    Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm) );
    //   TwoPhase Averages(Dm);
    Dm.CommInit(comm);

    // Mask that excludes the immobile phases
    Domain Mask(domain_db);
    MPI_Barrier(comm);
    
    Nx+=2; Ny+=2; Nz += 2;
    N = Nx*Ny*Nz;
    
}

ScaLBL_ColorModel::~ScaLBL_ColorModel(){
	
}

void ScaLBL_ColorModel::ReadInput(){
    //.......................................................................
    if (rank == 0)    printf("Read input media... \n");
    //.......................................................................
    sprintf(LocalRankString,"%05d",rank);
    sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
    sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
    // .......... READ THE INPUT FILE .......................................
    id = new char[N];
    double sum, sum_local;
    double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
    if (BoundaryCondition > 0) iVol_global = 1.0/(1.0*(Nx-2)*nprocx*(Ny-2)*nprocy*((Nz-2)*nprocz-6));
    //...........................................................................
    if (rank == 0) cout << "Reading in domain from signed distance function..." << endl;
    //.......................................................................
    sprintf(LocalRankString,"%05d",rank);
    sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
    ReadBinaryFile(LocalRankFilename, Averages->SDs.data(), N);
    MPI_Barrier(comm);
    if (rank == 0) cout << "Domain set." << endl;

    if (rank==0) printf("Initialize from segmented data: solid=0, NWP=1, WP=2 \n");
    sprintf(LocalRankFilename,"ID.%05i",rank);
    size_t readID;
    FILE *IDFILE = fopen(LocalRankFilename,"rb");
    if (IDFILE==NULL) ERROR("lbpm_color_simulator: Error opening file: ID.xxxxx");
    readID=fread(id,1,N,IDFILE);
    if (readID != size_t(N)) printf("lbpm_color_simulator: Error reading ID (rank=%i) \n",rank);
    fclose(IDFILE);

    // Read restart file
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
         FILE *RESTART = fopen(LocalRestartFile,"rb");
         if (IDFILE==NULL) ERROR("lbpm_color_simulator: Error opening file: Restart.xxxxx");
         readID=fread(id,1,N,RESTART);
         if (readID != size_t(N)) printf("lbpm_color_simulator: Error reading Restart (rank=%i) \n",rank);
         fclose(RESTART);

         MPI_Barrier(comm);
     }
}
void ScaLBL_ColorModel::Create(){
	/*
	 *  This function creates the variables needed to run a LBM 
	 */
     //.......................................................................
     // Compute the media porosity, assign phase labels and solid composition
     //.......................................................................
	double porosity;
	double sum_local=0.0;
     Np=0;  // number of local pore nodes
     //.......................................................................
     for (int k=1;k<Nz-1;k++){
         for (int j=1;j<Ny-1;j++){
             for (int i=1;i<Nx-1;i++){
                 int n = k*Nx*Ny+j*Nx+i;
                 if (id[n] > 0){
                     sum_local+=1.0;
                     Np++;
                 }
             }
         }
     }
     MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
     porosity = sum*iVol_global;
     if (rank==0) printf("Media porosity = %f \n",porosity);
     //.........................................................
     // If external boundary conditions are applied remove solid
     if (BoundaryCondition >  0  && Dm.kproc == 0){
         for (int k=0; k<3; k++){
             for (int j=0;j<Ny;j++){
                 for (int i=0;i<Nx;i++){
                     int n = k*Nx*Ny+j*Nx+i;
                     //id[n] = 1;
                     Averages->SDs(n) = max(Averages->SDs(n),1.0*(2.5-k));
                 }                    
             }
         }
     }
     if (BoundaryCondition >  0  && Dm.kproc == nprocz-1){
         for (int k=Nz-3; k<Nz; k++){
             for (int j=0;j<Ny;j++){
                 for (int i=0;i<Nx;i++){
                     int n = k*Nx*Ny+j*Nx+i;
                     //id[n] = 2;
                     Averages->SDs(n) = max(Averages->SDs(n),1.0*(k-Nz+2.5));
                 }                    
             }
         }
     }
     //.........................................................
     // don't perform computations at the eight corners
     id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
     id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;
     //.........................................................

     // Initialize communication structures in averaging domain
     for (int i=0; i<Mask.Nx*Mask.Ny*Mask.Nz; i++) Mask.id[i] = id[i];
     Mask.CommInit(comm);
     double *PhaseLabel;
     PhaseLabel = new double[N];
     Mask.AssignComponentLabels(PhaseLabel);
     
     //...........................................................................
        if (rank==0)    printf ("Create ScaLBL_Communicator \n");
        // Create a communicator for the device (will use optimized layout)
        ScaLBL_Communicator ScaLBL_Comm(Mask);
        //Create a second communicator based on the regular data layout
        //ScaLBL_Communicator ScaLBL_Comm_Regular(Mask);
        
        int Npad=(Np/16 + 2)*16;
        if (rank==0)    printf ("Set up memory efficient layout \n");
    	Map.resize(Nx,Ny,Nz);       Map.fill(0);
        auto neighborList= new int[18*Npad];
        Np = ScaLBL_Comm.MemoryOptimizedLayoutAA(Map,neighborList,Mask.id,Np);
        MPI_Barrier(comm);

        //...........................................................................
        //                MAIN  VARIABLES ALLOCATED HERE
        //...........................................................................
        // LBM variables
        if (rank==0)    printf ("Allocating distributions \n");
        //......................device distributions.................................
        int dist_mem_size = Np*sizeof(double);
        int neighborSize=18*(Np*sizeof(int));

        //...........................................................................
        ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
        ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);
        ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);
        ScaLBL_AllocateDeviceMemory((void **) &Aq, 7*dist_mem_size);
        ScaLBL_AllocateDeviceMemory((void **) &Bq, 7*dist_mem_size);
        ScaLBL_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
        ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*Np);        
        ScaLBL_AllocateDeviceMemory((void **) &Pressure, sizeof(double)*Np);
        ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
        ScaLBL_AllocateDeviceMemory((void **) &Gradient, 3*sizeof(double)*Np);
        ScaLBL_AllocateDeviceMemory((void **) &SolidPotential, 3*sizeof(double)*Np);
        
        //...........................................................................
        // Update GPU data structures
        if (rank==0)    printf ("Setting up device map and neighbor list \n");
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
        ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
        ScaLBL_DeviceBarrier();
        delete [] TmpMap;
}        
void ScaLBL_ColorModel::Initialize(){
	/*
	 * This function initializes model (includes both mobile and immobile components)
	 */

	// Compute the solid interaction potential and copy result to device
	if (rank==0) printf("Computing solid interaction potential \n");
	double *Tmp;
	Tmp=new double[3*Np];
	//Averages->UpdateMeshValues(); // this computes the gradient of distance field (among other things)
	// Create the distance stencil
	// Compute solid forces based on mean field approximation
	double *Dst;
	Dst = new double [5*5*5];
	for (int kk=0; kk<5; kk++){
		for (int jj=0; jj<5; jj++){
			for (int ii=0; ii<5; ii++){
				int index = kk*25+jj*5+ii;
				Dst[index] = sqrt(double(ii-2)*double(ii-2) + double(jj-2)*double(jj-2)+ double(kk-2)*double(kk-2));
			}
		}
	}
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				int idx=Map(i,j,k);
				if (!(idx < 0)){

					double phi_x = 0.f;
					double phi_y = 0.f;
					double phi_z = 0.f;
					for (int kk=0; kk<5; kk++){
						for (int jj=0; jj<5; jj++){
							for (int ii=0; ii<5; ii++){

								int index = kk*25+jj*5+ii;
								double distval= Dst[index];

								int idi=i+ii-2;
								int idj=j+jj-2;
								int idk=k+kk-2;

								if (idi < 0) idi=0;
								if (idj < 0) idj=0;
								if (idk < 0) idk=0;
								if (!(idi < Nx)) idi=Nx-1;
								if (!(idj < Ny)) idj=Ny-1;
								if (!(idk < Nz)) idk=Nz-1;

								int nn = idk*Nx*Ny + idj*Nx + idi;
								if (!(Mask.id[nn] > 0)){
									double vec_x = double(ii-2);
									double vec_y = double(jj-2);
									double vec_z = double(kk-2);

									double ALPHA=PhaseLabel[nn];
									double GAMMA=-2.f;
									if (distval > 2.f) ALPHA=0.f; // symmetric cutoff distance                                    
									phi_x += ALPHA*exp(GAMMA*distval)*vec_x/distval;
									phi_y += ALPHA*exp(GAMMA*distval)*vec_y/distval;
									phi_z += ALPHA*exp(GAMMA*distval)*vec_z/distval;
								}
							}
						}
					}
					Tmp[idx] = phi_x;
					Tmp[idx+Np] = phi_y;
					Tmp[idx+2*Np] = phi_z;

					/*                        double d = Averages->SDs(n);
                                         double dx = Averages->SDs_x(n);
                                         double dy = Averages->SDs_y(n);
                                         double dz = Averages->SDs_z(n);
                                         double value=cns*exp(-bns*fabs(d))-cws*exp(-bns*fabs(d));

                 Tmp[idx] = value*dx;
                 Tmp[idx+Np] = value*dy;
                 Tmp[idx+2*Np] = value*dz;
					 */
				}
			}
		}
	}
	ScaLBL_CopyToDevice(SolidPotential, Tmp, 3*sizeof(double)*Np);
	ScaLBL_DeviceBarrier();
	delete [] Tmp;
	delete [] Dst;

	DoubleArray Psx(Nx,Ny,Nz);
	DoubleArray Psy(Nx,Ny,Nz);
	DoubleArray Psz(Nx,Ny,Nz);
	DoubleArray Psnorm(Nx,Ny,Nz);
	ScaLBL_Comm.RegularLayout(Map,&SolidPotential[0],Psx);
	ScaLBL_Comm.RegularLayout(Map,&SolidPotential[Np],Psy);
	ScaLBL_Comm.RegularLayout(Map,&SolidPotential[2*Np],Psz);
	for (int n=0; n<N; n++) Psnorm(n) = Psx(n)*Psx(n)+Psy(n)*Psy(n)+Psz(n)*Psz(n);
	FILE *PFILE;
	sprintf(LocalRankFilename,"Potential.%05i.raw",rank);
	PFILE = fopen(LocalRankFilename,"wb");
	fwrite(Psnorm.data(),8,N,PFILE);
	fclose(PFILE);

	double count_wet=0.f;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				int idx=Map(i,j,k);
				int n = k*Nx*Ny+j*Nx+i;
				if (!(idx < 0)){
					if (Mask.id[n] == 1)
						PhaseLabel[idx] = 1.0;
					else {
						PhaseLabel[idx] = -1.0;
						count_wet+=1.f;
					}
				}
			}
		}
	}
	//printf("sw=%f \n",count_wet/double(Np));
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	// initialize phi based on PhaseLabel (include solid component labels)
	ScaLBL_CopyToDevice(Phi, PhaseLabel, Np*sizeof(double));
	//...........................................................................

    if (rank==0)    printf ("Initializing distributions \n");
    ScaLBL_D3Q19_Init(fq, Np);
    if (rank==0)    printf ("Initializing phase field \n");
    ScaLBL_DFH_Init(Phi, Den, Aq, Bq, 0, ScaLBL_Comm.last_interior, Np);
    
}

void ScaLBL_ColorModel::Run(){
    if (rank==0) printf("********************************************************\n");
    if (rank==0)    printf("No. of timesteps: %i \n", timestepMax);
    //.......create and start timer............
    double starttime,stoptime,cputime;
    ScaLBL_DeviceBarrier();
    MPI_Barrier(comm);
    starttime = MPI_Wtime();
    //.........................................
    //************ MAIN ITERATION LOOP ***************************************/
    PROFILE_START("Loop");
    runAnalysis analysis( analysis_db, rank_info, ScaLBL_Comm, Dm, Np, pBC, beta, Map );
    while (timestep < timestepMax ) {
        //if ( rank==0 ) { printf("Running timestep %i (%i MB)\n",timestep+1,(int)(Utilities::getMemoryUsage()/1048576)); }
        PROFILE_START("Update");
        // *************ODD TIMESTEP*************
        timestep++;
        // Compute the Phase indicator field
        // Read for Aq, Bq happens in this routine (requires communication)
        ScaLBL_Comm.BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAodd_DFH(NeighborList, Aq, Bq, Den, Phi, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);
        ScaLBL_Comm.BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
        ScaLBL_D3Q7_AAodd_DFH(NeighborList, Aq, Bq, Den, Phi, 0, ScaLBL_Comm.next, Np);
        
        // compute the gradient 
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);
        ScaLBL_Comm.SendHalo(Phi);
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, 0, ScaLBL_Comm.next, Np);
        ScaLBL_Comm.RecvGrad(Phi,Gradient);
        
        // Perform the collision operation
        ScaLBL_Comm.SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAodd_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                alpha, beta, Fx, Fy, Fz, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);
        ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        // Set BCs
        if (BoundaryCondition > 0){
            ScaLBL_Comm.Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
            ScaLBL_Comm.Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
        }
        if (BoundaryCondition == 3){
            ScaLBL_Comm.D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }
        if (BoundaryCondition == 4){
            din = ScaLBL_Comm.D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }
        ScaLBL_D3Q19_AAodd_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                alpha, beta, Fx, Fy, Fz, 0, ScaLBL_Comm.next, Np);
        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

        // *************EVEN TIMESTEP*************
        timestep++;
        // Compute the Phase indicator field
        ScaLBL_Comm.BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAeven_DFH(Aq, Bq, Den, Phi, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);
        ScaLBL_Comm.BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
        ScaLBL_D3Q7_AAeven_DFH(Aq, Bq, Den, Phi, 0, ScaLBL_Comm.next, Np);
        
        // compute the gradient 
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);
        ScaLBL_Comm.SendHalo(Phi);
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, 0, ScaLBL_Comm.next, Np);
        ScaLBL_Comm.RecvGrad(Phi,Gradient);

        // Perform the collision operation
        ScaLBL_Comm.SendD3Q19AA(fq); //READ FORM NORMAL
        ScaLBL_D3Q19_AAeven_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                alpha, beta, Fx, Fy, Fz, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);
        ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
        // Set boundary conditions
        if (BoundaryCondition > 0){
            ScaLBL_Comm.Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
            ScaLBL_Comm.Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
        }
        if (BoundaryCondition == 3){
            ScaLBL_Comm.D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }
        else if (BoundaryCondition == 4){
            din = ScaLBL_Comm.D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }
        ScaLBL_D3Q19_AAeven_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                alpha, beta, Fx, Fy, Fz,  0, ScaLBL_Comm.next, Np);
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
