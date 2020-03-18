
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/MPI_Helpers.h"

using namespace std;


//***************************************************************************************
int main(int argc, char **argv)
{
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	int check=0;
	{
		// parallel domain size (# of sub-domains)
		int i,j,k,n,Npad;
        auto filename = argv[1];
        auto db = std::make_shared<Database>( filename );
        auto domain_db = db->getDatabase( "Domain" );
        auto color_db = db->getDatabase( "Color" );
        auto analysis_db = db->getDatabase( "Analysis" );

        if (rank == 0){
            printf("********************************************************\n");
            printf("TestColorMassBounceback \n");
            printf("********************************************************\n");
        }
        // Initialize compute device
        //        int device=ScaLBL_SetDevice(rank);
        ScaLBL_DeviceBarrier();
        MPI_Barrier(comm);
        Utilities::setErrorHandlers();

        // Variables that specify the computational domain  
        string FILENAME;

        // Color Model parameters
        double tauA = color_db->getScalar<double>( "tauA" );
        double tauB = color_db->getScalar<double>( "tauB" );
        double rhoA = color_db->getScalar<double>( "rhoA" );
        double rhoB = color_db->getScalar<double>( "rhoB" );
        double Fx = color_db->getVector<double>( "F" )[0];
        double Fy = color_db->getVector<double>( "F" )[1];
        double Fz = color_db->getVector<double>( "F" )[2];
        double alpha = color_db->getScalar<double>( "alpha" );
        double beta = color_db->getScalar<double>( "beta" );
        bool Restart = color_db->getScalar<bool>( "Restart" );
        double din = color_db->getScalar<double>( "din" );

        // Read domain values
        auto size = domain_db->getVector<int>( "n" );
        auto nproc = domain_db->getVector<int>( "nproc" );
        int BoundaryCondition = domain_db->getScalar<int>( "BC" );
        int Nx = size[0];
        int Ny = size[1];
        int Nz = size[2];
        int nprocx = nproc[0];
        int nprocy = nproc[1];
        int nprocz = nproc[2];

        int timestep = 6;

        // Get the rank info
        const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);

        MPI_Barrier(comm);

        if (nprocs != nprocx*nprocy*nprocz){
            printf("nprocx =  %i \n",nprocx);
            printf("nprocy =  %i \n",nprocy);
            printf("nprocz =  %i \n",nprocz);
            INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
        }

        if (rank==0){
            printf("********************************************************\n");
            printf("tau (non-wetting) = %f \n", tauA);
            printf("tau (wetting) = %f \n", tauB);
            printf("density (non-wetting) = %f \n", rhoA);
            printf("density (wetting) = %f \n", rhoB);
            printf("alpha = %f \n", alpha);        
            printf("beta = %f \n", beta);
            printf("gamma_{wn} = %f \n", 5.796*alpha);
            printf("Force(x) = %f \n", Fx);
            printf("Force(y) = %f \n", Fy);
            printf("Force(z) = %f \n", Fz);
            printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
            printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
            if (BoundaryCondition==0) printf("Periodic boundary conditions will applied \n");
            if (BoundaryCondition==1) printf("Pressure boundary conditions will be applied \n");
            if (BoundaryCondition==2) printf("Velocity boundary conditions will be applied \n");
            if (BoundaryCondition==3) printf("Dynamic pressure boundary conditions will be applied \n");
            if (BoundaryCondition==4) printf("Average flux boundary conditions will be applied \n");
            if (!Restart) printf("Initial conditions assigned from phase ID file \n");
            if (Restart) printf("Initial conditions assigned from restart file \n");
            printf("********************************************************\n");
        }

        // Initialized domain and averaging framework for Two-Phase Flow
        bool pBC;
        if (BoundaryCondition==1 || BoundaryCondition==3 || BoundaryCondition == 4)
            pBC=true;
        else
            pBC=false;

        // Full domain used for averaging (do not use mask for analysis)
        std::shared_ptr<Domain> Dm(new Domain(domain_db,comm));
        for (int i=0; i<Dm->Nx*Dm->Ny*Dm->Nz; i++) Dm->id[i] = 1;
        Dm->CommInit();
        MPI_Barrier(comm);

        Nx+=2; Ny+=2; Nz += 2;
        int N = Nx*Ny*Nz;
        //.......................................................................
        if (rank == 0)    printf("Read input media... \n");
        //.......................................................................

		int Np=0;  // number of local pore nodes
		double *PhaseLabel;
		PhaseLabel = new double[N];
		//.......................................................................
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					Dm->id[n]=0;
				}
			}
		}
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					Dm->id[n]=1;
					Np++;
					// constant color
					PhaseLabel[n]= -1.0;
				}
			}
		}
		Dm->CommInit();
		MPI_Barrier(comm);
		if (rank == 0) cout << "Domain set." << endl;
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");

		//Create a second communicator based on the regular data layout
		std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm_Regular(new ScaLBL_Communicator(Dm));
		std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm(new ScaLBL_Communicator(Dm));

		// LBM variables
		if (rank==0)	printf ("Set up the neighborlist \n");

		int neighborSize=18*Np*sizeof(int);
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		Npad=Np+32;
		neighborList= new int[18*Npad];
		Np=ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Dm->id,Np);
		MPI_Barrier(comm);

		//......................device distributions.................................
		int dist_mem_size = Np*sizeof(double);
		if (rank==0)	printf ("Allocating distributions \n");
		int *NeighborList;
		int *dvcMap;
		double *fq, *Aq, *Bq;
		double *Den, *Phi;
		double *Gradient;
		double *SolidPotential;
		double *Vel;
		double *Pressure;

		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Aq, 7*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Bq, 7*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*Nx*Ny*Nz);		
		ScaLBL_AllocateDeviceMemory((void **) &Pressure, sizeof(double)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &Vel, 3*sizeof(double)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &Gradient, 3*sizeof(double)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &SolidPotential, 3*sizeof(double)*Np);
		
		//...........................................................................
		// Update GPU data structures
		if (rank==0)	printf ("Setting up device map and neighbor list \n");
		int *TmpMap;
		TmpMap=new int[Np*sizeof(int)];
		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int idx=Map(i,j,k);
					if (!(idx < 0))
						TmpMap[idx] = k*Nx*Ny+j*Nx+i;
				}
			}
		}
		ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
		ScaLBL_DeviceBarrier();
		delete [] TmpMap;		

		// copy the neighbor list 
		ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
		//...........................................................................

		// Distributions / densities for checking
		double nA,nB;
		double *DIST;
		DIST= new double [7*Np];
		double *DENSITY;
		DENSITY= new double [2*Np];
		int SIZE;
		
		int errc_odd_a=0;
		int errc_even_a=0;
		int errc_odd_b=0;
		int errc_even_b=0;
		
		//*******************Component A*******************
		// initialize phi based on PhaseLabel (include solid component labels)
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					// constant color
					PhaseLabel[n]= 1.0;
				}
			}
		}
		ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));

        if (rank==0)    printf ("Initializing distributions \n");
        ScaLBL_D3Q19_Init(fq, Np);
        if (rank==0)    printf ("Initializing phase field \n");
        ScaLBL_DFH_Init(Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastInterior(), Np);
		// *************ODD TIMESTEP*************
		// Compute the Phase indicator field
		// Read for Aq, Bq happens in this routine (requires communication)
        // Read for Aq, Bq happens in this routine (requires communication)
        ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAodd_DFH(NeighborList, Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
        ScaLBL_D3Q7_AAodd_DFH(NeighborList, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);
        
        // compute the gradient 
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->SendHalo(Phi);
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->RecvGrad(Phi,Gradient);
        
        // Perform the collision operation
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAodd_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                alpha, beta, Fx, Fy, Fz, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE

        ScaLBL_D3Q19_AAodd_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                alpha, beta, Fx, Fy, Fz, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

		timestep++;

		printf("Check after odd time \n");
		SIZE=2*Np*sizeof(double);
		ScaLBL_CopyToHost(&DENSITY[0],&Den[0],SIZE);

    	// Check the distributions
		SIZE=7*Np*sizeof(double);
		ScaLBL_CopyToHost(&DIST[0],&Aq[0],SIZE);

		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (Dm->id[n] > 0){
						int idx = Map(i,j,k);
    					nA=DENSITY[idx];
    					nB=DENSITY[Np+idx];
    					//printf("i,j,k=%i,%i,%i \n",i,j,k);
    					//printf("   nA=%f, nB=%f \n",nA,nB);
    					double val=DIST[idx];
    					double error = fabs(val - 0.3333333333333333*nA);
    					if (error > 1.0e-12) {
    						printf("   q=0, Aq=%f \n",val);
    						errc_odd_b++;
    					}
    					for (int q=1; q<7; q++){
    						val=DIST[q*Np+idx];
    						error = fabs(val - 0.1111111111111111*nA);
    						if (error > 1.0e-12) {
    							printf("   q=%i, Aq=%f \n",q,val);
    							errc_odd_b++;
    						}
    					}
					}
				}
			}
		}
		// EVEN TIMESTEP
        // Compute the Phase indicator field
         ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
         ScaLBL_D3Q7_AAeven_DFH(Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
         ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
         ScaLBL_D3Q7_AAeven_DFH(Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);
         
         // compute the gradient 
         ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
         ScaLBL_Comm->SendHalo(Phi);
         ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, 0, ScaLBL_Comm->LastExterior(), Np);
         ScaLBL_Comm->RecvGrad(Phi,Gradient);

         // Perform the collision operation
         ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
         ScaLBL_D3Q19_AAeven_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                 alpha, beta, Fx, Fy, Fz, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
         ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
         ScaLBL_D3Q19_AAeven_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                 alpha, beta, Fx, Fy, Fz,  0, ScaLBL_Comm->LastExterior(), Np);
         ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
         timestep++;
         //************************************************************************
		printf("Check after even time \n");

		SIZE=2*Np*sizeof(double);
		ScaLBL_CopyToHost(&DENSITY[0],&Den[0],SIZE);

    	// Check the distributions
		SIZE=7*Np*sizeof(double);
		ScaLBL_CopyToHost(&DIST[0],&Aq[0],SIZE);

		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (Dm->id[n] > 0){
						int idx = Map(i,j,k);
    					nA=DENSITY[idx];
    					nB=DENSITY[Np+idx];
    					//printf("i,j,k=%i,%i,%i \n",i,j,k);
    					//printf("   nA=%f, nB=%f \n",nA,nB);
    					double val=DIST[idx];
    					double error = fabs(val - 0.3333333333333333*nA);
    					if (error > 1.0e-12) {
    						printf("   q=0, Aq=%f \n",val);
    						errc_even_b++;
    					}
    					for (int q=1; q<7; q++){
    						val=DIST[q*Np+idx];
    						error = fabs(val - 0.1111111111111111*nA);
    						if (error > 1.0e-12) {
    							printf("   q=%i, Aq=%f \n",q,val);
    							errc_even_b++;
    						}
    					}
					}
				}
			}
		}

		//*******************Component B*******************
		// initialize phi based on PhaseLabel (include solid component labels)
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					// constant color
					PhaseLabel[n]= -1.0;
				}
			}
		}
		ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));


        if (rank==0)    printf ("Initializing distributions \n");
        ScaLBL_D3Q19_Init(fq, Np);
        if (rank==0)    printf ("Initializing phase field \n");
        ScaLBL_DFH_Init(Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastInterior(), Np);

		// *************ODD TIMESTEP*************
		// Compute the Phase indicator field
		// Read for Aq, Bq happens in this routine (requires communication)
        // Read for Aq, Bq happens in this routine (requires communication)
        ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
        ScaLBL_D3Q7_AAodd_DFH(NeighborList, Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
        ScaLBL_D3Q7_AAodd_DFH(NeighborList, Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);
        
        // compute the gradient 
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->SendHalo(Phi);
        ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_Comm->RecvGrad(Phi,Gradient);
        
        // Perform the collision operation
        ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
        ScaLBL_D3Q19_AAodd_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                alpha, beta, Fx, Fy, Fz, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
        ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE

        ScaLBL_D3Q19_AAodd_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                alpha, beta, Fx, Fy, Fz, 0, ScaLBL_Comm->LastExterior(), Np);
        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

		timestep++;

		printf("Check after odd time \n");
		SIZE=2*Np*sizeof(double);
		ScaLBL_CopyToHost(&DENSITY[0],&Den[0],SIZE);

    	// Check the distributions
		SIZE=7*Np*sizeof(double);
		ScaLBL_CopyToHost(&DIST[0],&Bq[0],SIZE);

		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (Dm->id[n] > 0){
						int idx = Map(i,j,k);
    					nA=DENSITY[idx];
    					nB=DENSITY[Np+idx];
    					//printf("i,j,k=%i,%i,%i \n",i,j,k);
    					//printf("   nA=%f, nB=%f \n",nA,nB);
    					double val=DIST[idx];
    					double error = fabs(val - 0.3333333333333333*nB);
    					if (error > 1.0e-12) {
    						printf("   q=0, Bq=%f \n",val);
    						errc_odd_b++;
    					}
    					for (int q=1; q<7; q++){
    						val=DIST[q*Np+idx];
    						error = fabs(val - 0.1111111111111111*nB);
    						if (error > 1.0e-12) {
    							printf("   q=%i, Bq=%f \n",q,val);
    							errc_odd_b++;
    						}
    					}
					}
				}
			}
		}

		// EVEN TIMESTEP
        // Compute the Phase indicator field
         ScaLBL_Comm->BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
         ScaLBL_D3Q7_AAeven_DFH(Aq, Bq, Den, Phi, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
         ScaLBL_Comm->BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
         ScaLBL_D3Q7_AAeven_DFH(Aq, Bq, Den, Phi, 0, ScaLBL_Comm->LastExterior(), Np);
         
         // compute the gradient 
         ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
         ScaLBL_Comm->SendHalo(Phi);
         ScaLBL_D3Q19_Gradient_DFH(NeighborList, Phi, Gradient, SolidPotential, 0, ScaLBL_Comm->LastExterior(), Np);
         ScaLBL_Comm->RecvGrad(Phi,Gradient);

         // Perform the collision operation
         ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
         ScaLBL_D3Q19_AAeven_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                 alpha, beta, Fx, Fy, Fz, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
         ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
         ScaLBL_D3Q19_AAeven_DFH(NeighborList, fq, Aq, Bq, Den, Phi, Gradient, rhoA, rhoB, tauA, tauB,
                 alpha, beta, Fx, Fy, Fz,  0, ScaLBL_Comm->LastExterior(), Np);
         ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
         timestep++;
         //************************************************************************
		printf("Check after even time \n");

		SIZE=2*Np*sizeof(double);
		ScaLBL_CopyToHost(&DENSITY[0],&Den[0],SIZE);

    	// Check the distributions
		SIZE=7*Np*sizeof(double);
		ScaLBL_CopyToHost(&DIST[0],&Bq[0],SIZE);

		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (Dm->id[n] > 0){
						int idx = Map(i,j,k);
    					nA=DENSITY[idx];
    					nB=DENSITY[Np+idx];
    					//printf("i,j,k=%i,%i,%i \n",i,j,k);
    					//printf("   nA=%f, nB=%f \n",nA,nB);
    					double val=DIST[idx];
    					double error = fabs(val - 0.3333333333333333*nB);
    					if (error > 1.0e-12) {
    						printf("   q=0, Bq=%f \n",val);
    						errc_even_b++;
    					}
    					for (int q=1; q<7; q++){
    						val=DIST[q*Np+idx];
    						error = fabs(val - 0.1111111111111111*nB);
    						if (error > 1.0e-12) {
    							printf("   q=%i, Bq=%f \n",q,val);
    							errc_even_b++;
    						}
    					}
					}
				}
			}
		}
		printf("Error counts: A even=%i, A odd=%i, B even=%i, B odd=%i \n",errc_even_a,errc_odd_a,errc_even_b,errc_odd_b);
		int errc_total=errc_even_a+errc_odd_a+errc_even_b+errc_odd_b;
		if (errc_total>0) check=1;
		else check=0;

	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
	return check;
}

