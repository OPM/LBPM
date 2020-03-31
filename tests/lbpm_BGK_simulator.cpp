#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "analysis/TwoPhase.h"
#include "common/MPI_Helpers.h"

//#define WRITE_SURFACES

/*
 * Simulator for two-phase flow in porous media
 * James E. McClure 2013-2014
 */

using namespace std;


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
	{
		// parallel domain size (# of sub-domains)
		int nprocx,nprocy,nprocz;

		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Single Phase Permeability Calculation \n");
			printf("********************************************************\n");
		}

		// Variables that specify the computational domain  
		string FILENAME;
		int Nx,Ny,Nz;		// local sub-domain size
		int nspheres;		// number of spheres in the packing
		double Lx,Ly,Lz;	// Domain length
		double D = 1.0;		// reference length for non-dimensionalization
		// Color Model parameters
		int timestepMax, interval;
		double tau,Fx,Fy,Fz,tol,err;
		double din,dout;
		bool pBC,Restart;
		int i,j,k,n;

		int RESTART_INTERVAL=20000;

		if (rank==0){
			//.............................................................
			//		READ SIMULATION PARMAETERS FROM INPUT FILE
			//.............................................................
			ifstream input("Permeability.in");
			// Line 1: model parameters (tau, alpha, beta, das, dbs)
			input >> tau;			// Viscosity parameter
			// Line 2: External force components (Fx,Fy, Fz)
			input >> Fx;
			input >> Fy;
			input >> Fz;
			// Line 3: Pressure Boundary conditions
			input >> Restart;
			input >> pBC;
			input >> din;
			input >> dout;
			// Line 4: time-stepping criteria
			input >> timestepMax;		// max no. of timesteps
			input >> interval;			// restart interval
			input >> tol;				// error tolerance
			//.............................................................

			//.......................................................................
			// Reading the domain information file
			//.......................................................................
			ifstream domain("Domain.in");
			domain >> nprocx;
			domain >> nprocy;
			domain >> nprocz;
			domain >> Nx;
			domain >> Ny;
			domain >> Nz;
			//domain >> nspheres;
			domain >> Lx;
			domain >> Ly;
			domain >> Lz;
			//.......................................................................

		}
		// **************************************************************
		// Broadcast simulation parameters from rank 0 to all other procs
		MPI_Barrier(comm);
		//.................................................
		MPI_Bcast(&tau,1,MPI_DOUBLE,0,comm);
		//MPI_Bcast(&pBC,1,MPI_LOGICAL,0,comm);
		//	MPI_Bcast(&Restart,1,MPI_LOGICAL,0,comm);
		MPI_Bcast(&din,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&dout,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Fx,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Fy,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Fz,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&timestepMax,1,MPI_INT,0,comm);
		MPI_Bcast(&interval,1,MPI_INT,0,comm);
		MPI_Bcast(&tol,1,MPI_DOUBLE,0,comm);
		// Computational domain
		MPI_Bcast(&Nx,1,MPI_INT,0,comm);
		MPI_Bcast(&Ny,1,MPI_INT,0,comm);
		MPI_Bcast(&Nz,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
		//MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
		MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
		//.................................................
		MPI_Barrier(comm);

		RESTART_INTERVAL=interval;
		// **************************************************************
		// **************************************************************
		double rlx = 1.f/tau;

		if (nprocs != nprocx*nprocy*nprocz){
			printf("nprocx =  %i \n",nprocx);
			printf("nprocy =  %i \n",nprocy);
			printf("nprocz =  %i \n",nprocz);
			INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
		}

		if (rank==0){
			printf("********************************************************\n");
			printf("tau = %f \n", tau);
			printf("Force(x) = %.5g \n", Fx);
			printf("Force(y) = %.5g \n", Fy);
			printf("Force(z) = %.5g \n", Fz);
			printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
			printf("Process grid = %i x %i x %i\n",nprocx,nprocy,nprocz);
			printf("********************************************************\n");
		}

		double viscosity=(tau-0.5)/3.0;
		// Initialized domain and averaging framework for Two-Phase Flow
		int BC=pBC;
		Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
		for (i=0; i<Dm.Nx*Dm.Ny*Dm.Nz; i++) Dm.id[i] = 1;
		Dm.CommInit();
		TwoPhase Averages(Dm);

		// Mask that excludes the solid phase
		Domain Mask(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
		MPI_Barrier(comm);

		Nx += 2;	Ny += 2;	Nz += 2;
		int N = Nx*Ny*Nz;

		//.......................................................................
		if (rank == 0)	printf("Read input media... \n");
		//.......................................................................

		//.......................................................................
		// Filenames used
		char LocalRankString[8];
		char LocalRankFilename[40];
		char LocalRestartFile[40];
		char tmpstr[10];
		sprintf(LocalRankString,"%05d",rank);
		sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
		sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);

		//	printf("Local File Name =  %s \n",LocalRankFilename);
		// .......... READ THE INPUT FILE .......................................
		//	char value;
		char *id;
		id = new char[N];
		double sum, sum_local;
		double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
		//if (BoundaryCondition > 0) iVol_global = 1.0/(1.0*(Nx-2)*nprocx*(Ny-2)*nprocy*((Nz-2)*nprocz-6));
		double porosity, pore_vol;
		//...........................................................................
		if (rank == 0) cout << "Reading in domain from signed distance function..." << endl;

		//.......................................................................
		// Read the signed distance
		sprintf(LocalRankString,"%05d",rank);
		sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
		ReadBinaryFile(LocalRankFilename, Averages.SDs.data(), N);
		MPI_Barrier(comm);
		if (rank == 0) cout << "Domain set." << endl;

		//.......................................................................
		// Assign the phase ID field based on the signed distance
		//.......................................................................

		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					int n = k*Nx*Ny+j*Nx+i;
					id[n] = 0;
				}
			}
		}
		sum=0.f;
		pore_vol = 0.0;
		for ( k=0;k<Nz;k++){
			for ( j=0;j<Ny;j++){
				for ( i=0;i<Nx;i++){
					int n = k*Nx*Ny+j*Nx+i;
					if (Averages.SDs(n) > 0.0){
						id[n] = 2;	
					}
					// compute the porosity (actual interface location used)
					if (Averages.SDs(n) > 0.0){
						sum++;	
					}
				}
			}
		}

		if (rank==0) printf("Initialize from segmented data: solid=0, NWP=1, WP=2 \n");
		sprintf(LocalRankFilename,"ID.%05i",rank);
		size_t readID;
		FILE *IDFILE = fopen(LocalRankFilename,"rb");
		if (IDFILE==NULL) ERROR("lbpm_permeability_simulator: Error opening file: ID.xxxxx");
		readID=fread(id,1,N,IDFILE);
		if (readID != size_t(N)) printf("lbpm_permeability_simulator: Error reading ID (rank=%i) \n",rank);
		fclose(IDFILE);
		
		//.......................................................................
		// Compute the media porosity, assign phase labels and solid composition
		//.......................................................................
		sum_local=0.0;
		int Np=0;  // number of local pore nodes
		//.......................................................................
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
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
		// don't perform computations at the eight corners
		id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
		id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;
		//.........................................................
		MPI_Barrier(comm);

		// Initialize communication structures in averaging domain
		for (i=0; i<Mask.Nx*Mask.Ny*Mask.Nz; i++) Mask.id[i] = id[i];
		Mask.CommInit(comm);

		//...........................................................................
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");
		// Create a communicator for the device

		int Npad=(Np/16 + 2)*16;
		ScaLBL_Communicator ScaLBL_Comm(Mask);
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Npad];
		Np = ScaLBL_Comm.MemoryOptimizedLayoutAA(Map,neighborList,Mask.id,Np);
		MPI_Barrier(comm);
		
		// LBM variables
		if (rank==0)	printf ("Allocating distributions \n");
		//......................device distributions.................................
		int dist_mem_size = Np*sizeof(double);
		int neighborSize=18*(Np*sizeof(int));

		int *NeighborList;
		//		double *f_even,*f_odd;
		double * dist;
		double * Velocity;
		double * Pressure;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &dist, 19*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &Pressure, 3*sizeof(double)*Np);
		ScaLBL_CopyToDevice(NeighborList,     neighborList, neighborSize);
		//...........................................................................

		//...........................................................................
		if (rank==0)	printf("Setting the distributions, size = %i\n", N);
		//...........................................................................

		// Finalize setup for averaging domain
		//Averages.SetupCubes(Dm);
		Averages.UpdateSolid();
		// Initialize two phase flow variables (all wetting phase)
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n=k*Nx*Ny+j*Nx+i;
					Averages.Phase(i,j,k) = -1.0;
					Averages.SDn(i,j,k) = Averages.Phase(i,j,k);
					Averages.Phase_tplus(i,j,k) = Averages.SDn(i,j,k);
					Averages.Phase_tminus(i,j,k) = Averages.SDn(i,j,k);
					Averages.DelPhi(i,j,k) = 0.0;
					Averages.Press(i,j,k) = 0.0;
					Averages.Vel_x(i,j,k) = 0.0;
					Averages.Vel_y(i,j,k) = 0.0;
					Averages.Vel_z(i,j,k) = 0.0;
				}
			}
		}

		//.......................................................................

		ScaLBL_D3Q19_Init(dist, Np);

		int timestep = 0;
		if (rank==0) printf("********************************************************\n");
		if (rank==0)	printf("No. of timesteps: %i \n", timestepMax);

		//.......create and start timer............
		double starttime,stoptime,cputime;
		MPI_Barrier(comm);
		starttime = MPI_Wtime();
		//.........................................

		double D32,Fo,Re,velocity,err1D,mag_force,vel_prev;
		err = vel_prev = 1.0;
		if (rank==0) printf("Begin timesteps: error tolerance is %f \n", tol);
		//************ MAIN ITERATION LOOP ***************************************/
		while (timestep < timestepMax && err > tol ){

			timestep++;
			ScaLBL_Comm.SendD3Q19AA(dist); //READ FROM NORMAL
			ScaLBL_D3Q19_AAodd_BGK(NeighborList, dist, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np, rlx, Fx, Fy, Fz);
			ScaLBL_Comm.RecvD3Q19AA(dist); //WRITE INTO OPPOSITE
			ScaLBL_D3Q19_AAodd_BGK(NeighborList, dist, 0, ScaLBL_Comm.next, Np, rlx, Fx, Fy, Fz);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

			timestep++;
			ScaLBL_Comm.SendD3Q19AA(dist); //READ FORM NORMAL
			ScaLBL_D3Q19_AAeven_BGK(dist, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np, rlx, Fx, Fy, Fz);
			ScaLBL_Comm.RecvD3Q19AA(dist); //WRITE INTO OPPOSITE
			ScaLBL_D3Q19_AAeven_BGK(dist, 0, ScaLBL_Comm.next, Np, rlx, Fx, Fy, Fz);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
			//************************************************************************/

			if (timestep%500 == 0){
				//...........................................................................
				// Copy the data for for the analysis timestep
				//...........................................................................
				// Copy the phase from the GPU -> CPU
				//...........................................................................
				ScaLBL_DeviceBarrier();
		        ScaLBL_D3Q19_Pressure(dist,Pressure,Np);
		        ScaLBL_D3Q19_Momentum(dist,Velocity,Np);
		        
		        ScaLBL_Comm.RegularLayout(Map,Pressure,Averages.Press);
		        ScaLBL_Comm.RegularLayout(Map,&Velocity[0],Averages.Vel_x);
		        ScaLBL_Comm.RegularLayout(Map,&Velocity[Np],Averages.Vel_y);
		        ScaLBL_Comm.RegularLayout(Map,&Velocity[2*Np],Averages.Vel_z);

				// Way more work than necessary -- this is just to get the solid interfacial area!!
				Averages.Initialize();
				Averages.UpdateMeshValues();
				Averages.ComputeLocal();
				Averages.Reduce();

				double vawx = Averages.vaw_global(0);
				double vawy = Averages.vaw_global(1);
				double vawz = Averages.vaw_global(2);
				if (rank==0){
					// ************* DIMENSIONLESS FORCHEIMER EQUATION *************************
					//  Dye, A.L., McClure, J.E., Gray, W.G. and C.T. Miller
					//  Description of Non-Darcy Flows in Porous Medium Systems
					//  Physical Review E 87 (3), 033012
					//  Fo := density*D32^3*(density*force) / (viscosity^2)
					//	Re := density*D32*velocity / viscosity
					//  Fo = a*Re + b*Re^2
					// *************************************************************************
					//viscosity = (tau-0.5)*0.333333333333333333;
					D32 = 6.0*(Dm.Volume-Averages.vol_w_global)/Averages.As_global;
					printf("Sauter Mean Diameter = %f \n",D32);
					mag_force = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
					Fo = D32*D32*D32*mag_force/viscosity/viscosity;
					// .... 1-D flow should be aligned with force ...
					velocity = vawx*Fx/mag_force + vawy*Fy/mag_force + vawz*Fz/mag_force;
					err1D = fabs(velocity-sqrt(vawx*vawx+vawy*vawy+vawz*vawz))/velocity;
					//.......... Computation of the Reynolds number Re ..............
					Re = D32*velocity/viscosity;
					printf("Force: %.5g,%.5g,%.5g \n",Fx,Fy,Fz);
					printf("Velocity: %.5g,%.5g,%.5g \n",vawx,vawy,vawz);
					printf("Relative error for 1D representation: %.5g \n",err1D);
					printf("Dimensionless force: %5g \n", Fo);
					printf("Reynolds number: %.5g \n", Re);
					printf("Dimensionless Permeability (k/D^2): %.5g \n", Re/Fo);
				}
			}
		}
		//************************************************************************/
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

		NULL_USE(RESTART_INTERVAL);
	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
