#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/TwoPhase.h"
#include "common/MPI_Helpers.h"

//#define WRITE_SURFACES

/*
 * Simulator for two-phase flow in porous media
 * James E. McClure 2013-2014
 */

using namespace std;

//*************************************************************************
// Steady State Single-Phase LBM to generate non-Darcy curves
//*************************************************************************

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
		int iproc,jproc,kproc;
		//*****************************************
		// MPI ranks for all 18 neighbors
		//**********************************
		int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
		int rank_xy,rank_XY,rank_xY,rank_Xy;
		int rank_xz,rank_XZ,rank_xZ,rank_Xz;
		int rank_yz,rank_YZ,rank_yZ,rank_Yz;
		//**********************************
		MPI_Request req1[18],req2[18];
		MPI_Status stat1[18],stat2[18];

		double REYNOLDS_NUMBER = 100.f;
		if (argc > 1){
			REYNOLDS_NUMBER=strtod(argv[1],NULL);
		}
		if (rank == 0){
			printf("********************************************************\n");
			printf("Simulating Single Phase Non-Darcy Curve, Re < %f \n",REYNOLDS_NUMBER);
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
			domain >> nspheres;
			domain >> Lx;
			domain >> Ly;
			domain >> Lz;
			//.......................................................................

			/*
			 * Set simulation parameters internally
			 */
			tau=1.f;
			Fx = 0.f;
			Fy = 0.f;
			Fz = 1.0e-7;
			pBC = 0;
			din = 1.0;
			dout = 1.0;
			timestepMax = nprocz*Nz*100;
			interval = 500;
			tol = 1.0e-4;


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
		MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
		MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
		//.................................................
		MPI_Barrier(comm);

		RESTART_INTERVAL=interval;
		// **************************************************************
		// **************************************************************
		double rlxA = 1.f/tau;
		double rlxB = 8.f*(2.f-rlxA)/(8.f-rlxA);


		if (nprocs != nprocx*nprocy*nprocz){
			printf("nprocx =  %i \n",nprocx);
			printf("nprocy =  %i \n",nprocy);
			printf("nprocz =  %i \n",nprocz);
			INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
		}

		if (rank==0){
			printf("********************************************************\n");
			printf("tau = %f \n", tau);
			printf("Force(x) = %f \n", Fx);
			printf("Force(y) = %f \n", Fy);
			printf("Force(z) = %f \n", Fz);
			printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
			printf("Process grid = %i x %i x %i\n",nprocx,nprocy,nprocz);
			printf("********************************************************\n");
		}

		double viscosity=(tau-0.5)/3.0;
		// Initialized domain and averaging framework for Two-Phase Flow
		int BC=pBC;
		Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
		TwoPhase Averages(Dm);

		InitializeRanks( rank, nprocx, nprocy, nprocz, iproc, jproc, kproc,
				rank_x, rank_y, rank_z, rank_X, rank_Y, rank_Z,
				rank_xy, rank_XY, rank_xY, rank_Xy, rank_xz, rank_XZ, rank_xZ, rank_Xz,
				rank_yz, rank_YZ, rank_yZ, rank_Yz );

		MPI_Barrier(comm);

		Nx += 2;	Ny += 2;	Nz += 2;

		int N = Nx*Ny*Nz;
		int dist_mem_size = N*sizeof(double);

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
		int sum = 0;
		double sum_local;
		double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
		if (pBC) iVol_global = 1.0/(1.0*(Nx-2)*nprocx*(Ny-2)*nprocy*((Nz-2)*nprocz-6));
		double porosity, pore_vol;
		//...........................................................................
		if (rank == 0) cout << "Reading in domain from signed distance function..." << endl;

		//.......................................................................
		sprintf(LocalRankString,"%05d",rank);
		//	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
		//	WriteLocalSolidID(LocalRankFilename, id, N);
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
					n = k*Nx*Ny+j*Nx+i;
					id[n] = 0;
				}
			}
		}
		sum=0;
		pore_vol = 0.0;
		for ( k=1;k<Nz-1;k++){
			for ( j=1;j<Ny-1;j++){
				for ( i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
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

		// Set up kstart, kfinish so that the reservoirs are excluded from averaging
		int kstart,kfinish;
		kstart = 1;
		kfinish = Nz-1;
		if (pBC && kproc==0)		kstart = 4;
		if (pBC && kproc==nprocz-1)	kfinish = Nz-4;

		// Compute the pore volume
		sum_local = 0.0;
		for ( k=kstart;k<kfinish;k++){
			for ( j=1;j<Ny-1;j++){
				for ( i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (id[n] > 0){
						sum_local += 1.0;
					}
				}
			}
		}
		MPI_Allreduce(&sum_local,&pore_vol,1,MPI_DOUBLE,MPI_SUM,comm);
		//	MPI_Allreduce(&sum_local,&porosity,1,MPI_DOUBLE,MPI_SUM,comm);
		porosity = pore_vol*iVol_global;
		if (rank==0) printf("Media porosity = %f \n",porosity);
		//.........................................................
		// If pressure boundary conditions are applied remove solid
		if (pBC && kproc == 0){
			for (k=0; k<3; k++){
				for (j=0;j<Ny;j++){
					for (i=0;i<Nx;i++){
						n = k*Nx*Ny+j*Nx+i;
						id[n] = 1;
						Averages.SDs(n) = max(Averages.SDs(n),1.0*(2.5-k));
					}
				}
			}
		}
		if (pBC && kproc == nprocz-1){
			for (k=Nz-3; k<Nz; k++){
				for (j=0;j<Ny;j++){
					for (i=0;i<Nx;i++){
						n = k*Nx*Ny+j*Nx+i;
						id[n] = 2;
						Averages.SDs(n) = max(Averages.SDs(n),1.0*(k-Nz+2.5));
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
		for (i=0; i<Dm.Nx*Dm.Ny*Dm.Nz; i++) Dm.id[i] = id[i];
		Dm.CommInit(comm);

		//...........................................................................
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");
		// Create a communicator for the device
		ScaLBL_Communicator ScaLBL_Comm(Dm);

		//...........device phase ID.................................................
		if (rank==0)	printf ("Copying phase ID to device \n");
		char *ID;
		AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
		// Copy to the device
		CopyToDevice(ID, id, N);
		//...........................................................................

		//...........................................................................
		//				MAIN  VARIABLES ALLOCATED HERE
		//...........................................................................
		// LBM variables
		if (rank==0)	printf ("Allocating distributions \n");
		//......................device distributions.................................
		double *f_even,*f_odd;
		//...........................................................................
		AllocateDeviceMemory((void **) &f_even, 10*dist_mem_size);	// Allocate device memory
		AllocateDeviceMemory((void **) &f_odd, 9*dist_mem_size);	// Allocate device memory
		//...........................................................................
		double *Velocity, *Pressure, *dvcSignDist;
		//...........................................................................
		AllocateDeviceMemory((void **) &Pressure, dist_mem_size);
		AllocateDeviceMemory((void **) &dvcSignDist, dist_mem_size);
		AllocateDeviceMemory((void **) &Velocity, 3*dist_mem_size);
		//...........................................................................

		// Copy signed distance for device initialization
		CopyToDevice(dvcSignDist, Averages.SDs.data(), dist_mem_size);
		//...........................................................................

		int logcount = 0; // number of surface write-outs

		//...........................................................................
		//				MAIN  VARIABLES INITIALIZED HERE
		//...........................................................................
		//...........................................................................
		if (rank==0)	printf("Setting the distributions, size = %i\n", N);
		//...........................................................................
		InitD3Q19(ID, f_even, f_odd, Nx, Ny, Nz);
		//......................................................................

		//.......................................................................
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

		if (rank==0 && pBC){
			printf("Setting inlet pressure = %f \n", din);
			printf("Setting outlet pressure = %f \n", dout);
		}
		if (pBC && kproc == 0)	{
			PressureBC_inlet(f_even,f_odd,din,Nx,Ny,Nz);
		}

		if (pBC && kproc == nprocz-1){
			PressureBC_outlet(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
		}

		int timestep = 0;
		if (rank==0) printf("********************************************************\n");
		if (rank==0)	printf("No. of timesteps: %i \n", timestepMax);

		//.......create and start timer............
		double starttime,stoptime,cputime;
		MPI_Barrier(comm);
		starttime = MPI_Wtime();
		//.........................................

		double D32,Fo,Re,velocity,err1D,mag_force,vel_prev;

		FILE * NONDARCY;
		if (rank == 0){
			NONDARCY = fopen("nondarcy.csv","a");
			fprintf(NONDARCY,"D32 Fx Fy Fz vx vy vz Re Fo\n");
		}

		Re = 0.f;
		// Generate a bunch of points until sufficiently high Re is obtained
		while (Re < REYNOLDS_NUMBER){
			// Increase the external force and simulate to steady state
			Fz = 2.0*Fz;

			err = vel_prev = 1.0;
			if (rank==0) printf("Begin timesteps: error tolerance is %f \n", tol);
			//************ MAIN ITERATION LOOP ***************************************/
			while (timestep < timestepMax && err > tol ){

				//*************************************************************************
				// Fused Color Gradient and Collision
				//*************************************************************************
				MRT( ID,f_even,f_odd,rlxA,rlxB,Fx,Fy,Fz,Nx,Ny,Nz);
				//*************************************************************************

				//*************************************************************************
				// Pack and send the D3Q19 distributions
				ScaLBL_Comm.SendD3Q19(f_even, f_odd);
				//*************************************************************************
				// 		Swap the distributions for momentum transport
				//*************************************************************************
				SwapD3Q19(ID, f_even, f_odd, Nx, Ny, Nz);
				//*************************************************************************
				// Wait for communications to complete and unpack the distributions
				ScaLBL_Comm.RecvD3Q19(f_even, f_odd);
				//*************************************************************************

				if (pBC && kproc == 0)	{
					PressureBC_inlet(f_even,f_odd,din,Nx,Ny,Nz);
				}

				if (pBC && kproc == nprocz-1){
					PressureBC_outlet(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
				}
				//...................................................................................
				DeviceBarrier();
				MPI_Barrier(comm);

				// Timestep completed!
				timestep++;


				if (rank==0){
					// write out csv file
					printf("D32 Fx Fy Fz vx vy vz err1d Fo Re K err\n");
				}

				if (timestep%500 == 0){
					//...........................................................................
					// Copy the data for for the analysis timestep
					//...........................................................................
					// Copy the phase from the GPU -> CPU
					//...........................................................................
					DeviceBarrier();
					ComputePressureD3Q19(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
					ComputeVelocityD3Q19(ID,f_even,f_odd,Velocity,Nx,Ny,Nz);
					CopyToHost(Averages.Press.data(),Pressure,N*sizeof(double));
					CopyToHost(Averages.Vel_x.data(),&Velocity[0],N*sizeof(double));
					CopyToHost(Averages.Vel_y.data(),&Velocity[N],N*sizeof(double));
					CopyToHost(Averages.Vel_z.data(),&Velocity[2*N],N*sizeof(double));

					// Way more work than necessary -- this is just to get the solid interfacial area!!
					Averages.Initialize();
					Averages.UpdateMeshValues();
					Averages.ComputeLocal();
					Averages.Reduce();

					double vawx = -Averages.vaw_global(0);
					double vawy = -Averages.vaw_global(1);
					double vawz = -Averages.vaw_global(2);

					// Compute local measures
					err = Re; // previous Reynolds number
					D32 = 6.0*(Dm.Volume-Averages.vol_w_global)/Averages.As_global;
					mag_force = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
					Fo = D32*D32*D32*mag_force/viscosity/viscosity;
					// .... 1-D flow should be aligned with force ...
					velocity = vawx*Fx/mag_force + vawy*Fy/mag_force + vawz*Fz/mag_force;
					err1D = fabs(velocity-sqrt(vawx*vawx+vawy*vawy+vawz*vawz))/velocity;
					//.......... Computation of the Reynolds number Re ..............
					Re = D32*velocity/viscosity;
					err = fabs(Re-err);

					if (rank==0){
						// ************* DIMENSIONLESS FORCHEIMER EQUATION *************************
						//  Dye, A.L., McClure, J.E., Gray, W.G. and C.T. Miller
						//  Description of Non-Darcy Flows in Porous Medium Systems
						//  Physical Review E 87 (3), 033012
						//  Fo := density*D32^3*(density*force) / (viscosity^2)
						//	Re := density*D32*velocity / viscosity
						//  Fo = a*Re + b*Re^2
						// *************************************************************************
						printf("%f ",D32);
						printf("%.5g,%.5g,%.5g ",Fx,Fy,Fz);
						printf("%.5g,%.5g,%.5g ",vawx,vawy,vawz);
						printf("%.5g ",err1D);
						printf("%5g ", Fo);
						printf("%.5g ", Re);
						printf("%.5g ", Re/Fo);
						printf("%.5g\n", err);
					}
				}
			}

			// Write steady state variables to csv file
			if (rank==0){
				fprintf(NONDARCY,"%.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g\n",D32,Fx,Fy,Fz,vx,vy,vz,Re,Fo);
				fflush(NONDARCY);
			}
		}
		//************************************************************************/
		fclose(NONDARCY);
		DeviceBarrier();
		MPI_Barrier(comm);
		stoptime = MPI_Wtime();
		if (rank==0) printf("-------------------------------------------------------------------\n");
		// Compute the walltime per timestep
		cputime = (stoptime - starttime)/timestep;
		// Performance obtained from each node
		double MLUPS = double(Nx*Ny*Nz)/cputime/1000000;

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
