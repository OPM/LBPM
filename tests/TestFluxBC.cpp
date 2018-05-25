#include <iostream>
#include "common/MPI_Helpers.h"
#include "common/Utilities.h"
#include "common/ScaLBL.h"

std::shared_ptr<Database> loadInputs( int nprocs )
{
  //auto db = std::make_shared<Database>( "Domain.in" );
    auto db = std::make_shared<Database>();
    db->putScalar<int>( "BC", 0 );
    db->putVector<int>( "nproc", { 1, 1, 1 } );
    db->putVector<int>( "n", { 16, 16, 16 } );
    db->putScalar<int>( "nspheres", 1 );
    db->putVector<double>( "L", { 1, 1, 1 } );
    return db;
}

int main (int argc, char **argv)
{
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	int rank = MPI_WORLD_RANK();
	int nprocs = MPI_WORLD_SIZE();

	// set the error code
	// Note: the error code should be consistent across all processors
	int error = 0;

	if (nprocs != 1){
		printf("FAIL: Unit test TestFluxBC requires 1 MPI process! \n");
		ASSERT(nprocs==1);
	}
	{
	  int i,j,k,n,Np;
		bool pBC=true;
		double Lx,Ly,Lz;
		Lx = Ly = Lz = 1.f;
		double din,dout;
		int BC=1;

	    // Load inputs
	    auto db = loadInputs( nprocs );
	    int Nx = db->getVector<int>( "n" )[0];
	    int Ny = db->getVector<int>( "n" )[1];
	    int Nz = db->getVector<int>( "n" )[2];
	    int nprocx = db->getVector<int>( "nproc" )[0];
	    int nprocy = db->getVector<int>( "nproc" )[1];
	    int nprocz = db->getVector<int>( "nproc" )[2];
		std::shared_ptr<Domain> Dm(new Domain(db,comm));
		
		Nx += 2;   Ny+=2;	Nz += 2;
		Nx = Ny = Nz;	// Cubic domain
		int N = Nx*Ny*Nz;

		//.......................................................................
		// Assign the phase ID
		//.......................................................................
		char *id;
		id = new char[N];
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					id[n] = 0;
				}
			}
		}
		// Set up parallel plates
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=2;i<Nx-2;i++){
					n = k*Nx*Ny+j*Nx+i;
					id[n] = 1;
				}
			}
		}

		//.........................................................
		// don't perform computations at the eight corners
		id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
		id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;
		//.........................................................
		// Initialize communication structures in averaging domain
		for (i=0; i<Dm->Nx*Dm->Ny*Dm->Nz; i++) Dm->id[i] = id[i];
		Dm->CommInit();
		Np=Dm->PoreCount();
		//................................................
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");
		// Create a communicator for the device
		std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm(new ScaLBL_Communicator(Dm));

		if (rank==0)	printf ("Set up the neighborlist \n");
		int Npad=Np+32;
		int neighborSize=18*Npad*sizeof(int);
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Npad];

		Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Dm->id,Np);
		MPI_Barrier(comm);

		//......................device distributions.................................
		int dist_mem_size = Np*sizeof(double);
		if (rank==0)	printf ("Allocating distributions \n");

		int *NeighborList;
		int *dvcMap;
		double *fq;
		
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);
		//...........................................................................
		// Update GPU data structures
		if (rank==0)	printf ("Setting up device map and neighbor list \n");
		int *TmpMap;
		TmpMap=new int[Np];
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

		//...........................................................................
		//			INITIALIZE DISTRIBUTIONS
		//...........................................................................
		//...........................................................................
		if (rank==0)	printf("Initializing the distributions, size = %i\n", Np);
		//...........................................................................
		ScaLBL_D3Q19_Init(fq, Np);
		//......................................................................
		double flux = 1.0;
		int timestep=0; 

		din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
		
		printf("Computed pressure for flux = %f\n",din);
		
		printf("Compute velocity \n");
		double *dvc_vel;
		ScaLBL_AllocateDeviceMemory((void **) &dvc_vel, 3*Np*sizeof(double));
		ScaLBL_D3Q19_Momentum(fq,dvc_vel,Np);

		printf("Copying velocity to host \n");
    	double *VEL;
    	VEL= new double [3*Np];
    	int SIZE=3*Np*sizeof(double);
    	ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    	ScaLBL_CopyToHost(&VEL[0],&dvc_vel[0],SIZE);

		double err,value,Q;

    	Q = 0.f;    	
    	k=1;
    	for (j=1;j<Ny-1;j++){
    		for (i=1;i<Nx-1;i++){
    			n = k*Nx*Ny+j*Nx+i;
    			if (Dm->id[n] > 0){
    				int idx = Map(i,j,k);
    				Q += VEL[2*Np+idx];
    			}
    		}
    	}

    	// respect backwards read / write!!!
		printf("Inlet Flux: input=%f, output=%f \n",flux,Q);
		err = fabs(flux + Q);
		if (err > 1e-12){
			error = 1;
			printf("  Inlet error %f \n",err);
		}		
		
		// Consider a larger number of timesteps and simulate flow
		double Fx, Fy, Fz;
		double tau = 1.0;
		double mu=(tau-0.5)/3.0;
		double rlx_setA=1.0/tau;
		double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
		dout=1.f;
		Fx = 0; Fy = 0; Fz = 0.f;
		ScaLBL_D3Q19_Init(fq, Np);
		timestep=1;
		printf("*** Running 2000 timesteps as a test *** \n");

		while (timestep < 2000) {

			ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
			ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq,  ScaLBL_Comm->first_interior, ScaLBL_Comm->last_interior, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
			ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
			din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq, 0, ScaLBL_Comm->next, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
			timestep++;

			ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
			ScaLBL_D3Q19_AAeven_MRT(fq,  ScaLBL_Comm->first_interior, ScaLBL_Comm->last_interior, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
			ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
			din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
			ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			ScaLBL_D3Q19_AAeven_MRT(fq, 0, ScaLBL_Comm->next, Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
			timestep++;
			//************************************************************************/

		}
		
		printf("Compute velocity \n");
		ScaLBL_D3Q19_Momentum(fq,dvc_vel,Np);

		printf("Copying velocity to host \n");
		ScaLBL_CopyToHost(&VEL[0],&dvc_vel[0],SIZE);

		printf("Printing velocity profile \n");
		j=4;
		for (k=1;k<Nz-1;k++){
			for (i=1;i<Nx-1;i++){ 
				n = k*Nx*Ny+j*Nx+i;
				if (Dm->id[n] > 0){
					int idx = Map(i,j,k);
					double vz = VEL[2*Np+idx];
					printf("%f ",vz);
				}
			}
			printf("\n");
		}
		


		printf("Printing INLET  velocity profile \n");
		k=1;
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){ 
				n = k*Nx*Ny+j*Nx+i;
				if (Dm->id[n] > 0){
					int idx = Map(i,j,k);
					double vz = VEL[2*Np+idx];
					printf("%f ",vz);
				}
			}
			printf("\n");
		}
		
    	Q = 0.f;    	
    	k=1;
    	for (j=1;j<Ny-1;j++){
    		for (i=1;i<Nx-1;i++){
    			n = k*Nx*Ny+j*Nx+i;
    			if (Dm->id[n] > 0){
    				int idx = Map(i,j,k);
    				Q += VEL[2*Np+idx];
    			}
    		}
    	}

		printf("Inlet Flux: input=%f, output=%f \n",flux,Q);
		err = fabs(flux - Q);
		if (err > 1e-12){
			error = 1;
			printf("  Inlet error %f \n",err);
		}
		

	}
	// Finished
	MPI_Barrier(comm);
	MPI_Finalize();
    return error; 
}
