#include <iostream>
#include "common/MPI_Helpers.h"
#include "common/Utilities.h"
#include "common/ScaLBL.h"


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
		int i,j,k,n;
		int Nx,Ny,Nz;
		bool pBC=true;
		int nprocx,nprocy,nprocz;
		double Lx,Ly,Lz;
		Nx = Ny = Nz = 16;
		Lx = Ly = Lz = 1.f;
		nprocx=nprocy=nprocz=1;
		double din,dout;
		int BC=1;

		Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);

		Nz += 2;
		Nx = Ny = Nz;	// Cubic domain

		int N = Nx*Ny*Nz;
		int dist_mem_size = N*sizeof(double);

		//.......................................................................
		// Assign the phase ID
		//.......................................................................
		char *id;
		id = new char[N];
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
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
		for (i=0; i<Dm.Nx*Dm.Ny*Dm.Nz; i++) Dm.id[i] = id[i];
		Dm.CommInit(comm);
		//...........................................................................
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");
		// Create a communicator for the device
		ScaLBL_Communicator ScaLBL_Comm(Dm);

		//...........device phase ID.................................................
		if (rank==0)	printf ("Copying phase ID to device \n");
		char *ID;
		ScaLBL_AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
		// Copy to the device
		ScaLBL_CopyToDevice(ID, id, N);
		//...........................................................................

		//...........................................................................
		//				MAIN  VARIABLES ALLOCATED HERE
		//...........................................................................
		// LBM variables
		if (rank==0)	printf ("Allocating distributions \n");
		//......................device distributions.................................
		double *f_even,*f_odd;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &f_even, 10*dist_mem_size);	// Allocate device memory
		ScaLBL_AllocateDeviceMemory((void **) &f_odd, 9*dist_mem_size);	// Allocate device memory
		//...........................................................................

		//...........................................................................
		//			INITIALIZE DISTRIBUTIONS
		//...........................................................................
		//...........................................................................
		if (rank==0)	printf("Setting the distributions, size = %i\n", N);
		//...........................................................................
		ScaLBL_D3Q19_Init(ID, f_even, f_odd, Nx, Ny, Nz);
		//......................................................................

		double flux = 0.1;

		//printf("kproc=%i \n",Dm.kproc);
		if (pBC && Dm.kproc == 0){
			din = ScaLBL_D3Q19_Flux_BC_z(f_even,f_odd,flux,Nx,Ny,Nz);
			printf("Computed inlet pressure: %f \n", din);
			ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
		}

		if (pBC && Dm.kproc == nprocz-1){
			dout = ScaLBL_D3Q19_Flux_BC_Z(f_even,f_odd,flux,Nx,Ny,Nz,Nx*Ny*(Nz-2));
			printf("Computed outlet pressure: %f \n", dout);
			ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
		}

		printf("Compute velocity \n");
		double *dvc_vel;
		ScaLBL_AllocateDeviceMemory((void **) &dvc_vel, 3*N*sizeof(double));
		ScaLBL_D3Q19_Velocity(ID,f_even,f_odd,dvc_vel,Nx,Ny,Nz);

		printf("Copying velocity to host \n");
		double * vel;
		vel = new double [3*N];
		ScaLBL_CopyToHost(vel,dvc_vel,3*N*sizeof(double));

		// Check the first layer
		int offset=2*N;
		double sum=0.f;
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				n = Nx*Ny+j*Nx + i;
				sum += vel[offset+n];
			}
		}
		double err;

		double value;
		value = sum/(Nx*Ny)/din;
		printf("Inlet Flux: input=%f, output=%f \n",flux,value);
		err = fabs(flux - value);
		if (err > 1e-14){
			error = 1;
			printf("  Inlet error %f \n",err);
		}

		// Check the last layer
		sum=0.f;
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				n = (Nz-2)*Nx*Ny+j*Nx + i;
				sum += vel[offset+n];
			}
		}
		value = sum/(Nx*Ny)/dout;
		err = fabs(flux - value);
		printf("Outlet Flux: input=%f, output=%f \n",flux,value);
		err = fabs(flux - value);
		if (err > 1e-14){
			error += 2;
			printf("   Outlet error %f \n",err);
		}

	}
	// Finished
	MPI_Barrier(comm);
	MPI_Finalize();
    return error; 
}
