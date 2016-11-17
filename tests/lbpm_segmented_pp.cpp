/*
 * Pre-processor to generate signed distance function from segmented data
 * segmented data should be stored in a raw binary file as 1-byte integer (type char)
 * will output distance functions for phases
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "common/Array.h"
#include "common/Domain.h"
#include "common/TwoPhase.h"

inline void MeanFilter(DoubleArray &Mesh){
	for (int k=1; k<(int)Mesh.size(2)-1; k++){
		for (int j=1; j<(int)Mesh.size(1)-1; j++){
			for (int i=1; i<(int)Mesh.size(0)-1; i++){
				double sum;
				sum=Mesh(i,j,k)+Mesh(i+1,j,k)+Mesh(i-1,j,k)+Mesh(i,j+1,k)+Mesh(i,j-1,k)+
						+Mesh(i,j,k+1)+Mesh(i,j,k-1);
				Mesh(i,j,k) = sum/7.0;
			}
		}
	}
}

inline double minmod(double &a, double &b){

	double value;

	value = a;
	if 	( a*b < 0.0)	    value=0.0;
	else if (fabs(a) > fabs(b)) value = b;

	return value;
}


inline double Eikonal(DoubleArray &Distance, char *ID, Domain &Dm, int timesteps){

	/*
	 * This routine converts the data in the Distance array to a signed distance
	 * by solving the equation df/dt = sign(1-|grad f|), where Distance provides
	 * the values of f on the mesh associated with domain Dm
	 * It has been tested with segmented data initialized to values [-1,1]
	 * and will converge toward the signed distance to the surface bounding the associated phases
	 *
	 * Reference:
	 * Min C (2010) On reinitializing level set functions, Journal of Computational Physics	229
	 */

	int i,j,k;
	double dt=0.1;
	double Dx,Dy,Dz;
	double Dxp,Dxm,Dyp,Dym,Dzp,Dzm;
	double Dxxp,Dxxm,Dyyp,Dyym,Dzzp,Dzzm;
	double sign,norm;
	double LocalVar,GlobalVar,LocalMax,GlobalMax;

	int xdim,ydim,zdim;
	xdim=Dm.Nx-2;
	ydim=Dm.Ny-2;
	zdim=Dm.Nz-2;
	fillHalo<double> fillData(Dm.Comm, Dm.rank_info,xdim,ydim,zdim,1,1,1,0,1);

	// Arrays to store the second derivatives
	DoubleArray Dxx(Dm.Nx,Dm.Ny,Dm.Nz);
	DoubleArray Dyy(Dm.Nx,Dm.Ny,Dm.Nz);
	DoubleArray Dzz(Dm.Nx,Dm.Ny,Dm.Nz);

	int count = 0;
	while (count < timesteps){

		// Communicate the halo of values
		fillData.fill(Distance);

		// Compute second order derivatives
		for (k=1;k<Dm.Nz-1;k++){
			for (j=1;j<Dm.Ny-1;j++){
				for (i=1;i<Dm.Nx-1;i++){
					Dxx(i,j,k) = Distance(i+1,j,k) + Distance(i-1,j,k) - 2*Distance(i,j,k);
					Dyy(i,j,k) = Distance(i,j+1,k) + Distance(i,j-1,k) - 2*Distance(i,j,k);
					Dzz(i,j,k) = Distance(i,j,k+1) + Distance(i,j,k-1) - 2*Distance(i,j,k);
				}
			}
		}
		fillData.fill(Dxx);
		fillData.fill(Dyy);
		fillData.fill(Dzz);

		LocalMax=LocalVar=0.0;
		// Execute the next timestep
		for (k=1;k<Dm.Nz-1;k++){
			for (j=1;j<Dm.Ny-1;j++){
				for (i=1;i<Dm.Nx-1;i++){

					int n = k*Dm.Nx*Dm.Ny + j*Dm.Nx + i;

					sign = -1;
					if (ID[n] == 1) sign = 1;

					// local second derivative terms
					Dxxp = minmod(Dxx(i,j,k),Dxx(i+1,j,k));
					Dyyp = minmod(Dyy(i,j,k),Dyy(i,j+1,k));
					Dzzp = minmod(Dzz(i,j,k),Dzz(i,j,k+1));
					Dxxm = minmod(Dxx(i,j,k),Dxx(i-1,j,k));
					Dyym = minmod(Dyy(i,j,k),Dyy(i,j-1,k));
					Dzzm = minmod(Dzz(i,j,k),Dzz(i,j,k-1));

					/* //............Compute upwind derivatives ...................
                    Dxp = Distance(i+1,j,k) - Distance(i,j,k) + 0.5*Dxxp;
                    Dyp = Distance(i,j+1,k) - Distance(i,j,k) + 0.5*Dyyp;
                    Dzp = Distance(i,j,k+1) - Distance(i,j,k) + 0.5*Dzzp;


                    Dxm = Distance(i,j,k) - Distance(i-1,j,k) + 0.5*Dxxm;
                    Dym = Distance(i,j,k) - Distance(i,j-1,k) + 0.5*Dyym;
                    Dzm = Distance(i,j,k) - Distance(i,j,k-1) + 0.5*Dzzm;
					 */
					Dxp = Distance(i+1,j,k)- Distance(i,j,k) - 0.5*Dxxp;
					Dyp = Distance(i,j+1,k)- Distance(i,j,k) - 0.5*Dyyp;
					Dzp = Distance(i,j,k+1)- Distance(i,j,k) - 0.5*Dzzp;

					Dxm = Distance(i,j,k) - Distance(i-1,j,k) + 0.5*Dxxm;
					Dym = Distance(i,j,k) - Distance(i,j-1,k) + 0.5*Dyym;
					Dzm = Distance(i,j,k) - Distance(i,j,k-1) + 0.5*Dzzm;

					// Compute upwind derivatives for Godunov Hamiltonian
					if (sign < 0.0){
						if (Dxp + Dxm > 0.f)  	Dx = Dxp*Dxp;
						else					Dx = Dxm*Dxm;

						if (Dyp + Dym > 0.f)  	Dy = Dyp*Dyp;
						else					Dy = Dym*Dym;

						if (Dzp + Dzm > 0.f)  	Dz = Dzp*Dzp;
						else					Dz = Dzm*Dzm;
					}
					else{

						if (Dxp + Dxm < 0.f)  	Dx = Dxp*Dxp;
						else					Dx = Dxm*Dxm;

						if (Dyp + Dym < 0.f)  	Dy = Dyp*Dyp;
						else					Dy = Dym*Dym;

						if (Dzp + Dzm < 0.f)  	Dz = Dzp*Dzp;
						else					Dz = Dzm*Dzm;
					}

					//Dx = max(Dxp*Dxp,Dxm*Dxm);
					//Dy = max(Dyp*Dyp,Dym*Dym);
					//Dz = max(Dzp*Dzp,Dzm*Dzm);

					norm=sqrt(Dx + Dy + Dz);
					if (norm > 1.0) norm=1.0;

					Distance(i,j,k) += dt*sign*(1.0 - norm);
					LocalVar +=  dt*sign*(1.0 - norm);

					if (fabs(dt*sign*(1.0 - norm)) > LocalMax)
						LocalMax = fabs(dt*sign*(1.0 - norm));
				}
			}
		}

		MPI_Allreduce(&LocalVar,&GlobalVar,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
		MPI_Allreduce(&LocalMax,&GlobalMax,1,MPI_DOUBLE,MPI_MAX,Dm.Comm);
		GlobalVar /= (Dm.Nx-2)*(Dm.Ny-2)*(Dm.Nz-2)*Dm.nprocx*Dm.nprocy*Dm.nprocz;
		count++;

		if (count%50 == 0 && Dm.rank==0 )
			printf("Time=%i, Max variation=%f, Global variation=%f \n",count,GlobalMax,GlobalVar);

		if (fabs(GlobalMax) < 1e-5){
			if (Dm.rank==0) printf("Exiting with max tolerance of 1e-5 \n");
			count=timesteps;
		}
	}
	return GlobalVar;
}


int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	{	
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
		double Lx, Ly, Lz;
		int Nx,Ny,Nz;
		int i,j,k,n;
		int BC=0;
		//  char fluidValue,solidValue;

		std::vector<char> solidValues;
		std::vector<char> nwpValues;
		std::string line;

		if (rank==0){
			ifstream domain("Domain.in");
			domain >> nprocx;
			domain >> nprocy;
			domain >> nprocz;
			domain >> nx;
			domain >> ny;
			domain >> nz;
			domain >> nspheres;
			domain >> Lx;
			domain >> Ly;
			domain >> Lz;

		}
		MPI_Barrier(comm);
		// Computational domain
		MPI_Bcast(&nx,1,MPI_INT,0,comm);
		MPI_Bcast(&ny,1,MPI_INT,0,comm);
		MPI_Bcast(&nz,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
		MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
		MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
		//.................................................
		MPI_Barrier(comm);

		// Check that the number of processors >= the number of ranks
		if ( rank==0 ) {
			printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
			printf("Number of MPI ranks used: %i \n", nprocs);
			printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
		}
		if ( nprocs < nprocx*nprocy*nprocz ){
			ERROR("Insufficient number of processors");
		}

		char LocalRankFilename[40];

		int N = (nx+2)*(ny+2)*(nz+2);
		Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
		for (n=0; n<N; n++) Dm.id[n]=1;
		Dm.CommInit(comm);

		// Read the phase ID
		size_t readID;
		sprintf(LocalRankFilename,"ID.%05i",rank);
		FILE *ID = fopen(LocalRankFilename,"rb");
		readID=fread(Dm.id,1,N,ID);
		if (readID != size_t(N)) printf("lbpm_segmented_pp: Error reading ID \n");
		fclose(ID);

		// Initialize the domain and communication
		nx+=2; ny+=2; nz+=2;

		char *id;
		id = new char [N];
		TwoPhase Averages(Dm);
		//	DoubleArray Distance(nx,ny,nz);
		//	DoubleArray Phase(nx,ny,nz);

		int count = 0;
		// Solve for the position of the solid phase
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n = k*nx*ny+j*nx+i;
					// Initialize the solid phase
					if (Dm.id[n] == 0)	id[n] = 0;
					else		      	id[n] = 1;
				}
			}
		}
		// Initialize the signed distance function
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n=k*nx*ny+j*nx+i;
					// Initialize distance to +/- 1
					Averages.SDs(i,j,k) = 2.0*double(id[n])-1.0;
				}
			}
		}
		MeanFilter(Averages.SDs);

		double LocalVar, TotalVar;
		if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
		int Maxtime=10*max(max(Dm.Nx*Dm.nprocx,Dm.Ny*Dm.nprocy),Dm.Nz*Dm.nprocz);
		LocalVar = Eikonal(Averages.SDs,id,Dm,Maxtime);

		MPI_Allreduce(&LocalVar,&TotalVar,1,MPI_DOUBLE,MPI_SUM,comm);
		TotalVar /= nprocs;
		if (rank==0) printf("Final variation in signed distance function %f \n",TotalVar);

		sprintf(LocalRankFilename,"SignDist.%05i",rank);
		FILE *DIST = fopen(LocalRankFilename,"wb");
		fwrite(Averages.SDs.data(),8,N,DIST);
		fclose(DIST);

	}
	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;

}
