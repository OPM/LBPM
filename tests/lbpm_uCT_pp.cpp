// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "common/Array.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Mesh.h"
#include "IO/Writer.h"
#include "IO/netcdf.h"

#include "ProfilerApp.h"

inline void Med3D(Array<float> &Input, Array<float> &Output){
	// Perform a 3D Median filter on Input array with specified width
	int i,j,k,ii,jj,kk;
	int imin,jmin,kmin,imax,jmax,kmax;

	float *List;
	List=new float[27];

	int Nx = int(Input.size(0));
	int Ny = int(Input.size(1));
	int Nz = int(Input.size(2));

	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){

				// Just use a 3x3x3 window (hit recursively if needed)
				imin = i-1;
				jmin = j-1;
				kmin = k-1;
				imax = i+2;
				jmax = j+2;
				kmax = k+2;

				// Populate the list with values in the window
				int Number=0;
				for (kk=kmin; kk<kmax; kk++){
					for (jj=jmin; jj<jmax; jj++){
						for (ii=imin; ii<imax; ii++){
							List[Number++] = Input(ii,jj,kk);
						}
					}
				}
				// Sort the first 5 entries and return the median
				for (ii=0; ii<14; ii++){
					for (jj=ii+1; jj<27; jj++){
						if (List[jj] < List[ii]){
							float tmp = List[ii];
							List[ii] = List[jj];
							List[jj] = tmp;
						}
					}
				}
				// Return the median
				Output(i,j,k) = List[13];
			}
		}
	}
}

inline void Sparsify(Array<float> &Fine, Array<float> &Coarse){

	// Create sparse version of Fine mesh to reduce filtering costs
	int i,j,k,ii,jj,kk;
	float x,y,z;

	// Fine mesh
	int Nx = int(Fine.size(0));
	int Ny = int(Fine.size(1));
	int Nz = int(Fine.size(2));

	// Coarse mesh
	int nx = int(Coarse.size(0));
	int ny = int(Coarse.size(1));
	int nz = int(Coarse.size(2));

	// compute the stride
	float hx = float(Nx-1) / float (nx-1);
	float hy = float(Ny-1) / float (ny-1);
	float hz = float(Nz-1) / float (nz-1);

	// Fill in the coarse mesh
	for (k=0; k<nz; k++){
		for (j=0; j<ny; j++){
			for (i=0; i<nx; i++){

				x = i*hx;
				y = j*hy;
				z = k*hz;

				ii = int(floor(x));
				jj = int(floor(y));
				kk = int(floor(z));

				// get the eight values in the cell
				float v1 = Fine(ii,jj,kk);
				float v2 = Fine(ii+1,jj,kk);
				float v3 = Fine(ii,jj+1,kk);
				float v4 = Fine(ii+1,jj+1,kk);
				float v5 = Fine(ii,jj,kk+1);
				float v6 = Fine(ii+1,jj,kk+1);
				float v7 = Fine(ii,jj+1,kk+1);
				float v8 = Fine(ii+1,jj+1,kk+1);

				Coarse(i,j,k)=0.125*(v1+v2+v3+v4+v5+v6+v7+v8);

				//Coarse(i,j,k) = Fine(ii,jj,kk);

			}
		}
	}
}

inline void InterpolateMesh(Array<float> &Coarse, Array<float> &Fine){

	// Interpolate values from a Coarse mesh to a fine one
	// This routine assumes that the mesh boundaries match
	int i,j,k,ii,jj,kk;
	float x,y,z;
	Array<float> Corners(2,2,2);
	float a,b,c,d,e,f,g,h;

	// Fine mesh
	int Nx = int(Fine.size(0));
	int Ny = int(Fine.size(1));
	int Nz = int(Fine.size(2));

	// Coarse mesh
	int nx = int(Coarse.size(0));
	int ny = int(Coarse.size(1));
	int nz = int(Coarse.size(2));

	// compute the stride
	float hx = float(Nx-1) / float (nx-1);
	float hy = float(Ny-1) / float (ny-1);
	float hz = float(Nz-1) / float (nz-1);

	// value to map distance between meshes (since distance is in voxels)
	//  usually hx=hy=hz (or something very close)
	//  the mapping is not exact
	//  however, it's assumed the coarse solution will be refined
	//  a good guess is the goal here!
	float mapvalue = sqrt(hx*hx+hy*hy+hz*hz);

	// Interpolate to the fine mesh
	for (k=0; k<nz-1; k++){
		for (j=0; j<ny-1; j++){
			for (i=0; i<nx-1; i++){

				// get the eight values in the cell
				Corners(0,0,0) = mapvalue*Coarse(i,j,k);
				Corners(1,0,0) = mapvalue*Coarse(i+1,j,k);
				Corners(0,1,0) = mapvalue*Coarse(i,j+1,k);
				Corners(1,1,0) = mapvalue*Coarse(i+1,j+1,k);
				Corners(0,0,1) = mapvalue*Coarse(i,j,k+1);
				Corners(1,0,1) = mapvalue*Coarse(i+1,j,k+1);
				Corners(0,1,1) = mapvalue*Coarse(i,j+1,k+1);
				Corners(1,1,1) = mapvalue*Coarse(i+1,j+1,k+1);

				// coefficients of the tri-linear approximation
				a = Corners(0,0,0);
				b = Corners(1,0,0)-a;
				c = Corners(0,1,0)-a;
				d = Corners(0,0,1)-a;
				e = Corners(1,1,0)-a-b-c;
				f = Corners(1,0,1)-a-b-d;
				g = Corners(0,1,1)-a-c-d;
				h = Corners(1,1,1)-a-b-c-d-e-f-g;

				// Interpolate to each point on the fine mesh
				for (kk=int(ceil(k*hz)); kk<int(ceil((k+1)*hz)); kk++){
					for (jj=int(ceil(j*hy)); jj<int(ceil((j+1)*hy)); jj++){
						for (ii=int(ceil(i*hx)); ii<int(ceil((i+1)*hx)); ii++){

							// get the value within the unit cube
							x = (ii-i*hx)/hx;
							y = (jj-j*hy)/hy;
							z = (kk-k*hz)/hz;

							if (ii<Nx && jj<Ny && kk<Nz)
								Fine(ii,jj,kk) = a + b*x + c*y+d*z + e*x*y + f*x*z + g*y*z + h*x*y*z;

						}
					}
				}

			}
		}
	}

}

inline float minmod(float &a, float &b){

	float value;

	value = a;
	if 	( a*b < 0.0)	    value=0.0;
	else if (fabs(a) > fabs(b)) value = b;

	return value;
}


inline float Eikonal3D(Array<float> &Distance, Array<char> &ID, Domain &Dm, int timesteps){

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
	float dt=0.1;
	float Dx,Dy,Dz;
	float Dxp,Dxm,Dyp,Dym,Dzp,Dzm;
	float Dxxp,Dxxm,Dyyp,Dyym,Dzzp,Dzzm;
	float sign,norm;
	float LocalVar,GlobalVar,LocalMax,GlobalMax;

	int xdim,ydim,zdim;
	xdim=Dm.Nx-2;
	ydim=Dm.Ny-2;
	zdim=Dm.Nz-2;
	fillHalo<float> fillData(Dm.Comm, Dm.rank_info,xdim,ydim,zdim,1,1,1,0,1);

	// Arrays to store the second derivatives
	Array<float> Dxx(Dm.Nx,Dm.Ny,Dm.Nz);
	Array<float> Dyy(Dm.Nx,Dm.Ny,Dm.Nz);
	Array<float> Dzz(Dm.Nx,Dm.Ny,Dm.Nz);

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
					if (ID(i,j,k) == 1) sign = 1;

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
					Dxp = Distance(i+1,j,k);
					Dyp = Distance(i,j+1,k);
					Dzp = Distance(i,j,k+1);

					Dxm = Distance(i-1,j,k);
					Dym = Distance(i,j-1,k);
					Dzm = Distance(i,j,k-1);

					// Compute upwind derivatives for Godunov Hamiltonian
					if (sign < 0.0){
						if (Dxp > Dxm)  Dx = Dxp - Distance(i,j,k) + 0.5*Dxxp;
						else			Dx = Distance(i,j,k) - Dxm + 0.5*Dxxm;

						if (Dyp > Dym)  Dy = Dyp - Distance(i,j,k) + 0.5*Dyyp;
						else			Dy = Distance(i,j,k) - Dym + 0.5*Dyym;

						if (Dzp > Dzm)  Dz = Dzp - Distance(i,j,k) + 0.5*Dzzp;
						else			Dz = Distance(i,j,k) - Dzm + 0.5*Dzzm;
					}
					else{
						if (Dxp < Dxm)  Dx = Dxp - Distance(i,j,k) + 0.5*Dxxp;
						else			Dx = Distance(i,j,k) - Dxm + 0.5*Dxxm;

						if (Dyp < Dym)  Dy = Dyp - Distance(i,j,k) + 0.5*Dyyp;
						else			Dy = Distance(i,j,k) - Dym + 0.5*Dyym;

						if (Dzp < Dzm)  Dz = Dzp - Distance(i,j,k) + 0.5*Dzzp;
						else			Dz = Distance(i,j,k) - Dzm + 0.5*Dzzm;
					}

					norm=sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
					if (norm > 1.0) norm=1.0;

					Distance(i,j,k) += dt*sign*(1.0 - norm);
					LocalVar +=  dt*sign*(1.0 - norm);

					if (fabs(dt*sign*(1.0 - norm)) > LocalMax)
						LocalMax = fabs(dt*sign*(1.0 - norm));
				}
			}
		}

		MPI_Allreduce(&LocalVar,&GlobalVar,1,MPI_FLOAT,MPI_SUM,Dm.Comm);
		MPI_Allreduce(&LocalMax,&GlobalMax,1,MPI_FLOAT,MPI_MAX,Dm.Comm);
		GlobalVar /= (Dm.Nx-2)*(Dm.Ny-2)*(Dm.Nz-2)*Dm.nprocx*Dm.nprocy*Dm.nprocz;
		count++;

		if (count%50 == 0 && Dm.rank==0 )
			printf("	Time=%i, Max variation=%f, Global variation=%f \n",count,GlobalMax,GlobalVar);

		if (fabs(GlobalMax) < 1e-5){
			if (Dm.rank==0) printf("	Exiting with max tolerance of 1e-5 \n");
			count=timesteps;
		}
	}
	return GlobalVar;
}


inline int NLM3D(Array<float> &Input, Array<float> &Mean, Array<float> &Distance, Array<float> &Output,
		const int d, const float h){
	// Implemenation of 3D non-local means filter
	// 		d determines the width of the search volume
	// 		h is a free parameter for non-local means (i.e. 1/sigma^2)
	// 		Distance is the signed distance function
	// 		If Distance(i,j,k) < THRESHOLD_DIST then don't compute NLM

	float THRESHOLD_DIST = float(d);
	float weight, sum;
	int i,j,k,ii,jj,kk;
	int imin,jmin,kmin,imax,jmax,kmax;
	int returnCount=0;

	int Nx = int(Input.size(0));
	int Ny = int(Input.size(1));
	int Nz = int(Input.size(2));

	// Compute the local means
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){

				imin = max(0,i-d);
				jmin = max(0,j-d);
				kmin = max(0,k-d);
				imax = min(Nx-1,i+d);
				jmax = min(Ny-1,j+d);
				kmax = min(Nz-1,k+d);

				// Populate the list with values in the window
				sum = 0; weight=0;
				for (kk=kmin; kk<kmax; kk++){
					for (jj=jmin; jj<jmax; jj++){
						for (ii=imin; ii<imax; ii++){
							sum += Input(ii,jj,kk);
							weight++;
						}
					}
				}

				Mean(i,j,k) = sum / weight;
			}
		}
	}

	// Compute the non-local means
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){


				if (fabs(Distance(i,j,k)) < THRESHOLD_DIST){
					// compute the expensive non-local means
					sum = 0; weight=0;

					imin = max(0,i-d);
					jmin = max(0,j-d);
					kmin = max(0,k-d);
					imax = min(Nx-1,i+d);
					jmax = min(Ny-1,j+d);
					kmax = min(Nz-1,k+d);

					for (kk=kmin; kk<kmax; kk++){
						for (jj=jmin; jj<jmax; jj++){
							for (ii=imin; ii<imax; ii++){
								float tmp = Mean(i,j,k) - Mean(ii,jj,kk);
								sum += exp(-tmp*tmp*h)*Input(ii,jj,kk);
								weight += exp(-tmp*tmp*h);
							}
						}
					}

					returnCount++;
					//Output(i,j,k) = Mean(i,j,k);
					Output(i,j,k) = sum / weight;
				}
				else{
					// Just return the mean
					Output(i,j,k) = Mean(i,j,k);
				}
			}
		}
	}
	// Return the number of sites where NLM was applied
	return returnCount;
}

int main(int argc, char **argv)
{

	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);

	//std::vector<std::string> filenames;
	std::string filename;
	if (rank==0){
		if ( argc==0 ) {
			printf("At least one filename must be specified\n");
			return 1;
		}
		else {
			filename=std::string(argv[1]);
			printf("Input data file: %s\n",filename.c_str());
		}
	}
	//.......................................................................
	// Reading the domain information file
	//.......................................................................
	int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
	double Lx, Ly, Lz;
	int Nx,Ny,Nz;
	int i,j,k,n;
	int BC=0;

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
	//.................................................
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

	// Allocate local arrays for every MPI rank
	Array<float> LOCVOL(nx+2,ny+2,nz+2);

	// Get the rank info
	int N = (nx+2)*(ny+2)*(nz+2);
	Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
	for (k=0;k<nz+2;k++){
		for (j=0;j<ny+2;j++){
			for (i=0;i<nx+2;i++){
				n = k*(nx+2)*(ny+2)+j*(nx+2)+i;
				Dm.id[n] = 1;
			}
		}
	}
	Dm.CommInit(comm);

	PROFILE_START("ReadVolume");
	{
		Array<float> VOLUME;

		// Read the input volume to rank 0 only, then distribute pieces to workers
		if (rank==0){
			// Open the netcdf file
			int fid = netcdf::open(filename);

			// Read all of the attributes
			std::vector<std::string> attr = netcdf::getAttNames( fid );
			for (size_t i=0; i<attr.size(); i++) {
				printf("Reading attribute %s\n",attr[i].c_str());
				netcdf::VariableType type = netcdf::getAttType( fid, attr[i] );
				if ( type == netcdf::STRING ){
					Array<std::string> tmp = netcdf::getAtt<std::string>( fid, attr[i] );
				}
				else{
					//Array<double> tmp = netcdf::getAtt<double>( fid, attr[i] );
				}
			}

			// Read the VOLUME data array
			std::string varname("VOLUME");
			printf("Reading %s\n",varname.c_str());
			VOLUME = netcdf::getVar<float>( fid, varname);
			Nx = int(VOLUME.size(0));
			Ny = int(VOLUME.size(1));
			Nz = int(VOLUME.size(2));
			printf("VOLUME dims =  %i x %i x %i \n",Nx,Ny,Nz);
			printf("Sucess!! \n");
			netcdf::close( fid );
		}
		PROFILE_SAVE("ReadVolume");

		MPI_Bcast(&Ny,1,MPI_INT,0,comm);
		MPI_Bcast(&Ny,1,MPI_INT,0,comm);
		MPI_Bcast(&Nz,1,MPI_INT,0,comm);

		MPI_Barrier(comm);

		// Set up the sub-domains
		int xStart,yStart,zStart;
		xStart=Nx/2;
		yStart=Ny/2;
		zStart=Nz/2;
		if (rank==0){
			printf("Distributing subdomains across %i processors \n",nprocs);
			printf("Process grid: %i x %i x %i \n",Dm.nprocx,Dm.nprocy,Dm.nprocz);
			printf("Subdomain size: %i \n",N);
			//	printf("Size of transition region: %i \n", z_transition_size);
			float *tmp;
			tmp = new float[N];
			for (int kp=0; kp<nprocz; kp++){
				for (int jp=0; jp<nprocy; jp++){
					for (int ip=0; ip<nprocx; ip++){
						// rank of the process that gets this subdomain
						int rnk = kp*Dm.nprocx*Dm.nprocy + jp*Dm.nprocx + ip;
						// Pack and send the subdomain for rnk
						for (k=0;k<nz+2;k++){
							for (j=0;j<ny+2;j++){
								for (i=0;i<nx+2;i++){
									int x = xStart + ip*nx + i-1;
									int y = yStart + jp*ny + j-1;
									int z = zStart + kp*nz + k-1;

									int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
									tmp[nlocal] = VOLUME(x,y,z);
								}
							}
						}
						if (rnk==0){
							for (k=0;k<nz+2;k++){
								for (j=0;j<ny+2;j++){
									for (i=0;i<nx+2;i++){
										int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
										LOCVOL(i,j,k) = tmp[nlocal];
									}
								}
							}
						}
						else{
							//printf("Sending data to process %i \n", rnk);
							MPI_Send(tmp,N,MPI_FLOAT,rnk,15,comm);
						}
					}
				}
			}
		}
		else{
			// Recieve the subdomain from rank = 0
			//printf("Ready to recieve data %i at process %i \n", N,rank);
			MPI_Recv(LOCVOL.get(),N,MPI_FLOAT,0,15,comm,MPI_STATUS_IGNORE);
		}
	}
	MPI_Barrier(comm);

	nx+=2; ny+=2; nz+=2;
	N=nx*ny*nz;

	if (rank==0) printf("All sub-domains recieved \n");

	// Initialize sparse domein
	int nsx,nsy,nsz;
	nsx=nx/8; nsy=ny/8; nsz=nz/8;

	Domain spDm(nsx-2,nsy-2,nsz-2,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
	for (k=0;k<nsz+2;k++){
		for (j=0;j<nsy+2;j++){
			for (i=0;i<nsx+2;i++){
				n = k*(nsx+2)*(nsy+2)+j*(nsx+2)+i;
				spDm.id[n] = 1;
			}
		}
	}
	spDm.CommInit(comm);

	fillHalo<float> fillFloat(Dm.Comm, Dm.rank_info,nx-2,ny-2,nz-2,1,1,1,0,1);
	fillHalo<char> fillChar(Dm.Comm, Dm.rank_info,nx-2,ny-2,nz-2,1,1,1,0,1);
	fillHalo<float> fillFloat_sp(spDm.Comm, spDm.rank_info,nsx-2,nsy-2,nsz-2,1,1,1,0,1);
	fillHalo<char>  fillChar_sp(spDm.Comm, spDm.rank_info,nsx-2,nsy-2,nsz-2,1,1,1,0,1);

	Array<float> spLOCVOL(nsx,nsy,nsz);	// this holds sparse original data
	Array<float> spM(nsx,nsy,nsz); 		// this holds sparse median filter
	Array<float> spSmooth(nsx,nsy,nsz); // this holds smoothed data

	Array<float> spDist(nsx,nsy,nsz);		// this holds sparse signed distance

	// sparse phase ID (segmented values)
	Array<char> spID(nsx,nsy,nsz);

	Array<char>  ID(nx,ny,nz);
	Array<float> Dist(nx,ny,nz);
	Array<float> MultiScaleSmooth(nx,ny,nz);
	Array<float> Mean(nx,ny,nz);
	Array<float> NonLocalMean(nx,ny,nz);

	if (rank==0) printf("Running segmentation workflow \n");
	if (rank==0) printf("Step 1. Sparsify space: \n");
	if (rank==0) printf("   Original Mesh: %ix%ix%i \n",nx,ny,nz);
	if (rank==0) printf("   Sparse Mesh: %ix%ix%i \n",nsx,nsy,nsz);

	// Sparsify the the mesh using a stride of 8
	Sparsify(LOCVOL,spLOCVOL);
	fillFloat_sp.fill(spLOCVOL);

	if (rank==0) printf("Step 2. Sparse median filter \n");
	// Compute the median filter on the sparse array
	Med3D(spLOCVOL,spM);
	fillFloat_sp.fill(spM);

	// quick & dirty sparse segmentation
	// this should be replaced
	//    (should use automated mixture model to approximate histograms)
	if (rank==0) printf("Step 3. Threshold for sparse segmentation \n");
	float THRESHOLD=50;
	for (k=0;k<nsz;k++){
		for (j=0;j<nsy;j++){
			for (i=0;i<nsx;i++){
				if (spM(i,j,k) > THRESHOLD) spID(i,j,k) = 0;
				else 						spID(i,j,k) = 1;
			}
		}
	}

	//..........................................
	// Compute the means for each region
	float mean_plus,mean_minus;
	int count_plus,count_minus;
	float mean_plus_global,mean_minus_global;
	int count_plus_global,count_minus_global;
	float *TmpMed;
	TmpMed = new float[nsx*nsy*nsz];

	// Compute median for regions of distance function
	count_plus=count_minus=0;
	mean_plus=mean_minus=0;
	for (k=1;k<nsz-1;k++){
		for (j=1;j<nsy-1;j++){
			for (i=1;i<nsx-1;i++){

				if (spDist(i,j,k) > 0.0)
					TmpMed[count_plus++]= spM(i,j,k);

			}
		}
	}
	for (int ii=0; ii<count_plus/2+1; ii++){
		for (int jj=ii+1; jj<count_plus; jj++){
			if (TmpMed[jj] < TmpMed[ii]){
				float tmp = TmpMed[ii];
				TmpMed[ii] = TmpMed[jj];
				TmpMed[jj] = tmp;
			}
		}
	}
	mean_plus = TmpMed[count_plus/2];

	for (k=1;k<nsz-1;k++){
		for (j=1;j<nsy-1;j++){
			for (i=1;i<nsx-1;i++){

				if (spDist(i,j,k) < 0.0)
					TmpMed[count_minus++]= spM(i,j,k);

			}
		}
	}
	for (int ii=0; ii<count_minus/2+1; ii++){
		for (int jj=ii+1; jj<count_minus; jj++){
			if (TmpMed[jj] < TmpMed[ii]){
				float tmp = TmpMed[ii];
				TmpMed[ii] = TmpMed[jj];
				TmpMed[jj] = tmp;
			}
		}
	}
	mean_minus = TmpMed[count_minus/2];

	MPI_Allreduce(&mean_plus,&mean_plus_global,1,MPI_FLOAT,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&mean_minus,&mean_minus_global,1,MPI_FLOAT,MPI_SUM,Dm.Comm);

	mean_plus_global /= nprocs;
	mean_minus_global /= nprocs;
	if (rank==0) printf("	Region 1 mean (+): %f, Region 2 mean (-): %f \n",mean_plus_global, mean_minus_global);
	delete[] TmpMed;
	//..........................................

	// intialize distance based on segmentation
	for (k=0;k<nsz;k++){
		for (j=0;j<nsy;j++){
			for (i=0;i<nsx;i++){
				spDist(i,j,k) = 2.0*spID(i,j,k)-1.0;
			}
		}
	}

	if (rank==0) printf("Step 4. Generate sparse distance function \n");
	// generate a sparse signed distance function
	Eikonal3D(spDist,spID,spDm,nsx*nprocx);

	if (rank==0) printf("Step 5. Interpolate to fine mesh \n");
	InterpolateMesh(spDist,Dist);

	float lambda = 0.5;
	for (k=0;k<nsz;k++){
		for (j=0;j<nsy;j++){
			for (i=0;i<nsx;i++){
				float dst = spDist(i,j,k);
				// use exponential weight based on the distance
				float temp = 1.0-exp(-lambda*fabs(dst)));
				float value;
				if (dst > 0){
					value = temp*mean_plus;
				}
				else{
					value = temp*mean_minus;
				}
				value += (1-temp)*spM(i,j,k);
				spSmooth(i,j,k) = value;
			}
		}
	}
	InterpolateMesh(spSmooth,MultiScaleSmooth);

	if (rank==0) printf("Step 6. Compute distance thresholded non-local mean \n");
	int depth = 5;
	float sigsq=0.1;
	int nlm_count=NLM3D(MultiScaleSmooth, Mean, Dist, NonLocalMean, depth, sigsq);

	if (rank==0) printf("Step 7. Threshold for segmentation \n");
	THRESHOLD=50;
	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				if (NonLocalMean(i,j,k) > THRESHOLD) ID(i,j,k) = 0;
				else 								 ID(i,j,k) = 1;

				// intialize distance based on segmentation
				Dist(i,j,k) = 2.0*ID(i,j,k)-1.0;
			}
		}
	}

	if (rank==0) printf("Step 8. Generate final distance function \n");
	// generate a sparse signed distance function
	Eikonal3D(Dist,ID,Dm,nx*nprocx);

	//printf("Non-local means count fraction = %f \n",float(nlm_count)/float(nx*ny*nz));
	/*    for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
			        n = k*nx*ny+j*nx+i;
			        if (Dm.id[n]==char(SOLID))     Dm.id[n] = 0;
			       	else if (Dm.id[n]==char(NWP))  Dm.id[n] = 1;
			       	else                           Dm.id[n] = 2;

			}
		}
	}
	if (rank==0) printf("Domain set \n");
	// Write the local volume files
	char LocalRankString[8];
	char LocalRankFilename[40];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"Seg.%s",LocalRankString);
	FILE * SEG;
	SEG=fopen(LocalRankFilename,"wb");
	fwrite(LOCVOL.get(),4,N,SEG);
	fclose(SEG);
	 */
	if (rank==0) printf("Writing output \n");

	std::vector<IO::MeshDataStruct> meshData(2);
	meshData[0].meshName = "Full domain";
	meshData[0].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(Dm.rank_info,nx-2,ny-2,nz-2,Lx,Ly,Lz) );
	meshData[1].meshName = "Sparse domain";
	meshData[1].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(Dm.rank_info,nsx-2,nsy-2,nsz-2,Lx,Ly,Lz) );

	std::shared_ptr<IO::Variable> OrigData( new IO::Variable() );
	std::shared_ptr<IO::Variable> spMedianData( new IO::Variable() );
	std::shared_ptr<IO::Variable> spSegData( new IO::Variable() );
	std::shared_ptr<IO::Variable> spDistData( new IO::Variable() );
	std::shared_ptr<IO::Variable> DistData( new IO::Variable() );
	std::shared_ptr<IO::Variable> MultiMean( new IO::Variable() );
	std::shared_ptr<IO::Variable> NonLocMean( new IO::Variable() );
	std::shared_ptr<IO::Variable> SegData( new IO::Variable() );

	// Full resolution data
	OrigData->name = "Source Data";
	OrigData->type = IO::VolumeVariable;
	OrigData->dim = 1;
	OrigData->data.resize(nx-2,ny-2,nz-2);
	meshData[0].vars.push_back(OrigData);

	MultiMean->name = "Multiscale Mean";
	MultiMean->type = IO::VolumeVariable;
	MultiMean->dim = 1;
	MultiMean->data.resize(nx-2,ny-2,nz-2);
	meshData[0].vars.push_back(MultiMean);

	NonLocMean->name = "Non-Local Mean";
	NonLocMean->type = IO::VolumeVariable;
	NonLocMean->dim = 1;
	NonLocMean->data.resize(nx-2,ny-2,nz-2);
	meshData[0].vars.push_back(NonLocMean);

	SegData->name = "Segmented Data";
	SegData->type = IO::VolumeVariable;
	SegData->dim = 1;
	SegData->data.resize(nx-2,ny-2,nz-2);
	meshData[0].vars.push_back(SegData);

	DistData->name = "Signed Distance";
	DistData->type = IO::VolumeVariable;
	DistData->dim = 1;
	DistData->data.resize(nx-2,ny-2,nz-2);
	meshData[0].vars.push_back(DistData);
	//..........................................

	// ....... Sparse resolution data .......
	spMedianData->name = "Sparse Median Filter";
	spMedianData->type = IO::VolumeVariable;
	spMedianData->dim = 1;
	spMedianData->data.resize(nsx-2,nsy-2,nsz-2);
	meshData[1].vars.push_back(spMedianData);

	spSegData->name = "Sparse Segmentation";
	spSegData->type = IO::VolumeVariable;
	spSegData->dim = 1;
	spSegData->data.resize(nsx-2,nsy-2,nsz-2);
	meshData[1].vars.push_back(spSegData);

	spDistData->name = "Sparse Distance";
	spDistData->type = IO::VolumeVariable;
	spDistData->dim = 1;
	spDistData->data.resize(nsx-2,nsy-2,nsz-2);
	meshData[1].vars.push_back(spDistData);
	//..........................................

	/*
	 * Only Array<double> works right now :(
	 *
    Array<float>& INPUT = meshData[0].vars[0]->data;
    Array<float>& spMEDIAN = meshData[1].vars[0]->data;
    Array<char>& spSEGMENTED = meshData[1].vars[1]->data;
    Array<float>& spDISTANCE = meshData[1].vars[2]->data;

    fillFloat.copy(LOCVOL,INPUT);
    fillFloat_sp.copy(spM,spMEDIAN);
    fillChar_sp.copy(spID,spSEGMENTED);
    fillFloat_sp.copy(spDist,spDISTANCE);
	 */

	Array<double>& INPUT = meshData[0].vars[0]->data;	
	Array<double>& MULTIMEAN = meshData[0].vars[1]->data;
	Array<double>& NONLOCALMEAN = meshData[0].vars[2]->data;
	Array<double>& SEGMENTED = meshData[0].vars[3]->data;
	Array<double>& DISTANCE = meshData[0].vars[4]->data;

	Array<double>& spMEDIAN = meshData[1].vars[0]->data;
	Array<double>& spSEGMENTED = meshData[1].vars[1]->data;
	Array<double>& spDISTANCE = meshData[1].vars[2]->data;


	// manually change to double and write
	for (k=1;k<nz-1;k++){
		for (j=1;j<ny-1;j++){
			for (i=1;i<nx-1;i++){
				INPUT(i-1,j-1,k-1) = double( LOCVOL(i,j,k));
				SEGMENTED(i-1,j-1,k-1) = double( ID(i,j,k));
				DISTANCE(i-1,j-1,k-1) = double( Dist(i,j,k));
				MULTIMEAN(i-1,j-1,k-1) = double( MultiScaleSmooth(i,j,k));
				NONLOCALMEAN(i-1,j-1,k-1) = double( NonLocalMean(i,j,k));
			}
		}
	}

	for (k=1;k<nsz-1;k++){
		for (j=1;j<nsy-1;j++){
			for (i=1;i<nsx-1;i++){
				spMEDIAN(i-1,j-1,k-1) = double( spM(i,j,k));
				spSEGMENTED(i-1,j-1,k-1) = double( spID(i,j,k));
				spDISTANCE(i-1,j-1,k-1) = double( spDist(i,j,k));
			}
		}
	}

	IO::writeData( 0, meshData, 2, comm );
	if (rank==0) printf("Finished. \n");

	/*    for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
			        n = k*nx*ny+j*nx+i;
			        if (Dm.id[n]==char(SOLID))     Dm.id[n] = 0;
			       	else if (Dm.id[n]==char(NWP))  Dm.id[n] = 1;
			       	else                           Dm.id[n] = 2;

			}
		}
	}
	if (rank==0) printf("Domain set \n");

	// Write the local volume files
	char LocalRankString[8];
	char LocalRankFilename[40];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"Seg.%s",LocalRankString);
	FILE * SEG;
	SEG=fopen(LocalRankFilename,"wb");
	fwrite(LOCVOL.get(),4,N,SEG);
	fclose(SEG);
	 */

	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;
}

