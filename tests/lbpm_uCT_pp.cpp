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
#include "IO/netcdf.h"

#include "ProfilerApp.h"

inline void Med3D(Array<float> &Input, Array<float> &Output){
	// Perform a 3D Median filter on Input array with specified width
	int i,j,k,ii,jj,kk;
	int imin,jmin,kmin,imax,jmax,kmax;

	float *List;
	List=new float[9];

	int Nx = int(Input.size(0));
	int Ny = int(Input.size(1));
	int Nz = int(Input.size(2));

	for (k=1; k<Nz-1; k++){
		for (j=1; Ny-1; j++){
			for (i=1; Nz-1; i++){

				// Just use a 3x3x3 window (hit recursively if needed)
				imin = i-1;
				jmin = j-1;
				kmin = k-1;
				imax = i+1;
				jmax = j+1;
				kmax = k+1;

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
				for (ii=0; ii<5; ii++){
					for (jj=ii+1; jj<9; jj++){
						if (List[jj] < List[ii]){
							float tmp = List[ii];
							List[ii] = List[jj];
							List[jj] = tmp;
						}
					}
				}
				// Return the median
				Output(i,j,k) = List[4];
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
			for (i=0; nz; i++){

				x = i*hx;
				y = j*hy;
				z = k*hz;

				ii = floor(x);
				jj = floor(y);
				kk = floor(z);

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
	  printf("Time=%i, Max variation=%f, Global variation=%f \n",count,GlobalMax,GlobalVar);

	if (fabs(GlobalMax) < 1e-5){
	  if (Dm.rank==0) printf("Exiting with max tolerance of 1e-5 \n");
	  count=timesteps;
	}
    }
    return GlobalVar;
}


inline void NLM3D(Array<float> &Input, Array<float> &Mean, Array<float> &Output,
					const int d, const float h){
	// Implemenation of 3D non-local means filter
	// 		d determines the width of the search volume
	// 		h is a free parameter for non-local means (i.e. 1/sigma^2)
	float weight, sum;
	int i,j,k,ii,jj,kk;
	int imin,jmin,kmin,imax,jmax,kmax;

	int Nx = int(Input.size(0));
	int Ny = int(Input.size(1));
	int Nz = int(Input.size(2));

	// Compute the local means
	for (k=1; k<Nz-1; k++){
		for (j=1; Ny-1; j++){
			for (i=1; Nz-1; i++){

				imin = max(0,i-d);
				jmin = max(0,j-d);
				kmin = max(0,k-d);
				imax = max(Nx-1,i+d);
				jmax = max(Ny-1,j+d);
				kmax = max(Nz-1,k+d);

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
		for (j=1; Ny-1; j++){
			for (i=1; Nz-1; i++){

				imin = max(0,i-d);
				jmin = max(0,j-d);
				kmin = max(0,k-d);
				imax = max(Nx-1,i+d);
				jmax = max(Ny-1,j+d);
				kmax = max(Nz-1,k+d);

				// compute the expensive non-local means
				sum = 0; weight=0;
				for (kk=kmin; kk<kmax; kk++){
					for (jj=jmin; jj<jmax; jj++){
						for (ii=imin; ii<imax; ii++){
							float tmp = Mean(i,j,k) - Mean(ii,jj,kk);
							sum += exp(-tmp*tmp*h)*Input(ii,jj,kk);
							weight += exp(-tmp*tmp*h);
						}
					}
				}
				Output(i,j,k) = sum / weight;
			}
		}
	}

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
	if ( argc==0 ) {
		printf("At least one filename must be specified\n");
		return 1;
	}
	else {
		filename=std::string(argv[1]);
		printf("Input data file: %s\n",filename.c_str());
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
/*	//.................................................	
	MPI_Bcast(&Ny,1,MPI_INT,0,comm);
	MPI_Bcast(&Ny,1,MPI_INT,0,comm);
	MPI_Bcast(&Nz,1,MPI_INT,0,comm);
	MPI_Bcast(&xStart,1,MPI_INT,0,comm);
	MPI_Bcast(&yStart,1,MPI_INT,0,comm);
	MPI_Bcast(&zStart,1,MPI_INT,0,comm);
*/	//.................................................
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

	// Allocate local arrays for every MPI rank
    Array<float> LOCVOL(nx+2,ny+2,nz+2);

	int xStart,yStart,zStart;
	xStart=Nx/2;
	yStart=Ny/2;
	zStart=Nz/2;

	// Read the volume file and distribute to all processes
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
							printf("Sending data to process %i \n", rnk);
							MPI_Send(tmp,N,MPI_FLOAT,rnk,15,comm);
						}
					}
				}
			}
		}
		else{
			// Recieve the subdomain from rank = 0
			printf("Ready to recieve data %i at process %i \n", N,rank);
			MPI_Recv(LOCVOL.get(),N,MPI_FLOAT,0,15,comm,MPI_STATUS_IGNORE);
		}
		MPI_Barrier(comm);
	}
	nx+=2; ny+=2; nz+=2;
	N=nx*ny*nz;


	if (rank==0) printf("All sub-domains recieved \n");

	// Filter the volume in distributed memory
	// Create sparse structures to make an initial guess, then refine

	int nsx,nsy,nsz;
	nsx=nx/8; nsy=ny/8; nsz=nz/8;

    fillHalo<float> fillFloat(Dm.Comm, Dm.rank_info,nx-2,ny-2,nz-2,1,1,1,0,1);
    fillHalo<char> fillChar(Dm.Comm, Dm.rank_info,nx-2,ny-2,nz-2,1,1,1,0,1);
    fillHalo<float> fillFloat_sp(Dm.Comm, Dm.rank_info,nsx-2,nsy-2,nsz-2,1,1,1,0,1);
    fillHalo<char>  fillChar_sp(Dm.Comm, Dm.rank_info,nsx-2,nsy-2,nsz-2,1,1,1,0,1);


	Array<float> spLOCVOL(nsx,nsy,nsz);	// this holds sparse original data
	Array<float> spM(nsx,nsy,nsz); 		// this holds sparse median filter
	Array<float> spDist(nsx,nsy,nsz);		// this holds sparse signed distance

	// sparse phase ID (segmented values)
	Array<char> spID(nsx,nsy,nsz);

	// Sparsify the the mesh using a stride of 8
	Sparsify(LOCVOL,spLOCVOL);

	// Compute the median filter on the sparse array
	Med3D(spLOCVOL,spM);

	// quick & dirty sparse segmentation
	// this should be replaced
	//    (should use automated mixture model to approximate histograms)
	float THRESHOLD=50;
	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				if (spM(i,j,k) > THRESHOLD) spID(i,j,k) = 0;
				else 						spID(i,j,k) = 1;

				// intialize distance based on segmentation
				spDist(i,j,k) = 2.0*spID(i,j,k)-1.0;
			}
		}
	}

	// generate a sparse signed distance function
	Eikonal3D(spDist,spID,Dm,nsx*nprocx);



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

    std::vector<IO::MeshDataStruct>& visData;
    ASSERT(visData[0].vars[0]->name=="Source Data");
    ASSERT(visData[0].vars[1]->name=="Sparse Med3");
    ASSERT(visData[0].vars[2]->name=="Sparse Segmentation");
    ASSERT(visData[0].vars[3]->name=="Sparse Distance");

    Array<float>& INPUT = visData[0].vars[0]->data;
    Array<float>& spMEDIAN = visData[0].vars[1]->data;
    Array<char>& spSEGMENTED  = visData[0].vars[2]->data;
    Array<float>& spDISTANCE = visData[0].vars[3]->data;

    fillFloat.copy(LOCVOL,INPUT);
    fillFloat_sp.copy(spM,spMEDIAN);
    fillChar_sp.copy(spID,spSEGMENTED);
    fillFloat_sp.copy(spDist,spDISTANCE);

    IO::writeData( 0, visData, 2, comm );
    
    MPI_Barrier(comm);
    MPI_Finalize();
	return 0;
}

