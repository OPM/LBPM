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
#include <Array.h>
#include <Domain.h>

int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
    double Lx, Ly, Lz;
    int Nx,Ny,Nz;
    int i,j,k,n;
	int BC=0;
    char Filename[40];
    int xStart,yStart,zStart;
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

    	ifstream image("Segmented.in");
    	image >> Filename; 	// Name of data file containing segmented data
    	image >> Nx;   		// size of the binary file
    	image >> Ny;
    	image >> Nz;
    	image >> xStart;	// offset for the starting voxel
    	image >> yStart;
    	image >> zStart;

    	getline(image,line);
    	std::istringstream solidLine(line);
    	while (solidLine >> n)	solidValues.push_back(n);
    	printf("Read %i solid values \n",solidValues.size());

    	getline(image,line);
    	std::istringstream nwpLine(line);
    	while (nwpLine >> n)	nwpValues.push_back(n);
    	printf("Read %i nwp values \n",nwpValues.size());


//    	image >> solidValue; // value assigned to the solid phase
//   	image >> fluidValue; // value assigned to the non-wetting phase
    }
	MPI_Barrier(MPI_COMM_WORLD);
	// Computational domain
	MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocy,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nspheres,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//.................................................
	MPI_Barrier(MPI_COMM_WORLD);

    // Check that the number of processors >= the number of ranks
    if ( rank==0 ) {
        printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
        printf("Number of MPI ranks used: %i \n", nprocs);
        printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
    }
    if ( nprocs < nprocx*nprocy*nprocz ){
        ERROR("Insufficient number of processors");
    }
    char *SegData;
    SegData = new char[Nx*Ny*Nz];
    // Rank=0 reads the entire segmented data and distributes to worker processes
    if (rank==0){
    	FILE *SEGDAT = fopen(Filename,"rb");
    	if (SEGDAT==NULL) ERROR("Error reading segmented data");
    	fread(SegData,1,Nx*Ny*Nz,SEGDAT);
    	fclose(SEGDAT);
        printf("Read segmented data from %s \n",Filename);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Get the rank info
    int N = (nx+2)*(ny+2)*(nz+2);
	Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				n = k*nx*ny+j*nx+i;
				Dm.id[n] = 1;
			}
		}
	}
	// Set up the sub-domains
	if (rank==0){
		char *tmp;
		tmp = new char[(nx+2)*(ny+2)*(nz+2)];
		for (int kp=0; kp<Dm.kproc; kp++){
			for (int jp=0; jp<Dm.jproc; jp++){
				for (int ip=0; ip<Dm.iproc; ip++){
					// rank of the process that gets this subdomain
					int rnk = kp*Dm.nprocx*Dm.nprocy + jp*Dm.nprocx + ip;
					// Pack and send the subdomain for rnk
					for (k=0;k<nz+2;k++){
						for (j=0;j<ny+2;j++){
							for (i=0;i<nx+2;j++){
								int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
								int nglobal = zStart*Nx*Ny + yStart*Nx + xStart +
												kp*nprocx*nprocy*nx*ny*nz+ jp*nprocx*nx*ny*nz
												 + ip*nx*ny*nz + k*nx*ny + j*nx + i;
								tmp[nlocal] = SegData[nglobal];
							}
						}
					}
					MPI_Send(&tmp,N,MPI_CHAR,rnk,15,MPI_COMM_WORLD);
				}
			}
		}
	}
	else{
		// Recieve the subdomain from rank = 0
		MPI_Recv(&Dm.id,N,MPI_CHAR,0,15,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	Dm.CommInit(MPI_COMM_WORLD);

	nx+=2; ny+=2; nz+=2;
	int count = 0;
	N=nx*ny*nz;

	char *id;
	id = new char [N];
	for (k=1;k<nz-1;k++){
		for (j=1;j<ny-1;j++){
			for (i=1;i<nx-1;i++){
				n = k*nx*ny+j*nx+i;
				bool solid=false;
				for (int idx=0; idx<solidValues.size(); idx++){
					if (Dm.id[n] == solidValues(idx)) solid=true;
				}
				if (solid==true)			id[n] = 0;
				else						id[n] = 1;
			}
		}
	}
	
	DoubleArray Distance(nx,ny,nz);
	// Initialize the signed distance function
	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				n=k*nx*ny+j*nx+i;
				// Initialize distance to +/- 1
				Distance(i,j,k) = 1.0*id[n]-0.5;
			}
		}
	}
	if (rank==0) printf("Nx = %i \n",(int)Distance.size(0));
	if (rank==0) printf("Ny = %i \n",(int)Distance.size(1));
	if (rank==0) printf("Nz = %i \n",(int)Distance.size(2));

	printf("Initialized! Converting to Signed Distance function \n");
	SSO(Distance,id,Dm,10);

	char LocalRankFilename[40];

    sprintf(LocalRankFilename,"Dist.%05i",rank);
    FILE *DIST = fopen(LocalRankFilename,"wb");
    fwrite(Distance.get(),8,Distance.length(),DIST);
    fclose(DIST);

    int symrank,sympz;
    sympz = 2*nprocz - Dm.kproc;
    symrank = sympz*nprocx*nprocy + Dm.jproc*nprocx + Dm.iproc;

    DoubleArray SymDist(nx,ny,nz);
	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				SymDist(i,j,nz-k-1)=Distance(i,j,k);
			}
		}
    }

    sprintf(LocalRankFilename,"Dist.%05i",symrank);
    FILE *SYMDIST = fopen(LocalRankFilename,"wb");
    fwrite(SymDist.get(),8,SymDist.length(),SYMDIST);
    fclose(SYMDIST);

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;

}
