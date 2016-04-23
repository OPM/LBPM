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
#include "IO/netcdf.h"

#include "ProfilerApp.h"

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


	PROFILE_START("ReadVolume");

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

	nx+=2; ny+=2; nz+=2;
	N=nx*ny*nz;

	if (rank==0) printf("All sub-domains recieved \n");

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

    std::vector<IO::MeshDataStruct> meshData(2);
    meshData[0].meshName = "Full domain";
    meshData[0].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(Dm.rank_info,nx-2,ny-2,nz-2,Lx,Ly,Lz) );
    meshData[1].meshName = "Sparse domain";
    meshData[1].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(Dm.rank_info,nsx-2,nsy-2,nsz-2,Lx,Ly,Lz) );

    std::shared_ptr<IO::Variable> OrigData( new IO::Variable() );
    std::shared_ptr<IO::Variable> spMedianData( new IO::Variable() );
    std::shared_ptr<IO::Variable> spSegData( new IO::Variable() );
    std::shared_ptr<IO::Variable> spDistData( new IO::Variable() );

    // Full resolution data
    OrigData->name = "Source Data";
    OrigData->type = IO::VolumeVariable;
    OrigData->dim = 1;
    OrigData->data.resize(nx-2,ny-2,nz-2);
    meshData[0].vars.push_back(OrigData);
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
    Array<double>& spMEDIAN = meshData[1].vars[0]->data;
    Array<double>& spSEGMENTED = meshData[1].vars[1]->data;
    Array<double>& spDISTANCE = meshData[1].vars[2]->data;

    // manually change to double and write
    for (k=0;k<nz;k++){
    	for (j=0;j<ny;j++){
    		for (i=0;i<nx;i++){
    			INPUT(i,j,k) = double( LOCVOL(i,j,k));
    		}
    	}
    }

    for (k=0;k<nsz;k++){
    	for (j=0;j<nsy;j++){
    		for (i=0;i<nsx;i++){
    			spMEDIAN(i,j,k) = double( spM(i,j,k));
    			spSEGMENTED(i,j,k) = double( spID(i,j,k));
    			spDISTANCE(i,j,k) = double( spDist(i,j,k));
    		}
    	}
    }

    IO::writeData( 0, meshData, 2, comm );
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
*/
	// Write the local volume files
	char LocalRankString[8];
	char LocalRankFilename[40];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"Seg.%s",LocalRankString);
	FILE * SEG;
	SEG=fopen(LocalRankFilename,"wb");
	fwrite(LOCVOL.get(),4,N,SEG);
	fclose(SEG);
    
    MPI_Barrier(comm);
    MPI_Finalize();
	return 0;
}

