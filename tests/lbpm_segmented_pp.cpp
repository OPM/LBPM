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
#include <TwoPhase.h>

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

	char LocalRankFilename[40];

    int N = (nx+2)*(ny+2)*(nz+2);
	Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
	// Read the phase ID
    sprintf(LocalRankFilename,"ID.%05i",rank);
    FILE *ID = fopen(LocalRankFilename,"rb");
    fread(Dm.id,1,N,ID);
    fclose(ID);
    // Initialize the domain and communication
    Dm.CommInit(MPI_COMM_WORLD);

	nx+=2; ny+=2; nz+=2;
	int count = 0;
	N=nx*ny*nz;

	char *id;
	id = new char [N];
	TwoPhase Averages(Dm);
//	DoubleArray Distance(nx,ny,nz);
//	DoubleArray Phase(nx,ny,nz);

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
				Averages.SDs(i,j,k) = 2.0*id[n]-1.0;
			}
		}
	}
	MeanFilter(Averages.SDs);

	if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
	SSO(Averages.SDs,id,Dm,25);

    sprintf(LocalRankFilename,"SignDist.%05i",rank);
    FILE *DIST = fopen(LocalRankFilename,"wb");
    fwrite(Averages.SDs.get(),8,Averages.SDs.length(),DIST);
    fclose(DIST);

	// Solve for the position of the non-wetting phase
	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				n = k*nx*ny+j*nx+i;
				// Initialize the non-wetting phase
				if (Dm.id[n] == 1)	id[n] = 1;
				else	       		id[n] = 0;
			}
		}
	}
	// Initialize the signed distance function
	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				n=k*nx*ny+j*nx+i;
				// Initialize distance to +/- 1
				Averages.Phase(i,j,k) = 2.0*id[n]-1.0;
			}
		}
	}
	MeanFilter(Averages.Phase);

	if (rank==0) printf("Initialized non-wetting phase -- Converting to Signed Distance function \n");
	SSO(Averages.Phase,id,Dm,25);

	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				n=k*nx*ny+j*nx+i;
				Averages.Phase(i,j,k) -= 1.0;
				// Initialize distance to +/- 1
				// Dilation of the non-wetting phase
				Averages.SDn(i,j,k) = -Averages.Phase(i,j,k);
				Averages.Phase(i,j,k) = Averages.SDn(i,j,k);
				Averages.Phase_tplus(i,j,k) = Averages.SDn(i,j,k);
				Averages.Phase_tminus(i,j,k) = Averages.SDn(i,j,k);
				Averages.DelPhi(i,j,k) = 0.0;
				Averages.Press(i,j,k) = 0.0;
				Averages.Vel_x(i,j,k) = 0.0;
				Averages.Vel_y(i,j,k) = 0.0;
				Averages.Vel_z(i,j,k) = 0.0;
				if (Averages.SDs(i,j,k) > 0.0){
					if (Averages.Phase(i,j,k) > 0.0){
						Dm.id[n] = 2;
					}
					else{
						Dm.id[n] = 1;
					}
				}
				else{
					Dm.id[n] = 0;
				}
			}
		}
	}

    // Create the MeshDataStruct
    fillHalo<double> fillData(Dm.rank_info,Nx-2,Ny-2,Nz-2,1,1,1,0,1);
    std::vector<IO::MeshDataStruct> meshData(1);
    meshData[0].meshName = "domain";
    meshData[0].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(Dm.rank_info,Nx-2,Ny-2,Nz-2,Lx,Ly,Lz) );
    std::shared_ptr<IO::Variable> PhaseVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> SignDistVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> BlobIDVar( new IO::Variable() );
    PhaseVar->name = "phase";
    PhaseVar->type = IO::VolumeVariable;
    PhaseVar->dim = 1;
    PhaseVar->data.resize(Nx-2,Ny-2,Nz-2);
    meshData[0].vars.push_back(PhaseVar);
    SignDistVar->name = "SignDist";
    SignDistVar->type = IO::VolumeVariable;
    SignDistVar->dim = 1;
    SignDistVar->data.resize(Nx-2,Ny-2,Nz-2);
    meshData[0].vars.push_back(SignDistVar);
    BlobIDVar->name = "BlobID";
    BlobIDVar->type = IO::VolumeVariable;
    BlobIDVar->dim = 1;
    BlobIDVar->data.resize(Nx-2,Ny-2,Nz-2);
    meshData[0].vars.push_back(BlobIDVar);
    
    fillData.copy(Averages.SDn,PhaseVar->data);
    fillData.copy(Averages.SDs,SignDistVar->data);
    fillData.copy(Averages.Label_NWP,BlobIDVar->data);
    IO::writeData( 0, meshData, 2 );
    
 //   sprintf(LocalRankFilename,"Phase.%05i",rank);
  //  FILE *PHASE = fopen(LocalRankFilename,"wb");
  //  fwrite(Averages.Phase.get(),8,Averages.Phase.length(),PHASE);
  //  fclose(PHASE);

	double vF,vS;
	vF = vS = 0.0;

    double beta = 0.95;
	if (rank==0) printf("initializing the system \n");
    Averages.UpdateSolid();
    Averages.UpdateMeshValues();
    Dm.CommunicateMeshHalo(Averages.Phase);
    Dm.CommunicateMeshHalo(Averages.SDn);
    Dm.CommunicateMeshHalo(Averages.SDs);

	int timestep=5;
	Averages.Initialize();
	if (rank==0) printf("computing phase components \n");
	Averages.ComponentAverages();
	if (rank==0) printf("sorting phase components \n");
	Averages.SortBlobs();
	Averages.PrintComponents(timestep);

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;

}
