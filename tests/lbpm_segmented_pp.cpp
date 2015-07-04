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

inline void WriteBlobs(TwoPhase Averages){
	printf("Writing the blob list \n");
    FILE *BLOBLOG;
  	BLOBLOG=fopen("blobs.tcat","w");
	fprintf(BLOBLOG,"%.5g %.5g %.5g\n",Averages.vol_w_global,Averages.paw_global,Averages.aws_global);
	for (int b=0; b<(int)Averages.BlobAverages.size(1); b++){
		if (Averages.BlobAverages(0,b) > 0.0){
			double Vn,pn,awn,ans,Jwn,Kwn,lwns,cwns;
			Vn = Averages.BlobAverages(1,b);
			pn = Averages.BlobAverages(2,b);
			awn = Averages.BlobAverages(3,b);
			ans = Averages.BlobAverages(4,b);
			Jwn = Averages.BlobAverages(5,b);
			Kwn = Averages.BlobAverages(6,b);
			lwns = Averages.BlobAverages(7,b);
			cwns = Averages.BlobAverages(8,b);

			fprintf(BLOBLOG,"%.5g ", Vn); //Vn
			fprintf(BLOBLOG,"%.5g ", pn); //pn
			fprintf(BLOBLOG,"%.5g ", awn); //awn
			fprintf(BLOBLOG,"%.5g ", ans); //ans
			fprintf(BLOBLOG,"%.5g ", Jwn); //Jwn
			fprintf(BLOBLOG,"%.5g ", Kwn); //Kwn
			fprintf(BLOBLOG,"%.5g ", lwns); //lwns
			fprintf(BLOBLOG,"%.5g\n",cwns); //cwns
		}
	}
	fclose(BLOBLOG);
}

inline void MeanFilter(DoubleArray &Mesh){
	for (int k=1; k<Mesh.size(2)-1; k++){
		for (int j=1; j<Mesh.size(1)-1; j++){
			for (int i=1; i<Mesh.size(0)-1; i++){
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
	SSO(Averages.SDs,id,Dm,20);

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
	SSO(Averages.Phase,id,Dm,20);

    sprintf(LocalRankFilename,"Phase.%05i",rank);
    FILE *PHASE = fopen(LocalRankFilename,"wb");
    fwrite(Averages.Phase.get(),8,Averages.Phase.length(),PHASE);
    fclose(PHASE);

	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				n=k*nx*ny+j*nx+i;
				Averages.Phase(i,j,k) += 1.0;
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
				// Initialize distance to +/- 1
				// Dilation of the non-wetting phase
				Averages.SDn(i,j,k) = Averages.Phase(i,j,k)+1.0;
				Averages.Phase(i,j,k) = Averages.SDn(i,j,k);
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

	double vF,vS;
	vF = vS = 0.0;

    double beta = 0.95;
	if (rank==0) printf("initializing the system \n");
    Averages.SetupCubes(Dm);
    Averages.UpdateSolid();
    Averages.Initialize();
    Averages.UpdateMeshValues();
    Dm.CommunicateMeshHalo(Averages.Phase);
    Dm.CommunicateMeshHalo(Averages.SDn);

	if (rank==0) printf("computing blobs \n");
 //   int nblobs_global = ComputeGlobalBlobIDs(Dm.Nx-2,Dm.Ny-2,Dm.Nz-2,Dm.rank_info,
  //  		Averages.Phase,Averages.SDs,vF,vS,Averages.BlobLabel);
//	if (Dm.rank==0) printf("Number of blobs is %i \n",nblobs_global);

//     int nblobs_global = ComputeGlobalBlobIDs(Dm.Nx-2,Dm.Ny-2,Dm.Nz-2,Dm.rank_info,
//					 Averages.SDn,Averages.SDs,vF,vS,Averages.BlobLabel);

	if (rank==0) printf("computing local averages  \n");
	Averages.ComputeLocalBlob();
	if (rank==0) printf("reducing averages  \n");
	Averages.Reduce();

	if (rank==0) printf("Writing blobs \n");
    // Write the local blob ids
    sprintf(LocalRankFilename,"BlobLabel.%05i",rank);
    FILE *BLOBLOCAL = fopen(LocalRankFilename,"wb");
    fwrite(Averages.BlobLabel.get(),4,Averages.BlobLabel.length(),BLOBLOCAL);
    fclose(BLOBLOCAL);
    printf("Wrote BlobLabel.%05i \n",rank);

	if (rank==0) printf("Sorting averages \n");
    //  Blobs.Set(Averages.BlobAverages.NBLOBS);
    int dimx = (int)Averages.BlobAverages.size(0);
    int dimy = (int)Averages.BlobAverages.size(1);
    int TotalBlobInfoSize=dimx*dimy;

    //      BlobContainer Blobs;
    DoubleArray RecvBuffer(dimx);
    //    MPI_Allreduce(&Averages.BlobAverages.get(),&Blobs.get(),1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("Number of components is %i \n",dimy);

    for (int b=0; b<dimy; b++){

    	MPI_Allreduce(&Averages.BlobAverages(0,b),&RecvBuffer(0),dimx,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    	for (int idx=0; idx<dimx-1; idx++) Averages.BlobAverages(idx,b)=RecvBuffer(idx);
    	MPI_Barrier(MPI_COMM_WORLD);

    	if (Averages.BlobAverages(0,b) > 0.0){
    		double Vn,pn,awn,ans,Jwn,Kwn,lwns,cwns,trawn,trJwn;
    		Vn = Averages.BlobAverages(1,b);
    		pn = Averages.BlobAverages(2,b)/Averages.BlobAverages(0,b);
    		awn = Averages.BlobAverages(3,b);
    		ans = Averages.BlobAverages(4,b);
    		if (awn != 0.0){
    			Jwn = Averages.BlobAverages(5,b)/Averages.BlobAverages(3,b);
    			Kwn = Averages.BlobAverages(6,b)/Averages.BlobAverages(3,b);
    		}
    		else Jwn=Kwn=0.0;

    		trawn = Averages.BlobAverages(12,b);
    		if (trawn != 0.0){
    			trJwn = Averages.BlobAverages(13,b)/trawn;
    		}
    		else trJwn=0.0;

    		lwns = Averages.BlobAverages(7,b);
    		if (lwns != 0.0) cwns = Averages.BlobAverages(8,b)/Averages.BlobAverages(7,b);
    		else  cwns=0.0;
    		Averages.BlobAverages(2,b) = pn;
    		Averages.BlobAverages(5,b) = trJwn;
    		Averages.BlobAverages(6,b) = Kwn;
    		Averages.BlobAverages(8,b) = cwns;
    		//	Averages.BlobAverages(13,b) = trJwn;
    	}
    }

    if (rank==0) printf("Sorting blobs by volume \n");
    Averages.SortBlobs();

    if (rank==0)   WriteBlobs(Averages);

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;

}
