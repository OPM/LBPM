// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "common/Communication.h"
#include "analysis/analysis.h"
#ifdef PROFILE
	#include "ProfilerApp.h"
#endif
#include "TwoPhase.h"

//#include "Domain.h"

using namespace std;

void readRankData( int proc, int nx, int ny, int nz, DoubleArray& Phase, DoubleArray& SignDist )
{
    Phase.resize(nx,ny,nz);
    SignDist.resize(nx,ny,nz);
    char file1[40], file2[40];
    sprintf(file1,"SignDist.%05d",proc);
    sprintf(file2,"Phase.%05d",proc);
    ReadBinaryFile(file1, Phase.get(), nx*ny*nz);
    ReadBinaryFile(file2, SignDist.get(), nx*ny*nz);
}

inline void  WriteBlobStates(TwoPhase TCAT, double D, double porosity){
	int a;
	double iVol=1.0/TCAT.Dm.Volume;
	double PoreVolume;
	double nwp_volume,vol_n,pan,pn,pw,pawn,pwn,awn,ans,aws,Jwn,Kwn,lwns,cwns,clwns;
	double sw,awnD,awsD,ansD,lwnsDD,JwnD,pc;
	nwp_volume=vol_n=pan=awn=ans=Jwn=Kwn=lwns=clwns=pawn=0.0;
	sw = TCAT.sat_w;
	pw = TCAT.paw_global;
	aws = TCAT.aws;
	// Compute the averages over the entire non-wetting phase
	printf("Writing blobstates.tcat for %i components \n",TCAT.nblobs_global);
	FILE *BLOBSTATES;
	BLOBSTATES = fopen("./blobstates.tcat","w");
	if (BLOBSTATES==NULL) ERROR("Cannot open blobstates.tcat for writing");
	for (a=0; a<TCAT.nblobs_global; a++){
		vol_n += TCAT.BlobAverages(0,a);
		pan += TCAT.BlobAverages(2,a)*TCAT.BlobAverages(0,a);
		awn += TCAT.BlobAverages(3,a);
		ans += TCAT.BlobAverages(4,a);
		Jwn += TCAT.BlobAverages(5,a)*TCAT.BlobAverages(3,a);
		Kwn += TCAT.BlobAverages(6,a)*TCAT.BlobAverages(3,a);
		lwns += TCAT.BlobAverages(7,a);
		clwns += TCAT.BlobAverages(8,a)*TCAT.BlobAverages(7,a);
		nwp_volume += TCAT.BlobAverages(1,a);
		pawn += TCAT.BlobAverages(2,a)*TCAT.BlobAverages(3,a);
	}

	// Compute the pore voume (sum of wetting an non-wetting phase volumes)
	PoreVolume=TCAT.wp_volume_global + nwp_volume;
	// Subtract off portions of non-wetting phase in order of size
	for (a=TCAT.nblobs_global-1; a>0; a--){
		// Subtract the features one-by-one
		vol_n -= TCAT.BlobAverages(0,a);
		pan -= TCAT.BlobAverages(2,a)*TCAT.BlobAverages(0,a);
		awn -= TCAT.BlobAverages(3,a);
		ans -= TCAT.BlobAverages(4,a);
		Jwn -= TCAT.BlobAverages(5,a)*TCAT.BlobAverages(3,a);
		Kwn -= TCAT.BlobAverages(6,a)*TCAT.BlobAverages(3,a);
		lwns -= TCAT.BlobAverages(7,a);
		clwns -= TCAT.BlobAverages(8,a)*TCAT.BlobAverages(7,a);
		nwp_volume -= TCAT.BlobAverages(1,a);
		pawn -= TCAT.BlobAverages(2,a)*TCAT.BlobAverages(3,a);

		// Update wetting phase averages
		aws += TCAT.BlobAverages(4,a);
		if (vol_n > 64){	// Only consider systems with "large enough" blobs -- 4^3
			if (fabs(1.0 - nwp_volume/PoreVolume - sw) > 0.005 || a == 1){
				sw = 1.0 - nwp_volume/PoreVolume;

				JwnD = Jwn*D/awn;
				//trJwnD = -trJwn*D/trawn;
				cwns = clwns / lwns;
				pwn = (pawn/awn-pw)*D/0.058;
				pn = pan/vol_n;
				awnD = awn*D*iVol;
				awsD = aws*D*iVol;
				ansD = ans*D*iVol;
				lwnsDD = lwns*D*D*iVol;
				pc = (pn-pw)*D/0.058; 	// hard-coded surface tension due to being lazy

				fprintf(BLOBSTATES,"%.5g %.5g %.5g ",sw,pn,pw);
				fprintf(BLOBSTATES,"%.5g %.5g %.5g %.5g ",awnD,awsD,ansD,lwnsDD);
				fprintf(BLOBSTATES,"%.5g %.5g %.5g %.5g %i\n",pc,pwn,JwnD,cwns,a);
			}
		}
	}
	fclose(BLOBSTATES);
}

int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    if ( rank==0 ) {
        printf("-----------------------------------------------------------\n");
        printf("Labeling Blobs from Two-Phase Lattice Boltzmann Simulation \n");
        printf("-----------------------------------------------------------\n");
    }

    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
    double Lx, Ly, Lz;
    int Nx,Ny,Nz;
    int i,j,k,n;

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
    if ( nprocs < nprocx*nprocy*nprocz )
        ERROR("Insufficient number of processors");

    // Set up the domain
	int BC=0;
    // Get the rank info
	Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
 //   const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
	TwoPhase Averages(Dm);
	int N = (nx+2)*(ny+2)*(nz+2);
	Nx = nx+2;
	Ny = ny+2;
	Nz = nz+2;
	// Read in sphere pack (initialize the non-wetting phase as inside of spheres)
	if (rank==1) printf("nspheres =%i \n",nspheres);
	//.......................................................................
	double *cx,*cy,*cz,*rad;
	cx = new double[nspheres];
	cy = new double[nspheres];
	cz = new double[nspheres];
	rad = new double[nspheres];
	//.......................................................................
	if (rank == 0)	printf("Reading the sphere packing \n");
	if (rank == 0)	ReadSpherePacking(nspheres,cx,cy,cz,rad);
	MPI_Barrier(MPI_COMM_WORLD);
	// Broadcast the sphere packing to all processes
	MPI_Bcast(cx,nspheres,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(cy,nspheres,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(cz,nspheres,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(rad,nspheres,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//...........................................................................
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) cout << "Domain set." << endl;
	//.......................................................................
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				Dm.id[n] = 1;
			}
		}
	}
	//.......................................................................
    Dm.CommInit(MPI_COMM_WORLD); // Initialize communications for domains
	//.......................................................................

	//.......................................................................
	SignedDistance(Averages.Phase.get(),nspheres,cx,cy,cz,rad,Lx,Ly,Lz,Nx,Ny,Nz,
					   Dm.iproc,Dm.jproc,Dm.kproc,Dm.nprocx,Dm.nprocy,Dm.nprocz);
	//.......................................................................
	// Assign the phase ID field based on the signed distance
	//.......................................................................
	if (rank==0) printf("Initializing the system \n");
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				// Shrink the sphere sizes by two voxels to make sure they don't touch
				Averages.SDs(i,j,k) = 100.0;
				Averages.Phase(i,j,k) -= 5.0;
				if (Averages.Phase(i,j,k) > 0.0){
					Dm.id[n] = 2;
				}
				else{
					Dm.id[n] = 1;
				}
				Averages.SDn(i,j,k) = Averages.Phase(i,j,k);
				Averages.Phase_tplus(i,j,k) = Averages.Phase(i,j,k);
				Averages.Phase_tminus(i,j,k) = Averages.Phase(i,j,k);
				Averages.Press(i,j,k) = 0.0;
				Averages.Vel_x(i,j,k) = 0.0;
				Averages.Vel_y(i,j,k) = 0.0;
				Averages.Vel_z(i,j,k) = 0.0;
			}
		}
	}

	if (rank==0) printf("Computing averages \n");
    double beta = 0.95;
    Averages.SetupCubes(Dm);
    Averages.UpdateSolid();
	if (rank==0) printf("initializing the system \n");
    Averages.Initialize();
	if (rank==0) printf("updating mesh halos \n");
    Averages.UpdateMeshValues();
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

    FILE *BLOBLOG;
    if (rank==0){
    	printf("Writing the blob list \n");
      	BLOBLOG=fopen("blobs.tcat","w");

    	//	printf("Reduced blob %i \n",b);
    	fprintf(BLOBLOG,"%.5g %.5g %.5g\n",Averages.vol_w_global,Averages.paw_global,Averages.aws_global);
    	for (int b=0; b<dimy; b++){
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

    if (rank==0) {
      int a;
      double D=1.0;
	double iVol=1.0/Averages.Dm.Volume;
	double PoreVolume;
	double nwp_volume,vol_n,pan,pn,pw,pawn,pwn,awn,ans,aws,Jwn,Kwn,lwns,cwns,clwns;
	double sw,awnD,awsD,ansD,lwnsDD,JwnD,pc;
	nwp_volume=vol_n=pan=awn=ans=Jwn=Kwn=lwns=clwns=pawn=0.0;
	sw = Averages.sat_w;
	pw = Averages.paw_global;
	aws = Averages.aws;
	// Compute the averages over the entire non-wetting phase
	printf("Writing blobstates.tcat for %i components \n",Averages.nblobs_global);
	FILE *BLOBSTATES;
	BLOBSTATES = fopen("./blobstates.tcat","w");
	if (BLOBSTATES==NULL) ERROR("Cannot open blobstates.tcat for writing");
	for (a=0; a<Averages.nblobs_global; a++){
		vol_n += Averages.BlobAverages(0,a);
		pan += Averages.BlobAverages(2,a)*Averages.BlobAverages(0,a);
		awn += Averages.BlobAverages(3,a);
		ans += Averages.BlobAverages(4,a);
		Jwn += Averages.BlobAverages(5,a)*Averages.BlobAverages(3,a);
		Kwn += Averages.BlobAverages(6,a)*Averages.BlobAverages(3,a);
		lwns += Averages.BlobAverages(7,a);
		clwns += Averages.BlobAverages(8,a)*Averages.BlobAverages(7,a);
		nwp_volume += Averages.BlobAverages(1,a);
		pawn += Averages.BlobAverages(2,a)*Averages.BlobAverages(3,a);
	}

	// Compute the pore voume (sum of wetting an non-wetting phase volumes)
	PoreVolume=Averages.wp_volume_global + nwp_volume;
	// Subtract off portions of non-wetting phase in order of size
	for (a=Averages.nblobs_global-1; a>0; a--){
		// Subtract the features one-by-one
		vol_n -= Averages.BlobAverages(0,a);
		pan -= Averages.BlobAverages(2,a)*Averages.BlobAverages(0,a);
		awn -= Averages.BlobAverages(3,a);
		ans -= Averages.BlobAverages(4,a);
		Jwn -= Averages.BlobAverages(5,a)*Averages.BlobAverages(3,a);
		Kwn -= Averages.BlobAverages(6,a)*Averages.BlobAverages(3,a);
		lwns -= Averages.BlobAverages(7,a);
		clwns -= Averages.BlobAverages(8,a)*Averages.BlobAverages(7,a);
		nwp_volume -= Averages.BlobAverages(1,a);
		pawn -= Averages.BlobAverages(2,a)*Averages.BlobAverages(3,a);

		// Update wetting phase averages
		aws += Averages.BlobAverages(4,a);
		if (vol_n > 64){	// Only consider systems with "large enough" blobs -- 4^3
			if (fabs(1.0 - nwp_volume/PoreVolume - sw) > 0.005 || a == 1){
				sw = 1.0 - nwp_volume/PoreVolume;

				JwnD = Jwn*D/awn;
				//trJwnD = -trJwn*D/trawn;
				cwns = clwns / lwns;
				pwn = (pawn/awn-pw)*D/0.058;
				pn = pan/vol_n;
				awnD = awn*D*iVol;
				awsD = aws*D*iVol;
				ansD = ans*D*iVol;
				lwnsDD = lwns*D*D*iVol;
				pc = (pn-pw)*D/0.058; 	// hard-coded surface tension due to being lazy

				fprintf(BLOBSTATES,"%.5g %.5g %.5g ",sw,pn,pw);
				fprintf(BLOBSTATES,"%.5g %.5g %.5g %.5g ",awnD,awsD,ansD,lwnsDD);
				fprintf(BLOBSTATES,"%.5g %.5g %.5g %.5g %i\n",pc,pwn,JwnD,cwns,a);
			}
		}
	}
	fclose(BLOBSTATES);

    }

    //WriteBlobStates(Averages,Length,porosity);

    /*FILE *BLOBS = fopen("Blobs.dat","wb");
    fwrite(GlobalBlobID.get(),4,Nx*Ny*Nz,BLOBS);
    fclose(BLOBS);*/

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;  
}

