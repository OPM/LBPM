// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "common/pmmc.h"
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "ProfilerApp.h"
#include "TwoPhase.h"

//#include "Domain.h"

using namespace std;


inline void ReadBinaryFile(char *FILENAME, double *Data, int N)
{
    int n;
    double value;
    ifstream File(FILENAME,ios::binary);
    for (n=0; n<N; n++){
        // Write the two density values
        File.read((char*) &value, sizeof(value));
        Data[n] = value;

    }
    File.close();
}


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
	FILE *BLOBSTATES= fopen("blobstates.tcat","w");
	int a;
	double iVol=1.0/TCAT.Dm.Volume;
	double PoreVolume;
	double nwp_volume,vol_n,pan,pn,pw,pawn,pwn,awn,ans,aws,Jwn,Kwn,lwns,cwns,clwns;
	double sw,awnD,awsD,ansD,lwnsDD,JwnD,pc;
	nwp_volume=vol_n=pan=awn=ans=Jwn=Kwn=lwns=clwns=pawn=0.0;
	pw = TCAT.paw_global;
	aws = TCAT.aws;
	// Compute the averages over the entire non-wetting phsae
	for (a=0; a<(int)TCAT.BlobAverages.size(1); a++){
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
	for (a=TCAT.BlobAverages.size(1)-1; a>0; a--){
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
    PROFILE_ENABLE(0);
    PROFILE_DISABLE_TRACE();
    PROFILE_SYNCHRONIZE();
    PROFILE_START("main");

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

    // Check that the number of processors >= the number of ranks
    if ( rank==0 ) {
        printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
        printf("Number of MPI ranks used: %i \n", nprocs);
        printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
    }
    if ( nprocs < nprocx*nprocy*nprocz )
        ERROR("Insufficient number of processors");

    // Get the rank info
	Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
    const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
	TwoPhase Averages(Dm);
	int N = (nx+2)*(ny+2)*(nz+2);

    // Read the local file
    DoubleArray Phase;
    DoubleArray SignDist;
    readRankData( rank, nx+2, ny+2, nz+2, Phase, SignDist );

    // Communication the halos
    fillHalo<double> fillData(rank_info,nx,ny,nz,1,1,1,0,1);
    fillData.fill(Phase);
    fillData.fill(SignDist);

    // Find blob domains
    if ( rank==0 ) { printf("Finding blob domains\n"); }
    double vF=0.0;
    double vS=0.0;
    IntArray GlobalBlobID;
    int nblobs = ComputeGlobalBlobIDs(nx,ny,nz,rank_info,
        Phase,SignDist,vF,vS,GlobalBlobID);
    if ( rank==0 ) { printf("Identified %i blobs\n",nblobs); }

    // Write the local blob ids
    char LocalRankFilename[100];
    sprintf(LocalRankFilename,"BlobLabel.%05i",rank);
    FILE *BLOBLOCAL = fopen(LocalRankFilename,"wb");
    fwrite(GlobalBlobID.get(),4,GlobalBlobID.length(),BLOBLOCAL);
    fclose(BLOBLOCAL);
    printf("Wrote BlobLabel.%05i \n",rank);

    //.......................................................................
	//copies of data needed to perform checkpointing from cpu
	double *Den, *DistEven, *DistOdd;
	Den = new double[2*N];
	DistEven = new double[10*N];
	DistOdd = new double[9*N];
	//.........................................................................
	if (rank==0) printf("Reading restart file! \n");
	// Read in the restart file to CPU buffers
	ReadCheckpoint(LocalRestartFile, Den, DistEven, DistOdd, N);
	MPI_Barrier(MPI_COMM_WORLD);
	//.........................................................................
	// Populate the arrays needed to perform averaging
	if (rank==0) printf("Populate arrays \n");
    // Compute porosity
	double porosity,sum,sum_global;
    sum=0.0;
    for (int n=0; n<(nx+2)*(ny+2)*(nz+2); n++){
        double phi,da,db,press,vx,vy,vz;
        double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
        da = Den[n];
        db = Den[N+n];
        f0 = DistEven[n];
        f2 = DistEven[N+n];
        f4 = DistEven[2*N+n];
        f6 = DistEven[3*N+n];
        f8 = DistEven[4*N+n];
        f10 = DistEven[5*N+n];
        f12 = DistEven[6*N+n];
        f14 = DistEven[7*N+n];
        f16 = DistEven[8*N+n];
        f18 = DistEven[9*N+n];
        //........................................................................
        f1 = DistOdd[n];
        f3 = DistOdd[1*N+n];
        f5 = DistOdd[2*N+n];
        f7 = DistOdd[3*N+n];
        f9 = DistOdd[4*N+n];
        f11 = DistOdd[5*N+n];
        f13 = DistOdd[6*N+n];
        f15 = DistOdd[7*N+n];
        f17 = DistOdd[8*N+n];
        //.................Compute the velocity...................................
        press = 0.3333333333333333*(f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+
                                          f9+f12+f11+f14+f13+f16+f15+f18+f17);
        vx = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
        vy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
        vz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
        Averages.Phase(n)=(da-db)/(da+db);
        Averages.Phase_tplus(n)=(da-db)/(da+db);
        Averages.Phase_tminus(n)=(da-db)/(da+db);
        Averages.Press(n)=press;
        Averages.Vel_x(n)=vx;
        Averages.Vel_y(n)=vy;
        Averages.Vel_z(n)=vz;
	if (Averages.SDs(n) > 0.0){
	  Dm.id[n]=1;
	  sum += 1.0;
	}
	else Dm.id[n]=0;
    }
    delete [] DistEven;
    delete [] DistOdd;

    MPI_Allreduce(&sum,&sum_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    porosity = sum_global/Dm.Volume;
    if (rank==0) printf("Porosity = %f \n",porosity);
    Dm.CommInit(MPI_COMM_WORLD);
    for (int i=0; i<N; i++) Averages.SDs(i) -= 1.0; // map the distance

    double beta = 0.95;

    Averages.SetupCubes(Dm);
    Averages.UpdateSolid();
    Averages.Initialize();
    Averages.ComputeDelPhi();
    Averages.ColorToSignedDistance(beta,Averages.Phase.get(),Averages.SDn.get());
    Averages.UpdateMeshValues();
    Averages.ComputeLocalBlob();
    Averages.Reduce();
    int b=0;

    //  Blobs.Set(Averages.BlobAverages.NBLOBS);
    int dimx = Averages.BlobAverages.size(0);
    int dimy = Averages.BlobAverages.size(1);
    int TotalBlobInfoSize=dimx*dimy;

    FILE *BLOBLOG;
    if (rank==0){
      	BLOBLOG=fopen("blobs.tcat","w");
        //printf("dimx=%i \n",dimx);
    }
    //      BlobContainer Blobs;
    DoubleArray RecvBuffer(dimx);
    //    MPI_Allreduce(&Averages.BlobAverages.get(),&Blobs.get(),1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("All ranks passed gate \n");

    for (int b=0; b<(int)Averages.BlobAverages.size(1); b++){

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

    Averages.SortBlobs();

      if (rank==0){
	//	printf("Reduced blob %i \n",b);
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
      }
    if (rank==0)  fclose(BLOBLOG);

    double Length=1.0;
    if (rank==0) WriteBlobStates(Averages,Length,porosity);

    /*FILE *BLOBS = fopen("Blobs.dat","wb");
    fwrite(GlobalBlobID.get(),4,Nx*Ny*Nz,BLOBS);
    fclose(BLOBS);*/

    PROFILE_STOP("main");
    PROFILE_SAVE("BlobIdentifyParallel",false);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;  
}

