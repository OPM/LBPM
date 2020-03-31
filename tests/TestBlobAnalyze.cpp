// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "analysis/TwoPhase.h"
#include "common/SpherePack.h"

//#include "Domain.h"

using namespace std;

void readRankData( int proc, int nx, int ny, int nz, DoubleArray& Phase, DoubleArray& SignDist )
{
    Phase.resize(nx,ny,nz);
    SignDist.resize(nx,ny,nz);
    char file1[40], file2[40];
    sprintf(file1,"SignDist.%05d",proc);
    sprintf(file2,"Phase.%05d",proc);
    ReadBinaryFile(file1, Phase.data(), nx*ny*nz);
    ReadBinaryFile(file2, SignDist.data(), nx*ny*nz);
}
inline void WriteBlobs(TwoPhase Averages){
	printf("Writing the blob list \n");
    FILE *BLOBLOG;
  	BLOBLOG=fopen("blobs.tcat","w");
	fprintf(BLOBLOG,"%.5g %.5g %.5g\n",Averages.vol_w_global,Averages.paw_global,Averages.aws_global);
	for (int b=0; b<(int)Averages.ComponentAverages_NWP.size(1); b++){
		if (Averages.ComponentAverages_NWP(0,b) > 0.0){
			double Vn,pn,awn,ans,Jwn,Kwn,lwns,cwns;
			Vn = Averages.ComponentAverages_NWP(1,b);
			pn = Averages.ComponentAverages_NWP(2,b);
			awn = Averages.ComponentAverages_NWP(3,b);
			ans = Averages.ComponentAverages_NWP(4,b);
			Jwn = Averages.ComponentAverages_NWP(5,b);
			Kwn = Averages.ComponentAverages_NWP(6,b);
			lwns = Averages.ComponentAverages_NWP(7,b);
			cwns = Averages.ComponentAverages_NWP(8,b);

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

inline void  WriteBlobStates(TwoPhase TCAT, double D, double porosity){
	int a;
	double iVol=1.0/TCAT.Dm->Volume;
	double PoreVolume;
	double nwp_volume,vol_n,pan,pn,pw,pawn,pwn,awn,ans,aws,Jwn,Kwn,lwns,cwns,clwns;
	double sw,awnD,awsD,ansD,lwnsDD,JwnD,pc;
	nwp_volume=vol_n=pan=awn=ans=Jwn=Kwn=lwns=clwns=pawn=0.0;
	sw = TCAT.sat_w;
	pw = TCAT.paw_global;
	aws = TCAT.aws;
	// Compute the averages over the entire non-wetting phase
	printf("Writing blobstates.tcat for %i components \n",TCAT.NumberComponents_NWP);
	FILE *BLOBSTATES;
	BLOBSTATES = fopen("./blobstates.tcat","w");
	if (BLOBSTATES==NULL) ERROR("Cannot open blobstates.tcat for writing");
	for (a=0; a<TCAT.NumberComponents_NWP; a++){
		vol_n += TCAT.ComponentAverages_NWP(0,a);
		pan += TCAT.ComponentAverages_NWP(2,a)*TCAT.ComponentAverages_NWP(0,a);
		awn += TCAT.ComponentAverages_NWP(3,a);
		ans += TCAT.ComponentAverages_NWP(4,a);
		Jwn += TCAT.ComponentAverages_NWP(5,a)*TCAT.ComponentAverages_NWP(3,a);
		Kwn += TCAT.ComponentAverages_NWP(6,a)*TCAT.ComponentAverages_NWP(3,a);
		lwns += TCAT.ComponentAverages_NWP(7,a);
		clwns += TCAT.ComponentAverages_NWP(8,a)*TCAT.ComponentAverages_NWP(7,a);
		nwp_volume += TCAT.ComponentAverages_NWP(1,a);
		pawn += TCAT.ComponentAverages_NWP(2,a)*TCAT.ComponentAverages_NWP(3,a);
	}

	// Compute the pore voume (sum of wetting an non-wetting phase volumes)
	PoreVolume=TCAT.wp_volume_global + nwp_volume;
	// Subtract off portions of non-wetting phase in order of size
	for (a=TCAT.NumberComponents_NWP-1; a>0; a--){
		// Subtract the features one-by-one
		vol_n -= TCAT.ComponentAverages_NWP(0,a);
		pan -= TCAT.ComponentAverages_NWP(2,a)*TCAT.ComponentAverages_NWP(0,a);
		awn -= TCAT.ComponentAverages_NWP(3,a);
		ans -= TCAT.ComponentAverages_NWP(4,a);
		Jwn -= TCAT.ComponentAverages_NWP(5,a)*TCAT.ComponentAverages_NWP(3,a);
		Kwn -= TCAT.ComponentAverages_NWP(6,a)*TCAT.ComponentAverages_NWP(3,a);
		lwns -= TCAT.ComponentAverages_NWP(7,a);
		clwns -= TCAT.ComponentAverages_NWP(8,a)*TCAT.ComponentAverages_NWP(7,a);
		nwp_volume -= TCAT.ComponentAverages_NWP(1,a);
		pawn -= TCAT.ComponentAverages_NWP(2,a)*TCAT.ComponentAverages_NWP(3,a);

		// Update wetting phase averages
		aws += TCAT.ComponentAverages_NWP(4,a);
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
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&nprocs);
  { // Limit scope so variables that contain communicators will free before MPI_Finialize

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
    	if (domain.good()){
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
    	else if (nprocs==1){
    		nprocx=nprocy=nprocz=1;
    		nx=ny=nz=50;
    		nspheres=0;
    		Lx=Ly=Lz=1;
    	}
    	else if (nprocs==2){
    		nprocx=nprocy=1;
    		nprocz=2;
    		nx=ny=nz=50;
    		nspheres=0;
    		Lx=Ly=Lz=1;
    	}
    	else if (nprocs==4){
    		nprocx=nprocy=2;
    		nprocz=1;
    		nx=ny=nz=50;
    		nspheres=0;
    		Lx=Ly=Lz=1;
    	}
    	else if (nprocs==8){
    		nprocx=nprocy=nprocz=2;
    		nx=ny=nz=50;
    		nspheres=0;
    		Lx=Ly=Lz=1;
    	}
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
    if ( nprocs < nprocx*nprocy*nprocz )
        ERROR("Insufficient number of processors");

    // Set up the domain
	int BC=0;
    // Get the rank info
	std::shared_ptr<Domain> Dm(new Domain(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC));
 //   const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
	TwoPhase Averages(Dm);

	Nx = nx+2;
	Ny = ny+2;
	Nz = nz+2;
	if (rank == 0) cout << "Domain set." << endl;
	//.......................................................................
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				Dm->id[n] = 1;
			}
		}
	}
	//.......................................................................
    Dm->CommInit(); // Initialize communications for domains
	//.......................................................................
	// Read in sphere pack (initialize the non-wetting phase as inside of spheres)
        //
    nspheres=4;
	if (rank==1) printf("nspheres =%i \n",nspheres);
	//.......................................................................
	double *cx = new double[nspheres];
	double *cy = new double[nspheres];
	double *cz = new double[nspheres];
	double *rad = new double[nspheres];
	//.......................................................................
	//if (rank == 0)	printf("Reading the sphere packing \n");
	//if (rank == 0)	ReadSpherePacking(nspheres,cx,cy,cz,rad);
	// Hard coding the list of four spheres
	cx[0]=0.25*Lx; cx[1]=0.5*Lx; cx[2]=0.5*Lx; cx[3]=0.75*Lx;
	cy[0]=0.5*Ly; cx[1]=0.25*Ly; cx[2]=0.75*Ly; cx[3]=0.5*Ly;
	cz[0]=0.25*Lz; cx[1]=0.75*Lz; cx[2]=0.25*Lz; cx[3]=0.25*Lz;
	rad[0]=rad[1]=rad[2]=rad[3]=0.1*Lx;

	MPI_Barrier(comm);
	// Broadcast the sphere packing to all processes
	MPI_Bcast(cx,nspheres,MPI_DOUBLE,0,comm);
	MPI_Bcast(cy,nspheres,MPI_DOUBLE,0,comm);
	MPI_Bcast(cz,nspheres,MPI_DOUBLE,0,comm);
	MPI_Bcast(rad,nspheres,MPI_DOUBLE,0,comm);
	//...........................................................................
	MPI_Barrier(comm);
	//.......................................................................
	SignedDistance(Averages.Phase.data(),nspheres,cx,cy,cz,rad,Lx,Ly,Lz,Nx,Ny,Nz,
		       Dm->iproc(),Dm->jproc(),Dm->kproc(),Dm->nprocx(),Dm->nprocy(),Dm->nprocz());
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
				Averages.Phase(i,j,k) += 2.0;
				if (Averages.Phase(i,j,k) > 0.0){
					Dm->id[n] = 2;
				}
				else{
					Dm->id[n] = 1;
				}
				Averages.SDn(i,j,k) = -Averages.Phase(i,j,k);
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
	if (rank==0) printf("initializing the system \n");

	Averages.UpdateSolid();
    Dm->CommunicateMeshHalo(Averages.Phase);
    Dm->CommunicateMeshHalo(Averages.SDn);

    Averages.Initialize();
    Averages.UpdateMeshValues();

	if (rank==0) printf("computing local averages  \n");
	Averages.AssignComponentLabels();
    Averages.ComponentAverages();
    Averages.PrintComponents(int(5));
	if (rank==0) printf("reducing averages  \n");
   // Averages.Reduce();

    // Free memory
	delete [] cx;
	delete [] cy;
	delete [] cz;
	delete [] rad;

  } // Limit scope so variables that contain communicators will free before MPI_Finialize
  MPI_Barrier(comm);
  MPI_Finalize();
  return 0;  
}

