/*
This code computes TCAT averages on a blob-by-blob basis in parallel
It requires that the blobs be labeled using BlobIdentify.cpp 
James E. McClure 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "Domain.h"
#include "TwoPhase.h"
#include "common/MPI_Helpers.h"
#include "common/Utilities.h"

inline void ReadBlobFile(char *FILENAME, int *Data, int N)
{
    int n;
    int value;
    ifstream File(FILENAME,ios::binary);
    if (File.good()){
        for (n=0; n<N; n++){
            // Write the two density values
            File.read((char*) &value, sizeof(value));
            Data[n] = value;
            
        }
    }
    else {
        for (n=0; n<N; n++) Data[n] = 1.2e-34;
    }
    File.close();
    
}

inline void  WriteBlobStates(TwoPhase TCAT, double D, double porosity);

int main(int argc, char **argv)
{
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	// parallel domain size (# of sub-domains)
	int nprocx,nprocy,nprocz;
	int iproc,jproc,kproc;
	int Nx,Ny,Nz,N,nspheres;
	double Lx,Ly,Lz;

	int BC=0;	// type of boundary condition applied: 0-periodic, 1-pressure/velocity
	int nblobs_global=0; 	// number of blobs in the global system
    
    // Get the global number of blobs from arguments
	/*    if (argc > 1){
        nblobs_global = atoi(argv[1]);
        if (rank==0)    printf("Number of global blobs is: %i \n",nblobs_global);
    }
    else{
        ERROR("Number of blobs was not specified");
    }
	*/
	
	int *CubeList;

	if (rank==0){
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		ifstream domain("Domain.in");
		domain >> nprocx;
		domain >> nprocy;
		domain >> nprocz;
		domain >> Nx;
		domain >> Ny;
		domain >> Nz;
		domain >> nspheres;
		domain >> Lx;
		domain >> Ly;
		domain >> Lz;
		//.......................................................................
	}
	//.................................................
	MPI_Barrier(MPI_COMM_WORLD);
	// Computational domain
	MPI_Bcast(&Nx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Ny,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Nz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocy,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nspheres,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	//.................................................
	Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
	TwoPhase Averages(Dm);
	//	BlobTwoPhase BlobAverages(nblobs_global);
	Nx+=2;Ny+=2;Nz+=2;
	N=Nx*Ny*Nz; // number of lattice points
	//.......................................................................
	// Filenames used
	char LocalRankString[8];
	char LocalRankFilename[40];
	char LocalRestartFile[40];
	char tmpstr[10];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
	//...........................................................................
	if (rank == 0) cout << "Reading in domain from signed distance function..." << endl;
	//.......................................................................
	sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
	ReadBinaryFile(LocalRankFilename, Averages.SDs.data, N);
	MPI_Barrier(MPI_COMM_WORLD);

	//	sprintf(LocalRankFilename,"%s%s","Pressure.",LocalRankString);
	//ReadBinaryFile(LocalRankFilename, Averages.Press.data, N);
	//MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) cout << "Domain set." << endl;
    //.......................................................................
    sprintf(LocalRankFilename,"%s%s","BlobLabel.",LocalRankString);
    ReadBlobFile(LocalRankFilename, Averages.BlobLabel.data, N);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) cout << "BlobLabel set." << endl;

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
    for (int n=0; n<Nx*Ny*Nz; n++){
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
        Averages.Phase.data[n]=(da-db)/(da+db);
        Averages.Phase_tplus.data[n]=(da-db)/(da+db);
        Averages.Phase_tminus.data[n]=(da-db)/(da+db);
        Averages.Press.data[n]=press;
        Averages.Vel_x.data[n]=vx;
        Averages.Vel_y.data[n]=vy;
        Averages.Vel_z.data[n]=vz;
	if (Averages.SDs.data[n] > 0.0){
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
    for (int i=0; i<N; i++) Averages.SDs.data[i] -= 1.0; // map the distance 
    
    double beta = 0.95;
     
    Averages.SetupCubes(Dm);
    Averages.UpdateSolid();
    Averages.Initialize();
    Averages.ComputeDelPhi();
    Averages.ColorToSignedDistance(beta,Averages.Phase.data,Averages.SDn.data);
    Averages.UpdateMeshValues();
    Averages.ComputeLocalBlob();
    Averages.Reduce();
    int b=0;

    //  Blobs.Set(Averages.BlobAverages.NBLOBS);
    int dimx = Averages.BlobAverages.m;
    int dimy = Averages.BlobAverages.n;
    int TotalBlobInfoSize=Averages.BlobAverages.m*Averages.BlobAverages.n;
   
    FILE *BLOBLOG;
    if (rank==0){
      	BLOBLOG=fopen("blobs.tcat","w");	
        //printf("dimx=%i \n",dimx);
    }
    //      BlobContainer Blobs;
    DoubleArray RecvBuffer(dimx);
    //    MPI_Allreduce(&Averages.BlobAverages.data,&Blobs.data,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("All ranks passed gate \n");

    for (int b=0; b<Averages.BlobAverages.n; b++){
      
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
    for (int b=0; b<Averages.BlobAverages.n; b++){
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

    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("Exit, rank=%i \n",rank);
	// ****************************************************
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	// ****************************************************
}

inline void  WriteBlobStates(TwoPhase TCAT, double D, double porosity){
	FILE *BLOBSTATES= fopen("blobstates.tcat","w");
	int a;
	double iVol=1.0/TCAT.Dm.Volume;
	double PoreVolume;
	double nwp_volume,vol_n,pan,pn,pw,pawn,pwn,awn,ans,aws,Jwn,Kwn,lwns,cwns,clwns;
	double sw,awnD,awsD,ansD,lwnsDD,JwnD,pc;
	nwp_volume=vol_n=pan=awn=ans=Jwn=Kwn=lwns=clwns=pawn=0.0;
	pw = TCAT.paw_global / TCAT.vol_w_global;
	aws = TCAT.aws;
	// Compute the averages over the entire non-wetting phsae
	for (a=0; a<TCAT.nblobs_global; a++){
		vol_n += TCAT.BlobAverages(0,a);
		pan += TCAT.BlobAverages(2,a)*TCAT.BlobAverages(1,a);
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
		pan -= TCAT.BlobAverages(2,a)*TCAT.BlobAverages(1,a);
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
	fclose(BLOBSTATES);
}
