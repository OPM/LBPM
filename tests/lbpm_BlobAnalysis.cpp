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

struct BlobTwoPhase{
    int COUNT; // number of averages to compute for each blob
    BlobTwoPhase(int size){
    	COUNT=26;
        NBLOBS=size;
        Data = new double [size*COUNT];
    }
    ~BlobTwoPhase(){
        delete [] Data;
    }
    int NBLOBS;
    double *Data;
    
    // if modified -- make sure to adjust COUNT so that
    // there is enough memory to save all the averages
    double Vn(int IDX){return Data[COUNT*IDX];}
    double pan(int IDX){return Data[COUNT*IDX+1];}
    double awn(int IDX){return Data[COUNT*IDX+2];}
    double ans(int IDX){return Data[COUNT*IDX+3];}
    double Jwn(int IDX){return Data[COUNT*IDX+4];}
    double Kwn(int IDX){return Data[COUNT*IDX+5];}
    double lwns(int IDX){return Data[COUNT*IDX+6];}
    double cwns(int IDX){return Data[COUNT*IDX+7];}
    double vanx(int IDX){return Data[COUNT*IDX+8];}
    double vany(int IDX){return Data[COUNT*IDX+9];}
    double vanz(int IDX){return Data[COUNT*IDX+10];}
    double vawnx(int IDX){return Data[COUNT*IDX+11];}
    double vawny(int IDX){return Data[COUNT*IDX+12];}
    double vawnz(int IDX){return Data[COUNT*IDX+13];}
    double Gwnxx(int IDX){return Data[COUNT*IDX+14];}
    double Gwnyy(int IDX){return Data[COUNT*IDX+15];}
    double Gwnzz(int IDX){return Data[COUNT*IDX+16];}
    double Gwnxy(int IDX){return Data[COUNT*IDX+17];}
    double Gwnxz(int IDX){return Data[COUNT*IDX+18];}
    double Gwnyz(int IDX){return Data[COUNT*IDX+19];}
    double Gnsxx(int IDX){return Data[COUNT*IDX+20];}
    double Gnsyy(int IDX){return Data[COUNT*IDX+22];}
    double Gnszz(int IDX){return Data[COUNT*IDX+23];}
    double Gnsxy(int IDX){return Data[COUNT*IDX+23];}
    double Gnsxz(int IDX){return Data[COUNT*IDX+24];}
    double Gnsyz(int IDX){return Data[COUNT*IDX+25];}

};

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

	int BC;	// type of boundary condition applied: 0-periodic, 1-pressure/velocity
	int nblobs_global; 	// number of blobs in the global system
    
    // Get the global number of blobs from arguments
    if (argc > 1){
        nblobs_global = atoi(argv[1]);
        if (rank==0)    printf("Number of global blobs is: %i \n",nblobs_global);
    }
    else{
        ERROR("Number of blobs was not specified");
    }

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
	BlobTwoPhase BlobAverages(nblobs_global);
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
    }
    int label;
	for (int k=1;k<Nz-1;k++){
		for (int j=1;j<Ny-1;j++){
			for (int i=1;i<Nx-1;i++){
			  // Assign the label for the cube
			  label = Averages.GetCubeLabel(i,j,k);
			  
			}
		}
	}

	/*	Averages.Initialize();
	Averages.ComputeDelPhi();
	Averages.ColorToSignedDistance(beta,Averages.Phase.data,Averages.SDn.data);
	Averages.UpdateMeshValues();
	Averages.ComputeLocal();
	Averages.Reduce();
	Averages.PrintAll(timestep);
*/
	// ****************************************************
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	// ****************************************************
}

