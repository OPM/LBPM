// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "analysis/pmmc.h"
#include "analysis/analysis.h"

//#include "Domain.h"

using namespace std;

inline void ReadCheckpoint(char *FILENAME, double *cDen, double *cDistEven, double *cDistOdd, int N)
{
    int q,n;
    double value;
    ifstream File(FILENAME,ios::binary);
    for (n=0; n<N; n++){
        // Write the two density values
        File.read((char*) &value, sizeof(value));
        cDen[n] = value;
    //    if (n== 66276)    printf("Density a  = %f \n",value);
        File.read((char*) &value, sizeof(value));
        cDen[N+n] = value;
    //    if (n== 66276)    printf("Density b  = %f \n",value);
        // Read the even distributions
        for (q=0; q<10; q++){
            File.read((char*) &value, sizeof(value));
            cDistEven[q*N+n] = value;
    //        if (n== 66276)    printf("dist even %i  = %f \n",q,value);
        }
        // Read the odd distributions
        for (q=0; q<9; q++){
            File.read((char*) &value, sizeof(value));
            cDistOdd[q*N+n] = value;
    //        if (n== 66276)    printf("dist even %i  = %f \n",q,value);
        }
    }
    File.close();
}

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

inline void SetPeriodicBC(DoubleArray &Scalar, int nx, int ny, int nz){
    
    int i,j,k,in,jn,kn;
    for (k=0; k<nz; k++){
        for (j=0; j<ny; j++){
            for (i=0; i<nx; i++){
                in = i; jn=j; kn=k;
                if (i==0) in = nx-2 ;
                else if (i==nx-1) in = 0;
                if (j==0) jn = ny-2;
                else if (j==ny-1) jn = 0;
                if (k==0) kn = nz-2;
                else if (k==nz-1) kn = 0;    
                Scalar(i,j,k) = Scalar(in,jn,kn);
            }
        }
    }
}
inline void ReadFromRank(char *FILENAME, DoubleArray &Phase, int nx, int ny, int nz, int iproc, int
                            jproc, int kproc)
{
    int i,j,k,q,n,N;
    int iglobal,jglobal,kglobal;
    double value;
    double denA,denB;

    N = nx*ny*nz;
    
    double *Den;
    
    Den = new double[2*N];

    ifstream File(FILENAME,ios::binary);
    for (n=0; n<N; n++){
        // Write the two density values
        File.read((char*) &value, sizeof(value));
        Den[2*n] = value;
        //    if (n== 66276)    printf("Density a  = %f \n",value);
        File.read((char*) &value, sizeof(value));
        Den[2*n+1] = value;

        //    if (n== 66276)    printf("Density b  = %f \n",value);
        // Read the even distributions
        for (q=0; q<10; q++){
            File.read((char*) &value, sizeof(value));
        }
        // Read the odd distributions
        for (q=0; q<9; q++){
            File.read((char*) &value, sizeof(value));
        }
    }
    File.close();
    
    // Compute the phase field
    for (k=1; k<nz-1; k++){
        for (j=1; j<ny-1; j++){
            for (i=1; i<nz-1; i++){
                //........................................................................
                n = k*nx*ny+j*nx+i;
                //........................................................................
                denA = Den[n];
                denB = Den[N+n];
                //........................................................................
                // save values in global arrays
                //........................................................................
                iglobal = iproc*(nx-2)+i;
                jglobal = jproc*(ny-2)+j;
                kglobal = kproc*(nz-2)+k;
                //........................................................................                
                Phase(iglobal,jglobal,kglobal) = (denA-denB)/(denA+denB);
                //........................................................................
            }
        }
    }
    
    delete Den;
}


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



int main(int argc, char **argv)
{
	// Initialize MPI
	MPI_Init(&argc,&argv);

    printf("-----------------------------------------------------------\n");
    printf("Labeling Blobs from Two-Phase Lattice Boltzmann Simulation \n");
    printf("-----------------------------------------------------------\n");

    //.......................................................................
    int nprocx,nprocy,nprocz,nprocs;
    int Nx, Ny, Nz;
    int nx,ny,nz;
    int nspheres;
    double Lx,Ly,Lz;

    //.......................................................................
    // Reading the domain information file
    //.......................................................................
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
    //.......................................................................

    nx+=2;
    ny+=2;
    nz+=2;
    
    nprocs = nprocx*nprocy*nprocz;
    printf("Number of MPI ranks: %i \n", nprocs);
    Nx = (nx-2)*nprocx+2;
    Ny = (ny-2)*nprocy+2;
    Nz = (nz-2)*nprocz+2;
    printf("Full domain size: %i x %i x %i  \n", Nx,Ny,Nz);

    DoubleArray Phase(Nx,Ny,Nz);
    DoubleArray SignDist(Nx,Ny,Nz);
    Phase.fill(0);
    SignDist.fill(-100.0);
    
    // read the files and populate main arrays
    for (int kproc=0; kproc<nprocz; kproc++){
        for (int jproc=0; jproc<nprocy; jproc++){
            for (int iproc=0; iproc<nprocx; iproc++){

                int proc = kproc*nprocx*nprocy + jproc*nprocx + iproc;
                DoubleArray PhaseTmp;
                DoubleArray SignDistTmp;
                readRankData( proc, nx, ny, nz, PhaseTmp, SignDistTmp );

                for (int k=1; k<nz-1; k++){
                    for (int j=1; j<ny-1; j++){
                        for (int i=1; i<nz-1; i++){
                            int iglobal = iproc*(nx-2)+i;
                            int jglobal = jproc*(ny-2)+j;
                            int kglobal = kproc*(nz-2)+k;
                            SignDist(iglobal,jglobal,kglobal) = SignDistTmp(i,j,k);
                        }
                    }
                }
                
                for (int k=1; k<nz-1; k++){
                    for (int j=1; j<ny-1; j++){
                        for (int i=1; i<nx-1; i++){
                            int iglobal = iproc*(nx-2)+i;
                            int jglobal = jproc*(ny-2)+j;
                            int kglobal = kproc*(nz-2)+k;
                            Phase(iglobal,jglobal,kglobal) = PhaseTmp(i,j,k);
                        }
                    }
                }
                
            }
        }
    }
    printf("Read %i ranks of Phase, SignDist \n",nprocs);

    SetPeriodicBC(SignDist, Nx, Ny, Nz);
    SetPeriodicBC(Phase, Nx, Ny, Nz);
    
    //FILE *PHASE;
    //PHASE = fopen("Phase.dat","wb");
    //fwrite(Phase.data,8,Nx*Ny*Nz,PHASE);
    //fclose(PHASE);
    
    
    // Compute the porosity
    double porosity=0.0;
    for (int k=0; k<Nz; k++){
        for (int j=0; j<Ny; j++){
            for (int i=0; i<Nx; i++){
                if (SignDist(i,j,k) > 0.0){ 
                    porosity += 1.0;
                }
            }
        }
    }
    //int N=int(porosity*1.25);
    //porosity /= (Nx*Ny*Nz*1.0);
    //printf("Media porosity is %f \n",porosity);


    /* ****************************************************************
                IDENTIFY ALL BLOBS: F > vF, S > vS
    ****************************************************************** */
    // Find blob domains, number of blobs
    double vF=0.0;
    double vS=0.0;
    printf("Execute blob identification algorithm... \n");
    IntArray GlobalBlobID;
    int nblobs = ComputeLocalBlobIDs( Phase, SignDist, vF, vS, GlobalBlobID );
    ReorderBlobIDs(GlobalBlobID,MPI_COMM_WORLD);       // This will reorder by blob size
    printf("Identified %i blobs. Writing per-process output files. \n",nblobs);

    int sizeLoc = nx*ny*nz;
    int *LocalBlobID;
    LocalBlobID = new int [sizeLoc];

    printf("File size (4 bytes per entry) %i, \n",sizeLoc);
    // read the files and populate main arrays
    for (int kproc=0; kproc<nprocz; kproc++){
        for (int jproc=0; jproc<nprocy; jproc++){
            for (int iproc=0; iproc<nprocx; iproc++){

                int proc = kproc*nprocx*nprocy + jproc*nprocx + iproc;
                char LocalRankFilename[40];
                sprintf(LocalRankFilename,"BlobLabel.%05i",proc);

                for (int k=0; k<nz; k++){
                    for (int j=0; j<ny; j++){
                        for (int i=0; i<nx; i++){
                            //........................................................................
                            int n = k*nx*ny+j*nx+i;
                            //........................................................................
                            int iglobal = iproc*(nx-2)+i;
                            int jglobal = jproc*(ny-2)+j;
                            int kglobal = kproc*(nz-2)+k;
                            // periodic BC
                            if (iglobal < 0 ) iglobal+=Nx;
                            if (jglobal < 0 ) jglobal+=Ny;
                            if (kglobal < 0 ) kglobal+=Nz;
                            if (!(iglobal < Nx) ) iglobal-=Nx;
                            if (!(jglobal < Ny) ) jglobal-=Ny;
                            if (!(kglobal < Nz) ) kglobal-=Nz;
                            //........................................................................
                            LocalBlobID[n] = GlobalBlobID(iglobal,jglobal,kglobal);
                            //........................................................................
                        }
                    }
                }

                FILE *BLOBLOCAL = fopen(LocalRankFilename,"wb");
                fwrite(&LocalBlobID[0],4,sizeLoc,BLOBLOCAL);
                fclose(BLOBLOCAL);
            }
        }
    }
    printf("Wrote %i ranks of BlobLabel.xxxxx \n",nprocs);


    FILE *BLOBS = fopen("Blobs.dat","wb");
    fwrite(GlobalBlobID.data(),4,Nx*Ny*Nz,BLOBS);
    fclose(BLOBS);
	MPI_Finalize();
    return 0;
}

