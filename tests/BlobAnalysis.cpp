// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "analysis/pmmc.h"
#include "analysis/analysis.h"
//#include "Domain.h"

#define NUM_AVERAGES 30

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
inline void ReadFromRank(char *FILENAME, DoubleArray &Phase, DoubleArray &Pressure, DoubleArray &Vel_x, 
                            DoubleArray &Vel_y, DoubleArray &Vel_z, int nx, int ny, int nz, int iproc, int
                            jproc, int kproc)
{
    int i,j,k,q,n,N;
    int iglobal,jglobal,kglobal;
    double value;
    double denA,denB;
    double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
    double f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double vx,vy,vz;
    
    N = nx*ny*nz;
    
    double *Den, *DistEven, *DistOdd;
    
    Den = new double[2*N];
    DistEven = new double[10*N];
    DistOdd = new double[9*N];

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
            DistEven[q*N+n] = value;
        }
        // Read the odd distributions
        for (q=0; q<9; q++){
            File.read((char*) &value, sizeof(value));
            DistOdd[q*N+n] = value;
        }
    }
    File.close();
    
    // Compute the phase field, pressure and velocity
    for (k=1; k<nz-1; k++){
        for (j=1; j<ny-1; j++){
            for (i=1; i<nz-1; i++){
                //........................................................................
                n = k*nx*ny+j*nx+i;
                //........................................................................
                denA = Den[n];
                denB = Den[N+n];
                //........................................................................
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
                //........................................................................
                //.................Compute the pressure....................................
                value = 0.3333333333333333*(f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17);
                //........................................................................
                //.................Compute the velocity...................................
                vx = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
                vy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
                vz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
                //........................................................................
                // save values in global arrays
                //........................................................................
                iglobal = iproc*(nx-2)+i;
                jglobal = jproc*(ny-2)+j;
                kglobal = kproc*(nz-2)+k;
                //........................................................................                
                Phase(iglobal,jglobal,kglobal) = (denA-denB)/(denA+denB);
                Pressure(iglobal,jglobal,kglobal) = value;
                Vel_x(iglobal,jglobal,kglobal) = vx;
                Vel_y(iglobal,jglobal,kglobal) = vy;
                Vel_z(iglobal,jglobal,kglobal) = vz;
                //........................................................................
            }
        }
    }
    
    delete Den;
    delete DistEven;
    delete DistOdd;
}

int main(int argc, char **argv)
{
    printf("----------------------------------------------------------\n");
    printf("COMPUTING TCAT ANALYSIS FOR NON-WETTING PHASE FEATURES \n");
    printf("----------------------------------------------------------\n");

    //.......................................................................
    int nprocx,nprocy,nprocz,nprocs;
    int Nx, Ny, Nz;
    int nx,ny,nz;
    int nspheres;
    double Lx,Ly,Lz;
    //.......................................................................
    int i,j,k,n,p,idx;
    int iproc,jproc,kproc;
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
    DoubleArray Press(Nx,Ny,Nz);
    DoubleArray Vel_x(Nx,Ny,Nz);            // Velocity
    DoubleArray Vel_y(Nx,Ny,Nz);
    DoubleArray Vel_z(Nx,Ny,Nz);
    DoubleArray dPdt(Nx,Ny,Nz);
    
    // Filenames used
    char LocalRankString[8];
    char LocalRankFilename[40];
    char BaseFilename[20];
    sprintf(BaseFilename,"%s","dPdt.");
                     
    int proc,iglobal,kglobal,jglobal;

    double * Temp;
    Temp = new double[nx*ny*nz];
    
    for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
            for (i=0; i<Nx; i++){
                SignDist(i,j,k) = -100.0;
            }
        }
    }
    
#ifdef GENTEST
    // Fill the arrays with test data
    double minValue;
    for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
            for (i=0; i<Nx; i++){
                SignDist(i,j,k) = 100.0;
                
                minValue = sqrt((1.0*i-0.3*Nx)*(1.0*i-0.3*Nx)+(1.0*j-0.3*Ny)*(1.0*j-0.3*Ny)
                        +(1.0*k-0.3*Nz)*(1.0*k-0.3*Nz))-0.1*Nx;
                
                if (sqrt((1.0*i-0.7*Nx)*(1.0*i-0.7*Nx)+(1.0*j-0.7*Ny)*(1.0*j-0.7*Ny)
                        +(1.0*k-0.7*Nz)*(1.0*k-0.7*Nz))-0.2*Nx < minValue){
                    minValue = sqrt((1.0*i-0.7*Nx)*(1.0*i-0.7*Nx)+(1.0*j-0.7*Ny)*(1.0*j-0.7*Ny)
                            +(1.0*k-0.7*Nz)*(1.0*k-0.7*Nz))-0.2*Nx;
                }
                
                if (sqrt((1.0*i-0.2*Nx)*(1.0*i-0.2*Nx)+(1.0*j-0.7*Ny)*(1.0*j-0.7*Ny)
                        +(1.0*k-0.6*Nz)*(1.0*k-0.6*Nz))-0.13*Nx < minValue){
                    minValue = sqrt((1.0*i-0.2*Nx)*(1.0*i-0.2*Nx)+(1.0*j-0.7*Ny)*(1.0*j-0.7*Ny)
                            +(1.0*k-0.6*Nz)*(1.0*k-0.6*Nz))-0.13*Nx;
                }
                
                if (sqrt((1.0*i-0.7*Nx)*(1.0*i-0.7*Nx)+(1.0*j-0.3*Ny)*(1.0*j-0.3*Ny)
                        +(1.0*k-0.7*Nz)*(1.0*k-0.7*Nz))-0.17*Nx < minValue){
                    minValue = sqrt((1.0*i-0.7*Nx)*(1.0*i-0.7*Nx)+(1.0*j-0.3*Ny)*(1.0*j-0.3*Ny)
                            +(1.0*k-0.7*Nz)*(1.0*k-0.7*Nz))-0.17*Nx;
                }
                
                if (minValue < -1.0)        Phase(i,j,k) = 1.0;
                else if (minValue < 1.0)    Phase(i,j,k) = -minValue;
                else                         Phase(i,j,k) = -1.0;
                
                if (Phase(i,j,k) > 0.0){
                    Press(i,j,k) = 0.34;
                }
                else {
                    Press(i,j,k) = 0.32;
                }
            }
        }
    }
#else
    // read the files and populate main arrays
    for ( kproc=0; kproc<nprocz; kproc++){
        for ( jproc=0; jproc<nprocy; jproc++){
            for ( iproc=0; iproc<nprocx; iproc++){
                
                proc = kproc*nprocx*nprocy + jproc*nprocx + iproc;

                sprintf(LocalRankString,"%05d",proc);
                sprintf(LocalRankFilename,"%s%s","dPdt.",LocalRankString);
    //            printf("Reading file %s \n",LocalRankFilename);
                ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);
                for (k=1; k<nz-1; k++){
                    for (j=1; j<ny-1; j++){
                        for (i=1; i<nz-1; i++){
                            //........................................................................
                            n = k*nx*ny+j*nx+i;
                            //........................................................................
                            iglobal = iproc*(nx-2)+i;
                            jglobal = jproc*(ny-2)+j;
                            kglobal = kproc*(nz-2)+k;
                            //........................................................................
                            dPdt(iglobal,jglobal,kglobal) = Temp[n]; 
                            //........................................................................
                        }
                    }
                }
                
                sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
        //        printf("Reading file %s \n",LocalRankFilename);
        //        printf("Sub-domain size: %i x %i x %i  \n", nx,ny,nz);
                ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);    
                for (k=1; k<nz-1; k++){
                    for (j=1; j<ny-1; j++){
                        for (i=1; i<nz-1; i++){

                            //........................................................................
                            n = k*nx*ny+j*nx+i;
                            //........................................................................
                            iglobal = iproc*(nx-2)+i;
                            jglobal = jproc*(ny-2)+j;
                            kglobal = kproc*(nz-2)+k;
                            //........................................................................
                            SignDist(iglobal,jglobal,kglobal) = Temp[n];
                            //........................................................................
                        }
                    }
                }
                
                sprintf(LocalRankFilename,"%s%s","Restart.",LocalRankString);

                ReadFromRank(LocalRankFilename,Phase,Press,Vel_x,Vel_y,Vel_z,nx,ny,nz,iproc,jproc,kproc);

                sprintf(LocalRankFilename,"%s%s","Pressure.",LocalRankString);
                
                ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);    
                for (k=1; k<nz-1; k++){
                    for (j=1; j<ny-1; j++){
                        for (i=1; i<nx-1; i++){

                            //........................................................................
                            n = k*nx*ny+j*nx+i;
                            //........................................................................
                            iglobal = iproc*(nx-2)+i;
                            jglobal = jproc*(ny-2)+j;
                            kglobal = kproc*(nz-2)+k;
                            //........................................................................
                            Press(iglobal,jglobal,kglobal) = Temp[n];
                            //........................................................................
                        }
                    }
                }

                sprintf(LocalRankFilename,"%s%s","Phase.",LocalRankString);
                ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);    
                for (k=1; k<nz-1; k++){
                    for (j=1; j<ny-1; j++){
                        for (i=1; i<nx-1; i++){

                            //........................................................................
                            n = k*nx*ny+j*nx+i;
                            //........................................................................
                            iglobal = iproc*(nx-2)+i;
                            jglobal = jproc*(ny-2)+j;
                            kglobal = kproc*(nz-2)+k;
                            //........................................................................
                            Phase(iglobal,jglobal,kglobal) = Temp[n];
                            //........................................................................
                        }
                    }
                }
                
                
            }
        }
    }
    printf("Read %i ranks of %s \n",nprocs,BaseFilename);
#endif
    
    delete Temp;
    
    DoubleArray MeanCurvature(Nx,Ny,Nz);
    DoubleArray GaussCurvature(Nx,Ny,Nz);
    DoubleArray SignDist_x(Nx,Ny,Nz);        // Gradient of the signed distance
    DoubleArray SignDist_y(Nx,Ny,Nz);
    DoubleArray SignDist_z(Nx,Ny,Nz);
    DoubleArray Phase_x(Nx,Ny,Nz);            // Gradient of the phase indicator field
    DoubleArray Phase_y(Nx,Ny,Nz);
    DoubleArray Phase_z(Nx,Ny,Nz);

    SetPeriodicBC(SignDist, Nx, Ny, Nz);
    SetPeriodicBC(Phase, Nx, Ny, Nz);
    SetPeriodicBC(Press, Nx, Ny, Nz);
    
    //...........................................................................
    // Compute the gradients of the phase indicator and signed distance fields
    //...........................................................................
    pmmc_MeshGradient(Phase,Phase_x,Phase_y,Phase_z,Nx,Ny,Nz);
    pmmc_MeshGradient(SignDist,SignDist_x,SignDist_y,SignDist_z,Nx,Ny,Nz);
    //...........................................................................
    // Compute the mesh curvature of the phase indicator field
    pmmc_MeshCurvature(Phase, MeanCurvature, GaussCurvature, Nx, Ny, Nz);
    //...........................................................................
    
    SetPeriodicBC(MeanCurvature, Nx, Ny, Nz);
    SetPeriodicBC(GaussCurvature, Nx, Ny, Nz);
    SetPeriodicBC(SignDist_x, Nx, Ny, Nz);
    SetPeriodicBC(SignDist_y, Nx, Ny, Nz);
    SetPeriodicBC(SignDist_z, Nx, Ny, Nz);
    SetPeriodicBC(Phase_x, Nx, Ny, Nz);
    SetPeriodicBC(Phase_y, Nx, Ny, Nz);
    SetPeriodicBC(Phase_z, Nx, Ny, Nz);
    
/*    for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
            for (i=0; i<Nx; i++){
                if (i==0)         SignDist(i,j,k) = -2.0;
                if (j==0)         SignDist(i,j,k) = -2.0; 
                if (k==0)         SignDist(i,j,k) = -2.0;
                if (i==Nx-1)     SignDist(i,j,k) = -2.0;
                if (j==Ny-1)     SignDist(i,j,k) = -2.0;
                if (k==Nz-1)     SignDist(i,j,k) = -2.0;    
            }
        }
    }
*/
    
    FILE *PHASE = fopen("Phase.dat","wb");
    fwrite(Phase.data(),8,Nx*Ny*Nz,PHASE);
    fclose(PHASE);
    
    // Compute the porosity
    double porosity=0.0;
    for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
            for (i=0; i<Nx; i++){
                if (SignDist(i,j,k) > 0.0){ 
                    porosity += 1.0;
                }
            }
        }
    }
    porosity /= (Nx*Ny*Nz*1.0);

    printf("Media porosity is %f \n",porosity);

    /* ****************************************************************
     VARIABLES FOR THE PMMC ALGORITHM
     ****************************************************************** */
    //...........................................................................
    // Averaging variables
    //...........................................................................
    double awn,ans,aws,lwns,nwp_volume;
    double sw,vol_n,vol_w,paw,pan;
    double efawns,Jwn,Kwn;
    double trawn,trJwn,trRwn;
    double As,dummy;
    //  double dEs,dAwn,dAns;    // Global surface energy (calculated by rank=0)
    //    bool add=1;            // Set to false if any corners contain nw-phase ( F > fluid_isovalue)
    
    int n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0;
    int n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
    
    //double s,s1,s2,s3;        // Triangle sides (lengths)
    Point A,B,C,P;
    //    double area;
    int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
    //    int count_in=0,count_out=0;
    //    int nodx,nody,nodz;
    // initialize lists for vertices for surfaces, common line
    DTMutableList<Point> nw_pts(20);
    DTMutableList<Point> ns_pts(20);
    DTMutableList<Point> ws_pts(20);
    DTMutableList<Point> nws_pts(20);
    // initialize triangle lists for surfaces
    IntArray nw_tris(3,20);
    IntArray ns_tris(3,20);
    IntArray ws_tris(3,20);
    // initialize list for line segments
    IntArray nws_seg(2,20);
    DTMutableList<Point> tmp(20);
    
    // Initialize arrays for local solid surface
    DTMutableList<Point> local_sol_pts(20);
    int n_local_sol_pts = 0;
    IntArray local_sol_tris(3,18);
    int n_local_sol_tris;
    DoubleArray values(20);
    DTMutableList<Point> local_nws_pts(20);
    int n_local_nws_pts;
    
    DoubleArray CubeValues(2,2,2);
    DoubleArray Values(20);
    DoubleArray ContactAngle(20);
    DoubleArray Curvature(20);
    DoubleArray DistValues(20);
    DoubleArray InterfaceSpeed(20);
    DoubleArray NormalVector(60);
    DoubleArray van(3);
    DoubleArray vaw(3);
    DoubleArray vawn(3);
    DoubleArray Gwn(6);
    DoubleArray Gns(6);
    DoubleArray Gws(6);
    //...........................................................................
    
    printf("Execute blob identification algorithm... \n");

    /* ****************************************************************
                IDENTIFY ALL BLOBS: F > vF, S > vS
    ****************************************************************** */
    double vF=0.0;
    double vS=0.0;
    double trimdist=1.0;
    // Compute the number of blobs and create the local blob ids
    IntArray LocalBlobID;
    int nblobs = ComputeLocalBlobIDs(Phase,SignDist,vF,vS,LocalBlobID,true); 
    // Compute the cells in each blob
    IntArray b(nblobs);        // number of nodes in each blob
    IntArray blobs(3,Nx*Ny*Nz); // The indicies for the cells in each blob
    b.fill(0);
    blobs.fill(0);
    int ncubes = 0;
    for (int bb=0; bb<nblobs; bb++) {
        for (k=0;k<Nz;k++){
            for (j=0;j<Ny;j++){
                for (i=0;i<Nx;i++){
                    if ( LocalBlobID(i,j,k)==bb ) {
                        b(bb)++;
                        blobs(0,ncubes) = i;
                        blobs(1,ncubes) = j;
                        blobs(2,ncubes) = k;
                        ncubes++;
                    }
                }
            }
        }
    }
    b.resize(ncubes);
    
    DoubleArray BlobAverages(NUM_AVERAGES,nblobs);
    
    // Map the signed distance for the analysis
    for (i=0; i<Nx*Ny*Nz; i++)    SignDist(i) -= (1.0); 
    
    // Compute the porosity
    porosity=0.0;
    for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
            for (i=0; i<Nx; i++){
                if (SignDist(i,j,k) > 0.0){ 
                    porosity += 1.0;
                }
            }
        }
    }
    porosity /= (Nx*Ny*Nz*1.0);

    printf("Media porosity is %f \n",porosity);

    /* ****************************************************************
            RUN TCAT AVERAGING ON EACH BLOB
    ****************************************************************** */
    int n_nw_tris_beg, n_ns_tris_beg, n_ws_tris_beg, n_nws_seg_beg;
    int start=0,finish;
    int a,c;
    
    printf("-----------------------------------------------\n");
    printf("Computing TCAT averages based on connectivity \n");
    printf("The number of non-wetting phase features is %i \n",nblobs-1);
    printf("-----------------------------------------------\n");

    // Wetting phase averages assume global connectivity
    As = 0.0;
    vol_w = 0.0;
    paw = 0.0;
    aws = 0.0;
    vaw(0) = vaw(1) = vaw(2) = 0.0;
    Gws(0) = Gws(1) = Gws(2) = 0.0;
    Gws(3) = Gws(4) = Gws(5) = 0.0;    

    // Don't compute the last blob unless specified
    // the last blob is the entire wetting phase
//    nblobs -=1;
#ifdef WP
    nblobs+=1;
#endif
    
    for (a=0;a<nblobs;a++){
        
        finish = start+b(a);

        /* ****************************************************************
         RUN PMMC ON EACH BLOB
         ****************************************************************** */

        // Store beginning points for surfaces for blob p
        n_nw_tris_beg = n_nw_tris;
        n_ns_tris_beg = n_ns_tris;
        n_ws_tris_beg = n_ws_tris;
        n_nws_seg_beg = n_nws_seg;
        // Loop over all cubes
        nwp_volume = 0.0;
        
        // Compute phase averages
        vol_n =0.0;
        pan = 0.0;
        awn  = ans = lwns = 0.0;
        van(0) = van(1) = van(2) = 0.0;
        vawn(0) = vawn(1) = vawn(2) = 0.0;
        Gwn(0) = Gwn(1) = Gwn(2) = 0.0;
        Gwn(3) = Gwn(4) = Gwn(5) = 0.0;
        Gns(0) = Gns(1) = Gns(2) = 0.0;
        Gns(3) = Gns(4) = Gns(5) = 0.0;
        Jwn = Kwn = efawns = 0.0;
        trJwn = trawn = trRwn = 0.0;
                
        for (c=start;c<finish;c++){
            // Get cube from the list
            i = blobs(0,c);
            j = blobs(1,c);
            k = blobs(2,c);
    
            // Use the cube to compute volume averages
            for (p=0;p<8;p++){
                if ( SignDist(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){
                    
                    n = i+cube[p][0] + Nx*(j+cube[p][1]) + Nx*Ny*(k+cube[p][2]);

                    // Compute the non-wetting phase volume contribution
                    if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 )
                        nwp_volume += 0.125;

                    // volume averages over the non-wetting phase
                    if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0.99 ){
                        // volume the excludes the interfacial region
                        vol_n += 0.125;
                        // pressure
                        pan += 0.125*Press(n);
                        // velocity
                        van(0) += 0.125*Vel_x(n);
                        van(1) += 0.125*Vel_y(n);
                        van(2) += 0.125*Vel_z(n);
                    }

                    // volume averages over the wetting phase
                    if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) < -0.99 ){
                        // volume the excludes the interfacial region
                        vol_w += 0.125;
                        // pressure
                        paw += 0.125*Press(n);
                        // velocity
                        vaw(0) += 0.125*Vel_x(n);
                        vaw(1) += 0.125*Vel_y(n);
                        vaw(2) += 0.125*Vel_z(n);
                    }
                }
            }

            // Interface and common curve averages
            n_local_sol_tris = 0;
            n_local_sol_pts = 0;
            n_local_nws_pts = 0;
    
            //...........................................................................
            // Construct the interfaces and common curve
            pmmc_ConstructLocalCube(SignDist, Phase, vS, vF,
                    nw_pts, nw_tris, values, ns_pts, ns_tris, ws_pts, ws_tris,
                    local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                    n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                    n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
                    i, j, k, Nx, Ny, Nz);

            // Integrate the contact angle
            efawns += pmmc_CubeContactAngle(CubeValues,Values,Phase_x,Phase_y,Phase_z,SignDist_x,SignDist_y,SignDist_z,
                                            local_nws_pts,i,j,k,n_local_nws_pts);

            // Integrate the mean curvature
            Jwn    += pmmc_CubeSurfaceInterpValue(CubeValues,MeanCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);
            Kwn    += pmmc_CubeSurfaceInterpValue(CubeValues,GaussCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);

            // Integrate the trimmed mean curvature (hard-coded to use a distance of 4 pixels)
            pmmc_CubeTrimSurfaceInterpValues(CubeValues,MeanCurvature,SignDist,nw_pts,nw_tris,Values,DistValues,
                        i,j,k,n_nw_pts,n_nw_tris,trimdist,trawn,trJwn);
            
            pmmc_CubeTrimSurfaceInterpInverseValues(CubeValues,MeanCurvature,SignDist,nw_pts,nw_tris,Values,DistValues,
                        i,j,k,n_nw_pts,n_nw_tris,trimdist,dummy,trRwn);
            
            // Compute the normal speed of the interface
            pmmc_InterfaceSpeed(dPdt, Phase_x, Phase_y, Phase_z, CubeValues, nw_pts, nw_tris,
                                NormalVector, InterfaceSpeed, vawn, i, j, k, n_nw_pts, n_nw_tris);
            
            As  += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);

            // Compute the surface orientation and the interfacial area
            awn += pmmc_CubeSurfaceOrientation(Gwn,nw_pts,nw_tris,n_nw_tris);
            ans += pmmc_CubeSurfaceOrientation(Gns,ns_pts,ns_tris,n_ns_tris);
            aws += pmmc_CubeSurfaceOrientation(Gws,ws_pts,ws_tris,n_ws_tris);
            lwns +=  pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
            //...........................................................................
    
            //*******************************************************************
            // Reset the triangle counts to zero
            n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0;
            n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
            
            n_nw_tris_beg = n_nw_tris;
            n_ns_tris_beg = n_ns_tris;
            n_ws_tris_beg = n_ws_tris;
            n_nws_seg_beg = n_nws_seg;
            //*******************************************************************
            
        }
        start = finish;

        if (a < nblobs-1){
        if (vol_n > 0.0){
            pan /= vol_n;
            for (i=0;i<3;i++)    van(i) /= vol_n;
        }
        if (awn > 0.0){
            Jwn /= awn;
            Kwn /= awn;
            for (i=0;i<3;i++)    vawn(i) /= awn;
            for (i=0;i<6;i++)    Gwn(i) /= awn;
        }
        if (lwns > 0.0){
            efawns /= lwns;
        }
        if (ans > 0.0){
            for (i=0;i<6;i++)    Gns(i) /= ans;
        }
        if (trawn > 0.0){
            trJwn /= trawn;
        }
        
        BlobAverages(0,a) = nwp_volume;
        BlobAverages(1,a) = pan;
        BlobAverages(2,a) = awn;
        BlobAverages(3,a) = ans;
        BlobAverages(4,a) = Jwn;
        BlobAverages(5,a) = Kwn;
        BlobAverages(6,a) = lwns;
        BlobAverages(7,a) = efawns;
        BlobAverages(8,a) = van(0);
        BlobAverages(9,a) = van(1);
        BlobAverages(10,a) = van(2);
        BlobAverages(11,a) = vawn(0);
        BlobAverages(12,a) = vawn(1);
        BlobAverages(13,a) = vawn(2);
        BlobAverages(14,a) = Gwn(0);
        BlobAverages(15,a) = Gwn(1);
        BlobAverages(16,a) = Gwn(2);
        BlobAverages(17,a) = Gwn(3);
        BlobAverages(18,a) = Gwn(4);
        BlobAverages(19,a) = Gwn(5);
        BlobAverages(20,a) = Gns(0);
        BlobAverages(21,a) = Gns(1);
        BlobAverages(22,a) = Gns(2);
        BlobAverages(23,a) = Gns(3);
        BlobAverages(24,a) = Gns(4);
        BlobAverages(25,a) = Gns(5);
        BlobAverages(26,a) = trawn;
        BlobAverages(27,a) = trJwn;
        BlobAverages(28,a) = vol_n;
        BlobAverages(29,a) = trRwn;

        printf("Computed TCAT averages for feature = %i \n", a);
        }
        
    }  // End of the blob loop
    NULL_USE(n_nw_tris_beg);
    NULL_USE(n_ns_tris_beg);
    NULL_USE(n_ws_tris_beg);
    NULL_USE(n_nws_seg_beg);
        
    nblobs -= 1;
    printf("-----------------------------------------------\n");
    printf("Sorting the blobs based on volume \n");
    printf("-----------------------------------------------\n");
    int TempLabel,aa,bb;
    double TempValue;
    IntArray OldLabel(nblobs);
    for (a=0; a<nblobs; a++)    OldLabel(a) = a;
    // Sort the blob averages based on volume
    for (aa=0; aa<nblobs-1; aa++){
        for ( bb=aa+1; bb<nblobs; bb++){
            if (BlobAverages(0,aa) < BlobAverages(0,bb)){
                // Exchange location of blobs aa and bb
                //printf("Switch blob %i with %i \n", OldLabel(aa),OldLabel(bb));
                // switch the label
                TempLabel = OldLabel(bb);
                OldLabel(bb) = OldLabel(aa);
                OldLabel(aa) = TempLabel;
                // switch the averages
                for (idx=0; idx<NUM_AVERAGES; idx++){
                    TempValue = BlobAverages(idx,bb);
                    BlobAverages(idx,bb) = BlobAverages(idx,aa);
                    BlobAverages(idx,aa) = TempValue;
                }
            }
        }        
    }
    
    IntArray NewLabel(nblobs);
    for (aa=0; aa<nblobs; aa++){
        // Match the new label for original blob aa
        bb=0;
        while (OldLabel(bb) != aa)    bb++;
        NewLabel(aa) = bb;
    }
    
    // Re-label the blob ID
    printf("Re-labeling the blobs, now indexed by volume \n");
    for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
            for (i=0; i<Nx; i++){
                if (LocalBlobID(i,j,k) > -1){
                    TempLabel = NewLabel(LocalBlobID(i,j,k));
                    LocalBlobID(i,j,k) = TempLabel;
                }
            }
        }
    }
        
    FILE *BLOBLOG= fopen("blobs.tcat","a");
    for (a=0; a<nblobs; a++){
        //printf("Blob id =%i \n",a);
        //printf("Original Blob id = %i \n",OldLabel(a));
        //printf("Blob volume (voxels) = %f \n", BlobAverages(0,a));
        for (idx=0; idx<28; idx++){
            fprintf(BLOBLOG,"%.8g ",BlobAverages(idx,a));
        }
        fprintf(BLOBLOG,"\n");
    }
    fclose(BLOBLOG);
    
    double iVol = 1.0/Nx/Ny/Nz;
    sw = 1.0;
    // Compute the Sauter mean grain diamter
    double D = 6.0*Nx*Ny*Nz*(1.0-porosity) / As;
    double pw,pn,pc,awnD,ansD,awsD,JwnD,trJwnD,lwnsDD,cwns;
    pw = paw/vol_w;
    printf("paw = %f \n", paw/vol_w);
    printf("vol_w = %f \n", vol_w);
    
    printf("-----------------------------------------------\n");
    double pwn=0.0;
    vol_n = nwp_volume = 0.0;
    pan = 0.0;
    awn  = ans = lwns = 0.0;
    van(0) = van(1) = van(2) = 0.0;
    vawn(0) = vawn(1) = vawn(2) = 0.0;
    Gwn(0) = Gwn(1) = Gwn(2) = 0.0;
    Gwn(3) = Gwn(4) = Gwn(5) = 0.0;
    Gns(0) = Gns(1) = Gns(2) = 0.0;
    Gns(3) = Gns(4) = Gns(5) = 0.0;
    Jwn = Kwn = efawns = 0.0;
    trJwn = trawn = trRwn = 0.0;    
    
    // Write out the "equilibrium" state with a 0.5 % change in saturation"
    // Always write the largest blob 
    // sw, pw, pn, pc*D/gamma, awn*D, aws*D, ans*D,lwns*D*D, Jwn*D, trJwn*D, cwns, nblobs
    printf("Computing equilibria from blobs \n");
    printf("Sauter mean diamter = %f \n",D);
    printf("WARNING: lazy coder hard-coded the surface tension as 0.058 \n");
    FILE *BLOBSTATES= fopen("blobstates.tcat","a");
    fprintf(BLOBSTATES,"%.5g %.5g %.5g\n",vol_w,pw,aws);
    // Compute the averages over the entire non-wetting phsae
    for (a=0; a<nblobs; a++){
        nwp_volume += BlobAverages(0,a);
        pwn += (BlobAverages(1,a)-pw)*BlobAverages(2,a);
        pan += BlobAverages(1,a)*BlobAverages(28,a);
        awn += BlobAverages(2,a);
        ans += BlobAverages(3,a);
        Jwn += BlobAverages(4,a)*BlobAverages(2,a);
        Kwn += BlobAverages(5,a)*BlobAverages(2,a);
        lwns += BlobAverages(6,a);
        efawns += BlobAverages(7,a)*BlobAverages(6,a);
        van(0) += BlobAverages(8,a)*BlobAverages(28,a);
        van(1) += BlobAverages(9,a)*BlobAverages(28,a);
        van(2) += BlobAverages(10,a)*BlobAverages(28,a);
        vawn(0) += BlobAverages(11,a)*BlobAverages(2,a);
        vawn(1) += BlobAverages(12,a)*BlobAverages(2,a);
        vawn(2) += BlobAverages(13,a)*BlobAverages(2,a);
        Gwn(0) += BlobAverages(14,a)*BlobAverages(2,a);
        Gwn(1) += BlobAverages(15,a)*BlobAverages(2,a);
        Gwn(2) += BlobAverages(16,a)*BlobAverages(2,a);
        Gwn(3) += BlobAverages(17,a)*BlobAverages(2,a);
        Gwn(4) += BlobAverages(18,a)*BlobAverages(2,a);
        Gwn(5) += BlobAverages(19,a)*BlobAverages(2,a);
        Gns(0) += BlobAverages(20,a)*BlobAverages(3,a);
        Gns(1) += BlobAverages(21,a)*BlobAverages(3,a);
        Gns(2) += BlobAverages(22,a)*BlobAverages(3,a);
        Gns(3) += BlobAverages(23,a)*BlobAverages(3,a);
        Gns(4) += BlobAverages(24,a)*BlobAverages(3,a);
        Gns(5) += BlobAverages(25,a)*BlobAverages(3,a);
        trawn += BlobAverages(26,a);
        trJwn += BlobAverages(27,a)*BlobAverages(26,a);
        vol_n += BlobAverages(28,a);
    }    
    
    // Subtract off portions of non-wetting phase in order of size
    for (a=nblobs-1; a>0; a--){
        // Subtract the features one-by-one
        nwp_volume -= BlobAverages(0,a);
        pan -= BlobAverages(1,a)*BlobAverages(28,a);
        pwn -= (BlobAverages(1,a)-pw)*BlobAverages(2,a);
        awn -= BlobAverages(2,a);
        ans -= BlobAverages(3,a);
        Jwn -= BlobAverages(4,a)*BlobAverages(2,a);
        Kwn -= BlobAverages(5,a)*BlobAverages(2,a);
        lwns -= BlobAverages(6,a);
        efawns -= BlobAverages(7,a)*BlobAverages(6,a);
        van(0) -= BlobAverages(8,a)*BlobAverages(28,a);
        van(1) -= BlobAverages(9,a)*BlobAverages(28,a);
        van(2) -= BlobAverages(10,a)*BlobAverages(28,a);
        vawn(0) -= BlobAverages(11,a)*BlobAverages(2,a);
        vawn(1) -= BlobAverages(12,a)*BlobAverages(2,a);
        vawn(2) -= BlobAverages(13,a)*BlobAverages(2,a);
        Gwn(0) -= BlobAverages(14,a)*BlobAverages(2,a);
        Gwn(1) -= BlobAverages(15,a)*BlobAverages(2,a);
        Gwn(2) -= BlobAverages(16,a)*BlobAverages(2,a);
        Gwn(3) -= BlobAverages(17,a)*BlobAverages(2,a);
        Gwn(4) -= BlobAverages(18,a)*BlobAverages(2,a);
        Gwn(5) -= BlobAverages(19,a)*BlobAverages(2,a);
        Gns(0) -= BlobAverages(20,a)*BlobAverages(3,a);
        Gns(1) -= BlobAverages(21,a)*BlobAverages(3,a);
        Gns(2) -= BlobAverages(22,a)*BlobAverages(3,a);
        Gns(3) -= BlobAverages(23,a)*BlobAverages(3,a);
        Gns(4) -= BlobAverages(24,a)*BlobAverages(3,a);
        Gns(5) -= BlobAverages(25,a)*BlobAverages(3,a);
        trawn -= BlobAverages(26,a);
        trJwn -= BlobAverages(27,a)*BlobAverages(26,a);
        vol_n -= BlobAverages(28,a);
        
        // Update wetting phase averages
        aws += BlobAverages(3,a);
        Gws(0) += BlobAverages(20,a)*BlobAverages(3,a);
        Gws(1) += BlobAverages(21,a)*BlobAverages(3,a);
        Gws(2) += BlobAverages(22,a)*BlobAverages(3,a);
        Gws(3) += BlobAverages(23,a)*BlobAverages(3,a);
        Gws(4) += BlobAverages(24,a)*BlobAverages(3,a);
        Gws(5) += BlobAverages(25,a)*BlobAverages(3,a);
        
        if (fabs(1.0 - nwp_volume*iVol/porosity - sw) > 0.0025 || a == 1){
            sw = 1.0 - nwp_volume*iVol/porosity;
            
            JwnD = -Jwn*D/awn;
            pn = pan/vol_n;
            trJwnD = -trJwn*D/trawn;
            cwns = -efawns / lwns;
            awnD = awn*D*iVol;
            awsD = aws*D*iVol;
            ansD = ans*D*iVol;
            lwnsDD = lwns*D*D*iVol;
            pc = pwn*D/0.058/awn;     // hard-coded surface tension due to being lazy
            
            fprintf(BLOBSTATES,"%.5g %.5g %.5g ",sw,pn,pw);
            fprintf(BLOBSTATES,"%.5g %.5g %.5g %.5g ",awnD,awsD,ansD,lwnsDD);
            fprintf(BLOBSTATES,"%.5g %.5g %.5g %.5g %i\n",pc,JwnD,trJwnD,cwns,a);
        }
    }
    fclose(BLOBSTATES);

    start = 0;
    for (a=0;a<nblobs;a++){

        finish = start+b(a);

        for (c=start;c<finish;c++){
            // Get cube from the list
            i = blobs(0,c);
            j = blobs(1,c);
            k = blobs(2,c);
            // Label the entire cube so that interfaces can be re-labled easily
            LocalBlobID(i,j,k) = NewLabel(a);
            LocalBlobID(i+1,j,k) = NewLabel(a);
            LocalBlobID(i,j+1,k) = NewLabel(a);
            LocalBlobID(i+1,j+1,k) = NewLabel(a);
            LocalBlobID(i,j,k+1) = NewLabel(a);
            LocalBlobID(i+1,j,k+1) = NewLabel(a);
            LocalBlobID(i,j+1,k+1) = NewLabel(a);
            LocalBlobID(i+1,j+1,k+1) = NewLabel(a);
        }
        start=finish;
    }

    FILE *BLOBS;
    BLOBS = fopen("Blobs.dat","wb");
    fwrite(LocalBlobID.data(),4,Nx*Ny*Nz,BLOBS);
    fclose(BLOBS);
    
    FILE *DISTANCE;
    DISTANCE = fopen("SignDist.dat","wb");
    fwrite(SignDist.data(),8,Nx*Ny*Nz,DISTANCE);
    fclose(DISTANCE);
    
}

