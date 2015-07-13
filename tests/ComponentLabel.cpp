// Sequential component labeling for two phase systems
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2015

#include <iostream>
#include <math.h>
#include "common/pmmc.h"
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
	//	if (n== 66276)	printf("Density a  = %f \n",value);
		File.read((char*) &value, sizeof(value));
		cDen[N+n] = value;
	//	if (n== 66276)	printf("Density b  = %f \n",value);
		// Read the even distributions
		for (q=0; q<10; q++){
			File.read((char*) &value, sizeof(value));
			cDistEven[q*N+n] = value;
	//		if (n== 66276)	printf("dist even %i  = %f \n",q,value);
		}
		// Read the odd distributions
		for (q=0; q<9; q++){
			File.read((char*) &value, sizeof(value));
			cDistOdd[q*N+n] = value;
	//		if (n== 66276)	printf("dist even %i  = %f \n",q,value);
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
		//	if (n== 66276)	printf("Density a  = %f \n",value);
		File.read((char*) &value, sizeof(value));
		Den[2*n+1] = value;

		//	if (n== 66276)	printf("Density b  = %f \n",value);
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
	DoubleArray Vel_x(Nx,Ny,Nz);			// Velocity
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
	// read the files and populate main arrays
	for ( kproc=0; kproc<nprocz; kproc++){
		for ( jproc=0; jproc<nprocy; jproc++){
			for ( iproc=0; iproc<nprocx; iproc++){
				
				proc = kproc*nprocx*nprocy + jproc*nprocx + iproc;

				sprintf(LocalRankString,"%05d",proc);
				sprintf(LocalRankFilename,"%s%s","dPdt.",LocalRankString);
	//			printf("Reading file %s \n",LocalRankFilename);
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
		//		printf("Reading file %s \n",LocalRankFilename);
		//		printf("Sub-domain size: %i x %i x %i  \n", nx,ny,nz);
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
	
	delete Temp;
	
	IntArray PhaseLabel(Nx,Ny,Nz); // label solid (0), non-wetting phase (1), wetting phase (2)
	IntArray WP(Nx,Ny,Nz);	// labels for wetting phase components
	IntArray NWP(Nx,Ny,Nz); // labels for non-wetting phase components
	DoubleArray MeanCurvature(Nx,Ny,Nz);
	DoubleArray GaussCurvature(Nx,Ny,Nz);
	DoubleArray SignDist_x(Nx,Ny,Nz);		// Gradient of the signed distance
	DoubleArray SignDist_y(Nx,Ny,Nz);
	DoubleArray SignDist_z(Nx,Ny,Nz);
	DoubleArray Phase_x(Nx,Ny,Nz);			// Gradient of the phase indicator field
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
	
	FILE *PHASE = fopen("Phase.dat","wb");
	fwrite(Phase.get(),8,Nx*Ny*Nz,PHASE);
	fclose(PHASE);
	
	// Initializing the blob ID
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				if (SignDist(i,j,k) < 0.0){
					// Solid phase 
					PhaseLabel(i,j,k) = 0;
					//WP(i,j,k) = -2;
					//NWP(i,j,k) = -2;
				}
				else if (Phase(i,j,k) < 0){
					// non-wetting phase
					PhaseLabel(i,j,k) = 2;
					//WP(i,j,k) = -2;
					//NWP(i,j,k) = -1;
				}
				else {
					// wetting phase
					PhaseLabel(i,j,k) = 1;
					//WP(i,j,k) = -1;
					//NWP(i,j,k) = -2;
				}
			}
		}
	}
	
	// Compute the porosity
	double porosity=0.0;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				if (PhaseLabel(i,j,k) == 0){
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
	//  double dEs,dAwn,dAns;	// Global surface energy (calculated by rank=0)
	//	bool add=1;			// Set to false if any corners contain nw-phase ( F > fluid_isovalue)
	
	int n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0;
	int n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
	
	//double s,s1,s2,s3;		// Triangle sides (lengths)
	Point A,B,C,P;
	//	double area;
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	//	int count_in=0,count_out=0;
	//	int nodx,nody,nodz;
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
				IDENTIFY ALL COMPONENTS FOR BOTH PHASES
	****************************************************************** */
    int number_NWP_components = ComputeLocalPhaseComponent(PhaseLabel,1,NWP,true);
    int number_WP_components = ComputeLocalPhaseComponent(PhaseLabel,2,WP,true);

    printf("Number of WP components = %i \n",number_WP_components);
    printf("Number of NWP components = %i \n",number_NWP_components);

	DoubleArray BlobAverages(NUM_AVERAGES,number_NWP_components);
	
	// Map the signed distance for the analysis
	for (i=0; i<Nx*Ny*Nz; i++)	SignDist(i) -= (1.0); 
	
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

	FILE *NWP_FILE;
	NWP_FILE = fopen("NWP.dat","wb");
	fwrite(NWP.get(),4,Nx*Ny*Nz,NWP_FILE);
	fclose(NWP_FILE);

	FILE *WP_FILE;
	WP_FILE = fopen("WP.dat","wb");
	fwrite(WP.get(),4,Nx*Ny*Nz,WP_FILE);
	fclose(WP_FILE);
	
	FILE *DISTANCE;
	DISTANCE = fopen("SignDist.dat","wb");
	fwrite(SignDist.get(),8,Nx*Ny*Nz,DISTANCE);
	fclose(DISTANCE);
	
}

