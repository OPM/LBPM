// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "pmmc.h"
#include "Domain.h"

using namespace std;

//--------------------------------------------------------------------------------------------------------
inline int ComputeBlob(IntArray &blobs, int &nblobs, int &ncubes, IntArray &indicator,
					   DoubleArray &F, DoubleArray &S, double vf, double vs, int startx, int starty,
					   int startz, IntArray &temp)
{
	// Compute the blob (F>vf|S>vs) starting from (i,j,k) - oil blob
	// F>vf => oil phase S>vs => in porespace
	// update the list of blobs, indicator mesh
	int m = F.m;  // maxima for the meshes
	int n = F.n;
	int o = F.o;

	int cubes_in_blob=0;
	int nrecent = 1;						// number of nodes added at most recent sweep
	temp(0,0) = startx;				// Set the initial point as a "seed" for the sweeps
	temp(1,0) = starty;
	temp(2,0) = startz;
	int ntotal = 1;					// total number of nodes in blob
	indicator(startx,starty,startz) = nblobs;

	int p,s,x,y,z,start,finish,nodx,nody,nodz;
	int imin=startx,imax=startx,jmin=starty,jmax=starty;	// initialize maxima / minima
	int kmin=startz,kmax=startz;
	int d[26][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
	{1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},
	{1,0,-1},{-1,0,-1},{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
	{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},
	{-1,-1,1},{-1,-1,-1}};   // directions to neighbors
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	bool status = 1;						// status == true => continue to look for points
	while (status == 1){
		start = ntotal - nrecent;
		finish = ntotal;
		nrecent = 0;						// set recent points back to zero for next sweep through
		for (s=start;s<finish;s++){
			// Loop over recent points; look for new points
			x = temp(0,s);
			y = temp(1,s);
			z = temp(2,s);
			// Looop over the directions
			for (p=0;p<26;p++){
				nodx=x+d[p][0];
				if (nodx < 0 ){ nodx = m-1; }		// Periodic BC for x
				if (nodx > m-1 ){ nodx = 0; }
				nody=y+d[p][1];
				if (nody < 0 ){ nody = n-1; }	// Periodic BC for y
				if (nody > n-1 ){ nody = 0; }
				nodz=z+d[p][2];
				if (nodz < 0 ){ nodz = 0; }		// No periodic BC for z
				if (nodz > o-1 ){ nodz = o-1; }
				if ( F(nodx,nody,nodz) > vf && S(nodx,nody,nodz) > vs
					 && indicator(nodx,nody,nodz) == -1 ){
					// Node is a part of the blob - add it to the list
					temp(0,ntotal) = nodx;
					temp(1,ntotal) = nody;
					temp(2,ntotal) = nodz;
					ntotal++;
					nrecent++;
					// Update the indicator map
					indicator(nodx,nody,nodz) = nblobs;
					// Update the min / max for the cube loop
					if ( nodx < imin ){ imin = nodx; }
					if ( nodx > imax ){ imax = nodx; }
					if ( nody < jmin ){ jmin = nody; }
					if ( nody > jmax ){ jmax = nody; }
					if ( nodz < kmin ){ kmin = nodz; }
					if ( nodz > kmax ){ kmax = nodz; }
				}
				else if (F(nodx,nody,nodz) > vf && S(nodx,nody,nodz) > vs
						 && indicator(nodx,nody,nodz) > -1 &&  indicator(nodx,nody,nodz) != nblobs){
					// Some kind of error in algorithm
					printf("Error in blob search algorithm!");
				}
			}

		}
		if ( nrecent == 0){
			status = 0;
		}
	}
	// Use points in temporary storage array to add cubes to the list of blobs
	if ( imin > 0) { imin = imin-1; }
//	if ( imax < m-1) { imax = imax+1; }
	if ( jmin > 0) { jmin = jmin-1; }
//	if ( jmax < n-1) { jmax = jmax+1; }
	if ( kmin > 0) { kmin = kmin-1; }
//	if ( kmax < o-1) { kmax = kmax+1; }
	int i,j,k;
	bool add;
	for (k=kmin;k<kmax;k++){
		for (j=jmin;j<jmax;j++){
			for (i=imin;i<imax;i++){
				// If cube(i,j,k) has any nodes in blob, add it to the list
				// Loop over cube edges
				add = 0;
				for (p=0;p<8;p++){
					nodx = i+cube[p][0];
					nody = j+cube[p][1];
					nodz = k+cube[p][2];
					if ( indicator(nodx,nody,nodz) == nblobs ){
						// Cube corner is in this blob
						add = 1;
					}
				}
				if (add == 1){
					// Add cube to the list
					blobs(0,ncubes) = i;
					blobs(1,ncubes) = j;
					blobs(2,ncubes) = k;
					ncubes++;
					cubes_in_blob++;
					// Loop again to check for overlap
					for (p=0;p<8;p++){
						nodx = i+cube[p][0];
						nody = j+cube[p][1];
						nodz = k+cube[p][2];
						if (indicator(nodx,nody,nodz) > -1 && indicator(nodx,nody,nodz) != nblobs){
							printf("Overlapping cube!");
							cout << i << ", " << j << ", " << k << endl;
						}
					}
				}
			}
		}
	}

	return cubes_in_blob;
}

inline void ReadFromAllRanks(char *FILENAME, DoubleArray &Phase, DoubleArray &Pressure, DoubleArray &Vel_x, 
							DoubleArray &Vel_y, DoubleArray &Vel_z, int nx, int ny, int nz)
{
	int q,n,N;
	int iglobal,jglobal,kglobal;
	double value;
	double denA,denB;
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double vx,vy,vz;
	
	N = nx*ny*nz;
	
	double *Den, double *DistEven, double *DistOdd,
	
	Den = new double[2*N];
	DistEven = new double[10*N];
	DistOdd = new double[9*N];

	ifstream File(FILENAME,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		File.read((char*) &value, sizeof(value));
		Den[n] = value;
		//	if (n== 66276)	printf("Density a  = %f \n",value);
		File.read((char*) &value, sizeof(value));
		Den[N+n] = value;
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
				f0 = disteven[n];
				f2 = disteven[N+n];
				f4 = disteven[2*N+n];
				f6 = disteven[3*N+n];
				f8 = disteven[4*N+n];
				f10 = disteven[5*N+n];
				f12 = disteven[6*N+n];
				f14 = disteven[7*N+n];
				f16 = disteven[8*N+n];
				f18 = disteven[9*N+n];
				//........................................................................
				f1 = distodd[n];
				f3 = distodd[1*N+n];
				f5 = distodd[2*N+n];
				f7 = distodd[3*N+n];
				f9 = distodd[4*N+n];
				f11 = distodd[5*N+n];
				f13 = distodd[6*N+n];
				f15 = distodd[7*N+n];
				f17 = distodd[8*N+n];
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
				iglobal = iproc*(nx-2)+i-1;
				jglobal = jproc*(ny-2)+j-1;
				kglobal = kproc*(nz-2)+k-1;
				//........................................................................				
				Phase(iglobal,jglobal,kglobal) = (denA+denB)/(denA-denB);
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
	//.......................................................................
	int nprocx,nprocy,nprocz,nprocs;
	int Nx, Ny, Nz;
	int nx,ny,nz;
	int nspheres;
	double Lx,Ly,Lz;
	//.......................................................................
	int i,j,k,n;
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

	nprocs = nprocx*nprocy*nprocz;
	printf("Number of MPI ranks: %i \n", nprocs);
	Nx = (nx-2)*nprocx;
	Ny = (ny-2)*nprocy;
	Nz = (nz-2)*nprocz;
	printf("Full domain size: %i x %i x %i  \n", Nx,Ny,Nz);
	
	// Filenames used
	char LocalRankString[8];
	char LocalRankFilename[40];
	char LocalRestartFile[40];
	char BaseFilename[20];
	char tmpstr[10];
	sprintf(BaseFilename,"%s","dPdt.");
	                   
	DoubleArray Phase(Nx,Ny,Nz);
	DoubleArray SignDist(Nx,Ny,Nz);
	DoubleArray Press(Nx,Ny,Nz);
	DoubleArray Vel_x(Nx,Ny,Nz);			// Velocity
	DoubleArray Vel_y(Nx,Ny,Nz);
	DoubleArray Vel_z(Nx,Ny,Nz);
	
	DoubleArray dPdt(Nx,Ny,Nz);
	
	double * Temp;
	Temp = new double[nx*ny*nz];

	// read the files and populate main arrays
	for ( kproc=0; kproc<nprocz; kproc++){
		for ( jproc=0; jproc<nprocy; jproc++){
			for ( iproc=0; iproc<nprocx; iproc++){
				
				proc = kproc*nprocx*nprocy + jproc*nprocx + iproc;

				sprintf(LocalRankString,"%05d",proc);
				sprintf(LocalRankFilename,"%s%s","dPdt.",LocalRankString);
				printf("Reading file %s \n",LocalRankFilename);
				ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);	
				for (k=1; k<nz-1; k++){
					for (j=1; j<ny-1; j++){
						for (i=1; i<nz-1; i++){
							//........................................................................
							n = k*nx*ny+j*nx+i;
							//........................................................................
							iglobal = iproc*(nx-2)+i-1;
							jglobal = jproc*(ny-2)+j-1;
							kglobal = kproc*(nz-2)+k-1;
							//........................................................................
							SignDist(iglobal,jglobal,kglobal) = Temp[n];
							//........................................................................
						}
					}
				}
				
				sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
				printf("Reading file %s \n",LocalRankFilename);
				ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);	
				for (k=1; k<nz-1; k++){
					for (j=1; j<ny-1; j++){
						for (i=1; i<nz-1; i++){
							//........................................................................
							n = k*nx*ny+j*nx+i;
							//........................................................................
							iglobal = iproc*(nx-2)+i-1;
							jglobal = jproc*(ny-2)+j-1;
							kglobal = kproc*(nz-2)+k-1;
							//........................................................................
							SignDist(iglobal,jglobal,kglobal) = Temp[n];
							//........................................................................
						}
					}
				}
				
				sprintf(LocalRankFilename,"%s%s","Restart.",LocalRankString);

				ReadFromAllRanks(LocalRankFilename,Phase,Pressure,Vel_x,Vel_y,Vel_z,nx,ny,nz);
			}
		}
	}
	printf("Read %i ranks of %s \n",nprocs,BaseFilename);
	delete Temp;
	
	IntArray LocalBlobID(Nx,Ny,Nz);
	DoubleArray MeanCurvature(Nx,Ny,Nz);
	DoubleArray GaussCurvature(Nx,Ny,Nz);
	DoubleArray SignDist_x(Nx,Ny,Nz);		// Gradient of the signed distance
	DoubleArray SignDist_y(Nx,Ny,Nz);
	DoubleArray SignDist_z(Nx,Ny,Nz);
	DoubleArray Phase_x(Nx,Ny,Nz);			// Gradient of the phase indicator field
	DoubleArray Phase_y(Nx,Ny,Nz);
	DoubleArray Phase_z(Nx,Ny,Nz);
	
	// Initialize the local blob ID
	// Initializing the blob ID
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				if (SignDist(i,j,k) < 0.0){
					// Solid phase 
					LocalBlobID(i,j,k) = -2;
				}
			}
		}
	}
	
	/* ****************************************************************
	 VARIABLES FOR THE PMMC ALGORITHM
	 ****************************************************************** */
	//...........................................................................
	// Averaging variables
	//...........................................................................
	double awn,ans,aws,lwns,nwp_volume;
	double sw,vol_n,vol_w,paw,pan;
	double efawns,Jwn;
	double As;
	double dEs,dAwn,dAns;			 // Global surface energy (calculated by rank=0)
	double awn_global,ans_global,aws_global,lwns_global,nwp_volume_global;	
	double As_global;
	//	bool add=1;			// Set to false if any corners contain nw-phase ( F > fluid_isovalue)
	
	int n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
	int n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
	
	double s,s1,s2,s3;		// Triangle sides (lengths)
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
	DoubleArray ContactAngle(20);
	DoubleArray Curvature(20);
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
	// Find blob domains, number of blobs
	int nblobs = 0;					// number of blobs
	int ncubes = 0;					// total number of nodes in any blob
	int N = (Nx-1)*(Ny-1)*(Nz-1);		// total number of nodes
	IntArray blobs(3,N);	// store indices for blobs (cubes)
	IntArray temp(3,N);	// temporary storage array
	IntArray  b(50);		// number of nodes in each blob
	
	std::vector<int> BlobList;
	BlobList.reserve[10000];

	std::vector<int> TempBlobList;
	TempBlobList.reserve[10000];
	
	// Loop over z=0 first -> blobs attached to this end considered "connected" for LB simulation
	i=0;
	int number=0;
	for (k=0;k<1;k++){
		for (j=0;j<Ny;j++){
			if ( F(i,j,k) > vF ){
				if ( S(i,j,k) > vS ){
					// node i,j,k is in the porespace
					number = number+ComputeBlob(blobs,nblobs,ncubes,indicator,F,S,vF,vS,i,j,k,temp);
				}
			}
		}
	}
	// Specify the blob on the z axis
	if (ncubes > 0){
		b(nblobs) = number;
		BlobList.push_back[number];
		printf("Number of blobs is: %i \n",nblobs);
		nblobs++;
	}
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=1;i<Nx;i++){
				if ( indicator(i,j,k) == -1 ){
					if ( F(i,j,k) > vF ){
						if ( S(i,j,k) > vS ){
							// node i,j,k is in the porespace
							b(nblobs) = ComputeBlob(blobs,nblobs,ncubes,indicator,F,S,vF,vS,i,j,k,temp);
							nblobs++;
						}
					}
				}
				// Otherwise, this point has already been assigned - ignore

				// Make sure list blob_nodes is large enough
				if ( nblobs > b.Length-1){
					printf("Increasing size of blob list \n");
					b = IncreaseSize(b,b.Length);
				}
			}
		}
	}
	// Go over all cubes again -> add any that do not contain nw phase
	bool add=1;			// Set to false if any corners contain nw-phase ( F > vF)
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	int count_in=0,count_out=0;
	int nodx,nody,nodz;
	for (k=0;k<Nz-1;k++){
		for (j=0;j<Ny-1;j++){
			for (i=0;i<Nx-1;i++){
				// Loop over cube corners
				add=1;				// initialize to true - add unless corner occupied by nw-phase
				for (p=0;p<8;p++){
					nodx=i+cube[p][0];
					nody=j+cube[p][1];
					nodz=k+cube[p][2];
					if ( indicator(nodx,nody,nodz) > -1 ){
						// corner occupied by nw-phase  -> do not add
						add = 0;
					}
				}
				if ( add == 1 ){
					blobs(0,ncubes) = i;
					blobs(1,ncubes) = j;
					blobs(2,ncubes) = k;
					ncubes++;
					count_in++;
				}
				else { count_out++; }
			}
		}
	}
	b(nblobs) = count_in;
	nblobs++;


	/* ****************************************************************
			RUN TCAT AVERAGING ON EACH BLOB
	****************************************************************** */
	int n_nw_tris_beg, n_ns_tris_beg, n_ws_tris_beg, n_nws_seg_beg;
	int start=0,finish;
	int a,c;
	int newton_steps = 0;
	double blob_volume;
	
	printf("Computing TCAT averages based on connectivity \n");
	printf("The number of blobs is %i \n",nblobs);
	
	// Wetting phase averages assume global connectivity
	vol_w = 0.0;
	vaw(0) = vaw(1) = vaw(2) = 0.0;
	Gws(0) = Gws(1) = Gws(2) = 0.0;
	Gws(3) = Gws(4) = Gws(5) = 0.0;
	
	BLOBLOG= fopen("finalstate.tcat","a");
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
		blob_volume = 0;	// Initialize the volume for blob a to zero			awn = aws = ans = lwns = 0.0;
		nwp_volume = 0.0;
		As = 0.0;
		
		// Compute phase averages
		vol_n =0.0;
		pan = paw = 0.0;
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
			i = cubeList(0,c);
			j = cubeList(1,c);
			k = cubeList(2,c);
	
			// Use the cube to compute volume averages
			for (p=0;p<8;p++){
				if ( SignDist(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){

					// 1-D index for this cube corner
					n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;

					// Compute the non-wetting phase volume contribution
					if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 )
						nwp_volume += 0.125;

					// volume averages over the non-wetting phase
					if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0.99 ){
						// volume the excludes the interfacial region
						vol_n += 0.125;
						// pressure
						pan += 0.125*Press.data[n];
						// velocity
						van(0) += 0.125*Vel_x.data[n];
						van(1) += 0.125*Vel_y.data[n];
						van(2) += 0.125*Vel_z.data[n];
					}

					// volume averages over the wetting phase
					if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) < -0.99 ){
						// volume the excludes the interfacial region
						vol_w += 0.125;
						// pressure
						paw += 0.125*Press.data[n];
						// velocity
						vaw(0) += 0.125*Vel_x.data[n];
						vaw(1) += 0.125*Vel_y.data[n];
						vaw(2) += 0.125*Vel_z.data[n];
					}
				}
			}

			// Interface and common curve averages
			n_local_sol_tris = 0;
			n_local_sol_pts = 0;
			n_local_nws_pts = 0;
	
			//...........................................................................
			// Construct the interfaces and common curve
			pmmc_ConstructLocalCube(SignDist, Phase, solid_isovalue, fluid_isovalue,
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
			n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
			n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
			
			n_nw_tris_beg = n_nw_tris;
			n_ns_tris_beg = n_ns_tris;
			n_ws_tris_beg = n_ws_tris;
			n_nws_seg_beg = n_nws_seg;
			//*******************************************************************
			
			
		}
		start = finish;

		volume(a) = blob_volume;
		ws_areas(a) = aws;
		nw_areas(a) = awn;
		ns_areas(a) = ans;
		
		// Last "blob" is just the ws interface
		if (a+1 < nblobs){
			printf("Blob id = %i \n", a);
			fprintf(BLOBLOG,"%i %.5g ",a);						// blob ID
			fprintf(BLOBLOG,"%.5g ",nwp_volume);				// blob volume
			fprintf(BLOBLOG,"%.5g ",pan);						// blob volume
			fprintf(BLOBLOG,"%.5g %.5g ",awn,ans);				// interfacial areas
			fprintf(BLOBLOG,"%.5g %5g ",Jwn, Kwn);				// curvature of wn interface
			fprintf(BLOBLOG,"%.5g ",lwns);								// common curve length
			fprintf(BLOBLOG,"%.5g ",efawns);											// average contact angle
			fprintf(BLOBLOG,"%.5g %.5g %.5g ",van(0),van(1),van(2));	// average velocity of n phase
			fprintf(BLOBLOG,"%.5g %.5g %.5g ",vawn(0),vawn(1),vawn(2));	// velocity of wn interface
			fprintf(BLOBLOG,"%.5g %.5g %.5g %.5g %.5g %.5g ",
					Gwn(0),Gwn(1),Gwn(2),Gwn(3),Gwn(4),Gwn(5));	// orientation of wn interface
			fprintf(BLOBLOG,"%.5g %.5g %.5g %.5g %.5g %.5g ",
					Gns(0),Gns(1),Gns(2),Gns(3),Gns(4),Gns(5));	// orientation of ns interface
			fprintf(BLOBLOG,"%.5g %5g %5g\n",trawn, trJwn, trRwn);						// Trimmed curvature
		}
		
	}  // End of the blob loop
	
	fprintf(BLOBLOG,"%.5g ", paw);							// blob volume
	fprintf(BLOBLOG,"%.5g %.5g %.5g %.5g %.5g %.5g ",
			Gws(0),Gws(1),Gws(2),Gws(3),Gws(4),Gws(5));	// orientation of ws interface
	fprintf(BLOBLOG,"%.5g %.5g %.5g ",vaw(0),vaw(1),vaw(2));	// average velocity of w phase
	fclose(BLOBLOG);

	
}

