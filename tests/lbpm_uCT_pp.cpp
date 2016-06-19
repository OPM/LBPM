// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>

#include "common/Array.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "common/MPI_Helpers.h"
#include "common/imfilter.h"
#include "IO/MeshDatabase.h"
#include "IO/Mesh.h"
#include "IO/Writer.h"
#include "IO/netcdf.h"
#include "analysis/analysis.h"
#include "analysis/eikonal.h"

#include "ProfilerApp.h"


template<class T>
inline int sign( T x )
{
    if ( x==0 )
        return 0;
    return x>0 ? 1:-1;
}


inline void Med3D( const Array<float> &Input, Array<float> &Output )
{
	PROFILE_START("Med3D");
	// Perform a 3D Median filter on Input array with specified width
	int i,j,k,ii,jj,kk;
	int imin,jmin,kmin,imax,jmax,kmax;

	float *List;
	List=new float[27];

	int Nx = int(Input.size(0));
	int Ny = int(Input.size(1));
	int Nz = int(Input.size(2));

	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){

				// Just use a 3x3x3 window (hit recursively if needed)
				imin = i-1;
				jmin = j-1;
				kmin = k-1;
				imax = i+2;
				jmax = j+2;
				kmax = k+2;

				// Populate the list with values in the window
				int Number=0;
				for (kk=kmin; kk<kmax; kk++){
					for (jj=jmin; jj<jmax; jj++){
						for (ii=imin; ii<imax; ii++){
							List[Number++] = Input(ii,jj,kk);
						}
					}
				}
				// Sort the first 5 entries and return the median
				for (ii=0; ii<14; ii++){
					for (jj=ii+1; jj<27; jj++){
						if (List[jj] < List[ii]){
							float tmp = List[ii];
							List[ii] = List[jj];
							List[jj] = tmp;
						}
					}
				}
				// Return the median
				Output(i,j,k) = List[13];
			}
		}
	}
	PROFILE_STOP("Med3D");
}


inline float trilinear( float dx, float dy, float dz, float f1, float f2,
    float f3, float f4, float f5, float f6, float f7, float f8 )
{
    double f, dx2, dy2, dz2, h0, h1;
    dx2 = 1.0 - dx;
    dy2 = 1.0 - dy;
    dz2 = 1.0 - dz;
    h0  = ( dx * f2 + dx2 * f1 ) * dy2 + ( dx * f4 + dx2 * f3 ) * dy;
    h1  = ( dx * f6 + dx2 * f5 ) * dy2 + ( dx * f8 + dx2 * f7 ) * dy;
    f   = h0 * dz2 + h1 * dz;
    return ( f );
}
inline void InterpolateMesh( const Array<float> &Coarse, Array<float> &Fine )
{
	PROFILE_START("InterpolateMesh");

	// Interpolate values from a Coarse mesh to a fine one
	// This routine assumes cell-centered meshes with 1 ghost cell

	// Fine mesh
	int Nx = int(Fine.size(0))-2;
	int Ny = int(Fine.size(1))-2;
	int Nz = int(Fine.size(2))-2;

	// Coarse mesh
	int nx = int(Coarse.size(0))-2;
	int ny = int(Coarse.size(1))-2;
	int nz = int(Coarse.size(2))-2;

	// compute the stride
	int hx = Nx/nx;
	int hy = Ny/ny;
	int hz = Nz/nz;
    ASSERT(nx*hx==Nx);
    ASSERT(ny*hy==Ny);
    ASSERT(nz*hz==Nz);

	// value to map distance between meshes (since distance is in voxels)
	//  usually hx=hy=hz (or something very close)
	//  the mapping is not exact
	//  however, it's assumed the coarse solution will be refined
	//  a good guess is the goal here!
	float mapvalue = sqrt(hx*hx+hy*hy+hz*hz);

	// Interpolate to the fine mesh
	for (int k=-1; k<Nz+1; k++){
        int k0 = floor((k-0.5*hz)/hz);
        int k1 = k0+1;
        int k2 = k0+2;
        float dz = ( (k+0.5) - (k0+0.5)*hz ) / hz;
        ASSERT(k0>=-1&&k0<nz+1&&dz>=0&&dz<=1);
		for (int j=-1; j<Ny+1; j++){
            int j0 = floor((j-0.5*hy)/hy);
            int j1 = j0+1;
            int j2 = j0+2;
            float dy = ( (j+0.5) - (j0+0.5)*hy ) / hy;
            ASSERT(j0>=-1&&j0<ny+1&&dy>=0&&dy<=1);
			for (int i=-1; i<Nx+1; i++){
                int i0 = floor((i-0.5*hx)/hx);
                int i1 = i0+1;
                int i2 = i0+2;
                float dx = ( (i+0.5) - (i0+0.5)*hx ) / hx;
                ASSERT(i0>=-1&&i0<nx+1&&dx>=0&&dx<=1);
                float val = trilinear( dx, dy, dz,
                    Coarse(i1,j1,k1), Coarse(i2,j1,k1), Coarse(i1,j2,k1), Coarse(i2,j2,k1),
                    Coarse(i1,j1,k2), Coarse(i2,j1,k2), Coarse(i1,j2,k2), Coarse(i2,j2,k2) );
                Fine(i+1,j+1,k+1) = mapvalue*val;
			}
		}
	}
	PROFILE_STOP("InterpolateMesh");
}


inline int NLM3D( const Array<float> &Input, Array<float> &Mean, 
    const Array<float> &Distance, Array<float> &Output, const int d, const float h)
{
	PROFILE_START("NLM3D");
	// Implemenation of 3D non-local means filter
	// 		d determines the width of the search volume
	// 		h is a free parameter for non-local means (i.e. 1/sigma^2)
	// 		Distance is the signed distance function
	// 		If Distance(i,j,k) > THRESHOLD_DIST then don't compute NLM

	float THRESHOLD_DIST = float(d);
	float weight, sum;
	int i,j,k,ii,jj,kk;
	int imin,jmin,kmin,imax,jmax,kmax;
	int returnCount=0;

	int Nx = int(Input.size(0));
	int Ny = int(Input.size(1));
	int Nz = int(Input.size(2));

	// Compute the local means
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){

				imin = max(0,i-d);
				jmin = max(0,j-d);
				kmin = max(0,k-d);
				imax = min(Nx-1,i+d);
				jmax = min(Ny-1,j+d);
				kmax = min(Nz-1,k+d);

				// Populate the list with values in the window
				sum = 0; weight=0;
				for (kk=kmin; kk<kmax; kk++){
					for (jj=jmin; jj<jmax; jj++){
						for (ii=imin; ii<imax; ii++){
							sum += Input(ii,jj,kk);
							weight++;
						}
					}
				}

				Mean(i,j,k) = sum / weight;
			}
		}
	}

	// Compute the non-local means
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){


				if (fabs(Distance(i,j,k)) < THRESHOLD_DIST){
					// compute the expensive non-local means
					sum = 0; weight=0;

					imin = max(0,i-d);
					jmin = max(0,j-d);
					kmin = max(0,k-d);
					imax = min(Nx-1,i+d);
					jmax = min(Ny-1,j+d);
					kmax = min(Nz-1,k+d);

					for (kk=kmin; kk<kmax; kk++){
						for (jj=jmin; jj<jmax; jj++){
							for (ii=imin; ii<imax; ii++){
								float tmp = Mean(i,j,k) - Mean(ii,jj,kk);
								sum += exp(-tmp*tmp*h)*Input(ii,jj,kk);
								weight += exp(-tmp*tmp*h);
							}
						}
					}

					returnCount++;
					//Output(i,j,k) = Mean(i,j,k);
					Output(i,j,k) = sum / weight;
				}
				else{
					// Just return the mean
					Output(i,j,k) = Mean(i,j,k);
				}
			}
		}
	}
	// Return the number of sites where NLM was applied
	PROFILE_STOP("NLM3D");
	return returnCount;
}


// Reading the domain information file
void read_domain( int rank, int nprocs, MPI_Comm comm, 
    int& nprocx, int& nprocy, int& nprocz, int& nx, int& ny, int& nz,
    int& nspheres, double& Lx, double& Ly, double& Lz )
{
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
	MPI_Barrier(comm);
	// Computational domain
	//.................................................
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
	MPI_Barrier(comm);
}



// Smooth the data using the distance
void smooth( const Array<float>& VOL, const Array<float>& Dist, float sigma, Array<float>& MultiScaleSmooth, fillHalo<float>& fillFloat )
{
    for (size_t i=0; i<VOL.length(); i++) {
		// use exponential weight based on the distance
		float dst = Dist(i);
		float tmp = exp(-(dst*dst)/(sigma*sigma));
		float value = dst>0 ? -1:1;
		MultiScaleSmooth(i) = tmp*VOL(i) + (1-tmp)*value;
	}
	fillFloat.fill(MultiScaleSmooth);
}


// Segment the data
void segment( const Array<float>& data, Array<char>& ID, float tol )
{
    ASSERT(data.size()==ID.size());
    for (size_t i=0; i<data.length(); i++) {
        if ( data(i) > tol )
            ID(i) = 0;
        else
            ID(i) = 1;
    }
}


// Remove disconnected phases
void removeDisconnected( Array<char>& ID, const Domain& Dm )
{
    // Run blob identification to remove disconnected volumes
    BlobIDArray GlobalBlobID;
    DoubleArray SignDist(ID.size());
    DoubleArray Phase(ID.size());
    for (size_t i=0; i<ID.length(); i++) {
        SignDist(i) = (2*ID(i)-1);
        Phase(i) = 1;
    }
    ComputeGlobalBlobIDs( ID.size(0)-2, ID.size(1)-2, ID.size(2)-2,
        Dm.rank_info, Phase, SignDist, 0, 0, GlobalBlobID, Dm.Comm );
    for (size_t i=0; i<ID.length(); i++) {
        if ( GlobalBlobID(i) > 0 )
            ID(i) = 0;
        ID(i) = GlobalBlobID(i);
    }
}


// Solve a level (without any coarse level information)
void solve( const Array<float>& VOL, Array<float>& Mean, Array<char>& ID,
    Array<float>& Dist, Array<float>& MultiScaleSmooth, Array<float>& NonLocalMean, 
    fillHalo<float>& fillFloat, const Domain& Dm, int nprocx )
{
    PROFILE_SCOPED(timer,"solve");
    // Compute the median filter on the sparse array
    Med3D( VOL, Mean );
    fillFloat.fill( Mean );
    segment( Mean, ID, 0.01 );
    // Compute the distance using the segmented volume
	Eikonal3D( Dist, ID, Dm, ID.size(0)*nprocx );
	fillFloat.fill(Dist);
    smooth( VOL, Dist, 2.0, MultiScaleSmooth, fillFloat );
    // Compute non-local mean
	int depth = 5;
	float sigsq=0.1;
	int nlm_count = NLM3D( MultiScaleSmooth, Mean, Dist, NonLocalMean, depth, sigsq);
	fillFloat.fill(NonLocalMean);
}


// Refine a solution from a coarse grid to a fine grid
void refine( const Array<float>& Dist_coarse, 
    const Array<float>& VOL, Array<float>& Mean, Array<char>& ID,
    Array<float>& Dist, Array<float>& MultiScaleSmooth, Array<float>& NonLocalMean, 
    fillHalo<float>& fillFloat, const Domain& Dm, int nprocx, int level )
{
    PROFILE_SCOPED(timer,"refine");
    int ratio[3] = { int(Dist.size(0)/Dist_coarse.size(0)),
                     int(Dist.size(1)/Dist_coarse.size(1)),
                     int(Dist.size(2)/Dist_coarse.size(2)) };
    // Interpolate the distance from the coarse to fine grid
	InterpolateMesh( Dist_coarse, Dist );
    // Compute the median filter on the array and segment
    Med3D( VOL, Mean );
    fillFloat.fill( Mean );
    segment( Mean, ID, 0.01 );
    // If the ID has the wrong distance, set the distance to 0 and run a simple filter to set neighbors to 0
    for (size_t i=0; i<ID.length(); i++) {
        char id = Dist(i)>0 ? 1:0;
        if ( id != ID(i) )
            Dist(i) = 0;
    }
    fillFloat.fill( Dist );
    std::function<float(int,const float*)> filter_1D = []( int N, const float* data )
    {
        bool zero = data[0]==0 || data[2]==0;
        return zero ? data[1]*1e-12 : data[1];
    };
    std::vector<imfilter::BC> BC(3,imfilter::BC::replicate);
    std::vector<std::function<float(int,const float*)>> filter_set(3,filter_1D);
    Dist = imfilter::imfilter_separable<float>( Dist, {1,1,1}, filter_set, BC );
    fillFloat.fill( Dist );
    // Smooth the volume data
    float lambda = 2*sqrt(double(ratio[0]*ratio[0]+ratio[1]*ratio[1]+ratio[2]*ratio[2]));
    smooth( VOL, Dist, lambda, MultiScaleSmooth, fillFloat );
    // Compute non-local mean
	int depth = 3;
	float sigsq = 0.1;
	int nlm_count = NLM3D( MultiScaleSmooth, Mean, Dist, NonLocalMean, depth, sigsq);
	fillFloat.fill(NonLocalMean);
    segment( NonLocalMean, ID, 0.001 );
    for (size_t i=0; i<ID.length(); i++) {
        char id = Dist(i)>0 ? 1:0;
        if ( id!=ID(i) || fabs(Dist(i))<1 )
    	    Dist(i) = 2.0*ID(i)-1.0;
    }
    // Remove disconnected domains
    //removeDisconnected( ID, Dm );
    // Compute the distance using the segmented volume
    if ( level > 0 ) {
    	//Eikonal3D( Dist, ID, Dm, ID.size(0)*nprocx );
        //CalcDist3D( Dist, ID, Dm );
        CalcDistMultiLevel( Dist, ID, Dm );
	    fillFloat.fill(Dist);
    }
}


// Remove regions that are likely noise by shrinking the volumes by dx,
// removing all values that are more than dx+delta from the surface, and then
// growing by dx+delta and intersecting with the original data
void filter_final( Array<char>& ID, Array<float>& Dist,
    fillHalo<float>& fillFloat, const Domain& Dm,
    Array<float>& Mean, Array<float>& Dist1, Array<float>& Dist2 )
{
    PROFILE_SCOPED(timer,"filter_final");
	int rank;
	MPI_Comm_rank(Dm.Comm,&rank);
    int Nx = Dm.Nx-2;
    int Ny = Dm.Ny-2;
    int Nz = Dm.Nz-2;
    // Calculate the distance
    CalcDistMultiLevel( Dist, ID, Dm );
    fillFloat.fill(Dist);
    // Compute the range to shrink the volume based on the L2 norm of the distance
    Array<float> Dist0(Nx,Ny,Nz);
    fillFloat.copy(Dist,Dist0);
    float tmp = 0;
    for (size_t i=0; i<Dist0.length(); i++)
        tmp += Dist0(i)*Dist0(i);
    tmp = sqrt( sumReduce(Dm.Comm,tmp) / sumReduce(Dm.Comm,(float)Dist0.length()) );
    const float dx1 = 0.3*tmp;
    const float dx2 = 1.05*dx1;
    if (rank==0)
        printf("   %0.1f %0.1f %0.1f\n",tmp,dx1,dx2);
    // Update the IDs/Distance removing regions that are < dx of the range
    Dist1 = Dist;
    Dist2 = Dist;
    Array<char> ID1 = ID;
    Array<char> ID2 = ID;
    for (size_t i=0; i<ID.length(); i++) {
        ID1(i) = Dist(i)<-dx1 ? 1:0;
        ID2(i) = Dist(i)> dx1 ? 1:0;
    }
    //Array<float> Dist1 = Dist;
    //Array<float> Dist2 = Dist;
    CalcDistMultiLevel( Dist1, ID1, Dm );
    CalcDistMultiLevel( Dist2, ID2, Dm );
    fillFloat.fill(Dist1);
    fillFloat.fill(Dist2);
    // Keep those regions that are within dx2 of the new volumes
    Mean = Dist;
    for (size_t i=0; i<ID.length(); i++) {
        if ( Dist1(i)+dx2>0 && ID(i)<=0 ) {
            Mean(i) = -1;
        } else if ( Dist2(i)+dx2>0 && ID(i)>0 ) {
            Mean(i) = 1;
        } else {
            Mean(i) = Dist(i)>0 ? 0.5:-0.5;
        }
    }
    // Find regions of uncertainty that are entirely contained within another region
    fillHalo<double> fillDouble(Dm.Comm,Dm.rank_info,Nx,Ny,Nz,1,1,1,0,1);
    fillHalo<BlobIDType> fillInt(Dm.Comm,Dm.rank_info,Nx,Ny,Nz,1,1,1,0,1);
    BlobIDArray GlobalBlobID;
    DoubleArray SignDist(ID.size());
    for (size_t i=0; i<ID.length(); i++)
        SignDist(i) = fabs(Mean(i))==1 ? -1:1;
    fillDouble.fill(SignDist);
    DoubleArray Phase(ID.size());
    Phase.fill(1);
    ComputeGlobalBlobIDs( Nx, Ny, Nz, Dm.rank_info, Phase, SignDist, 0, 0, GlobalBlobID, Dm.Comm );
    fillInt.fill(GlobalBlobID);
    int N_blobs = maxReduce(Dm.Comm,GlobalBlobID.max()+1);
    std::vector<float> mean(N_blobs,0);
    std::vector<int> count(N_blobs,0);
    for (int k=1; k<=Nz; k++) {
        for (int j=1; j<=Ny; j++) {
            for (int i=1; i<=Nx; i++) {
                int id = GlobalBlobID(i,j,k);
                if ( id >= 0 ) {
                    if ( GlobalBlobID(i-1,j,k)<0 ) {
                        mean[id] += Mean(i-1,j,k);
                        count[id]++;
                    }
                    if ( GlobalBlobID(i+1,j,k)<0 ) {
                        mean[id] += Mean(i+1,j,k);
                        count[id]++;
                    }
                    if ( GlobalBlobID(i,j-1,k)<0 ) {
                        mean[id] += Mean(i,j-1,k);
                        count[id]++;
                    }
                    if ( GlobalBlobID(i,j+1,k)<0 ) {
                        mean[id] += Mean(i,j+1,k);
                        count[id]++;
                    }
                    if ( GlobalBlobID(i,j,k-1)<0 ) {
                        mean[id] += Mean(i,j,k-1);
                        count[id]++;
                    }
                    if ( GlobalBlobID(i,j,k+1)<0 ) {
                        mean[id] += Mean(i,j,k+1);
                        count[id]++;
                    }
                }
            }
        }
    }
    mean = sumReduce(Dm.Comm,mean);
    count = sumReduce(Dm.Comm,count);
    for (size_t i=0; i<mean.size(); i++)
        mean[i] /= count[i];
    /*if (rank==0) {
        for (size_t i=0; i<mean.size(); i++)
            printf("%i %0.4f\n",i,mean[i]);
    }*/
    for (size_t i=0; i<Mean.length(); i++) {
        int id = GlobalBlobID(i);
        if ( id >= 0 ) {
            if ( fabs(mean[id]) > 0.95 ) {
                // Isolated domain surrounded by one domain
                GlobalBlobID(i) = -2;
                Mean(i) = sign(mean[id]);
            } else {
                // Boarder volume, set to liquid
                Mean(i) = 1;
            }
        }
    }
    // Perform the final segmentation and update the distance
    fillFloat.fill(Mean);
    segment( Mean, ID, 0.01 );
    CalcDistMultiLevel( Dist, ID, Dm );
    fillFloat.fill(Dist);
}


int main(int argc, char **argv)
{

	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
    Utilities::setErrorHandlers();
	PROFILE_START("Main");

	//std::vector<std::string> filenames;
	if ( argc<2 ) {
        if ( rank == 0 )
    		printf("At least one filename must be specified\n");
		return 1;
	}
	std::string filename = std::string(argv[1]);
    if ( rank == 0 )
		printf("Input data file: %s\n",filename.c_str());

	//.......................................................................
	// Reading the domain information file
	//.......................................................................
	int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
	double Lx, Ly, Lz;
    read_domain( rank, nprocs, comm, nprocx, nprocy, nprocz, nx, ny, nz, nspheres, Lx, Ly, Lz );
	int BC=0;


	// Check that the number of processors >= the number of ranks
	if ( rank==0 ) {
		printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
		printf("Number of MPI ranks used: %i \n", nprocs);
		printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
	}
	if ( nprocs < nprocx*nprocy*nprocz ){
		ERROR("Insufficient number of processors");
	}

    // Determine the maximum number of levels for the desired coarsen ratio
    int ratio[3] = {4,4,4};
    std::vector<int> Nx(1,nx), Ny(1,ny), Nz(1,nz);
    while ( Nx.back()%ratio[0]==0 && Nx.back()>8 &&
            Ny.back()%ratio[1]==0 && Ny.back()>8 &&
            Nz.back()%ratio[2]==0 && Nz.back()>8 )
    {
        Nx.push_back( Nx.back()/ratio[0] );
        Ny.push_back( Ny.back()/ratio[1] );
        Nz.push_back( Nz.back()/ratio[2] );
    }
    int N_levels = Nx.size();

	// Initialize the domain
	std::vector<std::shared_ptr<Domain>> Dm(N_levels);
    for (int i=0; i<N_levels; i++) {
        Dm[i].reset( new Domain(Nx[i],Ny[i],Nz[i],rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC) );
        int N = (Nx[i]+2)*(Ny[i]+2)*(Nz[i]+2);
    	for (int n=0; n<N; n++)
			Dm[i]->id[n] = 1;
    	Dm[i]->CommInit(comm);
	}

    // array containing a distance mask
    Array<float> MASK(Nx[0]+2,Ny[0]+2,Nz[0]+2);
    
    // Create the level data
    std::vector<Array<char>>  ID(N_levels);
    std::vector<Array<float>> LOCVOL(N_levels);
    std::vector<Array<float>> Dist(N_levels);
    std::vector<Array<float>> MultiScaleSmooth(N_levels);
    std::vector<Array<float>> Mean(N_levels);
    std::vector<Array<float>> NonLocalMean(N_levels);
	std::vector<std::shared_ptr<fillHalo<double>>> fillDouble(N_levels);
	std::vector<std::shared_ptr<fillHalo<float>>>  fillFloat(N_levels);
	std::vector<std::shared_ptr<fillHalo<char>>>   fillChar(N_levels);
    for (int i=0; i<N_levels; i++) {
        ID[i] = Array<char>(Nx[i]+2,Ny[i]+2,Nz[i]+2);
        LOCVOL[i] = Array<float>(Nx[i]+2,Ny[i]+2,Nz[i]+2);
        Dist[i] = Array<float>(Nx[i]+2,Ny[i]+2,Nz[i]+2);
        MultiScaleSmooth[i] = Array<float>(Nx[i]+2,Ny[i]+2,Nz[i]+2);
        Mean[i] = Array<float>(Nx[i]+2,Ny[i]+2,Nz[i]+2);
        NonLocalMean[i] = Array<float>(Nx[i]+2,Ny[i]+2,Nz[i]+2);
        ID[i].fill(0);
        LOCVOL[i].fill(0);
        Dist[i].fill(0);
        MultiScaleSmooth[i].fill(0);
        Mean[i].fill(0);
        NonLocalMean[i].fill(0);
	    fillDouble[i].reset(new fillHalo<double>(Dm[i]->Comm,Dm[i]->rank_info,Nx[i],Ny[i],Nz[i],1,1,1,0,1) );
	    fillFloat[i].reset(new fillHalo<float>(Dm[i]->Comm,Dm[i]->rank_info,Nx[i],Ny[i],Nz[i],1,1,1,0,1) );
	    fillChar[i].reset(new fillHalo<char>(Dm[i]->Comm,Dm[i]->rank_info,Nx[i],Ny[i],Nz[i],1,1,1,0,1) );
	}


    // Read the subvolume of interest on each processor
	PROFILE_START("ReadVolume");
    int fid = netcdf::open(filename);
    std::string varname("VOLUME");
    netcdf::VariableType type = netcdf::getVarType( fid, varname );
    std::vector<size_t> dim = netcdf::getVarDim( fid, varname );
    if ( rank == 0 ) {
    	printf("Reading %s (%s)\n",varname.c_str(),netcdf::VariableTypeName(type).c_str());
	    printf("   dims =  %i x %i x %i \n",int(dim[0]),int(dim[1]),int(dim[2]));
    }
    {
        RankInfoStruct info( rank, nprocx, nprocy, nprocz );
	    int x = info.ix*nx;
	    int y = info.jy*ny;
	    int z = info.kz*nz;
        // Read the local data
		Array<short> VOLUME = netcdf::getVar<short>( fid, varname, {x,y,z}, {nx,ny,nz}, {1,1,1} );
        // Copy the data and fill the halos
        LOCVOL[0].fill(0);
        fillFloat[0]->copy( VOLUME, LOCVOL[0] );
        fillFloat[0]->fill( LOCVOL[0] );
    }
    netcdf::close( fid );
	MPI_Barrier(comm);
	PROFILE_STOP("ReadVolume");
	if (rank==0) printf("Read complete\n");


    // Filter the original data
	PROFILE_START("Filter source data");
    {
        // Perform a hot-spot filter on the data
        std::vector<imfilter::BC> BC = { imfilter::BC::replicate, imfilter::BC::replicate, imfilter::BC::replicate };
        std::function<float(const Array<float>&)> filter_3D = []( const Array<float>& data )
        {
            float min1 = std::min(data(0,1,1),data(2,1,1));
            float min2 = std::min(data(1,0,1),data(1,2,1));
            float min3 = std::min(data(1,1,0),data(1,1,2));
            float max1 = std::max(data(0,1,1),data(2,1,1));
            float max2 = std::max(data(1,0,1),data(1,2,1));
            float max3 = std::max(data(1,1,0),data(1,1,2));
            float min = std::min(min1,std::min(min2,min3));
            float max = std::max(max1,std::max(max2,max3));
            return std::max(std::min(data(1,1,1),max),min);
        };
        std::function<float(const Array<float>&)> filter_1D = []( const Array<float>& data )
        {
            float min = std::min(data(0),data(2));
            float max = std::max(data(0),data(2));
            return std::max(std::min(data(1),max),min);
        };
        //LOCVOL[0] = imfilter::imfilter<float>( LOCVOL[0], {1,1,1}, filter_3D, BC );
        std::vector<std::function<float(const Array<float>&)>> filter_set(3,filter_1D);
        LOCVOL[0] = imfilter::imfilter_separable<float>( LOCVOL[0], {1,1,1}, filter_set, BC );
        fillFloat[0]->fill( LOCVOL[0] );
        // Perform a gaussian filter on the data
        int Nh[3] = { 2, 2, 2 };
        float sigma[3] = { 1.0, 1.0, 1.0 };
        std::vector<Array<float>> H(3);
        H[0] = imfilter::create_filter<float>( { Nh[0] }, "gaussian", &sigma[0] );
        H[1] = imfilter::create_filter<float>( { Nh[1] }, "gaussian", &sigma[1] );
        H[2] = imfilter::create_filter<float>( { Nh[2] }, "gaussian", &sigma[2] );
        LOCVOL[0] = imfilter::imfilter_separable( LOCVOL[0], H, BC );
        fillFloat[0]->fill( LOCVOL[0] );
	}
	PROFILE_STOP("Filter source data");

	
	// Set up the mask to be distance to cylinder (crop outside cylinder)
	float CylRad=900;
	for (int k=0;k<Nz[0]+2;k++) {
		for (int j=0;j<Ny[0]+2;j++) {
			for (int i=0;i<Nx[0]+2;i++) {
				
				int iproc = Dm[0]->iproc;
				int jproc = Dm[0]->jproc;
				int kproc = Dm[0]->kproc;
				
				int x=iproc*Nx[0]+i-1;
				int y=jproc*Ny[0]+j-1;
				int z=kproc*Nz[0]+k-1;
				
				int cx = 0.5*nprocx*Nx[0];
				int cy = 0.5*nprocy*Ny[0];
				int cz = 0.5*nprocz*Nz[0];

				// distance from the center line 
				MASK(i,j,k) = CylRad - sqrt(float((z-cz)*(z-cz) + (y-cy)*(y-cy)) );
				
			}
		}
	}

	// Compute the means for the high/low regions
	// (should use automated mixture model to approximate histograms)
	float THRESHOLD = 0.05*maxReduce( Dm[0]->Comm, std::max( LOCVOL[0].max(), fabs(LOCVOL[0].min()) ) );
	float mean_plus=0;
	float mean_minus=0;
	int count_plus=0;
	int count_minus=0;
	for (int k=1;k<Nz[0]+1;k++) {
		for (int j=1;j<Ny[0]+1;j++) {
			for (int i=1;i<Nx[0]+1;i++) {
				if (MASK(i,j,k) > 0.f ){
					auto tmp = LOCVOL[0](i,j,k);
					if ( tmp > THRESHOLD ) {
						mean_plus += tmp;
						count_plus++;
					} else if ( tmp < -THRESHOLD ) {
						mean_minus += tmp;
						count_minus++;
					}
				}
			}
		}
	}
    mean_plus = sumReduce( Dm[0]->Comm, mean_plus ) / sumReduce( Dm[0]->Comm, count_plus );
    mean_minus = sumReduce( Dm[0]->Comm, mean_minus ) / sumReduce( Dm[0]->Comm, count_minus );
	if (rank==0) printf("	Region 1 mean (+): %f, Region 2 mean (-): %f \n",mean_plus, mean_minus);


    // Scale the source data to +-1.0
    for (size_t i=0; i<LOCVOL[0].length(); i++) {
    	if (MASK(i) < 0.f){
    		LOCVOL[0](i) = 1.0;
    	}
    	else if ( LOCVOL[0](i) >= 0 ) {
            LOCVOL[0](i) /= mean_plus;
            LOCVOL[0](i) = std::min( LOCVOL[0](i), 1.0f );
        } else {
            LOCVOL[0](i) /= -mean_minus;
            LOCVOL[0](i) = std::max( LOCVOL[0](i), -1.0f );
        }
    }


	// Fill the source data for the coarse meshes
    PROFILE_START("CoarsenMesh");
    for (int i=1; i<N_levels; i++) {
        Array<float> filter(ratio[0],ratio[1],ratio[2]);
        filter.fill(1.0f/filter.length());
        Array<float> tmp(Nx[i-1],Ny[i-1],Nz[i-1]);
        fillFloat[i-1]->copy( LOCVOL[i-1], tmp );
        Array<float> coarse = tmp.coarsen( filter );
        fillFloat[i]->copy( coarse, LOCVOL[i] );
        fillFloat[i]->fill( LOCVOL[i] );
    }
    PROFILE_STOP("CoarsenMesh");


    // Initialize the coarse level
    PROFILE_START("Solve coarse mesh");
    if (rank==0)
        printf("Initialize coarse mesh\n");
    solve( LOCVOL.back(), Mean.back(), ID.back(), Dist.back(), MultiScaleSmooth.back(),
        NonLocalMean.back(), *fillFloat.back(), *Dm.back(), nprocx );
    PROFILE_STOP("Solve coarse mesh");

    // Refine the solution
    PROFILE_START("Refine distance");
    if (rank==0)
        printf("Refine mesh\n");
    for (int i=int(Nx.size())-2; i>=0; i--) {
        if (rank==0)
            printf("   Refining to level %i\n",int(i));
        refine( Dist[i+1], LOCVOL[i], Mean[i], ID[i], Dist[i], MultiScaleSmooth[i],
            NonLocalMean[i], *fillFloat[i], *Dm[i], nprocx, i );
    }
    PROFILE_STOP("Refine distance");

    // Perform a final filter
    PROFILE_START("Filtering final domains");
    if (rank==0)
        printf("Filtering final domains\n");
    Array<float> filter_Mean, filter_Dist1, filter_Dist2;
    filter_final( ID[0], Dist[0], *fillFloat[0], *Dm[0], filter_Mean, filter_Dist1, filter_Dist2 );
    PROFILE_STOP("Filtering final domains");

//removeDisconnected( ID[0], *Dm[0] );

    // Write the results to visit
	if (rank==0) printf("Writing output \n");
	std::vector<IO::MeshDataStruct> meshData(N_levels);
    for (size_t i=0; i<Nx.size(); i++) {
        // Mesh
    	meshData[i].meshName = "Level " + std::to_string(i+1);
    	meshData[i].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(Dm[i]->rank_info,Nx[i],Ny[i],Nz[i],Lx,Ly,Lz) );
        // Source data
        std::shared_ptr<IO::Variable> OrigData( new IO::Variable() );
	    OrigData->name = "Source Data";
	    OrigData->type = IO::VolumeVariable;
	    OrigData->dim = 1;
	    OrigData->data.resize(Nx[i],Ny[i],Nz[i]);
	    meshData[i].vars.push_back(OrigData);
        fillDouble[i]->copy( LOCVOL[i], OrigData->data );
        // Non-Local Mean
        std::shared_ptr<IO::Variable> NonLocMean( new IO::Variable() );
	    NonLocMean->name = "Non-Local Mean";
	    NonLocMean->type = IO::VolumeVariable;
	    NonLocMean->dim = 1;
	    NonLocMean->data.resize(Nx[i],Ny[i],Nz[i]);
	    meshData[i].vars.push_back(NonLocMean);
        fillDouble[i]->copy( NonLocalMean[i], NonLocMean->data );
        // Segmented Data
        std::shared_ptr<IO::Variable> SegData( new IO::Variable() );
	    SegData->name = "Segmented Data";
	    SegData->type = IO::VolumeVariable;
	    SegData->dim = 1;
	    SegData->data.resize(Nx[i],Ny[i],Nz[i]);
	    meshData[i].vars.push_back(SegData);
        fillDouble[i]->copy( ID[i], SegData->data );
        // Signed Distance
        std::shared_ptr<IO::Variable> DistData( new IO::Variable() );
	    DistData->name = "Signed Distance";
	    DistData->type = IO::VolumeVariable;
	    DistData->dim = 1;
	    DistData->data.resize(Nx[i],Ny[i],Nz[i]);
	    meshData[i].vars.push_back(DistData);
        fillDouble[i]->copy( Dist[i], DistData->data );
        // Smoothed Data
        std::shared_ptr<IO::Variable> SmoothData( new IO::Variable() );
	    SmoothData->name = "Smoothed Data";
	    SmoothData->type = IO::VolumeVariable;
	    SmoothData->dim = 1;
	    SmoothData->data.resize(Nx[i],Ny[i],Nz[i]);
	    meshData[i].vars.push_back(SmoothData);
        fillDouble[i]->copy( MultiScaleSmooth[i], SmoothData->data );
    }
    #if 0
        std::shared_ptr<IO::Variable> filter_Mean_var( new IO::Variable() );
	    filter_Mean_var->name = "Mean";
	    filter_Mean_var->type = IO::VolumeVariable;
	    filter_Mean_var->dim = 1;
	    filter_Mean_var->data.resize(Nx[0],Ny[0],Nz[0]);
	    meshData[0].vars.push_back(filter_Mean_var);
        fillDouble[0]->copy( filter_Mean, filter_Mean_var->data );
        std::shared_ptr<IO::Variable> filter_Dist1_var( new IO::Variable() );
	    filter_Dist1_var->name = "Dist1";
	    filter_Dist1_var->type = IO::VolumeVariable;
	    filter_Dist1_var->dim = 1;
	    filter_Dist1_var->data.resize(Nx[0],Ny[0],Nz[0]);
	    meshData[0].vars.push_back(filter_Dist1_var);
        fillDouble[0]->copy( filter_Dist1, filter_Dist1_var->data );
        std::shared_ptr<IO::Variable> filter_Dist2_var( new IO::Variable() );
	    filter_Dist2_var->name = "Dist2";
	    filter_Dist2_var->type = IO::VolumeVariable;
	    filter_Dist2_var->dim = 1;
	    filter_Dist2_var->data.resize(Nx[0],Ny[0],Nz[0]);
	    meshData[0].vars.push_back(filter_Dist2_var);
        fillDouble[0]->copy( filter_Dist2, filter_Dist2_var->data );
    #endif

	// Write visulization data
	IO::writeData( 0, meshData, 2, comm );
	if (rank==0) printf("Finished. \n");

	/* for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
			        n = k*nx*ny+j*nx+i;
			        if (Dm.id[n]==char(SOLID))     Dm.id[n] = 0;
			       	else if (Dm.id[n]==char(NWP))  Dm.id[n] = 1;
			       	else                           Dm.id[n] = 2;

			}
		}
	}
	if (rank==0) printf("Domain set \n");

	// Write the local volume files
	char LocalRankString[8];
	char LocalRankFilename[40];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"Seg.%s",LocalRankString);
	FILE * SEG;
	SEG=fopen(LocalRankFilename,"wb");
	fwrite(LOCVOL.get(),4,N,SEG);
	fclose(SEG);
	 */

	PROFILE_STOP("Main");
    PROFILE_SAVE("lbpm_uCT_pp",true);
	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;
}

