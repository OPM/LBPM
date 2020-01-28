#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

//#include "common/pmmc.h"
#include "common/Domain.h"
#include "common/SpherePack.h"
#include "common/MPI.h"
#include "common/Communication.h"

/*
 * Pre-Processor to generate signed distance function from sphere packing
 * to use as an input domain for lattice Boltzmann simulator
 * James E. McClure 2014
 */

using namespace std;

inline void PackID(int *list, int count, char *sendbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This packs up the values that need to be sent from one processor to another
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = ID[n];
	}
}
//***************************************************************************************

inline void UnpackID(int *list, int count, char *recvbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This unpacks the values once they have been recieved from neighbors
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		ID[n] = recvbuf[idx];
	}
}


inline void MorphOpen(DoubleArray SignDist, char *id, Domain &Dm, int nx, int ny, int nz, int rank, double SW){
 
	int i,j,k,n;
	double count,countGlobal,totalGlobal;
	count = 0.f;
	double maxdist=0.f;
	double maxdistGlobal;
        int nprocx=Dm.nprocx();
        int nprocy=Dm.nprocy();
        int nprocz=Dm.nprocz();
	for (int k=0; k<nz; k++){
		for (int j=0; j<ny; j++){
			for (int i=0; i<nx; i++){
				n = k*nx*ny+j*nx+i;
				if (SignDist(i,j,k) < 0.0)  id[n] = 0;
				else{
					// initially saturated with wetting phase
					id[n] = 2;
					count+=1.0;
					// extract maximum distance for critical radius
					if ( SignDist(i,j,k) > maxdist) maxdist=SignDist(i,j,k);
				}
			}
		}
	}
	// total Global is the number of nodes in the pore-space
	totalGlobal = Dm.Comm.sumReduce( count );
	maxdistGlobal = Dm.Comm.sumReduce( maxdist );
	double volume=double(nprocx*nprocy*nprocz)*double(nx-2)*double(ny-2)*double(nz-2);
	double porosity=totalGlobal/volume;
	if (rank==0) printf("Media Porosity: %f \n",porosity);
	if (rank==0) printf("Maximum pore size: %f \n",maxdistGlobal);


	// Generate the NWP configuration
	//if (rank==0) printf("Initializing morphological distribution with critical radius %f \n", Rcrit);
	if (rank==0) printf("Performing morphological opening with target saturation %f \n", SW);
	//	GenerateResidual(id,nx,ny,nz,Saturation);

	// Communication buffers
	char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
	char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
	char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
	char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
	char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
	char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
	// send buffers
	sendID_x = new char [Dm.sendCount_x];
	sendID_y = new char [Dm.sendCount_y];
	sendID_z = new char [Dm.sendCount_z];
	sendID_X = new char [Dm.sendCount_X];
	sendID_Y = new char [Dm.sendCount_Y];
	sendID_Z = new char [Dm.sendCount_Z];
	sendID_xy = new char [Dm.sendCount_xy];
	sendID_yz = new char [Dm.sendCount_yz];
	sendID_xz = new char [Dm.sendCount_xz];
	sendID_Xy = new char [Dm.sendCount_Xy];
	sendID_Yz = new char [Dm.sendCount_Yz];
	sendID_xZ = new char [Dm.sendCount_xZ];
	sendID_xY = new char [Dm.sendCount_xY];
	sendID_yZ = new char [Dm.sendCount_yZ];
	sendID_Xz = new char [Dm.sendCount_Xz];
	sendID_XY = new char [Dm.sendCount_XY];
	sendID_YZ = new char [Dm.sendCount_YZ];
	sendID_XZ = new char [Dm.sendCount_XZ];
	//......................................................................................
	// recv buffers
	recvID_x = new char [Dm.recvCount_x];
	recvID_y = new char [Dm.recvCount_y];
	recvID_z = new char [Dm.recvCount_z];
	recvID_X = new char [Dm.recvCount_X];
	recvID_Y = new char [Dm.recvCount_Y];
	recvID_Z = new char [Dm.recvCount_Z];
	recvID_xy = new char [Dm.recvCount_xy];
	recvID_yz = new char [Dm.recvCount_yz];
	recvID_xz = new char [Dm.recvCount_xz];
	recvID_Xy = new char [Dm.recvCount_Xy];
	recvID_xZ = new char [Dm.recvCount_xZ];
	recvID_xY = new char [Dm.recvCount_xY];
	recvID_yZ = new char [Dm.recvCount_yZ];
	recvID_Yz = new char [Dm.recvCount_Yz];
	recvID_Xz = new char [Dm.recvCount_Xz];
	recvID_XY = new char [Dm.recvCount_XY];
	recvID_YZ = new char [Dm.recvCount_YZ];
	recvID_XZ = new char [Dm.recvCount_XZ];
	//......................................................................................
	int sendtag,recvtag;
	sendtag = recvtag = 7;

	int ii,jj,kk;
	int Nx = nx;
	int Ny = ny;
	int Nz = nz;
	
	double sw_old=1.0;
	double sw_new=1.0; 
    double sw_diff_old = 1.0;
    double sw_diff_new = 1.0;

	// Increase the critical radius until the target saturation is met
	double deltaR=0.05; // amount to change the radius in voxel units
	double Rcrit_old;
	double Rcrit_new;

	int imin,jmin,kmin,imax,jmax,kmax;
    
	Rcrit_new = maxdistGlobal;
    while (sw_new > SW)
    {
        sw_diff_old = sw_diff_new;
        sw_old = sw_new;
        Rcrit_old = Rcrit_new;
		Rcrit_new -= deltaR*Rcrit_old;
	    int Window=round(Rcrit_new);
        if (Window == 0) Window = 1; // If Window = 0 at the begining, after the following process will have sw=1.0
                                     // and sw<Sw will be immediately broken
        double LocalNumber=0.f;
        for(k=0; k<Nz; k++){
            for(j=0; j<Ny; j++){
                for(i=0; i<Nx; i++){
                    n = k*nx*ny + j*nx+i;
                    if (SignDist(i,j,k) > Rcrit_new){
                        // loop over the window and update
                        imin=max(1,i-Window);
                        jmin=max(1,j-Window);
                        kmin=max(1,k-Window);
                        imax=min(Nx-1,i+Window);
                        jmax=min(Ny-1,j+Window);
                        kmax=min(Nz-1,k+Window);
                        for (kk=kmin; kk<kmax; kk++){
                            for (jj=jmin; jj<jmax; jj++){
                                for (ii=imin; ii<imax; ii++){
                                    int nn = kk*nx*ny+jj*nx+ii;
                                    double dsq = double((ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k));
                                    if (id[nn] == 2 && dsq <= Rcrit_new*Rcrit_new){
				      LocalNumber+=1.0;
                                        id[nn]=1;
                                    }
                                }
                            }
                        }

                    }
                    // move on
                }
            }
        }
     

        // Pack and send the updated ID values
        PackID(Dm.sendList_x, Dm.sendCount_x ,sendID_x, id);
        PackID(Dm.sendList_X, Dm.sendCount_X ,sendID_X, id);
        PackID(Dm.sendList_y, Dm.sendCount_y ,sendID_y, id);
        PackID(Dm.sendList_Y, Dm.sendCount_Y ,sendID_Y, id);
        PackID(Dm.sendList_z, Dm.sendCount_z ,sendID_z, id);
        PackID(Dm.sendList_Z, Dm.sendCount_Z ,sendID_Z, id);
        PackID(Dm.sendList_xy, Dm.sendCount_xy ,sendID_xy, id);
        PackID(Dm.sendList_Xy, Dm.sendCount_Xy ,sendID_Xy, id);
        PackID(Dm.sendList_xY, Dm.sendCount_xY ,sendID_xY, id);
        PackID(Dm.sendList_XY, Dm.sendCount_XY ,sendID_XY, id);
        PackID(Dm.sendList_xz, Dm.sendCount_xz ,sendID_xz, id);
        PackID(Dm.sendList_Xz, Dm.sendCount_Xz ,sendID_Xz, id);
        PackID(Dm.sendList_xZ, Dm.sendCount_xZ ,sendID_xZ, id);
        PackID(Dm.sendList_XZ, Dm.sendCount_XZ ,sendID_XZ, id);
        PackID(Dm.sendList_yz, Dm.sendCount_yz ,sendID_yz, id);
        PackID(Dm.sendList_Yz, Dm.sendCount_Yz ,sendID_Yz, id);
        PackID(Dm.sendList_yZ, Dm.sendCount_yZ ,sendID_yZ, id);
        PackID(Dm.sendList_YZ, Dm.sendCount_YZ ,sendID_YZ, id);
        //......................................................................................
        MPI_Sendrecv(sendID_x,Dm.sendCount_x,MPI_CHAR,Dm.rank_x(),sendtag,
		     recvID_X,Dm.recvCount_X,MPI_CHAR,Dm.rank_X(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_X,Dm.sendCount_X,MPI_CHAR,Dm.rank_X(),sendtag,
		     recvID_x,Dm.recvCount_x,MPI_CHAR,Dm.rank_x(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_y,Dm.sendCount_y,MPI_CHAR,Dm.rank_y(),sendtag,
		     recvID_Y,Dm.recvCount_Y,MPI_CHAR,Dm.rank_Y(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Y,Dm.sendCount_Y,MPI_CHAR,Dm.rank_Y(),sendtag,
		     recvID_y,Dm.recvCount_y,MPI_CHAR,Dm.rank_y(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_z,Dm.sendCount_z,MPI_CHAR,Dm.rank_z(),sendtag,
		     recvID_Z,Dm.recvCount_Z,MPI_CHAR,Dm.rank_Z(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Z,Dm.sendCount_Z,MPI_CHAR,Dm.rank_Z(),sendtag,
		     recvID_z,Dm.recvCount_z,MPI_CHAR,Dm.rank_z(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xy,Dm.sendCount_xy,MPI_CHAR,Dm.rank_xy(),sendtag,
		     recvID_XY,Dm.recvCount_XY,MPI_CHAR,Dm.rank_XY(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_XY,Dm.sendCount_XY,MPI_CHAR,Dm.rank_XY(),sendtag,
		     recvID_xy,Dm.recvCount_xy,MPI_CHAR,Dm.rank_xy(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Xy,Dm.sendCount_Xy,MPI_CHAR,Dm.rank_Xy(),sendtag,
		     recvID_xY,Dm.recvCount_xY,MPI_CHAR,Dm.rank_xY(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xY,Dm.sendCount_xY,MPI_CHAR,Dm.rank_xY(),sendtag,
		     recvID_Xy,Dm.recvCount_Xy,MPI_CHAR,Dm.rank_Xy(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xz,Dm.sendCount_xz,MPI_CHAR,Dm.rank_xz(),sendtag,
		     recvID_XZ,Dm.recvCount_XZ,MPI_CHAR,Dm.rank_XZ(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_XZ,Dm.sendCount_XZ,MPI_CHAR,Dm.rank_XZ(),sendtag,
		     recvID_xz,Dm.recvCount_xz,MPI_CHAR,Dm.rank_xz(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Xz,Dm.sendCount_Xz,MPI_CHAR,Dm.rank_Xz(),sendtag,
		     recvID_xZ,Dm.recvCount_xZ,MPI_CHAR,Dm.rank_xZ(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_xZ,Dm.sendCount_xZ,MPI_CHAR,Dm.rank_xZ(),sendtag,
		     recvID_Xz,Dm.recvCount_Xz,MPI_CHAR,Dm.rank_Xz(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_yz,Dm.sendCount_yz,MPI_CHAR,Dm.rank_yz(),sendtag,
		     recvID_YZ,Dm.recvCount_YZ,MPI_CHAR,Dm.rank_YZ(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_YZ,Dm.sendCount_YZ,MPI_CHAR,Dm.rank_YZ(),sendtag,
		     recvID_yz,Dm.recvCount_yz,MPI_CHAR,Dm.rank_yz(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_Yz,Dm.sendCount_Yz,MPI_CHAR,Dm.rank_Yz(),sendtag,
		     recvID_yZ,Dm.recvCount_yZ,MPI_CHAR,Dm.rank_yZ(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendID_yZ,Dm.sendCount_yZ,MPI_CHAR,Dm.rank_yZ(),sendtag,
		     recvID_Yz,Dm.recvCount_Yz,MPI_CHAR,Dm.rank_Yz(),recvtag,Dm.Comm.getCommunicator(),MPI_STATUS_IGNORE);
        //......................................................................................
        UnpackID(Dm.recvList_x, Dm.recvCount_x ,recvID_x, id);
        UnpackID(Dm.recvList_X, Dm.recvCount_X ,recvID_X, id);
        UnpackID(Dm.recvList_y, Dm.recvCount_y ,recvID_y, id);
        UnpackID(Dm.recvList_Y, Dm.recvCount_Y ,recvID_Y, id);
        UnpackID(Dm.recvList_z, Dm.recvCount_z ,recvID_z, id);
        UnpackID(Dm.recvList_Z, Dm.recvCount_Z ,recvID_Z, id);
        UnpackID(Dm.recvList_xy, Dm.recvCount_xy ,recvID_xy, id);
        UnpackID(Dm.recvList_Xy, Dm.recvCount_Xy ,recvID_Xy, id);
        UnpackID(Dm.recvList_xY, Dm.recvCount_xY ,recvID_xY, id);
        UnpackID(Dm.recvList_XY, Dm.recvCount_XY ,recvID_XY, id);
        UnpackID(Dm.recvList_xz, Dm.recvCount_xz ,recvID_xz, id);
        UnpackID(Dm.recvList_Xz, Dm.recvCount_Xz ,recvID_Xz, id);
        UnpackID(Dm.recvList_xZ, Dm.recvCount_xZ ,recvID_xZ, id);
        UnpackID(Dm.recvList_XZ, Dm.recvCount_XZ ,recvID_XZ, id);
        UnpackID(Dm.recvList_yz, Dm.recvCount_yz ,recvID_yz, id);
        UnpackID(Dm.recvList_Yz, Dm.recvCount_Yz ,recvID_Yz, id);
        UnpackID(Dm.recvList_yZ, Dm.recvCount_yZ ,recvID_yZ, id);
        UnpackID(Dm.recvList_YZ, Dm.recvCount_YZ ,recvID_YZ, id);
        //......................................................................................

        //double GlobalNumber = Dm.Comm.sumReduce( LocalNumber );

        count = 0.f;
        for (int k=1; k<Nz-1; k++){
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++){
                    n=k*Nx*Ny+j*Nx+i;
                    if (id[n] == 2){
                        count+=1.0;
                    }
                }
            }
        }
        countGlobal = Dm.Comm.sumReduce( count );
        sw_new = countGlobal/totalGlobal;
        sw_diff_new = abs(sw_new-SW);
        // for test only
        if (rank==0){
	  //printf("Final saturation=%f\n",sw_new);
	  // printf("Final critical radius=%f\n",Rcrit_new);
        }
    }

    if (sw_diff_new<sw_diff_old){
        if (rank==0){
            printf("Final saturation=%f\n",sw_new);
            printf("Final critical radius=%f\n",Rcrit_new);

        }
    }
    else{
        if (rank==0){
            printf("Final saturation=%f\n",sw_old);
            printf("Final critical radius=%f\n",Rcrit_old);

        }
    }
	
}


int main(int argc, char **argv)
{
	// Initialize MPI
	MPI_Init(&argc,&argv);
	Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();
	{
		// parallel domain size (# of sub-domains)
		int nprocx,nprocy,nprocz;
		int iproc,jproc,kproc;
		//*****************************************

		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Sphere Packing pre-processor for LBPM-WIA	\n");
			printf("********************************************************\n");
		}

		// Variables that specify the computational domain  
		string filename;
		int Nx,Ny,Nz;		// local sub-domain size
		int nspheres;		// number of spheres in the packing
		double Lx,Ly,Lz;	// Domain length
		double D = 1.0;		// reference length for non-dimensionalization

		int i,j,k,n;
		filename = argv[1];
		auto db = std::make_shared<Database>( filename );
		auto domain_db = db->getDatabase( "Domain" );
		
		auto Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
		Nx = Dm->Nx;
		Ny = Dm->Ny;
		Nz = Dm->Nz;
		Lx = Dm->Lx;
		Ly = Dm->Ly;
		Lz = Dm->Lz;
		iproc = Dm->iproc();
		jproc = Dm->jproc();
		kproc = Dm->kproc();
		nprocx = Dm->nprocx();
		nprocy = Dm->nprocy();
		nprocz = Dm->nprocz();
	        nspheres = domain_db->getScalar<int>( "nspheres");

		//printf("Set domain \n");
		//int BoundaryCondition=1;
		//Nz += 2;
		//Nx = Ny = Nz;	// Cubic domain
		int N = Nx*Ny*Nz;

		// Define Dm.Communication sub-domain -- everywhere
		for (int k=0; k<Nz; k++){
			for (int j=0; j<Ny; j++){
				for (int i=0; i<Nx; i++){
					n = k*Nx*Ny+j*Nx+i;
					Dm->id[n] = 1;
				}
			}
		}
		Dm->CommInit();

		//.......................................................................
		// Filenames used
		char LocalRankString[8];
		char LocalRankFilename[40];
		char LocalRestartFile[40];
		sprintf(LocalRankString,"%05d",rank);
		sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
		sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);

		//	printf("Local File Name =  %s \n",LocalRankFilename);
		// .......... READ THE INPUT FILE .......................................
		//	char value;
		char *id;
		id = new char[N];
		int sum = 0;
		double sum_local;
		double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
		double porosity;
		//...........................................................................
		DoubleArray SignDist(Nx,Ny,Nz);
		//.......................................................................

		// Read in sphere pack
		//if (rank==1) printf("nspheres =%i \n",nspheres);
		//.......................................................................
		double *cx,*cy,*cz,*rad;
		cx = new double[nspheres];
		cy = new double[nspheres];
		cz = new double[nspheres];
		rad = new double[nspheres];
		//.......................................................................
		if (rank == 0)	printf("Reading the sphere packing \n");
		if (rank == 0)	ReadSpherePacking(nspheres,cx,cy,cz,rad);
		comm.barrier();
		// Broadcast the sphere packing to all processes
		comm.bcast(cx,nspheres,0);
		comm.bcast(cy,nspheres,0);
		comm.bcast(cz,nspheres,0);
		comm.bcast(rad,nspheres,0);
		//...........................................................................
		comm.barrier();
		if (rank == 0) cout << "Domain set." << endl;
		if (rank == 0){
			// Compute the Sauter mean diameter
			double totVol = 0.0;
			double totArea = 0.0;
			// Compute the total volume and area of all spheres
			for (i=0; i<nspheres; i++){
				totVol += 1.3333333333333*3.14159265359*rad[i]*rad[i]*rad[i];
				totArea += 4.0*3.14159265359*rad[i]*rad[i];
			}
			D = 6.0*(Nx-2)*nprocx*totVol / totArea / Lx;
			printf("Sauter Mean Diameter (computed from sphere packing) = %f \n",D);
		}
		comm.bcast(&D,1,0);

		//.......................................................................
		SignedDistance(SignDist.data(),nspheres,cx,cy,cz,rad,Lx,Ly,Lz,Nx,Ny,Nz,
				iproc,jproc,kproc,nprocx,nprocy,nprocz);
		//.......................................................................
		// Assign the phase ID field based on the signed distance
		//.......................................................................
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					id[n] = 0;
				}
			}
		}
		sum=0;
		for ( k=1;k<Nz-1;k++){
			for ( j=1;j<Ny-1;j++){
				for ( i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (SignDist(n) > 0.0){ 
						id[n] = 2;	
					}
					// compute the porosity (actual interface location used)
					if (SignDist(n) > 0.0){ 
						sum++;	
					}
				}
			}
		}
		sum_local = 1.0*sum;
		porosity = comm.sumReduce(sum_local);
		porosity = porosity*iVol_global;
		if (rank==0) printf("Media porosity = %f \n",porosity);

		// Run Morphological opening to initialize 50% saturation
		double SW=0.50;
		if (rank==0) printf("MorphOpen: Initializing with saturation %f \n",SW);
		MorphOpen(SignDist, id, *Dm, Nx, Ny, Nz, rank, SW);

		//.........................................................
		// don't perform computations at the eight corners
		id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
		id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;
		//.........................................................

		/*		//.......................................................................
		sprintf(LocalRankString,"%05d",rank);
		sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
		FILE *DIST = fopen(LocalRankFilename,"wb");
		if (DIST==NULL) ERROR("Error opening file: ID.xxxxx");
		fwrite(SignDist.data(),1,N*sizeof(double),DIST);
		fclose(DIST);
		//......................................................................
		*/
		//.......................................................................
		sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);	
		FILE *IDFILE = fopen(LocalRankFilename,"wb");
		if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
		fwrite(id,1,N,IDFILE);
		fclose(IDFILE);
		//......................................................................
	}
	// ****************************************************
	comm.barrier();
	MPI_Finalize();
	// ****************************************************
}
