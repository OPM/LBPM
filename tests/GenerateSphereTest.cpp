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

inline void PackID(const int *list, int count, char *sendbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This packs up the values that need to be sent from one processor to another
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = ID[n];
	}
}
//***************************************************************************************

inline void UnpackID(const int *list, int count, char *recvbuf, char *ID){
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
	sendID_x = new char [Dm.sendCount("x")];
	sendID_y = new char [Dm.sendCount("y")];
	sendID_z = new char [Dm.sendCount("z")];
	sendID_X = new char [Dm.sendCount("X")];
	sendID_Y = new char [Dm.sendCount("Y")];
	sendID_Z = new char [Dm.sendCount("Z")];
	sendID_xy = new char [Dm.sendCount("xy")];
	sendID_yz = new char [Dm.sendCount("yz")];
	sendID_xz = new char [Dm.sendCount("xz")];
	sendID_Xy = new char [Dm.sendCount("Xy")];
	sendID_Yz = new char [Dm.sendCount("Yz")];
	sendID_xZ = new char [Dm.sendCount("xZ")];
	sendID_xY = new char [Dm.sendCount("xY")];
	sendID_yZ = new char [Dm.sendCount("yZ")];
	sendID_Xz = new char [Dm.sendCount("Xz")];
	sendID_XY = new char [Dm.sendCount("XY")];
	sendID_YZ = new char [Dm.sendCount("YZ")];
	sendID_XZ = new char [Dm.sendCount("XZ")];
	//......................................................................................
	// recv buffers
	recvID_x = new char [Dm.recvCount("x")];
	recvID_y = new char [Dm.recvCount("y")];
	recvID_z = new char [Dm.recvCount("z")];
	recvID_X = new char [Dm.recvCount("X")];
	recvID_Y = new char [Dm.recvCount("Y")];
	recvID_Z = new char [Dm.recvCount("Z")];
	recvID_xy = new char [Dm.recvCount("xy")];
	recvID_yz = new char [Dm.recvCount("yz")];
	recvID_xz = new char [Dm.recvCount("xz")];
	recvID_Xy = new char [Dm.recvCount("Xy")];
	recvID_xZ = new char [Dm.recvCount("xZ")];
	recvID_xY = new char [Dm.recvCount("xY")];
	recvID_yZ = new char [Dm.recvCount("yZ")];
	recvID_Yz = new char [Dm.recvCount("Yz")];
	recvID_Xz = new char [Dm.recvCount("Xz")];
	recvID_XY = new char [Dm.recvCount("XY")];
	recvID_YZ = new char [Dm.recvCount("YZ")];
	recvID_XZ = new char [Dm.recvCount("XZ")];
	//......................................................................................
	int sendtag,recvtag;
	sendtag = recvtag = 7;

	int ii,jj,kk;
	int Nx = nx;
	int Ny = ny;
	int Nz = nz;
	int imin,jmin,kmin,imax,jmax,kmax;
    
	
	double sw_old=1.0;
	double sw_new=1.0; 
	double sw_diff_old = 1.0;
	double sw_diff_new = 1.0;

	// Increase the critical radius until the target saturation is met
	double deltaR=0.05; // amount to change the radius in voxel units
	double Rcrit_new = maxdistGlobal;
	double Rcrit_old = maxdistGlobal;
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
        PackID(Dm.sendList("x"), Dm.sendCount("x") ,sendID_x, id);
        PackID(Dm.sendList("X"), Dm.sendCount("X") ,sendID_X, id);
        PackID(Dm.sendList("y"), Dm.sendCount("y") ,sendID_y, id);
        PackID(Dm.sendList("Y"), Dm.sendCount("Y") ,sendID_Y, id);
        PackID(Dm.sendList("z"), Dm.sendCount("z") ,sendID_z, id);
        PackID(Dm.sendList("Z"), Dm.sendCount("Z") ,sendID_Z, id);
        PackID(Dm.sendList("xy"), Dm.sendCount("xy") ,sendID_xy, id);
        PackID(Dm.sendList("Xy"), Dm.sendCount("Xy") ,sendID_Xy, id);
        PackID(Dm.sendList("xY"), Dm.sendCount("xY") ,sendID_xY, id);
        PackID(Dm.sendList("XY"), Dm.sendCount("XY") ,sendID_XY, id);
        PackID(Dm.sendList("xz"), Dm.sendCount("xz") ,sendID_xz, id);
        PackID(Dm.sendList("Xz"), Dm.sendCount("Xz") ,sendID_Xz, id);
        PackID(Dm.sendList("xZ"), Dm.sendCount("xZ") ,sendID_xZ, id);
        PackID(Dm.sendList("XZ"), Dm.sendCount("XZ") ,sendID_XZ, id);
        PackID(Dm.sendList("yz"), Dm.sendCount("yz") ,sendID_yz, id);
        PackID(Dm.sendList("Yz"), Dm.sendCount("Yz") ,sendID_Yz, id);
        PackID(Dm.sendList("yZ"), Dm.sendCount("yZ") ,sendID_yZ, id);
        PackID(Dm.sendList("YZ"), Dm.sendCount("YZ") ,sendID_YZ, id);
        //......................................................................................
        Dm.Comm.sendrecv(sendID_x,Dm.sendCount("x"),Dm.rank_x(),sendtag,recvID_X,Dm.recvCount("X"),Dm.rank_X(),recvtag);
        Dm.Comm.sendrecv(sendID_X,Dm.sendCount("X"),Dm.rank_X(),sendtag,recvID_x,Dm.recvCount("x"),Dm.rank_x(),recvtag);
        Dm.Comm.sendrecv(sendID_y,Dm.sendCount("y"),Dm.rank_y(),sendtag,recvID_Y,Dm.recvCount("Y"),Dm.rank_Y(),recvtag);
        Dm.Comm.sendrecv(sendID_Y,Dm.sendCount("Y"),Dm.rank_Y(),sendtag,recvID_y,Dm.recvCount("y"),Dm.rank_y(),recvtag);
        Dm.Comm.sendrecv(sendID_z,Dm.sendCount("z"),Dm.rank_z(),sendtag,recvID_Z,Dm.recvCount("Z"),Dm.rank_Z(),recvtag);
        Dm.Comm.sendrecv(sendID_Z,Dm.sendCount("Z"),Dm.rank_Z(),sendtag,recvID_z,Dm.recvCount("z"),Dm.rank_z(),recvtag);
        Dm.Comm.sendrecv(sendID_xy,Dm.sendCount("xy"),Dm.rank_xy(),sendtag,recvID_XY,Dm.recvCount("XY"),Dm.rank_XY(),recvtag);
        Dm.Comm.sendrecv(sendID_XY,Dm.sendCount("XY"),Dm.rank_XY(),sendtag,recvID_xy,Dm.recvCount("xy"),Dm.rank_xy(),recvtag);
        Dm.Comm.sendrecv(sendID_Xy,Dm.sendCount("Xy"),Dm.rank_Xy(),sendtag,recvID_xY,Dm.recvCount("xY"),Dm.rank_xY(),recvtag);
        Dm.Comm.sendrecv(sendID_xY,Dm.sendCount("xY"),Dm.rank_xY(),sendtag,recvID_Xy,Dm.recvCount("Xy"),Dm.rank_Xy(),recvtag);
        Dm.Comm.sendrecv(sendID_xz,Dm.sendCount("xz"),Dm.rank_xz(),sendtag,recvID_XZ,Dm.recvCount("XZ"),Dm.rank_XZ(),recvtag);
        Dm.Comm.sendrecv(sendID_XZ,Dm.sendCount("XZ"),Dm.rank_XZ(),sendtag,recvID_xz,Dm.recvCount("xz"),Dm.rank_xz(),recvtag);
        Dm.Comm.sendrecv(sendID_Xz,Dm.sendCount("Xz"),Dm.rank_Xz(),sendtag,recvID_xZ,Dm.recvCount("xZ"),Dm.rank_xZ(),recvtag);
        Dm.Comm.sendrecv(sendID_xZ,Dm.sendCount("xZ"),Dm.rank_xZ(),sendtag,recvID_Xz,Dm.recvCount("Xz"),Dm.rank_Xz(),recvtag);
        Dm.Comm.sendrecv(sendID_yz,Dm.sendCount("yz"),Dm.rank_yz(),sendtag,recvID_YZ,Dm.recvCount("YZ"),Dm.rank_YZ(),recvtag);
        Dm.Comm.sendrecv(sendID_YZ,Dm.sendCount("YZ"),Dm.rank_YZ(),sendtag,recvID_yz,Dm.recvCount("yz"),Dm.rank_yz(),recvtag);
        Dm.Comm.sendrecv(sendID_Yz,Dm.sendCount("Yz"),Dm.rank_Yz(),sendtag,recvID_yZ,Dm.recvCount("yZ"),Dm.rank_yZ(),recvtag);
        Dm.Comm.sendrecv(sendID_yZ,Dm.sendCount("yZ"),Dm.rank_yZ(),sendtag,recvID_Yz,Dm.recvCount("Yz"),Dm.rank_Yz(),recvtag);
        //......................................................................................
        UnpackID(Dm.recvList("x"), Dm.recvCount("x") ,recvID_x, id);
        UnpackID(Dm.recvList("X"), Dm.recvCount("X") ,recvID_X, id);
        UnpackID(Dm.recvList("y"), Dm.recvCount("y") ,recvID_y, id);
        UnpackID(Dm.recvList("Y"), Dm.recvCount("Y") ,recvID_Y, id);
        UnpackID(Dm.recvList("z"), Dm.recvCount("z") ,recvID_z, id);
        UnpackID(Dm.recvList("Z"), Dm.recvCount("Z") ,recvID_Z, id);
        UnpackID(Dm.recvList("xy"), Dm.recvCount("xy") ,recvID_xy, id);
        UnpackID(Dm.recvList("Xy"), Dm.recvCount("Xy") ,recvID_Xy, id);
        UnpackID(Dm.recvList("xY"), Dm.recvCount("xY") ,recvID_xY, id);
        UnpackID(Dm.recvList("XY"), Dm.recvCount("XY") ,recvID_XY, id);
        UnpackID(Dm.recvList("xz"), Dm.recvCount("xz") ,recvID_xz, id);
        UnpackID(Dm.recvList("Xz"), Dm.recvCount("Xz") ,recvID_Xz, id);
        UnpackID(Dm.recvList("xZ"), Dm.recvCount("xZ") ,recvID_xZ, id);
        UnpackID(Dm.recvList("XZ"), Dm.recvCount("XZ") ,recvID_XZ, id);
        UnpackID(Dm.recvList("yz"), Dm.recvCount("yz") ,recvID_yz, id);
        UnpackID(Dm.recvList("Yz"), Dm.recvCount("Yz") ,recvID_Yz, id);
        UnpackID(Dm.recvList("yZ"), Dm.recvCount("yZ") ,recvID_yZ, id);
        UnpackID(Dm.recvList("YZ"), Dm.recvCount("YZ") ,recvID_YZ, id);
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
    Utilities::startup( argc, argv );
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
    Utilities::shutdown();
}
