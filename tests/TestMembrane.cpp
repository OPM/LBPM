
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/MPI.h"
#include "common/Membrane.h"
#include "common/ScaLBL.h"

using namespace std;

std::shared_ptr<Database> loadInputs( int nprocs )
{
  //auto db = std::make_shared<Database>( "Domain.in" );
    auto db = std::make_shared<Database>();
    db->putScalar<int>( "BC", 0 );
    db->putVector<int>( "nproc", { 1, 1, 1 } );
    db->putVector<int>( "n", { 32, 32, 32 } );
    db->putScalar<int>( "nspheres", 1 );
    db->putVector<double>( "L", { 1, 1, 1 } );
    return db;
}

//***************************************************************************************
int main(int argc, char **argv)
{
	// Initialize MPI
    Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
	int check=0;
	{

		int i,j,k,n;
		bool Bounceback = false;
		
        int rank = comm.getRank();
		if (rank == 0){
			printf("********************************************************\n");
			printf("Running unit test: TestMembrane	\n");
			printf("********************************************************\n");
		}
		
	    // Load inputs
	    auto db = loadInputs( comm.getSize() );
	    int Nx = db->getVector<int>( "n" )[0];
	    int Ny = db->getVector<int>( "n" )[1];
	    int Nz = db->getVector<int>( "n" )[2];
		auto Dm = std::make_shared<Domain>(db,comm);

		Nx += 2;
		Ny += 2;
		Nz += 2;
		int N = Nx*Ny*Nz;
		//.......................................................................
		int Np = 0;
		double distance,radius;
		DoubleArray Distance(Nx,Ny,Nz);
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					Dm->id[n] = 1;
					radius = double(Nx)/4;
					distance = sqrt(double((i-0.5*Nx)*(i-0.5*Nx)+ (j-0.5*Ny)*(j-0.5*Ny)+ (k-0.5*Nz)*(k-0.5*Nz)))-radius;
					if (distance < 0.0 ){
						Dm->id[n] = 1;
					}
					Distance(i,j,k) = distance;
					Np++;
				}
			}
		}
		Dm->CommInit();

		// Create a communicator for the device (will use optimized layout)
		std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm(new ScaLBL_Communicator(Dm));
		//Create a second communicator based on the regular data layout
		std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm_Regular(new ScaLBL_Communicator(Dm));

		if (rank==0){
			printf("Total domain size = %i \n",N);
			printf("Reduced domain size = %i \n",Np);
		}

		// LBM variables
		if (rank==0)	printf ("Set up the neighborlist \n");
		int Npad=Np+32;
		int neighborSize=18*Npad*sizeof(int);
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Npad];
		
		//......................device distributions.................................
		int *NeighborList;
		int *dvcMap;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Npad);

		Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Dm->id.data(),Np,1);
		comm.barrier();
	    ScaLBL_CopyToDevice(NeighborList, neighborList, 18*Np*sizeof(int));
		
		double *dist;
		dist = new double [19*Np];
		
		// Check the neighborlist
		printf("Check neighborlist: exterior %i, first interior %i last interior %i \n",ScaLBL_Comm->LastExterior(),ScaLBL_Comm->FirstInterior(),ScaLBL_Comm->LastInterior());
		for (int idx=0; idx<ScaLBL_Comm->LastExterior(); idx++){
			for (int q=0; q<18; q++){
				int nn = neighborList[q*Np+idx]%Np;
				if (nn>Np) printf("neighborlist error (exterior) at q=%i, idx=%i \n",q,idx);
				dist[q*Np + idx] = 0.0;
			}
		}
		for (int idx=ScaLBL_Comm->FirstInterior(); idx<ScaLBL_Comm->LastInterior(); idx++){
			for (int q=0; q<18; q++){
				int nn = neighborList[q*Np+idx]%Np;
				if (nn>Np) printf("neighborlist error (exterior) at q=%i, idx=%i \n",q,idx);
				dist[q*Np + idx] = 0.0;
			}
		}
		
		/* create a membrane data structure */
	    Membrane M(ScaLBL_Comm, NeighborList, Np);

	    int MembraneCount = M.Create(Distance, Map);
		if (rank==0)	printf (" Number of membrane links: %i \n", MembraneCount);

		/* create a tagged array to show where the mebrane is*/
		double *MembraneLinks;
		MembraneLinks = new double [Nx*Ny*Nz];
		for (int n=0; n<Nx*Ny*Nz; n++) {
			MembraneLinks[n] = 0.0;
		}
		for (int mlink=0; mlink<MembraneCount; mlink++){
			int iq = M.membraneLinks[2*mlink];
			int jq = M.membraneLinks[2*mlink+1];
			dist[iq] = -1.0; // set these distributions to non-zero
			dist[jq] = 1.0;
		}
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					int idx = Map(i,j,k);
					double sum = 0.0;
					for (int q=0; q<19; q++){
						sum += dist[q*Np + idx];
					}
					int n = k*Nx*Ny + j*Nx + i;
					MembraneLinks[n] = sum;
					if (sum > 0.f){
						Dm->id[n] = 127;
					}
					if (sum < 0.f){
						Dm->id[n] = 64;
					}
				}
			}
		}
		if (argc > 1)
			Dm->AggregateLabels("membrane.raw");
		
		
		/* create a pair of distributions to test membrane mass transport routine */
		double *fq, *gq, *Ci, *Cj, *Psi, *Ci_host;
		Ci_host = new double [Np];
		
	    ScaLBL_AllocateDeviceMemory((void **)&fq,  19 * sizeof(double) * Np);
	    ScaLBL_AllocateDeviceMemory((void **)&gq,  19 * sizeof(double) * Np);
	    ScaLBL_AllocateDeviceMemory((void **)&Ci, sizeof(double) * Np);
	    ScaLBL_AllocateDeviceMemory((void **)&Cj, sizeof(double) * Np);
	    ScaLBL_AllocateDeviceMemory((void **)&Psi, sizeof(double) * Np);
		
		/* initialize concentration inside membrane */
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					int idx = Map(i,j,k);
					if (Distance(i,j,k) > 0.0)
						Ci_host[idx] = 1.0;
					else 
						Ci_host[idx] = 0.0;
				}
			}
		}
        ScaLBL_CopyToDevice(Ci, Ci_host,  sizeof(double) * Np);
        
        /* initialize the distributions */
        ScaLBL_D3Q7_Ion_Init_FromFile(fq, Ci, Np);
        ScaLBL_D3Q7_Ion_Init_FromFile(gq, Ci, Np);
        
        /*  Streaming with the usual neighborlist */
        ScaLBL_D3Q19_AAodd_Compact(NeighborList, fq, Np);
        
        /* Streaming with the membrane neighborlist*/
        ScaLBL_D3Q19_AAodd_Compact(M.NeighborList, gq, Np);

        /* explicit mass transfer step with the membrane*/
        M.AssignCoefficients(dvcMap, Psi, 0.0, 1.0, 1.0, 1.0, 1.0);
        M.IonTransport(gq, Cj);
        ScaLBL_CopyToHost(Ci_host, Cj, sizeof(double) * Np);
        
        double ionError = 0.0;
        for (int n=0; n<Np; n++){
        	ionError += Ci_host[n];
        }
        if (fabs(ionError) > 1e-12) {
        	printf(" Failed error tolerance in membrane ion transport routine! \n");
        	check = 2;
        }
        
        DoubleArray Ions(Nx,Ny,Nz);
        ScaLBL_Comm->RegularLayout(Map, Cj, Ions);
		if (argc > 1)
			Dm->AggregateLabels("membrane2.raw",Ions);
		
        /* now compare streaming */
        ScaLBL_D3Q7_Ion_Init_FromFile(gq, Ci, Np);
        M.IonTransport(gq, Cj);
        ScaLBL_D3Q19_AAodd_Compact(M.NeighborList, gq, Np);
        M.IonTransport(gq, Cj);
        
        /* now check that the two results agree*/
        double *fq_h, *gq_h;
        fq_h = new double [7*Np];
        gq_h = new double [7*Np];
        ScaLBL_CopyToHost(fq_h, fq, 7*sizeof(double) * Np);
        ScaLBL_CopyToHost(gq_h, gq, 7*sizeof(double) * Np);
        for (int n = 0; n<Np; n++){
        	for (int q=0; q<7; q++){
        		double gval = gq_h[q*Np + n];
        		double fval = fq_h[q*Np + n];
        		if (gval != fval ){
        			printf("  Membrane streaming mismatch at q=%i, n=%i \n",q,n);
        			printf("  .... gq = %f, fq = %f \n",gval, fval);
        			printf(" (unit test will fail) \n");
        			check = 3;
        		}
        	}
        }

        DoubleArray MembraneErrors(Nx,Ny,Nz);
        for (k=1;k<Nz-1;k++){
        	for (j=1;j<Ny-1;j++){
        		for (i=1;i<Nx-1;i++){
        			n = k*Nx*Ny+j*Nx+i;
        			int idx = Map(i,j,k);
        			MembraneErrors(i,j,k) = 0.0;
        			for (int q=0; q<7; q++){
        				double gval = gq_h[q*Np + idx];
        				double fval = fq_h[q*Np + idx];			
        				MembraneErrors(i,j,k) += gval - fval;
        			}
        		}
        	}
        }
		Dm->AggregateLabels("membrane3.raw",MembraneErrors);

        
		//...........................................................................
		// Update GPU data structures
		if (rank==0)	printf ("Setting up device map and neighbor list \n");
		int *TmpMap;
		TmpMap=new int[Np*sizeof(int)];
		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int idx=Map(i,j,k);
					if (!(idx < 0))
						TmpMap[idx] = k*Nx*Ny+j*Nx+i;
				}
			}
		}
		ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
		ScaLBL_DeviceBarrier();
		
		// Create a dummy distribution data structure
		double *fq_host;
		fq_host = new double[19*Np];
		if (rank==0)	printf ("Setting up Np=%i distributions \n",Np);
		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int idx=Map(i,j,k);
					if (!(idx<0)){
					  for (int q=0; q<19; q++){
					    fq_host[q*Np+idx]=(k*Nx*Ny+j*Nx+i)+0.01*q;
					  }
					}
				}
			}
		}

		/* Run dummy communications */
		/*initialize fq from host data */
		ScaLBL_CopyToDevice(fq, fq_host, sizeof(double)*7*Np);
		
		M.SendD3Q7AA(&fq[0]);
		M.RecvD3Q7AA(&gq[0],Bounceback);
		// this has only the communicated values
		//ScaLBL_CopyToHost(fq_host, gq, sizeof(double)*7*Np);
		if (rank==0)	printf ("Sum result \n");
		
        ScaLBL_D3Q7_AAeven_IonConcentration(&gq[0 * Np * 7], &Ci[0 * Np],
                                            0, ScaLBL_Comm->LastExterior(),
                                            Np);
        DoubleArray Result(Nx,Ny,Nz);
        
        ScaLBL_Comm->RegularLayout(Map, Ci, Result);

/*		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int idx=Map(i,j,k);
					double sum = 0.0;
					if (!(idx<0)){
					  for (int q=1; q<3; q++){
					    sum += fq_host[q*Np+idx];
					  }
					  Result[k*Nx*Ny+j*Nx+i] = sum;
					}
				}
			}
		}
*/
		FILE *OUTFILE;
		OUTFILE = fopen("D3Q7.raw","wb");
		fwrite(Result.data(),8,Nx*Ny*Nz,OUTFILE);
		fclose(OUTFILE);
		
		FILE *MAPFILE;
		MAPFILE = fopen("Map.raw","wb");
		fwrite(Map.data(),4,Nx*Ny*Nz,MAPFILE);
		fclose(MAPFILE);

		delete [] TmpMap;
		delete [] fq_host;

	}
    Utilities::shutdown();

	return check;
}

