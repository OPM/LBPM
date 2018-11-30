// Unit test to test mass conservation for D3Q7 mass transport LBE
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/MPI_Helpers.h"
#include "models/ColorModel.h"

inline void InitializeBubble(ScaLBL_ColorModel &ColorModel, double BubbleRadius){
	// initialize a bubble
	int i,j,k,n;
	int rank = ColorModel.Mask->rank();
	int nprocx = ColorModel.Mask->nprocx();
	int nprocy = ColorModel.Mask->nprocy();
	int nprocz = ColorModel.Mask->nprocz();
	int Nx = ColorModel.Mask->Nx;
	int Ny = ColorModel.Mask->Ny;
	int Nz = ColorModel.Mask->Nz;
	if (rank == 0) cout << "Setting up bubble..." << endl;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;
				ColorModel.Averages->SDs(i,j,k) = 100.f;
			}
		}
	}
	// Initialize the bubble
	int count_in_bubble=0;
	int count_out_bubble=0;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;
				int iglobal= i+(Nx-2)*ColorModel.Mask->iproc();
				int jglobal= j+(Ny-2)*ColorModel.Mask->jproc();
				int kglobal= k+(Nz-2)*ColorModel.Mask->kproc();
				// Initialize phase position field for parallel bubble test
				if (jglobal < 40){
					ColorModel.Mask->id[n] = 0;
				}
				else if ((iglobal-0.5*(Nx-2)*nprocx)*(iglobal-0.5*(Nx-2)*nprocx)
						+(jglobal-0.5*(Ny-2)*nprocy)*(jglobal-0.5*(Ny-2)*nprocy)
						+(kglobal-0.5*(Nz-2)*nprocz)*(kglobal-0.5*(Nz-2)*nprocz) < BubbleRadius*BubbleRadius){
					ColorModel.Mask->id[n] = 2;
					ColorModel.Mask->id[n] = 2;
					count_in_bubble++;
				}
				else{
					ColorModel.Mask->id[n]=1;
					ColorModel.Mask->id[n]=1;
					count_out_bubble++;
				}
				ColorModel.id[n] = ColorModel.Mask->id[n];
			}
		}
	}
	printf("Count in %i, out %i\n",count_in_bubble,count_out_bubble);
	// initialize the phase indicator field
}

int main(int argc, char **argv)
{
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	// parallel domain size (# of sub-domains)
	int nprocx,nprocy,nprocz;
	int iproc,jproc,kproc;
	int sendtag,recvtag;
	//*****************************************
	// MPI ranks for all 18 neighbors
	//**********************************
	int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
	int rank_xy,rank_XY,rank_xY,rank_Xy;
	int rank_xz,rank_XZ,rank_xZ,rank_Xz;
	int rank_yz,rank_YZ,rank_yZ,rank_Yz;
	//**********************************
	MPI_Request req1[18],req2[18];
	MPI_Status stat1[18],stat2[18];

	if (rank == 0){
		printf("********************************************************\n");
		printf("Running Unit Test for D3Q7 Mass Conservation	\n");
		printf("********************************************************\n");
		if ( argc < 2 ) {
			std::cerr << "Invalid number of arguments, no input file specified\n";
			return -1;
		}
	}
	{		
	auto filename = argv[1];
	ScaLBL_ColorModel CM(rank,nprocs,comm);
	CM.ReadParams(filename);
	CM.SetDomain();    
	
	int i,j,k,n;
	int Nx,Ny,Nz,N,Np;
	Nx = CM.Nx;
	Ny = CM.Ny;
	Nz = CM.Nz;
	Np = CM.Np;
	N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);
	double *DenOriginal, *DenFinal;
	DenOriginal = new double [2*Np];
	DenFinal = new double [2*Np];
	
	//CM.ReadInput(); 
	double radius=0.4*double(Nx);
	InitializeBubble(CM,radius);
 	CM.Create();       // creating the model will create data structure to match the pore structure and allocate variables
	CM.Initialize();   // initializing the model will set initial conditions for variables
	//CM.Run();	       
	CM.WriteDebug();

	// Run the odd timestep 
 	ScaLBL_CopyToHost(DenOriginal,CM.Den,2*Np*sizeof(double));

	CM.ScaLBL_Comm->BiSendD3Q7AA(CM.Aq,CM.Bq); //READ FROM NORMAL
	ScaLBL_D3Q7_AAodd_PhaseField(CM.NeighborList, CM.dvcMap, CM.Aq, CM.Bq, CM.Den, CM.Phi, CM.ScaLBL_Comm->FirstInterior(), CM.ScaLBL_Comm->LastInterior(), CM.Np);
	CM.ScaLBL_Comm->BiRecvD3Q7AA(CM.Aq,CM.Bq); //WRITE INTO OPPOSITE
	ScaLBL_DeviceBarrier();
	ScaLBL_D3Q7_AAodd_PhaseField(CM.NeighborList, CM.dvcMap, CM.Aq, CM.Bq, CM.Den, CM.Phi, 0, CM.ScaLBL_Comm->LastExterior(), CM.Np);

	FILE *OUTFILE;
	OUTFILE = fopen("phi-out.raw","wb");
	fwrite(CM.Phi,8,N,OUTFILE);
	fclose(OUTFILE);	

	// Compare and make sure mass is conserved at every lattice site
	double *Error;
	Error = new double [N];
	bool CleanCheck = true;
	double original,final;
	ScaLBL_CopyToHost(DenFinal,CM.Den,2*Np*sizeof(double));
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;

				int idx = CM.Map(i,j,k);
				printf("idx=%i\n",idx);
				if (idx < Np && !(idx<0)){
					if (fabs(DenFinal[idx] - DenOriginal[idx]) > 1e-15){
						final = DenFinal[idx];
						original = DenOriginal[idx];
						//if (CM.Dm->id[n] == 0) printf("Solid phase! \n");
						//if (CM.Dm->id[n] == 1) printf("Wetting phase! \n");
						//if (CM.Dm->id[n] == 2) printf("Non-wetting phase! \n");							
						printf("Mass not conserved: WP density, site=%i,%i,%i, original = %f, final = %f \n",i,j,k,original,final);
						CleanCheck=false;
					}
					if (fabs(DenFinal[Np+idx] - DenOriginal[Np+idx]) > 1e-15){
						//if (CM.Dm->id[n] == 0) printf("Solid phase! \n");
						//if (CM.Dm->id[n] == 1) printf("Wetting phase! \n");
						//if (CM.Dm->id[n] == 2) printf("Non-wetting phase! \n");
						final = DenFinal[Np+idx];
						original = DenOriginal[Np+idx];
						printf("Mass not conserved: NWP density, site=%i,%i,%i, original = %f, final = %f \n",i,j,k,original,final);
						CleanCheck=false;
					}
				}
			}
		}
	}
	/*	if (rank==0) printf("Checking that the correct velocity is retained \n");
	// Swap convention is observed -- velocity is negative
	double *Aeven,*Aodd,*Beven,*Bodd;
	Aeven = new double[4*N];
	Aodd = new double[3*N];
	Beven = new double[4*N];
	Bodd = new double[3*N];
	ScaLBL_CopyToHost(Aeven,A_even,4*dist_mem_size);
	ScaLBL_CopyToHost(Aodd,A_odd,3*dist_mem_size);
	ScaLBL_CopyToHost(Beven,B_even,4*dist_mem_size);
	ScaLBL_CopyToHost(Bodd,B_odd,3*dist_mem_size);
	double rho,ux,uy,uz;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				rho = ux = uy = uz = 0.0;
				n = k*Nx*Ny + j*Nx + i;
				if (id[n] != 0){
					rho = Aeven[n]+Aeven[N+n]+Aeven[2*N+n]+Aeven[3*N+n]+Aodd[n]+Aodd[N+n]+Aodd[2*N+n]+
							Beven[n]+Beven[N+n]+Beven[2*N+n]+Beven[3*N+n]+Bodd[n]+Bodd[N+n]+Bodd[2*N+n];
					ux = Aeven[N+n] - Aodd[n] + Beven[N+n] - Bodd[n];
					uy = Aeven[2*N+n] - Aodd[N+n] + Beven[2*N+n] - Bodd[N+n];
					uz = Aeven[3*N+n] - Aodd[2*N+n] + Beven[3*N+n] - Bodd[2*N+n];
					if ( fabs(0.1+ux / rho) > 1e-13 ){
							if (id[n] == 1) printf("Wetting phase! \n");
							if (id[n] == 2) printf("Non-wetting phase! \n");
							final = ux/rho;
							printf("Momentum (x) not conserved, site=%i,%i,%i, final = %f \n",i,j,k,final);
							CleanCheck=false;
					}
					if ( fabs(0.1+uy / rho) > 1e-13 ){
							if (id[n] == 1) printf("Wetting phase! \n");
							if (id[n] == 2) printf("Non-wetting phase! \n");
							final = uy/rho;
							printf("Momentum (y) not conserved, site=%i,%i,%i, final = %f \n",i,j,k,final);
							CleanCheck=false;
					}
					if ( fabs(0.1+uz / rho) > 1e-13 ){
							if (id[n] == 1) printf("Wetting phase! \n");
							if (id[n] == 2) printf("Non-wetting phase! \n");
							final = uz/rho;
							printf("Momentum (z) not conserved, site=%i,%i,%i, final = %f \n",i,j,k,final);
							CleanCheck=false;
					}
				}
			}
		}
	}
*/
	if (CleanCheck){
		if (rank==0) printf("Test passed: mass conservation for D3Q7 \n");
	}
	else {
		if (rank==0) printf("Test failed!: mass conservation for D3Q7 \n");

	}
}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
