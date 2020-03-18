/*
 * Pre-processor to generate signed distance function from segmented data
 * segmented data should be stored in a raw binary file as 1-byte integer (type char)
 * will output distance functions for phases
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "common/Array.h"
#include "common/Domain.h"



int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);

	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	{


		bool MULTINPUT=false;

		int NWP,SOLID,rank_offset;
		int NLABELS=atoi(argv[1]);
		//SOLID=atoi(argv[1]);
		//NWP=atoi(argv[2]);
		if (rank==0){
			//printf("Solid Label: %i \n",SOLID);
			//printf("NWP Label: %i \n",NWP);
		}
		if (argc > 2){
			MULTINPUT=true;
			rank_offset = atoi(argv[3]);
		}
		else{
			rank_offset=0;
		}

		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
		double Lx, Ly, Lz;
		int Nx,Ny,Nz;
		int i,j,k,n;
		int BC=0;
		char Filename[40];
		int xStart,yStart,zStart;
		//  char fluidValue,solidValue;

		std::vector<char> solidValues;
		std::vector<char> nwpValues;
		std::string line;

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

			ifstream image("Segmented.in");
			image >> Filename; 	// Name of data file containing segmented data
			image >> Nx;   		// size of the binary file
			image >> Ny;
			image >> Nz;
			image >> xStart;	// offset for the starting voxel
			image >> yStart;
			image >> zStart;

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
		//.................................................
		MPI_Bcast(&Nx,1,MPI_INT,0,comm);
		MPI_Bcast(&Ny,1,MPI_INT,0,comm);
		MPI_Bcast(&Nz,1,MPI_INT,0,comm);
		MPI_Bcast(&xStart,1,MPI_INT,0,comm);
		MPI_Bcast(&yStart,1,MPI_INT,0,comm);
		MPI_Bcast(&zStart,1,MPI_INT,0,comm);
		//.................................................
		MPI_Barrier(comm);

		// Check that the number of processors >= the number of ranks
		if ( rank==0 ) {
			printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
			printf("Number of MPI ranks used: %i \n", nprocs);
			printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
		}
		if ( nprocs < nprocx*nprocy*nprocz ){
			ERROR("Insufficient number of processors");
		}
		char *SegData = NULL;
		// Rank=0 reads the entire segmented data and distributes to worker processes
		if (rank==0){
			printf("Dimensions of segmented image: %i x %i x %i \n",Nx,Ny,Nz);
			SegData = new char[Nx*Ny*Nz];
			FILE *SEGDAT = fopen(Filename,"rb");
			if (SEGDAT==NULL) ERROR("Error reading segmented data");
			size_t ReadSeg;
			ReadSeg=fread(SegData,1,Nx*Ny*Nz,SEGDAT);
			if (ReadSeg != size_t(Nx*Ny*Nz)) printf("lbpm_segmented_decomp: Error reading segmented data (rank=%i)\n",rank);
			fclose(SEGDAT);
			printf("Read segmented data from %s \n",Filename);
		}
		MPI_Barrier(comm);

		// Get the rank info
		int N = (nx+2)*(ny+2)*(nz+2);
		Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
		for (k=0;k<nz+2;k++){
			for (j=0;j<ny+2;j++){
				for (i=0;i<nx+2;i++){
					n = k*(nx+2)*(ny+2)+j*(nx+2)+i;
					Dm.id[n] = 1;
				}
			}
		}
		Dm.CommInit();

		// number of sites to use for periodic boundary condition transition zone
		int z_transition_size = (nprocz*nz - (Nz - zStart))/2;
		if (z_transition_size < 0) z_transition_size=0;

		// Set up the sub-domains
		if (rank==0){
			printf("Distributing subdomains across %i processors \n",nprocs);
			printf("Process grid: %i x %i x %i \n",Dm.nprocx(),Dm.nprocy(),Dm.nprocz());
			printf("Subdomain size: %i \n",N);
			printf("Size of transition region: %i \n", z_transition_size);
			char *tmp;
			tmp = new char[N];
			for (int kp=0; kp<nprocz; kp++){
				for (int jp=0; jp<nprocy; jp++){
					for (int ip=0; ip<nprocx; ip++){
						// rank of the process that gets this subdomain
						int rnk = kp*Dm.nprocx()*Dm.nprocy() + jp*Dm.nprocx() + ip;
						// Pack and send the subdomain for rnk
						for (k=0;k<nz+2;k++){
							for (j=0;j<ny+2;j++){
								for (i=0;i<nx+2;i++){
									int x = xStart + ip*nx + i-1;
									int y = yStart + jp*ny + j-1;
									//		int z = zStart + kp*nz + k-1;
									int z = zStart + kp*nz + k-1 - z_transition_size;
									if (x<xStart) 	x=xStart;
									if (!(x<Nx))	x=Nx-1;
									if (y<yStart) 	y=yStart;
									if (!(y<Ny))	y=Ny-1;
									if (z<zStart) 	z=zStart;
									if (!(z<Nz))	z=Nz-1;
									int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
									int nglobal = z*Nx*Ny+y*Nx+x;
									tmp[nlocal] = SegData[nglobal];
								}
							}
						}
						if (rnk==0){
							for (k=0;k<nz+2;k++){
								for (j=0;j<ny+2;j++){
									for (i=0;i<nx+2;i++){
										int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
										Dm.id[nlocal] = tmp[nlocal];
									}
								}
							}
						}
						else{
							printf("Sending data to process %i \n", rnk);
							MPI_Send(tmp,N,MPI_CHAR,rnk,15,comm);
						}
					}
				}
			}
		}
		else{
			// Recieve the subdomain from rank = 0
			printf("Ready to recieve data %i at process %i \n", N,rank);
			MPI_Recv(Dm.id,N,MPI_CHAR,0,15,comm,MPI_STATUS_IGNORE);
		}
		MPI_Barrier(comm);

		nx+=2; ny+=2; nz+=2;
		N=nx*ny*nz;
		if (rank==0) printf("All sub-domains recieved \n");

		//  Assign New Labels
		int *LabelList;
		LabelList=new int[2*NLABELS];
		if (rank==0){
			printf("Assigning new lablels \n");
			if (rank==0){
				printf("Component labels:\n");
				ifstream iFILE("ComponentLabels.csv");
				if (iFILE.good()){
					int oldlabel, newlabel;
					int label=0;
					while (!iFILE.eof()){
						iFILE>>oldlabel;
						iFILE>>newlabel;
						LabelList[2*label] = (oldlabel);
						LabelList[2*label+1] = (newlabel);
						label++;
					}
				}
				else{
					printf("Using default labels: Solid (0 --> -1), NWP (1 --> 1), WP (2 --> 2)\n");
					// Set default values
					NLABELS=3;
					for (int label=0; label<NLABELS; label++){
						LabelList[2*label] = (label);
						LabelList[2*label+1] = (label);
					}
				}
			}
			for (int label=0; label<NLABELS; label++){
				int oldlabel=LabelList[2*label];
				int newlabel=LabelList[2*label+1];
				printf("Original label=%i, New label=%i \n",oldlabel,newlabel);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(LabelList,2*NLABELS,MPI_INT,0,MPI_COMM_WORLD);
		
		char *newIDs;
		newIDs= new char [nx*ny*nz];
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n = k*nx*ny+j*nx+i;
					for (int label=0; label<NLABELS; label++){
						int oldlabel=LabelList[2*label];
						int newlabel=LabelList[2*label+1];
						if (Dm.id[n]==char(oldlabel))  newIDs[n] = char(newlabel);
					}
				}
			}
		}
		if (rank==0) printf("Domain set \n");

		int count = 0;
		int total = 0;
		int countGlobal = 0;
		int totalGlobal = 0;
		for (k=1;k<nz-1;k++){
			for (j=1;j<ny-1;j++){
				for (i=1;i<nx-1;i++){
					n=k*nx*ny+j*nx+i;
					total++;
					Dm.id[n] = newIDs[n];
					if (Dm.id[n] == 0){
						count++;
					}
				}
			}
		}
		MPI_Allreduce(&count,&countGlobal,1,MPI_INT,MPI_SUM,comm);
		MPI_Allreduce(&total,&totalGlobal,1,MPI_INT,MPI_SUM,comm);


		float porosity = float(totalGlobal-countGlobal)/totalGlobal;
		if (rank==0) printf("Porosity=%f\n",porosity);

		if (rank==0){
			int xstart = xStart;    // Is this correct?
			int ystart = yStart;
			int zstart = zStart;
			//totalGlobal=(Nx-xstart)*(Ny-ystart)*(Nz-zstart);
			countGlobal = 0;
			for (k=zstart; k<zstart+nprocz*(nz-2); k++){
				for (j=ystart; j<ystart+nprocy*(ny-2); j++){
					for (i=xstart; i<xstart+nprocx*(nx-2); i++){

						n=k*Nx*Ny+j*Nx+i;
						if (n < Nx*Ny*Nz){
							if (SegData[n] == char(SOLID)){
								countGlobal++;
							}
						}
					}
				}
			}
			float porosity = float(totalGlobal-countGlobal)/totalGlobal;
			printf("Original Porosity=%f\n",porosity);
		}

		count = 0;
		total = 0;
		countGlobal = 0;
		totalGlobal = 0;
		for (k=1;k<nz-1;k++){
			for (j=1;j<ny-1;j++){
				for (i=1;i<nx-1;i++){
					n=k*nx*ny+j*nx+i;
					if (Dm.id[n] != 0)      total++;
					if (Dm.id[n] == 2)	count++;
				}
			}
		}
		MPI_Allreduce(&count,&countGlobal,1,MPI_INT,MPI_SUM,comm);
		MPI_Allreduce(&total,&totalGlobal,1,MPI_INT,MPI_SUM,comm);
		float saturation = float(countGlobal)/totalGlobal;
		if (rank==0) printf("wetting phase saturation=%f\n",saturation);


		char LocalRankFilename[40];

		sprintf(LocalRankFilename,"ID.%05i",rank+rank_offset);
		FILE *ID = fopen(LocalRankFilename,"wb");
		fwrite(Dm.id,1,N,ID);
		//    fwrite(Distance.get(),8,Distance.length(),ID);
		fclose(ID);

		if (!MULTINPUT){

			if (rank==0) printf("Writing symmetric domain reflection\n");
			MPI_Barrier(comm);
			int symrank,sympz;
			sympz = 2*nprocz - Dm.kproc() -1;
			symrank = sympz*nprocx*nprocy + Dm.jproc()*nprocx + Dm.iproc();

			//    DoubleArray SymDist(nx,ny,nz);
			char *symid;
			symid = new char [N];
			for (k=0;k<nz;k++){
				for (j=0;j<ny;j++){
					for (i=0;i<nx;i++){
						n=k*nx*ny+j*nx+i;
						int nsym=(nz-k-1)*nx*ny+j*nx+i;
						symid[nsym] = Dm.id[n];
						//SymDist(i,j,nz-k-1)=Distance(i,j,k);
					}
				}
			}

			sprintf(LocalRankFilename,"ID.%05i",symrank);
			FILE *SYMID = fopen(LocalRankFilename,"wb");
			//    fwrite(SymDist.get(),8,SymDist.length(),SYMDIST);
			fwrite(symid,1,N,SYMID);
			fclose(SYMID);
		}
	}
	MPI_Barrier(comm);
	MPI_Finalize();
}
