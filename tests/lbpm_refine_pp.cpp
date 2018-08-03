/*
 * Pre-processor to refine signed distance mesh
 * this is a good way to increase the resolution 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "common/Array.h"
#include "common/Communication.h"
#include "common/Domain.h"
#include "analysis/pmmc.h"

int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);

	{
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
		double Lx, Ly, Lz;
		int i,j,k,n;
		int BC=0;

		string filename;
		if (argc > 1){
			filename=argv[1];
		}
		else ERROR("No input database provided\n");
		// read the input database 
		auto db = std::make_shared<Database>( filename );
		auto domain_db = db->getDatabase( "Domain" );

		// Read domain parameters
		auto L = domain_db->getVector<double>( "L" );
		auto size = domain_db->getVector<int>( "n" );
		auto nproc = domain_db->getVector<int>( "nproc" );
		auto ReadValues = domain_db->getVector<char>( "ReadValues" );
		auto WriteValues = domain_db->getVector<char>( "WriteValues" );
		
		nx = size[0];
		ny = size[1];
		nz = size[2];
		nprocx = nproc[0];
		nprocy = nproc[1];
		nprocz = nproc[2];

		// Check that the number of processors >= the number of ranks
		if ( rank==0 ) {
			printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
			printf("Number of MPI ranks used: %i \n", nprocs);
			printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
		}
		if ( nprocs < nprocx*nprocy*nprocz ){
			ERROR("Insufficient number of processors");
		}

		char LocalRankFilename[40];

		int rnx,rny,rnz;
		rnx=2*nx;
		rny=2*ny;
		rnz=2*nz;

		if (rank==0) printf("Refining mesh to %i x %i x %i \n",rnx,rny,rnz);

		int BoundaryCondition=0;
		Domain Dm(rnx,rny,rnz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);

		// Communication the halos
		const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
		fillHalo<double> fillData(comm,rank_info,{rnx,rny,rnz},{1,1,1},0,1);

		nx+=2; ny+=2; nz+=2;
		rnx+=2; rny+=2; rnz+=2;
		int N = nx*ny*nz;

		// Define communication sub-domain -- everywhere
		for (int k=0; k<rnz; k++){
			for (int j=0; j<rny; j++){
				for (int i=0; i<rnx; i++){
					n = k*rnx*rny+j*rnx+i;
					Dm.id[n] = 1;
				}
			}
		}
		Dm.CommInit();

		DoubleArray SignDist(nx,ny,nz);
		// Read the signed distance from file
		sprintf(LocalRankFilename,"SignDist.%05i",rank);
		FILE *DIST = fopen(LocalRankFilename,"rb");
		size_t ReadSignDist;
		ReadSignDist=fread(SignDist.data(),8,N,DIST);
		if (ReadSignDist != size_t(N)) printf("lbpm_refine_pp: Error reading signed distance function (rank=%i)\n",rank);
		fclose(DIST);
		
		char *Labels;
		Labels = new char[N];
		sprintf(LocalRankFilename,"ID.%05i",rank);
		FILE *LABELS = fopen(LocalRankFilename,"rb");
		size_t ReadLabels;
		ReadLabels=fread(Labels,1,N,LABELS);
		if (ReadLabels != size_t(N)) printf("lbpm_refine_pp: Error reading ID  (rank=%i)\n",rank);
		fclose(LABELS);

		if ( rank==0 )   printf("Set up Domain, read input distance \n");

		DoubleArray RefinedSignDist(rnx,rny,rnz);
		TriLinPoly LocalApprox;
		Point pt;

		// double the distance value (since it is evaluated in pixels)
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					SignDist(i,j,k) =  2.f*SignDist(i,j,k);
				}
			}
		}

		int ri,rj,rk,rn; //refined mesh indices
		//char *RefineLabel;
		//RefineLabel = new char [rnx*rny*rnz];
		Array <char> RefineLabel(rnx,rny,rnz);
		for (rk=1; rk<rnz-1; rk++){
			for (rj=1; rj<rny-1; rj++){
				for (ri=1; ri<rnx-1; ri++){
					n = rk*rnx*rny+rj*rnx+ri;
					// starting node for each processor matches exactly
					i = (ri-1)/2+1;
					j = (rj-1)/2+1;
					k = (rk-1)/2+1;

					//printf("(%i,%i,%i -> %i,%i,%i) \n",ri,rj,rk,i,j,k);
					// Assign local tri-linear polynomial
					LocalApprox.assign(SignDist,i,j,k);
					pt.x=0.5*(ri-1)+1.f;
					pt.y=0.5*(rj-1)+1.f;
					pt.z=0.5*(rk-1)+1.f;
					RefinedSignDist(ri,rj,rk) = LocalApprox.eval(pt);
					RefineLabel(ri,rj,rk) = Labels[k*nx*ny+j*nx+i]; 
				}
			}
		}
		fillData.fill(RefinedSignDist);
		//	sprintf(LocalRankFilename,"ID.%05i",rank);
		//FILE *ID = fopen(LocalRankFilename,"wb");
		//fwrite(id,1,N,ID);
		//fclose(ID);
/*
		sprintf(LocalRankFilename,"RefineDist.%05i",rank);
		FILE *REFINEDIST = fopen(LocalRankFilename,"wb");
		fwrite(RefinedSignDist.data(),8,rnx*rny*rnz,REFINEDIST);
		fclose(REFINEDIST);
*/
		if ( rank==0 )   printf("Write output \n");

		DoubleArray BlockDist(nx,ny,nz);
		FILE *WRITEID, *REFINEDIST;
		char * id;
		id = new char [N];
		int writerank;

		// Write output blocks with the same sub-domain size as origina
		// refinement increases the size of the process grid
 		writerank = 8*Dm.kproc()*nprocx*nprocy + 4*Dm.jproc()*nprocx + 2*Dm.iproc();
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					BlockDist(i,j,k) = RefinedSignDist(i,j,k);
					if (BlockDist(i,j,k) > 0) 	id[k*nx*ny + j*nx + i]=2;
					else 						id[k*nx*ny + j*nx + i] = RefineLabel(i,j,k);
				}
			}
		}
		sprintf(LocalRankFilename,"RefineDist.%05i",writerank);
		REFINEDIST = fopen(LocalRankFilename,"wb");
		fwrite(BlockDist.data(),8,nx*ny*nz,REFINEDIST);
		fclose(REFINEDIST);

/*		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					if (BlockDist(i,j,k) > 0.f)
						id[k*nx*ny + j*nx + i]=2;
					else
						id[k*nx*ny + j*nx + i]= 0;
				}
			}
		}
		*/
		sprintf(LocalRankFilename,"RefineID.%05i",writerank);
		WRITEID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,nx*ny*nz,WRITEID);
		fclose(WRITEID);

		writerank = 8*Dm.kproc()*nprocx*nprocy + 4*Dm.jproc()*nprocx + 2*Dm.iproc()+1;
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					BlockDist(i,j,k) = RefinedSignDist(i+nx-2,j,k);
					if (BlockDist(i,j,k) > 0) 	id[k*nx*ny + j*nx + i]=2;
					else 						id[k*nx*ny + j*nx + i] = RefineLabel(i+nx-2,j,k);
				}
			}
		}
		sprintf(LocalRankFilename,"RefineDist.%05i",writerank);
		REFINEDIST = fopen(LocalRankFilename,"wb");
		fwrite(BlockDist.data(),8,nx*ny*nz,REFINEDIST);
		fclose(REFINEDIST);

/*		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					if (BlockDist(i,j,k) > 0.f)
						id[k*nx*ny + j*nx + i]=2;
					else
						id[k*nx*ny + j*nx + i]=0;
				}
			}
		}
		*/
		sprintf(LocalRankFilename,"RefineID.%05i",writerank);
		WRITEID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,nx*ny*nz,WRITEID);
		fclose(WRITEID);


		writerank = (2*Dm.kproc())*4*nprocx*nprocy + (2*Dm.jproc()+1)*2*nprocx + 2*Dm.iproc()+1;
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					BlockDist(i,j,k) = RefinedSignDist(i+nx-2,j+ny-2,k);
					if (BlockDist(i,j,k) > 0) 	id[k*nx*ny + j*nx + i]=2;
					else 						id[k*nx*ny + j*nx + i] = RefineLabel(i+nx-2,j+ny-2,k);
				}
			}
		}
		sprintf(LocalRankFilename,"RefineDist.%05i",writerank);
		REFINEDIST = fopen(LocalRankFilename,"wb");
		fwrite(BlockDist.data(),8,nx*ny*nz,REFINEDIST);
		fclose(REFINEDIST);

/*		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					if (BlockDist(i,j,k) > 0.f)
						id[k*nx*ny + j*nx + i]=2;
					else
						id[k*nx*ny + j*nx + i]=0;
				}
			}
		}
		*/
		sprintf(LocalRankFilename,"RefineID.%05i",writerank);
		WRITEID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,nx*ny*nz,WRITEID);
		fclose(WRITEID);

		writerank = (2*Dm.kproc())*4*nprocx*nprocy + (2*Dm.jproc()+1)*2*nprocx + 2*Dm.iproc();
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					BlockDist(i,j,k) = RefinedSignDist(i,j+ny-2,k);
					if (BlockDist(i,j,k) > 0) 	id[k*nx*ny + j*nx + i]=2;
					else 						id[k*nx*ny + j*nx + i] = RefineLabel(i,j+ny-2,k);
				}
			}
		}
		sprintf(LocalRankFilename,"RefineDist.%05i",writerank);
		REFINEDIST = fopen(LocalRankFilename,"wb");
		fwrite(BlockDist.data(),8,nx*ny*nz,REFINEDIST);
		fclose(REFINEDIST);
/*
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					if (BlockDist(i,j,k) > 0.f)
						id[k*nx*ny + j*nx + i]=2;
					else
						id[k*nx*ny + j*nx + i]=0;
				}
			}
		}
		*/
		sprintf(LocalRankFilename,"RefineID.%05i",writerank);
		WRITEID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,nx*ny*nz,WRITEID);
		fclose(WRITEID);

		writerank = (2*Dm.kproc()+1)*4*nprocx*nprocy + (2*Dm.jproc())*2*nprocx + 2*Dm.iproc();
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					BlockDist(i,j,k) = RefinedSignDist(i,j,k+nz-2);
					if (BlockDist(i,j,k) > 0) 	id[k*nx*ny + j*nx + i]=2;
					else 						id[k*nx*ny + j*nx + i] = RefineLabel(i,j,k+nz-2);
				}
			}
		}
		sprintf(LocalRankFilename,"RefineDist.%05i",writerank);
		REFINEDIST = fopen(LocalRankFilename,"wb");
		fwrite(BlockDist.data(),8,nx*ny*nz,REFINEDIST);
		fclose(REFINEDIST);
/*
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					if (BlockDist(i,j,k) > 0.f)
						id[k*nx*ny + j*nx + i]=2;
					else
						id[k*nx*ny + j*nx + i]=0;
				}
			}
		}
		*/
		sprintf(LocalRankFilename,"RefineID.%05i",writerank);
		WRITEID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,nx*ny*nz,WRITEID);
		fclose(WRITEID);

		writerank = (2*Dm.kproc()+1)*4*nprocx*nprocy + (2*Dm.jproc())*2*nprocx + 2*Dm.iproc()+1;
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					BlockDist(i,j,k) = RefinedSignDist(i+nx-2,j,k+nz-2);
					if (BlockDist(i,j,k) > 0) 	id[k*nx*ny + j*nx + i]=2;
					else 						id[k*nx*ny + j*nx + i] = RefineLabel(i+nx-2,j,k+nz-2);
				}
			}
		}
		sprintf(LocalRankFilename,"RefineDist.%05i",writerank);
		REFINEDIST = fopen(LocalRankFilename,"wb");
		fwrite(BlockDist.data(),8,nx*ny*nz,REFINEDIST);
		fclose(REFINEDIST);

/*		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					if (BlockDist(i,j,k) > 0.f)
						id[k*nx*ny + j*nx + i]=2;
					else
						id[k*nx*ny + j*nx + i]=0;
				}
			}
		}
		*/
		sprintf(LocalRankFilename,"RefineID.%05i",writerank);
		WRITEID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,nx*ny*nz,WRITEID);
		fclose(WRITEID);

		writerank = (2*Dm.kproc()+1)*4*nprocx*nprocy + (2*Dm.jproc()+1)*2*nprocx + 2*Dm.iproc();
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					BlockDist(i,j,k) = RefinedSignDist(i,j+ny-2,k+nz-2);
					if (BlockDist(i,j,k) > 0) 	id[k*nx*ny + j*nx + i]=2;
					else 						id[k*nx*ny + j*nx + i] = RefineLabel(i,j+ny-2,k+nz-2);
				}
			}
		}
		sprintf(LocalRankFilename,"RefineDist.%05i",writerank);
		REFINEDIST = fopen(LocalRankFilename,"wb");
		fwrite(BlockDist.data(),8,nx*ny*nz,REFINEDIST);
		fclose(REFINEDIST);
/*
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					if (BlockDist(i,j,k) > 0.f)
						id[k*nx*ny + j*nx + i]=2;
					else
						id[k*nx*ny + j*nx + i]=0;
				}
			}
		}
		*/
		sprintf(LocalRankFilename,"RefineID.%05i",writerank);
		WRITEID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,nx*ny*nz,WRITEID);
		fclose(WRITEID);

		writerank = (2*Dm.kproc()+1)*4*nprocx*nprocy + (2*Dm.jproc()+1)*2*nprocx + 2*Dm.iproc()+1;
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					BlockDist(i,j,k) = RefinedSignDist(i+nx-2,j+ny-2,k+nz-2);
					if (BlockDist(i,j,k) > 0) 	id[k*nx*ny + j*nx + i]=2;
					else 						id[k*nx*ny + j*nx + i] = RefineLabel(i+nx-2,j+ny-2,k+nz-2);
				}
			}
		}

		sprintf(LocalRankFilename,"RefineDist.%05i",writerank);
		REFINEDIST = fopen(LocalRankFilename,"wb");
		fwrite(BlockDist.data(),8,nx*ny*nz,REFINEDIST);
		fclose(REFINEDIST);
		
/*		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					if (BlockDist(i,j,k) > 0.f)
						id[k*nx*ny + j*nx + i]=2;
					else
						id[k*nx*ny + j*nx + i]=0;
				}
			}
		}
		*/
		sprintf(LocalRankFilename,"RefineID.%05i",writerank);
		WRITEID = fopen(LocalRankFilename,"wb");
		fwrite(id,1,nx*ny*nz,WRITEID);
		fclose(WRITEID);


	}
	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;
}
