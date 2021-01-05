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
#include "analysis/distance.h"

int main(int argc, char **argv)
{
	// Initialize MPI
    Utilities::startup( argc, argv );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();
	{
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		double Lx, Ly, Lz;
        Lx = Ly = Lz = 1.0;
		int i,j,k,n;

		string filename;
		if (argc > 1){
			filename=argv[1];
		}
		else ERROR("No input database provided\n");
		// read the input database 
		auto db = std::make_shared<Database>( filename );
		auto domain_db = db->getDatabase( "Domain" );

		// Read domain parameters
		auto size = domain_db->getVector<int>( "n" );
		auto nproc = domain_db->getVector<int>( "nproc" );
		auto ReadValues = domain_db->getVector<char>( "ReadValues" );
		auto WriteValues = domain_db->getVector<char>( "WriteValues" );
		
		int nx = size[0];
		int ny = size[1];
		int nz = size[2];
		int nprocx = nproc[0];
		int nprocy = nproc[1];
		int nprocz = nproc[2];
		int BoundaryCondition=0;

		// Check that the number of processors >= the number of ranks
		if ( rank==0 ) {
			printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
			printf("Number of MPI ranks used: %i \n", nprocs);
			printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
		}
		if ( nprocs < nprocx*nprocy*nprocz ){
			ERROR("Insufficient number of processors");
		}

		//Domain Mask(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);
		Domain Mask(domain_db,MPI_COMM_WORLD);
		if (domain_db->keyExists( "Filename" )){
			auto Filename = domain_db->getScalar<std::string>( "Filename" );
		        if (rank==0) printf("Reading domain from %s \n",Filename.c_str());
			Mask.Decomp(Filename);
			if (rank==0) printf("Complete. \n");
		}
		else{
			Mask.ReadIDs();
		}
		Mask.CommInit();

		char LocalRankFilename[40];
		int rnx=2*nx;
		int rny=2*ny;
		int rnz=2*nz;

		if (rank==0) printf("Refining mesh to %i x %i x %i \n",rnx,rny,rnz);

		Domain Dm(rnx,rny,rnz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);

		// Communication the halos
		const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
		fillHalo<double> fillData(comm,rank_info,{rnx,rny,rnz},{1,1,1},0,1);

		nx+=2; ny+=2; nz+=2;
		rnx+=2; rny+=2; rnz+=2;
		int N = nx*ny*nz;

		// Define communication sub-domain -- everywhere
		if (rank==0) printf("Initialize refined domain \n");
		for (int k=0; k<rnz; k++){
			for (int j=0; j<rny; j++){
				for (int i=0; i<rnx; i++){
					n = k*rnx*rny+j*rnx+i;
					Dm.id[n] = 1;
				}
			}
		}
		Dm.CommInit();
		
		// Generate the signed distance map
		// Initialize the domain and communication
		Array<char> Labels(nx,ny,nz);
		DoubleArray SignDist(nx,ny,nz);

		// Solve for the position of the solid phase
		for (int k=0;k<nz;k++){
			for (int j=0;j<ny;j++){
				for (int i=0;i<nx;i++){
					int n = k*nx*ny+j*nx+i;
					// Initialize the solid phase
					signed char label = Mask.id[n];
					if (label > 0)		Labels(i,j,k) = 1;
					else	     		Labels(i,j,k) = 0;
				}
			}
		}
		// Initialize the signed distance function
		for (int k=0;k<nz;k++){
			for (int j=0;j<ny;j++){
				for (int i=0;i<nx;i++){
					// Initialize distance to +/- 1
					SignDist(i,j,k) = 2.0*double(Labels(i,j,k))-1.0;
				}
			}
		}
	//	MeanFilter(Averages->SDs);
		if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
		CalcDist(SignDist,Labels,Mask);
		
		/*		// Read the signed distance from file
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
*/
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

		//char *RefineLabel;
		//RefineLabel = new char [rnx*rny*rnz];
		Array <char> RefineLabel(rnx,rny,rnz);
		for (int rk=1; rk<rnz-1; rk++){
			for (int rj=1; rj<rny-1; rj++){
				for (int ri=1; ri<rnx-1; ri++){
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
					RefineLabel(ri,rj,rk) = Labels(i,j,k);
					Dm.id[n] = Labels(i,j,k);
				}
			}
		}
		fillData.fill(RefinedSignDist);


		if (domain_db->keyExists( "Filename" )){
			auto Filename = domain_db->getScalar<std::string>( "Filename" );
			if ( rank==0 )   printf("Write output \n");
			Dm.AggregateLabels("id_2x.raw");
			Mask.AggregateLabels("id.raw");
			//FILE *WRITEID = fopen("refine.raw","wb");
			//fwrite(RefineLabel.data(),1,rnx*rny*rnz,WRITEID);
			//fclose(WRITEID);
		}
		else{
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

			sprintf(LocalRankFilename,"RefineID.%05i",writerank);
			WRITEID = fopen(LocalRankFilename,"wb");
			fwrite(id,1,nx*ny*nz,WRITEID);
			fclose(WRITEID);
		}
	}
        Utilities::shutdown();
	return 0;
}
