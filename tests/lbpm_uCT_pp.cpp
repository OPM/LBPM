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
#include "IO/MeshDatabase.h"
#include "IO/Mesh.h"
#include "IO/Writer.h"
#include "IO/netcdf.h"
#include "analysis/analysis.h"
#include "analysis/filters.h"
#include "analysis/uCT.h"

#include "ProfilerApp.h"


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
		if ( rank == 0 ){
			printf("At least one filename must be specified\n");
		}
		return 1;
	}
	std::string filename = std::string(argv[1]);
	if ( rank == 0 ){
		printf("Input data file: %s\n",filename.c_str());
	}

	auto db = std::make_shared<Database>( filename );
	auto domain_db = db->getDatabase( "Domain" );
	auto uct_db = db->getDatabase( "uCT" );
	auto analysis_db = db->getDatabase( "Analysis" );

	// Read domain values
	auto L = domain_db->getVector<double>( "L" );
	auto size = domain_db->getVector<int>( "n" );
	auto nproc = domain_db->getVector<int>( "nproc" );
	int BoundaryCondition = domain_db->getScalar<int>( "BC" );
	int nx = size[0];
	int ny = size[1];
	int nz = size[2];
	double Lx = L[0];
	double Ly = L[1];
	double Lz = L[2];
	int nprocx = nproc[0];
	int nprocy = nproc[1];
	int nprocz = nproc[2];

	auto InputFile=uct_db->getScalar<std::string>( "InputFile" );    
	auto rough_cutoff=uct_db->getScalar<float>( "rought_cutoff" );    
	auto lamda=uct_db->getScalar<float>( "lamda" );    
	auto nlm_sigsq=uct_db->getScalar<float>( "nlm_sigsq" );    
	auto nlm_depth=uct_db->getScalar<int>( "nlm_depth" );    

	//.......................................................................
	// Reading the domain information file
	//.......................................................................
	//    std::shared_ptr<Domain> Dm ();
	//for (int i=0; i<Dm->Nx*Dm->Ny*Dm->Nz; i++) Dm->id[i] = 1;
	//Dm->CommInit();

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
	int ratio[3] = {2,2,2};
	//std::vector<size_t> ratio = {4,4,4};
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
		Dm[i].reset( new Domain(domain_db, comm) );
		int N = (Nx[i]+2)*(Ny[i]+2)*(Nz[i]+2);
		for (int n=0; n<N; n++){
			Dm[i]->id[n] = 1;
		}
		Dm[i]->CommInit();
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
		fillDouble[i].reset(new fillHalo<double>(Dm[i]->Comm,Dm[i]->rank_info,{Nx[i],Ny[i],Nz[i]},{1,1,1},0,1) );
		fillFloat[i].reset(new fillHalo<float>(Dm[i]->Comm,Dm[i]->rank_info,{Nx[i],Ny[i],Nz[i]},{1,1,1},0,1) );
		fillChar[i].reset(new fillHalo<char>(Dm[i]->Comm,Dm[i]->rank_info,{Nx[i],Ny[i],Nz[i]},{1,1,1},0,1) );
	}

	// Read the subvolume of interest on each processor
	PROFILE_START("ReadVolume");
	int fid = netcdf::open(InputFile,netcdf::READ);
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
	filter_src( *Dm[0], LOCVOL[0] );

	// Set up the mask to be distance to cylinder (crop outside cylinder)
	float CylRad=900;
	for (int k=0;k<Nz[0]+2;k++) {
		for (int j=0;j<Ny[0]+2;j++) {
			for (int i=0;i<Nx[0]+2;i++) {
				int iproc = Dm[0]->iproc();
				int jproc = Dm[0]->jproc();
				int kproc = Dm[0]->kproc();

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


	MPI_Barrier(comm);
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
		if (rank==0){
			printf("Coarsen level %i \n",i);
			printf("   Nx=%i, Ny=%i, Nz=%i \n",int(tmp.size(0)),int(tmp.size(1)),int(tmp.size(2))  );
			printf("   filter_x=%i, filter_y=%i, filter_z=%i \n",int(filter.size(0)),int(filter.size(1)),int(filter.size(2))  );
			printf("   ratio= %i,%i,%i \n",int(ratio[0]),int(ratio[1]),int(ratio[2])  );
		}
		MPI_Barrier(comm);
	}
	PROFILE_STOP("CoarsenMesh");

	// Initialize the coarse level
	PROFILE_START("Solve full mesh");
	if (rank==0)
		printf("Initialize full mesh\n");
	solve( LOCVOL[0], Mean[0], ID[0], Dist[0], MultiScaleSmooth[0],
			NonLocalMean[0], *fillFloat[0], *Dm[0], nprocx, 
			rough_cutoff, lamda, nlm_sigsq, nlm_depth);
	PROFILE_STOP("Solve fill mesh");
	MPI_Barrier(comm);

/*	// Initialize the coarse level
	PROFILE_START("Solve coarse mesh");
	if (rank==0)
		printf("Initialize coarse mesh\n");
	solve( LOCVOL.back(), Mean.back(), ID.back(), Dist.back(), MultiScaleSmooth.back(),
			NonLocalMean.back(), *fillFloat.back(), *Dm.back(), nprocx );
	PROFILE_STOP("Solve coarse mesh");
	MPI_Barrier(comm);

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
	MPI_Barrier(comm);    

	// Perform a final filter
	PROFILE_START("Filtering final domains");
	if (rank==0)
		printf("Filtering final domains\n");
	Array<float> filter_Mean, filter_Dist1, filter_Dist2;
	filter_final( ID[0], Dist[0], *fillFloat[0], *Dm[0], filter_Mean, filter_Dist1, filter_Dist2 );
	PROFILE_STOP("Filtering final domains");
*/
	//removeDisconnected( ID[0], *Dm[0] );

	/*
        // Write the distance function to a netcdf file
    const char* netcdf_filename = "Distance.nc";
    {
        RankInfoStruct info( rank, nprocx, nprocy, nprocz );
        std::vector<int> dim = { Nx[0]*nprocx, Ny[0]*nprocy, Nz[0]*nprocz };
        int fid = netcdf::open( netcdf_filename, netcdf::CREATE, MPI_COMM_WORLD );
        auto dims =  netcdf::defDim( fid, {"X", "Y", "Z"}, dim );
        Array<float> data(Nx[0],Ny[0],Nz[0]);
        fillFloat[0]->copy( Dist[0], data );
        netcdf::write( fid, "Distance", dims, data, info );
        netcdf::close( fid );
    }
	 */

	// Write the results to visit
	if (rank==0) printf("Setting up visualization structure \n");
	//	std::vector<IO::MeshDataStruct> meshData(N_levels);
	std::vector<IO::MeshDataStruct> meshData(1);
	//    for (size_t i=0; i<Nx.size(); i++) {
	// Mesh
	meshData[0].meshName = "image";
	meshData[0].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(Dm[0]->rank_info,Nx[0],Ny[0],Nz[0],Lx,Ly,Lz) );
	// Source data
	std::shared_ptr<IO::Variable> OrigData( new IO::Variable() );
	OrigData->name = "Source Data";
	OrigData->type = IO::VariableType::VolumeVariable;
	OrigData->dim = 1;
	OrigData->data.resize(Nx[0],Ny[0],Nz[0]);
	meshData[0].vars.push_back(OrigData);
	fillDouble[0]->copy( LOCVOL[0], OrigData->data );
	// Non-Local Mean
	std::shared_ptr<IO::Variable> NonLocMean( new IO::Variable() );
	NonLocMean->name = "Non-Local Mean";
	NonLocMean->type = IO::VariableType::VolumeVariable;
	NonLocMean->dim = 1;
	NonLocMean->data.resize(Nx[0],Ny[0],Nz[0]);
	meshData[0].vars.push_back(NonLocMean);
	fillDouble[0]->copy( NonLocalMean[0], NonLocMean->data );
	std::shared_ptr<IO::Variable> SegData( new IO::Variable() );
	SegData->name = "Segmented Data";
	SegData->type = IO::VariableType::VolumeVariable;
	SegData->dim = 1;
	SegData->data.resize(Nx[0],Ny[0],Nz[0]);
	meshData[0].vars.push_back(SegData);
	fillDouble[0]->copy( ID[0], SegData->data );
	// Signed Distance
	std::shared_ptr<IO::Variable> DistData( new IO::Variable() );
	DistData->name = "Signed Distance";
	DistData->type = IO::VariableType::VolumeVariable;
	DistData->dim = 1;
	DistData->data.resize(Nx[0],Ny[0],Nz[0]);
	meshData[0].vars.push_back(DistData);
	fillDouble[0]->copy( Dist[0], DistData->data );
	// Smoothed Data
	std::shared_ptr<IO::Variable> SmoothData( new IO::Variable() );
	SmoothData->name = "Smoothed Data";
	SmoothData->type = IO::VariableType::VolumeVariable;
	SmoothData->dim = 1;
	SmoothData->data.resize(Nx[0],Ny[0],Nz[0]);
	meshData[0].vars.push_back(SmoothData);
	fillDouble[0]->copy( MultiScaleSmooth[0], SmoothData->data );
	
	/*// Segmented Data
	       std::shared_ptr<IO::Variable> SegData( new IO::Variable() );
	    SegData->name = "Segmented Data";
	    SegData->type = IO::VariableType::VolumeVariable;
	    SegData->dim = 1;
	    SegData->data.resize(Nx[i],Ny[i],Nz[i]);
	    meshData[i].vars.push_back(SegData);
        fillDouble[i]->copy( ID[i], SegData->data );
        // Signed Distance
        std::shared_ptr<IO::Variable> DistData( new IO::Variable() );
	    DistData->name = "Signed Distance";
	    DistData->type = IO::VariableType::VolumeVariable;
	    DistData->dim = 1;
	    DistData->data.resize(Nx[i],Ny[i],Nz[i]);
	    meshData[i].vars.push_back(DistData);
        fillDouble[i]->copy( Dist[i], DistData->data );
        // Smoothed Data
        std::shared_ptr<IO::Variable> SmoothData( new IO::Variable() );
	    SmoothData->name = "Smoothed Data";
	    SmoothData->type = IO::VariableType::VolumeVariable;
	    SmoothData->dim = 1;
	    SmoothData->data.resize(Nx[i],Ny[i],Nz[i]);
	    meshData[i].vars.push_back(SmoothData);
        fillDouble[i]->copy( MultiScaleSmooth[i], SmoothData->data );

    }
      #if 0
        std::shared_ptr<IO::Variable> filter_Mean_var( new IO::Variable() );
	    filter_Mean_var->name = "Mean";
	    filter_Mean_var->type = IO::VariableType::VolumeVariable;
	    filter_Mean_var->dim = 1;
	    filter_Mean_var->data.resize(Nx[0],Ny[0],Nz[0]);
	    meshData[0].vars.push_back(filter_Mean_var);
        fillDouble[0]->copy( filter_Mean, filter_Mean_var->data );
        std::shared_ptr<IO::Variable> filter_Dist1_var( new IO::Variable() );
	    filter_Dist1_var->name = "Dist1";
	    filter_Dist1_var->type = IO::VariableType::VolumeVariable;
	    filter_Dist1_var->dim = 1;
	    filter_Dist1_var->data.resize(Nx[0],Ny[0],Nz[0]);
	    meshData[0].vars.push_back(filter_Dist1_var);
        fillDouble[0]->copy( filter_Dist1, filter_Dist1_var->data );
        std::shared_ptr<IO::Variable> filter_Dist2_var( new IO::Variable() );
	    filter_Dist2_var->name = "Dist2";
	    filter_Dist2_var->type = IO::VariableType::VolumeVariable;
	    filter_Dist2_var->dim = 1;
	    filter_Dist2_var->data.resize(Nx[0],Ny[0],Nz[0]);
	    meshData[0].vars.push_back(filter_Dist2_var);
        fillDouble[0]->copy( filter_Dist2, filter_Dist2_var->data );
    #endif
	 */
	MPI_Barrier(comm);
	if (rank==0) printf("Writing output \n");
	// Write visulization data
	IO::writeData( 0, meshData, comm );
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

