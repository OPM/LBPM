/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * Multi-relaxation time LBM Model
 */
#include "models/MRTModel.h"

ScaLBL_MRTModel::ScaLBL_MRTModel(int RANK, int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tau(0),
Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),mu(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{

}
ScaLBL_MRTModel::~ScaLBL_MRTModel(){

}

void ScaLBL_MRTModel::ReadParams(string filename){
	// read the input database 
	db = std::make_shared<Database>( filename );
	domain_db = db->getDatabase( "Domain" );
	mrt_db = db->getDatabase( "MRT" );

	// Color Model parameters
	timestepMax = mrt_db->getScalar<int>( "timestepMax" );
	tau = mrt_db->getScalar<double>( "tau" );
	Fx = mrt_db->getVector<double>( "F" )[0];
	Fy = mrt_db->getVector<double>( "F" )[1];
	Fz = mrt_db->getVector<double>( "F" )[2];
	Restart = mrt_db->getScalar<bool>( "Restart" );
	din = mrt_db->getScalar<double>( "din" );
	dout = mrt_db->getScalar<double>( "dout" );
	flux = mrt_db->getScalar<double>( "flux" );
	
	// Read domain parameters
	auto L = domain_db->getVector<double>( "L" );
	auto size = domain_db->getVector<int>( "n" );
	auto nproc = domain_db->getVector<int>( "nproc" );
	BoundaryCondition = domain_db->getScalar<int>( "BC" );
	Nx = size[0];
	Ny = size[1];
	Nz = size[2];
	Lx = L[0];
	Ly = L[1];
	Lz = L[2];
	nprocx = nproc[0];
	nprocy = nproc[1];
	nprocz = nproc[2];
	mu=(tau-0.5)/3.0;
}
void ScaLBL_MRTModel::SetDomain(){
	Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
	Mask  = std::shared_ptr<Domain>(new Domain(domain_db,comm));    // mask domain removes immobile phases
	Nx+=2; Ny+=2; Nz += 2;
	N = Nx*Ny*Nz;
	Distance.resize(Nx,Ny,Nz);
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;               // initialize this way
	//Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm) ); // TwoPhase analysis object
	MPI_Barrier(comm);
	Dm->CommInit();
	MPI_Barrier(comm);
}

void ScaLBL_MRTModel::ReadInput(){
    int rank=Dm->rank();
    size_t readID;
    //.......................................................................
    if (rank == 0)    printf("Read input media... \n");
    //.......................................................................
    Mask->ReadIDs();
    
    sprintf(LocalRankString,"%05d",Dm->rank());
    sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
    sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);

    // .......... READ THE INPUT FILE .......................................
    //...........................................................................
    if (rank == 0) cout << "Reading in signed distance function..." << endl;
    //.......................................................................
    sprintf(LocalRankString,"%05d",rank);
    sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
    ReadBinaryFile(LocalRankFilename, Distance.data(), N);
    MPI_Barrier(comm);
    if (rank == 0) cout << "Domain set." << endl;
}

void ScaLBL_MRTModel::Create(){
	/*
	 *  This function creates the variables needed to run a LBM 
	 */
	int rank=Mask->rank();
	//.........................................................
	// Initialize communication structures in averaging domain
	for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i];
	Mask->CommInit();
	Np=Mask->PoreCount();
	//...........................................................................
	if (rank==0)    printf ("Create ScaLBL_Communicator \n");
	// Create a communicator for the device (will use optimized layout)
	// ScaLBL_Communicator ScaLBL_Comm(Mask); // original
	ScaLBL_Comm  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));

	int Npad=(Np/16 + 2)*16;
	if (rank==0)    printf ("Set up memory efficient layout \n");
	Map.resize(Nx,Ny,Nz);       Map.fill(-2);
	auto neighborList= new int[18*Npad];
	Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Mask->id,Np);
	MPI_Barrier(comm);
	//...........................................................................
	//                MAIN  VARIABLES ALLOCATED HERE
	//...........................................................................
	// LBM variables
	if (rank==0)    printf ("Allocating distributions \n");
	//......................device distributions.................................
	int dist_mem_size = Np*sizeof(double);
	int neighborSize=18*(Np*sizeof(int));
	//...........................................................................
	ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
	ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);  
	ScaLBL_AllocateDeviceMemory((void **) &Pressure, sizeof(double)*Np);
	ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
	//...........................................................................
	// Update GPU data structures
	if (rank==0)    printf ("Setting up device map and neighbor list \n");
	// copy the neighbor list 
	ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
	MPI_Barrier(comm);
}        

void ScaLBL_MRTModel::Initialize(){
	/*
	 * This function initializes model
	 */
    if (rank==0)    printf ("Initializing distributions \n");
    ScaLBL_D3Q19_Init(fq, Np);
}

void ScaLBL_MRTModel::Run(){
	double rlx_setA=1.0/tau;
	double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
	//.......create and start timer............
	double starttime,stoptime,cputime;
	ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
	starttime = MPI_Wtime();
	if (rank==0) printf("Beginning AA timesteps...\n");
	if (rank==0) printf("********************************************************\n");
	timestep=0;
	while (timestep < timestepMax) {
		//************************************************************************/
		timestep++;
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FROM NORMAL
		ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq,  ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		timestep++;
		ScaLBL_Comm->SendD3Q19AA(fq); //READ FORM NORMAL
		ScaLBL_D3Q19_AAeven_MRT(fq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
		ScaLBL_Comm->RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
		ScaLBL_D3Q19_AAeven_MRT(fq, 0, ScaLBL_Comm->LastExterior(), Np, rlx_setA, rlx_setB, Fx, Fy, Fz);
		ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
		//************************************************************************/
	}
	//************************************************************************/
	stoptime = MPI_Wtime();
	//	cout << "CPU time: " << (stoptime - starttime) << " seconds" << endl;
	cputime = stoptime - starttime;
	//	cout << "Lattice update rate: "<< double(Nx*Ny*Nz*timestep)/cputime/1000000 <<  " MLUPS" << endl;
	double MLUPS = double(Np*timestep)/cputime/1000000;
//	if (rank==0) printf("********************************************************\n");
//	if (rank==0) printf("CPU time = %f \n", cputime);
//	if (rank==0) printf("Lattice update rate (per process)= %f MLUPS \n", MLUPS);
	MLUPS *= nprocs;
}

void ScaLBL_MRTModel::VelocityField(double *Vz){

	int SIZE=Np*sizeof(double);
	ScaLBL_D3Q19_Momentum(fq,Velocity, Np);
	ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
	ScaLBL_CopyToHost(&Vz[0],&Velocity[2*Np],SIZE);
	
}
