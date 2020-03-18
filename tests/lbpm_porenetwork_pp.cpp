#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "analysis/TwoPhase.h"
#include "common/MPI_Helpers.h"

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
	{
	// parallel domain size (# of sub-domains)
	int nprocx,nprocy,nprocz;
	int iproc,jproc,kproc;
	int sendtag,recvtag;
	//*****************************************
	MPI_Request req1[18],req2[18];
	MPI_Status stat1[18],stat2[18];
	//**********************************

	int nsph,ncyl, BC;
	nsph = atoi(argv[1]);
	ncyl = atoi(argv[2]);
	BC = 0;
	if (rank == 0){
		printf("********************************************************\n");
		printf("Generate LBM input geometry from simple pore network");
		printf("********************************************************\n");
	}

	// Variables that specify the computational domain  
	string FILENAME;
	int Nx,Ny,Nz;		// local sub-domain size
	int nspheres;		// number of spheres in the packing
	double Lx,Ly,Lz;	// Domain length
	int i,j,k,n;

	// pmmc threshold values

	if (rank==0){
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		ifstream domain("Domain.in");
		domain >> nprocx;
		domain >> nprocy;
		domain >> nprocz;
		domain >> Nx;
		domain >> Ny;
		domain >> Nz;
		domain >> Lx;
		domain >> Ly;
		domain >> Lz;
		//.......................................................................
	}
	// **************************************************************
	// Broadcast simulation parameters from rank 0 to all other procs
	MPI_Barrier(comm);
	// Computational domain
	MPI_Bcast(&Nx,1,MPI_INT,0,comm);
	MPI_Bcast(&Ny,1,MPI_INT,0,comm);
	MPI_Bcast(&Nz,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
	MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
	//.................................................
	MPI_Barrier(comm);
	
	// **************************************************************
	if (nprocs != nprocx*nprocy*nprocz){
		printf("nprocx =  %i \n",nprocx);
		printf("nprocy =  %i \n",nprocy);
		printf("nprocz =  %i \n",nprocz);
		INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
	}

	if (rank==0){
		printf("********************************************************\n");
		printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
		printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
		printf("********************************************************\n");
	}

	// Initialized domain and averaging framework for Two-Phase Flow
	//	Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
	//Dm.CommInit();
	//TwoPhase Averages(Dm);
	
        std::shared_ptr<Domain> Dm (new Domain(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC));
	Dm->CommInit();
        std::shared_ptr<TwoPhase> Averages( new TwoPhase(Dm) );

	MPI_Barrier(comm);

	Nx += 2; Ny += 2; Nz += 2;

	int N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);

	//.......................................................................
	// Filenames used
	char LocalRankString[8];
	char LocalRankFilename[40];
	char LocalRestartFile[40];
	char tmpstr[10];
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
	//if (pBC) iVol_global = 1.0/(1.0*(Nx-2)*nprocx*(Ny-2)*nprocy*((Nz-2)*nprocz-6));
	double pore_vol;
	
	DoubleArray cylinders(7,ncyl); // (x, y, z, X, Y, Z, radius)
	DoubleArray spheres(4,nsph); // ( x, y, z, radius)

	// Read from text file
	printf("Reading spheres \n");
	ifstream SPHERES ( "spheres.csv" ); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
	string value;
	int index=0;
	while ( index < nsph )	{
		getline ( SPHERES, value, ' ' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		spheres(0,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		getline ( SPHERES, value, ' ' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		spheres(1,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		getline ( SPHERES, value, ' ' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		spheres(2,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		getline ( SPHERES, value ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		spheres(3,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		printf("cx=%f,cy=%f,cz=%f,r=%f\n",spheres(0,index),spheres(1,index),spheres(2,index),spheres(3,index));
		index++;
		//cout << string( value, 1, value.length()-2 ); // display value removing the first and the last character from it
	}	

	printf("Reading cylinders \n");
	ifstream CYLINDERS ( "cylinders.csv" ); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
	index=0;
	while ( index < ncyl )	{
		getline ( CYLINDERS, value, ' ' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		cylinders(0,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		getline ( CYLINDERS, value, ' ' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		cylinders(1,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		getline ( CYLINDERS, value, ' ' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		cylinders(2,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		getline ( CYLINDERS, value, ' ' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		cylinders(3,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		getline ( CYLINDERS, value, ' ' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		cylinders(4,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		getline ( CYLINDERS, value, ' ' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		cylinders(5,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		getline ( CYLINDERS, value ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		cylinders(6,index)=strtod(value.c_str(),NULL)*double(Nx-1)/Lx;
		printf("x=%f,y=%f,z=%f,r=%f\n",cylinders(0,index),cylinders(1,index),cylinders(2,index),cylinders(6,index));
		index++;
		//cout << string( value, 1, value.length()-2 ); // display value removing the first and the last character from it
	}	

	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				Averages->SDs(i,j,k) = -2.f*Nx;
			}
		}
	}
	sum=0;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;

				// Compute the distance to each cylinder
				for (int p=0; p<ncyl; p++){
				  double x = cylinders(0,p);
					double y = cylinders(1,p);
					double z = cylinders(2,p);
					double X = cylinders(3,p);
					double Y = cylinders(4,p);
					double Z = cylinders(5,p);
					double radius = cylinders(6,p);
					double length = sqrt(x*x+y*y+z*z);
					double alpha = (X - x)/length;
					double beta = (Y - y)/length;
					double gamma = (Z - z)/length;
					double xi = double(i) - x;
					double yj = double(j) - y;
					double zk = double(k) - z;
					// value of s along center line {x=alpha*s, y = beta*s, z=gamma*s} that is closest to i,j,k
					double s = (alpha*xi + beta*yj + gamma*zk)/(alpha*alpha + beta*beta + gamma*gamma);
					double distance=Averages->SDs(i,j,k);
					if (s > length){
						//distance = radius - sqrt((xi-X)*(xi-X) + (yj-Y)*(yj-Y) + (zk-Z)*(zk-Z));
					}
					else if (s<0){
						//distance = radius - sqrt((xi-x)*(xi-x) + (yj-y)*(yj-y) + (zk-z)*(zk-z));
					}
					else{
						double xs = alpha*s;
						double ys = beta*s;
						double zs = gamma*s;
						distance = radius - sqrt((xi-xs)*(xi-xs) + (yj-ys)*(yj-ys) + (zk-zs)*(zk-zs));
						//if (distance>0)printf("s=%f,alpha=%f,beta=%f,gamma=%f,distance=%f\n",s,alpha,beta,gamma,distance);
					}
					if (distance > Averages->SDs(i,j,k))		Averages->SDs(i,j,k) = distance;
				}
				
				// Compute the distance to each sphere
				for (int p=0; p<nsph; p++){
					double x = spheres(0,p);
					double y = spheres(1,p);
					double z = spheres(2,p);
					double radius = spheres(3,p);
					double xi = double(i);
					double yj = double(j);
					double zk = double(k);

					double distance = radius - sqrt((xi-x)*(xi-x) + (yj-y)*(yj-y) + (zk-z)*(zk-z));

					if (distance > Averages->SDs(i,j,k))		Averages->SDs(i,j,k) = distance;
				}
			}				
		}
	}
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;	
				// Initialize phase positions
				if (Averages->SDs(i,j,k) < 0.0){
					id[n] = 0;
				}
				else{
					id[n] = 2;
					sum++;
				}
			}
		}
	}
	// Compute the pore volume
	sum_local = 0.0;
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				if (id[n] > 0){
					sum_local += 1.0;
				}
			}
		}
	}
	MPI_Allreduce(&sum_local,&pore_vol,1,MPI_DOUBLE,MPI_SUM,comm);
	if (rank==0) printf("Pore volume = %f \n",pore_vol/double(Nx*Ny*Nz));
	//.........................................................
	// don't perform computations at the eight corners
	id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
	id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;
	//.........................................................

    sprintf(LocalRankFilename,"SignDist.%05i",rank);
    FILE *DIST = fopen(LocalRankFilename,"wb");
    fwrite(Averages->SDs.data(),8,Averages->SDs.length(),DIST);
    fclose(DIST);

	sprintf(LocalRankFilename,"ID.%05i",rank);
	FILE *ID = fopen(LocalRankFilename,"wb");
	fwrite(id,1,N,ID);
	fclose(ID);

	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
