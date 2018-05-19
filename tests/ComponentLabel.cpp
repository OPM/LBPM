// Sequential component labeling for two phase systems
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2015

#include <iostream>
#include <math.h>
#include "analysis/analysis.h"
#include "analysis/TwoPhase.h"

#define NUM_AVERAGES 30

using namespace std;


inline void ReadFromRank(char *FILENAME, DoubleArray &Phase, DoubleArray &Pressure, DoubleArray &Vel_x, 
							DoubleArray &Vel_y, DoubleArray &Vel_z, int nx, int ny, int nz, int iproc, int
							jproc, int kproc)
{
	int i,j,k,q,n,N;
	int iglobal,jglobal,kglobal;
	double value;
	double denA,denB;
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double vx,vy,vz;
	
	N = nx*ny*nz;
	
	double *Den, *DistEven, *DistOdd;
	
	Den = new double[2*N];
	DistEven = new double[10*N];
	DistOdd = new double[9*N];

	ifstream File(FILENAME,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		File.read((char*) &value, sizeof(value));
		Den[n] = value;
		//	if (n== 66276)	printf("Density a  = %f \n",value);
		File.read((char*) &value, sizeof(value));
		Den[N+n] = value;

		//	if (n== 66276)	printf("Density b  = %f \n",value);
		// Read the even distributions
		for (q=0; q<10; q++){
			File.read((char*) &value, sizeof(value));
			DistEven[q*N+n] = value;
		}
		// Read the odd distributions
		for (q=0; q<9; q++){
			File.read((char*) &value, sizeof(value));
			DistOdd[q*N+n] = value;
		}
	}
	File.close();
	
	// Compute the phase field, pressure and velocity
	for (k=1; k<nz-1; k++){
		for (j=1; j<ny-1; j++){
			for (i=1; i<nz-1; i++){
				//........................................................................
				n = k*nx*ny+j*nx+i;
				//........................................................................
				denA = Den[n];
				denB = Den[N+n];
				//........................................................................
				f0 = DistEven[n];
				f2 = DistEven[N+n];
				f4 = DistEven[2*N+n];
				f6 = DistEven[3*N+n];
				f8 = DistEven[4*N+n];
				f10 = DistEven[5*N+n];
				f12 = DistEven[6*N+n];
				f14 = DistEven[7*N+n];
				f16 = DistEven[8*N+n];
				f18 = DistEven[9*N+n];
				//........................................................................
				f1 = DistOdd[n];
				f3 = DistOdd[1*N+n];
				f5 = DistOdd[2*N+n];
				f7 = DistOdd[3*N+n];
				f9 = DistOdd[4*N+n];
				f11 = DistOdd[5*N+n];
				f13 = DistOdd[6*N+n];
				f15 = DistOdd[7*N+n];
				f17 = DistOdd[8*N+n];
				//........................................................................
				//.................Compute the pressure....................................
				value = 0.3333333333333333*(f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17);
				//........................................................................
				//.................Compute the velocity...................................
				vx = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
				vy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
				vz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
				//........................................................................
				// save values in global arrays
				//........................................................................
				iglobal = iproc*(nx-2)+i;
				jglobal = jproc*(ny-2)+j;
				kglobal = kproc*(nz-2)+k;
				//........................................................................				
				Phase(iglobal,jglobal,kglobal) = (denA-denB)/(denA+denB);
				Pressure(iglobal,jglobal,kglobal) = value;
				Vel_x(iglobal,jglobal,kglobal) = vx;
				Vel_y(iglobal,jglobal,kglobal) = vy;
				Vel_z(iglobal,jglobal,kglobal) = vz;
				//........................................................................
			}
		}
	}
	
	delete Den;
	delete DistEven;
	delete DistOdd;
}

int main(int argc, char **argv)
{
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);

	printf("----------------------------------------------------------\n");
	printf("COMPUTING TCAT ANALYSIS FOR NON-WETTING PHASE FEATURES \n");
	printf("----------------------------------------------------------\n");

	if (nprocs != 1) INSIST(nprocs == 1,"Error: ComponentLabel --serial case!");

	//.......................................................................
	int nprocx,nprocy,nprocz;
	int Nx, Ny, Nz;
	int nx,ny,nz;
	int nspheres;
	double Lx,Ly,Lz;
	//.......................................................................
	int i,j,k,n;
	int iproc,jproc,kproc;
	//.......................................................................
	// Reading the domain information file
	//.......................................................................
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
	//.......................................................................

	nx+=2;
	ny+=2;
	nz+=2;
	
	nprocs = nprocx*nprocy*nprocz;
	printf("Number of MPI ranks: %i \n", nprocs);
	int BoundaryCondition=0;
	Nx = (nx-2)*nprocx;
	Ny = (ny-2)*nprocy;
	Nz = (nz-2)*nprocz;
	Domain Dm(Nx,Ny,Nz,rank,1,1,1,Lx,Ly,Lz,BoundaryCondition);
	Nx+=2; Ny+=2; Nz+=2;
	printf("Full domain size: %i x %i x %i  \n", Nx,Ny,Nz);
	
	
	TwoPhase Averages(Dm);

	// Filenames used
	char LocalRankString[8];
	char LocalRankFilename[40];
	char BaseFilename[20];
	sprintf(BaseFilename,"%s","dPdt.");
	                 
	int proc,iglobal,kglobal,jglobal;

	double * Temp;
	Temp = new double[nx*ny*nz];

	// read the files and populate main arrays
	for ( kproc=0; kproc<nprocz; kproc++){
		for ( jproc=0; jproc<nprocy; jproc++){
			for ( iproc=0; iproc<nprocx; iproc++){
				
				proc = kproc*nprocx*nprocy + jproc*nprocx + iproc;

				sprintf(LocalRankString,"%05d",proc);
	//		sprintf(LocalRankFilename,"%s%s","dPdt.",LocalRankString);
	//			printf("Reading file %s \n",LocalRankFilename);
				ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);
				for (k=1; k<nz-1; k++){
					for (j=1; j<ny-1; j++){
						for (i=1; i<nz-1; i++){
							//........................................................................
							n = k*nx*ny+j*nx+i;
							//........................................................................
							iglobal = iproc*(nx-2)+i;
							jglobal = jproc*(ny-2)+j;
							kglobal = kproc*(nz-2)+k;
							//........................................................................
							Averages.dPdt(iglobal,jglobal,kglobal) = 0.0;
							Averages.Phase_tplus(iglobal,jglobal,kglobal) = 0.0;
							Averages.Phase_tminus(iglobal,jglobal,kglobal) = 0.0;
							//........................................................................
						}
					}
				}
				
				sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
		//		printf("Reading file %s \n",LocalRankFilename);
		//		printf("Sub-domain size: %i x %i x %i  \n", nx,ny,nz);
				ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);	
				for (k=1; k<nz-1; k++){
					for (j=1; j<ny-1; j++){
						for (i=1; i<nz-1; i++){

							//........................................................................
							n = k*nx*ny+j*nx+i;
							//........................................................................
							iglobal = iproc*(nx-2)+i;
							jglobal = jproc*(ny-2)+j;
							kglobal = kproc*(nz-2)+k;
							//........................................................................
							Averages.SDs(iglobal,jglobal,kglobal) = Temp[n];
							//........................................................................
						}
					}
				}
				
				sprintf(LocalRankFilename,"%s%s","Restart.",LocalRankString);

				ReadFromRank(LocalRankFilename,Averages.Phase,Averages.Press,
						Averages.Vel_x,Averages.Vel_y,Averages.Vel_z,
						nx,ny,nz,iproc,jproc,kproc);

/*				sprintf(LocalRankFilename,"%s%s","Pressure.",LocalRankString);
				
				ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);	
				for (k=1; k<nz-1; k++){
					for (j=1; j<ny-1; j++){
						for (i=1; i<nx-1; i++){

							//........................................................................
							n = k*nx*ny+j*nx+i;
							//........................................................................
							iglobal = iproc*(nx-2)+i;
							jglobal = jproc*(ny-2)+j;
							kglobal = kproc*(nz-2)+k;
							//........................................................................
							Averages.Press(iglobal,jglobal,kglobal) = Temp[n];
							//........................................................................
						}
					}
				}
*/
				//sprintf(LocalRankFilename,"%s%s","Phase.",LocalRankString);
				//ReadBinaryFile(LocalRankFilename, Temp, nx*ny*nz);
				for (k=1; k<nz-1; k++){
					for (j=1; j<ny-1; j++){
						for (i=1; i<nx-1; i++){

							//........................................................................
							n = k*nx*ny+j*nx+i;
							//........................................................................
							iglobal = iproc*(nx-2)+i;
							jglobal = jproc*(ny-2)+j;
							kglobal = kproc*(nz-2)+k;
							//........................................................................
							//Averages.Phase(iglobal,jglobal,kglobal) = Temp[n];
							Averages.SDn(iglobal,jglobal,kglobal) = Averages.Phase(iglobal,jglobal,kglobal);//Temp[n];
							//........................................................................
						}
					}
				}
				
				
			}
		}
	}
	printf("Read %i ranks of %s \n",nprocs,BaseFilename);
	
	delete Temp;
	
	// Initializing the blob ID
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				n = k*Nx*Ny+j*Nx+i;
				if (Averages.SDs(i,j,k) < 0.0){
					// Solid phase 
					Dm.id[n] = 0;
					//PhaseLabel(i,j,k) = 0;
					//WP(i,j,k) = -2;
					//NWP(i,j,k) = -2;
				}
				else if (Averages.Phase(i,j,k) < 0.0){
					// wetting phase
					Dm.id[n] = 2;
					//PhaseLabel(i,j,k) = 2;
					//WP(i,j,k) = -2;
					//NWP(i,j,k) = -1;
				}
				else {
					// non-wetting phase
					Dm.id[n] = 1;
					//PhaseLabel(i,j,k) = 1;
					//WP(i,j,k) = -1;
					//NWP(i,j,k) = -2;
				}
			}
		}
	}

	// Compute the porosity
	double porosity=0.0;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				if (Dm.id[k*Nx*Ny+j*Nx+i] != 0){
					porosity += 1.0;
				}
			}
		}
	}
	porosity /= (Nx*Ny*Nz*1.0);
	printf("Media porosity is %f \n",porosity);
	Dm.CommInit();

	/* ****************************************************************
				IDENTIFY ALL COMPONENTS FOR BOTH PHASES
	****************************************************************** */
   // int number_NWP_components = ComputeLocalPhaseComponent(PhaseLabel,1,NWP,true);
    //int number_WP_components = ComputeLocalPhaseComponent(PhaseLabel,2,WP,true);

    //printf("Number of WP components = %i \n",number_WP_components);
    //printf("Number of NWP components = %i \n",number_NWP_components);
	
	// Map the signed distance for the analysis
	
	// Compute the porosity
	porosity=0.0;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				if (Averages.SDs(i,j,k) > 0.0){
					porosity += 1.0;
				}
				Averages.SDs(i,j,k) -= (1.0);
			}
		}
	}
	porosity /= (Nx*Ny*Nz*1.0);
	//printf("Media porosity is %f \n",porosity);

	double beta=0.95;
	int timestep=5;
	Averages.Initialize();
	Averages.ComputeDelPhi();
	Averages.ColorToSignedDistance(beta,Averages.Phase,Averages.SDn);
	Averages.UpdateMeshValues();
	Averages.ComponentAverages();
	Averages.SortBlobs();
	Averages.PrintComponents(timestep);

    // Create the MeshDataStruct
    fillHalo<double> fillData(Dm.Comm,Dm.rank_info,{Nx-2,Ny-2,Nz-2},{1,1,1},0,1);
    std::vector<IO::MeshDataStruct> meshData(1);
    meshData[0].meshName = "domain";
    meshData[0].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(Dm.rank_info,Nx-2,Ny-2,Nz-2,Lx,Ly,Lz) );
    std::shared_ptr<IO::Variable> PhaseVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> SignDistVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> LabelWPVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> LabelNWPVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> PhaseIDVar( new IO::Variable() );

    PhaseVar->name = "phase";
    PhaseVar->type = IO::VariableType::VolumeVariable;
    PhaseVar->dim = 1;
    PhaseVar->data.resize(Nx-2,Ny-2,Nz-2);
    meshData[0].vars.push_back(PhaseVar);

    SignDistVar->name = "SignDist";
    SignDistVar->type = IO::VariableType::VolumeVariable;
    SignDistVar->dim = 1;
    SignDistVar->data.resize(Nx-2,Ny-2,Nz-2);
    meshData[0].vars.push_back(SignDistVar);

    LabelNWPVar->name = "LabelNWP";
    LabelNWPVar->type = IO::VariableType::VolumeVariable;
    LabelNWPVar->dim = 1;
    LabelNWPVar->data.resize(Nx-2,Ny-2,Nz-2);
    meshData[0].vars.push_back(LabelNWPVar);

    LabelWPVar->name = "LabelWP";
    LabelWPVar->type = IO::VariableType::VolumeVariable;
    LabelWPVar->dim = 1;
    LabelWPVar->data.resize(Nx-2,Ny-2,Nz-2);
    meshData[0].vars.push_back(LabelWPVar);

    PhaseIDVar->name = "PhaseID";
    PhaseIDVar->type = IO::VariableType::VolumeVariable;
    PhaseIDVar->dim = 1;
    PhaseIDVar->data.resize(Nx-2,Ny-2,Nz-2);
    meshData[0].vars.push_back(PhaseIDVar);

    fillData.copy(Averages.SDn,PhaseVar->data);
    fillData.copy(Averages.SDs,SignDistVar->data);
    fillData.copy(Averages.Label_WP,LabelWPVar->data);
    fillData.copy(Averages.Label_NWP,LabelNWPVar->data);
    fillData.copy(Averages.PhaseID,PhaseIDVar->data);
    IO::writeData( 0, meshData, comm );
/*
	FILE *NWP_FILE;
	NWP_FILE = fopen("NWP.dat","wb");
	fwrite(Averages.Label_NWP.get(),4,Nx*Ny*Nz,NWP_FILE);
	fclose(NWP_FILE);

	FILE *WP_FILE;
	WP_FILE = fopen("WP.dat","wb");
	fwrite(Averages.Label_WP.get(),4,Nx*Ny*Nz,WP_FILE);
	fclose(WP_FILE);
	
	FILE *DISTANCE;
	DISTANCE = fopen("SignDist.dat","wb");
	fwrite(Averages.SDs.get(),8,Nx*Ny*Nz,DISTANCE);
	fclose(DISTANCE);
	*/
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}

