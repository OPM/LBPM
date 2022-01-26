#include <iostream>
#include <math.h>
#include "analysis/Minkowski.h"
#include "analysis/morphology.h"
#include "common/Domain.h"
#include "common/SpherePack.h"

#include "ProfilerApp.h"


/*
 *  Compare the measured and analytical curvature for a sphere
 *
 */
std::shared_ptr<Database> loadInputs( )
{
	//auto db = std::make_shared<Database>( "Domain.in" );
	auto db = std::make_shared<Database>();
	db->putScalar<int>( "BC", 0 );
	db->putVector<int>( "nproc", { 1, 1, 1 } );
	db->putVector<int>( "n", { 64, 64, 64 } );
	db->putScalar<int>( "nspheres", 1 );
	db->putVector<double>( "L", { 1, 1, 1 } );
	return db;
}

int main(int argc, char **argv)
{
	Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
	int toReturn = 0;
	{
		int i,j,k;
		// Load inputs
		auto db = loadInputs( );
		int Nx = db->getVector<int>( "n" )[0];
		int Ny = db->getVector<int>( "n" )[1];
		int Nz = db->getVector<int>( "n" )[2];
		auto Dm = std::make_shared<Domain>( db, comm );
		Dm->CommInit();
		
		Nx+=2; Ny+=2; Nz+=2;
		int N = Nx*Ny*Nz;
		signed char *ID;
		ID = new signed char[N];
		DoubleArray Phase(Nx,Ny,Nz);
		DoubleArray Distance(Nx,Ny,Nz);

		Minkowski sphere(Dm);

		printf("Set distance map for a sphere \n");
		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int n = k*Nx*Ny + j*Nx + i;
					// put sphere with center at (0.8*Nx, 0.5*Ny, 0.5*Nz)
					if (i+1 > 0.3*(Nx-2)){
						Distance(i,j,k) = 0.4*(Nx-2) - sqrt((1.0*(i-1)-0.8*(Nx-2))*(1.0*i-0.8*(Nx-2))+(1.0*(j-1)-
								0.5*(Ny-2))*(1.0*(j-1)-0.5*(Ny-2))+(1.0*(k-1)-0.5*(Nz-2))*(1.0*(k-1)-0.5*(Nz-2)));
					}
					// reflection sphere center at (-0.2*Nx, 0.5*Ny, 0.5*Nz)
					else {
						Distance(i,j,k) = 0.4*(Nx-2) - sqrt((1.0*(i-1)+0.2*(Nx-2))*(1.0*(i-1)+0.2*(Nx-2))+(1.0*(j-1)-
								0.5*(Ny-2))*(1.0*(j-1)-0.5*(Ny-2))+(1.0*(k-1)-0.5*(Nz-2))*(1.0*(k-1)-0.5*(Nz-2)));
					}
					/* set the labels */
					if (Distance(i,j,k) > 0.0){
						ID[n] = 2;
					}
					else {
						ID[n] = 0;
					}
				}
			}
		}	
	   
		printf("Perform morphological opening \n");
		signed char ErodeLabel = 2;
		signed char NewLabel = 1;
		double Saturation = 0.9;
		MorphOpen(Distance,ID,Dm,Saturation,ErodeLabel,NewLabel);
		
		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int n = k*Nx*Ny + j*Nx + i;
					sphere.id(i,j,k) = 1;
					if (ID[n] == 1)
						sphere.id(i,j,k) = 0;
				}
			}
		}
		
		printf("   Measure the opening \n");
		sphere.MeasureObject();
		//sphere.ComputeScalar(Distance,0.f);
		/* Note 0.856 = (0.95)^3 */
		printf("   Volume  = %f (analytical = %f) \n", sphere.Vi,0.256*0.33333333333333*0.856*3.14159*double((Nx-2)*(Nx-2)*(Nx-2)));
		double error = fabs(sphere.Vi - 0.256*0.856*0.33333333333333*3.14159*double((Nx-2)*(Nx-2)*(Nx-2)))/ (0.256*0.33333333333333*3.14159*double((Nx-2)*(Nx-2)*(Nx-2)));
		printf("   Relative error %f \n", error);
		if (error  > 0.05){
			toReturn = 10;
			printf("ERROR (test_morph): difference from analytical volume is too large\n");
			
			FILE *OUTFILE;
			OUTFILE = fopen("test_morph_ID.raw","wb");
			fwrite(ID,1,N,OUTFILE);
			fclose(OUTFILE);

		}
		else {
			printf("SUCCESS\n");
		}
	}
	PROFILE_SAVE("test_morph");
	Utilities::shutdown();
	return toReturn;
}
