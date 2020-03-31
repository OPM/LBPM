#include <iostream>
#include <math.h>
#include "analysis/Minkowski.h"
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
    db->putVector<int>( "n", { 32, 32, 32 } );
    db->putScalar<int>( "nspheres", 1 );
    db->putVector<double>( "L", { 1, 1, 1 } );
    return db;
}

int main(int argc, char **argv)
{
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	//int rank = MPI_WORLD_RANK();
	//int nprocs = MPI_WORLD_SIZE();
	int toReturn = 0;
	{
		int i,j,k;

		// Load inputs
		auto db = loadInputs( );
		int Nx = db->getVector<int>( "n" )[0];
		int Ny = db->getVector<int>( "n" )[1];
		int Nz = db->getVector<int>( "n" )[2];
		auto Dm = std::make_shared<Domain>( db, comm );
		
		Nx+=2; Ny+=2; Nz+=2;
		DoubleArray Phase(Nx,Ny,Nz);
		
		Minkowski plane(Dm);

		printf("Set distance map for plane \n");
		for (k=0; k<Nz; k++){
			for (j=0; j<Ny; j++){
				for (i=0; i<Nx; i++){
				  Phase(i,j,k) = k-0.5*double(Nz+1);
				}
			}
		}

		printf("   Construct local isosurface \n");
		plane.ComputeScalar(Phase,0.f);

		printf("   Surface area  = %f (analytical = %f) \n", plane.Ai,double((Nx-2)*(Ny-2)));
		printf("   Mean curvature  = %f (analytical =0) \n", plane.Ji);
		printf("   Euler characteristic  = %f (analytical = 0) \n",plane.Xi);
		
		Minkowski cylinder(Dm);

		printf("Set distance map for cylinder \n");
		for (k=0; k<Nz; k++){
			for (j=0; j<Ny; j++){
				for (i=0; i<Nx; i++){
					Phase(i,j,k) = sqrt((1.0*i-0.5*Nx)*(1.0*i-0.5*Nx)+(1.0*k-0.5*Nz)*(1.0*k-0.5*Nz))-0.4*Nx;
				}
			}
		}

		printf("   Construct local isosurface \n");
		cylinder.ComputeScalar(Phase,0.f);

		printf("   Surface area  = %f (analytical = %f) \n", cylinder.Ai,2*3.14159*0.4*double((Nx-2)*Nx));
		printf("   Mean curvature  = %f (analytical = %f) \n", cylinder.Ji,2*3.14159*double((Nx-2)));
		printf("   Euler characteristic  = %f (analytical = 0) \n",cylinder.Xi);
	      		
		Minkowski sphere(Dm);

		printf("Set distance map for a sphere \n");
		for (k=0; k<Nz; k++){
			for (j=0; j<Ny; j++){
				for (i=0; i<Nx; i++){
					Phase(i,j,k) = sqrt((1.0*i-0.5*Nx)*(1.0*i-0.5*Nx)+(1.0*j-0.5*Ny)*(1.0*j-0.5*Ny)+(1.0*k-0.5*Nz)*(1.0*k-0.5*Nz))-0.4*Nx;
	      			}
			}
		}

		printf("   Construct local isosurface \n");
		sphere.ComputeScalar(Phase,0.f);

		printf("   Surface area  = %f (analytical = %f) \n", sphere.Ai,4*3.14159*0.16*double(Nx*Nx));
		printf("   Mean curvature  = %f (analytical = %f) \n", sphere.Ji,8*3.14159*0.4*double(Nx));
		printf("   Euler characteristic  = %f (analytical = 2.0) \n",sphere.Xi);
		
	}
    PROFILE_SAVE("test_dcel_minkowski");
	MPI_Barrier(comm);
	MPI_Finalize();
	return toReturn;
}
