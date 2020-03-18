#include "analysis/Minkowski.h"
#include "analysis/pmmc.h"
#include "analysis/analysis.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"

#include "ProfilerApp.h"

#include <memory>


#define PI 3.14159265359

// Constructor
Minkowski::Minkowski(std::shared_ptr <Domain> dm):
	kstart(0), kfinish(0), isovalue(0), Volume(0),
    LOGFILE(NULL), Dm(dm), Vi(0), Vi_global(0)
{
	Nx=dm->Nx; Ny=dm->Ny; Nz=dm->Nz;
	Volume=double((Nx-2)*(Ny-2)*(Nz-2))*double(Dm->nprocx()*Dm->nprocy()*Dm->nprocz());
	
	id.resize(Nx,Ny,Nz);         	id.fill(0);
	label.resize(Nx,Ny,Nz);         label.fill(0);
	distance.resize(Nx,Ny,Nz);      distance.fill(0);

	if (Dm->rank()==0){
		LOGFILE = fopen("minkowski.csv","a+");
		if (fseek(LOGFILE,0,SEEK_SET) == fseek(LOGFILE,0,SEEK_CUR))
		{
			// If LOGFILE is empty, write a short header to list the averages
			//fprintf(LOGFILE,"--------------------------------------------------------------------------------------\n");
			fprintf(LOGFILE,"Vn An Jn Xn\n"); 			//miknowski measures,
		}
	}
}


// Destructor
Minkowski::~Minkowski()
{
    if ( LOGFILE!=NULL ) { fclose(LOGFILE); }
}

void Minkowski::ComputeScalar(const DoubleArray& Field, const double isovalue)
{
    PROFILE_START("ComputeScalar");
	Xi = Ji = Ai = 0.0;
	DECL object;
	int e1,e2,e3;
	double s,s1,s2,s3;
	double a1,a2,a3;
	//double Vx,Vy,Vz,Wx,Wy,Wz,nx,ny,nz,norm;
	//int Nx = Field.size(0);
	//int Ny = Field.size(1);
	//int Nz = Field.size(2);
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				object.LocalIsosurface(Field,isovalue,i,j,k);
				for (int idx=0; idx<object.TriangleCount; idx++){
					e1 = object.Face(idx); 
					e2 = object.halfedge.next(e1);
					e3 = object.halfedge.next(e2);
					auto P1 = object.vertex.coords(object.halfedge.v1(e1));
					auto P2 = object.vertex.coords(object.halfedge.v1(e2));
					auto P3 = object.vertex.coords(object.halfedge.v1(e3));
					// Surface area
					s1 = Distance( P1, P2 );
					s2 = Distance( P2, P3 );
					s3 = Distance( P1, P3 );
					s = 0.5*(s1+s2+s3);
					Ai += sqrt(s*(s-s1)*(s-s2)*(s-s3));
					// Mean curvature based on half edge angle
					a1 = object.EdgeAngle(e1);
					a2 = object.EdgeAngle(e2);
					a3 = object.EdgeAngle(e3);
					Ji += (a1*s1+a2*s2+a3*s3);
					//if (0.08333333333333*(a1*s1+a2*s2+a3*s3) < 0.f){
					//double intcurv=0.08333333333333*(a1*s1+a2*s2+a3*s3);
					//double surfarea=sqrt(s*(s-s1)*(s-s2)*(s-s3));
					//printf("   (%i,%i,%i) PQ(%i,%i)={%f,%f,%f} {%f,%f,%f} a=%f l=%f \n",i,j,k,e1,object.halfedge.twin(e1),P1.x,P1.y,P1.z,P2.x,P2.y,P2.z,a1,s1);
					// printf("   (%i,%i,%i) QR(%i,%i)={%f,%f,%f} {%f,%f,%f} a=%f l=%f \n",i,j,k,e2,object.halfedge.twin(e2),P2.x,P2.y,P2.z,P3.x,P3.y,P3.z,a2,s2);
					// printf("   (%i,%i,%i) RP(%i,%i)={%f,%f,%f} {%f,%f,%f} a=%f l=%f \n",i,j,k,e3,object.halfedge.twin(e3),P3.x,P3.y,P3.z,P1.x,P1.y,P1.z,a3,s3);
					  //}
					// Euler characteristic (half edge rule: one face - 0.5*(three edges))
					Xi -= 0.5;
				}
				// Euler characteristic -- each vertex shared by four cubes
				Xi += 0.25*double(object.VertexCount);
			}
		}
	}
	// Voxel counting for volume fraction
	Vi = 0.f;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				if (Field(i,j,k) < isovalue){
					Vi += 1.0;
				}
			}
		}
	}
	// convert X for 2D manifold to 3D object
	Xi *= 0.5;
	
	MPI_Barrier(Dm->Comm);
	// Phase averages
	MPI_Allreduce(&Vi,&Vi_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	MPI_Allreduce(&Xi,&Xi_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	MPI_Allreduce(&Ai,&Ai_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	MPI_Allreduce(&Ji,&Ji_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	MPI_Barrier(Dm->Comm);
    PROFILE_STOP("ComputeScalar");
}


void Minkowski::MeasureObject(){
	/*
	 *  compute the distance to an object 
	 * 
	 * THIS ALGORITHM ASSUMES THAT id() is populated with phase id to distinguish objects
	 *    0 - labels the object
	 *    1 - labels the rest of the 
	 */

	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				distance(i,j,k) =2.0*double(id(i,j,k))-1.0;
			}
		}
	}	
	CalcDist(distance,id,*Dm);
	ComputeScalar(distance,0.0);

}


int Minkowski::MeasureConnectedPathway(){
	/*
	 * compute the connected pathway for object with LABEL in id field
	 * compute the labels for connected components
	 * compute the distance to the connected pathway
	 * 
	 * THIS ALGORITHM ASSUMES THAT id() is populated with phase id to distinguish objects
	 */
	
	char LABEL = 0;
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				if (id(i,j,k) == LABEL){
					distance(i,j,k) = 1.0;
				}
				else
					distance(i,j,k) = -1.0;
			}
		}
	}
	
	// Extract only the connected part of NWP
	double vF=0.0; 
	n_connected_components = ComputeGlobalBlobIDs(Nx-2,Ny-2,Nz-2,Dm->rank_info,distance,distance,vF,vF,label,Dm->Comm);
//	int n_connected_components = ComputeGlobalPhaseComponent(Nx-2,Ny-2,Nz-2,Dm->rank_info,const IntArray &PhaseID, int &VALUE, BlobIDArray &GlobalBlobID, Dm->Comm )
	MPI_Barrier(Dm->Comm);
	
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				if ( label(i,j,k) == 0){
					id(i,j,k) = 0;
				}
				else{
					id(i,j,k) = 1;
				}
			}
		}
	}
	MeasureObject();
	return n_connected_components; 
}


void Minkowski::PrintAll()
{
	if (Dm->rank()==0){
		fprintf(LOGFILE,"%.5g %.5g %.5g %.5g\n",Vi_global, Ai_global, Ji_global, Xi_global);			// minkowski measures
		fflush(LOGFILE);
	}
}

