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
#include "analysis/Minkowski.h"
#include "analysis/pmmc.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "analysis/analysis.h"

#include "shared_ptr.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"

#define PI 3.14159265359

// Constructor
Minkowski::Minkowski(std::shared_ptr <Domain> dm):
	kstart(0), kfinish(0), isovalue(0), Volume(0),
    LOGFILE(NULL), Dm(dm), Vi(0), Vi_global(0)
{
	Nx=dm->Nx; Ny=dm->Ny; Nz=dm->Nz;
	Volume=double((Nx-2)*(Ny-2)*(Nz-2))*double(Dm->nprocx()*Dm->nprocy()*Dm->nprocz());
	
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

double Minkowski::V(){
	return Vi_global;
}

double Minkowski::A(){
	return Ai_global;
}

double Minkowski::J(){
	return Ji_global;
}

double Minkowski::X(){
	return Xi_global;
}

void Minkowski::ComputeScalar(const DoubleArray Field, const double isovalue)
{
	Xi = Ji = Ai = 0.0;
	DECL object;
	Point P1,P2,P3;
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
					P1 = object.vertex.coords(object.halfedge.v1(e1));
					P2 = object.vertex.coords(object.halfedge.v1(e2));
					P3 = object.vertex.coords(object.halfedge.v1(e3));
					// Surface area
					s1 = sqrt((P1.x-P2.x)*(P1.x-P2.x)+(P1.y-P2.y)*(P1.y-P2.y)+(P1.z-P2.z)*(P1.z-P2.z));
					s2 = sqrt((P2.x-P3.x)*(P2.x-P3.x)+(P2.y-P3.y)*(P2.y-P3.y)+(P2.z-P3.z)*(P2.z-P3.z));
					s3 = sqrt((P1.x-P3.x)*(P1.x-P3.x)+(P1.y-P3.y)*(P1.y-P3.y)+(P1.z-P3.z)*(P1.z-P3.z));
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
	MPI_Barrier(Dm->Comm);
	// Phase averages
	MPI_Allreduce(&Vi,&Vi_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	MPI_Allreduce(&Xi,&Xi_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	MPI_Allreduce(&Ai,&Ai_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	MPI_Allreduce(&Ji,&Ji_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	MPI_Barrier(Dm->Comm);

}

void Minkowski::PrintAll()
{
	if (Dm->rank()==0){
		fprintf(LOGFILE,"%.5g %.5g %.5g %.5g\n",Vi_global, Ai_global, Ji_global, Xi_global);			// minkowski measures
		fflush(LOGFILE);
	}
}

