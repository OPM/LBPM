/* Flow adaptor class for multiphase flow methods */

#include "analysis/FlowAdaptor.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"

Membrane::Membrane(ScaLBL_Communicator &ScaLBL, int *neighborList) {

}

Membrane::~Membrane() {
	
}

int Membrane::Create(DoubleArray &Distance, IntArray &Map, int *neighborList, int *membrane, int Np){
	int mlink = 0;
	int n, idx, neighbor;
	double dist, locdist;
	/* go through the neighborlist structure */
	for (idx=0; idx<last_interior; idx++){
		if (idx == next) idx=first_interior; // skip to interior
		
		
	}
	/* count & cut the links */
	for (k=1;k<Nz-1;k++){
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				n=k*Nx*Ny+j*Nx+i;
				idx=Map(i,j,k);
				locdist=Distance(i,j,k);

				else if (!(idx<0)){
					
					int neighbor;    // cycle through the neighbors of lattice site idx
					neighbor=Map(i-1,j,k);
					dist=Distance(i-1,j,k);
					if (dist*locdist < 0.0){
						neighborList[idx]=idx + 2*Np;
					}

					neighbor=Map(i+1,j,k);
					dist=Distance(i+1,j,k);
					if (dist*locdist < 0.0){
						neighborList[Np+idx] = idx + 1*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k);
					dist=Distance(i,j-1,k);
					if (dist*locdist < 0.0){
						neighborList[2*Np+idx]=idx + 4*Np;
					}

					neighbor=Map(i,j+1,k);
					dist=Distance(i,j+1,k);
					if (dist*locdist < 0.0){
						neighborList[3*Np+idx]=idx + 3*Np;
						mlink++;
					}

					neighbor=Map(i,j,k-1);
					dist=Distance(i,j,k-1);
					if (dist*locdist < 0.0){
						neighborList[4*Np+idx]=idx + 6*Np;
					}

					neighbor=Map(i,j,k+1);
					dist=Distance(i,j,k+1);
					if (dist*locdist < 0.0){
						neighborList[5*Np+idx]=idx + 5*Np;
						mlink++;
					}

					neighbor=Map(i-1,j-1,k);
					dist=Distance(i-1,j-1,k);
					if (dist*locdist < 0.0){
						neighborList[6*Np+idx]=idx + 8*Np;
					}

					neighbor=Map(i+1,j+1,k);
					dist=Distance(i+1,j+1,k);
					if (dist*locdist < 0.0){
						neighborList[7*Np+idx]=idx + 7*Np;
						mlink++;
					}

					neighbor=Map(i-1,j+1,k);
					dist=Distance(i-1,j+1,k);
					if (dist*locdist < 0.0){
						neighborList[8*Np+idx]=idx + 10*Np;
					}

					neighbor=Map(i+1,j-1,k);
					dist=Distance(i+1,j-1,k);
					if (dist*locdist < 0.0){
						neighborList[9*Np+idx]=idx + 9*Np;
						mlink++;
					}

					neighbor=Map(i-1,j,k-1);
					dist=Distance(i-1,j,k-1);
					if (dist*locdist < 0.0){
						neighborList[10*Np+idx]=idx + 12*Np;
					}

					neighbor=Map(i+1,j,k+1);
					dist=Distance(i+1,j,k+1);
					if (dist*locdist < 0.0){
						neighborList[11*Np+idx]=idx + 11*Np;
						mlink++;
					}

					neighbor=Map(i-1,j,k+1);
					dist=Distance(i-1,j,k+1);
					if (dist*locdist < 0.0){
						neighborList[12*Np+idx]=idx + 14*Np;
					}

					neighbor=Map(i+1,j,k-1);
					dist=Distance(i+1,j,k-1);
					if (dist*locdist < 0.0){
						neighborList[13*Np+idx]=idx + 13*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k-1);
					dist=Distance(i,j-1,k-1);
					if (dist*locdist < 0.0){
						neighborList[14*Np+idx]=idx + 16*Np;
					}

					neighbor=Map(i,j+1,k+1);
					dist=Distance(i,j+1,k+1);
					if (dist*locdist < 0.0){
						neighborList[15*Np+idx]=idx + 15*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k+1);
					dist=Distance(i,j-1,k+1);
					if (dist*locdist < 0.0){
						neighborList[16*Np+idx]=idx + 18*Np;
					}

					neighbor=Map(i,j+1,k-1);
					dist=Distance(i,j+1,k-1);
					if (dist*locdist < 0.0){
						neighborList[17*Np+idx]=idx + 17*Np;
						mlink++;
					}
				}
			}
		}
	}
	
	/* allocate memory */
	membrane = new int [2*mlink];
	
	/* construct the membrane*/
	mlink = 0;
	for (k=1;k<Nz-1;k++){
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				n=k*Nx*Ny+j*Nx+i;
				idx=Map(i,j,k);
				locdist=Distance(i,j,k);

				else if (!(idx<0)){
					
					int neighbor;    // cycle through the neighbors of lattice site idx
					neighbor=Map(i-1,j,k);
					dist=Distance(i-1,j,k);
					if (dist*locdist < 0.0){
						neighborList[idx]=idx + 2*Np;
						//membrane[2*mlink] = idx + 2*Np;
						//membrane[2*mlink+1] = neighbor + 1*Np;
						//mlink++;
					}

					neighbor=Map(i+1,j,k);
					dist=Distance(i+1,j,k);
					if (dist*locdist < 0.0){
						neighborList[Np+idx] = idx + 1*Np;
						membrane[2*mlink] = idx + 1*Np;
						membrane[2*mlink+1] = neighbor + 2*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k);
					dist=Distance(i,j-1,k);
					if (dist*locdist < 0.0){
						neighborList[2*Np+idx]=idx + 4*Np;
						//membrane[2*mlink] = idx + 4*Np;
						//membrane[2*mlink+1] = neighbor + 3*Np;
						//mlink++;
					}

					neighbor=Map(i,j+1,k);
					dist=Distance(i,j+1,k);
					if (dist*locdist < 0.0){
						neighborList[3*Np+idx]=idx + 3*Np;
						membrane[2*mlink] = idx + 3*Np;
						membrane[2*mlink+1] = neighbor + 4*Np;
						mlink++;
					}

					neighbor=Map(i,j,k-1);
					dist=Distance(i,j,k-1);
					if (dist*locdist < 0.0){
						neighborList[4*Np+idx]=idx + 6*Np;
						//membrane[2*mlink] = idx + 6*Np;
						//membrane[2*mlink+1] = neighbor + 5*Np;
						//mlink++;
					}

					neighbor=Map(i,j,k+1);
					dist=Distance(i,j,k+1);
					if (dist*locdist < 0.0){
						neighborList[5*Np+idx]=idx + 5*Np;
						membrane[2*mlink] = idx + 5*Np;
						membrane[2*mlink+1] = neighbor + 6*Np;
						mlink++;
					}

					neighbor=Map(i-1,j-1,k);
					dist=Distance(i-1,j-1,k);
					if (dist*locdist < 0.0){
						neighborList[6*Np+idx]=idx + 8*Np;
						//membrane[2*mlink] = idx + 8*Np;
						//membrane[2*mlink+1] = neighbor + 7*Np
						//mlink++;
					}

					neighbor=Map(i+1,j+1,k);
					dist=Distance(i+1,j+1,k);
					if (dist*locdist < 0.0){
						neighborList[7*Np+idx]=idx + 7*Np;
						membrane[2*mlink] = idx + 7*Np;
						membrane[2*mlink+1] = neighbor+8*Np;
						mlink++;
					}

					neighbor=Map(i-1,j+1,k);
					dist=Distance(i-1,j+1,k);
					if (dist*locdist < 0.0){
						neighborList[8*Np+idx]=idx + 10*Np;
						//membrane[2*mlink] = idx + 10*Np;
						//membrane[2*mlink+1] = neighbor + 9*Np;
						//mlink++;
					}

					neighbor=Map(i+1,j-1,k);
					dist=Distance(i+1,j-1,k);
					if (dist*locdist < 0.0){
						neighborList[9*Np+idx]=idx + 9*Np;
						membrane[2*mlink] = idx + 9*Np;
						membrane[2*mlink+1] = neighbor + 10*Np;
						mlink++;
					}

					neighbor=Map(i-1,j,k-1);
					dist=Distance(i-1,j,k-1);
					if (dist*locdist < 0.0){
						neighborList[10*Np+idx]=idx + 12*Np;
						//membrane[2*mlink] = idx + 12*Np;
						//membrane[2*mlink+1] = neighbor + 11*Np;;
						//mlink++;
					}

					neighbor=Map(i+1,j,k+1);
					dist=Distance(i+1,j,k+1);
					if (dist*locdist < 0.0){
						neighborList[11*Np+idx]=idx + 11*Np;
						membrane[2*mlink] = idx + 11*Np;
						membrane[2*mlink+1] = neighbor + 12*Np;
						mlink++;
					}

					neighbor=Map(i-1,j,k+1);
					dist=Distance(i-1,j,k+1);
					if (dist*locdist < 0.0){
						neighborList[12*Np+idx]=idx + 14*Np;
						//membrane[2*mlink] = idx + 14*Np;
						//membrane[2*mlink+1] = neighbor + 13*Np;
						//mlink++;
					}

					neighbor=Map(i+1,j,k-1);
					dist=Distance(i+1,j,k-1);
					if (dist*locdist < 0.0){
						neighborList[13*Np+idx]=idx + 13*Np;
						membrane[2*mlink] = idx + 13*Np;
						membrane[2*mlink+1] = neighbor + 14*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k-1);
					dist=Distance(i,j-1,k-1);
					if (dist*locdist < 0.0){
						neighborList[14*Np+idx]=idx + 16*Np;
						//membrane[2*mlink] = idx + 16*Np;
						//membrane[2*mlink+1] = neighbor + 15*Np;
						//mlink++;
					}

					neighbor=Map(i,j+1,k+1);
					dist=Distance(i,j+1,k+1);
					if (dist*locdist < 0.0){
						neighborList[15*Np+idx]=idx + 15*Np;
						membrane[2*mlink] = idx + 15*Np;
						membrane[2*mlink+1] =neighbor + 16*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k+1);
					dist=Distance(i,j-1,k+1);
					if (dist*locdist < 0.0){
						neighborList[16*Np+idx]=idx + 18*Np;
						//membrane[2*mlink] = idx + 18*Np;
						//membrane[2*mlink+1] = neighbor + 17*Np;
						//mlink++;
					}

					neighbor=Map(i,j+1,k-1);
					dist=Distance(i,j+1,k-1);
					if (dist*locdist < 0.0){
						neighborList[17*Np+idx]=idx + 17*Np;
						membrane[2*mlink] = idx + 17*Np;
						membrane[2*mlink+1] = neighbor + 18*Np;
						mlink++;
					}
				}
			}
		}
	}
	
	return mlink;
}
