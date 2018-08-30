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
#include <math.h>
#include <stdio.h>
#include <cuda_profiler_api.h>

#define NBLOCKS 1024
#define NTHREADS 256

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* address, double val) { 
   unsigned long long int* address_as_ull = (unsigned long long int*)address;
   unsigned long long int old = *address_as_ull, assumed;

   do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val+__longlong_as_double(assumed)));
   } while (assumed != old);
   return __longlong_as_double(old);
}
#endif

__global__ void dvc_ScaLBL_Gradient_Unpack(double weight, double Cqx, double Cqy, double Cqz, 
		int *list, int start, int count, double *recvbuf, double *phi, double *grad, int N){
	//....................................................................................
	// Unpack distribution from the recv buffer
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int n,idx;
	double value, tmp;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
		// Get the index from the list
		n = list[start+idx];
		// unpack the distribution to the proper location
		if (!(n<0)){
			value=weight*(recvbuf[idx] - phi[n]);
			// PARALLEL UPDATE MUST BE DONE ATOMICALLY
			tmp = Cqx*value;
			atomicAdd(&grad[n],tmp);
			tmp = Cqy*value;
			atomicAdd(&grad[N+n],tmp);
			tmp = Cqz*value;
			atomicAdd(&grad[2*N+n],tmp);
		}
	}
}

__global__ void dvc_ScaLBL_DFH_Init(double *Phi, double *Den, double *Aq, double *Bq, int start, int finish, int Np){
	int idx;
	double phi,nA,nB;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		idx =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x+start;
		if (idx<finish) {

			phi = Phi[idx];
			if (phi > 0.f){
				nA = 1.0; nB = 0.f;
			}
			else{
				nB = 1.0; nA = 0.f;
			}
			Den[idx] = nA;
			Den[Np+idx] = nB;

			Aq[idx]=0.3333333333333333*nA;
			Aq[Np+idx]=0.1111111111111111*nA;
			Aq[2*Np+idx]=0.1111111111111111*nA;
			Aq[3*Np+idx]=0.1111111111111111*nA;
			Aq[4*Np+idx]=0.1111111111111111*nA;
			Aq[5*Np+idx]=0.1111111111111111*nA;
			Aq[6*Np+idx]=0.1111111111111111*nA;

			Bq[idx]=0.3333333333333333*nB;
			Bq[Np+idx]=0.1111111111111111*nB;
			Bq[2*Np+idx]=0.1111111111111111*nB;
			Bq[3*Np+idx]=0.1111111111111111*nB;
			Bq[4*Np+idx]=0.1111111111111111*nB;
			Bq[5*Np+idx]=0.1111111111111111*nB;
			Bq[6*Np+idx]=0.1111111111111111*nB;
		}
	}
}


// LBM based on density functional hydrodynamics
__global__ void dvc_ScaLBL_D3Q19_AAeven_DFH(int *neighborList, double *dist, double *Aq, double *Bq, double *Den, double *Phi,
			double *Gradient, double *SolidForce, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
			double Fx, double Fy, double Fz, int start, int finish, int Np){
	int nn,n;
	double fq;
	// conserved momemnts
	double rho,jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m3,m5,m7;
	double nA,nB; // number density
	double a1,b1,a2,b2,nAB,delta;
	double C,nx,ny,nz; //color gradient magnitude and direction
	double ux,uy,uz;
	double phi,tau,rho0,rlx_setA,rlx_setB;
	double force_x,force_y,force_z;

	const double mrt_V1=0.05263157894736842;
	const double mrt_V2=0.012531328320802;
	const double mrt_V3=0.04761904761904762;
	const double mrt_V4=0.004594820384294068;
	const double mrt_V5=0.01587301587301587;
	const double mrt_V6=0.0555555555555555555555555;
	const double mrt_V7=0.02777777777777778;
	const double mrt_V8=0.08333333333333333;
	const double mrt_V9=0.003341687552213868;
	const double mrt_V10=0.003968253968253968;
	const double mrt_V11=0.01388888888888889;
	const double mrt_V12=0.04166666666666666;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			// read the component number densities
			nA = Den[n];
			nB = Den[Np + n];

			// compute phase indicator field
			phi=(nA-nB)/(nA+nB);

			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
			// local relaxation time
			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
			rlx_setA = 1.f/tau;
			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);

			//...........Read the Color Gradient.................................
			nx = Gradient[n];
			ny = Gradient[n+Np];
			nz = Gradient[n+2*Np];
			C = sqrt(nx*nx+ny*ny+nz*nz);
			if (C==0.0) C=1.0;
			nx = nx/C;
			ny = ny/C;
			nz = nz/C;		

			// q=0
			fq = dist[n];
			rho = fq;
			m1  = -30.0*fq;
			m2  = 12.0*fq;

			// q=1
			fq = dist[2*Np+n];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jx = fq;
			m4 = -4.0*fq;
			m9 = 2.0*fq;
			m10 = -4.0*fq;

			// f2 = dist[10*Np+n];
			fq = dist[1*Np+n];
			rho += fq;
			m1 -= 11.0*(fq);
			m2 -= 4.0*(fq);
			jx -= fq;
			m4 += 4.0*(fq);
			m9 += 2.0*(fq);
			m10 -= 4.0*(fq);

			// q=3
			fq = dist[4*Np+n];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jy = fq;
			m6 = -4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 = fq;
			m12 = -2.0*fq;

			// q = 4
			fq = dist[3*Np+n];
			rho+= fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jy -= fq;
			m6 += 4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 += fq;
			m12 -= 2.0*fq;

			// q=5
			fq = dist[6*Np+n];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jz = fq;
			m8 = -4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 -= fq;
			m12 += 2.0*fq;

			// q = 6
			fq = dist[5*Np+n];
			rho+= fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jz -= fq;
			m8 += 4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 -= fq;
			m12 += 2.0*fq;

			// q=7
			fq = dist[8*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jy += fq;
			m6 += fq;
			m9  += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 = fq;
			m16 = fq;
			m17 = -fq;

			// q = 8
			fq = dist[7*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jy -= fq;
			m6 -= fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 += fq;
			m16 -= fq;
			m17 += fq;

			// q=9
			fq = dist[10*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jy -= fq;
			m6 -= fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 -= fq;
			m16 += fq;
			m17 += fq;

			// q = 10
			fq = dist[9*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jy += fq;
			m6 += fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 -= fq;
			m16 -= fq;
			m17 -= fq;

			// q=11
			fq = dist[12*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jz += fq;
			m8 += fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 = fq;
			m16 -= fq;
			m18 = fq;

			// q=12
			fq = dist[11*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jz -= fq;
			m8 -= fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 += fq;
			m16 += fq;
			m18 -= fq;

			// q=13
			fq = dist[14*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jz -= fq;
			m8 -= fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 -= fq;
			m16 -= fq;
			m18 -= fq;

			// q=14
			fq = dist[13*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jz += fq;
			m8 += fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 -= fq;
			m16 += fq;
			m18 += fq;

			// q=15
			fq = dist[16*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy += fq;
			m6 += fq;
			jz += fq;
			m8 += fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 = fq;
			m17 += fq;
			m18 -= fq;

			// q=16
			fq = dist[15*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy -= fq;
			m6 -= fq;
			jz -= fq;
			m8 -= fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 += fq;
			m17 -= fq;
			m18 += fq;

			// q=17
			fq = dist[18*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy += fq;
			m6 += fq;
			jz -= fq;
			m8 -= fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 -= fq;
			m17 += fq;
			m18 += fq;

			// q=18
			fq = dist[17*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy -= fq;
			m6 -= fq;
			jz += fq;
			m8 += fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 -= fq;
			m17 -= fq;
			m18 -= fq;

			//........................................................................
			//..............carry out relaxation process..............................
			//..........Toelke, Fruediger et. al. 2006................................
			if (C == 0.0)	nx = ny = nz = 0.0;
			m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho0 - 11*rho) -alpha*C - m1);
			m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho0)- m2);
			m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
			m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
			m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
			m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho0) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
			m10 = m10 + rlx_setA*( - m10);
			m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho0) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
			m12 = m12 + rlx_setA*( - m12);
			m13 = m13 + rlx_setA*( (jx*jy/rho0) + 0.5*alpha*C*nx*ny - m13);
			m14 = m14 + rlx_setA*( (jy*jz/rho0) + 0.5*alpha*C*ny*nz - m14);
			m15 = m15 + rlx_setA*( (jx*jz/rho0) + 0.5*alpha*C*nx*nz - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);

			
			//.......................................................................................................
			// assign force with wetting BC
			force_x = alpha*(nA-nB)*SolidForce[n] + Fx;
			force_y = alpha*(nA-nB)*SolidForce[n+Np] + Fy;
			force_z = alpha*(nA-nB)*SolidForce[n+2*Np] + Fz;
			//.................inverse transformation......................................................
			
			// q=0
			fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
			dist[n] = fq;

			// q = 1
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10) + 0.16666666*force_x;
			dist[1*Np+n] = fq;

			// q=2
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10) -  0.16666666*force_x;
			dist[2*Np+n] = fq;

			// q = 3
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) + 0.16666666*force_y;
			dist[3*Np+n] = fq;

			// q = 4
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) - 0.16666666*force_y;
			dist[4*Np+n] = fq;

			// q = 5
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) + 0.16666666*force_z;
			dist[5*Np+n] = fq;

			// q = 6
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) - 0.16666666*force_z;
			dist[6*Np+n] = fq;

			// q = 7
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17) + 0.08333333333*(force_x+force_y);
			dist[7*Np+n] = fq;


			// q = 8
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
					+mrt_V12*m12+0.25*m13+0.125*(m17-m16) - 0.08333333333*(force_x+force_y);
			dist[8*Np+n] = fq;

			// q = 9
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17) + 0.08333333333*(force_x-force_y);
			dist[9*Np+n] = fq;

			// q = 10
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17)- 0.08333333333*(force_x-force_y);
			dist[10*Np+n] = fq;


			// q = 11
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12+0.25*m15+0.125*(m18-m16) + 0.08333333333*(force_x+force_z);
			dist[11*Np+n] = fq;

			// q = 12
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
					mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18)-0.08333333333*(force_x+force_z);
			dist[12*Np+n] = fq;

			// q = 13
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15-0.125*(m16+m18) + 0.08333333333*(force_x-force_z);
			dist[13*Np+n] = fq;

			// q= 14
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15+0.125*(m16+m18) - 0.08333333333*(force_x-force_z);

			dist[14*Np+n] = fq;

			// q = 15
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18) + 0.08333333333*(force_y+force_z);
			dist[15*Np+n] = fq;

			// q = 16
			fq =  mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17)- 0.08333333333*(force_y+force_z);
			dist[16*Np+n] = fq;


			// q = 17
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18) + 0.08333333333*(force_y-force_z);
			dist[17*Np+n] = fq;

			// q = 18
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18) - 0.08333333333*(force_y-force_z);
			dist[18*Np+n] = fq;

			//........................................................................

			// write the velocity 
			ux = (jx + force_x) / rho0;
			uy = (jy + force_y) / rho0;
			uz = (jz + force_z) / rho0;
			//Velocity[n] = ux;
			//Velocity[Np+n] = uy;
			//Velocity[2*Np+n] = uz;

			// Instantiate mass transport distributions
			// Stationary value - distribution 0

			nAB = 1.0/(nA+nB);
			Aq[n] = 0.3333333333333333*nA;
			Bq[n] = 0.3333333333333333*nB;

			//...............................................
			// q = 0,2,4
			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
			delta = beta*nA*nB*nAB*0.1111111111111111*nx;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;

			Aq[1*Np+n] = a1;
			Bq[1*Np+n] = b1;
			Aq[2*Np+n] = a2;
			Bq[2*Np+n] = b2;

			//...............................................
			// q = 2
			// Cq = {0,1,0}
			delta = beta*nA*nB*nAB*0.1111111111111111*ny;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*uy))+delta;

			Aq[3*Np+n] = a1;
			Bq[3*Np+n] = b1;
			Aq[4*Np+n] = a2;
			Bq[4*Np+n] = b2;
			//...............................................
			// q = 4
			// Cq = {0,0,1}
			delta = beta*nA*nB*nAB*0.1111111111111111*nz;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*uz))+delta;

			Aq[5*Np+n] = a1;
			Bq[5*Np+n] = b1;
			Aq[6*Np+n] = a2;
			Bq[6*Np+n] = b2;
			//...............................................
		}
	}
}


__global__ void dvc_ScaLBL_D3Q19_AAodd_DFH(int *neighborList, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *Gradient, double *SolidForce, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int start, int finish, int Np){

	int n,nn,nread;
	int nr1,nr2,nr3,nr4,nr5,nr6;
	int nr7,nr8,nr9,nr10;
	int nr11,nr12,nr13,nr14;
	//int nr15,nr16,nr17,nr18;
	double fq;
	// conserved momemnts
	double rho,jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m3,m5,m7;
	double nA,nB; // number density
	double a1,b1,a2,b2,nAB,delta;
	double C,nx,ny,nz; //color gradient magnitude and direction
	double ux,uy,uz;
	double phi,tau,rho0,rlx_setA,rlx_setB;
	double force_x,force_y,force_z;

	const double mrt_V1=0.05263157894736842;
	const double mrt_V2=0.012531328320802;
	const double mrt_V3=0.04761904761904762;
	const double mrt_V4=0.004594820384294068;
	const double mrt_V5=0.01587301587301587;
	const double mrt_V6=0.0555555555555555555555555;
	const double mrt_V7=0.02777777777777778;
	const double mrt_V8=0.08333333333333333;
	const double mrt_V9=0.003341687552213868;
	const double mrt_V10=0.003968253968253968;
	const double mrt_V11=0.01388888888888889;
	const double mrt_V12=0.04166666666666666;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			// read the component number densities
			nA = Den[n];
			nB = Den[Np + n];

			// compute phase indicator field
			phi=(nA-nB)/(nA+nB);

			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
			// local relaxation time
			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
			rlx_setA = 1.f/tau;
			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
			
			//...........Read the Color Gradient.................................
			nx = Gradient[n];
			ny = Gradient[n+Np];
			nz = Gradient[n+2*Np];
			C = sqrt(nx*nx+ny*ny+nz*nz);
			if (C==0.0) C=1.0;
			nx = nx/C;
			ny = ny/C;
			nz = nz/C;			

			// q=0
			fq = dist[n];
			rho = fq;
			m1  = -30.0*fq;
			m2  = 12.0*fq;

			// q=1
			//nread = neighborList[n]; // neighbor 2 
			//fq = dist[nread]; // reading the f1 data into register fq		
			nr1 = neighborList[n]; 
			fq = dist[nr1]; // reading the f1 data into register fq
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jx = fq;
			m4 = -4.0*fq;
			m9 = 2.0*fq;
			m10 = -4.0*fq;

			// f2 = dist[10*Np+n];
			//nread = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
			//fq = dist[nread];  // reading the f2 data into register fq
			nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
			fq = dist[nr2];  // reading the f2 data into register fq
			rho += fq;
			m1 -= 11.0*(fq);
			m2 -= 4.0*(fq);
			jx -= fq;
			m4 += 4.0*(fq);
			m9 += 2.0*(fq);
			m10 -= 4.0*(fq);

			// q=3
			//nread = neighborList[n+2*Np]; // neighbor 4
			//fq = dist[nread];
			nr3 = neighborList[n+2*Np]; // neighbor 4
			fq = dist[nr3];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jy = fq;
			m6 = -4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 = fq;
			m12 = -2.0*fq;

			// q = 4
			//nread = neighborList[n+3*Np]; // neighbor 3
			//fq = dist[nread];
			nr4 = neighborList[n+3*Np]; // neighbor 3
			fq = dist[nr4];
			rho+= fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jy -= fq;
			m6 += 4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 += fq;
			m12 -= 2.0*fq;

			// q=5
			//nread = neighborList[n+4*Np];
			//fq = dist[nread];
			nr5 = neighborList[n+4*Np];
			fq = dist[nr5];
			rho += fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jz = fq;
			m8 = -4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 -= fq;
			m12 += 2.0*fq;


			// q = 6
			//nread = neighborList[n+5*Np];
			//fq = dist[nread];
			nr6 = neighborList[n+5*Np];
			fq = dist[nr6];
			rho+= fq;
			m1 -= 11.0*fq;
			m2 -= 4.0*fq;
			jz -= fq;
			m8 += 4.0*fq;
			m9 -= fq;
			m10 += 2.0*fq;
			m11 -= fq;
			m12 += 2.0*fq;

			// q=7
			//nread = neighborList[n+6*Np];
			//fq = dist[nread];
			nr7 = neighborList[n+6*Np];
			fq = dist[nr7];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jy += fq;
			m6 += fq;
			m9  += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 = fq;
			m16 = fq;
			m17 = -fq;

			// q = 8
			//nread = neighborList[n+7*Np];
			//fq = dist[nread];
			nr8 = neighborList[n+7*Np];
			fq = dist[nr8];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jy -= fq;
			m6 -= fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 += fq;
			m16 -= fq;
			m17 += fq;

			// q=9
			//nread = neighborList[n+8*Np];
			//fq = dist[nread];
			nr9 = neighborList[n+8*Np];
			fq = dist[nr9];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jy -= fq;
			m6 -= fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 -= fq;
			m16 += fq;
			m17 += fq;

			// q = 10
			//nread = neighborList[n+9*Np];
			//fq = dist[nread];
			nr10 = neighborList[n+9*Np];
			fq = dist[nr10];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jy += fq;
			m6 += fq;
			m9 += fq;
			m10 += fq;
			m11 += fq;
			m12 += fq;
			m13 -= fq;
			m16 -= fq;
			m17 -= fq;

			// q=11
			//nread = neighborList[n+10*Np];
			//fq = dist[nread];
			nr11 = neighborList[n+10*Np];
			fq = dist[nr11];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jz += fq;
			m8 += fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 = fq;
			m16 -= fq;
			m18 = fq;

			// q=12
			//nread = neighborList[n+11*Np];
			//fq = dist[nread];
			nr12 = neighborList[n+11*Np];
			fq = dist[nr12];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jz -= fq;
			m8 -= fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 += fq;
			m16 += fq;
			m18 -= fq;

			// q=13
			//nread = neighborList[n+12*Np];
			//fq = dist[nread];
			nr13 = neighborList[n+12*Np];
			fq = dist[nr13];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx += fq;
			m4 += fq;
			jz -= fq;
			m8 -= fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 -= fq;
			m16 -= fq;
			m18 -= fq;

			// q=14
			//nread = neighborList[n+13*Np];
			//fq = dist[nread];
			nr14 = neighborList[n+13*Np];
			fq = dist[nr14];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jx -= fq;
			m4 -= fq;
			jz += fq;
			m8 += fq;
			m9 += fq;
			m10 += fq;
			m11 -= fq;
			m12 -= fq;
			m15 -= fq;
			m16 += fq;
			m18 += fq;

			// q=15
			nread = neighborList[n+14*Np];
			fq = dist[nread];
			//fq = dist[17*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy += fq;
			m6 += fq;
			jz += fq;
			m8 += fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 = fq;
			m17 += fq;
			m18 -= fq;

			// q=16
			nread = neighborList[n+15*Np];
			fq = dist[nread];
			//fq = dist[8*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy -= fq;
			m6 -= fq;
			jz -= fq;
			m8 -= fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 += fq;
			m17 -= fq;
			m18 += fq;

			// q=17
			//fq = dist[18*Np+n];
			nread = neighborList[n+16*Np];
			fq = dist[nread];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy += fq;
			m6 += fq;
			jz -= fq;
			m8 -= fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 -= fq;
			m17 += fq;
			m18 += fq;

			// q=18
			nread = neighborList[n+17*Np];
			fq = dist[nread];
			//fq = dist[9*Np+n];
			rho += fq;
			m1 += 8.0*fq;
			m2 += fq;
			jy -= fq;
			m6 -= fq;
			jz += fq;
			m8 += fq;
			m9 -= 2.0*fq;
			m10 -= 2.0*fq;
			m14 -= fq;
			m17 -= fq;
			m18 -= fq;
			
			//........................................................................
			//..............carry out relaxation process..............................
			//..........Toelke, Fruediger et. al. 2006................................
			if (C == 0.0)	nx = ny = nz = 0.0;
			m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho0 - 11*rho) -alpha*C - m1);
			m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho0)- m2);
			m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
			m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
			m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
			m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho0) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
			m10 = m10 + rlx_setA*( - m10);
			m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho0) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
			m12 = m12 + rlx_setA*( - m12);
			m13 = m13 + rlx_setA*( (jx*jy/rho0) + 0.5*alpha*C*nx*ny - m13);
			m14 = m14 + rlx_setA*( (jy*jz/rho0) + 0.5*alpha*C*ny*nz - m14);
			m15 = m15 + rlx_setA*( (jx*jz/rho0) + 0.5*alpha*C*nx*nz - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);

			// assign force with wetting BC
			force_x = alpha*(nA-nB)*SolidForce[n] + Fx;
			force_y = alpha*(nA-nB)*SolidForce[n+Np] + Fy;
			force_z = alpha*(nA-nB)*SolidForce[n+2*Np] + Fz;
			
			//.................inverse transformation......................................................
			// q=0
			fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
			dist[n] = fq;

			// q = 1
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10)+0.16666666*force_x;
			//nread = neighborList[n+Np];
			dist[nr2] = fq;

			// q=2
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10) -  0.16666666*force_x;
			//nread = neighborList[n];
			dist[nr1] = fq;

			// q = 3
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) + 0.16666666*force_y;
			//nread = neighborList[n+3*Np];
			dist[nr4] = fq;

			// q = 4
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) - 0.16666666*force_y;
			//nread = neighborList[n+2*Np];
			dist[nr3] = fq;

			// q = 5
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) + 0.16666666*force_z;
			//nread = neighborList[n+5*Np];
			dist[nr6] = fq;

			// q = 6
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) - 0.16666666*force_z;
			//nread = neighborList[n+4*Np];
			dist[nr5] = fq;

			// q = 7
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17) + 0.08333333333*(force_x+force_y);
			//nread = neighborList[n+7*Np];
			dist[nr8] = fq;

			// q = 8
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
					+mrt_V12*m12+0.25*m13+0.125*(m17-m16) - 0.08333333333*(force_x+force_y);
			//nread = neighborList[n+6*Np];
			dist[nr7] = fq;

			// q = 9
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17) + 0.08333333333*(force_x-force_y);
			//nread = neighborList[n+9*Np];
			dist[nr10] = fq;

			// q = 10
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17)- 0.08333333333*(force_x-force_y);
			//nread = neighborList[n+8*Np];
			dist[nr9] = fq;

			// q = 11
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12+0.25*m15+0.125*(m18-m16) + 0.08333333333*(force_x+force_z);
			//nread = neighborList[n+11*Np];
			dist[nr12] = fq;

			// q = 12
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
					mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18) - 0.08333333333*(force_x+force_z);
			//nread = neighborList[n+10*Np];
			dist[nr11]= fq;

			// q = 13
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15-0.125*(m16+m18) + 0.08333333333*(force_x-force_z);
			//nread = neighborList[n+13*Np];
			dist[nr14] = fq;

			// q= 14
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15+0.125*(m16+m18) - 0.08333333333*(force_x-force_z);
			//nread = neighborList[n+12*Np];
			dist[nr13] = fq;


			// q = 15
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18) + 0.08333333333*(force_y+force_z);
			nread = neighborList[n+15*Np];
			dist[nread] = fq;

			// q = 16
			fq =  mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17)- 0.08333333333*(force_y+force_z);
			nread = neighborList[n+14*Np];
			dist[nread] = fq;


			// q = 17
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18) + 0.08333333333*(force_y-force_z);
			nread = neighborList[n+17*Np];
			dist[nread] = fq;

			// q = 18
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18) - 0.08333333333*(force_y-force_z);
			nread = neighborList[n+16*Np];
			dist[nread] = fq;

			// write the velocity 
			ux = (jx + force_x) / rho0;
			uy = (jy + force_y) / rho0;
			uz = (jz + force_z) / rho0;
			//Velocity[n] = ux;
			//Velocity[Np+n] = uy;
			//Velocity[2*Np+n] = uz;

			// Instantiate mass transport distributions
			// Stationary value - distribution 0
			nAB = 1.0/(nA+nB);
			Aq[n] = 0.3333333333333333*nA;
			Bq[n] = 0.3333333333333333*nB;

			//...............................................
			// q = 0,2,4
			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
			delta = beta*nA*nB*nAB*0.1111111111111111*nx;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;

			// q = 1
			//nread = neighborList[n+Np];
			Aq[nr2] = a1;
			Bq[nr2] = b1;
			// q=2
			//nread = neighborList[n];
			Aq[nr1] = a2;
			Bq[nr1] = b2;

			//...............................................
			// Cq = {0,1,0}
			delta = beta*nA*nB*nAB*0.1111111111111111*ny;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*uy))+delta;

			// q = 3
			//nread = neighborList[n+3*Np];
			Aq[nr4] = a1;
			Bq[nr4] = b1;
			// q = 4
			//nread = neighborList[n+2*Np];
			Aq[nr3] = a2;
			Bq[nr3] = b2;

			//...............................................
			// q = 4
			// Cq = {0,0,1}
			delta = beta*nA*nB*nAB*0.1111111111111111*nz;
			if (!(nA*nB*nAB>0)) delta=0;
			a1 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
			b1 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
			a2 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
			b2 = nB*(0.1111111111111111*(1-4.5*uz))+delta;

			// q = 5
			//nread = neighborList[n+5*Np];
			Aq[nr6] = a1;
			Bq[nr6] = b1;
			// q = 6
			//nread = neighborList[n+4*Np];
			Aq[nr5] = a2;
			Bq[nr5] = b2;
			//...............................................
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAodd_DFH(int *neighborList, double *Aq, double *Bq, 
		double *Den, double *Phi, int start, int finish, int Np){
	int n,nread;
	double fq,nA,nB;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			//..........Compute the number density for each component ............
			// q=0
			fq = Aq[n];
			nA = fq;
			fq = Bq[n];
			nB = fq;
			
			// q=1
			nread = neighborList[n]; 
			fq = Aq[nread];
			nA += fq;
			fq = Bq[nread]; 
			nB += fq;
			
			// q=2
			nread = neighborList[n+Np]; 
			fq = Aq[nread];  
			nA += fq;
			fq = Bq[nread]; 
			nB += fq;
			
			// q=3
			nread = neighborList[n+2*Np]; 
			fq = Aq[nread];
			nA += fq;
			fq = Bq[nread];
			nB += fq;
			
			// q = 4
			nread = neighborList[n+3*Np]; 
			fq = Aq[nread];
			nA += fq;
			fq = Bq[nread];
			nB += fq;

			// q=5
			nread = neighborList[n+4*Np];
			fq = Aq[nread];
			nA += fq;
			fq = Bq[nread];
			nB += fq;
			
			// q = 6
			nread = neighborList[n+5*Np];
			fq = Aq[nread];
			nA += fq;
			fq = Bq[nread];
			nB += fq;

			// save the number densities
			Den[n] = nA;
			Den[Np+n] = nB;

			// save the phase indicator field
			//idx = Map[n];
			Phi[n] = (nA-nB)/(nA+nB); 
		}
	}
}

__global__  void dvc_ScaLBL_D3Q7_AAeven_DFH(double *Aq, double *Bq, double *Den, double *Phi, 
		int start, int finish, int Np){
	int idx,n;
	double fq,nA,nB;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			// compute number density for each component
			// q=0
			fq = Aq[n];
			nA = fq;
			fq = Bq[n];
			nB = fq;
			
			// q=1
			fq = Aq[2*Np+n];
			nA += fq;
			fq = Bq[2*Np+n];
			nB += fq;

			// q=2
			fq = Aq[1*Np+n];
			nA += fq;
			fq = Bq[1*Np+n];
			nB += fq;

			// q=3
			fq = Aq[4*Np+n];
			nA += fq;
			fq = Bq[4*Np+n];
			nB += fq;

			// q = 4
			fq = Aq[3*Np+n];
			nA += fq;
			fq = Bq[3*Np+n];
			nB += fq;
			
			// q=5
			fq = Aq[6*Np+n];
			nA += fq;
			fq = Bq[6*Np+n];
			nB += fq;
			
			// q = 6
			fq = Aq[5*Np+n];
			nA += fq;
			fq = Bq[5*Np+n];
			nB += fq;

			// save the number densities
			Den[n] = nA;
			Den[Np+n] = nB;

			// save the phase indicator field
			//idx = Map[n];
			Phi[n] = (nA-nB)/(nA+nB); 	
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_Gradient_DFH(int *neighborList, double *Phi, double *ColorGrad, int start, int finish, int Np){

	int n,nn;
	// distributions
	double m1,m2,m3,m4,m5,m6,m7,m8,m9;
	double m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double nx,ny,nz;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			nn = neighborList[n+Np]%Np;
			m1 = Phi[nn];
			nn = neighborList[n]%Np;
			m2 = Phi[nn];
			nn = neighborList[n+3*Np]%Np;
			m3 = Phi[nn];
			nn = neighborList[n+2*Np]%Np;
			m4 = Phi[nn];		
			nn = neighborList[n+5*Np]%Np;
			m5 = Phi[nn];
			nn = neighborList[n+4*Np]%Np;
			m6 = Phi[nn];		
			nn = neighborList[n+7*Np]%Np;
			m7 = Phi[nn];
			nn = neighborList[n+6*Np]%Np;
			m8 = Phi[nn];		
			nn = neighborList[n+9*Np]%Np;
			m9 = Phi[nn];
			nn = neighborList[n+8*Np]%Np;
			m10 = Phi[nn];		
			nn = neighborList[n+11*Np]%Np;
			m11 = Phi[nn];
			nn = neighborList[n+10*Np]%Np;
			m12 = Phi[nn];		
			nn = neighborList[n+13*Np]%Np;
			m13 = Phi[nn];
			nn = neighborList[n+12*Np]%Np;
			m14 = Phi[nn];		
			nn = neighborList[n+15*Np]%Np;
			m15 = Phi[nn];
			nn = neighborList[n+14*Np]%Np;
			m16 = Phi[nn];		
			nn = neighborList[n+17*Np]%Np;
			m17 = Phi[nn];
			nn = neighborList[n+16*Np]%Np;
			m18 = Phi[nn];					
			
			//............Compute the Color Gradient...................................
			//............Compute the wn fluid Gradient...................................
			nx = (m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
			ny = (m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
			nz = (m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
			
	/*		// .... read the solid potential gradient.....................
			m1 = SolidPotential[n];
			m2 = SolidPotential[n+Np];
			m3 = SolidPotential[n+2*Np];
			nx += m1;
			ny += m2;
			nz += m3;
		*/	
			//...........Normalize the Color Gradient.................................
			//	C = sqrt(nx*nx+ny*ny+nz*nz);
			//	nx = nx/C;
			//	ny = ny/C;
			//	nz = nz/C;
			//...Store the Color Gradient....................
			ColorGrad[n] = nx;
			ColorGrad[Np+n] = ny;
			ColorGrad[2*Np+n] = nz;
			//...............................................
		}
	}
}

extern "C" void ScaLBL_Gradient_Unpack(double weight, double Cqx, double Cqy, double Cqz, 
		int *list, int start, int count, double *recvbuf, double *phi, double *grad, int N){
	int GRID = count / 512 + 1;
	dvc_ScaLBL_Gradient_Unpack<<<GRID,512 >>>(weight, Cqx, Cqy, Cqz, list, start, count, recvbuf, phi, grad, N);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_Gradient_Unpack: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_DFH_Init(double *Phi, double *Den, double *Aq, double *Bq, int start, int finish, int Np){
	dvc_ScaLBL_DFH_Init<<<NBLOCKS,NTHREADS >>>(Phi, Den, Aq, Bq, start, finish, Np); 
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_DFH_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_DFH(int *neighborList, double *dist, double *Aq, double *Bq, double *Den, double *Phi,
		double *Gradient, double *SolidForce, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int start, int finish, int Np){

	cudaProfilerStart();
	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAeven_DFH, cudaFuncCachePreferL1);

	dvc_ScaLBL_D3Q19_AAeven_DFH<<<NBLOCKS,NTHREADS >>>(neighborList, dist, Aq, Bq, Den, Phi, Gradient, SolidForce, rhoA, rhoB, tauA, tauB, 
			alpha, beta, Fx, Fy, Fz, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_DFH: %s \n",cudaGetErrorString(err));
	}
	cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q19_AAodd_DFH(int *neighborList, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *Gradient, double *SolidForce, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int start, int finish, int Np){

	cudaProfilerStart();
	cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAodd_DFH, cudaFuncCachePreferL1);
	
	dvc_ScaLBL_D3Q19_AAodd_DFH<<<NBLOCKS,NTHREADS >>>(neighborList,dist, Aq, Bq, Den, Phi, Gradient, 
			SolidForce, rhoA, rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz,  start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_DFH: %s \n",cudaGetErrorString(err));
	}
	cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAodd_DFH(int *NeighborList, double *Aq, double *Bq, 
		double *Den, double *Phi, int start, int finish, int Np){

	cudaProfilerStart();
	dvc_ScaLBL_D3Q7_AAodd_DFH<<<NBLOCKS,NTHREADS >>>(NeighborList, Aq, Bq, Den, Phi, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAodd_DFH: %s \n",cudaGetErrorString(err));
	}
	cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAeven_DFH(double *Aq, double *Bq, double *Den, double *Phi, 
		int start, int finish, int Np){

	cudaProfilerStart();
	dvc_ScaLBL_D3Q7_AAeven_DFH<<<NBLOCKS,NTHREADS >>>(Aq, Bq, Den, Phi, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q7_AAeven_DFH: %s \n",cudaGetErrorString(err));
	}
	cudaProfilerStop();

}

extern "C" void ScaLBL_D3Q19_Gradient_DFH(int *neighborList, double *Phi, double *ColorGrad, int start, int finish, int Np){

	dvc_ScaLBL_D3Q19_Gradient_DFH<<<NBLOCKS,NTHREADS >>>(neighborList, Phi, ColorGrad, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_Gradient_DFH: %s \n",cudaGetErrorString(err));
	}

}
