#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <stdio.h>
#include <math.h>

#define NBLOCKS 1024
#define NTHREADS 256

//Model-1 & 4
void dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor(int *neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den,
		 double *Phi, double *GreySolidGrad, double *Poros,double *Perm, double *Velocity, double *Pressure,
         double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff,double alpha, double beta,
		double Gx, double Gy, double Gz, int strideY, int strideZ, int start, int finish, int Np,
		sycl::nd_item<3> item_ct1){

	int n,nn,ijk,nread;
	int nr1,nr2,nr3,nr4,nr5,nr6;
	int nr7,nr8,nr9,nr10;
	int nr11,nr12,nr13,nr14;
	//int nr15,nr16,nr17,nr18;
	double fq;
	// conserved momemnts
	double rho,jx,jy,jz;
	double vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m3,m5,m7;
	double nA,nB; // number density
	double a1,b1,a2,b2,nAB,delta;
	double C,nx,ny,nz; //color gradient magnitude and direction
	double phi,tau,rho0,rlx_setA,rlx_setB;

    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double tau_eff;
    double mu_eff;//kinematic viscosity
    double nx_gs,ny_gs,nz_gs;//grey-solid color gradient
    double nx_phase,ny_phase,nz_phase,C_phase;
    double Fx,Fy,Fz;

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
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {
			// read the component number densities
			nA = Den[n];
			nB = Den[Np + n];

            porosity = Poros[n];
            perm = Perm[n];
            nx_gs = GreySolidGrad[n+0*Np];
            ny_gs = GreySolidGrad[n+1*Np];
            nz_gs = GreySolidGrad[n+2*Np];

			// compute phase indicator field
			phi=(nA-nB)/(nA+nB);

			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
			// local relaxation time
			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
			rlx_setA = 1.f/tau;
			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
            mu_eff = (tau_eff-0.5)/3.0;//kinematic viscosity
			
			// Get the 1D index based on regular data layout
			ijk = Map[n];
			//					COMPUTE THE COLOR GRADIENT
			//........................................................................
			//.................Read Phase Indicator Values............................
			//........................................................................
			nn = ijk-1;							// neighbor index (get convention)
			m1 = Phi[nn];						// get neighbor for phi - 1
			//........................................................................
			nn = ijk+1;							// neighbor index (get convention)
			m2 = Phi[nn];						// get neighbor for phi - 2
			//........................................................................
			nn = ijk-strideY;							// neighbor index (get convention)
			m3 = Phi[nn];					// get neighbor for phi - 3
			//........................................................................
			nn = ijk+strideY;							// neighbor index (get convention)
			m4 = Phi[nn];					// get neighbor for phi - 4
			//........................................................................
			nn = ijk-strideZ;						// neighbor index (get convention)
			m5 = Phi[nn];					// get neighbor for phi - 5
			//........................................................................
			nn = ijk+strideZ;						// neighbor index (get convention)
			m6 = Phi[nn];					// get neighbor for phi - 6
			//........................................................................
			nn = ijk-strideY-1;						// neighbor index (get convention)
			m7 = Phi[nn];					// get neighbor for phi - 7
			//........................................................................
			nn = ijk+strideY+1;						// neighbor index (get convention)
			m8 = Phi[nn];					// get neighbor for phi - 8
			//........................................................................
			nn = ijk+strideY-1;						// neighbor index (get convention)
			m9 = Phi[nn];					// get neighbor for phi - 9
			//........................................................................
			nn = ijk-strideY+1;						// neighbor index (get convention)
			m10 = Phi[nn];					// get neighbor for phi - 10
			//........................................................................
			nn = ijk-strideZ-1;						// neighbor index (get convention)
			m11 = Phi[nn];					// get neighbor for phi - 11
			//........................................................................
			nn = ijk+strideZ+1;						// neighbor index (get convention)
			m12 = Phi[nn];					// get neighbor for phi - 12
			//........................................................................
			nn = ijk+strideZ-1;						// neighbor index (get convention)
			m13 = Phi[nn];					// get neighbor for phi - 13
			//........................................................................
			nn = ijk-strideZ+1;						// neighbor index (get convention)
			m14 = Phi[nn];					// get neighbor for phi - 14
			//........................................................................
			nn = ijk-strideZ-strideY;					// neighbor index (get convention)
			m15 = Phi[nn];					// get neighbor for phi - 15
			//........................................................................
			nn = ijk+strideZ+strideY;					// neighbor index (get convention)
			m16 = Phi[nn];					// get neighbor for phi - 16
			//........................................................................
			nn = ijk+strideZ-strideY;					// neighbor index (get convention)
			m17 = Phi[nn];					// get neighbor for phi - 17
			//........................................................................
			nn = ijk-strideZ+strideY;					// neighbor index (get convention)
			m18 = Phi[nn];					// get neighbor for phi - 18
			//............Compute the Color Gradient...................................
			nx_phase = -(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
			ny_phase = -(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
			nz_phase = -(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
                        C_phase = sycl::sqrt(nx_phase * nx_phase +
                                             ny_phase * ny_phase +
                                             nz_phase * nz_phase);

            //correct the normal color gradient by considering the effect of grey solid
            nx = nx_phase + (1.0-porosity)*nx_gs; 
            ny = ny_phase + (1.0-porosity)*ny_gs; 
            nz = nz_phase + (1.0-porosity)*nz_gs; 
            if (C_phase==0.0){
                nx = nx_phase; 
                ny = ny_phase;
                nz = nz_phase;
            }

			//...........Normalize the Color Gradient.................................
                        C = sycl::sqrt(nx * nx + ny * ny + nz * nz);
                        double ColorMag = C;
			if (C==0.0) ColorMag=1.0;
			nx = nx/ColorMag;
			ny = ny/ColorMag;
			nz = nz/ColorMag;		

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
			
            // Compute greyscale related parameters
            c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
            c1 = porosity * 0.5 * GeoFun / sycl::sqrt(perm);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jx/rho0+0.5*(porosity*Gx);
            vy = jy/rho0+0.5*(porosity*Gy);
            vz = jz/rho0+0.5*(porosity*Gz);
            v_mag = sycl::sqrt(vx * vx + vy * vy + vz * vz);
            ux = vx / (c0 + sycl::sqrt(c0 * c0 + c1 * v_mag));
            uy = vy / (c0 + sycl::sqrt(c0 * c0 + c1 * v_mag));
            uz = vz / (c0 + sycl::sqrt(c0 * c0 + c1 * v_mag));
            u_mag = sycl::sqrt(ux * ux + uy * uy + uz * uz);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            Fx = rho0 * (-porosity * mu_eff / perm * ux -
                         porosity * GeoFun / sycl::sqrt(perm) * u_mag * ux +
                         porosity * Gx);
            Fy = rho0 * (-porosity * mu_eff / perm * uy -
                         porosity * GeoFun / sycl::sqrt(perm) * u_mag * uy +
                         porosity * Gy);
            Fz = rho0 * (-porosity * mu_eff / perm * uz -
                         porosity * GeoFun / sycl::sqrt(perm) * u_mag * uz +
                         porosity * Gz);
            if (porosity==1.0){
                Fx=rho0*(Gx);
                Fy=rho0*(Gy);
                Fz=rho0*(Gz);
            }

			// write the velocity 
			Velocity[n] = ux;
			Velocity[Np+n] = uy;
			Velocity[2*Np+n] = uz;
            //Pressure[n] = rho/3.f/porosity;
            Pressure[n] = rho/3.f;

			//........................................................................
			//..............carry out relaxation process..............................
			//..........Toelke, Fruediger et. al. 2006................................
            //---------------- NO higher-order force -------------------------------//
			if (C == 0.0)	nx = ny = nz = 0.0;
			m1 = m1 + rlx_setA*((19*(ux*ux+uy*uy+uz*uz)*rho0/porosity - 11*rho) -19*alpha*C - m1);
			m2 = m2 + rlx_setA*((3*rho - 5.5*(ux*ux+uy*uy+uz*uz)*rho0/porosity)- m2);
            jx = jx + Fx;
			m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0)- m4)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
            jy = jy + Fy;
			m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0)- m6)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
            jz = jz + Fz;
			m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0)- m8)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
			m9 = m9 + rlx_setA*(((2*ux*ux-uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
			m10 = m10 + rlx_setA*( - m10);
            //m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
			m11 = m11 + rlx_setA*(((uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
			m12 = m12 + rlx_setA*( - m12);
            //m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
			m13 = m13 + rlx_setA*( (ux*uy*rho0/porosity) + 0.5*alpha*C*nx*ny - m13);
			m14 = m14 + rlx_setA*( (uy*uz*rho0/porosity) + 0.5*alpha*C*ny*nz - m14);
			m15 = m15 + rlx_setA*( (ux*uz*rho0/porosity) + 0.5*alpha*C*nx*nz - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);
            //----------------------------------------------------------------------//

            //----------------With higher-order force ------------------------------//
			//if (C == 0.0)	nx = ny = nz = 0.0;
			//m1 = m1 + rlx_setA*((19*(ux*ux+uy*uy+uz*uz)*rho0/porosity - 11*rho) -19*alpha*C - m1)
            //        + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
			//m2 = m2 + rlx_setA*((3*rho - 5.5*(ux*ux+uy*uy+uz*uz)*rho0/porosity)- m2)
            //        + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
            //jx = jx + Fx;
			//m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0)- m4)
            //        + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
            //jy = jy + Fy;
			//m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0)- m6)
            //        + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
            //jz = jz + Fz;
			//m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0)- m8)
            //        + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
			//m9 = m9 + rlx_setA*(((2*ux*ux-uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9)
            //        + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
			////m10 = m10 + rlx_setA*( - m10);
            //m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
            //          + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
			//m11 = m11 + rlx_setA*(((uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(ny*ny-nz*nz)- m11)
            //          + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
			////m12 = m12 + rlx_setA*( - m12);
            //m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12)
            //          + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
			//m13 = m13 + rlx_setA*( (ux*uy*rho0/porosity) + 0.5*alpha*C*nx*ny - m13);
            //          + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
			//m14 = m14 + rlx_setA*( (uy*uz*rho0/porosity) + 0.5*alpha*C*ny*nz - m14);
            //          + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
			//m15 = m15 + rlx_setA*( (ux*uz*rho0/porosity) + 0.5*alpha*C*nx*nz - m15);
            //          + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
			//m16 = m16 + rlx_setB*( - m16);
			//m17 = m17 + rlx_setB*( - m17);
			//m18 = m18 + rlx_setB*( - m18);
            //----------------------------------------------------------------------//

			//.................inverse transformation......................................................
			// q=0
			fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
			dist[n] = fq;

			// q = 1
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
			//nread = neighborList[n+Np];
			dist[nr2] = fq;

			// q=2
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
			//nread = neighborList[n];
			dist[nr1] = fq;

			// q = 3
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
			//nread = neighborList[n+3*Np];
			dist[nr4] = fq;

			// q = 4
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
			//nread = neighborList[n+2*Np];
			dist[nr3] = fq;

			// q = 5
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
			//nread = neighborList[n+5*Np];
			dist[nr6] = fq;

			// q = 6
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
			//nread = neighborList[n+4*Np];
			dist[nr5] = fq;

			// q = 7
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
			//nread = neighborList[n+7*Np];
			dist[nr8] = fq;

			// q = 8
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
					+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
			//nread = neighborList[n+6*Np];
			dist[nr7] = fq;

			// q = 9
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
			//nread = neighborList[n+9*Np];
			dist[nr10] = fq;

			// q = 10
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
			//nread = neighborList[n+8*Np];
			dist[nr9] = fq;

			// q = 11
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
			//nread = neighborList[n+11*Np];
			dist[nr12] = fq;

			// q = 12
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
					mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
			//nread = neighborList[n+10*Np];
			dist[nr11]= fq;

			// q = 13
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
			//nread = neighborList[n+13*Np];
			dist[nr14] = fq;

			// q= 14
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
			//nread = neighborList[n+12*Np];
			dist[nr13] = fq;


			// q = 15
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
			nread = neighborList[n+15*Np];
			dist[nread] = fq;

			// q = 16
			fq =  mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
			nread = neighborList[n+14*Np];
			dist[nread] = fq;


			// q = 17
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
			nread = neighborList[n+17*Np];
			dist[nread] = fq;

			// q = 18
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
			nread = neighborList[n+16*Np];
			dist[nread] = fq;
			//........................................................................

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

//Model-1 & 4
void dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
        double *Phi, double *GreySolidGrad, double *Poros,double *Perm, double *Velocity, double *Pressure, 
        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
		double Gx, double Gy, double Gz, int strideY, int strideZ, int start, int finish, int Np,
		sycl::nd_item<3> item_ct1){
	int ijk,nn,n;
	double fq;
	// conserved momemnts
	double rho,jx,jy,jz;
	double vx,vy,vz,v_mag;
    double ux,uy,uz,u_mag;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m3,m5,m7;
	double nA,nB; // number density
	double a1,b1,a2,b2,nAB,delta;
	double C,nx,ny,nz; //color gradient magnitude and direction
	double phi,tau,rho0,rlx_setA,rlx_setB;

    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double c0, c1; //Guo's model parameters
    double tau_eff;
    double mu_eff;//kinematic viscosity
    double nx_gs,ny_gs,nz_gs;//grey-solid color gradient
    double nx_phase,ny_phase,nz_phase,C_phase;
    double Fx,Fy,Fz;

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
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {

			// read the component number densities
			nA = Den[n];
			nB = Den[Np + n];
            porosity = Poros[n];
            perm = Perm[n];
            nx_gs = GreySolidGrad[n+0*Np];
            ny_gs = GreySolidGrad[n+1*Np];
            nz_gs = GreySolidGrad[n+2*Np];

			// compute phase indicator field
			phi=(nA-nB)/(nA+nB);

			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
			// local relaxation time
			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
			rlx_setA = 1.f/tau;
			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
            mu_eff = (tau_eff-0.5)/3.0;//kinematic viscosity

			// Get the 1D index based on regular data layout
			ijk = Map[n];
			//					COMPUTE THE COLOR GRADIENT
			//........................................................................
			//.................Read Phase Indicator Values............................
			//........................................................................
			nn = ijk-1;							// neighbor index (get convention)
			m1 = Phi[nn];						// get neighbor for phi - 1
			//........................................................................
			nn = ijk+1;							// neighbor index (get convention)
			m2 = Phi[nn];						// get neighbor for phi - 2
			//........................................................................
			nn = ijk-strideY;							// neighbor index (get convention)
			m3 = Phi[nn];					// get neighbor for phi - 3
			//........................................................................
			nn = ijk+strideY;							// neighbor index (get convention)
			m4 = Phi[nn];					// get neighbor for phi - 4
			//........................................................................
			nn = ijk-strideZ;						// neighbor index (get convention)
			m5 = Phi[nn];					// get neighbor for phi - 5
			//........................................................................
			nn = ijk+strideZ;						// neighbor index (get convention)
			m6 = Phi[nn];					// get neighbor for phi - 6
			//........................................................................
			nn = ijk-strideY-1;						// neighbor index (get convention)
			m7 = Phi[nn];					// get neighbor for phi - 7
			//........................................................................
			nn = ijk+strideY+1;						// neighbor index (get convention)
			m8 = Phi[nn];					// get neighbor for phi - 8
			//........................................................................
			nn = ijk+strideY-1;						// neighbor index (get convention)
			m9 = Phi[nn];					// get neighbor for phi - 9
			//........................................................................
			nn = ijk-strideY+1;						// neighbor index (get convention)
			m10 = Phi[nn];					// get neighbor for phi - 10
			//........................................................................
			nn = ijk-strideZ-1;						// neighbor index (get convention)
			m11 = Phi[nn];					// get neighbor for phi - 11
			//........................................................................
			nn = ijk+strideZ+1;						// neighbor index (get convention)
			m12 = Phi[nn];					// get neighbor for phi - 12
			//........................................................................
			nn = ijk+strideZ-1;						// neighbor index (get convention)
			m13 = Phi[nn];					// get neighbor for phi - 13
			//........................................................................
			nn = ijk-strideZ+1;						// neighbor index (get convention)
			m14 = Phi[nn];					// get neighbor for phi - 14
			//........................................................................
			nn = ijk-strideZ-strideY;					// neighbor index (get convention)
			m15 = Phi[nn];					// get neighbor for phi - 15
			//........................................................................
			nn = ijk+strideZ+strideY;					// neighbor index (get convention)
			m16 = Phi[nn];					// get neighbor for phi - 16
			//........................................................................
			nn = ijk+strideZ-strideY;					// neighbor index (get convention)
			m17 = Phi[nn];					// get neighbor for phi - 17
			//........................................................................
			nn = ijk-strideZ+strideY;					// neighbor index (get convention)
			m18 = Phi[nn];					// get neighbor for phi - 18
			//............Compute the Color Gradient...................................
			nx_phase = -(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
			ny_phase = -(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
			nz_phase = -(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
                        C_phase = sycl::sqrt(nx_phase * nx_phase +
                                             ny_phase * ny_phase +
                                             nz_phase * nz_phase);

            //correct the normal color gradient by considering the effect of grey solid
            nx = nx_phase + (1.0-porosity)*nx_gs; 
            ny = ny_phase + (1.0-porosity)*ny_gs; 
            nz = nz_phase + (1.0-porosity)*nz_gs; 
            if (C_phase==0.0){
                nx = nx_phase; 
                ny = ny_phase;
                nz = nz_phase;
            }

			//...........Normalize the Color Gradient.................................
                        C = sycl::sqrt(nx * nx + ny * ny + nz * nz);
                        double ColorMag = C;
			if (C==0.0) ColorMag=1.0;
			nx = nx/ColorMag;
			ny = ny/ColorMag;
			nz = nz/ColorMag;		

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

            // Compute greyscale related parameters
            c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
            c1 = porosity * 0.5 * GeoFun / sycl::sqrt(perm);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jx/rho0+0.5*(porosity*Gx);
            vy = jy/rho0+0.5*(porosity*Gy);
            vz = jz/rho0+0.5*(porosity*Gz);
            v_mag = sycl::sqrt(vx * vx + vy * vy + vz * vz);
            ux = vx / (c0 + sycl::sqrt(c0 * c0 + c1 * v_mag));
            uy = vy / (c0 + sycl::sqrt(c0 * c0 + c1 * v_mag));
            uz = vz / (c0 + sycl::sqrt(c0 * c0 + c1 * v_mag));
            u_mag = sycl::sqrt(ux * ux + uy * uy + uz * uz);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            Fx = rho0 * (-porosity * mu_eff / perm * ux -
                         porosity * GeoFun / sycl::sqrt(perm) * u_mag * ux +
                         porosity * Gx);
            Fy = rho0 * (-porosity * mu_eff / perm * uy -
                         porosity * GeoFun / sycl::sqrt(perm) * u_mag * uy +
                         porosity * Gy);
            Fz = rho0 * (-porosity * mu_eff / perm * uz -
                         porosity * GeoFun / sycl::sqrt(perm) * u_mag * uz +
                         porosity * Gz);
            if (porosity==1.0){
                Fx=rho0*(Gx);
                Fy=rho0*(Gy);
                Fz=rho0*(Gz);
            }

			// write the velocity 
			Velocity[n] = ux;
			Velocity[Np+n] = uy;
			Velocity[2*Np+n] = uz;
            //Pressure[n] = rho/3.f/porosity;
            Pressure[n] = rho/3.f;

			//........................................................................
			//..............carry out relaxation process..............................
			//..........Toelke, Fruediger et. al. 2006................................
            //---------------- NO higher-order force -------------------------------//
			if (C == 0.0)	nx = ny = nz = 0.0;
			m1 = m1 + rlx_setA*((19*(ux*ux+uy*uy+uz*uz)*rho0/porosity - 11*rho) -19*alpha*C - m1);
			m2 = m2 + rlx_setA*((3*rho - 5.5*(ux*ux+uy*uy+uz*uz)*rho0/porosity)- m2);
            jx = jx + Fx;
			m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0)- m4)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
            jy = jy + Fy;
			m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0)- m6)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
            jz = jz + Fz;
			m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0)- m8)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
			m9 = m9 + rlx_setA*(((2*ux*ux-uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
			m10 = m10 + rlx_setA*( - m10);
            //m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
			m11 = m11 + rlx_setA*(((uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
			m12 = m12 + rlx_setA*( - m12);
            //m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
			m13 = m13 + rlx_setA*( (ux*uy*rho0/porosity) + 0.5*alpha*C*nx*ny - m13);
			m14 = m14 + rlx_setA*( (uy*uz*rho0/porosity) + 0.5*alpha*C*ny*nz - m14);
			m15 = m15 + rlx_setA*( (ux*uz*rho0/porosity) + 0.5*alpha*C*nx*nz - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);
            //----------------------------------------------------------------------//

            //----------------With higher-order force ------------------------------//
			//if (C == 0.0)	nx = ny = nz = 0.0;
			//m1 = m1 + rlx_setA*((19*(ux*ux+uy*uy+uz*uz)*rho0/porosity - 11*rho) -19*alpha*C - m1)
            //        + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
			//m2 = m2 + rlx_setA*((3*rho - 5.5*(ux*ux+uy*uy+uz*uz)*rho0/porosity)- m2)
            //        + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
            //jx = jx + Fx;
			//m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0)- m4)
            //        + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
            //jy = jy + Fy;
			//m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0)- m6)
            //        + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
            //jz = jz + Fz;
			//m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0)- m8)
            //        + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
			//m9 = m9 + rlx_setA*(((2*ux*ux-uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9)
            //        + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
			////m10 = m10 + rlx_setA*( - m10);
            //m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
            //          + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
			//m11 = m11 + rlx_setA*(((uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(ny*ny-nz*nz)- m11)
            //          + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
			////m12 = m12 + rlx_setA*( - m12);
            //m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12)
            //          + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
			//m13 = m13 + rlx_setA*( (ux*uy*rho0/porosity) + 0.5*alpha*C*nx*ny - m13);
            //          + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
			//m14 = m14 + rlx_setA*( (uy*uz*rho0/porosity) + 0.5*alpha*C*ny*nz - m14);
            //          + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
			//m15 = m15 + rlx_setA*( (ux*uz*rho0/porosity) + 0.5*alpha*C*nx*nz - m15);
            //          + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
			//m16 = m16 + rlx_setB*( - m16);
			//m17 = m17 + rlx_setB*( - m17);
			//m18 = m18 + rlx_setB*( - m18);
            //----------------------------------------------------------------------//

			//.................inverse transformation......................................................
			// q=0
			fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
			dist[n] = fq;

			// q = 1
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
			dist[1*Np+n] = fq;

			// q=2
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
			dist[2*Np+n] = fq;

			// q = 3
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
			dist[3*Np+n] = fq;

			// q = 4
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
			dist[4*Np+n] = fq;

			// q = 5
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
			dist[5*Np+n] = fq;

			// q = 6
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
			dist[6*Np+n] = fq;

			// q = 7
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
			dist[7*Np+n] = fq;


			// q = 8
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
					+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
			dist[8*Np+n] = fq;

			// q = 9
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
			dist[9*Np+n] = fq;

			// q = 10
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
			dist[10*Np+n] = fq;


			// q = 11
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
			dist[11*Np+n] = fq;

			// q = 12
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
					mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
			dist[12*Np+n] = fq;

			// q = 13
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
			dist[13*Np+n] = fq;

			// q= 14
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15+0.125*(m16+m18);

			dist[14*Np+n] = fq;

			// q = 15
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
			dist[15*Np+n] = fq;

			// q = 16
			fq =  mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
			dist[16*Np+n] = fq;


			// q = 17
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
			dist[17*Np+n] = fq;

			// q = 18
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
			dist[18*Np+n] = fq;
			//........................................................................

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

//CP: capillary penalty
// also turn off recoloring for grey nodes
void dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor_CP(int *neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den,
		 double *Phi, double *GreySolidW, double *GreySn, double *GreySw, double *GreyKn, double *GreyKw, double *Poros,
		 double *Perm, double *Velocity, double *MobilityRatio, double *Pressure,
         double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff,double alpha, double beta,
		double Gx, double Gy, double Gz, bool RecoloringOff, int strideY, int strideZ, int start, int finish, int Np,
		sycl::nd_item<3> item_ct1){

	int n,nn,ijk,nread;
	int nr1,nr2,nr3,nr4,nr5,nr6;
	int nr7,nr8,nr9,nr10;
	int nr11,nr12,nr13,nr14;
	//int nr15,nr16,nr17,nr18;
	double fq;
	// conserved momemnts
	double rho,jx,jy,jz;
    double ux,uy,uz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m3,m5,m7;
	double nA,nB; // number density
	double a1,b1,a2,b2,nAB,delta;
	double C,nx,ny,nz; //color gradient magnitude and direction
	double phi,tau,rho0,rlx_setA,rlx_setB;

    double porosity;
    double perm;//voxel permeability
    double tau_eff;
    double mu_eff;//kinematic viscosity
    double Fx,Fy,Fz;
    double Fcpx,Fcpy,Fcpz;//capillary penalty force
    double W;//greyscale wetting strength
    double Sn_grey,Sw_grey;
    double GreyDiff=0.0e-4;
    
    /* Corey model parameters */
    double Kn_grey,Kw_grey;    
    double Swn,Krn_grey,Krw_grey,mobility_ratio,jA,jB;
    
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
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {
			// read the component number densities
			nA = Den[n];
			nB = Den[Np + n];

	        porosity = Poros[n];
	        //GreyDiff = Perm[n];
	        perm = 1.0;
	        W = GreySolidW[n];
	        Sn_grey = GreySn[n];
	        Sw_grey = GreySw[n];
	        Kn_grey = GreyKn[n];
	        Kw_grey = GreyKw[n];
	        
			// compute phase indicator field
			phi=(nA-nB)/(nA+nB);
			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);

			// local relaxation time
			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
			rlx_setA = 1.f/tau;
			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
	        mu_eff = (tau_eff-0.5)/3.0;//kinematic viscosity
			
	        mobility_ratio = 1.0;
	        Krn_grey = 0.0;
	        Krw_grey = 0.0;
	        if (nA/(nA+nB)<Sn_grey && porosity !=1.0){
	        	perm = Kw_grey;
	        	Krw_grey = Kw_grey;
	        	Swn = 0.0;
	        }
	        else if (nA/(nA+nB)>=Sn_grey && nA/(nA+nB) <= Sw_grey && porosity !=1.0){ 
	        	Swn = (nA/(nA+nB) - Sn_grey) /(Sw_grey - Sn_grey);
	        	Krn_grey = Kn_grey*Swn*Swn; // Corey model with exponent = 2, make sure that W cannot shift to zero
	        	Krw_grey = Kw_grey*(1.0-Swn)*(1.0-Swn); // Corey model with exponent = 4, make sure that W cannot shift to zero
	        	// recompute the effective permeability
	        	perm = mu_eff*(Krn_grey*3.0/(tauA-0.5) + Krw_grey*3.0/(tauB-0.5));
	        	//mobility_ratio =(nA*Krn_grey*3.0/(tauA-0.5) - nB*Krw_grey*3.0/(tauB-0.5))/(nA*Krn_grey*3.0/(tauA-0.5) + nB*Krw_grey*3.0/(tauB-0.5));
	        }
	        else if (nA/(nA+nB)>Sw_grey && porosity !=1.0){
	        	perm = Kn_grey;
	        	Krn_grey = Kn_grey;
	        	Swn = 1.0;
	        }
	        /* compute the mobility ratio */
	        if (porosity != 1.0){
	        	mobility_ratio =(Krn_grey/(tauA-0.5) - Krw_grey/(tauB-0.5))/(Krn_grey/(tauA-0.5) + Krw_grey/(tauB-0.5));
	        }
	        else if (phi > 0.0){
	        	mobility_ratio = 1.0;
	        }
	        else {
	        	mobility_ratio = -1.0;
	        }
	        MobilityRatio[n] = mobility_ratio;

			// Get the 1D index based on regular data layout
			ijk = Map[n];
			//					COMPUTE THE COLOR GRADIENT
			//........................................................................
			//.................Read Phase Indicator Values............................
			//........................................................................
			nn = ijk-1;							// neighbor index (get convention)
			m1 = Phi[nn];						// get neighbor for phi - 1
			//........................................................................
			nn = ijk+1;							// neighbor index (get convention)
			m2 = Phi[nn];						// get neighbor for phi - 2
			//........................................................................
			nn = ijk-strideY;							// neighbor index (get convention)
			m3 = Phi[nn];					// get neighbor for phi - 3
			//........................................................................
			nn = ijk+strideY;							// neighbor index (get convention)
			m4 = Phi[nn];					// get neighbor for phi - 4
			//........................................................................
			nn = ijk-strideZ;						// neighbor index (get convention)
			m5 = Phi[nn];					// get neighbor for phi - 5
			//........................................................................
			nn = ijk+strideZ;						// neighbor index (get convention)
			m6 = Phi[nn];					// get neighbor for phi - 6
			//........................................................................
			nn = ijk-strideY-1;						// neighbor index (get convention)
			m7 = Phi[nn];					// get neighbor for phi - 7
			//........................................................................
			nn = ijk+strideY+1;						// neighbor index (get convention)
			m8 = Phi[nn];					// get neighbor for phi - 8
			//........................................................................
			nn = ijk+strideY-1;						// neighbor index (get convention)
			m9 = Phi[nn];					// get neighbor for phi - 9
			//........................................................................
			nn = ijk-strideY+1;						// neighbor index (get convention)
			m10 = Phi[nn];					// get neighbor for phi - 10
			//........................................................................
			nn = ijk-strideZ-1;						// neighbor index (get convention)
			m11 = Phi[nn];					// get neighbor for phi - 11
			//........................................................................
			nn = ijk+strideZ+1;						// neighbor index (get convention)
			m12 = Phi[nn];					// get neighbor for phi - 12
			//........................................................................
			nn = ijk+strideZ-1;						// neighbor index (get convention)
			m13 = Phi[nn];					// get neighbor for phi - 13
			//........................................................................
			nn = ijk-strideZ+1;						// neighbor index (get convention)
			m14 = Phi[nn];					// get neighbor for phi - 14
			//........................................................................
			nn = ijk-strideZ-strideY;					// neighbor index (get convention)
			m15 = Phi[nn];					// get neighbor for phi - 15
			//........................................................................
			nn = ijk+strideZ+strideY;					// neighbor index (get convention)
			m16 = Phi[nn];					// get neighbor for phi - 16
			//........................................................................
			nn = ijk+strideZ-strideY;					// neighbor index (get convention)
			m17 = Phi[nn];					// get neighbor for phi - 17
			//........................................................................
			nn = ijk-strideZ+strideY;					// neighbor index (get convention)
			m18 = Phi[nn];					// get neighbor for phi - 18
			//............Compute the Color Gradient...................................
			nx = -3.0/18.0*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
			ny = -3.0/18.0*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
			nz = -3.0/18.0*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));

	        Fcpx = nx;
	        Fcpy = ny;
	        Fcpz = nz;
                double Fcp_mag =
                    sycl::sqrt(Fcpx * Fcpx + Fcpy * Fcpy + Fcpz * Fcpz);
                if (Fcp_mag==0.0) Fcpx=Fcpy=Fcpz=0.0;
	        //NOTE for open node (porosity=1.0),Fcp=0.0
                Fcpx *= alpha * W * (1.0 - porosity) / sycl::sqrt(perm);
                Fcpy *= alpha * W * (1.0 - porosity) / sycl::sqrt(perm);
                Fcpz *= alpha * W * (1.0 - porosity) / sycl::sqrt(perm);

                        //...........Normalize the Color Gradient.................................
                        C = sycl::sqrt(nx * nx + ny * ny + nz * nz);
                        double ColorMag = C;
			if (C==0.0) ColorMag=1.0;
			nx = nx/ColorMag;
			ny = ny/ColorMag;
			nz = nz/ColorMag;		

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
			
	        // Compute greyscale related parameters
	        ux = (jx/rho0+0.5*porosity*Gx+0.5*Fcpx/rho0)/(1.0+0.5*porosity*mu_eff/perm);
	        uy = (jy/rho0+0.5*porosity*Gy+0.5*Fcpy/rho0)/(1.0+0.5*porosity*mu_eff/perm);
	        uz = (jz/rho0+0.5*porosity*Gz+0.5*Fcpz/rho0)/(1.0+0.5*porosity*mu_eff/perm);
	        if (porosity==1.0){//i.e. open nodes
	            ux = (jx/rho0+0.5*porosity*Gx);
	            uy = (jy/rho0+0.5*porosity*Gy);
	            uz = (jz/rho0+0.5*porosity*Gz);
	        }

	        //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
	        Fx = rho0*(-porosity*mu_eff/perm*ux + porosity*Gx)+Fcpx;
	        Fy = rho0*(-porosity*mu_eff/perm*uy + porosity*Gy)+Fcpy;
	        Fz = rho0*(-porosity*mu_eff/perm*uz + porosity*Gz)+Fcpz;
	        if (porosity==1.0){
	            Fx=rho0*(porosity*Gx);
	            Fy=rho0*(porosity*Gy);
	            Fz=rho0*(porosity*Gz);
	        }

			// write the velocity 
			Velocity[n] = ux;
			Velocity[Np+n] = uy;
			Velocity[2*Np+n] = uz;
	        //Pressure[n] = rho/3.f/porosity;
	        Pressure[n] = rho/3.f;

			//........................................................................
			//..............carry out relaxation process..............................
			//..........Toelke, Fruediger et. al. 2006................................
	        //---------------- NO higher-order force -------------------------------//
			if (C == 0.0)	nx = ny = nz = 0.0;
			m1 = m1 + rlx_setA*((19*(ux*ux+uy*uy+uz*uz)*rho0/porosity - 11*rho) -19*alpha*C - m1);
			m2 = m2 + rlx_setA*((3*rho - 5.5*(ux*ux+uy*uy+uz*uz)*rho0/porosity)- m2);
	        jx = jx + Fx;
			m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0)- m4)
	                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
	        jy = jy + Fy;
			m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0)- m6)
	                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
	        jz = jz + Fz;
			m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0)- m8)
	                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
			m9 = m9 + rlx_setA*(((2*ux*ux-uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
			m10 = m10 + rlx_setA*( - m10);
	        //m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
			m11 = m11 + rlx_setA*(((uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
			m12 = m12 + rlx_setA*( - m12);
	        //m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
			m13 = m13 + rlx_setA*( (ux*uy*rho0/porosity) + 0.5*alpha*C*nx*ny - m13);
			m14 = m14 + rlx_setA*( (uy*uz*rho0/porosity) + 0.5*alpha*C*ny*nz - m14);
			m15 = m15 + rlx_setA*( (ux*uz*rho0/porosity) + 0.5*alpha*C*nx*nz - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);
	        //----------------------------------------------------------------------//

			//.................inverse transformation......................................................
			// q=0
			fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
			dist[n] = fq;

			// q = 1
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
			//nread = neighborList[n+Np];
			dist[nr2] = fq;

			// q=2
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
			//nread = neighborList[n];
			dist[nr1] = fq;

			// q = 3
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
			//nread = neighborList[n+3*Np];
			dist[nr4] = fq;

			// q = 4
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
			//nread = neighborList[n+2*Np];
			dist[nr3] = fq;

			// q = 5
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
			//nread = neighborList[n+5*Np];
			dist[nr6] = fq;

			// q = 6
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
			//nread = neighborList[n+4*Np];
			dist[nr5] = fq;

			// q = 7
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
			//nread = neighborList[n+7*Np];
			dist[nr8] = fq;

			// q = 8
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
					+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
			//nread = neighborList[n+6*Np];
			dist[nr7] = fq;

			// q = 9
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
			//nread = neighborList[n+9*Np];
			dist[nr10] = fq;

			// q = 10
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
			//nread = neighborList[n+8*Np];
			dist[nr9] = fq;

			// q = 11
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
			//nread = neighborList[n+11*Np];
			dist[nr12] = fq;

			// q = 12
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
					mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
			//nread = neighborList[n+10*Np];
			dist[nr11]= fq;

			// q = 13
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
			//nread = neighborList[n+13*Np];
			dist[nr14] = fq;

			// q= 14
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
			//nread = neighborList[n+12*Np];
			dist[nr13] = fq;


			// q = 15
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
			nread = neighborList[n+15*Np];
			dist[nread] = fq;

			// q = 16
			fq =  mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
			nread = neighborList[n+14*Np];
			dist[nread] = fq;


			// q = 17
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
			nread = neighborList[n+17*Np];
			dist[nread] = fq;

			// q = 18
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
			nread = neighborList[n+16*Np];
			dist[nread] = fq;
			//........................................................................

			// Instantiate mass transport distributions
			// Stationary value - distribution 0
			nAB = 1.0/(nA+nB);
			Aq[n] = 0.3333333333333333*nA;
			Bq[n] = 0.3333333333333333*nB;
			
			//...............................................
			// q = 0,2,4
			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
			jA = nA*ux;
			jB = nB*ux;		
			delta = beta*nA*nB*nAB*0.1111111111111111*nx;
			if (!(nA*nB*nAB>0)) delta=0;
	        //----------------newly added for better control of recoloring---------------//
	        if (nA/(nA+nB)>=Sn_grey && nA/(nA+nB) <= Sw_grey && porosity !=1.0){
	        	//delta = 0.0; 
	        	delta = -0.111111111111111*C*W*GreyDiff*nA*nB*nAB*nx;
	        	jA = 0.5*ux*(nA+nB)*(1.0+mobility_ratio);
	    		jB = 0.5*ux*(nA+nB)*(1.0-mobility_ratio);
	        }
	        if (nA/(nA+nB)>Sw_grey && porosity !=1.0) delta = -1.0*delta; 
	        //---------------------------------------------------------------------------//
	        if (RecoloringOff==true && porosity !=1.0) delta=0;
			a1 = (0.1111111111111111*(nA+4.5*jA))+delta;
			b1 = (0.1111111111111111*(nB+4.5*jB))-delta;
			a2 = (0.1111111111111111*(nA-4.5*jA))-delta;
			b2 = (0.1111111111111111*(nB-4.5*jB))+delta;

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
			jA = nA*uy;
			jB = nB*uy;		
			delta = beta*nA*nB*nAB*0.1111111111111111*ny;
			if (!(nA*nB*nAB>0)) delta=0;
	        //----------------newly added for better control of recoloring---------------//
	        if (nA/(nA+nB)>=Sn_grey && nA/(nA+nB) <= Sw_grey && porosity !=1.0){
	        	//delta = 0.0; 
	        	delta = -0.111111111111111*C*W*GreyDiff*nA*nB*nAB*ny;
	    		jA = 0.5*uy*(nA+nB)*(1.0+mobility_ratio);
	    		jB = 0.5*uy*(nA+nB)*(1.0-mobility_ratio);
	        }
	        if (nA/(nA+nB)>Sw_grey && porosity !=1.0) delta = -1.0*delta; 
	        //---------------------------------------------------------------------------//
	        if (RecoloringOff==true && porosity !=1.0) delta=0;
			a1 = (0.1111111111111111*(nA+4.5*jA))+delta;
			b1 = (0.1111111111111111*(nB+4.5*jB))-delta;
			a2 = (0.1111111111111111*(nA-4.5*jA))-delta;
			b2 = (0.1111111111111111*(nB-4.5*jB))+delta;

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
			jA = nA*uz;
			jB = nB*uz;		
			delta = beta*nA*nB*nAB*0.1111111111111111*nz;
			if (!(nA*nB*nAB>0)) delta=0;
	        //----------------newly added for better control of recoloring---------------//
	        if (nA/(nA+nB)>=Sn_grey && nA/(nA+nB) <= Sw_grey && porosity !=1.0){
	        	//delta = 0.0; 
	        	delta = -0.111111111111111*C*W*GreyDiff*nA*nB*nAB*nz;
	    		jA = 0.5*uz*(nA+nB)*(1.0+mobility_ratio);
	    		jB = 0.5*uz*(nA+nB)*(1.0-mobility_ratio);
	        }
	        if (nA/(nA+nB)>Sw_grey && porosity !=1.0) delta = -1.0*delta; 
	        //---------------------------------------------------------------------------//
	        if (RecoloringOff==true && porosity !=1.0) delta=0;

			a1 = (0.1111111111111111*(nA+4.5*jA))+delta;
			b1 = (0.1111111111111111*(nB+4.5*jB))-delta;
			a2 = (0.1111111111111111*(nA-4.5*jA))-delta;
			b2 = (0.1111111111111111*(nB-4.5*jB))+delta;

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

//CP: capillary penalty
// also turn off recoloring for grey nodes
void dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor_CP(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
        double *Phi, double *GreySolidW, double *GreySn, double *GreySw, double *GreyKn, double *GreyKw, double *Poros,
        double *Perm, double *Velocity, double *MobilityRatio, double *Pressure, 
        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
		double Gx, double Gy, double Gz, bool RecoloringOff, int strideY, int strideZ, int start, int finish, int Np,
		sycl::nd_item<3> item_ct1){
	int ijk,nn,n;
	double fq;
	// conserved momemnts
	double rho,jx,jy,jz;
    double ux,uy,uz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double m3,m5,m7;
	double nA,nB; // number density
	double a1,b1,a2,b2,nAB,delta;
	double C,nx,ny,nz; //color gradient magnitude and direction
	double phi,tau,rho0,rlx_setA,rlx_setB;

    double porosity;
    double perm;//voxel permeability
    double tau_eff;
    double mu_eff;//kinematic viscosity
    double Fx,Fy,Fz;
    double Fcpx,Fcpy,Fcpz;//capillary penalty force
    double W;//greyscale wetting strength
    double Sn_grey,Sw_grey;
    double GreyDiff=0.0e-4;

    /* Corey model parameters */
    double Kn_grey,Kw_grey;    
    double Swn,Krn_grey,Krw_grey,mobility_ratio,jA,jB;

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
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {
			// read the component number densities
			nA = Den[n];
			nB = Den[Np + n];

	        porosity = Poros[n];
	        //GreyDiff = Perm[n];
	        perm = 1.0;
	        W = GreySolidW[n];
	        Sn_grey = GreySn[n];
	        Sw_grey = GreySw[n];
	        Kn_grey = GreyKn[n];
	        Kw_grey = GreyKw[n];
	        
			// compute phase indicator field
			phi=(nA-nB)/(nA+nB);

			// local density
			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);

			// local relaxation time
			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
			rlx_setA = 1.f/tau;
			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
	        mu_eff = (tau_eff-0.5)/3.0;//kinematic viscosity

	        mobility_ratio = 1.0;
	        Krn_grey = 0.0;
	        Krw_grey = 0.0;
	        if (nA/(nA+nB)<Sn_grey && porosity !=1.0){
	        	perm = Kw_grey;
	        	Krw_grey = Kw_grey;
	        	Swn = 0.0;
	        }
	        else if (nA/(nA+nB)>=Sn_grey && nA/(nA+nB) <= Sw_grey && porosity !=1.0){ 
	        	Swn = (nA/(nA+nB) - Sn_grey) /(Sw_grey - Sn_grey);
	        	Krn_grey = Kn_grey*Swn*Swn; // Corey model with exponent = 2, make sure that W cannot shift to zero
	        	Krw_grey = Kw_grey*(1.0-Swn)*(1.0-Swn); // Corey model with exponent = 4, make sure that W cannot shift to zero
	        	// recompute the effective permeability
	        	perm = mu_eff*(Krn_grey*3.0/(tauA-0.5) + Krw_grey*3.0/(tauB-0.5));
	        	//mobility_ratio =(nA*Krn_grey*3.0/(tauA-0.5) - nB*Krw_grey*3.0/(tauB-0.5))/(nA*Krn_grey*3.0/(tauA-0.5) + nB*Krw_grey*3.0/(tauB-0.5));
	        }
	        else if (nA/(nA+nB)>Sw_grey && porosity !=1.0){
	        	perm = Kn_grey;
	        	Krn_grey = Kn_grey;
	        	Swn = 1.0;
	        }
	        /* compute the mobility ratio */
	        if (porosity != 1.0){
	        	mobility_ratio =(Krn_grey/(tauA-0.5) - Krw_grey/(tauB-0.5))/(Krn_grey/(tauA-0.5) + Krw_grey/(tauB-0.5));
	        }
	        else if (phi > 0.0){
	        	mobility_ratio = 1.0;
	        }
	        else {
	        	mobility_ratio = -1.0;
	        }
	        MobilityRatio[n] = mobility_ratio;
	        
			// Get the 1D index based on regular data layout
			ijk = Map[n];
			//					COMPUTE THE COLOR GRADIENT
			//........................................................................
			//.................Read Phase Indicator Values............................
			//........................................................................
			nn = ijk-1;							// neighbor index (get convention)
			m1 = Phi[nn];						// get neighbor for phi - 1
			//........................................................................
			nn = ijk+1;							// neighbor index (get convention)
			m2 = Phi[nn];						// get neighbor for phi - 2
			//........................................................................
			nn = ijk-strideY;							// neighbor index (get convention)
			m3 = Phi[nn];					// get neighbor for phi - 3
			//........................................................................
			nn = ijk+strideY;							// neighbor index (get convention)
			m4 = Phi[nn];					// get neighbor for phi - 4
			//........................................................................
			nn = ijk-strideZ;						// neighbor index (get convention)
			m5 = Phi[nn];					// get neighbor for phi - 5
			//........................................................................
			nn = ijk+strideZ;						// neighbor index (get convention)
			m6 = Phi[nn];					// get neighbor for phi - 6
			//........................................................................
			nn = ijk-strideY-1;						// neighbor index (get convention)
			m7 = Phi[nn];					// get neighbor for phi - 7
			//........................................................................
			nn = ijk+strideY+1;						// neighbor index (get convention)
			m8 = Phi[nn];					// get neighbor for phi - 8
			//........................................................................
			nn = ijk+strideY-1;						// neighbor index (get convention)
			m9 = Phi[nn];					// get neighbor for phi - 9
			//........................................................................
			nn = ijk-strideY+1;						// neighbor index (get convention)
			m10 = Phi[nn];					// get neighbor for phi - 10
			//........................................................................
			nn = ijk-strideZ-1;						// neighbor index (get convention)
			m11 = Phi[nn];					// get neighbor for phi - 11
			//........................................................................
			nn = ijk+strideZ+1;						// neighbor index (get convention)
			m12 = Phi[nn];					// get neighbor for phi - 12
			//........................................................................
			nn = ijk+strideZ-1;						// neighbor index (get convention)
			m13 = Phi[nn];					// get neighbor for phi - 13
			//........................................................................
			nn = ijk-strideZ+1;						// neighbor index (get convention)
			m14 = Phi[nn];					// get neighbor for phi - 14
			//........................................................................
			nn = ijk-strideZ-strideY;					// neighbor index (get convention)
			m15 = Phi[nn];					// get neighbor for phi - 15
			//........................................................................
			nn = ijk+strideZ+strideY;					// neighbor index (get convention)
			m16 = Phi[nn];					// get neighbor for phi - 16
			//........................................................................
			nn = ijk+strideZ-strideY;					// neighbor index (get convention)
			m17 = Phi[nn];					// get neighbor for phi - 17
			//........................................................................
			nn = ijk-strideZ+strideY;					// neighbor index (get convention)
			m18 = Phi[nn];					// get neighbor for phi - 18
			//............Compute the Color Gradient...................................
			nx = -3.0/18.0*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
			ny = -3.0/18.0*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
			nz = -3.0/18.0*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));

	        Fcpx = nx; 
	        Fcpy = ny; 
	        Fcpz = nz;
                double Fcp_mag =
                    sycl::sqrt(Fcpx * Fcpx + Fcpy * Fcpy + Fcpz * Fcpz);
                if (Fcp_mag==0.0) Fcpx=Fcpy=Fcpz=0.0;
	        //NOTE for open node (porosity=1.0),Fcp=0.0
                Fcpx *= alpha * W * (1.0 - porosity) / sycl::sqrt(perm);
                Fcpy *= alpha * W * (1.0 - porosity) / sycl::sqrt(perm);
                Fcpz *= alpha * W * (1.0 - porosity) / sycl::sqrt(perm);

                        //...........Normalize the Color Gradient.................................
                        C = sycl::sqrt(nx * nx + ny * ny + nz * nz);
                        double ColorMag = C;
			if (C==0.0) ColorMag=1.0;
			nx = nx/ColorMag;
			ny = ny/ColorMag;
			nz = nz/ColorMag;		

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

	        // Compute greyscale related parameters
	        ux = (jx/rho0+0.5*porosity*Gx+0.5*Fcpx/rho0)/(1.0+0.5*porosity*mu_eff/perm);
	        uy = (jy/rho0+0.5*porosity*Gy+0.5*Fcpy/rho0)/(1.0+0.5*porosity*mu_eff/perm);
	        uz = (jz/rho0+0.5*porosity*Gz+0.5*Fcpz/rho0)/(1.0+0.5*porosity*mu_eff/perm);
	        if (porosity==1.0){//i.e. open nodes
	            ux = (jx/rho0+0.5*porosity*Gx);
	            uy = (jy/rho0+0.5*porosity*Gy);
	            uz = (jz/rho0+0.5*porosity*Gz);
	        }

	        //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
	        Fx = rho0*(-porosity*mu_eff/perm*ux + porosity*Gx)+Fcpx;
	        Fy = rho0*(-porosity*mu_eff/perm*uy + porosity*Gy)+Fcpy;
	        Fz = rho0*(-porosity*mu_eff/perm*uz + porosity*Gz)+Fcpz;
	        if (porosity==1.0){
	            Fx=rho0*(porosity*Gx);
	            Fy=rho0*(porosity*Gy);
	            Fz=rho0*(porosity*Gz);
	        }

			// write the velocity 
			Velocity[n] = ux;
			Velocity[Np+n] = uy;
			Velocity[2*Np+n] = uz;
	        //Pressure[n] = rho/3.f/porosity;
	        Pressure[n] = rho/3.f;

			//........................................................................
			//..............carry out relaxation process..............................
			//..........Toelke, Fruediger et. al. 2006................................
	        //---------------- NO higher-order force -------------------------------//
			if (C == 0.0)	nx = ny = nz = 0.0;
			m1 = m1 + rlx_setA*((19*(ux*ux+uy*uy+uz*uz)*rho0/porosity - 11*rho) -19*alpha*C - m1);
			m2 = m2 + rlx_setA*((3*rho - 5.5*(ux*ux+uy*uy+uz*uz)*rho0/porosity)- m2);
	        jx = jx + Fx;
			m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0)- m4)
	                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
	        jy = jy + Fy;
			m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0)- m6)
	                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
	        jz = jz + Fz;
			m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0)- m8)
	                + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
			m9 = m9 + rlx_setA*(((2*ux*ux-uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
			m10 = m10 + rlx_setA*( - m10);
	        //m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
			m11 = m11 + rlx_setA*(((uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
			m12 = m12 + rlx_setA*( - m12);
	        //m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
			m13 = m13 + rlx_setA*( (ux*uy*rho0/porosity) + 0.5*alpha*C*nx*ny - m13);
			m14 = m14 + rlx_setA*( (uy*uz*rho0/porosity) + 0.5*alpha*C*ny*nz - m14);
			m15 = m15 + rlx_setA*( (ux*uz*rho0/porosity) + 0.5*alpha*C*nx*nz - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);
	        //----------------------------------------------------------------------//

			//.................inverse transformation......................................................
			// q=0
			fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
			dist[n] = fq;

			// q = 1
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
			dist[1*Np+n] = fq;

			// q=2
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
			dist[2*Np+n] = fq;

			// q = 3
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
			dist[3*Np+n] = fq;

			// q = 4
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
			dist[4*Np+n] = fq;

			// q = 5
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
			dist[5*Np+n] = fq;

			// q = 6
			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
			dist[6*Np+n] = fq;

			// q = 7
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
			dist[7*Np+n] = fq;


			// q = 8
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
					+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
			dist[8*Np+n] = fq;

			// q = 9
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
			dist[9*Np+n] = fq;

			// q = 10
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
			dist[10*Np+n] = fq;


			// q = 11
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
			dist[11*Np+n] = fq;

			// q = 12
			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
					mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
			dist[12*Np+n] = fq;

			// q = 13
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
			dist[13*Np+n] = fq;

			// q= 14
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
					-mrt_V12*m12-0.25*m15+0.125*(m16+m18);

			dist[14*Np+n] = fq;

			// q = 15
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
			dist[15*Np+n] = fq;

			// q = 16
			fq =  mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
			dist[16*Np+n] = fq;


			// q = 17
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
			dist[17*Np+n] = fq;

			// q = 18
			fq = mrt_V1*rho+mrt_V9*m1
					+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
					-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
			dist[18*Np+n] = fq;
			//........................................................................

			// Instantiate mass transport distributions
			// Stationary value - distribution 0
			nAB = 1.0/(nA+nB);
			Aq[n] = 0.3333333333333333*nA;
			Bq[n] = 0.3333333333333333*nB;
			
			//...............................................
			// q = 0,2,4
			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
			jA = nA*ux;
			jB = nB*ux;		
			delta = beta*nA*nB*nAB*0.1111111111111111*nx;
			if (!(nA*nB*nAB>0)) delta=0;
	        //----------------newly added for better control of recoloring---------------//
	        if (nA/(nA+nB)>=Sn_grey && nA/(nA+nB) <= Sw_grey && porosity !=1.0){
	        	//delta = 0.0; 
	        	delta = -0.111111111111111*C*W*GreyDiff*nA*nB*nAB*nx;
	    		jA = 0.5*ux*(nA+nB)*(1.0+mobility_ratio);
	    		jB = 0.5*ux*(nA+nB)*(1.0-mobility_ratio);
	        }
	        if (nA/(nA+nB)>Sw_grey && porosity !=1.0) delta = -1.0*delta; 
	        //---------------------------------------------------------------------------//
	        if (RecoloringOff==true && porosity !=1.0) delta=0;
			a1 = (0.1111111111111111*(nA+4.5*jA))+delta;
			b1 = (0.1111111111111111*(nB+4.5*jB))-delta;
			a2 = (0.1111111111111111*(nA-4.5*jA))-delta;
			b2 = (0.1111111111111111*(nB-4.5*jB))+delta;

			Aq[1*Np+n] = a1;
			Bq[1*Np+n] = b1;
			Aq[2*Np+n] = a2;
			Bq[2*Np+n] = b2;

			//...............................................
			// Cq = {0,1,0}
			jA = nA*uy;
			jB = nB*uy;		
			delta = beta*nA*nB*nAB*0.1111111111111111*ny;
			if (!(nA*nB*nAB>0)) delta=0;
	        //----------------newly added for better control of recoloring---------------//
	        if (nA/(nA+nB)>=Sn_grey && nA/(nA+nB) <= Sw_grey && porosity !=1.0){
	        	//delta = 0.0; 
	        	delta = -0.111111111111111*C*W*GreyDiff*nA*nB*nAB*ny;
	    		jA = 0.5*uy*(nA+nB)*(1.0+mobility_ratio);
	    		jB = 0.5*uy*(nA+nB)*(1.0-mobility_ratio);
	        }
	        if (nA/(nA+nB)>Sw_grey && porosity !=1.0) delta = -1.0*delta; 
	        //---------------------------------------------------------------------------//
	        if (RecoloringOff==true && porosity !=1.0) delta=0;
			a1 = (0.1111111111111111*(nA+4.5*jA))+delta;
			b1 = (0.1111111111111111*(nB+4.5*jB))-delta;
			a2 = (0.1111111111111111*(nA-4.5*jA))-delta;
			b2 = (0.1111111111111111*(nB-4.5*jB))+delta;

			Aq[3*Np+n] = a1;
			Bq[3*Np+n] = b1;
			Aq[4*Np+n] = a2;
			Bq[4*Np+n] = b2;

			//...............................................
			// q = 4
			// Cq = {0,0,1}
			jA = nA*uz;
			jB = nB*uz;		
			delta = beta*nA*nB*nAB*0.1111111111111111*nz;
			if (!(nA*nB*nAB>0)) delta=0;
	        //----------------newly added for better control of recoloring---------------//
	        if (nA/(nA+nB)>=Sn_grey && nA/(nA+nB) <= Sw_grey && porosity !=1.0){
	        	//delta = 0.0; 
	        	delta = -0.111111111111111*C*W*GreyDiff*nA*nB*nAB*nz;
	    		jA = 0.5*uz*(nA+nB)*(1.0+mobility_ratio);
	    		jB = 0.5*uz*(nA+nB)*(1.0-mobility_ratio);
	        }
	        if (nA/(nA+nB)>Sw_grey && porosity !=1.0) delta = -1.0*delta; 
	        //---------------------------------------------------------------------------//
	        if (RecoloringOff==true && porosity !=1.0) delta=0;

			a1 = (0.1111111111111111*(nA+4.5*jA))+delta;
			b1 = (0.1111111111111111*(nB+4.5*jB))-delta;
			a2 = (0.1111111111111111*(nA-4.5*jA))-delta;
			b2 = (0.1111111111111111*(nB-4.5*jB))+delta;

			Aq[5*Np+n] = a1;
			Bq[5*Np+n] = b1;
			Aq[6*Np+n] = a2;
			Bq[6*Np+n] = b2;
			//...............................................

		}
	}
}

void dvc_ScaLBL_PhaseField_InitFromRestart(double *Den, double *Aq, double *Bq, int start, int finish, int Np,
                                           sycl::nd_item<3> item_ct1){
	int idx;
	double nA,nB;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                idx = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                      s * item_ct1.get_local_range(2) +
                      item_ct1.get_local_id(2) + start;
                if (idx<finish) {

			nA = Den[idx];
			nB = Den[Np+idx];

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

//Model-1 & 4
extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColor(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
        double *Phi,double *GreySolidGrad, double *Poros,double *Perm,double *Vel, double *Pressure,
        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np){

	//cudaProfilerStart();
	//cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor, cudaFuncCachePreferL1);

        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor(
                        Map, dist, Aq, Bq, Den, Phi, GreySolidGrad, Poros, Perm,
                        Vel, Pressure, rhoA, rhoB, tauA, tauB, tauA_eff,
                        tauB_eff, alpha, beta, Fx, Fy, Fz, strideY, strideZ,
                        start, finish, Np, item_ct1);
            });
        /*
        DPCT1010:301: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;

        //cudaProfilerStop();

}

//Model-1 & 4
extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColor(int *d_neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *GreySolidGrad, double *Poros,double *Perm,double *Vel,double *Pressure, 
        double rhoA, double rhoB, double tauA, double tauB, double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np){

	//cudaProfilerStart();
	//cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor, cudaFuncCachePreferL1);

        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor(
                        d_neighborList, Map, dist, Aq, Bq, Den, Phi,
                        GreySolidGrad, Poros, Perm, Vel, Pressure, rhoA, rhoB,
                        tauA, tauB, tauA_eff, tauB_eff, alpha, beta, Fx, Fy, Fz,
                        strideY, strideZ, start, finish, Np, item_ct1);
            });

        /*
        DPCT1010:303: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;

        //cudaProfilerStop();
}

extern "C" void ScaLBL_PhaseField_InitFromRestart(double *Den, double *Aq, double *Bq, int start, int finish, int Np){
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_PhaseField_InitFromRestart(Den, Aq, Bq, start,
                                                          finish, Np, item_ct1);
            });
        /*
        DPCT1010:305: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}

//Model-1 & 4 with capillary pressure penalty
extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColor_CP(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
        double *Phi, double *GreySolidW, double *GreySn, double *GreySw, double *GreyKn, double *GreyKw,  double *Poros,
        double *Perm,double *Vel,double *MobilityRatio, double *Pressure,
        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, bool RecoloringOff, int strideY, int strideZ, int start, int finish, int Np){

        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor_CP(
                        Map, dist, Aq, Bq, Den, Phi, GreySolidW, GreySn, GreySw,
                        GreyKn, GreyKw, Poros, Perm, Vel, MobilityRatio,
                        Pressure, rhoA, rhoB, tauA, tauB, tauA_eff, tauB_eff,
                        alpha, beta, Fx, Fy, Fz, RecoloringOff, strideY,
                        strideZ, start, finish, Np, item_ct1);
            });

        /*
        DPCT1010:307: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}

//Model-1 & 4 with capillary pressure penalty
extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColor_CP(int *d_neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *GreySolidW, double *GreySn, double *GreySw, double *GreyKn, double *GreyKw, double *Poros, 
		double *Perm,double *Vel, double *MobilityRatio, double *Pressure, 
        double rhoA, double rhoB, double tauA, double tauB, double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, bool RecoloringOff, int strideY, int strideZ, int start, int finish, int Np){

        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor_CP(
                        d_neighborList, Map, dist, Aq, Bq, Den, Phi, GreySolidW,
                        GreySn, GreySw, GreyKn, GreyKw, Poros, Perm, Vel,
                        MobilityRatio, Pressure, rhoA, rhoB, tauA, tauB,
                        tauA_eff, tauB_eff, alpha, beta, Fx, Fy, Fz,
                        RecoloringOff, strideY, strideZ, start, finish, Np,
                        item_ct1);
            });

        /*
        DPCT1010:309: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}
