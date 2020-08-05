#include <stdio.h>
#include <math.h>

#define NBLOCKS 1024
#define NTHREADS 256

//Model-1 & 4
__global__ void dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor(int *neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den,
		 double *Phi, double *GreySolidGrad, double *Poros,double *Perm, double *Velocity, 
         double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff,double alpha, double beta,
		double Gx, double Gy, double Gz, int strideY, int strideZ, int start, int finish, int Np){

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
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
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
			C_phase = sqrt(nx_phase*nx_phase+ny_phase*ny_phase+nz_phase*nz_phase);

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
			C = sqrt(nx*nx+ny*ny+nz*nz);
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
            c1 = porosity*0.5*GeoFun/sqrt(perm);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jx/rho0+0.5*(porosity*Gx);
            vy = jy/rho0+0.5*(porosity*Gy);
            vz = jz/rho0+0.5*(porosity*Gz);
            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
            ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
            uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
            uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
            u_mag=sqrt(ux*ux+uy*uy+uz*uz);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx);
            Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy);
            Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz);
            if (porosity==1.0){
                Fx=rho0*(Gx);
                Fy=rho0*(Gy);
                Fz=rho0*(Gz);
            }

			// write the velocity 
			Velocity[n] = ux;
			Velocity[Np+n] = uy;
			Velocity[2*Np+n] = uz;

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
__global__  void dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
        double *Phi, double *GreySolidGrad, double *Poros,double *Perm, double *Velocity, 
        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
		double Gx, double Gy, double Gz, int strideY, int strideZ, int start, int finish, int Np){
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
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
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
			C_phase = sqrt(nx_phase*nx_phase+ny_phase*ny_phase+nz_phase*nz_phase);

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
			C = sqrt(nx*nx+ny*ny+nz*nz);
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
            c1 = porosity*0.5*GeoFun/sqrt(perm);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jx/rho0+0.5*(porosity*Gx);
            vy = jy/rho0+0.5*(porosity*Gy);
            vz = jz/rho0+0.5*(porosity*Gz);
            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
            ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
            uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
            uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
            u_mag=sqrt(ux*ux+uy*uy+uz*uz);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx);
            Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy);
            Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz);
            if (porosity==1.0){
                Fx=rho0*(Gx);
                Fy=rho0*(Gy);
                Fz=rho0*(Gz);
            }

			// write the velocity 
			Velocity[n] = ux;
			Velocity[Np+n] = uy;
			Velocity[2*Np+n] = uz;

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

////Model-2&3
//__global__ void dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor(int *neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den,
//		 double *Phi, double *GreySolidGrad, double *Poros,double *Perm, double *Velocity, 
//         double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff,double alpha, double beta,
//		 double Gx, double Gy, double Gz, int strideY, int strideZ, int start, int finish, int Np){
//
//	int n,nn,ijk,nread;
//	int nr1,nr2,nr3,nr4,nr5,nr6;
//	int nr7,nr8,nr9,nr10;
//	int nr11,nr12,nr13,nr14;
//	//int nr15,nr16,nr17,nr18;
//	double fq;
//	// conserved momemnts
//	double rho,jx,jy,jz;
//	double vx,vy,vz,v_mag;
//    double ux,uy,uz,u_mag;
//	// non-conserved moments
//	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
//	double m3,m5,m7;
//	double t1,t2,t4,t6,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18;
//	double t3,t5,t7;
//	double nA,nB; // number density
//	double a1,b1,a2,b2,nAB,delta;
//	double C,nx,ny,nz; //color gradient magnitude and direction
//	double phi,tau,rho0,rlx_setA,rlx_setB;
//
//    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
//    double porosity;
//    double perm;//voxel permeability
//    double c0, c1; //Guo's model parameters
//    double tau_eff;
//    double mu_eff;//kinematic viscosity
//    double nx_phase,ny_phase,nz_phase,C_phase;
//    double Fx,Fy,Fz;
//
//	const double mrt_V1=0.05263157894736842;
//	const double mrt_V2=0.012531328320802;
//	const double mrt_V3=0.04761904761904762;
//	const double mrt_V4=0.004594820384294068;
//	const double mrt_V5=0.01587301587301587;
//	const double mrt_V6=0.0555555555555555555555555;
//	const double mrt_V7=0.02777777777777778;
//	const double mrt_V8=0.08333333333333333;
//	const double mrt_V9=0.003341687552213868;
//	const double mrt_V10=0.003968253968253968;
//	const double mrt_V11=0.01388888888888889;
//	const double mrt_V12=0.04166666666666666;
//
//	int S = Np/NBLOCKS/NTHREADS + 1;
//	for (int s=0; s<S; s++){
//		//........Get 1-D index for this thread....................
//		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
//		if (n<finish) {
//			// read the component number densities
//			nA = Den[n];
//			nB = Den[Np + n];
//            porosity = Poros[n];
//            perm = Perm[n];
//
//			// compute phase indicator field
//			phi=(nA-nB)/(nA+nB);
//
//			// local density
//			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
//			// local relaxation time
//			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
//			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
//			rlx_setA = 1.f/tau;
//			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
//            mu_eff = (tau_eff-0.5)/3.0;//kinematic viscosity
//			
//			// Get the 1D index based on regular data layout
//			ijk = Map[n];
//			//					COMPUTE THE COLOR GRADIENT
//			//........................................................................
//			//.................Read Phase Indicator Values............................
//			//........................................................................
//			nn = ijk-1;							// neighbor index (get convention)
//			m1 = Phi[nn];						// get neighbor for phi - 1
//			t1 = m1+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t1)>1.0) t1 =((t1>0.0)-(t1<0.0))*(1.0-fabs(t1))+t1;
//			//........................................................................
//			nn = ijk+1;							// neighbor index (get convention)
//			m2 = Phi[nn];						// get neighbor for phi - 2
//			t2 = m2+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t2)>1.0) t2 =((t2>0.0)-(t2<0.0))*(1.0-fabs(t2))+t2;
//			//........................................................................
//			nn = ijk-strideY;							// neighbor index (get convention)
//			m3 = Phi[nn];					// get neighbor for phi - 3
//			t3 = m3+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t3)>1.0) t3 =((t3>0.0)-(t3<0.0))*(1.0-fabs(t3))+t3;
//			//........................................................................
//			nn = ijk+strideY;							// neighbor index (get convention)
//			m4 = Phi[nn];					// get neighbor for phi - 4
//			t4 = m4+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t4)>1.0) t4 =((t4>0.0)-(t4<0.0))*(1.0-fabs(t4))+t4;
//			//........................................................................
//			nn = ijk-strideZ;						// neighbor index (get convention)
//			m5 = Phi[nn];					// get neighbor for phi - 5
//			t5 = m5+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t5)>1.0) t5 =((t5>0.0)-(t5<0.0))*(1.0-fabs(t5))+t5;
//			//........................................................................
//			nn = ijk+strideZ;						// neighbor index (get convention)
//			m6 = Phi[nn];					// get neighbor for phi - 6
//			t6 = m6+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t6)>1.0) t6 =((t6>0.0)-(t6<0.0))*(1.0-fabs(t6))+t6;
//			//........................................................................
//			nn = ijk-strideY-1;						// neighbor index (get convention)
//			m7 = Phi[nn];					// get neighbor for phi - 7
//			t7 = m7+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t7)>1.0) t7 =((t7>0.0)-(t7<0.0))*(1.0-fabs(t7))+t7;
//			//........................................................................
//			nn = ijk+strideY+1;						// neighbor index (get convention)
//			m8 = Phi[nn];					// get neighbor for phi - 8
//			t8 = m8+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t8)>1.0) t8 =((t8>0.0)-(t8<0.0))*(1.0-fabs(t8))+t8;
//			//........................................................................
//			nn = ijk+strideY-1;						// neighbor index (get convention)
//			m9 = Phi[nn];					// get neighbor for phi - 9
//			t9 = m9+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t9)>1.0) t9 =((t9>0.0)-(t9<0.0))*(1.0-fabs(t9))+t9;
//			//........................................................................
//			nn = ijk-strideY+1;						// neighbor index (get convention)
//			m10 = Phi[nn];					// get neighbor for phi - 10
//			t10 = m10+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t10)>1.0) t10 =((t10>0.0)-(t10<0.0))*(1.0-fabs(t10))+t10;
//			//........................................................................
//			nn = ijk-strideZ-1;						// neighbor index (get convention)
//			m11 = Phi[nn];					// get neighbor for phi - 11
//			t11 = m11+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t11)>1.0) t11 =((t11>0.0)-(t11<0.0))*(1.0-fabs(t11))+t11;
//			//........................................................................
//			nn = ijk+strideZ+1;						// neighbor index (get convention)
//			m12 = Phi[nn];					// get neighbor for phi - 12
//			t12 = m12+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t12)>1.0) t12 =((t12>0.0)-(t12<0.0))*(1.0-fabs(t12))+t12;
//			//........................................................................
//			nn = ijk+strideZ-1;						// neighbor index (get convention)
//			m13 = Phi[nn];					// get neighbor for phi - 13
//			t13 = m13+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t13)>1.0) t13 =((t13>0.0)-(t13<0.0))*(1.0-fabs(t13))+t13;
//			//........................................................................
//			nn = ijk-strideZ+1;						// neighbor index (get convention)
//			m14 = Phi[nn];					// get neighbor for phi - 14
//			t14 = m14+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t14)>1.0) t14 =((t14>0.0)-(t14<0.0))*(1.0-fabs(t14))+t14;
//			//........................................................................
//			nn = ijk-strideZ-strideY;					// neighbor index (get convention)
//			m15 = Phi[nn];					// get neighbor for phi - 15
//			t15 = m15+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t15)>1.0) t15 =((t15>0.0)-(t15<0.0))*(1.0-fabs(t15))+t15;
//			//........................................................................
//			nn = ijk+strideZ+strideY;					// neighbor index (get convention)
//			m16 = Phi[nn];					// get neighbor for phi - 16
//			t16 = m16+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t16)>1.0) t16 =((t16>0.0)-(t16<0.0))*(1.0-fabs(t16))+t16;
//			//........................................................................
//			nn = ijk+strideZ-strideY;					// neighbor index (get convention)
//			m17 = Phi[nn];					// get neighbor for phi - 17
//			t17 = m17+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t17)>1.0) t17 =((t17>0.0)-(t17<0.0))*(1.0-fabs(t17))+t17;
//			//........................................................................
//			nn = ijk-strideZ+strideY;					// neighbor index (get convention)
//			m18 = Phi[nn];					// get neighbor for phi - 18
//			t18 = m18+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t18)>1.0) t18 =((t18>0.0)-(t18<0.0))*(1.0-fabs(t18))+t18;
//			//............Compute the Color Gradient...................................
//			nx_phase = -(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
//			ny_phase = -(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
//			nz_phase = -(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
//			C_phase = sqrt(nx_phase*nx_phase+ny_phase*ny_phase+nz_phase*nz_phase);
//            //correct the normal color gradient by considering the effect of grey solid
//			nx = -(t1-t2+0.5*(t7-t8+t9-t10+t11-t12+t13-t14));
//			ny = -(t3-t4+0.5*(t7-t8-t9+t10+t15-t16+t17-t18));
//			nz = -(t5-t6+0.5*(t11-t12-t13+t14+t15-t16-t17+t18));
//
//            if (C_phase==0.0){//i.e. if in a bulk phase, there is no need for grey-solid correction
//                nx = nx_phase; 
//                ny = ny_phase;
//                nz = nz_phase;
//            }
//
//			//...........Normalize the Color Gradient.................................
//			C = sqrt(nx*nx+ny*ny+nz*nz);
//			double ColorMag = C;
//			if (C==0.0) ColorMag=1.0;
//			nx = nx/ColorMag;
//			ny = ny/ColorMag;
//			nz = nz/ColorMag;		
//
//			// q=0
//			fq = dist[n];
//			rho = fq;
//			m1  = -30.0*fq;
//			m2  = 12.0*fq;
//
//			// q=1
//			//nread = neighborList[n]; // neighbor 2 
//			//fq = dist[nread]; // reading the f1 data into register fq		
//			nr1 = neighborList[n]; 
//			fq = dist[nr1]; // reading the f1 data into register fq
//			rho += fq;
//			m1 -= 11.0*fq;
//			m2 -= 4.0*fq;
//			jx = fq;
//			m4 = -4.0*fq;
//			m9 = 2.0*fq;
//			m10 = -4.0*fq;
//
//			// f2 = dist[10*Np+n];
//			//nread = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
//			//fq = dist[nread];  // reading the f2 data into register fq
//			nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
//			fq = dist[nr2];  // reading the f2 data into register fq
//			rho += fq;
//			m1 -= 11.0*(fq);
//			m2 -= 4.0*(fq);
//			jx -= fq;
//			m4 += 4.0*(fq);
//			m9 += 2.0*(fq);
//			m10 -= 4.0*(fq);
//
//			// q=3
//			//nread = neighborList[n+2*Np]; // neighbor 4
//			//fq = dist[nread];
//			nr3 = neighborList[n+2*Np]; // neighbor 4
//			fq = dist[nr3];
//			rho += fq;
//			m1 -= 11.0*fq;
//			m2 -= 4.0*fq;
//			jy = fq;
//			m6 = -4.0*fq;
//			m9 -= fq;
//			m10 += 2.0*fq;
//			m11 = fq;
//			m12 = -2.0*fq;
//
//			// q = 4
//			//nread = neighborList[n+3*Np]; // neighbor 3
//			//fq = dist[nread];
//			nr4 = neighborList[n+3*Np]; // neighbor 3
//			fq = dist[nr4];
//			rho+= fq;
//			m1 -= 11.0*fq;
//			m2 -= 4.0*fq;
//			jy -= fq;
//			m6 += 4.0*fq;
//			m9 -= fq;
//			m10 += 2.0*fq;
//			m11 += fq;
//			m12 -= 2.0*fq;
//
//			// q=5
//			//nread = neighborList[n+4*Np];
//			//fq = dist[nread];
//			nr5 = neighborList[n+4*Np];
//			fq = dist[nr5];
//			rho += fq;
//			m1 -= 11.0*fq;
//			m2 -= 4.0*fq;
//			jz = fq;
//			m8 = -4.0*fq;
//			m9 -= fq;
//			m10 += 2.0*fq;
//			m11 -= fq;
//			m12 += 2.0*fq;
//
//
//			// q = 6
//			//nread = neighborList[n+5*Np];
//			//fq = dist[nread];
//			nr6 = neighborList[n+5*Np];
//			fq = dist[nr6];
//			rho+= fq;
//			m1 -= 11.0*fq;
//			m2 -= 4.0*fq;
//			jz -= fq;
//			m8 += 4.0*fq;
//			m9 -= fq;
//			m10 += 2.0*fq;
//			m11 -= fq;
//			m12 += 2.0*fq;
//
//			// q=7
//			//nread = neighborList[n+6*Np];
//			//fq = dist[nread];
//			nr7 = neighborList[n+6*Np];
//			fq = dist[nr7];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx += fq;
//			m4 += fq;
//			jy += fq;
//			m6 += fq;
//			m9  += fq;
//			m10 += fq;
//			m11 += fq;
//			m12 += fq;
//			m13 = fq;
//			m16 = fq;
//			m17 = -fq;
//
//			// q = 8
//			//nread = neighborList[n+7*Np];
//			//fq = dist[nread];
//			nr8 = neighborList[n+7*Np];
//			fq = dist[nr8];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx -= fq;
//			m4 -= fq;
//			jy -= fq;
//			m6 -= fq;
//			m9 += fq;
//			m10 += fq;
//			m11 += fq;
//			m12 += fq;
//			m13 += fq;
//			m16 -= fq;
//			m17 += fq;
//
//			// q=9
//			//nread = neighborList[n+8*Np];
//			//fq = dist[nread];
//			nr9 = neighborList[n+8*Np];
//			fq = dist[nr9];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx += fq;
//			m4 += fq;
//			jy -= fq;
//			m6 -= fq;
//			m9 += fq;
//			m10 += fq;
//			m11 += fq;
//			m12 += fq;
//			m13 -= fq;
//			m16 += fq;
//			m17 += fq;
//
//			// q = 10
//			//nread = neighborList[n+9*Np];
//			//fq = dist[nread];
//			nr10 = neighborList[n+9*Np];
//			fq = dist[nr10];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx -= fq;
//			m4 -= fq;
//			jy += fq;
//			m6 += fq;
//			m9 += fq;
//			m10 += fq;
//			m11 += fq;
//			m12 += fq;
//			m13 -= fq;
//			m16 -= fq;
//			m17 -= fq;
//
//			// q=11
//			//nread = neighborList[n+10*Np];
//			//fq = dist[nread];
//			nr11 = neighborList[n+10*Np];
//			fq = dist[nr11];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx += fq;
//			m4 += fq;
//			jz += fq;
//			m8 += fq;
//			m9 += fq;
//			m10 += fq;
//			m11 -= fq;
//			m12 -= fq;
//			m15 = fq;
//			m16 -= fq;
//			m18 = fq;
//
//			// q=12
//			//nread = neighborList[n+11*Np];
//			//fq = dist[nread];
//			nr12 = neighborList[n+11*Np];
//			fq = dist[nr12];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx -= fq;
//			m4 -= fq;
//			jz -= fq;
//			m8 -= fq;
//			m9 += fq;
//			m10 += fq;
//			m11 -= fq;
//			m12 -= fq;
//			m15 += fq;
//			m16 += fq;
//			m18 -= fq;
//
//			// q=13
//			//nread = neighborList[n+12*Np];
//			//fq = dist[nread];
//			nr13 = neighborList[n+12*Np];
//			fq = dist[nr13];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx += fq;
//			m4 += fq;
//			jz -= fq;
//			m8 -= fq;
//			m9 += fq;
//			m10 += fq;
//			m11 -= fq;
//			m12 -= fq;
//			m15 -= fq;
//			m16 -= fq;
//			m18 -= fq;
//
//			// q=14
//			//nread = neighborList[n+13*Np];
//			//fq = dist[nread];
//			nr14 = neighborList[n+13*Np];
//			fq = dist[nr14];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx -= fq;
//			m4 -= fq;
//			jz += fq;
//			m8 += fq;
//			m9 += fq;
//			m10 += fq;
//			m11 -= fq;
//			m12 -= fq;
//			m15 -= fq;
//			m16 += fq;
//			m18 += fq;
//
//			// q=15
//			nread = neighborList[n+14*Np];
//			fq = dist[nread];
//			//fq = dist[17*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jy += fq;
//			m6 += fq;
//			jz += fq;
//			m8 += fq;
//			m9 -= 2.0*fq;
//			m10 -= 2.0*fq;
//			m14 = fq;
//			m17 += fq;
//			m18 -= fq;
//
//			// q=16
//			nread = neighborList[n+15*Np];
//			fq = dist[nread];
//			//fq = dist[8*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jy -= fq;
//			m6 -= fq;
//			jz -= fq;
//			m8 -= fq;
//			m9 -= 2.0*fq;
//			m10 -= 2.0*fq;
//			m14 += fq;
//			m17 -= fq;
//			m18 += fq;
//
//			// q=17
//			//fq = dist[18*Np+n];
//			nread = neighborList[n+16*Np];
//			fq = dist[nread];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jy += fq;
//			m6 += fq;
//			jz -= fq;
//			m8 -= fq;
//			m9 -= 2.0*fq;
//			m10 -= 2.0*fq;
//			m14 -= fq;
//			m17 += fq;
//			m18 += fq;
//
//			// q=18
//			nread = neighborList[n+17*Np];
//			fq = dist[nread];
//			//fq = dist[9*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jy -= fq;
//			m6 -= fq;
//			jz += fq;
//			m8 += fq;
//			m9 -= 2.0*fq;
//			m10 -= 2.0*fq;
//			m14 -= fq;
//			m17 -= fq;
//			m18 -= fq;
//			
//            // Compute greyscale related parameters
//            c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
//            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
//            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
//            c1 = porosity*0.5*GeoFun/sqrt(perm);
//            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes
//
//            vx = jx/rho0+0.5*(porosity*Gx);
//            vy = jy/rho0+0.5*(porosity*Gy);
//            vz = jz/rho0+0.5*(porosity*Gz);
//            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
//            ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
//            uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
//            uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
//            u_mag=sqrt(ux*ux+uy*uy+uz*uz);
//
//            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
//            Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx);
//            Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy);
//            Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz);
//            if (porosity==1.0){
//                Fx=rho0*(Gx);
//                Fy=rho0*(Gy);
//                Fz=rho0*(Gz);
//            }
//
//			// write the velocity 
//			Velocity[n] = ux;
//			Velocity[Np+n] = uy;
//			Velocity[2*Np+n] = uz;
//
//			//........................................................................
//			//..............carry out relaxation process..............................
//			//..........Toelke, Fruediger et. al. 2006................................
//			if (C == 0.0)	nx = ny = nz = 0.0;
//			m1 = m1 + rlx_setA*((19*(ux*ux+uy*uy+uz*uz)*rho0/porosity - 11*rho) -19*alpha*C - m1);
//			m2 = m2 + rlx_setA*((3*rho - 5.5*(ux*ux+uy*uy+uz*uz)*rho0/porosity)- m2);
//            jx = jx + Fx;
//			m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0)- m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//			m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0)- m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//			m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0)- m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//			m9 = m9 + rlx_setA*(((2*ux*ux-uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
//			m10 = m10 + rlx_setA*( - m10);
//            //m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
//			m11 = m11 + rlx_setA*(((uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
//			m12 = m12 + rlx_setA*( - m12);
//            //m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
//			m13 = m13 + rlx_setA*( (ux*uy*rho0/porosity) + 0.5*alpha*C*nx*ny - m13);
//			m14 = m14 + rlx_setA*( (uy*uz*rho0/porosity) + 0.5*alpha*C*ny*nz - m14);
//			m15 = m15 + rlx_setA*( (ux*uz*rho0/porosity) + 0.5*alpha*C*nx*nz - m15);
//			m16 = m16 + rlx_setB*( - m16);
//			m17 = m17 + rlx_setB*( - m17);
//			m18 = m18 + rlx_setB*( - m18);
//
//			//.................inverse transformation......................................................
//			// q=0
//			fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
//			dist[n] = fq;
//
//			// q = 1
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
//			//nread = neighborList[n+Np];
//			dist[nr2] = fq;
//
//			// q=2
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
//			//nread = neighborList[n];
//			dist[nr1] = fq;
//
//			// q = 3
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
//			//nread = neighborList[n+3*Np];
//			dist[nr4] = fq;
//
//			// q = 4
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
//			//nread = neighborList[n+2*Np];
//			dist[nr3] = fq;
//
//			// q = 5
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
//			//nread = neighborList[n+5*Np];
//			dist[nr6] = fq;
//
//			// q = 6
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
//			//nread = neighborList[n+4*Np];
//			dist[nr5] = fq;
//
//			// q = 7
//			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
//					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
//			//nread = neighborList[n+7*Np];
//			dist[nr8] = fq;
//
//			// q = 8
//			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
//					+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
//			//nread = neighborList[n+6*Np];
//			dist[nr7] = fq;
//
//			// q = 9
//			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
//					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
//			//nread = neighborList[n+9*Np];
//			dist[nr10] = fq;
//
//			// q = 10
//			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
//					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
//			//nread = neighborList[n+8*Np];
//			dist[nr9] = fq;
//
//			// q = 11
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
//					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
//					-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
//			//nread = neighborList[n+11*Np];
//			dist[nr12] = fq;
//
//			// q = 12
//			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
//					mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
//			//nread = neighborList[n+10*Np];
//			dist[nr11]= fq;
//
//			// q = 13
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
//					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
//					-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
//			//nread = neighborList[n+13*Np];
//			dist[nr14] = fq;
//
//			// q= 14
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
//					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
//					-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
//			//nread = neighborList[n+12*Np];
//			dist[nr13] = fq;
//
//
//			// q = 15
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
//					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
//			nread = neighborList[n+15*Np];
//			dist[nread] = fq;
//
//			// q = 16
//			fq =  mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
//					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
//			nread = neighborList[n+14*Np];
//			dist[nread] = fq;
//
//
//			// q = 17
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
//					-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
//			nread = neighborList[n+17*Np];
//			dist[nread] = fq;
//
//			// q = 18
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
//					-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
//			nread = neighborList[n+16*Np];
//			dist[nread] = fq;
//			//........................................................................
//
//			// Instantiate mass transport distributions
//			// Stationary value - distribution 0
//			nAB = 1.0/(nA+nB);
//			Aq[n] = 0.3333333333333333*nA;
//			Bq[n] = 0.3333333333333333*nB;
//
//			//...............................................
//			// q = 0,2,4
//			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
//			delta = beta*nA*nB*nAB*0.1111111111111111*nx;
//			if (!(nA*nB*nAB>0)) delta=0;
//			a1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
//			b1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
//			a2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
//			b2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;
//
//			// q = 1
//			//nread = neighborList[n+Np];
//			Aq[nr2] = a1;
//			Bq[nr2] = b1;
//			// q=2
//			//nread = neighborList[n];
//			Aq[nr1] = a2;
//			Bq[nr1] = b2;
//
//			//...............................................
//			// Cq = {0,1,0}
//			delta = beta*nA*nB*nAB*0.1111111111111111*ny;
//			if (!(nA*nB*nAB>0)) delta=0;
//			a1 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
//			b1 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
//			a2 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
//			b2 = nB*(0.1111111111111111*(1-4.5*uy))+delta;
//
//			// q = 3
//			//nread = neighborList[n+3*Np];
//			Aq[nr4] = a1;
//			Bq[nr4] = b1;
//			// q = 4
//			//nread = neighborList[n+2*Np];
//			Aq[nr3] = a2;
//			Bq[nr3] = b2;
//
//			//...............................................
//			// q = 4
//			// Cq = {0,0,1}
//			delta = beta*nA*nB*nAB*0.1111111111111111*nz;
//			if (!(nA*nB*nAB>0)) delta=0;
//			a1 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
//			b1 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
//			a2 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
//			b2 = nB*(0.1111111111111111*(1-4.5*uz))+delta;
//
//			// q = 5
//			//nread = neighborList[n+5*Np];
//			Aq[nr6] = a1;
//			Bq[nr6] = b1;
//			// q = 6
//			//nread = neighborList[n+4*Np];
//			Aq[nr5] = a2;
//			Bq[nr5] = b2;
//			//...............................................
//		}
//	}
//}
//
////Model-2&3
//__global__  void dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
//        double *Phi, double *GreySolidGrad, double *Poros,double *Perm, double *Velocity, 
//        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
//		double Gx, double Gy, double Gz, int strideY, int strideZ, int start, int finish, int Np){
//	int ijk,nn,n;
//	double fq;
//	// conserved momemnts
//	double rho,jx,jy,jz;
//	double vx,vy,vz,v_mag;
//    double ux,uy,uz,u_mag;
//	// non-conserved moments
//	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
//	double m3,m5,m7;
//	double t1,t2,t4,t6,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18;
//	double t3,t5,t7;
//	double nA,nB; // number density
//	double a1,b1,a2,b2,nAB,delta;
//	double C,nx,ny,nz; //color gradient magnitude and direction
//	double phi,tau,rho0,rlx_setA,rlx_setB;
//
//    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
//    double porosity;
//    double perm;//voxel permeability
//    double c0, c1; //Guo's model parameters
//    double tau_eff;
//    double mu_eff;//kinematic viscosity
//    double nx_phase,ny_phase,nz_phase,C_phase;
//    double Fx,Fy,Fz;
//
//	const double mrt_V1=0.05263157894736842;
//	const double mrt_V2=0.012531328320802;
//	const double mrt_V3=0.04761904761904762;
//	const double mrt_V4=0.004594820384294068;
//	const double mrt_V5=0.01587301587301587;
//	const double mrt_V6=0.0555555555555555555555555;
//	const double mrt_V7=0.02777777777777778;
//	const double mrt_V8=0.08333333333333333;
//	const double mrt_V9=0.003341687552213868;
//	const double mrt_V10=0.003968253968253968;
//	const double mrt_V11=0.01388888888888889;
//	const double mrt_V12=0.04166666666666666;
//
//	int S = Np/NBLOCKS/NTHREADS + 1;
//	for (int s=0; s<S; s++){
//		//........Get 1-D index for this thread....................
//		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
//		if (n<finish) {
//
//			// read the component number densities
//			nA = Den[n];
//			nB = Den[Np + n];
//            porosity = Poros[n];
//            perm = Perm[n];
//
//			// compute phase indicator field
//			phi=(nA-nB)/(nA+nB);
//
//			// local density
//			rho0=rhoA + 0.5*(1.0-phi)*(rhoB-rhoA);
//			// local relaxation time
//			tau=tauA + 0.5*(1.0-phi)*(tauB-tauA);
//			tau_eff=tauA_eff + 0.5*(1.0-phi)*(tauB_eff-tauA_eff);
//			rlx_setA = 1.f/tau;
//			rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
//            mu_eff = (tau_eff-0.5)/3.0;//kinematic viscosity
//
//			// Get the 1D index based on regular data layout
//			ijk = Map[n];
//			//					COMPUTE THE COLOR GRADIENT
//			//........................................................................
//			//.................Read Phase Indicator Values............................
//			//........................................................................
//			nn = ijk-1;							// neighbor index (get convention)
//			m1 = Phi[nn];						// get neighbor for phi - 1
//			t1 = m1+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t1)>1.0) t1 =((t1>0.0)-(t1<0.0))*(1.0-fabs(t1))+t1;
//			//........................................................................
//			nn = ijk+1;							// neighbor index (get convention)
//			m2 = Phi[nn];						// get neighbor for phi - 2
//			t2 = m2+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t2)>1.0) t2 =((t2>0.0)-(t2<0.0))*(1.0-fabs(t2))+t2;
//			//........................................................................
//			nn = ijk-strideY;							// neighbor index (get convention)
//			m3 = Phi[nn];					// get neighbor for phi - 3
//			t3 = m3+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t3)>1.0) t3 =((t3>0.0)-(t3<0.0))*(1.0-fabs(t3))+t3;
//			//........................................................................
//			nn = ijk+strideY;							// neighbor index (get convention)
//			m4 = Phi[nn];					// get neighbor for phi - 4
//			t4 = m4+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t4)>1.0) t4 =((t4>0.0)-(t4<0.0))*(1.0-fabs(t4))+t4;
//			//........................................................................
//			nn = ijk-strideZ;						// neighbor index (get convention)
//			m5 = Phi[nn];					// get neighbor for phi - 5
//			t5 = m5+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t5)>1.0) t5 =((t5>0.0)-(t5<0.0))*(1.0-fabs(t5))+t5;
//			//........................................................................
//			nn = ijk+strideZ;						// neighbor index (get convention)
//			m6 = Phi[nn];					// get neighbor for phi - 6
//			t6 = m6+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t6)>1.0) t6 =((t6>0.0)-(t6<0.0))*(1.0-fabs(t6))+t6;
//			//........................................................................
//			nn = ijk-strideY-1;						// neighbor index (get convention)
//			m7 = Phi[nn];					// get neighbor for phi - 7
//			t7 = m7+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t7)>1.0) t7 =((t7>0.0)-(t7<0.0))*(1.0-fabs(t7))+t7;
//			//........................................................................
//			nn = ijk+strideY+1;						// neighbor index (get convention)
//			m8 = Phi[nn];					// get neighbor for phi - 8
//			t8 = m8+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t8)>1.0) t8 =((t8>0.0)-(t8<0.0))*(1.0-fabs(t8))+t8;
//			//........................................................................
//			nn = ijk+strideY-1;						// neighbor index (get convention)
//			m9 = Phi[nn];					// get neighbor for phi - 9
//			t9 = m9+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t9)>1.0) t9 =((t9>0.0)-(t9<0.0))*(1.0-fabs(t9))+t9;
//			//........................................................................
//			nn = ijk-strideY+1;						// neighbor index (get convention)
//			m10 = Phi[nn];					// get neighbor for phi - 10
//			t10 = m10+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t10)>1.0) t10 =((t10>0.0)-(t10<0.0))*(1.0-fabs(t10))+t10;
//			//........................................................................
//			nn = ijk-strideZ-1;						// neighbor index (get convention)
//			m11 = Phi[nn];					// get neighbor for phi - 11
//			t11 = m11+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t11)>1.0) t11 =((t11>0.0)-(t11<0.0))*(1.0-fabs(t11))+t11;
//			//........................................................................
//			nn = ijk+strideZ+1;						// neighbor index (get convention)
//			m12 = Phi[nn];					// get neighbor for phi - 12
//			t12 = m12+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t12)>1.0) t12 =((t12>0.0)-(t12<0.0))*(1.0-fabs(t12))+t12;
//			//........................................................................
//			nn = ijk+strideZ-1;						// neighbor index (get convention)
//			m13 = Phi[nn];					// get neighbor for phi - 13
//			t13 = m13+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t13)>1.0) t13 =((t13>0.0)-(t13<0.0))*(1.0-fabs(t13))+t13;
//			//........................................................................
//			nn = ijk-strideZ+1;						// neighbor index (get convention)
//			m14 = Phi[nn];					// get neighbor for phi - 14
//			t14 = m14+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t14)>1.0) t14 =((t14>0.0)-(t14<0.0))*(1.0-fabs(t14))+t14;
//			//........................................................................
//			nn = ijk-strideZ-strideY;					// neighbor index (get convention)
//			m15 = Phi[nn];					// get neighbor for phi - 15
//			t15 = m15+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t15)>1.0) t15 =((t15>0.0)-(t15<0.0))*(1.0-fabs(t15))+t15;
//			//........................................................................
//			nn = ijk+strideZ+strideY;					// neighbor index (get convention)
//			m16 = Phi[nn];					// get neighbor for phi - 16
//			t16 = m16+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t16)>1.0) t16 =((t16>0.0)-(t16<0.0))*(1.0-fabs(t16))+t16;
//			//........................................................................
//			nn = ijk+strideZ-strideY;					// neighbor index (get convention)
//			m17 = Phi[nn];					// get neighbor for phi - 17
//			t17 = m17+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t17)>1.0) t17 =((t17>0.0)-(t17<0.0))*(1.0-fabs(t17))+t17;
//			//........................................................................
//			nn = ijk-strideZ+strideY;					// neighbor index (get convention)
//			m18 = Phi[nn];					// get neighbor for phi - 18
//			t18 = m18+(1.0-porosity)*GreySolidGrad[nn];				
//            if (fabs(t18)>1.0) t18 =((t18>0.0)-(t18<0.0))*(1.0-fabs(t18))+t18;
//			//............Compute the Color Gradient...................................
//			nx_phase = -(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
//			ny_phase = -(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
//			nz_phase = -(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
//			C_phase = sqrt(nx_phase*nx_phase+ny_phase*ny_phase+nz_phase*nz_phase);
//            //correct the normal color gradient by considering the effect of grey solid
//			nx = -(t1-t2+0.5*(t7-t8+t9-t10+t11-t12+t13-t14));
//			ny = -(t3-t4+0.5*(t7-t8-t9+t10+t15-t16+t17-t18));
//			nz = -(t5-t6+0.5*(t11-t12-t13+t14+t15-t16-t17+t18));
//
//            if (C_phase==0.0){
//                nx = nx_phase; 
//                ny = ny_phase;
//                nz = nz_phase;
//            }
//
//			//...........Normalize the Color Gradient.................................
//			C = sqrt(nx*nx+ny*ny+nz*nz);
//			double ColorMag = C;
//			if (C==0.0) ColorMag=1.0;
//			nx = nx/ColorMag;
//			ny = ny/ColorMag;
//			nz = nz/ColorMag;		
//
//			// q=0
//			fq = dist[n];
//			rho = fq;
//			m1  = -30.0*fq;
//			m2  = 12.0*fq;
//
//			// q=1
//			fq = dist[2*Np+n];
//			rho += fq;
//			m1 -= 11.0*fq;
//			m2 -= 4.0*fq;
//			jx = fq;
//			m4 = -4.0*fq;
//			m9 = 2.0*fq;
//			m10 = -4.0*fq;
//
//			// f2 = dist[10*Np+n];
//			fq = dist[1*Np+n];
//			rho += fq;
//			m1 -= 11.0*(fq);
//			m2 -= 4.0*(fq);
//			jx -= fq;
//			m4 += 4.0*(fq);
//			m9 += 2.0*(fq);
//			m10 -= 4.0*(fq);
//
//			// q=3
//			fq = dist[4*Np+n];
//			rho += fq;
//			m1 -= 11.0*fq;
//			m2 -= 4.0*fq;
//			jy = fq;
//			m6 = -4.0*fq;
//			m9 -= fq;
//			m10 += 2.0*fq;
//			m11 = fq;
//			m12 = -2.0*fq;
//
//			// q = 4
//			fq = dist[3*Np+n];
//			rho+= fq;
//			m1 -= 11.0*fq;
//			m2 -= 4.0*fq;
//			jy -= fq;
//			m6 += 4.0*fq;
//			m9 -= fq;
//			m10 += 2.0*fq;
//			m11 += fq;
//			m12 -= 2.0*fq;
//
//			// q=5
//			fq = dist[6*Np+n];
//			rho += fq;
//			m1 -= 11.0*fq;
//			m2 -= 4.0*fq;
//			jz = fq;
//			m8 = -4.0*fq;
//			m9 -= fq;
//			m10 += 2.0*fq;
//			m11 -= fq;
//			m12 += 2.0*fq;
//
//			// q = 6
//			fq = dist[5*Np+n];
//			rho+= fq;
//			m1 -= 11.0*fq;
//			m2 -= 4.0*fq;
//			jz -= fq;
//			m8 += 4.0*fq;
//			m9 -= fq;
//			m10 += 2.0*fq;
//			m11 -= fq;
//			m12 += 2.0*fq;
//
//			// q=7
//			fq = dist[8*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx += fq;
//			m4 += fq;
//			jy += fq;
//			m6 += fq;
//			m9  += fq;
//			m10 += fq;
//			m11 += fq;
//			m12 += fq;
//			m13 = fq;
//			m16 = fq;
//			m17 = -fq;
//
//			// q = 8
//			fq = dist[7*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx -= fq;
//			m4 -= fq;
//			jy -= fq;
//			m6 -= fq;
//			m9 += fq;
//			m10 += fq;
//			m11 += fq;
//			m12 += fq;
//			m13 += fq;
//			m16 -= fq;
//			m17 += fq;
//
//			// q=9
//			fq = dist[10*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx += fq;
//			m4 += fq;
//			jy -= fq;
//			m6 -= fq;
//			m9 += fq;
//			m10 += fq;
//			m11 += fq;
//			m12 += fq;
//			m13 -= fq;
//			m16 += fq;
//			m17 += fq;
//
//			// q = 10
//			fq = dist[9*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx -= fq;
//			m4 -= fq;
//			jy += fq;
//			m6 += fq;
//			m9 += fq;
//			m10 += fq;
//			m11 += fq;
//			m12 += fq;
//			m13 -= fq;
//			m16 -= fq;
//			m17 -= fq;
//
//			// q=11
//			fq = dist[12*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx += fq;
//			m4 += fq;
//			jz += fq;
//			m8 += fq;
//			m9 += fq;
//			m10 += fq;
//			m11 -= fq;
//			m12 -= fq;
//			m15 = fq;
//			m16 -= fq;
//			m18 = fq;
//
//			// q=12
//			fq = dist[11*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx -= fq;
//			m4 -= fq;
//			jz -= fq;
//			m8 -= fq;
//			m9 += fq;
//			m10 += fq;
//			m11 -= fq;
//			m12 -= fq;
//			m15 += fq;
//			m16 += fq;
//			m18 -= fq;
//
//			// q=13
//			fq = dist[14*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx += fq;
//			m4 += fq;
//			jz -= fq;
//			m8 -= fq;
//			m9 += fq;
//			m10 += fq;
//			m11 -= fq;
//			m12 -= fq;
//			m15 -= fq;
//			m16 -= fq;
//			m18 -= fq;
//
//			// q=14
//			fq = dist[13*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jx -= fq;
//			m4 -= fq;
//			jz += fq;
//			m8 += fq;
//			m9 += fq;
//			m10 += fq;
//			m11 -= fq;
//			m12 -= fq;
//			m15 -= fq;
//			m16 += fq;
//			m18 += fq;
//
//			// q=15
//			fq = dist[16*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jy += fq;
//			m6 += fq;
//			jz += fq;
//			m8 += fq;
//			m9 -= 2.0*fq;
//			m10 -= 2.0*fq;
//			m14 = fq;
//			m17 += fq;
//			m18 -= fq;
//
//			// q=16
//			fq = dist[15*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jy -= fq;
//			m6 -= fq;
//			jz -= fq;
//			m8 -= fq;
//			m9 -= 2.0*fq;
//			m10 -= 2.0*fq;
//			m14 += fq;
//			m17 -= fq;
//			m18 += fq;
//
//			// q=17
//			fq = dist[18*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jy += fq;
//			m6 += fq;
//			jz -= fq;
//			m8 -= fq;
//			m9 -= 2.0*fq;
//			m10 -= 2.0*fq;
//			m14 -= fq;
//			m17 += fq;
//			m18 += fq;
//
//			// q=18
//			fq = dist[17*Np+n];
//			rho += fq;
//			m1 += 8.0*fq;
//			m2 += fq;
//			jy -= fq;
//			m6 -= fq;
//			jz += fq;
//			m8 += fq;
//			m9 -= 2.0*fq;
//			m10 -= 2.0*fq;
//			m14 -= fq;
//			m17 -= fq;
//			m18 -= fq;
//
//            // Compute greyscale related parameters
//            c0 = 0.5*(1.0+porosity*0.5*mu_eff/perm);
//            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
//            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
//            c1 = porosity*0.5*GeoFun/sqrt(perm);
//            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes
//
//            vx = jx/rho0+0.5*(porosity*Gx);
//            vy = jy/rho0+0.5*(porosity*Gy);
//            vz = jz/rho0+0.5*(porosity*Gz);
//            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
//            ux = vx/(c0+sqrt(c0*c0+c1*v_mag));
//            uy = vy/(c0+sqrt(c0*c0+c1*v_mag));
//            uz = vz/(c0+sqrt(c0*c0+c1*v_mag));
//            u_mag=sqrt(ux*ux+uy*uy+uz*uz);
//
//            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
//            Fx = rho0*(-porosity*mu_eff/perm*ux - porosity*GeoFun/sqrt(perm)*u_mag*ux + porosity*Gx);
//            Fy = rho0*(-porosity*mu_eff/perm*uy - porosity*GeoFun/sqrt(perm)*u_mag*uy + porosity*Gy);
//            Fz = rho0*(-porosity*mu_eff/perm*uz - porosity*GeoFun/sqrt(perm)*u_mag*uz + porosity*Gz);
//            if (porosity==1.0){
//                Fx=rho0*(Gx);
//                Fy=rho0*(Gy);
//                Fz=rho0*(Gz);
//            }
//
//			// write the velocity 
//			Velocity[n] = ux;
//			Velocity[Np+n] = uy;
//			Velocity[2*Np+n] = uz;
//
//			//........................................................................
//			//..............carry out relaxation process..............................
//			//..........Toelke, Fruediger et. al. 2006................................
//			if (C == 0.0)	nx = ny = nz = 0.0;
//			m1 = m1 + rlx_setA*((19*(ux*ux+uy*uy+uz*uz)*rho0/porosity - 11*rho) -19*alpha*C - m1);
//			m2 = m2 + rlx_setA*((3*rho - 5.5*(ux*ux+uy*uy+uz*uz)*rho0/porosity)- m2);
//            jx = jx + Fx;
//			m4 = m4 + rlx_setB*((-0.6666666666666666*ux*rho0)- m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//			m6 = m6 + rlx_setB*((-0.6666666666666666*uy*rho0)- m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//			m8 = m8 + rlx_setB*((-0.6666666666666666*uz*rho0)- m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//			m9 = m9 + rlx_setA*(((2*ux*ux-uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
//			m10 = m10 + rlx_setA*( - m10);
//            //m10 = m10 + rlx_setA*(-0.5*rho0*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10);
//			m11 = m11 + rlx_setA*(((uy*uy-uz*uz)*rho0/porosity) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
//			m12 = m12 + rlx_setA*( - m12);
//            //m12 = m12 + rlx_setA*(-0.5*(rho0*(uy*uy-uz*uz)/porosity)- m12);
//			m13 = m13 + rlx_setA*( (ux*uy*rho0/porosity) + 0.5*alpha*C*nx*ny - m13);
//			m14 = m14 + rlx_setA*( (uy*uz*rho0/porosity) + 0.5*alpha*C*ny*nz - m14);
//			m15 = m15 + rlx_setA*( (ux*uz*rho0/porosity) + 0.5*alpha*C*nx*nz - m15);
//			m16 = m16 + rlx_setB*( - m16);
//			m17 = m17 + rlx_setB*( - m17);
//			m18 = m18 + rlx_setB*( - m18);
//
//			//.................inverse transformation......................................................
//			// q=0
//			fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
//			dist[n] = fq;
//
//			// q = 1
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10);
//			dist[1*Np+n] = fq;
//
//			// q=2
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10);
//			dist[2*Np+n] = fq;
//
//			// q = 3
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
//			dist[3*Np+n] = fq;
//
//			// q = 4
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12);
//			dist[4*Np+n] = fq;
//
//			// q = 5
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
//			dist[5*Np+n] = fq;
//
//			// q = 6
//			fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11);
//			dist[6*Np+n] = fq;
//
//			// q = 7
//			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)+
//					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12+0.25*m13+0.125*(m16-m17);
//			dist[7*Np+n] = fq;
//
//
//			// q = 8
//			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
//					+mrt_V12*m12+0.25*m13+0.125*(m17-m16);
//			dist[8*Np+n] = fq;
//
//			// q = 9
//			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)+
//					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13+0.125*(m16+m17);
//			dist[9*Np+n] = fq;
//
//			// q = 10
//			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)+
//					mrt_V7*m9+mrt_V11*m10+mrt_V8*m11+mrt_V12*m12-0.25*m13-0.125*(m16+m17);
//			dist[10*Np+n] = fq;
//
//
//			// q = 11
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
//					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
//					-mrt_V12*m12+0.25*m15+0.125*(m18-m16);
//			dist[11*Np+n] = fq;
//
//			// q = 12
//			fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)+
//					mrt_V7*m9+mrt_V11*m10-mrt_V8*m11-mrt_V12*m12+0.25*m15+0.125*(m16-m18);
//			dist[12*Np+n] = fq;
//
//			// q = 13
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
//					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
//					-mrt_V12*m12-0.25*m15-0.125*(m16+m18);
//			dist[13*Np+n] = fq;
//
//			// q= 14
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
//					+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
//					-mrt_V12*m12-0.25*m15+0.125*(m16+m18);
//
//			dist[14*Np+n] = fq;
//
//			// q = 15
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
//					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18);
//			dist[15*Np+n] = fq;
//
//			// q = 16
//			fq =  mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
//					-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17);
//			dist[16*Np+n] = fq;
//
//
//			// q = 17
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
//					-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18);
//			dist[17*Np+n] = fq;
//
//			// q = 18
//			fq = mrt_V1*rho+mrt_V9*m1
//					+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
//					-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18);
//			dist[18*Np+n] = fq;
//			//........................................................................
//
//			// Instantiate mass transport distributions
//			// Stationary value - distribution 0
//			nAB = 1.0/(nA+nB);
//			Aq[n] = 0.3333333333333333*nA;
//			Bq[n] = 0.3333333333333333*nB;
//
//			//...............................................
//			// q = 0,2,4
//			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
//			delta = beta*nA*nB*nAB*0.1111111111111111*nx;
//			if (!(nA*nB*nAB>0)) delta=0;
//			a1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
//			b1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
//			a2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
//			b2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;
//
//			Aq[1*Np+n] = a1;
//			Bq[1*Np+n] = b1;
//			Aq[2*Np+n] = a2;
//			Bq[2*Np+n] = b2;
//
//			//...............................................
//			// q = 2
//			// Cq = {0,1,0}
//			delta = beta*nA*nB*nAB*0.1111111111111111*ny;
//			if (!(nA*nB*nAB>0)) delta=0;
//			a1 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
//			b1 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
//			a2 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
//			b2 = nB*(0.1111111111111111*(1-4.5*uy))+delta;
//
//			Aq[3*Np+n] = a1;
//			Bq[3*Np+n] = b1;
//			Aq[4*Np+n] = a2;
//			Bq[4*Np+n] = b2;
//			//...............................................
//			// q = 4
//			// Cq = {0,0,1}
//			delta = beta*nA*nB*nAB*0.1111111111111111*nz;
//			if (!(nA*nB*nAB>0)) delta=0;
//			a1 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
//			b1 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
//			a2 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
//			b2 = nB*(0.1111111111111111*(1-4.5*uz))+delta;
//
//			Aq[5*Np+n] = a1;
//			Bq[5*Np+n] = b1;
//			Aq[6*Np+n] = a2;
//			Bq[6*Np+n] = b2;
//			//...............................................
//
//		}
//	}
//}

//__global__ void dvc_ScaLBL_D3Q19_GreyscaleColor_Init(double *dist, double *Porosity, int Np)
//{
//	int n;
//	int S = Np/NBLOCKS/NTHREADS + 1;
//    double porosity;
//	for (int s=0; s<S; s++){
//		//........Get 1-D index for this thread....................
//		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
//		if (n<Np ){
//            porosity = Porosity[n];
//            if (porosity==0.0) porosity=1.f;
//			dist[n] = 0.3333333333333333/porosity;
//			dist[Np+n] = 0.055555555555555555/porosity;		//double(100*n)+1.f;
//			dist[2*Np+n] = 0.055555555555555555/porosity;	//double(100*n)+2.f;
//			dist[3*Np+n] = 0.055555555555555555/porosity;	//double(100*n)+3.f;
//			dist[4*Np+n] = 0.055555555555555555/porosity;	//double(100*n)+4.f;
//			dist[5*Np+n] = 0.055555555555555555/porosity;	//double(100*n)+5.f;
//			dist[6*Np+n] = 0.055555555555555555/porosity;	//double(100*n)+6.f;
//			dist[7*Np+n] = 0.0277777777777778/porosity;   //double(100*n)+7.f;
//			dist[8*Np+n] = 0.0277777777777778/porosity;   //double(100*n)+8.f;
//			dist[9*Np+n] = 0.0277777777777778/porosity;   //double(100*n)+9.f;
//			dist[10*Np+n] = 0.0277777777777778/porosity;  //double(100*n)+10.f;
//			dist[11*Np+n] = 0.0277777777777778/porosity;  //double(100*n)+11.f;
//			dist[12*Np+n] = 0.0277777777777778/porosity;  //double(100*n)+12.f;
//			dist[13*Np+n] = 0.0277777777777778/porosity;  //double(100*n)+13.f;
//			dist[14*Np+n] = 0.0277777777777778/porosity;  //double(100*n)+14.f;
//			dist[15*Np+n] = 0.0277777777777778/porosity;  //double(100*n)+15.f;
//			dist[16*Np+n] = 0.0277777777777778/porosity;  //double(100*n)+16.f;
//			dist[17*Np+n] = 0.0277777777777778/porosity;  //double(100*n)+17.f;
//			dist[18*Np+n] = 0.0277777777777778/porosity;  //double(100*n)+18.f;
//		}
//	}
//}


//extern "C" void ScaLBL_D3Q19_GreyscaleColor_Init(double *dist, double *Porosity, int Np){
//	dvc_ScaLBL_D3Q19_GreyscaleColor_Init<<<NBLOCKS,NTHREADS >>>(dist,Porosity,Np);
//	cudaError_t err = cudaGetLastError();
//	if (cudaSuccess != err){
//		printf("CUDA error in ScaLBL_D3Q19_GreyscaleColor_Init: %s \n",cudaGetErrorString(err));
//	}
//}

//Model-1 & 4
extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColor(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
        double *Phi,double *GreySolidGrad, double *Poros,double *Perm,double *Vel, 
        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np){

	//cudaProfilerStart();
	//cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor, cudaFuncCachePreferL1);

	dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor<<<NBLOCKS,NTHREADS >>>(Map, dist, Aq, Bq, Den, Phi, GreySolidGrad, Poros, Perm, Vel, 
            rhoA, rhoB, tauA, tauB, tauA_eff, tauB_eff, alpha, beta, Fx, Fy, Fz, strideY, strideZ, start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_GreyscaleColor: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();

}

//Model-1 & 4
extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColor(int *d_neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *GreySolidGrad, double *Poros,double *Perm,double *Vel, 
        double rhoA, double rhoB, double tauA, double tauB, double tauA_eff,double tauB_eff, double alpha, double beta,
		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np){

	//cudaProfilerStart();
	//cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor, cudaFuncCachePreferL1);
	
	dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor<<<NBLOCKS,NTHREADS >>>(d_neighborList, Map, dist, Aq, Bq, Den, Phi,  GreySolidGrad, Poros, Perm,Vel, 
			rhoA, rhoB, tauA, tauB, tauA_eff, tauB_eff,alpha, beta, Fx, Fy, Fz, strideY, strideZ, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_GreyscaleColor: %s \n",cudaGetErrorString(err));
	}
	//cudaProfilerStop();
}

////Model-2&3
//extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleColor(int *Map, double *dist, double *Aq, double *Bq, double *Den, 
//        double *Phi,double *GreySolidGrad, double *Poros,double *Perm,double *Vel, 
//        double rhoA, double rhoB, double tauA, double tauB,double tauA_eff,double tauB_eff, double alpha, double beta,
//		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np){
//
//	//cudaProfilerStart();
//	//cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor, cudaFuncCachePreferL1);
//
//	dvc_ScaLBL_D3Q19_AAeven_GreyscaleColor<<<NBLOCKS,NTHREADS >>>(Map, dist, Aq, Bq, Den, Phi, GreySolidGrad, Poros, Perm, Vel, 
//            rhoA, rhoB, tauA, tauB, tauA_eff, tauB_eff, alpha, beta, Fx, Fy, Fz, strideY, strideZ, start, finish, Np);
//	cudaError_t err = cudaGetLastError();
//	if (cudaSuccess != err){
//		printf("CUDA error in ScaLBL_D3Q19_AAeven_GreyscaleColor: %s \n",cudaGetErrorString(err));
//	}
//	//cudaProfilerStop();
//
//}
//
////Model-2&3
//extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleColor(int *d_neighborList, int *Map, double *dist, double *Aq, double *Bq, double *Den, 
//		double *Phi, double *GreySolidGrad, double *Poros,double *Perm,double *Vel, 
//        double rhoA, double rhoB, double tauA, double tauB, double tauA_eff,double tauB_eff, double alpha, double beta,
//		double Fx, double Fy, double Fz, int strideY, int strideZ, int start, int finish, int Np){
//
//	//cudaProfilerStart();
//	//cudaFuncSetCacheConfig(dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor, cudaFuncCachePreferL1);
//	
//	dvc_ScaLBL_D3Q19_AAodd_GreyscaleColor<<<NBLOCKS,NTHREADS >>>(d_neighborList, Map, dist, Aq, Bq, Den, Phi,  GreySolidGrad, Poros, Perm,Vel, 
//			rhoA, rhoB, tauA, tauB, tauA_eff, tauB_eff,alpha, beta, Fx, Fy, Fz, strideY, strideZ, start, finish, Np);
//
//	cudaError_t err = cudaGetLastError();
//	if (cudaSuccess != err){
//		printf("CUDA error in ScaLBL_D3Q19_AAodd_GreyscaleColor: %s \n",cudaGetErrorString(err));
//	}
//	//cudaProfilerStop();
//}
