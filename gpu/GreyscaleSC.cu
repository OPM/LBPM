#include <stdio.h>

#define NBLOCKS 1024
#define NTHREADS 256

__global__ void dvc_ScaLBL_D3Q19_AAodd_GreyscaleSC(int *neighborList, double *distA, double *distB, double *Den, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

	int n, nread;
	double vx,vy,vz,v_mag;
    double ux_A,uy_A,uz_A,ux_B,uy_B,uz_B,u_mag;
    double ux,uy,uz;
	// conserved momemnts
	double jxA,jyA,jzA;
	double jxB,jyB,jzB;
    double rhoA,rhoB;
    double rhoA_next,rhoB_next;
	// non-conserved moments
	double m1A,m2A,m4A,m6A,m8A,m9A,m10A,m11A,m12A,m13A,m14A,m15A,m16A,m17A,m18A;
	double m1B,m2B,m4B,m6B,m8B,m9B,m10B,m11B,m12B,m13B,m14B,m15B,m16B,m17B,m18B;
    double fq;
	//double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double permA,permB;//effective relative perm
    double c0, c1; //Guo's model parameters
    double muA_eff = (tauA_eff-0.5)/3.0;//kinematic viscosity
    double muB_eff = (tauB_eff-0.5)/3.0;//kinematic viscosity
    double FxA, FyA, FzA;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
    double FxB, FyB, FzB;
    double rlx_setA,rlx_setB;
    double rhoA_gradx,rhoA_grady,rhoA_gradz;
    double rhoB_gradx,rhoB_grady,rhoB_gradz;
    double GffA_x,GffA_y,GffA_z;
    double GfsA_x,GfsA_y,GfsA_z;
    double GffB_x,GffB_y,GffB_z;
    double GfsB_x,GfsB_y,GfsB_z;

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
	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){		

            // Load common parameters shared by two fluid components
            rhoA = Den[n];
            rhoB = Den[n+Np];
            porosity = Poros[n];
            perm = Perm[n];
            permA = perm*rhoA/(rhoA+rhoB);//effective relative perm
            permB = perm*rhoB/(rhoA+rhoB);
			rhoA_gradx = DenGradA[n+0*Np];
			rhoA_grady = DenGradA[n+1*Np];
			rhoA_gradz = DenGradA[n+2*Np];
			rhoB_gradx = DenGradB[n+0*Np];
			rhoB_grady = DenGradB[n+1*Np];
			rhoB_gradz = DenGradB[n+2*Np];

            // ------------------- Fluid component A ---------------------------------//
            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
			// q=0
			fq = distA[n];
            rhoA_next = fq;
			m1A  = -30.0*fq;
			m2A  = 12.0*fq;

			// q=1
			nread = neighborList[n]; // neighbor 2 
			fq = distA[nread]; // reading the f1 data into register fq		
            rhoA_next += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jxA = fq;
			m4A = -4.0*fq;
			m9A = 2.0*fq;
			m10A = -4.0*fq;

			// q=2
			nread = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
			fq = distA[nread];  // reading the f2 data into register fq
            rhoA_next += fq;
			m1A -= 11.0*(fq);
			m2A -= 4.0*(fq);
			jxA -= fq;
			m4A += 4.0*(fq);
			m9A += 2.0*(fq);
			m10A -= 4.0*(fq);

			// q=3
			nread = neighborList[n+2*Np]; // neighbor 4
			fq = distA[nread];
            rhoA_next += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jyA = fq;
			m6A = -4.0*fq;
			m9A -= fq;
			m10A += 2.0*fq;
			m11A = fq;
			m12A = -2.0*fq;

			// q = 4
			nread = neighborList[n+3*Np]; // neighbor 3
			fq = distA[nread];
            rhoA_next += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jyA -= fq;
			m6A += 4.0*fq;
			m9A -= fq;
			m10A += 2.0*fq;
			m11A += fq;
			m12A -= 2.0*fq;

			// q=5
			nread = neighborList[n+4*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jzA = fq;
			m8A = -4.0*fq;
			m9A -= fq;
			m10A += 2.0*fq;
			m11A -= fq;
			m12A += 2.0*fq;


			// q = 6
			nread = neighborList[n+5*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jzA -= fq;
			m8A += 4.0*fq;
			m9A -= fq;
			m10A += 2.0*fq;
			m11A -= fq;
			m12A += 2.0*fq;

			// q=7
			nread = neighborList[n+6*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA += fq;
			m4A += fq;
			jyA += fq;
			m6A += fq;
			m9A  += fq;
			m10A += fq;
			m11A += fq;
			m12A += fq;
			m13A = fq;
			m16A = fq;
			m17A = -fq;

			// q = 8
			nread = neighborList[n+7*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA -= fq;
			m4A -= fq;
			jyA -= fq;
			m6A -= fq;
			m9A += fq;
			m10A += fq;
			m11A += fq;
			m12A += fq;
			m13A += fq;
			m16A -= fq;
			m17A += fq;

			// q=9
			nread = neighborList[n+8*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA += fq;
			m4A += fq;
			jyA -= fq;
			m6A -= fq;
			m9A += fq;
			m10A += fq;
			m11A += fq;
			m12A += fq;
			m13A -= fq;
			m16A += fq;
			m17A += fq;

			// q = 10
			nread = neighborList[n+9*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA -= fq;
			m4A -= fq;
			jyA += fq;
			m6A += fq;
			m9A += fq;
			m10A += fq;
			m11A += fq;
			m12A += fq;
			m13A -= fq;
			m16A -= fq;
			m17A -= fq;

			// q=11
			nread = neighborList[n+10*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA += fq;
			m4A += fq;
			jzA += fq;
			m8A += fq;
			m9A += fq;
			m10A += fq;
			m11A -= fq;
			m12A -= fq;
			m15A = fq;
			m16A -= fq;
			m18A = fq;

			// q=12
			nread = neighborList[n+11*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA -= fq;
			m4A -= fq;
			jzA -= fq;
			m8A -= fq;
			m9A += fq;
			m10A += fq;
			m11A -= fq;
			m12A -= fq;
			m15A += fq;
			m16A += fq;
			m18A -= fq;

			// q=13
			nread = neighborList[n+12*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA += fq;
			m4A += fq;
			jzA -= fq;
			m8A -= fq;
			m9A += fq;
			m10A += fq;
			m11A -= fq;
			m12A -= fq;
			m15A -= fq;
			m16A -= fq;
			m18A -= fq;

			// q=14
			nread = neighborList[n+13*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA -= fq;
			m4A -= fq;
			jzA += fq;
			m8A += fq;
			m9A += fq;
			m10A += fq;
			m11A -= fq;
			m12A -= fq;
			m15A -= fq;
			m16A += fq;
			m18A += fq;

			// q=15
			nread = neighborList[n+14*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jyA += fq;
			m6A += fq;
			jzA += fq;
			m8A += fq;
			m9A -= 2.0*fq;
			m10A -= 2.0*fq;
			m14A = fq;
			m17A += fq;
			m18A -= fq;

			// q=16
			nread = neighborList[n+15*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jyA -= fq;
			m6A -= fq;
			jzA -= fq;
			m8A -= fq;
			m9A -= 2.0*fq;
			m10A -= 2.0*fq;
			m14A += fq;
			m17A -= fq;
			m18A += fq;

			// q=17
			nread = neighborList[n+16*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jyA += fq;
			m6A += fq;
			jzA -= fq;
			m8A -= fq;
			m9A -= 2.0*fq;
			m10A -= 2.0*fq;
			m14A -= fq;
			m17A += fq;
			m18A += fq;

			// q=18
			nread = neighborList[n+17*Np];
			fq = distA[nread];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jyA -= fq;
			m6A -= fq;
			jzA += fq;
			m8A += fq;
			m9A -= 2.0*fq;
			m10A -= 2.0*fq;
			m14A -= fq;
			m17A -= fq;
			m18A -= fq;
            //---------------------------------------------------------------------//

            // ------------------- Fluid component B ---------------------------------//
            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
			// q=0
			fq = distB[n];
            rhoB_next = fq;
			m1B  = -30.0*fq;
			m2B  = 12.0*fq;

			// q=1
			nread = neighborList[n]; // neighbor 2 
			fq = distB[nread]; // reading the f1 data into register fq		
            rhoB_next += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jxB = fq;
			m4B = -4.0*fq;
			m9B = 2.0*fq;
			m10B = -4.0*fq;

			// q=2
			nread = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
			fq = distB[nread];  // reading the f2 data into register fq
            rhoB_next += fq;
			m1B -= 11.0*(fq);
			m2B -= 4.0*(fq);
			jxB -= fq;
			m4B += 4.0*(fq);
			m9B += 2.0*(fq);
			m10B -= 4.0*(fq);

			// q=3
			nread = neighborList[n+2*Np]; // neighbor 4
			fq = distB[nread];
            rhoB_next += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jyB = fq;
			m6B = -4.0*fq;
			m9B -= fq;
			m10B += 2.0*fq;
			m11B = fq;
			m12B = -2.0*fq;

			// q = 4
			nread = neighborList[n+3*Np]; // neighbor 3
			fq = distB[nread];
            rhoB_next += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jyB -= fq;
			m6B += 4.0*fq;
			m9B -= fq;
			m10B += 2.0*fq;
			m11B += fq;
			m12B -= 2.0*fq;

			// q=5
			nread = neighborList[n+4*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jzB = fq;
			m8B = -4.0*fq;
			m9B -= fq;
			m10B += 2.0*fq;
			m11B -= fq;
			m12B += 2.0*fq;


			// q = 6
			nread = neighborList[n+5*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jzB -= fq;
			m8B += 4.0*fq;
			m9B -= fq;
			m10B += 2.0*fq;
			m11B -= fq;
			m12B += 2.0*fq;

			// q=7
			nread = neighborList[n+6*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB += fq;
			m4B += fq;
			jyB += fq;
			m6B += fq;
			m9B  += fq;
			m10B += fq;
			m11B += fq;
			m12B += fq;
			m13B = fq;
			m16B = fq;
			m17B = -fq;

			// q = 8
			nread = neighborList[n+7*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB -= fq;
			m4B -= fq;
			jyB -= fq;
			m6B -= fq;
			m9B += fq;
			m10B += fq;
			m11B += fq;
			m12B += fq;
			m13B += fq;
			m16B -= fq;
			m17B += fq;

			// q=9
			nread = neighborList[n+8*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB += fq;
			m4B += fq;
			jyB -= fq;
			m6B -= fq;
			m9B += fq;
			m10B += fq;
			m11B += fq;
			m12B += fq;
			m13B -= fq;
			m16B += fq;
			m17B += fq;

			// q = 10
			nread = neighborList[n+9*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB -= fq;
			m4B -= fq;
			jyB += fq;
			m6B += fq;
			m9B += fq;
			m10B += fq;
			m11B += fq;
			m12B += fq;
			m13B -= fq;
			m16B -= fq;
			m17B -= fq;

			// q=11
			nread = neighborList[n+10*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB += fq;
			m4B += fq;
			jzB += fq;
			m8B += fq;
			m9B += fq;
			m10B += fq;
			m11B -= fq;
			m12B -= fq;
			m15B = fq;
			m16B -= fq;
			m18B = fq;

			// q=12
			nread = neighborList[n+11*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB -= fq;
			m4B -= fq;
			jzB -= fq;
			m8B -= fq;
			m9B += fq;
			m10B += fq;
			m11B -= fq;
			m12B -= fq;
			m15B += fq;
			m16B += fq;
			m18B -= fq;

			// q=13
			nread = neighborList[n+12*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB += fq;
			m4B += fq;
			jzB -= fq;
			m8B -= fq;
			m9B += fq;
			m10B += fq;
			m11B -= fq;
			m12B -= fq;
			m15B -= fq;
			m16B -= fq;
			m18B -= fq;

			// q=14
			nread = neighborList[n+13*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB -= fq;
			m4B -= fq;
			jzB += fq;
			m8B += fq;
			m9B += fq;
			m10B += fq;
			m11B -= fq;
			m12B -= fq;
			m15B -= fq;
			m16B += fq;
			m18B += fq;

			// q=15
			nread = neighborList[n+14*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jyB += fq;
			m6B += fq;
			jzB += fq;
			m8B += fq;
			m9B -= 2.0*fq;
			m10B -= 2.0*fq;
			m14B = fq;
			m17B += fq;
			m18B -= fq;

			// q=16
			nread = neighborList[n+15*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jyB -= fq;
			m6B -= fq;
			jzB -= fq;
			m8B -= fq;
			m9B -= 2.0*fq;
			m10B -= 2.0*fq;
			m14B += fq;
			m17B -= fq;
			m18B += fq;

			// q=17
			nread = neighborList[n+16*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jyB += fq;
			m6B += fq;
			jzB -= fq;
			m8B -= fq;
			m9B -= 2.0*fq;
			m10B -= 2.0*fq;
			m14B -= fq;
			m17B += fq;
			m18B += fq;

			// q=18
			nread = neighborList[n+17*Np];
			fq = distB[nread];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jyB -= fq;
			m6B -= fq;
			jzB += fq;
			m8B += fq;
			m9B -= 2.0*fq;
			m10B -= 2.0*fq;
			m14B -= fq;
			m17B -= fq;
			m18B -= fq;
            //---------------------------------------------------------------------//


            // Compute SC fluid-fluid interaction force
            GffA_x = -Gsc*rhoB_gradx;
            GffA_y = -Gsc*rhoB_grady;
            GffA_z = -Gsc*rhoB_gradz;
            GffB_x = -Gsc*rhoA_gradx;
            GffB_y = -Gsc*rhoA_grady;
            GffB_z = -Gsc*rhoA_gradz;
            // Compute SC fluid-solid force
            GfsA_x = SolidForceA[n+0*Np];    
            GfsA_y = SolidForceA[n+1*Np];    
            GfsA_z = SolidForceA[n+2*Np];    
            GfsB_x = SolidForceB[n+0*Np];    
            GfsB_y = SolidForceB[n+1*Np];    
            GfsB_z = SolidForceB[n+2*Np];    

            // Compute greyscale related parameters
            // ------------------- Fluid Component A -----------------------//
            c0 = 0.5*(1.0+porosity*0.5*muA_eff/permA);
            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
            c1 = porosity*0.5*GeoFun/sqrt(permA);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jxA/rhoA+0.5*(porosity*Gx+GffA_x+GfsA_x);
            vy = jyA/rhoA+0.5*(porosity*Gy+GffA_y+GfsA_y);
            vz = jzA/rhoA+0.5*(porosity*Gz+GffA_z+GfsA_z);
            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
            ux_A = vx/(c0+sqrt(c0*c0+c1*v_mag));
            uy_A = vy/(c0+sqrt(c0*c0+c1*v_mag));
            uz_A = vz/(c0+sqrt(c0*c0+c1*v_mag));
            u_mag=sqrt(ux_A*ux_A+uy_A*uy_A+uz_A*uz_A);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            FxA = rhoA*(-porosity*muA_eff/permA*ux_A - porosity*GeoFun/sqrt(permA)*u_mag*ux_A + porosity*Gx + GffA_x + GfsA_x);
            FyA = rhoA*(-porosity*muA_eff/permA*uy_A - porosity*GeoFun/sqrt(permA)*u_mag*uy_A + porosity*Gy + GffA_y + GfsA_y);
            FzA = rhoA*(-porosity*muA_eff/permA*uz_A - porosity*GeoFun/sqrt(permA)*u_mag*uz_A + porosity*Gz + GffA_z + GfsA_z);
            if (porosity==1.0){
                FxA=rhoA*(Gx + GffA_x + GfsA_x);
                FyA=rhoA*(Gy + GffA_y + GfsA_y);
                FzA=rhoA*(Gz + GffA_z + GfsA_z);
            }
            // ------------------- Fluid Component B -----------------------//
            // Compute greyscale related parameters
            c0 = 0.5*(1.0+porosity*0.5*muB_eff/permB);
            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
            c1 = porosity*0.5*GeoFun/sqrt(permB);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jxB/rhoB+0.5*(porosity*Gx+GffB_x+GfsB_x);
            vy = jyB/rhoB+0.5*(porosity*Gy+GffB_y+GfsB_y);
            vz = jzB/rhoB+0.5*(porosity*Gz+GffB_z+GfsB_z);
            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
            ux_B = vx/(c0+sqrt(c0*c0+c1*v_mag));
            uy_B = vy/(c0+sqrt(c0*c0+c1*v_mag));
            uz_B = vz/(c0+sqrt(c0*c0+c1*v_mag));
            u_mag=sqrt(ux_B*ux_B+uy_B*uy_B+uz_B*uz_B);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            FxB = rhoB*(-porosity*muB_eff/permB*ux_B - porosity*GeoFun/sqrt(permB)*u_mag*ux_B + porosity*Gx + GffB_x + GfsB_x);
            FyB = rhoB*(-porosity*muB_eff/permB*uy_B - porosity*GeoFun/sqrt(permB)*u_mag*uy_B + porosity*Gy + GffB_y + GfsB_y);
            FzB = rhoB*(-porosity*muB_eff/permB*uz_B - porosity*GeoFun/sqrt(permB)*u_mag*uz_B + porosity*Gz + GffB_z + GfsB_z);
            if (porosity==1.0){
                FxB=rhoB*(Gx + GffB_x + GfsB_x);
                FyB=rhoB*(Gy + GffB_y + GfsB_y);
                FzB=rhoB*(Gz + GffB_z + GfsB_z);
            }

            // Calculate barycentric velocity of the fluid mixture
            ux = (rhoA*ux_A+rhoB*ux_B)/(rhoA+rhoB);
            uy = (rhoA*uy_A+rhoB*uy_B)/(rhoA+rhoB);
            uz = (rhoA*uz_A+rhoB*uz_B)/(rhoA+rhoB);

//            //..............carry out relaxation process...............................................
//            m1 = m1 + rlx_setA*((-30*Den+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
//                    + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m2 = m2 + rlx_setA*((12*Den - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
//                    + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
//            jx = jx + Fx;
//            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*Den) - m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*Den) - m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*Den) - m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//            m9 = m9 + rlx_setA*((Den*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
//                    + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
//            m10 = m10 + rlx_setA*(-0.5*Den*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
//                      + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m11 = m11 + rlx_setA*((Den*(uy*uy-uz*uz)/porosity) - m11)
//                      + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
//            m12 = m12 + rlx_setA*(-0.5*(Den*(uy*uy-uz*uz)/porosity)- m12)
//                      + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
//            m13 = m13 + rlx_setA*((Den*ux*uy/porosity) - m13)
//                      + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
//            m14 = m14 + rlx_setA*((Den*uy*uz/porosity) - m14)
//                      + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
//            m15 = m15 + rlx_setA*((Den*ux*uz/porosity) - m15)
//                      + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
//            m16 = m16 + rlx_setB*( - m16);
//            m17 = m17 + rlx_setB*( - m17);
//            m18 = m18 + rlx_setB*( - m18);
//            //.......................................................................................................
           
            // ------------------- Fluid Component A -----------------------//
            rlx_setA = 1.0/tauA;
            rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
            //-------------------- MRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            //TODO need to incoporate porosity
			m1A = m1A + rlx_setA*((19*rhoA*(ux*ux+uy*uy+uz*uz) - 11*rhoA_next) - m1A)
                      + (1-0.5*rlx_setA)*38*(FxA*ux+FyA*uy+FzA*uz);
			m2A = m2A + rlx_setA*((3*rhoA_next - 5.5*rhoA*(ux*ux+uy*uy+uz*uz))- m2A)
                      + (1-0.5*rlx_setA)*11*(-FxA*ux-FyA*uy-FzA*uz);
            jxA = jxA + FxA;
			m4A = m4A + rlx_setB*((-0.6666666666666666*ux*rhoA)- m4A)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FxA);
            jyA = jyA + FyA;
			m6A = m6A + rlx_setB*((-0.6666666666666666*uy*rhoA)- m6A)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FyA);
            jzA = jzA + FzA;
			m8A = m8A + rlx_setB*((-0.6666666666666666*uz*rhoA)- m8A)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FzA);
			m9A = m9A + rlx_setA*((rhoA*(2*ux*ux-uy*uy-uz*uz)) - m9A)
                      + (1-0.5*rlx_setA)*(4*FxA*ux-2*FyA*uy-2*FzA*uz);
			m10A = m10A + rlx_setA*( - m10A)
                        + (1-0.5*rlx_setA)*(-2*FxA*ux+FyA*uy+FzA*uz);
			m11A = m11A + rlx_setA*((rhoA*(uy*uy-uz*uz)) - m11A)
                        + (1-0.5*rlx_setA)*(2*FyA*uy-2*FzA*uz);
			m12A = m12A + rlx_setA*( - m12A)
                        + (1-0.5*rlx_setA)*(-FyA*uy+FzA*uz);
			m13A = m13A + rlx_setA*( rhoA*(ux*uy) - m13A)
                        + (1-0.5*rlx_setA)*(FyA*ux+FxA*uy);
			m14A = m14A + rlx_setA*( rhoA*(uy*uz) - m14A)
                        + (1-0.5*rlx_setA)*(FzA*uy+FyA*uz);
			m15A = m15A + rlx_setA*( rhoA*(ux*uz) - m15A)
                        + (1-0.5*rlx_setA)*(FzA*ux+FxA*uz);
			m16A = m16A + rlx_setB*( - m16A);
			m17A = m17A + rlx_setB*( - m17A);
			m18A = m18A + rlx_setB*( - m18A);
            //.......................................................................................................


            // ------------------- Fluid Component A -----------------------//
            //.................inverse transformation......................................................
            // q=0
            fq = mrt_V1*rhoA_next-mrt_V2*m1A+mrt_V3*m2A;
            //f0 = mrt_V1*rhoA_next-mrt_V2*m1A+mrt_V3*m2A;
            distA[n] = fq;

            // q = 1
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jxA-m4A)+mrt_V6*(m9A-m10A);
            //f1 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jxA-m4A)+mrt_V6*(m9A-m10A);
            nread = neighborList[n+Np];
            distA[nread] = fq;

            // q=2
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m4A-jxA)+mrt_V6*(m9A-m10A);
            //f2 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m4A-jxA)+mrt_V6*(m9A-m10A);
            nread = neighborList[n];
            distA[nread] = fq;

            // q = 3
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jyA-m6A)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            //f3 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jyA-m6A)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            nread = neighborList[n+3*Np];
            distA[nread] = fq;

            // q = 4
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m6A-jyA)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            //f4 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m6A-jyA)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            nread = neighborList[n+2*Np];
            distA[nread] = fq;

            // q = 5
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jzA-m8A)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            //f5 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jzA-m8A)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            nread = neighborList[n+5*Np];
            distA[nread] = fq;

            // q = 6
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m8A-jzA)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            //f6 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m8A-jzA)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            nread = neighborList[n+4*Np];
            distA[nread] = fq;

            // q = 7
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jyA)+0.025*(m4A+m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m16A-m17A);
            //f7 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jyA)+0.025*(m4A+m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m16A-m17A);
            nread = neighborList[n+7*Np];
            distA[nread] = fq;

            // q = 8
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jyA)-0.025*(m4A+m6A) +mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m17A-m16A);
            //f8 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jyA)-0.025*(m4A+m6A) +mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m17A-m16A);
            nread = neighborList[n+6*Np];
            distA[nread] = fq;

            // q = 9
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jyA)+0.025*(m4A-m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A+0.125*(m16A+m17A);
            //f9 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jyA)+0.025*(m4A-m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A+0.125*(m16A+m17A);
            nread = neighborList[n+9*Np];
            distA[nread] = fq;

            // q = 10
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jxA)+0.025*(m6A-m4A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A-0.125*(m16A+m17A);
            //f10 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jxA)+0.025*(m6A-m4A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A-0.125*(m16A+m17A);
            nread = neighborList[n+8*Np];
            distA[nread] = fq;

            // q = 11
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jzA)+0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m18A-m16A);
            //f11 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jzA)+0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m18A-m16A);
            nread = neighborList[n+11*Np];
            distA[nread] = fq;

            // q = 12
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jzA)-0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m16A-m18A);
            //f12 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jzA)-0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m16A-m18A);
            nread = neighborList[n+10*Np];
            distA[nread]= fq;

            // q = 13
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jzA)+0.025*(m4A-m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A-0.125*(m16A+m18A);
            //f13 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jzA)+0.025*(m4A-m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A-0.125*(m16A+m18A);
            nread = neighborList[n+13*Np];
            distA[nread] = fq;

            // q= 14
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jxA)+0.025*(m8A-m4A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A+0.125*(m16A+m18A);
            //f14 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jxA)+0.025*(m8A-m4A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A+0.125*(m16A+m18A);
            nread = neighborList[n+12*Np];
            distA[nread] = fq;

            // q = 15
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA+jzA)+0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m17A-m18A);
            //f15 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA+jzA)+0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m17A-m18A);
            nread = neighborList[n+15*Np];
            distA[nread] = fq;

            // q = 16
            fq =  mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jyA+jzA)-0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m18A-m17A);
            //f16 =  mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jyA+jzA)-0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m18A-m17A);
            nread = neighborList[n+14*Np];
            distA[nread] = fq;

            // q = 17
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jzA)+0.025*(m6A-m8A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A+0.125*(m17A+m18A);
            //f17 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jzA)+0.025*(m6A-m8A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A+0.125*(m17A+m18A);
            nread = neighborList[n+17*Np];
            distA[nread] = fq;

            // q = 18
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jyA)+0.025*(m8A-m6A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A-0.125*(m17A+m18A);
            //f18 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jyA)+0.025*(m8A-m6A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A-0.125*(m17A+m18A);
            nread = neighborList[n+16*Np];
            distA[nread] = fq;
            //........................................................................


            //Den[n] = f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;


//            //..............carry out relaxation process...............................................
//            m1 = m1 + rlx_setA*((-30*Den+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
//                    + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m2 = m2 + rlx_setA*((12*Den - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
//                    + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
//            jx = jx + Fx;
//            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*Den) - m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*Den) - m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*Den) - m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//            m9 = m9 + rlx_setA*((Den*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
//                    + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
//            m10 = m10 + rlx_setA*(-0.5*Den*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
//                      + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m11 = m11 + rlx_setA*((Den*(uy*uy-uz*uz)/porosity) - m11)
//                      + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
//            m12 = m12 + rlx_setA*(-0.5*(Den*(uy*uy-uz*uz)/porosity)- m12)
//                      + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
//            m13 = m13 + rlx_setA*((Den*ux*uy/porosity) - m13)
//                      + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
//            m14 = m14 + rlx_setA*((Den*uy*uz/porosity) - m14)
//                      + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
//            m15 = m15 + rlx_setA*((Den*ux*uz/porosity) - m15)
//                      + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
//            m16 = m16 + rlx_setB*( - m16);
//            m17 = m17 + rlx_setB*( - m17);
//            m18 = m18 + rlx_setB*( - m18);
//            //.......................................................................................................
           
            // ------------------- Fluid Component B -----------------------//
            rlx_setA = 1.0/tauB;
            rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
            //-------------------- MRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            //TODO need to incoporate porosity
			m1B = m1B + rlx_setA*((19*rhoB*(ux*ux+uy*uy+uz*uz) - 11*rhoB_next) - m1B)
                      + (1-0.5*rlx_setA)*38*(FxB*ux+FyB*uy+FzB*uz);
			m2B = m2B + rlx_setA*((3*rhoB_next - 5.5*rhoB*(ux*ux+uy*uy+uz*uz))- m2B)
                      + (1-0.5*rlx_setA)*11*(-FxB*ux-FyB*uy-FzB*uz);
            jxB = jxB + FxB;
			m4B = m4B + rlx_setB*((-0.6666666666666666*ux*rhoB)- m4B)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FxB);
            jyB = jyB + FyB;
			m6B = m6B + rlx_setB*((-0.6666666666666666*uy*rhoB)- m6B)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FyB);
            jzB = jzB + FzB;
			m8B = m8B + rlx_setB*((-0.6666666666666666*uz*rhoB)- m8B)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FzB);
			m9B = m9B + rlx_setA*((rhoB*(2*ux*ux-uy*uy-uz*uz)) - m9B)
                      + (1-0.5*rlx_setA)*(4*FxB*ux-2*FyB*uy-2*FzB*uz);
			m10B = m10B + rlx_setA*( - m10B)
                        + (1-0.5*rlx_setA)*(-2*FxB*ux+FyB*uy+FzB*uz);
			m11B = m11B + rlx_setA*((rhoB*(uy*uy-uz*uz)) - m11B)
                        + (1-0.5*rlx_setA)*(2*FyB*uy-2*FzB*uz);
			m12B = m12B + rlx_setA*( - m12B)
                        + (1-0.5*rlx_setA)*(-FyB*uy+FzB*uz);
			m13B = m13B + rlx_setA*( rhoB*(ux*uy) - m13B)
                        + (1-0.5*rlx_setA)*(FyB*ux+FxB*uy);
			m14B = m14B + rlx_setA*( rhoB*(uy*uz) - m14B)
                        + (1-0.5*rlx_setA)*(FzB*uy+FyB*uz);
			m15B = m15B + rlx_setA*( rhoB*(ux*uz) - m15B)
                        + (1-0.5*rlx_setA)*(FzB*ux+FxB*uz);
			m16B = m16B + rlx_setB*( - m16B);
			m17B = m17B + rlx_setB*( - m17B);
			m18B = m18B + rlx_setB*( - m18B);
            //.......................................................................................................


            // ------------------- Fluid Component B -----------------------//
            //.................inverse transformation......................................................
            // q=0
            fq = mrt_V1*rhoB_next-mrt_V2*m1B+mrt_V3*m2B;
            //f0 = mrt_V1*rhoB_next-mrt_V2*m1B+mrt_V3*m2B;
            distB[n] = fq;

            // q = 1
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jxB-m4B)+mrt_V6*(m9B-m10B);
            //f1 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jxB-m4B)+mrt_V6*(m9B-m10B);
            nread = neighborList[n+Np];
            distB[nread] = fq;

            // q=2
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m4B-jxB)+mrt_V6*(m9B-m10B);
            //f2 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m4B-jxB)+mrt_V6*(m9B-m10B);
            nread = neighborList[n];
            distB[nread] = fq;

            // q = 3
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jyB-m6B)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            //f3 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jyB-m6B)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            nread = neighborList[n+3*Np];
            distB[nread] = fq;

            // q = 4
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m6B-jyB)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            //f4 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m6B-jyB)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            nread = neighborList[n+2*Np];
            distB[nread] = fq;

            // q = 5
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jzB-m8B)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            //f5 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jzB-m8B)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            nread = neighborList[n+5*Np];
            distB[nread] = fq;

            // q = 6
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m8B-jzB)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            //f6 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m8B-jzB)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            nread = neighborList[n+4*Np];
            distB[nread] = fq;

            // q = 7
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jyB)+0.025*(m4B+m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m16B-m17B);
            //f7 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jyB)+0.025*(m4B+m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m16B-m17B);
            nread = neighborList[n+7*Np];
            distB[nread] = fq;

            // q = 8
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jyB)-0.025*(m4B+m6B) +mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m17B-m16B);
            //f8 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jyB)-0.025*(m4B+m6B) +mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m17B-m16B);
            nread = neighborList[n+6*Np];
            distB[nread] = fq;

            // q = 9
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jyB)+0.025*(m4B-m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B+0.125*(m16B+m17B);
            //f9 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jyB)+0.025*(m4B-m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B+0.125*(m16B+m17B);
            nread = neighborList[n+9*Np];
            distB[nread] = fq;

            // q = 10
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jxB)+0.025*(m6B-m4B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B-0.125*(m16B+m17B);
            //f10 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jxB)+0.025*(m6B-m4B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B-0.125*(m16B+m17B);
            nread = neighborList[n+8*Np];
            distB[nread] = fq;

            // q = 11
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jzB)+0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m18B-m16B);
            //f11 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jzB)+0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m18B-m16B);
            nread = neighborList[n+11*Np];
            distB[nread] = fq;

            // q = 12
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jzB)-0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m16B-m18B);
            //f12 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jzB)-0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m16B-m18B);
            nread = neighborList[n+10*Np];
            distB[nread]= fq;

            // q = 13
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jzB)+0.025*(m4B-m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B-0.125*(m16B+m18B);
            //f13 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jzB)+0.025*(m4B-m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B-0.125*(m16B+m18B);
            nread = neighborList[n+13*Np];
            distB[nread] = fq;

            // q= 14
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jxB)+0.025*(m8B-m4B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B+0.125*(m16B+m18B);
            //f14 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jxB)+0.025*(m8B-m4B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B+0.125*(m16B+m18B);
            nread = neighborList[n+12*Np];
            distB[nread] = fq;

            // q = 15
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB+jzB)+0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m17B-m18B);
            //f15 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB+jzB)+0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m17B-m18B);
            nread = neighborList[n+15*Np];
            distB[nread] = fq;

            // q = 16
            fq =  mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jyB+jzB)-0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m18B-m17B);
            //f16 =  mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jyB+jzB)-0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m18B-m17B);
            nread = neighborList[n+14*Np];
            distB[nread] = fq;

            // q = 17
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jzB)+0.025*(m6B-m8B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B+0.125*(m17B+m18B);
            //f17 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jzB)+0.025*(m6B-m8B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B+0.125*(m17B+m18B);
            nread = neighborList[n+17*Np];
            distB[nread] = fq;

            // q = 18
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jyB)+0.025*(m8B-m6B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B-0.125*(m17B+m18B);
            //f18 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jyB)+0.025*(m8B-m6B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B-0.125*(m17B+m18B);
            nread = neighborList[n+16*Np];
            distB[nread] = fq;
            //........................................................................

            //Den[n+Np] = f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;

            //Update velocity on device
            Velocity[0*Np+n] = ux;
            Velocity[1*Np+n] = uy;
            Velocity[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = (rhoA+rhoB+Gsc*rhoA*rhoB)/3.0;
            //Update density
            //Den[n] = rhoA_next;
            //Den[n+Np] = rhoB_next;

		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_AAeven_GreyscaleSC(double *distA, double *distB, double *Den, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

	int n;
	double vx,vy,vz,v_mag;
    double ux_A,uy_A,uz_A,ux_B,uy_B,uz_B,u_mag;
    double ux,uy,uz;
	// conserved momemnts
	double jxA,jyA,jzA;
	double jxB,jyB,jzB;
    double rhoA,rhoB;
    double rhoA_next,rhoB_next;
	// non-conserved moments
	double m1A,m2A,m4A,m6A,m8A,m9A,m10A,m11A,m12A,m13A,m14A,m15A,m16A,m17A,m18A;
	double m1B,m2B,m4B,m6B,m8B,m9B,m10B,m11B,m12B,m13B,m14B,m15B,m16B,m17B,m18B;
    double fq;
	//double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double permA,permB;//effective relative perm
    double c0, c1; //Guo's model parameters
    double muA_eff = (tauA_eff-0.5)/3.0;//kinematic viscosity
    double muB_eff = (tauB_eff-0.5)/3.0;//kinematic viscosity
    double FxA, FyA, FzA;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
    double FxB, FyB, FzB;
    double rlx_setA,rlx_setB;
    double rhoA_gradx,rhoA_grady,rhoA_gradz;
    double rhoB_gradx,rhoB_grady,rhoB_gradz;
    double GffA_x,GffA_y,GffA_z;
    double GfsA_x,GfsA_y,GfsA_z;
    double GffB_x,GffB_y,GffB_z;
    double GfsB_x,GfsB_y,GfsB_z;

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
	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){		

            // Load common parameters shared by two fluid components
            rhoA = Den[n];
            rhoB = Den[n+Np];
            porosity = Poros[n];
            perm = Perm[n];
            permA = perm*rhoA/(rhoA+rhoB);//effective relative perm
            permB = perm*rhoB/(rhoA+rhoB);
			rhoA_gradx = DenGradA[n+0*Np];
			rhoA_grady = DenGradA[n+1*Np];
			rhoA_gradz = DenGradA[n+2*Np];
			rhoB_gradx = DenGradB[n+0*Np];
			rhoB_grady = DenGradB[n+1*Np];
			rhoB_gradz = DenGradB[n+2*Np];

            // ------------------- Fluid component A ---------------------------------//
            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
			// q=0
			fq = distA[n];
            rhoA_next = fq;
			m1A  = -30.0*fq;
			m2A  = 12.0*fq;

			// q=1
            fq = distA[2*Np+n];
            rhoA_next += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jxA = fq;
			m4A = -4.0*fq;
			m9A = 2.0*fq;
			m10A = -4.0*fq;

			// q=2
            fq = distA[1*Np+n];
            rhoA_next += fq;
			m1A -= 11.0*(fq);
			m2A -= 4.0*(fq);
			jxA -= fq;
			m4A += 4.0*(fq);
			m9A += 2.0*(fq);
			m10A -= 4.0*(fq);

			// q=3
            fq = distA[4*Np+n];
            rhoA_next += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jyA = fq;
			m6A = -4.0*fq;
			m9A -= fq;
			m10A += 2.0*fq;
			m11A = fq;
			m12A = -2.0*fq;

			// q = 4
            fq = distA[3*Np+n];
            rhoA_next += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jyA -= fq;
			m6A += 4.0*fq;
			m9A -= fq;
			m10A += 2.0*fq;
			m11A += fq;
			m12A -= 2.0*fq;

			// q=5
            fq = distA[6*Np+n];
            rhoA_next += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jzA = fq;
			m8A = -4.0*fq;
			m9A -= fq;
			m10A += 2.0*fq;
			m11A -= fq;
			m12A += 2.0*fq;


			// q = 6
            fq = distA[5*Np+n];
            rhoA_next += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jzA -= fq;
			m8A += 4.0*fq;
			m9A -= fq;
			m10A += 2.0*fq;
			m11A -= fq;
			m12A += 2.0*fq;

			// q=7
            fq = distA[8*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA += fq;
			m4A += fq;
			jyA += fq;
			m6A += fq;
			m9A  += fq;
			m10A += fq;
			m11A += fq;
			m12A += fq;
			m13A = fq;
			m16A = fq;
			m17A = -fq;

			// q = 8
            fq = distA[7*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA -= fq;
			m4A -= fq;
			jyA -= fq;
			m6A -= fq;
			m9A += fq;
			m10A += fq;
			m11A += fq;
			m12A += fq;
			m13A += fq;
			m16A -= fq;
			m17A += fq;

			// q=9
            fq = distA[10*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA += fq;
			m4A += fq;
			jyA -= fq;
			m6A -= fq;
			m9A += fq;
			m10A += fq;
			m11A += fq;
			m12A += fq;
			m13A -= fq;
			m16A += fq;
			m17A += fq;

			// q = 10
            fq = distA[9*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA -= fq;
			m4A -= fq;
			jyA += fq;
			m6A += fq;
			m9A += fq;
			m10A += fq;
			m11A += fq;
			m12A += fq;
			m13A -= fq;
			m16A -= fq;
			m17A -= fq;

			// q=11
            fq = distA[12*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA += fq;
			m4A += fq;
			jzA += fq;
			m8A += fq;
			m9A += fq;
			m10A += fq;
			m11A -= fq;
			m12A -= fq;
			m15A = fq;
			m16A -= fq;
			m18A = fq;

			// q=12
            fq = distA[11*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA -= fq;
			m4A -= fq;
			jzA -= fq;
			m8A -= fq;
			m9A += fq;
			m10A += fq;
			m11A -= fq;
			m12A -= fq;
			m15A += fq;
			m16A += fq;
			m18A -= fq;

			// q=13
            fq = distA[14*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA += fq;
			m4A += fq;
			jzA -= fq;
			m8A -= fq;
			m9A += fq;
			m10A += fq;
			m11A -= fq;
			m12A -= fq;
			m15A -= fq;
			m16A -= fq;
			m18A -= fq;

			// q=14
            fq = distA[13*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jxA -= fq;
			m4A -= fq;
			jzA += fq;
			m8A += fq;
			m9A += fq;
			m10A += fq;
			m11A -= fq;
			m12A -= fq;
			m15A -= fq;
			m16A += fq;
			m18A += fq;

			// q=15
            fq = distA[16*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jyA += fq;
			m6A += fq;
			jzA += fq;
			m8A += fq;
			m9A -= 2.0*fq;
			m10A -= 2.0*fq;
			m14A = fq;
			m17A += fq;
			m18A -= fq;

			// q=16
            fq = distA[15*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jyA -= fq;
			m6A -= fq;
			jzA -= fq;
			m8A -= fq;
			m9A -= 2.0*fq;
			m10A -= 2.0*fq;
			m14A += fq;
			m17A -= fq;
			m18A += fq;

			// q=17
            fq = distA[18*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jyA += fq;
			m6A += fq;
			jzA -= fq;
			m8A -= fq;
			m9A -= 2.0*fq;
			m10A -= 2.0*fq;
			m14A -= fq;
			m17A += fq;
			m18A += fq;

			// q=18
            fq = distA[17*Np+n];
            rhoA_next += fq;
			m1A += 8.0*fq;
			m2A += fq;
			jyA -= fq;
			m6A -= fq;
			jzA += fq;
			m8A += fq;
			m9A -= 2.0*fq;
			m10A -= 2.0*fq;
			m14A -= fq;
			m17A -= fq;
			m18A -= fq;
            //---------------------------------------------------------------------//

            // ------------------- Fluid component B ---------------------------------//
            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
			// q=0
			fq = distB[n];
            rhoB_next = fq;
			m1B  = -30.0*fq;
			m2B  = 12.0*fq;

			// q=1
            fq = distB[2*Np+n];
            rhoB_next += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jxB = fq;
			m4B = -4.0*fq;
			m9B = 2.0*fq;
			m10B = -4.0*fq;

			// q=2
            fq = distB[1*Np+n];
            rhoB_next += fq;
			m1B -= 11.0*(fq);
			m2B -= 4.0*(fq);
			jxB -= fq;
			m4B += 4.0*(fq);
			m9B += 2.0*(fq);
			m10B -= 4.0*(fq);

			// q=3
            fq = distB[4*Np+n];
            rhoB_next += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jyB = fq;
			m6B = -4.0*fq;
			m9B -= fq;
			m10B += 2.0*fq;
			m11B = fq;
			m12B = -2.0*fq;

			// q = 4
            fq = distB[3*Np+n];
            rhoB_next += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jyB -= fq;
			m6B += 4.0*fq;
			m9B -= fq;
			m10B += 2.0*fq;
			m11B += fq;
			m12B -= 2.0*fq;

			// q=5
            fq = distB[6*Np+n];
            rhoB_next += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jzB = fq;
			m8B = -4.0*fq;
			m9B -= fq;
			m10B += 2.0*fq;
			m11B -= fq;
			m12B += 2.0*fq;


			// q = 6
            fq = distB[5*Np+n];
            rhoB_next += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jzB -= fq;
			m8B += 4.0*fq;
			m9B -= fq;
			m10B += 2.0*fq;
			m11B -= fq;
			m12B += 2.0*fq;

			// q=7
            fq = distB[8*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB += fq;
			m4B += fq;
			jyB += fq;
			m6B += fq;
			m9B  += fq;
			m10B += fq;
			m11B += fq;
			m12B += fq;
			m13B = fq;
			m16B = fq;
			m17B = -fq;

			// q = 8
            fq = distB[7*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB -= fq;
			m4B -= fq;
			jyB -= fq;
			m6B -= fq;
			m9B += fq;
			m10B += fq;
			m11B += fq;
			m12B += fq;
			m13B += fq;
			m16B -= fq;
			m17B += fq;

			// q=9
            fq = distB[10*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB += fq;
			m4B += fq;
			jyB -= fq;
			m6B -= fq;
			m9B += fq;
			m10B += fq;
			m11B += fq;
			m12B += fq;
			m13B -= fq;
			m16B += fq;
			m17B += fq;

			// q = 10
            fq = distB[9*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB -= fq;
			m4B -= fq;
			jyB += fq;
			m6B += fq;
			m9B += fq;
			m10B += fq;
			m11B += fq;
			m12B += fq;
			m13B -= fq;
			m16B -= fq;
			m17B -= fq;

			// q=11
            fq = distB[12*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB += fq;
			m4B += fq;
			jzB += fq;
			m8B += fq;
			m9B += fq;
			m10B += fq;
			m11B -= fq;
			m12B -= fq;
			m15B = fq;
			m16B -= fq;
			m18B = fq;

			// q=12
            fq = distB[11*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB -= fq;
			m4B -= fq;
			jzB -= fq;
			m8B -= fq;
			m9B += fq;
			m10B += fq;
			m11B -= fq;
			m12B -= fq;
			m15B += fq;
			m16B += fq;
			m18B -= fq;

			// q=13
            fq = distB[14*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB += fq;
			m4B += fq;
			jzB -= fq;
			m8B -= fq;
			m9B += fq;
			m10B += fq;
			m11B -= fq;
			m12B -= fq;
			m15B -= fq;
			m16B -= fq;
			m18B -= fq;

			// q=14
            fq = distB[13*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jxB -= fq;
			m4B -= fq;
			jzB += fq;
			m8B += fq;
			m9B += fq;
			m10B += fq;
			m11B -= fq;
			m12B -= fq;
			m15B -= fq;
			m16B += fq;
			m18B += fq;

			// q=15
            fq = distB[16*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jyB += fq;
			m6B += fq;
			jzB += fq;
			m8B += fq;
			m9B -= 2.0*fq;
			m10B -= 2.0*fq;
			m14B = fq;
			m17B += fq;
			m18B -= fq;

			// q=16
            fq = distB[15*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jyB -= fq;
			m6B -= fq;
			jzB -= fq;
			m8B -= fq;
			m9B -= 2.0*fq;
			m10B -= 2.0*fq;
			m14B += fq;
			m17B -= fq;
			m18B += fq;

			// q=17
            fq = distB[18*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jyB += fq;
			m6B += fq;
			jzB -= fq;
			m8B -= fq;
			m9B -= 2.0*fq;
			m10B -= 2.0*fq;
			m14B -= fq;
			m17B += fq;
			m18B += fq;

			// q=18
            fq = distB[17*Np+n];
            rhoB_next += fq;
			m1B += 8.0*fq;
			m2B += fq;
			jyB -= fq;
			m6B -= fq;
			jzB += fq;
			m8B += fq;
			m9B -= 2.0*fq;
			m10B -= 2.0*fq;
			m14B -= fq;
			m17B -= fq;
			m18B -= fq;
            //---------------------------------------------------------------------//

            
            // Compute SC fluid-fluid interaction force
            GffA_x = -Gsc*rhoB_gradx;
            GffA_y = -Gsc*rhoB_grady;
            GffA_z = -Gsc*rhoB_gradz;
            GffB_x = -Gsc*rhoA_gradx;
            GffB_y = -Gsc*rhoA_grady;
            GffB_z = -Gsc*rhoA_gradz;
            // Compute SC fluid-solid force
            GfsA_x = SolidForceA[n+0*Np];    
            GfsA_y = SolidForceA[n+1*Np];    
            GfsA_z = SolidForceA[n+2*Np];    
            GfsB_x = SolidForceB[n+0*Np];    
            GfsB_y = SolidForceB[n+1*Np];    
            GfsB_z = SolidForceB[n+2*Np];    

            // Compute greyscale related parameters
            // ------------------- Fluid Component A -----------------------//
            c0 = 0.5*(1.0+porosity*0.5*muA_eff/permA);
            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
            c1 = porosity*0.5*GeoFun/sqrt(permA);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jxA/rhoA+0.5*(porosity*Gx+GffA_x+GfsA_x);
            vy = jyA/rhoA+0.5*(porosity*Gy+GffA_y+GfsA_y);
            vz = jzA/rhoA+0.5*(porosity*Gz+GffA_z+GfsA_z);
            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
            ux_A = vx/(c0+sqrt(c0*c0+c1*v_mag));
            uy_A = vy/(c0+sqrt(c0*c0+c1*v_mag));
            uz_A = vz/(c0+sqrt(c0*c0+c1*v_mag));
            u_mag=sqrt(ux_A*ux_A+uy_A*uy_A+uz_A*uz_A);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            FxA = rhoA*(-porosity*muA_eff/permA*ux_A - porosity*GeoFun/sqrt(permA)*u_mag*ux_A + porosity*Gx + GffA_x + GfsA_x);
            FyA = rhoA*(-porosity*muA_eff/permA*uy_A - porosity*GeoFun/sqrt(permA)*u_mag*uy_A + porosity*Gy + GffA_y + GfsA_y);
            FzA = rhoA*(-porosity*muA_eff/permA*uz_A - porosity*GeoFun/sqrt(permA)*u_mag*uz_A + porosity*Gz + GffA_z + GfsA_z);
            if (porosity==1.0){
                FxA=rhoA*(Gx + GffA_x + GfsA_x);
                FyA=rhoA*(Gy + GffA_y + GfsA_y);
                FzA=rhoA*(Gz + GffA_z + GfsA_z);
            }
            // ------------------- Fluid Component B -----------------------//
            // Compute greyscale related parameters
            c0 = 0.5*(1.0+porosity*0.5*muB_eff/permB);
            if (porosity==1.0) c0 = 0.5;//i.e. apparent pore nodes
            //GeoFun = 1.75/sqrt(150.0*porosity*porosity*porosity);
            c1 = porosity*0.5*GeoFun/sqrt(permB);
            if (porosity==1.0) c1 = 0.0;//i.e. apparent pore nodes

            vx = jxB/rhoB+0.5*(porosity*Gx+GffB_x+GfsB_x);
            vy = jyB/rhoB+0.5*(porosity*Gy+GffB_y+GfsB_y);
            vz = jzB/rhoB+0.5*(porosity*Gz+GffB_z+GfsB_z);
            v_mag=sqrt(vx*vx+vy*vy+vz*vz);
            ux_B = vx/(c0+sqrt(c0*c0+c1*v_mag));
            uy_B = vy/(c0+sqrt(c0*c0+c1*v_mag));
            uz_B = vz/(c0+sqrt(c0*c0+c1*v_mag));
            u_mag=sqrt(ux_B*ux_B+uy_B*uy_B+uz_B*uz_B);

            //Update the total force to include linear (Darcy) and nonlinear (Forchheimer) drags due to the porous medium
            FxB = rhoB*(-porosity*muB_eff/permB*ux_B - porosity*GeoFun/sqrt(permB)*u_mag*ux_B + porosity*Gx + GffB_x + GfsB_x);
            FyB = rhoB*(-porosity*muB_eff/permB*uy_B - porosity*GeoFun/sqrt(permB)*u_mag*uy_B + porosity*Gy + GffB_y + GfsB_y);
            FzB = rhoB*(-porosity*muB_eff/permB*uz_B - porosity*GeoFun/sqrt(permB)*u_mag*uz_B + porosity*Gz + GffB_z + GfsB_z);
            if (porosity==1.0){
                FxB=rhoB*(Gx + GffB_x + GfsB_x);
                FyB=rhoB*(Gy + GffB_y + GfsB_y);
                FzB=rhoB*(Gz + GffB_z + GfsB_z);
            }

            // Calculate barycentric velocity of the fluid mixture
            ux = (rhoA*ux_A+rhoB*ux_B)/(rhoA+rhoB);
            uy = (rhoA*uy_A+rhoB*uy_B)/(rhoA+rhoB);
            uz = (rhoA*uz_A+rhoB*uz_B)/(rhoA+rhoB);


//            //..............carry out relaxation process...............................................
//            m1 = m1 + rlx_setA*((-30*Den+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
//                    + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m2 = m2 + rlx_setA*((12*Den - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
//                    + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
//            jx = jx + Fx;
//            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*Den) - m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*Den) - m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*Den) - m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//            m9 = m9 + rlx_setA*((Den*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
//                    + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
//            m10 = m10 + rlx_setA*(-0.5*Den*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
//                      + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m11 = m11 + rlx_setA*((Den*(uy*uy-uz*uz)/porosity) - m11)
//                      + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
//            m12 = m12 + rlx_setA*(-0.5*(Den*(uy*uy-uz*uz)/porosity)- m12)
//                      + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
//            m13 = m13 + rlx_setA*((Den*ux*uy/porosity) - m13)
//                      + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
//            m14 = m14 + rlx_setA*((Den*uy*uz/porosity) - m14)
//                      + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
//            m15 = m15 + rlx_setA*((Den*ux*uz/porosity) - m15)
//                      + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
//            m16 = m16 + rlx_setB*( - m16);
//            m17 = m17 + rlx_setB*( - m17);
//            m18 = m18 + rlx_setB*( - m18);
//            //.......................................................................................................
           

            // ------------------- Fluid Component A -----------------------//
            rlx_setA = 1.0/tauA;
            rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
            //-------------------- MRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            //TODO need to incoporate porosity
			m1A = m1A + rlx_setA*((19*rhoA*(ux*ux+uy*uy+uz*uz) - 11*rhoA_next) - m1A)
                      + (1-0.5*rlx_setA)*38*(FxA*ux+FyA*uy+FzA*uz);
			m2A = m2A + rlx_setA*((3*rhoA_next - 5.5*rhoA*(ux*ux+uy*uy+uz*uz))- m2A)
                      + (1-0.5*rlx_setA)*11*(-FxA*ux-FyA*uy-FzA*uz);
            jxA = jxA + FxA;
			m4A = m4A + rlx_setB*((-0.6666666666666666*ux*rhoA)- m4A)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FxA);
            jyA = jyA + FyA;
			m6A = m6A + rlx_setB*((-0.6666666666666666*uy*rhoA)- m6A)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FyA);
            jzA = jzA + FzA;
			m8A = m8A + rlx_setB*((-0.6666666666666666*uz*rhoA)- m8A)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FzA);
			m9A = m9A + rlx_setA*((rhoA*(2*ux*ux-uy*uy-uz*uz)) - m9A)
                      + (1-0.5*rlx_setA)*(4*FxA*ux-2*FyA*uy-2*FzA*uz);
			m10A = m10A + rlx_setA*( - m10A)
                        + (1-0.5*rlx_setA)*(-2*FxA*ux+FyA*uy+FzA*uz);
			m11A = m11A + rlx_setA*((rhoA*(uy*uy-uz*uz)) - m11A)
                        + (1-0.5*rlx_setA)*(2*FyA*uy-2*FzA*uz);
			m12A = m12A + rlx_setA*( - m12A)
                        + (1-0.5*rlx_setA)*(-FyA*uy+FzA*uz);
			m13A = m13A + rlx_setA*( rhoA*(ux*uy) - m13A)
                        + (1-0.5*rlx_setA)*(FyA*ux+FxA*uy);
			m14A = m14A + rlx_setA*( rhoA*(uy*uz) - m14A)
                        + (1-0.5*rlx_setA)*(FzA*uy+FyA*uz);
			m15A = m15A + rlx_setA*( rhoA*(ux*uz) - m15A)
                        + (1-0.5*rlx_setA)*(FzA*ux+FxA*uz);
			m16A = m16A + rlx_setB*( - m16A);
			m17A = m17A + rlx_setB*( - m17A);
			m18A = m18A + rlx_setB*( - m18A);
            //.......................................................................................................


            // ------------------- Fluid Component A -----------------------//
            //.................inverse transformation......................................................
            // q=0
            fq = mrt_V1*rhoA_next-mrt_V2*m1A+mrt_V3*m2A;
            //f0 = mrt_V1*rhoA_next-mrt_V2*m1A+mrt_V3*m2A;
            distA[n] = fq;

            // q = 1
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jxA-m4A)+mrt_V6*(m9A-m10A);
            //f1 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jxA-m4A)+mrt_V6*(m9A-m10A);
            distA[1*Np+n] = fq;

            // q=2
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m4A-jxA)+mrt_V6*(m9A-m10A);
            //f2 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m4A-jxA)+mrt_V6*(m9A-m10A);
            distA[2*Np+n] = fq;

            // q = 3
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jyA-m6A)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            //f3 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jyA-m6A)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            distA[3*Np+n] = fq;

            // q = 4
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m6A-jyA)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            //f4 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m6A-jyA)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            distA[4*Np+n] = fq;

            // q = 5
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jzA-m8A)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            //f5 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(jzA-m8A)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            distA[5*Np+n] = fq;

            // q = 6
            fq = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m8A-jzA)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            //f6 = mrt_V1*rhoA_next-mrt_V4*m1A-mrt_V5*m2A+0.1*(m8A-jzA)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            distA[6*Np+n] = fq;

            // q = 7
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jyA)+0.025*(m4A+m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m16A-m17A);
            //f7 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jyA)+0.025*(m4A+m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m16A-m17A);
            distA[7*Np+n] = fq;

            // q = 8
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jyA)-0.025*(m4A+m6A) +mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m17A-m16A);
            //f8 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jyA)-0.025*(m4A+m6A) +mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m17A-m16A);
            distA[8*Np+n] = fq;

            // q = 9
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jyA)+0.025*(m4A-m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A+0.125*(m16A+m17A);
            //f9 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jyA)+0.025*(m4A-m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A+0.125*(m16A+m17A);
            distA[9*Np+n] = fq;

            // q = 10
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jxA)+0.025*(m6A-m4A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A-0.125*(m16A+m17A);
            //f10 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jxA)+0.025*(m6A-m4A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A-0.125*(m16A+m17A);
            distA[10*Np+n] = fq;

            // q = 11
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jzA)+0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m18A-m16A);
            //f11 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jzA)+0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m18A-m16A);
            distA[11*Np+n] = fq;

            // q = 12
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jzA)-0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m16A-m18A);
            //f12 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jzA)-0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m16A-m18A);
            distA[12*Np+n] = fq;

            // q = 13
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jzA)+0.025*(m4A-m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A-0.125*(m16A+m18A);
            //f13 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jzA)+0.025*(m4A-m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A-0.125*(m16A+m18A);
            distA[13*Np+n] = fq;

            // q= 14
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jxA)+0.025*(m8A-m4A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A+0.125*(m16A+m18A);
            //f14 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jxA)+0.025*(m8A-m4A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A+0.125*(m16A+m18A);
            distA[14*Np+n] = fq;

            // q = 15
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA+jzA)+0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m17A-m18A);
            //f15 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA+jzA)+0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m17A-m18A);
            distA[15*Np+n] = fq;

            // q = 16
            fq =  mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jyA+jzA)-0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m18A-m17A);
            //f16 =  mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A-0.1*(jyA+jzA)-0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m18A-m17A);
            distA[16*Np+n] = fq;

            // q = 17
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jzA)+0.025*(m6A-m8A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A+0.125*(m17A+m18A);
            //f17 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jzA)+0.025*(m6A-m8A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A+0.125*(m17A+m18A);
            distA[17*Np+n] = fq;

            // q = 18
            fq = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jyA)+0.025*(m8A-m6A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A-0.125*(m17A+m18A);
            //f18 = mrt_V1*rhoA_next+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jyA)+0.025*(m8A-m6A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A-0.125*(m17A+m18A);
            distA[18*Np+n] = fq;
            //........................................................................

            //Den[n] = f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;

//            //..............carry out relaxation process...............................................
//            m1 = m1 + rlx_setA*((-30*Den+19*(ux*ux+uy*uy+uz*uz)/porosity + 57*pressure*porosity) - m1) 
//                    + (1-0.5*rlx_setA)*38*(Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m2 = m2 + rlx_setA*((12*Den - 5.5*(ux*ux+uy*uy+uz*uz)/porosity-27*pressure*porosity) - m2)
//                    + (1-0.5*rlx_setA)*11*(-Fx*ux-Fy*uy-Fz*uz)/porosity;
//            jx = jx + Fx;
//            m4 = m4 + rlx_setB*((-0.6666666666666666*ux*Den) - m4)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fx);
//            jy = jy + Fy;
//            m6 = m6 + rlx_setB*((-0.6666666666666666*uy*Den) - m6)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fy);
//            jz = jz + Fz;
//            m8 = m8 + rlx_setB*((-0.6666666666666666*uz*Den) - m8)
//                    + (1-0.5*rlx_setB)*(-0.6666666666666666*Fz);
//            m9 = m9 + rlx_setA*((Den*(2*ux*ux-uy*uy-uz*uz)/porosity) - m9)
//                    + (1-0.5*rlx_setA)*(4*Fx*ux-2*Fy*uy-2*Fz*uz)/porosity;
//            m10 = m10 + rlx_setA*(-0.5*Den*((2*ux*ux-uy*uy-uz*uz)/porosity)- m10)
//                      + (1-0.5*rlx_setA)*(-2*Fx*ux+Fy*uy+Fz*uz)/porosity;
//            m11 = m11 + rlx_setA*((Den*(uy*uy-uz*uz)/porosity) - m11)
//                      + (1-0.5*rlx_setA)*(2*Fy*uy-2*Fz*uz)/porosity;
//            m12 = m12 + rlx_setA*(-0.5*(Den*(uy*uy-uz*uz)/porosity)- m12)
//                      + (1-0.5*rlx_setA)*(-Fy*uy+Fz*uz)/porosity;
//            m13 = m13 + rlx_setA*((Den*ux*uy/porosity) - m13)
//                      + (1-0.5*rlx_setA)*(Fy*ux+Fx*uy)/porosity;
//            m14 = m14 + rlx_setA*((Den*uy*uz/porosity) - m14)
//                      + (1-0.5*rlx_setA)*(Fz*uy+Fy*uz)/porosity;
//            m15 = m15 + rlx_setA*((Den*ux*uz/porosity) - m15)
//                      + (1-0.5*rlx_setA)*(Fz*ux+Fx*uz)/porosity;
//            m16 = m16 + rlx_setB*( - m16);
//            m17 = m17 + rlx_setB*( - m17);
//            m18 = m18 + rlx_setB*( - m18);
//            //.......................................................................................................
           

            // ------------------- Fluid Component B -----------------------//
            rlx_setA = 1.0/tauB;
            rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
            //-------------------- MRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            //TODO need to incoporate porosity
			m1B = m1B + rlx_setA*((19*rhoB*(ux*ux+uy*uy+uz*uz) - 11*rhoB_next) - m1B)
                      + (1-0.5*rlx_setA)*38*(FxB*ux+FyB*uy+FzB*uz);
			m2B = m2B + rlx_setA*((3*rhoB_next - 5.5*rhoB*(ux*ux+uy*uy+uz*uz))- m2B)
                      + (1-0.5*rlx_setA)*11*(-FxB*ux-FyB*uy-FzB*uz);
            jxB = jxB + FxB;
			m4B = m4B + rlx_setB*((-0.6666666666666666*ux*rhoB)- m4B)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FxB);
            jyB = jyB + FyB;
			m6B = m6B + rlx_setB*((-0.6666666666666666*uy*rhoB)- m6B)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FyB);
            jzB = jzB + FzB;
			m8B = m8B + rlx_setB*((-0.6666666666666666*uz*rhoB)- m8B)
                    + (1-0.5*rlx_setB)*(-0.6666666666666666*FzB);
			m9B = m9B + rlx_setA*((rhoB*(2*ux*ux-uy*uy-uz*uz)) - m9B)
                      + (1-0.5*rlx_setA)*(4*FxB*ux-2*FyB*uy-2*FzB*uz);
			m10B = m10B + rlx_setA*( - m10B)
                        + (1-0.5*rlx_setA)*(-2*FxB*ux+FyB*uy+FzB*uz);
			m11B = m11B + rlx_setA*((rhoB*(uy*uy-uz*uz)) - m11B)
                        + (1-0.5*rlx_setA)*(2*FyB*uy-2*FzB*uz);
			m12B = m12B + rlx_setA*( - m12B)
                        + (1-0.5*rlx_setA)*(-FyB*uy+FzB*uz);
			m13B = m13B + rlx_setA*( rhoB*(ux*uy) - m13B)
                        + (1-0.5*rlx_setA)*(FyB*ux+FxB*uy);
			m14B = m14B + rlx_setA*( rhoB*(uy*uz) - m14B)
                        + (1-0.5*rlx_setA)*(FzB*uy+FyB*uz);
			m15B = m15B + rlx_setA*( rhoB*(ux*uz) - m15B)
                        + (1-0.5*rlx_setA)*(FzB*ux+FxB*uz);
			m16B = m16B + rlx_setB*( - m16B);
			m17B = m17B + rlx_setB*( - m17B);
			m18B = m18B + rlx_setB*( - m18B);
            //.......................................................................................................


            // ------------------- Fluid Component B -----------------------//
            //.................inverse transformation......................................................
            // q=0
            fq = mrt_V1*rhoB_next-mrt_V2*m1B+mrt_V3*m2B;
            //f0 = mrt_V1*rhoB_next-mrt_V2*m1B+mrt_V3*m2B;
            distB[n] = fq;

            // q = 1
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jxB-m4B)+mrt_V6*(m9B-m10B);
            //f1 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jxB-m4B)+mrt_V6*(m9B-m10B);
            distB[1*Np+n] = fq;

            // q=2
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m4B-jxB)+mrt_V6*(m9B-m10B);
            //f2 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m4B-jxB)+mrt_V6*(m9B-m10B);
            distB[2*Np+n] = fq;

            // q = 3
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jyB-m6B)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            //f3 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jyB-m6B)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            distB[3*Np+n] = fq;

            // q = 4
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m6B-jyB)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            //f4 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m6B-jyB)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            distB[4*Np+n] = fq;

            // q = 5
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jzB-m8B)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            //f5 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(jzB-m8B)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            distB[5*Np+n] = fq;

            // q = 6
            fq = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m8B-jzB)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            //f6 = mrt_V1*rhoB_next-mrt_V4*m1B-mrt_V5*m2B+0.1*(m8B-jzB)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            distB[6*Np+n] = fq;

            // q = 7
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jyB)+0.025*(m4B+m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m16B-m17B);
            //f7 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jyB)+0.025*(m4B+m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m16B-m17B);
            distB[7*Np+n] = fq;

            // q = 8
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jyB)-0.025*(m4B+m6B) +mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m17B-m16B);
            //f8 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jyB)-0.025*(m4B+m6B) +mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m17B-m16B);
            distB[8*Np+n] = fq;

            // q = 9
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jyB)+0.025*(m4B-m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B+0.125*(m16B+m17B);
            //f9 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jyB)+0.025*(m4B-m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B+0.125*(m16B+m17B);
            distB[9*Np+n] = fq;

            // q = 10
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jxB)+0.025*(m6B-m4B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B-0.125*(m16B+m17B);
            //f10 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jxB)+0.025*(m6B-m4B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B-0.125*(m16B+m17B);
            distB[10*Np+n] = fq;

            // q = 11
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jzB)+0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m18B-m16B);
            //f11 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jzB)+0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m18B-m16B);
            distB[11*Np+n] = fq;

            // q = 12
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jzB)-0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m16B-m18B);
            //f12 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jzB)-0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m16B-m18B);
            distB[12*Np+n] = fq;

            // q = 13
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jzB)+0.025*(m4B-m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B-0.125*(m16B+m18B);
            //f13 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jzB)+0.025*(m4B-m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B-0.125*(m16B+m18B);
            distB[13*Np+n] = fq;

            // q= 14
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jxB)+0.025*(m8B-m4B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B+0.125*(m16B+m18B);
            //f14 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jxB)+0.025*(m8B-m4B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B+0.125*(m16B+m18B);
            distB[14*Np+n] = fq;

            // q = 15
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB+jzB)+0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m17B-m18B);
            //f15 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB+jzB)+0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m17B-m18B);
            distB[15*Np+n] = fq;

            // q = 16
            fq =  mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jyB+jzB)-0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m18B-m17B);
            //f16 =  mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B-0.1*(jyB+jzB)-0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m18B-m17B);
            distB[16*Np+n] = fq;

            // q = 17
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jzB)+0.025*(m6B-m8B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B+0.125*(m17B+m18B);
            //f17 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jzB)+0.025*(m6B-m8B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B+0.125*(m17B+m18B);
            distB[17*Np+n] = fq;

            // q = 18
            fq = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jyB)+0.025*(m8B-m6B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B-0.125*(m17B+m18B);
            //f18 = mrt_V1*rhoB_next+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jyB)+0.025*(m8B-m6B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B-0.125*(m17B+m18B);
            distB[18*Np+n] = fq;
            //........................................................................

            //Den[n+Np] = f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
            //Update velocity on device
            Velocity[0*Np+n] = ux;
            Velocity[1*Np+n] = uy;
            Velocity[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = (rhoA+rhoB+Gsc*rhoA*rhoB)/3.0;
            //Update density
            //Den[n] = rhoA_next;
            //Den[n+Np] = rhoB_next;

		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_GreyscaleSC_Init(double *distA, double *distB, double *Den, int Np)
{
	int n;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<Np ){
			distA[0*Np+n]  = Den[0*Np+n]*0.3333333333333333;
			distA[1*Np+n]  = Den[0*Np+n]*0.055555555555555555;		//double(100*n)+1.f;
			distA[2*Np+n]  = Den[0*Np+n]*0.055555555555555555;	//double(100*n)+2.f;
			distA[3*Np+n]  = Den[0*Np+n]*0.055555555555555555;	//double(100*n)+3.f;
			distA[4*Np+n]  = Den[0*Np+n]*0.055555555555555555;	//double(100*n)+4.f;
			distA[5*Np+n]  = Den[0*Np+n]*0.055555555555555555;	//double(100*n)+5.f;
			distA[6*Np+n]  = Den[0*Np+n]*0.055555555555555555;	//double(100*n)+6.f;
			distA[7*Np+n]  = Den[0*Np+n]*0.0277777777777778;   //double(100*n)+7.f;
			distA[8*Np+n]  = Den[0*Np+n]*0.0277777777777778;   //double(100*n)+8.f;
			distA[9*Np+n]  = Den[0*Np+n]*0.0277777777777778;   //double(100*n)+9.f;
			distA[10*Np+n] = Den[0*Np+n]*0.0277777777777778;  //double(100*n)+10.f;
			distA[11*Np+n] = Den[0*Np+n]*0.0277777777777778;  //double(100*n)+11.f;
			distA[12*Np+n] = Den[0*Np+n]*0.0277777777777778;  //double(100*n)+12.f;
			distA[13*Np+n] = Den[0*Np+n]*0.0277777777777778;  //double(100*n)+13.f;
			distA[14*Np+n] = Den[0*Np+n]*0.0277777777777778;  //double(100*n)+14.f;
			distA[15*Np+n] = Den[0*Np+n]*0.0277777777777778;  //double(100*n)+15.f;
			distA[16*Np+n] = Den[0*Np+n]*0.0277777777777778;  //double(100*n)+16.f;
			distA[17*Np+n] = Den[0*Np+n]*0.0277777777777778;  //double(100*n)+17.f;
			distA[18*Np+n] = Den[0*Np+n]*0.0277777777777778;  //double(100*n)+18.f;

			distB[0*Np+n]  = Den[1*Np+n]*0.3333333333333333;
			distB[1*Np+n]  = Den[1*Np+n]*0.055555555555555555;		//double(100*n)+1.f;
			distB[2*Np+n]  = Den[1*Np+n]*0.055555555555555555;	//double(100*n)+2.f;
			distB[3*Np+n]  = Den[1*Np+n]*0.055555555555555555;	//double(100*n)+3.f;
			distB[4*Np+n]  = Den[1*Np+n]*0.055555555555555555;	//double(100*n)+4.f;
			distB[5*Np+n]  = Den[1*Np+n]*0.055555555555555555;	//double(100*n)+5.f;
			distB[6*Np+n]  = Den[1*Np+n]*0.055555555555555555;	//double(100*n)+6.f;
			distB[7*Np+n]  = Den[1*Np+n]*0.0277777777777778;   //double(100*n)+7.f;
			distB[8*Np+n]  = Den[1*Np+n]*0.0277777777777778;   //double(100*n)+8.f;
			distB[9*Np+n]  = Den[1*Np+n]*0.0277777777777778;   //double(100*n)+9.f;
			distB[10*Np+n] = Den[1*Np+n]*0.0277777777777778;  //double(100*n)+10.f;
			distB[11*Np+n] = Den[1*Np+n]*0.0277777777777778;  //double(100*n)+11.f;
			distB[12*Np+n] = Den[1*Np+n]*0.0277777777777778;  //double(100*n)+12.f;
			distB[13*Np+n] = Den[1*Np+n]*0.0277777777777778;  //double(100*n)+13.f;
			distB[14*Np+n] = Den[1*Np+n]*0.0277777777777778;  //double(100*n)+14.f;
			distB[15*Np+n] = Den[1*Np+n]*0.0277777777777778;  //double(100*n)+15.f;
			distB[16*Np+n] = Den[1*Np+n]*0.0277777777777778;  //double(100*n)+16.f;
			distB[17*Np+n] = Den[1*Np+n]*0.0277777777777778;  //double(100*n)+17.f;
			distB[18*Np+n] = Den[1*Np+n]*0.0277777777777778;  //double(100*n)+18.f;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAodd_GreyscaleSC_Density(int *neighborList, double *distA, double *distB, double *Den, int start, int finish, int Np){
	int n,nread;
	double fq,nA,nB;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			//..........Compute the number density for each component ............
			// q=0
			fq = distA[n];
			nA = fq;
			fq = distB[n];
			nB = fq;
			
			// q=1
			nread = neighborList[n]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;
			
			// q=2
			nread = neighborList[n+Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=3
			nread = neighborList[n+2*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=4
			nread = neighborList[n+3*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=5
			nread = neighborList[n+4*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=6
			nread = neighborList[n+5*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=7
			nread = neighborList[n+6*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=8
			nread = neighborList[n+7*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=9
			nread = neighborList[n+8*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=10
			nread = neighborList[n+9*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=11
			nread = neighborList[n+10*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=12
			nread = neighborList[n+11*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=13
			nread = neighborList[n+12*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=14
			nread = neighborList[n+13*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=15
			nread = neighborList[n+14*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=16
			nread = neighborList[n+15*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=17
			nread = neighborList[n+16*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// q=18
			nread = neighborList[n+17*Np]; 
			fq = distA[nread];
			nA += fq;
			fq = distB[nread]; 
			nB += fq;

			// save the number densities
			Den[n] = nA;
			Den[Np+n] = nB;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAeven_GreyscaleSC_Density(double *distA, double *distB, double *Den, int start, int finish, int Np){
	int n,nread;
	double fq,nA,nB;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {
			//..........Compute the number density for each component ............
			// q=0
			fq = distA[n];
			nA = fq;
			fq = distB[n];
			nB = fq;
			
            // q=1
            fq = distA[2*Np+n];
			nA += fq;
            fq = distB[2*Np+n];
			nB += fq;

            // q=2
            fq = distA[1*Np+n];
			nA += fq;
            fq = distB[1*Np+n];
			nB += fq;

            // q=3
            fq = distA[4*Np+n];
			nA += fq;
            fq = distB[4*Np+n];
			nB += fq;

            // q = 4
            fq = distA[3*Np+n];
			nA += fq;
            fq = distB[3*Np+n];
			nB += fq;

            // q=5
            fq = distA[6*Np+n];
			nA += fq;
            fq = distB[6*Np+n];
			nB += fq;

            // q = 6
            fq = distA[5*Np+n];
			nA += fq;
            fq = distB[5*Np+n];
			nB += fq;

            // q=7
            fq = distA[8*Np+n];
			nA += fq;
            fq = distB[8*Np+n];
			nB += fq;

            // q = 8
            fq = distA[7*Np+n];
			nA += fq;
            fq = distB[7*Np+n];
			nB += fq;

            // q=9
            fq = distA[10*Np+n];
			nA += fq;
            fq = distB[10*Np+n];
			nB += fq;

            // q = 10
            fq = distA[9*Np+n];
			nA += fq;
            fq = distB[9*Np+n];
			nB += fq;

            // q=11
            fq = distA[12*Np+n];
			nA += fq;
            fq = distB[12*Np+n];
			nB += fq;

            // q=12
            fq = distA[11*Np+n];
			nA += fq;
            fq = distB[11*Np+n];
			nB += fq;

            // q=13
            fq = distA[14*Np+n];
			nA += fq;
            fq = distB[14*Np+n];
			nB += fq;

            // q=14
            fq = distA[13*Np+n];
			nA += fq;
            fq = distB[13*Np+n];
			nB += fq;

            // q=15
            fq = distA[16*Np+n];
			nA += fq;
            fq = distB[16*Np+n];
			nB += fq;

            // q=16
            fq = distA[15*Np+n];
			nA += fq;
            fq = distB[15*Np+n];
			nB += fq;

            // q=17
            fq = distA[18*Np+n];
			nA += fq;
            fq = distB[18*Np+n];
			nB += fq;

            // q=18
            fq = distA[17*Np+n];
			nA += fq;
            fq = distB[17*Np+n];
			nB += fq;

			// save the number densities
			Den[n] = nA;
			Den[Np+n] = nB;
		}
	}
}

extern "C" void ScaLBL_D3Q19_GreyscaleSC_Init(double *distA,double *distB, double *Den, int Np){
	dvc_ScaLBL_D3Q19_GreyscaleSC_Init<<<NBLOCKS,NTHREADS >>>(distA,distB,Den,Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_GreyscaleSC_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleSC_Density(int *NeighborList, double *distA, double *distB, double *Den, int start, int finish, int Np){

	dvc_ScaLBL_D3Q19_AAodd_GreyscaleSC_Density<<<NBLOCKS,NTHREADS >>>(NeighborList, distA, distB, Den, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_GreyscaleSC_Density: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleSC_Density(double *distA, double *distB, double *Den, int start, int finish, int Np){

	dvc_ScaLBL_D3Q19_AAeven_GreyscaleSC_Density<<<NBLOCKS,NTHREADS >>>(distA, distB, Den, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_GreyscaleSC_Density: %s \n",cudaGetErrorString(err));
	}
}


extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleSC(int *neighborList, double *distA, double *distB, double *Den, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

    dvc_ScaLBL_D3Q19_AAodd_GreyscaleSC<<<NBLOCKS,NTHREADS >>>(neighborList,distA,distB,Den,DenGradA,DenGradB,SolidForceA,SolidForceB,Poros,Perm,Velocity,Pressure, 
                tauA,tauB,tauA_eff,tauB_eff,Gsc,Gx,Gy,Gz,start,finish,Np);

    cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_GreyscaleSC: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleSC(double *distA, double *distB, double *Den, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

    dvc_ScaLBL_D3Q19_AAeven_GreyscaleSC<<<NBLOCKS,NTHREADS >>>(distA,distB,Den,DenGradA,DenGradB,SolidForceA,SolidForceB,Poros,Perm,Velocity,Pressure, 
                tauA,tauB,tauA_eff,tauB_eff,Gsc,Gx,Gy,Gz,start,finish,Np);

    cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_GreyscaleSC: %s \n",cudaGetErrorString(err));
	}
}
