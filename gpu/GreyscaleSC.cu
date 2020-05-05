#include <stdio.h>

#define NBLOCKS 1024
#define NTHREADS 256

__global__ void dvc_ScaLBL_D3Q19_AAodd_GreyscaleSC_MRT(int *neighborList,int *Map, double *distA, double *distB, double *DenA,double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

    int ijk;
	int n, nread;
	double vx,vy,vz,v_mag;
    double ux_A,uy_A,uz_A,ux_B,uy_B,uz_B,u_mag;
    double ux,uy,uz;
	// conserved momemnts
	double jxA,jyA,jzA;
	double jxB,jyB,jzB;
    double rhoA,rhoB;
    double nA,nB;
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
    double nA_gradx,nA_grady,nA_gradz;
    double nB_gradx,nB_grady,nB_gradz;
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
            ijk = Map[n];
            nA = DenA[ijk];
            nB = DenB[ijk];
            porosity = Poros[n];
            perm = Perm[n];
            permA = perm*nA/(nA+nB);//effective relative perm
            permB = perm*nB/(nA+nB);
			nA_gradx = DenGradA[n+0*Np];
			nA_grady = DenGradA[n+1*Np];
			nA_gradz = DenGradA[n+2*Np];
			nB_gradx = DenGradB[n+0*Np];
			nB_grady = DenGradB[n+1*Np];
			nB_gradz = DenGradB[n+2*Np];

            // ------------------- Fluid component A ---------------------------------//
            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
			// q=0
			fq = distA[n];
            rhoA = fq;
			m1A  = -30.0*fq;
			m2A  = 12.0*fq;

			// q=1
			nread = neighborList[n]; // neighbor 2 
			fq = distA[nread]; // reading the f1 data into register fq		
            rhoA += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jxA = fq;
			m4A = -4.0*fq;
			m9A = 2.0*fq;
			m10A = -4.0*fq;

			// q=2
			nread = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
			fq = distA[nread];  // reading the f2 data into register fq
            rhoA += fq;
			m1A -= 11.0*(fq);
			m2A -= 4.0*(fq);
			jxA -= fq;
			m4A += 4.0*(fq);
			m9A += 2.0*(fq);
			m10A -= 4.0*(fq);

			// q=3
			nread = neighborList[n+2*Np]; // neighbor 4
			fq = distA[nread];
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoB = fq;
			m1B  = -30.0*fq;
			m2B  = 12.0*fq;

			// q=1
			nread = neighborList[n]; // neighbor 2 
			fq = distB[nread]; // reading the f1 data into register fq		
            rhoB += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jxB = fq;
			m4B = -4.0*fq;
			m9B = 2.0*fq;
			m10B = -4.0*fq;

			// q=2
			nread = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
			fq = distB[nread];  // reading the f2 data into register fq
            rhoB += fq;
			m1B -= 11.0*(fq);
			m2B -= 4.0*(fq);
			jxB -= fq;
			m4B += 4.0*(fq);
			m9B += 2.0*(fq);
			m10B -= 4.0*(fq);

			// q=3
			nread = neighborList[n+2*Np]; // neighbor 4
			fq = distB[nread];
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            GffA_x = -Gsc*nB_gradx;
            GffA_y = -Gsc*nB_grady;
            GffA_z = -Gsc*nB_gradz;
            GffB_x = -Gsc*nA_gradx;
            GffB_y = -Gsc*nA_grady;
            GffB_z = -Gsc*nA_gradz;
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
            FxA = nA*(-porosity*muA_eff/permA*ux_A - porosity*GeoFun/sqrt(permA)*u_mag*ux_A + porosity*Gx + GffA_x + GfsA_x);
            FyA = nA*(-porosity*muA_eff/permA*uy_A - porosity*GeoFun/sqrt(permA)*u_mag*uy_A + porosity*Gy + GffA_y + GfsA_y);
            FzA = nA*(-porosity*muA_eff/permA*uz_A - porosity*GeoFun/sqrt(permA)*u_mag*uz_A + porosity*Gz + GffA_z + GfsA_z);
            if (porosity==1.0){
                FxA=nA*(Gx + GffA_x + GfsA_x);
                FyA=nA*(Gy + GffA_y + GfsA_y);
                FzA=nA*(Gz + GffA_z + GfsA_z);
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
            FxB = nB*(-porosity*muB_eff/permB*ux_B - porosity*GeoFun/sqrt(permB)*u_mag*ux_B + porosity*Gx + GffB_x + GfsB_x);
            FyB = nB*(-porosity*muB_eff/permB*uy_B - porosity*GeoFun/sqrt(permB)*u_mag*uy_B + porosity*Gy + GffB_y + GfsB_y);
            FzB = nB*(-porosity*muB_eff/permB*uz_B - porosity*GeoFun/sqrt(permB)*u_mag*uz_B + porosity*Gz + GffB_z + GfsB_z);
            if (porosity==1.0){
                FxB=nB*(Gx + GffB_x + GfsB_x);
                FyB=nB*(Gy + GffB_y + GfsB_y);
                FzB=nB*(Gz + GffB_z + GfsB_z);
            }

            // Calculate barycentric velocity of the fluid mixture
            ux = (nA*ux_A+nB*ux_B)/(nA+nB);
            uy = (nA*uy_A+nB*uy_B)/(nA+nB);
            uz = (nA*uz_A+nB*uz_B)/(nA+nB);

            // ------------------- Fluid Component A -----------------------//
            rlx_setA = 1.0/tauA;
            rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
            //-------------------- MRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            //TODO need to incoporate porosity
			m1A = m1A + rlx_setA*((19*rhoA*(ux*ux+uy*uy+uz*uz) - 11*rhoA) - m1A)
                      + (1-0.5*rlx_setA)*38*(FxA*ux+FyA*uy+FzA*uz);
			m2A = m2A + rlx_setA*((3*rhoA - 5.5*rhoA*(ux*ux+uy*uy+uz*uz))- m2A)
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
			//m10A = m10A + rlx_setA*( - m10A)
            //            + (1-0.5*rlx_setA)*(-2*FxA*ux+FyA*uy+FzA*uz);
			m10A = m10A + rlx_setA*( -0.5*(rhoA*(2*ux*ux-uy*uy-uz*uz))- m10A)
                        + (1-0.5*rlx_setA)*(-2*FxA*ux+FyA*uy+FzA*uz);
			m11A = m11A + rlx_setA*((rhoA*(uy*uy-uz*uz)) - m11A)
                        + (1-0.5*rlx_setA)*(2*FyA*uy-2*FzA*uz);
			//m12A = m12A + rlx_setA*( - m12A)
            //            + (1-0.5*rlx_setA)*(-FyA*uy+FzA*uz);
			m12A = m12A + rlx_setA*( -0.5*(rhoA*(uy*uy-uz*uz))- m12A)
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
            fq = mrt_V1*rhoA-mrt_V2*m1A+mrt_V3*m2A;
            distA[n] = fq;

            // q = 1
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(jxA-m4A)+mrt_V6*(m9A-m10A);
            nread = neighborList[n+Np];
            distA[nread] = fq;

            // q=2
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(m4A-jxA)+mrt_V6*(m9A-m10A);
            nread = neighborList[n];
            distA[nread] = fq;

            // q = 3
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(jyA-m6A)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            nread = neighborList[n+3*Np];
            distA[nread] = fq;

            // q = 4
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(m6A-jyA)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            nread = neighborList[n+2*Np];
            distA[nread] = fq;

            // q = 5
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(jzA-m8A)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            nread = neighborList[n+5*Np];
            distA[nread] = fq;

            // q = 6
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(m8A-jzA)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            nread = neighborList[n+4*Np];
            distA[nread] = fq;

            // q = 7
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jyA)+0.025*(m4A+m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m16A-m17A);
            nread = neighborList[n+7*Np];
            distA[nread] = fq;

            // q = 8
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jyA)-0.025*(m4A+m6A) +mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m17A-m16A);
            nread = neighborList[n+6*Np];
            distA[nread] = fq;

            // q = 9
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jyA)+0.025*(m4A-m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A+0.125*(m16A+m17A);
            nread = neighborList[n+9*Np];
            distA[nread] = fq;

            // q = 10
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jxA)+0.025*(m6A-m4A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A-0.125*(m16A+m17A);
            nread = neighborList[n+8*Np];
            distA[nread] = fq;

            // q = 11
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jzA)+0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m18A-m16A);
            nread = neighborList[n+11*Np];
            distA[nread] = fq;

            // q = 12
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jzA)-0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m16A-m18A);
            nread = neighborList[n+10*Np];
            distA[nread]= fq;

            // q = 13
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jzA)+0.025*(m4A-m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A-0.125*(m16A+m18A);
            nread = neighborList[n+13*Np];
            distA[nread] = fq;

            // q= 14
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jxA)+0.025*(m8A-m4A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A+0.125*(m16A+m18A);
            nread = neighborList[n+12*Np];
            distA[nread] = fq;

            // q = 15
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA+jzA)+0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m17A-m18A);
            nread = neighborList[n+15*Np];
            distA[nread] = fq;

            // q = 16
            fq =  mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A-0.1*(jyA+jzA)-0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m18A-m17A);
            nread = neighborList[n+14*Np];
            distA[nread] = fq;

            // q = 17
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jzA)+0.025*(m6A-m8A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A+0.125*(m17A+m18A);
            nread = neighborList[n+17*Np];
            distA[nread] = fq;

            // q = 18
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jyA)+0.025*(m8A-m6A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A-0.125*(m17A+m18A);
            nread = neighborList[n+16*Np];
            distA[nread] = fq;
            //........................................................................

            // ------------------- Fluid Component B -----------------------//
            rlx_setA = 1.0/tauB;
            rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
            //-------------------- MRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            //TODO need to incoporate porosity
			m1B = m1B + rlx_setA*((19*rhoB*(ux*ux+uy*uy+uz*uz) - 11*rhoB) - m1B)
                      + (1-0.5*rlx_setA)*38*(FxB*ux+FyB*uy+FzB*uz);
			m2B = m2B + rlx_setA*((3*rhoB - 5.5*rhoB*(ux*ux+uy*uy+uz*uz))- m2B)
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
			//m10B = m10B + rlx_setA*( - m10B)
            //            + (1-0.5*rlx_setA)*(-2*FxB*ux+FyB*uy+FzB*uz);
			m10B = m10B + rlx_setA*( -0.5*(rhoB*(2*ux*ux-uy*uy-uz*uz))- m10B)
                        + (1-0.5*rlx_setA)*(-2*FxB*ux+FyB*uy+FzB*uz);
			m11B = m11B + rlx_setA*((rhoB*(uy*uy-uz*uz)) - m11B)
                        + (1-0.5*rlx_setA)*(2*FyB*uy-2*FzB*uz);
			//m12B = m12B + rlx_setA*( - m12B)
            //            + (1-0.5*rlx_setA)*(-FyB*uy+FzB*uz);
			m12B = m12B + rlx_setA*( -0.5*(rhoB*(uy*uy-uz*uz))- m12B)
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
            fq = mrt_V1*rhoB-mrt_V2*m1B+mrt_V3*m2B;
            distB[n] = fq;

            // q = 1
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(jxB-m4B)+mrt_V6*(m9B-m10B);
            nread = neighborList[n+Np];
            distB[nread] = fq;

            // q=2
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(m4B-jxB)+mrt_V6*(m9B-m10B);
            nread = neighborList[n];
            distB[nread] = fq;

            // q = 3
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(jyB-m6B)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            nread = neighborList[n+3*Np];
            distB[nread] = fq;

            // q = 4
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(m6B-jyB)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            nread = neighborList[n+2*Np];
            distB[nread] = fq;

            // q = 5
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(jzB-m8B)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            nread = neighborList[n+5*Np];
            distB[nread] = fq;

            // q = 6
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(m8B-jzB)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            nread = neighborList[n+4*Np];
            distB[nread] = fq;

            // q = 7
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jyB)+0.025*(m4B+m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m16B-m17B);
            nread = neighborList[n+7*Np];
            distB[nread] = fq;

            // q = 8
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jyB)-0.025*(m4B+m6B) +mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m17B-m16B);
            nread = neighborList[n+6*Np];
            distB[nread] = fq;

            // q = 9
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jyB)+0.025*(m4B-m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B+0.125*(m16B+m17B);
            nread = neighborList[n+9*Np];
            distB[nread] = fq;

            // q = 10
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jxB)+0.025*(m6B-m4B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B-0.125*(m16B+m17B);
            nread = neighborList[n+8*Np];
            distB[nread] = fq;

            // q = 11
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jzB)+0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m18B-m16B);
            nread = neighborList[n+11*Np];
            distB[nread] = fq;

            // q = 12
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jzB)-0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m16B-m18B);
            nread = neighborList[n+10*Np];
            distB[nread]= fq;

            // q = 13
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jzB)+0.025*(m4B-m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B-0.125*(m16B+m18B);
            nread = neighborList[n+13*Np];
            distB[nread] = fq;

            // q= 14
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jxB)+0.025*(m8B-m4B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B+0.125*(m16B+m18B);
            nread = neighborList[n+12*Np];
            distB[nread] = fq;

            // q = 15
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB+jzB)+0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m17B-m18B);
            nread = neighborList[n+15*Np];
            distB[nread] = fq;

            // q = 16
            fq =  mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B-0.1*(jyB+jzB)-0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m18B-m17B);
            nread = neighborList[n+14*Np];
            distB[nread] = fq;

            // q = 17
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jzB)+0.025*(m6B-m8B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B+0.125*(m17B+m18B);
            nread = neighborList[n+17*Np];
            distB[nread] = fq;

            // q = 18
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jyB)+0.025*(m8B-m6B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B-0.125*(m17B+m18B);
            nread = neighborList[n+16*Np];
            distB[nread] = fq;
            //........................................................................

            //Update velocity on device
            Velocity[0*Np+n] = ux;
            Velocity[1*Np+n] = uy;
            Velocity[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = (nA+nB+Gsc*nA*nB)/3.0;
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_AAeven_GreyscaleSC_MRT(int *Map,double *distA, double *distB, double *DenA,double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

    int ijk;
	int n;
	double vx,vy,vz,v_mag;
    double ux_A,uy_A,uz_A,ux_B,uy_B,uz_B,u_mag;
    double ux,uy,uz;
	// conserved momemnts
	double jxA,jyA,jzA;
	double jxB,jyB,jzB;
    double rhoA,rhoB;
    double nA,nB;
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
    double nA_gradx,nA_grady,nA_gradz;
    double nB_gradx,nB_grady,nB_gradz;
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
            ijk = Map[n];
            nA = DenA[ijk];
            nB = DenB[ijk];
            porosity = Poros[n];
            perm = Perm[n];
            permA = perm*nA/(nA+nB);//effective relative perm
            permB = perm*nB/(nA+nB);
			nA_gradx = DenGradA[n+0*Np];
			nA_grady = DenGradA[n+1*Np];
			nA_gradz = DenGradA[n+2*Np];
			nB_gradx = DenGradB[n+0*Np];
			nB_grady = DenGradB[n+1*Np];
			nB_gradz = DenGradB[n+2*Np];

            // ------------------- Fluid component A ---------------------------------//
            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
			// q=0
			fq = distA[n];
            rhoA = fq;
			m1A  = -30.0*fq;
			m2A  = 12.0*fq;

			// q=1
            fq = distA[2*Np+n];
            rhoA += fq;
			m1A -= 11.0*fq;
			m2A -= 4.0*fq;
			jxA = fq;
			m4A = -4.0*fq;
			m9A = 2.0*fq;
			m10A = -4.0*fq;

			// q=2
            fq = distA[1*Np+n];
            rhoA += fq;
			m1A -= 11.0*(fq);
			m2A -= 4.0*(fq);
			jxA -= fq;
			m4A += 4.0*(fq);
			m9A += 2.0*(fq);
			m10A -= 4.0*(fq);

			// q=3
            fq = distA[4*Np+n];
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoA += fq;
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
            rhoB = fq;
			m1B  = -30.0*fq;
			m2B  = 12.0*fq;

			// q=1
            fq = distB[2*Np+n];
            rhoB += fq;
			m1B -= 11.0*fq;
			m2B -= 4.0*fq;
			jxB = fq;
			m4B = -4.0*fq;
			m9B = 2.0*fq;
			m10B = -4.0*fq;

			// q=2
            fq = distB[1*Np+n];
            rhoB += fq;
			m1B -= 11.0*(fq);
			m2B -= 4.0*(fq);
			jxB -= fq;
			m4B += 4.0*(fq);
			m9B += 2.0*(fq);
			m10B -= 4.0*(fq);

			// q=3
            fq = distB[4*Np+n];
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            rhoB += fq;
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
            GffA_x = -Gsc*nB_gradx;
            GffA_y = -Gsc*nB_grady;
            GffA_z = -Gsc*nB_gradz;
            GffB_x = -Gsc*nA_gradx;
            GffB_y = -Gsc*nA_grady;
            GffB_z = -Gsc*nA_gradz;
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
            FxA = nA*(-porosity*muA_eff/permA*ux_A - porosity*GeoFun/sqrt(permA)*u_mag*ux_A + porosity*Gx + GffA_x + GfsA_x);
            FyA = nA*(-porosity*muA_eff/permA*uy_A - porosity*GeoFun/sqrt(permA)*u_mag*uy_A + porosity*Gy + GffA_y + GfsA_y);
            FzA = nA*(-porosity*muA_eff/permA*uz_A - porosity*GeoFun/sqrt(permA)*u_mag*uz_A + porosity*Gz + GffA_z + GfsA_z);
            if (porosity==1.0){
                FxA=nA*(Gx + GffA_x + GfsA_x);
                FyA=nA*(Gy + GffA_y + GfsA_y);
                FzA=nA*(Gz + GffA_z + GfsA_z);
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
            FxB = nB*(-porosity*muB_eff/permB*ux_B - porosity*GeoFun/sqrt(permB)*u_mag*ux_B + porosity*Gx + GffB_x + GfsB_x);
            FyB = nB*(-porosity*muB_eff/permB*uy_B - porosity*GeoFun/sqrt(permB)*u_mag*uy_B + porosity*Gy + GffB_y + GfsB_y);
            FzB = nB*(-porosity*muB_eff/permB*uz_B - porosity*GeoFun/sqrt(permB)*u_mag*uz_B + porosity*Gz + GffB_z + GfsB_z);
            if (porosity==1.0){
                FxB=nB*(Gx + GffB_x + GfsB_x);
                FyB=nB*(Gy + GffB_y + GfsB_y);
                FzB=nB*(Gz + GffB_z + GfsB_z);
            }

            // Calculate barycentric velocity of the fluid mixture
            ux = (nA*ux_A+nB*ux_B)/(nA+nB);
            uy = (nA*uy_A+nB*uy_B)/(nA+nB);
            uz = (nA*uz_A+nB*uz_B)/(nA+nB);


            // ------------------- Fluid Component A -----------------------//
            rlx_setA = 1.0/tauA;
            rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
            //-------------------- MRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            //TODO need to incoporate porosity
			m1A = m1A + rlx_setA*((19*rhoA*(ux*ux+uy*uy+uz*uz) - 11*rhoA) - m1A)
                      + (1-0.5*rlx_setA)*38*(FxA*ux+FyA*uy+FzA*uz);
			m2A = m2A + rlx_setA*((3*rhoA - 5.5*rhoA*(ux*ux+uy*uy+uz*uz))- m2A)
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
			//m10A = m10A + rlx_setA*( - m10A)
            //            + (1-0.5*rlx_setA)*(-2*FxA*ux+FyA*uy+FzA*uz);
			m10A = m10A + rlx_setA*( -0.5*(rhoA*(2*ux*ux-uy*uy-uz*uz))- m10A)
                        + (1-0.5*rlx_setA)*(-2*FxA*ux+FyA*uy+FzA*uz);
			m11A = m11A + rlx_setA*((rhoA*(uy*uy-uz*uz)) - m11A)
                        + (1-0.5*rlx_setA)*(2*FyA*uy-2*FzA*uz);
			//m12A = m12A + rlx_setA*( - m12A)
            //            + (1-0.5*rlx_setA)*(-FyA*uy+FzA*uz);
			m12A = m12A + rlx_setA*( -0.5*(rhoA*(uy*uy-uz*uz))- m12A)
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
            fq = mrt_V1*rhoA-mrt_V2*m1A+mrt_V3*m2A;
            distA[n] = fq;

            // q = 1
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(jxA-m4A)+mrt_V6*(m9A-m10A);
            distA[1*Np+n] = fq;

            // q=2
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(m4A-jxA)+mrt_V6*(m9A-m10A);
            distA[2*Np+n] = fq;

            // q = 3
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(jyA-m6A)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            distA[3*Np+n] = fq;

            // q = 4
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(m6A-jyA)+mrt_V7*(m10A-m9A)+mrt_V8*(m11A-m12A);
            distA[4*Np+n] = fq;

            // q = 5
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(jzA-m8A)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            distA[5*Np+n] = fq;

            // q = 6
            fq = mrt_V1*rhoA-mrt_V4*m1A-mrt_V5*m2A+0.1*(m8A-jzA)+mrt_V7*(m10A-m9A)+mrt_V8*(m12A-m11A);
            distA[6*Np+n] = fq;

            // q = 7
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jyA)+0.025*(m4A+m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m16A-m17A);
            distA[7*Np+n] = fq;

            // q = 8
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jyA)-0.025*(m4A+m6A) +mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A+0.25*m13A+0.125*(m17A-m16A);
            distA[8*Np+n] = fq;

            // q = 9
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jyA)+0.025*(m4A-m6A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A+0.125*(m16A+m17A);
            distA[9*Np+n] = fq;

            // q = 10
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jxA)+0.025*(m6A-m4A)+mrt_V7*m9A+mrt_V11*m10A+mrt_V8*m11A+mrt_V12*m12A-0.25*m13A-0.125*(m16A+m17A);
            distA[10*Np+n] = fq;

            // q = 11
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA+jzA)+0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m18A-m16A);
            distA[11*Np+n] = fq;

            // q = 12
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A-0.1*(jxA+jzA)-0.025*(m4A+m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A+0.25*m15A+0.125*(m16A-m18A);
            distA[12*Np+n] = fq;

            // q = 13
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jxA-jzA)+0.025*(m4A-m8A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A-0.125*(m16A+m18A);
            distA[13*Np+n] = fq;

            // q= 14
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jxA)+0.025*(m8A-m4A)+mrt_V7*m9A+mrt_V11*m10A-mrt_V8*m11A-mrt_V12*m12A-0.25*m15A+0.125*(m16A+m18A);
            distA[14*Np+n] = fq;

            // q = 15
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA+jzA)+0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m17A-m18A);
            distA[15*Np+n] = fq;

            // q = 16
            fq =  mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A-0.1*(jyA+jzA)-0.025*(m6A+m8A)-mrt_V6*m9A-mrt_V7*m10A+0.25*m14A+0.125*(m18A-m17A);
            distA[16*Np+n] = fq;

            // q = 17
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jyA-jzA)+0.025*(m6A-m8A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A+0.125*(m17A+m18A);
            distA[17*Np+n] = fq;

            // q = 18
            fq = mrt_V1*rhoA+mrt_V9*m1A+mrt_V10*m2A+0.1*(jzA-jyA)+0.025*(m8A-m6A)-mrt_V6*m9A-mrt_V7*m10A-0.25*m14A-0.125*(m17A+m18A);
            distA[18*Np+n] = fq;
            //........................................................................
           

            // ------------------- Fluid Component B -----------------------//
            rlx_setA = 1.0/tauB;
            rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
            //-------------------- MRT collison where body force has NO higher-order terms -------------//
            //..............carry out relaxation process...............................................
            //TODO need to incoporate porosity
			m1B = m1B + rlx_setA*((19*rhoB*(ux*ux+uy*uy+uz*uz) - 11*rhoB) - m1B)
                      + (1-0.5*rlx_setA)*38*(FxB*ux+FyB*uy+FzB*uz);
			m2B = m2B + rlx_setA*((3*rhoB - 5.5*rhoB*(ux*ux+uy*uy+uz*uz))- m2B)
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
			//m10B = m10B + rlx_setA*( - m10B)
            //            + (1-0.5*rlx_setA)*(-2*FxB*ux+FyB*uy+FzB*uz);
			m10B = m10B + rlx_setA*( -0.5*(rhoB*(2*ux*ux-uy*uy-uz*uz))- m10B)
                        + (1-0.5*rlx_setA)*(-2*FxB*ux+FyB*uy+FzB*uz);
			m11B = m11B + rlx_setA*((rhoB*(uy*uy-uz*uz)) - m11B)
                        + (1-0.5*rlx_setA)*(2*FyB*uy-2*FzB*uz);
			//m12B = m12B + rlx_setA*( - m12B)
            //            + (1-0.5*rlx_setA)*(-FyB*uy+FzB*uz);
			m12B = m12B + rlx_setA*( -0.5*(rhoB*(uy*uy-uz*uz))- m12B)
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
            fq = mrt_V1*rhoB-mrt_V2*m1B+mrt_V3*m2B;
            distB[n] = fq;

            // q = 1
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(jxB-m4B)+mrt_V6*(m9B-m10B);
            distB[1*Np+n] = fq;

            // q=2
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(m4B-jxB)+mrt_V6*(m9B-m10B);
            distB[2*Np+n] = fq;

            // q = 3
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(jyB-m6B)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            distB[3*Np+n] = fq;

            // q = 4
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(m6B-jyB)+mrt_V7*(m10B-m9B)+mrt_V8*(m11B-m12B);
            distB[4*Np+n] = fq;

            // q = 5
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(jzB-m8B)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            distB[5*Np+n] = fq;

            // q = 6
            fq = mrt_V1*rhoB-mrt_V4*m1B-mrt_V5*m2B+0.1*(m8B-jzB)+mrt_V7*(m10B-m9B)+mrt_V8*(m12B-m11B);
            distB[6*Np+n] = fq;

            // q = 7
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jyB)+0.025*(m4B+m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m16B-m17B);
            distB[7*Np+n] = fq;

            // q = 8
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jyB)-0.025*(m4B+m6B) +mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B+0.25*m13B+0.125*(m17B-m16B);
            distB[8*Np+n] = fq;

            // q = 9
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jyB)+0.025*(m4B-m6B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B+0.125*(m16B+m17B);
            distB[9*Np+n] = fq;

            // q = 10
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jxB)+0.025*(m6B-m4B)+mrt_V7*m9B+mrt_V11*m10B+mrt_V8*m11B+mrt_V12*m12B-0.25*m13B-0.125*(m16B+m17B);
            distB[10*Np+n] = fq;

            // q = 11
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB+jzB)+0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m18B-m16B);
            distB[11*Np+n] = fq;

            // q = 12
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B-0.1*(jxB+jzB)-0.025*(m4B+m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B+0.25*m15B+0.125*(m16B-m18B);
            distB[12*Np+n] = fq;

            // q = 13
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jxB-jzB)+0.025*(m4B-m8B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B-0.125*(m16B+m18B);
            distB[13*Np+n] = fq;

            // q= 14
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jxB)+0.025*(m8B-m4B)+mrt_V7*m9B+mrt_V11*m10B-mrt_V8*m11B-mrt_V12*m12B-0.25*m15B+0.125*(m16B+m18B);
            distB[14*Np+n] = fq;

            // q = 15
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB+jzB)+0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m17B-m18B);
            distB[15*Np+n] = fq;

            // q = 16
            fq =  mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B-0.1*(jyB+jzB)-0.025*(m6B+m8B)-mrt_V6*m9B-mrt_V7*m10B+0.25*m14B+0.125*(m18B-m17B);
            distB[16*Np+n] = fq;

            // q = 17
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jyB-jzB)+0.025*(m6B-m8B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B+0.125*(m17B+m18B);
            distB[17*Np+n] = fq;

            // q = 18
            fq = mrt_V1*rhoB+mrt_V9*m1B+mrt_V10*m2B+0.1*(jzB-jyB)+0.025*(m8B-m6B)-mrt_V6*m9B-mrt_V7*m10B-0.25*m14B-0.125*(m17B+m18B);
            distB[18*Np+n] = fq;
            //........................................................................

            //Update velocity on device
            Velocity[0*Np+n] = ux;
            Velocity[1*Np+n] = uy;
            Velocity[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = (nA+nB+Gsc*nA*nB)/3.0;
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_AAodd_GreyscaleSC_BGK(int *neighborList, int *Map, double *distA, double *distB, double *DenA, double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

	int n;
    int ijk;
	int nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,nr15,nr16,nr17,nr18;
	double vx,vy,vz,v_mag;
    double ux_A,uy_A,uz_A,ux_B,uy_B,uz_B,u_mag;
    double ux,uy,uz;
    double rhoA,rhoB;
	double jxA,jyA,jzA;
	double jxB,jyB,jzB;
	// distribution functions
	double f0A,f1A,f2A,f3A,f4A,f5A,f6A,f7A,f8A,f9A,f10A,f11A,f12A,f13A,f14A,f15A,f16A,f17A,f18A;
	double f0B,f1B,f2B,f3B,f4B,f5B,f6B,f7B,f8B,f9B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B;
    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double permA,permB;//effective relative perm
    double c0, c1; //Guo's model parameters
    double muA_eff = (tauA_eff-0.5)/3.0;//kinematic viscosity
    double muB_eff = (tauB_eff-0.5)/3.0;//kinematic viscosity
    double FxA, FyA, FzA;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
    double FxB, FyB, FzB;
    double tau,rlx;
    double phi;//phase field indicator
    double rhoA_gradx,rhoA_grady,rhoA_gradz;
    double rhoB_gradx,rhoB_grady,rhoB_gradz;
    double GffA_x,GffA_y,GffA_z;
    double GfsA_x,GfsA_y,GfsA_z;
    double GffB_x,GffB_y,GffB_z;
    double GfsB_x,GfsB_y,GfsB_z;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
	    //........Get 1-D index for this thread....................
	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){		

            // Load common parameters shared by two fluid components
            ijk = Map[n];
            rhoA = DenA[ijk];
            rhoB = DenB[ijk];
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
            phi = (rhoA-rhoB)/(rhoA+rhoB);
            tau = tauA+0.5*(1.0-phi)*(tauB-tauA); 
            rlx = 1.0/tau;
            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
            // q=0
            f0A = distA[n];
            f0B = distB[n];
            // q=1
            nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
            f1A = distA[nr1]; // reading the f1 data into register fq
            f1B = distB[nr1]; // reading the f1 data into register fq

            nr2 = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
            f2A = distA[nr2];  // reading the f2 data into register fq
            f2B = distB[nr2];  // reading the f2 data into register fq

            // q=3
            nr3 = neighborList[n+2*Np]; // neighbor 4
            f3A = distA[nr3];
            f3B = distB[nr3];

            // q = 4
            nr4 = neighborList[n+3*Np]; // neighbor 3
            f4A = distA[nr4];
            f4B = distB[nr4];

            // q=5
            nr5 = neighborList[n+4*Np];
            f5A = distA[nr5];
            f5B = distB[nr5];

            // q = 6
            nr6 = neighborList[n+5*Np];
            f6A = distA[nr6];
            f6B = distB[nr6];
            
            // q=7
            nr7 = neighborList[n+6*Np];
            f7A = distA[nr7];
            f7B = distB[nr7];

            // q = 8
            nr8 = neighborList[n+7*Np];
            f8A = distA[nr8];
            f8B = distB[nr8];

            // q=9
            nr9 = neighborList[n+8*Np];
            f9A = distA[nr9];
            f9B = distB[nr9];

            // q = 10
            nr10 = neighborList[n+9*Np];
            f10A = distA[nr10];
            f10B = distB[nr10];

            // q=11
            nr11 = neighborList[n+10*Np];
            f11A = distA[nr11];
            f11B = distB[nr11];

            // q=12
            nr12 = neighborList[n+11*Np];
            f12A = distA[nr12];
            f12B = distB[nr12];

            // q=13
            nr13 = neighborList[n+12*Np];
            f13A = distA[nr13];
            f13B = distB[nr13];

            // q=14
            nr14 = neighborList[n+13*Np];
            f14A = distA[nr14];
            f14B = distB[nr14];

            // q=15
            nr15 = neighborList[n+14*Np];
            f15A = distA[nr15];
            f15B = distB[nr15];

            // q=16
            nr16 = neighborList[n+15*Np];
            f16A = distA[nr16];
            f16B = distB[nr16];

            // q=17
            //fq = dist[18*Np+n];
            nr17 = neighborList[n+16*Np];
            f17A = distA[nr17];
            f17B = distB[nr17];

            // q=18
            nr18 = neighborList[n+17*Np];
            f18A = distA[nr18];
            f18B = distB[nr18];
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
            jxA = f1A-f2A+f7A-f8A+f9A-f10A+f11A-f12A+f13A-f14A;
            jyA = f3A-f4A+f7A-f8A-f9A+f10A+f15A-f16A+f17A-f18A;
            jzA = f5A-f6A+f11A-f12A-f13A+f14A+f15A-f16A-f17A+f18A;

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
            jxB = f1B-f2B+f7B-f8B+f9B-f10B+f11B-f12B+f13B-f14B;
            jyB = f3B-f4B+f7B-f8B-f9B+f10B+f15B-f16B+f17B-f18B;
            jzB = f5B-f6B+f11B-f12B-f13B+f14B+f15B-f16B-f17B+f18B;

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

            //..............carry out relaxation process...............................................
            // ------------------- Fluid Component A -----------------------//
            // q=0
            distA[n] = f0A*(1.0-rlx) + rlx*0.3333333333333333*rhoA*(1. - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                      + 0.3333333333333333*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(0. - (3.*uy)/porosity) + FzA*(0. - (3.*uz)/porosity));

            // q = 1
            distA[nr2] = f1A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 + 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(3. + (6.*ux)/porosity) + FyA*(0. - (3.*uy)/porosity) + FzA*(0. - (3.*uz)/porosity));

            // q=2
            distA[nr1] = f2A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 - 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(-3. + (6.*ux)/porosity) + FyA*(0. - (3.*uy)/porosity) + FzA*(0. - (3.*uz)/porosity));

            // q = 3
            distA[nr4] = f3A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 + 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(3. + (6.*uy)/porosity) + FzA*(0. - (3.*uz)/porosity));

            // q = 4
            distA[nr3] = f4A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 - 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)  
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(-3. + (6.*uy)/porosity) + FzA*(0. - (3.*uz)/porosity));

            // q = 5
            distA[nr6] = f5A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 + 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(0. - (3.*uy)/porosity) + FzA*(3. + (6.*uz)/porosity));

            // q = 6
            distA[nr5] = f6A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 - 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux+ uy*uy + uz*uz))/porosity) 
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(0. - (3.*uy)/porosity) + FzA*(-3. + (6.*uz)/porosity));

            // q = 7
            distA[nr8] = f7A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(ux + uy) + (4.5*(ux + uy)*(ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(3. - (3.*ux)/porosity + (9.*(ux + uy))/porosity) + FyA*(3. - (3.*uy)/porosity + (9.*(ux + uy))/porosity) + 
      FzA*(0. - (3.*uz)/porosity));

            // q = 8
            distA[nr7] = f8A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-ux - uy) + (4.5*(-ux - uy)*(-ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(-3. - (3.*ux)/porosity - (9.*(-ux - uy))/porosity) + FyA*(-3. - (9.*(-ux - uy))/porosity - (3.*uy)/porosity) + 
      FzA*(0. - (3.*uz)/porosity));

            // q = 9
            distA[nr10] = f9A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(ux - uy) + (4.5*(ux - uy)*(ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(3. - (3.*ux)/porosity + (9.*(ux - uy))/porosity) + FyA*(-3. - (9.*(ux - uy))/porosity - (3.*uy)/porosity) + 
      FzA*(0. - (3.*uz)/porosity));

            // q = 10
            distA[nr9] = f10A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-ux + uy) + (4.5*(-ux + uy)*(-ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(-3. - (3.*ux)/porosity - (9.*(-ux + uy))/porosity) + FyA*(3. - (3.*uy)/porosity + (9.*(-ux + uy))/porosity) + 
      FzA*(0. - (3.*uz)/porosity));

            // q = 11
            distA[nr12] = f11A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(ux + uz) + (4.5*(ux + uz)*(ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyA*(0. - (3.*uy)/porosity) + FxA*(3. - (3.*ux)/porosity + (9.*(ux + uz))/porosity) + 
      FzA*(3. - (3.*uz)/porosity + (9.*(ux + uz))/porosity));

            // q = 12
            distA[nr11] = f12A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-ux - uz) + (4.5*(-ux - uz)*(-ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyA*(0. - (3.*uy)/porosity) + FxA*(-3. - (3.*ux)/porosity - (9.*(-ux - uz))/porosity) + 
      FzA*(-3. - (9.*(-ux - uz))/porosity - (3.*uz)/porosity));

            // q = 13
            distA[nr14] = f13A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(ux - uz) + (4.5*(ux - uz)*(ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyA*(0. - (3.*uy)/porosity) + FxA*(3. - (3.*ux)/porosity + (9.*(ux - uz))/porosity) + 
      FzA*(-3. - (9.*(ux - uz))/porosity - (3.*uz)/porosity));

            // q= 14
            distA[nr13] = f14A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-ux + uz) + (4.5*(-ux + uz)*(-ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyA*(0. - (3.*uy)/porosity) + FxA*(-3. - (3.*ux)/porosity - (9.*(-ux + uz))/porosity) + 
      FzA*(3. - (3.*uz)/porosity + (9.*(-ux + uz))/porosity));

            // q = 15
            distA[nr16] = f15A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(uy + uz) + (4.5*(uy + uz)*(uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(3. - (3.*uy)/porosity + (9.*(uy + uz))/porosity) + 
      FzA*(3. - (3.*uz)/porosity + (9.*(uy + uz))/porosity));

            // q = 16
            distA[nr15] = f16A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-uy - uz) + (4.5*(-uy - uz)*(-uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(-3. - (3.*uy)/porosity - (9.*(-uy - uz))/porosity) + 
      FzA*(-3. - (9.*(-uy - uz))/porosity - (3.*uz)/porosity));

            // q = 17
            distA[nr18] = f17A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(uy - uz) + (4.5*(uy - uz)*(uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(3. - (3.*uy)/porosity + (9.*(uy - uz))/porosity) + 
      FzA*(-3. - (9.*(uy - uz))/porosity - (3.*uz)/porosity));

            // q = 18
            distA[nr17] = f18A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-uy + uz) + (4.5*(-uy + uz)*(-uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(-3. - (3.*uy)/porosity - (9.*(-uy + uz))/porosity) + 
      FzA*(3. - (3.*uz)/porosity + (9.*(-uy + uz))/porosity));


            // ------------------- Fluid Component B -----------------------//
            // q=0
            distB[n] = f0B*(1.0-rlx) + rlx*0.3333333333333333*rhoB*(1. - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                      + 0.3333333333333333*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(0. - (3.*uy)/porosity) + FzB*(0. - (3.*uz)/porosity));

            // q = 1
            distB[nr2] = f1B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 + 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(3. + (6.*ux)/porosity) + FyB*(0. - (3.*uy)/porosity) + FzB*(0. - (3.*uz)/porosity));

            // q=2
            distB[nr1] = f2B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 - 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(-3. + (6.*ux)/porosity) + FyB*(0. - (3.*uy)/porosity) + FzB*(0. - (3.*uz)/porosity));

            // q = 3
            distB[nr4] = f3B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 + 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(3. + (6.*uy)/porosity) + FzB*(0. - (3.*uz)/porosity));

            // q = 4
            distB[nr3] = f4B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 - 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)  
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(-3. + (6.*uy)/porosity) + FzB*(0. - (3.*uz)/porosity));

            // q = 5
            distB[nr6] = f5B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 + 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(0. - (3.*uy)/porosity) + FzB*(3. + (6.*uz)/porosity));

            // q = 6
            distB[nr5] = f6B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 - 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux+ uy*uy + uz*uz))/porosity) 
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(0. - (3.*uy)/porosity) + FzB*(-3. + (6.*uz)/porosity));

            // q = 7
            distB[nr8] = f7B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(ux + uy) + (4.5*(ux + uy)*(ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(3. - (3.*ux)/porosity + (9.*(ux + uy))/porosity) + FyB*(3. - (3.*uy)/porosity + (9.*(ux + uy))/porosity) + 
      FzB*(0. - (3.*uz)/porosity));

            // q = 8
            distB[nr7] = f8B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-ux - uy) + (4.5*(-ux - uy)*(-ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(-3. - (3.*ux)/porosity - (9.*(-ux - uy))/porosity) + FyB*(-3. - (9.*(-ux - uy))/porosity - (3.*uy)/porosity) + 
      FzB*(0. - (3.*uz)/porosity));

            // q = 9
            distB[nr10] = f9B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(ux - uy) + (4.5*(ux - uy)*(ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(3. - (3.*ux)/porosity + (9.*(ux - uy))/porosity) + FyB*(-3. - (9.*(ux - uy))/porosity - (3.*uy)/porosity) + 
      FzB*(0. - (3.*uz)/porosity));

            // q = 10
            distB[nr9] = f10B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-ux + uy) + (4.5*(-ux + uy)*(-ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(-3. - (3.*ux)/porosity - (9.*(-ux + uy))/porosity) + FyB*(3. - (3.*uy)/porosity + (9.*(-ux + uy))/porosity) + 
      FzB*(0. - (3.*uz)/porosity));

            // q = 11
            distB[nr12] = f11B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(ux + uz) + (4.5*(ux + uz)*(ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyB*(0. - (3.*uy)/porosity) + FxB*(3. - (3.*ux)/porosity + (9.*(ux + uz))/porosity) + 
      FzB*(3. - (3.*uz)/porosity + (9.*(ux + uz))/porosity));

            // q = 12
            distB[nr11] = f12B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-ux - uz) + (4.5*(-ux - uz)*(-ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyB*(0. - (3.*uy)/porosity) + FxB*(-3. - (3.*ux)/porosity - (9.*(-ux - uz))/porosity) + 
      FzB*(-3. - (9.*(-ux - uz))/porosity - (3.*uz)/porosity));

            // q = 13
            distB[nr14] = f13B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(ux - uz) + (4.5*(ux - uz)*(ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyB*(0. - (3.*uy)/porosity) + FxB*(3. - (3.*ux)/porosity + (9.*(ux - uz))/porosity) + 
      FzB*(-3. - (9.*(ux - uz))/porosity - (3.*uz)/porosity));

            // q= 14
            distB[nr13] = f14B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-ux + uz) + (4.5*(-ux + uz)*(-ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyB*(0. - (3.*uy)/porosity) + FxB*(-3. - (3.*ux)/porosity - (9.*(-ux + uz))/porosity) + 
      FzB*(3. - (3.*uz)/porosity + (9.*(-ux + uz))/porosity));

            // q = 15
            distB[nr16] = f15B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(uy + uz) + (4.5*(uy + uz)*(uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(3. - (3.*uy)/porosity + (9.*(uy + uz))/porosity) + 
      FzB*(3. - (3.*uz)/porosity + (9.*(uy + uz))/porosity));

            // q = 16
            distB[nr15] = f16B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-uy - uz) + (4.5*(-uy - uz)*(-uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(-3. - (3.*uy)/porosity - (9.*(-uy - uz))/porosity) + 
      FzB*(-3. - (9.*(-uy - uz))/porosity - (3.*uz)/porosity));

            // q = 17
            distB[nr18] = f17B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(uy - uz) + (4.5*(uy - uz)*(uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(3. - (3.*uy)/porosity + (9.*(uy - uz))/porosity) + 
      FzB*(-3. - (9.*(uy - uz))/porosity - (3.*uz)/porosity));

            // q = 18
            distB[nr17] = f18B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-uy + uz) + (4.5*(-uy + uz)*(-uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(-3. - (3.*uy)/porosity - (9.*(-uy + uz))/porosity) + 
      FzB*(3. - (3.*uz)/porosity + (9.*(-uy + uz))/porosity));


            //Update velocity on device
            Velocity[0*Np+n] = ux;
            Velocity[1*Np+n] = uy;
            Velocity[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = (rhoA+rhoB+Gsc*rhoA*rhoB)/3.0;
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_AAeven_GreyscaleSC_BGK(int *Map, double *distA, double *distB, double *DenA, double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

	int n;
    int ijk;
	double vx,vy,vz,v_mag;
    double ux_A,uy_A,uz_A,ux_B,uy_B,uz_B,u_mag;
    double ux,uy,uz;
    double rhoA,rhoB;
	double jxA,jyA,jzA;
	double jxB,jyB,jzB;
	// distribution functions
	double f0A,f1A,f2A,f3A,f4A,f5A,f6A,f7A,f8A,f9A,f10A,f11A,f12A,f13A,f14A,f15A,f16A,f17A,f18A;
	double f0B,f1B,f2B,f3B,f4B,f5B,f6B,f7B,f8B,f9B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B;
    double GeoFun=0.0;//geometric function from Guo's PRE 66, 036304 (2002)
    double porosity;
    double perm;//voxel permeability
    double permA,permB;//effective relative perm
    double c0, c1; //Guo's model parameters
    double muA_eff = (tauA_eff-0.5)/3.0;//kinematic viscosity
    double muB_eff = (tauB_eff-0.5)/3.0;//kinematic viscosity
    double FxA, FyA, FzA;//The total body force including Brinkman force and user-specified (Gx,Gy,Gz)
    double FxB, FyB, FzB;
    double tau,rlx;
    double phi;//phase field indicator
    double rhoA_gradx,rhoA_grady,rhoA_gradz;
    double rhoB_gradx,rhoB_grady,rhoB_gradz;
    double GffA_x,GffA_y,GffA_z;
    double GfsA_x,GfsA_y,GfsA_z;
    double GffB_x,GffB_y,GffB_z;
    double GfsB_x,GfsB_y,GfsB_z;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
	    //........Get 1-D index for this thread....................
	    n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;

		if ( n<finish ){		

            // Load common parameters shared by two fluid components
            ijk = Map[n];
            rhoA = DenA[ijk];
            rhoB = DenB[ijk];
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
            phi = (rhoA-rhoB)/(rhoA+rhoB);
            tau = tauA+0.5*(1.0-phi)*(tauB-tauA);
            rlx = 1.0/tau;

            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................
            // fluid component A
            f0A  = distA[n];
            f1A  = distA[2*Np+n];
            f2A  = distA[1*Np+n];
            f3A  = distA[4*Np+n];
            f4A  = distA[3*Np+n];
            f5A  = distA[6*Np+n];
            f6A  = distA[5*Np+n];
            f7A  = distA[8*Np+n];
            f8A  = distA[7*Np+n];
            f9A  = distA[10*Np+n];
            f10A = distA[9*Np+n];
            f11A = distA[12*Np+n];
            f12A = distA[11*Np+n];
            f13A = distA[14*Np+n];
            f14A = distA[13*Np+n];
            f15A = distA[16*Np+n];
            f16A = distA[15*Np+n];
            f17A = distA[18*Np+n];
            f18A = distA[17*Np+n];
            // fluid component B
            f0B  = distB[n];
            f1B  = distB[2*Np+n];
            f2B  = distB[1*Np+n];
            f3B  = distB[4*Np+n];
            f4B  = distB[3*Np+n];
            f5B  = distB[6*Np+n];
            f6B  = distB[5*Np+n];
            f7B  = distB[8*Np+n];
            f8B  = distB[7*Np+n];
            f9B  = distB[10*Np+n];
            f10B = distB[9*Np+n];
            f11B = distB[12*Np+n];
            f12B = distB[11*Np+n];
            f13B = distB[14*Np+n];
            f14B = distB[13*Np+n];
            f15B = distB[16*Np+n];
            f16B = distB[15*Np+n];
            f17B = distB[18*Np+n];
            f18B = distB[17*Np+n];
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
            jxA = f1A-f2A+f7A-f8A+f9A-f10A+f11A-f12A+f13A-f14A;
            jyA = f3A-f4A+f7A-f8A-f9A+f10A+f15A-f16A+f17A-f18A;
            jzA = f5A-f6A+f11A-f12A-f13A+f14A+f15A-f16A-f17A+f18A;

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
            jxB = f1B-f2B+f7B-f8B+f9B-f10B+f11B-f12B+f13B-f14B;
            jyB = f3B-f4B+f7B-f8B-f9B+f10B+f15B-f16B+f17B-f18B;
            jzB = f5B-f6B+f11B-f12B-f13B+f14B+f15B-f16B-f17B+f18B;

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

            //..............carry out relaxation process...............................................
            // ------------------- Fluid Component A -----------------------//
            // q=0
            distA[n] = f0A*(1.0-rlx) + rlx*0.3333333333333333*rhoA*(1. - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                      + 0.3333333333333333*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(0. - (3.*uy)/porosity) + FzA*(0. - (3.*uz)/porosity));

            // q = 1
            distA[1*Np+n] = f1A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 + 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(3. + (6.*ux)/porosity) + FyA*(0. - (3.*uy)/porosity) + FzA*(0. - (3.*uz)/porosity));

            // q=2
            distA[2*Np+n] = f2A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 - 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(-3. + (6.*ux)/porosity) + FyA*(0. - (3.*uy)/porosity) + FzA*(0. - (3.*uz)/porosity));

            // q = 3
            distA[3*Np+n] = f3A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 + 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(3. + (6.*uy)/porosity) + FzA*(0. - (3.*uz)/porosity));

            // q = 4
            distA[4*Np+n] = f4A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 - 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)  
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(-3. + (6.*uy)/porosity) + FzA*(0. - (3.*uz)/porosity));

            // q = 5
            distA[5*Np+n] = f5A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 + 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(0. - (3.*uy)/porosity) + FzA*(3. + (6.*uz)/porosity));

            // q = 6
            distA[6*Np+n] = f6A*(1.0-rlx) + rlx*0.05555555555555555*rhoA*(1 - 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux+ uy*uy + uz*uz))/porosity) 
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(0. - (3.*uy)/porosity) + FzA*(-3. + (6.*uz)/porosity));

            // q = 7
            distA[7*Np+n] = f7A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(ux + uy) + (4.5*(ux + uy)*(ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(3. - (3.*ux)/porosity + (9.*(ux + uy))/porosity) + FyA*(3. - (3.*uy)/porosity + (9.*(ux + uy))/porosity) + 
      FzA*(0. - (3.*uz)/porosity));

            // q = 8
            distA[8*Np+n] = f8A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-ux - uy) + (4.5*(-ux - uy)*(-ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(-3. - (3.*ux)/porosity - (9.*(-ux - uy))/porosity) + FyA*(-3. - (9.*(-ux - uy))/porosity - (3.*uy)/porosity) + 
      FzA*(0. - (3.*uz)/porosity));

            // q = 9
            distA[9*Np+n] = f9A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(ux - uy) + (4.5*(ux - uy)*(ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(3. - (3.*ux)/porosity + (9.*(ux - uy))/porosity) + FyA*(-3. - (9.*(ux - uy))/porosity - (3.*uy)/porosity) + 
      FzA*(0. - (3.*uz)/porosity));

            // q = 10
            distA[10*Np+n] = f10A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-ux + uy) + (4.5*(-ux + uy)*(-ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(-3. - (3.*ux)/porosity - (9.*(-ux + uy))/porosity) + FyA*(3. - (3.*uy)/porosity + (9.*(-ux + uy))/porosity) + 
      FzA*(0. - (3.*uz)/porosity));

            // q = 11
            distA[11*Np+n] = f11A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(ux + uz) + (4.5*(ux + uz)*(ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyA*(0. - (3.*uy)/porosity) + FxA*(3. - (3.*ux)/porosity + (9.*(ux + uz))/porosity) + 
      FzA*(3. - (3.*uz)/porosity + (9.*(ux + uz))/porosity));

            // q = 12
            distA[12*Np+n] = f12A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-ux - uz) + (4.5*(-ux - uz)*(-ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyA*(0. - (3.*uy)/porosity) + FxA*(-3. - (3.*ux)/porosity - (9.*(-ux - uz))/porosity) + 
      FzA*(-3. - (9.*(-ux - uz))/porosity - (3.*uz)/porosity));

            // q = 13
            distA[13*Np+n] = f13A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(ux - uz) + (4.5*(ux - uz)*(ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyA*(0. - (3.*uy)/porosity) + FxA*(3. - (3.*ux)/porosity + (9.*(ux - uz))/porosity) + 
      FzA*(-3. - (9.*(ux - uz))/porosity - (3.*uz)/porosity));

            // q= 14
            distA[14*Np+n] = f14A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-ux + uz) + (4.5*(-ux + uz)*(-ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyA*(0. - (3.*uy)/porosity) + FxA*(-3. - (3.*ux)/porosity - (9.*(-ux + uz))/porosity) + 
      FzA*(3. - (3.*uz)/porosity + (9.*(-ux + uz))/porosity));

            // q = 15
            distA[15*Np+n] = f15A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(uy + uz) + (4.5*(uy + uz)*(uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(3. - (3.*uy)/porosity + (9.*(uy + uz))/porosity) + 
      FzA*(3. - (3.*uz)/porosity + (9.*(uy + uz))/porosity));

            // q = 16
            distA[16*Np+n] = f16A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-uy - uz) + (4.5*(-uy - uz)*(-uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(-3. - (3.*uy)/porosity - (9.*(-uy - uz))/porosity) + 
      FzA*(-3. - (9.*(-uy - uz))/porosity - (3.*uz)/porosity));

            // q = 17
            distA[17*Np+n] = f17A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(uy - uz) + (4.5*(uy - uz)*(uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(3. - (3.*uy)/porosity + (9.*(uy - uz))/porosity) + 
      FzA*(-3. - (9.*(uy - uz))/porosity - (3.*uz)/porosity));

            // q = 18
            distA[18*Np+n] = f18A*(1.0-rlx) + rlx*0.027777777777777776*rhoA*(1 + 3.*(-uy + uz) + (4.5*(-uy + uz)*(-uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxA*(0. - (3.*ux)/porosity) + FyA*(-3. - (3.*uy)/porosity - (9.*(-uy + uz))/porosity) + 
      FzA*(3. - (3.*uz)/porosity + (9.*(-uy + uz))/porosity));


            // ------------------- Fluid Component B -----------------------//
            // q=0
            distB[n] = f0B*(1.0-rlx) + rlx*0.3333333333333333*rhoB*(1. - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                      + 0.3333333333333333*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(0. - (3.*uy)/porosity) + FzB*(0. - (3.*uz)/porosity));

            // q = 1
            distB[1*Np+n] = f1B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 + 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(3. + (6.*ux)/porosity) + FyB*(0. - (3.*uy)/porosity) + FzB*(0. - (3.*uz)/porosity));

            // q=2
            distB[2*Np+n] = f2B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 - 3.*ux + (4.5*ux*ux)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(-3. + (6.*ux)/porosity) + FyB*(0. - (3.*uy)/porosity) + FzB*(0. - (3.*uz)/porosity));

            // q = 3
            distB[3*Np+n] = f3B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 + 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(3. + (6.*uy)/porosity) + FzB*(0. - (3.*uz)/porosity));

            // q = 4
            distB[4*Np+n] = f4B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 - 3.*uy + (4.5*uy*uy)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)  
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(-3. + (6.*uy)/porosity) + FzB*(0. - (3.*uz)/porosity));

            // q = 5
            distB[5*Np+n] = f5B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 + 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(0. - (3.*uy)/porosity) + FzB*(3. + (6.*uz)/porosity));

            // q = 6
            distB[6*Np+n] = f6B*(1.0-rlx) + rlx*0.05555555555555555*rhoB*(1 - 3.*uz + (4.5*uz*uz)/porosity - (1.5*(ux*ux+ uy*uy + uz*uz))/porosity) 
                    +0.05555555555555555*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(0. - (3.*uy)/porosity) + FzB*(-3. + (6.*uz)/porosity));

            // q = 7
            distB[7*Np+n] = f7B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(ux + uy) + (4.5*(ux + uy)*(ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(3. - (3.*ux)/porosity + (9.*(ux + uy))/porosity) + FyB*(3. - (3.*uy)/porosity + (9.*(ux + uy))/porosity) + 
      FzB*(0. - (3.*uz)/porosity));

            // q = 8
            distB[8*Np+n] = f8B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-ux - uy) + (4.5*(-ux - uy)*(-ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(-3. - (3.*ux)/porosity - (9.*(-ux - uy))/porosity) + FyB*(-3. - (9.*(-ux - uy))/porosity - (3.*uy)/porosity) + 
      FzB*(0. - (3.*uz)/porosity));

            // q = 9
            distB[9*Np+n] = f9B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(ux - uy) + (4.5*(ux - uy)*(ux - uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(3. - (3.*ux)/porosity + (9.*(ux - uy))/porosity) + FyB*(-3. - (9.*(ux - uy))/porosity - (3.*uy)/porosity) + 
      FzB*(0. - (3.*uz)/porosity));

            // q = 10
            distB[10*Np+n] = f10B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-ux + uy) + (4.5*(-ux + uy)*(-ux + uy))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(-3. - (3.*ux)/porosity - (9.*(-ux + uy))/porosity) + FyB*(3. - (3.*uy)/porosity + (9.*(-ux + uy))/porosity) + 
      FzB*(0. - (3.*uz)/porosity));

            // q = 11
            distB[11*Np+n] = f11B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(ux + uz) + (4.5*(ux + uz)*(ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyB*(0. - (3.*uy)/porosity) + FxB*(3. - (3.*ux)/porosity + (9.*(ux + uz))/porosity) + 
      FzB*(3. - (3.*uz)/porosity + (9.*(ux + uz))/porosity));

            // q = 12
            distB[12*Np+n] = f12B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-ux - uz) + (4.5*(-ux - uz)*(-ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyB*(0. - (3.*uy)/porosity) + FxB*(-3. - (3.*ux)/porosity - (9.*(-ux - uz))/porosity) + 
      FzB*(-3. - (9.*(-ux - uz))/porosity - (3.*uz)/porosity));

            // q = 13
            distB[13*Np+n] = f13B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(ux - uz) + (4.5*(ux - uz)*(ux - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyB*(0. - (3.*uy)/porosity) + FxB*(3. - (3.*ux)/porosity + (9.*(ux - uz))/porosity) + 
      FzB*(-3. - (9.*(ux - uz))/porosity - (3.*uz)/porosity));

            // q= 14
            distB[14*Np+n] = f14B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-ux + uz) + (4.5*(-ux + uz)*(-ux + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FyB*(0. - (3.*uy)/porosity) + FxB*(-3. - (3.*ux)/porosity - (9.*(-ux + uz))/porosity) + 
      FzB*(3. - (3.*uz)/porosity + (9.*(-ux + uz))/porosity));

            // q = 15
            distB[15*Np+n] = f15B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(uy + uz) + (4.5*(uy + uz)*(uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(3. - (3.*uy)/porosity + (9.*(uy + uz))/porosity) + 
      FzB*(3. - (3.*uz)/porosity + (9.*(uy + uz))/porosity));

            // q = 16
            distB[16*Np+n] = f16B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-uy - uz) + (4.5*(-uy - uz)*(-uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(-3. - (3.*uy)/porosity - (9.*(-uy - uz))/porosity) + 
      FzB*(-3. - (9.*(-uy - uz))/porosity - (3.*uz)/porosity));

            // q = 17
            distB[17*Np+n] = f17B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(uy - uz) + (4.5*(uy - uz)*(uy - uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity) 
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(3. - (3.*uy)/porosity + (9.*(uy - uz))/porosity) + 
      FzB*(-3. - (9.*(uy - uz))/porosity - (3.*uz)/porosity));

            // q = 18
            distB[18*Np+n] = f18B*(1.0-rlx) + rlx*0.027777777777777776*rhoB*(1 + 3.*(-uy + uz) + (4.5*(-uy + uz)*(-uy + uz))/porosity - (1.5*(ux*ux + uy*uy + uz*uz))/porosity)
                    +0.027777777777777776*(1. - 0.5*rlx)*(FxB*(0. - (3.*ux)/porosity) + FyB*(-3. - (3.*uy)/porosity - (9.*(-uy + uz))/porosity) + 
      FzB*(3. - (3.*uz)/porosity + (9.*(-uy + uz))/porosity));


            //Update velocity on device
            Velocity[0*Np+n] = ux;
            Velocity[1*Np+n] = uy;
            Velocity[2*Np+n] = uz;
            //Update pressure on device
            Pressure[n] = (rhoA+rhoB+Gsc*rhoA*rhoB)/3.0;
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_GreyscaleSC_Init(int *Map, double *distA, double *distB, double *DenA, double *DenB, int Np)
{
	int n;
    int ijk;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		if (n<Np ){
            ijk = Map[n];

			distA[0*Np+n]  = DenA[ijk]*0.3333333333333333;
			distA[1*Np+n]  = DenA[ijk]*0.055555555555555555;		//double(100*n)+1.f;
			distA[2*Np+n]  = DenA[ijk]*0.055555555555555555;	//double(100*n)+2.f;
			distA[3*Np+n]  = DenA[ijk]*0.055555555555555555;	//double(100*n)+3.f;
			distA[4*Np+n]  = DenA[ijk]*0.055555555555555555;	//double(100*n)+4.f;
			distA[5*Np+n]  = DenA[ijk]*0.055555555555555555;	//double(100*n)+5.f;
			distA[6*Np+n]  = DenA[ijk]*0.055555555555555555;	//double(100*n)+6.f;
			distA[7*Np+n]  = DenA[ijk]*0.0277777777777778;   //double(100*n)+7.f;
			distA[8*Np+n]  = DenA[ijk]*0.0277777777777778;   //double(100*n)+8.f;
			distA[9*Np+n]  = DenA[ijk]*0.0277777777777778;   //double(100*n)+9.f;
			distA[10*Np+n] = DenA[ijk]*0.0277777777777778;  //double(100*n)+10.f;
			distA[11*Np+n] = DenA[ijk]*0.0277777777777778;  //double(100*n)+11.f;
			distA[12*Np+n] = DenA[ijk]*0.0277777777777778;  //double(100*n)+12.f;
			distA[13*Np+n] = DenA[ijk]*0.0277777777777778;  //double(100*n)+13.f;
			distA[14*Np+n] = DenA[ijk]*0.0277777777777778;  //double(100*n)+14.f;
			distA[15*Np+n] = DenA[ijk]*0.0277777777777778;  //double(100*n)+15.f;
			distA[16*Np+n] = DenA[ijk]*0.0277777777777778;  //double(100*n)+16.f;
			distA[17*Np+n] = DenA[ijk]*0.0277777777777778;  //double(100*n)+17.f;
			distA[18*Np+n] = DenA[ijk]*0.0277777777777778;  //double(100*n)+18.f;

			distB[0*Np+n]  = DenB[ijk]*0.3333333333333333;
			distB[1*Np+n]  = DenB[ijk]*0.055555555555555555;		//double(100*n)+1.f;
			distB[2*Np+n]  = DenB[ijk]*0.055555555555555555;	//double(100*n)+2.f;
			distB[3*Np+n]  = DenB[ijk]*0.055555555555555555;	//double(100*n)+3.f;
			distB[4*Np+n]  = DenB[ijk]*0.055555555555555555;	//double(100*n)+4.f;
			distB[5*Np+n]  = DenB[ijk]*0.055555555555555555;	//double(100*n)+5.f;
			distB[6*Np+n]  = DenB[ijk]*0.055555555555555555;	//double(100*n)+6.f;
			distB[7*Np+n]  = DenB[ijk]*0.0277777777777778;   //double(100*n)+7.f;
			distB[8*Np+n]  = DenB[ijk]*0.0277777777777778;   //double(100*n)+8.f;
			distB[9*Np+n]  = DenB[ijk]*0.0277777777777778;   //double(100*n)+9.f;
			distB[10*Np+n] = DenB[ijk]*0.0277777777777778;  //double(100*n)+10.f;
			distB[11*Np+n] = DenB[ijk]*0.0277777777777778;  //double(100*n)+11.f;
			distB[12*Np+n] = DenB[ijk]*0.0277777777777778;  //double(100*n)+12.f;
			distB[13*Np+n] = DenB[ijk]*0.0277777777777778;  //double(100*n)+13.f;
			distB[14*Np+n] = DenB[ijk]*0.0277777777777778;  //double(100*n)+14.f;
			distB[15*Np+n] = DenB[ijk]*0.0277777777777778;  //double(100*n)+15.f;
			distB[16*Np+n] = DenB[ijk]*0.0277777777777778;  //double(100*n)+16.f;
			distB[17*Np+n] = DenB[ijk]*0.0277777777777778;  //double(100*n)+17.f;
			distB[18*Np+n] = DenB[ijk]*0.0277777777777778;  //double(100*n)+18.f;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAodd_GreyscaleSC_Density(int *neighborList, int *Map, double *distA, double *distB, double *DenA, double *DenB, int start, int finish, int Np){
	int n,nread;
	double fq,nA,nB;
    int idx;

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
            idx = Map[n];
			DenA[idx] = nA;
			DenB[idx] = nB;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAeven_GreyscaleSC_Density(int *Map, double *distA, double *distB, double *DenA, double *DenB, int start, int finish, int Np){
	int n,nread;
	double fq,nA,nB;
    int idx;

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
            idx = Map[n];
			DenA[idx] = nA;
			DenB[idx] = nB;
		}
	}
}

__global__ void dvc_ScaLBL_D3Q19_GreyscaleSC_Gradient(int *neighborList, int *Map, double *Den, double *DenGrad, int strideY, int strideZ,int start, int finish, int Np){

	int n,nn;
    int ijk;
	// distributions
	double m1,m2,m3,m4,m5,m6,m7,m8,m9;
	double m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double nx,ny,nz;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
		n =  S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x + start;
		if (n<finish) {

			// Get the 1D index based on regular data layout
			ijk = Map[n];
			//					COMPUTE THE COLOR GRADIENT
			//........................................................................
			//.................Read Phase Indicator Values............................
			//........................................................................
			nn = ijk-1;							// neighbor index (get convention)
			m1 = Den[nn];						// get neighbor for phi - 1
			//........................................................................
			nn = ijk+1;							// neighbor index (get convention)
			m2 = Den[nn];						// get neighbor for phi - 2
			//........................................................................
			nn = ijk-strideY;							// neighbor index (get convention)
			m3 = Den[nn];					// get neighbor for phi - 3
			//........................................................................
			nn = ijk+strideY;							// neighbor index (get convention)
			m4 = Den[nn];					// get neighbor for phi - 4
			//........................................................................
			nn = ijk-strideZ;						// neighbor index (get convention)
			m5 = Den[nn];					// get neighbor for phi - 5
			//........................................................................
			nn = ijk+strideZ;						// neighbor index (get convention)
			m6 = Den[nn];					// get neighbor for phi - 6
			//........................................................................
			nn = ijk-strideY-1;						// neighbor index (get convention)
			m7 = Den[nn];					// get neighbor for phi - 7
			//........................................................................
			nn = ijk+strideY+1;						// neighbor index (get convention)
			m8 = Den[nn];					// get neighbor for phi - 8
			//........................................................................
			nn = ijk+strideY-1;						// neighbor index (get convention)
			m9 = Den[nn];					// get neighbor for phi - 9
			//........................................................................
			nn = ijk-strideY+1;						// neighbor index (get convention)
			m10 = Den[nn];					// get neighbor for phi - 10
			//........................................................................
			nn = ijk-strideZ-1;						// neighbor index (get convention)
			m11 = Den[nn];					// get neighbor for phi - 11
			//........................................................................
			nn = ijk+strideZ+1;						// neighbor index (get convention)
			m12 = Den[nn];					// get neighbor for phi - 12
			//........................................................................
			nn = ijk+strideZ-1;						// neighbor index (get convention)
			m13 = Den[nn];					// get neighbor for phi - 13
			//........................................................................
			nn = ijk-strideZ+1;						// neighbor index (get convention)
			m14 = Den[nn];					// get neighbor for phi - 14
			//........................................................................
			nn = ijk-strideZ-strideY;					// neighbor index (get convention)
			m15 = Den[nn];					// get neighbor for phi - 15
			//........................................................................
			nn = ijk+strideZ+strideY;					// neighbor index (get convention)
			m16 = Den[nn];					// get neighbor for phi - 16
			//........................................................................
			nn = ijk+strideZ-strideY;					// neighbor index (get convention)
			m17 = Den[nn];					// get neighbor for phi - 17
			//........................................................................
			nn = ijk-strideZ+strideY;					// neighbor index (get convention)
			m18 = Den[nn];					// get neighbor for phi - 18
			//............Compute the Color Gradient...................................
			nx = -1.f/18.f*(m1-m2+0.5*(m7-m8+m9-m10+m11-m12+m13-m14));
			ny = -1.f/18.f*(m3-m4+0.5*(m7-m8-m9+m10+m15-m16+m17-m18));
			nz = -1.f/18.f*(m5-m6+0.5*(m11-m12-m13+m14+m15-m16-m17+m18));
			
			DenGrad[n] = nx;
			DenGrad[Np+n] = ny;
			DenGrad[2*Np+n] = nz;
		}
	}
}

extern "C" void ScaLBL_D3Q19_GreyscaleSC_Init(int *Map, double *distA,double *distB, double *DenA, double *DenB, int Np){
	dvc_ScaLBL_D3Q19_GreyscaleSC_Init<<<NBLOCKS,NTHREADS >>>(Map,distA,distB,DenA,DenB,Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_GreyscaleSC_Init: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleSC_Density(int *NeighborList, int *Map, double *distA, double *distB, double *DenA, double *DenB, int start, int finish, int Np){

	dvc_ScaLBL_D3Q19_AAodd_GreyscaleSC_Density<<<NBLOCKS,NTHREADS >>>(NeighborList, Map, distA, distB, DenA, DenB, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_GreyscaleSC_Density: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleSC_Density(int *Map, double *distA, double *distB, double *DenA, double *DenB, int start, int finish, int Np){

	dvc_ScaLBL_D3Q19_AAeven_GreyscaleSC_Density<<<NBLOCKS,NTHREADS >>>(Map,distA, distB, DenA, DenB, start, finish, Np);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_GreyscaleSC_Density: %s \n",cudaGetErrorString(err));
	}
}


extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleSC_MRT(int *neighborList, int *Map, double *distA, double *distB, double *DenA,double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

    dvc_ScaLBL_D3Q19_AAodd_GreyscaleSC_MRT<<<NBLOCKS,NTHREADS >>>(neighborList,Map,distA,distB,DenA,DenB,DenGradA,DenGradB,SolidForceA,SolidForceB,Poros,Perm,Velocity,Pressure, 
                tauA,tauB,tauA_eff,tauB_eff,Gsc,Gx,Gy,Gz,start,finish,Np);

    cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_GreyscaleSC_MRT: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleSC_MRT(int *Map,double *distA, double *distB, double *DenA,double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

    dvc_ScaLBL_D3Q19_AAeven_GreyscaleSC_MRT<<<NBLOCKS,NTHREADS >>>(Map,distA,distB,DenA,DenB,DenGradA,DenGradB,SolidForceA,SolidForceB,Poros,Perm,Velocity,Pressure, 
                tauA,tauB,tauA_eff,tauB_eff,Gsc,Gx,Gy,Gz,start,finish,Np);

    cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_GreyscaleSC_MRT: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_GreyscaleSC_BGK(int *neighborList, int *Map, double *distA, double *distB, double *DenA, double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

    dvc_ScaLBL_D3Q19_AAodd_GreyscaleSC_BGK<<<NBLOCKS,NTHREADS >>>(neighborList,Map,distA,distB,DenA,DenB,DenGradA,DenGradB,SolidForceA,SolidForceB,Poros,Perm,Velocity,Pressure, 
                tauA,tauB,tauA_eff,tauB_eff,Gsc,Gx,Gy,Gz,start,finish,Np);

    cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAodd_GreyscaleSC_BGK: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_GreyscaleSC_BGK(int *Map,double *distA, double *distB, double *DenA, double *DenB, double *DenGradA, double *DenGradB, 
                double *SolidForceA, double *SolidForceB, double *Poros,double *Perm, double *Velocity,double *Pressure, 
                double tauA,double tauB,double tauA_eff,double tauB_eff, double Gsc, double Gx, double Gy, double Gz,                                                 
                int start, int finish, int Np){

    dvc_ScaLBL_D3Q19_AAeven_GreyscaleSC_BGK<<<NBLOCKS,NTHREADS >>>(Map,distA,distB,DenA,DenB,DenGradA,DenGradB,SolidForceA,SolidForceB,Poros,Perm,Velocity,Pressure, 
                tauA,tauB,tauA_eff,tauB_eff,Gsc,Gx,Gy,Gz,start,finish,Np);

    cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_AAeven_GreyscaleSC_BGK: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_D3Q19_GreyscaleSC_Gradient(int *neighborList, int *Map, double *Den, double *DenGrad, int strideY, int strideZ,int start, int finish, int Np){

	dvc_ScaLBL_D3Q19_GreyscaleSC_Gradient<<<NBLOCKS,NTHREADS >>>(neighborList, Map, Den, DenGrad, strideY, strideZ,start, finish, Np);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("CUDA error in ScaLBL_D3Q19_GreyscaleSC_Gradient: %s \n",cudaGetErrorString(err));
	}
}
