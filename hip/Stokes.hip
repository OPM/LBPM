#include <stdio.h>
#include <math.h>
#include "hip/hip_runtime.h"

#define NBLOCKS 1024
#define NTHREADS 256


__global__  void dvc_ScaLBL_D3Q19_AAodd_StokesMRT(int *neighborList, double *dist, double *Velocity, double *ChargeDensity, double *ElectricField, double rlx_setA, double rlx_setB, double Gx, double Gy, double Gz, double rho0, double den_scale, double h, double time_conv,bool UseSlippingVelBC,int start, int finish, int Np){

    int n;
    double fq;
	// conserved momemnts
	double rho,jx,jy,jz;
    double ux,uy,uz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	int nread;
    // body force due to electric field
    double rhoE;//charge density
    double Ex,Ey,Ez;
    // total body force
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

            //Load data
            rhoE = ChargeDensity[n];
            Ex = ElectricField[n+0*Np];
            Ey = ElectricField[n+1*Np];
            Ez = ElectricField[n+2*Np];
            //compute total body force, including input body force (Gx,Gy,Gz)
            Fx = (UseSlippingVelBC==1) ? Gx : Gx + rhoE * Ex * (time_conv * time_conv) / (h * h * 1.0e-12) /
                          den_scale; //the extra factors at the end necessarily convert unit from phys to LB
            Fy = (UseSlippingVelBC==1) ? Gy : Gy + rhoE * Ey * (time_conv * time_conv) / (h * h * 1.0e-12) /
                          den_scale;
            Fz = (UseSlippingVelBC==1) ? Gz : Gz + rhoE * Ez * (time_conv * time_conv) / (h * h * 1.0e-12) /
                          den_scale;

            // q=0
            fq = dist[n];
            rho = fq;
            m1  = -30.0*fq;
            m2  = 12.0*fq;

            // q=1
            nread = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
            fq = dist[nread]; // reading the f1 data into register fq
            //fp = dist[10*Np+n];
            rho += fq;
            m1 -= 11.0*fq;
            m2 -= 4.0*fq;
            jx = fq;
            m4 = -4.0*fq;
            m9 = 2.0*fq;
            m10 = -4.0*fq;

            // f2 = dist[10*Np+n];
            nread = neighborList[n+Np]; // neighbor 1 ( < 10Np => even part of dist)
            fq = dist[nread];  // reading the f2 data into register fq
            //fq = dist[Np+n];
            rho += fq;
            m1 -= 11.0*(fq);
            m2 -= 4.0*(fq);
            jx -= fq;
            m4 += 4.0*(fq);
            m9 += 2.0*(fq);
            m10 -= 4.0*(fq);

            // q=3
            nread = neighborList[n+2*Np]; // neighbor 4
            fq = dist[nread];
            //fq = dist[11*Np+n];
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
            nread = neighborList[n+3*Np]; // neighbor 3
            fq = dist[nread];
            //fq = dist[2*Np+n];
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
            nread = neighborList[n+4*Np];
            fq = dist[nread];
            //fq = dist[12*Np+n];
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
            nread = neighborList[n+5*Np];
            fq = dist[nread];
            //fq = dist[3*Np+n];
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
            nread = neighborList[n+6*Np];
            fq = dist[nread];
            //fq = dist[13*Np+n];
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
            nread = neighborList[n+7*Np];
            fq = dist[nread];
            //fq = dist[4*Np+n];
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
            nread = neighborList[n+8*Np];
            fq = dist[nread];
            //fq = dist[14*Np+n];
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
            nread = neighborList[n+9*Np];
            fq = dist[nread];
            //fq = dist[5*Np+n];
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
            nread = neighborList[n+10*Np];
            fq = dist[nread];
            //fq = dist[15*Np+n];
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
            nread = neighborList[n+11*Np];
            fq = dist[nread];
            //fq = dist[6*Np+n];
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
            nread = neighborList[n+12*Np];
            fq = dist[nread];
            //fq = dist[16*Np+n];
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
            nread = neighborList[n+13*Np];
            fq = dist[nread];
            //fq = dist[7*Np+n];
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

            // write the velocity 
            ux = jx / rho0;
            uy = jy / rho0;
            uz = jz / rho0;
            Velocity[n] = ux;
            Velocity[Np+n] = uy;
            Velocity[2*Np+n] = uz;

            //..............incorporate external force................................................
            //..............carry out relaxation process...............................................
            m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho0 - 11*rho) - m1);
            m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho0) - m2);
            m4 = m4 + rlx_setB*((-0.6666666666666666*jx) - m4);
            m6 = m6 + rlx_setB*((-0.6666666666666666*jy) - m6);
            m8 = m8 + rlx_setB*((-0.6666666666666666*jz) - m8);
            m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho0) - m9);
            m10 = m10 + rlx_setA*(-0.5*((2*jx*jx-jy*jy-jz*jz)/rho) - m10);
            m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho0) - m11);
            m12 = m12 + rlx_setA*(-0.5*((jy*jy-jz*jz)/rho0) - m12);
            m13 = m13 + rlx_setA*((jx*jy/rho0) - m13);
            m14 = m14 + rlx_setA*((jy*jz/rho0) - m14);
            m15 = m15 + rlx_setA*((jx*jz/rho0) - m15);
            m16 = m16 + rlx_setB*( - m16);
            m17 = m17 + rlx_setB*( - m17);
            m18 = m18 + rlx_setB*( - m18);
            //.......................................................................................................
            //.................inverse transformation......................................................

            // q=0
            fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
            dist[n] = fq;

            // q = 1
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10)+0.16666666*Fx;
            nread = neighborList[n+Np];
            dist[nread] = fq;

            // q=2
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10) -  0.16666666*Fx;
            nread = neighborList[n];
            dist[nread] = fq;

            // q = 3
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) + 0.16666666*Fy;
            nread = neighborList[n+3*Np];
            dist[nread] = fq;

            // q = 4
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) - 0.16666666*Fy;
            nread = neighborList[n+2*Np];
            dist[nread] = fq;

            // q = 5
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) + 0.16666666*Fz;
            nread = neighborList[n+5*Np];
            dist[nread] = fq;

            // q = 6
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) - 0.16666666*Fz;
            nread = neighborList[n+4*Np];
            dist[nread] = fq;

            // q = 7
            fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)
                                                                                    +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                                                                                    +mrt_V12*m12+0.25*m13+0.125*(m16-m17) + 0.08333333333*(Fx+Fy);
            nread = neighborList[n+7*Np];
            dist[nread] = fq;

            // q = 8
            fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                    +mrt_V12*m12+0.25*m13+0.125*(m17-m16) - 0.08333333333*(Fx+Fy);
            nread = neighborList[n+6*Np];
            dist[nread] = fq;

            // q = 9
            fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)
                                                                                    +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                                                                                    +mrt_V12*m12-0.25*m13+0.125*(m16+m17) + 0.08333333333*(Fx-Fy);
            nread = neighborList[n+9*Np];
            dist[nread] = fq;

            // q = 10
            fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)
                                                                                    +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                                                                                    +mrt_V12*m12-0.25*m13-0.125*(m16+m17)- 0.08333333333*(Fx-Fy);
            nread = neighborList[n+8*Np];
            dist[nread] = fq;

            // q = 11
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
                    +mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
                    -mrt_V12*m12+0.25*m15+0.125*(m18-m16) + 0.08333333333*(Fx+Fz);
            nread = neighborList[n+11*Np];
            dist[nread] = fq;

            // q = 12
            fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)
                                                                            +mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
                                                                            -mrt_V12*m12+0.25*m15+0.125*(m16-m18) - 0.08333333333*(Fx+Fz);
            nread = neighborList[n+10*Np];
            dist[nread]= fq;

            // q = 13
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
                    +mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
                    -mrt_V12*m12-0.25*m15-0.125*(m16+m18) + 0.08333333333*(Fx-Fz);
            nread = neighborList[n+13*Np];
            dist[nread] = fq;

            // q= 14
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
                    +mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
                    -mrt_V12*m12-0.25*m15+0.125*(m16+m18) - 0.08333333333*(Fx-Fz);
            nread = neighborList[n+12*Np];
            dist[nread] = fq;


            // q = 15
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
                    -mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18) + 0.08333333333*(Fy+Fz);
            nread = neighborList[n+15*Np];
            dist[nread] = fq;

            // q = 16
            fq =  mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
                    -mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17)- 0.08333333333*(Fy+Fz);
            nread = neighborList[n+14*Np];
            dist[nread] = fq;


            // q = 17
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
                    -mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18) + 0.08333333333*(Fy-Fz);
            nread = neighborList[n+17*Np];
            dist[nread] = fq;

            // q = 18
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
                    -mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18) - 0.08333333333*(Fy-Fz);
            nread = neighborList[n+16*Np];
            dist[nread] = fq;
		}
	}
}

__global__  void dvc_ScaLBL_D3Q19_AAeven_StokesMRT(double *dist, double *Velocity, double *ChargeDensity, double *ElectricField, double rlx_setA, double rlx_setB, double Gx, double Gy, double Gz,double rho0, double den_scale, double h, double time_conv, bool UseSlippingVelBC,int start, int finish, int Np){

    int n;
    double fq;
	// conserved momemnts
	double rho,jx,jy,jz;
    double ux,uy,uz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    // body force due to electric field
    double rhoE;//charge density
    double Ex,Ey,Ez;
    // total body force
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

            //Load data
            rhoE = ChargeDensity[n];
            Ex = ElectricField[n+0*Np];
            Ey = ElectricField[n+1*Np];
            Ez = ElectricField[n+2*Np];
            //compute total body force, including input body force (Gx,Gy,Gz)
            Fx = (UseSlippingVelBC==1) ? Gx : Gx + rhoE * Ex * (time_conv * time_conv) / (h * h * 1.0e-12) /
                          den_scale; //the extra factors at the end necessarily convert unit from phys to LB
            Fy = (UseSlippingVelBC==1) ? Gy : Gy + rhoE * Ey * (time_conv * time_conv) / (h * h * 1.0e-12) /
                          den_scale;
            Fz = (UseSlippingVelBC==1) ? Gz : Gz + rhoE * Ez * (time_conv * time_conv) / (h * h * 1.0e-12) /
                          den_scale;

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

            // write the velocity 
            ux = jx / rho0;
            uy = jy / rho0;
            uz = jz / rho0;
            Velocity[n] = ux;
            Velocity[Np+n] = uy;
            Velocity[2*Np+n] = uz;


            //........................................................................
            //					READ THE DISTRIBUTIONS
            //		(read from opposite array due to previous swap operation)
            //........................................................................

            //..............incorporate external force................................................
            //..............carry out relaxation process...............................................
            m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho0 - 11*rho) - m1);
            m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho0) - m2);
            m4 = m4 + rlx_setB*((-0.6666666666666666*jx) - m4);
            m6 = m6 + rlx_setB*((-0.6666666666666666*jy) - m6);
            m8 = m8 + rlx_setB*((-0.6666666666666666*jz) - m8);
            m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho0) - m9);
            m10 = m10 + rlx_setA*(-0.5*((2*jx*jx-jy*jy-jz*jz)/rho) - m10);
            m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho0) - m11);
            m12 = m12 + rlx_setA*(-0.5*((jy*jy-jz*jz)/rho0) - m12);
            m13 = m13 + rlx_setA*((jx*jy/rho0) - m13);
            m14 = m14 + rlx_setA*((jy*jz/rho0) - m14);
            m15 = m15 + rlx_setA*((jx*jz/rho0) - m15);
            m16 = m16 + rlx_setB*( - m16);
            m17 = m17 + rlx_setB*( - m17);
            m18 = m18 + rlx_setB*( - m18);
            //.......................................................................................................
            //.................inverse transformation......................................................

            // q=0
            fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
            dist[n] = fq;

            // q = 1
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10) + 0.16666666*Fx;
            dist[1*Np+n] = fq;

            // q=2
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10) -  0.16666666*Fx;
            dist[2*Np+n] = fq;

            // q = 3
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) + 0.16666666*Fy;
            dist[3*Np+n] = fq;

            // q = 4
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) - 0.16666666*Fy;
            dist[4*Np+n] = fq;

            // q = 5
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) + 0.16666666*Fz;
            dist[5*Np+n] = fq;

            // q = 6
            fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) - 0.16666666*Fz;
            dist[6*Np+n] = fq;

            // q = 7
            fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)
                                                                                    +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                                                                                    +mrt_V12*m12+0.25*m13+0.125*(m16-m17) + 0.08333333333*(Fx+Fy);
            dist[7*Np+n] = fq;


            // q = 8
            fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                    +mrt_V12*m12+0.25*m13+0.125*(m17-m16) - 0.08333333333*(Fx+Fy);
            dist[8*Np+n] = fq;

            // q = 9
            fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)
                                                                                    +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                                                                                    +mrt_V12*m12-0.25*m13+0.125*(m16+m17) + 0.08333333333*(Fx-Fy);
            dist[9*Np+n] = fq;

            // q = 10
            fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)
                                                                                    +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                                                                                    +mrt_V12*m12-0.25*m13-0.125*(m16+m17)- 0.08333333333*(Fx-Fy);
            dist[10*Np+n] = fq;


            // q = 11
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
                    +mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
                    -mrt_V12*m12+0.25*m15+0.125*(m18-m16) + 0.08333333333*(Fx+Fz);
            dist[11*Np+n] = fq;

            // q = 12
            fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)
                                                                            +mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
                                                                            -mrt_V12*m12+0.25*m15+0.125*(m16-m18) - 0.08333333333*(Fx+Fz);
            dist[12*Np+n] = fq;

            // q = 13
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
                    +mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
                    -mrt_V12*m12-0.25*m15-0.125*(m16+m18) + 0.08333333333*(Fx-Fz);
            dist[13*Np+n] = fq;

            // q= 14
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
                    +mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
                    -mrt_V12*m12-0.25*m15+0.125*(m16+m18) - 0.08333333333*(Fx-Fz);

            dist[14*Np+n] = fq;

            // q = 15
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
                    -mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18) + 0.08333333333*(Fy+Fz);
            dist[15*Np+n] = fq;

            // q = 16
            fq =  mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
                    -mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17)- 0.08333333333*(Fy+Fz);
            dist[16*Np+n] = fq;


            // q = 17
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
                    -mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18) + 0.08333333333*(Fy-Fz);
            dist[17*Np+n] = fq;

            // q = 18
            fq = mrt_V1*rho+mrt_V9*m1
                    +mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
                    -mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18) - 0.08333333333*(Fy-Fz);
            dist[18*Np+n] = fq;

            //........................................................................
		}
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_StokesMRT(int *neighborList, double *dist, double *Velocity, double *ChargeDensity, double *ElectricField, double rlx_setA, double rlx_setB, double Gx, double Gy, double Gz,double rho0, double den_scale, double h, double time_conv, bool UseSlippingVelBC, int start, int finish, int Np){
                
	//cudaProfilerStart();
	dvc_ScaLBL_D3Q19_AAodd_StokesMRT<<<NBLOCKS,NTHREADS >>>(neighborList,dist,Velocity,ChargeDensity,ElectricField,rlx_setA,rlx_setB,Gx,Gy,Gz,rho0,den_scale,h,time_conv,UseSlippingVelBC, start,finish,Np);

	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q19_AAodd_StokesMRT: %s \n",hipGetErrorString(err));
	}
	//cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q19_AAeven_StokesMRT(double *dist, double *Velocity, double *ChargeDensity, double *ElectricField, double rlx_setA, double rlx_setB, double Gx, double Gy, double Gz,double rho0, double den_scale, double h, double time_conv, bool UseSlippingVelBC, int start, int finish, int Np){
                
	//cudaProfilerStart();
	dvc_ScaLBL_D3Q19_AAeven_StokesMRT<<<NBLOCKS,NTHREADS >>>(dist,Velocity,ChargeDensity,ElectricField,rlx_setA,rlx_setB,Gx,Gy,Gz,rho0,den_scale,h,time_conv,UseSlippingVelBC, start,finish,Np);

	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("hip error in ScaLBL_D3Q19_AAeven_StokesMRT: %s \n",hipGetErrorString(err));
	}
	//cudaProfilerStop();
}

