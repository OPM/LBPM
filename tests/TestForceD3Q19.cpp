#include <iostream>
#include "common/MPI_Helpers.h"
#include "common/Utilities.h"
#include <math.h>


double mrt_V1=0.05263157894736842;
double mrt_V2=0.012531328320802;
double mrt_V3=0.04761904761904762;
double mrt_V4=0.004594820384294068;
double mrt_V5=0.01587301587301587;
double mrt_V6=0.0555555555555555555555555;
double mrt_V7=0.02777777777777778;
double mrt_V8=0.08333333333333333;
double mrt_V9=0.003341687552213868;
double mrt_V10=0.003968253968253968;
double mrt_V11=0.01388888888888889;
double mrt_V12=0.04166666666666666;

/*
# Rcode to check the moments

f=c(1.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18)

rho=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
jx=c(0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0)
jy=c(0,0,0,1,-1,0,0,1,-1,-1,1,0,0,0,0,1,-1,1,-1)
jz=c(0,0,0,0,0,1,-1,0,0,0,0,1,-1,-1,1,1,-1,-1,1)

M1=c(-30,-11,-11,-11,-11,-11,-11,8,8,8,8,8,8,8,8,8,8,8,8)
M2=c(12,-4,-4,-4,-4,-4,-4,1,1,1,1,1,1,1,1,1,1,1,1)
M4=c(0,-4,4,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0)
M6=c(0,0,0,-4,4,0,0,1,-1,-1,1,0,0,0,0,1,-1,1,-1)
M8=c(0,0,0,0,0,-4,4,0,0,0,0,1,-1,-1,1,1,-1,-1,1)
M9=c(0,2,2,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-2,-2,-2,-2)
M10=c(0,-4,-4,2,2,2,2,1,1,1,1,1,1,1,1,-2,-2,-2,-2)
M11=c(0,0,0,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,0,0,0,0)
M12=c(0,0,0,-2,-2,2,2,1,1,1,1,-1,-1,-1,-1,0,0,0,0)
M13=c(0,0,0,0,0,0,0,1,1,-1,-1,0,0,0,0,0,0,0,0)
M14=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,-1,-1)
M15=c(0,0,0,0,0,0,0,0,0,0,0,1,1,-1,-1,0,0,0,0)
M16=c(0,0,0,0,0,0,0,1,-1,1,-1,-1,1,-1,1,0,0,0,0)
M17=c(0,0,0,0,0,0,0,-1,1,1,-1,0,0,0,0,1,-1,1,-1)
M18=c(0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,-1,1,1,-1)
 */
inline void MRT_Transform(double *dist, int Np, double Fx, double Fy, double Fz) {

	double fq,fp;
	// conserved momemnts
	double rho,jx,jy,jz;

	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	for (int n=0; n<Np; n++){

		//........................................................................
		//					READ THE DISTRIBUTIONS
		//		(read from opposite array due to previous swap operation)
		//........................................................................
		fq = dist[n];
		//printf("q=0: %f\n",fq);
		rho = fq;
		m1  = -30.0*fq;
		m2  = 12.0*fq;

		// q=1
		fp = dist[10*Np+n];
		//printf("q=1: %f\n",fp);
		rho += fp;
		m1 -= 11.0*fp;
		m2 -= 4.0*fp;
		jx = fp;
		m4 = -4.0*fp;
		m9 = 2.0*fp;
		m10 = -4.0*fp;

		// f2 = dist[10*Np+n];
		fq = dist[Np+n];
		//printf("q=2: %f\n",fq);
		rho += fq;
		m1 -= 11.0*(fq);
		m2 -= 4.0*(fq);
		jx -= fq;
		m4 += 4.0*(fq);
		m9 += 2.0*(fq);
		m10 -= 4.0*(fq);

		// q=3
		fq = dist[11*Np+n];
		//printf("q=3: %f\n",fq);
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
		fq = dist[2*Np+n];
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
		fq = dist[12*Np+n];
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
		fq = dist[3*Np+n];
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
		fq = dist[13*Np+n];
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
		fq = dist[4*Np+n];
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
		fq = dist[14*Np+n];
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
		fq = dist[5*Np+n];
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
		fq = dist[15*Np+n];
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
		fq = dist[6*Np+n];
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
		fq = dist[16*Np+n];
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
		fq = dist[7*Np+n];
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
		fq = dist[17*Np+n];
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
		fq = dist[8*Np+n];
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

		// q=18
		fq = dist[9*Np+n];
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

		//.................inverse transformation......................................................
		// q=0
		fq = mrt_V1*rho-mrt_V2*m1+mrt_V3*m2;
		dist[n] = fq;

		// q = 1
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jx-m4)+mrt_V6*(m9-m10)+0.16666666*Fx;
		dist[10*Np+n] = fq;

		// q=2
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m4-jx)+mrt_V6*(m9-m10) -  0.16666666*Fx;
		dist[Np+n] = fq;

		// q = 3
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jy-m6)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) + 0.16666666*Fy;
		dist[11*Np+n] = fq;

		// q = 4
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m6-jy)+mrt_V7*(m10-m9)+mrt_V8*(m11-m12) - 0.16666666*Fy;
		dist[2*Np+n] = fq;

		// q = 5
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(jz-m8)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) + 0.16666666*Fz;
		dist[12*Np+n] = fq;

		// q = 6
		fq = mrt_V1*rho-mrt_V4*m1-mrt_V5*m2+0.1*(m8-jz)+mrt_V7*(m10-m9)+mrt_V8*(m12-m11) - 0.16666666*Fz;
		dist[3*Np+n] = fq;

		// q = 7
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx+jy)+0.025*(m4+m6)
                                                																														+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                                                																														+mrt_V12*m12+0.25*m13+0.125*(m16-m17) + 0.08333333333*(Fx+Fy);
		dist[13*Np+n] = fq;


		// q = 8
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jy)-0.025*(m4+m6) +mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
				+mrt_V12*m12+0.25*m13+0.125*(m17-m16) - 0.08333333333*(Fx+Fy);
		dist[4*Np+n] = fq;

		// q = 9
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jx-jy)+0.025*(m4-m6)
                                                																														+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                                                																														+mrt_V12*m12-0.25*m13+0.125*(m16+m17) + 0.08333333333*(Fx-Fy);
		dist[14*Np+n] = fq;

		// q = 10
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2+0.1*(jy-jx)+0.025*(m6-m4)
                                                																														+mrt_V7*m9+mrt_V11*m10+mrt_V8*m11
                                                																														+mrt_V12*m12-0.25*m13-0.125*(m16+m17)- 0.08333333333*(Fx-Fy);
		dist[5*Np+n] = fq;

		// q = 11
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jx+jz)+0.025*(m4+m8)
				+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
				-mrt_V12*m12+0.25*m15+0.125*(m18-m16) + 0.08333333333*(Fx+Fz);
		dist[15*Np+n] = fq;

		// q = 12
		fq = mrt_V1*rho+mrt_V9*m1+mrt_V10*m2-0.1*(jx+jz)-0.025*(m4+m8)
                                        																														+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
                                        																														-mrt_V12*m12+0.25*m15+0.125*(m16-m18) - 0.08333333333*(Fx+Fz);
		dist[6*Np+n] = fq;

		// q = 13
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jx-jz)+0.025*(m4-m8)
				+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
				-mrt_V12*m12-0.25*m15-0.125*(m16+m18) + 0.08333333333*(Fx-Fz);
		dist[16*Np+n] = fq;

		// q= 14
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jz-jx)+0.025*(m8-m4)
				+mrt_V7*m9+mrt_V11*m10-mrt_V8*m11
				-mrt_V12*m12-0.25*m15+0.125*(m16+m18) - 0.08333333333*(Fx-Fz);
		dist[7*Np+n] = fq;

		// q = 15
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jy+jz)+0.025*(m6+m8)
				-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m17-m18) + 0.08333333333*(Fy+Fz);
		dist[17*Np+n] = fq;

		// q = 16
		fq =  mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2-0.1*(jy+jz)-0.025*(m6+m8)
				-mrt_V6*m9-mrt_V7*m10+0.25*m14+0.125*(m18-m17)- 0.08333333333*(Fy+Fz);
		dist[8*Np+n] = fq;

		// q = 17
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jy-jz)+0.025*(m6-m8)
				-mrt_V6*m9-mrt_V7*m10-0.25*m14+0.125*(m17+m18) + 0.08333333333*(Fy-Fz);
		dist[18*Np+n] = fq;

		// q = 18
		fq = mrt_V1*rho+mrt_V9*m1
				+mrt_V10*m2+0.1*(jz-jy)+0.025*(m8-m6)
				-mrt_V6*m9-mrt_V7*m10-0.25*m14-0.125*(m17+m18) - 0.08333333333*(Fy-Fz);
		dist[9*Np+n] = fq;

	}
}


int main (int argc, char **argv)
{
	MPI_Init(&argc,&argv);
	int rank = MPI_WORLD_RANK();
	int nprocs = MPI_WORLD_SIZE();

	for (int i=0; i<nprocs; i++) {
		if ( rank==i )
			printf("%i of %i: Testing force term \n",rank,nprocs);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	// Create a memory leak for valgrind to find
	if ( nprocs==1 ) {
		double *x = new double[1];
		ASSERT(x!=NULL);
	}

    // set the error code
    // Note: the error code should be consistent across all processors
    int error = 0;
    
    int Np = 1;
    int Q = 9;

    double Fx = 1.0;
    double Fy = 1.0;
    double Fz = 1.0;
    
    double *dist;    
    double * Velocity;
    
    dist = new double [19*Np];
    Velocity = new double [3*Np];


    for (int n=0; n<Np; n++){
    	dist[n] = 0.3333333333333333;
    	dist[10*Np+n] = 0.055555555555555555;                //double(100*n)+1.f;
    	dist[Np+n] = 0.055555555555555555;    //double(100*n)+2.f;
    	dist[11*Np+n] = 0.055555555555555555;     //double(100*n)+3.f;
    	dist[2*Np+n] = 0.055555555555555555;  //double(100*n)+4.f;
    	dist[12*Np+n] = 0.055555555555555555;   //double(100*n)+5.f;
    	dist[3*Np+n] = 0.055555555555555555;  //double(100*n)+6.f;
    	dist[13*Np+n] = 0.0277777777777778;   //double(100*n)+7.f;
    	dist[4*Np+n] = 0.0277777777777778;   //double(100*n)+8.f;
    	dist[14*Np+n] = 0.0277777777777778;   //double(100*n)+9.f;
    	dist[5*Np+n] = 0.0277777777777778;  //double(100*n)+10.f;
    	dist[15*Np+n] = 0.0277777777777778;  //double(100*n)+11.f;
    	dist[6*Np+n] = 0.0277777777777778;  //double(100*n)+12.f;
    	dist[16*Np+n] = 0.0277777777777778;  //double(100*n)+13.f;
    	dist[7*Np+n] = 0.0277777777777778;  //double(100*n)+14.f;
    	dist[17*Np+n] = 0.0277777777777778;  //double(100*n)+15.f;
    	dist[8*Np+n] = 0.0277777777777778;  //double(100*n)+16.f;
    	dist[18*Np+n] = 0.0277777777777778;  //double(100*n)+17.f;
    	dist[9*Np+n] = 0.0277777777777778;  //double(100*n)+18.f;
    }

    MRT_Transform(dist,Np,Fx,Fy,Fz);

    double *vel;
    vel= new double [3*Np];
    // distributions
    double f1,f2,f3,f4,f5,f6,f7,f8,f9;
    double f10,f11,f12,f13,f14,f15,f16,f17,f18;
    double vx,vy,vz;

    for (int n=0; n<Np; n++){
    	//........................................................................
    	// Registers to store the distributions
    	//........................................................................
    	f2 = dist[Np+n];
    	f4 = dist[2*Np+n];
    	f6 = dist[3*Np+n];
    	f8 = dist[4*Np+n];
    	f10 = dist[5*Np+n];
    	f12 = dist[6*Np+n];
    	f14 = dist[7*Np+n];
    	f16 = dist[8*Np+n];
    	f18 = dist[9*Np+n];
    	//........................................................................
		f1 = dist[10*Np+n];
    	f3 = dist[11*Np+n];
    	f5 = dist[12*Np+n];
    	f7 = dist[13*Np+n];
    	f9 = dist[14*Np+n];
    	f11 = dist[15*Np+n];
    	f13 = dist[16*Np+n];
    	f15 = dist[17*Np+n];
    	f17 = dist[18*Np+n];
    	//.................Compute the velocity...................................
		vx = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
		vy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		vz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
		//..................Write the velocity.....................................
		vel[n] = vx;
		vel[Np+n] = vy;
		vel[2*Np+n] = vz;
		printf("vx=%f, vy=%f, vz=%f \n",vx,vy,vz);
		//........................................................................
    }
    

    printf("Fx = %f; Computed vx=%f \n",Fx,vel[0]);
    printf("Fy = %f; Computed vy=%f \n",Fy,vel[Np+0]);
    printf("Fz = %f; Computed vz=%f \n",Fz,vel[2*Np+0]);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
