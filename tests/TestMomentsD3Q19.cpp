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
inline void MRT_Transform(double *dist, int Np) {

	double fq,fp;
	// conserved momemnts
	double rho,jx,jy,jz;
	double Fx,Fy,Fz;
	Fx=Fy=Fz=0;
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

				printf("rho: %f\n",rho);
				printf("jx: %f\n",jx);
				printf("jy: %f\n",jy);
				printf("jz: %f\n",jz);
				printf("m1: %f\n",m1);
				printf("m2: %f\n",m2);
				printf("m4: %f\n",m4);
				printf("m6: %f\n",m6);
				printf("m8: %f\n",m8);
				printf("m9: %f\n",m9);
				printf("m10: %f\n",m10);
				printf("m11: %f\n",m11);
				printf("m12: %f\n",m12);
				printf("m13: %f\n",m13);
				printf("m14: %f\n",m14);
				printf("m15: %f\n",m15);
				printf("m16: %f\n",m16);
				printf("m17: %f\n",m17);
				printf("m18: %f\n",m18);

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
            printf("%i of %i: TestMoments\n",rank,nprocs);
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
    
    double *dist;
    dist = new double[(2*Q+1)*Np];
    
    // create the test system
    for (int n=0; n<Np; n++){
	dist[n] = 1.f*Np;
    }      
    for (int q=0; q<Q; q++){
      // set up the odd distributions
      int qodd=2*q+1;
      int qeven=2*(q+1);
      for (int n=0; n<Np; n++){
	//dist[(q+1)*Np + n] = 1.f*n + 0.01*qodd;
	dist[(q+1)*Np + n] = 1.f*n + 0.01*qeven;
      }      
      // set up the even distributions
      for (int n=0; n<Np; n++){
	//dist[(q+10)*Np + n] = 1.f*n + 0.01*qeven;
	dist[(q+10)*Np + n] = 1.f*n + 0.01*qodd;
      }      
    }

    MRT_Transform(dist,Np);

    // Check the result
    double *diff;
    diff = new double [(2*Q+1)*Np];
    
    for (int n=0; n<Np; n++){
	diff[n] = dist[n] - 1.f*Np;
    }      
    for (int q=0; q<Q; q++){
      int qodd=2*q+1;
      int qeven=2*(q+1);
      for (int n=0; n<Np; n++){
	diff[(q+1)*Np + n] = dist[(q+1)*Np + n] - (1.f*n + 0.01*qeven);
      }      
      for (int n=0; n<Np; n++){
	diff[(q+10)*Np + n] = dist[(q+10)*Np + n] - (1.f*n + 0.01*qodd);
      }      
    }

    double tol = 1e-13;
    int count=0;
    for (int idx = 0; idx<(2*Q+1)*Np; idx++){
      if(fabs(diff[idx]) > tol){
	printf("Error at %f, value=%f\n",dist[idx]-diff[idx],dist[idx]);
	count++;
      }
    }

    error=count;
    // Finished
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return error; 
}
