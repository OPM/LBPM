#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cuda.h>
#include <mpi.h>

inline void PackID(int *list, int count, char *sendbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This packs up the values that need to be sent from one processor to another
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = ID[n];
	}
}

inline void UnpackID(int *list, int count, char *recvbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This unpacks the values once they have been recieved from neighbors
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		ID[n] = recvbuf[idx];
	}
}


__global__ void PackDist(int q, int *list, int start, int count, double *sendbuf, double *dist, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
//	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[start+idx] = dist[q*N+n];
	}
}


__global__ void MapRecvDist(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
					   double *recvbuf, double *dist, int Nx, int Ny, int Nz){
	//....................................................................................
	// Unack distribution from the recv buffer
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int i,j,k,n,nn,idx;
	int N = Nx*Ny*Nz;
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx<count){
//	for (idx=0; idx<count; idx++){
		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = list[idx];
		// Get the 3-D indices
		k = n/(Nx*Ny);
		j = (n-Nx*Ny*k)/Nx;
		i = n-Nx*Ny*k-Nz*j;
		// Streaming for the non-local distribution
		i += Cqx;
		j += Cqy;
		k += Cqz;
/*		if (i < 0) i += Nx;
		if (j < 0) j += Ny;
		if (k < 0) k += Nz;
		if (!(i<Nx)) i -= Nx;
		if (!(j<Ny)) j -= Ny;
		if (!(k<Nz)) k -= Nz;
*/
		nn = k*Nx*Ny+j*Nx+i;
		// unpack the distribution to the proper location
	//	if (recvbuf[start+idx] != dist[q*N+nn]){
	//		printf("Stopping to check error \n");
	//		printf("recvbuf[start+idx] = %f \n",recvbuf[start+idx]);
	//		printf("dist[q*N+nn] = %f \n",dist[q*N+nn]);
	//		printf("A bug! Again? \n");
	//		idx = count;
	//	}
//		list[idx] = nn;
		dist[q*N+nn] = recvbuf[start+idx];
	}
}


//************************************************************************* 
__global__ void INITIALIZE(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz, int S)
{
	int n,N;
	N = Nx*Ny*Nz;
	
	for (int s=0; s<S; s++){

		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		
		if (n<N){
			if (ID[n] > 0){
				f_even[n] = 0.3333333333333333;	
				f_odd[n] = 0.055555555555555555;		//double(100*n)+1.f;
				f_even[N+n] = 0.055555555555555555;	//double(100*n)+2.f;
				f_odd[N+n] = 0.055555555555555555;	//double(100*n)+3.f;
				f_even[2*N+n] = 0.055555555555555555;	//double(100*n)+4.f;
				f_odd[2*N+n] = 0.055555555555555555;	//double(100*n)+5.f;
				f_even[3*N+n] = 0.055555555555555555;	//double(100*n)+6.f;
				f_odd[3*N+n] = 0.0277777777777778;   //double(100*n)+7.f;
				f_even[4*N+n] = 0.0277777777777778;   //double(100*n)+8.f;
				f_odd[4*N+n] = 0.0277777777777778;   //double(100*n)+9.f;
				f_even[5*N+n] = 0.0277777777777778;  //double(100*n)+10.f;
				f_odd[5*N+n] = 0.0277777777777778;  //double(100*n)+11.f;
				f_even[6*N+n] = 0.0277777777777778;  //double(100*n)+12.f;
				f_odd[6*N+n] = 0.0277777777777778;  //double(100*n)+13.f;
				f_even[7*N+n] = 0.0277777777777778;  //double(100*n)+14.f;
				f_odd[7*N+n] = 0.0277777777777778;  //double(100*n)+15.f;
				f_even[8*N+n] = 0.0277777777777778;  //double(100*n)+16.f;
				f_odd[8*N+n] = 0.0277777777777778;  //double(100*n)+17.f;
				f_even[9*N+n] = 0.0277777777777778;  //double(100*n)+18.f;
			}
			else{
				for(int q=0; q<9; q++){
					f_even[q*N+n] = -1.0;
					f_odd[q*N+n] = -1.0;
				}
				f_even[9*N+n] = -1.0;
			}
		}
	}
}

__global__ void Compute_VELOCITY(char *ID, double *disteven, double *distodd, double *vel, int Nx, int Ny, int Nz, int S)
{
	int n,N;
	// distributions 
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double vx,vy,vz;
		
	N = Nx*Ny*Nz;

	// S - number of threadblocks per grid block
	for (int s=0; s<S; s++){

		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		
		if (n<N){
			if (ID[n] > 0){
				//........................................................................
				// Registers to store the distributions
				//........................................................................
				f2 = disteven[N+n];
				f4 = disteven[2*N+n];
				f6 = disteven[3*N+n];
				f8 = disteven[4*N+n];
				f10 = disteven[5*N+n];
				f12 = disteven[6*N+n];
				f14 = disteven[7*N+n];
				f16 = disteven[8*N+n];
				f18 = disteven[9*N+n];
				//........................................................................			
				f1 = distodd[n];
				f3 = distodd[1*N+n];
				f5 = distodd[2*N+n];
				f7 = distodd[3*N+n];
				f9 = distodd[4*N+n];
				f11 = distodd[5*N+n];
				f13 = distodd[6*N+n];
				f15 = distodd[7*N+n];
				f17 = distodd[8*N+n];
				//.................Compute the velocity...................................	
				vx = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
				vy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
				vz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
				//..................Write the velocity.....................................	
				vel[n] = vx;
				vel[N+n] = vy;
				vel[2*N+n] = vz;
				//........................................................................			

			}
		}
	}
}


//************************************************************************* 
__global__ void SWAP(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz, int S)
{
	int n,nn,N;
	// distributions 
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	
	N = Nx*Ny*Nz;
	
	// S - number of threadblocks per grid block
	for (int s=0; s<S; s++){
		
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;

//	for (n=0; n<N; n++){
		//.......Back out the 3-D indices for node n..............
		int	k = n/(Nx*Ny);
		int j = (n-Nx*Ny*k)/Nx;
		int i = n-Nx*Ny*k-Nz*j;
		
		if (n<N){
			if (ID[n] > 0){				
				//........................................................................
				// Retrieve even distributions from the local node (swap convention)
				//		f0 = disteven[n];  // Does not particupate in streaming
				f1 = distodd[n];
				f3 = distodd[N+n];
				f5 = distodd[2*N+n];
				f7 = distodd[3*N+n];
				f9 = distodd[4*N+n];
				f11 = distodd[5*N+n];
				f13 = distodd[6*N+n];
				f15 = distodd[7*N+n];
				f17 = distodd[8*N+n];
				//........................................................................
				
				//........................................................................
				// Retrieve odd distributions from neighboring nodes (swap convention)
				//........................................................................	
				nn = n+1;							// neighbor index (pull convention)
					if (!(i+1<Nx))	nn -= Nx;			// periodic BC along the x-boundary
				//if (i+1<Nx){
					f2 = disteven[N+nn];					// pull neighbor for distribution 2
					if (f2 > 0){
						distodd[n] = f2;
						disteven[N+nn] = f1;
					}
				//}
				//........................................................................	
				nn = n+Nx;							// neighbor index (pull convention)
				if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
				//if (j+1<Ny){
					f4 = disteven[2*N+nn];				// pull neighbor for distribution 4
					if (f4 > 0){
						distodd[N+n] = f4;
						disteven[2*N+nn] = f3;
				//	}
				}
				//........................................................................	
				nn = n+Nx*Ny;						// neighbor index (pull convention)
					if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
				//if (k+1<Nz){
					f6 = disteven[3*N+nn];				// pull neighbor for distribution 6
					if (f6 > 0){
						distodd[2*N+n] = f6;
						disteven[3*N+nn] = f5;
				//	}
				}
				//........................................................................	
				nn = n+Nx+1;						// neighbor index (pull convention)
					if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
					if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
				//if ((i+1<Nx) && (j+1<Ny)){
					f8 = disteven[4*N+nn];				// pull neighbor for distribution 8
					if (f8 > 0){
						distodd[3*N+n] = f8;
						disteven[4*N+nn] = f7;
				//	}
				}
				//........................................................................			
				nn = n-Nx+1;						// neighbor index (pull convention)
					if (!(i+1<Nx))	nn -= Nx;		// periodic BC along the x-boundary
					if (j-1<0)		nn += Nx*Ny;	// Perioidic BC along the y-boundary
				//if (!(i-1<0) && (j+1<Ny)){
					f10 = disteven[5*N+nn];					// pull neighbor for distribution 9
					if (f10 > 0){
						distodd[4*N+n] = f10;
						disteven[5*N+nn] = f9;
				//	}
				}
				//........................................................................	
				nn = n+Nx*Ny+1;						// neighbor index (pull convention)
					if (!(i+1<Nx))	nn -= Nx;		// periodic BC along the x-boundary
					if (!(k+1<Nz))	nn -= Nx*Ny*Nz;	// Perioidic BC along the z-boundary
				//if ( !(i-1<0) && !(k-1<0)){
					f12 = disteven[6*N+nn];				// pull distribution 11
					if (f12 > 0){
						distodd[5*N+n] = f12;
						disteven[6*N+nn] = f11;
				//	}
				}
				//........................................................................			
				nn = n-Nx*Ny+1;						// neighbor index (pull convention)
					if (!(i+1<Nx))	nn -= Nx;		// periodic BC along the x-boundary
					if (k-1<0)		nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
				//if (!(i-1<0) && (k+1<Nz)){
					f14 = disteven[7*N+nn];				// pull neighbor for distribution 13
					if (f14 > 0){
						distodd[6*N+n] = f14;
						disteven[7*N+nn] = f13;
				//	}
				}
				//........................................................................							
				nn = n+Nx*Ny+Nx;					// neighbor index (pull convention)
					if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
					if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary	
				//if (!(j-1<0) && !(k-1<0)){
					f16 = disteven[8*N+nn];				// pull neighbor for distribution 15
					if (f16 > 0){
						distodd[7*N+n] = f16;
						disteven[8*N+nn] = f15;
				//	}
				}
				//........................................................................											
				nn = n-Nx*Ny+Nx;					// neighbor index (pull convention)
					if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
					if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary	
				//if (!(j-1<0) && (k+1<Nz)){
					f18 = disteven[9*N+nn];				// pull neighbor for distribution 17
					if (f18 > 0){
						distodd[8*N+n] = f18;
						disteven[9*N+nn] = f17;
				//	}
				}
				//........................................................................		
				
			}
		}
	}
}
//************************************************************************* 

//************************************************************************* 
__global__ void MRT(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz, int S,
					double rlx_setA, double rlx_setB, double Fx, double Fy, double Fz)
{
		
	int n,N;
	// distributions 
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	
	// conserved momemnts
	double rho,jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;	
	
	N = Nx*Ny*Nz;
	
	char id;
	
	// S - number of threadblocks per grid block
	for (int s=0; s<S; s++){
//	for (int n=0; n<N; n++){
		//........Get 1-D index for this thread....................
		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		
		id = ID[n];
		
		if (n<N){
			if (id > 0){				
				//........................................................................
				// Registers to store the distributions - read based on swap convention
				//........................................................................
				f2 = distodd[n];
				f4 = distodd[N+n];
				f6 = distodd[2*N+n];
				f8 = distodd[3*N+n];
				f10 = distodd[4*N+n];
				f12 = distodd[5*N+n];
				f14 = distodd[6*N+n];
				f16 = distodd[7*N+n];
				f18 = distodd[8*N+n];
				//........................................................................			
				f0 = disteven[n];
				f1 = disteven[N+n];
				f3 = disteven[2*N+n];
				f5 = disteven[3*N+n];
				f7 = disteven[4*N+n];
				f9 = disteven[5*N+n];
				f11 = disteven[6*N+n];
				f13 = disteven[7*N+n];
				f15 = disteven[8*N+n];
				f17 = disteven[9*N+n];
				//........................................................................		
				//....................compute the moments...............................................		
				rho = f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
				m1 = -30*f0-11*(f2+f1+f4+f3+f6+f5)+8*(f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18 +f17);
				m2 = 12*f0-4*(f2+f1 +f4+f3+f6 +f5)+f8+f7+f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
				jx = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
				m4 = 4*(-f1+f2)+f7-f8+f9-f10+f11-f12+f13-f14;
				jy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
				m6 = -4*(f3-f4)+f7-f8-f9+f10+f15-f16+f17-f18;
				jz = f5-f6+f11-f12-f13+f14+f15-f16-f17+f18;
				m8 = -4*(f5-f6)+f11-f12-f13+f14+f15-f16-f17+f18;
				m9 = 2*(f1+f2)-f3-f4-f5-f6+f7+f8+f9+f10+f11+f12+f13+f14-2*(f15+f16+f17+f18);
				m10 = -4*(f1+f2)+2*(f4+f3+f6+f5)+f8+f7+f10+f9+f12+f11+f14+f13-2*(f16+f15+f18+f17);
				m11 = f4+f3-f6-f5+f8+f7+f10+f9-f12-f11-f14-f13;
				m12 = -2*(f4+f3-f6-f5)+f8+f7+f10+f9-f12-f11-f14-f13;
				m13 = f8+f7-f10-f9;
				m14 = f16+f15-f18-f17;
				m15 = f12+f11-f14-f13;
				m16 = f7-f8+f9-f10-f11+f12-f13+f14;
				m17 = -f7+f8+f9-f10+f15-f16+f17-f18;
				m18 = f11-f12-f13+f14-f15+f16+f17-f18;
				//..............incorporate external force................................................
				//jx += 0.5*Fx;
				//jy += 0.5*Fy;
				//jz += 0.5*Fz;
				//..............carry out relaxation process...............................................
				m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho - 11*rho) - m1);	
				m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho) - m2);
				m4 = m4 + rlx_setB*((-0.6666666666666666*jx) - m4);
				m6 = m6 + rlx_setB*((-0.6666666666666666*jy) - m6);
				m8 = m8 + rlx_setB*((-0.6666666666666666*jz) - m8);
				m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho) - m9);
				m10 = m10 + rlx_setA*(-0.5*((2*jx*jx-jy*jy-jz*jz)/rho) - m10);
				m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho) - m11);
				m12 = m12 + rlx_setA*(-0.5*((jy*jy-jz*jz)/rho) - m12);
				m13 = m13 + rlx_setA*((jx*jy/rho) - m13);
				m14 = m14 + rlx_setA*((jy*jz/rho) - m14);
				m15 = m15 + rlx_setA*((jx*jz/rho) - m15);
				m16 = m16 + rlx_setB*( - m16);
				m17 = m17 + rlx_setB*( - m17);
				m18 = m18 + rlx_setB*( - m18);
				//.................inverse transformation......................................................
				f0 = 0.05263157894736842*rho-0.012531328320802*m1+0.04761904761904762*m2;
				f1 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
				+0.1*(jx-m4)+0.05555555555555555*(m9-m10);		
				f2 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
				+0.1*(m4-jx)+0.05555555555555555*(m9-m10);	
				f3 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
				+0.1*(jy-m6)+0.02777777777777778*(m10-m9)+0.08333333333333333*(m11-m12);		
				f4 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
				+0.1*(m6-jy)+0.02777777777777778*(m10-m9)+0.08333333333333333*(m11-m12);
				f5 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
				+0.1*(jz-m8)+0.02777777777777778*(m10-m9)+0.08333333333333333*(m12-m11);		
				f6 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
				+0.1*(m8-jz)+0.02777777777777778*(m10-m9)+0.08333333333333333*(m12-m11);
				f7 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2+0.1*(jx+jy)+0.025*(m4+m6)
				+0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333*m11
				+0.04166666666666666*m12+0.25*m13+0.125*(m16-m17);
				f8 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2-0.1*(jx+jy)-0.025*(m4+m6)
				+0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333*m11
				+0.04166666666666666*m12+0.25*m13+0.125*(m17-m16);
				f9 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2+0.1*(jx-jy)+0.025*(m4-m6)
				+0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333*m11
				+0.04166666666666666*m12-0.25*m13+0.125*(m16+m17);
				f10 = 0.05263157894736842*rho+0.003341687552213868*m1+0.003968253968253968*m2+0.1*(jy-jx)+0.025*(m6-m4)
				+0.02777777777777778*m9+0.01388888888888889*m10+0.08333333333333333*m11
				+0.04166666666666666*m12-0.25*m13-0.125*(m16+m17);
				f11 = 0.05263157894736842*rho+0.003341687552213868*m1
				+0.003968253968253968*m2+0.1*(jx+jz)+0.025*(m4+m8)
				+0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333*m11
				-0.04166666666666666*m12+0.25*m15+0.125*(m18-m16);
				f12 = 0.05263157894736842*rho+0.003341687552213868*m1
				+0.003968253968253968*m2-0.1*(jx+jz)-0.025*(m4+m8)
				+0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333*m11
				-0.04166666666666666*m12+0.25*m15+0.125*(m16-m18);
				f13 = 0.05263157894736842*rho+0.003341687552213868*m1
				+0.003968253968253968*m2+0.1*(jx-jz)+0.025*(m4-m8)
				+0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333*m11
				-0.04166666666666666*m12-0.25*m15-0.125*(m16+m18);
				f14 = 0.05263157894736842*rho+0.003341687552213868*m1
				+0.003968253968253968*m2+0.1*(jz-jx)+0.025*(m8-m4)
				+0.02777777777777778*m9+0.01388888888888889*m10-0.08333333333333333*m11
				-0.04166666666666666*m12-0.25*m15+0.125*(m16+m18);
				f15 = 0.05263157894736842*rho+0.003341687552213868*m1
				+0.003968253968253968*m2+0.1*(jy+jz)+0.025*(m6+m8)
				-0.05555555555555555*m9-0.02777777777777778*m10+0.25*m14+0.125*(m17-m18);
				f16 =  0.05263157894736842*rho+0.003341687552213868*m1
				+0.003968253968253968*m2-0.1*(jy+jz)-0.025*(m6+m8)
				-0.05555555555555555*m9-0.02777777777777778*m10+0.25*m14+0.125*(m18-m17);
				f17 = 0.05263157894736842*rho+0.003341687552213868*m1
				+0.003968253968253968*m2+0.1*(jy-jz)+0.025*(m6-m8)
				-0.05555555555555555*m9-0.02777777777777778*m10-0.25*m14+0.125*(m17+m18);
				f18 = 0.05263157894736842*rho+0.003341687552213868*m1
				+0.003968253968253968*m2+0.1*(jz-jy)+0.025*(m8-m6)
				-0.05555555555555555*m9-0.02777777777777778*m10-0.25*m14-0.125*(m17+m18);
				//.......................................................................................................
				// incorporate external force
				f1 += 0.16666666*Fx;		
				f2 -= 0.16666666*Fx;		
				f3 += 0.16666666*Fy;		
				f4 -= 0.16666666*Fy;		
				f5 += 0.16666666*Fz;		
				f6 -= 0.16666666*Fz;	
				f7 += 0.08333333333*(Fx+Fy);	
				f8 -= 0.08333333333*(Fx+Fy);	
				f9 += 0.08333333333*(Fx-Fy);	
				f10 -= 0.08333333333*(Fx-Fy);	
				f11 += 0.08333333333*(Fx+Fz);	
				f12 -= 0.08333333333*(Fx+Fz);	
				f13 += 0.08333333333*(Fx-Fz);	
				f14 -= 0.08333333333*(Fx-Fz);	
				f15 += 0.08333333333*(Fy+Fz);	
				f16 -= 0.08333333333*(Fy+Fz);	
				f17 += 0.08333333333*(Fy-Fz);	
				f18 -= 0.08333333333*(Fy-Fz);	
				//.......................................................................................................
				// Write data based on un-swapped convention				
				disteven[n] = f0;
				disteven[N+n] = f2;
				disteven[2*N+n] = f4;
				disteven[3*N+n] = f6;
				disteven[4*N+n] = f8;
				disteven[5*N+n] = f10;
				disteven[6*N+n] = f12;
				disteven[7*N+n] = f14;
				disteven[8*N+n] = f16;
				disteven[9*N+n] = f18;
				
				distodd[n] = f1;
				distodd[N+n] = f3;
				distodd[2*N+n] = f5;
				distodd[3*N+n] = f7;				
				distodd[4*N+n] = f9;				
				distodd[5*N+n] = f11;
				distodd[6*N+n] = f13;
				distodd[7*N+n] = f15;
				distodd[8*N+n] = f17;
				//.......................................................................................................
			}
		}
	}
}
//************************************************************************* 

using namespace std;

void Write_Out(double *array, int Nx, int Ny, int Nz){
	int value;
	FILE *output;
	output = fopen("dist.list","w");
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				int index = k*Nx*Ny+j*Nx+i;
				value = int(array[index]);
				fprintf(output, "| %i",value);
			}
			fprintf(output, " | \n");
		}
		fprintf(output,"************************************** \n");	
	}
	fclose(output);
}

//************************************************************************* 
// MRT implementation of the LBM using CUDA
//************************************************************************* 

int main(int argc, char **argv)
{
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	// parallel domain size (# of sub-domains)
	int nprocx,nprocy,nprocz;
	int iproc,jproc,kproc;
	int sendtag,recvtag;
	//*****************************************
	// MPI ranks for all 18 neighbors
	//**********************************
	int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
	int rank_xy,rank_XY,rank_xY,rank_Xy;
	int rank_xz,rank_XZ,rank_xZ,rank_Xz;
	int rank_yz,rank_YZ,rank_yZ,rank_Yz;
	//**********************************
	MPI_Request req1[18],req2[18];
	MPI_Status stat1[18],stat2[18];
	//**********************************
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!!!!!!! Random debugging communications!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//	rank_X = rank+1;
//	if (!(rank_X < nprocs)) rank_X-=nprocs;
//	rank_x = rank-1;
//	if (rank_x < 0) rank_x +=nprocs;
//	rank_y = rank_z = rank_xy = rank_Xy = rank_xz = rank_Xz = rank_yz = rank_Yz = rank_x;
//	rank_Y = rank_Z = rank_XY = rank_xY = rank_XZ = rank_xZ = rank_YZ = rank_yZ = rank_X;
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	int device = 1;
	if (rank==0)	printf("Number of devices = %i \n", deviceCount);
	if (rank==0)	printf("Current device is = %i \n", device);
	cudaSetDevice(device);
	
	// BGK Model parameters
	string FILENAME;	
	unsigned int nBlocks, nthreads;
	int iterMax, interval;
	double tau,Fx,Fy,Fz,tol;
	// Domain variables
	int Nx,Ny,Nz;
	int i,j,k,n;

	if (rank==0){
		ifstream input("MRT.in");
		input >> FILENAME;		// name of the input file
		input >> Nz;			// number of nodes (x,y,z)
		input >> nBlocks;
		input >> nthreads;
		input >> tau;				// relaxation time
		input >> Fx;			// External force components (x,y,z)
		input >> Fy;
		input >> Fz;
		input >> iterMax;			// max no. of iterations
		input >> interval;			// error interval
		input >> tol;				// error tolerance

		ifstream domain("Domain.in");
		domain >> nprocx;
		domain >> nprocy;
		domain >> nprocz;
	}

	// **************************************************************
	// Broadcast simulation parameters from rank 0 to all other procs
	MPI_Barrier(comm);
	//.................................................
	MPI_Bcast(&Nz,1,MPI_INT,0,comm);
	MPI_Bcast(&nBlocks,1,MPI_INT,0,comm);
	MPI_Bcast(&nthreads,1,MPI_INT,0,comm);
	MPI_Bcast(&tau,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Fx,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Fy,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Fz,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&iterMax,1,MPI_INT,0,comm);
	MPI_Bcast(&interval,1,MPI_INT,0,comm);
	MPI_Bcast(&tol,1,MPI_DOUBLE,0,comm);

	MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
	//.................................................
	MPI_Barrier(comm);
	// **************************************************************

	double rlx_setA = 1.f/tau;
	double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);

	if (nprocs != nprocx*nprocy*nprocz){
		printf("Fatal error in processor number! \n");
		printf("nprocx =  %i \n",nprocx);
		printf("nprocy =  %i \n",nprocy);
		printf("nprocz =  %i \n",nprocz);
	}

	if (rank==0){
		printf("tau = %f \n", tau);
		printf("Set A = %f \n", rlx_setA);
		printf("Set B = %f \n", rlx_setB);
		printf("Force(x) = %f \n", Fx);
		printf("Force(y) = %f \n", Fy);
		printf("Force(z) = %f \n", Fz);
		printf("Sub-domain size = %i x %i x %i\n",Nz,Nz,Nz);
	}

	MPI_Barrier(comm);
	kproc = rank/(nprocx*nprocy);
	jproc = (rank-nprocx*nprocy*kproc)/nprocx;
	iproc = rank-nprocx*nprocy*kproc-nprocz*jproc;

	//..........................................
	// set up the neighbor ranks
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=1;
	j+=0;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_X = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i-=1;
	j+=0;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_x = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j+=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_Y = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j-=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_y = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j+=0;
	k+=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-= nprocy;
	if (!(k<nprocz)) k-= nprocz;
	rank_Z = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j+=0;
	k-=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-= nprocy;
	if (!(k<nprocz)) k-= nprocz;
	rank_z = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=1;
	j+=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_XY = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i-=1;
	j-=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_xy = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=1;
	j-=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_Xy = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i-=1;
	j+=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_xY = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=1;
	j+=0;
	k+=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_XZ = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i-=1;
	j+=0;
	k-=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_xz = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i-=1;
	j+=0;
	k+=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_xZ = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=1;
	j+=0;
	k-=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_Xz = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j+=1;
	k+=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_YZ = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j-=1;
	k-=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_yz = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j-=1;
	k+=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_yZ = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k=kproc;
	i+=0;
	j+=1;
	k-=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_Yz = k*nprocx*nprocy+j*nprocx+i;
	//..........................................

	Nz += 2;
	Nx = Ny = Nz;	// Cubic domain
	
	int N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);
	
//	unsigned int nBlocks = 32;
//	int nthreads = 128;
	int S = N/nthreads/nBlocks;
	
//	unsigned int nBlocks = N/nthreads + (N%nthreads == 0?0:1);
	dim3 grid(nBlocks,1,1);
		
	if (rank==0) printf("Number of blocks = %i \n", nBlocks);
	if (rank==0) printf("Threads per block = %i \n", nthreads);
	if (rank==0) printf("Sweeps per thread = %i \n", S);
	if (rank==0) printf("Number of nodes per side = %i \n", Nx);
	if (rank==0) printf("Total Number of nodes = %i \n", N);
	
	//.......................................................................
	if (rank == 0)	printf("Read input media... \n");
	//.......................................................................
	char LocalRankString[8];
	char LocalRankFilename[40];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	printf("Local File Name =  %s \n",LocalRankFilename);
	// .......... READ THE INPUT FILE .......................................
	char value;
	char *id;
	id = new char[N];	
	int sum = 0;
	double porosity;
	//.......................................................................
	ifstream PM(LocalRankFilename,ios::binary);
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				id[n] = 0;
			}
		}
	}
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				PM.read((char *) (&value), sizeof(value));
				n = k*Nx*Ny+j*Nx+i;
				id[n] = value;
				if (value > 0) sum++;
			}
		}
	}
	PM.close();
//	printf("File porosity = %f\n", double(sum)/N);
	//...........................................................................
	MPI_Barrier(comm);
	if (rank == 0) cout << "Domain set." << endl;
	//...........................................................................
	// Write the communcation structure into a file for debugging
	char LocalCommFile[40];
	sprintf(LocalCommFile,"%s%s","Comm.",LocalRankString);
	FILE *CommFile;
	CommFile = fopen(LocalCommFile,"w");
	fprintf(CommFile,"rank=%d, ",rank);
	fprintf(CommFile,"i=%d,j=%d,k=%d :",iproc,jproc,kproc);
	fprintf(CommFile,"x=%d, ",rank_x);
	fprintf(CommFile,"X=%d, ",rank_X);
	fprintf(CommFile,"y=%d, ",rank_y);
	fprintf(CommFile,"Y=%d, ",rank_Y);
	fprintf(CommFile,"z=%d, ",rank_z);
	fprintf(CommFile,"Z=%d, ",rank_Z);
	fprintf(CommFile,"xy=%d, ",rank_xy);
	fprintf(CommFile,"XY=%d, ",rank_XY);
	fprintf(CommFile,"xY=%d, ",rank_xY);
	fprintf(CommFile,"Xy=%d, ",rank_Xy);
	fprintf(CommFile,"xz=%d, ",rank_xz);
	fprintf(CommFile,"XZ=%d, ",rank_XZ);
	fprintf(CommFile,"xZ=%d, ",rank_xZ);
	fprintf(CommFile,"Xz=%d, ",rank_Xz);
	fprintf(CommFile,"yz=%d, ",rank_yz);
	fprintf(CommFile,"YZ=%d, ",rank_YZ);
	fprintf(CommFile,"yZ=%d, ",rank_yZ);
	fprintf(CommFile,"Yz=%d, ",rank_Yz);
	fprintf(CommFile,"\n");
	fclose(CommFile);
	//...........................................................................

	// Set up MPI communication structures
	if (rank==0)	printf ("Setting up communication control structures \n");
	//......................................................................................
	// Get the actual D3Q19 communication counts (based on location of solid phase)
	// Discrete velocity set symmetry implies the sendcount = recvcount
	int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z;
	int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz, sendCount_xZ;
	int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ, sendCount_XZ;
	sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
	sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
	sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;
	//......................................................................................
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				// Check the phase ID
				if (id[k*Nx*Ny+j*Nx+i] != 0){
					// Counts for the six faces
					if (i==1)	sendCount_x++;
					if (j==1)	sendCount_y++;
					if (k==1)	sendCount_z++;
					if (i==Nx-2)	sendCount_X++;
					if (j==Ny-2)	sendCount_Y++;
					if (k==Nz-2)	sendCount_Z++;
					// Counts for the twelve edges
					if (i==1 && j==1)	sendCount_xy++;
					if (i==1 && j==Ny-2)	sendCount_xY++;
					if (i==Nx-2 && j==1)	sendCount_Xy++;
					if (i==Nx-2 && j==Ny-2)	sendCount_XY++;

					if (i==1 && k==1)	sendCount_xz++;
					if (i==1 && k==Nz-2)	sendCount_xZ++;
					if (i==Nx-2 && k==1)	sendCount_Xz++;
					if (i==Nx-2 && k==Nz-2)	sendCount_XZ++;

					if (j==1 && k==1)	sendCount_yz++;
					if (j==1 && k==Nz-2)	sendCount_yZ++;
					if (j==Ny-2 && k==1)	sendCount_Yz++;
					if (j==Ny-2 && k==Nz-2)	sendCount_YZ++;
				}
			}
		}
	}
	//......................................................................................
	int *sendList_x, *sendList_y, *sendList_z, *sendList_X, *sendList_Y, *sendList_Z;
	int *sendList_xy, *sendList_yz, *sendList_xz, *sendList_Xy, *sendList_Yz, *sendList_xZ;
	int *sendList_xY, *sendList_yZ, *sendList_Xz, *sendList_XY, *sendList_YZ, *sendList_XZ;
	//......................................................................................
	// send buffers
	sendList_x = new int [sendCount_x];
	sendList_y = new int [sendCount_y];
	sendList_z = new int [sendCount_z];
	sendList_X = new int [sendCount_X];
	sendList_Y = new int [sendCount_Y];
	sendList_Z = new int [sendCount_Z];
	sendList_xy = new int [sendCount_xy];
	sendList_yz = new int [sendCount_yz];
	sendList_xz = new int [sendCount_xz];
	sendList_Xy = new int [sendCount_Xy];
	sendList_Yz = new int [sendCount_Yz];
	sendList_xZ = new int [sendCount_xZ];
	sendList_xY = new int [sendCount_xY];
	sendList_yZ = new int [sendCount_yZ];
	sendList_Xz = new int [sendCount_Xz];
	sendList_XY = new int [sendCount_XY];
	sendList_YZ = new int [sendCount_YZ];
	sendList_XZ = new int [sendCount_XZ];
	if (rank==0)	printf ("Preparing the sendlists \n");
	//......................................................................................
	// Populate the send list
	sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
	sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
	sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				// Local value to send
				n = k*Nx*Ny+j*Nx+i;
				if (id[n] != 0){
					// Counts for the six faces
					if (i==1)		sendList_x[sendCount_x++]=n;
					if (j==1)		sendList_y[sendCount_y++]=n;
					if (k==1)		sendList_z[sendCount_z++]=n;
					if (i==Nx-2)	sendList_X[sendCount_X++]=n;
					if (j==Ny-2)	sendList_Y[sendCount_Y++]=n;
					if (k==Nz-2)	sendList_Z[sendCount_Z++]=n;
					// Counts for the twelve edges
					if (i==1 && j==1)		sendList_xy[sendCount_xy++]=n;
					if (i==1 && j==Ny-2)	sendList_xY[sendCount_xY++]=n;
					if (i==Nx-2 && j==1)	sendList_Xy[sendCount_Xy++]=n;
					if (i==Nx-2 && j==Ny-2)	sendList_XY[sendCount_XY++]=n;

					if (i==1 && k==1)		sendList_xz[sendCount_xz++]=n;
					if (i==1 && k==Nz-2)	sendList_xZ[sendCount_xZ++]=n;
					if (i==Nx-2 && k==1)	sendList_Xz[sendCount_Xz++]=n;
					if (i==Nx-2 && k==Nz-2)	sendList_XZ[sendCount_XZ++]=n;

					if (j==1 && k==1)		sendList_yz[sendCount_yz++]=n;
					if (j==1 && k==Nz-2)	sendList_yZ[sendCount_yZ++]=n;
					if (j==Ny-2 && k==1)	sendList_Yz[sendCount_Yz++]=n;
					if (j==Ny-2 && k==Nz-2)	sendList_YZ[sendCount_YZ++]=n;
				}
			}
		}
	}
	MPI_Barrier(comm);
	if (rank==0)	printf ("SendLists are ready on host\n");
	//......................................................................................
	// Use MPI to fill in the recvCounts form the associated processes
	int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z;
	int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz, recvCount_xZ;
	int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ, recvCount_XZ;
	//......................................................................................
	//**********************************************************************************
	// Fill in the recieve counts using MPI
	sendtag = recvtag = 3;
	MPI_Send(&sendCount_x,1,MPI_INT,rank_X,sendtag,comm);
	MPI_Recv(&recvCount_X,1,MPI_INT,rank_x,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_X,1,MPI_INT,rank_x,sendtag,comm);
	MPI_Recv(&recvCount_x,1,MPI_INT,rank_X,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_y,1,MPI_INT,rank_Y,sendtag,comm);
	MPI_Recv(&recvCount_Y,1,MPI_INT,rank_y,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_Y,1,MPI_INT,rank_y,sendtag,comm);
	MPI_Recv(&recvCount_y,1,MPI_INT,rank_Y,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_z,1,MPI_INT,rank_Z,sendtag,comm);
	MPI_Recv(&recvCount_Z,1,MPI_INT,rank_z,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_Z,1,MPI_INT,rank_z,sendtag,comm);
	MPI_Recv(&recvCount_z,1,MPI_INT,rank_Z,recvtag,comm,MPI_STATUS_IGNORE);

	MPI_Send(&sendCount_xy,1,MPI_INT,rank_XY,sendtag,comm);
	MPI_Recv(&recvCount_XY,1,MPI_INT,rank_xy,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_XY,1,MPI_INT,rank_xy,sendtag,comm);
	MPI_Recv(&recvCount_xy,1,MPI_INT,rank_XY,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_Xy,1,MPI_INT,rank_xY,sendtag,comm);
	MPI_Recv(&recvCount_xY,1,MPI_INT,rank_Xy,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_xY,1,MPI_INT,rank_Xy,sendtag,comm);
	MPI_Recv(&recvCount_Xy,1,MPI_INT,rank_xY,recvtag,comm,MPI_STATUS_IGNORE);

	MPI_Send(&sendCount_xz,1,MPI_INT,rank_XZ,sendtag,comm);
	MPI_Recv(&recvCount_XZ,1,MPI_INT,rank_xz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_XZ,1,MPI_INT,rank_xz,sendtag,comm);
	MPI_Recv(&recvCount_xz,1,MPI_INT,rank_XZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_Xz,1,MPI_INT,rank_xZ,sendtag,comm);
	MPI_Recv(&recvCount_xZ,1,MPI_INT,rank_Xz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_xZ,1,MPI_INT,rank_Xz,sendtag,comm);
	MPI_Recv(&recvCount_Xz,1,MPI_INT,rank_xZ,recvtag,comm,MPI_STATUS_IGNORE);

	MPI_Send(&sendCount_yz,1,MPI_INT,rank_YZ,sendtag,comm);
	MPI_Recv(&recvCount_YZ,1,MPI_INT,rank_yz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_YZ,1,MPI_INT,rank_yz,sendtag,comm);
	MPI_Recv(&recvCount_yz,1,MPI_INT,rank_YZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_Yz,1,MPI_INT,rank_yZ,sendtag,comm);
	MPI_Recv(&recvCount_yZ,1,MPI_INT,rank_Yz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_yZ,1,MPI_INT,rank_Yz,sendtag,comm);
	MPI_Recv(&recvCount_Yz,1,MPI_INT,rank_yZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Barrier(comm);
	//**********************************************************************************
	//recvCount_x = sendCount_x;
	//recvCount_X = sendCount_X;
	//recvCount_y = sendCount_y;
	//recvCount_Y = sendCount_Y;
	//recvCount_z = sendCount_z;
	//recvCount_Z = sendCount_Z;
	//recvCount_xy = sendCount_xy;
	//recvCount_xY = sendCount_xY;
	//recvCount_Xy = sendCount_Xy;
	//recvCount_XY = sendCount_XY;
	//recvCount_xz = sendCount_xz;
	//recvCount_xZ = sendCount_xZ;
	//recvCount_Xz = sendCount_XZ;
	//recvCount_XZ = sendCount_XZ;
	//recvCount_yz = sendCount_yz;
	//recvCount_Yz = sendCount_Yz;
	//recvCount_yZ = sendCount_yZ;
	//recvCount_YZ = sendCount_YZ;
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//......................................................................................
	// Use MPI to fill in the appropriate values
	//	int tag = 5;
	//	MPI_Sendrecv(sendCount_x,1,MPI_INT,rank_x,tag,sendCount_X,1,MPI_INT,comm,req);
	//......................................................................................
	int *recvList_x, *recvList_y, *recvList_z, *recvList_X, *recvList_Y, *recvList_Z;
	int *recvList_xy, *recvList_yz, *recvList_xz, *recvList_Xy, *recvList_Yz, *recvList_xZ;
	int *recvList_xY, *recvList_yZ, *recvList_Xz, *recvList_XY, *recvList_YZ, *recvList_XZ;
	//......................................................................................
	// recv buffers
	recvList_x = new int [recvCount_x];
	recvList_y = new int [recvCount_y];
	recvList_z = new int [recvCount_z];
	recvList_X = new int [recvCount_X];
	recvList_Y = new int [recvCount_Y];
	recvList_Z = new int [recvCount_Z];
	recvList_xy = new int [recvCount_xy];
	recvList_yz = new int [recvCount_yz];
	recvList_xz = new int [recvCount_xz];
	recvList_Xy = new int [recvCount_Xy];
	recvList_Yz = new int [recvCount_Yz];
	recvList_xZ = new int [recvCount_xZ];
	recvList_xY = new int [recvCount_xY];
	recvList_yZ = new int [recvCount_yZ];
	recvList_Xz = new int [recvCount_Xz];
	recvList_XY = new int [recvCount_XY];
	recvList_YZ = new int [recvCount_YZ];
	recvList_XZ = new int [recvCount_XZ];
	//......................................................................................
	//......................................................................................
	// Use MPI to fill in the appropriate values for recvList
	// Fill in the recieve lists using MPI
	sendtag = recvtag = 4;
	MPI_Isend(sendList_x, sendCount_x,MPI_INT,rank_X,sendtag,comm,&req1[0]);
	MPI_Irecv(recvList_X, recvCount_X,MPI_INT,rank_x,recvtag,comm,&req2[0]);
	MPI_Isend(sendList_X, sendCount_X,MPI_INT,rank_x,sendtag,comm,&req1[1]);
	MPI_Irecv(recvList_x, recvCount_x,MPI_INT,rank_X,recvtag,comm,&req2[1]);
	MPI_Isend(sendList_y, sendCount_y,MPI_INT,rank_Y,sendtag,comm,&req1[2]);
	MPI_Irecv(recvList_Y, recvCount_Y,MPI_INT,rank_y,recvtag,comm,&req2[2]);
	MPI_Isend(sendList_Y, sendCount_Y,MPI_INT,rank_y,sendtag,comm,&req1[3]);
	MPI_Irecv(recvList_y, recvCount_y,MPI_INT,rank_Y,recvtag,comm,&req2[3]);
	MPI_Isend(sendList_z, sendCount_z,MPI_INT,rank_Z,sendtag,comm,&req1[4]);
	MPI_Irecv(recvList_Z, recvCount_Z,MPI_INT,rank_z,recvtag,comm,&req2[4]);
	MPI_Isend(sendList_Z, sendCount_Z,MPI_INT,rank_z,sendtag,comm,&req1[5]);
	MPI_Irecv(recvList_z, recvCount_z,MPI_INT,rank_Z,recvtag,comm,&req2[5]);

	MPI_Isend(sendList_xy, sendCount_xy,MPI_INT,rank_XY,sendtag,comm,&req1[6]);
	MPI_Irecv(recvList_XY, recvCount_XY,MPI_INT,rank_xy,recvtag,comm,&req2[6]);
	MPI_Isend(sendList_XY, sendCount_XY,MPI_INT,rank_xy,sendtag,comm,&req1[7]);
	MPI_Irecv(recvList_xy, recvCount_xy,MPI_INT,rank_XY,recvtag,comm,&req2[7]);
	MPI_Isend(sendList_Xy, sendCount_Xy,MPI_INT,rank_xY,sendtag,comm,&req1[8]);
	MPI_Irecv(recvList_xY, recvCount_xY,MPI_INT,rank_Xy,recvtag,comm,&req2[8]);
	MPI_Isend(sendList_xY, sendCount_xY,MPI_INT,rank_Xy,sendtag,comm,&req1[9]);
	MPI_Irecv(recvList_Xy, recvCount_Xy,MPI_INT,rank_xY,recvtag,comm,&req2[9]);

	MPI_Isend(sendList_xz, sendCount_xz,MPI_INT,rank_XZ,sendtag,comm,&req1[10]);
	MPI_Irecv(recvList_XZ, recvCount_XZ,MPI_INT,rank_xz,recvtag,comm,&req2[10]);
	MPI_Isend(sendList_XZ, sendCount_XZ,MPI_INT,rank_xz,sendtag,comm,&req1[11]);
	MPI_Irecv(recvList_xz, recvCount_xz,MPI_INT,rank_XZ,recvtag,comm,&req2[11]);
	MPI_Isend(sendList_Xz, sendCount_Xz,MPI_INT,rank_xZ,sendtag,comm,&req1[12]);
	MPI_Irecv(recvList_xZ, recvCount_xZ,MPI_INT,rank_Xz,recvtag,comm,&req2[12]);
	MPI_Isend(sendList_xZ, sendCount_xZ,MPI_INT,rank_Xz,sendtag,comm,&req1[13]);
	MPI_Irecv(recvList_Xz, recvCount_Xz,MPI_INT,rank_xZ,recvtag,comm,&req2[13]);

	MPI_Isend(sendList_yz, sendCount_yz,MPI_INT,rank_YZ,sendtag,comm,&req1[14]);
	MPI_Irecv(recvList_YZ, recvCount_YZ,MPI_INT,rank_yz,recvtag,comm,&req2[14]);
	MPI_Isend(sendList_YZ, sendCount_YZ,MPI_INT,rank_yz,sendtag,comm,&req1[15]);
	MPI_Irecv(recvList_yz, recvCount_yz,MPI_INT,rank_YZ,recvtag,comm,&req2[15]);
	MPI_Isend(sendList_Yz, sendCount_Yz,MPI_INT,rank_yZ,sendtag,comm,&req1[16]);
	MPI_Irecv(recvList_yZ, recvCount_yZ,MPI_INT,rank_Yz,recvtag,comm,&req2[16]);
	MPI_Isend(sendList_yZ, sendCount_yZ,MPI_INT,rank_Yz,sendtag,comm,&req1[17]);
	MPI_Irecv(recvList_Yz, recvCount_Yz,MPI_INT,rank_yZ,recvtag,comm,&req2[17]);
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	MPI_Barrier(comm);
	//......................................................................................
	double *sendbuf_x, *sendbuf_y, *sendbuf_z, *sendbuf_X, *sendbuf_Y, *sendbuf_Z;
	double *sendbuf_xy, *sendbuf_yz, *sendbuf_xz, *sendbuf_Xy, *sendbuf_Yz, *sendbuf_xZ;
	double *sendbuf_xY, *sendbuf_yZ, *sendbuf_Xz, *sendbuf_XY, *sendbuf_YZ, *sendbuf_XZ;
	double *recvbuf_x, *recvbuf_y, *recvbuf_z, *recvbuf_X, *recvbuf_Y, *recvbuf_Z;
	double *recvbuf_xy, *recvbuf_yz, *recvbuf_xz, *recvbuf_Xy, *recvbuf_Yz, *recvbuf_xZ;
	double *recvbuf_xY, *recvbuf_yZ, *recvbuf_Xz, *recvbuf_XY, *recvbuf_YZ, *recvbuf_XZ;
	//......................................................................................
	cudaMalloc((void **) &sendbuf_x, 5*sendCount_x*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_X, 5*sendCount_X*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_y, 5*sendCount_y*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_Y, 5*sendCount_Y*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_z, 5*sendCount_z*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_Z, 5*sendCount_Z*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_xy, sendCount_xy*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_xY, sendCount_xY*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_Xy, sendCount_Xy*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_XY, sendCount_XY*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_xz, sendCount_xz*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_xZ, sendCount_xZ*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_Xz, sendCount_Xz*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_XZ, sendCount_XZ*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_yz, sendCount_yz*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_yZ, sendCount_yZ*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_Yz, sendCount_Yz*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &sendbuf_YZ, sendCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	cudaMalloc((void **) &recvbuf_x, 5*recvCount_x*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_X, 5*recvCount_X*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_y, 5*recvCount_y*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_Y, 5*recvCount_Y*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_z, 5*recvCount_z*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_Z, 5*recvCount_Z*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_xy, recvCount_xy*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_xY, recvCount_xY*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_Xy, recvCount_Xy*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_XY, recvCount_XY*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_xz, recvCount_xz*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_xZ, recvCount_xZ*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_Xz, recvCount_Xz*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_XZ, recvCount_XZ*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_yz, recvCount_yz*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_yZ, recvCount_yZ*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_Yz, recvCount_Yz*sizeof(double));	// Allocate device memory
	cudaMalloc((void **) &recvbuf_YZ, recvCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	int *dvcSendList_x, *dvcSendList_y, *dvcSendList_z, *dvcSendList_X, *dvcSendList_Y, *dvcSendList_Z;
	int *dvcSendList_xy, *dvcSendList_yz, *dvcSendList_xz, *dvcSendList_Xy, *dvcSendList_Yz, *dvcSendList_xZ;
	int *dvcSendList_xY, *dvcSendList_yZ, *dvcSendList_Xz, *dvcSendList_XY, *dvcSendList_YZ, *dvcSendList_XZ;
	//......................................................................................
	int *dvcRecvList_x, *dvcRecvList_y, *dvcRecvList_z, *dvcRecvList_X, *dvcRecvList_Y, *dvcRecvList_Z;
	int *dvcRecvList_xy, *dvcRecvList_yz, *dvcRecvList_xz, *dvcRecvList_Xy, *dvcRecvList_Yz, *dvcRecvList_xZ;
	int *dvcRecvList_xY, *dvcRecvList_yZ, *dvcRecvList_Xz, *dvcRecvList_XY, *dvcRecvList_YZ, *dvcRecvList_XZ;
	//......................................................................................
	cudaMalloc((void **) &dvcSendList_x, sendCount_x*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_X, sendCount_X*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_y, sendCount_y*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_Y, sendCount_Y*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_z, sendCount_z*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_Z, sendCount_Z*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_xy, sendCount_xy*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_xY, sendCount_xY*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_Xy, sendCount_Xy*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_XY, sendCount_XY*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_xz, sendCount_xz*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_xZ, sendCount_xZ*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_Xz, sendCount_Xz*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_XZ, sendCount_XZ*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_yz, sendCount_yz*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_yZ, sendCount_yZ*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_Yz, sendCount_Yz*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcSendList_YZ, sendCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	cudaMalloc((void **) &dvcRecvList_x, recvCount_x*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_X, recvCount_X*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_y, recvCount_y*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_Y, recvCount_Y*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_z, recvCount_z*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_Z, recvCount_Z*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_xy, recvCount_xy*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_xY, recvCount_xY*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_Xy, recvCount_Xy*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_XY, recvCount_XY*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_xz, recvCount_xz*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_xZ, recvCount_xZ*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_Xz, recvCount_Xz*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_XZ, recvCount_XZ*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_yz, recvCount_yz*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_yZ, recvCount_yZ*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_Yz, recvCount_Yz*sizeof(int));	// Allocate device memory
	cudaMalloc((void **) &dvcRecvList_YZ, recvCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	if (rank==0)	printf ("Prepare to copy send/recv Lists to device \n");
	cudaMemcpy(dvcSendList_x,sendList_x,sendCount_x*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_X,sendList_X,sendCount_X*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_y,sendList_y,sendCount_y*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_Y,sendList_Y,sendCount_Y*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_z,sendList_z,sendCount_z*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_Z,sendList_Z,sendCount_Z*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_xy,sendList_xy,sendCount_xy*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_XY,sendList_XY,sendCount_XY*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_xY,sendList_xY,sendCount_xY*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_Xy,sendList_Xy,sendCount_Xy*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_xz,sendList_xz,sendCount_xz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_XZ,sendList_XZ,sendCount_XZ*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_xZ,sendList_xZ,sendCount_xZ*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_Xz,sendList_Xz,sendCount_Xz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_yz,sendList_yz,sendCount_yz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_YZ,sendList_YZ,sendCount_YZ*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_yZ,sendList_yZ,sendCount_yZ*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSendList_Yz,sendList_Yz,sendCount_Yz*sizeof(int),cudaMemcpyHostToDevice);
	//......................................................................................
	cudaMemcpy(dvcRecvList_x,recvList_x,recvCount_x*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_X,recvList_X,recvCount_X*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_y,recvList_y,recvCount_y*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_Y,recvList_Y,recvCount_Y*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_z,recvList_z,recvCount_z*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_Z,recvList_Z,recvCount_Z*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_xy,recvList_xy,recvCount_xy*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_XY,recvList_XY,recvCount_XY*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_xY,recvList_xY,recvCount_xY*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_Xy,recvList_Xy,recvCount_Xy*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_xz,recvList_xz,recvCount_xz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_XZ,recvList_XZ,recvCount_XZ*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_xZ,recvList_xZ,recvCount_xZ*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_Xz,recvList_Xz,recvCount_Xz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_yz,recvList_yz,recvCount_yz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_YZ,recvList_YZ,recvCount_YZ*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_yZ,recvList_yZ,recvCount_yZ*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dvcRecvList_Yz,recvList_Yz,recvCount_Yz*sizeof(int),cudaMemcpyHostToDevice);
	//......................................................................................
	// Fill in the phase ID from neighboring processors
	char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
	char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
	char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
	char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
	char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
	char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
	// send buffers
	sendID_x = new char [sendCount_x];
	sendID_y = new char [sendCount_y];
	sendID_z = new char [sendCount_z];
	sendID_X = new char [sendCount_X];
	sendID_Y = new char [sendCount_Y];
	sendID_Z = new char [sendCount_Z];
	sendID_xy = new char [sendCount_xy];
	sendID_yz = new char [sendCount_yz];
	sendID_xz = new char [sendCount_xz];
	sendID_Xy = new char [sendCount_Xy];
	sendID_Yz = new char [sendCount_Yz];
	sendID_xZ = new char [sendCount_xZ];
	sendID_xY = new char [sendCount_xY];
	sendID_yZ = new char [sendCount_yZ];
	sendID_Xz = new char [sendCount_Xz];
	sendID_XY = new char [sendCount_XY];
	sendID_YZ = new char [sendCount_YZ];
	sendID_XZ = new char [sendCount_XZ];
	//......................................................................................
	// recv buffers
	recvID_x = new char [recvCount_x];
	recvID_y = new char [recvCount_y];
	recvID_z = new char [recvCount_z];
	recvID_X = new char [recvCount_X];
	recvID_Y = new char [recvCount_Y];
	recvID_Z = new char [recvCount_Z];
	recvID_xy = new char [recvCount_xy];
	recvID_yz = new char [recvCount_yz];
	recvID_xz = new char [recvCount_xz];
	recvID_Xy = new char [recvCount_Xy];
	recvID_xZ = new char [recvCount_xZ];
	recvID_xY = new char [recvCount_xY];
	recvID_yZ = new char [recvCount_yZ];
	recvID_Yz = new char [recvCount_Yz];
	recvID_Xz = new char [recvCount_Xz];
	recvID_XY = new char [recvCount_XY];
	recvID_YZ = new char [recvCount_YZ];
	recvID_XZ = new char [recvCount_XZ];
	//......................................................................................
	sendtag = recvtag = 7;
	PackID(sendList_x, sendCount_x ,sendID_x, id);
	PackID(sendList_X, sendCount_X ,sendID_X, id);
	PackID(sendList_y, sendCount_y ,sendID_y, id);
	PackID(sendList_Y, sendCount_Y ,sendID_Y, id);
	PackID(sendList_z, sendCount_z ,sendID_z, id);
	PackID(sendList_Z, sendCount_Z ,sendID_Z, id);
	PackID(sendList_xy, sendCount_xy ,sendID_xy, id);
	PackID(sendList_Xy, sendCount_Xy ,sendID_Xy, id);
	PackID(sendList_xY, sendCount_xY ,sendID_xY, id);
	PackID(sendList_XY, sendCount_XY ,sendID_XY, id);
	PackID(sendList_xz, sendCount_xz ,sendID_xz, id);
	PackID(sendList_Xz, sendCount_Xz ,sendID_Xz, id);
	PackID(sendList_xZ, sendCount_xZ ,sendID_xZ, id);
	PackID(sendList_XZ, sendCount_XZ ,sendID_XZ, id);
	PackID(sendList_yz, sendCount_yz ,sendID_yz, id);
	PackID(sendList_Yz, sendCount_Yz ,sendID_Yz, id);
	PackID(sendList_yZ, sendCount_yZ ,sendID_yZ, id);
	PackID(sendList_YZ, sendCount_YZ ,sendID_YZ, id);
	//......................................................................................
	MPI_Sendrecv(sendID_x,sendCount_x,MPI_CHAR,rank_X,sendtag,
			recvID_X,recvCount_X,MPI_CHAR,rank_x,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_X,sendCount_X,MPI_CHAR,rank_x,sendtag,
			recvID_x,recvCount_x,MPI_CHAR,rank_X,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_y,sendCount_y,MPI_CHAR,rank_Y,sendtag,
			recvID_Y,recvCount_Y,MPI_CHAR,rank_y,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Y,sendCount_Y,MPI_CHAR,rank_y,sendtag,
			recvID_y,recvCount_y,MPI_CHAR,rank_Y,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_z,sendCount_z,MPI_CHAR,rank_Z,sendtag,
			recvID_Z,recvCount_Z,MPI_CHAR,rank_z,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Z,sendCount_Z,MPI_CHAR,rank_z,sendtag,
			recvID_z,recvCount_z,MPI_CHAR,rank_Z,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xy,sendCount_xy,MPI_CHAR,rank_XY,sendtag,
			recvID_XY,recvCount_XY,MPI_CHAR,rank_xy,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XY,sendCount_XY,MPI_CHAR,rank_xy,sendtag,
			recvID_xy,recvCount_xy,MPI_CHAR,rank_XY,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xy,sendCount_Xy,MPI_CHAR,rank_xY,sendtag,
			recvID_xY,recvCount_xY,MPI_CHAR,rank_Xy,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xY,sendCount_xY,MPI_CHAR,rank_Xy,sendtag,
			recvID_Xy,recvCount_Xy,MPI_CHAR,rank_xY,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xz,sendCount_xz,MPI_CHAR,rank_XZ,sendtag,
			recvID_XZ,recvCount_XZ,MPI_CHAR,rank_xz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XZ,sendCount_XZ,MPI_CHAR,rank_xz,sendtag,
			recvID_xz,recvCount_xz,MPI_CHAR,rank_XZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xz,sendCount_Xz,MPI_CHAR,rank_xZ,sendtag,
			recvID_xZ,recvCount_xZ,MPI_CHAR,rank_Xz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xZ,sendCount_xZ,MPI_CHAR,rank_Xz,sendtag,
			recvID_Xz,recvCount_Xz,MPI_CHAR,rank_xZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yz,sendCount_yz,MPI_CHAR,rank_YZ,sendtag,
			recvID_YZ,recvCount_YZ,MPI_CHAR,rank_yz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_YZ,sendCount_YZ,MPI_CHAR,rank_yz,sendtag,
			recvID_yz,recvCount_yz,MPI_CHAR,rank_YZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Yz,sendCount_Yz,MPI_CHAR,rank_yZ,sendtag,
			recvID_yZ,recvCount_yZ,MPI_CHAR,rank_Yz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yZ,sendCount_yZ,MPI_CHAR,rank_Yz,sendtag,
			recvID_Yz,recvCount_Yz,MPI_CHAR,rank_yZ,recvtag,comm,MPI_STATUS_IGNORE);
	//......................................................................................
	UnpackID(recvList_x, recvCount_x ,recvID_x, id);
	UnpackID(recvList_X, recvCount_X ,recvID_X, id);
	UnpackID(recvList_y, recvCount_y ,recvID_y, id);
	UnpackID(recvList_Y, recvCount_Y ,recvID_Y, id);
	UnpackID(recvList_z, recvCount_z ,recvID_z, id);
	UnpackID(recvList_Z, recvCount_Z ,recvID_Z, id);
	UnpackID(recvList_xy, recvCount_xy ,recvID_xy, id);
	UnpackID(recvList_Xy, recvCount_Xy ,recvID_Xy, id);
	UnpackID(recvList_xY, recvCount_xY ,recvID_xY, id);
	UnpackID(recvList_XY, recvCount_XY ,recvID_XY, id);
	UnpackID(recvList_xz, recvCount_xz ,recvID_xz, id);
	UnpackID(recvList_Xz, recvCount_Xz ,recvID_Xz, id);
	UnpackID(recvList_xZ, recvCount_xZ ,recvID_xZ, id);
	UnpackID(recvList_XZ, recvCount_XZ ,recvID_XZ, id);
	UnpackID(recvList_yz, recvCount_yz ,recvID_yz, id);
	UnpackID(recvList_Yz, recvCount_Yz ,recvID_Yz, id);
	UnpackID(recvList_yZ, recvCount_yZ ,recvID_yZ, id);
	UnpackID(recvList_YZ, recvCount_YZ ,recvID_YZ, id);
	//.....................................................................................
	// Once the ID is saved, free memory allocated to the buffers (no longer needed)
	//......................................................................................
	free(sendID_x); free(sendID_X); free(sendID_y); free(sendID_Y); free(sendID_z); free(sendID_Z);
	free(sendID_xy); free(sendID_XY); free(sendID_xY); free(sendID_Xy);
	free(sendID_xz); free(sendID_XZ); free(sendID_xZ); free(sendID_Xz);
	free(sendID_yz); free(sendID_YZ); free(sendID_yZ); free(sendID_Yz);
	free(recvID_x); free(recvID_X); free(recvID_y); free(recvID_Y); free(recvID_z); free(recvID_Z);
	free(recvID_xy); free(recvID_XY); free(recvID_xY); free(recvID_Xy);
	free(recvID_xz); free(recvID_XZ); free(recvID_xZ); free(recvID_Xz);
	free(recvID_yz); free(recvID_YZ); free(recvID_yZ); free(recvID_Yz);
	//......................................................................................
	if (rank==0)	printf ("Devices are ready to communicate. \n");
	MPI_Barrier(comm);

	//...........device phase ID.................................................
	if (rank==0)	printf ("Copying phase ID to device \n");
	char *ID;
	cudaMalloc((void **) &ID, N);						// Allocate device memory
	// Copy to the device
	cudaMemcpy(ID, id, N, cudaMemcpyHostToDevice);
	//...........................................................................

	if (rank==0)	printf ("Allocating distributions \n");
	//......................device distributions.................................
	double *f_even,*f_odd;
	//...........................................................................
	cudaMalloc((void **) &f_even, 10*dist_mem_size);	// Allocate device memory
	cudaMalloc((void **) &f_odd, 9*dist_mem_size);		// Allocate device memory
	//...........................................................................

	if (rank==0)	printf("Setting the distributions, size = : %i\n", N);
	//...........................................................................
	INITIALIZE <<< grid, nthreads >>>  (ID, f_even, f_odd, Nx, Ny, Nz, S);
	//...........................................................................

	//...........................................................................
	// Grids used to pack faces on the GPU for MPI
	int faceGrid,edgeGrid,packThreads;
	packThreads=512;
	edgeGrid=1;
	faceGrid=Nx*Ny/packThreads;
	//...........................................................................

	int iter = 0;
	if (rank==0)	printf("No. of iterations: %i \n", iterMax);
	
	//.......create a stream for the LB calculation.......
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	
	//.......create and start timer............
	double starttime,stoptime,cputime;
	MPI_Barrier(comm);
	starttime = MPI_Wtime();
	// Old cuda timer is below
//	cudaEvent_t start, stop;
//	float time;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);
//	cudaEventRecord( start, 0 );
	//.........................................

	sendtag = recvtag = 5;

	//************ MAIN ITERATION LOOP ***************************************/
	while (iter < iterMax){
		//...................................................................................
		PackDist<<<faceGrid,packThreads>>>(1,dvcSendList_x,0,sendCount_x,sendbuf_x,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(4,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(5,dvcSendList_x,2*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(6,dvcSendList_x,3*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(7,dvcSendList_x,4*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		//...Packing for X face<<<faceGrid,packThreads>>>(1,7,9,11,13)................................
		PackDist<<<faceGrid,packThreads>>>(0,dvcSendList_X,0,sendCount_X,sendbuf_X,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(3,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(4,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(5,dvcSendList_X,3*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(6,dvcSendList_X,4*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		//...Packing for y face<<<faceGrid,packThreads>>>(4,8,9,16,18).................................
		PackDist<<<faceGrid,packThreads>>>(2,dvcSendList_y,0,sendCount_y,sendbuf_y,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(4,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(4,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(8,dvcSendList_y,3*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(9,dvcSendList_y,4*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		//...Packing for Y face<<<faceGrid,packThreads>>>(3,7,10,15,17).................................
		PackDist<<<faceGrid,packThreads>>>(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(3,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(5,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(7,dvcSendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(8,dvcSendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		//...Packing for z face<<<faceGrid,packThreads>>>(6,12,13,16,17)................................
		PackDist<<<faceGrid,packThreads>>>(3,dvcSendList_z,0,sendCount_z,sendbuf_z,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(6,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(6,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(8,dvcSendList_z,3*sendCount_z,sendCount_z,sendbuf_z,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(8,dvcSendList_z,4*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
		//...Packing for Z face<<<faceGrid,packThreads>>>(5,11,14,15,18)................................
		PackDist<<<faceGrid,packThreads>>>(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(5,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(7,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
		PackDist<<<faceGrid,packThreads>>>(7,dvcSendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
		PackDist<<<faceGrid,packThreads>>>(9,dvcSendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
		//...Pack the xy edge <<<edgeGrid,packThreads>>>(8)................................
		PackDist<<<edgeGrid,packThreads>>>(4,dvcSendList_xy,0,sendCount_xy,sendbuf_xy,f_even,N);
		//...Pack the Xy edge <<<edgeGrid,packThreads>>>(9)................................
		PackDist<<<edgeGrid,packThreads>>>(4,dvcSendList_Xy,0,sendCount_Xy,sendbuf_Xy,f_odd,N);
		//...Pack the xY edge <<<edgeGrid,packThreads>>>(10)................................
		PackDist<<<edgeGrid,packThreads>>>(5,dvcSendList_xY,0,sendCount_xY,sendbuf_xY,f_even,N);
		//...Pack the XY edge <<<edgeGrid,packThreads>>>(7)................................
		PackDist<<<edgeGrid,packThreads>>>(3,dvcSendList_XY,0,sendCount_XY,sendbuf_XY,f_odd,N);
		//...Pack the xz edge <<<edgeGrid,packThreads>>>(12)................................
		PackDist<<<edgeGrid,packThreads>>>(6,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
		//...Pack the xZ edge <<<edgeGrid,packThreads>>>(14)................................
		PackDist<<<edgeGrid,packThreads>>>(7,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
		//...Pack the Xz edge <<<edgeGrid,packThreads>>>(13)................................
		PackDist<<<edgeGrid,packThreads>>>(6,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
		//...Pack the XZ edge <<<edgeGrid,packThreads>>>(11)................................
		PackDist<<<edgeGrid,packThreads>>>(5,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
		//...Pack the xz edge <<<edgeGrid,packThreads>>>(12)................................
		PackDist<<<edgeGrid,packThreads>>>(6,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
		//...Pack the xZ edge <<<edgeGrid,packThreads>>>(14)................................
		PackDist<<<edgeGrid,packThreads>>>(7,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
		//...Pack the Xz edge <<<edgeGrid,packThreads>>>(13)................................
		PackDist<<<edgeGrid,packThreads>>>(6,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
		//...Pack the XZ edge <<<edgeGrid,packThreads>>>(11)................................
		PackDist<<<edgeGrid,packThreads>>>(5,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
		//...Pack the yz edge <<<edgeGrid,packThreads>>>(16)................................
		PackDist<<<edgeGrid,packThreads>>>(8,dvcSendList_yz,0,sendCount_yz,sendbuf_yz,f_even,N);
		//...Pack the yZ edge <<<edgeGrid,packThreads>>>(18)................................
		PackDist<<<edgeGrid,packThreads>>>(9,dvcSendList_yZ,0,sendCount_yZ,sendbuf_yZ,f_even,N);
		//...Pack the Yz edge <<<edgeGrid,packThreads>>>(17)................................
		PackDist<<<edgeGrid,packThreads>>>(8,dvcSendList_Yz,0,sendCount_Yz,sendbuf_Yz,f_odd,N);
		//...Pack the YZ edge <<<edgeGrid,packThreads>>>(15)................................
		PackDist<<<edgeGrid,packThreads>>>(7,dvcSendList_YZ,0,sendCount_YZ,sendbuf_YZ,f_odd,N);
		//...................................................................................

		//*****************************************************************************
		//*****************************************************************************
		//*****************************************************************************
		//........ Execute the swap kernel (device) .........................
		//*****************************************************************************
		//*****************************************************************************
		SWAP <<< grid, nthreads >>> (ID, f_even, f_odd, Nx, Ny, Nz, S);
		//*****************************************************************************
		//*****************************************************************************
		//*****************************************************************************
		//*****************************************************************************

		//...................................................................................
		// Send all the distributions
		MPI_Isend(sendbuf_x, 5*sendCount_x,MPI_DOUBLE,rank_X,sendtag,comm,&req1[0]);
		MPI_Irecv(recvbuf_X, 5*recvCount_X,MPI_DOUBLE,rank_x,recvtag,comm,&req2[0]);
		MPI_Isend(sendbuf_X, 5*sendCount_X,MPI_DOUBLE,rank_x,sendtag,comm,&req1[1]);
		MPI_Irecv(recvbuf_x, 5*recvCount_x,MPI_DOUBLE,rank_X,recvtag,comm,&req2[1]);
		MPI_Isend(sendbuf_y, 5*sendCount_y,MPI_DOUBLE,rank_Y,sendtag,comm,&req1[2]);
		MPI_Irecv(recvbuf_Y, 5*recvCount_Y,MPI_DOUBLE,rank_y,recvtag,comm,&req2[2]);
		MPI_Isend(sendbuf_Y, 5*sendCount_Y,MPI_DOUBLE,rank_y,sendtag,comm,&req1[3]);
		MPI_Irecv(recvbuf_y, 5*recvCount_y,MPI_DOUBLE,rank_Y,recvtag,comm,&req2[3]);
		MPI_Isend(sendbuf_z, 5*sendCount_z,MPI_DOUBLE,rank_Z,sendtag,comm,&req1[4]);
		MPI_Irecv(recvbuf_Z, 5*recvCount_Z,MPI_DOUBLE,rank_z,recvtag,comm,&req2[4]);
		MPI_Isend(sendbuf_Z, 5*sendCount_Z,MPI_DOUBLE,rank_z,sendtag,comm,&req1[5]);
		MPI_Irecv(recvbuf_z, 5*recvCount_z,MPI_DOUBLE,rank_Z,recvtag,comm,&req2[5]);
		MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_XY,sendtag,comm,&req1[6]);
		MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_xy,recvtag,comm,&req2[6]);
		MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_xy,sendtag,comm,&req1[7]);
		MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_XY,recvtag,comm,&req2[7]);
		MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_xY,sendtag,comm,&req1[8]);
		MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_Xy,recvtag,comm,&req2[8]);
		MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_Xy,sendtag,comm,&req1[9]);
		MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_xY,recvtag,comm,&req2[9]);
		MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_XZ,sendtag,comm,&req1[10]);
		MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_xz,recvtag,comm,&req2[10]);
		MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_xz,sendtag,comm,&req1[11]);
		MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_XZ,recvtag,comm,&req2[11]);
		MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_xZ,sendtag,comm,&req1[12]);
		MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_Xz,recvtag,comm,&req2[12]);
		MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_Xz,sendtag,comm,&req1[13]);
		MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_xZ,recvtag,comm,&req2[13]);
		MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_YZ,sendtag,comm,&req1[14]);
		MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_yz,recvtag,comm,&req2[14]);
		MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_yz,sendtag,comm,&req1[15]);
		MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_YZ,recvtag,comm,&req2[15]);
		MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_yZ,sendtag,comm,&req1[16]);
		MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_Yz,recvtag,comm,&req2[16]);
		MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_Yz,sendtag,comm,&req1[17]);
		MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_yZ,recvtag,comm,&req2[17]);
		//...................................................................................

		//...................................................................................
		// Wait for completion of D3Q19 communication
		MPI_Waitall(18,req1,stat1);
		MPI_Waitall(18,req2,stat2);
		//...................................................................................
		// Unpack the distributions on the device
		//...................................................................................
		//...Map recieve list for the X face: q=2,8,10,12,13 .................................
		MapRecvDist<<<faceGrid,packThreads>>>(0,-1,0,0,dvcRecvList_X,0,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(3,-1,-1,0,dvcRecvList_X,recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(4,-1,1,0,dvcRecvList_X,2*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(5,-1,0,-1,dvcRecvList_X,3*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(6,-1,0,1,dvcRecvList_X,4*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the x face: q=1,7,9,11,13..................................
		MapRecvDist<<<faceGrid,packThreads>>>(1,1,0,0,dvcRecvList_x,0,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(4,1,1,0,dvcRecvList_x,recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(5,1,-1,0,dvcRecvList_x,2*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(6,1,0,1,dvcRecvList_x,3*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(7,1,0,-1,dvcRecvList_x,4*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
		MapRecvDist<<<faceGrid,packThreads>>>(1,0,-1,0,dvcRecvList_Y,0,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(3,-1,-1,0,dvcRecvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(5,1,-1,0,dvcRecvList_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(7,0,-1,-1,dvcRecvList_Y,3*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(8,0,-1,1,dvcRecvList_Y,4*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
		MapRecvDist<<<faceGrid,packThreads>>>(2,0,1,0,dvcRecvList_y,0,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(4,1,1,0,dvcRecvList_y,recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(4,-1,1,0,dvcRecvList_y,2*recvCount_y,recvCount_y,recvbuf_y,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(8,0,1,1,dvcRecvList_y,3*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(9,0,1,-1,dvcRecvList_y,4*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the z face<<<faceGrid,packThreads>>>(6,12,13,16,17)..............................................
		MapRecvDist<<<faceGrid,packThreads>>>(2,0,0,-1,dvcRecvList_Z,0,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(5,-1,0,-1,dvcRecvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(7,1,0,-1,dvcRecvList_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(7,0,-1,-1,dvcRecvList_Z,3*recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(9,0,1,-1,dvcRecvList_Z,4*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
		//...Map recieve list for the Z face<<<faceGrid,packThreads>>>(5,11,14,15,18)..............................................
		MapRecvDist<<<faceGrid,packThreads>>>(3,0,0,1,dvcRecvList_z,0,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(6,1,0,1,dvcRecvList_z,recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(6,-1,0,1,dvcRecvList_z,2*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(8,0,1,1,dvcRecvList_z,3*recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		MapRecvDist<<<faceGrid,packThreads>>>(8,0,-1,1,dvcRecvList_z,4*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
		//..................................................................................
		//...Map recieve list for the xy edge <<<edgeGrid,packThreads>>>(8)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(3,-1,-1,0,dvcRecvList_XY,0,recvCount_XY,recvbuf_XY,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Xy edge <<<edgeGrid,packThreads>>>(9)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(5,1,-1,0,dvcRecvList_xY,0,recvCount_xY,recvbuf_xY,f_even,Nx,Ny,Nz);
		//...Map recieve list for the xY edge <<<edgeGrid,packThreads>>>(10)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(4,-1,1,0,dvcRecvList_Xy,0,recvCount_Xy,recvbuf_Xy,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the XY edge <<<edgeGrid,packThreads>>>(7)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(4,1,1,0,dvcRecvList_xy,0,recvCount_xy,recvbuf_xy,f_even,Nx,Ny,Nz);
		//...Map recieve list for the xz edge <<<edgeGrid,packThreads>>>(12)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(5,-1,0,-1,dvcRecvList_XZ,0,recvCount_XZ,recvbuf_XZ,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the xZ edge <<<edgeGrid,packThreads>>>(14)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(6,-1,0,1,dvcRecvList_Xz,0,recvCount_Xz,recvbuf_Xz,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Xz edge <<<edgeGrid,packThreads>>>(13)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(7,1,0,-1,dvcRecvList_xZ,0,recvCount_xZ,recvbuf_xZ,f_even,Nx,Ny,Nz);
		//...Map recieve list for the XZ edge <<<edgeGrid,packThreads>>>(11)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(6,1,0,1,dvcRecvList_xz,0,recvCount_xz,recvbuf_xz,f_even,Nx,Ny,Nz);
		//...Map recieve list for the yz edge <<<edgeGrid,packThreads>>>(16)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(7,0,-1,-1,dvcRecvList_YZ,0,recvCount_YZ,recvbuf_YZ,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the yZ edge <<<edgeGrid,packThreads>>>(18)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(8,0,-1,1,dvcRecvList_Yz,0,recvCount_Yz,recvbuf_Yz,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Yz edge <<<edgeGrid,packThreads>>>(17)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(9,0,1,-1,dvcRecvList_yZ,0,recvCount_yZ,recvbuf_yZ,f_even,Nx,Ny,Nz);
		//...Map recieve list for the YZ edge <<<edgeGrid,packThreads>>>(15)................................
		MapRecvDist<<<edgeGrid,packThreads>>>(8,0,1,1,dvcRecvList_yz,0,recvCount_yz,recvbuf_yz,f_even,Nx,Ny,Nz);
		//...................................................................................

		//*****************************************************************************
		//*****************************************************************************
		//*****************************************************************************
		//........ Execute the collision kernel (device) ....................
		//*****************************************************************************
		//*****************************************************************************
		MRT <<< grid, nthreads >>> (ID, f_even, f_odd, Nx, Ny, Nz, S,
									rlx_setA, rlx_setB, Fx, Fy, Fz);
		//*****************************************************************************
		//*****************************************************************************
		//*****************************************************************************
		//*****************************************************************************

		MPI_Barrier(comm);
		// Iteration completed!
		iter++;
		//...................................................................
	}
	//************************************************************************/
	
	cudaThreadSynchronize();
	MPI_Barrier(comm);
	stoptime = MPI_Wtime();
//	cout << "CPU time: " << (stoptime - starttime) << " seconds" << endl;
	cputime = stoptime - starttime;
//	cout << "Lattice update rate: "<< double(Nx*Ny*Nz*iter)/cputime/1000000 <<  " MLUPS" << endl;
	double MLUPS = double(Nx*Ny*Nz*iter)/cputime/1000000;

	if (rank==0) printf("CPU time = %f \n", cputime);
	if (rank==0) printf("Lattice update rate = %f MLUPS \n", MLUPS);
	//.......... stop and destroy timer.............................
//	cudaEventRecord( stop, stream);
//	cudaEventSynchronize( stop );
//	cudaEventElapsedTime( &time, start, stop );
//	printf("CPU time = %f \n", time);
//	float MLUPS = 0.001*float(Nx*Ny*Nz)*iter/time;
//	printf("MLUPS = %f \n", MLUPS);

	cudaStreamDestroy(stream);
//	cudaEventDestroy( start );
//	cudaEventDestroy( stop );
	//..............................................................
	
	//..............................................................
	//.........Compute the velocity and copy result to host ........
	double *velocity;
	velocity = new double[3*N];
	//......................device distributions....................................
	double *vel;
	//..............................................................................
	cudaMalloc((void **) &vel, 3*dist_mem_size);	// Allocate device memory
	//..............................................................................
	Compute_VELOCITY <<< grid, nthreads >>>  (ID, f_even, f_odd, vel, Nx, Ny, Nz, S);
	//..............................................................................
	cudaMemcpy(velocity, vel, 3*dist_mem_size, cudaMemcpyDeviceToHost);
	//..............................................................................
	cudaThreadSynchronize();
	MPI_Barrier(comm);
	//............................................................	
	//....Write the z-velocity to test poiseuille flow............
	double vz,vz_avg;	
	vz_avg = 0.0;

	FILE *output;
	output = fopen("velocity.out","w");
	for (int k=0; k<1; k++){
		for (int j=0; j<1; j++){
			for (int i=0; i<Nx; i++){
				int n = k*Nx*Ny+j*Nx+i;
				//.....print value........
				vz = velocity[2*N+n];
				vz_avg += vz;
				fprintf(output, " %e",vz);
			}
		}
	}
	fclose(output);
	
	vz = vz_avg/double(sum);
	printf("Average Velocity = %e\n", vz);

	// cleanup	
	cudaFree(f_even);	cudaFree(f_odd);	cudaFree(vel);	cudaFree(ID);
	free (velocity);	free(id);
	
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
