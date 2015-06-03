// Compute the signed distance from a digitized image 
// Two phases are present
// Phase 1 has value -1
// Phase 2 has value 1
// this code uses the segmented image to generate the signed distance 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <Array.h>

inline void SSO(DoubleArray &Distance, char *ID, int timesteps){

	int Q=26;
	int q,i,j,k,n;
	int Nx = Distance.m;
	int Ny = Distance.n;
	int Nz = Distance.o;
	const static int D3Q27[26][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
	{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},
	{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1},{1,1,1},{-1,-1,-1},{1,1,-1},{-1,-1,1},
	{-1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1}};

	double weights[26];
	// Compute the weights from the finite differences
	for (q=0; q<Q; q++){
		weights[q] = sqrt(1.0*(D3Q27[q][0]*D3Q27[q][0]) + 1.0*(D3Q27[q][1]*D3Q27[q][1]) + 1.0*(D3Q27[q][2]*D3Q27[q][2]));
	}

	// Initialize the Distance from ID
//	for (i=0; i<Nx*Ny*Nz; i++) Distance.data[i] = -0.5;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n=k*Nx*Ny+j*Nx+i;
				// Initialize distance to +/- 1
				Distance(i,j,k) = 1.0*ID[n]-0.5;
			}
		}
	}

	int count = 0;
	double dt=0.1;
	int in,jn,kn,nn;
	double Dqx,Dqy,Dqz,Dx,Dy,Dz,W;
	double nx,ny,nz,Cqx,Cqy,Cqz,sign,norm;
	double f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;

	printf("Number of timesteps is %i \n",timesteps);
	printf("Mesh is %i,%i,%i \n",Nx,Ny,Nz);

	while (count < timesteps){

		printf("count=%i \n",count);
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){

					n = k*Nx*Ny + j*Nx + i;
					sign = Distance.data[n] / fabs(Distance.data[n]);

					//............Compute the Gradient...................................
					if (!(i+1<Nx)) 	nx=0.5*Distance(i,j,k);
					else 			nx=0.5*Distance(i+1,j,k);;
					if (!(j+1<Ny)) 	ny=0.5*Distance(i,j,k);
					else 			ny=0.5*Distance(i,j+1,k);
					if (!(k+1<Nz)) 	nz=0.5*Distance(i,j,k);
					else 			nz=0.5*Distance(i,j,k+1);
					if (i<1)  		nx-=0.5*Distance(i,j,k);
					else 			nx-=0.5*Distance(i-1,j,k);
					if (j<1) 		ny-=0.5*Distance(i,j,k);
					else 			ny-=0.5*Distance(i,j-1,k);
					if (k<1)  		nz-=0.5*Distance(i,j,k);
					else 			nz-=0.5*Distance(i,j,k-1);

//					nx = 0.5*(Distance(i+1,j,k) - Distance(i-1,j,k));
	//				ny = 0.5*(Distance(i,j+1,k) - Distance(i,j-1,k));
		//			nz = 0.5*(Distance(i,j,k+1) - Distance(i,j,k-1));

					W = 0.0;	Dx = Dy = Dz = 0.0;
					if (nx*nx+ny*ny+nz*nz > 0.0){
						for (q=0; q<26; q++){
							Cqx = 1.0*D3Q27[q][0];
							Cqy = 1.0*D3Q27[q][1];
							Cqz = 1.0*D3Q27[q][2];
							// get the associated neighbor
							in = i + D3Q27[q][0];
							jn = j + D3Q27[q][1];
							kn = k + D3Q27[q][2];
							// make sure the neighbor is in the domain (periodic BC)
			/*				if (in < 0 ) in +=Nx;
							if (jn < 0 ) jn +=Ny;
							if (kn < 0 ) kn +=Nz;
							if (!(in < Nx) ) in -=Nx;
							if (!(jn < Ny) ) jn -=Ny;
							if (!(kn < Nz) ) kn -=Nz;
			*/				// symmetric boundary
							if (in < 0 ) in = i;
							if (jn < 0 ) jn = j;
							if (kn < 0 ) kn = k;
							if (!(in < Nx) ) in = i;
							if (!(jn < Ny) ) jn = k;
							if (!(kn < Nz) ) kn = k;
							// 1-D index
							nn = kn*Nx*Ny + jn*Nx + in;

							// Compute the gradient using upwind finite differences
							Dqx = weights[q]*(Distance.data[n] - Distance.data[nn])*Cqx;
							Dqy = weights[q]*(Distance.data[n] - Distance.data[nn])*Cqy;
							Dqz = weights[q]*(Distance.data[n] - Distance.data[nn])*Cqz;

							// Only include upwind derivatives
							if (sign*(nx*Cqx + ny*Cqy + nz*Cqz) < 0.0 ){

								Dx += Dqx;
								Dy += Dqy;
								Dz += Dqz;
								W += weights[q];
							}
						}
						// Normalize by the weight to get the approximation to the gradient
						Dx /= W;
						Dy /= W;
						Dz /= W;

						norm = sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
					}
					else{
						norm = 0.0;
					}
					Distance.data[n] += dt*sign*(1.0 - norm);

					// Disallow any change in phase
					if (Distance.data[n]*2.0*(ID[n]-1.0) < 0) Distance.data[n] = -Distance.data[n];
				}
			}
		}

		count++;
	}
}

int main(){
	
	int i,j,k,n,nn;
	int Nx, Ny, Nz, N;
	Nx = Ny = Nz = 50;
	N = Nx*Ny*Nz;
	
	double err = 1.0;
	double err_prev=1.0;
	double tol = 1e-6;
	double f_x, f_y, f_z, norm;
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	double fxm,fym,fzm,fxp,fyp,fzp;	
	double fxy,fXy,fxY,fXY,fxz,fXz,fxZ,fXZ,fyz,fYz,fyZ,fYZ;
	double nx,ny,nz;
	
	int ip,im,jp,jm,kp,km;
	
	int count = 0;
	double sign;
	double dt=1.0;
	printf("Nx=%i, Ny=%i, Nz= %i, \n",Nx,Ny,Nz);
	
	char *id;

#ifdef READMEDIA
	Nx = 347;
	Ny = 347;
	Nz = 235;
	Nx = 512;
	Ny = 512;
	Nz = 512;
	N = Nx*Ny*Nz;
	id = new char [N];

	FILE *INPUT = fopen("Solid.dat","rb");
	fread(id,1,Nx*Ny*Nz,INPUT);
	fclose(INPUT);
#else
	id = new char [N];
	double BubbleRadius = 5;
	// Initialize the bubble
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;
				// Initialize phase positions field
				if ((i-0.5*Nx)*(i-0.5*Nx)+(j-0.5*Ny)*(j-0.5*Ny)+(k-0.5*Nz)*(k-0.5*Nz) < BubbleRadius*BubbleRadius){
					id[n] = 0;
				}
				else{
					id[n]=1;
				}
			}
		}
	}
#endif
	
	DoubleArray Distance(Nx,Ny,Nz);
	// Initialize the signed distance function
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n=k*Nx*Ny+j*Nx+i;
				// Initialize distance to +/- 1
				Distance(i,j,k) = 2.0*ID[n]-1.0;
			}
		}
	}

	printf("Initialized! Converting to Signed Distance function \n");
	SSO(Distance,id,20);

/*	double *f,*f_old,*f_new;
	f = new double[N];
	f_old = new double[N];
	f_new = new double[N];

	for (int n=0; n<N; n++)	Distance.data[n] = 0.5*(id[n]-1);
	for (int n=0; n<N; n++)	f_old[n] = Distance.data[n];
	for (int n=0; n<N; n++)	f_new[n] = Distance.data[n];
	for (int n=0; n<N; n++)	f[n] = Distance.data[n];

	count=0;
	dt=1.0;
	while (count < 10 && dt > 1.0e-6){

		err = 0.0;
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){

					n = k*Nx*Ny + j*Nx + i;

					ip = i+1;
					im = i-1;
					jp = j+1;
					jm = j-1;
					kp = k+1;
					km = k-1;
					if (!(ip<Nx))	ip-=Nx;
					if (im < 0)		im+=Nx;
					if (!(jp<Ny))	jp-=Ny;
					if (jm < 0)		jm+=Ny;
					if (!(kp<Nz))	kp-=Nz;
					if (km < 0)		km+=Nz;

					//........................................................................
					//					COMPUTE THE COLOR GRADIENT
					//........................................................................
					//.................Read Phase Indicator Values............................
					//........................................................................
					nn = n-1;							// neighbor index (get convention)
					if (i-1<0)		nn += Nx;			// periodic BC along the x-boundary
					f1 = f[nn];						// get neighbor for phi - 1
					//........................................................................
					nn = n+1;							// neighbor index (get convention)
					if (!(i+1<Nx))	nn -= Nx;			// periodic BC along the x-boundary
					f2 = f[nn];						// get neighbor for phi - 2
					//........................................................................
					nn = n-Nx;							// neighbor index (get convention)
					if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
					f3 = f[nn];					// get neighbor for phi - 3
					//........................................................................
					nn = n+Nx;							// neighbor index (get convention)
					if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
					f4 = f[nn];					// get neighbor for phi - 4
					//........................................................................
					nn = n-Nx*Ny;						// neighbor index (get convention)
					if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
					f5 = f[nn];					// get neighbor for phi - 5
					//........................................................................
					nn = n+Nx*Ny;						// neighbor index (get convention)
					if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
					f6 = f[nn];					// get neighbor for phi - 6
					//........................................................................
					nn = n-Nx-1;						// neighbor index (get convention)
					if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
					if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
					f7 = f[nn];					// get neighbor for phi - 7
					//........................................................................
					nn = n+Nx+1;						// neighbor index (get convention)
					if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
					if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
					f8 = f[nn];					// get neighbor for phi - 8
					//........................................................................
					nn = n+Nx-1;						// neighbor index (get convention)
					if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
					if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
					f9 = f[nn];					// get neighbor for phi - 9
					//........................................................................
					nn = n-Nx+1;						// neighbor index (get convention)
					if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
					if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
					f10 = f[nn];					// get neighbor for phi - 10
					//........................................................................
					nn = n-Nx*Ny-1;						// neighbor index (get convention)
					if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
					if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
					f11 = f[nn];					// get neighbor for phi - 11
					//........................................................................
					nn = n+Nx*Ny+1;						// neighbor index (get convention)
					if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
					if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
					f12 = f[nn];					// get neighbor for phi - 12
					//........................................................................
					nn = n+Nx*Ny-1;						// neighbor index (get convention)
					if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
					if (!(k+1<Nz))		nn -= Nx*Ny*Nz;	// Perioidic BC along the z-boundary
					f13 = f[nn];					// get neighbor for phi - 13
					//........................................................................
					nn = n-Nx*Ny+1;						// neighbor index (get convention)
					if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
					if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
					f14 = f[nn];					// get neighbor for phi - 14
					//........................................................................
					nn = n-Nx*Ny-Nx;					// neighbor index (get convention)
					if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
					if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
					f15 = f[nn];					// get neighbor for phi - 15
					//........................................................................
					nn = n+Nx*Ny+Nx;					// neighbor index (get convention)
					if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
					if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
					f16 = f[nn];					// get neighbor for phi - 16
					//........................................................................
					nn = n+Nx*Ny-Nx;					// neighbor index (get convention)
					if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
					if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
					f17 = f[nn];					// get neighbor for phi - 17
					//........................................................................
					nn = n-Nx*Ny+Nx;					// neighbor index (get convention)
					if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
					if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
					f18 = f[nn];					// get neighbor for phi - 18
					//............Compute the Color Gradient...................................
					nx = -(f1-f2+0.5*(f7-f8+f9-f10+f11-f12+f13-f14));
					ny = -(f3-f4+0.5*(f7-f8-f9+f10+f15-f16+f17-f18));
					nz = -(f5-f6+0.5*(f11-f12-f13+f14+f15-f16-f17+f18));
					//...........Normalize the Color Gradient.................................
					
					f_x = 0.5*(f[k*Nx*Ny + j*Nx + ip] - f[k*Nx*Ny + j*Nx + im]);
					f_y = 0.5*(f[k*Nx*Ny + jp*Nx + i] - f[k*Nx*Ny + jm*Nx + i]);
					f_z = 0.5*(f[kp*Nx*Ny + j*Nx + i] - f[km*Nx*Ny + j*Nx + i]);					

					if (id[n] < 0){
						if ( nx > 0.0)	f_x = fxp;
						else			f_x = fxm;
						if ( ny > 0.0)	f_y = fyp;
						else			f_y = fym;
						if ( nz > 0.0)	f_z = fzp;
						else			f_z = fzm;
					}
					else{
						if ( nx > 0.0)	f_x = fxm;
						else			f_x = fxp;
						if ( ny > 0.0)	f_y = fym;
						else			f_y = fyp;
						if ( nz > 0.0)	f_z = fzm;
						else			f_z = fzp;
					}
										
					norm = sqrt(f_x*f_x+f_y*f_y+f_z*f_z);

					if (err < (1.0 - norm))	err = (1.0 - norm);

					if (fabs(1.0-norm) > err) 	err =  fabs(1.0-norm);

					sign =1.0;
					if (!(id[n] > 0)) sign = -1.0;

					f_new[n] = f_old[n] + dt*sign*(1.0 - norm);

				}
			}
		}

		for (int n=0; n<N; n++)	f[n] = f_new[n];
		for (int n=0; n<N; n++)	f_old[n] = f[n];

		printf("Error %i = %f \n",count,err);

		if (err > err_prev)	dt *= 0.2;
		err_prev = err;
		
		count++;			
	}
*/

	FILE *DIST;
	DIST = fopen("SignDist","wb");
	fwrite(Distance.data,8,N,DIST);
	fclose(DIST);
	
}
