#include <math.h>

#define STOKES

extern "C" void InitDenColor(char *ID, double *Den, double *Phi, double das, double dbs, int Nx, int Ny, int Nz)
{
	int i,j,k,n,N;

	N = Nx*Ny*Nz;

	for (n=0; n<N; n++){

		//.......Back out the 3-D indices for node n..............
		k = n/(Nx*Ny);
		j = (n-Nx*Ny*k)/Nx;
		i = n-Nx*Ny*k-Nx*j;

		if ( ID[n] == 1){
			Den[n] = 1.0;
			Den[N+n] = 0.0;
			Phi[n] = 1.0;
		}
		else if ( ID[n] == 2){
			Den[n] = 0.0;
			Den[N+n] = 1.0;
			Phi[n] = -1.0;
		}
		else{
			Den[n] = das;
			Den[N+n] = dbs;
			Phi[n] = (das-dbs)/(das+dbs);
		}
	}
}
extern "C" void InitDenColorDistancePacked(char *ID, double *Den, double *Phi, double *Distance,
								double das, double dbs, double beta, double xp, int Nx, int Ny, int Nz)
{
	int i,j,k,n,N;
	double d;

	N = Nx*Ny*Nz;

	for (n=0; n<N; n++){

		//.......Back out the 3-D indices for node n..............
		k = n/(Nx*Ny);
		j = (n-Nx*Ny*k)/Nx;
		i = n-Nx*Ny*k-Nx*j;

		if ( ID[n] == 1){
			Den[2*n] = 1.0;
			Den[2*n+1] = 0.0;
			Phi[n] = 1.0;
		}
		if (i == 0 || j == 0 || k == 0 || i == Nx-1 || j == Ny-1 || k == Nz-1){
			Den[2*n] = 0.0;
			Den[2*n+1] = 0.0;
		}
		else if ( ID[n] == 1){
			Den[2*n] = 1.0;
			Den[2*n+1] = 0.0;
			Phi[n] = 1.0;
		}
		else if ( ID[n] == 2){
			Den[2*n] = 0.0;
			Den[2*n+1] = 1.0;
			Phi[n] = -1.0;
		}
		else{
			Den[2*n] = das;
			Den[2*n+1] = dbs;
			Phi[n] = (das-dbs)/(das+dbs);
			d = fabs(Distance[n]);
			Phi[n] = (2.f*(exp(-2.f*beta*(d+xp)))/(1.f+exp(-2.f*beta*(d+xp))) - 1.f);
		}
	}
}

extern "C" void InitDenColorDistance(char *ID, double *Den, double *Phi, double *Distance,
								double das, double dbs, double beta, double xp, int Nx, int Ny, int Nz)
{
	int i,j,k,n,N;
	double d;

	N = Nx*Ny*Nz;

	for (n=0; n<N; n++){

		//.......Back out the 3-D indices for node n..............
		k = n/(Nx*Ny);
		j = (n-Nx*Ny*k)/Nx;
		i = n-Nx*Ny*k-Nx*j;

		if ( ID[n] == 1){
			Den[n] = 1.0;
			Den[N+n] = 0.0;
			Phi[n] = 1.0;
		}
		else if ( ID[n] == 2){
			Den[n] = 0.0;
			Den[N+n] = 1.0;
			Phi[n] = -1.0;
		}
		else{
			Den[n] = das;
			Den[N+n] = dbs;
			Phi[n] = (das-dbs)/(das+dbs);
			d = fabs(Distance[n]);
			Phi[n] = (2.f*(exp(-2.f*beta*(d+xp)))/(1.f+exp(-2.f*beta*(d+xp))) - 1.f);
		}
	}
}



//*************************************************************************
extern "C" void ColorBC_inlet(double *Phi, double *Den, double *A_even, double *A_odd,
								  double *B_even, double *B_odd, int Nx, int Ny, int Nz)
{
	int i,j,k,n,N;
	N = Nx*Ny*Nz;
	// Fill the inlet with component a
	for (k=0; k<1; k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				Phi[n] = 1.0;
			}					
		}
	}
	
	for (k=1; k<3; k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				Phi[n] = 1.0;
				Den[n] = 1.0;
				Den[N+n] = 0.0;

				A_even[n] = 0.3333333333333333;
				A_odd[n] = 0.1111111111111111;	
				A_even[N+n] = 0.1111111111111111;
				A_odd[N+n] = 0.1111111111111111;	
				A_even[2*N+n] = 0.1111111111111111;
				A_odd[2*N+n] = 0.1111111111111111;
				A_even[3*N+n] = 0.1111111111111111;
				
				B_even[n] = 0.0;
				B_odd[n] = 0.0;	
				B_even[N+n] = 0.0;
				B_odd[N+n] = 0.0;	
				B_even[2*N+n] = 0.0;
				B_odd[2*N+n] = 0.0;
				B_even[3*N+n] = 0.0;
			}					
		}
	}
}
//*************************************************************************
extern "C" void ColorBC_outlet(double *Phi, double *Den, double *A_even, double *A_odd,
								  double *B_even, double *B_odd, int Nx, int Ny, int Nz)
{
	int i,j,k,n,N;
	N = Nx*Ny*Nz;
	// Fill the outlet with component b
	for (k=Nz-3; k<Nz-1; k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){

				n = k*Nx*Ny+j*Nx+i;
				Phi[n] = -1.0;
				Den[n] = 0.0;
				Den[N+n] = 1.0;
				
				A_even[n] = 0.0;
				A_odd[n] = 0.0;	
				A_even[N+n] = 0.0;
				A_odd[N+n] = 0.0;	
				A_even[2*N+n] = 0.0;
				A_odd[2*N+n] = 0.0;
				A_even[3*N+n] = 0.0;
				
				B_even[n] = 0.3333333333333333;
				B_odd[n] = 0.1111111111111111;	
				B_even[N+n] = 0.1111111111111111;
				B_odd[N+n] = 0.1111111111111111;	
				B_even[2*N+n] = 0.1111111111111111;
				B_odd[2*N+n] = 0.1111111111111111;
				B_even[3*N+n] = 0.1111111111111111;
				
			}					
		}
	}
	for (k=Nz-1; k<Nz; k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				Phi[n] = -1.0;
			}					
		}
	}
}
//*************************************************************************
extern "C" void PressureBC_inlet(double *disteven, double *distodd, double din,
								  int Nx, int Ny, int Nz)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double uz;

	N = Nx*Ny*Nz;

	for (n=Nx*Ny; n<2*Nx*Ny; n++){

		//........................................................................
		// Read distributions from "opposite" memory convention
		//........................................................................
		//........................................................................
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
		f0 = disteven[n];
		f2 = disteven[N+n];
		f4 = disteven[2*N+n];
		f6 = disteven[3*N+n];
		f8 = disteven[4*N+n];
		f10 = disteven[5*N+n];
		f12 = disteven[6*N+n];
		f14 = disteven[7*N+n];
		f16 = disteven[8*N+n];
		f18 = disteven[9*N+n];
		//...................................................
		//........Determine the intlet flow velocity.........
		//			uz = -1 + (f0+f3+f4+f1+f2+f7+f8+f10+f9
		//					   + 2*(f5+f15+f18+f11+f14))/din;
		//........Set the unknown distributions..............
		//			f6 = f5 - 0.3333333333333333*din*uz;
		//			f16 = f15 - 0.1666666666666667*din*uz;
		//			f17 = f16 - f3 + f4-f15+f18-f7+f8-f10+f9;
		//			f12= 0.5*(-din*uz+f5+f15+f18+f11+f14-f6-f16-
		//					  f17+f1-f2-f14+f11+f7-f8-f10+f9);
		//			f13= -din*uz+f5+f15+f18+f11+f14-f6-f16-f17-f12;

		// Determine the outlet flow velocity
		uz = 1.0 - (f0+f4+f3+f2+f1+f8+f7+f9+ f10 +
				2*(f5+ f15+f18+f11+f14))/din;
		// Set the unknown distributions:
		f6 = f5 + 0.3333333333333333*din*uz;
		f16 = f15 + 0.1666666666666667*din*uz;
		f17 = f16 + f4 - f3-f15+f18+f8-f7	+f9-f10;
		f12= (din*uz+f5+ f15+f18+f11+f14-f6-f16-f17-f2+f1-f14+f11-f8+f7+f9-f10)*0.5;
		f13= din*uz+f5+ f15+f18+f11+f14-f6-f16-f17-f12;

		//........Store in "opposite" memory location..........
		disteven[3*N+n] = f6;
		disteven[6*N+n] = f12;
		distodd[6*N+n] = f13;
		disteven[8*N+n] = f16;
		distodd[8*N+n] = f17;
		//...................................................
	}
}

extern "C" void PressureBC_outlet(double *disteven, double *distodd, double dout,
								   int Nx, int Ny, int Nz, int outlet)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double uz;

	N = Nx*Ny*Nz;

	// Loop over the boundary - threadblocks delineated by start...finish
	for (n=outlet; n<N-Nx*Ny; n++){

		//........................................................................
		// Read distributions from "opposite" memory convention
		//........................................................................
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
		f0 = disteven[n];
		f2 = disteven[N+n];
		f4 = disteven[2*N+n];
		f6 = disteven[3*N+n];
		f8 = disteven[4*N+n];
		f10 = disteven[5*N+n];
		f12 = disteven[6*N+n];
		f14 = disteven[7*N+n];
		f16 = disteven[8*N+n];
		f18 = disteven[9*N+n];
		//........Determine the outlet flow velocity.........
		//			uz = 1 - (f0+f3+f4+f1+f2+f7+f8+f10+f9+
		//					  2*(f6+f16+f17+f12+f13))/dout;
		//...................................................
		//........Set the Unknown Distributions..............
		//			f5 = f6 + 0.33333333333333338*dout*uz;
		//			f15 = f16 + 0.16666666666666678*dout*uz;
		//			f18 = f15+f3-f4-f16+f17+f7-f8+f10-f9;
		//			f11= 0.5*(dout*uz+f6+ f16+f17+f12+f13-f5
		//				  -f15-f18-f1+f2-f13+f12-f7+f8+f10-f9);
		//			f14= dout*uz+f6+ f16+f17+f12+f13-f5-f15-f18-f11;

		uz = -1.0 + (f0+f4+f3+f2+f1+f8+f7+f9+f10 + 2*(f6+f16+f17+f12+f13))/dout;

		f5 = f6 - 0.33333333333333338*dout* uz;
		f15 = f16 - 0.16666666666666678*dout* uz;
		f18 = f15 - f4 + f3-f16+f17-f8+f7-f9+f10;
		f11 = (-dout*uz+f6+ f16+f17+f12+f13-f5-f15-f18+f2-f1-f13+f12+f8-f7-f9+f10)*0.5;
		f14 = -dout*uz+f6+ f16+f17+f12+f13-f5-f15-f18-f11;
		//........Store in "opposite" memory location..........
		distodd[2*N+n] = f5;
		distodd[5*N+n] = f11;
		disteven[7*N+n] = f14;
		distodd[7*N+n] = f15;
		disteven[9*N+n] = f18;
		//...................................................

	}
}
//*************************************************************************
extern "C" void ComputeColorGradient(char *ID, double *phi, double *ColorGrad, int Nx, int Ny, int Nz)
{
	int n,N,i,j,k,nn;
	// distributions
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double nx,ny,nz;

	// non-conserved moments
	// additional variables needed for computations

	N = Nx*Ny*Nz;

	for (n=0; n<N; n++){

		//.......Back out the 3-D indices for node n..............
		k = n/(Nx*Ny);
		j = (n-Nx*Ny*k)/Nx;
		i = n-Nx*Ny*k-Nx*j;
		//........................................................................
		//........Get 1-D index for this thread....................
		//		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
		//........................................................................
		//					COMPUTE THE COLOR GRADIENT
		//........................................................................
		//.................Read Phase Indicator Values............................
		//........................................................................
		nn = n-1;							// neighbor index (get convention)
		if (i-1<0)		nn += Nx;			// periodic BC along the x-boundary
		f1 = phi[nn];						// get neighbor for phi - 1
		//........................................................................
		nn = n+1;							// neighbor index (get convention)
		if (!(i+1<Nx))	nn -= Nx;			// periodic BC along the x-boundary
		f2 = phi[nn];						// get neighbor for phi - 2
		//........................................................................
		nn = n-Nx;							// neighbor index (get convention)
		if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
		f3 = phi[nn];					// get neighbor for phi - 3
		//........................................................................
		nn = n+Nx;							// neighbor index (get convention)
		if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
		f4 = phi[nn];					// get neighbor for phi - 4
		//........................................................................
		nn = n-Nx*Ny;						// neighbor index (get convention)
		if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
		f5 = phi[nn];					// get neighbor for phi - 5
		//........................................................................
		nn = n+Nx*Ny;						// neighbor index (get convention)
		if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
		f6 = phi[nn];					// get neighbor for phi - 6
		//........................................................................
		nn = n-Nx-1;						// neighbor index (get convention)
		if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
		if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
		f7 = phi[nn];					// get neighbor for phi - 7
		//........................................................................
		nn = n+Nx+1;						// neighbor index (get convention)
		if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
		if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
		f8 = phi[nn];					// get neighbor for phi - 8
		//........................................................................
		nn = n+Nx-1;						// neighbor index (get convention)
		if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
		if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
		f9 = phi[nn];					// get neighbor for phi - 9
		//........................................................................
		nn = n-Nx+1;						// neighbor index (get convention)
		if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
		if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
		f10 = phi[nn];					// get neighbor for phi - 10
		//........................................................................
		nn = n-Nx*Ny-1;						// neighbor index (get convention)
		if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
		if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
		f11 = phi[nn];					// get neighbor for phi - 11
		//........................................................................
		nn = n+Nx*Ny+1;						// neighbor index (get convention)
		if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
		if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
		f12 = phi[nn];					// get neighbor for phi - 12
		//........................................................................
		nn = n+Nx*Ny-1;						// neighbor index (get convention)
		if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
		if (!(k+1<Nz))		nn -= Nx*Ny*Nz;	// Perioidic BC along the z-boundary
		f13 = phi[nn];					// get neighbor for phi - 13
		//........................................................................
		nn = n-Nx*Ny+1;						// neighbor index (get convention)
		if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
		if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
		f14 = phi[nn];					// get neighbor for phi - 14
		//........................................................................
		nn = n-Nx*Ny-Nx;					// neighbor index (get convention)
		if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
		if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
		f15 = phi[nn];					// get neighbor for phi - 15
		//........................................................................
		nn = n+Nx*Ny+Nx;					// neighbor index (get convention)
		if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
		if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
		f16 = phi[nn];					// get neighbor for phi - 16
		//........................................................................
		nn = n+Nx*Ny-Nx;					// neighbor index (get convention)
		if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
		if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
		f17 = phi[nn];					// get neighbor for phi - 17
		//........................................................................
		nn = n-Nx*Ny+Nx;					// neighbor index (get convention)
		if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
		if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
		f18 = phi[nn];					// get neighbor for phi - 18
		//............Compute the Color Gradient...................................
		nx = -(f1-f2+0.5*(f7-f8+f9-f10+f11-f12+f13-f14));
		ny = -(f3-f4+0.5*(f7-f8-f9+f10+f15-f16+f17-f18));
		nz = -(f5-f6+0.5*(f11-f12-f13+f14+f15-f16-f17+f18));
		//...........Normalize the Color Gradient.................................
		//	C = sqrt(nx*nx+ny*ny+nz*nz);
		//	nx = nx/C;
		//	ny = ny/C;
		//	nz = nz/C;
		//...Store the Color Gradient....................
		ColorGrad[n] = nx;
		ColorGrad[N+n] = ny;
		ColorGrad[2*N+n] = nz;
		//...............................................
	}
}
//*************************************************************************
extern "C" void ColorCollide( char *ID, double *disteven, double *distodd, double *ColorGrad,
								double *Velocity, int Nx, int Ny, int Nz, double rlx_setA, double rlx_setB,
								double alpha, double beta, double Fx, double Fy, double Fz, bool pBC)
{

	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	// additional variables needed for computations
	double rho,jx,jy,jz,C,nx,ny,nz;

	N = Nx*Ny*Nz;
	char id;

	for (n=0; n<N; n++){

		id = ID[n];

		if (id > 0){

			// Retrieve the color gradient
			nx = ColorGrad[n];
			ny = ColorGrad[N+n];
			nz = ColorGrad[2*N+n];
			//...........Normalize the Color Gradient.................................
			C = sqrt(nx*nx+ny*ny+nz*nz);
			if (C==0.0) C=1.0;
			nx = nx/C;
			ny = ny/C;
			nz = nz/C;
			//......No color gradient at z-boundary if pressure BC are set.............
			//	if (pBC && k==0) nx = ny = nz = 0.f;
			//	if (pBC && k==Nz-1) nx = ny = nz = 0.f;
			//........................................................................
			//					READ THE DISTRIBUTIONS
			//		(read from opposite array due to previous swap operation)
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
			//					PERFORM RELAXATION PROCESS
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
			//..........Toelke, Fruediger et. al. 2006...............
			if (C == 0.0)	nx = ny = nz = 1.0;
#ifdef STOKES
			m1 = m1 + rlx_setA*(- 11*rho -alpha*C - m1);
			m2 = m2 + rlx_setA*(3*rho - m2);
			m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
			m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
			m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
			m9 = m9 + rlx_setA*( 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
			m10 = m10 + rlx_setA*( - m10);
			m11 = m11 + rlx_setA*( 0.5*alpha*C*(ny*ny-nz*nz)- m11);
			m12 = m12 + rlx_setA*( - m12);
			m13 = m13 + rlx_setA*(  0.5*alpha*C*nx*ny - m13);
			m14 = m14 + rlx_setA*(  0.5*alpha*C*ny*nz - m14);
			m15 = m15 + rlx_setA*(  0.5*alpha*C*nx*nz - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);
#else
			m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho - 11*rho) -alpha*C - m1);
			m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho)- m2);
			m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
			m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
			m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
			m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
			m10 = m10 + rlx_setA*( - m10);
			m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
			m12 = m12 + rlx_setA*( - m12);
			m13 = m13 + rlx_setA*( (jx*jy/rho) + 0.5*alpha*C*nx*ny - m13);
			m14 = m14 + rlx_setA*( (jy*jz/rho) + 0.5*alpha*C*ny*nz - m14);
			m15 = m15 + rlx_setA*( (jx*jz/rho) + 0.5*alpha*C*nx*nz - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);
#endif
			//.................inverse transformation......................................................
			f0 = 0.05263157894736842*rho-0.012531328320802*m1+0.04761904761904762*m2;
			f1 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
					+0.1*(jx-m4)+0.0555555555555555555555555*(m9-m10);
			f2 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
					+0.1*(m4-jx)+0.0555555555555555555555555*(m9-m10);
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
					-0.0555555555555555555555555*m9-0.02777777777777778*m10+0.25*m14+0.125*(m17-m18);
			f16 =  0.05263157894736842*rho+0.003341687552213868*m1
					+0.003968253968253968*m2-0.1*(jy+jz)-0.025*(m6+m8)
					-0.0555555555555555555555555*m9-0.02777777777777778*m10+0.25*m14+0.125*(m18-m17);
			f17 = 0.05263157894736842*rho+0.003341687552213868*m1
					+0.003968253968253968*m2+0.1*(jy-jz)+0.025*(m6-m8)
					-0.0555555555555555555555555*m9-0.02777777777777778*m10-0.25*m14+0.125*(m17+m18);
			f18 = 0.05263157894736842*rho+0.003341687552213868*m1
					+0.003968253968253968*m2+0.1*(jz-jy)+0.025*(m8-m6)
					-0.0555555555555555555555555*m9-0.02777777777777778*m10-0.25*m14-0.125*(m17+m18);
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
			//*********** WRITE UPDATED VALUES TO MEMORY ******************
			// Write the updated distributions
			//....EVEN.....................................
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
			//....ODD......................................
			distodd[n] = f1;
			distodd[N+n] = f3;
			distodd[2*N+n] = f5;
			distodd[3*N+n] = f7;
			distodd[4*N+n] = f9;
			distodd[5*N+n] = f11;
			distodd[6*N+n] = f13;
			distodd[7*N+n] = f15;
			distodd[8*N+n] = f17;

			//...Store the Velocity..........................
			Velocity[n] = jx;
			Velocity[N+n] = jy;
			Velocity[2*N+n] = jz;
		/*	Velocity[3*n] = jx;
			Velocity[3*n+1] = jy;
			Velocity[3*n+2] = jz;
		*/	//...Store the Color Gradient....................
			//			ColorGrad[3*n] = nx*C;
			//			ColorGrad[3*n+1] = ny*C;
			//			ColorGrad[3*n+2] = nz*C;
			//...............................................
			//***************************************************************
		}	// check if n is in the solid
	} // loop over n
}

extern "C" void ColorCollideOpt( char *ID, double *disteven, double *distodd, double *phi, double *ColorGrad,
								double *Velocity, int Nx, int Ny, int Nz, double rlx_setA, double rlx_setB, 
								double alpha, double beta, double Fx, double Fy, double Fz)
{
		
	int i,j,k,n,nn,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	// additional variables needed for computations
	double rho,jx,jy,jz,C,nx,ny,nz;

	N = Nx*Ny*Nz;
	char id;

	for (n=0; n<N; n++){

		id = ID[n];

		if (id > 0){

			//.......Back out the 3-D indices for node n..............
			k = n/(Nx*Ny);
			j = (n-Nx*Ny*k)/Nx;
			i = n-Nx*Ny*k-Nx*j;
			//........................................................................
			//........Get 1-D index for this thread....................
			//		n = S*blockIdx.x*blockDim.x + s*blockDim.x + threadIdx.x;
			//........................................................................
			//					COMPUTE THE COLOR GRADIENT
			//........................................................................
			//.................Read Phase Indicator Values............................
			//........................................................................
			nn = n-1;							// neighbor index (get convention)
			if (i-1<0)		nn += Nx;			// periodic BC along the x-boundary
			f1 = phi[nn];						// get neighbor for phi - 1
			//........................................................................
			nn = n+1;							// neighbor index (get convention)
			if (!(i+1<Nx))	nn -= Nx;			// periodic BC along the x-boundary
			f2 = phi[nn];						// get neighbor for phi - 2
			//........................................................................
			nn = n-Nx;							// neighbor index (get convention)
			if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
			f3 = phi[nn];					// get neighbor for phi - 3
			//........................................................................
			nn = n+Nx;							// neighbor index (get convention)
			if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
			f4 = phi[nn];					// get neighbor for phi - 4
			//........................................................................
			nn = n-Nx*Ny;						// neighbor index (get convention)
			if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
			f5 = phi[nn];					// get neighbor for phi - 5
			//........................................................................
			nn = n+Nx*Ny;						// neighbor index (get convention)
			if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
			f6 = phi[nn];					// get neighbor for phi - 6
			//........................................................................
			nn = n-Nx-1;						// neighbor index (get convention)
			if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
			if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
			f7 = phi[nn];					// get neighbor for phi - 7
			//........................................................................
			nn = n+Nx+1;						// neighbor index (get convention)
			if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
			if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
			f8 = phi[nn];					// get neighbor for phi - 8
			//........................................................................
			nn = n+Nx-1;						// neighbor index (get convention)
			if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
			if (!(j+1<Ny))		nn -= Nx*Ny;	// Perioidic BC along the y-boundary
			f9 = phi[nn];					// get neighbor for phi - 9
			//........................................................................
			nn = n-Nx+1;						// neighbor index (get convention)
			if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
			if (j-1<0)			nn += Nx*Ny;	// Perioidic BC along the y-boundary
			f10 = phi[nn];					// get neighbor for phi - 10
			//........................................................................
			nn = n-Nx*Ny-1;						// neighbor index (get convention)
			if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
			if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
			f11 = phi[nn];					// get neighbor for phi - 11
			//........................................................................
			nn = n+Nx*Ny+1;						// neighbor index (get convention)
			if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
			if (!(k+1<Nz))		nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
			f12 = phi[nn];					// get neighbor for phi - 12
			//........................................................................
			nn = n+Nx*Ny-1;						// neighbor index (get convention)
			if (i-1<0)			nn += Nx;		// periodic BC along the x-boundary
			if (!(k+1<Nz))		nn -= Nx*Ny*Nz;	// Perioidic BC along the z-boundary
			f13 = phi[nn];					// get neighbor for phi - 13
			//........................................................................
			nn = n-Nx*Ny+1;						// neighbor index (get convention)
			if (!(i+1<Nx))		nn -= Nx;		// periodic BC along the x-boundary
			if (k-1<0)			nn += Nx*Ny*Nz;	// Perioidic BC along the z-boundary
			f14 = phi[nn];					// get neighbor for phi - 14
			//........................................................................
			nn = n-Nx*Ny-Nx;					// neighbor index (get convention)
			if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
			if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
			f15 = phi[nn];					// get neighbor for phi - 15
			//........................................................................
			nn = n+Nx*Ny+Nx;					// neighbor index (get convention)
			if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
			if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
			f16 = phi[nn];					// get neighbor for phi - 16
			//........................................................................
			nn = n+Nx*Ny-Nx;					// neighbor index (get convention)
			if (j-1<0)		nn += Nx*Ny;		// Perioidic BC along the y-boundary
			if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
			f17 = phi[nn];					// get neighbor for phi - 17
			//........................................................................
			nn = n-Nx*Ny+Nx;					// neighbor index (get convention)
			if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
			if (k-1<0)		nn += Nx*Ny*Nz;		// Perioidic BC along the z-boundary
			f18 = phi[nn];					// get neighbor for phi - 18
			//............Compute the Color Gradient...................................
			nx = -(f1-f2+0.5*(f7-f8+f9-f10+f11-f12+f13-f14));
			ny = -(f3-f4+0.5*(f7-f8-f9+f10+f15-f16+f17-f18));
			nz = -(f5-f6+0.5*(f11-f12-f13+f14+f15-f16-f17+f18));
			//...Store the Color Gradient....................
			ColorGrad[n] = nx;
			ColorGrad[N+n] = ny;
			ColorGrad[2*N+n] = nz;
			//...............................................
			//...........Normalize the Color Gradient.................................
			C = sqrt(nx*nx+ny*ny+nz*nz);
			if (C==0.0) C=1.0;
			nx = nx/C;
			ny = ny/C;
			nz = nz/C;
			//......No color gradient at z-boundary if pressure BC are set.............
			//	if (pBC && k==0) nx = ny = nz = 0.f;
			//	if (pBC && k==Nz-1) nx = ny = nz = 0.f;
			//........................................................................
			//					READ THE DISTRIBUTIONS
			//		(read from opposite array due to previous swap operation)
			//........................................................................
			f2 = distodd[n];
			f4 = distodd[N+n];
			f6 = distodd[2*N+n];
			f0 = disteven[n];
			f1 = disteven[N+n];
			f3 = disteven[2*N+n];
			f5 = disteven[3*N+n];
			//........................................................................
			//....................compute the moments...............................................
			rho = f0+f2+f1+f4+f3+f6+f5;
			m1 = -30*f0-11*(f2+f1+f4+f3+f6+f5);
			m2 = 12*f0-4*(f2+f1 +f4+f3+f6 +f5);
			jx = f1-f2;
			m4 = 4*(-f1+f2);
			jy = f3-f4;
			m6 = -4*(f3-f4);
			jz = f5-f6;
			m8 = -4*(f5-f6);
			m9 = 2*(f1+f2)-f3-f4-f5-f6;
			m10 = -4*(f1+f2)+2*(f4+f3+f6+f5);
			m11 = f4+f3-f6-f5;
			m12 = -2*(f4+f3-f6-f5);
			//........................................................................
			f8 = distodd[3*N+n];
			f10 = distodd[4*N+n];
			f7 = disteven[4*N+n];
			f9 = disteven[5*N+n];
			//........................................................................
			rho += f8+f7+f10+f9;
			m1 += 8*(f8+f7+f10+f9);
			m2 += f8+f7+f10+f9;
			jx += f7-f8+f9-f10;
			m4 += f7-f8+f9-f10;
			jy += f7-f8-f9+f10;
			m6 += f7-f8-f9+f10;
			m9 += f7+f8+f9+f10;
			m10 += f8+f7+f10+f9;
			m11 += f8+f7+f10+f9;
			m12 += f8+f7+f10+f9;
			m13 = f8+f7-f10-f9;
			m16 = f7-f8+f9-f10;
			m17 = -f7+f8+f9-f10;
			//........................................................................
			f11 = disteven[6*N+n];
			f13 = disteven[7*N+n];
			f12 = distodd[5*N+n];
			f14 = distodd[6*N+n];
			//........................................................................
			//........................................................................
			f15 = disteven[8*N+n];
			f17 = disteven[9*N+n];
			f16 = distodd[7*N+n];
			f18 = distodd[8*N+n];
			//........................................................................
			//....................compute the moments...............................................
			rho += f12+f11+f14+f13+f16+f15+f18+f17;
			m1 += 8*(f12+f11+f14+f13+f16+f15+f18+f17);
			m2 += f12+f11+f14+f13+f16+f15+f18+f17;
			jx += f11-f12+f13-f14;
			m4 += f11-f12+f13-f14;
			jy += f15-f16+f17-f18;
			m6 += f15-f16+f17-f18;
			jz += f11-f12-f13+f14+f15-f16-f17+f18;
			m8 += f11-f12-f13+f14+f15-f16-f17+f18;
			m9 += f11+f12+f13+f14-2*(f15+f16+f17+f18);
			m10 += f12+f11+f14+f13-2*(f16+f15+f18+f17);
			m11 += -f12-f11-f14-f13;
			m12 += -f12-f11-f14-f13;
			m14 = f16+f15-f18-f17;
			m15 = f12+f11-f14-f13;
			m16 += -f11+f12-f13+f14;
			m17 += f15-f16+f17-f18;
			m18 = f11-f12-f13+f14-f15+f16+f17-f18;
			//........................................................................

			/*				f2 = distodd[n];
				f4 = distodd[N+n];
				f6 = distodd[2*N+n];
				f8 = distodd[3*N+n];
				//........................................................................
				f0 = disteven[n];
				f1 = disteven[N+n];
				f3 = disteven[2*N+n];
				f5 = disteven[3*N+n];
				f7 = disteven[4*N+n];
				//........................................................................
				//........................................................................
				//....................compute the moments...............................................
				rho = f0+f2+f1+f4+f3+f6+f5+f8+f7;
				m1 = -30*f0-11*(f2+f1+f4+f3+f6+f5)+8*(f8+f7);
				m2 = 12*f0-4*(f2+f1 +f4+f3+f6 +f5)+f8+f7;
				jx = f1-f2+f7-f8;
				m4 = 4*(-f1+f2)+f7-f8;
				jy = f3-f4+f7-f8;
				m6 = -4*(f3-f4)+f7-f8;
				jz = f5-f6;
				m8 = -4*(f5-f6);
				m9 = 2*(f1+f2)-f3-f4-f5-f6+f7+f8;
				m10 = -4*(f1+f2)+2*(f4+f3+f6+f5)+f8+f7;
				m11 = f4+f3-f6-f5+f8+f7;
				m12 = -2*(f4+f3-f6-f5)+f8+f7;
				m13 = f8+f7;
				m16 = f7-f8;
				m17 = -f7+f8;
				//........................................................................
				f9 = disteven[5*N+n];
				f11 = disteven[6*N+n];
				f13 = disteven[7*N+n];
				f15 = disteven[8*N+n];
				f17 = disteven[9*N+n];
				f10 = distodd[4*N+n];
				f12 = distodd[5*N+n];
				f14 = distodd[6*N+n];
				f16 = distodd[7*N+n];
				f18 = distodd[8*N+n];
				//........................................................................
				rho += f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
				m1 += 8*(f10+f9+f12+f11+f14+f13+f16+f15+f18 +f17);
				m2 += f10+f9+f12+f11+f14+f13+f16+f15+f18+f17;
				jx += f9-f10+f11-f12+f13-f14;
				m4 += f9-f10+f11-f12+f13-f14;
				jy += -f9+f10+f15-f16+f17-f18;
				m6 += -f9+f10+f15-f16+f17-f18;
				jz += f11-f12-f13+f14+f15-f16-f17+f18;
				m8 += f11-f12-f13+f14+f15-f16-f17+f18;
				m9 += f9+f10+f11+f12+f13+f14-2*(f15+f16+f17+f18);
				m10 += f10+f9+f12+f11+f14+f13-2*(f16+f15+f18+f17);
				m11 += f10+f9-f12-f11-f14-f13;
				m12 += f10+f9-f12-f11-f14-f13;
				m13 += -f10-f9;
				m14 = f16+f15-f18-f17;
				m15 = f12+f11-f14-f13;
				m16 += f9-f10-f11+f12-f13+f14;
				m17 += f9-f10+f15-f16+f17-f18;
				m18 = f11-f12-f13+f14-f15+f16+f17-f18;
			 */			//........................................................................
			//					PERFORM RELAXATION PROCESS
			//........................................................................
			//..........Toelke, Fruediger et. al. 2006...............
			if (C == 0.0)	nx = ny = nz = 1.0;
			m1 = m1 + rlx_setA*((19*(jx*jx+jy*jy+jz*jz)/rho - 11*rho) -alpha*C - m1);
			m2 = m2 + rlx_setA*((3*rho - 5.5*(jx*jx+jy*jy+jz*jz)/rho)- m2);
			m4 = m4 + rlx_setB*((-0.6666666666666666*jx)- m4);
			m6 = m6 + rlx_setB*((-0.6666666666666666*jy)- m6);
			m8 = m8 + rlx_setB*((-0.6666666666666666*jz)- m8);
			m9 = m9 + rlx_setA*(((2*jx*jx-jy*jy-jz*jz)/rho) + 0.5*alpha*C*(2*nx*nx-ny*ny-nz*nz) - m9);
			m10 = m10 + rlx_setA*( - m10);
			m11 = m11 + rlx_setA*(((jy*jy-jz*jz)/rho) + 0.5*alpha*C*(ny*ny-nz*nz)- m11);
			m12 = m12 + rlx_setA*( - m12);
			m13 = m13 + rlx_setA*( (jx*jy/rho) + 0.5*alpha*C*nx*ny - m13);
			m14 = m14 + rlx_setA*( (jy*jz/rho) + 0.5*alpha*C*ny*nz - m14);
			m15 = m15 + rlx_setA*( (jx*jz/rho) + 0.5*alpha*C*nx*nz - m15);
			m16 = m16 + rlx_setB*( - m16);
			m17 = m17 + rlx_setB*( - m17);
			m18 = m18 + rlx_setB*( - m18);
			//.................inverse transformation......................................................
			f0 = 0.05263157894736842*rho-0.012531328320802*m1+0.04761904761904762*m2;
			f1 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
					+0.1*(jx-m4)+0.0555555555555555555555555*(m9-m10);
			f2 = 0.05263157894736842*rho-0.004594820384294068*m1-0.01587301587301587*m2
					+0.1*(m4-jx)+0.0555555555555555555555555*(m9-m10);
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
					-0.0555555555555555555555555*m9-0.02777777777777778*m10+0.25*m14+0.125*(m17-m18);
			f16 =  0.05263157894736842*rho+0.003341687552213868*m1
					+0.003968253968253968*m2-0.1*(jy+jz)-0.025*(m6+m8)
					-0.0555555555555555555555555*m9-0.02777777777777778*m10+0.25*m14+0.125*(m18-m17);
			f17 = 0.05263157894736842*rho+0.003341687552213868*m1
					+0.003968253968253968*m2+0.1*(jy-jz)+0.025*(m6-m8)
					-0.0555555555555555555555555*m9-0.02777777777777778*m10-0.25*m14+0.125*(m17+m18);
			f18 = 0.05263157894736842*rho+0.003341687552213868*m1
					+0.003968253968253968*m2+0.1*(jz-jy)+0.025*(m8-m6)
					-0.0555555555555555555555555*m9-0.02777777777777778*m10-0.25*m14-0.125*(m17+m18);
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
			//*********** WRITE UPDATED VALUES TO MEMORY ******************
			// Write the updated distributions
			//....EVEN.....................................
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
			//....ODD......................................
			distodd[n] = f1;
			distodd[N+n] = f3;
			distodd[2*N+n] = f5;
			distodd[3*N+n] = f7;
			distodd[4*N+n] = f9;
			distodd[5*N+n] = f11;
			distodd[6*N+n] = f13;
			distodd[7*N+n] = f15;
			distodd[8*N+n] = f17;
			//...Store the Velocity..........................
			Velocity[n] = jx;
			Velocity[N+n] = jy;
			Velocity[2*N+n] = jz;
			//***************************************************************
		}	// check if n is in the solid
	} // loop over n
}

extern "C" void MassColorCollideD3Q7(char *ID, double *A_even, double *A_odd, double *B_even, double *B_odd, 
		double *Den, double *Phi, double *ColorGrad, double *Velocity, double beta, int N, bool pBC)
{
	char id;

	int idx,n,q,Cqx,Cqy,Cqz;
	//	int sendLoc;

	double f0,f1,f2,f3,f4,f5,f6;
	double na,nb;		// density values
	double ux,uy,uz;	// flow velocity
	double nx,ny,nz,C;	// color gradient components
	double a1,a2,b1,b2;
	double sp,delta;
	//double feq[6];		// equilibrium distributions
	// Set of Discrete velocities for the D3Q19 Model
	//int D3Q7[3][3]={{1,0,0},{0,1,0},{0,0,1}};

	for (n=0; n<N; n++){
		id = ID[n];
		if (id != 0 ){

			//.....Load the Color gradient.........
			nx = ColorGrad[n];
			ny = ColorGrad[N+n];
			nz = ColorGrad[2*N+n];
			C = sqrt(nx*nx+ny*ny+nz*nz);
			if (C==0.0) C=1.0;
			nx = nx/C;
			ny = ny/C;
			nz = nz/C;
			//....Load the flow velocity...........
			ux = Velocity[n];
			uy = Velocity[N+n];
			uz = Velocity[2*N+n];
			//........................................................................
			//					READ THE DISTRIBUTIONS
			//		(read from opposite array due to previous swap operation)
			//........................................................................
			f2 = A_odd[n];
			f4 = A_odd[N+n];
			f6 = A_odd[2*N+n];
			f0 = A_even[n];
			f1 = A_even[N+n];
			f3 = A_even[2*N+n];
			f5 = A_even[3*N+n];
			na = f0+f1+f2+f3+f4+f5+f6;
			//........................................................................
			f2 = B_odd[n];
			f4 = B_odd[N+n];
			f6 = B_odd[2*N+n];
			f0 = B_even[n];
			f1 = B_even[N+n];
			f3 = B_even[2*N+n];
			f5 = B_even[3*N+n];
			nb = f0+f1+f2+f3+f4+f5+f6;
			//........................................................................
			//....Instantiate the density distributions
			// Generate Equilibrium Distributions and stream
			// Stationary value - distribution 0
			A_even[n] = 0.3333333333333333*na;
			B_even[n] = 0.3333333333333333*nb;
			// Non-Stationary equilibrium distributions
			//feq[0] = 0.1111111111111111*(1+4.5*ux);
			//feq[1] = 0.1111111111111111*(1-4.5*ux);
			//feq[2] = 0.1111111111111111*(1+4.5*uy);
			//feq[3] = 0.1111111111111111*(1-4.5*uy);
			//feq[4] = 0.1111111111111111*(1+4.5*uz);
			//feq[5] = 0.1111111111111111*(1-4.5*uz);
			
			//...............................................
			// q = 0,2,4
			// Cq = {1,0,0}, {0,1,0}, {0,0,1}
			a1 = na*(0.1111111111111111*(1+4.5*ux));
			b1 = nb*(0.1111111111111111*(1+4.5*ux));
			a2 = na*(0.1111111111111111*(1-4.5*ux));
			b2 = nb*(0.1111111111111111*(1-4.5*ux));
			if (C > 0.0){
				delta = na*nb/(na+nb)*0.1111111111111111*nx;
				a1 += beta*delta;
				a2 -= beta*delta;
				b1 -= beta*delta;
				b2 += beta*delta;
			}
			A_odd[n] 	= a1;
			A_even[N+n] = a2;
			B_odd[n] 	= b1;
			B_even[N+n] = b2;
			//...............................................
			// q = 2
			// Cq = {0,1,0}
			a1 = na*(0.1111111111111111*(1+4.5*uy));
			b1 = nb*(0.1111111111111111*(1+4.5*uy));
			a2 = na*(0.1111111111111111*(1-4.5*uy));
			b2 = nb*(0.1111111111111111*(1-4.5*uy));
			if (C > 0.0){
				delta = na*nb/(na+nb)*0.1111111111111111*ny;
				a1 += beta*delta;
				a2 -= beta*delta;
				b1 -= beta*delta;
				b2 += beta*delta;
			}
			A_odd[N+n] 	= a1;
			A_even[2*N+n] = a2;
			B_odd[N+n] 	= b1;
			B_even[2*N+n] = b2;
			//...............................................
			// q = 4
			// Cq = {0,0,1}
			a1 = na*(0.1111111111111111*(1+4.5*uz));
			b1 = nb*(0.1111111111111111*(1+4.5*uz));
			a2 = na*(0.1111111111111111*(1-4.5*uz));
			b2 = nb*(0.1111111111111111*(1-4.5*uz));
			if (C > 0.0){
				delta = na*nb/(na+nb)*0.1111111111111111*nz;
				a1 += beta*delta;
				a2 -= beta*delta;
				b1 -= beta*delta;
				b2 += beta*delta;
			}
			A_odd[2*N+n] = a1;
			A_even[3*N+n] = a2;
			B_odd[2*N+n] = b1;
			B_even[3*N+n] = b2;
			//...............................................

	/*		// Construction and streaming for the components
			for (idx=0; idx<3; idx++){
				//...............................................
				// Distribution index
				q = 2*idx;
				// Associated discrete velocity
				Cqx = D3Q7[idx][0];
				Cqy = D3Q7[idx][1];
				Cqz = D3Q7[idx][2];
				// Generate the Equilibrium Distribution
				a1 = na*feq[q];
				b1 = nb*feq[q];
				a2 = na*feq[q+1];
				b2 = nb*feq[q+1];
				// Recolor the distributions
				if (C > 0.0){
					sp = nx*double(Cqx)+ny*double(Cqy)+nz*double(Cqz);
					//if (idx > 2)	sp = 0.7071067811865475*sp;
					//delta = sp*min( min(a1,a2), min(b1,b2) );
					delta = na*nb/(na+nb)*0.1111111111111111*sp;
					//if (a1>0 && b1>0){
					a1 += beta*delta;
					a2 -= beta*delta;
					b1 -= beta*delta;
					b2 += beta*delta;
				}
				// Save the re-colored distributions
				A_odd[N*idx+n] 		= a1;
				A_even[N*(idx+1)+n] = a2;
				B_odd[N*idx+n] 		= b1;
				B_even[N*(idx+1)+n] = b2;
				//...............................................
			}
	*/
		}
	}
}

//*************************************************************************
extern "C" void DensityStreamD3Q7(char *ID, double *Den, double *Copy, double *Phi, double *ColorGrad, double *Velocity,
		double beta, int Nx, int Ny, int Nz, bool pBC, int S)
{
	char id;

	int idx;
	int in,jn,kn,n,nn,N;
	int q,Cqx,Cqy,Cqz;
	//	int sendLoc;

	double na,nb;		// density values
	double ux,uy,uz;	// flow velocity
	double nx,ny,nz,C;	// color gradient components
	double a1,a2,b1,b2;
	double sp,delta;
	double feq[6];		// equilibrium distributions
	// Set of Discrete velocities for the D3Q19 Model
	int D3Q7[3][3]={{1,0,0},{0,1,0},{0,0,1}};
	N = Nx*Ny*Nz;

	for (n=0; n<N; n++){
		id = ID[n];
		// Local Density Values
		na = Copy[2*n];
		nb = Copy[2*n+1];
		if (id > 0 && na+nb > 0.0){
			//.......Back out the 3-D indices for node n..............
			int	k = n/(Nx*Ny);
			int j = (n-Nx*Ny*k)/Nx;
			int i = n-Nx*Ny*k-Nx*j;
			//.....Load the Color gradient.........
			nx = ColorGrad[n];
			ny = ColorGrad[N+n];
			nz = ColorGrad[2*N+n];
			C = sqrt(nx*nx+ny*ny+nz*nz);
			nx = nx/C;
			ny = ny/C;
			nz = nz/C;
			//....Load the flow velocity...........
			ux = Velocity[n];
			uy = Velocity[N+n];
			uz = Velocity[2*N+n];
			//....Instantiate the density distributions
			// Generate Equilibrium Distributions and stream
			// Stationary value - distribution 0
			//			Den[2*n] += 0.3333333333333333*na;
			//			Den[2*n+1] += 0.3333333333333333*nb;
			Den[2*n] += 0.3333333333333333*na;
			Den[2*n+1] += 0.3333333333333333*nb;
			// Non-Stationary equilibrium distributions
			feq[0] = 0.1111111111111111*(1+3*ux);
			feq[1] = 0.1111111111111111*(1-3*ux);
			feq[2] = 0.1111111111111111*(1+3*uy);
			feq[3] = 0.1111111111111111*(1-3*uy);
			feq[4] = 0.1111111111111111*(1+3*uz);
			feq[5] = 0.1111111111111111*(1-3*uz);
			// Construction and streaming for the components
			for (idx=0; idx<3; idx++){
				// Distribution index
				q = 2*idx;
				// Associated discrete velocity
				Cqx = D3Q7[idx][0];
				Cqy = D3Q7[idx][1];
				Cqz = D3Q7[idx][2];
				// Generate the Equilibrium Distribution
				a1 = na*feq[q];
				b1 = nb*feq[q];
				a2 = na*feq[q+1];
				b2 = nb*feq[q+1];
				// Recolor the distributions
				if (C > 0.0){
					sp = nx*double(Cqx)+ny*double(Cqy)+nz*double(Cqz);
					//if (idx > 2)	sp = 0.7071067811865475*sp;
					//delta = sp*min( min(a1,a2), min(b1,b2) );
					delta = na*nb/(na+nb)*0.1111111111111111*sp;
					//if (a1>0 && b1>0){
					a1 += beta*delta;
					a2 -= beta*delta;
					b1 -= beta*delta;
					b2 += beta*delta;
				}

				// .......Get the neighbor node..............
				//nn = n + Stride[idx];
				in = i+Cqx;
				jn = j+Cqy;
				kn = k+Cqz;

				// Adjust for periodic BC, if necessary
				//				if (in<0) in+= Nx;
				//				if (jn<0) jn+= Ny;
				//				if (kn<0) kn+= Nz;
				//				if (!(in<Nx)) in-= Nx;
				//				if (!(jn<Ny)) jn-= Ny;
				//				if (!(kn<Nz)) kn-= Nz;
				// Perform streaming or bounce-back as needed
				id = ID[kn*Nx*Ny+jn*Nx+in];
				if (id == 0){							//.....Bounce-back Rule...........
					//						Den[2*n] += a1;
					//						Den[2*n+1] += b1;
					Den[2*n] += a1;
					Den[2*n+1] += b1;
				}
				else{
					//......Push the "distribution" to neighboring node...........
					// Index of the neighbor in the local process
					//nn = (kn-zmin[rank]+1)*Nxp*Nyp + (jn-ymin[rank]+1)*Nxp + (in-xmin[rank]+1);
					nn = kn*Nx*Ny+jn*Nx+in;
					// Push to neighboring node
					//						Den[2*nn] += a1;
					//						Den[2*nn+1] += b1;
					Den[2*nn] += a1;
					Den[2*nn+1] += b1;
				}

				// .......Get the neighbor node..............
				q = 2*idx+1;
				in = i-Cqx;
				jn = j-Cqy;
				kn = k-Cqz;
				// Adjust for periodic BC, if necessary
				//				if (in<0) in+= Nx;
				//				if (jn<0) jn+= Ny;
				//				if (kn<0) kn+= Nz;
				//				if (!(in<Nx)) in-= Nx;
				//				if (!(jn<Ny)) jn-= Ny;
				//				if (!(kn<Nz)) kn-= Nz;
				// Perform streaming or bounce-back as needed
				id = ID[kn*Nx*Ny+jn*Nx+in];
				if (id == 0){
					//.....Bounce-back Rule...........
					//						Den[2*n] += a2;
					//					Den[2*n+1] += b2;
					Den[2*n] += a2;
					Den[2*n+1] += b2;
				}
				else{
					//......Push the "distribution" to neighboring node...........
					// Index of the neighbor in the local process
					//nn = (kn-zmin[rank]+1)*Nxp*Nyp + (jn-ymin[rank]+1)*Nxp + (in-xmin[rank]+1);
					nn = kn*Nx*Ny+jn*Nx+in;
					// Push to neighboring node
					//					Den[2*nn] += a2;
					//					Den[2*nn+1] += b2;
					Den[2*nn] += a2;
					Den[2*nn+1] += b2;
				}
			}
		}
	}
}
/*
extern "C" void ComputePhi(char *ID, double *Phi, double *Copy, double *Den, int N, int S)
{
	int n;
	double Na,Nb;
	//...................................................................
	// Update Phi
	for (n=0; n<N; n++){

		if (ID[n] > 0 && n<N){
			// Get the density value (Streaming already performed)
			Na = Den[2*n];
			Nb = Den[2*n+1];
			Phi[n] = (Na-Nb)/(Na+Nb);
			// Store the copy of the current density
			Copy[2*n] = Na;
			Copy[2*n+1] = Nb;
			// Zero the Density value to get ready for the next streaming
			Den[2*n] = 0.0;
			Den[2*n+1] = 0.0;
		}
	}
	//...................................................................
}
*/
extern "C" void ComputePhi(char *ID, double *Phi, double *Den, int N)
{
	int n;
	double Na,Nb;
	//...................................................................
	// Update Phi
	for (n=0; n<N; n++){

		if (ID[n] > 0 && n<N){
			// Get the density value (Streaming already performed)
			Na = Den[n];
			Nb = Den[N+n];
			Phi[n] = (Na-Nb)/(Na+Nb);
		}
	}
	//...................................................................
}

/*
//*************************************************************************
extern "C" void InitDenColor( int nblocks, int nthreads, int S,
		char *ID, double *Den, double *Phi, double das, double dbs, int Nx, int Ny, int Nz)
{
	InitDenColor <<<nblocks, nthreads>>>  (ID, Den, Phi, das, dbs, Nx, Ny, Nz, S);
}
//*************************************************************************
extern "C" void ComputeColorGradient(int nBlocks, int nthreads, int S,
		char *ID, double *Phi, double *ColorGrad, int Nx, int Ny, int Nz)
{
	ComputeColorGradient<<<nBlocks,nthreads>>>(ID, Phi, ColorGrad, Nx, Ny, Nz, S);
}
//*************************************************************************
extern "C" void ColorCollide(int nBlocks, int nthreads, int S,
		char *ID, double *f_even, double *f_odd, double *ColorGrad, double *Velocity,
		double rlxA, double rlxB,double alpha, double beta, double Fx, double Fy, double Fz,
		int Nx, int Ny, int Nz, bool pBC)
{
	ColorCollide<<<nBlocks, nthreads>>>(ID, f_even, f_odd, ColorGrad, Velocity, Nx, Ny, Nz, S,
							 rlxA, rlxB, alpha, beta, Fx, Fy, Fz, pBC);
}
//*************************************************************************
extern "C" void ColorCollideOpt(int nBlocks, int nthreads, int S,
								char *ID, double *f_even, double *f_odd, double *Phi, double *ColorGrad,
								double *Velocity, int Nx, int Ny, int Nz,double rlxA, double rlxB, 
								double alpha, double beta, double Fx, double Fy, double Fz)
{
	ColorCollideOpt<<<nBlocks, nthreads>>>(ID, f_even, f_odd, Phi, ColorGrad, Velocity, Nx, Ny, Nz, S,
							 rlxA, rlxB, alpha, beta, Fx, Fy, Fz);
//	bool pBC = false;
//	ColorCollide<<<nBlocks, nthreads>>>(ID, f_even, f_odd, ColorGrad, Velocity, Nx, Ny, Nz, S,
//							 rlxA, rlxB, alpha, beta, Fx, Fy, Fz, pBC);
}

//*************************************************************************
extern "C" void DensityStreamD3Q7(int nBlocks, int nthreads, int S,
		char *ID, double *Den, double *Copy, double *Phi, double *ColorGrad, double *Velocity,
		double beta, int Nx, int Ny, int Nz, bool pBC)
{
	DensityStreamD3Q7<<<nBlocks, nthreads>>>(ID,Den,Copy,Phi,ColorGrad,Velocity,beta,Nx,Ny,Nz,pBC,S);
}
//*************************************************************************
extern "C" void ComputePhi(int nBlocks, int nthreads, int S,
		char *ID, double *Phi, double *Copy, double *Den, int N)
{
	ComputePhi<<<nBlocks, nthreads>>>(ID,Phi,Copy,Den,N,S);
}
//*************************************************************************
extern "C" void ComputePressure(int nBlocks, int nthreads, int S,
									char *ID, double *disteven, double *distodd, 
									double *Pressure, int Nx, int Ny, int Nz)
{

	ComputePressureD3Q19<<<nBlocks, nthreads>>>(ID,disteven,distodd,Pressure,Nx,Ny,Nz,S);

}
*/
