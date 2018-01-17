#include <stdio.h>

extern "C" void ScaLBL_D3Q19_Pack(int q, int *list, int start, int count, double *sendbuf, double *dist, int N){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[start+idx] = dist[q*N+n];
	}
}

/*extern "C" void ScaLBL_D3Q19_Unpack(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
					   double *recvbuf, double *dist, int Nx, int Ny, int Nz){
	//....................................................................................
	// Unack distribution from the recv buffer
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int i,j,k,n,nn,idx;
	int N = Nx*Ny*Nz;
	for (idx=0; idx<count; idx++){
		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = list[idx];
		// Get the 3-D indices
		k = n/(Nx*Ny);
		j = (n-Nx*Ny*k)/Nx;
		i = n-Nx*Ny*k-Nx*j;
		// Streaming for the non-local distribution
		i += Cqx;
		j += Cqy;
		k += Cqz;

		nn = k*Nx*Ny+j*Nx+i;

		// unpack the distribution to the proper location
		dist[q*N+nn] = recvbuf[start+idx];
	}
}*/


extern "C" void ScaLBL_D3Q19_Unpack(int q, int *list,  int start, int count,
					   double *recvbuf, double *dist, int N){
	//....................................................................................
	// Unack distribution from the recv buffer
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int n,idx;
	for (idx=0; idx<count; idx++){
		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = list[start+idx];
		// unpack the distribution to the proper location
		if (!(n<0)) dist[q*N+n] = recvbuf[start+idx];
		//dist[q*N+n] = recvbuf[start+idx];

	}
}

/*
extern "C" void ScaLBL_D3Q19_MapRecv(int q, int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
					   int *d3q19_recvlist, int Nx, int Ny, int Nz){
	//....................................................................................
	// Map the recieve distributions to
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................

	int i,j,k,n,nn,idx;
	int N = Nx*Ny*Nz;
	for (idx=0; idx<count; idx++){
		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = list[idx];
		// Get the 3-D indices
		k = n/(Nx*Ny);
		j = (n-Nx*Ny*k)/Nx;
		i = n-Nx*Ny*k-Nx*j;
		// Streaming for the non-local distribution
		i += Cqx;
		j += Cqy;
		k += Cqz;
		// compute 1D index for the neighbor and save
		nn = k*Nx*Ny+j*Nx+i;
		d3q19_recvlist[start+idx] = nn;
	}
}

*/
extern "C" void ScaLBL_D3Q19_Init(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz)
{
	int n,N;
	N = Nx*Ny*Nz;

	for (n=0; n<N; n++){

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

//*************************************************************************
extern "C" void ScaLBL_D3Q19_Swap(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz)
{
	int i,j,k,n,nn,N;
	// distributions
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	
	N = Nx*Ny*Nz;
	
	for (n=0; n<N; n++){
		//.......Back out the 3-D indices for node n..............
		k = n/(Nx*Ny);
		j = (n-Nx*Ny*k)/Nx;
		i = n-Nx*Ny*k-Nx*j;
		
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

extern "C" void ScaLBL_D3Q19_Swap_Compact(int *neighborList, double *disteven, double *distodd, int Np)
{
	int q,n,nn;
	double f1,f2;
	for (q=0; q<9; q++){
		for (n=0; n<Np; n++){
			nn = neighborList[q*Np+n];
			if (!(nn<0)){
				f1 = distodd[q*Np+n];
				f2 = disteven[(q+1)*Np+nn];
				disteven[(q+1)*Np+nn] = f1;
				distodd[q*Np+n] = f2;
			}
		}
	}
}


extern "C" double ScaLBL_D3Q19_Flux_BC_z(char *ID,  double *disteven, double *distodd, double Q, double area,
								  int Nx, int Ny, int Nz){
	// Note that this routine assumes the distributions are stored "opposite"
	// odd distributions in disteven and even distributions in distodd.
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double din = 0.f;
	N = Nx*Ny*Nz;

	double sum = 0.f;
	char id;
	for (n=Nx*Ny; n<2*Nx*Ny; n++){
		id = ID[n];
		if (id > 0){
			//........................................................................
			// Read distributions from "opposite" memory convention
			//........................................................................
			//........................................................................
			f2 = distodd[n];
			f4 = distodd[N+n];
			f6 = distodd[2*N+n];
			f8 = distodd[3*N+n];
			f10 = distodd[4*N+n];
			f12 = distodd[5*N+n];
			//f14 = distodd[6*N+n];
			f16 = distodd[7*N+n];
			//f18 = distodd[8*N+n];
			//........................................................................
			f0 = disteven[n];
			f1 = disteven[N+n];
			f3 = disteven[2*N+n];
			//f5 = disteven[3*N+n];
			f7 = disteven[4*N+n];
			f9 = disteven[5*N+n];
			//f11 = disteven[6*N+n];
			f13 = disteven[7*N+n];
			//f15 = disteven[8*N+n];
			f17 = disteven[9*N+n];
			//...................................................

			// Determine the outlet flow velocity
			//sum += 1.0 - (f0+f4+f3+f2+f1+f8+f7+f9+ f10 + 2*(f5+ f15+f18+f11+f14))/din;
			//sum += (f0+f4+f3+f2+f1+f8+f7+f9+ f10 + 2*(f5+f15+f18+f11+f14));
			sum += (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));
		}
	}
	din = sum/area;

	return din;
}

extern "C" double ScaLBL_D3Q19_Flux_BC_Z(char *ID, double *disteven, double *distodd, double Q, double area,
		int Nx, int Ny, int Nz, int outlet){
	// Note that this routine assumes the distributions are stored "opposite"
	// odd distributions in disteven and even distributions in distodd.
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double dout = 0.f;
	N = Nx*Ny*Nz;

	// Loop over the boundary - threadblocks delineated by start...finish
	double sum = 0.f;
    char id;
	for (n=outlet; n<N-Nx*Ny; n++){
        id = ID[n];
        if (id>0){
            //........................................................................
            // Read distributions from "opposite" memory convention
            //........................................................................
            f2 = distodd[n];
            f4 = distodd[N+n];
            //f6 = distodd[2*N+n];
            f8 = distodd[3*N+n];
            f10 = distodd[4*N+n];
            //f12 = distodd[5*N+n];
            f14 = distodd[6*N+n];
            //f16 = distodd[7*N+n];
            f18 = distodd[8*N+n];
            //........................................................................
            f0 = disteven[n];
            f1 = disteven[N+n];
            f3 = disteven[2*N+n];
            f5 = disteven[3*N+n];
            f7 = disteven[4*N+n];
            f9 = disteven[5*N+n];
            f11 = disteven[6*N+n];
            //f13 = disteven[7*N+n];
            f15 = disteven[8*N+n];
            //f17 = disteven[9*N+n];

            sum += (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));
        }
	}
	dout = sum/area;
	return dout;
}

extern "C" void ScaLBL_D3Q19_Pressure_BC_z(double *disteven, double *distodd, double din,
								  int Nx, int Ny, int Nz)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz;
	ux = uy = 0.0;

	N = Nx*Ny*Nz;

	double Cxz,Cyz;

	for (n=Nx*Ny; n<2*Nx*Ny; n++){

		//........................................................................
		// Read distributions from "opposite" memory convention
		//........................................................................
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
		//...................................................
		//........Determine the inlet flow velocity.........
		//			uz = -1 + (f0+f3+f4+f1+f2+f7+f8+f10+f9
		//					   + 2*(f5+f15+f18+f11+f14))/din;
		//........Set the unknown distributions..............
		//			f6 = f5 - 0.3333333333333333*din*uz;
		//			f16 = f15 - 0.1666666666666667*din*uz;
		//			f17 = f16 - f3 + f4-f15+f18-f7+f8-f10+f9;
		//			f12= 0.5*(-din*uz+f5+f15+f18+f11+f14-f6-f16-
		//					  f17+f1-f2-f14+f11+f7-f8-f10+f9);
		//			f13= -din*uz+f5+f15+f18+f11+f14-f6-f16-f17-f12;
		// Determine the inlet flow velocity
		//ux = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14);
		//uy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18);
		uz = din - (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f5 = f6 + 0.3333333333333333*uz;
		f11 = f12 + 0.1666666666666667*(uz+ux)-Cxz;
		f14 = f13 + 0.1666666666666667*(uz-ux)+Cxz;
		f15 = f16 + 0.1666666666666667*(uy+uz)-Cyz;
		f18 = f17 + 0.1666666666666667*(uz-uy)+Cyz;
		//........Store in "opposite" memory location..........
		/*		distodd[2*N+n] = f5;
		distodd[5*N+n] = f11;
		disteven[7*N+n] = f14;
		distodd[7*N+n] = f15;
		disteven[9*N+n] = f18;

		*/

		disteven[3*N+n] = f5;
		disteven[6*N+n] = f11;
		distodd[6*N+n] = f14;
		disteven[8*N+n] = f15;
		distodd[8*N+n] = f18;
		/*
		printf("Site=%i\n",n);
		printf("ux=%f, uy=%f, uz=%f\n",ux,uy,uz);
		printf("Cxz=%f, Cyz=%f\n",Cxz,Cyz);
		n = N;
		*/
		//...................................................
	}
}

extern "C" void ScaLBL_D3Q19_Pressure_BC_Z(double *disteven, double *distodd, double dout,
								   int Nx, int Ny, int Nz, int outlet)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz;
	ux = uy = 0.0;

	double Cxz,Cyz;
	N = Nx*Ny*Nz;

	// Loop over the boundary - threadblocks delineated by start...finish
	for (n=outlet; n<N-Nx*Ny; n++){

		//........................................................................
		// Read distributions from "opposite" memory convention
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
		// Determine the outlet flow velocity
		//ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
		//uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		//uz = -1.0 + (f0+f4+f3+f2+f1+f8+f7+f9+f10 + 2*(f6+f16+f17+f12+f13))/dout;

		// Determine the inlet flow velocity
		//ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
		//uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		uz = -dout + (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f6 = f5 - 0.3333333333333333*uz;
		f12 = f11 - 0.1666666666666667*(uz+ux)+Cxz;
		f13 = f14 - 0.1666666666666667*(uz-ux)-Cxz;
		f16 = f15 - 0.1666666666666667*(uy+uz)+Cyz;
		f17 = f18 - 0.1666666666666667*(uz-uy)-Cyz;

		//........Store in "opposite" memory location..........
		distodd[2*N+n] = f6;
		distodd[5*N+n] = f12;
		disteven[7*N+n] = f13;
		distodd[7*N+n] = f16;
		disteven[9*N+n] = f17;
		//...................................................

	}
}

extern "C" void ScaLBL_D3Q19_Velocity_BC_z(double *disteven, double *distodd, double uz,
								   int Nx, int Ny, int Nz)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double din;
    double Cxz,Cyz;
	N = Nx*Ny*Nz;

	for (n=Nx*Ny; n<2*Nx*Ny; n++){

		//........................................................................
		// Read distributions from "opposite" memory convention
		//........................................................................
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
		//...................................................

		// Determine 'din' based on the inlet flow velocity
        // NOTE: Default: ux = uy = 0.0, we only specify 'uz'
		//ux = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14);
		//uy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18);
		din = (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17))/(1.0-uz);

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8); // ux = 0.0
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8); // uy = 0.0

		f5 = f6 + 0.3333333333333333*din*uz;
		f11 = f12 + 0.1666666666666667*din*uz-Cxz;
		f14 = f13 + 0.1666666666666667*din*uz+Cxz;
		f15 = f16 + 0.1666666666666667*din*uz-Cyz;
		f18 = f17 + 0.1666666666666667*din*uz+Cyz;


		//........Store in "opposite" memory location..........
		disteven[3*N+n] = f5;
		disteven[6*N+n] = f11;
		distodd[6*N+n] = f14;
		disteven[8*N+n] = f15;
		distodd[8*N+n] = f18;
		//...................................................
	}
}

extern "C" void ScaLBL_D3Q19_Velocity_BC_Z(double *disteven, double *distodd, double uz,
								   int Nx, int Ny, int Nz, int outlet)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double dout;
    double Cxz,Cyz;

	N = Nx*Ny*Nz;

	// Loop over the boundary - threadblocks delineated by start...finish
	for (n=outlet; n<N-Nx*Ny; n++){
		//........................................................................
		// Read distributions from "opposite" memory convention
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

		// Determine the 'dout' based on the outlet flow velocity
        // Default: ux = uy = 0.0, we only specify 'uz'
//		ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
//		uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		dout = (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18))/(1.0+uz);

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8); // ux = 0.0
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8); // uy = 0.0

		f6 = f5 - 0.3333333333333333*dout*uz;
		f12 = f11 - 0.16666666666666678*dout*uz+Cxz;
		f13 = f14 - 0.16666666666666678*dout*uz-Cxz;
		f16 = f15 - 0.16666666666666678*dout*uz+Cyz;
		f17 = f18 - 0.16666666666666678*dout*uz-Cyz;

		//........Store in "opposite" memory location..........
		distodd[2*N+n] = f6;
		distodd[5*N+n] = f12;
		disteven[7*N+n] = f13;
		distodd[7*N+n] = f16;
		disteven[9*N+n] = f17;
		//...................................................
	}
}

extern "C" void ScaLBL_D3Q19_Velocity(char *ID, double *disteven, double *distodd, double *vel, int Nx, int Ny, int Nz)
{
	int n,N;
	// distributions
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double vx,vy,vz;

	N = Nx*Ny*Nz;

	for (n=0; n<N; n++){
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
		else{
			for(int q=0; q<9; q++){
				disteven[q*N+n] = -1.0;
				distodd[q*N+n] = -1.0;
			}
			disteven[9*N+n] = -1.0;
            
            //For ID[n]<0 - solid nodes
			vel[n] = 0.0;
			vel[N+n] = 0.0;
			vel[2*N+n] = 0.0;
		}
	}
}

extern "C" void ScaLBL_D3Q19_Pressure(const char *ID, const double *disteven, const double *distodd,
    double *Pressure, int Nx, int Ny, int Nz)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	N = Nx*Ny*Nz;

	for (n=0; n<N; n++){

		if (ID[n] > 0){
			//........................................................................
			// Registers to store the distributions
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
			Pressure[n] = 0.3333333333333333*(f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+
					f9+f12+f11+f14+f13+f16+f15+f18+f17);
		}
	}
}

