/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
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

extern "C" void ScaLBL_D3Q19_AA_Init(double *f_even, double *f_odd, int Np)
{
	int n;
	for (n=0; n<Np; n++){
		f_even[n] = 0.3333333333333333;
		f_odd[n] = 0.055555555555555555;		//double(100*n)+1.f;
		f_even[Np+n] = 0.055555555555555555;	//double(100*n)+2.f;
		f_odd[Np+n] = 0.055555555555555555;	//double(100*n)+3.f;
		f_even[2*Np+n] = 0.055555555555555555;	//double(100*n)+4.f;
		f_odd[2*Np+n] = 0.055555555555555555;	//double(100*n)+5.f;
		f_even[3*Np+n] = 0.055555555555555555;	//double(100*n)+6.f;
		f_odd[3*Np+n] = 0.0277777777777778;   //double(100*n)+7.f;
		f_even[4*Np+n] = 0.0277777777777778;   //double(100*n)+8.f;
		f_odd[4*Np+n] = 0.0277777777777778;   //double(100*n)+9.f;
		f_even[5*Np+n] = 0.0277777777777778;  //double(100*n)+10.f;
		f_odd[5*Np+n] = 0.0277777777777778;  //double(100*n)+11.f;
		f_even[6*Np+n] = 0.0277777777777778;  //double(100*n)+12.f;
		f_odd[6*Np+n] = 0.0277777777777778;  //double(100*n)+13.f;
		f_even[7*Np+n] = 0.0277777777777778;  //double(100*n)+14.f;
		f_odd[7*Np+n] = 0.0277777777777778;  //double(100*n)+15.f;
		f_even[8*Np+n] = 0.0277777777777778;  //double(100*n)+16.f;
		f_odd[8*Np+n] = 0.0277777777777778;  //double(100*n)+17.f;
		f_even[9*Np+n] = 0.0277777777777778;  //double(100*n)+18.f;

	}
}

extern "C" void ScaLBL_D3Q19_Init(double *dist, int Np)
{
	int n;
	for (n=0; n<Np; n++){
		dist[n] = 0.3333333333333333;
		dist[Np+n] = 0.055555555555555555;		//double(100*n)+1.f;
		dist[2*Np+n] = 0.055555555555555555;	//double(100*n)+2.f;
		dist[3*Np+n] = 0.055555555555555555;	//double(100*n)+3.f;
		dist[4*Np+n] = 0.055555555555555555;	//double(100*n)+4.f;
		dist[5*Np+n] = 0.055555555555555555;	//double(100*n)+5.f;
		dist[6*Np+n] = 0.055555555555555555;	//double(100*n)+6.f;
		dist[7*Np+n] = 0.0277777777777778;   //double(100*n)+7.f;
		dist[8*Np+n] = 0.0277777777777778;   //double(100*n)+8.f;
		dist[9*Np+n] = 0.0277777777777778;   //double(100*n)+9.f;
		dist[10*Np+n] = 0.0277777777777778;  //double(100*n)+10.f;
		dist[11*Np+n] = 0.0277777777777778;  //double(100*n)+11.f;
		dist[12*Np+n] = 0.0277777777777778;  //double(100*n)+12.f;
		dist[13*Np+n] = 0.0277777777777778;  //double(100*n)+13.f;
		dist[14*Np+n] = 0.0277777777777778;  //double(100*n)+14.f;
		dist[15*Np+n] = 0.0277777777777778;  //double(100*n)+15.f;
		dist[16*Np+n] = 0.0277777777777778;  //double(100*n)+16.f;
		dist[17*Np+n] = 0.0277777777777778;  //double(100*n)+17.f;
		dist[18*Np+n] = 0.0277777777777778;  //double(100*n)+18.f;
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


extern "C" double ScaLBL_D3Q19_Flux_BC_z(double *disteven, double *distodd, double flux,
		int Nx, int Ny, int Nz){
	// Note that this routine assumes the distributions are stored "opposite"
	// odd distributions in disteven and even distributions in distodd.
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double din = 0.f;
	N = Nx*Ny*Nz;

	double A = 1.f*double(Nx*Ny);
	double sum = 0.f;
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

		// Determine the outlet flow velocity
		//sum += 1.0 - (f0+f4+f3+f2+f1+f8+f7+f9+ f10 + 2*(f5+ f15+f18+f11+f14))/din;
		//sum += (f0+f4+f3+f2+f1+f8+f7+f9+ f10 + 2*(f5+f15+f18+f11+f14));
		sum += (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));
	}
	din = sum/(A*(1.0-flux));
	return din;
}

extern "C" double ScaLBL_D3Q19_AAodd_Flux_BC_z(int *d_neighborList, int *list, double *dist, double flux, 
		double area, int count, int Np){
	int idx, n;
	int nread;

	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double factor = 1.f/(area);
	double sum = 0.f;
	
	for (idx=0; idx<count; idx++){
		n = list[idx];

		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+6*Np];
		f7 = dist[nread];

		nread = d_neighborList[n+8*Np];
		f9 = dist[nread];

		nread = d_neighborList[n+12*Np];
		f13 = dist[nread];

		nread = d_neighborList[n+16*Np];
		f17 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

		nread = d_neighborList[n+7*Np];
		f8 = dist[nread];

		nread = d_neighborList[n+9*Np];
		f10 = dist[nread];

		nread = d_neighborList[n+11*Np];
		f12 = dist[nread];

		nread = d_neighborList[n+15*Np];
		f16 = dist[nread];

		sum += factor*(f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));
	}
	
	return sum;
}

extern "C" double ScaLBL_D3Q19_AAeven_Flux_BC_z(int *list, double *dist, double flux, double area, 
		 int count, int Np){
	int idx, n;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double factor = 1.f/(area);
	double sum = 0.f;
	
	for (idx=0; idx<count; idx++){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
		f7 = dist[8*Np+n];
		f8 = dist[7*Np+n];
		f9 = dist[10*Np+n];
		f10 = dist[9*Np+n];
		f12 = dist[11*Np+n];
		f13 = dist[14*Np+n];
		f16 = dist[15*Np+n];
		f17 = dist[18*Np+n];
		sum += factor*(f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));
	}
	return sum;
}


extern "C" double ScaLBL_D3Q19_Flux_BC_Z(double *disteven, double *distodd, double flux,
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
	double A = 1.f*double(Nx*Ny);
	double sum = 0.f;
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

		sum += (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));

	}
	dout = sum/(A*(1.0+flux));
	return dout;
}

extern "C" void ScaLBL_D3Q19_AAeven_Pressure_BC_z(int *list, double *dist, double din, int count, int Np)
{
	int idx, n;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz,Cyz,Cxz;
	ux = uy = 0.0;
	for (int idx=0; idx<count; idx++){
		n = list[idx];
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f6 = dist[5*Np+n];
		f7 = dist[8*Np+n];
		f8 = dist[7*Np+n];
		f9 = dist[10*Np+n];
		f10 = dist[9*Np+n];
		f12 = dist[11*Np+n];
		f13 = dist[14*Np+n];
		f16 = dist[15*Np+n];
		f17 = dist[18*Np+n];
		//...................................................
		// Determine the inlet flow velocity
		//ux = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14);
		//uy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18);
		uz = din - (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f5 = f6 + 0.33333333333333338*uz;
		f11 = f12 + 0.16666666666666678*(uz+ux)-Cxz;
		f14 = f13 + 0.16666666666666678*(uz-ux)+Cxz;
		f15 = f16 + 0.16666666666666678*(uy+uz)-Cyz;
		f18 = f17 + 0.16666666666666678*(uz-uy)+Cyz;

		dist[6*Np+n] = f5;
		dist[12*Np+n] = f11;
		dist[13*Np+n] = f14;
		dist[16*Np+n] = f15;
		dist[17*Np+n] = f18;
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_Pressure_BC_Z(int *list, double *dist, double dout, int count, int Np)
{
	int idx, n;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz,Cyz,Cxz;
	ux = uy = 0.0;
	for (int idx=0; idx<count; idx++){
		n = list[idx];
		//........................................................................
		// Read distributions 
		//........................................................................
		f0 = dist[n];
		f1 = dist[2*Np+n];
		f2 = dist[1*Np+n];
		f3 = dist[4*Np+n];
		f4 = dist[3*Np+n];
		f5 = dist[6*Np+n];
		f7 = dist[8*Np+n];
		f8 = dist[7*Np+n];
		f9 = dist[10*Np+n];
		f10 = dist[9*Np+n];
		f11 = dist[12*Np+n];
		f14 = dist[13*Np+n];
		f15 = dist[16*Np+n];
		f18 = dist[17*Np+n];
		
		// Determine the outlet flow velocity
		//ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
		//uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		uz = -dout + (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f6 = f5 - 0.33333333333333338*uz;
		f12 = f11 - 0.16666666666666678*(uz+ux)+Cxz;
		f13 = f14 - 0.16666666666666678*(uz-ux)-Cxz;
		f16 = f15 - 0.16666666666666678*(uy+uz)+Cyz;
		f17 = f18 - 0.16666666666666678*(uz-uy)-Cyz;

		dist[5*Np+n] = f6;
		dist[11*Np+n] = f12;
		dist[14*Np+n] = f13;
		dist[15*Np+n] = f16;
		dist[18*Np+n] = f17;
		//...................................................
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Pressure_BC_z(int *d_neighborList, int *list, double *dist, double din, int count, int Np)
{
	int idx, n;
	int nread;
	int nr5,nr11,nr14,nr15,nr18;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz,Cyz,Cxz;
	ux = uy = 0.0;

	for (int idx=0; idx<count; idx++){
		n = list[idx];
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+6*Np];
		f7 = dist[nread];

		nread = d_neighborList[n+8*Np];
		f9 = dist[nread];

		nread = d_neighborList[n+12*Np];
		f13 = dist[nread];

		nread = d_neighborList[n+16*Np];
		f17 = dist[nread];

		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+5*Np];
		f6 = dist[nread];

		nread = d_neighborList[n+7*Np];
		f8 = dist[nread];

		nread = d_neighborList[n+9*Np];
		f10 = dist[nread];

		nread = d_neighborList[n+11*Np];
		f12 = dist[nread];

		nread = d_neighborList[n+15*Np];
		f16 = dist[nread];

		// Unknown distributions
		nr5 = d_neighborList[n+4*Np];
		nr11 = d_neighborList[n+10*Np];
		nr15 = d_neighborList[n+14*Np];
		nr14 = d_neighborList[n+13*Np];
		nr18 = d_neighborList[n+17*Np];
		
		//...................................................
		// Determine the inlet flow velocity
		//ux = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14);
		//uy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18);
		uz = din - (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f5 = f6 + 0.33333333333333338*uz;
		f11 = f12 + 0.16666666666666678*(uz+ux)-Cxz;
		f14 = f13 + 0.16666666666666678*(uz-ux)+Cxz;
		f15 = f16 + 0.16666666666666678*(uy+uz)-Cyz;
		f18 = f17 + 0.16666666666666678*(uz-uy)+Cyz;

		dist[nr5] = f5;
		dist[nr11] = f11;
		dist[nr14] = f14;
		dist[nr15] = f15;
		dist[nr18] = f18;
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Pressure_BC_Z(int *d_neighborList, int *list, double *dist, double dout, int count, int Np)
{
	int idx,n,nread;
	int nr6,nr12,nr13,nr16,nr17;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz,Cyz,Cxz;
	ux = uy = 0.0;

	for (int idx=0; idx<count; idx++){
		n = list[idx];
		//........................................................................
		// Read distributions 
		//........................................................................
		f0 = dist[n];

		nread = d_neighborList[n];
		f1 = dist[nread];

		nread = d_neighborList[n+2*Np];
		f3 = dist[nread];

		nread = d_neighborList[n+4*Np];
		f5 = dist[nread];

		nread = d_neighborList[n+6*Np];
		f7 = dist[nread];

		nread = d_neighborList[n+8*Np];
		f9 = dist[nread];

		nread = d_neighborList[n+10*Np];
		f11 = dist[nread];

		nread = d_neighborList[n+14*Np];
		f15 = dist[nread];


		nread = d_neighborList[n+Np];
		f2 = dist[nread];

		nread = d_neighborList[n+3*Np];
		f4 = dist[nread];

		nread = d_neighborList[n+7*Np];
		f8 = dist[nread];

		nread = d_neighborList[n+9*Np];
		f10 = dist[nread];

		nread = d_neighborList[n+13*Np];
		f14 = dist[nread];

		nread = d_neighborList[n+17*Np];
		f18 = dist[nread];
		
		// unknown distributions
		nr6 = d_neighborList[n+5*Np];
		nr12 = d_neighborList[n+11*Np];
		nr16 = d_neighborList[n+15*Np];
		nr17 = d_neighborList[n+16*Np];
		nr13 = d_neighborList[n+12*Np];

		// Determine the inlet flow velocity
		//ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
		//uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		uz = -dout + (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f6 = f5 - 0.33333333333333338*uz;
		f12 = f11 - 0.16666666666666678*(uz+ux)+Cxz;
		f13 = f14 - 0.16666666666666678*(uz-ux)-Cxz;
		f16 = f15 - 0.16666666666666678*(uy+uz)+Cyz;
		f17 = f18 - 0.16666666666666678*(uz-uy)-Cyz;

		//........Store in "opposite" memory location..........
		dist[nr6] = f6;
		dist[nr12] = f12;
		dist[nr13] = f13;
		dist[nr16] = f16;
		dist[nr17] = f17;
		//...................................................
	}
}

extern "C" void ScaLBL_D3Q19_Pressure_BC_z(int *list, double *dist, double din, int count, int Np)
{
	int n;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz;
	double Cxz,Cyz;

	for (int idx=0; idx<count; idx++){
		n = list[idx];
		//........................................................................
		// Read distributions from "opposite" memory convention
		//........................................................................
		//........................................................................
		f0 = dist[n];
		f1 = dist[Np+n];
		f2 = dist[2*Np+n];
		f3 = dist[3*Np+n];
		f4 = dist[4*Np+n];
		f6 = dist[6*Np+n];
		f7 = dist[7*Np+n];
		f8 = dist[8*Np+n];
		f9 = dist[9*Np+n];
		f10 = dist[10*Np+n];
		f12 = dist[12*Np+n];
		f13 = dist[13*Np+n];
		f16 = dist[16*Np+n];
		f17 = dist[17*Np+n];
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
		ux = (f1-f2+f7-f8+f9-f10+f11-f12+f13-f14);
		uy = (f3-f4+f7-f8-f9+f10+f15-f16+f17-f18);
		uz = din - (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f6+f12+f13+f16+f17));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f5 = f6 + 0.33333333333333338*uz;
		f11 = f12 + 0.16666666666666678*(uz+ux)-Cxz;
		f14 = f13 + 0.16666666666666678*(uz-ux)+Cxz;
		f15 = f16 + 0.16666666666666678*(uy+uz)-Cyz;
		f18 = f17 + 0.16666666666666678*(uz-uy)+Cyz;
		//........Store in "opposite" memory location..........
		dist[5*Np+n] = f5;
		dist[11*Np+n] = f11;
		dist[14*Np+n] = f14;
		dist[15*Np+n] = f15;
		dist[18*Np+n] = f18;
		
		/*
		printf("Site=%i\n",n);
		printf("ux=%f, uy=%f, uz=%f\n",ux,uy,uz);
		printf("Cxz=%f, Cyz=%f\n",Cxz,Cyz);
		n = N;
		 */
		//...................................................
	}
}

extern "C" void ScaLBL_D3Q19_Pressure_BC_Z(int *list, double *dist, double dout, int count, int Np)
{
	int n;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double ux,uy,uz;
	double Cxz,Cyz;

	for (int idx=0; idx<count; idx++){
		n = list[idx];
		
		//........................................................................
		// Read distributions 
		//........................................................................
		f0 = dist[n];
		f1 = dist[Np+n];
		f2 = dist[2*Np+n];
		f3 = dist[3*Np+n];
		f4 = dist[4*Np+n];
		f5 = dist[5*Np+n];
		f7 = dist[7*Np+n];
		f8 = dist[8*Np+n];
		f9 = dist[9*Np+n];
		f10 = dist[10*Np+n];
		f11 = dist[11*Np+n];
		f14 = dist[14*Np+n];
		f15 = dist[15*Np+n];
		f18 = dist[18*Np+n];
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
		ux = f1-f2+f7-f8+f9-f10+f11-f12+f13-f14;
		uy = f3-f4+f7-f8-f9+f10+f15-f16+f17-f18;
		uz = -dout + (f0+f1+f2+f3+f4+f7+f8+f9+f10 + 2*(f5+f11+f14+f15+f18));

		Cxz = 0.5*(f1+f7+f9-f2-f10-f8) - 0.3333333333333333*ux;
		Cyz = 0.5*(f3+f7+f10-f4-f9-f8) - 0.3333333333333333*uy;

		f6 = f5 - 0.33333333333333338*uz;
		f12 = f11 - 0.16666666666666678*(uz+ux)+Cxz;
		f13 = f14 - 0.16666666666666678*(uz-ux)-Cxz;
		f16 = f15 - 0.16666666666666678*(uy+uz)+Cyz;
		f17 = f18 - 0.16666666666666678*(uz-uy)-Cyz;

		//........Store in "opposite" memory location..........
		dist[6*Np+n] = f6;
		dist[12*Np+n] = f12;
		dist[13*Np+n] = f13;
		dist[16*Np+n] = f16;
		dist[17*Np+n] = f17;
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

		// Determine the outlet flow velocity
		//	uz = 1.0 - (f0+f4+f3+f2+f1+f8+f7+f9+f10 +
		//			2*(f5+f15+f18+f11+f14))/din;
		din = (f0+f4+f3+f2+f1+f8+f7+f9+f10+2*(f5+f15+f18+f11+f14))/(1.0-uz);
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

extern "C" void ScaLBL_D3Q19_Velocity_BC_Z(double *disteven, double *distodd, double uz,
		int Nx, int Ny, int Nz, int outlet)
{
	int n,N;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double dout;

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
		//uz = -1.0 + (f0+f4+f3+f2+f1+f8+f7+f9+f10 + 2*(f6+f16+f17+f12+f13))/dout;
		dout = (f0+f4+f3+f2+f1+f8+f7+f9+f10 + 2*(f6+f16+f17+f12+f13))/(1.0+uz);
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

extern "C" void ScaLBL_D3Q19_Momentum(double *dist, double *vel, int Np)
{
	int n;
	int N =Np;
	// distributions
	double f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	double vx,vy,vz;

	for (n=0; n<N; n++){
		//........................................................................
		// Registers to store the distributions
		//........................................................................
		f2 = dist[2*N+n];
		f4 = dist[4*N+n];
		f6 = dist[6*N+n];
		f8 = dist[8*N+n];
		f10 = dist[10*N+n];
		f12 = dist[12*N+n];
		f14 = dist[14*N+n];
		f16 = dist[16*N+n];
		f18 = dist[18*N+n];
		//........................................................................
		f1 = dist[N+n];
		f3 = dist[3*N+n];
		f5 = dist[5*N+n];
		f7 = dist[7*N+n];
		f9 = dist[9*N+n];
		f11 = dist[11*N+n];
		f13 = dist[13*N+n];
		f15 = dist[15*N+n];
		f17 = dist[17*N+n];
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

extern "C" void ScaLBL_D3Q19_Pressure(double *dist, double *Pressure, int N)
{
	int n;
	// distributions
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	for (n=0; n<N; n++){
		//........................................................................
		// Registers to store the distributions
		//........................................................................
		f0 = dist[n];
		f2 = dist[2*N+n];
		f4 = dist[4*N+n];
		f6 = dist[6*N+n];
		f8 = dist[8*N+n];
		f10 = dist[10*N+n];
		f12 = dist[12*N+n];
		f14 = dist[14*N+n];
		f16 = dist[16*N+n];
		f18 = dist[18*N+n];
		//........................................................................
		f1 = dist[N+n];
		f3 = dist[3*N+n];
		f5 = dist[5*N+n];
		f7 = dist[7*N+n];
		f9 = dist[9*N+n];
		f11 = dist[11*N+n];
		f13 = dist[13*N+n];
		f15 = dist[15*N+n];
		f17 = dist[17*N+n];
		//.................Compute the velocity...................................
		Pressure[n] = 0.3333333333333333*(f0+f2+f1+f4+f3+f6+f5+f8+f7+f10+
				f9+f12+f11+f14+f13+f16+f15+f18+f17);
	}
}

extern "C" void ScaLBL_D3Q19_AAeven_MRT(double *dist, int start, int finish, int Np, double rlx_setA, double rlx_setB, double Fx,
		double Fy, double Fz){
	int n;
	double fq,fp;
	// conserved momemnts
	double rho,jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;

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


	for (int n=start; n<finish; n++){
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



		//........................................................................
		//					READ THE DISTRIBUTIONS
		//		(read from opposite array due to previous swap operation)
		//........................................................................

		//..............incorporate external force................................................
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

extern "C" void ScaLBL_D3Q19_AAodd_MRT(int *neighborList, double *dist, int start, int finish, int Np, double rlx_setA, double rlx_setB, double Fx,
		double Fy, double Fz){
	int n;
	double fq,fp;
	// conserved momemnts
	double rho,jx,jy,jz;
	// non-conserved moments
	double m1,m2,m4,m6,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
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


	int nread;
	for (int n=start; n<finish; n++){
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

		//..............incorporate external force................................................
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

extern "C" void ScaLBL_D3Q19_AAeven_Compact(char * ID, double *dist,  int Np) {

	int n;
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;

	for (int n=0; n<Np; n++){

		//........................................................................
		//                                      READ THE DISTRIBUTIONS
		//              (read from opposite array due to previous swap operation)
		//........................................................................
		// even
		f2 = dist[10*Np+n];
		f4 = dist[11*Np+n];
		f6 = dist[12*Np+n];
		f8 = dist[13*Np+n];
		f10 = dist[14*Np+n];
		f12 = dist[15*Np+n];
		f14 = dist[16*Np+n];
		f16 = dist[17*Np+n];
		f18 = dist[18*Np+n];

		f0 = dist[n];
		// odd
		f1 = dist[Np+n];
		f3 = dist[2*Np+n];
		f5 = dist[3*Np+n];
		f7 = dist[4*Np+n];
		f9 = dist[5*Np+n];
		f11 = dist[6*Np+n];
		f13 = dist[7*Np+n];
		f15 = dist[8*Np+n];
		f17 = dist[9*Np+n];

		//........................................................................
		//                                      WRITE THE DISTRIBUTIONS
		// even
		//disteven[n] = f0;
		dist[Np+n] = f2;
		dist[2*Np+n] = f4;
		dist[3*Np+n] = f6;
		dist[4*Np+n] = f8;
		dist[5*Np+n] = f10;
		dist[6*Np+n] = f12;
		dist[7*Np+n] = f14;
		dist[8*Np+n] = f16;
		dist[9*Np+n] = f18;

		// odd
		dist[10*Np+n] = f1;
		dist[11*Np+n] = f3;
		dist[12*Np+n] = f5;
		dist[13*Np+n] = f7;
		dist[14*Np+n] = f9;
		dist[15*Np+n] = f11;
		dist[16*Np+n] = f13;
		dist[17*Np+n] = f15;
		dist[18*Np+n] = f17;
		//........................................................................
	}
}

extern "C" void ScaLBL_D3Q19_AAodd_Compact(char * ID, int *neighborList, double *dist, int Np) {
	int n;
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
	double f10,f11,f12,f13,f14,f15,f16,f17,f18;
	int nread;

	for (int n=0; n<Np; n++){
		//........Get 1-D index for this thread....................

		f0 = dist[n];

		nread = neighborList[n]; // + 0*Np
		f2 = dist[nread];

		nread = neighborList[n+2*Np];
		f4 = dist[nread];

		nread = neighborList[n+4*Np];
		f6 = dist[nread];

		nread = neighborList[n+6*Np];
		f8 = dist[nread];

		nread = neighborList[n+8*Np];
		f10 = dist[nread];

		nread = neighborList[n+10*Np];
		f12 = dist[nread];

		nread = neighborList[n+12*Np];
		f14 = dist[nread];

		nread = neighborList[n+14*Np];
		f16 = dist[nread];

		nread = neighborList[n+16*Np];
		f18 = dist[nread];


		nread = neighborList[n+Np];
		f1 = dist[nread];

		nread = neighborList[n+3*Np];
		f3 = dist[nread];

		nread = neighborList[n+5*Np];
		f5 = dist[nread];

		nread = neighborList[n+7*Np];
		f7 = dist[nread];

		nread = neighborList[n+9*Np];
		f9 = dist[nread];

		nread = neighborList[n+11*Np];
		f11 = dist[nread];

		nread = neighborList[n+13*Np];
		f13 = dist[nread];

		nread = neighborList[n+15*Np];
		f15 = dist[nread];

		nread = neighborList[n+17*Np];
		f17 = dist[nread];


		nread = neighborList[n];
		dist[nread] = f1;

		nread = neighborList[n+2*Np];
		dist[nread] = f3;

		nread = neighborList[n+4*Np];
		dist[nread] = f5;

		nread = neighborList[n+6*Np];
		dist[nread] = f7;

		nread = neighborList[n+8*Np];
		dist[nread] = f9;

		nread = neighborList[n+10*Np];
		dist[nread] = f11;

		nread = neighborList[n+12*Np];
		dist[nread] = f13;

		nread = neighborList[n+14*Np];
		dist[nread] = f15;

		nread = neighborList[n+16*Np];
		dist[nread] = f17;


		nread = neighborList[n+Np];
		dist[nread] = f2;

		nread = neighborList[n+3*Np];
		dist[nread] = f4;

		nread = neighborList[n+5*Np];
		dist[nread] = f6;

		nread = neighborList[n+7*Np];
		dist[nread] = f8;

		nread = neighborList[n+9*Np];
		dist[nread] = f10;

		nread = neighborList[n+11*Np];
		dist[nread] = f12;

		nread = neighborList[n+13*Np];
		dist[nread] = f14;

		nread = neighborList[n+15*Np];
		dist[nread] = f16;

		nread = neighborList[n+17*Np];
		dist[nread] = f18;
	}
}
